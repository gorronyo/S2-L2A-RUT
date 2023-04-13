#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sept 20 14:03:11 2016

@author: Jan
"""

import os, math
import datetime
import numpy as np
import snappy
import glob

from snappy import GPF
from snappy import HashMap as hash

# from osgeo import gdal
import xml.etree.ElementTree as ET
import zipfile

import matplotlib

matplotlib.use('Agg')  # this does not show the plot on the screen

S2_MSI_L1C_TYPE_STRING = 'S2_MSI_Level-1C'
S2_MSI_L2A_TYPE_STRING = ['S2_MSI_Level-2Ap', 'S2_MSI_Level-2A']
S2_BAND_NAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
S2_BAND_SAMPLING = {'B1': 60, 'B2': 10, 'B3': 10, 'B4': 10, 'B5': 20, 'B6': 20, 'B7': 20, 'B8': 10, 'B8A': 20, 'B9': 60,
                    'B10': 60, 'B11': 20, 'B12': 20}  # product before October 2016 is not automatic
S2_L2AQUALITY_SAMPLING = {'quality_aot': 10, 'quality_wvp': 10, 'quality_cloud_confidence': 20,
                          'quality_snow_confidence': 20, 'quality_scene_classification': 20}


class S2Processor:
    def __init__(self):
        # INPUTS -------------------------------------------------------------------------------------------------------
        self.prod_path = None  # path to S2 L1C or S2 L2A product
        self.selected_bands = None  # array with the tags for the bands to be processed
        # http://step.esa.int/docs/v5.0/apidoc/engine/org/esa/snap/core/datamodel/GeoPos.html
        self.lat_centre = None  # ROI centre. The geographical latitude in decimal degree, valid range is -90 to +90.
        self.lon_centre = None  # ROI centre. The geographical longitude in decimal degree, valid range is -180 to +180.

        self.sensing_time = None  # Timestamp at the S2 L1C acquisition

        self.w = None  # width of the desired data array (meters)
        self.h = None  # height of the desired data array (meters)

        # OUTPUTS ------------------------------------------------------------------------------------------------------
        self.L1C_ref = []  # list of arrays with ROI pixels per band. They provide TOA reflectance values as float.
        self.L1C_rad = []
        self.view_azimuth = []  # same for each angle
        self.view_zenith = []
        self.sun_azimuth = []
        self.sun_zenith = []
        self.cloud = None
        self.cirrus = None

        self.atmosdict = {}  # this dictionary contains the ECMWF data for the data ROI
        self.metadatadict = {'spacecraft': None, 'quant': None, 'A': [], 'alpha': [],
                             'beta': [], 'Esun': []}  # this dictionary contains the relevant metadata parameters

        # Possible scene classification values
        # 0 : SC_NODATA, 1 : SC_SATURATED_DEFECTIVE, 2 : SC_DARK_FEATURE_SHADOW, 3 : SC_CLOUD_SHADOW, 4 : SC_VEGETATION
        # 5 : SC_NOT_VEGETATED, 6 : SC_WATER, 7 : SC_UNCLASSIFIED, 8 : SC_CLOUD_MEDIUM_PROBA, 9 : SC_CLOUD_HIGH_PROBA
        # 10 : SC_THIN_CIRRUS, 11 : SC_SNOW_ICE
        # atm_season can be "h" for midlatitude summer and "w" for midlatitude winter.
        # aot_type values are: rura, mari, urba, dese
        self.L2A_ancillary = {'quality_aot': None, 'quality_wvp': None, 'quality_cloud_confidence': None,
                              'quality_snow_confidence': None, 'quality_scene_classification': None, 'ozone': None,
                              'ozone_source': None, 'aot_method': None, 'altitude': None, 'atm_season': None,
                              'aot_type': None}

    def get_roidata(self):
        '''
        Executes a set a routines for each of the selected bands
        '''
        # select the specific image and masks for the band
        source_product = snappy.ProductIO.readProduct(self.prod_path)  # S2 L1C or L2A product
        # checks that the product is the correct one. There are two types of L2A tag names and one of L1C.
        if source_product.getProductType() == S2_MSI_L1C_TYPE_STRING:
            product_meta, datastrip_meta, granules_meta = self.parse_product(source_product, 'Level-1C')
        elif any(source_product.getProductType() == t for t in S2_MSI_L2A_TYPE_STRING):
            product_meta, datastrip_meta, granules_meta = self.parse_product(source_product, 'Level-2A')
            # S2 L2A routines for each band. WE ONLY NEED BASIC AOT, WV, AND SCENE!!!
            for band in ['quality_aot', 'quality_wvp', 'quality_cloud_confidence', 'quality_snow_confidence',
                         'quality_scene_classification']:
                source_band = source_product.getBand(band)
                wpix = int(round(self.w / S2_L2AQUALITY_SAMPLING[band]))
                hpix = int(round(self.h / S2_L2AQUALITY_SAMPLING[band]))
                if any([wpix, hpix]) == 0:  # minimum pixel size is 1
                    wpix = 1
                    hpix = 1
                x_off, y_off = self.get_corners(wpix, hpix, source_band)
                roi_data = np.zeros(wpix * hpix, np.float32)
                source_band.readPixels(x_off, y_off, wpix, hpix, roi_data)
                roi_data.shape = wpix, hpix
                self.L2A_ancillary[band] = roi_data
            self.L2A_ancillary['ozone'] = np.float(
                product_meta.getElement('Quality_Indicators_Info').getElement('Image_Content_QI').getAttribute(
                    'OZONE_VALUE').getData().getElemString())
            self.L2A_ancillary['ozone_source'] = product_meta.getElement('Quality_Indicators_Info').getElement(
                'Image_Content_QI').getAttribute('OZONE_SOURCE').getData().getElemString()  # AUX_ECMWFT or CONFIG
            self.L2A_ancillary['aot_method'] = product_meta.getElement('Quality_Indicators_Info').getElement(
                'Image_Content_QI').getAttribute(
                'AOT_RETRIEVAL_METHOD').getData().getElemString()  # SEN2COR_DDV, CAMS, DEFAULT
            # For altitude, SNAP cannot parse L2A_QUALITY.xml, We get only a mean altitude...limited for mountains!
            zf = zipfile.ZipFile(self.prod_path, 'r')
            l2aquality_path = [i for i in zf.namelist() if 'L2A_QUALITY.xml' in i][0]
            l2aquality = zf.read(l2aquality_path)
            root = ET.fromstring(l2aquality)
            h = np.float([i.text for i in root[1].findall('.//') if i.attrib == {'name': 'DEM_MEAN_ALTITUDE_KM'}][0])
            self.L2A_ancillary['altitude'] = h
            # [S2-PDGS-MPC-L2A-IODD-V2.8.pdf] LUT convention: 16 characters/numbers followed by the extension ‘.atm’.
            # The first character defines the atmospheric temperature profile (h=summer, w=winter) and ozone content,
            # followed by ‘99000’ (indicating the symbolic satellite height of 99,000 m), followed by ‘_’, then ‘wvxy’
            # where xy is the sea-level water vapour column, followed by ‘_’ and a 4 letter aerosol identifier ‘_rura’.
            lut = [i.text for i in root[1].findall('.//') if i.attrib == {'name': 'LUT_DATA_FILES'}][0]
            self.L2A_ancillary['atm_season'] = lut[2]
            self.L2A_ancillary['atm_type'] = lut[14:18]
            return
        else:
            raise RuntimeError('Source product must be of type ' + S2_MSI_L1C_TYPE_STRING + ' or ' +
                               S2_MSI_L2A_TYPE_STRING[0] + 'or' + S2_MSI_L2A_TYPE_STRING[1])

        # sensing time of S2 is not correctly implemented. PSD definition is correct but it's equal to "GENERATION_TIME"
        # This is the first scan line of the image. (This is a UTC time)
        sensing_time = product_meta.getElement('General_Info').getElement('Product_Info').getAttribute(
            'PRODUCT_START_TIME').getData().getElemString()
        self.sensing_time = datetime.datetime.strptime(sensing_time, '%Y-%m-%dT%H:%M:%S.%fZ')

        # S2 L1C routines for each band
        for band in self.selected_bands:
            if not band in S2_BAND_NAMES:  # The band name is checked to confirm it is valid band.
                raise RuntimeError(
                    'Source band ' + band + ' is not a S2 band. Valid S2 bands are:\n' + ', '.join(S2_BAND_NAMES))

            # execute the ROI of reflectance and masks for each one of the bands
            source_band = source_product.getBand(band)
            wpix = int(round(self.w / S2_BAND_SAMPLING[band]))
            hpix = int(round(self.h / S2_BAND_SAMPLING[band]))
            if any([wpix, hpix]) == 0:  # minimum pixel size is 1
                wpix = 1
                hpix = 1
            x_off, y_off = self.get_corners(wpix, hpix, source_band)
            roi_data = np.zeros(wpix * hpix, np.float32)
            source_band.readPixels(x_off, y_off, wpix, hpix, roi_data)
            roi_data.shape = wpix, hpix
            self.L1C_ref.append(roi_data)
        self.cloud = self.get_cloud('B_opaque_clouds', source_product)
        self.cirrus = self.get_cloud('B_cirrus_clouds', source_product)
        self.get_roiangles()
        self.radiance_convert(product_meta)
        # self.ecmwf_extract()
        self.get_metadata(product_meta, datastrip_meta)

    def get_roiangles(self):
        '''
        NOTE: must be call only after get_roidata() has been processed!!!
        Obtains the angular source band by resampling at the image band spatial resolution the original 22x22 angles
        Interpolation is bilinear since 1) the angular change can be easily linearised and 2) no
        kernel edge effect will exist for VZA and VAA:
        http://forum.step.esa.int/t/sentinel-2-viewing-angles-interpolation-on-per-pixel-basis/2776
        The last row/colum pixels might produce a NaN since bilinear does not extrapolate values yet.
        In that case you might need to specify the interpolation as: parameters.put('upsampling', 'Nearest')
        '''
        # select the specific image and masks for the band
        source_product = snappy.ProductIO.readProduct(self.prod_path)  # S2 L1C or L2A product
        # checks that the product is the correct one. There are two types of L1C tag names and one of L1C.
        if source_product.getProductType() != S2_MSI_L1C_TYPE_STRING:
            raise RuntimeError('Source product must be of type ' + S2_MSI_L1C_TYPE_STRING)

        # S2 L1C routines for each band
        for band in self.selected_bands:
            if not band in S2_BAND_NAMES:  # The band name is checked to confirm it is valid band.
                raise RuntimeError(
                    'Source band ' + band + ' is not a S2 band. Valid S2 bands are:\n' + ', '.join(S2_BAND_NAMES))
            source_band = source_product.getBand(band)
            # execute the ROI of reflectance and masks for each one of the bands
            wpix = int(round(self.w / S2_BAND_SAMPLING[band]))
            hpix = int(round(self.h / S2_BAND_SAMPLING[band]))
            if any([wpix, hpix]) == 0:  # minimim pixel size is 1
                wpix = 1
                hpix = 1
            x_off, y_off = self.get_corners(wpix, hpix, source_band)

            parameters = hash()
            parameters.put('targetResolution', S2_BAND_SAMPLING[band])
            parameters.put('upsampling', 'Bilinear')
            parameters.put('downsampling', 'Mean')  # indiferent since angles will be always upsampled
            parameters.put('flagDownsampling', 'FlagMedianAnd')
            parameters.put('resampleOnPyramidLevels', True)

            product = GPF.createProduct('Resample', parameters, source_product)
            for band_angle in ['sun_zenith', 'sun_azimuth', 'view_zenith_' + band, 'view_azimuth_' + band]:
                source_angle = product.getBand(band_angle)

                roi_data = np.zeros(wpix * hpix, np.float32)
                source_angle.readPixels(x_off, y_off, wpix, hpix, roi_data)
                roi_data.shape = wpix, hpix

                if 'view_azimuth' in source_angle.getName():
                    self.view_azimuth.append(roi_data)
                elif 'view_zenith' in source_angle.getName():
                    self.view_zenith.append(roi_data)
                elif source_angle.getName() == 'sun_azimuth':
                    self.sun_azimuth.append(roi_data)
                elif source_angle.getName() == 'sun_zenith':
                    self.sun_zenith.append(roi_data)

    def parse_product(self, source_product, leveltag):
        '''
        Obtains the S2 L1C/L2A product structure and brings the metadata
        '''

        # S2 L1C/L2A metadata
        metadata_root = source_product.getMetadataRoot()
        product_meta = metadata_root.getElement(leveltag + '_User_Product')
        datastrip_meta = metadata_root.getElement(leveltag + '_DataStrip_ID')
        granules = metadata_root.getElement('Granules')
        # todo - check if there is a granule

        # with the >PSD14 format only one granule will exist. Initially contain several granules, a list is required.
        granules_meta = []
        for granule in granules.getElements():
            granules_meta.append(granule)
        return (product_meta, datastrip_meta, granules_meta)

    def get_corners(self, wpix, hpix, source_band):
        '''
        Convert lat/lon centre coordinates of a ROI in upper corner pixel coordinates
        '''

        # print (scene_width,scene_height)
        geo_pos = snappy.GeoPos()
        geo_pos.lat = self.lat_centre
        geo_pos.lon = self.lon_centre
        pix_pos = snappy.PixelPos()
        geo_code = source_band.getGeoCoding()
        geo_code.getPixelPos(geo_pos, pix_pos)

        # upper corner of the ROI need to subtract half the ROI size.
        x_off = int(pix_pos.getX()) - int(wpix / 2)
        y_off = int(pix_pos.getY()) - int(hpix / 2)

        return (x_off, y_off)

    def get_cloud(self, band, source_product):
        source_band = source_product.getBand(band)
        wpix = int(round(self.w / 20))
        hpix = int(round(self.h / 20))
        if any([wpix, hpix]) == 0:  # minimum pixel size is 1
            wpix = 1
            hpix = 1
        x_off, y_off = self.get_corners(wpix, hpix, source_band)
        roi_data = np.zeros(wpix * hpix, np.float32)
        source_band.readPixels(x_off, y_off, wpix, hpix, roi_data)
        roi_data.shape = wpix, hpix
        return (roi_data)

    def radiance_convert(self, product_meta):
        # Get metadata values to convert to radiance
        u_sun = product_meta.getElement('General_Info').getElement(
            'Product_Image_Characteristics').getElement('Reflectance_Conversion').getAttributeDouble('U')
        for band in range(len(self.selected_bands)):
            bandind = np.where(np.array(S2_BAND_NAMES) == self.selected_bands[band])[0][0]
            e_sun = float([i for i in product_meta.getElement('General_Info').getElement(
                'Product_Image_Characteristics').getElement('Reflectance_Conversion').getElement(
                'Solar_Irradiance_list').getAttributes() if i.getName() == 'SOLAR_IRRADIANCE'][
                              bandind].getData().getElemString())
            # Convert TOA reflectance to TOA radiance
            self.L1C_rad.append(
                (e_sun * u_sun * np.cos(np.radians(self.sun_zenith[band])) / math.pi) * np.mean(self.L1C_ref[band]))

    def get_metadata(self, product_meta, datastrip_meta):
        self.metadatadict['spacecraft'] = datastrip_meta.getElement('General_Info').getElement(
            'Datatake_Info').getAttributeString('SPACECRAFT_NAME')
        self.metadatadict['quant'] = product_meta.getElement('General_info').getElement(
            'Product_Image_Characteristics').getAttributeDouble('QUANTIFICATION_VALUE')
        for band in range(len(self.selected_bands)):
            bandind = np.where(np.array(S2_BAND_NAMES) == self.selected_bands[band])[0][0]
            self.metadatadict['beta'].append([i for i in
                                              datastrip_meta.getElement('Quality_Indicators_Info').getElement(
                                                  'Radiometric_Info').getElement(
                                                  'Radiometric_Quality_list').getElements() if
                                              i.getName() == 'Radiometric_Quality'][bandind].getElement(
                'Noise_Model').getAttributeDouble('BETA'))

            self.metadatadict['alpha'].append(
                [i for i in datastrip_meta.getElement('Quality_Indicators_Info').getElement('Radiometric_Info').
                getElement('Radiometric_Quality_list').getElements() if i.getName() == 'Radiometric_Quality'][bandind]
                .getElement('Noise_Model').getAttributeDouble('ALPHA'))

            self.metadatadict['A'].append(
                [i for i in datastrip_meta.getElement('Image_Data_Info').getElement('Sensor_Configuration').
                getElement('Acquisition_Configuration').getElement('Spectral_Band_Info').getElements()
                 if i.getName() == 'Spectral_Band_Information'][bandind].getAttributeDouble('PHYSICAL_GAINS'))

            self.metadatadict['Esun'].append(float([i for i in product_meta.getElement('General_Info').getElement(
                'Product_Image_Characteristics').getElement('Reflectance_Conversion').getElement(
                'Solar_Irradiance_list').getAttributes() if i.getName() == 'SOLAR_IRRADIANCE'][
                                                       bandind].getData().getElemString()))

    # def ecmwf_extract(self):
    #     '''
    #     Obtains nformation of water content, pressure and ozone in ECMWF data embedded in S2 L1C product
    #     NOTE: it uses GDAL library. http://www.gdal.org/frmt_grib.html
    #
    #     S2 PSD 14.9
    #     The ECMWF auxiliary data embedded in the Level-1C at Tile level includes the following parameters:
    #         - Total column ozone (TCO3) [Kg/m2];
    #         - Total column water vapour (TCWV) [Kg/m2];
    #         - Mean sea level pressure (MSL) [hPa];.
    #         - 10-meter u/v wind components (10 u and 10 v) [m s-1];
    #         - Relative Humidity (at surface pressure) [%]
    #
    #     The CAMS auxiliary data embedded in the Level-1C at Tile level includes the following parameters:
    #         - Total Aerosol Optical Depth at 550nm (aod550): mandatory parameter
    #         - Surface geopotential (z): mandatory parameter
    #         - Black Carbon Aerosol Optical Depth at 550nm (bcaod550)
    #         - Dust Aerosol Optical Depth at 550nm (duaod550)
    #         - Organic Matter Aerosol Optical Depth at 550nm (omaod550)
    #         - Sea Salt Aerosol Optical Depth at 550nm (ssaod550)
    #         - Sulphate Aerosol Optical Depth at 550nm (suaod550
    #         - Total Aerosol Optical Depth at 469nm (aod469)
    #         - Total Aerosol Optical Depth at 670nm (aod670)
    #         - Total Aerosol Optical Depth at 865nm (aod865)
    #         - Total Aerosol Optical Depth at 1240nm (aod1240)
    #
    #     Resulting from a temporal and spatial interpolation of the raw ECMWF global forecast dataset, this
    #     data will be provided as part of the Level-1C auxiliary data resampled and distributed in grid
    #     information tiles with the same dimensions as the Level-1C Tiles. Grid points are provided in
    #     latitude/longitude using WGS84 reference system.
    #     They are interpolated from original ECMWF data to match L1C Tiles both temporally (linear) and
    #     geometrically (bilinear with a Ground Sample Distance of 12.5km) and provided in GRIB V1 format
    #
    #     Just used the PYGRIB to check the correct order of ECMWF file:
    #         1:Total column water vapour:kg m**-2
    #         2:Mean sea level pressure:Pa
    #         3:Total column ozone:kg m**-2
    #         4:10 metre U wind component:m s**-1
    #         5:10 metre V wind component:m s**-1
    #         6:Relative humidity:% (instant)
    #     and CAMS file:
    #         1:Geopotential:m**2 s**-2
    #         2:Total Aerosol Optical Depth at 469nm:~
    #         3:Total Aerosol Optical Depth at 550nm:~
    #         4:Total Aerosol Optical Depth at 670nm:~
    #         5:Total Aerosol Optical Depth at 865nm:~
    #         6:Total Aerosol Optical Depth at 1240nm:~
    #         7:Black Carbon Aerosol Optical Depth at 550nm:~
    #         8:Dust Aerosol Optical Depth at 550nm:~
    #         9:Organic Matter Aerosol Optical Depth at 550nm:~
    #         10:Sea Salt Aerosol Optical Depth at 550nm:~
    #         11:Sulphate Aerosol Optical Depth at 550nm:~
    #     '''
    #     # STEP 1. Setup the path and read the ECMWF and CAMS file
    #     # TODO - include a decompression of zip beforehand
    #     aux_path = glob.glob(os.path.join(self.prod_path[:-4] + '.SAFE', 'GRANULE', 'L1C_*', 'AUX_DATA'))[0]
    #     ecmwf = gdal.Open(os.path.join(aux_path, 'AUX_ECMWFT'))
    #     for lab in ['AUX_CAMSFO', 'AUX_CAMSAN', 'AUX_CAMSRE']:
    #         try:
    #             cams = gdal.Open(os.path.join(aux_path, lab))
    #         except:
    #             pass
    #         if cams is not None:
    #             break  # if exists, exits the loop
    #     # STEP 2. Find the nearest pixel. It is the same for all bands and both ECMWF and CAMS.
    #     # TODO - include an interpolation at S2 pixel level (might come in future SNAP S2 toolbox)
    #     geo_matrix = ecmwf.GetGeoTransform()
    #     ul_x = geo_matrix[0]
    #     ul_y = geo_matrix[3]
    #     x_dist = geo_matrix[1]
    #     y_dist = geo_matrix[5]
    #     x = int((self.lon_centre - ul_x) / x_dist)
    #     y = -int((ul_y - self.lat_centre) / y_dist)
    #
    #     # STEP 3. Populate the dictionary with each band
    #     ecmwf_tags = ['TCWV', 'MSL', 'TCO3', 'U10', 'V10', 'humidity']
    #     for i in range(6):
    #         band = ecmwf.GetRasterBand(i + 1)
    #         val = band.ReadAsArray()
    #         self.atmosdict[ecmwf_tags[i]] = val[x, y]
    #
    #     cams_tags = ['z', 'aod469', 'aod550', 'aod670', 'aod865', 'aod1240', 'bcaod550', 'duaod550', 'omaod550',
    #                  'ssaod550', 'suaod550', ]
    #     for i in range(11):
    #         band = cams.GetRasterBand(i + 1)
    #         val = band.ReadAsArray()
    #         self.atmosdict[cams_tags[i]] = val[x, y]
