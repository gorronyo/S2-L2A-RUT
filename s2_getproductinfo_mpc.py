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

from snappy import GPF
from snappy import HashMap as hash

import matplotlib

matplotlib.use('Agg')  # this does not show the plot on the screen

S2_MSI_L1C_TYPE_STRING = 'S2_MSI_Level-1C'
S2_BAND_NAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
S2_BAND_SAMPLING = {'B1': 60, 'B2': 10, 'B3': 10, 'B4': 10, 'B5': 20, 'B6': 20, 'B7': 20, 'B8': 10, 'B8A': 20, 'B9': 60,
                    'B10': 60, 'B11': 20, 'B12': 20}  # product before October 2016 is not automatic


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

        self.metadatadict = {'spacecraft': None, 'quant': None, 'A': [], 'alpha': [],
                             'beta': [], 'Esun': [], 'Usun': None}  # this dictionary contains the relevant metadata parameters

    def get_roidata(self):
        '''
        Executes a set a routines for each of the selected bands
        '''
        # select the specific image and masks for the band
        source_product = snappy.ProductIO.readProduct(self.prod_path)  # S2 L1C or L2A product
        # checks that the product is the correct one. There are two types of L2A tag names and one of L1C.
        if source_product.getProductType() == S2_MSI_L1C_TYPE_STRING:
            product_meta, datastrip_meta, granules_meta = self.parse_product(source_product, 'Level-1C')
        else:
            raise RuntimeError('Source product must be of type ' + S2_MSI_L1C_TYPE_STRING)

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
        self.metadatadict['Usun'] = product_meta.getElement('General_Info').getElement(
            'Product_Image_Characteristics').getElement('Reflectance_Conversion').getAttributeDouble('U')
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

