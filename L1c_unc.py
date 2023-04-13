# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:48:33 2016
@author: jg9
"""

import math
import numpy as np
import utils

S2_BAND_NAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']


class S2RutOp:
    def __init__(self):
        """
        Contains Sentinel-2 Level-1 radiometric configuration and uncertainty contributions placeholders.
        """

        self.Lref = [129.11, 128, 128, 108, 74.6, 68.23, 66.70, 103, 52.39, 8.77, 6, 4, 1.70]

        self.u_stray_rand_all = {'Sentinel-2A': [0.1, 0.1, 0.08, 0.12, 0.44, 0.16, 0.2, 0.2, 0.04, 0.8, 0, 0, 0],
                                 'Sentinel-2B': [0.1, 0.1, 0.08, 0.12, 0.44, 0.16, 0.2, 0.2, 0.04, 0.8, 0, 0, 0]}

        # These are the estimated error of optical crosstalk before the correction.
        # self.u_xtalk_all = {'Sentinel-2A': [0.05, 0.01, 0.01, 0.01, 0.04, 0.03, 0.04, 0.02, 0.03, 0.02, 0.19, 0.15, 0.02],
        #                     'Sentinel-2B': [0.05, 0.01, 0.01, 0.01, 0.04, 0.03, 0.04, 0.02, 0.03, 0.02, 0.19, 0.15, 0.02]}

        # units in W.m-2.sr-1.μm-1
        self.u_xtalk_all = {
            'Sentinel-2A': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            'Sentinel-2B': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]}

        self.u_DS_all = {'Sentinel-2A': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.24, 0.12, 0.16],
                         'Sentinel-2B': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.24, 0.12, 0.16]}

        # values from ICCDB (S2A at S2_OPER_MSI_DIFF___20150519T000000_0009.xml and
        # S2B at S2_OPER_MSI_DIFF___20160415T000000_0001.xml)
        self.u_diff_absarray = {
            'Sentinel-2A': [1.09, 1.08, 0.84, 0.73, 0.68, 0.97, 0.83, 0.81, 0.88, 0.97, 1.39, 1.39, 1.58],
            'Sentinel-2B': [1.16, 1.00, 0.79, 0.70, 0.85, 0.77, 0.80, 0.80, 0.85, 0.66, 1.70, 1.46, 2.13]}

        # self.u_diff_temp_rate = {'Sentinel-2A': [0.15, 0.09, 0.04, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        #                     'Sentinel-2B': [0.15, 0.09, 0.04, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
        self.u_diff_cos = 0.4  # [%]from 0.13° diffuser planarity/micro as in (AIRBUS 2015). Assumed same for S2A/S2B.
        self.u_diff_k = 0.3  # [%] as a conservative residual (AIRBUS 2015). Assumed same for S2A/S2B.
        self.u_diff_temp = 1.0  # This value is correctly redefined for specific satellite at the S2RutOp.
        self.u_ADC = 0.5  # [DN](rectangular distribution, see combination)
        self.u_gamma = 0.4
        self.k = 1  # This value is correctly redefined for specific satellite at the S2RutOp.
        self.unc_select = [True, True, True, True, True, True, True, True, False, True, True,
                           True]  # list of booleans with user selected uncertainty sources(order as in interface)

    def unc_calculation(self, band_data, band_rad, bandname, bandind, metadatadict):
        """
        This function represents the core of the RUTv1. It takes as an input the pixel data of a specific band and tile
        in a S2-L1C product and produces an image with the same dimensions that contains the radiometric uncertainty of
        each pixel reflectance factor. The steps and its numbering is equivalent to the RUT-DPM. This document can be
        found in the tool github. Also there a more detailed explanation of the theoretical background can be found.
        :param band_data: list with the quantized L1C reflectance pixels of a band
        :param band_rad: list with the quantized L1C radiance pixels of a band
        :param bandname: tag name of the S2 L1C band
        :param bandind: zero-based index of the band
        :param metadatadict: dictionary contains the relevant metadata parameters
        :return: u_int8 with uncertainty associated to each pixel.
        """
        band_id = np.where(np.array(S2_BAND_NAMES) == bandname)[0][0]  # check for index in full L1C band names!

        #######################################################################
        # 1.	Undo reflectance conversion
        #######################################################################
        # a.	No action required
        # b.	[product metadata] #issue: missing one band
        #    General_Info/Product_Image_Characteristics/PHYSICAL_GAINS [bandId]
        #    [datastrip metadata]
        #    Image_Data_Info/Sensor_Configuration/Acquisition_Configuration/
        #    Spectral_Band_Info/Spectral_Band_Information [bandId]/ PHYSICAL_GAINS

        # Replace the reflectance factors by CN values
        # cn = (self.a * self.e_sun * self.u_sun * math.cos(math.radians(self.tecta)) / math.pi) * band_data
        # cn = (self.a * self.e_sun * self.u_sun * np.cos(np.radians(self.tecta)) / math.pi) * band_data
        cn = metadatadict['A'][bandind] * band_rad

        #######################################################################
        # 2.	Orthorectification process
        #######################################################################

        # TBD. Here both terms will be used with no distinction.

        #######################################################################
        # 3.	L1B uncertainty contributors: raw and dark signal
        #######################################################################

        if self.unc_select[0]:
            u_noise = 100 * np.sqrt(metadatadict['alpha'][bandind] ** 2 + metadatadict['beta'][bandind] * cn) / cn
        else:
            u_noise = 0

        # [W.m-2.sr-1.μm-1] 0.3%*Lref all bands (AIRBUS 2015) and (AIRBUS 2014)
        if self.unc_select[1]:
            u_stray_sys = 0.3 * self.Lref[band_id] / 100
        else:
            u_stray_sys = 0

        if self.unc_select[2]:
            u_stray_rand = self.u_stray_rand_all[metadatadict['spacecraft']][
                band_id]  # [%](AIRBUS 2015) and (AIRBUS 2012)
        else:
            u_stray_rand = 0

        if self.unc_select[3]:
            u_xtalk = self.u_xtalk_all[metadatadict['spacecraft']][band_id]  # [W.m-2.sr-1.μm-1](AIRBUS 2015)
        else:
            u_xtalk = 0

        if not self.unc_select[4]:
            self.u_ADC = 0  # predefined but updated to 0 if deselected by user

        if self.unc_select[5]:
            u_DS = self.u_DS_all[metadatadict['spacecraft']][band_id]
        else:
            u_DS = 0

        #######################################################################
        # 4.	L1B uncertainty contributors: gamma correction
        #######################################################################

        if self.unc_select[6]:
            self.u_gamma = 0.4  # [%] (AIRBUS 2015)
        else:
            self.u_gamma = 0

        #######################################################################
        # 5.	L1C uncertainty contributors: absolute calibration coefficient
        #######################################################################

        if self.unc_select[7]:
            u_diff_abs = self.u_diff_absarray[metadatadict['spacecraft']][bandind]
        else:
            u_diff_abs = 0

        if not self.unc_select[8]:
            self.u_diff_temp = 0  # calculated in s2_rut.py. Updated to 0 if deselected by user

        if not self.unc_select[9]:
            self.u_diff_cos = 0  # predefined but updated to 0 if deselected by user

        if not self.unc_select[10]:
            self.u_diff_k = 0  # predefined but updated to 0 if deselected by user

        #######################################################################
        # 6.	L1C uncertainty contributors: reflectance conversion
        #######################################################################

        if self.unc_select[11]:
            u_ref_quant = 100 * (0.5 / math.sqrt(3)) / (
                    metadatadict['quant'] * band_data)  # [%]scaling 0-1 in steps number=quant
        else:
            u_ref_quant = 0

        #######################################################################
        # 7.	Combine uncertainty contributors
        #######################################################################
        # NOTE: no gamma propagation for RUTv1!!!
        # values given as percentages. Multiplied by 10 and saved to 1 byte(uint8)
        # Clips values to 0-250 --> uncertainty >=25%  assigns a value 250.
        # Uncertainty <=0 represents a processing error (uncertainty is positive)
        u_adc = (100 * self.u_ADC / math.sqrt(3)) / cn
        u_ds = (100 * u_DS) / cn
        u_stray = np.sqrt(u_stray_rand ** 2 + ((100 * metadatadict['A'][bandind] * u_xtalk) / cn) ** 2)
        u_diff = math.sqrt(u_diff_abs ** 2 + self.u_diff_cos ** 2 + self.u_diff_k ** 2)
        u_1sigma = np.sqrt(u_ref_quant ** 2 + self.u_gamma ** 2 + u_stray ** 2 + u_diff ** 2 +
                           u_noise ** 2 + u_adc ** 2 + u_ds ** 2)
        u_expand = np.round(
            10 * (self.u_diff_temp + ((100 * metadatadict['A'][bandind] * u_stray_sys) / cn) + self.k * u_1sigma))
        u_ref = np.uint8(np.clip(u_expand, 0, 250))

        return u_ref

    def unc_spectralcorrelation(self, band_rad, metadatadict, rep):
        """
        :param band_data: list with the quantized L1C reflectance values for each band (e.g. at Lref)
        :param band_rad: list with the quantized L1C radiance values for each band (e.g. Lref)
        :param metadatadict: dictionary contains the relevant metadata parameters for typically S2 L2A bands
        :return: samples out of the potential uncertainty L1C for each band. list of arrays [band, #samples]
        """

        #######################################################################
        # 1.	Undo reflectance conversion
        #######################################################################
        cn = np.array(metadatadict['A']) * band_rad
        #######################################################################
        # 2.	Orthorectification process
        #######################################################################
        # TBD. Here both terms will be used with no distinction.
        #######################################################################
        # 3.	L1B uncertainty contributors: raw and dark signal
        #######################################################################
        u_expand = []
        noise = 100 * np.sqrt(np.array(metadatadict['alpha']) ** 2 + np.array(metadatadict['beta']) * cn) / cn
        Lref = np.hstack((self.Lref[0:10], self.Lref[11:13]))  # We remove B10
        u_stray_sys = 0.3 * np.array(Lref) / band_rad  # adapted to L2A bands 0.3% of Lref.
        # [W.m-2.sr-1.μm-1] 0.3%*Lref all bands (AIRBUS 2015) and (AIRBUS 2014)
        # This a known bias added linearly. when a e.g. band ratio is set, the term Lref_i/Lref_j should be kept.
        stray_rand = np.hstack((self.u_stray_rand_all[metadatadict['spacecraft']][0:10],
                                self.u_stray_rand_all[metadatadict['spacecraft']][
                                11:13]))  # [%](AIRBUS 2015) and (AIRBUS 2012)
        xtalk = np.hstack((self.u_xtalk_all[metadatadict['spacecraft']][0:10],
                           self.u_xtalk_all[metadatadict['spacecraft']][11:13]))  # [W.m-2.sr-1.μm-1](AIRBUS 2015)
        u_DSunc = np.hstack((self.u_DS_all[metadatadict['spacecraft']][0:10],
                             self.u_DS_all[metadatadict['spacecraft']][11:13]))

        corr_matrix = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B1 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B2 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B3 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B4 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B5 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B6 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B7 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B8 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B8A vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B9 vs B1-12
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],  # B11 vs B1-12
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])  # B12 vs B1-12
        U = np.array(
            self.u_diff_absarray[metadatadict['spacecraft']][0:10] + self.u_diff_absarray[metadatadict['spacecraft']][
                                                                     11:13])
        cov_diff_abs = np.dot(U.reshape(-1,1), U.reshape(-1,1).T) * corr_matrix

        for i in range(rep):
            u_noise = np.random.normal(0, noise)  # uncorrelated samples for each band
            u_stray_rand = np.random.normal(0, stray_rand)  # uncorrelated samples for each band
            u_xtalk = np.random.normal(0, xtalk)  # uncorrelated samples for each band
            u_DS = np.random.normal(0, u_DSunc)  # correlated samples for each band
            #######################################################################
            # 4.	L1B uncertainty contributors: gamma correction
            #######################################################################
            u_gamma = utils.correlated_samples(np.zeros(12),  self.u_gamma**2 * corr_matrix, 1)  # [%] (AIRBUS 2015)
            #######################################################################
            # 5.	L1C uncertainty contributors: absolute calibration coefficient
            #######################################################################
            u_diff_abs = utils.correlated_samples(np.zeros(12), cov_diff_abs, 1)
            #######################################################################
            # 6.	L1C uncertainty contributors: reflectance conversion
            #######################################################################
            # not included as a simplification. Small impact under really low reflectance e.g. 0.01 reflectance, 0.3% unc.
            # u_ref_quant = 100 * (0.5 / math.sqrt(3)) / (metadatadict['quant'] * band_data)  # [%]scaling 0-1 in steps number=quant
            #######################################################################
            # 7.	Combine uncertainty contributors
            #######################################################################
            # NOTE: no gamma propagation for RUTv1!!!
            # values given as percentages. Multiplied by 10 and saved to 1 byte(uint8)
            # Clips values to 0-250 --> uncertainty >=25%  assigns a value 250.
            # Uncertainty <=0 represents a processing error (uncertainty is positive)
            u_adc = 100 * np.random.uniform(-self.u_ADC, self.u_ADC, 12) / cn
            u_ds = (100 * u_DS) / cn

            u_stray = u_stray_rand * ((100 * np.array(metadatadict['A']) * u_xtalk) / cn)
            u_diff = u_diff_abs * np.random.normal(0, self.u_diff_cos, 12) + np.random.normal(0, self.u_diff_k, 12)
            u_1sigma = u_gamma + u_stray + u_diff + u_noise + u_adc + u_ds
            u_expand.append(u_stray_sys + self.k * u_1sigma)
        return [np.array(u_expand)[:, i] for i in range(np.array(u_expand).shape[1])]


