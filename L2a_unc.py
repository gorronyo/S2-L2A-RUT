#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 14:03:11 2020

@author: Javier Gorroño
"""

import subprocess
import numpy as np
import pandas as pd
import os, math
from scipy import interpolate
import multiprocessing

import libradtran_wrapper
import L1c_unc
import utils
import plots_L2aunc

IRRADIANCE_FILE = 'Solar_irradiance_Thuillier_2002'
DATA_FOLDER = os.path.join(os.getcwd(), 'data')
S2SRFFILE = os.path.join(DATA_FOLDER, 'S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.0.xlsx')
LIBRAD_REPTRANCHANNEL = {'B1': 'sentinel2a_msi_b01', 'B2': 'sentinel2a_msi_b02', 'B3': 'sentinel2a_msi_b03',
                         'B4': 'sentinel2a_msi_b04', 'B5': 'sentinel2a_msi_b05', 'B6': 'sentinel2a_msi_b06',
                         'B7': 'sentinel2a_msi_b07', 'B8': 'sentinel2a_msi_b08', 'B8A': 'sentinel2a_msi_b08a',
                         'B9': 'sentinel2a_msi_b09', 'B11': 'sentinel2a_msi_b11', 'B12': 'sentinel2a_msi_b12'}
LREF = {'B1': 129.11, 'B2': 128, 'B3': 128, 'B4': 108, 'B5': 74.6, 'B6': 68.23, 'B7': 66.70, 'B8': 103, 'B8A': 52.39,
        'B9': 8.77, 'B11': 4, 'B12': 1.70}  # reference radiance in (W/m²/sr/μm), same as (mW/m²/sr/nm)
S2BAND_RANGE = {'B1': [429, 458], 'B2': [439, 539], 'B3': [536, 583], 'B4': [645, 685], 'B5': [693, 714],
                'B6': [730, 750], 'B7': [768, 798], 'B8': [759, 909], 'B8A': [847, 882], 'B9': [931, 959],
                'B11': [1538, 1683], 'B12': [2078, 2320]}  # limiting the spectral range to most of the bandpass
S2_CW = {'B1': 443, 'B2': 490, 'B3': 560, 'B4': 665, 'B5': 705, 'B6': 740, 'B7': 783, 'B8': 842, 'B8A': 865,
         'B9': 945, 'B11': 1610, 'B12': 2190}


class L2aUnc:
    def __init__(self, libradfolder, rep):
        '''

        :param libradfolder: (string) folder where Libradtran is installed
        :param rep: number of surface reflectance samples generated
        '''
        self.rep = rep
        self.libradbin = os.path.join(libradfolder, 'bin')  # path to libradtran binary file
        # this defines the sampling of Thuillier irradiance (also the SRF) and the output grid step in Libradtran.
        # The wavelength range is delimited by "wavelength" in the input file.
        self.libradstep = 0.1
        # Input LibRadTran parameterisation. Tuple with value and standard deviation.
        self.ch4 = (1.8, 0.1)
        self.co2 = (400, 40)
        self.ozoneunc = 3  # dobson units. Default 3% uncertainty based on TCO comparison in Rossana, 2011
        self.hunc = 10  # 10% given as a placeholder. Limited info in L2A but mean slope could be used.
        # Seven radiative transfer terms are stored in the database: path radiance, diffuse solar flux (at
        # sensor), direct (beam) irradiance (at sensor), direct and diffuse ground-to-sensor transmittance,
        # spherical albedo, and direct sun-to-ground transmittance.
        self.lpath = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.edir = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.ediff = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.sph_albedo = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.tx_g2toa = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.tx_g2toa_dir = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.tx_g2toa_diff = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))
        self.L2Arho = np.zeros((len(LIBRAD_REPTRANCHANNEL), self.rep))

        self.L1C_rad = None  # this is not predefined but directly managed by external function L1c_unc.py

        self.r = []  # Correlation matrix for different atmospheric functions at each specific band
        self.rbis = []  # this is the same correlation matrix but with combined downwelling irradiance.
        self.r_etoa = []  # same as self.rbis but Ltoa is substituted by Etoa to understand its effect.

        # Correlation matrix between bands for a specific atmospheric function
        self.r_L1Crad = None
        self.r_lpath = None
        self.r_edir = None
        self.r_ediff = None
        self.r_sph_albedo = None
        self.r_trans_g2toa = None
        self.r_L2Arho = None

        self.l2a_unc = []  # this is the theoretical uncertainty for a flat surface that will be compared to MCM

        self.paramdict = None  # simulation parameters

    def get_libradunc(self, path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                      libradunc_flag):
        '''
        :param path_l1c: (string) absolute path to the L1C ZIP file
        :param path_l2a: (string) absolute path to the L2A ZIP file
        :param latlon: (tuple) latitude and longitude in decimal degrees of the area to be processed
        :param roisize: (tuple) length in x and y direction of the rectangle area to be processed
        :param toairrad_flag: String. If True, includes TOA irradiance uncertainty.(can be set to zero internally) and
        False the precompiled S2 reptran band model is used (no irradiance uncertainty).
        Recommended False mode since the other option will take much longer
        :param adjacency_flag: Boolean. If TRUE, includes adjacency correction uncertainty.
        :param lambertian_flag: Boolean. If TRUE, includes Lambertian assumption uncertainty.
        :param libradunc_flag: Boolean. If TRUE, includes Libradtran error.
        :return:
        '''

        # ----------------------------------  STEP 1. L1C TOA radiance + L1/L2 info  -----------------------------------
        l1c_rad, vza, vaa, sza, saa, l1cmetadatadict, aot, wv, ozone, aot_method, h, atm_season, aot_type = self.get_productinfo(
            path_l1c, path_l2a, latlon, roisize)
        # # The 3 lines below where added to force the code to work with maximum value of WV and VIS (AOT) on the LUT.
        # vis = 5 # maximum VIS is 5km
        # aot = np.exp(1.467 - 0.830 * np.log(vis))  # conversion from vis to AOT at 0km height from Thesis Guanter Fig 3.7
        # wv = 5 # maximum water vapour in summer for the LUT
        # Ozone is based on metadata for tile and h is the mean altitude in the tile (for mountain regions limited)
        # Visibility/AOT @550nm uncertainty in S2 L2A data quality report Ref.: OMPC.CS.DQR.002.12-2022, Issue: 57.0
        # the formula is just the uncertainty of the requirement + the offset found during validation.
        if aot_method == 'CAMS':
            aotunc = np.abs(-0.46 * aot + 0.09) + (0.1*aot+0.03)  # cams fallbacksolution
        elif aot_method == 'SEN2COR_DDV':
            aotunc = np.abs(-0.56 * aot + 0.07) + (0.1*aot+0.03)  # ddv uncertainty
        elif aot_method == 'DEFAULT':
            aotunc = np.abs(-0.56 * aot + 0.07) + (0.1*aot+0.03)  # ddv uncertainty
        else:
            print('Not recognised AOT method. Given default value of 0.15')
            aotunc = 0.15

        # wv uncertainty in S2 L2A data quality report Ref.: OMPC.CS.DQR.002.12-2022, Issue: 57.0
        wvunc = np.abs(-0.1 * wv + 0.03) + (0.1*wv+0.2)

        if atm_season == 'h':
            midlat = 'midlatitude_summer'
        elif atm_season == 'w':
            midlat = 'midlatitude_winter'
        else:
            # PDGS products LUT file are from "f" to "u". Seems a level of ozone and automatically given summer profile.
            print('Atm season code is: ' + atm_season + '. Automatically assigned midlatitude summer profile')
            midlat = 'midlatitude_summer'

        librad = libradtran_wrapper.libradwrapper()
        if aot_type == 'rura':
            pass  # we keep the default aerosol in input file that is the rural profile
        elif aot_type == 'mari':
            librad.input_libradtran['aerosol_species_file'] = ('maritime_clean',)
        elif aot_type == 'urba':
            librad.input_libradtran['aerosol_species_file'] = ('urban',)
        elif aot_type == 'dese':
            librad.input_libradtran['aerosol_species_file'] = ('desert',)

        # ---------------------------------------  STEP 2. L1C TOA uncertainty  ----------------------------------------
        l1c_unc_rho = self.get_l1unc(l1c_rad, l1cmetadatadict)  # dispersion % for value of radiance.

        # ------------------------------  STEP 3. Spectral irradiance (Etoa) uncertainty  ------------------------------
        cw_etoa = [370, 500, 650, 850, 950, 1500, 2000, 2500]
        s2cw = np.array([i for i in S2_CW.values()])
        if toairrad_flag:
            # TOA irradiance uncertainty based on values quoted in VNIR (first paper) and SWIR (second paper):
            # Thuillier, G. et al. (1998). The Visible Solar Spectral Irradiance from 350 to 850 nm as Measured by the Solspec Spectrometer During the Atlas I Mission. In: Pap, J.M., Fröhlich, C., Ulrich, R.K. (eds) Solar Electromagnetic Radiation Study for Solar Cycle 22. Springer, Dordrecht. https://doi.org/10.1007/978-94-011-5000-2_4
            # G. Thuillier, et al., “The Solar Spectral Irradiance from 200 to 2400 nm as Measured by the SOLSPEC Spectrometer from the ATLAS and EURECA Missions,” Solar Physics, Vol. 214, No. 1, 2003, pp. 1-22. doi:10.1023/A:1024048429145
            # normalisation of 1.4%  mainly originating from the pyrometer calibration, and not from the PTB scale
            # However, we keep the same values as a conservative approach.
            unc_etoa = [2.15, 1.05, 0.85, 1.35, 1.1, 0.8, 0.65, 0.6]  # divided by 2 from original 2sigma in paper.
            # unc_etoa = [0, 0, 0, 0, 0, 0, 0, 0]  # we can set it to zero to show low impact at surface reflectance.
            # interpolate etoa uncertainty and rescale from %etoa to %l1c_rad
            s2etoaunc = np.interp(s2cw, cw_etoa, unc_etoa)
            # No info spectral correlation. RBF kernel (var=1, gamma=400) high corr depending on wavelength distance.
            corr_etoa = utils.corr_rbf(1, 400, s2cw)
            cov_etoa = np.dot(s2etoaunc.reshape(-1, 1), s2etoaunc.reshape(-1, 1).T) * corr_etoa
            etoa_unc = np.moveaxis(utils.correlated_samples(np.zeros(12), cov_etoa, self.rep), 1, 0)

            # read Thuillier file and generate an interpolated version that we can modify at each iteration
            wav = []
            irrad_interp = []
            for bandname in LIBRAD_REPTRANCHANNEL.keys():
                irrads = np.loadtxt(os.path.join(DATA_FOLDER, IRRADIANCE_FILE))
                func = interpolate.interp1d(irrads[:, 0], irrads[:, 1], kind='linear')
                interpwav = np.arange(S2BAND_RANGE[bandname][0], S2BAND_RANGE[bandname][1] + self.libradstep,
                                      self.libradstep)
                wav.append(interpwav)
                irrads = func(interpwav)
                irrad_interp.append(irrads)  # input sampling of solar irradiance will define the exit range
                # the method resamples the irradiance to SRF at self.libradstep. Faster than SRF resampled to irradiance
                # We tested that alternative approach and shows no remarkable diff. (code at the end of script)
                with open(os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname]), 'w') as f:
                    f.write('# nm \t mW/m2/nm\n')
                    for j in range(irrads.shape[0]):
                        row = str(np.round(interpwav[j], 2)) + '\t' + str(np.round( irrads[j], 6)) + '\n'
                        f.write(row)

                srf = np.loadtxt(os.path.join(DATA_FOLDER, LIBRAD_REPTRANCHANNEL[bandname]))
                # last value is out of range but the SRF is negligible so we set to zero
                func = interpolate.interp1d(srf[:, 0], srf[:, 1], kind='cubic', bounds_error=False, fill_value=0)
                srfinterp = func(interpwav)
                with open(os.path.join(DATA_FOLDER, LIBRAD_REPTRANCHANNEL[bandname] + 'interp'), 'w') as f:
                    f.write('# nm \t srf\n')
                    for j in range(srfinterp.shape[0]):
                        row = str(np.round(interpwav[j], 2)) + '\t' + str(np.round(srfinterp[j], 6)) + '\n'
                        f.write(row)

            # We convert the dispersion in % of radiance to absolute values.
            self.L1C_rad = l1c_rad[:, None] + etoa_unc * l1c_rad[:, None] / 100 + l1c_unc_rho * l1c_rad[:, None] / 100
        else:
            # We convert the dispersion in % of radiance to absolute values.
            self.L1C_rad = l1c_rad[:, None] + l1c_unc_rho * l1c_rad[:, None] / 100

        # ---------------------------------  STEP 4. Generate L2A uncertainty samples  ---------------------------------
        for i in range(self.rep):
            print(i)
            # --------------------------------  STEP 4.1 Define libradtran input file  ---------------------------------
            librad.input_libradtran['atmosphere_file'] = (midlat,)
            aotsample = np.clip(np.round(np.random.normal(aot, aotunc), 4), 0.0001, None)  # could be negative
            librad.input_libradtran['aerosol_set_tau_at_wvl'] = (550.0, aotsample)
            wvsample = np.clip(np.round(np.random.normal(10 * wv, 10 * wvunc), 4), 0.0001, None)  # could be negative
            librad.input_libradtran['mol_modify h2o'] = (wvsample, 'mm')  # S2 unit 'cm', converted to 'mm'.

            librad.input_libradtran['mol_modify O3'] = (
                np.round(np.random.normal(ozone, self.ozoneunc), 2), 'DU')  # Set ozone column
            librad.input_libradtran['mixing_ratio ch4'] = (np.round(np.random.normal(self.ch4[0], self.ch4[1]), 4),)
            librad.input_libradtran['mixing_ratio co2'] = (np.round(np.random.normal(self.co2[0], self.co2[1]), 1),)
            hsample = np.clip(np.abs(np.round(np.random.normal(h, h * (1 + self.hunc / 100)), 3)), 0.0001, None)
            librad.input_libradtran['altitude'] = (hsample,)  # altitude

            # Angular setup.
            librad.input_libradtran['sza'] = (sza,)
            librad.input_libradtran['umu'] = (np.cos(np.radians(vza)),)
            # reference North for S2, South for Libradtran. Clockwise or both cases.
            librad.input_libradtran['phi'] = (vaa,)
            librad.input_libradtran['phi0'] = ((saa + 180) % 360,)

            for bandname, band in zip(LIBRAD_REPTRANCHANNEL.keys(), range(len(LIBRAD_REPTRANCHANNEL))):
                if toairrad_flag:   # For Thuillier case, we redefine the irradiance at each iteration
                    with open(os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname]), 'w') as f:
                        f.write('# nm \t mW/m2/nm\n')
                        for j in range(irrad_interp[band].shape[0]):
                            row = str(np.round(wav[band][j], 2)) + '\t' + str(
                                np.round((1 + etoa_unc[band, i] / 100) * irrad_interp[band][j], 6)) + '\n'
                            f.write(row)

                    librad.input_libradtran['wavelength'] = (S2BAND_RANGE[bandname][0], S2BAND_RANGE[bandname][1])
                    # per_nm libRadtran assumes that the unit of the spectrum is W/(m2 nm)
                    librad.input_libradtran['source'] = (
                        'solar', os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname]), 'per_nm')
                    # default "coarse". Small test showed minor changes with "medium". Use latter if hi-spec computer.
                    librad.input_libradtran['mol_abs_param'] = ('reptran', 'medium')  # fine|medium|coarse
                    librad.input_libradtran['filter_function_file'] = (
                        os.path.join(DATA_FOLDER, LIBRAD_REPTRANCHANNEL[bandname] + 'interp'), 'normalize')
                    librad.input_libradtran['output_process'] = ('integrate',)
                else:
                    librad.input_libradtran['mol_abs_param'] = ('reptran_channel', LIBRAD_REPTRANCHANNEL[bandname])
                librad.libradtran_input(bandname)  # writes the updated input file. For each band generates input files
            # ----------------------------  STEP 4.2 Execute libradtran and extract output  ----------------------------
            # PARALLEL OPTION. Create a pool that runs all bands asynchronous with the number of processes in cpu.
            pool = multiprocessing.Pool(multiprocessing.cpu_count())
            results = [pool.apply_async(librad.libradtran_call, (band, self.libradbin,)) for band in
                       LIBRAD_REPTRANCHANNEL.keys()]
            pool.close()
            pool.join()
            # predefine array of atmospheric variables
            lambd = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            lpath = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            edir = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            ediff = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            sph_albedo = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            tx_g2toa = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            tx_g2toa_dir = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            tx_g2toa_diff = np.zeros(len(LIBRAD_REPTRANCHANNEL))
            values = [(res.get()[0], res.get()[1], res.get()[2], res.get()[3], res.get()[4], res.get()[5], res.get()[6],
                       res.get()[7]) for res in results]

            # First value is the bandname because asynchronous processes might not return in order.
            for val in values:
                bandindex = list(LIBRAD_REPTRANCHANNEL.keys()).index(val[0])
                lambd[bandindex] = val[1][0]
                lpath[bandindex] = val[2][0]
                edir[bandindex] = val[3][0]
                ediff[bandindex] = val[4][0]
                if toairrad_flag == 'Thuillier' or toairrad_flag == 'zero':
                    sph_albedo[bandindex] = self.convolve_sphalbedo(val[5][0], val[0], wav[bandindex])
                else:
                    sph_albedo[bandindex] = val[5][0]
                tx_g2toa_dir[bandindex] = val[6][0]
                tx_g2toa_diff[bandindex] = val[7][0]
                tx_g2toa[bandindex] = val[6][0] + val[7][0]  # tx direct + tx diffuse

            # -----------------------------  STEP 4.3 RHO Lambertian flat uniform surface  -----------------------------
            # calculates rho as lambertian flat and uniform case
            rho = self.calculate_rho(self.L1C_rad[:, i], np.array(lpath), sph_albedo, np.array(tx_g2toa),
                                     np.array(edir), np.array(ediff))
            # --------------------------------  STEP 4.4 Add adjacency correction error --------------------------------
            adj_unc = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3])  # a default 3%
            # Diff with adjacency modelled as 3% with same correlation as in etoa. Mean correction is 0!Not a problem!
            if adjacency_flag:
                corr_adj = utils.corr_rbf(1, 400, s2cw)
                cov_diff_adj = np.dot(adj_unc.reshape(-1, 1), adj_unc.reshape(-1, 1).T) * corr_adj
                rho_diff_adj = utils.correlated_samples(np.zeros(12), cov_diff_adj, 1)
                rho += (tx_g2toa_diff / tx_g2toa_dir) * (rho - (rho + (rho_diff_adj / 100)))
            # ---------------------------------  STEP 4.5 Lambertian assumption error ----------------------------------
            # Hu et al. 1999 and Franch et al. 2013 show values 2-5% depending on the angular setup and optical depth.
            lambertian_unc = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3])  # a default 3%
            # Absence of spectral correlation info. Use same spectral correlation as in etoa
            if lambertian_flag:
                corr_lamb = utils.corr_rbf(1, 400, s2cw)
                cov_diff_lamb = np.dot(lambertian_unc.reshape(-1, 1), lambertian_unc.reshape(-1, 1).T) * corr_lamb
                rho_diff_lamb = utils.correlated_samples(np.zeros(12), cov_diff_lamb, 1)
                rho += (rho * rho_diff_lamb) / 100

            # --------------------------------------  STEP 4.6 Libradtran error ----------------------------------------
            # Based on values in Govaerts 2022. B8> 3%. Here 1.2% excluding 6SV due to limited spectral resolution.
            # B5,B6,B7 not studied but given a 1.3%. B9 not studied but given a 2% due to water vapour effect.
            # Gven as a value range. Thus, we take half of this range to represent uncertainty.
            libradunc = np.array([1.7, 1.8, 1.9, 1.3, 1.3, 1.3, 1.3, 1.2, 0.8, 2, 2.9, 2]) / 2
            if libradunc_flag:
                corr_librad = utils.corr_rbf(1, 400, s2cw)
                cov_diff_librad = np.dot(libradunc.reshape(-1, 1), libradunc.reshape(-1, 1).T) * corr_librad
                rho_diff_librad = utils.correlated_samples(np.zeros(12), cov_diff_librad, 1)
                rho += (rho * rho_diff_librad) / 100

            # -----------------------------  STEP 4.7 Store the samples and next iteration -----------------------------
            self.lpath[:, i] = lpath
            self.edir[:, i] = edir
            self.ediff[:, i] = ediff
            self.sph_albedo[:, i] = sph_albedo
            self.tx_g2toa[:, i] = tx_g2toa
            self.tx_g2toa_dir[:, i] = tx_g2toa_dir
            self.tx_g2toa_diff[:, i] = tx_g2toa_diff
            self.L2Arho[:, i] = rho

        # --------------------------------  STEP 5. Generate atm.function correlation  ---------------------------------
        for band in range(len(LIBRAD_REPTRANCHANNEL)):
            df = pd.DataFrame(np.array([self.L1C_rad[band, :], self.lpath[band, :], self.tx_g2toa[band, :],
                                        self.edir[band, :], self.ediff[band, :], self.sph_albedo[band, :]]).T,
                              columns=[r'$L_{toa}$', r'$L_{path}$', r'$\tau_{up}$', r'$E_{g(dir)}$', r'$E_{g(diff)}$',
                                       r'$s_{atm}$'])
            self.r.append(df.corr())

            df = pd.DataFrame(np.array([self.L1C_rad[band, :], self.lpath[band, :], self.tx_g2toa[band, :],
                                        self.edir[band, :] + self.ediff[band, :], self.sph_albedo[band, :]]).T,
                              columns=[r'$L_{toa}$', r'$L_{path}$', r'$\tau_{up}$', r'$E_{g}$', r'$s_{atm}$'])
            self.rbis.append(df.corr())
            if toairrad_flag == 'Thuillier':
                df = pd.DataFrame(np.array([etoa_unc[band, :], self.lpath[band, :], self.tx_g2toa[band, :],
                                            self.edir[band, :], self.ediff[band, :], self.sph_albedo[band, :]]).T,
                                  columns=[r'$E_{toa}$', r'$L_{path}$', r'$\tau_{up}$', r'$E_{g(dir)}$',
                                           r'$E_{g(diff)}$', r'$s_{atm}$'])
                self.r_etoa.append(df.corr())

        # ----------------------------------  STEP 6. Generate spectral correlation  -----------------------------------
        self.r_L1Crad = self.get_spectralcorrelation(self.L1C_rad[:, :])
        self.r_lpath = self.get_spectralcorrelation(self.lpath[:, :])
        self.r_edir = self.get_spectralcorrelation(self.edir[:, :])
        self.r_ediff = self.get_spectralcorrelation(self.ediff[:, :])
        self.r_sph_albedo = self.get_spectralcorrelation(self.sph_albedo[:, :])
        self.r_trans_g2toa = self.get_spectralcorrelation(self.tx_g2toa[:, :])
        self.r_L2Arho = self.get_spectralcorrelation(self.L2Arho[:, :])

        # this will be useful to store the parameterisation of the simulation.
        self.paramdict = {'tagfile': os.path.basename(path_l1c[0:-4]), 'sza': sza, 'vza': vza, 'saa': saa, 'vaa': saa,
                          'aot': aot, 'wv': wv, 'h': h, 'ozone': ozone, 'atm': midlat, 'toa_irrad': toairrad_flag,
                          'aot_method': aot_method, 'adjacency': adjacency_flag, 'lamb_correction_unc': lambertian_flag,
                          'libradtran_unc': libradunc_flag}

        # generate an L2A uncertainty theory to compare. Due to complex modelling only for Lambertian uniform case.
        if not adjacency_flag and not lambertian_flag and not libradunc_flag:
            self.get_l2aunc_theory()

    def get_spectralcorrelation(self, values):
        '''

        :param values:
        :return:
        '''
        df = pd.DataFrame(values.T, columns=LIBRAD_REPTRANCHANNEL.keys())
        return (df.corr())

    def convolve_sphalbedo(self, values, band, wav):
        '''
        For the spherical albedo, we define the convolution a-posteriori since it seems a bug in Libradtran does not
        allow for a correct convolution.
        :param values: values of spherical albedo from Libradtran. Grid given by "wav", input grid defined by TOA irrad
        :param band:
        :param wav:
        :return:
        '''
        srf = np.loadtxt(os.path.join(DATA_FOLDER, LIBRAD_REPTRANCHANNEL[band]))
        srf_interp = np.interp(wav, srf[:, 0], srf[:, 1])
        return np.sum(values * srf_interp) / np.sum(srf_interp)

    def calculate_rho(self, L1C_rad, lpath, sph_albedo, trans_g2toa, edir, ediff):
        # we also calculate the rho toa dispersion to understand the combined impact on the distribution
        rho = (L1C_rad - lpath) * math.pi * (1 - 0.15 * sph_albedo) / (trans_g2toa * (edir + ediff))
        return rho

    def get_l2aunc_theory(self):
        '''
        Get a theoretical uncertainty based on jacobian matrices that will be used as a validation.
        Here no uncertainty of the TOA radiance is set and tests for Lref and assumed L1c uncertainty = 3%
        :param band: string with selected S2 band
        :return:
        '''
        for bandname, band in zip(LIBRAD_REPTRANCHANNEL.keys(), range(len(LIBRAD_REPTRANCHANNEL))):
            eglob = self.edir[band, :] + self.ediff[band, :]
            s2ltoa_jacob = np.mean(
                math.pi * (1 - 0.15 * (self.sph_albedo[band, :])) / (self.tx_g2toa[band, :] * np.mean(eglob)))
            s2trans_g2toa_jacob = np.mean(-self.L2Arho[band, :] / self.tx_g2toa[band, :])
            eglob_jacob = np.mean(-self.L2Arho[band, :] / np.mean(eglob))
            s2sph_albedo_jacob = np.mean(-0.15 * math.pi * (self.L1C_rad - self.lpath) / (
                    self.tx_g2toa[band, :] * np.mean(eglob)))

            c = np.array([s2ltoa_jacob, -s2ltoa_jacob, s2trans_g2toa_jacob, eglob_jacob, s2sph_albedo_jacob])
            v = np.array(
                [np.std(self.L1C_rad[band, :]), np.std(self.lpath[band, :]), np.std(self.tx_g2toa[band, :]),
                 np.std(self.edir[band, :] + self.ediff[band, :]), np.std(self.sph_albedo[band, :])])

            self.l2a_unc.append(
                np.sqrt(np.sum(np.diag(c) * np.diag(v) * self.rbis[band].values * np.diag(v).T * np.diag(c).T)))

    def get_productinfo(self, path_L1C, path_l2A, latlon, roisize):
        '''

        :param path_L1C:
        :param latlon:
        :param roisize:
        :return:
        '''
        s2proc = utils.get_l1c([i for i in LIBRAD_REPTRANCHANNEL.keys()], latlon, roisize, path_L1C, path_l2A)
        radval = np.array([np.mean(rad) for rad in s2proc.L1C_rad])  # we will calculate over the mean radiance ROI
        # Mean ROI (or pixel) angle. Reference B8. Minimum angle change other bands and uncertainty impact negligible.
        vza = np.mean(s2proc.view_zenith[7])
        vaa = np.mean(s2proc.view_azimuth[7])
        sza = np.mean(s2proc.sun_zenith[7])
        saa = np.mean(s2proc.sun_azimuth[7])
        aot = np.mean(s2proc.L2A_ancillary['quality_aot'])  # aerosol optical thickness
        wv = np.mean(s2proc.L2A_ancillary['quality_wvp'])  # water vapour in cm

        return (radval, vza, vaa, sza, saa, s2proc.metadatadict, aot, wv, s2proc.L2A_ancillary['ozone'],
                s2proc.L2A_ancillary['aot_method'], s2proc.L2A_ancillary['altitude'],
                s2proc.L2A_ancillary['atm_season'], s2proc.L2A_ancillary['aot_type'])

    def get_l1unc(self, radval, l1cmetadatadict):
        '''
        Wrapper that calls the L1C processor and the L1C montecarlo samples generator
        :return:
        '''
        # run an adapted version of the RUT to generate spectral correlated samples.
        # it returns a distribution of percentage error. u_stray_sys is added linear. At Lref mean distribution is 0.3%
        rutl1 = L1c_unc.S2RutOp()
        l1c_uncsamples = rutl1.unc_spectralcorrelation(radval, l1cmetadatadict, self.rep)
        return np.array(l1c_uncsamples)

    def plot_results(self, toairrad_flag):
        plots_L2aunc.plot_spectralcorrelation(self.r_L2Arho, r'$\rho _{surf}$', 'L2Arho', self.paramdict['tagfile'])
        plots_L2aunc.plot_spectralcorrelation(self.r_L1Crad, r'$L_{TOA}$', 'ltoa', self.paramdict['tagfile'])
        for bandname, band in zip(LIBRAD_REPTRANCHANNEL.keys(), range(len(LIBRAD_REPTRANCHANNEL))):
            if not self.l2a_unc:  # no theoretical uncertainty is calculated
                plots_L2aunc.plot_histograms(self.L1C_rad[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$L_{TOA}[mWm^{−2}sr^{−1}nm^{−1}]$', 'ltoa', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.edir[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$E_{g(dir)}[mWm^{−2}nm^{−1}]$', 'edir', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.ediff[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$E_{g(diff)}[mWm^{−2}nm^{−1}]$', 'ediff', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.lpath[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$L_{path}[mWm^{−2}sr^{−1}nm^{−1}]$', 'lpath', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.sph_albedo[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$s_{atm}$', 'sph_albedo', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.tx_g2toa[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                             r'$\tau_{up}$', 'trans_g2toa', self.paramdict['tagfile'])
                plots_L2aunc.plot_histograms(self.L2Arho[band, :], bandname, r'$\rho _{surf}$', 'L2Arho',
                                             self.paramdict['tagfile'])
            else:  # theoretical uncertainty is calculated. Plot includes normal comparison
                plots_L2aunc.plot_histograms_unctheory(self.L1C_rad[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$L_{TOA}[mWm^{−2}sr^{−1}nm^{−1}]$', 'ltoa',
                                                       self.paramdict['tagfile'], self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.edir[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$E_{g(dir)}[mWm^{−2}nm^{−1}]$', 'edir',
                                                       self.paramdict['tagfile'], self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.ediff[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$E_{g(diff)}[mWm^{−2}nm^{−1}]$', 'ediff',
                                                       self.paramdict['tagfile'], self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.lpath[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$L_{path}[mWm^{−2}sr^{−1}nm^{−1}]$', 'lpath',
                                                       self.paramdict['tagfile'], self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.sph_albedo[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$s_{atm}$', 'sph_albedo', self.paramdict['tagfile'],
                                                       self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.tx_g2toa[band, :], LIBRAD_REPTRANCHANNEL[bandname],
                                                       r'$\tau_{up}$', 'trans_g2toa', self.paramdict['tagfile'],
                                                       self.l2a_unc[band])
                plots_L2aunc.plot_histograms_unctheory(self.L2Arho[band, :], bandname, r'$\rho _{surf}$', 'L2Arho',
                                                       self.paramdict['tagfile'], self.l2a_unc[band])

            plots_L2aunc.plot_functioncorrelation(self.r[band], bandname, 'corr', self.paramdict['tagfile'])
            plots_L2aunc.plot_functioncorrelation(self.rbis[band], bandname, 'corrbis', self.paramdict['tagfile'])
            # if toairrad_flag:
            #     # NOTE: empty plot if irradiance uncertainty is zero. Not considered if reptran channel is used
            #     plots_L2aunc.plot_etoacorrelation(
            #         np.array([self.r_etoa[band].values[0, 1:] for band in range(len(LIBRAD_REPTRANCHANNEL))]),
            #         r'$E_{g(TOA)}$', self.paramdict['tagfile'], LIBRAD_REPTRANCHANNEL.keys())

# irrads = np.loadtxt(
#     os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname] + ' (1)'))
# wav.append(irrads[:,0])
# with open(os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname] + 'modified'), 'w') as f:
#     f.write('# nm \t mW/m2/nm\n')
#     for j in range(irrads.shape[0]):
#         row = str(np.round(irrads[j,0], 2)) + '\t' + str(
#             np.round((1 + etoa_unc[band, i] / 100) * irrads[j,1], 6)) + '\n'
#         f.write(row)
# librad.input_libradtran['wavelength'] = (irrads[0,0], irrads[-1,0])
# per_nm libRadtran assumes that the unit of the spectrum is W/(m2 nm)
# librad.input_libradtran['source'] = (
#     'solar', os.path.join(DATA_FOLDER, IRRADIANCE_FILE + LIBRAD_REPTRANCHANNEL[bandname] + 'modified'),
#     'per_nm')
# librad.input_libradtran['filter_function_file'] = (
#     os.path.join(DATA_FOLDER, LIBRAD_REPTRANCHANNEL[bandname] + 'resampled'), 'normalize')
