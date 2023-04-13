#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 10:23:48 2021

@author: Javier Gorroño
"""

import os
import numpy as np
import pandas as pd
import subprocess
import sys

TOAFILE = 'sen2corlutuncertainty_toa'
PATHFILE = 'sen2corlutuncertainty_path'
SPHERFILE = 'sen2corlutuncertainty_toasph'
SURFFILE = 'sen2corlutuncertainty_surface'
TXFILE = 'sen2corlutuncertainty_tx'

PATH_INP_OUT_FILES = os.path.join(os.getcwd(), 'Libradtran_inp_out')

class libradwrapper:
    def __init__(self):
        # Variables used as input and output of the radiative transfer.
        self.librad_wav = None
        self.lpath = None
        self.surf_albedo = None
        self.sph_albedo = None
        self.edir = None
        self.ediff = None
        self.trans_g2toa = None
        self.ltoa = None

        self.input_libradtran = {
            'source': ('solar',),
            'number_of_streams': (8,),
            'output_user': ('lambda', 'uu'),
            'sza': (30.0,),
            'atmosphere_file': ('midlatitude_summer',),
            'mol_abs_param': ('reptran_channel', 'sentinel2a_msi_b08'),
            'mixing_ratio ch4': (1.8,),
            'mixing_ratio co2': (415.0,),
            'mol_modify h2o': (40.0, 'mm'),
            'altitude': (0.2,),
            'albedo': (0.15,),
            'umu': (1.0,),
            'phi': (60.0,),
            'aerosol_default': (),
            'zout': ('TOA',),
            'phi0': (10.0,),
            'rte_solver': ('disort',)
        }

    def libradtran_input(self, band):
        # pipe = open(os.path.join(PATH_INP_OUT_FILES, TOAFILE + band + '.inp'), 'w')
        # for key, values in self.input_libradtran.items():
        #     input_line = '{} {}\n'.format(key, ' '.join([str(v) for v in values]))
        #     pipe.write(input_line)
        # pipe.write('quiet\n')
        # pipe.close()

        pipe = open(os.path.join(PATH_INP_OUT_FILES, PATHFILE + band + '.inp'), 'w')
        for key, values in self.input_libradtran.items():
            if key == 'albedo':
                continue  # for path radiance, we must set surface albedo to zero.
            input_line = '{} {}\n'.format(key, ' '.join([str(v) for v in values]))
            pipe.write(input_line)
        pipe.write('quiet\n')
        pipe.close()

        pipe = open(os.path.join(PATH_INP_OUT_FILES, SPHERFILE + band + '.inp'), 'w')
        for key, values in self.input_libradtran.items():
            # for spherical albedo, surface albedo and viewing angles should be omitted.
            # also 'filter_function_file' and 'output_process'.A bug in Libradtran seems not to let the convolution here
            if key == 'phi' or key == 'phi0' or key == 'albedo' or key == 'filter_function_file' or key == 'output_process':
                continue
            if key == 'output_user':
                values = ['lambda', 'spher_alb']  # only spherical albedo is a relevant output.
            input_line = '{} {}\n'.format(key, ' '.join([str(v) for v in values]))
            pipe.write(input_line)
        pipe.write('disort_spherical_albedo\n')  # This line sets DISORT just to calculate spherical albedo.
        pipe.write('quiet\n')
        pipe.close()

        pipe = open(os.path.join(PATH_INP_OUT_FILES, SURFFILE + band + '.inp'), 'w')
        for key, values in self.input_libradtran.items():
            if key == 'output_user':
                values = ['lambda', 'uu', 'edir', 'edn']  # relevant downwelling (direc and diffuse) irradiance.
            if key == 'zout':
                values = ['sur']  # we set altitude to surface.
            input_line = '{} {}\n'.format(key, ' '.join([str(v) for v in values]))
            pipe.write(input_line)
        pipe.write('quiet\n')
        pipe.close()

        pipe = open(os.path.join(PATH_INP_OUT_FILES, TXFILE + band + '.inp'), 'w')
        for key, values in self.input_libradtran.items():
            if key == 'output_user':
                values = ['lambda', 'uu', 'edir', 'edn']  # also relevant downwelling (direc and diffuse) irradiance.
            if key == 'zout':
                values = ['sur']  # we set altitude to surface.
            if key == 'sza':
                values = np.rad2deg(np.arccos(self.input_libradtran['umu']))  # VZA set to SZA (reciprocity for tx)
            input_line = '{} {}\n'.format(key, ' '.join([str(v) for v in values]))
            pipe.write(input_line)
        pipe.write('output_quantity transmittance\n')  # Sets Etoa=1 and obtain Tx(dir) and Tx(diff) from 'edir','edn'.
        pipe.write('quiet\n')
        pipe.close()

    def libradtran_call(self, band, libradbin):
        # proc1 = subprocess.Popen(
        #     './uvspec < ' + os.path.join(PATH_INP_OUT_FILES, TOAFILE + band + '.inp') + ' > "' + os.path.join(
        #         PATH_INP_OUT_FILES,
        #         TOAFILE + band + '.out') + '"',
        #     shell=True, cwd=UVSPEC_PATH)
        proc2 = subprocess.Popen(
            './uvspec < ' + os.path.join(PATH_INP_OUT_FILES, PATHFILE + band + '.inp') + ' > "' + os.path.join(
                PATH_INP_OUT_FILES,
                PATHFILE + band + '.out') + '"',
            shell=True, cwd=libradbin)
        proc3 = subprocess.Popen(
            './uvspec < ' + os.path.join(PATH_INP_OUT_FILES, SPHERFILE + band + '.inp') + ' > "' + os.path.join(
                PATH_INP_OUT_FILES,
                SPHERFILE + band + '.out') + '"',
            shell=True, cwd=libradbin)
        proc4 = subprocess.Popen(
            './uvspec < ' + os.path.join(PATH_INP_OUT_FILES, SURFFILE + band + '.inp') + ' > "' + os.path.join(
                PATH_INP_OUT_FILES,
                SURFFILE + band + '.out') + '"',
            shell=True, cwd=libradbin)
        proc5 = subprocess.Popen(
            './uvspec < ' + os.path.join(PATH_INP_OUT_FILES, TXFILE + band + '.inp') + ' > "' + os.path.join(
                PATH_INP_OUT_FILES,
                TXFILE + band + '.out') + '"',
            shell=True, cwd=libradbin)
        # proc1.wait()
        proc2.wait()
        proc3.wait()
        proc4.wait()
        proc5.wait()
        # Default output for disort is output_user lambda, edir, edn, eup, uavgdir, uavgdn, uavgup
        # libradtoa = pd.read_csv(os.path.join(PATH_INP_OUT_FILES, TOAFILE + band + '.out'), delimiter='\s+',
        #                         names=['lambd', 'uu'])
        libradpath = pd.read_csv(os.path.join(PATH_INP_OUT_FILES, PATHFILE + band + '.out'), delimiter='\s+',
                                 names=['lambd', 'uu'])
        libradspher = pd.read_csv(os.path.join(PATH_INP_OUT_FILES, SPHERFILE + band + '.out'), delimiter='\s+',
                                  names=['lambd', 'spher_alb'])
        libradsurf = pd.read_csv(os.path.join(PATH_INP_OUT_FILES, SURFFILE + band + '.out'), delimiter='\s+',
                                 names=['lambd', 'uu', 'edir', 'edn'])
        libradtx = pd.read_csv(os.path.join(PATH_INP_OUT_FILES, TXFILE + band + '.out'), delimiter='\s+',
                                 names=['lambd', 'uu', 'edir', 'edn'])

        # # Libradtran output is divided to obtain the transmittance. If 0 value in denominator, convert Nan to 0.
        # trans_g2toa = (libradtoa['uu'].to_numpy().astype('float32') - libradpath['uu'].to_numpy().astype('float32')) / \
        #               libradsurf['uu'].to_numpy().astype('float32')
        # # The ratio and substraction could produce negative or -inf results. In thar case, set them to zero.
        # trans_g2toa[np.isnan(trans_g2toa)] = 0
        # trans_g2toa[trans_g2toa == float('-inf')] = 0
        # trans_g2toa[trans_g2toa < 0] = 0

        return (band, libradpath['lambd'].to_numpy().astype('float32'), libradpath['uu'].to_numpy().astype('float32'),
                libradsurf['edir'].to_numpy().astype('float32'), libradsurf['edn'].to_numpy().astype('float32'),
                libradspher['spher_alb'].to_numpy().astype('float32'), libradtx['edir'].to_numpy().astype('float32'),
                libradtx['edn'].to_numpy().astype('float32'))
