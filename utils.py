# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:48:33 2016
@author: jg9
"""

import numpy as np
from scipy.stats import multivariate_normal
import L2a_unc
import s2_getproductinfo


def correlated_samples(means, cov_matrix, numbersamples):
    '''
    E.g. for the case of two variables it would be
    cov_matrix = [[(SIGMA1 ** 2), CORR * SIGMA1 * SIGMA2], [CORR * SIGMA1 * SIGMA2, (SIGMA2 ** 2)]]
    distribution = multivariate_normal(mean=[MU1, MU2], cov=cov_matrix)
    samples = distribution.rvs(1000000)
    samples1 = samples[:, 0]
    sample2 = samples[:, 1]
    # A nice example can be found here.
    https://towardsdatascience.com/correlated-variables-in-monte-carlo-simulations-19266fb1cf29
    :param means: list of means e.g. [MU1, MU2].
    :param cov_matrix: covariance matrix (covariance is equal to CORR * SIGMA1 * SIGMA2)
    :param numbersamples: integer with the number of samples that are required.
    :return: array of samples with shape [numbersamples, number variables]
    '''
    # TODO - a multivariate normal is a first approach but there are many other options that will need to consider.
    # see more here https://docs.scipy.org/doc/scipy/reference/stats.html#multivariate-distributions
    distribution = multivariate_normal(mean=means, cov=cov_matrix, allow_singular=True)

    return distribution.rvs(numbersamples)


def samples_check():
    L2Arhos = []
    L2Arhosmean = []
    L2Arhosstd = []
    for sampnum in [10, 25, 50, 100, 200, 300, 400, 500, 750, 1000]:
        lut = L2a_unc.L2aUnc(sampnum)
        lut.get_libradunc(30, 0, 100, 40, 2.9, 0, 331, 'midlatitude_summer', True, 'cams')
        L2Arhos.append(lut.L2Arho)
        L2Arhosmean.append(np.mean(lut.L2Arho, axis=1))
        L2Arhosstd.append(np.std(lut.L2Arho, axis=1))
    return (L2Arhosmean, L2Arhosstd)


def corr_rbf(var, gamma, s2cw):
    '''

    :param var: (float) scaling factor. Typically set to 1 (maximum correlation).
    :param gamma: (float) gamma value of distribution
    :param s2cw: (numpy array) S2 bands CW
    :return:
    '''
    # Set a RBF kernel simulating large corr but depends on wavelength dist.
    X = np.matmul(s2cw.reshape(-1, 1), s2cw.reshape(-1, 1).T)
    X = X / np.abs(X).sum()  # normalise matrix proportional to wavelength distance
    X_norm = np.sum(X ** 2, axis=-1)
    return var * np.exp(-gamma * (X_norm[:, None] + X_norm[None, :] - 2 * np.dot(X, X.T)))


def get_l1c(L1Cbands, latlon, roisize, path_L1C, path_l2A):
    '''
    You can use this function to obtain all parameters from a S2 L1C ROI. For collection 1 products!!!
    NOTE: L99 "self.ecmwf_extract()" in s2_getproductinfo.py is commented due to limitation in ECMWF parsing.
    If ever want to run it (not needed L2a_unc.py), uncompress the ZIP file before calling this function!
    :param L2Abands: (list of strings) set of L1C bands to be processed. Valid values are:
    ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
    :param latlon:(tuple) latitude and longitude in decimal degrees of the area to be processed
    :param roisize: (tuple) length in x and y direction of the rectangle area to be processed
    :param path_L1C: (string) absolute path to the L1C ZIP file
    :param path_L2A: (string) absolute path to the L2A ZIP file. Ypu also get AOT, WV and classification from L2A prod.
    :return: class s2proc. Contains L1C info ncluding atm.param., metadata, radiance conversion... see s2_getL1C.init()
    e.g. how to call this:
    from utils import get_l1c
    path_l1A = 'ABSOLUTE_PATH_TO_S2L1C_PRODUCT/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916.zip'
    path_l2A = 'ABSOLUTE_PATH_TO_S2L1C_PRODUCT/S2B_MSIL2A_20221107T144729_N0400_R139_T19MGP_20221107T183144.zip'
    L1Cbands = ['B12']
    latlon = (42.792954, -105.354420)
    roisize = (500, 500)
    s2proc = get_l1c(L1Cbands, latlon, roisize, path_L1C, pathl2a)
    '''

    s2proc = s2_getproductinfo.S2Processor()
    s2proc.selected_bands = L1Cbands
    s2proc.lat_centre = latlon[0]
    s2proc.lon_centre = latlon[1]
    s2proc.w = roisize[0]
    s2proc.h = roisize[1]
    s2proc.prod_path = path_L1C
    s2proc.get_roidata()
    s2proc.prod_path = path_l2A
    s2proc.get_roidata()
    return s2proc

# def get_l2a(L2Abands, latlon, roisize, path_L2A):

def unc_montecarlo(pdf, bins, meanval):
    up = [t for t, x in enumerate(bins) if x > meanval][0]
    area = sum(np.diff(bins[up - 1:up + 1]) * pdf[up - 1])
    inc = 1
    while area < 0.68268 and inc <= up and up + inc < np.size(bins) - 1 and up - 1 - inc > 0:
        area = sum(np.diff(bins[up - 1 - inc:up + inc + 1]) * pdf[up - 1 - inc:up + inc])
        inc += 1
    bin_down = bins[up - 1 - inc]
    bin_up = bins[up + inc]

    return bin_down, bin_up, 100 * area