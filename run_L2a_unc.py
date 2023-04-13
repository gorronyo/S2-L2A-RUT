from datetime import datetime
import numpy as np
import L2a_unc
import plots_L2aunc

samples = 1000
librad_bin = None # path to folder with libradtran binary. E.g. '/home/gorrono/libRadtran-2.0.4'
lut = L2a_unc.L2aUnc(librad_bin, samples)
case = 'user'  # 'user', 'libya4' , 'amazon', 'winterwheat' or 'maize'
subcase = 'standard' # 'standard', 'curuc' and 'reptrantest' options for cases 'libya4' and 'amazon'

if case == 'user':  # user defined case.
    path_l1c = None # FULLPATH_TO_L1C.zip'
    path_l2a = None # FULLPATH_TO_L2A.zip'
    latlon = (None, None) # LATITUDE AND LONGITUDE
    roisize = (None, None) # AREA SIZE IN METERS
    toairrad_flag = False
    adjacency_flag = True
    lambertian_flag = True
    libradunc_flag = True
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                      libradunc_flag)
    lut.plot_results(toairrad_flag)
elif case == 'libya4':  # Libya4, PICS with sand dunes in June northern hemisphere (very high radiance).
    path_l1c = '/home/gorrono/Desktop/UPV/Methane/NASdata/S2B_MSIL1C_20220621T085559_N0400_R007_T34RGS_20220621T093642.zip'
    path_l2a = '/home/gorrono/Desktop/UPV/Methane/NASdata/S2B_MSIL2A_20220621T085559_N0400_R007_T34RGS_20220621T105401.zip'
    latlon = (28.55, 23.39)
    roisize = (500, 500)
    toairrad_flag = False
    adjacency_flag = True
    lambertian_flag = True
    libradunc_flag = True
    # Definition of uncertainty setup based on subcase
    if subcase=='standard':
        toairrad_flag = False
        adjacency_flag = True
        lambertian_flag = True
        libradunc_flag = True
    elif subcase=='curuc':
        toairrad_flag = False
        adjacency_flag = False # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
        lambertian_flag = False
        libradunc_flag = False
    elif subcase=='reptrantest':
        toairrad_flag = True  # in True runs special configurations for test purposes.
        adjacency_flag = False
        lambertian_flag = False
        libradunc_flag = False
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                      libradunc_flag)
    lut.plot_results(toairrad_flag)
elif case == 'amazon':  # Amazon forest.
    path_l1c = '/home/gorrono/Desktop/UPV/Methane/NASdata/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916.zip'
    path_l2a = '/home/gorrono/Desktop/UPV/Methane/NASdata/S2B_MSIL2A_20221107T144729_N0400_R139_T19MGP_20221107T183144.zip'
    latlon = (-6.2869, -66.8652)
    roisize = (500, 500)
    # Definition of uncertainty setup based on subcase
    if subcase == 'standard':
        toairrad_flag = False
        adjacency_flag = True
        lambertian_flag = True
        libradunc_flag = True
    elif subcase == 'curuc':
        toairrad_flag = False
        adjacency_flag = False  # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
        lambertian_flag = False
        libradunc_flag = False
    elif subcase == 'reptrantest':
        toairrad_flag = False  # in True runs special configurations for test purposes.
        adjacency_flag = False
        lambertian_flag = False
        libradunc_flag = False
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                      libradunc_flag)
    lut.plot_results(toairrad_flag)
elif case == 'winterwheat':  # Winterwheat field in Winterthur (Switzerland).
    paths_l1c = [
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220322T101711_N0400_R065_T32TMT_20220322T123834.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220421T101601_N0400_R065_T32TMT_20220421T123339.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220620T102041_N0400_R065_T32TMT_20220620T140020.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2B_MSIL1C_20220725T101559_N0400_R065_T32TMT_20220725T122447.zip'
        ]
    paths_l2a = [
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220322T101711_N0400_R065_T32TMT_20220322T141030.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220421T101601_N0400_R065_T32TMT_20220421T134744.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220620T102041_N0400_R065_T32TMT_20220620T162319.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2B_MSIL2A_20220725T101559_N0400_R065_T32TMT_20220725T132135.zip'
        ]
    latlon = (47.452876567, 8.696025314)
    roisize = (60, 60)
    # Definition of uncertainty setup
    toairrad_flag = False  # this should be typically False unless you want to run special configurations for test purposes.
    adjacency_flag = True
    lambertian_flag = True
    libradunc_flag = True
    ndvi_samp = []  # contains the samples of NDVI for each of the four acquisitions
    evi_samp = []  # contains the samples of EVI for each of the four acquisitions
    ndvi_sampcorr = []  # contains the samples of NDVI for each of the four acquisitions with full spectral correlation
    evi_sampcorr = []  # contains the samples of EVI for each of the four acquisitions with full spectral correlation
    ndvi_sampuncorr = []  # contains the samples of NDVI for each of the four acquisitions with no spectral correlation
    evi_sampuncorr = []  # contains the samples of EVI for each of the four acquisitions with no spectral correlation
    meas_dates = []  # contains the dates for each of the four acquisitions
    for path_l1c, path_l2a in zip(paths_l1c, paths_l2a):
        lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                          libradunc_flag)
        lut.plot_results(toairrad_flag)
        # ndvi = (B8-B4)/(B8+B4); evi = 2.5 * (B8-B4) / (B8 + 6*B4 − 7.5*B2 + 1)
        ndvi_samp.append((lut.L2Arho[7, :] - lut.L2Arho[3, :]) / (lut.L2Arho[7, :] + lut.L2Arho[3, :]))
        evi_samp.append(2.5 * (lut.L2Arho[7, :] - lut.L2Arho[3, :]) / (
                    lut.L2Arho[7, :] + 6 * lut.L2Arho[3, :] - 7.5 * lut.L2Arho[1, :] + 1))
        b2uncorr = np.random.normal(loc=np.mean(lut.L2Arho[1, :]), scale=np.std(lut.L2Arho[1, :]), size=samples)
        b4uncorr = np.random.normal(loc=np.mean(lut.L2Arho[3, :]), scale=np.std(lut.L2Arho[3, :]), size=samples)
        b8uncorr = np.random.normal(loc=np.mean(lut.L2Arho[7, :]), scale=np.std(lut.L2Arho[7, :]), size=samples)
        sampscorr = np.random.normal(loc=0, scale=1, size=samples)
        b2corr = sampscorr * np.std(lut.L2Arho[1, :]) + np.mean(lut.L2Arho[1, :])
        b4corr = sampscorr * np.std(lut.L2Arho[3, :]) + np.mean(lut.L2Arho[3, :])
        b8corr = sampscorr * np.std(lut.L2Arho[7, :]) + np.mean(lut.L2Arho[7, :])

        ndvi_sampcorr.append((b8corr - b4corr) / (b8corr + b4corr))
        evi_sampcorr.append(2.5 * (b8corr - b4corr) / (b8corr + 6 * b4corr - 7.5 * b2corr + 1))
        ndvi_sampuncorr.append((b8uncorr - b4uncorr) / (b8uncorr + b4uncorr))
        evi_sampuncorr.append(2.5 * (b8uncorr - b4uncorr) / (b8uncorr + 6 * b4uncorr - 7.5 * b2uncorr + 1))

        meas_dates.append(datetime.strptime(
            lut.paramdict['tagfile'][11:15] + '-' + lut.paramdict['tagfile'][15:17] + '-' + lut.paramdict['tagfile'][
                                                                                            17:19], '%Y-%m-%d'))

    plots_L2aunc.evi_ndvi_trenduncertainty(ndvi_samp, evi_samp, meas_dates, ndvi_sampcorr, evi_sampcorr,
                                           ndvi_sampuncorr, evi_sampuncorr, 'ndvi_evi_uncertaintytrendwinterwheat')

elif case == 'maize':  # Maize field in Winterthur (Switzerland).
    paths_l1c = [
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220322T101711_N0400_R065_T32TMT_20220322T123834.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220421T101601_N0400_R065_T32TMT_20220421T123339.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL1C_20220620T102041_N0400_R065_T32TMT_20220620T140020.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2B_MSIL1C_20220725T101559_N0400_R065_T32TMT_20220725T122447.zip'
    ]
    paths_l2a = [
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220322T101711_N0400_R065_T32TMT_20220322T141030.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220421T101601_N0400_R065_T32TMT_20220421T134744.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2A_MSIL2A_20220620T102041_N0400_R065_T32TMT_20220620T162319.zip',
        '/home/gorrono/Desktop/UPV/L2ARUT/casestudy/S2B_MSIL2A_20220725T101559_N0400_R065_T32TMT_20220725T132135.zip'
    ]
    latlon = (47.451120279, 8.697838380)
    roisize = (60, 60)
    # Definition of uncertainty setup
    toairrad_flag = False  # this should be typically False unless you want to run special configurations for test purposes.
    adjacency_flag = True
    lambertian_flag = True
    libradunc_flag = True
    ndvi_samp = []  # contains the samples of NDVI for each of the four acquisitions
    evi_samp = []  # contains the samples of EVI for each of the four acquisitions
    ndvi_sampcorr = []  # contains the samples of NDVI for each of the four acquisitions with full spectral correlation
    evi_sampcorr = []  # contains the samples of EVI for each of the four acquisitions with full spectral correlation
    ndvi_sampuncorr = []  # contains the samples of NDVI for each of the four acquisitions with no spectral correlation
    evi_sampuncorr = []  # contains the samples of EVI for each of the four acquisitions with no spectral correlation
    meas_dates = []  # contains the dates for each of the four acquisitions
    for path_l1c, path_l2a in zip(paths_l1c, paths_l2a):
        lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, toairrad_flag, adjacency_flag, lambertian_flag,
                          libradunc_flag)
        lut.plot_results(toairrad_flag)
        # ndvi = (B8-B4)/(B8+B4); evi = 2.5 * (B8-B4) / (B8 + 6*B4 − 7.5*B2 + 1)
        ndvi_samp.append((lut.L2Arho[7, :] - lut.L2Arho[3, :]) / (lut.L2Arho[7, :] + lut.L2Arho[3, :]))
        evi_samp.append(2.5 * (lut.L2Arho[7, :] - lut.L2Arho[3, :]) / (
                lut.L2Arho[7, :] + 6 * lut.L2Arho[3, :] - 7.5 * lut.L2Arho[1, :] + 1))
        b2uncorr = np.random.normal(loc=np.mean(lut.L2Arho[1, :]), scale=np.std(lut.L2Arho[1, :]), size=samples)
        b4uncorr = np.random.normal(loc=np.mean(lut.L2Arho[3, :]), scale=np.std(lut.L2Arho[3, :]), size=samples)
        b8uncorr = np.random.normal(loc=np.mean(lut.L2Arho[7, :]), scale=np.std(lut.L2Arho[7, :]), size=samples)
        sampscorr = np.random.normal(loc=0, scale=1, size=samples)
        b2corr = sampscorr * np.std(lut.L2Arho[1, :]) + np.mean(lut.L2Arho[1, :])
        b4corr = sampscorr * np.std(lut.L2Arho[3, :]) + np.mean(lut.L2Arho[3, :])
        b8corr = sampscorr * np.std(lut.L2Arho[7, :]) + np.mean(lut.L2Arho[7, :])

        ndvi_sampcorr.append((b8corr - b4corr) / (b8corr + b4corr))
        evi_sampcorr.append(2.5 * (b8corr - b4corr) / (b8corr + 6 * b4corr - 7.5 * b2corr + 1))
        ndvi_sampuncorr.append((b8uncorr - b4uncorr) / (b8uncorr + b4uncorr))
        evi_sampuncorr.append(2.5 * (b8uncorr - b4uncorr) / (b8uncorr + 6 * b4uncorr - 7.5 * b2uncorr + 1))

        meas_dates.append(datetime.strptime(
            lut.paramdict['tagfile'][11:15] + '-' + lut.paramdict['tagfile'][15:17] + '-' + lut.paramdict['tagfile'][
                                                                                            17:19], '%Y-%m-%d'))

    plots_L2aunc.evi_ndvi_trenduncertainty(ndvi_samp, evi_samp, meas_dates, ndvi_sampcorr, evi_sampcorr,
                                           ndvi_sampuncorr, evi_sampuncorr, 'ndvi_evi_uncertaintytrendmaize')

# typical values in the Sen2Cor LUT that can be used as a reference.
# szalut = [0, 10, 20, 30, 40, 50, 60, 70]  # in degrees
# vzalut = [0, 10]  # in degrees
# raalut = [0, 30, 60, 90, 120, 150, 180]  # in degrees
# hlut = [0, 0.5, 1, 1.5, 2, 2.5]  # in km
# ozonesummerlut = [250, 290, 331, 370, 410, 450]  # dobson units
# ozonewinterlut = [250, 290, 330, 377, 420, 460]  # dobson units
# wvsummerlut = [0.4, 1, 2, 2.9, 4, 5]  # in cm
# wvwinterlut = [0.2, 0.4, 0.8, 1.1]  # in cm
# vislut = [5, 7, 10, 15, 23, 40, 80, 120]  # in km
