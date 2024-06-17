import numpy as np
import L2a_unc
import plots_L2aunc
import utils

samples = 500
librad_bin = '/home/larsshared/libRadtran-2.0.5/bin' # path to folder with libradtran binary. E.g. '/home/gorrono/libRadtran-2.0.4'
case = 'amazon'  # 'user', 'libya4' , 'amazon',
subcase = 'curuc'  # 'standard' or 'curuc' options

gum_unc = []
mcm_std = []
mcm_unc = []
gum_unc_percent = []
mcm_std_percent = []
mcm_unc_percent = []
aots = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
wvs = [0, 0.5, 1, 1.5, 2, 2.5, 3]

# aots = [0, 0.25, 0.5, 0.75, 1]
# wvs = [0, 1.5, 3, 4.5, 6]

# aots = [0, 0.25, 0.5, 0.75, 1]
# wvs = [0, 3, 6]
bandnames = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8a', 'B9', 'B11', 'B12']
# for aot in aots:
#     for wv in wvs:
#         lut = L2a_unc.L2aUnc(librad_bin, samples)  # reset the class at each loop to better manage memory
#         lut.aotuser = aot
#         lut.wvuser = wv
#         if case == 'libya4':  # Libya4, PICS with sand dunes in June northern hemisphere (very high radiance).
#             path_l1c = '/home/gorrono/Desktop/UPV/L2ARUT/S2B_MSIL1C_20220621T085559_N0400_R007_T34RGS_20220621T093642.SAFE.zip'
#             path_l2a = '/home/gorrono/Desktop/UPV/L2ARUT/S2B_MSIL2A_20220621T085559_N0400_R007_T34RGS_20220621T105401.SAFE.zip'
#             latlon = (28.55, 23.39)
#             roisize = (500, 500)
#             adjacency_flag = True
#             lambertian_flag = True
#             libradunc_flag = True
#             # Definition of uncertainty setup based on subcase
#             if subcase == 'standard':
#                 adjacency_flag = True
#                 lambertian_flag = True
#                 libradunc_flag = True
#             elif subcase == 'curuc':
#                 adjacency_flag = False  # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
#                 lambertian_flag = False
#                 libradunc_flag = False
#             lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag, libradunc_flag)
#             gum_unc.append(lut.l2a_unc)
#             mcm_std.append(np.std(lut.L2Arho[:, :], axis=1))
#             gum_unc_percent.append(100 * np.abs(np.array(lut.l2a_unc) / np.nanmean(lut.L2Arho[:, :], axis=1)))
#             mcm_std_percent.append(
#                 100 * np.abs(np.std(lut.L2Arho[:, :], axis=1) / np.nanmean(lut.L2Arho[:, :], axis=1)))
#             # we repeat the histogram for fine bins in uncertainty area
#             mcm_uncbands = []
#             for i in range(len(bandnames)):
#                 [pdf, bins] = np.histogram(lut.L2Arho[i, :], bins=10000, normed=True)
#                 [bin_down, bin_up, area] = utils.unc_montecarlo(pdf, bins, np.nanmedian(lut.L2Arho[i, :]))
#                 # print(bin_down, bin_up, area)
#                 mcm_uncbands.append((bin_up - bin_down) / 2)
#             mcm_unc.append(mcm_uncbands)
#             mcm_unc_percent.append(100 * np.abs(np.array(mcm_uncbands) / np.nanmean(lut.L2Arho[:, :], axis=1)))
#
#         elif case == 'amazon':  # Amazon forest.
#             path_l1c = '/home/gorrono/Desktop/UPV/L2ARUT/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916.SAFE.zip'
#             path_l2a = '/home/gorrono/Desktop/UPV/L2ARUT/S2B_MSIL2A_20221107T144729_N0400_R139_T19MGP_20221107T183144.SAFE.zip'
#             latlon = (-6.2869, -66.8652)
#             roisize = (500, 500)
#             # Definition of uncertainty setup based on subcase
#             if subcase == 'standard':
#                 toairrad_flag = False
#                 adjacency_flag = True
#                 lambertian_flag = True
#                 libradunc_flag = True
#             elif subcase == 'curuc':
#                 toairrad_flag = False
#                 adjacency_flag = False  # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
#                 lambertian_flag = False
#                 libradunc_flag = False
#             lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag, libradunc_flag)
#             gum_unc.append(lut.l2a_unc)
#             mcm_std.append(np.std(lut.L2Arho[:, :], axis=1))
#             gum_unc_percent.append(100*np.abs(np.array(lut.l2a_unc)/np.nanmean(lut.L2Arho[:, :], axis=1)))
#             mcm_std_percent.append(100*np.abs(np.std(lut.L2Arho[:, :], axis=1)/np.nanmean(lut.L2Arho[:, :], axis=1)))
#             # we repeat the histogram for fine bins in uncertainty area
#             mcm_uncbands = []
#             for i in range(len(bandnames)):
#                 [pdf, bins] = np.histogram(lut.L2Arho[i, :], bins=10000, normed=True)
#                 [bin_down, bin_up, area] = utils.unc_montecarlo(pdf, bins, np.nanmedian(lut.L2Arho[i,:]))
#                 # print(bin_down, bin_up, area)
#                 mcm_uncbands.append((bin_up - bin_down) / 2)
#             mcm_unc.append(mcm_uncbands)
#             mcm_unc_percent.append(100 * np.abs(np.array(mcm_uncbands) / np.nanmean(lut.L2Arho[:, :], axis=1)))
#
# with open('gum_unc.npy', 'wb') as f:
#     np.save(f, np.array(gum_unc))
# with open('mcm_std.npy', 'wb') as f:
#     np.save(f, np.array(mcm_std))
# with open('mcm_unc.npy', 'wb') as f:
#     np.save(f, np.array(mcm_unc))
# with open('gum_unc_percent.npy', 'wb') as f:
#     np.save(f, np.array(gum_unc_percent))
# with open('mcm_std_percent.npy', 'wb') as f:
#     np.save(f, np.array(mcm_std_percent))
# with open('mcm_unc_percent.npy', 'wb') as f:
#     np.save(f, np.array(mcm_unc_percent))
#
#
#
# for i, bandname in zip(range(len(bandnames)), bandnames):
#     plots_L2aunc.plot_validationmap(bandname, aots, wvs, np.array(gum_unc)[:, i], np.array(mcm_std)[:, i],
#                                     np.array(mcm_unc)[:, i], np.array(gum_unc_percent)[:, i],
#                                     np.array(mcm_std_percent)[:, i], np.array(mcm_unc_percent)[:, i],
#                                     lut.paramdict['tagfile'])

import os
path_l1c = '/home/gorrono/Desktop/UPV/L2ARUT/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916.SAFE.zip'
# aots = [0, 0.25, 0.5, 0.75, 1]
# wvs = [0, 1.5, 3, 4.5, 6]
# mcm_std_percent = np.load('./L2Aunc_results/amazon_validation/mcm_std_percent.npy')
# mcm_unc_percent = np.load('./L2Aunc_results/amazon_validation/mcm_unc_percent.npy')
# gum_unc_percent = np.load('./L2Aunc_results/amazon_validation/gum_unc_percent.npy')
# mcm_std = np.load('./L2Aunc_results/amazon_validation/mcm_std.npy')
# mcm_unc = np.load('./L2Aunc_results/amazon_validation/mcm_unc.npy')
# gum_unc = np.load('./L2Aunc_results/amazon_validation/gum_unc.npy')

aots = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
wvs = [0, 0.5, 1, 1.5, 2, 2.5, 3]
mcm_std_percent = np.load('./L2Aunc_results/amazon_validation_zoom/mcm_std_percent.npy')
mcm_unc_percent = np.load('./L2Aunc_results/amazon_validation_zoom/mcm_unc_percent.npy')
gum_unc_percent = np.load('./L2Aunc_results/amazon_validation_zoom/gum_unc_percent.npy')
mcm_std = np.load('./L2Aunc_results/amazon_validation_zoom/mcm_std.npy')
mcm_unc = np.load('./L2Aunc_results/amazon_validation_zoom/mcm_unc.npy')
gum_unc = np.load('./L2Aunc_results/amazon_validation_zoom/gum_unc.npy')
for i, bandname in zip(range(len(bandnames)), bandnames):
    plots_L2aunc.plot_validationmap(bandname, aots, wvs, np.array(gum_unc)[:, i], np.array(mcm_std)[:, i],
                                    np.array(mcm_unc)[:, i], np.array(gum_unc_percent)[:, i],
                                    np.array(mcm_std_percent)[:, i], np.array(mcm_unc_percent)[:, i],
                                    os.path.basename(path_l1c[0:-9]))

