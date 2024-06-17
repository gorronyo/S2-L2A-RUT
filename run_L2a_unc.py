import L2a_unc

samples = 1000
n_proc = 1  # check max processes with multiprocessing.cpu_count()
librad_bin = None # path to folder with libradtran binary. E.g. '/home/larsshared/libRadtran-2.0.5'
lut = L2a_unc.L2aUnc(librad_bin, samples)
case = 'user'  # 'user', 'libya4' , 'amazon', 'winterwheat' or 'maize'
subcase = 'standard'  # 'standard' and 'curuc' options for cases 'libya4' and 'amazon'

if case == 'user':  # user defined case.
    path_l1c = None  # FULLPATH_TO_L1C.zip'
    path_l2a = None  # FULLPATH_TO_L2A.zip'
    latlon = (None, None)  # LATITUDE AND LONGITUDE
    roisize = (None, None)  # AREA SIZE IN METERS
    adjacency_flag = True
    lambertian_flag = True
    libradunc_flag = True
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag, libradunc_flag, n_proc)
    lut.plot_results()
elif case == 'libya4':  # Libya4, PICS with sand dunes in June northern hemisphere (very high radiance).
    path_l1c = '/home/larsshared/PROJECTS/L2A_RUT/S2products/S2B_MSIL1C_20220621T085559_N0400_R007_T34RGS_20220621T093642.SAFE.zip'
    path_l2a = '/home/larsshared/PROJECTS/L2A_RUT/S2products/S2B_MSIL2A_20220621T085559_N0400_R007_T34RGS_20220621T105401.SAFE.zip'
    latlon = (28.55, 23.39)
    roisize = (500, 500)
    # Definition of uncertainty setup based on subcase
    if subcase=='standard':
        adjacency_flag = True
        lambertian_flag = True
        libradunc_flag = True
    elif subcase=='curuc':
        adjacency_flag = False # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
        lambertian_flag = False
        libradunc_flag = False
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag, libradunc_flag, n_proc)
    lut.plot_results()
elif case == 'amazon':  # Amazon forest.
    path_l1c = '/home/larsshared/PROJECTS/L2A_RUT/S2products/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916.SAFE.zip'
    path_l2a = '/home/larsshared/PROJECTS/L2A_RUT/S2products/S2B_MSIL2A_20221107T144729_N0400_R139_T19MGP_20221107T183144.SAFE.zip'
    latlon = (-6.2869, -66.8652)
    roisize = (500, 500)
    # Definition of uncertainty setup based on subcase
    if subcase == 'standard':
        adjacency_flag = True
        lambertian_flag = True
        libradunc_flag = True
    elif subcase == 'curuc':
        adjacency_flag = False  # If adjacency, lambertian and libradtran flags set to False, CURUC uncertainty is added
        lambertian_flag = False
        libradunc_flag = False
    lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag, libradunc_flag, n_proc)
    lut.plot_results()

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
