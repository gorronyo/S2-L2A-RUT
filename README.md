# S2-L2A-RUT
The software generates radiometric uncertainty estimates and spectral error correlation
for the Copernicus Sentinel-2 L2A products (i.e. surface reflectance) and L1C products (i.e. TOA reflectance).

| ![L2A-RUT workflow](./L2Aunc_results//l2arut_workflow.png) |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|                                         <b>Simplified workflow of the L2A-RUT software.</b>                                          |

The software first defines a class named *L2A_unc* and calls a function named *get_libradunc*.
This function reads the L1C and L2A data products to extract all possible metadata that defined the Libradtran parameterisation (e.g. estimated value of AOT).
Furthermore, it also extracts the L1C reflectance, converts to radiance and generates a set of TOA radiance distribution values.
Once all this information has been obtained, the Libradtran input can be defined, run the radiative transfer code and extract the atmospheric parameters.
Finally, the surface reflectance for all bands is calculated and other effects such as adjacency are also included.
This process is repeated (defined by the user value) and results in a set of surface reflectance samples that define the surface reflectance distribution for all bands out of which a covariance matrix can be obtained.

From the output, here we offer the uncertainty (*k*=1) based on the standard deviation of the surface reflectance samples in a normal distribution.
If the distribution is considerably asymmetric, the reference [GUM-MCM2008](https://www.bipm.org/documents/20126/2071204/JCGM_101_2008_E.pdf/325dcaad-c15a-407c-1105-8b7f322d651c) recommends the use of the shortest coverage interval that corresponds to a specific coverage factor (e.g. 68.27% for *k*=1).


## Get the code
You will need to have installed a [Git client](https://git-scm.com) for fetching the source code

Clone or fork the repository at [GitHub](https://github.com/senbox-org/snap-rut)
```
> git clone https://github.com/gorronyo/S2-L2A-RUT.git
> cd snap-rut
```
You can update your checked-out sources from the remote repository by running 
```
> git pull --rebase
```

## Requirements
The code uses *common* libraries from Python. You will need to have installed those listed in the *requirements.tx*.
You can directly installed them as follows:
```
> pip install -r requirements.txt
```

You will also need to install two separate software programmes:
- [SNAP](https://step.esa.int/main/download/snap-download/) and configure the [Python snappy wrapper](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface) during (or after) the installation. Recommended SNAP latest version 9.0.0.
- [Libradtran](http://www.libradtran.org/doku.php). Recommended latest version 2.0.5.

The code has been tested with Python 3.6. Higher versions of Python might not work correctly with the Python snappy wrapper (see [here](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface)).

## Running the code
The code includes a script *run_L2a_unc.py* that wraps the processing workflow and runs several examples.
The lines 3-20 in this script are the only ones of interest to the user:

```python class:"lineNo"
3   samples = 1000
4   n_proc = 20  # check max processes with multiprocessing.cpu_count()
5   librad_bin = None # path to folder with libradtran binary. E.g. '/home/gorrono/libRadtran-2.0.4'
6   lut = L2a_unc.L2aUnc(librad_bin, samples)
7   case = 'user'  # 'user', 'libya4' , 'amazon', 'winterwheat' or 'maize'
8   subcase = 'standard' # 'standard', 'curuc' and 'reptrantest' options for cases 'libya4' and 'amazon'
9 
10  if case == 'user':  # user defined case.
11      path_l1c = None # FULLPATH_TO_L1C.zip'
12      path_l2a = None # FULLPATH_TO_L2A.zip'
13      latlon = (None, None) # LATITUDE AND LONGITUDE
14      roisize = (None, None) # AREA SIZE IN METERS
15      adjacency_flag = True
16      lambertian_flag = True
17      libradunc_flag = True
18      lut.get_libradunc(path_l1c, path_l2a, latlon, roisize, adjacency_flag, lambertian_flag,
19                        libradunc_flag)
20      lut.plot_results()
```
The variable ```samples``` in line 3 must be set to a sufficient number for convergence (1000 is default) whereas the user must define
the location of the libradtran binary ```librad_bin``` in line 5.
The number of parallel processes is dedined with ```nproc``` in line 4. Make sure how make processes your system can deal with. The python command ```multiprocessing.cpu_count()``` can tell you how many processes you have and/or directly set to this value.

You must download the zip products for both L1C and L2A products and define the full path in lines 13-16 (```path_L1C``` and ```path_L2A```).
The Sentinel-2 products must be from [Collection 1 reprocessing](https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-2-msi/copernicus-sentinel-2-collection-1-availability-status). 

The software calculates the uncertainty and spectral correlation at a specific location defined in ```latitude``` and ```longitude``` for an area size defined in ```roisize```. If the user sets ```roisize``` at or below 10 meters, the calculation will be generated for a single pixel. Otherwise, it will be performed for
the mean value of the rectangular area defined by ```roisize```.

Now you can simply run/import the script! 🚀

## Results
All the results will be available in class variable ```lut```. For example ```lut.L2Arho``` will contain an array
with the 12 L2A bands and, for each one of them, the distribution samples.

All the results will be automatically plotted and stored in the folder *./L2Aunc_results*.
If you prefer not to obtain these plots, please comment line 23!. 


| ![B5 surface reflectance_Amazonforest](./L2Aunc_results/amazon_reptranband1000samp/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916_hist_L2Arho_B5.png) |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|                                         <b>Surface reflectance distribution for Sentinel-2 L2A B5 at the Amazon forest</b>                                          |

| ![Spectral_correlation_Amazonforest](./L2Aunc_results/amazon_reptranband1000samp/S2B_MSIL1C_20221107T144729_N0400_R139_T19MGP_20221107T174916_speccorrL2Arho.png) |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|                                          <b>Spectral error correlation for Sentinel-2 L2A bands at the Amazon forest</b>                                          |



