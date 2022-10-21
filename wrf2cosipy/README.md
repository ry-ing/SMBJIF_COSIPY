  This folder contains the code to generate the COSIPY input data from the raw WRF data.
  
  wrf_cfsr_prep.py and wrf_gfdl-ccsm_prep.py takes the raw WRF data (with corrupt coordinates), subsets to the JIF region and fixes any time and space coordinate issues.
  
  wrf2cosipyinput.py takes the fixed WRF data and generates the main COSIPY input files. This includes the downscaling and interpolation to the previously generated 600m static data grid.
  
  wrf_bias_corr.py takes the COSIPY input files and bias corrects the future projections.
