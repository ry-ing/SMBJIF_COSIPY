#from COSIPY import model_file

import pickle

model_file_file = open('model_file_file', 'rb')
model_file = pickle.load(model_file_file)
model_file_file.close()

model_input_file_file = open('model_input_file_file', 'rb')
model_input_file = pickle.load(model_input_file_file)
model_input_file_file.close()

"""
 This is the COSIPY configuration (init) file.
 Please make your changes here.
"""
#-----------------------------------
# SIMULATION PERIOD 
#-----------------------------------
# Taku test run
time_start = '2003-07-01'
time_end   = '2008-12-30'

#-----------------------------------
# FILENAMES AND PATHS 
#-----------------------------------
time_start_str=(time_start[0:10]).replace('-','')
time_end_str=(time_end[0:10]).replace('-','')

data_path = './data/'

# lapse rate of -5oc per km
input_netcdf= model_input_file #'/fastdata/ggp21rni/cosipy_input/cosipy_input_GFDLdownscaled_LC-C_2000-2010_v1.nc' 
output_netcdf = model_file

#-----------------------------------
# RESTART 
#-----------------------------------
restart = False                                             # set to true if you want to start from restart file

#-----------------------------------
# STAKE DATA 
#-----------------------------------
stake_evaluation = False 
stakes_loc_file = '/data/ggp21rni/cosipy_input/taku_stakes_loc_sharc_v2.csv'         # path to stake location file
stakes_data_file = '/data/ggp21rni/cosipy_input/taku_stakes_data_sharc_v3.csv'   # path to stake data file
eval_method = 'rmse'                                        # how to evaluate the simulations ('rmse')
obs_type = 'mb'                                     # What kind of stake data is used 'mb' or 'snowheight'

#-----------------------------------
# STANDARD LAT/LON or WRF INPUT 
#-----------------------------------
# Dimensions
WRF = False                                                 # Set to True if you use WRF as input

northing = 'lat'	                                    # name of dimension	in in- and -output
easting = 'lon'					                        # name of dimension in in- and -output
if WRF:
    northing = 'south_north'                                # name of dimension in WRF in- and output
    easting = 'west_east'                                   # name of dimension in WRF in- and output

# Interactive simulation with WRF
WRF_X_CSPY = False

#-----------------------------------
# COMPRESSION of output netCDF
#-----------------------------------
compression_level = 2                                       # Choose value between 1 and 9 (highest compression)
                                                            # Recommendation: choose 1, 2 or 3 (higher not worthwhile, because of needed time for writing output)
#-----------------------------------
# PARALLELIZATION 
#-----------------------------------
slurm_use = False                                           # use SLURM
workers = None                                              # number of workers, if local cluster is used
local_port = 8786                                           # port for local cluster

#-----------------------------------
# WRITE FULL FIELDS 
#-----------------------------------    
full_field = False                                          # write full fields (2D data) to file
if WRF_X_CSPY:
    full_field = True
    
#-----------------------------------
# TOTAL PRECIPITATION  
#-----------------------------------
force_use_TP = True                                        # If total precipitation and snowfall in input data;
                                                            # use total precipitation

#-----------------------------------
# CLOUD COVER FRACTION  
#-----------------------------------
force_use_N = False #gets lwc_melted_layers index out of range error              # If cloud cover fraction and incoming longwave radiation
                                                                       # in input data use cloud cover fraction

#-----------------------------------
# SUBSET  (provide pixel values) 
#-----------------------------------
tile = False
xstart = 20
xend = 40
ystart = 20
yend = 40
