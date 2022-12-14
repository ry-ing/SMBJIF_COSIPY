#from COSIPY import albsnow, albfirn, albice, albtime, albdepth, rsnow, rfirn, rice, surf_em

import pickle
import numpy as np

param_file = open('param_values', 'rb')
param_array = pickle.load(param_file)
param_file.close()

param_array = param_array.astype(np.float)
param_array = np.round(param_array, 2)

albsnow = param_array[0]
albfirn = param_array[1]
albice = param_array[2]
albtime = param_array[3]
albdepth = param_array[4]
rsnow = param_array[5]
rfirn = param_array[6]
rice = param_array[7]
surf_em = param_array[8]
mult_rrr = param_array[9]

WRF_X_CSPY=False 

"""
    Declaration of constants
    Do not modify unless you are absolutely sure what you are doing.
"""

' GENERAL INFORMATION ' 
dt = 86400                                       # Time step in the input files [s]
max_layers = 300    #200                            # Max. number of layers, just for the restart file
z = 2.0                                         # Measurement height [m]

' PARAMETERIZATIONS '
stability_correction = 'Ri'                     # possibilities: 'Ri','MO'
albedo_method = 'Oerlemans98'                   # possibilities: 'Oerlemans98'
densification_method = 'Boone'                  # possibilities: 'Boone','empirical','constant' TODO: solve error Vionnet
penetrating_method = 'Bintanja95'               # possibilities: 'Bintanja95'
roughness_method = 'Moelg12'                    # possibilities: 'Moelg12'
saturation_water_vapour_method = 'Sonntag90'    # possibilities: 'Sonntag90'
thermal_conductivity_method = 'bulk'		    # possibilities: 'bulk', 'empirical'
sfc_temperature_method = 'SLSQP'                # possibilities: 'L-BFGS-B', 'SLSQP'(faster), 'Newton' (Secant, fastest)'

# WRF_X_CSPY: for efficiency and consistency
if WRF_X_CSPY:
    stability_correction = 'MO'
    sfc_temperature_method = 'Newton'


' INITIAL CONDITIONS '
initial_snowheight_constant = 0.2               # Initial snowheight
initial_snow_layer_heights = 0.10               # Initial thickness of snow layers
initial_glacier_height = 100.0                   # Initial glacier height without snowlayers ORIGINALLY 40
initial_glacier_layer_heights = 0.5             # Initial thickness of glacier ice layers

initial_top_density_snowpack = 300.0            # Top density for initial snowpack
initial_bottom_density_snowpack = 600.0         # Bottom density for initial snowpack

temperature_bottom = 270.16 #270.16                     # Lower boundary condition for initial temperature profile (K)
const_init_temp = 0.1                           # constant for init temperature profile used in exponential function (exponential decay)

zlt1 = 0.06					                    # First depth for temperature interpolation which is used for calculation of ground heat flux
zlt2 = 0.1					                    # Second depth for temperature interpolation which is used for calculation of ground heat flux

' MODEL CONSTANTS '
center_snow_transfer_function = 1.0             # center (50/50) when total precipitation is transferred to snow and rain
spread_snow_transfer_function = 1               # 1: +-2.5
mult_factor_RRR = mult_rrr #1.0                           # multiplication factor for RRR

minimum_snow_layer_height = 0.001               # minimum layer height
minimum_snowfall = 0.001                        # minimum snowfall per time step in m which is added as new snow


' REMESHING OPTIONS'
remesh_method = 'log_profile'                   # Remeshing (log_profile or adaptive_profile)
first_layer_height = 0.01                       # The first layer will always have the defined height (m)
layer_stretching = 1.20                         # Stretching factor used by the log_profile method (e.g. 1.1 mean the subsequent layer is 10% greater than the previous

merge_max = 1                                   # How many mergings are allowed per time step
density_threshold_merging = 5                   # If merging is true threshold for layer densities difference two layer try: 5-10 (kg m^-3)
temperature_threshold_merging = 0.01            # If mering is true threshold for layer temperatures to merge  try: 0.05-0.1 (K)


' PHYSICAL CONSTANTS '
constant_density = 300.                         # constant density of freshly fallen snow [kg m-3], if densification_method is set to 'constant'

albedo_fresh_snow = albsnow                    # albedo of fresh snow [-] (Moelg et al. 2012, TC)
albedo_firn = albfirn                           # albedo of firn [-] (Moelg et al. 2012, TC)
albedo_ice = albice                            # albedo of ice [-] (Moelg et al. 2012, TC)
albedo_mod_snow_aging = albtime               # effect of ageing on snow albedo [days] (Oerlemans and Knap 1998, J. Glaciol.)
albedo_mod_snow_depth = albdepth                 # effect of snow depth on albedo [cm] (Oerlemans and Knap 1998, J. Glaciol.)

### For tropical glaciers or High Mountain Asia summer-accumulation glaciers (low latitude), the Moelg et al. 2012, TC should be tested for a possible better albedo fit 
#albedo_mod_snow_aging = 6                      # effect of ageing on snow albedo [days] (Moelg et al. 2012, TC)
#albedo_mod_snow_depth = 8                      # effect of snow depth on albedo [cm] (Moelg et al. 2012, TC)

roughness_fresh_snow = rsnow               # surface roughness length for fresh snow [mm] (Moelg et al. 2012, TC)
roughness_firn = rfirn                      # surface roughness length for aged snow [mm] (Moelg et al. 2012, TC)
roughness_ice = rice               # surface roughness length for ice [mm] (Moelg et al. 2012, TC)
aging_factor_roughness = 0.0026                 # effect of ageing on roughness lenght (hours) 60 days from 0.24 to 4.0 => 0.0026

snow_ice_threshold = 900.0                      # pore close of density [kg m^(-3)]

lat_heat_melting = 3.34e5                       # latent heat for melting [J kg-1]
lat_heat_vaporize = 2.5e6                       # latent heat for vaporization [J kg-1]
lat_heat_sublimation = 2.834e6                  # latent heat for sublimation [J kg-1]

spec_heat_air = 1004.67                         # specific heat of air [J kg-1 K-1]
spec_heat_ice = 2050.00                         # specific heat of ice [J Kg-1 K-1]
spec_heat_water = 4217.00                       # specific heat of water [J Kg-1 K-1]

k_i = 2.22                                      # thermal conductivity ice [W m^-1 K^-1]
k_w = 0.55                                      # thermal conductivity water [W m^-1 K^-1]
k_a = 0.024                                     # thermal conductivity air [W m^-1 K^-1]

water_density = 1000.0                          # density of water [kg m^(-3)]
ice_density = 917.                              # density of ice [kg m^(-3)]
air_density = 1.1                               # density of air [kg m^(-3)]

sigma = 5.67e-8                                 # Stefan-Bolzmann constant [W m-2 K-4]
zero_temperature = 273.16                       # Melting temperature [K]
surface_emission_coeff = surf_em                  # surface emission coefficient [-]
