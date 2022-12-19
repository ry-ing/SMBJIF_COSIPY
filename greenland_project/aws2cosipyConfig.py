"""
 This is the configuration (init) file for the utility aws2cosipy.
 Please make your changes here.
"""

#------------------------
# Declare variable names 
#------------------------

# Pressure
PRES_var = 'PSFC'

# Temperature
T2_var = 'T2'
in_K = False

# Relative humidity
RH2_var = 'RH2'

# Incoming shortwave radiation
G_var = 'SWDNB'

# Precipitation
RRR_var = 'PCPT'

# Wind velocity
U2_var = 'U2'

# Incoming longwave radiation
LWin_var = 'LWDNB'

# Snowfall
SNOWFALL_var = 'ACSNOW'

# Cloud cover fraction
N_var = 'N'

#------------------------
# Aggregation to hourly data
#------------------------
aggregate = False
aggregation_step = 'H'

# Delimiter in csv file
delimiter = ','

# WRF non uniform grid
WRF = False

#------------------------
# Radiation module 
#------------------------
radiationModule = 'Wohlfahrt2016' # 'Moelg2009', 'Wohlfahrt2016', 'none'
LUT = False                   # If there is already a Look-up-table for topographic shading and sky-view-factor built for this area, set to True

dtstep = 3600               # time step (s)
stationLat =  67.09           # Latitude of station
tcart = -7                   # Station time correction in hour angle units (1 is 4 min)
timezone_lon = -49.96	      # Longitude of station

# Zenit threshold (>threshold == zenit): maximum potential solar zenith angle during the whole year, specific for each location
zeni_thld = 43.8              # If you do not know the exact value for your location, set value to 89.0

#------------------------
# Point model 
#------------------------

#C161:  -134.42, 58.63, 1488
#DG1:   -134.13, 58.62, 1018
#TKG4:  -134.24, 58.63, 1116
#DG3:   -134.09, 58.72, 1352

point_model = False
plon = -49.96
plat = 67.09
hgt = 675.0

#------------------------
# Interpolation arguments 
#------------------------
stationName = 'KAN_L'
stationAlt = 675.0

lapse_T         = -0.0055    # Temp K per  m
lapse_RH        =  -0.000    # RH % per  m (0 to 1)
lapse_RRR       =  0.000   # mm per m
lapse_SNOWFALL  =  0.000   # Snowfall % per m (0 to 1)
