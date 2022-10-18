print('Importing modules')
import xarray as xr
import pandas as pd
import numpy as np
import netCDF4 as nc
import time
import dateutil
from calendar import isleap
import os
import sys

def add_variable_along_timelatlon_wrf(ds, var, name, units, long_name):
    """ This function adds missing variables to the DATA class """
    ds[name] = (('time','south_north','west_east'), var)	
    return ds
    
def add_variable_along_latlon_wrf(ds, var, name, units, long_name):
    """ This function self.adds missing variables to the self.DATA class """
    ds[name] = (('south_north','west_east'), var)
    return ds

#----------------------------------------------------------------------------------------------
# Read data              CHANGE DATE TO: 1980_1990
#------------------------------------------------------------------------------------------------
# Make sure to change year in all file names and TIMESTAMP CODE

period = '2001_2010'
start_date = '2001-01-01'
end_date = '2010-12-31'

data_dir = '/fastdata/ggp21rni/wrf_data/ccsm_4km/merged_files/%s/' % period #CHANGE DATE

#temporary files for processing
wrf_file = data_dir + '%s_merge_temp.nc' % period                           #CHANGE DATE
reformat_file = data_dir + 'reformat_temp.nc' 
grid_file = '/fastdata/ggp21rni/wrf_data/cfsr_4km/merged_files/grid.txt' #grid file to copy over to output file
temp_file = data_dir + 'regrid_temp.nc'

output_file = '/fastdata/ggp21rni/wrf_data/ccsm_4km/merged_files/ccsm_final/CCSM_%s_output.nc' % period #CHANGE DATE

coord_file = '/fastdata/ggp21rni/wrf_data/cfsr_4km/2019/WRFDS_2019-01-01.nc' # coord file for regridding coordinates


print('Loading in datasets with dask')
df_coord = xr.open_dataset(coord_file) #WRF file with the correct coordinates

# Read CFSR file
df = xr.open_dataset(wrf_file, chunks='auto') #.chunk({'time': 146}) #loading in using dask array chunks
print(df)
print('Generating time coordinates')
#---------Creating and assigning time coordinates----------------------------------------
timestamp = pd.date_range(start=start_date, end=end_date, freq='D')            # CHANGE DATE
timestamp = timestamp[(timestamp.day != 29) | (timestamp.month != 2)] # removing leap year dates

hour_array = np.full((len(timestamp.values), 1), 12)
dt = pd.DataFrame(hour_array, columns = ['Hour']) #creating and adding hour column

timestamp += pd.to_timedelta(dt.Hour, unit='h')
timestamp = pd.to_datetime(timestamp).dt.normalize() #reseting the hour time to midnight

df = df.rename({'Time':'Time'})
# Re-format timestamp (only hour and minutes, no seconds)
df = df.assign_coords(Time=pd.to_datetime(timestamp.values).strftime('%Y-%m-%d'))

print('Assigning coordinates')
# Create intermediate input file
dso = xr.Dataset(df_coord.attrs) #adding attributes
dso.coords['time'] = (('time'), df.Time.values)
dso = dso.assign_coords(time=pd.to_datetime(dso['time'].values)) #Needed for radiation module
dso.coords['XLAT'] = (('south_north', 'west_east'), df_coord.XLAT[:].values)
dso.coords['XLON'] = (('south_north', 'west_east'), df_coord.XLONG[:].values)
	    
# Calculating a rough elevation field needed for the downscaling of T2	    
p_0 = df.SLP[0].values
p = df.PSFC[0].values
p = p/100.0
T2 = df.T2[0].values
HGT = ( (((p_0/p)**(1/5.257)) - 1) * (T2) )/0.0065

print('Adding variables to the dataset')
dso = add_variable_along_latlon_wrf(dso, HGT,'HGT', 'm', 'Elevation')	 

# Add variables to file 
print('Adding T2 (2/16)')
dso = add_variable_along_timelatlon_wrf(dso, df.T2.values, 'T2', 'm', 'Temperature at 2 m')
dso = add_variable_along_timelatlon_wrf(dso, df.Q2.values, 'Q2', 'm', 'Specific humidity at 2m')
dso = add_variable_along_timelatlon_wrf(dso, df.PSFC.values, 'PSFC', 'hPa', 'Atmospheric Pressure')
print('Adding SNOW (6/16)')
dso = add_variable_along_timelatlon_wrf(dso, df.U10.values, 'U10', 'm s^-1', 'Wind velocity at 10 m')
dso = add_variable_along_timelatlon_wrf(dso, df.V10.values, 'V10', 'm s^-1', 'Wind velocity at 10 m')
print('Adding LWDNB (11/16)')
dso = add_variable_along_timelatlon_wrf(dso, df.LWDNB.values, 'LWDNB', 'W m\u207b\xb2', 'Incoming longwave radiation')
dso = add_variable_along_timelatlon_wrf(dso, df.SWDNB.values, 'SWDNB', 'W m\u207b\xb2', 'Incoming shortwave radiation')
print('Adding PCPT (14/16)')
dso = add_variable_along_timelatlon_wrf(dso, df.PCPT.values, 'PCPT', 'mm', 'Total precipitation')
dso = add_variable_along_timelatlon_wrf(dso, df.ACSNOW, 'ACSNOW', 'kg m^-2', 'Snowfall')

print(dso)  

dso.to_netcdf(reformat_file)

print('Performing CDO grid routines')
os.system('cdo setgrid,' + grid_file + ' ' + reformat_file + ' ' + temp_file)
os.system('cdo sellonlatbox,-135.5,-133.3,58.325,59.75 ' + temp_file + ' ' + output_file)


print('\n Input file created \n')
print('-------------------------------------------')
