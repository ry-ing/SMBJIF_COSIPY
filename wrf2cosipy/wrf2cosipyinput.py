
print('Importing modules')
import sys
import xarray as xr
import pandas as pd
import numpy as np
import netCDF4 as nc
import time
import dateutil
from itertools import product
import metpy.calc
from metpy.units import units
from numba import jit
import os
import sys

sys.path.append('../../')

from utilities.aws2cosipy.aws2cosipyConfig import *
from cosipy.modules.radCor import *

import argparse

def create_2D_input(wrf_file, cosipy_file, static_file, start_date, end_date, x0=None, x1=None, y0=None, y1=None):
    """ This function creates an input for COSIPY by interpolating fields from 
        the WRF data to the grid of the static file (aggregated DEM).
	"""

    print('-------------------------------------------')
    print('Remapping select variables using CDO remapbil \n')
   
    temp_file1 = '/home/ggp21rni/cosipy/utilities/aws2cosipy/temp_file1.nc'
    temp_file2 = '/home/ggp21rni/cosipy/utilities/aws2cosipy/temp_file2.nc'
    
    #os.system('rm ' + temp_file1 + ' ' + temp_file2)
    
    os.system('cdo selvar,Q2,SWDNB,LWDNB,PCPT,V10,U10,ACSNOW ' + wrf_file + ' ' + temp_file1)
    os.system('cdo remapbil,' + static_file + ' ' + temp_file1 + ' ' + temp_file2)
    remap_ds = xr.open_dataset(temp_file2)

    #----------------------------------------------------------------------------------------------
    # Loading in WRF file
    #------------------------------------------------------------------------------------------------
    print('Read input file %s' % (wrf_file))
    # Read WRF file
    df = xr.open_dataset(wrf_file).chunk({"time":20})
    
    # Rename the time coordinate
    df = df.rename({'time':'Time'})
    
    # Re-format timestamp (only hour and minutes, no seconds)
    df = df.assign_coords(Time=pd.to_datetime(df['Time'].values).strftime('%Y-%m-%d'))

    # Select the specified period
    if ((start_date!=None) & (end_date!=None)):
        df = df.sel(Time=slice(start_date,end_date))
   
    #----------------------------------------------------------------------------------------------
    # Creating intermediate file for downscaling of T2
    #------------------------------------------------------------------------------------------------
    print('\n Creating intermediate file for downscaling of T2')         
	    
    # Create intermediate input file
    dso = xr.Dataset()
    dso.coords['time'] = (('time'), df.Time.values)
    dso = dso.assign_coords(time=pd.to_datetime(dso['time'].values)) #Needed for radiation module
    dso.coords['lat'] = (('south_north', 'west_east'), df.XLAT[:].values)
    dso.coords['lon'] = (('south_north', 'west_east'), df.XLONG[:].values)
    print(dso.time.values)
	      	   
    # Add variables to file 
    dso = add_variable_along_latlon_wrf(dso, df.HGT.values, 'HGT', 'm', 'Elevation')
    dso = add_variable_along_timelatlon_wrf(dso, df.T2.values, 'T2', 'm', 'Temperature at 2 m')
    dso = add_variable_along_timelatlon_wrf(dso, df.PSFC.values/100.0, 'PRES', 'hPa', 'Atmospheric Pressure')


    # Wind velocity at 2 m (assuming neutral stratification)
    z  = 10.0     # Height of measurement
    z0 = 0.00212 #0.0040 # Roughness length for momentum
    umag = np.sqrt(remap_ds.V10.values**2+remap_ds.U10.values**2)   # Mean wind velocity
    U2 = umag * (np.log(2 / z0) / np.log(10 / z0))

    #----------------------------------------------------------------------------------------------
    # Load static data
    #----------------------------------------------------------------------------------------------
    # Static file defines the grid values are interpolated to and the glacier mask
    
    print('\n Read static file %s \n' % (static_file))
    ds = xr.open_dataset(static_file).chunk({"lat": 100, "lon": 100})

    #-----------------------------------
    # Create subset
    #-----------------------------------
    ds = ds.sel(lat=slice(y0,y1), lon=slice(x0,x1))
    ds.coords['time'] = (('time'), dso.time.values) #set time coords for static file

    #-----------------------------------
    # Get values from file
    #-----------------------------------
    PRES = dso.PRES.values     # Pressure
    T2 = dso.T2.values         # Temperature
    mask = ds.MASK.values
    
    #-----------------------------------
    # Create numpy arrays for the 2D fields
    #-----------------------------------
    T_interp = np.full([len(ds.time), len(ds.lat), len(ds.lon)], np.nan)
    P_interp = np.full([len(ds.time), len(ds.lat), len(ds.lon)], np.nan)
    G_interp = np.full([len(ds.time), len(ds.lat), len(ds.lon)], np.nan)

    #-----------------------------------
    # Interpolate point data to grid 
    #-----------------------------------
    print('Interpolate wrf file to finer grid')
    
    # Interpolate data (T, RH, RRR, U)  to grid using lapse rates
    # Xarray datasets set to numpy arrays to be run using Numba (JIT)
    ds_lat = ds.lat.values
    ds_lon = ds.lon.values
    len_time = len(ds.time.values)
    len_lat = len(ds.lat.values)
    len_lon = len(ds.lon.values)
    HGT = ds.HGT.values
    
    dso_lat = dso.lat.values
    dso_lon = dso.lon.values
    hgt_wrf = dso.HGT.values
      
    # Function to interpolate wrf data to the aggregated DEM grid
    # Function run through Numba (JIT) to speed up the process  
    @jit(nopython=True)
    def interpolate_loop(T_interp, P_interp, ds_lat, hgt_wrf, ds_lon, dso_lat, dso_lon, HGT, mask):
    
      for t in range(len_time):
          print('TIMESTEP:', t, '/', len_time)
        
          for i in range(len_lat):            
              for j in range(len_lon):
                  if (mask[i,j] == 1):
		     
	          	#Finding nearest lat,lon coords from WRF data to the new finer grid
                  	abslat = np.abs(dso_lat - ds_lat[i])
                  	abslon = np.abs(dso_lon - ds_lon[j])
                  	c = np.maximum(abslon, abslat)
                  	([xloc], [yloc]) = np.where(c == np.min(c)) #x,y index for coords

                  	T_interp[t,i,j] = (T2[t,xloc,yloc]) + (HGT[i,j] - hgt_wrf[xloc,yloc])*lapse_T		  
		  
                  	# Interpolate pressure using the barometric equation
                  	SLP = PRES[t,xloc,yloc]/np.power((1-(0.0065*hgt_wrf[xloc,yloc])/(288.15)), 5.255)
                  	P_interp[t,i,j] = SLP * np.power((1-(0.0065*HGT[i,j])/(288.15)), 5.255)   
	     	         
		  
    interpolate_loop(T_interp, P_interp, ds_lat, hgt_wrf, ds_lon, dso_lat, dso_lon, HGT, mask)
    
    print(('\n Number of grid cells: %i') % (np.count_nonzero(~np.isnan(ds['MASK'].values))))
    print(('\n Number of glacier cells: %i') % (np.nansum(ds['MASK'].values)))

    # Auxiliary variables
    
    hgt = ds.HGT.values
    slope = ds.SLOPE.values
    aspect = ds.ASPECT.values
    lats = ds.lat.values
    lons = ds.lon.values
    sw = remap_ds.SWDNB.values

    #-----------------------------------
    # Run radiation module 
    #-----------------------------------
    if radiationModule == 'Wohlfahrt2016':
        print('Run the Radiation Module: Wohlfahrt2016')

        # Change aspect to south==0, east==negative, west==positive
        aspect = ds['ASPECT'].values - 180.0
        ds['ASPECT'] = (('lat', 'lon'), aspect)
	
        doy = np.zeros(len_time)
        hour = np.zeros(len_time)
        for t in range(len_time):
           doy[t] = dso.time[t].dt.dayofyear
           hour[t] = dso.time[t].dt.hour
  
        @jit
        def radiation_module(G_interp, sw, lats, lons, slope, aspect, doy, hour):
          for t in range(len_time):
              #doy = dso.time[t].dt.dayofyear
              #hour = dso.time[t].dt.hour
              print('TIMESTEP RADIATION:', t, '/', len_time)
	    
              for i in range(len_lat):
                  for j in range(len_lon):
                      if (mask[i, j] == 1):
                          if radiationModule == 'Wohlfahrt2016':			    
                              G_interp[t, i, j] = np.maximum(0.0, correctRadiation(lats[i], lons[j], timezone_lon, doy[t], hour[t], slope[i, j], aspect[i, j], sw[t,i,j], zeni_thld))
                          else:
                              G_interp[t, i, j] = sw[t,i,j]
			    
        radiation_module(G_interp, sw, lats, lons, slope, aspect, doy, hour)
    #-----------------------------------
    # Add variables to file 
    #-----------------------------------
    
    add_variable_along_timelatlon(ds, T_interp, 'T2', 'K', 'Temperature at 2 m')
    add_variable_along_timelatlon(ds, wrf_rh(T_interp, remap_ds.Q2.values, P_interp*100.0), 'RH2', '%', 'Relative humidity at 2 m')
    add_variable_along_timelatlon(ds, U2, 'U2', 'm s\u207b\xb9', 'Wind velocity at 2 m')
    add_variable_along_timelatlon(ds, P_interp, 'PRES', 'hPa', 'Atmospheric Pressure')
    add_variable_along_timelatlon(ds, remap_ds.PCPT.values, 'RRR', 'mm', 'Total precipitation (liquid+solid)')
    add_variable_along_timelatlon(ds, remap_ds.ACSNOW.values/1000.0 , 'SNOWFALL', 'm', 'Snowfall')
    add_variable_along_timelatlon(ds, remap_ds.LWDNB.values, 'LWin', 'W m\u207b\xb2', 'Incoming longwave radiation') 
    add_variable_along_timelatlon(ds, G_interp, 'G', 'W m\u207b\xb2', 'Incoming shortwave radiation')

    #-----------------------------------
    # Write file to disc 
    #-----------------------------------
    check_for_nan(ds)
    ds.to_netcdf(cosipy_file)

    print('\n Input file created \n')
    print('-------------------------------------------')
    
    os.system('rm ' + temp_file1 + ' ' + temp_file2)

#----------- Additional functions for adding variables to XARRAY dataset ----------------------
def add_variable_along_timelatlon_wrf(ds, var, name, units, long_name): #Add variable along time and latlon dimensions for WRF data
    """ This function adds missing variables to the DATA class """
    ds[name] = (('time','south_north','west_east'), var)	
    return ds
    
def add_variable_along_timelatlon(ds, var, name, units, long_name): #Add variable along time and latlon dimensions
    """ This function adds missing variables to the DATA class """
    ds[name] = (('time','lat','lon'), var)
    ds[name].attrs['units'] = units
    ds[name].attrs['long_name'] = long_name
    return ds

def add_variable_along_latlon_wrf(ds, var, name, units, long_name): #Add variable along JUST latlon dimensions for WRF data
    """ This function self.adds missing variables to the self.DATA class """
    ds[name] = (('south_north','west_east'), var)
    return ds
    
def add_variable_along_latlon(ds, var, name, units, long_name): #Add variable along JUST latlon dimensions
    """ This function self.adds missing variables to the self.DATA class """
    ds[name] = (('lat','lon'), var)
    ds[name].attrs['units'] = units
    ds[name].attrs['long_name'] = long_name
    ds[name].encoding['_FillValue'] = -9999
    return ds
    
def wrf_rh(T2, Q2, PSFC): #For calculating relative humidity 
    pq0 = 379.90516
    a2 = 17.2693882
    a3 = 273.16
    a4 = 35.86
    rh = Q2 * 100 / ( (pq0 / PSFC) * np.exp(a2 * (T2 - a3) / (T2 - a4)) )
    
    rh[rh>100.0] = 100.0
    rh[rh<0.0] = 0.0
    return rh
    
def check_for_nan(ds): #Checking for NaNs in dataset
        for y,x in product(range(ds.dims['lat']),range(ds.dims['lon'])):
            mask = ds.MASK.isel(lat=y, lon=x)
            if mask==1:
                if np.isnan(ds.isel(lat=y, lon=x).to_array()).any():
                    print('ERROR!!!!!!!!!!! There are NaNs in the dataset')

#Defining files, dates and running function
wrf_file = '/fastdata/ggp21rni/wrf_data/ccsm_4km/merged_files/ccsm_final/WRF_CCSM_1980_2010_final_output.nc' #CFSR_1980-2012_reformat.nc'
cosipy_file ='/fastdata/ggp21rni/cosipy_input/cosipy_input_CCSM_JIF_1980_2010.nc' 
static_file = '/data/ggp21rni/static_data/JIF_static.nc' #JIF_static, lemoncreek_static, taku_static
start_date ='1980-01-01'
end_date='2010-12-31'

create_2D_input(wrf_file, cosipy_file, static_file, start_date, end_date)
