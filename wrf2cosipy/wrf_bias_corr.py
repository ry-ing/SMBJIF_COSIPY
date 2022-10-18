import sys
import xarray as xr
import pandas as pd
import numpy as np
import sys
import scipy
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numba_scipy
from numba import jit

def bias_correction(obs, p, s, method='delta', nbins=10, extrapolate=None):
    """Bias Correction techniques for correcting simulated output based on differences between the CDFs of
    observed and simulated output for a training period.

    three different methods are available
    'delta'   This is the simplest bias correction method, which consists on adding the mean change signal
              to the observations (delta method). This method corresponds to case g=1 and f=0 in Amengual
              et al. (2012). This method is applicable to any kind of variable but it is preferable not to
              apply it to bounded variables (e.g. precipitation, wind speed, etc.) because values out of
              range could be obtained.
    'scaling' This method is very similar to the delta method but, in this case, the correction consist on
              scaling the simulation with the difference (additive: 'scaling_add') or quotient
              (multiplicative: 'scaling_multi') between the mean of the observations and the simulation in
              the train period.
    'eqm'     Empirical Quantile Mapping (eQM) This is the most popular bias correction method which consists
              on calibrating the simulated Cumulative Distribution Function (CDF) by adding to the observed
              quantiles both the mean delta change and the individual delta changes in the corresponding
              quantiles. This is equivalent to f=g=1 in Amengual et al. (2012). This method is applicable to
              any kind of variable.

    input args
    obs:      observed climate data for the training period
    p:        simulated climate by the model for the same variable obs for the training period.
    s:        simulated climate for the variables used in p, but considering the test/projection period.
    method:   'delta', 'scaling_add', 'scaling_multi', 'eqm', see explenation above
    nbins:    for 'eqm' method only: number of quantile bins in case of 'eqm' method (default = 10)
    extrapolate: for 'eqm' method only: None (default) or 'constant' indicating the extrapolation method to
              be applied to correct values in 's' that are out of the range of lowest and highest quantile of 'p'

    output
    c:        bias corrected series for s


    ref:
    Amengual, A., Homar, V., Romero, R., Alonso, S., & Ramis, C. (2012). A statistical adjustment of regional
    climate model outputs to local scales: application to Platja de Palma, Spain. Journal of Climate, 25(3), 939-957.
    http://journals.ametsoc.org/doi/pdf/10.1175/JCLI-D-10-05024.1

    """

    if (method == 'eqm') and (nbins > 1):
        binmid = np.arange((1./nbins)*0.5, 1., 1./nbins)
        qo = mquantiles(obs[np.isfinite(obs)], prob=binmid)
        qp = mquantiles(p[np.isfinite(p)], prob=binmid)
        p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
        c = p2o(s)
        if extrapolate is None:
            c[s > np.max(qp)] = qo[-1]
            c[s < np.min(qp)] = qo[0]
        elif extrapolate == 'constant':
            c[s > np.max(qp)] = s[s > np.max(qp)] + qo[-1] - qp[-1]
            c[s < np.min(qp)] = s[s < np.min(qp)] + qo[0] - qp[0]

    elif method == 'delta':
        c = obs + (np.nanmean(s) - np.nanmean(p))

    elif method == 'scaling_add':
        c = s - np.nanmean(p) + np.nanmean(obs)

    elif method == 'scaling_multi':
        c = (s/np.nanmean(p)) * np.nanmean(obs)

    else:
        raise ValueError("incorrect method, choose from 'delta', 'scaling_add', 'scaling_multi' or 'eqm'")

    return c

def COSIPY_projection_file(obs_file, hist_file, proj_file, cosipy_file):
    """ This function creates an input for COSIPY by interpolating fields from 
        the WRF data to the grid of the static file (aggregated DEM).
	"""
    #----------------------------------------------------------------------------------------------
    # Loading in WRF file
    #------------------------------------------------------------------------------------------------
    print('Reading input files %s,\n%s,\n%s' % (obs_file, hist_file, proj_file))
    # Read WRF file
    CFSR_DATA = xr.open_dataset(obs_file, chunks='auto') #acts as observations
    HIST_DATA = xr.open_dataset(hist_file, chunks='auto')
    PROJ_DATA = xr.open_dataset(proj_file, chunks='auto')
    BC_DATA = PROJ_DATA.copy()
    
    mask = CFSR_DATA.MASK.values
    len_lat = len(HIST_DATA.lat.values)
    len_lon = len(HIST_DATA.lon.values)
    len_time = len(HIST_DATA.time.values)
    fields = ['T2', 'RH2', 'PRES', 'U2', 'LWin', 'G', 'RRR', 'SNOWFALL']
    
    for field in fields:
        print('\n---------------\n')
        print('Bias correcting %s:' % field)
        field_data = np.full((len_time, len_lat, len_lon), np.nan)
	
        print('HIST DATA, MIN IS:', HIST_DATA[field].min().values, 'MAX IS:', HIST_DATA[field].max().values)
        print('PROJ DATA, MIN IS:', PROJ_DATA[field].min().values, 'MAX IS:', PROJ_DATA[field].max().values)
        print('CFSR DATA, MIN IS:', CFSR_DATA[field].min().values, 'MAX IS:', CFSR_DATA[field].max().values)

        for i in range(len_lat):
           print(i, '/', len_lat)            
           for j in range(len_lon):
               if (mask[i,j] == 1):
                  #print('bias correcting for single point')
                  CFSR_DATA_var = CFSR_DATA[field][0:len_time,i,j].values
                  HIST_DATA_var = HIST_DATA[field][:,i,j].values
                  PROJ_DATA_var = PROJ_DATA[field][:,i,j].values
			
                  #@jit
                  def numba_routine(CFSR_DATA_var, HIST_DATA_var, PROJ_DATA_var):
                      #print('Performing EQM routines')
                      corrected_field = bias_correction(CFSR_DATA_var, HIST_DATA_var, PROJ_DATA_var, method='eqm', extrapolate='constant', nbins=100)
                  
                      if field in ['RRR', 'LWin', 'G', 'U2']:
                            corrected_field[corrected_field < 0] = 0
		
                      if field in ['RH2']:
                            corrected_field[corrected_field < 0] = 0
                            corrected_field[corrected_field > 100] = 100
			
                      field_data[:,i,j] = corrected_field
		      
                  numba_routine(CFSR_DATA_var, HIST_DATA_var, PROJ_DATA_var)
		  
                  #x_series = np.linspace(0, len_time, len_time)
                  #plt.plot(x_series, CFSR_DATA_var, label='OBS')
                  #plt.plot(x_series, HIST_DATA_var, label='HIST')
                  #plt.plot(x_series, PROJ_DATA_var, label='PROJ')
                  #plt.plot(x_series, corrected_field, label='CORR')
                  #plt.legend()
                  #plt.savefig('bias_corr_test_%s.png' % field)
                  #plt.close()
		  		  
                  #obs_data = np.sort(CFSR_DATA_var.flatten())	
                  #hist_data = np.sort(HIST_DATA_var.flatten())
                  #proj_data = np.sort(PROJ_DATA_var.flatten())
                  #corr_data = np.sort(corrected_field.flatten())
		  	  
                  #obs = scipy.stats.norm.cdf(obs_data)
                  #hist = scipy.stats.norm.cdf(hist_data) 
                  #proj = scipy.stats.norm.cdf(proj_data) 
                  #corr = scipy.stats.norm.cdf(corr_data) 
		  
                  #plt.plot(obs_data, obs, label='OBS') 
                  #plt.plot(hist_data, hist, label='HIST')
                  #plt.plot(proj_data, proj, label='PROJ') 
                  #plt.plot(corr_data, corr, label='CORR')
                  #plt.xlim(0,10)
                  #plt.xscale('log')
                  #plt.legend()
                  #plt.savefig('Q_Q_plot_%s.png' % field)  
                  #sys.exit()
		
        BC_DATA[field] = (('time', 'lat', 'lon'), field_data)
        print('Bias corrected:', field, ':\n', BC_DATA[field], '\n')
	    
    #-----------------------------------
    # Write file to disc 
    #-----------------------------------
    #check_for_nan(BC_DATA)
    BC_DATA.to_netcdf(cosipy_file)
    print(BC_DATA)

    print('\n Input file created \n')
    print('-------------------------------------------')
    
def check_for_nan(ds): #Checking for NaNs in dataset
        for y,x in product(range(ds.dims['lat']),range(ds.dims['lon'])):
            mask = ds.MASK.isel(lat=y, lon=x)
            if mask==1:
                if np.isnan(ds.isel(lat=y, lon=x).to_array()).any():
                    print('ERROR!!!!!!!!!!! There are NaNs in the dataset')

#Defining files, dates and running function
obs_file =  '/fastdata/ggp21rni/cosipy_input/cosipy_input_CFSR_JIF_1980_2019.nc'
hist_file = '/fastdata/ggp21rni/cosipy_input/cosipy_input_CCSM_JIF_1980_2010.nc'
proj_file = '/fastdata/ggp21rni/cosipy_input/cosipy_input_CCSM_JIF_2030_2060.nc'
cosipy_file ='/fastdata/ggp21rni/cosipy_input/cosipy_input_CCSM_JIF_2030_2060_bias-corr.nc' 


COSIPY_projection_file(obs_file, hist_file, proj_file, cosipy_file)
