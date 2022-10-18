#!/usr/bin/env python

#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#import netCDF4 as nc
import datetime
from dateutil.relativedelta import relativedelta
#import matplotlib.dates as mdates
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import math

import argparse

parser = argparse.ArgumentParser(description='COSIPY evaluate')
parser.add_argument('site', help='site name or glacier')
parser.add_argument('run_number', help='number of model run')
parser.add_argument('model_file', help='model output file')
parser.add_argument('obs_file', help='total obs smb measurement file')
    
args = parser.parse_args()

site = args.site
run_number = args.run_number
run_number = float(run_number)
model_file = args.model_file
obs_file = args.obs_file
glacier = site

#print('RUN NUMBER IS', run_number)

#Importing stake data from JIRP
#stake_file = stake_file
obs_data = pd.read_csv(obs_file, delimiter=',', index_col='Ba_Date')

#Importing COSIPY output file 
DATA = xr.open_dataset(model_file)
DATA['time'] = pd.to_datetime(DATA['time'].values).strftime('%Y-%m-%d')

def stake_evaluation(obs_data, DATA, start_date, end_date, site, run_number, outfile_smb, outfile_stats):

    obs_dates = pd.to_datetime(obs_data.index, format='%d/%m/%Y') # convert time to datetime and set as dataframe index
        
    obs_data = obs_data.loc[(obs_dates > start_date)  &   (obs_dates < end_date)] #select dates for evaluation
    obs_dates = pd.to_datetime(obs_data.index, format='%d/%m/%Y') # convert to datetime object
    
    obs_data = obs_data['Ba']
           
    model_smb = np.zeros(len(obs_data)) #create array for model SMB
    for i in range(len(obs_dates)):
    
          lon = DATA.lon[:].values
          lat = DATA.lat[:].values
            
          edate = obs_dates[i]
          sdate = obs_dates[i] - relativedelta(years=1) #1 year SMB period to match observations
                
          edate = pd.to_datetime(edate).strftime('%Y-%m-%d') # convert to string type
          sdate = pd.to_datetime(sdate).strftime('%Y-%m-%d')
                        
          mask_sel = (DATA['time'] > sdate) & (DATA['time'] <= edate) # select time period from model for evaluation against observations
          smb_period = DATA.surfMB.loc[mask_sel]

          glacier_cells = np.nansum(DATA.MASK.values)
          sum_smb = smb_period.sum() / glacier_cells
          model_smb[i] = sum_smb.values
                
    rmse_smb = math.sqrt(mean_squared_error(obs_data.values, model_smb)) #calculate RMSE
    rmse_smb = round(rmse_smb, 2)
    
    r_smb = r2_score(obs_data.values, model_smb) #calculate R squared score
    r_smb = round(r_smb, 2)
    
    
    
    #print('RUN NUMBER IS:', run_number)
    if run_number == 1.0:  
       print('RUN NUMBER IS 1, creating output file')  
       obs_data = obs_data.values
       df_array = np.array([obs_data, model_smb])
       df_array = np.transpose(df_array)
       smb_dataframe = pd.DataFrame(df_array, columns=[glacier, run_number], index=obs_dates)
       smb_dataframe.to_csv(outfile_smb)
       
       stats = [run_number, rmse_smb, r_smb]
       stats_dataframe = pd.DataFrame(columns=['MODEL_RUN', 'RMSE', 'R2'])
       
       row_index = run_number - 1.0
       stats_dataframe.loc[row_index] = stats
       stats_dataframe.to_csv(outfile_stats)
       
       
    else: #load in dataframe, add column, save dataframe
       smb_dataframe = pd.read_csv(outfile_smb, delimiter=',', index_col=0)
       smb_dataframe[run_number] = model_smb.tolist() 
       smb_dataframe.to_csv(outfile_smb)
       
       stats_dataframe = pd.read_csv(outfile_stats, delimiter=',', index_col=0)
       #stats = [run_number, rmse_smb, r_smb]
       
       stats = {'MODEL_RUN':run_number, 'RMSE':rmse_smb, 'R2':r_smb}
       stats_dataframe = stats_dataframe.append(stats, ignore_index=True)
       
       #row_index = run_number - 1.0
       #stats_dataframe.loc[row_index] = stats
       stats_dataframe.to_csv(outfile_stats)
       

outfile_smb='/fastdata/ggp21rni/cosipy_output/optimization/%s/output/%s_2000-2010_SMB_results_v1.csv' % (site, site)
outfile_stats='/fastdata/ggp21rni/cosipy_output/optimization/%s/output/%s_2000-2010_model_score_v1.csv' % (site, site)

stake_evaluation(obs_data, DATA, '2003-08-01', '2008-12-31', site, run_number, outfile_smb, outfile_stats)
