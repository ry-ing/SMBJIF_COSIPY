#!/bin/bash -i
#$ -cwd
#$ -M rning3@sheffield.ac.uk
#$ -m abe
#$ -l h_rt=01:00:00
#$ -N cdo_merge
# Request 16 cores in an OpenMP environment
#$ -pe openmp 4
# Request 64 gigabytes of real memory (RAM) 16 cores *4G = 16
#$ -l rmem=4G

YEAR=2010

echo ${YEAR}

# Load the modules required by our program
module load apps/python/anaconda3-4.2.0

source activate COSIPYenv

cdo mergetime *.nc WRF_CCSM_${YEAR}.nc

cdo selvar,lat,lon,SLP,PSFC,U10,V10,T2,Q2,PCPT,ACSNOW,LWDNB,SWDNB WRF_CCSM_${YEAR}.nc /fastdata/ggp21rni/wrf_data/ccsm_4km/merged_files/2001_2010/WRF_CCSM_${YEAR}_selvar_v1.nc



