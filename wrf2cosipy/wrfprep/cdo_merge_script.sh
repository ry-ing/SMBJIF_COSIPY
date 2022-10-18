#!/bin/bash
# Request 16 cores in an OpenMP environment
#$ -pe openmp 4
# Request 64 gigabytes of real memory (RAM) 16 cores *4G = 16
#$ -l rmem=4G
#$ -V


# Load the modules required by our program
module load apps/python/anaconda3-4.2.0

source activate COSIPYenv2

cdo mergetime *.nc WRF_CFSR_2019_v2.nc
