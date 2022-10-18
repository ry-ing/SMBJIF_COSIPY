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

# Load the modules required by our program
module load apps/python/anaconda3-4.2.0

source activate COSIPYenv

cdo --sortname mergetime *.nc 2001_2010_merge_temp.nc



