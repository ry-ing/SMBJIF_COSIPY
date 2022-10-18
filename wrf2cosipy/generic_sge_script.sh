#!/bin/bash
#$ -cwd
#$ -M rning3@sheffield.ac.uk
#$ -m abe
#$ -l h_rt=09:30:00 #01:00:00
#$ -N bias_correction
#$ -pe openmp 5 #4
#$ -l rmem=6G #6
module load apps/python/anaconda3-4.2.0

source activate COSIPYenv

python wrf_bias_corr.py 
