#!/bin/bash
#$ -cwd
#$ -M rning3@sheffield.ac.uk
#$ -m abe
#$ -l h_rt=01:00:00 
#$ -N KAN_L_COSIPY
#$ -pe smp 7
#$ -l rmem=10G  

#Memory and time allocations:
#TAKU: 3 HRS, (6 cores, 10GB) for 30 year run
#JIF: 28 HRS, (8 cores, 24GB) for 30 year run
#SSC run: 5 HRS, (6 cores, 10GB) for 24 month run?
#greenland run: 30 minutes, 6 cores 10 GB
#LEMONCREEK: 

module load apps/python/anaconda3-4.2.0

source activate COSIPYenv

python COSIPY.py
