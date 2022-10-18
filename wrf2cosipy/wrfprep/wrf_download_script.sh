#!/bin/bash
# Request 5 gigabytes of real memory (mem)
#$ -l rmem=5G
# Email notifications 
#$ -M ryan.ing@outlook.com
# Email notifications if the job aborts
#$ -m a


aws s3 cp s3://http://wrf-se-ak-ar5/cfsr/4km/daily/ /fastdata/ggp21rni/wrf_data/cfsr_4km --no-sign-request
