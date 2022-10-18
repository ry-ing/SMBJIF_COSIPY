#!/bin/bash
#$ -cwd
#$ -M rning3@sheffield.ac.uk
#$ -m abe
#$ -l h_rt=04:00:00
#$ -N cosipy_model_opt_C
#$ -pe smp 4 #10
#$ -l rmem=8G  #6

# load and activate conda environment
module load apps/python/anaconda3-4.2.0

source activate COSIPYenv

LHS_FILE=model_params_lhc_v2_rrr.csv  #model_params_lhc_v1.csv # model parameters file

sitename=C #C161  still to do: C

# Define paramter arrays
albsnow_arr=()
albfirn_arr=()
albice_arr=()
albtime_arr=()
albdepth_arr=()
rsnow_arr=()
rfirn_arr=()
rice_arr=()
surf_em_arr=()
mult_rrr_arr=() 
vals=($(seq 1 1 100)) #first, step, last (array for file name and run number, should be same size as parameter arrays)
#printf "%s\n" "${vals[@]}"

while IFS="," read -r index albsnow albfirn albice albtime albdepth rsnow rfirn rice surf_em mult_rrr
do
 albsnow_arr+=("$albsnow") #setting csv column to array
 albfirn_arr+=("$albfirn")
 albice_arr+=("$albice")
 albtime_arr+=("$albtime")
 albdepth_arr+=("$albdepth")
 rsnow_arr+=("$rsnow")
 rfirn_arr+=("$rfirn")
 rice_arr+=("$rice")
 surf_em_arr+=("$surf_em")
 mult_rrr_arr+=("$mult_rrr")
done < <(tail -n +2 $LHS_FILE) # only read from the second row


#loop through array parameters
#for i in "${vals[@]}"
for (( i=0; i<${#vals[@]}; i++ ));
do
 printf "Run number: %s\n" "${vals[$i]}"
 albedo_snow=${albsnow_arr[$i]}
 albedo_firn=${albfirn_arr[$i]} 
 albedo_ice=${albice_arr[$i]}
 albedo_time=${albtime_arr[$i]}
 albedo_depth=${albdepth_arr[$i]}
 rlength_snow=${rsnow_arr[$i]}
 rlength_firn=${rfirn_arr[$i]}
 rlength_ice=${rice_arr[$i]}
 surf_em=${surf_em_arr[$i]}
 mult_rrr=${mult_rrr_arr[$i]}
 run_number=${vals[$i]}
 model_file=/fastdata/ggp21rni/cosipy_output/optimization/${sitename}/${sitename}_run_${run_number}_v2.nc

 python COSIPY.py ${albedo_snow} ${albedo_firn} ${albedo_ice} ${albedo_time} ${albedo_depth} ${rlength_snow} ${rlength_firn} ${rlength_ice} ${surf_em} ${mult_rrr} ${model_file}
 #cosipy_script_id=$!
 wait #$cosipy_script_id
 python optimization_eval.py ${sitename} ${run_number} ${model_file}
 #eval_script_id=$!
 wait #$eval_script_id
done
