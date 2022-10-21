COSIPY.py, config.py and constants.py contain the COSIPY parameters and options used in the main Juneau icefield wide runs

COSIPY_opt.py, config_opt.py and constants_opt.py contains the COSIPY options for the optimisation model ensembles run by the cosipy_sge_optimisaztion.sh job shell script. This script runs the model in a loop with the model parameters contained within 'model_params_lhc_v2_rrr.csv', and then sends the results to optimization_eval_dist.py to evaluate the modelled results and calculate RMSE and R2 for each run.
