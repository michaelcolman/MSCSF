#!/bin/sh

# Control pacing pacing at different cycle lengths

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

BCL=$(printf "%04d" 750)
./model_single_native Model $model BCL $BCL Beats 200 Reference Control_pacing_${model} Results_Reference BCL_${BCL}

BCL=$(printf "%04d" 350)
./model_single_native Model $model BCL $BCL Beats 200 Reference Control_pacing_${model} Results_Reference BCL_${BCL}

# This will create an Outputs directory named for the model:  "Outputs_single_native_Control_pacing_model"
# Within which results directory for the two BCLs "Results_BCL"
# Then "Outputs_single_native_Control_pacing_model/Properties_log.dat" will have 2 entries containing BCL, APD, dvdt for the two simulations
# And "Results_BCL/Currents.dat" will contain the AP and currents for that simulation BCL
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


