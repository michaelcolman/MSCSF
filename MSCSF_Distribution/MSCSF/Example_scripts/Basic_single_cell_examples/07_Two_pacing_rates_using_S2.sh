#!/bin/sh

# Pacing using two diffent cycle lengths using S2 and NS2

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

BCL=$(printf "%04d" 1000)  # cycle length 1
beats=$(printf "%04d" 50) # number of beats at cycle length 1
S2=$(printf "%04d" 500)  # cycle length 2
NS2=$(printf "%04d" 100) # number of beats at cycle length 2

./model_single_native Model $model BCL $BCL Beats $beats S2 $S2 NS2 $NS2 Reference Two_pacing_rates_using_S2

# This will create an Outputs directory named for the model:  "Outputs_single_native_Two_pacing_rates_using_S2"
# No Results_Reference so directory will simply be "Results"
# "Results/Currents.dat" plotting columns 1 (time) vs 2 (AP) for example - should clearly see change in CL
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


