#!/bin/sh

# Example of protocol to produce multiple APDr relationships (models and modulation)
# Assuming we want to store raw data, and have unique APDr data for all conditions
# Can organise in different ways, but need to ensure:
#   An Outputs directory name ("Reference X") unique to each condition
#   Raw data within each Outputs for every BCL
# The example here will create a new simulations folder to store all of the data, within
# which model and conditon specific Outputs directories will be created by the simulation
# Instead of creating the new parent directory, you could simply add "APDr_" to the Reference
# The new directory just seems cleaner, especially in the context of running multiple of these
# script examples from within the same directory.

# Create parent simulation directory, copy executable and PATH to it, then nagivate to it
mkdir APD_multiple_examples
cp model_single_native PATH.txt APD_multiple_examples
cd APD_multiple_examples

# run simulations
for model in hAM_CRN hAM_GB # will run for both models
do
    for remodelling in none AF_Col_4
    do
        for i in 250 300 350 400 500 750 1000 1500 2000
        do
            BCL=$(printf "%0d" $i)
            ./model_single_native Model $model BCL $BCL Remodelling ${remodelling} Beats 20 Reference ${model}_${remodelling} Results_Reference BCL_${BCL}    #Don't use 20 beats for real data!!!!
        done
    done
done

# Now, all within "APD_multiple_examples"
# This will create Outputs directories "Outputs_single_native_model_remodelling"
# Within which results directory for the all BCLs "Results_BCL_X"
# Then "Outputs_single_native_model_remodelling/Properties_log.dat" will contain a list of data 
# from all simulations, BCL vs APD, dvdt etc, for that model and remodelling condition
# And "Results_BCL/Currents.dat" will contain the raw AP and currents for each BCL
# If you do not need to store the raw data, then simply don't pass a Results_Reference; 
# the Results directory will then contain most recent simulation only 
# Now can plot Outputs_single_native_model_remodelling/Properties_log.dat using 1 vs 11 
# for all model and remodeling to overlay APD_90 of final beat for each condition, for example
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
