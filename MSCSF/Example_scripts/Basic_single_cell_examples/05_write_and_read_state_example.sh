#!/bin/sh

# Control pacing pacing at a single cycle length (with a modulation example), pacing to stable conditions, 
# Writing then state variables to file
# Then running a new simulation reading in the state files

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

BCL=$(printf "%04d" 350)

# Pace to stable conditions and write state
./model_single_native Model $model BCL $BCL Beats 200 ISO 1 Write_state On Reference State_write_read_example Results_Reference Write

# Read state and pace for just a few beats (now already at stable conditions)
./model_single_native Model $model BCL $BCL Beats 5 ISO 1 Read_state On Reference State_write_read_example Results_Reference Read

# This will create an Outputs directory:  "Outputs_single_native_State_write_read_example""
# Within which results directory for the two simulations ("Results_Write", "Results_Read")
# You can check the measured parameters at the end of each simulation match (i.e. stable, as 200 vs 5 beats)
# using "Outputs_single_native_State_write_read_example/Properties_log.dat" and comparing the 2 lines
# NOTE: If using the above settings, the AP exhibits alternans and so the entries in Properties_log.dat (which contains data
# for the final 2 beats) are swapped for adjacent paired columns as the simulations ended on different alternate. 
# so while the entries are not identical, the data and simulations are and the new simulation has started in this alternan
# stable condition.
# And can plot the APs also to check the same (offset the time axis of one to align final beats of initial simulation
# with the 5 single beats of read simulation)
# "Results_Write/Currents.dat" using column 1 (data minus 69300) vs 2, and "Results_Read/Currents.dat" using columns 1 vs 2

# Note that state files are saved independently of the results data (being placed in the PATH directory), and so you need not 
# worry about accidentally deleting state files when deleting redundant or old data.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


