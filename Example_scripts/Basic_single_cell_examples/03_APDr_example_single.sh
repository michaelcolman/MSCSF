#!/bin/sh

# Example of protocol to produce a single APDr relationship

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

for i in 250 300 350 400 500 750 1000 1500 2000
do
    BCL=$(printf "%0d" $i)
    ./model_single_native Model $model BCL $BCL Beats 200 Reference APDr_example Results_Reference BCL_${BCL}
done

# This will create an Outputs directory named:  "Outputs_single_native_APDr_example""
# Within which results directory for the all BCLs "Results_BCL_X"
# Then "Outputs_single_native_APDr_example/Properties_log.dat" will contain a list of data from all simulations, BCL vs APD, dvdt etc
# And "Results_BCL/Currents.dat" will contain the raw AP and currents for each BCL
# If you do not need to store the raw data, then simply don't pass a Results_Reference; the Results directory will then contain most recent simulation only
# Plot Outputs_single_native_APDr_example/Properties_log.dat using colums 1 vs 11 for BCL vs APD_90 and 1 vs 15 for BCL vs minimum voltage for final beat
# properties
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
