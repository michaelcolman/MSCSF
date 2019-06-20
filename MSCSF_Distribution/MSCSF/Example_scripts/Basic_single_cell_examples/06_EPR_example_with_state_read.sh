#!/bin/sh

# Example of an ERP protocol using S1-S2 pacing for control and remodelling

# As we need to run multiple simulations (S2) for a single cycle length (BCL), 
# we want to write the state file for the BCL at stable state to reduce the number
# of beats required for each S2 simulation.
# We may also want to do this for multiple BCLs
# Example here just does two to save on time
# Similarly, the S2 interval is 5 ms to save time - you may want to first run a coarse
# S2 interval to identify the ERP range, then run with a finer interval (1ms) within
# just this range.
# As with APDr, we need to think about data storage so that Properties_log contains a list of all
# S2 for each individual condition. Thus we want:
#   An Outputs directory name ("Reference X") unique to each condition and BCL
#   Raw data within each Outputs for every S2 (optional)


model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

# First, pace to stability and write state for each BCL and condition
for i in 400 1000
do
    BCL=$(printf "%0d" $i)
    for remodelling in none AF_Col_4
    do
        ./model_single_native Model $model BCL $BCL Remodelling $remodelling Beats 200 Write_state On
    done
done

# As with multiple APD example, we will create a new simulations folder to tidy data here
# (above data is not stored in this example)
mkdir ERP_example
cp model_single_native PATH.txt ERP_example
cd ERP_example

# Now, cycle through S2 intervals for each BCL and condition
for i in 400 1000
do
    BCL=$(printf "%d" $i)
    for remodelling in none AF_Col_4
    do
        for j in `seq 400 -5 80`
        do
            S2=$(printf "%04d" $j)
            ./model_single_native Model $model BCL $BCL Remodelling $remodelling Beats 5 Read_state On S2 $S2 Reference BCL_${BCL}_rem_${remodelling} Results_Reference ${S2} 
        done
    done
done


# Within "ERP_example"
# We will have 4 Output directorties "Outputs_single_native_BCL_X_rem_Y"
# With "Outputs_single_native_BCL_X_rem_Y/Properties_log.dat" containing the measured data (relevantly Vmin, Vmax, Vamp) for each S2.
# Plot: column 2 (S2) vs 15 (Vmin) 17 (Vmax) or 19 (Vamp) for these properties in the final beat (i.e. S2)
# Can compare to columns 16, 18 and 20 for the same properties on the penultimate beat (i.e. S1 and thus should be a flat line)
# And/Or inspect / run a script to determine at what S2 the Vamp (column 19) reaches below e.g. 80% of that for S1 (column 20)
# Most efficient way would be to run coarse S2, inspect by plotting to determine actual S2 range, then run script within that for 80%.
# If doing this, use a new Reference for the finer-S2 interval simulation, or delete the original coarse one, such that Properties_log is 
# specific to the finer simulation.
# May also want dvdt_max: column 13 for final, S2 beat and 14 for penultimate, S1 beat.
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
