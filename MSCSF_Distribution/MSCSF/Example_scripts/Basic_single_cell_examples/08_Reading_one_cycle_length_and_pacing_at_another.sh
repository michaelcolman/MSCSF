#!/bin/sh

# Pacing using two diffent cycle lengths
# Here, we pace at CL 1 to stability and write the state
# Then, we read the state and pace at a different cycle length
# E.g. effect of change in CL once stable state has been reached, 
# wanting multiple new CLs and to save simulation time
# Because the state files are named according to BCL, if we want to read one
# BCL and pace at another, we must use S2 and NS2 to make the S2 our regular pacing

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc

BCL=$(printf "%04d" 1000)  # cycle length 1

# Do store write data, so we can plot both cycle lengths if we want (will have to offset time)
./model_single_native Model $model BCL $BCL Beats 200 Write_state On Reference Reading_and_pacing_at_different_CLs Results_Reference CL1

# Now read in state at that BCL, but sets beats to 1 such that only the S2 train is applied
# (first stimulus, beat 1 at BCL is at t = 0 ms and independent of BCL; second stimulus is thus
# first S2 beat if Beats is set to 1, at S2 CL)

for S2 in 300 700
do 
    ./model_single_native Model $model Read_state On BCL $BCL Beats 1 S2 $S2 NS2 50 Reference Reading_and_pacing_at_different_CLs Results_Reference CL2_$S2
done

# This will create an Outputs directory named for the model:  "Outputs_single_native_Reading_and_pacing_at_different_CLs"
# With three results directories, one for the first CL and one for each second cycle length (S2)
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


