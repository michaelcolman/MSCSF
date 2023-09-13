#!/bin/sh

# Apply rapid pacing protocol to load the CaSR, then run multiple simulations to produce
# a dataset on which spontaneous relase can be analysed
# This will also explore and use the State_write "ave" or "On" functionality (similar to tissue models)
# In this example, we are not rigorously testing to ensure we have reached a stable-state, 
# but you should do this for real investigation!

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hVM_ORD_s

# Rapid pace testing cell and write state as whole-cell averages; apply ISO and further enhance Jup in order to promote CaSR loading
# we'll also set tau_ss_type to medium_fast to promote spontaneous release, for this example.
# do this to the testing cell size in order to produce state file more quickly (or full size if too noisy)
# Note: The spatial and 0D models, by default, apply pacing and then leave a 2 second quiescent period to analyse
# spontaneous release; in order to correctly write the state file, we need to explicitly set the Total_time
# State write to ave, as we want to read the state into the full sized model later - and clearly cannot do this with a spatial state file

BCL=$(printf "350")
Beats=$(printf "50")
./model_single_3D Model $model Beats ${Beats} BCL ${BCL} ISO 1 Jup_scale 2 tau_ss_type medium_fast Sim_cell_size testing Total_time $(( ${BCL} * ${Beats} )) Write_state ave Reference Spontaneous_release Results_Reference Pre_pace

# We can check in "Outputs_3Dcell_Spontaneous_release/Results_Pre_pace/CRU.dat" using columns 1 and 7 to plot CaSR against time;
# in this example, it loads to ~1100 mM by the end of the simulation (just enough to induce some spontaneous release!)

# Now read in that state (which is written as whole-cell averages) into the full-size cell model.
# Pace this model for a couple of beats, and write the whole state file which includes the spatial
# state (local concentrations and RyR/LTCC state). Here you'd probably want more than 2 beats for real work.
Beats=$(printf "2")
./model_single_3D Model $model Beats 2 BCL ${BCL} Read_state On ISO 1 Jup_scale 2 tau_ss_type medium_fast Sim_cell_size full Read_state ave Write_state On Total_time $(( ${BCL} * ${Beats} )) Reference Spontaneous_release Results_Reference Pre_pace_full

# Then run multiple simulations to produce a dataset which can be analysed statistically
# Apply just 1 beat, as this is why we wrote the state file
# Produce spatial data outputs for the first 2 simulations (to visualise spontnaeous calcium activity).
# We now want the Total_time to include a 2 second quiescent period, so no need to pass the argument.
for i in 0 1
do
    run=$(printf "%04d" $i)
    ./model_single_3D Model $model Beats 1 BCL ${BCL} Read_state On ISO 1 Jup_scale 2 tau_ss_type medium_fast Sim_cell_size full Reference Spontaneous_release Results_Reference Run_${run}_spatial_vis
done

# but not for the others in order to save disk space. 10 total here - need more for actual data.
for i in `seq 2 1 9`
do
    run=$(printf "%04d" $i)
    ./model_single_3D Model $model Beats 1 BCL ${BCL} Read_state On ISO 1 Jup_scale 2 tau_ss_type medium_fast Sim_cell_size full Spatial_output_interval_data 0 Reference Spontaneous_release Results_Reference Run_${run}
done

# We have a look at "Outputs_3Dcell_Spontaneous_release/Results_Run_0000_spatial_vis/CRU.dat" looking at columns 5 (Cai) and 9 (open RyR state occpancy)
# and indeed see there is spontaneous activity after the 5th applied beat (likley for at least one of the two visualisation simulations under these example conditions).
# So we may want to convert the binary spatial data into vtks for visualisation.
./bin_to_vtk_3Dcell Sim_cell_size portion Reference Spontaneous_release Results_Reference Run_0000_spatial_vis start_time 300 end_time 2350 interval 10 Variable Ca

# And use open RyR state (CRU.dat column 9) to analyse spontaneous release timing statistics, using some threshold of open RyR occpancy
# occuring after 300 ms to determine if a sptonaneous release has occured and its timing, or simply plot the overlays of all
# simulations' open RyR or Cai (CRU.dat columns 9 or 5) to inspect spontaneous release variability.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
