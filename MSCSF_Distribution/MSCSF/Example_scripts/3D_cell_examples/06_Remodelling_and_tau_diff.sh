#!/bin/sh

# Illustrating the effect of the two "toy" remodelling models on SR-loading spontaneous release

# We are essentially going to repeat the previous example, with the additional modulations:

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hVM_ORD_s

# Pre pace testing cell 
BCL=$(printf "350")
Beats=$(printf "50")
for rem in none RSERCA_NCX RCRU
do
    ./model_single_3D Model $model Beats ${Beats} BCL ${BCL} ISO 1 Jup_scale 1.5 Sim_cell_size testing Remodelling ${rem} Total_time $(( ${BCL} * ${Beats} )) Write_state ave Reference Spontaneous_release_remodelling Results_Reference Pre_pace_${rem} State_Reference_write pre_pace_${rem}_ave Spatial_output_interval_data 0 Spatial_output_interval_vtk 0
done

# Notice that RSERCA_NCX has a significant effect on SR Ca2+ loading, whereas RCRU does not.

# Now read in state and pace full cell for a couple of beats
Beats=$(printf "2")
for rem in none RSERCA_NCX RCRU
do
    ./model_single_3D Model $model Beats 2 BCL ${BCL} ISO 1 Jup_scale 1.5 Sim_cell_size full Remodelling ${rem} Read_state ave State_Reference_read pre_pace_${rem}_ave Write_state On State_Reference_write pre_pace_${rem}_full Total_time $(( ${BCL} * ${Beats} )) Reference Spontaneous_release_remodelling Results_Reference Pre_pace_${rem}_full Spatial_output_interval_data 0 Spatial_output_interval_vtk 0
done

# And finally, a couple of sims to look at spontaneous release
for rem in none RSERCA_NCX RCRU
do
    for i in 0 1
    do
        run=$(printf "%04d" $i)
        ./model_single_3D Model $model Beats 1 BCL ${BCL} Remodelling ${rem} Read_state On State_Reference_read pre_pace_${rem}_full ISO 1 Jup_scale 1.5 Sim_cell_size full Reference Spontaneous_release_remodelling Results_Reference Spont_${rem}
    done
done


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
