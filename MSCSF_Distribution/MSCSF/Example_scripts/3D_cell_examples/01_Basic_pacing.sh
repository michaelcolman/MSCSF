#!/bin/sh

# Control pacing using the 3D cell model

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hVM_ORD_s

# testing portion of cell for rapid testing
./model_single_3D Model $model Beats 20 Sim_cell_size testing Reference Control_pacing_${model} Results_Reference testing

# full 3D cell
./model_single_3D Model $model Beats 1 Total_time 1000 Sim_cell_size full Reference Control_pacing_${model} Results_Reference full


# This will create an Outputs directory named for the model:  "Outputs_3Dcell_Control_pacing_model"
# Within which results directory for the two model sizes "Results_{testing/full}"
# Further to "Results_{testing/full}/Currents.dat", which has time courses for the AP and ion currents, 
# "Results_{testing/full}/CRU.dat" contains the intracellular Ca2+ concentrations and RyR, LTCC and SERCA dynamics.
# Linescan data (cross section through centre of cell) "Results_{testing/full}/Ca_linescan_{X/Y_Z}.dat"
# Spatial data for calicum concentrations in the 3D cell (Ca, CaSR, Cads) is "Spatial_Results_{testing/full}"
# Under the default settings, these spatial data will be output in binary format every 5 ms
# See next script for examples of how to convert binary spatial data to plain text and vtk.
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
