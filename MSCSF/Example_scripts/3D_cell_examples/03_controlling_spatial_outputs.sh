#!/bin/sh

# Control pacing using the 3D cell model; controlling spatial outputs
# By default, the simulations output binary spatial data for Ca, CaDS and CaSR every 5 ms
# We may want to do the following:
#   Output vtk data from the simulation directly for qucik visual checking of spatial data (Cai only)
#   Output binary data at a finer or coarser spatial resolution
#   Not output spatial data at all (to save disc space)
# Which we do using "Spatial_output_interval_{data/vtk}" (set to 0 in order not to output)

# Output binary data every ms, and a vtk Cai datafile every 20 ms
./model_single_3D Model minimal Beats 1 Total_time 100 Sim_cell_size full Spatial_output_interval_data 1 Spatial_output_interval_vtk 20 Reference Determining_spatial_outputs Results_Reference 1ms_binary_plus_vtk

# Output no spatial data
./model_single_3D Model minimal Beats 1 Total_time 100 Sim_cell_size full Spatial_output_interval_data 0 Spatial_output_interval_vtk 0 Reference Determining_spatial_outputs Results_Reference no_spatial

# This will create an Outputs directory named for the model:  "Outputs_3Dcell_Control_pacing_model"
# Within which results directory for the two model sizes "Results_{testing/full}"
# Further to "Results_{testing/full}/Currents.dat", which has time courses for the AP and ion currents, 
# "Results_{testing/full}/CRU.dat" contains the intracellular Ca2+ concentrations and RyR, LTCC and SERCA dynamics.
# Linescan data (cross section through centre of cell) "Results_{testing/full}/Ca_linescan_{X/Y_Z}.dat"
# Spatial data for calicum concentrations in the 3D cell (Ca, CaSR, Cads) is "Spatial_Results_{testing/full}"
# Under the default settings, these spatial data will be output in binary format every 5 ms
# See next script for examples of how to convert binary spatial data to plain text and vtk.
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
