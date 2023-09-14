#!/bin/sh

# Visualise spatial data associated with Control pacing using the basic tissue model
# The previous script produced 3 simulations:

# 1D - 3D, tissue model basic, default settings
#./model_tissue_native Model minimal BCL 500 Beats 20 Tissue_model basic Tissue_order 1D Reference Basic_model_pacing_1D
#./model_tissue_native Model minimal BCL 500 Beats 2 Tissue_model basic Tissue_order 2D Reference Basic_model_pacing_2D
#./model_tissue_native Model minimal BCL 500 Beats 1 Tissue_model basic Tissue_order 3D Reference Basic_model_pacing_3D

# In order to convert the binary data to plain text and/or vtk data, we must run:

# ./bin_to_vtk_tissue

# This must be executed in the same directory as the simulations were performed, or a parent directory which contains the 
# “Outputs_” directories required.

# In order to access the data correctly and with the correct settings, we must pass arguments:
#   1 - Any Reference and/or Results_Reference passed (to look in correct Outputs and Spatial_Results directories)
#   2 - The variable we wish to convert (Vm only for native tissue models by default)
#   3 - The start, end and interval time over which convert data (binaries written at interval of 5 ms as default)
#   4 - The Tissue_order and Tissue_model used to perform the simulation (necessary for array sizes etc)
#   5 - The data you want to write: Write_vtk, Write_data, On or Off
#   6 - The model type (integrated or native)

# So, for the 2D and 3D simulatons performed in the previous example, we will run:
# (note – we did not pass a Results_Reference for the simulations, so don’t pass one here!)

# Convert the 2D data into vtk and plain text, for the second beat (t = 500 ms to 1000 ms in 20 ms intervals)
./bin_to_vtk_tissue Reference Basic_model_pacing_2D Tissue_model basic Tissue_order 2D start_time 500 end_time 1000 interval 20 Variable Vm Write_vtk On Write_data On Model_type native

# Convert the 3D data into vtk only, for the first beat (t = 500 ms to 1000 ms in 10 ms intervals)
./bin_to_vtk_tissue Reference Basic_model_pacing_3D Tissue_model basic Tissue_order 3D start_time 0 end_time 500 interval 10 Variable Vm Write_vtk On Write_data Off Model_type native
