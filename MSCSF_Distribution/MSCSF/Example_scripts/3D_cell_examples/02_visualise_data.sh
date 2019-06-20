#!/bin/sh

# Visualise spatial data associated with Control pacing using the 3D cell model
# The previous script produced 2 simulations:

# portion of cell for rapid testing
#./model_single_3D Model $model Beats 20 Sim_cell_size testing Reference Control_pacing_${model} Results_Reference testing

# full 3D cell
#./model_single_3D Model $model Beats 1 Total_time 100 Sim_cell_size full Reference Control_pacing_${model} Results_Reference full

# Let's also run a third simulation with different baseline cell size option (no Results Reference passed here):
./model_single_3D Model minimal Beats 1 Total_time 100 Cell_size thin Sim_cell_size full Reference Thin_cell_pacing_minimal

# In order to convert the binary data to plain text and/or vtk data, we must run:

# ./bin_to_vtk_3Dcell

# In order to access the data correctly and with the correct settings, we must pass arguments:
#   1 - Any Reference and/or Results_Reference passed (to look in correct Outputs and Spatial_Results directories)
#   2 - The variable we wish to convert (Ca, CaSR, CaDS are all written as binary)
#   3 - The start, end and interval time over which convert data (binaries written at interval of 5 ms as default)
#   4 - The Sim_cell_size to match the simulation
#   5 - The actual cell size/type (standard, thin)
#   6 - The data you want to write: Write_vtk, Write_data, Write_slices  - On or Off
#   7 - If writing slices, you can also specify the coordinate of the normal to slice at which the slice is taken (XY_slice_z x etc)


# Visualise intracellular calcium at 10 ms intervals in testing cell control simulation
./bin_to_vtk_3Dcell Reference Control_pacing_minimal Results_Reference testing Sim_cell_size testing start_time 0 end_time 50 interval 10 Variable Ca

# Visualise dyadic calicum at 5 ms intervals; write vtk and text data; full cell control simulation
./bin_to_vtk_3Dcell Reference Control_pacing_minimal Results_Reference full Sim_cell_size full start_time 0 end_time 50 interval 5 Variable Ca Write_data On

# Visualise intracellular calcium slices (plain text only) for the thin cell full size simulation (ran above); output xz slices at surface and centre
# Don't produce vtk files or 3D text data for this case
./bin_to_vtk_3Dcell Reference Thin_cell_pacing_minimal Sim_cell_size full start_time 0 end_time 50 interval 5 Variable Ca Write_data Off Write_vtk Off Write_slices On XZ_slice_y 0
./bin_to_vtk_3Dcell Reference Thin_cell_pacing_minimal Sim_cell_size full start_time 0 end_time 50 interval 5 Variable Ca Write_data Off Write_vtk Off Write_slices On XZ_slice_y 7
./bin_to_vtk_3Dcell Reference Thin_cell_pacing_minimal Sim_cell_size full start_time 0 end_time 50 interval 5 Variable Ca Write_data Off Write_vtk Off Write_slices On XZ_slice_y 15




