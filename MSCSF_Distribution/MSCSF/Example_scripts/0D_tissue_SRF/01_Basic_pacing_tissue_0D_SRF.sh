#!/bin/sh

# We have all the same options as the native tissue models, but can now combine with the SRF.
# SRF are implemented in exactly the same way as 0D cell models.
# We'll increase dt as this is a simple case and it will speed it up

# This example will pre pace a tissue model to load the CaSR, then apply a beat and with SRF
# to see what happens, then convert the spatial data.

# We'll start by reading in the state files from the high CaSR load single cell 0D example, 
# and writing the whole tissue state
./model_tissue_0D Tissue_order 2D Tissue_model basic dt 0.025 BCL 400 ISO 1 Jup_scale 2 Read_state single_cell state_reference_read 50beats_Jup Write_state On State_Reference_write 2D_tissue_CaSR_load Beats 2 Total_time 800 Reference Tissue_SRF_pre_pace

# Now we'll read that state file and apply the RCRU SRF
./model_tissue_0D Tissue_order 2D Tissue_model basic dt 0.025 BCL 400 ISO 1 Jup_scale 2 Beats 1 Read_state On State_Reference_read 2D_tissue_CaSR_load  SRF_mode Dynamic SRF_model 3D_cell SRF_Pset RCRU Reference Tissue_SRF

# And convert Vm and Cai binary spatial data to vtk. Note the "Moidel_type integrated" argument
./bin_to_vtk_tissue Tissue_order 2D Tissue_model basic Model_type integrated Reference Tissue_SRF start_time 0 end_time 2000 interval 10 Variable Vm
./bin_to_vtk_tissue Tissue_order 2D Tissue_model basic Model_type integrated Reference Tissue_SRF start_time 0 end_time 2000 interval 10 Variable Ca

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


