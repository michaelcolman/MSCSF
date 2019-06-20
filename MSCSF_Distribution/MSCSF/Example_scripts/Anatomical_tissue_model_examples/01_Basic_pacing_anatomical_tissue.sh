#!/bin/sh

# Control pacing using an anatomical tissue model
# Simple example using an anatomically detailed model (Tissue_order geo). We’ll use Human_vent_wedge  for this first illustration. 
./model_tissue_native Tissue_order geo Tissue_model Human_vent_wedge Reference Human_vent_control Beats 1 BCL 500

# And visualise data (note now it is important PATH.txt is in the directory from which “./bin_to_vtk_tissue” is being performed, 
# as it now needs to be able to load the correct tissue geometry file.
./bin_to_vtk_tissue Tissue_order geo Tissue_model Human_vent_wedge Reference Human_vent_control start_time 0 end_time 500 interval 20

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


