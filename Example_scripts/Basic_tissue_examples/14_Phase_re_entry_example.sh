#!/bin/sh

# Inducing re-entry using the phase-distribution method

# An alternative method to induce re-entry is to use the phase-distribution approach 
# – where state files for a coupled cell are written throughout the duration of an AP 
# (representing 0 - 2π) and then read in as initial-conditions applied to a phase map. 

# First, we must create the phase state files. This must be performed in 1D, with a 
# BCL of 400 ms or longer. This is because the state files are output over an interval 
# of 400 ms, hard coded in the model; you may change this in the appropriate main tissue file. 
./model_tissue_native Tissue_order 1D Tissue_model re-entry Beats 2 BCL 600 Write_state phase Reference Phase_re-entry Results_Reference phase_write

# And now we read in this phase file to a 2D or 3D tissue model (we are using Tissue_mdoel re-entry, 
# but this is not required), and ensure we set beats to 0.
./model_tissue_native Tissue_order 3D Tissue_model re-entry Beats 0 BCL 600 Read_state phase Total_time 500 Reference Phase_re-entry Results_Reference phase_read Spatial_output_interval_vtk 20

# We can visualise the phase map itself in “Outputs_tissue_native_Phase_re-entry/Phase_Map_3D.vtk”.

# Note: we can impose the phase method on a heterogeneous tissue model; doing so requires the phase 
# files to be created for every celltype in the tissue model (the code then automatically reads the 
# correct files for each region).


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


