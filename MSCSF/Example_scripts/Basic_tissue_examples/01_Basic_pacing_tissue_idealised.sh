#!/bin/sh

# Control pacing in a tissue model using different dimensions

    ./model_tissue_native Model minimal BCL 500 Beats 20 Tissue_model basic Tissue_order 1D Reference Basic_model_pacing_1D
    ./model_tissue_native Model minimal BCL 500 Beats 2 Tissue_model basic Tissue_order 2D Reference Basic_model_pacing_2D
    ./model_tissue_native Model minimal BCL 500 Beats 1 Tissue_model basic Tissue_order 3D Reference Basic_model_pacing_3D


# This will create an Outputs directory named for the model:  "Outputs_tissue_native_Basic_model_pacing_{1D/2D/3D}"
# In the Outputs directories are the files "Map_{S1/S2}.vtk" and "Map_{D1/D2}.vtk".
# Notice for the 1D models, these contain just one line of data, with either the value of D1 or, for the sitm map, 
# 1s for stimlus regions and 0s for non-stimulus regions.
# For the 2D or 3D models, we can visualise this using visualisation software (such as the free Paraview (tm) program).
# Our results are now in "Results", which contains files "Currents_cell{1-3}.dat" and "Properties_cell{1-3}.dat" -
# simply the single-cell outputs, as produced using model_single_native, for three different cells in the tissue.
# "Results/Vm_linescan_x" contains a linescan of voltage in the x-direction, each line for 1 ms interals.
# "Spatial_Results" then contains outputs of the voltage in all space at specific time points (by default, binary data
# with an interval of 5 ms).
# As well as, on completion of the simulation, the activation map (the most recent excitation time for each cell)
# The next script gives an example of how to convert these binary data into plain text and/or vtk.
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


