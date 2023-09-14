#!/bin/sh

# Controlling stimulus files and/or cooridinates:

# We can specify different files for stimulus (S1 or S2), or set the stimulus using coordinates 
# (just as with the idealised tissue models).

# Let’s first see what happens under default settings for this model – we’ll set total time to 0 
# as we are only interested in the stimulus locations, which we can see from the map files produced. 
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 Reference Default_stimulus

# You may notice the on-screen output during setup: 
# “File loaded /…/MSCSF_state_and_geometry_files/Tissue_geometries/functional_model_test_stim.dat”, which is the 
# default file set in lib/Tissue.cpp, and we can inspect what this stimulus location will be in 
# “Outputs../Default_stimulus/Map_S1.vtk”. 

# We can also pass an S2 argument to see what the default S2 stimulus file is and its location. 
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 S2 50 Reference Default_stimulus_S2

# And we can see that it loads file “functional_model_test_S2.dat”, which is output in “../Map_S2.vtk”. 
# There is also another S2 stimulus file available to us – “functional_model_test_S2_2.dat”. Let’s now 
# set the S1 stimulus file to being the default S2 stimulus file, and the S2 file to be this alternative 
# (there is nothing specific about S1 or S2 stimulus files – we can just name them as is most useful for us). 
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 S2 50 Reference Different_stimulus_files stim_file functional_model_test_S2.dat S2_stim_file functional_model_test_S2_2.dat

# And we now notice that “Map_S1.vtk” looks like the previous S2 (as we would expect), while “Map_S2.vtk” is clearly a new file. 
# We can also set either or both stimuli by coordinates by passing “{S2_}Stimulus_loc_type coords”. 
# There are default coordinate settings in the tissue model, so in the first instance we don’t need to speficy any coordinate parameters directly:
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 Stimulus_location_type coords S2_Stimulus_location_type coords S2 50 Reference Default_coords_stimulus

# And we can see that the S2 stimulus is perfectly setup for cross-field stimulation. We can also control any of the stimulus location 
# parameters directly, just as with example 4 for idealised models. Let’s adjust some S2 parameters:
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 Stimulus_location_type coords S2_Stimulus_location_type coords S2 50 S2_shape cuboid S2_y_loc 25 S2_x_size 15 S2_y_size 5 S2_z_size 1 Reference run_time_S2_coords_stimulus

# Note: if the tissue model does not have default settings for the stimulus coordinates, all parameters must be
# set at run-time if "coords" set.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


