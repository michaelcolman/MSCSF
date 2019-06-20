#!/bin/sh

# There are different ways we can control the stimulus spatial properties from the command-line 
# (i.e. to differ from any default settings). For the Tissue_model basic, which we used for the 
# above examples, the default settings for both S1 and S2 stimulus “edge” with a shape “cuboid” 
# – which simply means stimulate from the small-x end of the tissue, using a cuboid stimulus shape.

# The first way we can change the stimulus spatial settings is to pass in some of the other 
# pre-defined types, and change the shape from cuboid to sphere. These pre-defined types set 
# the location and sizes of the stimulus in all directions, based on the tissue dimensions:
# First, let’s make the S1 stimulus from the centre and a sphere/circular shape 
# (we limit total time as this is sufficient to inspect the stimulus):
./model_tissue_native Tissue_model basic Tissue_order 2D Stimulus_location_type centre S1_shape sphere Reference Stimulus_S1_centre_sphere Total_time 20 Spatial_output_interval_vtk 5

# We can now inspect the stimulus site uisng "Outputs_tissue_native_Stimulus_S1_centre_sphere/Map_S1.vtk"
# and can check the vtk files for voltage evolve as we would expect.

# We can set the S1 stimulus to be “edge”, “centre” or “whole-tissue” and shape to be “cuboid” or “sphere”. 
# Using “S2_Stimulus_location_type”, we can also do the same for the S2 stimulus, where now our options 
# are: “S1” – to be identical to the S1 stimulus; “edge”, “centre” or “cross_field”.
# Let’s try setting an S2 of 300 ms and set the stimulus to be “cross_field”:
./model_tissue_native Tissue_model basic Tissue_order 2D Reference Stimulus_S2_cross_field Total_time 500 Beats 1 S2 300 S2_Stimulus_location_type cross_field Spatial_output_interval_vtk 5

# Finally, we can also control explicitly the location (centre-points) and 
# sizes (distance from the centre to the edge) of both stimuli. We can define 
# every single value ourselves, or start from a pre-defined setting and overwrite just some parameters.

# Let’s start with Stimulus_location_type edge, but move it in the x-direction; let’s also set the S2 stimulus to be centre, move its centre in y and adjust its size in x and y:

./model_tissue_native Tissue_model basic Tissue_order 2D Reference Explicit_stimulus_settings Stimulus_location_type edge S1_shape cuboid S1_x_loc 50 S2_Stimulus_location_type centre S2_shape cuboid S2_y_loc 25 S2_x_size 15 S2_y_size 5 S2 300 Total_time 500 

# Note: if we set the shape to “sphere”, then the only size parameter is the radius (as elliptical stimulus 
# functionality has not been included yet) – the only size argument we can use is “..x_size”, to define the 
# radius; y_size and z_size will have no effect if shape is sphere. 

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


