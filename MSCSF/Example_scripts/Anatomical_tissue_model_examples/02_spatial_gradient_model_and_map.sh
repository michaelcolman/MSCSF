#!/bin/sh

# Control pacing using a spatial gradient model with underlying regional heterogeneity

# Unlike idealised models, we can now also impose a spatial-gradient model. This may be apico-basal 
# ventricular heterogeneity, or a gradual shortening of the AP with distance from SAN, or similar. 
# i.e. a spatial gradient imposed on top of possible underlying cellular heterogeneity.

# All we need is:
#   1) a map file containing the spatial gradient (continuous 0 – 1; must be in PATH -> Tissue_geometries) 
#   2) an electrophysiology gradient model to apply (we have two testing examples implemented, “apico_basal_example” 
#   and “SAN_distance_example”, neither of which are based on real data or observations). 
# Functional model test has an associated spatial-gradient file. Even though this filename is set in the tissue 
# model settings (in lib/Tissue.cpp), we will specify it at run-time as an example of how to do this (we could select 
# between multiple available files this way) – note we just specify the file name, which is stored in PATH. In general, 
# we do not need to pass in the filename directly; any default filenames set in the tissue model will be used, unless overwritten.
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Spatial_gradient apico_basal_example Spatial_gradient_map_file functional_model_test_spatial_gradient.dat BCL 500 Beats 1 Reference Spatial_gradient_example

# We can then have a look at the spatial gradient map, as well as the geometry itself, in “Outputs../Map_Spatial_gradient.vtk” 
# – note that you may have to rescale the colourmap to 0-1 (rather than scaling for empty space of -100). Note also that the 
# geometry is heterogeneous, and the spatial gradient applies on top of that
./bin_to_vtk_tissue Tissue_order geo Tissue_model functional_model_test Reference Spatial_gradient_example start_time 0 end_time 350 interval 10

# And now we can see the effect of the spatial gradient and heterogeneity on repolarisation; the default D_uniformity 
# is also regional, so there is variation in conduction velocity also.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


