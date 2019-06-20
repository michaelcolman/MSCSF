#!/bin/sh

# We will now show examples of turning on or off heterogeneity, anisotropy and D_uniformity. To turn any of these on, 
# first, the tissue model must have appropriate settings (in lib/Tissue.cpp) – we will use the Tissue model “functional_model_test” 
# now, as this has settings for all of these functionalies, while the defaults are all off (homogeneous, isotropic and uniform).

# First, let’s run the model with defaults (so we can make comparisons):
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Reference functional_test_defaults Total_time 500

#Note that the on-screen settings show: Tissue type = homogeneous, Orientation type = isotropic and D_uniformity = uniform. 

# Now, let’s make the model heterogeneous, anisotropic, and vary D according to region:
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Reference functional_test_het_aniso_regional_D Total_time 500 Tissue_type heterogeneous Orientation_type anisotropic D_uniformity regional

# First, we can compare “Outputs…Defaults/Geometry_idealised.vtk” with the file in “Outputs…het_aniso_regional_D/” to 
# inspect the geometry and the celltype segmentation. 
# We can also compare “Outputs_tissue_native_functional_test_defaults/Map_{D1/D2}.vtk” with those in 
# “Outputs_tissue_native_functional_test_het_aniso_regional_D” and confirm that they both do, indeed, 
# vary in space. Moreover, “Outputs_tissue_native_functional_test_het_aniso_regional_D” also contains a “
# Orientation.vtk” which contains the fibre orientation (note this is global with them all pointing in the same direction).

# Let’s also inspect the evolution of voltage in the tissue:
./bin_to_vtk_tissue Tissue_model functional_model_test Tissue_order 2D start_time 0 end_time 500 interval 10 Reference functional_test_defaults
./bin_to_vtk_tissue Tissue_model functional_model_test Tissue_order 2D start_time 0 end_time 500 interval 10 Reference functional_test_het_aniso_regional_D

# And we can clearly see the influence of the myocyte orientation, regionally dependent D1 and D2, and the heterogeneity (and thus repolarisation time) on the activation and repolarisation patterns. 

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


