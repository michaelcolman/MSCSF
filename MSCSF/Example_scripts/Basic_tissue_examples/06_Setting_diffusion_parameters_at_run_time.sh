#!/bin/sh

# You can explicitly set the values of the diffusion properties (D1 and D_AR) and the spatial integration step (dx) 
# through command-line arguments. Note: If D_uniformity is not uniform, then D1 and D_AR set the baseline values, 
# from which regional variation is scaled. Let’s run the heterogeneous, anisotropic and regional D simulation from 
# the previous example, except with a larger D1, larger anisotropy ratio, and larger space step. 
# You may have noticed the on-screen outputs for the previous simulation 
# (also in “Outputs_tissue_native_functional_test_het_aniso_regional_D/Settings.dat”) that this model has default 
# settings: D1 = 0.2 mm/ms; D_AR = 4; dx = 0.2 mm.

./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Reference functional_Diffusion_params Total_time 500 Tissue_type heterogeneous Orientation_type anisotropic D_uniformity regional D1 0.3 D_AR 10 dx 0.25

# Now, inspecting the D1 and D2 maps compared to the previous simulation should clearly show the larger D1 and anisotropy ratio 
# (smaller D2 relative to D1). And if we analyse the voltage evolution data, we can see the impact this has (also reflected in the Activation_output.vtk data):
./bin_to_vtk_tissue Tissue_model functional_model_test Tissue_order 2D start_time 0 end_time 500 interval 10 Reference functional_Diffusion_params

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents

