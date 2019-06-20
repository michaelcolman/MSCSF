#!/bin/sh

# Idealised tissue models can only have global myocyte orientation (i.e. the same orientation for every cell in the tissue). 
# You can directly control what this orientation is, through either explicitly setting each component (OX, OY, OZ) 
# or by selecting a pre-set global orientation direction using the Global_orientation_direction argument. 
# If setting directly, you only need to specify N-1 components, as the code will automatically set the final component from 
# the normalisation constraint. Obviously, Orientation_type anisotropic must be set in order for the orientation to have any effect.

# Let’s set the stimulus to be centre and sphere, and change the default components explicitly, in both 2D and 3D:
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Stimulus_location_type centre S1_shape sphere Total_time 200 Orientation_type anisotropic OX 1 Reference 2D_OX_1 Spatial_output_interval_vtk 10

./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Stimulus_location_type centre S1_shape sphere Total_time 200 Orientation_type anisotropic OY 0.2656 Reference 2D_OY_0p2656 Spatial_output_interval_vtk 10

./model_tissue_native Tissue_order 3D Tissue_model functional_model_test Stimulus_location_type centre S1_shape sphere Total_time 200 Orientation_type anisotropic OX 0.7 OY 0.1 Reference 3D_OX_0p7_OY_0p1 Spatial_output_interval_vtk 10

./model_tissue_native Tissue_order 3D Tissue_model functional_model_test Stimulus_location_type centre S1_shape sphere  Total_time  200 Orientation_type anisotropic OX 0.2 OZ 0.5 Reference 3D_OX_0p2_OZ_0p5 Spatial_output_interval_vtk 10

# Note that the on-screen settings (and in Settings.dat) prints the values of all the orientation components. Visualise the 
# orientations using “Outputs_X/Orientation.vtk” and check they align with what is expected from the settings, and do indeed 
# control the preferential direction of propagation correctly (which should be apparent from the voltage spatial outputs and the activation map).

# We can also use Global_orientation_direction to set the fibres pointing in any of the principal or diagonal directions. 
# Let’s set Total_time to zero so we can quickly output the Orientation for different pre-set directions.
for orientation in X Y Z XY_plus XY_minus XZ_plus XZ_minus YZ_plus YZ_minus XYZ_ppp XYZ_ppm XYZ_pmp XYZ_mpp
do
    ./model_tissue_native Tissue_order 3D Tissue_model functional_model_test Total_time  0 Orientation_type anisotropic Global_orientation_direction ${orientation} Reference ${orientation}
done

# And again we can inspect the orientation in "Outputs_X/Orientation.vtk"

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


