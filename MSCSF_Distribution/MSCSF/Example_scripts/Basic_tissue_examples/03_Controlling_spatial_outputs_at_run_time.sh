#!/bin/sh
# We may want to change the interval over which binary data are written, or to output vtk files directly from the 
# simulation for fast inspection of results, or output no spatial data at all if the excitation data is sufficient 
# and we want to save disk space. To do this, we use “Spatial_output_interval_{data/vtk} [x]”, where “0” will 
# prevent data being output.

# We will run the same simulation as example 1, for the 2D case, except now we will output binary spatial data
# every ms, and a vtk for inspection every 10 ms

./model_tissue_native Model minimal BCL 500 Beats 2 Tissue_model basic Tissue_order 2D Reference Basic_model_pacing_2D Results_Reference Controlling_output_intervals Spatial_output_interval_data 1 Spatial_output_interval_vtk 10

# Or output no spatial data:
./model_tissue_native Model minimal BCL 500 Beats 2 Tissue_model basic Tissue_order 2D Reference Basic_model_pacing_2D Results_Reference Preventing_spatial_data_output Spatial_output_interval_data 0 Spatial_output_interval_vtk 0

# Notice now, in "Outputs_tissue_native_Basic_model_pacing_2D", we have new Results and Spatial_Results directories, 
# with the Results_References appended. Note that in "../Spatial_Results_Controlling_output_intervals", there are binary
# data at intervals of 1 ms rather than 5 ms, and a vtk at intervals of 10 ms, while in "../Spatial_Results/Preventing_spatial_data_output" there are no spatial data.
# Also note that the Currents data in "../Results_X/" are the same for all simulations, as the conditions were identical.

# We may now want to convert the 1ms interval binary data into vtk files, for which we use the "bin_to_vtk_tissue" tool, as before, 
# but now need to also pass the Results Reference
./bin_to_vtk_tissue Reference Basic_model_pacing_2D Tissue_model basic Tissue_order 2D start_time 500 end_time 700 interval 1 Variable Vm Write_vtk On Write_data On Model_type native Results_Reference Controlling_output_intervals
# Note how now, vtk files before 500 ms have intervals of 10 ms (as output by the simulation), whereas between 500 and 700 the
# vtk interval is 1ms. 

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


