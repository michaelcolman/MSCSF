#!/bin/sh

# Control stimulus protocol with multiple sites and timings

# Using files rather than creating idealised stimulus regions easily allows for the stimulus to apply to multiple 
# different locations. We can also apply a stimulus protocol where different regions are stimulated at different 
# times (for example, simulating the Purkinje-fibre breakthrough sites in the ventricle during sinus rhythm). To 
# do this, we need to ensure we have the right stimulus file and settings in lib/Tissue.cpp for this functionality 
# (which just requires a stimulus file with different regions numbered sequentially in the order to apply the stimuli, 
# and settings as to what the relative delay of each region is – see the functional_model_test setup for an example) 
# and then to pass “Multi_stim On”. As multi_stim is not the default for this tissue model, we also need to tell it 
# to read the multi_stim stimulus file rather than the standard one. In general, if a model requires multi-stim functionality, 
# it will likely be for control and so set the default stimulus file to the multi-stim file. 
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 100 Multi_stim On stim_file functional_model_test_multi_stim.dat Reference Multi_stim

# First, we can inspect “Map_S1.vtk” and notice it now has 5 different regions (rather than just one), and applying the post-processing:
./bin_to_vtk_tissue Tissue_order geo Tissue_model functional_model_test Reference Multi_stim start_time 0 end_time 100 interval

# We can clearly see the effect of the delayed stimulus points. 


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


