#!/bin/sh

# Functionality to use two different cell models for different regions has been included. For example, 
# to use a ventricular and atrial cell model for different parts of a whole-heart model. The tissue 
# model must have the appropriate settings, and the celltypes must exist for each region which that 
# cell model is applied to.

# We’ll again use functional_model_test, where the heterogeneity celltypes are set to ENDO, M and RA. 
# We’ll turn Multiple_models On, use the minimal model for the ventricular cell types (model 1) and 
# set the second model (applied for region RA) to an atrial cell model. 
./model_tissue_native Tissue_model functional_model_test Tissue_type heterogeneous Multiple_models On Model minimal Tissue_model_2 hAM_GB Beats 2 Reference Multiple_model_example

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
