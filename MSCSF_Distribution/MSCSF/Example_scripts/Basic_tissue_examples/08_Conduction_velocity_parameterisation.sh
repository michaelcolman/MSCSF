#!/bin/sh

# Using the conduction_velocity model to parameterise a tissue model

# The Tissue_model conduction_velocity has been provided to calculate and parameterise conduction 
# velocities relative to the myocyte orientation. For this specific model only, certain cells are 
# set for the conduction velocity to be automatically calculated in all principal and diagonal 
# directions (rather than requiring post-processing of the results data to make these calculations). 
# Combining this tissue model with direct settings of D1, D_AR and dx allows other tissue model properties 
# to be matched, to provide a clean method to calculate and parameterise these parameters to match 
# experimentally measured conduction velocities. The model must be run using homogeneous, where the influence 
# of celltype on conduction velocity can be investigated by passing different celltypes. 

# By default, the model is 2D and sets the global orientation in the x-direction, and applies a spherical 
# stimulus to the centre; the conduction velocity in the x-direction thus gives the longitudinal velocity, 
# and in the y-direction gives the transverse conduction velocity. The velocities and diffusion parameters are 
# stored in the file “Outputs_tissue_native_X/Conduction_velocity_log.dat” and output to screen. 

# Let’s run the model with default parameters, then set a different dx (to correspond to the space step of an 
# anatomical reconstruction) and cycle through a couple of diffusion parameter combinations to see the effect 
# on the conduction velocity. We’ll use just 1 beat here to save time; ensure you are at stable state for real 
# parameterisation. We’ll use Results_Reference to distinguish individual simulations, such that a list of our 
# conduction velocity data is stored in Conduction_velocity_log.dat.
./model_tissue_native Tissue_model conduction_velocity Beats 1 BCL 500 Reference Conduction_velocity Results_Reference default_params
./model_tissue_native Tissue_model conduction_velocity Beats 1 BCL 500 dx 0.3 Reference Conduction_velocity Results_Reference dx
./model_tissue_native Tissue_model conduction_velocity Beats 1 BCL 500 dx 0.3 D1 0.25 D_AR 8 Reference Conduction_velocity Results_Reference dx_D1_D_AR_1
./model_tissue_native Tissue_model conduction_velocity Beats 1 BCL 500 dx 0.3 D1 0.35 D_AR 10 Reference Conduction_velocity Results_Reference dx_D1_D_AR_2

# We are then interested in Conduction_velocity_log.dat which lists D1 (column 2), D2 (3), dx (4) and the 
# conduction velocities in the x-direction (longitudinal; column 8), y-direction (transverse; column 9) and 
# the xy-diagonal-directions (11 and 12).

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


