#!/bin/sh

# Example re-entry induction and state save protocol

# Tissue_model re-entry has been designed for simulations of re-entry in idealised 2D and 3D sheets. 
# It is a relatively large tissue, with S1 stimulus from the x=0 edge and the default 
# S2_Stimulus_location_type set to “cross_field”. Thus, passing an S2 argument will apply cross-field 
# stimulation in an attempt to induce re-entry. By default, the model is isotropic, but anisotropy can 
# be included by setting Orientation_type to anisotropic and setting the orientation directions as in 
# example 7. You may also want to define the diffusion parameters at run-time (example 6) as the defaults 
# may not be desired for your investigation.

# In this example, we’ll set D1 and dx directly, and set Dscale to 0.75 to promote sustained re-entry. 
# We’ll use the state write and read functionality to pace to stable state, and save the stable re-entry.

# First, we’ll pace a 1D model with matched D1 and dx to stable state at a relative rapid pacing rate. 
# We’ll write the state as an ave coupled cell, and pass a state reference to identify the specific tissue simulation.
./model_tissue_native Tissue_model re-entry Tissue_order 1D D1 0.25 dx 0.3 Dscale 0.75 BCL 300 Beats 100 Write_state ave State_Reference_write Re_entry_1D_pre_pace Reference Re_entry_example Results_Reference 1D_pre_pace

# We’ll then pace the 2D re-entry tissue model for a couple beats to settle it into a stable state 
# (you should ensure this is actually the case for any real investigations) and write the whole-tissue 
# state file (with another reference – this is now very important as we will also perform a whole-tissue 
# state write after inducing re-entry, and we want to preserve the separate state files):
./model_tissue_native Tissue_model re-entry D1 0.25 dx 0.3 Dscale 0.75 BCL 300 Beats 2 Read_state ave State_Reference_read Re_entry_1D_pre_pace Write_state On State_Reference_write Re_entry_2D_pre_pace Reference Re_entry_example Results_Reference 2D_pre_pace

# And finally, apply just one S1 beat, an S2 of 155 (or cycle through different intervals – which is why we wrote the 
# tissue state at stable S1 pacing). These parameters don’t actually lead to sustained re-entry – it self-terminates 
# within one second – so we will set the Total_time to 450, where a re-entrant circuit still persists, to illustrate 
# the state-read protocol next. 
./model_tissue_native Tissue_model re-entry D1 0.25 dx 0.3 Dscale 0.75 BCL 300 Beats 1 Total_time 450 S2 155 Read_state On State_Reference_read Re_entry_2D_pre_pace Write_state On State_Reference_write Re_entry_2D_stable_at_3secs Reference Re_entry_example Results_Reference S2_155

# And create vtk files of the data
./bin_to_vtk_tissue Reference Re_entry_example Results_Reference S2_155 Tissue_model re-entry start_time 0 end_time 450 interval 10

# We can then read this (semi) stable re-entrant state to start any number of simulations from a stable re-entry point (for 
# example to apply some Direct_modulation to test potential intervention to terminate the re-entry). Set beats to 0 as we 
# don’t want to apply a stimulus (we must still specify a BCL in order to read in the correct state file, which has the BCL 
# in its filename).
./model_tissue_native Tissue_model re-entry D1 0.25 dx 0.3 Dscale 0.75 BCL 300 Beats 0 Total_time 500 Read_state On State_Reference_read Re_entry_2D_stable_at_3secs IKr_scale 0.25 IKs_scale 0.25 Reference Re_entry_example Results_Reference Termination_attempt

# And create vtk files of the data
./bin_to_vtk_tissue Reference Re_entry_example Results_Reference Termination_attempt Tissue_model re-entry start_time 0 end_time 500 interval 10

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


