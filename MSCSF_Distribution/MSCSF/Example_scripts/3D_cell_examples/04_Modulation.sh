#!/bin/sh

# Pacing with different modulations.
# One key aspect here is the different ways to control RyR and LTCC activity:
#   Change in channel expression: ICaL_scale and Jrel_scale will scale NRyR and NLTCC (ave per dyad)
#   Change in channel activity: LTCC_Po and RyR_Po will scale the transition rate to the open, active channel state
#       (LTCC this is a voltage independent transition from the voltage activation gate; for RyRs, it scales the transition
#       rate from C to O i.e. sensitivity to CaDS). 
#       These control the model variables GLTCC_kva1_va2 and GRyR_kCO.
# We can do both types of modification simultaneously (e.g. remodelling reduces channel numbers, but increases RyR sensitivity).

# As an example, let's compare what happens when we increase NLTCC and the ICaL open rate - both increase channel activity, different ways
# Check the screen settings outputs for NLTCC to see how it changes according to GCaL
# We'll turn voltage clamp on (which for the 3D cell model will clamp ICaL stochastically using spatial CICR) so we can compare peak current
# (Note that with the testing Sim_cell_size, the ICaL trace will be very noisy)
./model_single_3D Model minimal Beats 10 Sim_cell_size testing Vclamp On Reference Modulation_example Results_Reference Control_ICAL
./model_single_3D Model minimal Beats 10 Sim_cell_size testing Vclamp On ICaL_scale 2 Reference Modulation_example Results_Reference NLTCC
./model_single_3D Model minimal Beats 10 Sim_cell_size testing Vclamp On LTCC_Po 2 Reference Modulation_example Results_Reference LTCC_Po

# In "Outputs_Modulation_example/Results_{Control_ICAL/NLTCC/LTCC_Po}/Parameters" are now data for ICaL voltage clamp
# To plot/analyse the raw current data, use ICaL_Vclamp_traces.dat - plot time (column 1) against voltage (2), current (3) or Cai (4)
# To plot IV relationships, use ICaL_IV.dat (column 1 vs 2 for voltage and current)
# Look at "Outputs_Modulation_example/Results_{Control_ICAL/NLTCC/LTCC_Po}/CRU.dat" using columns:
#   1 - time; 1 - AP; 5 - Cai; 16 - JCaL
# And "../Currents.dat" using 1 (time) and 12 (ICaL)
# Note that a scale factor of 2 for ICaL_scale and LTCC_Po do not have identical effects; if you are matching data based on
# current measurements, you will need to set the scale to give the correct magnitude.

# Apply multiple:
# Increase SERCA (Jup_scale), increase activity of ICaL but not channel numbers (LTCC_Po), reduce numbers of RyRs (Jrel_scale) but increase open rate (RyR_Po)
./model_single_3D Model minimal Beats 10 Sim_cell_size testing Jup_scale 2 LTCC_Po 1.5 Jrel_scale 0.5 RyR_Po 2.5 Reference Modulation_example Results_Reference Multiple


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
