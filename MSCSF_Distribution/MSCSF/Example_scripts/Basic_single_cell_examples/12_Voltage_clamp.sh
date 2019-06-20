#!/bin/sh

# Apply voltage clamp (ICaL, Ito, IKur, IK1 currently implemented)

# Apply 0 beats as do not need pacing, or apply voltage clamp and then run simulation
./model_single_native Model hAM_CRN Beats 0 Vclamp On Reference Voltage_clamp_example

# In "Outputs_single_native_Voltage_clamp_example/Parameters" are now data for ICaL, Ito, IKur and IK1 following voltage clamp
# To plot/analyse the raw current data, use I{CaL/to}_Vclamp_trace{s/_Z}.dat
# To plot IV relationships, use {CaL/to/K1}_IV.dat
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


