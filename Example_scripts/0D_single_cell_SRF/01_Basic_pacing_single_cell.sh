#!/bin/sh

# Control pacing using the 0D moddel

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hAM_ORD_s

BCL=$(printf "%04d" 500)
./model_single_0D Model $model BCL $BCL Beats 200 Reference Control_pacing_${model}

# Output data will be in "Outputs_0Dcell"
# Further to "Results/Currents.dat", which has time courses for the AP and ion currents,
# "Results/CRU.dat" contains the intracellular Ca2+ concentrations and RyR, LTCC and SERCA dynamics.
# (exactly as with 3D cell models)

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


