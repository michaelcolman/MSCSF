#!/bin/sh

# Simplest example of the spontaneous release functions, using direct control presets
# "Direct_Control" simply means the SRF parameters will be sampled from a single, defined
# distribution for each beat. We will use the preset model "General_1" for this first example

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hAM_ORD_s

# First, just one beat so see some SRF
BCL=$(printf "%04d" 500)
./model_single_0D Model $model BCL $BCL Beats 1 SRF_mode Direct_Control SRF_Pset General_1 Reference Direct_control Results_Reference 1_beat

# And if we look at "Outputs_0Dcell_Direct_control/Results_1_beat/CRU.dat", we can see:
#   A DAD (colums 1 vs 2 - voltage)
#   The spontaneous calcium transient underlying the DAD (columns 1 vs 5 - Cai)
#   And the open RyR SRF waveform (columns 1 vs 9 - NRyR_O)

# This implementation can be used with multiple beats - the distribution is the same each time.
# Lets apply a few beats, one at a BCL faster than the likley fastest ti, and one slower
BCL=$(printf "%04d" 350)
./model_single_0D Model $model BCL $BCL Beats 3 SRF_mode Direct_Control SRF_Pset General_1 Reference Direct_control Results_Reference 3_beats_BCL_${BCL}

BCL=$(printf "%04d" 1500)
./model_single_0D Model $model BCL $BCL Beats 3 SRF_mode Direct_Control SRF_Pset General_1 Reference Direct_control Results_Reference 3_beats_BCL_${BCL}

# And we can see that when the BCL is slower, we get SCRE on each beat (and it has variability, even with the same distribution)

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


