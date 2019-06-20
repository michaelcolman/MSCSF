#!/bin/sh

# Now using a dynamic model, fit to the 3D cell model behaviour
# As this model depends on CaSR, we will prepace the model to 3 different CaSRs so we can see
# the effect. We will apply ISO, and maybe some Jup_scale, to help load the CaSR
# We need to pass "SRF_model 3D_cell" to select from models fit to the 3D cell model (rather than general)

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hAM_ORD_s
BCL=$(printf "%04d" 400)

# recall that we need to set total time direclty to avoid leaving a 2sec quiescent period
# Below settings load CaSR to ~ 800, 1000 and 1200 mM
./model_single_0D Model minimal BCL ${BCL} ISO 1 NBeats 5 Total_time 2000 Write_state   On state_reference_write   5beats Reference Dynamic_3D_cell Results_Reference Pre_pace_5
./model_single_0D Model minimal BCL ${BCL} ISO 1 NBeats 50 Total_time 20000 Write_state   On state_reference_write   50beats Reference Dynamic_3D_cell Results_Reference Pre_pace_50
./model_single_0D Model minimal BCL ${BCL} ISO 1 NBeats 50 Total_time 20000 Jup_scale 2 Write_state   On state_reference_write   50beats_Jup Reference Dynamic_3D_cell Results_Reference Pre_pace_50_Jup

# Now, read in each and apply 1 beat, then apply the Dynamic SRF mode and different Psets
for pset in Control RSERCA_NCX RCRU
do
    ./model_single_0D Model minimal BCL ${BCL} ISO 1 Beats 1 Read_state On state_reference_read 5beats SRF_mode Dynamic SRF_model 3D_cell SRF_Pset ${pset} Reference Dynamic_3D_cell Results_Reference ${pset}_800
    ./model_single_0D Model minimal BCL ${BCL} ISO 1 Beats 1 Read_state On state_reference_read 50beats SRF_mode Dynamic SRF_model 3D_cell SRF_Pset ${pset} Reference Dynamic_3D_cell Results_Reference ${pset}_1000
    ./model_single_0D Model minimal BCL ${BCL} ISO 1 Beats 1 Jup_scale 2 Read_state On state_reference_read 50beats_Jup SRF_mode Dynamic SRF_model 3D_cell SRF_Pset ${pset} Reference Dynamic_3D_cell Results_Reference ${pset}_1200
done

# Notice how the SRF models depend on CaSR differently
# If we also passed "Remodelling RSERCA_NCX" into that case, we would also see different CaSR loading
# (would have to create new state files).
# Notice therefore that the SRF model can be set independent of cell model conditions

# We can also see what the SRF distributions at each CaSR for each model are
# in "../SRF_distributions/X_distribution_YmM.dat" and how the parameter vary with CaSR in "../CaSR_dependence.dat"


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


