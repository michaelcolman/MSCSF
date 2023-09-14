#!/bin/sh

# We can also use a controllable dynamic SRF model
# Some presets Psets have been included, so we will first run one of these (set in lib/Spontaneous_release_functions.cpp)
# And then demonstrate how to set the CaSR dependence using the User_control approach

# We will use the CaSR loaded state files from the previous example
./model_single_0D Model minimal BCL 400 ISO 1 Beats 1 Read_state On state_reference_read 5beats SRF_mode Dynamic SRF_model General SRF_Pset General_1p10 Reference Dynamic_controllable Results_Reference General_800
./model_single_0D Model minimal BCL 400 ISO 1 Beats 1 Read_state On state_reference_read 50beats SRF_mode Dynamic SRF_model General SRF_Pset General_1p10 Reference Dynamic_controllable Results_Reference General_1000
./model_single_0D Model minimal BCL 400 ISO 1 Beats 1 Jup_scale 2 Read_state On state_reference_read 50beats_Jup SRF_mode Dynamic SRF_model General SRF_Pset General_1p10 Reference Dynamic_controllable Results_Reference General_1200

# Now we will control the CaSR dependence directly at run-time.
# Again, we will use a settings file to make this tidier
# We must pass:
#   SRF_Dyn_PSCR_threshold      CaSR threhsold above which SCRE occurs
#   SRF_Dyn_CaSR_max            CaSR value at which SCRE dynamcis converge
#   SRF_Dyn_CaSR_Prange         Range of CaSR for which probability goes from 0 to 1
#   SRF_Dyn_ti_sep_max          ti_sep at CaSR_max
#   SRF_Dyn_ti_sep_min          ti_sep at CaSR_threshold
#   SRF_Dyn_ti_width_max        ti dist width at CaSR_max
#   SRF_Dyn_ti_width_min        ti dist width at CaSR_threshold
#   SRF_Dyn_MD_max              MD at CaSR_max
#   SRF_Dyn_MD_min              MD at CaSR_threshold
#   SRF_Dyn_duration_width_max  duration dist width at CaSR_max
#   SRF_Dyn_duration_width_min  duration dist width at CaSR_threhold
#   SRF_Dyn_H                   power of CaSR dependent functions

./model_single_0D Settings_file SF_Dyn_User_control.txt Model minimal BCL 400 ISO 1 Beats 1 Read_state On state_reference_read 5beats SRF_mode Dynamic SRF_model General SRF_Pset User_control Reference Dynamic_controllable Results_Reference User_control_800
./model_single_0D Settings_file SF_Dyn_User_control.txt Model minimal BCL 400 ISO 1 Beats 1 Read_state On state_reference_read 50beats SRF_mode Dynamic SRF_model General SRF_Pset User_control Reference Dynamic_controllable Results_Reference User_control_1000 
./model_single_0D Settings_file SF_Dyn_User_control.txt Model minimal BCL 400 ISO 1 Beats 1 Jup_scale 1 Read_state On state_reference_read 50beats_Jup SRF_mode Dynamic SRF_model General SRF_Pset User_control Reference Dynamic_controllable Results_Reference User_control_1200

# And we can inspect the results and distributions exactly as before
# Have a play with changing the settings and seeing what happens at different CaSRs etc

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


