#!/bin/sh

# Example of using state references to identify specific state files where the automatic filenaming
# won't be unique to every simulation.

# The following settings are automatically encoded into the name of any single-cell state file:
#   Model, BCL, celltype, ISO, ACh, Remodelling, Agent, Mutation, environment
# This means that any other conditions/settings (ISO_model, Remodelling_prop, any direct control modulation)
# are not encoded automatically - state files with different settings here will overwrite each other.
# In order to distinguish and identify specific state files which are not differentiated by the automatically
# encoded conditions, you can use a State_reference (for reading or writing) - this appends the reference
# to the state files.

# example - write
./model_single_native Model hAM_CRN Beats 100 ISO 1 Write_state On Reference State_reference_example Results_Reference hAM_CRN_ISO_control_write
./model_single_native Model hAM_CRN Beats 100 ISO 1 IK1_scale 1.5 Write_state On State_Reference_write IK1_scale_1p5 Reference State_reference_example Results_Reference hAM_CRN_ISO_IK1_write

# now, running:
./model_single_native Model hAM_CRN Beats 5 ISO 1 IK1_scale 1.5 Read_state On Reference State_reference_example Results_Reference hAM_CRN_ISO_control_read
# will read in the state written by ./model_single_native Model hAM_CRN Beats 100 ISO 1 Write_state On, despite the presence
# of "IK1_scale 1.5" - which is NOT what is desired.
# Whereas running
./model_single_native Model hAM_CRN Beats 5 ISO 1 IK1_scale 1.5 Read_state On State_Reference_read IK1_scale_1p5 Reference State_reference_example Results_Reference hAM_CRN_ISO_IK1_read
# will read in the state file created when running IK1_scale 1.5, via looking for the file appended with "IK1_scale_1p5"
# Plot: "Results_hAM_CRN_ISO_control_write/Currents.dat" ($1-100000) vs 2,  "Results_hAM_CRN_ISO_IK1_write/Currents.dat" ($1-100000) vs 2,
# "Results_hAM_CRN_ISO_control_read/Currents.dat" 1 vs 2, "Results_hAM_CRN_ISO_IK1_read/Currents.dat" 1 vs 2
# Note the reference you pass is completely up to you and doesn't have to link to the additional modulation - use
# sensible references to idnetify and distnguish your simulations.

# Another reason to use state references is to output the state at different times during a simulation
# As the state files are named by the condition (which is constant in this example),
# and written only at the end of the simulation, in order to do this you need to use
# state references. These provide specific references for the state files to distinguish
# them within files associated with the same conditions.
# An example may be if you want to know what happens over time when a modifier is applied (e.g. ISO)
# during which it is possible pro-arrhythmic behaviour is transient or temporally dependent (e.g. EADs/DADs),
# and which intervention may have the best effect and at what time (i.e. different behaviours over the evolutio).
# The example here has no such interesting behaviour, but shows how to output state files at different points to
# read in and apply a test intervention applied at that time.

model=$(printf "hAM_CRN") # choose model here: mininmal, hAM_CRN, hAM_GB, dAM_VA etc
BCL=$(printf "%04d" 500)

# Pace with ISO, writing state files at different time points and identifyng with references
for beats in 5 10 30 100
do
    ./model_single_native Model $model BCL $BCL Beats $beats ISO 1 Write_state On State_Reference_write beats_$beats Reference State_reference_example Results_Reference Pre_pace_beats_${beats} 
done

# Read state and apply modulation test
# So this will read in the state with reference "beats_$beats", which was written at that number of beats into the sim.
for beats_start in 5 10 30 100
do
    for i in `seq 0.5 0.25 1.5`
    do
        mod=$(printf "%0.2f" $i)
        ./model_single_native Model $model BCL $BCL Beats 5 ISO 1 Read_state On State_Reference_read beats_${beats_start} Jup_scale $mod ICaL_scale $mod Reference State_reference_example Results_Reference Beats_start_${beats_start}_intervention_Jup_ICaL_${mod}
    done
done



# Output Cai and using 1.0 as control, show beats 5 vs 10 different etc

# plot "Outputs_single_native_State_reference_example/Results_Beats_start_5_intervention_Jup_ICaL_1.25/Currents.dat" u 1:34 w l, 
# "Outputs_single_native_State_reference_example/Results_Beats_start_10_intervention_Jup_ICaL_1.25/Currents.dat" u 1:34 w l, 
# "Outputs_single_native_State_reference_example/Results_Beats_start_10_intervention_Jup_ICaL_1.00/Currents.dat" u 1:34 w l, 
# "Outputs_single_native_State_reference_example/Results_Beats_start_5_intervention_Jup_ICaL_1.00/Currents.dat" u 1:34 w l
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


