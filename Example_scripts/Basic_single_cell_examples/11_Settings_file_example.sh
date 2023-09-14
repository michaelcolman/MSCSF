#!/bin/sh

# The example with multiple modulations including direct control involved many arguments to be passed
# In order to make our lives simpler regarding repeatedly running simulations with multiple constant and
# a few varied settings, we can use a settings file.
# A settings file can be used in combination with further command-line arguments to provide additional
# settings and/or overwrite setttings. This example simply provides additional settings (the Model, Results_Reference
# and ICaL_scale), but any option passed after the settings file will overwrite that option if in the file.
# Settings_file argument must be passed first.

# We basically copy the commands from script example 9 into a settings file "Settings_file_for_example_9_and_11.txt"
# including the number of options at the top of the file. Have a look. Note it has all options from example 9 other
# than Model, Results_Reference and ICaL_scale


for model in hAM_CRN hAM_WL_CRN
do
    # Produce data without modifications
    ./model_single_native Model ${model} BCL 1000 Beats 10 Celltype PV Reference Direct_modulation_settings_file Results_Reference Model_${model}_control
    for i in `seq 0 0.5 2`
    do
        ICAL=$(printf "%0.2f" $i)
        ./model_single_native Settings_file Settings_file_for_example_9_and_11.txt Model ${model} ICaL_scale $ICAL Results_Reference Model_${model}_ICAL_${ICAL}
    done
done

# Exaclty the same outputs as example 9 except with Outputs appended with "_settings_file"
# This will create an Outputs directory associated with each model: "Outputs_single_native_Direct_modulation"
# And Results directorties for each model and ICaL scale: "Results_Model_Y_ICAL_X", "Results_Model_Y_control"
# There is no need to provide a refernce noting the other modulation variables; because these are constant, 
# the most recent simulation is valid for all modified variables other than ICaL scale.
# Therefore, Log.dat and "Outputs_single_native_Direct_modulation/Settings.dat" contains a record of the other modifications.
# Furthermore, all of the parameters for the simulations are in FUNCTIONS OUT including plotting new 
# ICaL steady-states and Ito tau
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
