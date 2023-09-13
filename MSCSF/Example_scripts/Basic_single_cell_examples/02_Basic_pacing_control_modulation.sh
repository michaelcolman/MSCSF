#!/bin/sh

# Pacing pacing at single cycle length using different models comparing control and different modulations (example has a remodelling and +/- ISO)

for model in hAM_CRN hAM_GB hAM_WL_CRN # will run for all three models
do
    for remodelling in none AF_Col_4 AF_GB
    do
        for i in `seq 0 0.5 1`
        do
            ISO=$(printf "%0.2f" $i)
            ./model_single_native Model ${model} BCL 1000 Beats 50 Remodelling ${remodelling} ISO ${ISO} Reference Control_modulation_${model} Results_Reference ISO_${ISO}_remodelling_${remodelling} 
        done
    done
done

# This will create an Outputs directory associated with each model: "Outputs_single_native_Control_modulation_model"
# And Results directorties for each condition: "Results_ISO_X_remodelling_Y"
# So, for example, for each model, can plot the AP/currents/CaT overlaid for each condition using "Results_ISO_X_remodelling_Y/Currents.dat"
# plot "Outputs_single_native_Control_modulation_model/Results_ISO_0_remodelling_none.dat" using columns 1 and 2, 
#           'Outputs_single_native_Control_modulation_model/Results_ISO_0_remodelling_AF_GB.dat" using columns 1 and 2 etc
# Can also compare outputs in "Outputs_single_native_Control_modulation_model/Parameters_ISO_X_remodelling_Y/"
# which contains parameter list, modification variable list and the voltage dependent variables (transition rates, steaty-states, time constants) 
# for each current as used in each simulation i.e. for plotting and/or confirming implementation of the modulation.
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
