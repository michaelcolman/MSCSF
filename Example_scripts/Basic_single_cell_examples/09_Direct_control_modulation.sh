#!/bin/sh

# Example to apply multiple modulations including direct control of modifier parameters
# Example has one varaible being varied, within the condition of multiple others applied constantly
# Celltype PV, ISO 0.5
# Apply direct control:
#   INa scale 0.5
#   SERCA / Jup scale 0.5
#   increase the time constant of the activation gate of Ito by factor 5
#   shift the voltage dependence of the inactivation gate for ICaL by -10 mV
#   scale the gradient parameter of IKur activation steady state by 2
#   shift the voltage-dependence of the time constant for Ito inactivation by +10 mV
#   scale ICaL between 0 and 2 

for model in hAM_CRN hAM_WL_CRN
do
    # Produce data without modifications
    ./model_single_native Model ${model} BCL 1000 Beats 10 Celltype PV Reference Direct_modulation Results_Reference Model_${model}_control
    for i in `seq 0 0.5 2`
    do
        ICAL=$(printf "%0.2f" $i)
        ./model_single_native Model ${model} BCL 1000 Beats 10 Celltype PV ISO 0.5 INa_scale 0.5 Jup_scale 0.5 Ito_va_tau_scale 5 ICaL_vi_ss_shift -10 IKur_va_ss_kscale 2 Ito_vi_tau_shift 10 ICaL_scale $ICAL Reference Direct_modulation Results_Reference Model_${model}_ICAL_${ICAL}
    done
done

# This will create an Outputs directory associated with each model: "Outputs_single_native_Direct_modulation"
# And Results directorties for each model and ICaL scale: "Results_Model_Y_ICAL_X", "Results_Model_Y_control"
# There is no need to provide a reference noting the other modulation variables because these are constant, 
# the most recent simulation is valid for all modified variables other than ICaL scale.
# Therefore, Log.dat and "Outputs_single_native_Direct_modulation/Settings.dat" contains a record of the other modifications.
# Further to the results data, compare the parameters list in “../Parameters_Without_direct_modulation” and 
# “../Parameters_With_direct_modulation” and the voltage-dependent steady-state of ICaL_vi, IKur_va and 
# time constants of Ito_va and vi in in “../Parameters_X/Ix.dat” to visualise the effect of the voltage shift, 
# time-constant scaling and kscale:
# diff Parameters_Model_hAM_CRN_control/Magnitude_parameters.dat Parameters_Model_hAM_CRN_ICAL_1.00/Magnitude_parameters.dat
# "../Parameters_Model_hAM_CRN_control/ICaL.dat" 1 vs 4, "../Parameters_Model_hAM_CRN_ICAL_1.00/ICaL.dat" 1 vs 4
# "../Parameters_Model_hAM_CRN_control/IKur.dat" 1 vs 2, "../Parameters_Model_hAM_CRN_ICAL_1.00/IKur.dat" 1 vs 2
# "../Parameters_Model_hAM_CRN_control/Ito.dat" 1 vs 3, "../Parameters_Model_hAM_CRN_ICAL_1.00/Ito.dat" 1 vs 3
# "../Parameters_Model_hAM_CRN_control/Ito.dat" 1 vs 5, "../Parameters_Model_hAM_CRN_ICAL_1.00/Ito.dat" 1 vs 5
# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
