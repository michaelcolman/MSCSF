#!/bin/sh

# We will still use the Direct_Control mode, but this time we will apply the User_Control
# Pset, which requires all the distribution parameters to be set at run-time.
# We will use settings files to make it neater; you certainly don't have to.

# Settings which need to be passed:
#   SRF_DC_PSCRE        - probability of release (per excitation)
#   SRF_DC_CF           - Cumulative frequency at ti separation point
#   SRF_DC_ti_sep       - ti at this point
#   SRF_DC_ti_W1        - width of ti disitritbion to left of ti_sep
#   SRF_DC_ti_W2        - width of ti distribution to right of ti_sep
#   SRF_DC_MD           - Median distribution
#   SRF_DC_duration_W   - width of duration distribution (*optional; we set from MD is not passed)SRF_DC_PSCRE

model=$(printf "minimal") # choose model here: mininmal, hAM_CAZ_s, hAM_ORD_s
BCL=$(printf "%04d" 1000)

# Let's run a couple of very different distributions to see their effect.
# Dist 1 settings file: "SF_DC_User_control_low_SCRE.txt" contains:
# SRF_DC_PSCRE    0.5
# SRF_DC_CF       0.4
# SRF_DC_ti_sep   1000
# SRF_DC_ti_W1    400
# SRF_DC_ti_W2    700
# SRF_DC_MD       500
# SRF_DC_duration_W   300

./model_single_0D Settings_file SF_DC_User_control_low_SCRE.txt Model $model BCL $BCL Beats 1 SRF_mode Direct_Control SRF_Pset User_control Reference DC_User_control Results_Reference low_SCRE

# Dist 2 settings file: "SF_DC_User_control_high_SCRE.txt" contains:
# SRF_DC_PSCRE    1
# SRF_DC_CF       0.4
# SRF_DC_ti_sep   500
# SRF_DC_ti_W1    100
# SRF_DC_ti_W2    200
# SRF_DC_MD       150
# SRF_DC_duration_W   50

./model_single_0D Settings_file SF_DC_User_control_high_SCRE.txt Model $model BCL $BCL Beats 1 SRF_mode Direct_Control SRF_Pset User_control Reference DC_User_control Results_Reference high_SCRE

# The effect of these settings is clear from the results: dist 2 leads to much larger spontaneous calcium transients, 
# with some eliciting a triggered AP (run a few simlations - append Results_Reference with a run number if you like).
# We can also view/plot the distributions directly:
#   plot: Outputs_0Dcell_DC_User_control/Parameters_low_SCRE/SRF_distributions/ti_distribution.dat using 1 and 2, 
#       and compare it to Outputs_0Dcell_DC_User_control/Parameters_high_SCRE/SRF_distributions/ti_distribution.dat
#   And similar for "duration_distribution" and "NRyRo_peak_distribution"


# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


