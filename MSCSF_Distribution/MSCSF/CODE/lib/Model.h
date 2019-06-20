// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Global Model file, which contains functions =  //
// which: 1) selects appropriate set and solve functions; =  //
// 2) sets global and selects appropriate model-specific ==  //
// heterogeneity and modulation; 3) determines stimulus and  //
// measurement variables; 4) contains the Luo-Rudy ========  //
// implementation of INa: =================================  //
// Luo-Rudy Circulation Research. 1991;68:1501-1526 =======  //
// ========================================================  //
// GNU 3 LICENSE TEXT =====================================  //
// COPYRIGHT (C) 2016-2019 MICHAEL A. COLMAN ==============  //
// THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT =  //
// AND/OR MODIFY IT UNDER THE TERMS OF THE GNU GENERAL ====  //
// PUBLIC LICENSE AS PUBLISHED BY THE FREE SOFTWARE =======  //
// FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT YOUR  //
// OPTION) ANY LATER VERSION. =============================  //
// THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE=  //
// USEFUL, BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE =====  //
// IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS FOR A ===  //
// PARTICULAR PURPOSE.  SEE THE GNU GENERAL PUBLIC LICENSE=  //
// FOR MORE DETAILS. ======================================  //
// YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL =====  //
// PUBLIC LICENSE ALONG WITH THIS PROGRAM.  IF NOT, SEE ===  //
// <https://www.gnu.org/licenses/>. =======================  //
// ========================================================  //
// ADDITIONAL LICENSE TEXT ================================  //
// THIS SOFTWARE IS PROVIDED OPEN SOURCE AND MAY BE FREELY=  //
// USED, DISTRIBUTED AND UPDATED, PROVIDED: ===============  //
//  (i) THE APPROPRIATE WORK(S) IS(ARE) CITED. THIS =======  //
//      PERTAINS TO THE CITATION OF COLMAN 2019 PLOS COMP =  //
//      BIOL (FOR THIS IMPLEMTATION) AND ALL WORKS ========  //
//      ASSOCIATED WITH THE SPECIFIC MODELS AND COMPONENTS=  //
//      USED IN PARTICULAR SIMULATIONS. IT IS THE USER'S ==  //
//      RESPONSIBILITY TO ENSURE ALL RELEVANT WORKS ARE ===  //
//      CITED. PLEASE SEE FULL DOCUMENTATION AND ON-SCREEN=  //
//      DISCLAIMER OUTPUTS FOR A GUIDE. ===================  //
//  (ii) ALL OF THIS TEXT IS RETAINED WITHIN OR ASSOCIATED=  //
//      WITH THE SOURCE CODE AND/OR BINARY FORM OF THE ====  //
//      SOFTWARE. =========================================  //
// ========================================================  //
// ANY INTENDED COMMERCIAL USE OF THIS SOFTWARE MUST BE BY   //
// EXPRESS PERMISSION OF MICHAEL A COLMAN ONLY. IN NO EVENT  //
// ARE THE COPYRIGHT HOLDERS LIABLE FOR ANY DIRECT, =======  //
// INDIRECT INCIDENTAL, SPECIAL, EXEMPLARY OR CONSEQUENTIAL  //
// DAMAGES ASSOCIATED WITH USE OF THIS SOFTWARE ===========  //
// ========================================================  //
// THIS SOFTWARE CONTAINS IMPLEMENTATIONS OF MODELS AND ===  //
// COMPONENTS WHICH I (MICHAEL COLMAN) DID NOT DEVELOP.====  //
// ALL OF THESE COMPONENTS HAVE BEEN CODED FROM PROVIDED ==  //
// SOURCE CODE OR INFORMATION IN THE PUBLICATIONS. ========  //
// I CLAIM NO RIGHTS OR INTELLECTUAL PROPERTY OWNERSHIP ===  //
// FOR THESE MODELS AND COMPONENTS, OTHER THAN THEIR ======  //
// SPECIFIC IMPLEMENTATION IN THIS CODE PACKAGE. FURTHER TO  //
// THE ABOVE STATEMENT, ANY INDTENDED COMMERCIAL USE OF ===  //
// THOSE COMPONENTS MUST BE BY EXPRESS PERMISSION OF THE ==  //
// ORIGINAL COPYRIGHT HOLDERS. ============================  //
// WHERE IMPLEMENTED FROM PROVIDED CODE, ANY DISCLAIMERS ==  //
// PRESENT IN THE ORIGINAL CODE HAVE BEEN RETAINED IN THE =  //
// RELEVANT FILE. =========================================  //
// ========================================================  //
// Contact: m.a.colman@leeds.ac.uk ========================  //
// For updates, corrections etc, please check: ============  //
// 1. http://physicsoftheheart.com/ =======================  //
// 2. https://github.com/michaelcolman ====================  //
// ========================================================  //

#ifndef MODEL_H
#define MODEL_H

#include "Structs.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Common functions  ========================================================\\|
// Set parameters functions - choses which set parameters function to call
void set_parameters_native(Cell_parameters *p, char const *Model); 
void set_parameters_spatial_Ca(Cell_parameters *p, char const *Model);

// Initial conditions - choses which initial conditions function to call
void initial_conditions_native(State_variables *s, Cell_parameters p, char const *Model);

// Compute and update ionic currents - choses which total current function set to call
void compute_model_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_model_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);

// Output functions for checking
void compute_and_output_current_functions(Cell_parameters p, Model_variables *var, char const *directory);

// Stimulus current functions
void stimulus_setup(Cell_parameters p, Model_variables *var, double dt, int BCL, int S2, int Paced_time);
void compute_Istim(Cell_parameters p, Model_variables *var, double Paced_time, double S2_time, double time, int time_int);

// Current modification variables - sets global modification and selects apporpriate specific functions 
void set_heterogeneity_and_modulation_native(Cell_parameters *p);
void update_heterogeneity_and_modulation_integrated(Cell_parameters *p);

// Global or common het and modulation to multiple models
void set_global_Agents(Cell_parameters *p);
void set_celltype_hAM(Cell_parameters *p);
void set_ISO_hAM(Cell_parameters *p);
void set_global_remodelling(Cell_parameters *p);
void set_remodelling_hAM(Cell_parameters *p);
void set_mutation_hAM(Cell_parameters *p);
void set_ACh_global(Cell_parameters *p);
void set_spatial_gradient(Cell_parameters *p);
void set_MODIFIER_X_Y(Cell_parameters *p);

// Reversal potentials
void compute_reversal_potentials(Cell_parameters p, Model_variables *var, State_variables *s);

// Excitation properties / measurements
void determine_excitation_state(Model_variables *var, double Vm, double time);
void determine_excitation_state_integrated_0D(Model_variables *var, double Vm, double time, double *Ca_JSR_t_ex, double Ca_JSR, double *dyad_SRF_prop_active, double srf_SRF_prop_active, int *srf_init, int *srf_set, const char *SRF_Mode);
void calculate_measurement_properties(Model_variables *var, double Vm1, double Vm2, double time, double dt, double APD_threshold, double CaT, double CaSR);
void calculate_flux_integrals(Cell_parameters p, Model_variables *var, double J_SERCA, double J_NCX, double J_rel, double J_LTCC);

// Voltage clamp
void run_voltage_clamp(Cell_parameters p, Model_variables *var, State_variables *s, char const *directory, double dt);

// Frequently used functions
double rush_larsen(double y, double ss, double tau, double dt);
double sigmoid(double V, double V_half, double k);

// Formulation of INa from Luo-Rudy 1991, used in multiple models
void set_INa_LR_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_INa_LR(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_INa_LR(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End Common functions  ====================================================//|

// Minimal model functions ==================================================\\|
// Parameters and specific settings
void set_parameters_native_minimal(Cell_parameters *p);
void initial_conditions_native_minimal(State_variables *s, Cell_parameters p);
void set_het_mod_minimal(Cell_parameters *p);
void set_celltype_native_minimal(Cell_parameters *p);
void set_modulation_ISO_native_minimal(Cell_parameters *p);
void set_modulation_Agent_native_minimal(Cell_parameters *p);
void set_modulation_Remodelling_native_minimal(Cell_parameters *p);
void set_modulation_Mutation_native_minimal(Cell_parameters *p);
void set_modulation_ACh_minimal(Cell_parameters *p);

// Solve model parent functions
void compute_model_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_model_minimal_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_minimal_native(Cell_parameters p, Model_variables *var, double Vm);
void update_gating_variables_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_Itot_minimal_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_minimal(Cell_parameters p, Model_variables *var, char const * directory);

// Specific current functions
// Ip0d
void set_Ip0d_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ip0d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ip0d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ip1r
void set_Ip1r_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ip1r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ip1r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ip2d
void set_Ip2d_rates(Cell_parameters p, Model_variables *var, double Vm);
void set_Ip2d_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_Ip2d_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_Ip2d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ip2d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ip2r
void set_Ip2r_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ip2r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ip2r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ip3r
void set_Ip3r_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ip3r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ip3r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ip4r
void set_Ip4r_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_Ip4r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End Minimal model functions ==============================================//|

// hAM_WL model functions ==================================================\\|
// Parameters and specific settings
void update_parameters_native_hAM_WL(Cell_parameters *p);					
void initial_conditions_native_hAM_WL(State_variables *s, Cell_parameters p);
void set_het_mod_hAM_WL(Cell_parameters *p);
void set_celltype_native_hAM_WL(Cell_parameters *p);
void set_modulation_ISO_native_hAM_WL(Cell_parameters *p);
void set_modulation_Agent_native_hAM_WL(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_WL(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_WL(Cell_parameters *p);
void set_modulation_ACh_hAM_WL(Cell_parameters *p);

// Solve model parent functions
void compute_ICaL_hAM_WL_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICaL_hAM_WL_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

void compute_model_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_WL_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

void compute_and_output_current_functions_hAM_WL(Cell_parameters p, Model_variables *var, char const *directory);

// Specific current functions
// INa
void set_INa_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm);

// Ito
void set_Ito_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_WL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_WL_vi_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);

void set_ICaL_hAM_CRN_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_CRN_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_CRN_mWL_vi_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);

void set_ICaL_hAM_GB_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_GB_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_GB_mWL_vi_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);

void set_ICaL_hAM_NG_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_NG_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_NG_mWL_vi_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);

void update_gates_ICaL_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void update_gates_ICaL_hAM_WL_CRN_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void update_gates_ICaL_hAM_WL_GB_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void update_gates_ICaL_hAM_WL_NG_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);

void compute_ICaL_hAM_CRN_mWL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICaL_hAM_WL_CRN_bar(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double Cai);
void compute_ICaL_hAM_WL_GB_bar(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICaL_hAM_NG_mWL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void compute_IK1_hAM_WL_isolated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_IK1_hAM_WL_intact(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End hAM_WL model functions ===============================================//|

// Maleckar et al Ito and IKur hAM currents =================================\\|
void update_parameters_native_hAM_MT(Cell_parameters *p);
void initial_conditions_native_hAM_MT(State_variables *s, Cell_parameters p);

void set_het_mod_hAM_MT(Cell_parameters *p);
void set_celltype_native_hAM_MT(Cell_parameters *p);
void set_modulation_ISO_native_hAM_MT(Cell_parameters *p);
void set_modulation_Agent_native_hAM_MT(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_MT(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_MT(Cell_parameters *p);
void set_modulation_ACh_hAM_MT(Cell_parameters *p);


// Solve model functions
void compute_model_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_MT_native(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Ko);
void update_gating_variables_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hAM_MT(Cell_parameters p, Model_variables *var, char const * directory);

// Ito
void set_Ito_hAM_MT_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_hAM_MT_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End Maleckar et al Ito and IKur hAM currents =============================//|

// hAM_GB model functions (Grandi-Bers lab model, 2011) =====================\\|
// Parameters and specific settings
void set_parameters_native_hAM_GB(Cell_parameters *p);
void initial_conditions_native_hAM_GB(State_variables *s, Cell_parameters p);
void set_het_mod_hAM_GB(Cell_parameters *p);
void set_modulation_ISO_native_hAM_GB(Cell_parameters *p);
void set_modulation_Agent_native_hAM_GB(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_GB(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_GB(Cell_parameters *p);
void set_modulation_ACh_hAM_GB(Cell_parameters *p);

// Solve model functions
void compute_model_hAM_GB_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_GB_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hAM_GB_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_GB_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hAM_GB(Cell_parameters p, Model_variables *var, char const * directory);
void set_celltype_native_hAM_GB(Cell_parameters *p);

// INa
void compute_INa_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaL
void set_INaL_hAM_GB_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_INaL_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_INaL_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ito  || MT formulation
void set_Ito_hAM_GB_rates(Cell_parameters p, Model_variables *var, double Vm);

// ICaL
void set_ICaL_hAM_GB_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_GB_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_GB_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
double compute_ICaL_bar_hAM_GB(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Cao);
void compute_ICaL_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
double compute_ICaL_bar_Na_hAM_GB(Cell_parameters p, Model_variables *var, double Vm, double Nai, double Nao);
double compute_ICaL_bar_K_hAM_GB(Cell_parameters p, Model_variables *var, double Vm, double Ki, double Ko);

// IKur || MT formulation

// IKr
void set_IKr_hAM_GB_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_hAM_GB_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_hAM_GB_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_IK1_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INCX
void compute_INCX_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaK
void compute_INaK_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IClCa  | IClb
void compute_IClCa_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_IClb_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// PMCA
void compute_ICaP_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background 
void compute_INab_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_IKb_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Homeostasis
void comp_homeostasis_hAM_GB(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);

// Modded
void initial_conditions_native_hAM_GB_modded(State_variables *s, Cell_parameters p);
void compute_model_hAM_GB_modded_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_GB_modded_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hAM_GB_modded_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void update_parameters_native_hAM_GB_modded(Cell_parameters *p);
void compute_Itot_hAM_GB_modded_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// 5 state Markov IcaL
void set_ICaL_5sm_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_states_ICaL_5sm(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_hAM_GB_5sm(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// end hAM_GB model functions ===============================================//|

// hAM_CRN model functions ==================================================\\|
// Parameters and specific settings
void set_parameters_native_hAM_CRN(Cell_parameters *p);
void initial_conditions_native_hAM_CRN(State_variables *s, Cell_parameters p);
void set_het_mod_hAM_CRN(Cell_parameters *p);
void set_modulation_ISO_native_hAM_CRN(Cell_parameters *p);
void set_modulation_Agent_native_hAM_CRN(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_CRN(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_CRN(Cell_parameters *p);
void set_modulation_ACh_hAM_CRN(Cell_parameters *p);

// Solve model parent functions
void compute_model_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_CRN_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hAM_CRN(Cell_parameters p, Model_variables *var, char const * directory);
void set_celltype_native_hAM_CRN(Cell_parameters *p);

// Specific current functions
// INa
void set_INa_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm);

// Ito
void set_Ito_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_CRN_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_CRN_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_hAM_CRN_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_IK1_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INCX
void compute_INCX_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaK
void compute_INaK_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// PMCA
void compute_ICaP_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background
void compute_INab_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Homeostasis
void comp_homeostasis_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
// End hAM_CRN model functions ===============================================//|

// hAM_NG model functions ====================================================\\|
// Parameters and specific settings
void set_parameters_native_hAM_NG(Cell_parameters *p);
void initial_conditions_native_hAM_NG(State_variables *s, Cell_parameters p);
void set_het_mod_hAM_NG(Cell_parameters *p);
void set_celltype_native_hAM_NG(Cell_parameters *p);
void set_modulation_ISO_native_hAM_NG(Cell_parameters *p);
void set_modulation_Agent_native_hAM_NG(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_NG(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_NG(Cell_parameters *p);
void set_modulation_ACh_hAM_NG(Cell_parameters *p);

// Solve model parent functions
void compute_model_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_NG_native(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Ko);
void update_gating_variables_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hAM_NG(Cell_parameters p, Model_variables *var, char const * directory);
void set_celltype_native_hAM_NG(Cell_parameters *p);

// Specific current functions
// INa
void set_INa_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm);
void compute_INa_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ito
void set_Ito_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hAM_NG_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_NG_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_hAM_NG_variables(Cell_parameters p, Model_variables *var, double Vm, double Ko);
void compute_IK1_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INCX
void compute_INCX_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaK
void compute_INaK_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// PMCA
void compute_ICaP_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background
void compute_INab_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Homeostasis
void comp_homeostasis_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
// End hAM_NG model functions ===============================================//|

// ratAM_CAL model functions ====================================================\\|
// Parameters and specific settings
void set_parameters_native_ratAM_CAL(Cell_parameters *p);
void initial_conditions_native_ratAM_CAL(State_variables *s, Cell_parameters p);
void set_het_mod_ratAM_CAL(Cell_parameters *p);
void set_celltype_native_ratAM_CAL(Cell_parameters *p);
void set_modulation_ISO_native_ratAM_CAL(Cell_parameters *p);
void set_modulation_Agent_native_ratAM_CAL(Cell_parameters *p);
void set_modulation_Remodelling_native_ratAM_CAL(Cell_parameters *p);
void set_modulation_Mutation_native_ratAM_CAL(Cell_parameters *p);

// Solve model parent functions
void compute_model_ratAM_CAL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_ratAM_CAL_native(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Ko);
void update_gating_variables_ratAM_CAL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_ratAM_CAL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_ratAM_CAL(Cell_parameters p, Model_variables *var, char const * directory);
void set_celltype_native_ratAM_CAL(Cell_parameters *p);

// Specific current functions
// INa
void set_INa_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm);
void compute_INa_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ito
void set_Ito_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_ratAM_CAL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_ratAM_CAL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_ratAM_CAL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_ratAM_CAL_variables(Cell_parameters p, Model_variables *var, double Vm, double Ko);
void compute_IK1_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INCX
void compute_INCX_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaK
void compute_INaK_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// PMCA
void compute_ICaP_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background
void compute_INab_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Homeostasis
void comp_homeostasis_ratAM_CAL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
// End ratAM_CAL model functions ============================================//|

// ORD hVM simplified model functions =======================================\\|
void set_parameters_native_hVM_ORD_simple(Cell_parameters *p);
void update_parameters_native_hVM_ORD_simple(Cell_parameters *p);
void initial_conditions_native_hVM_ORD_simple(State_variables *s, Cell_parameters p);
void set_het_mod_hVM_ORD_simple(Cell_parameters *p);
void set_modulation_ISO_native_hVM_ORD_simple(Cell_parameters *p);
void set_modulation_Agent_native_hVM_ORD_simple(Cell_parameters *p);
void set_modulation_Remodelling_native_hVM_ORD_simple(Cell_parameters *p);
void set_modulation_Mutation_native_hVM_ORD_simple(Cell_parameters *p);
void set_celltype_native_hVM_ORD_simple(Cell_parameters *p);
void set_modulation_ACh_hVM_ORD_simple(Cell_parameters *p);

void compute_model_hVM_ORD_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hVM_ORD_simple_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hVM_ORD_simple_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hVM_ORD_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hVM_ORD_simple(Cell_parameters p, Model_variables *var, char const * directory);

// INaL
void set_INaL_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_INaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_INaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ito
void set_Ito_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_hVM_ORD_simple_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hVM_ORD_simple_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

void set_IKs_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

//IK1
void update_gates_IK1_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IK1_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background
void compute_IKb_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_INab_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End ORD hVM simplified model functions ===================================//|

// CAZ hAM simplified model =================================================\\|
void set_parameters_native_hAM_CAZ_simple(Cell_parameters *p);
void update_parameters_integrated_hAM_CAZ_simple(Cell_parameters *p);
void initial_conditions_native_hAM_CAZ_simple(State_variables *s, Cell_parameters p);
void set_het_mod_hAM_CAZ_simple(Cell_parameters *p);
void set_modulation_ISO_native_hAM_CAZ_simple(Cell_parameters *p);
void set_modulation_Agent_native_hAM_CAZ_simple(Cell_parameters *p);
void set_modulation_Remodelling_native_hAM_CAZ_simple(Cell_parameters *p);
void set_modulation_Mutation_native_hAM_CAZ_simple(Cell_parameters *p);
void set_celltype_native_hAM_CAZ_simple(Cell_parameters *p);
void set_modulation_ACh_hAM_CAZ_simple(Cell_parameters *p);

void compute_model_hAM_CAZ_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_hAM_CAZ_simple_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_hAM_CAZ_simple_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_hAM_CAZ_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, char const * directory);

// Ito
void set_Ito_hAM_CAZ_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_hAM_CAZ_simple_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_hAM_CAZ_simple_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);

// IKr
void set_IKr_hAM_CAZ_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_hAM_CAZ_simple_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_hAM_CAZ_simple_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_IK1_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INab
void compute_INab_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
// End CAZ hAM simplified model =============================================//|

// Varela-Aslanidi dog atrial myocyte (dAM_VA) ==============================\\|
void set_parameters_native_dAM_VA(Cell_parameters *p);
void update_parameters_integrated_dAM_VA(Cell_parameters *p);
void initial_conditions_native_dAM_VA(State_variables *s, Cell_parameters p);

void set_het_mod_dAM_VA(Cell_parameters *p);
void update_het_and_mod_dAM_VA_integrated(Cell_parameters *p);

void set_celltype_native_dAM_VA(Cell_parameters *p);
void update_celltype_integrated_dAM_VA(Cell_parameters *p);
void set_modulation_ISO_native_dAM_VA(Cell_parameters *p);
void set_modulation_Agent_native_dAM_VA(Cell_parameters *p);
void set_modulation_Remodelling_native_dAM_VA(Cell_parameters *p);
void set_modulation_Mutation_native_dAM_VA(Cell_parameters *p);
void set_modulation_ACh_dAM_VA(Cell_parameters *p);

void compute_model_dAM_VA_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_model_dAM_VA_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_dAM_VA_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_dAM_VA_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_dAM_VA_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_Itot_dAM_VA_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_dAM_VA(Cell_parameters p, Model_variables *var, char const *directory);

// Ito
void set_Ito_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_dAM_VA_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_dAM_VA_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKACh
void set_IKACh_dAM_VA_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKACh_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKACh_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_dAM_VA_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_IK1_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background and Ca2+ handling
void compute_INCX_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_INaK_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICaP_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_INab_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_IClb_Varela_dAM(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void comp_homeostasis_dAM_VA(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
// End Varela-Aslanidi dog atrial myocyte (dAM_VA) ==========================//|

// TEMPLATE FOR NEW MODEL ===================================================\\|
// Copy all of these + add new ones for a new model, renaming to appropriate
// Parameters and specific settings
void set_parameters_native_speciesCELL_MODEL(Cell_parameters *p);
void update_parameters_native_speciesCELL_MODEL(Cell_parameters *p);
void update_parameters_integrated_speciesCELL_MODEL(Cell_parameters *p);
void initial_conditions_native_speciesCELL_MODEL(State_variables *s, Cell_parameters p);
void set_het_mod_speciesCELL_MODEL(Cell_parameters *p);
void update_het_and_mod_speciesCELL_MODEL_integrated(Cell_parameters *p);
void set_modulation_ISO_native_speciesCELL_MODEL(Cell_parameters *p);
void set_modulation_Agent_native_speciesCELL_MODEL(Cell_parameters *p);
void set_modulation_Remodelling_native_speciesCELL_MODEL(Cell_parameters *p);
void set_modulation_Mutation_native_speciesCELL_MODEL(Cell_parameters *p);
void set_celltype_native_speciesCELL_MODEL(Cell_parameters *p);
void set_modulation_ACh_speciesCELL_MODEL(Cell_parameters *p);
void update_celltype_integrated_speciesCELL_MODEL(Cell_parameters *p);

// Solve model parent functions
void compute_model_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_model_speciesCELL_MODEL_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void set_gate_rates_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void update_gating_variables_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Itot_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_Itot_speciesCELL_MODEL_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_and_output_current_functions_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, char const * directory);

// Specific current functions
// INa
void set_INa_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_INa_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_INa_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaL
void set_INaL_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_INaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_INaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Ito
void set_Ito_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_Ito_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_Ito_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// ICaL
void set_ICaL_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai);
void set_ICaL_speciesCELL_MODEL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale);
void set_ICaL_speciesCELL_MODEL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale);
void update_gates_ICaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_ICaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKur
void set_IKur_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKur_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKur_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKr
void set_IKr_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKr_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKr_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IKs
void set_IKs_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm);
void update_gates_IKs_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
void compute_IKs_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// IK1
void set_IK1_speciesCELL_MODEL_variables(Cell_parameters p, Model_variables *var, double Vm);
void compute_IK1_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INCX
void compute_INCX_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// INaK
void compute_INaK_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// PMCA
void compute_ICaP_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Background
void compute_INab_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);
void compute_ICab_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm);

// Homeostasis
void comp_homeostasis_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt);
// EMD TEMPLATE FOR NEW MODEL ===============================================//|

#endif
