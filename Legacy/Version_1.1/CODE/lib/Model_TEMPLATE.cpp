// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Template for adding a new model =============  //
// see below for specific references etc ==================  //
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
// If relevant, ensure the below is filled in/replaced ====  //
// Model-specific disclaimer ==============================  //
// Implementation of model XXXX ===========================  //
// Coded from publication YYY =============================  //
// or from distributed source code downloaded from XYZ ====  //
// by YOUR NAME 20XX ======================================  //
// If coded using other code which contains a disclaimer/ =  //
// liscence statement, copy that here also. ===============  //
// ========================================================  //

#include "Model.h"
#include "Structs.h"

// INTSRUCTIONS ============================================  //
// 1) determine an identifier for your model. ==============  //
//  This can have any form you like, with the convention: ==  //
//  speciesCELL_MODEL. For example, Grandi-Bers lab human ==  //
//  atrial myocyte cell model is "hAM_GB" ==================  //
//  This will be used for function and file naming as well =  //
//  as the string used in IF statements to select your =====  //
//  model ==================================================  //
// =========================================================  //
// 2) Copy this file to a new file called lib/Model_X.cpp, =  //
//  where "X" is your model identifier =====================  //
// =========================================================  //
// 3) Search and replace all "speciesCELL_MODEL" in the ====  //
//  new model file with your identifier ====================  //
// =========================================================  //
// 4) In each main file for which your model is valid ======  //
//  (single_cell/tissue; native/integrated) ================  //
//      Single_cell_native_main.cc =========================  //
//      Tissue_native_main.cc ==============================  //
//      Single_cell_3D_main.cc =============================  //
//      Single_cell_0D_main.cc =============================  //
//      Tissue_integrated_main.cc ==========================  //
//  Add a clause exception for your model in the section ===  //
//  "Check Model is appropriate for this code version" =====  //
// =========================================================  //
// 5) Now go through all the functions below and input =====  //
//  your model. Instructions are given for each function ===  //
//  where relevant. ========================================  //
//  Delete unused parameters/functions as neccessary =======  //
//  Functions cover: =======================================  //
//      set_parameters =====================================  //
//      initial_conditions =================================  //
//      compute_model (including all the actual model fns) =  //
//      heterogeneity and modulation =======================  //
// =========================================================  //
// 6) Add the relevant function calls to "lib/Model.c": ====  //
//  Compulsory: ============================================  //
//      set_parameters_native_X() ==========================  //
//      initial_conditions_native_X() ======================  //
//      compute_model_X_native() ===========================  //
//      set_het_mod_X() ====================================  //
//  And also the following if implemented in the model: ====  //
//      update_parameters_integrated_X() ===================  //
//      update_het_and_mod_speciesCELL_MODEL_integrated() ==  //
//      compute_model_X_integrated() =======================  //
// =========================================================  //
// You are now ready to run your model! Compile the code ===  //
// and access your model by typing "Model [identifier]" ====  //
// =========================================================  //

// NOTE on adding new parameters/variables =================  //
// Parameters and variables already exist for the majority =  //
// of what would be required by a new cell model ===========  //
// If exact parameters don't exist, there may be ones which=  //
// are suitable to repurpose (e.g. Ca_j could be used for ==  //
// Ca_ds, CanSR/CajSR as CaSR_1 and 2 etc)==================  //
// However, it may be necessary to add new variables if ====  //
// suitably similar ones do not already exist. =============  //
// Below lists all major types of variables - if adding a ==  //
// whole new current, then most/all of the below will need =  //
// to be added. ============================================  //
// =========================================================  //
// If the variable is only needed by one function, delcare =  //
// it locally. =============================================  //
// =========================================================  //
// If it is a parameter (e.g. conductance, affinity etc) ===  //
// then simply add it to the "Cell_parameters" struct in ===  //
// lib/Structs.h. It is now ready to use - simply set it in=  //
// "set_parameters_native()" ===============================  //
// =========================================================  //
// If it is a computed variable (transition rate, SS..) ====  //
// then add it to "Model_variables" struct. It can now be ==  //
// set and accessed from any relevant function. ============  //
// This includes a current variable itself (e.g. IKACh) ====  //
// State and modifiers will also need to be added if a whole  //
// new current is being added, as well as relevant functions  //
// =========================================================  //
// If it is a new state variable (e.g. gating variable) ====  //
// Then it must be added to the struct "State_variables" ===  //
// Ensure that you set an initial value in =================  //
// "initial_conditions_native_X()", and that it is updated =  //
// within the code. Finally, go to "lib/Read_write_state.c"=  //
// and include an additional IF statement at the bottom of == //
// the functions "{Write/Read}_state_variables_native()" ==== //
// so that it can be saved in state reading/writing. ======== //
// Instructions are presented in those functions ============ //
// ========================================================== //
// If it is a new modification variable / modifier ========== //
// (voltage shift, tau scale, current scale) ================ //
// Then it must be added to the "Cell_parameters" struct. === //
// It MUST ALSO BE DEFAULTED in ============================= //
// "set_modification_defaults_native()" in lib/Model.c ====== //
// Then simply set as needed in set_het_and_mod ============= //
// If the new modifier applies to currents or components ==== //
// present in the models already implemented, then add it to  //
// all appropriate models if you want to be able to apply the //
// new modification to all relevant models. ================= //
// Note: modulation voltage shift variables are implemented = //
// as a shift of the V1/2: V1/2 + vshift which is equivilant= //
// to Vm - vshift (which is how it is implemented) ========== // 
// ========================================================== //

// Function order in this file ========================================================\\|
//  set_parameters_native_speciesCELL_MODEL()
//  update_parameters_native_speciesCELL_MODEL()
//  update_parameters_integrated_speciesCELL_MODEL()
//  initial_conditions_native_speciesCELL_MODEL()
//  compute_model_speciesCELL_MODEL_native()
//  compute_model_speciesCELL_MODEL_integrated()
//  set_gate_rates_speciesCELL_MODEL_native()
//  update_gating_variables_speciesCELL_MODEL_native()
//  compute_Itot_speciesCELL_MODEL_native()
//  compute_Itot_speciesCELL_MODEL_integrated()

//  set_INa_speciesCELL_MODEL_rates() ; update_gates_INa_speciesCELL_MODEL() ; compute_INa_speciesCELL_MODEL()
//  set_INaL_speciesCELL_MODEL_rates() ; update_gates_INaL_speciesCELL_MODEL() ; compute_INaL_speciesCELL_MODEL()
//  set_Ito_speciesCELL_MODEL_rates() ; update_gates_Ito_speciesCELL_MODEL() ; compute_Ito_speciesCELL_MODEL()
//  set_ICaL_speciesCELL_MODEL_rates() ; update_gates_ICaL_speciesCELL_MODEL() ; compute_ICaL_speciesCELL_MODEL()
//  set_IKur_speciesCELL_MODEL_rates() ; update_gates_IKur_speciesCELL_MODEL() ; compute_IKur_speciesCELL_MODEL()
//  set_IKr_speciesCELL_MODEL_rates() ; update_gates_IKr_speciesCELL_MODEL() ; compute_IKr_speciesCELL_MODEL()
//  set_IKs_speciesCELL_MODEL_rates() ; update_gates_IKs_speciesCELL_MODEL() ; compute_IKs_speciesCELL_MODEL()
//  set_IKACh_speciesCELL_MODEL_rates() ; update_gates_IKACh_speciesCELL_MODEL() ; compute_IKACh_speciesCELL_MODEL()
//  set_IK1_speciesCELL_MODEL_variables() ; compute_IK1_speciesCELL_MODEL
//  set_If_speciesCELL_MODEL_rates() ; update_gates_If_speciesCELL_MODEL() ; compute_If_speciesCELL_MODEL()
//  compute_INCX_speciesCELL_MODEL()
//  compute_INaK_speciesCELL_MODEL()
//  compute_ICaP_speciesCELL_MODEL()
//  compute_INab_speciesCELL_MODEL()
//  compute_ICab_speciesCELL_MODEL()
//  compute_IKb_speciesCELL_MODEL()
//  compute_IClCa_speciesCELL_MODEL()
//  compute_IClb_speciesCELL_MODEL()

//  comp_homeostasis_speciesCELL_MODEL()

//  set_het_mod_speciesCELL_MODEL()
//  update_het_and_mod_speciesCELL_MODEL_integrated()
//  set_celltype_native_speciesCELL_MODEL()
//  update_celltype_integrated_speciesCELL_MODEL()
//  set_modulation_ISO_native_speciesCELL_MODEL()
//  set_modulation_Agent_native_speciesCELL_MODEL()
//  set_modulation_Remodelling_native_speciesCELL_MODEL()
//  set_modulation_Mutation_native_speciesCELL_MODEL()
//  set_modulation_ACh_speciesCELL_MODEL()
// End Function order in this file =====================================================//|

// Parameters and specific settings =============================================================\\|
// Set model dependent parameters 
// If your model does not inheret from a cell model already included, define all parameters:
// If your model is similar to an already implemented one, use those same variable names

// !! MUST BE CALLED in lib/Model.c -> set_parameters_native()
void set_parameters_native_speciesCELL_MODEL(Cell_parameters *p)
{
	// Note on units: There is flexibility in the units being used, as long as:
	// 1 - your are self-consistent
	// 2 - the output of the currents is in current density (A/F); Vm is in mV and time-constants/rates are in ms/ms^-1

	// Full parameter list is too large to cleanly put here;
	// some basic parameters below
	// Select which params you need and set to correct values; delete the rest
	// Please see lib/Structs.h -> Cell_parameters struct for full parameter list

	// Default dt if required (if model is stiff then want smaller than global default of 0.02)
	//p->dt               = 0.01; // ms

	// Capacitance and cell structure =========\\|
	//p->Cm               = 100;      // pF   
	//p->Cm_F             = 1.0e-10;  // F    || may need one or both depending on model
	//p->Vcell            = 20100.0;  // um^3
	//p->Vcyto            = 0.68      * p->Vcell; // um^3
	//p->VjSR             = 0.0048    * p->Vcell; // um^3
	//p->VnSR             = 0.0552    * p->Vcell; // um^3

	//p->Fjunc            = 0.11;     // Fraction in junction, where relevant
	//p->Fjunc_ICaL       = 0.9;      // Fraction of ICaL in junction, where relevant
	// End Capacitance and cell structure =====//|

	// Concentrations =========================\\|
	// (constant OR initial condition)
	// Note: CaSR concs are not set as parameters, 
	// only state variables
	//p->Nai              = 11.2;     // mM
	//p->Nao              = 140;      // mM
	//p->Ki               = 139;      // mM
	//p->Ko               = 5.4;      // mM
	//p->Cai              = 0.000102; // mM
	//p->Cao              = 1.8;      // mM

	//p->Nai_sl           = 11.2;     // mM
	//p->Nai_j            = 11.2;     // mM
	//p->Cai_sl           = 0.0001;     // mM
	//p->Cai_j            = 0.0001;     // mM
	// End Concentrations =====================//|

	// Current parameters =====================\\|
	// Magnitudes
	//p->gNa              = 7.8;      // s/mF
	//p->gNaL             = 0.02;     // s/mF
	//p->gto              = 0.1652;   // s/mF
	//p->gKur             = 0.0115;   // s/mF
	//p->gKs              = 0.129;    // s/mF
	//p->gClCa            = 0.1005;   // s/mF
	//p->gKr              = 0.0294*sqrt(p->Ko/5.4);   // s/mF
	//p->gK1              = 0.09*pow(p->Ko/5.4,0.4);  // s/mF
	//p->gf				  = x;`		// s/mF

	//p->gCaL             = 0.1238;   // s/mF  || current density
	//p->pCaL             = 2.7e-4;   // cm/s  || channel flux
	//p->pCaL_K           = 1.35e-7;  // cm/s
	//p->pCaL_Na          = 0.75e-8;  // cm/s
	//p->ICaL_vi_Fs       = 0.4;      // fraction slow voltage inactivation gate

	//p->INCX_bar         = 1600;     // pA/pF
	//p->INaK_bar         = 0.59933874;   // pA/pF
	//p->ICaP_bar         = 0.275;        // pA/p
	//p->gNab             = 0.000674; // s/mF
	//p->gCab             = 0.001131; // s/mF

	//p->gKACh_max	      = x; // s/mF  (total conductance, or maximum to multiply by conc dependency)

	// Affinity constants etc
	//p->INCX_kNao        = 87.5;     // mM
	//p->INCX_kNao        = 87.5;     // mM
	//p->INCX_kCao        = 1.3;      // mM
	//p->INCX_kNai        = 12.3;     // mM
	//p->INCX_kCai        = 3.59e-3;  // mM
	//p->INCX_k           = 0.1;
	//p->INCX_gamma       = 0.35;
	//p->INCX_kda         = 0.384e-3; // mM
	//p->INCX_ksat        = 0.27;

	//p->INaK_kK          = 1.5;          // mM
	//p->INaK_kNa         = 10.0;         // mM

	// p->ICaP_kCa         = 0.0005;       // mM
	// End Current parameters =================//|

	// Stimulus parameters ====================\\|
	// (if different to default)
	//p->stimduration 	= 5.0;      // ms
	//p->stimmag      	= -13.5;    // pA/pF
	// End Stimulus parameters ================//|

	// Celltype and modulation defaults =======\\|
	// Default celltype
	//p->Celltype			= "some_region";

	// Default ISO/ACh model (if multiple exist; if not, will be set to "default" already and no need to add here)
	//p->ISO_model	    = "Version_2";
	//p->ACh_model		= "Version_B";
	// End Celltype and modulation defaults ===//|

	// Ca-handling parameters =================\\|
	//p->cmdnbar          = 0.050;    // mM
	//p->trpnbar          = 0.070;    // mM
	//p->csqnbar          = 10.0;     // mM
	//p->cmdn_k           = 0.00238;  // mM
	//p->trpn_k           = 0.0005;   // mM
	//p->csqn_k           = 0.8;      // mM

	//p->J_rel_max        = 30.0;     // mM/ms
	//p->J_SERCA_max      = 0.005;    // mM/ms
	//p->J_SERCA_kCa      = 0.00092;  // mM
	//p->J_leak_max       = 0.005;    // mM/ms
	//p->J_leak_kCaSR     = 15.0;     // mM
	//p->J_jsr_nsr_tau    = 180.0;    // ms
	// End Ca-handling parameters =============//|
}

// ============= OR ==============\\|
// If your model does inheret from a cell model already included, define only different parameters here
// Anything here will overwrite those set by original function
// Remember to call the first function in Model.c !!

// !! MUST BE CALLED in lib/Model.c -> set_parameters_native()
void update_parameters_native_speciesCELL_MODEL(Cell_parameters *p)
{
	// Updated parameters
	//p->gto                = X;        // s/mF
}

// AND for integrated, may want to change the Ca parameters (which have same variable names for native or integrated, 
// and are overwritten from native with integrated defaults)
// This is only if your model has different settings to default integrated
// Defaults are set in lib/Initialisation.c -> set_parameters_spatial_Ca_defaults()

// !! MUST BE CALLED in lib/Model.c -> set_parameters_spatial_Ca()
// OPTIONAL
void update_parameters_integrated_speciesCELL_MODEL(Cell_parameters *p)
{
	// Look in ib/Initialisation.c -> set_parameters_spatial_Ca_defaults()
	// and place here any parameter which is different
}

// Initial conditions
// !! MUST BE CALLED in lib/Model.c -> initial_conditions_native()
void initial_conditions_native_speciesCELL_MODEL(State_variables *s, Cell_parameters p)
{
	// Define ICs for ALL state variables used in the model
	// These are just baseline (steady-state) ICs: specific, conditional
	// paced ICs controlled through read/write state.
	// All state variables listed below - delete and change as appropriate

	// Vm
	s->Vm                   = -85; // mV

	// Ion currents
	//s->INa_va               = 0;              // voltage activation
	//s->INa_vi_1             = 0;            // voltage inactivation
	//s->INa_vi_2             = 0;            // voltage inactivation 2
	//s->INaL_va              = 0;             // voltage activation
	//s->INaL_vi              = 0;             // voltage inactivation
	//s->Ito_va               = 0;              // voltage activation
	//s->Ito_vi               = 0;              // voltage inactivation (fast if slow also used)
	//s->Ito_vi_s             = 0;            // voltage inactivation, slow
	//s->Ito_vi_3             = 0;            // third Ito voltage inactivation gate
	//s->ICaL_va              = 0;             // voltage activation
	//s->ICaL_vi              = 0;             // voltage inactivation (fast if slow also used)
	//s->ICaL_vi_s            = 0;           // voltage inactivation, slow
	//s->ICaL_ci              = 0;             // calcium inactivation
	//s->ICaL_ci_j            = 0;           // calcium inactivation, junctional
	//s->IKur_va              = 0;             // voltage activation
	//s->IKur_vi              = 0;             // voltage inactivation (fast if slow also used)
	//s->IKr_va               = 0;              // voltage activation
	//s->IKr_vi               = 0;              // used where IKr is used for IKf in Rabbit atria || also used for va_slow
	//s->IKs_va               = 0;              // voltage activation
	//s->IKs_va_2             = 0;            // voltage activation 2
	//s->IK1_va               = 0;
	//s->IKACh_va             = 0;            // voltage activation
	//s->IKACh_vi             = 0;            // voltage inactivation
	//s->If_va                = 0;               // voltage activation

	// Ca2+ handling
	//s->cmdn = 0;                // calmodulin concentration     (mM)
	//s->trpn = 0;                // troponin concentration       (mM)
	//s->csqn = 0;                // calsequestrin concenrration  (mM)
	//s->CajSR = 0;               // junctional SR / release compartment Ca concentration (mM)
	//s->CanSR = 0;               // network SR / uptake compartment Ca concentration     (mM)
	//s->RyRo = 0;                // proportion open RyRs / activation gate
	//s->RyRr = 0;                // proportion refactory RyRs / inactivation gate 1
	//s->RyRi = 0;                // proportion inactivated RyRs / inactvation gate 2
	//s->Tn_CHm = 0;
	//s->Tn_CHc = 0;
	//s->Myo_m = 0;
	//s->Myo_c = 0;
	//s->Tn_CL = 0;
	//s->CaCalse = 0;
	//s->CaCal = 0;
	//s->Catrop = 0;
	//s->Camg = 0;
	//s->Mgmg = 0;

	// Concentrations
	// Assign state from param (even if constant)
	// The state variable is what will actually be used in calculations
	// If constant, it simply won't be updated
	// Ensure the ICs for all dynamic concentrations are set as a parameter
	//s->Nai              = p.Nai;
	//s->Nai_j            = p.Nai_j;
	//s->Nai_sl           = p.Nai_sl;
	//s->Nao              = p.Nao;
	//s->Ki               = p.Ki;
	//s->Ko               = p.Ko;
	//s->Cai              = p.Cai;
	//s->Cai_j            = p.Cai_j;
	//s->Cai_sl           = p.Cai_sl;
	//s->Cao              = p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Compute model functions ======================================================================\\|
// Your model may have more or fewer currents than this template - just follow the procedure and add/delete as appropriate

// !! MUST BE CALLED in lib/Model.c -> compute_model_native()
void compute_model_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);     // lib/Model.c || replace with model-specific function if different/more complex
	set_gate_rates_speciesCELL_MODEL_native(p, var, Vm, s->Cai);
	update_gating_variables_speciesCELL_MODEL_native(p, var, s, Vm, dt);
	compute_Itot_speciesCELL_MODEL_native(p, var, s, Vm);
	comp_homeostasis_speciesCELL_MODEL(p, var, s, Vm, dt);
}

// !! MUST BE CALLED in lib/Model.c -> compute_model_integrated()
// OPTIONAL (not needed if not integrating model with Ca2+ handling system!)
void compute_model_speciesCELL_MODEL_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);     // lib/Model.c || replace with model-specific function if different/more complex
	set_gate_rates_speciesCELL_MODEL_native(p, var, Vm, s->Cai);
	update_gating_variables_speciesCELL_MODEL_native(p, var, s, Vm, dt);
	compute_Itot_speciesCELL_MODEL_integrated(p, var, s, Vm);
	// NO homeostasis here, as done in Ca handling model
	// Can add a function which does K+ and Na+ cycling if required
	// Or just define those explicitly here
}

// !! MUST BE CALLED in above compute_model_species() function
void set_gate_rates_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	// Call only currents you need
	// you can call existing functions from other models here also - doesn't have to be new model specific
	// just ensure the appropriate parameters are set and the variables are consistent between set, update and compute functions
	// if mixing models.

	// INa - common to use LR
	//set_INa_LR_rates(p, var, Vm);                   // lib/Model.c

	// OR model-specific
	//set_INa_speciesCELL_MODEL_rates(p, var, Vm);				

	//set_INaL_speciesCELL_MODEL_rates(p, var, Vm);				
	//set_Ito_speciesCELL_MODEL_rates(p, var, Vm);
	//set_IKur_speciesCELL_MODEL_rates(p, var, Vm);
	//set_IKs_speciesCELL_MODEL_rates(p, var, Vm);
	//set_IKr_speciesCELL_MODEL_rates(p, var, Vm);
	//set_IK1_speciesCELL_MODEL_variables(p, var, Vm);
	//set_ICaL_speciesCELL_MODEL_rates(p, var, Vm, Cai);
	//set_If_speciesCELL_MODEL_rates(p, var, Vm, Cai);
}

// !! MUST BE CALLED in above compute_model_species() function
void update_gating_variables_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// INa - common to use LR
	//update_gates_INa_LR(p, var, s, Vm, dt);         // lib/Model.c

	// OR model-specific
	//update_gates_INa_speciesCELL_MODEL(p, var, s, Vm, dt); 	

	//update_gates_INaL_speciesCELL_MODEL(p, var, s, Vm, dt); 	
	//update_gates_IKs_speciesCELL_MODEL(p, var, s, Vm, dt);
	//update_gates_IKr_speciesCELL_MODEL(p, var, s, Vm, dt);
	//update_gates_ICaL_speciesCELL_MODEL(p, var, s, Vm, dt);
	//update_gates_Ito_speciesCELL_MODEL(p, var, s, Vm, dt);
	//update_gates_IKur_speciesCELL_MODEL(p, var, s, Vm, dt);
	//update_gates_If_speciesCELL_MODEL(p, var, s, Vm, dt);
}

// !! MUST BE CALLED in above compute_model_species() function
void compute_Itot_speciesCELL_MODEL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	// INa - common to use LR
	// compute_INa_LR(p, var, s, Vm);                  // lib/Model.c

	// OR model-specific || Note: rates can be LR and compute is model-specific e.g. SL/cyto comparts
	//compute_INa_speciesCELL_MODEL(p, var, s, Vm);					

	//compute_INaL_speciesCELL_MODEL(p, var, s, Vm);					
	//compute_Ito_speciesCELL_MODEL(p, var, s, Vm);
	//compute_ICaL_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKur_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKr_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKs_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IK1_speciesCELL_MODEL(p, var, s, Vm);
	//compute_If_speciesCELL_MODEL(p, var, s, Vm);
	//compute_INCX_speciesCELL_MODEL(p, var, s, Vm);
	//compute_INaK_speciesCELL_MODEL(p, var, s, Vm);
	//compute_ICaP_speciesCELL_MODEL(p, var, s, Vm);
	//compute_INab_speciesCELL_MODEL(p, var, s, Vm);
	//compute_ICab_speciesCELL_MODEL(p, var, s, Vm);

	// Add all the computed currents here
    // DO NOT include Istim, as this is added in main
	//var->Itot   = var->INa + var->Ito + var->IK1 + var->ICaL + var->IKur + var->INCX + var->INaK + var->ICaP + var->INab + var->ICab + var->IKr + var->IKs;
}

// !! MUST BE CALLED in above compute_model_species() function
// OPTIONAL (integrated models only)
void compute_Itot_speciesCELL_MODEL_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	// INa - common to use LR
	// compute_INa_LR(p, var, s, Vm);                  // lib/Model.c

	// OR model-specific || Note: rates can be LR and compute is model-specific e.g. SL/cyto comparts
	//compute_INa_speciesCELL_MODEL(p, var, s, Vm);					

	//compute_INaL_speciesCELL_MODEL(p, var, s, Vm);					
	//compute_INa_speciesCELL_MODEL(p, var, s, Vm);                 
	//compute_Ito_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKur_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKr_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IKs_speciesCELL_MODEL(p, var, s, Vm);
	//compute_IK1_speciesCELL_MODEL(p, var, s, Vm);
	//compute_INaK_speciesCELL_MODEL(p, var, s, Vm);
	//compute_INab_speciesCELL_MODEL(p, var, s, Vm);

	// Add all the non-Ca computed currents here
	//var->Itot   = var->INa + var->Ito + var->IK1 + var->ICaL + var->IKur + var->INCX + var->INaK + var->ICaP + var->INab + var->ICab + var->IKr + var->IKs;
	// Then add the Ca currents, computed in CRU
	var->Itot   += var->INCX + var->ICaP + var->ICab + var->ICaL;   // Adds Ca currents computed in CRU
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|

// NOTES: ===================================================================\\|
// All available and relevant state variables, computed variables, current parameters
// and modulation parameters are listed in each current's function; use only as needed.
// The modulation variables (shitfs/scales) used are what need to be modified in het and modulation.
// You can use none, any or all of them.
// And they can be used however you like.
// However, use as suggested to ensure consistency between different models.
// End NOTES: ===============================================================//|

// INa ======================================================================\\|
// Identical to that of LR model, found in lib/Model.c 
// or 
// new formulation here
void set_INa_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	// Implementation for one activation and two inactivation gating variables
	// Implementation for global INa and compartments INa_sl and INa_j (single voltage kinetics)
	// Can be alpha/beta or steady-state/tau for any/all gates

	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables: (not used in this function, but the gating
	// variables for which rates/time-constants/steady-states must be set)
	//s->INa_va;
	//s->INa_vi_1;
	//s->INa_vi_2;

	// Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
	//p.INa_va_shift;         // volttage dep activation -> V1/2 shift of SS and V of tau/rates, mV
	//p.INa_vi_shift;         // volttage dep inactivation -> V1/2 shift of SS and V of tau/rates, mV
	//p.INa_va_tau_scale;    // Scales time constant for INa voltage activation
	//p.INa_vi_1_tau_scale;  // Scales time constant for INa voltage inactivation 1
	//p.INa_vi_2_tau_scale;  // Scales time constant for INa voltage inactivation 2

	// Relevant computed variables:
	//var->INa_va_ss;           // voltage activation, steady state
	//var->INa_va_tau;          // voltage activation, time constant
	//var->INa_vi_1_ss;         // voltage inactivation, steady state
	//var->INa_vi_1_tau;        // voltage inactivation, time constant
	//var->INa_vi_2_ss;         // voltage inactivation, steady state
	//var->INa_vi_2_tau;        // voltage inactivation, time constant
	//var->INa_va_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->INa_va_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->INa_vi_1_al;         // voltage activation, alpha transition rate (1-y -> y)
	//var->INa_vi_1_bet;        // voltage activation, beta transition rate  (y -> 1-y)
	//var->INa_vi_2_al;         // voltage activation, alpha transition rate (1-y -> y)
	//var->INa_vi_2_bet;        // voltage activation, beta transition rate  (y -> 1-y)
	// End Available Variables and parameters for reference =======//|

	// Set gate rates =============================================\\|
	// Define voltage to be used in calulcations || Note here, ss and tau/rate not distinguished
	// -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
	// Note: modulation voltage shift variables are implemented = //
	// as a shift of the V1/2: V1/2 + vshift which is equivilant= //
	// to Vm - vshift (which is how it is implemented) ========== // 
	double Vm_ac        = Vm - p.INa_va_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac      = Vm - p.INa_vi_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

	// Set Activation gate alpha and beta
	//var->INa_va_al                 = f(Vm_ac)
	//var->INa_va_bet                = f(Vm_ac)

	// Set inactivation gates alphas and betas
	//var->INa_vi_1_al           = f(Vm_inac)
	//var->INa_vi_1_bet          = f(Vm_inac)
	//var->INa_vi_2_al           = f(Vm_inac)
	//var->INa_vi_2_bet          = f(Vm_inac)

	// Set tau and SS from alpha and beta
	//var->INa_va_tau                = 1.0/(var->INa_va_al + var->INa_va_bet); // 1/(a+b)
	//var->INa_vi_1_tau              = 1.0/(var->INa_vi_1_al + var->INa_vi_1_bet);
	//var->INa_vi_2_tau              = 1.0/(var->INa_vi_2_al + var->INa_vi_2_bet);
	//var->INa_va_ss                 = var->INa_va_al * var->INa_va_tau; // a*tau
	//var->INa_vi_1_ss               = var->INa_vi_1_al * var->INa_vi_1_tau;
	//var->INa_vi_2_ss               = var->INa_vi_2_al * var->INa_vi_2_tau;

	// OR define steady-state and tau directly
	//var->INa_va_ss      = f(Vm_ac) e.g:
	//  var->INa_va_ss      = sigmoid(Vm_ac, V1/2, -k);  OR:
	//  var->INa_va_ss      = 1/(1 + exp((Vm_ac - V1/2)/-k) ); or more complex
	//var->INa_va_tau    = f(Vm_ac) or constant

	//var->INa_vi_1_ss      = f(Vm_inac) e.g:
	//  var->INa_vi_1_ss      = sigmoid(Vm_inac, V1/2, k);  OR:
	//  var->INa_vi_1_ss      = 1/(1 + exp((Vm_inac - V1/2)/k) ); or more complex
	//var->INa_vi_1_tau    = f(Vm_inac) or constant
	// and again for vi_2 = f(Vm_inac) 

	// Don't forget to scale time constants, or tau_scale will do nothing!
	var->INa_va_tau                 *= p.INa_va_tau_scale;
	var->INa_vi_1_tau               *= p.INa_vi_1_tau_scale;
	var->INa_vi_2_tau               *= p.INa_vi_2_tau_scale;
	// End Set gate rates =========================================//|
}

void update_gates_INa_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables - update ALL which are used!
	//s->INa_va;
	//s->INa_vi_1;            
	//s->INa_vi_2;         

	// Relevant computed variables:
	//var->INa_va_ss;           // voltage activation, steady state
	//var->INa_va_tau;          // voltage activation, time constant
	//var->INa_vi_1_ss;         // voltage inactivation, steady state
	//var->INa_vi_1_tau;        // voltage inactivation, time constant
	//var->INa_vi_2_ss;         // voltage inactivation, steady state
	//var->INa_vi_2_tau;        // voltage inactivation, time constant
	//var->INa_va_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->INa_va_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->INa_vi_1_al;         // voltage activation, alpha transition rate (1-y -> y)
	//var->INa_vi_1_bet;        // voltage activation, beta transition rate  (y -> 1-y)
	//var->INa_vi_2_al;         // voltage activation, alpha transition rate (1-y -> y)
	//var->INa_vi_2_bet;        // voltage activation, beta transition rate  (y -> 1-y)
	// End Available Variables and parameters for reference =======//|

	// Update by rush-larsen                      variable     steady-state  time-constant   dt
	//s->INa_va                      = rush_larsen(s->INa_va, var->INa_va_ss, var->INa_va_tau, dt); // lib/Membrane.c
	//s->INa_vi_1                    = rush_larsen(s->INa_vi_1, var->INa_vi_1_ss, var->INa_vi_1_tau, dt);
	//s->INa_vi_2                    = rush_larsen(s->INa_vi_2, var->INa_vi_2_ss, var->INa_vi_2_tau, dt);

	// Or explicitly, FE or other
	//s->INa_va                       += dt*(differential = f(ss, tau, alpha, beta..))
}

void compute_INa_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	// Relevant state variables, reversal potential, conductance and scale factors, 
	// junction/fast-slor
	//s->INa_va;
	//s->INa_vi_1;
	//s->INa_vi_2;
	//var->ENa, var->ENa_sl, var->ENa_j;
	//p.gNa;  - conductance
	//p.GNa;  - scale factor
	//p.Fjunc;

	// e.g, simple:
	//var->INa        = p.gNa * pow(s->INa_va, 3) * s->INa_vi_1 * s->INa_vi_2 * (Vm - var->ENa);
	var->INa        *= p.GNa; // SCALE!!  Must be here for GNa to have any effect

	// e.g, cyto/sl:
	//var->INa_sl         = (1 - p.Fjunc) * p.gNa * pow(s->INa_va, 3) * s->INa_vi_1 * s->INa_vi_2 * (Vm - var->ENa_sl);
	//var->INa_j          = (    p.Fjunc) * p.gNa * pow(s->INa_va, 3) * s->INa_vi_1 * s->INa_vi_2 * (Vm - var->ENa_j);
	var->INa_sl         *= p.GNa; // SCALE!!  Must be here for GNa to have any effect
	var->INa_j          *= p.GNa;
	var->INa            = var->INa_sl + var->INa_j;
}
// End INa ==================================================================//|

// INaL =====================================================================\\|
void set_INaL_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	// Implementation for one activtion and one inactivation gating variable
	// Implementation for global INaL and compartments INaL_sl and INaL_j (single voltage kinetics)
	// Can be alpha/beta or steady-state/tau for either/both gates

	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables: (not used in this function, but the gating
	// variables for which rates/time-constants/steady-states must be set)
	// s->INaL_va;
	// s->INaL_vi;

	// Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
	// Alpha beta:
	//  p.INaL_va_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	//  p.INaL_vi_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation
	//
	// OR steady-state, time constant
	//  p.INaL_va_ss_shift; // shift of the activation steady state
	//  p.INaL_vi_ss_shift; // shift of the inactivation steady state
	//  p.INaL_va_tau_shift; // Shift of voltage dependence of activation tau
	//  p.INaL_vi_tau_shift; // Shift of voltage dependence of inactivation tau
	//
	// tau scale
	//  p.INaL_va_tau_scale;
	//  p.INaL_vi_tau_scale;
	//
	// steady-state gradient scale
	//  p.INaL_va_ss_kscale
	//  p.INaL_vi_ss_kscale

	// Relevant computed variables:   
	//var->INaL_va_ss;          // voltage activation, steady state
	//var->INaL_va_tau;         // voltage activation, time constant
	//var->INaL_vi_ss;          // voltage inactivation, steady state
	//var->INaL_vi_tau;         // voltage inactivation, time constant
	//var->INaL_va_al;          // voltage activation, alpha transition rate (1-y-> y)
	//var->INaL_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
	//var->INaL_vi_al;          // voltage inactivation, alpha transition rate (1-y-> y)
	//var->INaL_vi_bet;         // voltage inactivation, beta transition rate  (y -> 1-y)
	// End Available Variables and parameters for reference =======//|

	// Set gate rates =============================================\\|
	// alpha-beta style
	// -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
	// Note: modulation voltage shift variables are implemented = //
	// as a shift of the V1/2: V1/2 + vshift which is equivilant= //
	// to Vm - vshift (which is how it is implemented) ========== // 
	double Vm_ac        = Vm - p.INaL_va_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac      = Vm - p.INaL_vi_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

	// Set Activation gate alpha and beta
	//var->INaL_va_al                 = f(Vm_ac)
	//var->INaL_va_bet                = f(Vm_ac)

	// Set inactivation gates alphas and betas
	//var->INaL_vi_al           = f(Vm_inac)
	//var->INaL_vi_bet          = f(Vm_inac)

	// Set tau and SS from alpha and beta
	//var->INaL_va_tau                = 1.0/(var->INaL_va_al + var->INaL_va_bet); // 1/(a+b)
	//var->INaL_vi_tau              = 1.0/(var->INaL_vi_al + var->INaL_vi_bet);
	//var->INaL_va_ss                 = var->INaL_va_al * var->INaL_va_tau; // a*tau
	//var->INaL_vi_ss               = var->INaL_vi_al * var->INaL_vi_tau;

	// steady-state, tau style
	// -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
	double Vm_ac_ss     = Vm - p.INaL_va_ss_shift;  // Shift of the voltage used to calculate activation steady state
	double Vm_inac_ss   = Vm - p.INaL_vi_ss_shift;  // Shift of the voltage used to calculate inactivation steady state
	double Vm_ac_tau    = Vm - p.INaL_va_tau_shift;  // Shift of the voltage used to calculate activation time constant
	double Vm_inac_tau  = Vm - p.INaL_vi_tau_shift;  // Shift of the voltage used to calculate inactivation time constant

	// simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

	//var->INaL_va_ss      = f(Vm_ac_ss, p.INaL_va_ss_kscale) e.g:
	//  var->INaL_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.INaL_va_ss_kscale);  OR:
	//  var->INaL_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.INaL_va_ss_kscale)) ); or more complex
	//var->INaL_va_tau    = f(Vm_ac_tau) or constant

	//var->INaL_vi_ss      = f(Vm_inac_ss, p.INaL_vi_ss_kscale) e.g:
	//  var->INaL_vi_ss      = sigmoid(Vm_inac_ss, V1/2, k*p.INaL_vi_ss_kscale);  OR:
	//  var->INaL_vi_ss      = 1/(1 + exp((Vm_inac_ss - V1/2)/(k*p.INaL_vi_ss_kscale)) ); or more complex
	//var->INaL_vi_tau    = f(Vm_inac_tau) or constant

	// Scale time constants for modulation
	var->INaL_va_tau    *= p.INaL_va_tau_scale;
	var->INaL_vi_tau    *= p.INaL_vi_tau_scale;
	// End Set gate rates =========================================//|
}

void update_gates_INaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables - update ALL which are used!
	//s->INaL_va;
	//s->INaL_vi;

	// Relevant computed variables:
	//var->INaL_va_ss;          // voltage activation, steady state
	//var->INaL_va_tau;         // voltage activation, time constant
	//var->INaL_vi_ss;          // voltage inactivation, steady state
	//var->INaL_vi_tau;         // voltage inactivation, time constant
	//var->INaL_va_al;          // voltage activation, alpha transition rate (1-y-> y)
	//var->INaL_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
	//var->INaL_vi_al;          // voltage inactivation, alpha transition rate (1-y-> y)
	//var->INaL_vi_bet;         // voltage inactivation, beta transition rate  (y -> 1-y)
	// End Available Variables and parameters for reference =======//|

	// Update by rush-larsen                variable     steady-state  time-constant   dt
	//s->INaL_va              = rush_larsen(s->INaL_va, var->INaL_va_ss, var->INaL_va_tau, dt);
	//s->INaL_vi              = rush_larsen(s->INaL_vi, var->INaL_vi_ss, var->INaL_vi_tau, dt);

	// Or explicitly, FE or other
	//s->INaL_va                       += dt*(differential = f(ss, tau, alpha, beta..))
}

void compute_INaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	// Relevant state variables, reversal potential, conductance and scale factor +  junction factor
	//s->INaL_va;
	//s->INaL_vi;
	//var->ENa, var->ENa_sl, var->ENa_j;
	//p.gNaL;  - conductance
	//p.GNaL;  - scale factor
	//p.Fjunc;

	// global current
	//var->INaL               = p.gNaL * s->INaL_va * s->INaL_vi * (Vm - var->ENa);
	var->INaL               *= p.GNaL;

	// compartments
	//var->INaL_sl        = (1 - p.Fjunc)  * p.gNaL * pow(s->INaL_va, 3) * s->INaL_vi * (Vm - var->ENa_sl);
	//var->INaL_j         = (    p.Fjunc)  * p.gNaL * pow(s->INaL_va, 3) * s->INaL_vi * (Vm - var->ENa_j);
	var->INaL_sl        *= p.GNaL;
	var->INaL_j         *= p.GNaL;
	var->INaL           = var->INaL_sl + var->INaL_j;
}
// End INaL =================================================================//|

// Ito ======================================================================\\|
void set_Ito_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	// Implementation for one activtion and three inactivation gating variables
	// Can be alpha/beta or steady-state/tau for ac and inac gates
	// Has functionality for fast and slow inactivation gates set by proportion

	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables: (not used in this function, but the gating
	// variables for which rates/time-constants/steady-states must be set)
	// s->Ito_va;              // voltage activation
	// s->Ito_vi;              // voltage inactivation (fast if slow also used)
	// s->Ito_vi_s;            // voltage inactivation, slow
	// s->Ito_vi_3;            // third Ito voltage inactivation gate

	// Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
	// Voltage shifts
	//  p.Ito_va_ss_shift; // shift of the activation steady state
	//  p.Ito_vi_ss_shift; // shift of the inactivation steady state
	//  p.Ito_va_tau_shift; // Shift of voltage dependence of activation tau
	//  p.Ito_vi_tau_shift; // Shift of voltage dependence of inactivation tau
	//
	// tau scale
	//  p.Ito_va_tau_scale;
	//  p.Ito_vi_tau_scale;
	//
	// steady-state gradient scale
	//  p.Ito_va_ss_kscale
	//  p.Ito_vi_ss_kscale

	// Relevant computed variables:
	//var->Ito_va_ss;           // voltage activation, steady state
	//var->Ito_va_tau;          // voltage activation, time constant
	//var->Ito_vi_ss;           // voltage inactivation, steady state
	//var->Ito_vi_3_ss;         // voltage inactivation, steady state, alt
	//var->Ito_vi_tau;          // voltage inactivation, time constant
	//var->Ito_vi_s_tau;        // voltage inactivation, time constant, slow
	//var->Ito_vi_3_tau;        // voltage inactivation, time constant, alt
	//var->Ito_va_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->Ito_va_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->Ito_vi_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->Ito_vi_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->Ito_vi_Fs;           // Fraction of slow voltage inactivation
	// End Available Variables and parameters for reference =======//|

	// Set gate rates =============================================\\|
	// steady-state, tau style
	// -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
	double Vm_ac_ss     = Vm - p.Ito_va_ss_shift;  // Shift of the voltage used to calculate activation steady state
	double Vm_inac_ss   = Vm - p.Ito_vi_ss_shift;  // Shift of the voltage used to calculate inactivation steady state
	double Vm_ac_tau    = Vm - p.Ito_va_tau_shift;  // Shift of the voltage used to calculate activation time constant
	double Vm_inac_tau  = Vm - p.Ito_vi_tau_shift;  // Shift of the voltage used to calculate inactivation time constant

	// simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

	//var->Ito_va_ss      = f(Vm_ac_ss, p.Ito_va_ss_kscale) e.g:
	//  var->Ito_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.Ito_va_ss_kscale);  OR:
	//  var->Ito_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.Ito_va_ss_kscale)) ); or more complex
	//var->Ito_va_tau    = f(Vm_ac_tau) or constant

	//var->Ito_vi_ss      = f(Vm_inac_ss, p.Ito_vi_ss_kscale) e.g:
	//  var->Ito_vi_ss      = sigmoid(Vm_inac_ss, V1/2, k*p.Ito_vi_ss_kscale);  OR:
	//  var->Ito_vi_ss      = 1/(1 + exp((Vm_inac_ss - V1/2)/(k*p.Ito_vi_ss_kscale)) ); or more complex
	//var->Ito_vi_tau    = f(Vm_inac_tau) or constant
	// Ito_vi_3_ss, Ito_vi_s_tau, Ito_vi_3_tau similar

	// alpha-beta style
	// -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
	double Vm_ac        = Vm - p.Ito_va_ss_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac      = Vm - p.Ito_vi_ss_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

	// Set Activation gate alpha and beta
	//var->Ito_va_al                 = f(Vm_ac)
	//var->Ito_va_bet                = f(Vm_ac)

	// Set inactivation gates alphas and betas
	//var->Ito_vi_al           = f(Vm_inac)
	//var->Ito_vi_bet          = f(Vm_inac)

	// Set tau and SS from alpha and beta
	//var->Ito_va_tau                = 1.0/(var->Ito_va_al + var->Ito_va_bet); // 1/(a+b)
	//var->Ito_vi_tau              = 1.0/(var->Ito_vi_al + var->Ito_vi_bet);
	//var->Ito_va_ss                 = var->Ito_va_al * var->Ito_va_tau; // a*tau
	//var->Ito_vi_ss               = var->Ito_vi_al * var->Ito_vi_tau;

	// Set fraction of fast and slow if relevant
	// var->Ito_vi_Fs = f(Vm) or constant or function of state or other variable

	// and don't forget to scale the time-constants!
	var->Ito_va_tau			*= p.Ito_va_tau_scale;
	var->Ito_vi_tau			*= p.Ito_vi_tau_scale;
	var->Ito_vi_s_tau		*= p.Ito_vi_tau_scale;
	var->Ito_vi_3_tau		*= p.Ito_vi_tau_scale;
	// End Set gate rates =========================================//|
}

void update_gates_Ito_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Available Variables and parameters for reference ===========\\|
	// Relevant state variables - update all which are used
	//s->Ito_va;              // voltage activation
	//s->Ito_vi;              // voltage inactivation (fast if slow also used)
	//s->Ito_vi_s;            // voltage inactivation, slow
	//s->Ito_vi_3;            // third Ito voltage inactivation gate

	// Relevant computed variables:
	//var->Ito_va_ss;           // voltage activation, steady state
	//var->Ito_va_tau;          // voltage activation, time constant
	//var->Ito_vi_ss;           // voltage inactivation, steady state
	//var->Ito_vi_3_ss;         // voltage inactivation, steady state, alt
	//var->Ito_vi_tau;          // voltage inactivation, time constant
	//var->Ito_vi_s_tau;        // voltage inactivation, time constant, slow
	//var->Ito_vi_3_tau;        // voltage inactivation, time constant, alt
	//var->Ito_va_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->Ito_va_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->Ito_vi_al;           // voltage activation, alpha transition rate (1-y-> y)
	//var->Ito_vi_bet;          // voltage activation, beta transition rate  (y -> 1-y)
	//var->Ito_vi_Fs;           // Fraction of slow voltage inactivation
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
	//s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);

    // Or explicitly, FE or other
	//s->Ito_vi              += dt*(differential= f(ss, tau, alpha, beta..))
}

void compute_Ito_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potential, conductance and scale factor 
    //s->Ito_va;              // voltage activation
    //s->Ito_vi;              // voltage inactivation (fast if slow also used)
    //s->Ito_vi_s;            // voltage inactivation, slow
    //s->Ito_vi_3;            // third Ito voltage inactivation gate
    //var->EK
    //p.gto;  - conductance
    //p.Gto;  - scale factor
    //var->Ito_vi_Fs - fraction of slow inactivation channels
	
    // Simple
	//var->Ito				= p.gto * pow(s->Ito_va, 3) * s->Ito_vi * (Vm - var->EK);

    // fast-slow fraction
	//var->Ito				= p.gto * s->Ito_va * ( var->Ito_vi_Fs * s->Ito_vi_s + (1 - var->Ito_vi_Fs) * Ito_vi ) * (Vm - var->EK);
    
    // etc for any combination of gates

	// Don't forget to multiply by scale factor!!
	var->Ito				*= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
// ICaL can be written exactly as above for Ito. However, if there are intentions
// to integrate the new model with the "integrated" calcium handling system
// e.g. for spatial cell models or spontaneous release functions, please 
// follow the procedure below: voltage-dependent gates have their own
// functions so can be called elsewhere (i.e. from CRU)
void set_ICaL_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	// Implementation for one voltage activtion, two voltage inactivation 
    // and one calcium inactivation gating variable.
    // Voltage gates steady-state/tau only
    // Calcium inactivation gate alpha/beta or stady-state tau
    // Implementation for global ICaL and compartments ICaL_sl and ICaL_j (single voltage kinetics)

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    //s->ICaL_va;             // voltage activation
    //s->ICaL_vi;             // voltage inactivation (fast if slow also used)
    //s->ICaL_vi_s;           // voltage inactivation, slow
    //s->ICaL_ci;             // calcium inactivation
    //s->ICaL_ci_j;           // calcium inactivation, junctional

    // Calcium concentrations 
    // s->Cai, s->Cai_sl, s->Cai_j

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //  p.ICaL_va_ss_shift; // shift of the activation steady state
    //  p.ICaL_vi_ss_shift; // shift of the inactivation steady state
    //  p.ICaL_va_tau_shift; // Shift of voltage dependence of activation tau
    //  p.ICaL_vi_tau_shift; // Shift of voltage dependence of inactivation tau
    //
    // tau scale
    //  p.ICaL_va_tau_scale;
    //  p.ICaL_vi_tau_scale;
    //
    // steady-state gradient scale
    //  p.ICaL_va_ss_kscale
    //  p.ICaL_vi_ss_kscale

    // Relevant computed variables:
    //var->ICaL_va_ss;          // voltage activation, steady state
    //var->ICaL_va_tau;         // voltage activation, time constant
    //var->ICaL_vi_ss;          // voltage inactivation, steady state
    //var->ICaL_vi_tau;         // voltage inactivation, time constant
    //var->ICaL_vi_s_tau;       // voltage inactivation, time constant, slow
    //var->ICaL_ci_ss;          // calcium inactivation, steady state
    //var->ICaL_ci_j_ss;        // calcium inactivation, steady state, junctional
    //var->ICaL_ci_tau;         // calcium inactivation, time constant
    //var->ICaL_ci_al;          // Calcium induced activation, alpha rate       (ms^-1)
    //var->ICaL_ci_bet;         // Calcium induced activation, beta rate        (ms^-1)

    // Relevant parameters (must be set in set_parameters function)
    // p.ICaL_vi_Fs               // fraction of slow voltage inactivation
    // End Available Variables and parameters for reference =======//|
    
    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	// voltage-activation gate
    set_ICaL_speciesCELL_MODEL_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

    // voltage inactivation gate
	set_ICaL_speciesCELL_MODEL_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

    //var->ICaL_vi_s_tau = f(Vm_inac_tau) or constant
	//var->ICaL_vi_s_tau        *= p.ICaL_vi_tau_scale;

	// calcium inactivation
	//var->ICaL_ci_ss		= f(Cai_x) or constant
	//var->ICaL_ci_j_ss		= f(Cai_x) or constant
	//var->ICaL_ci_tau      = f(Cai_x) or constant
    // or
	//var->ICaL_ci_al		= f(Cai_x) or constant
	//var->ICaL_ci_bet      = f(Cai_x) or constant
    // End Set gate rates =========================================//|
}

// voltage activation
void set_ICaL_speciesCELL_MODEL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	//*va_ss        = sigmoid(Vm_ss, V1/2, -gradient*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k) OR
    //*va_ss     = 1/(1 + exp((Vm_ss - V1/2)/(-k*kscale)) ); or more complex
	//*va_tau  		= f(Vm_tau)
}

// voltage inactivation
void set_ICaL_speciesCELL_MODEL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	//*vi_ss        = sigmoid(Vm_ss, V1/2, gradient*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k) OR:
    //*vi_ss        = 1/(1 + exp((Vm_ss - V1/2)/(-k*kscale)) ); or more complex
	//*vi_tau  		= f(Vm_tau)
}

// Everything else follows the same format
void update_gates_ICaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    //s->ICaL_va;             // voltage activation
    //s->ICaL_vi;             // voltage inactivation (fast if slow also used)
    //s->ICaL_vi_s;           // voltage inactivation, slow
    //s->ICaL_ci;             // calcium inactivation
    //s->ICaL_ci_j;           // calcium inactivation, junctional
    
    // Relevant computed variables:
    //var->ICaL_va_ss;          // voltage activation, steady state
    //var->ICaL_va_tau;         // voltage activation, time constant
    //var->ICaL_vi_ss;          // voltage inactivation, steady state
    //var->ICaL_vi_tau;         // voltage inactivation, time constant
    //var->ICaL_vi_s_tau;       // voltage inactivation, time constant, slow
    //var->ICaL_ci_ss;          // calcium inactivation, steady state
    //var->ICaL_ci_j_ss;        // calcium inactivation, steady state, junctional
    //var->ICaL_ci_tau;         // calcium inactivation, time constant
    //var->ICaL_ci_al;          // Calcium induced activation, alpha rate       (ms^-1)
    //var->ICaL_ci_bet;         // Calcium induced activation, beta rate        (ms^-1)

    // Calcium concentrations
    // s->Cai, s->Cai_sl, s->Cai_j
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->ICaL_va       		= rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
	//s->ICaL_vi    			= rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
	//s->ICaL_vi_s    			= rush_larsen(s->ICaL_vi_s, var->ICaL_vi_ss, var->ICaL_vi_s_tau, dt);
	//s->ICaL_ci    			= rush_larsen(s->ICaL_ci, var->ICaL_ci_ss, var->ICaL_ci_tau, dt);

    // Or explicitly, FE or other
    //s->ICaL_ci_{j/sl}              += dt*(differential= f( {ss, tau}/{alpha, beta},, Cai_{j/sl} ..))
}

void compute_ICaL_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, conductance/ICaL_bar and scale factor
    //s->ICaL_va;             // voltage activation
    //s->ICaL_vi;             // voltage inactivation (fast if slow also used)
    //s->ICaL_vi_s;           // voltage inactivation, slow
    //s->ICaL_ci;             // calcium inactivation
    //s->ICaL_ci_j;           // calcium inactivation, junctional
    //p.gCaL;  - conductance or:
    //var->ICaL_bar_sl, var->ICaL_bar_j
    //p.GCaL;  - scale factor
    //p.ICaL_vi_Fs - fraction slow voltage inactivation

    // Simple, HH-like, one inactivation gate and one compartment:
    //var->ICaL               = p.gCaL * s->ICaL_va * s->ICaL_vi * s->ICaL_ci * (Vm - 55);
	var->ICaL 				*= p.GCaL;

    // ICaL_bar, multiple inactivation gates and sl/j compartments:
    //var->ICaL_bar_sl        = compute_ICaL_bar_hAM_GB(p, var, Vm, s->Cai_sl, s->Cao); // function from lib/Model_hAM_GB and must be cited as such
    //var->ICaL_bar_j         = compute_ICaL_bar_hAM_GB(p, var, Vm, s->Cai_j, s->Cao);

    //double Op               = s->ICaL_va * ((1-p.ICaL_vi_Fs)*s->ICaL_vi + p.ICaL_vi_Fs*s->ICaL_vi_s);
    //var->ICaL_Ca_sl         =  (1 - p.Fjunc_ICaL) * p.pCaL * Op * s->ICaL_ci   * var->ICaL_bar_sl;
    //var->ICaL_Ca_j          =  (    p.Fjunc_ICaL) * p.pCaL * Op * s->ICaL_ci_j * var->ICaL_bar_j;
    var->ICaL_Ca_sl         *= p.GCaL;
    var->ICaL_Ca_j          *= p.GCaL;
    var->ICaL               = var->ICaL_Ca_sl + var->ICaL_Ca_j;
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
void set_IKur_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // Implementation for one activtion and one inactivation gating variable
    // Can be alpha/beta or steady-state/tau for ac and inac gates
    // Has functionality for a dynamically determined conductance parameter

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    // s->IKur_va;             // voltage activation
    // s->IKur_vi;             // voltage inactivation 

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //  p.IKur_va_ss_shift;    // Shift of the V1/2 of acitvation steady state     mV
    //  p.IKur_vi_ss_shift;    // Shift of the V1/2 of inactivation steady state   mV
    //  p.IKur_va_tau_shift;   // Shift of the voltage dependence of time constant mV
    //  p.IKur_vi_tau_shift;   // Shift of the voltage dependence of time constant mV
    //
    // tau scale
    //  p.IKur_va_tau_scale;   // Scales time constant for IKur voltage activation
    //  p.IKur_vi_tau_scale;   // Scales time constant for IKur voltage inactivation
    //
    // steady-state gradient scale
    //  p.IKur_va_ss_kscale;   // Scales the gradient parameter of activation steady state
    //  p.IKur_vi_ss_kscale;   // Scales the gradient parameter of inactivation steady state

    // Relevant computed variables:
    //var->IKur_va_ss;          // voltage activation, steady state
    //var->IKur_va_tau;         // voltage activation, time constant
    //var->IKur_vi_ss;          // voltage inactivation, steady state
    //var->IKur_vi_tau;         // voltage inactivation, time constant
    //var->IKur_va_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKur_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    //var->IKur_vi_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKur_vi_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    //var->IKur_dynamic_g;      // Dynamic conductance factor (s/mF)
    // End Available Variables and parameters for reference =======//|

    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Dynamic conductance factor
    //var->IKur_dynamic_g		= f(Vm) or other

    // simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

    //var->IKur_va_ss      = f(Vm_ac_ss, p.IKur_va_ss_kscale) e.g:
    //  var->IKur_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.IKur_va_ss_kscale);  OR:
    //  var->IKur_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.IKur_va_ss_kscale)) ); or more complex
    //var->IKur_va_tau    = f(Vm_ac_tau) or constant

    //var->IKur_vi_ss      = f(Vm_inac_ss, p.IKur_vi_ss_kscale) e.g:
    //  var->IKur_vi_ss      = sigmoid(Vm_inac_ss, V1/2, k*p.IKur_vi_ss_kscale);  OR:
    //  var->IKur_vi_ss      = 1/(1 + exp((Vm_inac_ss - V1/2)/(k*p.IKur_vi_ss_kscale)) ); or more complex
    //var->IKur_vi_tau    = f(Vm_inac_tau) or constant

    // alpha-beta style
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac        = Vm - p.IKur_va_ss_shift;  // Shift of the voltage used to calculate alpha and beta, activation
    double Vm_inac      = Vm - p.IKur_vi_ss_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

    // Set Activation gate alpha and beta
    //var->IKur_va_al                 = f(Vm_ac)
    //var->IKur_va_bet                = f(Vm_ac)

    // Set inactivation gates alphas and betas
    //var->IKur_vi_al           = f(Vm_inac)
    //var->IKur_vi_bet          = f(Vm_inac)

    // Set tau and SS from alpha and beta
    //var->IKur_va_tau                = 1.0/(var->IKur_va_al + var->IKur_va_bet); // 1/(a+b)
    //var->IKur_vi_tau              = 1.0/(var->IKur_vi_al + var->IKur_vi_bet);
    //var->IKur_va_ss                 = var->IKur_va_al * var->IKur_va_tau; // a*tau
    //var->IKur_vi_ss               = var->IKur_vi_al * var->IKur_vi_tau;
	
    // and don't forget to scale the time-constants!
	var->IKur_va_tau			*= p.IKur_va_tau_scale;
	var->IKur_vi_tau			*= p.IKur_vi_tau_scale;
    // End Set gate rates =========================================//|
}

void update_gates_IKur_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    // s->IKur_va;             // voltage activation
    // s->IKur_vi;             // voltage inactivation 

    // Relevant computed variables:
    //var->IKur_va_ss;          // voltage activation, steady state
    //var->IKur_va_tau;         // voltage activation, time constant
    //var->IKur_vi_ss;          // voltage inactivation, steady state
    //var->IKur_vi_tau;         // voltage inactivation, time constant
    //var->IKur_va_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKur_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    //var->IKur_vi_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKur_vi_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->IKur_va              = rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);

    // Or explicitly, FE or other
    //s->IKur_vi              += dt*(differential= f(ss, tau, alpha, beta..))
}

void compute_IKur_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potential, conductance and scale factors
    // s->IKur_va;             // voltage activation
    // s->IKur_vi;             // voltage inactivation 
    // var->EK
    // p.gKur, var->IKur_dynamic_g
    // p.GKur

    // Example with static conductance parameter
    // var->IKur                 = p.gKur * s->IKur_va * s->IKur_vi * (Vm - var->EK);

    // Example with dynamic conductance
	//var->IKur 				= var->IKur_dynamic_g * pow(s->IKur_va, 3) * s->IKur_vi * (Vm - var->EK);
    // Scale current
    var->IKur				*= p.GKur;
}
// End IKur =================================================================//|

// IKr ======================================================================\\|
void set_IKr_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // Implementation for one activtion and one inactivation dynamic gating variables, 
    // with one inactivation time-independent gating variable
    // Can be alpha/beta or steady-state/tau for ac and inac gates
    // Has functionality for fast and slow inactivation gates set by proportion

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    // s->IKr_va;              // voltage activation
    // s->IKr_vi;              // used where IKr is used for IKf in Rabbit atria || also used for va_slow

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //p.IKr_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.IKr_va_tau_shift;    // Shift of the voltage dependence of time constant mV
    //p.IKr_vi_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.IKr_vi_tau_shift;    // Shift of the voltage dependence of time constant mV
    
    // tau scale
    //p.IKr_va_tau_scale;    // Scales time constant for IKr voltage activation
    //p.IKr_vi_tau_scale;    // Scales time constant for IKr voltage inactivation
    
    // steady-state gradient scale
    //p.IKr_va_ss_kscale;    // Scales the gradient parameter of activation steady state
    //p.IKr_vi_ss_kscale;    // Scales the gradient parameter of activation steady state

    // Relevant computed variables:
    //var->IKr_va_ss;           // voltage activation, steady state
    //var->IKr_va_tau;          // voltage activation, time constant
    //var->IKr_va_al;           // alpha transition rate
    //var->IKr_va_bet;          // beta transition rate

    //var->IKr_vi_ss;           // voltage inactivation, steady state
    //var->IKr_vi_tau;          // voltage inactivation, time constant
    //var->IKr_vi_al;           // alpha transition rate
    //var->IKr_vi_bet;          // beta transition rate

    //var->IKr_vi_ti;           // Time-independent inactivation gate

    //var->IKr_va_Fs;           // Fraction of gates if using dynamic and static
    // End Available Variables and parameters for reference =======//|

    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
    double Vm_inac_ss       = Vm - p.IKr_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
    double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant
    double Vm_inac_tau      = Vm - p.IKr_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Time-independant gate
	//var->IKr_vi_ti          = f(Vm) or f(Vm_inac_ss)

    // simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

    //var->IKr_va_ss      = f(Vm_ac_ss, p.IKr_va_ss_kscale) e.g:
    //  var->IKr_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.IKr_va_ss_kscale);  OR:
    //  var->IKr_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.IKr_va_ss_kscale)) ); or more complex
    //var->IKr_va_tau    = f(Vm_ac_tau) or constant

    //var->IKr_vi_ss      = f(Vm_inac_ss, p.IKr_vi_ss_kscale) e.g:
    //  var->IKr_vi_ss      = sigmoid(Vm_inac_ss, V1/2, k*p.IKr_vi_ss_kscale);  OR:
    //  var->IKr_vi_ss      = 1/(1 + exp((Vm_inac_ss - V1/2)/(k*p.IKr_vi_ss_kscale)) ); or more complex
    //var->IKr_vi_tau    = f(Vm_inac_tau) or constant

    // alpha-beta style
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac        = Vm - p.IKr_va_ss_shift;  // Shift of the voltage used to calculate alpha and beta, activation
    double Vm_inac      = Vm - p.IKr_vi_ss_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

    // Set Activation gate alpha and beta
    //var->IKr_va_al                 = f(Vm_ac)
    //var->IKr_va_bet                = f(Vm_ac)

    // Set inactivation gates alphas and betas
    //var->IKr_vi_al           = f(Vm_inac)
    //var->IKr_vi_bet          = f(Vm_inac)

    // Set tau and SS from alpha and beta
    //var->IKr_va_tau                = 1.0/(var->IKr_va_al + var->IKr_va_bet); // 1/(a+b)
    //var->IKr_vi_tau              = 1.0/(var->IKr_vi_al + var->IKr_vi_bet);
    //var->IKr_va_ss                 = var->IKr_va_al * var->IKr_va_tau; // a*tau
    //var->IKr_vi_ss               = var->IKr_vi_al * var->IKr_vi_tau;
    
    // and don't forget to scale the time-constants!
	var->IKr_va_tau			*= p.IKr_va_tau_scale;
	var->IKr_vi_tau			*= p.IKr_vi_tau_scale;
    // End Set gate rates =========================================//|
}

void update_gates_IKr_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    // s->IKr_va;             // voltage activation
    // s->IKr_vi;             // voltage inactivation

    // Relevant computed variables:
    //var->IKr_va_ss;          // voltage activation, steady state
    //var->IKr_va_tau;         // voltage activation, time constant
    //var->IKr_vi_ss;          // voltage inactivation, steady state
    //var->IKr_vi_tau;         // voltage inactivation, time constant
    //var->IKr_va_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKr_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    //var->IKr_vi_al;          // voltage activation, alpha transition rate (1-y-> y)
    //var->IKr_vi_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->IKr_va              = rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);

    // Or explicitly, FE or other
    //s->IKr_vi              += dt*(differential= f(ss, tau, alpha, beta..))
}

void compute_IKr_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potential, conductance and scale factors
    // s->IKr_va;             // voltage activation
    // s->IKr_vi;             // voltage inactivation
    // var->EK
    // p.gKr
    // p.GKr
    // var->IKr_vi_ti;           // Time-independent inactivation gate
    // var->IKr_va_Fs;           // Fraction of gates if using dynamic and static

    // Example with dynamic activation and static gate
    // var->IKr                  = p.gKr * s->IKr_va * var->IKr_vi_ti * (Vm - var->EK);

    // Using fraction of channels (as in ORd)
    // var->IKr = p.gKr	* ((1-var->IKr_va_Fs)*s->IKr_va + var->IKr_va_Fs*s->IKr_vi) * var->IKr_vi_ti * (Vm - var->EK);

    // Scale current
    var->IKr *= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // Implementation for two activation variables (one ss, 2 taus; one tau scale)
    // Can be alpha/beta or steady-state/tau
    // implementation for sl/j compartementalisation 

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    // s->IKs_va;              // voltage activation
    // s->IKs_va_2;            // voltage activation

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //p.IKs_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.IKs_va_tau_shift;    // Shift of the voltage dependence of time constant mV
    
    // tau scale
    //p.IKs_va_tau_scale;    // Scales time constant for IKs voltage activation
    
    // steady-state gradient scale
    //p.IKs_va_ss_kscale;    // Scales the gradient parameter of activation steady state

    // Relevant computed variables:
    //var->IKs_va_ss;           // voltage activation, steady state
    //var->IKs_va_tau;          // voltage activation, time constant
    //var->IKs_va_2_tau;        // voltage activation, time constant
    //var->IKs_va_al;           // alpha transition rate
    //var->IKs_va_bet;          // beta transition rate
    // End Available Variables and parameters for reference =======//|

    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
    double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

    // simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

    //var->IKs_va_ss      = f(Vm_ac_ss, p.IKs_va_ss_kscale) e.g:
    //  var->IKs_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.IKs_va_ss_kscale);  OR:
    //  var->IKs_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.IKs_va_ss_kscale)) ); or more complex
    //var->IKs_va_tau    = f(Vm_ac_tau) or constant
    //var->IKs_va_2_tau    = f(Vm_ac_tau) or constant

    // alpha-beta style
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac        = Vm - p.IKs_va_ss_shift;  // Shift of the voltage used to calculate alpha and beta, activation

    // Set Activation gate alpha and beta
    //var->IKs_va_al                 = f(Vm_ac)
    //var->IKs_va_bet                = f(Vm_ac)

    // Set tau and SS from alpha and beta
    //var->IKs_va_tau                = 1.0/(var->IKs_va_al + var->IKs_va_bet); // 1/(a+b)
    //var->IKs_va_ss                 = var->IKs_va_al * var->IKs_va_tau; // a*tau
    
    // and don't forget to scale the time-constants!
	var->IKs_va_tau			    *= p.IKs_va_tau_scale;
	var->IKs_va_2_tau			*= p.IKs_va_tau_scale;
    // End Set gate rates =========================================//|

}

void update_gates_IKs_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    // s->IKs_va;             // voltage activation
    // s->IKs_va_2;           // voltage inactivation

    // Relevant computed variables:
    //var->IKs_va_ss;           // voltage activation, steady state
    //var->IKs_va_tau;          // voltage activation, time constant
    //var->IKs_va_2_tau;        // voltage activation, time constant
    //var->IKs_va_al;           // alpha transition rate
    //var->IKs_va_bet;          // beta transition rate
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->IKs_va              = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
    //s->IKs_va_2            = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_2_tau, dt);

    // Or explicitly, FE or other
    //s->IKs_va_2              += dt*(differential= f(ss, tau, alpha, beta..))
}

void compute_IKs_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potentials, conductance and scale factors
    // s->IKs_va;             // voltage activation
    // s->IKs_va_2;             // voltage inactivation
    // Multiple reversal potential options (all computed by default in Model.c so use any)
    // var->EK - potassium only
    // var->EKs - potassium + sodium: ((p.R * p.T)/p.F)*log((s->Ko + 0.01833*s->Nao)/(s->Ki + 0.01833*s->Nai));
    // var->EKs_ORD - potassium + sodium alt: ((p.R * p.T)/p.F)*log((s->Ko + s->Nao)/(s->Ki + 0.01833*s->Nai));
    // p.gKs
    // p.GKs

    // Single compartment
    // var->IKs                  = p.gKs * s->IKs_va * s->IKs_va_2 * (Vm - var->EKs);
    // Scale current
    var->IKs *= p.GKs;

    // Multi
    // var->IKs                =  p.gKs * s->IKs_va * s->IKs_va * (Vm - var->EK);
    var->IKs                *= p.GKs;
    var->IKs_sl             = (1 - p.Fjunc)*var->IKs;
    var->IKs_j              = (    p.Fjunc)*var->IKs;
}
// End IKs ==================================================================//|

// IKACh ====================================================================\\|
void set_IKACh_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // Implementation for one activtion and one inactivation dynamic gating variables, 
    // with one time-independent gating variable
    // steady-state/tau only

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    // s->IKACh_va;              // voltage activation
    // s->IKACh_vi;              // 

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //p.IKACh_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.IKACh_va_tau_shift;    // Shift of the voltage dependence of time constant mV
    //p.IKACh_vi_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.IKACh_vi_tau_shift;    // Shift of the voltage dependence of time constant mV
    
    // tau scale
    //p.IKACh_va_tau_scale;    // Scales time constant for IKACh voltage activation
    //p.IKACh_vi_tau_scale;    // Scales time constant for IKACh voltage inactivation
    
    // steady-state gradient scale
    //p.IKACh_va_ss_kscale;    // Scales the gradient parameter of activation steady state
    //p.IKACh_vi_ss_kscale;    // Scales the gradient parameter of activation steady state

    // Relevant computed variables:
    //var->IKACh_va_ss;           // voltage activation, steady state
    //var->IKACh_va_tau;          // voltage activation, time constant

    //var->IKACh_vi_ss;           // voltage inactivation, steady state
    //var->IKACh_vi_tau;          // voltage inactivation, time constant

    //var->IKACh_v_ti;           // Time-independent gate
    // End Available Variables and parameters for reference =======//|

    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.IKACh_va_ss_shift;   // Voltage modified by shift applied to activation steady state
    double Vm_inac_ss       = Vm - p.IKACh_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
    double Vm_ac_tau        = Vm - p.IKACh_va_tau_shift;  // Voltage modified by shift applied to activation time constant
    double Vm_inac_tau      = Vm - p.IKACh_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Time-independant gate
	//var->IKACh_v_ti          = f(Vm)

    // simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

    //var->IKACh_va_ss      = f(Vm_ac_ss, p.IKACh_va_ss_kscale) e.g:
    //  var->IKACh_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.IKACh_va_ss_kscale);  OR:
    //  var->IKACh_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.IKACh_va_ss_kscale)) ); or more complex
    //var->IKACh_va_tau    = f(Vm_ac_tau) or constant

    //var->IKACh_vi_ss      = f(Vm_inac_ss, p.IKACh_vi_ss_kscale) e.g:
    //  var->IKACh_vi_ss      = sigmoid(Vm_inac_ss, V1/2, k*p.IKACh_vi_ss_kscale);  OR:
    //  var->IKACh_vi_ss      = 1/(1 + exp((Vm_inac_ss - V1/2)/(k*p.IKACh_vi_ss_kscale)) ); or more complex
    //var->IKACh_vi_tau    = f(Vm_inac_tau) or constant
    
    // and don't forget to scale the time-constants!
	var->IKACh_va_tau			*= p.IKACh_va_tau_scale;
	var->IKACh_vi_tau			*= p.IKACh_vi_tau_scale;
    // End Set gate rates =========================================//|
}

void update_gates_IKACh_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    // s->IKACh_va;             // voltage activation
    // s->IKACh_vi;             // voltage inactivation

    // Relevant computed variables:
    //var->IKACh_va_ss;          // voltage activation, steady state
    //var->IKACh_va_tau;         // voltage activation, time constant
    //var->IKACh_vi_ss;          // voltage inactivation, steady state
    //var->IKACh_vi_tau;         // voltage inactivation, time constant
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->IKACh_va              = rush_larsen(s->IKACh_va, var->IKACh_va_ss, var->IKACh_va_tau, dt);

    // Or explicitly, FE or other
    //s->IKACh_vi              += dt*(differential= f(ss, tau..))
}

void compute_IKACh_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potential, conductance and scale factors
    // s->IKACh_va;             // voltage activation
    // s->IKACh_vi;             // voltage inactivation
    // var->EK
    // p.gKACh, p.gKACh_max
    // p.GKACh
    // var->IKACh_v_ti;           // Time-independent inactivation gate

    // Example with dynamic activation and static gate
    // var->IKACh                  = p.gKACh * s->IKACh_va * s->IKACh_vi_ti * (Vm - var->EK);

    // Scale current
    var->IKACh *= p.GKACh;
}
// End IKACh ================================================================//|

// IK1 ======================================================================\\|
// IK1 can be implemented in multiple ways:
// 1 - single, compute IK1 function
// 2 - set variables + compute Ik1 ; this allows the computed part to be computed for lookup table separately to the current
// 3 - with time-depdendent state variable (not added to template; see ORd implementation (lib/Model_hVM_ORd_simple.cc) for example)

void set_IK1_speciesCELL_MODEL_variables(Cell_parameters p, Model_variables *var, double Vm)
{
    // Relevant modulation variables
    // p.IK1_va_shift // voltage shift
    
    // Relevant computed variables
    //var->IK1_va_ti    // time independent voltage dependant factor
	
    // Voltage + shift
    double Vm_in        = Vm - p.IK1_va_shift;
	//var->IK1_va_ti      = f(Vm_in)
}

void compute_IK1_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant modulation variables
    // p.IK1_Erev_shift       // shift of the V term in V-Ek

    double Vm_in         = Vm - p.IK1_Erev_shift; // if using Erev v shift

    // Simmple, signle function
    // var->IK1               = p.gK1 * (Vm_in - var->EK)*var->IK1_va_ti;

    // using pre computed variable:
	//var->IK1               = p.gK1 * (Vm_in - var->EK)/var->IK1_va_ti;

    // Scale!
	var->IK1               *= p.GK1;
}
// End IK1 ==================================================================//|

// If ======================================================================\\|
void set_If_speciesCELL_MODEL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // Implementation for one activation gate
    // steady-state/tau

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (not used in this function, but the gating
    // variables for which rates/time-constants/steady-states must be set)
    // s->If_va;              // voltage activation

    // Relevant modulation parameters (which are defaulted to 1 (scale) or 0 (shift) and set in het and modulation functions):
    // Voltage shifts
    //p.If_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
    //p.If_va_tau_shift;    // Shift of the voltage dependence of time constant mV
    
    // tau scale
    //p.If_va_tau_scale;    // Scales time constant for If voltage activation
    
    // steady-state gradient scale
    //p.If_va_ss_kscale;    // Scales the gradient parameter of activation steady state

    // Relevant computed variables:
    //var->If_va_ss;           // voltage activation, steady state
    //var->If_va_tau;          // voltage activation, time constant
    // End Available Variables and parameters for reference =======//|

    // Set gate rates =============================================\\|
    // -> assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    double Vm_ac_ss         = Vm - p.If_va_ss_shift;   // Voltage modified by shift applied to activation steady state
    double Vm_ac_tau        = Vm - p.If_va_tau_shift;  // Voltage modified by shift applied to activation time constant

    // simgoid function for use:  sigmoid(V, V1/2, k*p._kscale)  =  1/(1 + exp( (V - V1/2)/K ) )

    //var->If_va_ss      = f(Vm_ac_ss, p.If_va_ss_kscale) e.g:
    //  var->If_va_ss      = sigmoid(Vm_ac_ss, V1/2, -k*p.If_va_ss_kscale);  OR:
    //  var->If_va_ss      = 1/(1 + exp((Vm_ac_ss - V1/2)/(-k*p.If_va_ss_kscale)) ); or more complex
    //var->If_va_tau    = f(Vm_ac_tau) or constant

    // and don't forget to scale the time-constants!
	var->If_va_tau			    *= p.If_va_tau_scale;
    // End Set gate rates =========================================//|
}

void update_gates_If_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables - update all which are used
    // s->If_va;             // voltage activation

    // Relevant computed variables:
    //var->If_va_ss;           // voltage activation, steady state
    //var->If_va_tau;          // voltage activation, time constant
    //var->If_va_2_tau;        // voltage activation, time constant
    // End Available Variables and parameters for reference =======//|

    // Update by rush-larsen                variable     steady-state  time-constant   dt
    //s->If_va              = rush_larsen(s->If_va, var->If_va_ss, var->If_va_tau, dt);

    // Or explicitly, FE or other
    //s->If_va              += dt*(differential= f(ss, tau..))
}

void compute_If_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Relevant state variables, reversal potentials, conductance and scale factors
    // s->If_va;             // voltage activation
    // var->EK, var->ENa
    // p.gf
    // p.Gf

    // Single compartment
    // var->If                  = p.gf * f(s->If_va, var->EK, var->ENa)
    // Scale current
    var->If *= p.Gf;
}
// End If ==================================================================//|

// Ca2+ handling, background and pump currents ==============================\\|
// INCX
void compute_INCX_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single current, or SL/j compartments

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (concentrations):
    // s->Cai, s->Cao, s->Nai, s->Nao
    // s->Cai_sl, s->Cai_j, s->Nai_sl, s->Nai_j

    // Relevant parameters (ensure which are used are set in set parameters!)
    //p.INCX_bar;            // Maximal current scale factor for NCX         (pA/pF) or (um^3.uM.ms^-1)
    //p.INCX_kNao;           // [Na+]o saturation constant for INCX          (mM)
    //p.INCX_kNai;           // [Na+]i saturation constant for INCX          (mM)
    //p.INCX_kCao;           // [Ca2+]o saturation constant for INCX         (mM)
    //p.INCX_kCai;           // [Ca2+]i saturation constant for INCX         (mM)    or (uM)
    //p.INCX_k;              // Saturation factor for INCX                   (dimensionless)
    //p.INCX_gamma;          // Voltage dependence factor for INCX           (dimensionless)
    //p.INCX_kda;            // Magnitude calcium saturation constant        (mM)    or (uM)
    //p.INCX_ksat;           // Saturation constant                          (mM)

    // p.Fjunc;

    // Modulation:
    // p.GNCX - current scale
    // End Available Variables and parameters for reference =======//|
	
    // Single compartment example
    // Create some local variables to make equation easier to read
	//double Cai 		 	= s->Cai;  
	//double Cao			= s->Cao;
	//double Nai			= s->Nai;
	//double Nao			= s->Nao;

	//var->INCX			= p.INCX_bar *  f(Cai, Cao, Nai, Nao, Vm, p.INCX_kNao, p.INCX_gamma, p.INCX_kda, ....)

    // Scale
	var->INCX			*= p.GNCX;

    // Multi compartment
    //double Cao            = s->Cao;
    //double Nao            = s->Nao;
    

    // Junctional
    //double Cai            = s->Cai_j;
    //double Nai            = s->Nai_j;
    //var->INCX_j         = p.Fjunc * p.INCX_bar * f(Cai, Cao, Nai, Nao, Vm, p.INCX_kNao, p.INCX_gamma, p.INCX_kda, ....)

    // sl
    //double Cai            = s->Cai_sl;
    //double Nai            = s->Nai_sl;
    //var->INCX_sl         = (1 - p.Fjunc) * p.INCX_bar * f(Cai, Cao, Nai, Nao, Vm, p.INCX_kNao, p.INCX_gamma, p.INCX_kda, ....)

    // Scale and total
    var->INCX_sl        *= p.GNCX;
    var->INCX_j         *= p.GNCX;
    var->INCX           = var->INCX_sl + var->INCX_j;
}

// INaK
void compute_INaK_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single current, or SL/j compartments

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (concentrations):
    // s->Ki, s->Ko, s->Nai, s->Nao
    // s->Nai_sl, s->Nai_j

    // Parameters:
    // p.INaK_bar;            // Maximal current scale factor for INaK        (pA/pF)
    // p.INaK_kK;             // Half-Saturation constant for Ko INaK         (mM)
    // p.INaK_kNa;            // Half-saturation constant for Nai INaK        (mM)
    // p.Fjunc
    
    // Modulation:
    // p.GNaK
    // End Available Variables and parameters for reference =======//|

    // Possible local variables for clarity
    // double sigma = f(Nao) e.g. (exp(s->Nao/67.3)-1.0)/7.0;
    // double FNaK  = f(Vm, F, R, T) e.g. pow(1.0+0.1245*exp(-0.1*p.F*Vm/(p.R*p.T))+0.0365*sigma*exp(-Vm*p.FoRT ), -1.0);

    // Single compartment
    // var->INaK    = p.INaK_bar * f(Vm, concentrations (Nai), sigma, FNaK, p.INaK_kK, p.INaK_kNa)

    // Scale
    var->INaK *= p.GNaK; 

    // Multi-compartment
    // var->INaK_j    = p.Fjunc * p.INaK_bar * f(Vm, concentrations (Nai_j), sigma, FNaK, p.INaK_kK, p.INaK_kNa)
    // var->INaK_sl    = (1-p.Fjunc) * p.INaK_bar * f(Vm, concentrations (Nai_sl), sigma, FNaK, p.INaK_kK, p.INaK_kNa)

    // Scale and total
    var->INaK_sl        *= p.GNaK;
    var->INaK_j         *= p.GNaK;
    var->INaK           = var->INaK_sl + var->INaK_j;
}

// ICaP
void compute_ICaP_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single current, or SL/j compartments

    // Available Variables and parameters for reference ===========\\|
    // Relevant state variables: (concentrations):
    // s->Cai, s->Cai_j, s->Cai_sl

    // Parameters:
    // p.ICaP_bar;            // Conductance factor for calcium pump (pA/pF)  or (um^3 . uM . ms^-1)
    // p.ICaP_kCa;            // Saturation constant for calcium pump (mM)    or (uM)
    // p.Fjunc

    // Modulation
    // p.GCaP   // scale
    // End Available Variables and parameters for reference =======//|

    // var->ICaP    = f(s->Cai, Vm)
    
    // Simple, single compartment
    // var->ICaP           = (p.ICaP_bar * s->Cai) / (s->Cai + p.ICaP_kCa);
    // Scale
    var->ICaP           *= p.GCaP;

    // Multi compartment example
    //var->ICaP_sl        = (1 - p.Fjunc) * p.ICaP_bar * pow(s->Cai_sl, 1.6)/(pow(p.ICaP_kCa, 1.6) + pow(s->Cai_sl, 1.6));
    //var->ICaP_j         = (    p.Fjunc) * p.ICaP_bar * pow(s->Cai_j, 1.6)/(pow(p.ICaP_kCa, 1.6) + pow(s->Cai_j, 1.6));
    var->ICaP_sl        *= p.GCaP;
    var->ICaP_j         *= p.GCaP;
    var->ICaP           = var->ICaP_sl + var->ICaP_j;
}

// INab
void compute_INab_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single current, or SL/j compartments

    // Available Variables and parameters for reference ===========\\|
    // Parameters:
    // p.gNab;                // Conductance of background sodium current     (s/mF)
    // var->ENa, var->ENa_sl, var->ENa_j
    // p.Fjunc

    // Modulation
    // p.GNab   // scale
    // End Available Variables and parameters for reference =======//|

    // Single compartment
    //var->INab           = p.gNab * (Vm - var->ENa);
    var->INab           *= p.GNab;

    // Multi compartments
    //var->INab_sl        = (1 - p.Fjunc) * p.gNab * (Vm - var->ENa_sl);
    //var->INab_j         = (    p.Fjunc) * p.gNab * (Vm - var->ENa_j);
    var->INab_sl        *= p.GNab;
    var->INab_j         *= p.GNab;
    var->INab           = var->INab_sl + var->INab_j;
}

// ICab
void compute_ICab_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single current, or SL/j compartments

    // Available Variables and parameters for reference ===========\\|
    // Parameters:
    // p.gCab;                // Conductance of background sodium current     (s/mF)
    // var->ECa, var->ECa_sl, var->ECa_j
    // p.Fjunc

    // Modulation
    // p.GCab   // scale
    // End Available Variables and parameters for reference =======//|

    // Single compartment
    //var->ICab           = p.gCab * (Vm - var->ECa);
    var->ICab           *= p.GCab;

    // Multi compartments
    //var->ICab_sl        = (1 - p.Fjunc) * p.gCab * (Vm - var->ECa_sl);
    //var->ICab_j         = (    p.Fjunc) * p.gCab * (Vm - var->ECa_j);
    var->ICab_sl        *= p.GCab;
    var->ICab_j         *= p.GCab;
    var->ICab           = var->ICab_sl + var->ICab_j;
}

// IKb
void compute_IKb_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single compartment
    //
    // Available Variables and parameters for reference ===========\\|
    // Parameters:
    // p.gKb;                // Conductance of background potassium current     (s/mF)
    // var->EK

    // Modulation
    // p.GKb // scale
    // End Available Variables and parameters for reference =======//|

    // Simple eg
    //var->IKb    = p.gKb*(Vm - var->EK);

    // or more complex
    //var->IKb    = p.gKb*(vm - var->EK)* f(Vm);

    // Scale
    var->IKb *= p.GKb;
}

// IClCa
void compute_IClCa_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single compartment or sl/j compartments
    //
    // Available Variables and parameters for reference ===========\\|
    // Concentrations:
    // s->Cai, s->Cai_slm s->Cai_j

    // Parameters:
    // p.gKb;                // Conductance of Ca2+-activated chloride       (s/mF)
    // var->ECl
    // p.IClCa_kd;            // Saturation constant for Ca2+-activated Cl current
    // p.Fjunc

    // Modulation
    // p.GClCa // scale
    // End Available Variables and parameters for reference =======//|

    // Single implementation
    // var->IClCa       = p.gClCa * (vm - var->ECl) * f(Vm, p.IClCa_kd, s->Cai)

    // Scale 
    var->IClCa *= p.GClCa;

    // Multi compartment
    //var->IClCa_sl       = (1 - p.Fjunc) * p.gClCa * f(Vm, p.IClCa_kd, s->Cai_sl)
    //var->IClCa_j        = (    p.Fjunc) * p.gClCa * f(Vm, p.IClCa_kd, s->Cai_j)
    var->IClCa_sl       *= p.GClCa;
    var->IClCa_j        *= p.GClCa;
    var->IClCa          = var->IClCa_sl + var->IClCa_j;
}

// IClb
void compute_IClb_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    // Single function; create variables which don't need to be used elsewhere locally
    // Implementation for single compartment
    //
    // Available Variables and parameters for reference ===========\\|
    // Parameters:
    // p.gClb
    // var->ECl
    
    // Modifiers:
    // p.GClb
    // End Available Variables and parameters for reference =======//|

    //var->IClb = p.gClb*(Vm 0 var->ECl);
    var->IClb   *= p.GClb;
}
// End Ca2+ handling, background and pump currents ==========================//|

// Homeostasis ==============================================================\\|
void comp_homeostasis_speciesCELL_MODEL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    // (Some) available variables and parameters ==================\\|
    // Please see the Ca2+ handling elements in the relevant structs (lib/Structs.h) for a full list, 
    // which includes variables and parameters for multiple different models. This should cover most
    // requirements. Use a similar already implemented model (if available) as a guide, and/or repurpose
    // variable and parameters where appropriate.
    // Here only basic parameters and variables are listed

	// Currents
    // var->ICaL, var->INCX or var->INCX_sl/j etc

	// Differentials
	//var->dNai 			
	//var->dKi 			
    //var->dCai

    // Concentration state variables
    // s->Cai, s->Cai_sl, s->Cai_j, s->Cao, s->Nai, s->Nai_sl, s->Nai_j, s->Nao, s->Ki, s->Ko
    // s->CajSR, s->CanSR

    //Modifiers
    //p.Gup;                 // Scale factor for intracellular, SERCA Ca2+ uptake
    //p.Gleak;               // Scale factor for intracellular Ca2+ leak
    //p.Grel;                // Scale factor for intracellualr Ca2+ release
    
    // Cell structure parameters:
    //p.Cm;                  // Membrane capacitance                         (pF)
    //p.Cm_F;                // Membrane capacitance                         (F)

    //p.Vcell;               // volume of whole cell (um^3)
    //p.Vcyto;               // volume of intracellular Ca2+ space
    //p.VjSR;                // volume of jSR / release Ca2+ space
    //p.VnSR;                // volume of nSR / uptake Ca2+ space
    //p.Vjunc;               // Volume of intracellular junctional release space
    //p.Vsl;                 // Volume of intracellular sub-sarcolemmal space
    //p.Vc;                  // Volume of extracellular cleft space

    //p.Fjunc;               // Proportion of currents in junctional vs non-junctional compartments
    //p.Fjunc_ICaL;          // Proportion of currents in junctional vs non-junctional compartments, LTCCs/ICaL

    // Ca2+ handling parameters
    //p.J_rel_max;           // Maximal flux rate of intracellular Ca2+ release  (mM/ms) or (um^3/ms)
    //p.J_SERCA_max;         // Maximal flux rate of intracellular Ca2+ uptake   (mM/ms) or (uM/ms)
    //p.J_SERCA_kCa;         // SERCA Cai constant                               (mM)    or (uM)
    //p.J_SERCA_kCaSR;       // SERCA CaSR constant                              (uM)
    //p.J_leak_max;          // Maximal flux rate of intracellular Ca2+ leak     (mM/ms) or ms^-1
    //p.J_leak_kCaSR;        // Jleak CaSR2+ constant                            (-)     or uM
    //p.J_jsr_nsr_tau;       // Time constant of transfer between jsr and nsr    (ms)
    // End (Some) available variables and parameters ===============//|

	// Ca2+ handling
    // As models have different structures, components and variables, there are multiple options here.
    // Many variables as named in this code are general and common to differnet implementations 
    // (e.g. the RyR gates RyRo, RyRi and RyRr can represent any RyR gating variables, 
    // SR Ca concentrations CajSR and CanSR could represent Caup, Carel, or single CanSR, 
    // Cai_j could represent dyadic cleft etc)
    // so these should be used where they can.
    // Create variables which don't need to be accessed elsewhere locally here

	// RyR / Jrel =============================\\|
    // Gating/state variables:
    // s->RyRo;                // proportion open RyRs / activation gate or other
    // s->RyRr;                // proportion refactory RyRs / inactivation gate 1 or other
    // s->RyRi;                // proportion inactivated RyRs / inactvation gate 2 or other
	//var->J_rel = p.J_rel_max * f(s->RyR, Ca)
    var->J_rel *= p.Grel; 	// Scale for MODIFIERS	

	//s->RyRo	= 
	//s->RyRr	= 
	//s->RyRi	= 
	//s->RyRo 	+= dt*(differential)
	//s->RyRr 	+= dt*(differential)
	//s->RyRi 	+= dt*(differential)
	// End RyR / Jrel =========================//|

	// Jup/leak ===============================\\|
	//var->J_SERCA = p.J_SERCA_max * f(e.g. s->Cai, s->CanSR, p.J_SERCA_kCa, p.J_SERCA_kCaSRk).
	var->J_SERCA	*= p.Gup; // Scale for MODIFIERS
	//var->J_leak 	= p.J_leak_max * f(e.g. s->CanSR, s->CajSR, p.J_leak_kCaSR)
	var->J_leak	*= p.Gleak; // Scale for MODIFIERS
	// End Jup/leak ===========================//|

	// Buffering ==============================\\| 
    // some parameters:
    //p.Kcam;                // uM
    //p.Bcam;                // uM
    //p.Kbsr;                // uM
    //p.Bbsr;                // uM
    //p.Kmca;                // uM
    //p.Bmca;                // uM
    //p.Kmmg;                // uM
    //p.Bmmg;                // uM

    //p.Kcsqn;               // mM
    //p.Bcsqn;               // mM

    //p.cmdnbar;             // Total calmodulin concentration
    //p.trpnbar;             // Total troponin concentration
    //p.csqnbar;             // Total csqn concentration
    //p.cmdn_k;              // half-saturation constant for calmodulin  (mM)
    //p.trpn_k;              // half-saturation constant for troponin    (mM)
    //p.csqn_k;              // half-saturation constant for calsequestrin (mM)

    // Some state variables:
    //s->cmdn;                // calmodulin concentration     (mM)
    //s->trpn;                // troponin concentration       (mM)
    //s->csqn;                // calsequestrin concenrration  (mM)
    //s->Tn_CHm;
    //s->Tn_CHc;
    //s->Myo_m;
    //s->Myo_c;
    //s->Tn_CL;

    //double beta_X = 
	// Buffering ==============================//| 

	// transfer equations etc etc
    // var->J_jsr_nsr          = (s->CanSR - s->CajSR)/p.J_jsr_nsr_tau;

	// Differentials
	//var->dCai 				= beta_X*f(Cm, Ca, Ix)
	//var->dKi		    		= f(Cm, Ix)
	//var->dNai		    		= f(Cm, Ix)

	// Update concentrations
	//s->Nai  				+= dt* var->dNai;
	//s->Ki  				+= dt* var->dKi;
	//s->Cai  				+= dt* var->dCai;
}
// End Homeostasis ==========================================================//|
// End Current formulations =====================================================================//|

// Heterogeneity and modulation =================================================================\\|
// Here, add model-specific heterogeneity and modulation implementations.
// If you implementation does, or is likley to apply to multiple models 
// (e.g. all human atrial models), then add it as a function in lib/Model.c
// "set_heterogeneity_and_modulation_native()" and follow instructions for global or common
// modulation

// !! MUST BE CALLED in lib/Model.c -> set_heterogeneity_and_modulation_native()
void set_het_mod_speciesCELL_MODEL(Cell_parameters *p)
{
    // Here must either call model-specific function or return an error if modulation is set but does not exist
    // Select per modulation type
    
    // Note: Need the "p->X_set_ref == 0" to check if modulation has already been set
    // globally in Model.c; if so, we don't want to overwrite or return error;
    // you could remove this clause to invalidate your model from any global moduation.

    // No settings exist
    if (strcmp(p->Celltype, "default") != 0 && p->Het_set_ref == 0) 		{ printf("ERROR: speciesCELL_MODEL has no function for heterogeneity\n"); exit(1); }
    if (p->ISO > 0.0                        && p->ISO_set_ref == 0) 		{ printf("ERROR: speciesCELL_MODEL has no function for ISO, but ISO is non-zero\n"); exit(1); }
    if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0) 		{ printf("ERROR: speciesCELL_MODEL has no function for Agent, but an Agent is set\n"); exit(1); }
    if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) { printf("ERROR: speciesCELL_MODEL has no function for Remodelling, but a Remodelling is set\n"); exit(1); }
    if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    { printf("ERROR: speciesCELL_MODEL has no function for Mutation, but a Mutation is set\n"); exit(1); }
    if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         { printf("ERROR: speciesCELL_MODEL has no function for ACh, but ACh is non-zero\n"); exit(1); }

    // OR settings do exist; call function	
    if (                                       p->Het_set_ref == 0)         set_celltype_native_speciesCELL_MODEL(p);
    if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_speciesCELL_MODEL(p);
    if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_speciesCELL_MODEL(p);
    if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_speciesCELL_MODEL(p);
    if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_speciesCELL_MODEL(p);
    if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_speciesCELL_MODEL(p);
}

// any updates for integrated model versions
// !! MUST BE CALLED in lib/Model.c -> update_heterogeneity_and_modulation_integrated()
void update_het_and_mod_speciesCELL_MODEL_integrated(Cell_parameters *p)
{
    // Add function calls here to "update_celltype/modulation" etc.
    // Ensure you add the new function
    // Function should follow the format of set_native, except only have entries
    // where thigns are updated, and doesn't need the error checking\
    
    // Note: DON'T want set refs here, as it will have already been set for this to be an update!!!

    // e.g.
    if (strcmp(p->Celltype, "default") != 0) update_celltype_integrated_speciesCELL_MODEL(p);
}

// Celltype/Heterogeneity
void set_celltype_native_speciesCELL_MODEL(Cell_parameters *p)
{
    // ================================================================================\\|
    // Can modify parameters (e.g. conductance) directly, or scaling factors.
    // (conductance denoted "g", scale factor "G"; current = f(g*G))
    // In general, use the scale factors as these are multiplicative/additive throughout
    // the code. 
    // Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // But overwrites any previous settings absolutely
    // ================================================================================//|

    // NOTE: relevant scale factors/shifts for each current are found in the current's function
    // If adding global style here, ensure all relevant models have the scale-factor/shift
    // and use it the same way!

    // NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

    if (strcmp(p->Celltype, "some_region") == 0); 	// Do nothing for baseline celltype; ensure this is what is set in set_parameters_native_species()

    // testing exmaple illustration of model-specific implementation
    //else if (strcmp(p->Celltype, "test") == 0) 
    //{
    //    p->Gto				*= 2; 		// Scale factor = Multiplies previous settings
    //    p->gKur				= 0.003;	// Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. Don't use this style unless intentional to absolutely set 
    //    p->IKr_va_ss_kscale	*= 1.25;	// Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
    //    p->ICaL_vi_ss_shift	+= 5;		// Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
    //    p->IKs_va_tau_scale	*= 0.75;	// Multiplies time constant of voltage activation. Multiplies previous settings.
    //    p->Gup				*= 2;		// Scales intracellular upatke rate. Multiplies previous settings.
	//    p->GCaL, p->Grel *= x ; // use GCaL and Grel here, as celltypes is about expression not activity
    //    p->ko                 = 5;        // explicitly set extracellular potassium (or other parameter)
    //}
    // Add new celltypes here: else if (strcmp(p->Celltype, "X") == 0) {   }
    else 
    {
        printf("ERROR: \"%s\" is not a valid Celltype for the speciesCELL_MODEL models. Please check Model.c and Model_speciesCELL_MODEL.cpp for options\n\n", p->Celltype);
        exit(1);
    }
}

// Any integrated updates here
void update_celltype_integrated_speciesCELL_MODEL(Cell_parameters *p)
{
    if (strcmp(p->Celltype, "test") == 0)
    {
       // p->Gto              *= 2; // Relative to native model settings!!!
    }
}

// ISO
void set_modulation_ISO_native_speciesCELL_MODEL(Cell_parameters *p)
{
    // ================================================================================\\|
    // Can modify parameters (e.g. conductance) directly, or scaling factors.
    // (conductance denoted "g", scale factor "G"; current = f(g*G))
    // In general, use the scale factors as these are multiplicative/additive throughout
    // the code. 
    // Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // But overwrites any previous settings absolutely
    // ===========================================
    // "p->ISO" is concentration of ISO; "p->ISO_model" is the ISO model to be applied
    // ===========================================
    // This approach assumes linear transition between ISO = 0 and full effect at ISO = 1. 
    // If non-linear effect is reqired, multiply by another function:
    // e.g. double ISO_scaled = 1/ ( 1 + exp((-4.362*(log(p->ISO) + 3.27)) ) )  and use "ISO_scaled" in all below functions
    // This could be done per target for different dose dependency
    // ===========================================
    // NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
    // ================================================================================//|
    
    // NOTE: relevant scale factors/shifts for each current are found in the current's function
    // If adding global style here, ensure all relevant models have the scale-factor/shift
    // and use it the same way!

    // testing exmaple illustration of model-specific implementation
    if (strcmp(p->ISO_model, "default") == 0)
    {
        //p->Gto				*= (1.0 + p->ISO*(3.0-1.0));	// At maximal effect, x3, 
		//p->IKr_va_ss_kscale	*= (1.0 + p->ISO*(0.5-1.0));  	// At maximal effect, x0.5
		//p->ICaL_vi_ss_shift	+= p->ISO*5;					// At maximal effect, + 5
		//p->Ito_va_tau_scale *= (1.0 + p->ISO*(2.3-1.0));	// At maximal effect x2.3
		//p->Gup				*= (1.0 + p->ISO*(2.5-1.0));  	// At maximal effect, x2.5
		// For ICaL and Jrel, use  GLTCC_kva1_va2 and GRyR_kCO rather than GCaL and Grel, 
		// as this scales the activity not the expression.
		// In native models, the variables do exactly the same thing as GCaL Grel; in integrated models, they
		// scale transition rates in Markov models rather than NLTCC and NRyR - ISO definitely affects activity, not expression
	}
	// Add new ISO_models here: else if (strcmp(p->ISO_model, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid ISO model for the speciesCELL_MODEL models. Please check Model.c and Model_speciesCELL_MODEL.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_speciesCELL_MODEL(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // But overwrites any previous settings absolutely
	// ===========================================
	// "p->Agent_prop" is proportion of maximal Agent effect; "p->Agent" is the Agent to be applied
	// ===========================================
	// This approach assumes linear transition between Agent_prop = 0 and full effect at Agent_prop = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double Agent_scaled = 1/ ( 1 + exp((-4.362*(log(p->Agent_prop) + 3.27)) ) )  and use "Agent_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ===========================================
	// NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
	// ================================================================================//|
    
    // NOTE: relevant scale factors/shifts for each current are found in the current's function
    // If adding global style here, ensure all relevant models have the scale-factor/shift
    // and use it the same way!

	if (strcmp(p->Agent, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	//else if (strcmp(p->Agent, "test") == 0) 
	//{
	//	p->Gto              *= (1.0 + p->Agent_prop*(3.0-1.0));    // At maximal effect, x3, 
	//    p->IKr_va_ss_kscale *= (1.0 + p->Agent_prop*(0.5-1.0));    // At maximal effect, x0.5
	//    p->ICaL_vi_ss_shift += p->Agent_prop*5;                    // At maximal effect, + 5
	//    p->Ito_va_tau_scale *= (1.0 + p->Agent_prop*(2.3-1.0));    // At maximal effect x2.3
	//    p->Gup              *= (1.0 + p->Agent_prop*(2.5-1.0));    // At maximal effect, x2.5
	// For ICaL and Jrel, use  GLTCC_kva1_va2 and GRyR_kCO rather than GCaL and Grel, 
	// as this scales the activity not the expression.
	// In native models, the variables do exactly the same thing as GCaL Grel; in integrated models, they
	// scale transition rates in Markov models rather than NLTCC and NRyR
	//}
	// Add new pharmacological agent here: else if (strcmp(p->Agent, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the speciesCELL_MODEL models. Please check Model.c and Model_speciesCELL_MODEL.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_speciesCELL_MODEL(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// But overwrites any previous settings absolutely
	// ===========================================
	// "p->Remodelling_prop" is proportion of maximal Remodelling effect; "p->Remodelling" is the Remodelling to be applied
	// ===========================================
	// This approach assumes linear transition between Remodelling_prop = 0 and full effect at Remodelling_prop = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double Remodelling_scaled = 1/ ( 1 + exp((-4.362*(log(p->Remodelling_prop) + 3.27)) ) )  and use "Remodelling_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ===========================================
	// NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
	// ================================================================================//|

	// NOTE: relevant scale factors/shifts for each current are found in the current's function
	// If adding global style here, ensure all relevant models have the scale-factor/shift
	// and use it the same way!

	if (strcmp(p->Remodelling, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	//else if (strcmp(p->Remodelling, "test") == 0)
	//{
	//	p->Gto              *= (1.0 + p->Remodelling_prop*(3.0-1.0));    // At maximal effect, x3, 
	//	p->IKr_va_ss_kscale *= (1.0 + p->Remodelling_prop*(0.5-1.0));    // At maximal effect, x0.5
	//	p->ICaL_vi_ss_shift += p->Remodelling_prop*5;                    // At maximal effect, + 5
	//	p->Ito_va_tau_scale *= (1.0 + p->Remodelling_prop*(2.3-1.0));    // At maximal effect x2.3
	//	p->Gup              *= (1.0 + p->Remodelling_prop*(2.5-1.0));    // At maximal effect, x2.5
	//	p->Ko				= 4; // explicitly set extracellular potassium or similar (e.g. ischemia)
	//    p->GCaL, p->Grel *= x ; // use GCaL and Grel here, as remodelling is about expression not activity
	//}
	// Add new Remodelling here: else if (strcmp(p->Remodelling, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Remodelling model for the speciesCELL_MODEL models. Please check Model.c and Model_speciesCELL_MODEL.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_speciesCELL_MODEL(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// But overwrites any previous settings absolutely
	// ===========================================
	// NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
	// ================================================================================//|

	// NOTE: relevant scale factors/shifts for each current are found in the current's function
	// If adding global style here, ensure all relevant models have the scale-factor/shift
	// and use it the same way!

	if (strcmp(p->Mutation, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	//else if (strcmp(p->Mutation, "test") == 0)
	//{
	//	p->Gto              *= 2;       // Scale factor = Multiplies previous settings
	//	p->gKur             = 0.003;    // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
	//	p->IKr_va_ss_kscale *= 1.25;    // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
	//	p->ICaL_vi_ss_shift += 5;       // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
	//	p->IKs_va_tau_scale *= 0.75;    // Multiplies time constant of voltage activation. Multiplies previous settings.
	//	p->Gup              *= 2;       // Scales intracellular upatke rate. Multiplies previous settings.
	//    p->GCaL, p->Grel *= x ; // use GCaL and Grel here, as mutation is about expression not activity
	//}
	// Add new mutation here: else if (strcmp(p->mutation, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Mutation model for the speciesCELL_MODEL models. Please check Model.c and Model_speciesCELL_MODEL.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}
// ACh ===================================\\|
void set_modulation_ACh_speciesCELL_MODEL(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // But overwrites any previous settings absolutely
	// ===========================================
	// "p->ACh" is concentration of ACh; "p->ACh_model" is the ACh model to be applied
	// ===========================================
	// This approach assumes linear transition between ACh = 0 and full effect at ACh = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double ACh_scaled = 1/ ( 1 + exp((-4.362*(log(p->ACh) + 3.27)) ) )  and use "ACh_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ===========================================
	// NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
	// ================================================================================//|
    
    // NOTE: relevant scale factors/shifts for each current are found in the current's function
    // If adding global style here, ensure all relevant models have the scale-factor/shift
    // and use it the same way!

	// testing exmaple illustration of model-specific implementation

	// NOTE: if the model does not have IKACh (or If) in, then modifying them won't make a difference!
	if (strcmp(p->ACh_model, "default") == 0)
	{
		// Put actual things in here
	}          
	//else if (strcmp(p->ACh_model, "test") == 0)  // just an example of how to implement: NOT an actual implementation
	//{
	//	p->gKACh    		= p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) );
	//	// OR
	//	p->GKACh			*= (1.0 + p->ACh*(2.5-1.0));
	//	p->If_va_ss_shift 	+= -7.2*p->ACh;
	//	p->Gto              *= (1.0 + p->ACh*(3.0-1.0));    // At maximal effect, x3,
	// For ICaL and Jrel, use  GLTCC_kva1_va2 and GRyR_kCO rather than GCaL and Grel, 
	// as this scales the activity not the expression.
	// In native models, the variables do exactly the same thing as GCaL Grel; in integrated models, they
	// scale transition rates in Markov models rather than NLTCC and NRyR

	//}
	else
	{
		printf("ERROR: \"%s\" is not a valid ACh for the speciesCELL_MODEL model. Please check Model_speciesCELL_MODEL.cpp for options\n\n", p->ACh_model);
		exit(1);
	}
}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

