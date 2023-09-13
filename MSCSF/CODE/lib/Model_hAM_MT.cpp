// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of the human atrial cell model  //
// by Maleckar-Trayanova:  M.M. Maleckar et al. Am. J. ====  //
// Physiol Heart Circ Physiol 2009;297(4):H1398-410 =======  //
// In turn based on: Nygren et al Circ Res 1998;82(1):63-81  //
// Coded from the publication, 2017-2018 ==================  //
// In-code identifier: hAM_MT =============================  //
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

#include "Model.h"
#include "Structs.h"

// Parameters and specific settings =============================================================\\|
// Set model dependent parameters 
void update_parameters_native_hAM_MT(Cell_parameters *p)
{
	// Updated, MT parameters
	p->gto              = 0.1652;		// s/mF
	p->gKur             = 0.045;		// s/mF

	p->stimduration     = 5.0;      // ms
	p->stimmag          = -13.5;    // pA/pF

	// Default celltype
	p->Celltype			= "RA";
	p->ISO_model		= "Col";
}

// Initial conditions
void initial_conditions_native_hAM_MT(State_variables *s, Cell_parameters p)
{
	s->Vm               = -74;
	s->INa_va           = 0.00291;
	s->INa_vi_1         = 0.9791;
	s->INa_vi_2         = 0.9869;
	s->Ito_va           = 0.07;
	s->Ito_vi           = 0.99;
	s->Ito_vi_s         = 0.361531;
	s->ICaL_va          = 0.0;
	s->ICaL_vi          = 1.0;
	s->ICaL_vi_s        = 1.0;
	s->ICaL_ci          = 1.0;
	s->IKur_va          = 0.0;
	s->IKur_vi          = 1.0;
	s->IKr_va           = 0.0;
	s->IKs_va           = 0.0;

	// Ca handling
	s->CajSR            = 0.6646;       // mM
	s->CanSR            = 0.6646;       // mM
	s->RyRo             = 0.0;
	s->RyRr             = 0.4;
	s->CaCal            = 0.0275;
	s->Catrop           = 0.0133;
	s->Camg             = 0.1961;
	s->Mgmg             = 0.7094;
	s->CaCalse          = 0.44369;

	// Assign state from param (even if constant)
	s->Nai              = p.Nai;
	s->Nao              = p.Nao;
	s->Ki               = p.Ki;
	s->Ko               = p.Ko;
	s->Cai              = p.Cai;
	s->Cai_j            = p.Cai;
	s->Cao              = p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Heterogeneity and modulation =================================================================\\|
// Parent function
void set_het_mod_hAM_MT(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hAM_MT(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hAM_MT(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hAM_MT(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hAM_MT(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hAM_MT(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hAM_MT(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE:: Most het has been set in Model.c as common to multiple human atrial cell models. Please see that function first.
	// set_celltype_hAM() in lib/Model.c

	if (strcmp(p->Celltype, "RA") == 0);    // Do nothing for baseline celltype
	else if (strcmp(p->Celltype, "WL_scaling") == 0)    // Scaling of current magntiudes to fit WL_2006 data
	{
		p->gto              = 0.72*0.1652;   // s/mF        
		p->gCaL             = 1.31*0.1238;   // s/mF
		p->GKur             *= 1.25;   // scaling
		p->gK1              = 1.4*0.09*pow(p->Ko/5.4,0.4);
	}
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Celltype, "test") == 0)
	{
		p->Gto              *= 2;       // Scale factor = Multiplies previous settings
		p->gKur             = 0.003;    // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
		p->IKr_va_ss_kscale *= 1.25;    // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->ICaL_vi_ss_shift += 5;       // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->IKs_va_tau_scale *= 0.75;    // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= 2;       // Scales intracellular upatke rate. Multiplies previous settings.
	}
	// Add new celltypes here: else if (strcmp(p->Celltype, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Celltype for the hAM_MT models. Please check Model.c and Model_hAM_MT.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

// ISO
void set_modulation_ISO_native_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// "ISO" is 0-1 representing 0 to saturating solution i.e. describes the proportion of full ISO scaling to apply
	// This is thus a linear and homogeneous application. If desired, rescale this according to a non-linear
	// function to apply conc dependency etc up to saturating solution. If so, rescale here. 
	// "ISO_model" may be used also to implement different implementations of ISO. If not, will be set to "Default"

	// NOTE: ISO has been set in Model.c set_ISO_hAM() as commong to many human atrial models
	// Place only a model-specific ISO here

	// testing exmaple illustration of model-specific implementation
	if (strcmp(p->ISO_model, "test") == 0)
	{
		p->GKr              *= (1.0 + p->ISO*(2.0 - 1.0)); // x2
		p->GLTCC_kva1_va2   *= (1.0 + p->ISO*(2.0 - 1.0)); // x2
		p->Gup              *= (1.0 + p->ISO*(1.75- 1.0)); // x1.75
	}

	// Add new ISO_models here: else if (strcmp(p->ISO_model, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid ISO model for the hAM_MT models. Please check Model.c and Model_hAM_MT.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// For Agents: "Agent_prop" is 0-1 which describes the proportion of full agent scaling to apply
	// This is thus a linear and homogeneous application. If desired, rescale this according to a non-linear
	// function to apply conc dependency etc up to saturating solution. If so, rescale here. 

	// Note: Many Agents will be common to multiple or all models, and thus should be defined in
	// lib/Model.c. Only model-specific drugs should be here

	if (strcmp(p->Agent, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Agent, "test") == 0)
	{
		p->Gto              *= (1.0 + p->Agent_prop*(2.0  -1));     // Scale factor = Multiplies previous settings  | 2 is scale factor, 2-1 is additive scale factor
		p->IKr_va_ss_kscale *= (1.0 + p->Agent_prop*(1.25 -1));     // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->IKs_va_tau_scale *= (1.0 + p->Agent_prop*(0.75 -1));     // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= (1.0 + p->Agent_prop*(2.0  -1));     // Scales intracellular upatke rate. Multiplies previous settings.
		p->ICaL_vi_ss_shift += 5*p->Agent_prop;                     // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->gKur             = 0.003;                                // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
	}

	// Add new pharmacological agent here: else if (strcmp(p->Agent, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hAM_MT models. Please check Model.c and Model_hAM_MT.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE: AF remodelling from both Grandi 2011 and Colman 2013, common
	// to many human atrial models, is in Model.c set_remodelling_hAM().
	// Add only a model-specific remodelling here

	// For remodelling: "Remodelling_prop" is 0-1 which describes the proportion of full agent scaling to apply
	// This is thus a linear and homogeneous application. If desired, rescale this according to a non-linear
	// function. This may have less/no meaning for remodelling (compared to what it means in Agent/ISO).
	// In general, this can be used for continuous gradients between healthy and remodelling regions

	if (strcmp(p->Remodelling, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Remodelling, "test") == 0)
	{
		p->Gto              *= (1.0 + p->Remodelling_prop*(2.0  -1));     // Scale factor = Multiplies previous settings  | 2 is scale factor, 2-1 is additive scale factor
		p->IKr_va_ss_kscale *= (1.0 + p->Remodelling_prop*(1.25 -1));     // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->IKs_va_tau_scale *= (1.0 + p->Remodelling_prop*(0.75 -1));     // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= (1.0 + p->Remodelling_prop*(2.0  -1));     // Scales intracellular upatke rate. Multiplies previous settings.
		p->ICaL_vi_ss_shift += 5*p->Remodelling_prop;                     // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->gKur             = 0.003;                                // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods.
	}

	// Add new Remodelling here: else if (strcmp(p->Remodelling, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hAM_MT models. Please check Model.c and Model_hAM_MT.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE: Many mutations may be common to multiple models and should be in Model.c
	// Add only a model-specific remodelling here

	if (strcmp(p->Mutation, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Mutation, "test") == 0)
	{
		p->Gto              *= 2;       // Scale factor = Multiplies previous settings
		p->gKur             = 0.003;    // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
		p->IKr_va_ss_kscale *= 1.25;    // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->ICaL_vi_ss_shift += 5;       // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->IKs_va_tau_scale *= 0.75;    // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= 2;       // Scales intracellular upatke rate. Multiplies previous settings.
	}
	// Add new mutation here: else if (strcmp(p->mutation, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Mutation model for the hAM_MT models. Please check Model.c and Model_hAM_MT.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}

// ACh ===================================\\|
void set_modulation_ACh_hAM_MT(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// Note: setting X_set_ref to 1 and then 0 if not set is important for error checking
	// and for being able to jump into global/common AND specific functions properly

	// NOTE: if the model does not have IKACh (or If) in, then modifying them won't make a difference!
	/*if (strcmp(p->ACh_model, "default") == 0)
	  {
	// Put actual things in here
	}          
	else*/ if (strcmp(p->ACh_model, "test") == 0)  // just an example of how to implement: NOT an actual implementation
	{
		p->gKACh    = p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) );
		p->If_va_ss_shift += -7.2*p->ACh;
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid ACh for the hAM_MT model. Please check Model_hAM_MT.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|

// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
void compute_model_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hAM_MT_native(p, var, Vm, s->Cai, s->Ko);
	update_gating_variables_hAM_MT_native(p, var, s, Vm, dt);
	compute_Itot_hAM_MT_native(p, var, s, Vm);
	comp_homeostasis_hAM_NG(p, var, s, Vm, dt);	// lib/Model_hAM_NG.cpp
}

void set_gate_rates_hAM_MT_native(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Ko)
{
	set_INa_hAM_NG_rates(p, var, Vm);
	set_Ito_hAM_MT_rates(p, var, Vm);
	set_ICaL_hAM_NG_rates(p, var, Vm, Cai);
	set_IKur_hAM_MT_rates(p, var, Vm);
	set_IKs_hAM_NG_rates(p, var, Vm);
	set_IKr_hAM_NG_rates(p, var, Vm);
	set_IK1_hAM_NG_variables(p, var, Vm, Ko);
}

void update_gating_variables_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_INa_LR(p, var, s, Vm, dt);         // lib/Model.c
	update_gates_IKs_hAM_NG(p, var, s, Vm, dt);
	update_gates_IKr_hAM_NG(p, var, s, Vm, dt);
	update_gates_Ito_hAM_MT(p, var, s, Vm, dt);
	update_gates_ICaL_hAM_NG(p, var, s, Vm, dt);
	update_gates_IKur_hAM_MT(p, var, s, Vm, dt);
}

void compute_Itot_hAM_MT_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_INa_hAM_NG(p, var, s, Vm);
	compute_Ito_hAM_MT(p, var, s, Vm);
	compute_ICaL_hAM_NG(p, var, s, Vm);
	compute_IKur_hAM_MT(p, var, s, Vm);
	compute_IKr_hAM_NG(p, var, s, Vm);
	compute_IKs_hAM_NG(p, var, s, Vm);
	compute_IK1_hAM_NG(p, var, s, Vm);
	compute_INCX_hAM_NG(p, var, s, Vm);
	compute_INaK_hAM_NG(p, var, s, Vm);
	compute_ICaP_hAM_NG(p, var, s, Vm);
	compute_INab_hAM_NG(p, var, s, Vm);
	compute_ICab_hAM_NG(p, var, s, Vm);	

	if (strcmp(p.environment, "isolated") == 0)
	{
		compute_IK1_hAM_WL_isolated(p, var, s, Vm);
		var->IKr = var->IKs = 0;
	}

	var->Itot   = var->INa + var->Ito + var->IK1 + var->ICaL + var->IKur + var->INCX + var->INaK + var->ICaP + var->INab + var->ICab + var->IKr + var->IKs;
	if (strcmp(p.environment, "isolated") == 0) var->Itot +=  p.AIhyp;
}
// End Compute model functions ==================================================================//|

// Ito ======================================================================\\|
void set_Ito_hAM_MT_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation   
	var->Ito_va_ss          = sigmoid(Vm_ac_ss, 1.0, -11.0*p.Ito_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_va_tau         = 1000*(0.0035 * exp(-(Vm_ac_tau/30.0)*(Vm_ac_tau/30.0)) + 0.0015);
	var->Ito_va_tau         *= p.Ito_va_tau_scale;

	// Voltage inactivation
	var->Ito_vi_ss          = sigmoid(Vm_inac_ss, -40.5, 11.5*p.Ito_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_vi_tau         = 1000*(0.025635 * exp (-((Vm_inac_tau+52.45)/15.8827)*((Vm_inac_tau+52.45)/15.8827)) + 0.01414);
	var->Ito_vi_tau         *= p.Ito_vi_tau_scale;
}

void update_gates_Ito_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi              = rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
}

void compute_Ito_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ito                = p.gto * s->Ito_va * s->Ito_vi * (Vm - var->EK);
	var->Ito                *= p.Gto;
}
// End Ito ==================================================================//|

// IKur =====================================================================\\|
void set_IKur_hAM_MT_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation   
	var->IKur_va_ss         = sigmoid(Vm_ac_ss, -6, -8.6*p.IKur_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKur_va_tau        = (0.009/(1.0 + exp((Vm_ac_tau - (-5.0))/12.0)) + 0.0005)*1000; // ms
	var->IKur_va_tau         *= p.IKur_va_tau_scale;

	// Voltage inactivation
	var->IKur_vi_ss         = sigmoid(Vm_inac_ss, -7.5, 10*p.IKur_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKur_vi_tau         =  (0.59/(1 + exp((Vm_inac_tau - (-60.0))/10.0)) + 3.05)*1000; // ms
	var->IKur_vi_tau         *= p.IKur_vi_tau_scale;
}

void update_gates_IKur_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKur_va              = rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);
	s->IKur_vi              = rush_larsen(s->IKur_vi, var->IKur_vi_ss, var->IKur_vi_tau, dt);
}

void compute_IKur_hAM_MT(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKur               =  p.gKur * s->IKur_va * s->IKur_vi * (Vm - var->EK);
	var->IKur               *= p.GKur;
}
// End IKur =================================================================//|

