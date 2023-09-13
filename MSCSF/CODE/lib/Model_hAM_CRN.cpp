// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of the human atrial cell model  //
// by Courtemanche, Ramirez and Nattel:  M. Courtemanche, =  //
// R.J. Ramirez and S. Nattel. Am J Physiol 1998; =========  //
// 275(1 Pt 2):H301-21. ===================================  // 
// In-code identifier: hAM_CRN ============================  //
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
void set_parameters_native_hAM_CRN(Cell_parameters *p)
{
	// Capacitance and cell structure
	p->Cm				= 100;		// pF
	p->Vcell			= 20100.0;	// um^3
	p->Vcyto			= 0.68 		* p->Vcell;	// um^3
	p->VjSR				= 0.0048 	* p->Vcell; // um^3
	p->VnSR				= 0.0552	* p->Vcell; // um^3

	// Concentrations (constant OR initial condition)
	p->Nai				= 11.2;		// mM
	p->Nao				= 140;		// mM
	p->Ki				= 139;		// mM
	p->Ko				= 5.4;		// mM
	p->Cai				= 0.000102; // mM 
	p->Cao				= 1.8;		// mM

	// Current parameters
	p->gNa				= 7.8; 		// s/mF
	p->gto				= 0.1652;	// s/mF
	p->gCaL				= 0.1238; 	// s/mF
	p->gKur				= 0.0115; 	// s/mF 
	p->gKr				= 0.0294*sqrt(p->Ko/5.4); 	// s/mF
	p->gKs				= 0.129;	// s/mF
	p->gK1				= 0.09*pow(p->Ko/5.4,0.4); 	// s/mF

	p->INCX_bar			= 1600;		// pA/pF
	p->INCX_kNao		= 87.5;		// mM
	p->INCX_kCao		= 1.38;		// mM
	p->INCX_k			= 0.1;
	p->INCX_gamma		= 0.35;

	p->INaK_bar			= 0.59933874; 	// pA/pF
	p->INaK_kK			= 1.5;			// mM
	p->INaK_kNa			= 10.0;			// mM

	p->ICaP_bar			= 0.275;		// pA/pF
	p->ICaP_kCa			= 0.5e-3; 		// mM
	p->gNab				= 0.000674; // s/mF
	p->gCab				= 0.001131; // s/mF

	p->ICaL_ci_tau      = 2.0;		// ms

	// Ca2+ handling
	p->cmdnbar			= 0.050;	// mM
	p->trpnbar			= 0.070;	// mM
	p->csqnbar			= 10.0;		// mM
	p->cmdn_k			= 0.00238;	// mM
	p->trpn_k			= 0.0005;	// mM
	p->csqn_k			= 0.8;		// mM

	p->J_rel_max		= 30.0;		// mM/ms
	p->J_SERCA_max		= 0.005;	// mM/ms
	p->J_SERCA_kCa		= 0.00092;	// mM
	p->J_leak_max		= 0.005;	// mM/ms
	p->J_leak_kCaSR		= 15.0;		// mM
	p->J_jsr_nsr_tau	= 180.0;	// ms

	// Stimulus parameters
	//p->stimduration 	= 0.5;   	// ms
	//p->stimmag      	= -200;    	// pA/pF
	//p->stimduration 	= 2.0;   	// ms
	//p->stimmag      	= -20;    	// pA/pF
	p->stimduration 	= 5.0;      // ms
	p->stimmag      	= -13.5;    // pA/pF

	// Default celltype
	p->Celltype			= "RA";
	p->ISO_model		= "Col";

	// Ihyp
	p->AIhyp				= 0.63; 	// pA/pF
}

// Initial conditions
void initial_conditions_native_hAM_CRN(State_variables *s, Cell_parameters p)
{
	s->Vm      			= -82;
	s->INa_va  			= 0.00291;
	s->INa_vi_1			= 0.9791;
	s->INa_vi_2			= 0.9869;
	s->Ito_va  			= 0.07;
	s->Ito_vi  			= 0.99;
	s->Ito_vi_s			= 0.361531;
	s->ICaL_va			= 0.0;
	s->ICaL_vi			= 1.0;
	s->ICaL_vi_s		= 1.0;
	s->ICaL_ci			= 1.0;
	s->IKur_va			= 0.0;
	s->IKur_vi			= 1.0;
	s->IKr_va			= 0.0;
	s->IKs_va			= 0.0;

	// Ca handling
	s->CajSR			= 1.49;		// mM
	s->CanSR			= 1.49;		// mM
	s->cmdn				= 0.00205; 	// mM
	s->trpn				= 0.0118;	// mM
	s->csqn				= 6.432;	// mM
	s->RyRo				= 0.0;
	s->RyRr				= 1.0;
	s->RyRi				= 0.9992;

	// Assign state from param (even if constant)
	s->Nai				= p.Nai;
	s->Nao				= p.Nao;
	s->Ki				= p.Ki;
	s->Ko				= p.Ko;
	s->Cai				= p.Cai;
	s->Cao				= p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Heterogeneity and modulation =================================================================\\|
// Parent function
void set_het_mod_hAM_CRN(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hAM_CRN(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hAM_CRN(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hAM_CRN(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hAM_CRN(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hAM_CRN(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hAM_CRN(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hAM_CRN(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE:: Most het has been set in Model.c as common to multiple human atrial cell models. Please see that function first.
	// set_celltype_hAM() in lib/Model.c

	if (strcmp(p->Celltype, "RA") == 0); 	// Do nothing for baseline celltype
	else if (strcmp(p->Celltype, "WL_scaling") == 0)	// Scaling of current magntiudes to fit WL_2006 data
	{
		p->gto              = 0.72*0.1652;   // s/mF		
		p->gCaL             = 1.31*0.1238;   // s/mF
		p->GKur             *= 1.25;   // scaling
		p->gK1              = 1.4*0.09*pow(p->Ko/5.4,0.4);	
	}
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Celltype, "test") == 0) 
	{
		p->Gto				*= 2; 		// Scale factor = Multiplies previous settings
		p->gKur				= 0.003;	// Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
		p->IKr_va_ss_kscale	*= 1.25;	// Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->ICaL_vi_ss_shift	+= 5;		// Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->IKs_va_tau_scale	*= 0.75;	// Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup				*= 2;		// Scales intracellular upatke rate. Multiplies previous settings.
	}
	// Add new celltypes here: else if (strcmp(p->Celltype, "X") == 0) {   }
	else 
	{
		printf("ERROR: \"%s\" is not a valid Celltype for the hAM_CRN models. Please check Model.c and Model_hAM_CRN.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

// ISO
void set_modulation_ISO_native_hAM_CRN(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ISO model for the hAM_CRN models. Please check Model.c and Model_hAM_CRN.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hAM_CRN(Cell_parameters *p)
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
		p->gKur             = 0.003;  
	}
	// Add new pharmacological agent here: else if (strcmp(p->Agent, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hAM_CRN models. Please check Model.c and Model_hAM_CRN.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hAM_CRN(Cell_parameters *p)
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

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

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
    else if (strcmp(p->Remodelling, "METS") == 0)
            {
                p->GCaL *=2;
            }

	// Add new Remodelling here: else if (strcmp(p->Remodelling, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hAM_CRN models. Please check Model.c and Model_hAM_CRN.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hAM_CRN(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Mutation model for the hAM_CRN models. Please check Model.c and Model_hAM_CRN.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}
// ACh ===================================\\|
void set_modulation_ACh_hAM_CRN(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ACh for the hAM_CRN model. Please check Model_hAM_CRN.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
void compute_model_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hAM_CRN_native(p, var, Vm, s->Cai);
	update_gating_variables_hAM_CRN_native(p, var, s, Vm, dt);
	compute_Itot_hAM_CRN_native(p, var, s, Vm);
	comp_homeostasis_hAM_CRN(p, var, s, Vm, dt);
}

void set_gate_rates_hAM_CRN_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	set_INa_LR_rates(p, var, Vm);					// lib/Model.c
	set_Ito_hAM_CRN_rates(p, var, Vm);
	set_IKur_hAM_CRN_rates(p, var, Vm);
	set_IKs_hAM_CRN_rates(p, var, Vm);
	set_IKr_hAM_CRN_rates(p, var, Vm);
	set_IK1_hAM_CRN_variables(p, var, Vm);
	set_ICaL_hAM_CRN_rates(p, var, Vm, Cai);
}

void update_gating_variables_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_INa_LR(p, var, s, Vm, dt); 		// lib/Model.c
	update_gates_IKs_hAM_CRN(p, var, s, Vm, dt);
	update_gates_IKr_hAM_CRN(p, var, s, Vm, dt);
	update_gates_ICaL_hAM_CRN(p, var, s, Vm, dt);
	update_gates_Ito_hAM_CRN(p, var, s, Vm, dt);
	update_gates_IKur_hAM_CRN(p, var, s, Vm, dt);
}

void compute_Itot_hAM_CRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_INa_LR(p, var, s, Vm);					// lib/Model.c
	compute_Ito_hAM_CRN(p, var, s, Vm);
	compute_ICaL_hAM_CRN(p, var, s, Vm);
	compute_IKur_hAM_CRN(p, var, s, Vm);
	compute_IKr_hAM_CRN(p, var, s, Vm);
	compute_IKs_hAM_CRN(p, var, s, Vm);
	compute_IK1_hAM_CRN(p, var, s, Vm);
	compute_INCX_hAM_CRN(p, var, s, Vm);
	compute_INaK_hAM_CRN(p, var, s, Vm);
	compute_ICaP_hAM_CRN(p, var, s, Vm);
	compute_INab_hAM_CRN(p, var, s, Vm);
	compute_ICab_hAM_CRN(p, var, s, Vm);

	// Overwrites for isolated
	if (strcmp(p.environment, "isolated") == 0)
	{   
		compute_IK1_hAM_WL_isolated(p, var, s, Vm);
		var->IKr = var->IKs = 0;
	}

	var->Itot   = var->INa + var->Ito + var->IK1 + var->ICaL + var->IKur + var->INCX + var->INaK + var->ICaP + var->INab + var->ICab + var->IKr + var->IKs;
	if (strcmp(p.environment, "isolated") == 0)	var->Itot +=  p.AIhyp;
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
// Identical to that of LR model, found in lib/Model.c
// End INa ==================================================================//|

// Ito ======================================================================\\|
void set_Ito_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;  // Voltage modified by shift applied to activation time constant (and rates as they only define tau)
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation	
	var->Ito_va_ss			= sigmoid(Vm_ac_ss, -20.47, -17.54*p.Ito_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_va_al			= 0.65/(exp(-(Vm_ac_tau+10)/8.5)+exp(-(Vm_ac_tau-30)/59));
	var->Ito_va_bet			= 0.65/(2.5+exp((Vm_ac_tau+82)/17));
	var->Ito_va_tau			= 1.0/(3*(var->Ito_va_al + var->Ito_va_bet));
	var->Ito_va_tau			*= p.Ito_va_tau_scale;

	// Voltage inactivation
	var->Ito_vi_ss			= sigmoid(Vm_inac_ss, -43.1, 5.3*p.Ito_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_vi_al			= 1.0/(18.53+exp((Vm_inac_tau+113.7)/10.95));
	var->Ito_vi_bet			= 1.0/(35.56+exp(-(Vm_inac_tau+1.26)/7.44));
	var->Ito_vi_tau			= 1.0/(3*(var->Ito_vi_al +  var->Ito_vi_bet)); 
	var->Ito_vi_tau			*= p.Ito_vi_tau_scale;
}

void update_gates_Ito_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi              = rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
}

void compute_Ito_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ito				= p.gto * pow(s->Ito_va, 3) * s->Ito_vi * (Vm - var->EK);
	var->Ito				*= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
void set_ICaL_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_CRN_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_CRN_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// calcium inactivation
	var->ICaL_ci_ss			= 1/(1+Cai/0.00035);
	var->ICaL_ci_tau		= p.ICaL_ci_tau; // ms
}

void set_ICaL_hAM_CRN_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss                  = sigmoid(Vm_ss, -10, -8*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	if (fabs(Vm_tau+10) < 1.0e-10)
	{
		*va_tau  			= 4.579/(1.0+exp((Vm_tau+10.0)/-6.24));
	}
	else *va_tau  			= (1-exp((Vm_tau+10)/-6.24))/(0.035*(Vm_tau+10)*(1+exp((Vm_tau+10)/-6.24)));
}

void set_ICaL_hAM_CRN_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss                  = exp(-(Vm_ss+28.0)/6.9)/(1.0+exp(-(Vm_ss+28.0)/6.9));
	*vi_tau                 = 9.0/(0.0197*exp(-pow(0.0337,2)*pow((Vm_tau+10),2))+0.02);
}

void update_gates_ICaL_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_va       		= rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
	s->ICaL_vi    			= rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
	s->ICaL_ci    			= rush_larsen(s->ICaL_ci, var->ICaL_ci_ss, var->ICaL_ci_tau, dt);
}

void compute_ICaL_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaL   			= p.gCaL * s->ICaL_va * s->ICaL_vi * s->ICaL_ci * (Vm - 65);
	var->ICaL 				*= p.GCaL;
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
void set_IKur_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation   
	var->IKur_va_ss         = sigmoid(Vm_ac_ss, -30.3, -9.6*p.IKur_va_ss_kscale);
	var->IKur_va_al			= 0.65/(exp(-(Vm_ac_tau+10)/8.5)+exp(-(Vm_ac_tau-30)/59.0));
	var->IKur_va_bet		= 0.65/(2.5+exp((Vm_ac_tau+82)/17.0));
	var->IKur_va_tau		= 1.0/(3*(var->IKur_va_al + var->IKur_va_bet));
	var->IKur_va_tau         *= p.IKur_va_tau_scale;

	// Voltage inactivation
	var->IKur_vi_ss    		= sigmoid(Vm_inac_ss, 99.5, 27.48*p.IKur_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKur_vi_al			= 1/(21+exp(-(Vm_inac_tau-185)/28));
	var->IKur_vi_bet		= exp((Vm_inac_tau-158)/16); 
	var->IKur_vi_tau		= 1.0/(3*(var->IKur_vi_al + var->IKur_vi_bet));
	var->IKur_vi_tau        *= p.IKur_vi_tau_scale;

	// Dynamic conductance factor
	var->IKur_dynamic_g		= 0.005+0.05/(1+exp(-(Vm-15)/13));
}

void update_gates_IKur_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKur_va  			= rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);
	s->IKur_vi  			= rush_larsen(s->IKur_vi, var->IKur_vi_ss, var->IKur_vi_tau, dt);
}

void compute_IKur_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKur 				= var->IKur_dynamic_g * pow(s->IKur_va, 3) * s->IKur_vi * (Vm - var->EK);
	var->IKur				*= p.GKur;
}
// End IKur =================================================================//|

// IKr ======================================================================\\|
void set_IKr_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	// Voltage activation   
	var->IKr_va_ss          = sigmoid(Vm_ac_ss, -14.1, -6.5*p.IKr_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKr_va_tau         = 1/(0.0003*(Vm_ac_tau+14.1)/(1-exp(-(Vm_ac_tau+14.1)/5))+0.000073898*(Vm_ac_tau-3.3328)/(exp((Vm_ac_tau-3.3328)/5.1237)-1));
	var->IKr_va_tau         *= p.IKr_va_tau_scale;

	// Time-independant gate
	var->IKr_vi_ti          = sigmoid(Vm_ac_ss, -15, 22.4);
}

void update_gates_IKr_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKr_va              = rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);
}

void compute_IKr_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKr               =  p.gKr * s->IKr_va * var->IKr_vi_ti * (Vm - var->EK);
	var->IKr               *= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_hAM_CRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	// Voltage activation   
	var->IKs_va_ss          = 1.0/pow((1+exp(-(Vm_ac_ss - (19.9))/12.7)),0.5);
	var->IKs_va_tau         = 0.5/(0.00004*(Vm_ac_tau-19.9)/(1-exp(-(Vm_ac_tau-19.9)/17))+0.000035*(Vm_ac_tau-19.9)/(exp((Vm_ac_tau-19.9)/9)-1));
	var->IKs_va_tau         *= p.IKs_va_tau_scale;
}

void update_gates_IKs_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKs_va              = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
}

void compute_IKs_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKs               = p.gKs * s->IKs_va * s->IKs_va * (Vm - var->EK);
	var->IKs               *= p.GKs;
}
// End IKs ==================================================================//|

// IK1 ======================================================================\\|
void set_IK1_hAM_CRN_variables(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_in        = Vm - p.IK1_va_shift;
	var->IK1_va_ti      = (1.0 + exp(0.07*(Vm_in-(-80))));
}

void compute_IK1_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IK1               = p.gK1 * (Vm - var->EK - p.IK1_Erev_shift)/var->IK1_va_ti;
	var->IK1               *= p.GK1;
}
// End IK1 ==================================================================//|

// Ca2+ handling, background and pump currents ==============================\\|
// INCX
void compute_INCX_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	// Create some local variables to make equation easier to read
	double Cai 		 	= s->Cai;
	double Cao			= s->Cao;
	double Nai			= s->Nai;
	double Nao			= s->Nao;

	var->INCX 			= p.INCX_bar *( exp(p.INCX_gamma* p.FoRT*Vm) * Nai*Nai*Nai * Cao - exp((p.INCX_gamma-1)* p.FoRT*Vm) * Nao*Nao*Nao*Cai)		\
						  / ((pow(p.INCX_kNao,3) + Nao*Nao*Nao)*(p.INCX_kCao + Cao) * (1 + p.INCX_k * exp((p.INCX_gamma-1)* p.FoRT*Vm)));
	var->INCX			*= p.GNCX;

}

// INaK
void compute_INaK_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double sigma		= (exp(s->Nao/67.3)-1.0)/7.0;
	double FNaK			= pow(1.0+0.1245*exp(-0.1*p.F*Vm/(p.R*p.T))+0.0365*sigma*exp(-Vm*p.FoRT ), -1.0);
	var->INaK			= p.INaK_bar * FNaK * (1.0/(1.0+pow((p.INaK_kNa /s->Nai),1.5))) * (s->Ko/(s->Ko + p.INaK_kK));
	var->INaK			*= p.GNaK;
}

// ICaP
void compute_ICaP_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaP  			= (p.ICaP_bar * s->Cai) / (s->Cai + p.ICaP_kCa);
	var->ICaP			*= p.GCaP;
}

// INab
void compute_INab_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->INab  			= p.gNab * (Vm - var->ENa);	
	var->INab			*= p.GNab;
}

// ICab
void compute_ICab_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICab  			= p.gCab * (Vm - var->ECa);	
	var->ICab			*= p.GCab;
}
// End Ca2+ handling currents ===============================================//|

// Homeostasis ==============================================================\\|
void comp_homeostasis_hAM_CRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Sodium and Potassium differentials
	var->dNai 				= p.Cm*(-3*var->INaK - 3*var->INCX - var->INab - var->INa)/(p.F*p.Vcyto);
	var->dKi  				= p.Cm*(2*var->INaK - var->IK1 - var->Ito - var->IKur - var->IKr - var->IKs)/(p.F*p.Vcyto);

	// Buffering
	s->cmdn                 = p.cmdnbar * (s->Cai/(s->Cai + p.cmdn_k)); // mM 
	s->trpn                 = p.trpnbar * (s->Cai/(s->Cai + p.trpn_k)); // mM 
	s->csqn                 = p.csqnbar * (s->CajSR/(s->CajSR + p.csqn_k));

	// Jrel
	var->J_rel 				= p.J_rel_max * s->RyRo*s->RyRo * s->RyRr * s->RyRi *(s->CajSR - s->Cai);
	var->J_rel				*= p.Grel;
	var->Cai_Fn 			= 1.0e3*( 1.0e-15*p.VjSR*var->J_rel - 1.0e-15/(2.0*p.F)*(0.5*var->ICaL-0.2*var->INCX)*p.Cm );
	s->RyRo                 = rush_larsen(s->RyRo, 1.0/(1.0+exp(-(var->Cai_Fn-3.4175e-13 )/13.67e-16)) , 8.0, dt);
	s->RyRr                 = rush_larsen(s->RyRr, 1.0-1.0/(1.0+exp(-(var->Cai_Fn-6.835e-14)/13.67e-16)) , 1.91+2.09/(1.0+exp(-(var->Cai_Fn-3.4175e-13)/13.67e-16)), dt);
	s->RyRi                 = rush_larsen(s->RyRi, 1.0-1.0/(1.0+exp(-(Vm -40.0)/17.0)) , 6.0*(1.0-exp(-(Vm -7.9)/5.0))/((1.0+0.3*exp(-(Vm -7.9)/5.0))*(Vm -7.9)), dt);

	// Jup
	var->J_SERCA   			= p.J_SERCA_max * s->Cai/(s->Cai + p.J_SERCA_kCa);
	var->J_leak		   		= p.J_leak_max  * s->CanSR / p.J_leak_kCaSR;
	var->J_SERCA			*= p.Gup;
	var->J_leak				*= p.Gleak;

	// SR compart transfer
	var->J_jsr_nsr			= (s->CanSR - s->CajSR)/p.J_jsr_nsr_tau;

	// Update Ca concentrations
	var->dCai 				= p.Cm*(2.0*var->INCX-(var->ICaP+var->ICaL+var->ICab))/(2.0*p.Vcyto*p.F)+(p.VnSR*(var->J_leak - var->J_SERCA)+var->J_rel*p.VjSR)/p.Vcyto;
	var->dCai 				*= 1.0/(1.0+p.trpnbar*p.trpn_k/pow(s->Cai + p.trpn_k, 2.0) + p.cmdnbar*p.cmdn_k/pow(s->Cai + p.cmdn_k, 2.0));
	s->CanSR   				+= dt* (var->J_SERCA - var->J_jsr_nsr * p.VjSR/p.VnSR - var->J_leak);
	s->CajSR  				+= dt* (var->J_jsr_nsr - var->J_rel)/((1+p.csqnbar*0.8/pow((s->CajSR+0.8),2)));
	s->Cai  				+= dt* var->dCai;

	// Sodium and Potassium concs  || commented out for stability
	//s->Nai                  += dt*var->dNai;
	//s->Ki                   += dt*var->dKi;
}
// End Homeostasis ==========================================================//|
// End Current formulations =====================================================================//|
