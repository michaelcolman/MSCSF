// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of the human atrial cell model  //
// by Nygren-Giles: Nygren et al Circ Res 1998;82(1):63-81=  //
// Coded from the publication, 2017-2018 ==================  //
// In-code identifier: hAM_NG =============================  //
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

// Specific notes: rates are multipled by 1e-3 due to =====  //
// conversion from s to ms ================================  //

#include "Model.h"
#include "Structs.h"

// Parameters and specific settings =============================================================\\|
// Set model dependent parameters 
void set_parameters_native_hAM_NG(Cell_parameters *p)
{
	// Capacitance and cell structure
	p->Cm				= 50;		// pF
	p->Vcyto			= 5.884E-6;	//uL
	p->VjSR				= 4.41E-8;
	p->VnSR				= 3.969E-7;
	p->Vc				= 0.136 * p->Vcyto;
	p->Vjunc			= 0.02	* p->Vcyto;

	// Concentrations (constant OR initial condition)
	p->Nai				= 8.5547;		// mM
	p->Nao				= 130;		// mM
	p->Ki				= 129.435;	// mM
	p->Ko				= 5.3581;	// mM
	p->Cai				= 0.00006729; // mM 
	p->Cao				= 1.8147;	// mM

	// Current parameters
	p->gNa				= 0.000032;	// s/mF
	p->gto				= 0.150000; // s/mF
	p->gCaL				= 0.135000; // s/mF
	p->gKur				= 0.055000; // s/mF
	p->gKr				= 0.010000; // s/mF
	p->gKs				= 0.020000; // s/mF
	p->gK1				= 0.060000; // s/mF

	p->INaK_bar			= 1.416506;	// pA/pF
	p->INCX_bar			= 0.000750;	// pA/((mmol/L)^4.pF)
	p->ICaP_bar			= 0.080000;	// pA/pF
	p->ICaP_kCa			= 0.000200;	// mmol/L
	p->gNab				= 0.001212; // s/mF
	p->gCab				= 0.001574; // s/mF

	// Environment params
	p->T            	= 306.15;      // K
	p->R            	= 8.31441;    // J mol^-1 K^-1
	p->FoRT         	= p->F/(p->R*p->T);

	p->stimduration 	= 5.0;      // ms
	//p->stimmag      	= -8;    	// pA/pF  || original value
	p->stimmag      	= -13.5;    // pA/pF  || for consistency with other hAM models

	// Default celltype
	p->Celltype			= "RA";
	p->ISO_model		= "Col";

	// Ihyp
	p->AIhyp			= 0.63; 	// pA/pF
}

// Initial conditions
void initial_conditions_native_hAM_NG(State_variables *s, Cell_parameters p)
{
	s->Vm      			= -74.2525;
	s->INa_va  			= 0.0032017;
	s->INa_vi_1			= 0.8814;
	s->INa_vi_2			= 0.8742;
	s->Ito_va  			= 0.0010678;
	s->Ito_vi  			= 0.9490;
	s->Ito_vi_s			= 0.361531;
	s->ICaL_va			= 0.000013;
	s->ICaL_vi			= 0.9986;
	s->ICaL_vi_s		= 0.9986;
	s->ICaL_ci			= 1.0;
	s->IKur_va			= 0.00015949;
	s->IKur_vi			= 0.9912;
	s->IKr_va			= 0.0001;
	s->IKs_va			= 0.0048357;

	// Ca handling
	s->CajSR			= 0.6465;		// mM
	s->CanSR			= 0.6646;		// mM
	s->RyRo				= 0.0028;
	s->RyRr				= 0.4284;
	s->CaCal     		= 0.0275;
	s->Catrop      		= 0.0133;
	s->Camg        		= 0.1961;
	s->Mgmg       		= 0.7094;
	s->CaCalse    		= 0.44369;

	// Assign state from param (even if constant)
	s->Nai				= p.Nai;
	s->Nao				= p.Nao;
	s->Ki				= p.Ki;
	s->Ko				= p.Ko;
	s->Cai				= p.Cai;
	s->Cai_j			= p.Cai;
	s->Cao				= p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Heterogeneity and modulation =================================================================\\|
// Parent function
void set_het_mod_hAM_NG(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hAM_NG(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hAM_NG(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hAM_NG(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hAM_NG(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hAM_NG(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hAM_NG(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hAM_NG(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	if (strcmp(p->Celltype, "RA") == 0); 	// Do nothing for baseline celltype
	else if (strcmp(p->Celltype, "WL_scaling") == 0)	// Scaling of current magntiudes to fit WL_2006 data
	{
		p->Gto              *= 0.88;
		p->GCaL             *= 1.2;
		p->GKur             *= 1.2;   // scaling
		p->GK1				*= 1.111;
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
		printf("ERROR: \"%s\" is not a valid Celltype for the hAM_NG models. Please check Model.c and Model_hAM_NG.cpp for options\n\n", p->Celltype);
		exit(1);
	}

}

// ISO
void set_modulation_ISO_native_hAM_NG(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ISO model for the hAM_NG models. Please check Model.c and Model_hAM_NG.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hAM_NG(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// Note: Many Agents will be commong to multiple or all models, and thus should be defined in
	// lib/Model.c. Only model-specific drugs should be here

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
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hAM_NG models. Please check Model.c and Model_hAM_NG.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hAM_NG(Cell_parameters *p)
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

	// Add new Remodelling here: else if (strcmp(p->Remodelling, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hAM_NG models. Please check Model.c and Model_hAM_NG.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hAM_NG(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Mutation model for the hAM_NG models. Please check Model.c and Model_hAM_NG.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}

// ACh ===================================\\|
void set_modulation_ACh_hAM_NG(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ACh for the hAM_NG model. Please check Model_hAM_NG.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
void compute_model_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hAM_NG_native(p, var, Vm, s->Cai, s->Ko);
	update_gating_variables_hAM_NG_native(p, var, s, Vm, dt);
	compute_Itot_hAM_NG_native(p, var, s, Vm);
	comp_homeostasis_hAM_NG(p, var, s, Vm, dt);
}

void set_gate_rates_hAM_NG_native(Cell_parameters p, Model_variables *var, double Vm, double Cai, double Ko)
{
	set_INa_hAM_NG_rates(p, var, Vm);
	set_Ito_hAM_NG_rates(p, var, Vm);
	set_ICaL_hAM_NG_rates(p, var, Vm, Cai);
	set_IKur_hAM_NG_rates(p, var, Vm);
	set_IKs_hAM_NG_rates(p, var, Vm);
	set_IKr_hAM_NG_rates(p, var, Vm);
	set_IK1_hAM_NG_variables(p, var, Vm, Ko);
}

void update_gating_variables_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_INa_LR(p, var, s, Vm, dt); 		// lib/Model.c
	update_gates_IKs_hAM_NG(p, var, s, Vm, dt);
	update_gates_IKr_hAM_NG(p, var, s, Vm, dt);
	update_gates_ICaL_hAM_NG(p, var, s, Vm, dt);
	update_gates_Ito_hAM_NG(p, var, s, Vm, dt);
	update_gates_IKur_hAM_NG(p, var, s, Vm, dt);
}

void compute_Itot_hAM_NG_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_INa_hAM_NG(p, var, s, Vm);		
	compute_Ito_hAM_NG(p, var, s, Vm);
	compute_ICaL_hAM_NG(p, var, s, Vm);
	compute_IKur_hAM_NG(p, var, s, Vm);
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
	if (strcmp(p.environment, "isolated") == 0) var->Itot +=  p.AIhyp; // Ihyp
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
void set_INa_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac        = Vm - p.INa_va_shift;
	double Vm_inac      = Vm - p.INa_vi_shift;

	// Set Activation gate alpha and beta
	var->INa_va_ss     			= 1.0 / (1 + exp(-(Vm_ac + 27.12) / 8.21));
	var->INa_va_tau  			= 1e3*((0.000042 * exp(-pow( ((Vm_ac + 25.57) / 28.8 ) , 2))) + 0.000024); // s -> ms gives 1e3 for ALL time-constants (1e-3 for all rates) throughout

	// Set inactivation gates alphas and betas
	var->INa_vi_1_ss   			= 1.0 / (1 + exp((Vm_inac + 63.6) / 5.3));
	var->INa_vi_2_ss  			= var->INa_vi_1_ss;
	var->INa_vi_1_tau  			= 1e3*(0.03 / (1 + exp((Vm_inac + 35.1) / 3.2)) + 0.0003);
	var->INa_vi_2_tau  			= 1e3*(0.12 / (1 + exp((Vm_inac + 35.1) / 3.2)) + 0.003);

	// Scale time constants
	var->INa_va_tau    			*= p.INa_va_tau_scale;
	var->INa_vi_1_tau  			*= p.INa_vi_1_tau_scale;
	var->INa_vi_2_tau     		*= p.INa_vi_2_tau_scale;	
}

void compute_INa_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	if (fabs(Vm) > 1.0E-4)
	{
		var->INa                      = p.gNa * pow(s->INa_va, 3) * (0.9 * s->INa_vi_1 + 0.1 * s->INa_vi_2) * s->Nao * Vm * 1e3*p.F * p.FoRT * (exp((Vm - var->ENa) * p.FoRT) - 1) /
			(exp(Vm * p.FoRT) - 1);
	}
	else
	{
		var->INa                      = p.gNa * pow(s->INa_va, 3) * (0.9 * s->INa_vi_1 + 0.1 * s->INa_vi_2) * s->Nao * 1e3*p.F * (exp((Vm - var->ENa) *p.FoRT) - 1);
	}
	var->INa        *= p.GNa;
}
// End INa ==================================================================//|

// Ito ======================================================================\\|
void set_Ito_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;  // Voltage modified by shift applied to activation time constant (and rates as they only define tau)
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation	
	var->Ito_va_ss			= sigmoid(Vm_ac_ss, 1.0, -11.0*p.Ito_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_va_tau			= ((0.0035 * exp(-pow((Vm / 30.0), 2))) + 0.0015)*1e3;	// 1e3 s->ms
	var->Ito_va_tau			*= p.Ito_va_tau_scale;

	// Voltage inactivation
	var->Ito_vi_ss			= sigmoid(Vm_inac_ss, -40.5, 11.5*p.Ito_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->Ito_vi_tau			= ((0.4812 * exp(-pow((Vm + 52.45) / 14.97, 2))) + 0.01414)*1e3; 
	var->Ito_vi_tau			*= p.Ito_vi_tau_scale;
}

void update_gates_Ito_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi              = rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
}

void compute_Ito_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ito				= p.gto * s->Ito_va * s->Ito_vi * (Vm - var->EK);
	var->Ito				*= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
void set_ICaL_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_NG_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_NG_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// Slow V inac
	var->ICaL_vi_s_tau		= ((0.3323 * exp(-pow(((Vm + 40.0) / 14.2), 2))) + 0.0626)*1e3;
}

void set_ICaL_hAM_NG_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss    			= sigmoid(Vm_ss, -9.0, -5.8*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*va_tau  			= ((0.0027 * exp(-pow(((Vm_tau + 35.0) / 30.0), 2))) + 0.002)*1e3;
}

void set_ICaL_hAM_NG_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss                  = sigmoid(Vm_ss, -27.4, 7.1*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*vi_tau                 = (0.161 * exp(-pow( ( (Vm_tau + 40.0) / 14.4 ) , 2)) + 0.01)*1e3;
}

void update_gates_ICaL_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_va       		= rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
	s->ICaL_vi    			= rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
	s->ICaL_vi_s            = rush_larsen(s->ICaL_vi_s, var->ICaL_vi_ss, var->ICaL_vi_s_tau, dt);
	s->ICaL_ci    			= s->Cai_j / (s->Cai_j + 0.025);
}

void compute_ICaL_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaL   			= p.gCaL * s->ICaL_va * (s->ICaL_ci*s->ICaL_vi + (1 - s->ICaL_ci)*s->ICaL_vi_s) * (Vm - 60);
	var->ICaL 				*= p.GCaL;
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
void set_IKur_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation   
	var->IKur_va_ss         = sigmoid(Vm_ac_ss, -4.3, -8.0*p.IKur_va_ss_kscale);		// Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKur_va_tau		= ((0.009 / (1 + exp((Vm + 5.0) / 12.0))) + 0.0005)*1e3;
	var->IKur_va_tau         *= p.IKur_va_tau_scale;

	// Voltage inactivation
	var->IKur_vi_ss    		= 0.4*sigmoid(Vm_inac_ss, -20.0, 10.0*p.IKur_vi_ss_kscale) + 0.6; // Vm, V1/2, k  0.4/(1+exp((V-V1/2)/k) + 0.6
	var->IKur_vi_tau		= ((0.047 / (1 + exp((Vm + 60.0) / 10.0))) + 0.3)*1e3;
	var->IKur_vi_tau        *= p.IKur_vi_tau_scale;
}

void update_gates_IKur_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKur_va  			= rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);
	s->IKur_vi  			= rush_larsen(s->IKur_vi, var->IKur_vi_ss, var->IKur_vi_tau, dt);
}

void compute_IKur_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKur 				= p.gKur * s->IKur_va * s->IKur_vi * (Vm - var->EK);
	var->IKur				*= p.GKur;
}
// End IKur =================================================================//|

// IKr ("IKf" in original NG model) =========================================\\|
void set_IKr_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	// Voltage activation   
	var->IKr_va_ss          = sigmoid(Vm_ac_ss, -15.0, -6.0*p.IKr_va_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKr_va_tau         = (0.03118 + 0.21718 * exp(-pow((Vm_ac_tau + 20.1376) / 22.1996, 2)))*1e3;
	var->IKr_va_tau         *= p.IKr_va_tau_scale;

	// Time-independant gate
	var->IKr_vi_ti          = sigmoid(Vm_ac_ss, -55, 24.0);
}

void update_gates_IKr_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKr_va              = rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);
}

void compute_IKr_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKr               =  p.gKr * s->IKr_va * var->IKr_vi_ti * (Vm - var->EK);
	var->IKr               *= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_hAM_NG_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	// Voltage activation   
	var->IKs_va_ss          = sigmoid(Vm_ac_ss, 19.9, -12.7*p.IKs_va_ss_kscale);
	var->IKs_va_tau         = (0.7 + 0.4 * exp(-pow((Vm_ac_tau - 20.0) / 20.0, 2)))*1e3;
	var->IKs_va_tau         *= p.IKs_va_tau_scale;
}

void update_gates_IKs_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKs_va              = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
}

void compute_IKs_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKs                =  p.gKs * s->IKs_va * (Vm - var->EK);
	var->IKs                *= p.GKs;
}
// End IKs ==================================================================//|

// IK1 ======================================================================\\|
void set_IK1_hAM_NG_variables(Cell_parameters p, Model_variables *var, double Vm, double Ko)
{
	double Vm_in        = Vm - p.IK1_va_shift;
	var->IK1_va_ti      = pow(Ko, 0.4457)/(1 + exp(1.5 * (Vm_in - var->EK + 3.6) * p.FoRT));
}

void compute_IK1_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IK1               = p.gK1 * (Vm - var->EK - p.IK1_Erev_shift)*var->IK1_va_ti;
	var->IK1               *= p.GK1;
}
// End IK1 ==================================================================//|

// Ca2+ handling, background and pump currents ==============================\\|
// INCX
void compute_INCX_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	// Create some local variables to make equation easier to read
	double Cai 		 	= s->Cai;
	double Cao			= s->Cao;
	double Nai			= s->Nai;
	double Nao			= s->Nao;

	var->INCX			= p.INCX_bar*(pow(Nai, 3) * Cao * exp(0.450 * Vm * p.FoRT) - pow(Nao, 3) * s->Cai * exp(Vm * (0.45 - 1) * p.FoRT)) /				\
						  (1 + 0.0003 * (s->Cai * pow(Nao, 3) + Cao * pow(Nai, 3)));
	var->INCX			*= p.GNCX;

}

// INaK
void compute_INaK_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double FNaK			= (Vm +150)/(Vm + 200);
	var->INaK           = p.INaK_bar * FNaK *  s->Ko / (s->Ko + 1) * (pow(s->Nai, 1.5) / (pow(s->Nai, 1.5) + pow(11.0, 1.5)));
	var->INaK			*= p.GNaK;
}

// ICaP
void compute_ICaP_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaP  			= (p.ICaP_bar * s->Cai) / (s->Cai + p.ICaP_kCa);
	var->ICaP			*= p.GCaP;
}

// INab
void compute_INab_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->INab  			= p.gNab * (Vm - var->ENa);	
	var->INab			*= p.GNab;
}

// ICab
void compute_ICab_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICab  			= p.gCab * (Vm - var->ECa);	
	var->ICab			*= p.GCab;
}
// End Ca2+ handling currents ===============================================//|

// Homeostasis ==============================================================\\|
void comp_homeostasis_hAM_NG(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Sodium and Potassium differentials
	var->dNai   			= 1e-3*p.Cm * (-(3 * var->INaK + 3 * var->INCX + var->INab + var->INa - 0.0336) / (1e3*p.F * p.Vcyto));		
	var->dKi           		= 1e-3*p.Cm * ((2.0 * var->INaK - var->IKr - var->IKs - var->Ito - var->IKur - var->IK1) / (1e3*p.F *p.Vcyto));

	// Jrel
	var->J_rel       	= 200.0 * pow((s->RyRo / (s->RyRo + 0.25)), 2) * (s->CajSR - s->Cai);
	var->J_rel			*= p.Grel;

	double ract, rinact;
	ract                    = 203.8 * pow((s->Cai_j / (s->Cai_j + 0.003)), 4) + 203.8 * pow((s->Cai / (s->Cai + 0.0003)), 4);
	rinact                  = 33.96 + 339.6 * pow((s->Cai / (s->Cai + 0.0003)), 4);

	double dRyRr, dRyRo;
	dRyRr                   = ((0.815 * (1 - s->RyRr - s->RyRo) - ract * s->RyRr));
	dRyRo                   = (ract * s->RyRr - rinact * s->RyRo);
	s->RyRr                 += 1e-3*dt*dRyRr;
	s->RyRo                 += 1e-3*dt*dRyRo;

	// Jup
	var->J_SERCA             = (2.8 * (((s->Cai / 0.0003) - (0.4 * 0.4 * s->CanSR / 0.5)) / (((s->Cai + 0.0003) / 0.0003) + (0.4 * (s->CanSR + 0.5) / 0.5))));
	var->J_SERCA            *= p.Gup;

	// SR compart transfer
	var->J_jsr_nsr          = (s->CanSR - s->CajSR)* 2.0 * 1e3*p.F * p.VjSR / 0.01;

	// Cai factors
	double fab;
	fab = 0.08 * (78400.0 * s->Cai * (1 - s->Catrop) - 392.0 * s->Catrop)
		+ 0.16 * (200000.0 * s->Cai * (1 - s->Camg - s->Mgmg) - 6.6 * s->Camg)
		+ 0.045 * (200000.0 * s->Cai * (1 - s->CaCal) - 476.0 * s->CaCal);

	double idi;
	idi                     = (s->Cai_j - s->Cai) * 2.0 * p.Vjunc * 1e3*p.F / 0.01;


	// Update Ca compartments and buffers
	double dCanSR, dCajSR, dCaCal, dCatrop, dCamg, dMgmg, dCaCalse;  // due to inter-dependence of variables, differentials must be set before updated
	dCanSR                   = ((var->J_SERCA - var->J_jsr_nsr) / (2.0 * 1e3*p.F * p.VnSR));
	dCajSR                  = ((var->J_jsr_nsr - var->J_rel) / (2.0 * 1e3*p.F * p.VjSR) - 31.0 * (480 * s->CajSR * (1 - s->CaCalse) - 400 * s->CaCalse));
	dCaCal                  = ((200000.0 * s->Cai * (1 - s->CaCal) - 476.0 * s->CaCal));
	dCatrop                 = ((78400.0 * s->Cai * (1 - s->Catrop) - 392.0 * s->Catrop));
	dCamg                   = ((200000.0 * s->Cai * (1 - s->Camg - s->Mgmg) - 6.6 * s->Camg));
	dMgmg                   = (2000 * 2.5 * (1 - s->Camg - s->Mgmg) - 666.0 * s->Mgmg);
	dCaCalse                = ((480 * s->CajSR * (1 - s->CaCalse) - 400.0*s->CaCalse));

	s->CanSR                += 1e-3*dt*dCanSR;
	s->CajSR                += 1e-3*dt*dCajSR;
	s->CaCal                += 1e-3*dt*dCaCal;
	s->Catrop               += 1e-3*dt*dCatrop;
	s->Camg                 += 1e-3*dt*dCamg;
	s->Mgmg                 += 1e-3*dt*dMgmg;
	s->CaCalse              += 1e-3*dt*dCaCalse;

	s->Cai                  += 1e-3*dt*((idi + 1e-3*p.Cm*(2.0 * var->INCX - var->ICaP - var->ICab) - var->J_SERCA + var->J_rel) / (2.0 * p.Vcyto * 1e3*p.F) - 1.0 * fab);
	s->Cai_j                += 1e-3*dt*(-(1e-3*p.Cm*var->ICaL + idi) / (2.0 * p.Vjunc * 1e3*p.F));

	// Sodium and Potassium concs
	s->Nai                  += 1e-3*dt*var->dNai;
	s->Ki                   += 1e-3*dt*var->dKi;
	s->Nao				  	+= 1e-3*dt*( (130.0 - s->Nao)/14.3 + 1e-3*p.Cm*(var->INa + var->INab + 3*var->INaK + 3*var->INCX - 0.00168)/(1e3*p.F * p.Vcyto));
	s->Ko					+= 1e-3*dt*( (5.4 - s->Ko)/10 + 1e-3*p.Cm*(-2*var->INaK + var->IKr + var->IKs + var->Ito + var->IK1 + var->IKur)/(1e3*p.F * p.Vcyto));
}
// End Homeostasis ==========================================================//|
// End Current formulations =====================================================================//|
