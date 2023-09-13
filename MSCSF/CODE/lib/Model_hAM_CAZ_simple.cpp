// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of a simplified version of the  // 
// human atrial cell model model associated with:  ========  //
// M.A. Colman et al. 2013 J. Physiol =====================  //
// Simplified model is for integration with spatial Ca2+ ==  //
// handling model. Original full model is not packaged with  //
// this code. =============================================  //
// In-code identifier: hAM_CAZ_s ==========================  //
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
// If your model does not inhereted from a cell model already included, define all parameters:
void set_parameters_native_hAM_CAZ_simple(Cell_parameters *p)
{
	// Capacitance and cell structure
	//p->Cm				= X;		// pF
	//p->Vcell			= Y;	// um^3

	// Concentrations (constant OR initial condition)
	p->Ko                   = 4.5;      // mM
	p->Nao                  = 140;      // mM
	p->Cao                  = 1.8;      // mM
	p->Ki                   = 143.103;  // mM
	p->Cai                  = 0.1;      // mM
	p->Nai                  = 7.95;     // mM


	// Current parameters
	p->gNa                  = 16.0;			// s/mF
	p->gto      			= 0.1652;		
	p->gKur     			= 0.05874;		
	p->gKs      			= 0.12941176;
	p->gKr      			= 1.2*0.029411765;
	p->gK1      			= 0.1;
	p->gNab     			= 0.004;

	// Stimulus parameters (if different to default)
	p->stimduration 	= 5.0;      // ms
	p->stimmag      	= -12.5;    // pA/pF

	// Default celltype
	p->ISO_model        = "Toy_methods_demonstration";
	p->Celltype			= "RA";
}

void update_parameters_integrated_hAM_CAZ_simple(Cell_parameters *p)
{
	// All of these are involved in the Ca handling model
	p->J_SERCA_max			= 0.91*0.3729;       // uM/ms
	p->J_leak_max       	= 0.91*1.41265*1e-5; // ms^-1
	p->INCX_bar         	= 1.15*0.9*0.3726;   // (um^3.uM.ms^-1)
	p->ICab_bar         	= 1.15*0.6*1.8237*1e-5; // (um^3 . uM . ms^-1)
	p->ICaP_bar         	= 1.15*0.9*0.137*1e-2;   // (um^3 . uM . ms^-1)
	p->NLTCC_mean   		= 10;               // Number of LTCCs per dyad 
}

// Initial conditions
void initial_conditions_native_hAM_CAZ_simple(State_variables *s, Cell_parameters p)
{
	// Define ICs for ALL state variables used in the model
	s->Vm      			= -85;
	s->INa_va  			= 0.001231;
	s->INa_vi_1			= 0.992842;
	s->INa_vi_2			= 0.988864;
	s->Ito_va			= 0;
	s->Ito_vi			= 0.5;
	s->IKur_va			= 0;
	s->IKur_vi			= 1;
	s->IKr_va			= 0.000175;
	s->IKs_va			= 0.397626;

	// Ca handling
	//s->CajSR				= X;		// mM
	//s->RyRo				= X;

	// Assign state from param (even if constant) for concentrations
	// Ensure the ICs for all dynamic concenstrations is set as a parameter
	s->Nai              = p.Nai; 
	s->Nao              = p.Nao;
	s->Ki               = p.Ki;
	s->Ko               = p.Ko;
	s->Cai              = p.Cai;
	s->Cao              = p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Heterogeneity and modulation =================================================================\\|
// Parent function
void set_het_mod_hAM_CAZ_simple(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hAM_CAZ_simple(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hAM_CAZ_simple(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hAM_CAZ_simple(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hAM_CAZ_simple(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hAM_CAZ_simple(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hAM_CAZ_simple(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hAM_CAZ_simple(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

	if (strcmp(p->Celltype, "RA") == 0); 	// Do nothing for baseline celltype
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
		printf("ERROR: \"%s\" is not a valid Celltype for the hAM_CAZ_simple models. Please check Model.c and Model_hAM_CAZ_simple.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

// ISO
void set_modulation_ISO_native_hAM_CAZ_simple(Cell_parameters *p)
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

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

	// testing exmaple illustration of model-specific implementation
	if (strcmp(p->ISO_model, "Toy_methods_demonstration") == 0)
	{
		p->GKur				*= (1.0 + p->ISO*0.5);
		p->Gto				*= (1.0 + p->ISO*0.5);		
		p->GKr				*= (1.0 + p->ISO*0.5);		
		p->Gup				*= (1.0 + p->ISO*0.7);    	
		p->GLTCC_kva1_va2   *= (1.0 + p->ISO*1.0);
	}
	// Add new ISO_models here: else if (strcmp(p->ISO_model, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid ISO model for the hAM_CAZ_simple models. Please check Model.c and Model_hAM_CAZ_simple.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hAM_CAZ_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hAM_CAZ_simple models. Please check Model.c and Model_hAM_CAZ_simple.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hAM_CAZ_simple(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

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
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hAM_CAZ_simple models. Please check Model.c and Model_hAM_CAZ_simple.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hAM_CAZ_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Mutation model for the hAM_CAZ_simple models. Please check Model.c and Model_hAM_CAZ_simple.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}
// ACh ===================================\\|
void set_modulation_ACh_hAM_CAZ_simple(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// Note: setting X_set_ref to 1 and then 0 if not set is important for error checking
	// and for being able to jump into global/common AND specific functions properly

	// NOTE: if the model does not have IKACh (or If) in, then modifying them won't make a difference!
	if (strcmp(p->ACh_model, "default") == 0)
	{
		// Put actual things in here
	}
	else if (strcmp(p->ACh_model, "test") == 0)  // just an example of how to implement: NOT an actual implementation
	{
		p->gKACh    = p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) );
		p->If_va_ss_shift += -7.2*p->ACh;
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid ACh for the hAM_CAZ model. Please check Model_hAM_CAZ_simple.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
// Your model may have more or fewer currents than this template - just follow the procedure and add/delete as appropriate
void compute_model_hAM_CAZ_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hAM_CAZ_simple_native(p, var, Vm, s->Cai);
	update_gating_variables_hAM_CAZ_simple_native(p, var, s, Vm, dt);
	compute_Itot_hAM_CAZ_simple_integrated(p, var, s, Vm);
}

void set_gate_rates_hAM_CAZ_simple_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	// Call only current you need
	set_INa_LR_rates(p, var, Vm);					// lib/Model.c
	set_Ito_hAM_MT_rates(p, var, Vm);
	set_IKur_hAM_MT_rates(p, var, Vm);
	set_IKs_hAM_CAZ_simple_rates(p, var, Vm);
	set_IKr_hAM_CAZ_simple_rates(p, var, Vm);
	set_IK1_hAM_CAZ_simple_variables(p, var, Vm);
}

void update_gating_variables_hAM_CAZ_simple_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_INa_LR(p, var, s, Vm, dt); 		// lib/Model.c	
	update_gates_Ito_hAM_MT(p, var, s, Vm, dt);
	update_gates_IKur_hAM_MT(p, var, s, Vm, dt);
	update_gates_IKs_hAM_CAZ_simple(p, var, s, Vm, dt);
	update_gates_IKr_hAM_CAZ_simple(p, var, s, Vm, dt);
}

void compute_Itot_hAM_CAZ_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_INa_LR(p, var, s, Vm);                  // lib/Model.c
	compute_Ito_hAM_MT(p, var, s, Vm);
	compute_IKur_hAM_MT(p, var, s, Vm);
	compute_IKr_hAM_CAZ_simple(p, var, s, Vm);
	compute_IKs_hAM_CAZ_simple(p, var, s, Vm);
	compute_IK1_hAM_CAZ_simple(p, var, s, Vm);

	var->Itot	= var->INa + var->Ito + var->IKur + var->IKr + var->IKs + var->IK1;
	var->Itot   += var->INCX + var->ICaP + var->ICab + var->ICaL;   // Adds Ca currents computed in CRU
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
// Identical to that of LR model, found in lib/Model.c or new formulation here
// End INa ==================================================================//|

// Ito ======================================================================\\|
// Identical to Maleckar et al formulation: lib/Model_hAM_MT.cpp
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
void set_ICaL_hAM_CAZ_simple_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss                          = sigmoid(Vm_ss, 5.0, -8.0*kscale);
	if (Vm_tau != 5.00)
	{
		*va_tau                     = (1/(1+exp((Vm_tau-5)/-6.24)))*(1-exp((Vm_tau-5)/-6.24))/(0.035*(Vm_tau-5));
	}
	else *va_tau                    = (1/(1+exp((Vm_tau-5)/-6.24))) * 4.579;
}

void set_ICaL_hAM_CAZ_simple_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss      = sigmoid(Vm_ss, -32.06, 8.6*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*vi_tau  	= 2.0/(0.0197*exp(- (0.0337*(Vm_tau - (10)))*(0.0337*(Vm_tau - (10))) ) + 0.02 );
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
// MT formulation lib/Model_hAM_MT.cpp
// End IKur =================================================================//|

// IKr ======================================================================\\|
void set_IKr_hAM_CAZ_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKr_va_ss		= sigmoid(Vm_ac_ss, -14.10, -6.5*p.IKr_va_ss_kscale);

	if (fabs(Vm_ac_tau+14.1)< 1e-10) var->IKr_va_al  = 0.0015; // Denominator = 0 clause
	else var->IKr_va_al                       = 0.0003*(Vm_ac_tau+14.1)/(1-exp((Vm_ac_tau+14.1)/-5));
	// Beta
	if (fabs(Vm_ac_tau-3.3328) < 1e-10) var->IKr_va_bet = 3.7836118e-4;
	else var->IKr_va_bet           			            = 0.000073898*(Vm_ac_tau-3.3328)/(exp((Vm_ac_tau-3.3328)/5.1237)-1);

	var->IKr_va_tau		= 1.0/(var->IKr_va_al + var->IKr_va_bet);
	var->IKr_va_tau		*= p.IKr_va_tau_scale;

	// Time-independant gate
	var->IKr_vi_ti          = sigmoid(Vm_ac_ss, -15, 22.4);
}

void update_gates_IKr_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKr_va             = rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);
}

void compute_IKr_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKr				= p.gKr * s->IKr_va * var->IKr_vi_ti * (Vm - var->EK);
	var->IKr				*= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_hAM_CAZ_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKs_va_ss                    = sqrt(1.0/(1+exp((Vm_ac_ss - (19.9))/(-12.7))));

	if (fabs(Vm_ac_tau-19.9) < 1e-10) /* denominator = 0 */
	{
		var->IKs_va_al                = 0.00068;
		var->IKs_va_bet                 = 0.000315;
	}
	else
	{
		var->IKs_va_al                = 0.00004*(Vm_ac_tau-19.9)/(1-exp((Vm_ac_tau - 19.9)/-17));
		var->IKs_va_bet               = 0.000035*(Vm_ac_tau-19.9)/(exp((Vm_ac_tau - 19.9)/9)-1);
	}

	var->IKs_va_tau                 = 0.5/(var->IKs_va_al + var->IKs_va_bet);
	var->IKs_va_tau                 *= p.IKs_va_tau_scale;
}

void update_gates_IKs_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKs_va             = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
}

void compute_IKs_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKs                = p.gKs * s->IKs_va * s->IKs_va * (Vm - var->EK);
	var->IKs                *= p.GKs;
}
// End IKs ==================================================================//|

// IK1 ======================================================================\\|
void set_IK1_hAM_CAZ_simple_variables(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_in        = Vm - p.IK1_va_shift;
	var->IK1_va_ti      = (1.0 + exp(0.07*(Vm_in-(-80))));
}

void compute_IK1_hAM_CAZ_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IK1               = p.gK1 * (Vm - var->EK - p.IK1_Erev_shift)/var->IK1_va_ti;
	var->IK1               *= p.GK1;
}
// End IK1 ==================================================================//|

