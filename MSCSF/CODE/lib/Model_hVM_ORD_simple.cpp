// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of a simplified version of ===  //
// the human ventricular cell model by O'Hara-Rudy: =======  //
// O'Hara et al. 2011 PLOS Comp. Biol. 7, e1002061. =======  //
// Simplified model is for integration with spatial Ca2+ ==  //
// handling model. Original full model is not packaged with  //
// this code. =============================================  //
// In-code identifier: hVM_ORD_s ==========================  //
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

// Disclaimer copied from  source code for ORD model ==================================\\|
// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy, Washington University in St. Louis.
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

// C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the
// undiseased human ventricular action potential and calcium transient
// The ORd model is described in the article "Simulation of the Undiseased
// Human Cardiac Ventricular Action Potential: Model Formulation and
// Experimental Validation"
// by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
// The article and supplemental materails are freely available in the
// Open Access jounal PLoS Computational Biology
// Link to Article:
// http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
// Email: tom.ohara@gmail.com / rudy@wustl.edu
// Web: http://rudylab.wustl.edu
// End Disclaimer in source code for ORD model ===================================//|

#include "Model.h"
#include "Structs.h"

// Parameters and specific settings =============================================================\\|
// Set model dependent parameters 
// If your model does not inhereted from a cell model already included, define all parameters:
void set_parameters_native_hVM_ORD_simple(Cell_parameters *p)
{
	// Capacitance and cell structure
	//p->Cm				= X;		// pF
	//p->Vcell			= Y;	// um^3

	// Concentrations (constant OR initial condition)
	p->Ko					= 4.5;		// mM
	p->Nao					= 140;		// mM
	p->Cao					= 1.8;		// mM
	p->Ki					= 143.103; 	// mM
	p->Cai					= 0.1;		// mM
	p->Nai					= 7.95;		// mM

	// Current parameters
	p->gNa					= 16.0;		// s/mF
	p->gNaL					= 0.0075;	// s/mF
	p->gto					= 0.02;		// s/mF
	p->gKr					= 1.5*0.0529; // s/mF || Mod from original as part of simplification
	p->gK1					= 0.3;		// s/mF
	p->gKs					= 1.5*0.00391; // s/mF || Mod from original as part of simplification
	p->gKb					= 0.003;	// s/mF
	p->gNab					= 3.75e-10;	// s/mF 

	// Stimulus parameters (if different to default)
	p->stimduration 	= 5.0;      // ms
	p->stimmag      	= -12.5;    // pA/pF

	// Default celltype
	p->Celltype			= "EPI";
}

// Initial conditions
void initial_conditions_native_hVM_ORD_simple(State_variables *s, Cell_parameters p)
{
	// Define ICs for ALL state variables used in the model
	s->Vm      			= -85;
	s->INa_va  			= 0.001231;
	s->INa_vi_1			= 0.992842; 
	s->INa_vi_2			= 0.988864;
	s->INaL_va			= 0.000195011;
	s->INaL_vi			= 0.438143;
	s->Ito_va			= 0;
	s->Ito_vi			= 0.999539;
	s->Ito_vi_s			= 0.361531;
	s->ICaL_va			= 0; // not actually used as integrated only
	s->ICaL_vi			= 1; // not actually used as integrated only
	s->IKr_va			= 0.000264226;
	s->IKr_vi			= 0.669; 	// va s
	s->IKs_va			= 0.397626;	
	s->IKs_va_2			= 0.000197522;
	s->IK1_va			= 0.99681;

	// Ca handling
	//s->CajSR				= X;		// mM
	//s->RyRo				= X;

	// Assign state from param (even if constant) for concentrations
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
void set_het_mod_hVM_ORD_simple(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hVM_ORD_simple(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hVM_ORD_simple(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hVM_ORD_simple(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hVM_ORD_simple(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hVM_ORD_simple(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hVM_ORD_simple(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hVM_ORD_simple(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

	if (strcmp(p->Celltype, "ENDO") == 0)
	{
		p->Gto				*= 4;
	}
	else if (strcmp(p->Celltype, "M") == 0)
	{   
		p->Gto              *= 4;
		p->GKr				*= 0.8;	
		p->GK1				*= 1.3;
	}
	else if (strcmp(p->Celltype, "EPI") == 0)
	{   
		p->Gto              *= 4;
		p->GKr				*= 1.3;
		p->GKs				*= 1.4;	
		p->GK1				*= 1.2;
		p->GNaL				*= 0.6;
	}
	else 
	{
		printf("ERROR: \"%s\" is not a valid Celltype for the hVM_ORD_simple models. Please check Model.c and Model_hVM_ORD_simple.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

// ISO
void set_modulation_ISO_native_hVM_ORD_simple(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// "p->ISO" is concentration of ISO; "p-ISO_model" is the ISO model to be applied

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c

	// testing exmaple illustration of model-specific implementation
	if (strcmp(p->ISO_model, "Toy_methods_demonstration") == 0)
	{
		p->GKr				*= (1.0 + p->ISO*(2.0 - 1.0)); // x2
		p->GLTCC_kva1_va2	*= (1.0 + p->ISO*(2.0 - 1.0)); // x2
		p->Gup				*= (1.0 + p->ISO*(1.75- 1.0)); // x1.75
	}
	// Add new ISO_models here: else if (strcmp(p->ISO_model, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid ISO model for the hVM_ORD_simple models. Please check Model.c and Model_hVM_ORD_simple.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hVM_ORD_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hVM_ORD_simple models. Please check Model.c and Model_hVM_ORD_simple.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hVM_ORD_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hVM_ORD_simple models. Please check Model.c and Model_hVM_ORD_simple.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hVM_ORD_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Mutation model for the hVM_ORD_simple models. Please check Model.c and Model_hVM_ORD_simple.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}

// ACh ===================================\\|
void set_modulation_ACh_hVM_ORD_simple(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ACh for the ORD_simple model. Please check Model_hVM_ORD_simple.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
// Your model may have more or fewer currents than this template - just follow the procedure and add/delete as appropriate
void compute_model_hVM_ORD_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hVM_ORD_simple_native(p, var, Vm, s->Cai);
	update_gating_variables_hVM_ORD_simple_native(p, var, s, Vm, dt);
	compute_Itot_hVM_ORD_simple_integrated(p, var, s, Vm);
}

void set_gate_rates_hVM_ORD_simple_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	// Call only current you need
	set_INa_LR_rates(p, var, Vm);					// lib/Model.c
	set_INaL_hVM_ORD_simple_rates(p, var, Vm);
	set_Ito_hVM_ORD_simple_rates(p, var, Vm);
	set_IKs_hVM_ORD_simple_rates(p, var, Vm);
	set_IKr_hVM_ORD_simple_rates(p, var, Vm);
}

void update_gating_variables_hVM_ORD_simple_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_INa_LR(p, var, s, Vm, dt); 		// lib/Model.c
	update_gates_INaL_hVM_ORD_simple(p, var, s, Vm, dt);
	update_gates_IKs_hVM_ORD_simple(p, var, s, Vm, dt);
	update_gates_IKr_hVM_ORD_simple(p, var, s, Vm, dt);
	update_gates_Ito_hVM_ORD_simple(p, var, s, Vm, dt);
	update_gates_IK1_hVM_ORD_simple(p, var, s, Vm, dt);
}

void compute_Itot_hVM_ORD_simple_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_INa_LR(p, var, s, Vm);					// lib/Model.c
	compute_INaL_hVM_ORD_simple(p, var, s, Vm);
	compute_Ito_hVM_ORD_simple(p, var, s, Vm);
	compute_IKr_hVM_ORD_simple(p, var, s, Vm);
	compute_IKs_hVM_ORD_simple(p, var, s, Vm);
	compute_IK1_hVM_ORD_simple(p, var, s, Vm);
	compute_IKb_hVM_ORD_simple(p, var, s, Vm);
	compute_INab_hVM_ORD_simple(p, var, s, Vm);

	var->Itot	= var->INa + var->INaL + var->Ito + var->IKr + var->IKs + var->IK1 + var->IKb + var->INab;
	var->Itot   += var->INCX + var->ICaP + var->ICab + var->ICaL;   // Adds Ca currents computed in CRU
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
// Identical to that of LR model, found in lib/Model.c or new formulation here
// End INa ==================================================================//|

// INaL =====================================================================\\|
void set_INaL_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss  	= Vm - p.INaL_va_ss_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac_ss 	= Vm - p.INaL_vi_ss_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation
	double Vm_ac_tau  	= Vm - p.INaL_va_tau_shift;  // Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac_tau 	= Vm - p.INaL_vi_tau_shift;  // Shift of the voltage used to calculate alpha and beta, inactivation

	var->INaL_va_ss		= sigmoid(Vm_ac_ss, -42.85, -5.264*p.INaL_va_ss_kscale);
	var->INaL_va_tau	= 1.0/(6.765*exp((Vm_ac_tau + 11.64)/34.77)+8.552*exp(-(Vm_ac_tau+77.42)/5.955));
	var->INaL_va_tau	*= p.INaL_va_tau_scale;

	var->INaL_vi_ss		= sigmoid(Vm_inac_ss, -87.61, 7.488*p.INaL_vi_ss_kscale);
	var->INaL_vi_tau	= 200.0;
	var->INaL_vi_tau	*= p.INaL_vi_tau_scale;
}

void update_gates_INaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->INaL_va              = rush_larsen(s->INaL_va, var->INaL_va_ss, var->INaL_va_tau, dt);
	s->INaL_vi              = rush_larsen(s->INaL_vi, var->INaL_vi_ss, var->INaL_vi_tau, dt);
}

void compute_INaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->INaL				= p.gNaL * s->INaL_va * s->INaL_vi * (Vm - var->ENa);
	var->INaL				*= p.GNaL;
}
// End INaL =================================================================//|

// Ito ======================================================================\\|
void set_Ito_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	var->Ito_va_ss			= sigmoid(Vm_ac_ss, 14.34, -14.82*p.Ito_va_ss_kscale);
	var->Ito_va_tau			= 1.0515/(1.0/(1.2089*(1.0+exp(-(Vm_ac_tau-18.4099)/29.3814)))+3.5/(1.0+exp((Vm_ac_tau+100.0)/29.3814))); 
	var->Ito_va_tau			*= p.Ito_va_tau_scale;

	var->Ito_vi_ss			= sigmoid(Vm_inac_ss, -43.94, 5.711*p.Ito_vi_ss_kscale);
	var->Ito_vi_tau			= 4.562+1/(0.3933*exp((-(Vm_inac_tau+100.0))/100.0)+0.08004*exp((Vm_inac_tau+50.0)/16.59));
	var->Ito_vi_tau			*= p.Ito_vi_tau_scale;
	var->Ito_vi_s_tau		= 23.62+1/(0.001416*exp((-(Vm_inac_tau+96.52))/59.05)+1.780e-8*exp((Vm_inac_tau+114.1)/8.079));
	var->Ito_vi_s_tau  		*= p.Ito_vi_tau_scale;

	var->Ito_vi_Fs			= 1 - (1.0/(1.0+exp((Vm_inac_ss-214.6)/151.2))); 
}

void update_gates_Ito_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Update the gates, using rush_larsen, forward Euler, or other
	s->Ito_va   			= rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi   			= rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
	s->Ito_vi_s 			= rush_larsen(s->Ito_vi_s, var->Ito_vi_ss, var->Ito_vi_s_tau, dt); // same steady state as fast
}

void compute_Ito_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ito                    =  p.gto * s->Ito_va * ((1-var->Ito_vi_Fs)*s->Ito_vi + var->Ito_vi_Fs*s->Ito_vi_s) * (Vm - var->EK);
	var->Ito                    *= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
// ICaL can be written exactly as above for Ito. However, if there are intentions
// to integrate the new model with the "integrated" calcium handling system
// e.g. for spatial cell models or spontaneous release functions, please 
// follow the procedure below; voltage-dependent gates have their own
// functions so can be called elsewhere
void set_ICaL_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hVM_ORD_simple_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hVM_ORD_simple_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// calcium inactivation || Note included as only integrated implementation
}

void set_ICaL_hVM_ORD_simple_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss  		= sigmoid(Vm_ss, -3.940, -4.230*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*va_tau			= 0.6 + 1.0/(exp(-0.05*(Vm_tau - (-6.0))) + exp(0.09*(Vm_tau - -(14.0))));
}

void set_ICaL_hVM_ORD_simple_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss        	= sigmoid(Vm_ss, -19.58, 3.696*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*vi_tau			= 10.0 + 1.0/(0.000035*exp(-(Vm_tau+5.0)/4.0)+0.000035*exp((Vm_tau+5.0)/6.0));
}

// Everything else follows the same format
// Below not used as only integrated implementation, which takes care of this elsewhere (in lib/CRU.cpp)
/*void update_gates_ICaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
  {
  s->ICaL_va       		= rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
  s->ICaL_vi    			= rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
  }

  void compute_ICaL_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
  {
  var->ICaL   			= p.gCaL * s->ICaL_va * s->ICaL_vi * (Vm - 22);
  var->ICaL 				*= p.GCaL;
  }*/
// End ICaL =================================================================//|

// IKr ======================================================================\\|
// NOTE: "IKr_vi" gate and associated variables used for IKr_va_slow parameters
void set_IKr_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	// va | fast va
	var->IKr_va_ss			= sigmoid(Vm_ac_ss, -8.337, -6.789*p.IKr_va_ss_kscale);
	var->IKr_va_tau			= 12.98+1.0/(0.3652*exp((Vm_ac_tau-31.66)/3.869)+4.123e-5*exp((-(Vm_ac_tau-47.78))/20.38)); 
	var->IKr_va_tau			*= p.IKr_va_tau_scale;

	// slow gate  (va slow using vi variables)
	var->IKr_vi_tau			= 1.865+1.0/(0.06629*exp((Vm_ac_tau-34.70)/7.355)+1.128e-5*exp((-(Vm_ac_tau-29.74))/25.94));
	var->IKr_vi_tau			*= p.IKr_vi_tau_scale;

	// Time-independant gate
	var->IKr_vi_ti          = 1.0/(1.0+exp((Vm_ac_ss-(-55.0))/75.0))*1.0/(1.0+exp((Vm_ac_ss-(10.0))/30.0));
	var->IKr_vi_ti			*= sqrt(p.Ko/5.4);

	// Proportion fast/slow
	var->IKr_va_Fs			= 1 - (1.0/(1.0+exp((Vm_ac_ss+54.81)/38.21)));
}

// NOTE: "IKr_vi" gate and associated variables used for IKr_va_slow parameters
void update_gates_IKr_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKr_va				= rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);
	s->IKr_vi				= rush_larsen(s->IKr_vi, var->IKr_va_ss, var->IKr_vi_tau, dt);
}

// NOTE: "IKr_vi" gate and associated variables used for IKr_va_slow parameters
void compute_IKr_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKr				= p.gKr	* ((1-var->IKr_va_Fs)*s->IKr_va + var->IKr_va_Fs*s->IKr_vi) * var->IKr_vi_ti * (Vm - var->EK);
	var->IKr				*= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_hVM_ORD_simple_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKs_va_ss			= sigmoid(Vm_ac_ss, -11.60, -8.932*p.IKs_va_ss_kscale);
	var->IKs_va_tau			= 817.3+1.0/(2.326e-4*exp((Vm_ac_tau+48.28)/17.80)+0.001292*exp((-(Vm_ac_tau+210.0))/230.0));
	var->IKs_va_tau			*= p.IKs_va_tau_scale;
	var->IKs_va_2_tau		= 1.0/(0.01*exp((Vm_ac_tau-50.0)/20.0)+0.0193*exp((-(Vm_ac_tau+66.54))/31.0));
	var->IKs_va_2_tau		*= p.IKs_va_tau_scale;
}

void update_gates_IKs_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKs_va             = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
	s->IKs_va_2           = rush_larsen(s->IKs_va_2, var->IKs_va_ss, var->IKs_va_2_tau, dt);
}

void compute_IKs_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double IKs_ci;
	IKs_ci					= 1.0+0.6/(1.0+pow(3.8e-5/s->Cai,1.4));

	var->IKs				= p.gKs * s->IKs_va * s->IKs_va_2 * IKs_ci * (Vm - var->EKs_ORD);
	var->IKs				*= p.GKs;
}
// End IKs ==================================================================//|

// IK1 ======================================================================\\|
//void set_IK1_hVM_ORD_simple_variables(Cell_parameters p, Model_variables *var, double Vm)
//{
//	double Vm_in        = Vm - p.IK1_va_shift;
//	
//	//var->IK1_va_ti      = (1.0 + exp(0.07*(Vm_in-(-80))));
//}

void update_gates_IK1_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	double Vm_in        = Vm - p.IK1_va_shift;
	var->IK1_va_ss		= 1.0/(1.0+exp(-(Vm_in+2.5538*s->Ko+144.59)/(1.5692*s->Ko+3.8115)));
	var->IK1_va_tau		= 122.2/(exp((-(Vm_in+127.2))/20.36)+exp((Vm_in+236.8)/69.33));
	var->IK1_va_ti		= 1.0/(1.0+exp((Vm_in+105.8-2.6*s->Ko)/9.493));
	var->IK1_va_ti		*= sqrt(s->Ko);
	s->IK1_va			= rush_larsen(s->IK1_va, var->IK1_va_ss, var->IK1_va_tau, dt);
}

void compute_IK1_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IK1               = p.gK1 * s->IK1_va * var->IK1_va_ti * (Vm - var->EK - p.IK1_Erev_shift);
	var->IK1               *= p.GK1;
}
// End IK1 ==================================================================//|

// IKb ======================================================================\\|
void compute_IKb_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKb	= p.gKb * 1.0/(1.0+exp(-(Vm-14.48)/18.34)) * (Vm - var->EK);
	var->IKb	*= p.GKb;
}
// End IKb ==================================================================//|

// INab
void compute_INab_hVM_ORD_simple(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double vfrt		= Vm * p.FoRT;
	double vffrt	= vfrt * p.F;

	var->INab	= p.gNab * vffrt*(s->Nai*exp(vfrt)-s->Nao)/(exp(vfrt)-1.0);
	var->INab	*= p.GNab;
}

