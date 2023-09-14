// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of the human atrial cell model  //
// by Colman-Workman: Colman et al. 2018 Front. Physiol. ==  //
// 9 1211, and incorporating components from: =============  //
// hAM_CRN ================================================  //
// "Ionic mechanisms underlying human atrial action =======  //
// potential properties: insights from a mathematical =====  //
// model." M. Courtemanche, R.J. Ramirez and S. Nattel. Am=  //
// J Physiol 1998; 275(1 Pt 2):H301-21. ===================  //
// hAM_NG =================================================  //
// "Mathematical model of an adult human atrial cell: the==  //
// role of K+ currents in repolarization." A. Nygren, C. ==  //
// Fiset, L. Firek, J.W. Clark, D.S. Lindblad, R.B. Clark =  //
// W.R. Giles. Circ Res 1998;82(1):63-81 ==================  //
// hAM_GB =================================================  //
// "Human atrial action potential and Ca2+ model: sinus ===  //
// rhythm and chronic atrial fibrillation." E. Grandi, ====  //
// S.V. Pandit, N. Voigt, A.J. Workman, D. Dobrev, ========  //
// J. Jalife, D.M. Bers. Circ Res 2011; 109(9):1055-66 ====  //
// In-code identifier: hAM_WL_X, hAM_X_mWL ================  //
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

// Models:
//  hAM_WL_CRN      - WL AP, CRN Ca2+
//  hAM_WL_GB       - WL AP, GB Ca2+
//  hAM_CRN_mWL     - CRN + WL INa, Ito, IKur, IK1 and ICaL_mod
//  hAM_GB_mWL      - GB + WL INa, Ito, IKur, IK1 and ICaL_mod
//  hAM_NG_mWL      - NG + WL INa, Ito, IKur, IK1 and ICaL_mod

#include "Model.h"
#include "Structs.h"

// Parameters and specific settings =============================================================\\|
// Set model dependent parameters (updates params set by baseline model (CRN, GB)
void update_parameters_native_hAM_WL(Cell_parameters *p)
{
	// Ca handling model selection
	if (strcmp(p->Model, "hAM_WL_CRN") == 0	|| strcmp(p->Model, "hAM_CRN_mWL") == 0) // if CRN Ca handling
	{
		p->Ca_handling	= "CRN";
	}
	else if (strcmp(p->Model, "hAM_WL_GB") == 0 || strcmp(p->Model, "hAM_GB_mWL") == 0) // if GB Ca handling
	{
		p->Ca_handling   = "GB";
	}
	else if (strcmp(p->Model, "hAM_NG_mWL") == 0)	// if NG Ca handling
	{	
		p->Ca_handling   = "NG";
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid model type, Ca handling system cannot be set (lib/Model_hAM_WL.cpp)\n\n", p->Model);
		exit(1);
	}

	// WL currents
	p->gNa              = 17.55;
	p->gto				= 0.102816;
	p->gKur				= 0.067626;
	p->ICaL_vi_Fs		= 0.0;	
	p->pCaL         	= 2.7e-4; 	// overwritten for each model  

	// Updates to maintain CICR
	if (strcmp(p->Ca_handling, "GB") == 0)
	{
		p->Gleak        *= 0.3; 
		p->ICaP_bar     *= 0.3;//0.5*0.6;
		p->INCX_bar     *= 0.64;//0.8*0.8;
	}

	// Model-dependant parameters
	if (strcmp(p->Model, "hAM_WL_CRN") == 0)  // hAM WL minimal CRN Ca handling
	{
		p->pCaL      	= 5.346e-4;	//1.25 * 1.1*1.6*0.9*2.7e-4;
		p->ICaL_vi_Fs	= 0.2;

		p->Nai       	= 11.2;     // mM
		p->Nao       	= 140;      // mM
		p->Ki        	= 139;      // mM
		p->Ko        	= 4;        // mM
		p->Cai       	= 0.000102; // mM 
		p->Cao      	= 1.8;      // mM	
	}
	else if (strcmp(p->Model, "hAM_WL_GB") == 0)  // hAM WL minimal GB Ca handling
	{
		p->pCaL         = 6.4152e-4;//1.1*1.5 * 1.6*0.9*2.7e-4;
		p->ICaL_vi_Fs   = 0.2;

		p->gKr			*= 4.5;
		p->gKs			*= 4.5;
	}
	// Modified WL models
	else if (strcmp(p->Model, "hAM_CRN_mWL") == 0)
	{
		p->gCaL         = 1.725*0.1238;	// scale * original	
	}
	else if (strcmp(p->Model, "hAM_GB_mWL") == 0)
	{
		p->pCaL         = 0.9*2.7e-4; 
	}
	else if (strcmp(p->Model, "hAM_NG_mWL") == 0)
	{
		p->gCaL         = 1.28*0.135000;	
		if (strcmp(p->environment, "intact") == 0)
		{
			p->ICaP_bar     *= 1.5;		// Update to maintain CICR stable in intact
			p->gCab         *= 0.3;
			p->gKr			*= 5;
			p->gKs			*= 5;
		}
	}

	// Stimulus
	p->stimduration     = 5.0;      // ms
	p->stimmag          = -13.5;    // pA/pF

	// Default celltype
	p->Celltype			= "RA";
	p->ISO_model		= "Col";

	// Ihyp
	p->AIhyp			= 0.63; 	// pA/pF || will be overwritten if argument specefier is passed
}

// Heterogeneity and modulation =================================================================\\|
// Parent function
void set_het_mod_hAM_WL(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_hAM_WL(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_hAM_WL(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_hAM_WL(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_hAM_WL(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_hAM_WL(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_hAM_WL(p);
}

// Celltype/Heterogeneity
void set_celltype_native_hAM_WL(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	// NOTE:: Most het has been set in Model.c as common to multiple human atrial cell models. Please see that function first.
	// set_celltype_hAM() in lib/Model.c

	if (strcmp(p->Celltype, "RA") == 0);	// Default, do nothing
	else if (strcmp(p->Celltype, "CRN_scaling") == 0)
	{
		p->GCaL				*= (1.0/1.31);
		p->Gto				*= (1.0/0.72);
		p->GKur				*= (1.0/1.24);
	}
	else if (strcmp(p->Celltype, "GB_scaling") == 0)
	{
		p->GCaL             *= (1.0/0.61);
		p->Gto              *= (1.0/0.93);
		p->GKur             *= (1.0/1.48);
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
		printf("ERROR: \"%s\" is not a valid Celltype for the hAM_WL models. Please check Model.c and Model_hAM_WL.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

// ISO
void set_modulation_ISO_native_hAM_WL(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ISO model for the hAM_WL models. Please check Model.c and Model_hAM_WL.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_hAM_WL(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the hAM_WL models. Please check Model.c and Model_hAM_WL.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_hAM_WL(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Remodelling model for the hAM_WL models. Please check Model.c and Model_hAM_WL.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

// Mutation
void set_modulation_Mutation_native_hAM_WL(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid Mutation model for the hAM_WL models. Please check Model.c and Model_hAM_WL.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}

// ACh ===================================\\|
void set_modulation_ACh_hAM_WL(Cell_parameters *p)
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
	else */if (strcmp(p->ACh_model, "test") == 0)  // just an example of how to implement: NOT an actual implementation
	{
		p->gKACh    = p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) );
		p->If_va_ss_shift += -7.2*p->ACh;
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid ACh for the hAM_WL model. Please check Model_hAM_WL.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End Heterogeneity and modulation =============================================================//|

// Initial conditions
void initial_conditions_native_hAM_WL(State_variables *s, Cell_parameters p)
{
	// Common ICs
	s->Vm               = -82;
	s->INa_va           = 0.00291;
	s->INa_vi_1         = 0.9791;
	s->INa_vi_2         = 0.9869;
	s->Ito_va           = 0.07;
	s->Ito_vi           = 0.99;
	s->Ito_vi_s         = 0.361531;
	s->ICaL_va          = 0.0;
	s->ICaL_vi          = 0.0;
	s->ICaL_vi_s        = 0.0;
	s->IKur_va          = 0.0;
	s->IKur_vi          = 1.0;
	s->IKr_va           = 0.0;
	s->IKs_va           = 0.0;
	s->INaL_va          = 0.000102;
	s->INaL_vi          = 0.050442;

	s->ICaL_ci          = 1.0;
	s->ICaL_ci_j        = 1.0;

	if (strcmp(p.Model, "hAM_WL_CRN") == 0 || strcmp(p.Model, "hAM_CRN_mWL") == 0)
	{
		s->Vm               = -80;
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
		s->CajSR            = 1.49;     // mM
		s->CanSR            = 1.49;     // mM
		s->cmdn             = 0.00205;  // mM
		s->trpn             = 0.0118;   // mM
		s->csqn             = 6.432;    // mM
		s->RyRo             = 0.0;
		s->RyRr             = 1.0;
		s->RyRi             = 0.9992;
	}

	if (strcmp(p.Model, "hAM_WL_GB") == 0 || strcmp(p.Model, "hAM_GB_mWL") == 0)
	{
		s->Vm               = -80;
		s->INa_va           = 0.010085;
		s->INa_vi_1         = 0.829053;
		s->INa_vi_2         = 0.849974;
		s->INaL_va          = 0.000102;
		s->INaL_vi          = 0.050442;
		s->Ito_va           = 0.07;
		s->Ito_vi           = 0.99;
		s->Ito_vi_s         = 0.361531;
		s->ICaL_va          = 0.0;
		s->ICaL_vi          = 1.0;
		s->ICaL_ci          = 0;//1.0;
		s->ICaL_ci_j        = 0;//1.0;
		s->IKur_va          = 0.0;
		s->IKur_vi          = 1.0;
		s->IKr_va           = 0.000024;
		s->IKs_va           = 0;

		// Ca handling
		s->CajSR            = 0.52;     // mM
		s->CanSR            = 0.52;     // mM
		s->RyRo             = 0.000002;
		s->RyRi             = 0.0;
		s->RyRr             = 0.841294;
		s->Myo_c            = 0.003737;
		s->Myo_m            = 0.135752;
		s->Tn_CHc           = 0.127315;
		s->Tn_CHm           = 0.005947;
		s->Tn_CL            = 0.017693;
	}

	if (strcmp(p.Model, "hAM_NG_mWL") == 0)
	{
		s->Vm               = -80;
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
	}

	// Assign state from param (even if constant)
	s->Nai              = p.Nai;
	s->Nai_sl           = p.Nai;
	s->Nai_j            = p.Nai;
	s->Nao              = p.Nao;
	s->Ki               = p.Ki;
	s->Ko               = p.Ko;
	s->Cai              = p.Cai;
	s->Cai_j            = p.Cai;
	s->Cai_sl           = p.Cai;
	s->Cao              = p.Cao;
}
// end Parameters and specific settings =========================================================//|

// Compute model functions ======================================================================\\|
void compute_model_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	if (strcmp(p.Ca_handling, "CRN") == 0) s->Cai_sl   =   s->Cai; // so functions can read Cai_sl for CRN or GB
	compute_reversal_potentials(p, var, s);
	set_gate_rates_hAM_WL_native(p, var, Vm, s->Cai);
	update_gating_variables_hAM_WL_native(p, var, s, Vm, dt);
	compute_Itot_hAM_WL_native(p, var, s, Vm);

	if (strcmp(p.Ca_handling, "CRN") == 0) 				comp_homeostasis_hAM_CRN(p, var, s, Vm, dt); 
	else if (strcmp(p.Ca_handling, "GB") == 0)			comp_homeostasis_hAM_GB(p, var, s, Vm, dt);
	else if (strcmp(p.Ca_handling, "NG") == 0)			comp_homeostasis_hAM_NG(p, var, s, Vm, dt);
}

void set_gate_rates_hAM_WL_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	// WL currents
	set_INa_hAM_WL_rates(p, var, Vm);
	set_Ito_hAM_WL_rates(p, var, Vm);
	set_IKur_hAM_WL_rates(p, var, Vm);

	// WL ICaL
	if (strcmp(p.Model, "hAM_WL_CRN") == 0 || strcmp(p.Model, "hAM_WL_GB") == 0)  set_ICaL_hAM_WL_rates(p, var, Vm, Cai);
	// mWL ICaL
	else if (strcmp(p.Model, "hAM_CRN_mWL") == 0)   set_ICaL_hAM_CRN_mWL_rates(p, var, Vm, Cai);
	else if (strcmp(p.Model, "hAM_GB_mWL") == 0)	set_ICaL_hAM_GB_mWL_rates(p, var, Vm, Cai);
	else if (strcmp(p.Model, "hAM_NG_mWL") == 0)    set_ICaL_hAM_NG_mWL_rates(p, var, Vm, Cai);

	// Ca handling model dependent currents (i.e., inherited from native models - not necessarily currents which are involved in Ca handling)
	if (strcmp(p.Ca_handling, "CRN") == 0)
	{
		set_IKs_hAM_CRN_rates(p, var, Vm);
		set_IKr_hAM_CRN_rates(p, var, Vm);	
	}
	else if (strcmp(p.Ca_handling, "GB") == 0)
	{
		set_INaL_hAM_GB_rates(p, var, Vm);
		set_IKr_hAM_GB_rates(p, var, Vm);
		set_IKs_hAM_GB_rates(p, var, Vm);
	}	
	else if (strcmp(p.Ca_handling, "NG") == 0)
	{
		set_IKs_hAM_NG_rates(p, var, Vm);
		set_IKr_hAM_NG_rates(p, var, Vm);
	}	
}

void update_gating_variables_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// WL specific currents
	update_gates_INa_LR(p, var, s, Vm, dt);         // lib/Model.c
	update_gates_Ito_hAM_WL(p, var, s, Vm, dt);
	update_gates_IKur_hAM_WL(p, var, s, Vm, dt);	

	update_gates_ICaL_hAM_WL(p, var, s, Vm, dt);    // Updates v gates | Ci gates below	
	if (strcmp(p.Ca_handling, "CRN") == 0) 		update_gates_ICaL_hAM_WL_CRN_ci(p, var, s, Vm, dt);
	else if (strcmp(p.Ca_handling, "GB") == 0) 	update_gates_ICaL_hAM_WL_GB_ci(p, var, s, Vm, dt);
	else if (strcmp(p.Ca_handling, "NG") == 0) 	update_gates_ICaL_hAM_WL_NG_ci(p, var, s, Vm, dt);

	if (strcmp(p.Ca_handling, "CRN") == 0)
	{	
		update_gates_IKs_hAM_CRN(p, var, s, Vm, dt);
		update_gates_IKr_hAM_CRN(p, var, s, Vm, dt);
	}
	else if (strcmp(p.Ca_handling, "GB") == 0)
	{	
		update_gates_INaL_hAM_GB(p, var, s, Vm, dt);
		update_gates_IKr_hAM_GB(p, var, s, Vm, dt);
		update_gates_IKs_hAM_GB(p, var, s, Vm, dt);		
	}
	else if (strcmp(p.Ca_handling, "NG") == 0)
	{
		update_gates_IKs_hAM_NG(p, var, s, Vm, dt);
		update_gates_IKr_hAM_NG(p, var, s, Vm, dt);
	}
}

void compute_Itot_hAM_WL_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	// WL common functions
	compute_Ito_hAM_WL(p, var, s, Vm);
	compute_IKur_hAM_WL(p, var, s, Vm);

	// Intact vs isolated IK1
	if (strcmp(p.environment, "isolated") == 0) compute_IK1_hAM_WL_isolated(p, var, s, Vm);
	else compute_IK1_hAM_WL_intact(p, var, s, Vm);

	if (strcmp(p.Ca_handling, "CRN") == 0)
	{
		compute_INa_LR(p, var, s, Vm);                  // lib/Model.c

		compute_INCX_hAM_CRN(p, var, s, Vm); 
		compute_INaK_hAM_CRN(p, var, s, Vm);
		compute_ICaP_hAM_CRN(p, var, s, Vm);
		compute_INab_hAM_CRN(p, var, s, Vm);
		compute_ICab_hAM_CRN(p, var, s, Vm);

		compute_IKr_hAM_CRN(p, var, s, Vm);
		compute_IKs_hAM_CRN(p, var, s, Vm);

		if (strcmp(p.Model, "hAM_WL_CRN") == 0)			compute_ICaL_hAM_WL_CRN_bar(p, var, s, Vm, s->Cai);
		else if (strcmp(p.Model, "hAM_CRN_mWL") == 0)	compute_ICaL_hAM_CRN_mWL(p, var, s, Vm);

		// GB currents not in CRN thus need to be zeroed
		var->IClCa = var->IClb = var->INaL = var->IKb = 0;
	}
	else if (strcmp(p.Ca_handling, "GB") == 0)
	{
		compute_INa_hAM_GB(p, var, s, Vm);
		compute_INab_hAM_GB(p, var, s, Vm);
		compute_ICab_hAM_GB(p, var, s, Vm);
		compute_ICaP_hAM_GB(p, var, s, Vm);
		compute_INCX_hAM_GB(p, var, s, Vm);

		compute_INaL_hAM_GB(p, var, s, Vm);
		compute_IKb_hAM_GB(p, var, s, Vm);
		compute_IKr_hAM_GB(p, var, s, Vm);
		compute_IKs_hAM_GB(p, var, s, Vm);
		compute_INaK_hAM_GB(p, var, s, Vm);
		compute_IClCa_hAM_GB(p, var, s, Vm);
		compute_IClb_hAM_GB(p, var, s, Vm);

		compute_ICaL_hAM_WL_GB_bar(p, var, s, Vm); // common to both ICaL types

		if (strcmp(p.environment, "isolated") == 0)
		{
			var->IClCa 	*= 0.5;
			var->IClb 	*= 0.5;
			var->IKb	*= 0.2;
		}
	}
	else if (strcmp(p.Ca_handling, "NG") == 0)
	{
		compute_INa_LR(p, var, s, Vm);                  // lib/Model.c
		compute_IKr_hAM_NG(p, var, s, Vm);
		compute_IKs_hAM_NG(p, var, s, Vm);
		compute_INCX_hAM_NG(p, var, s, Vm);
		compute_INaK_hAM_NG(p, var, s, Vm);
		compute_ICaP_hAM_NG(p, var, s, Vm);
		compute_INab_hAM_NG(p, var, s, Vm);
		compute_ICab_hAM_NG(p, var, s, Vm);

		compute_ICaL_hAM_NG(p, var, s, Vm);

		// GB currents not in NG thus need to be zeroed
		var->IClCa = var->IClb = var->INaL = var->IKb = 0;
	}

	if (strcmp(p.environment, "isolated") == 0) var->IKs = var->IKr = 0;

	var->Itot   = var->INa + var->IK1 + var->INab + var->ICab + var->ICaP + var->INCX + var->Ito + var->ICaL + var->IKur + var->INaK + var->IKr + var->IKs + var->IClCa + var->IClb + var->INaL + var->IKb;

	if (strcmp(p.environment, "isolated") == 0) var->Itot 	+=  p.AIhyp; 	// add hyperpolarizing current if isolated conditions
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
void set_INa_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac        = Vm - 8 - p.INa_va_shift;
	double Vm_inac      = Vm + 5 - p.INa_vi_shift;

	// Set Activation gate alpha and beta
	var->INa_va_al                 = 1.5 * 0.32*(Vm_ac+47.13)/(1-exp(-0.09*(Vm_ac+47.13))); // 1.1 - 1.25 x added for tissue
	var->INa_va_bet                = 0.08*exp(-Vm_ac/11.0);

	// Set inactivation gates alphas and betas
	if (Vm_inac < -40.0)
	{
		var->INa_vi_1_al           = 0.135*exp((80+Vm_inac)/-6.8);
		var->INa_vi_1_bet          = 0.9*3.56*exp(0.079*Vm_inac)+310000*exp(0.35*Vm_inac);
		var->INa_vi_2_al           = (-127140*exp(0.2444*Vm_inac)-0.00003474*exp(-0.04391*Vm_inac))*((Vm_inac+37.78)/(1+exp(0.311*(Vm_inac+79.23))));
		var->INa_vi_2_bet          = 0.9*(0.1212*exp(-0.01052*Vm_inac))/(1+exp(-0.1378*(Vm_inac+40.14)));
	}
	else
	{
		var->INa_vi_1_al           = 0;
		var->INa_vi_1_bet          = 1.0/(0.13*(1+exp((Vm_inac+10.66)/-11.1)));
		var->INa_vi_2_al           = 0;
		var->INa_vi_2_bet          = (0.3*exp(-0.0000002535*Vm_inac))/(1+exp(-0.1*(Vm_inac+32)));
	}

	// Set tau and SS from alpha and beta
	var->INa_va_tau                = 1/(var->INa_va_al + var->INa_va_bet); // 1/(a+b)
	var->INa_vi_1_tau              = 1/(var->INa_vi_1_al + var->INa_vi_1_bet);
	var->INa_vi_2_tau              = 1/(var->INa_vi_2_al + var->INa_vi_2_bet);
	var->INa_va_ss                 = var->INa_va_al * var->INa_va_tau; // a*tau
	var->INa_vi_1_ss               = var->INa_vi_1_al * var->INa_vi_1_tau;
	var->INa_vi_2_ss               = var->INa_vi_2_al * var->INa_vi_2_tau;

	var->INa_va_tau                 *= p.INa_va_tau_scale;
	var->INa_vi_1_tau               *= p.INa_vi_1_tau_scale;
	var->INa_vi_2_tau               *= p.INa_vi_2_tau_scale;
}

// Update gates and compute current identical to that of LR model, found in lib/Model.c

// End INa ==================================================================//|

// Ito ======================================================================\\|
void set_Ito_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;//-7;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;//-7;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;//-7;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;//-7;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation	
	var->Ito_va_ss          = sigmoid(Vm_ac_ss, 15, -7*p.Ito_va_ss_kscale);
	var->Ito_va_tau			= 0.5*(0.79 + 36.2*exp(-pow((Vm_ac_tau + 40)/45,2))); 
	var->Ito_va_tau         *= p.Ito_va_tau_scale;

	// Voltage inactivation	
	var->Ito_vi_ss          = 1.0/(1.0 + exp((Vm_inac_ss - (-23))/5.3));
	var->Ito_vi_tau         = (8.6 + 62.3*exp(-pow((Vm_inac_tau + 32)/27,2))); 
	var->Ito_vi_tau         *= p.Ito_vi_tau_scale;
	var->Ito_vi_s_tau       = 15.0 + 27.93/(1.0 + exp(0.0696*(Vm_inac_tau - 2.72)));
	var->Ito_vi_Fs          = 0.2/(1+exp((Vm - 35)/-5));
}

void update_gates_Ito_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi              = rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
	s->Ito_vi_s            = rush_larsen(s->Ito_vi_s, var->Ito_vi_ss, var->Ito_vi_s_tau, dt);
}

void compute_Ito_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ito					=  p.gto * s->Ito_va * ((1-var->Ito_vi_Fs)*s->Ito_vi + var->Ito_vi_Fs*s->Ito_vi_s) * (Vm - var->EK);
	var->Ito					*= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
// CRN mWL version ========\\|
void set_ICaL_hAM_CRN_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift-3;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift-3;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift-3;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift-3;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_CRN_mWL_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_CRN_mWL_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// calcium inactivation
	var->ICaL_ci_tau      = 2;	
}

void set_ICaL_hAM_CRN_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss                  = sigmoid(Vm_ss, -10, -0.95*7.45*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)  K = 8 -> 0.95*7.45 is WL MOD
	if (fabs(Vm_tau+10) < 1.0e-10)
	{
		*va_tau             = 4.579/(1.0+exp((Vm_tau+10.0)/-6.24));
	}
	else *va_tau            = (1-exp((Vm_tau+10)/-6.24))/(0.035*(Vm_tau+10)*(1+exp((Vm_tau+10)/-6.24)));
}

void set_ICaL_hAM_CRN_mWL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss                  = exp(-(Vm_ss+28.0)/6.9)/(1.0+exp(-(Vm_ss+28.0)/6.9));
	*vi_tau                 = 9.0/(0.0197*exp(-pow(0.0337,2)*pow((Vm_tau+10),2))+0.02);
}

void compute_ICaL_hAM_CRN_mWL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double ci;
	ci = s->ICaL_ci;

	var->ICaL               = p.gCaL * s->ICaL_va * s->ICaL_vi * ci * (Vm - 55);  // Note, ci = 1 is fully inactivated
	var->ICaL               *= p.GCaL;
}
// End CRN mWL version ====//|

// GB mWL version =========\\|
void set_ICaL_hAM_GB_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift   -10;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift   -10;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift  -10;    // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift  -10;    // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_GB_mWL_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_GB_mWL_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	var->ICaL_vi_s_tau      = 0; // No slow gate, but can cause issues if not set to 0

	// calcium inactivation
	var->ICaL_ci_al       = 1.7;
	var->ICaL_ci_bet      = 11.9e-3;
}

void set_ICaL_hAM_GB_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss                  = sigmoid(Vm_ss, -9.0, -0.95*7.2*kscale);  // V, V1/2, k
	*va_tau                 = 1.0 * (*va_ss) * (1.0 - exp(-(Vm_tau - (-9.0)) / 6.0)) / (0.035 * (Vm_tau - (-9.0)));
}

void set_ICaL_hAM_GB_mWL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss                  = 1.0 / (1.0 + exp((Vm_ss - (-30.0)) / 7.0)) + 0.2 / (1.0 + exp((50.0 - Vm_ss) / 20.0));
	*vi_tau                 = 1.0 / (0.0197 * exp(-pow(0.0337 * (Vm_tau - (-25.0)), 2.0)) + 0.02);
}
// End GB mWL version =====//|

// NG mWL version =========//|
void set_ICaL_hAM_NG_mWL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift - 4;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift  - 4;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift - 4;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift - 4;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_NG_mWL_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_NG_mWL_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// Slow V inac
	var->ICaL_vi_s_tau      = ((0.3323 * exp(-pow(((Vm + 40.0) / 14.2), 2))) + 0.0626)*1e3;
}

void set_ICaL_hAM_NG_mWL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss              = sigmoid(Vm_ss, -9.0, -0.95*0.975*1.2*5.8*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*va_tau             = ((0.0027 * exp(-pow(((Vm_tau + 35.0) / 30.0), 2))) + 0.002)*1e3;
}

void set_ICaL_hAM_NG_mWL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss                  = sigmoid(Vm_ss, -27.4, 7.1*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*vi_tau                 = (0.161 * exp(-pow( ( (Vm_tau + 40.0) / 14.4 ) , 2)) + 0.01)*1e3;
}
// End NG mWL version =====//|

// Full WL version ========\\|
void set_ICaL_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	set_ICaL_hAM_WL_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_hAM_WL_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	var->ICaL_vi_s_tau      = 12424 - 12027*exp(- ((Vm_inac_tau - 13)/83)*((Vm_inac_tau - 13)/83) );

	// calcium inactivation
	var->ICaL_ci_tau      = 50; // consider making faster
	var->ICaL_ci_al       = 5.1;//3*1.7;
	var->ICaL_ci_bet      = 8.33e-3;//0.7*11.9e-3;
}

void set_ICaL_hAM_WL_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
    // Original/as published
	*va_ss                  = sigmoid(Vm_ss, 0.5, -5.967*kscale); 
	*va_tau                 = 7.02 - 2.37*exp(- ((Vm_tau - 14.45)/52.33)*((Vm_tau - 14.45)/52.33) );
}

void set_ICaL_hAM_WL_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	 // Original/as published
    *vi_ss                  = sigmoid(Vm_ss, -18, 3.8*kscale);
	*vi_tau                 = 16.48 - 10.72*exp(- ((Vm_tau - -2.22)/22.64)*((Vm_tau - -2.22)/22.64) );  // v2
}
// End Full WL version ====//|

// Grandi-style ICaL bar, single compartment
void compute_ICaL_hAM_WL_CRN_bar(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double Cai)
{
	double ci;
	ci = s->ICaL_ci;

	var->ICaL_bar_sl        = compute_ICaL_bar_hAM_GB(p, var, Vm, Cai, s->Cao); 

	var->ICaL               = p.pCaL * s->ICaL_va * ((1-p.ICaL_vi_Fs)*s->ICaL_vi + p.ICaL_vi_Fs*s->ICaL_vi_s)  * ci * var->ICaL_bar_sl; 
	var->ICaL               *= p.GCaL;
}

void compute_ICaL_hAM_WL_GB_bar(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double ci, ci_j;
	ci = 1 - s->ICaL_ci; 
	ci_j = 1 - s->ICaL_ci_j;

	var->ICaL_bar_sl        = compute_ICaL_bar_hAM_GB(p, var, Vm, s->Cai_sl, s->Cao);
	var->ICaL_bar_j         = compute_ICaL_bar_hAM_GB(p, var, Vm, s->Cai_j, s->Cao);

	double Op               = s->ICaL_va * ((1-p.ICaL_vi_Fs)*s->ICaL_vi + p.ICaL_vi_Fs*s->ICaL_vi_s);
	var->ICaL_Ca_sl         =  (1 - p.Fjunc_ICaL) * p.pCaL * Op * ci   * var->ICaL_bar_sl;
	var->ICaL_Ca_j          =  (    p.Fjunc_ICaL) * p.pCaL * Op * ci_j * var->ICaL_bar_j;
	var->ICaL_Ca_sl         *= p.GCaL;
	var->ICaL_Ca_j          *= p.GCaL;
	var->ICaL               = var->ICaL_Ca_sl + var->ICaL_Ca_j;

	// To ensure K and Na components don't contribute to homeostasis as not in this model
	var->ICaL_Na_sl = var->ICaL_Na_j = var->ICaL_K = 0.0;
}

// Update gating variables functions
void update_gates_ICaL_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_va              = rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
	s->ICaL_vi              = rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
	s->ICaL_vi_s            = rush_larsen(s->ICaL_vi_s, var->ICaL_vi_ss, var->ICaL_vi_s_tau, dt);
}

void update_gates_ICaL_hAM_WL_CRN_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_ci              = rush_larsen(s->ICaL_ci, 1.0/(1.0+s->Cai_sl/0.00035), var->ICaL_ci_tau, dt);
	s->ICaL_ci_j            = rush_larsen(s->ICaL_ci_j, 1.0/(1.0+s->Cai_j/0.00035), var->ICaL_ci_tau, dt);
}

void update_gates_ICaL_hAM_WL_NG_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_ci              = s->Cai_j / (s->Cai_j + 0.025);
}

void update_gates_ICaL_hAM_WL_GB_ci(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_ci              += dt*(var->ICaL_ci_al * s->Cai_sl * (1.0 - s->ICaL_ci) -  var->ICaL_ci_bet * s->ICaL_ci);
	s->ICaL_ci_j            += dt*(var->ICaL_ci_al * s->Cai_j * (1.0 - s->ICaL_ci_j) - var->ICaL_ci_bet * s->ICaL_ci_j);
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
void set_IKur_hAM_WL_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	// Voltage activation   
	var->IKur_va_ss         = 1.0/(1 + exp((Vm_ac_ss-(-4.2516))/5.61))*4.1503*exp(0.1828*Vm_ac_ss-0.9849);
	if (var->IKur_va_ss > 1.0) var->IKur_va_ss = 1.0;
	var->IKur_va_tau        = (0.009/(1.0 + exp((Vm_ac_tau - (-5.0))/12.0)) + 0.0005)*1000; // ms    // 009 003
	var->IKur_va_tau         *= p.IKur_va_tau_scale;

	// Voltage inactivation
	var->IKur_vi_ss    		= sigmoid(Vm_inac_ss, -7.5, 10*p.IKur_vi_ss_kscale); // Vm, V1/2, k  1/(1+exp((V-V1/2)/k)
	var->IKur_vi_tau         =  (0.59/(1 + exp((Vm_inac_tau - (-60.0))/10.0)) + 3.05)*1000; // ms
	var->IKur_vi_tau         *= p.IKur_vi_tau_scale;
}

void update_gates_IKur_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKur_va  			= rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);
	s->IKur_vi  			= rush_larsen(s->IKur_vi, var->IKur_vi_ss, var->IKur_vi_tau, dt);
}

void compute_IKur_hAM_WL(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKur 				= p.gKur * s->IKur_va * s->IKur_vi * (Vm - var->EK);
	var->IKur				*= p.GKur;
}
// End IKur =================================================================//|

// IK1 ======================================================================\\|
void compute_IK1_hAM_WL_isolated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double ca, cb, cc, cd, ce;
	ca = 0.0955775201140965;
	cb = 0.0071299443560193;
	cc = 0.0000895076773421;
	cd = -0.0000007131973142;
	ce = -0.0000000136553877;

	double Vm_in    = Vm - p.IK1_va_shift;

	var->IK1   		= 4 * (ca + cb*Vm_in + cc*Vm_in*Vm_in + cd*Vm_in*Vm_in*Vm_in + ce*Vm_in*Vm_in*Vm_in*Vm_in);
	var->IK1  		*= p.GK1;
}

void compute_IK1_hAM_WL_intact(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double ca, cb, cc, cd, ce;
	ca = 0.0029122730527995;
	cb = 0.0012982932986851;
	cc = 0.0000351110152611;
	cd = -0.0000009764299296;
	ce = -0.0000000147248251;

	double Vm_in	= Vm - p.IK1_va_shift;

	var->IK1        = 4 * (ca + cb*Vm_in + cc*Vm_in*Vm_in + cd*Vm_in*Vm_in*Vm_in + ce*Vm_in*Vm_in*Vm_in*Vm_in);
	var->IK1        *= p.GK1;
}

// End IK1 =================================================================//|
// NaK, NCX, ICaP, INab, ICab, Ca handling all from CRN, NG or GB models
