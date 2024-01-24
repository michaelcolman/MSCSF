// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Modified ver. of the human atrial cell model=  //
// by Courtemanche, Ramirez and Nattel:  M. Courtemanche, =  //
// R.J. Ramirez and S. Nattel. Am J Physiol 1998; =========  //
// 275(1 Pt 2):H301-21. ===================================  //
// In-code identifier: mCRN ===============================  //
// Presented in Colman et al. 2023 Interface Focus=========  // 
// (In press). doi: 10.1098/rsfs.2023.0041=================  //
// ========================================================  //
// GNU 3 LICENSE TEXT =====================================  //
// COPYRIGHT (C) 2016-2023 MICHAEL A. COLMAN ==============  //
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
// If your model does not inheret from a cell model already included, define all parameters:
void set_parameters_native_mCRN(Cell_parameters *p)
{
	// Capacitance and cell structure
	p->Cm		= 100;		        // pF
	p->Vcell	= 20100.0;	        // um^3
	p->Vcyto    = 0.68*p->Vcell;    // Vi in original model
	p->VjSR     = 0.0048*p->Vcell;  // Vrel in original model
	p->VnSR     = 0.0552*p->Vcell;  // Vup in original model:

	// Concentrations (constant OR initial condition)
	p->Nao		= 140;		// mM
	p->Cao		= 1.8;
	p->Ko		= 5.4;
	p->Nai		= 7.5;//11.75; // This is what it appears to tend towards
	p->Ki		= 138.4;
	p->Cai		= 0.0001024;

	// Current parameters
	p->gNa          = 7.8; // s/mF
	p->gto          = 0.096;
	p->gKur         = 0.0115;
	p->gKr          = 0.00899;//0.0145;
	p->gKs          = 0.052;
	p->gCaL         = 0.34;
	p->gK1          = 0.1;
	p->gKACh        = 0.0045;

	p->gNab         = 1.0e-5;
	p->gCab         = 1.0e-5;

	p->INaK_bar		= 0.6;
	p->INaK_kK		= 1.5;
	p->INaK_kNa		= 10.0;

	p->ICaP_bar		= 0.275;
	p->ICaP_kCa		= 0.5e-3;

	p->INCX_bar 	= 1600;     // pA/pF
	p->INCX_kNao 	= 87.5;     // mM
	p->INCX_kCao	= 1.38;     // mM
	p->INCX_k		= 0.1;
	p->INCX_gamma 	= 0.35;

	// Ca2+ handling
	p->cmdnbar      = 0.045;
	p->trpnbar      = 0.35;
	p->csqnbar      = 10.0;
	p->J_rel_max    = 8.0;  
	p->J_SERCA_max  = 0.0035;  
	p->J_SERCA_kCa  = 6e-4;	
	p->J_leak_max	= 0.0035;
	p->J_leak_kCaSR	= 27.0;
	p->J_jsr_nsr_tau = 180.0; // ms

	// Stimulus parameters (if different to default)
	p->stimduration 	= 2.0;      // ms
	p->stimmag      	= -40.0;    // pA/pF

	// Default celltype
	p->Celltype			= "RA";

	// ISO model
	p->ISO_model        = "Toy_CaSR_load";
}

void update_parameters_integrated_mCRN(Cell_parameters *p)
{
	// All of these are involved in the Ca handling model
	p->LTCC_Ca_bar  		= 1.0;
	p->LTCC_kva1_va2		= 0.6*0.075;          // ms^-1
	p->J_SERCA_max          = 1.45*0.3729;       // uM/ms
	p->J_leak_max           = 1.45*1.41265*1e-5; // ms^-1
	p->INCX_bar             = 0.45*0.9*0.3726;   // (um^3.uM.ms^-1)
	p->ICab_bar             = 0.4*0.6*1.8237*1e-5; // (um^3 . uM . ms^-1)
	p->ICaP_bar             = 0.4*0.9*0.137*1e-2;   // (um^3 . uM . ms^-1)
	p->NLTCC_mean           = 2*12;               // Number of LTCCs per dyad 
}

// Initial conditions
void initial_conditions_native_mCRN(State_variables *s, Cell_parameters p)
{
	// Define ICs for ALL state variables used in the model
	s->Vm      		= -85;
	s->INa_va       = 0.00291;
	s->INa_vi_1     = 0.9791;
	s->INa_vi_2     = 0.9869;
	s->Ito_va       = 0.07;
	s->Ito_vi       = 0.99;
	s->ICaL_va      = 4.757e-6;
	s->ICaL_vi      = 0.9999;
	s->ICaL_ci      = 0.7484;
	s->IKur_va      = 0.06;
	s->IKur_vi      = 0.99;
	s->IKr_va       = 0.000072;
	s->IKs_va       = 0.022846;
	s->IKACh_va     = 1e-5;
	s->cmdn         = 1.856e-3;
	s->trpn         = 7.022e-3;
	s->csqn         = 6.432;
	s->RyRo         = 0;            // "uu" in Varela model
	s->RyRr         = 1.0;          // "vv" in Varela model
	s->RyRi         = 0.9992;       // "w" in Varela model

	s->CajSR        = 1.502;        // CaSR
	s->CanSR        = 1.502;

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
void set_het_mod_mCRN(Cell_parameters *p)
{
	if (                                       p->Het_set_ref == 0)         set_celltype_native_mCRN(p);
	if (p->ISO > 0.0                        && p->ISO_set_ref == 0)         set_modulation_ISO_native_mCRN(p);
	if (strcmp(p->Agent, "none") != 0       && p->Agent_set_ref == 0)       set_modulation_Agent_native_mCRN(p);
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_mCRN(p);
	if (strcmp(p->Mutation, "none") != 0    && p->Mutation_set_ref == 0)    set_modulation_Mutation_native_mCRN(p);
	if (p->ACh > 0.0                        && p->ACh_set_ref == 0)         set_modulation_ACh_mCRN(p);
}

void update_het_and_mod_mCRN_integrated(Cell_parameters *p)
{
    update_celltype_integrated_mCRN(p);
}

// Celltype/Heterogeneity
void set_celltype_native_mCRN(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ================================================================================//|

	// NOTE:: Ensure setting only model-specific conditions here - others defined in Model.c
	if (strcmp(p->Celltype, "RA") == 0); 	// Do nothing for baseline celltype; ensure this is what is set in set_parameters_native_species()
	else if (strcmp(p->Celltype, "LA") == 0)
	{
		p->IKr_rect_vshift      += -34;
		p->IKr_rect_grad_scale  *= (24.4/20.4);
		p->IKr_rect_exp_factor  *= 0.5;
		p->GKr					*= (0.0145/0.00899);
	}
	else if (strcmp(p->Celltype, "PM") == 0)
	{
		//p->GKr                  *= (0.00899/0.0145);
	}
	else if (strcmp(p->Celltype, "BB") == 0)
	{
		p->GCaL                 *= (0.58/0.34);
		//p->GKr                  *= (0.00899/0.0145);
	}
	else if (strcmp(p->Celltype, "PV") == 0)
	{
		p->IK1_VEexp_scale      *= (0.078/0.072);
		p->IK1_den_add_factor   *= 0.66;
		p->IK1_va_shift         += 18; // Note: + 18 as my shift is inside V bracket
		p->IK1_Erev_shift       += 9; // Note to Oleg: term added to give less negative RP

		p->IKr_rect_vshift      += -34;
		p->IKr_rect_grad_scale  *= (24.4/20.4);
		p->IKr_rect_exp_factor  *= 0.5;

		p->GCaL                 *= (0.255/0.34);
		p->GK1                  *= (0.036/0.1);
		p->Gto                  *= (0.07104/0.096);
		p->GKr                  *= (0.022475/0.00899);
		p->GKs                  *= (0.0832/0.052);
		p->GKACh                *= (0.0065/0.0045);
		p->GClb                 *= (0.0055/8.0e-4);
	}
	else 
	{
		printf("ERROR: \"%s\" is not a valid Celltype for the mCRN models. Please check Model.c and Model_mCRN.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

void update_celltype_integrated_mCRN(Cell_parameters *p)
{
    // Basically, here are any het differences in integrated which are not present in native, which result from the diff Ca handling system behaviour
    if (strcmp(p->Celltype, "BB") == 0)
    {
        p->Gup                  *= 0.75;
        p->Gleak                *= 2;
        p->GCab                 *= 0.25;
        p->GCaP                 *= 4;
    }
    else if (strcmp(p->Celltype, "PV") == 0)
    {
        p->Gup                 *= 1.5;
    }
    // no need to error check here, as any invalid celltype will have been caught by native functions
}

// ISO
void set_modulation_ISO_native_mCRN(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
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

	// testing exmaple illustration of model-specific implementation
	// Illustrative ISO model; gain in CICR, little change in APD
	if (strcmp(p->ISO_model, "Toy_CaSR_load") == 0)
	{
		p->GLTCC_kva1_va2   *= (1.0 + p->ISO*(2.5   -1));       // x 2.5 at maximum ISO -> this is GCaL in native and open rate scale in integrated
		p->GKur             *= (1.0 + p->ISO*(2.5   -1));
		p->GKs              *= (1.0 + p->ISO*(2.5   -1));
		p->Gup              *= (1.0 + p->ISO*(2.5   -1));
		//p->GNa              *= (1.0 + p->ISO*(1.25  -1));
		//p->GNCX				*= (1.0 + p->ISO*(1.5   -1));
		//p->GRyR_kCO			*= (1.0 + p->ISO*(2.0   -1));		// this is Grel in native and open rate scale in integrated
	}
	// Add new ISO_models here: else if (strcmp(p->ISO_model, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid ISO model for the mCRN models. Please check Model.c and Model_mCRN.cpp for options\n\n", p->ISO_model);
		exit(1);
	}
}

// Pharmacological modulation
void set_modulation_Agent_native_mCRN(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
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

	if (strcmp(p->Agent, "none") == 0); // do nothing
	// testing exmaple illustration of model-specific implementation
	else if (strcmp(p->Agent, "test") == 0) 
	{
		p->Gto              *= (1.0 + p->Agent_prop*(3.0-1.0));    // At maximal effect, x3, 
		p->IKr_va_ss_kscale *= (1.0 + p->Agent_prop*(0.5-1.0));    // At maximal effect, x0.5
		p->ICaL_vi_ss_shift += p->Agent_prop*5;                    // At maximal effect, + 5
		p->Ito_va_tau_scale *= (1.0 + p->Agent_prop*(2.3-1.0));    // At maximal effect x2.3
		p->Gup              *= (1.0 + p->Agent_prop*(2.5-1.0));    // At maximal effect, x2.5
	}
	// Add new pharmacological agent here: else if (strcmp(p->Agent, "X") == 0) {   }
	else
	{
		printf("ERROR: \"%s\" is not a valid pharmacological agent for the mCRN models. Please check Model.c and Model_mCRN.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

// Remodelling
void set_modulation_Remodelling_native_mCRN(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
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

	if (strcmp(p->Remodelling, "none") == 0); // do nothing
	//else if (strcmp(p->Remodelling, "short_AP_CICR") == 0)
	//{
	//	p->GNa				*= (1.0 + p->Remodelling_prop*(0.9  - 1.0));
	//	p->GCaL				*= (1.0 + p->Remodelling_prop*(0.7  - 1.0));
	//	p->GK1				*= (1.0 + p->Remodelling_prop*(1.5  - 1.0));
	//	p->Gto				*= (1.0 + p->Remodelling_prop*(0.54 - 1.0));
	//	p->GKr				*= (1.0 + p->Remodelling_prop*(1.5  - 1.0));
	//	p->GKs				*= (1.0 + p->Remodelling_prop*(1.6  - 1.0));
	//	p->GNCX				*= (1.0 + p->Remodelling_prop*(0.75 - 1.0));
	//	p->GKACh			*= (1.0 + p->Remodelling_prop*(1.8  - 1.0));
	//	p->Gup				*= (1.0 + p->Remodelling_prop*(2.0  - 1.0));
	//	p->Gleak            *= (1.0 + p->Remodelling_prop*(0.5  - 1.0));
	//	p->GRyR_kCO			*= (1.0 + p->Remodelling_prop*(2.0  - 1.0));
	//}
    else if (strcmp(p->Remodelling, "AF") == 0)
    {
        p->GNa                  *= (1.0 + p->Remodelling_prop*(0.5  -1));
        p->GCaL                 *= (1.0 + p->Remodelling_prop*(0.35  -1));
        p->Gto                  *= (1.0 + p->Remodelling_prop*(0.2  -1));
        p->GKur                 *= (1.0 + p->Remodelling_prop*(0.75  -1));
        p->GKr                  *= (1.0 + p->Remodelling_prop*(3  -1));
        p->GKs                  *= (1.0 + p->Remodelling_prop*(3  -1));
        p->GKACh                *= (1.0 + p->Remodelling_prop*(1.8 - 1.0));
    }
    else if (strcmp(p->Remodelling, "test") == 0)
    {
        p->Gto              *= (1.0 + p->Remodelling_prop*(3.0-1.0));    // At maximal effect, x3, 
        p->IKr_va_ss_kscale *= (1.0 + p->Remodelling_prop*(0.5-1.0));    // At maximal effect, x0.5
        p->ICaL_vi_ss_shift += p->Remodelling_prop*5;                    // At maximal effect, + 5
        p->Ito_va_tau_scale *= (1.0 + p->Remodelling_prop*(2.3-1.0));    // At maximal effect x2.3
        p->Gup              *= (1.0 + p->Remodelling_prop*(2.5-1.0));    // At maximal effect, x2.5
    }
    // Add new Remodelling here: else if (strcmp(p->Remodelling, "X") == 0) {   }
    else
    {
        printf("ERROR: \"%s\" is not a valid Remodelling model for the mCRN models. Please check Model.c and Model_mCRN.cpp for options\n\n", p->Remodelling);
        exit(1);
    }
}

// Mutation
void set_modulation_Mutation_native_mCRN(Cell_parameters *p)
{
    // ================================================================================\\|
    // Can modify parameters (e.g. conductance) directly, or scaling factors.
    // (conductance denoted "g", scale factor "G"; current = f(g*G))
    // In general, use the scale factors as these are multiplicative/additive throughout
    // the code. 
    // Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // ===========================================
    // NOTE:: Ensure setting only model-specific conditions here - global defined in Model.c
    // ================================================================================//|

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
        printf("ERROR: \"%s\" is not a valid Mutation model for the mCRN models. Please check Model.c and Model_mCRN.cpp for options\n\n", p->Mutation);
        exit(1);
    }
}
// ACh ===================================\\|
void set_modulation_ACh_mCRN(Cell_parameters *p)
{
    // ================================================================================\\|
    // Can modify parameters (e.g. conductance) directly, or scaling factors.
    // (conductance denoted "g", scale factor "G"; current = f(g*G))
    // In general, use the scale factors as these are multiplicative/additive throughout
    // the code. 
    // Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
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

    // testing exmaple illustration of model-specific implementation

    // NOTE: if the model does not have IKACh (or If) in, then modifying them won't make a difference!
    /*if (strcmp(p->ACh_model, "default") == 0)
      {
    // Put actual things in here
    }          
    else*/ if (strcmp(p->ACh_model, "test") == 0)  // just an example of how to implement: NOT an actual implementation
    {
        p->gKACh    		= p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) );
        // OR
        p->gKACh			*= (1.0 + p->ACh*(2.5-1.0));
        p->If_va_ss_shift 	+= -7.2*p->ACh;
        p->Gto              *= (1.0 + p->ACh*(3.0-1.0));    // At maximal effect, x3,
    }
    else
    {
        printf("ERROR: \"%s\" is not a valid ACh for the mCRN model. Please check Model_mCRN.cpp for options\n\n", p->ACh_model);
        exit(1);
    }
}
// End ACh ===============================//|
// End heterogeneity and modulation =============================================================//|

// Compute model functions ======================================================================\\|
// Your model may have more or fewer currents than this template - just follow the procedure and add/delete as appropriate
void compute_model_mCRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	compute_reversal_potentials(p, var, s);
	set_gate_rates_mCRN_native(p, var, Vm, s->Cai);
	update_gating_variables_mCRN_native(p, var, s, Vm, dt);
	compute_Itot_mCRN_native(p, var, s, Vm);
	comp_homeostasis_mCRN(p, var, s, Vm, dt);
}

void compute_model_mCRN_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    compute_reversal_potentials(p, var, s);
    set_gate_rates_mCRN_native(p, var, Vm, s->Cai);
    update_gating_variables_mCRN_native(p, var, s, Vm, dt);
    compute_Itot_mCRN_integrated(p, var, s, Vm);
    // NO homeostasis here, as done in Ca handling model
    // Can add a function which does K+ and Na+ cycling if required
}

void set_gate_rates_mCRN_native(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
    // Call only currents you need
    set_INa_LR_rates(p, var, Vm);				
    set_Ito_mCRN_rates(p, var, Vm);
    set_IKur_mCRN_rates(p, var, Vm);
    set_IKs_mCRN_rates(p, var, Vm);
    set_IKr_mCRN_rates(p, var, Vm);
    set_IKACh_mCRN_rates(p, var, Vm);
    set_IK1_mCRN_variables(p, var, Vm);
    set_ICaL_mCRN_rates(p, var, Vm, Cai);
}

void update_gating_variables_mCRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
    update_gates_INa_LR(p, var, s, Vm, dt); 	
    update_gates_IKs_mCRN(p, var, s, Vm, dt);
    update_gates_IKr_mCRN(p, var, s, Vm, dt);
    update_gates_IKACh_mCRN(p, var, s, Vm, dt);
    update_gates_ICaL_mCRN(p, var, s, Vm, dt);
    update_gates_Ito_mCRN(p, var, s, Vm, dt);
    update_gates_IKur_mCRN(p, var, s, Vm, dt);
}

void compute_Itot_mCRN_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    var->Itot   = 0;

    compute_INa_LR(p, var, s, Vm);					
    compute_Ito_mCRN(p, var, s, Vm);
    compute_ICaL_mCRN(p, var, s, Vm);
    compute_IKur_mCRN(p, var, s, Vm);
    compute_IKr_mCRN(p, var, s, Vm);
    compute_IKs_mCRN(p, var, s, Vm);
    compute_IKACh_mCRN(p, var, s, Vm);
    compute_IK1_mCRN(p, var, s, Vm);
    compute_INCX_mCRN(p, var, s, Vm);
    compute_INaK_mCRN(p, var, s, Vm);
    compute_ICaP_mCRN(p, var, s, Vm);
    compute_INab_mCRN(p, var, s, Vm);
    compute_ICab_mCRN(p, var, s, Vm);

    // Add all the computed currents here
    var->Itot   = var->INa + var->Ito + var->ICaL + var->IKur + var->IKr + var->IKs + var->IKACh + var->IK1 + var->INCX + var->INaK + var->ICaP + var->INab + var->ICab + var->IClb;
}

void compute_Itot_mCRN_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
    var->Itot   = 0;

    compute_INa_LR(p, var, s, Vm);
    compute_Ito_mCRN(p, var, s, Vm);
    compute_IKur_mCRN(p, var, s, Vm);
    compute_IKr_mCRN(p, var, s, Vm);
    compute_IKs_mCRN(p, var, s, Vm);
    compute_IKACh_mCRN(p, var, s, Vm);
    compute_IK1_mCRN(p, var, s, Vm);
    compute_INaK_mCRN(p, var, s, Vm);
    compute_INab_mCRN(p, var, s, Vm);

    // Add all the non-Ca computed currents here
    var->Itot   = var->INa + var->Ito + var->IKur + var->IKr + var->IKs + var->IKACh + var->IK1 + var->INaK + var->INab + var->IClb;
    // Then add the Ca currents, computed in CRU
    var->Itot   += var->INCX + var->ICaP + var->ICab + var->ICaL;   // Adds Ca currents computed in CRU
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// INa ======================================================================\\|
// Identical to that of LR model, found in lib/Model.c 
// End INa ==================================================================//|

// Ito ======================================================================\\|
void set_Ito_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
    // First, assign local voltages for each gate and type such that shifts can be applied by MODIFIERS
    //e.g.:
    double Vm_ac_ss         = Vm - p.Ito_va_ss_shift;   // Voltage modified by shift applied to activation steady state
    double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
    double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;  // Voltage modified by shift applied to activation time constant (and rates as they only define tau)
    double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant
    // Be sure to use the relevant local voltage in the correct functions

    // Voltage activation (va), voltage inactivation (vi)
    var->Ito_va_ss		= pow((1.0+exp((Vm_ac_ss - (12.0))/(-11.5*p.Ito_va_ss_kscale))),(-1.0/3.0));
    var->Ito_va_tau		= 0.4/(exp((Vm_ac_tau -15.0)/20.0)); 

    var->Ito_vi_ss                    = 1.0/(1.0+exp((Vm_inac_ss - (-31.0))/(6.45*p.Ito_vi_ss_kscale))); // Steady-state
    var->Ito_vi_al                    = 1.0/(1.2+exp((Vm_inac_tau  - (-95.2))/5.85));
	var->Ito_vi_bet                   = 1.0/(9.54+exp((Vm_inac_tau  - 48.0)/-20.0));
	var->Ito_vi_tau                   = 1.0/(var->Ito_vi_al +  var->Ito_vi_bet);

	// and don't forget to scale the time-constants!
	var->Ito_va_tau			*= p.Ito_va_tau_scale;
	var->Ito_vi_tau			*= p.Ito_vi_tau_scale;
}

void update_gates_Ito_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// Update the gates, using rush_larsen, forward Euler, or other
	s->Ito_va              = rush_larsen(s->Ito_va, var->Ito_va_ss, var->Ito_va_tau, dt);
	s->Ito_vi              = rush_larsen(s->Ito_vi, var->Ito_vi_ss, var->Ito_vi_tau, dt);
}

void compute_Ito_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	// And finally, compute the actual current
	var->Ito				= p.gto * pow(s->Ito_va, 3) * s->Ito_vi * (Vm - var->EK);
	var->Ito				*= p.Gto;
}
// End Ito ==================================================================//|

// ICaL =====================================================================\\|
// ICaL can be written exactly as above for Ito. However, if there are intentions
// to integrate the new model with the "integrated" calcium handling system
// e.g. for spatial cell models or spontaneous release functions, please 
// follow the procedure below; voltage-dependent gates have their own
// functions so can be called elsewhere (i.e. from CRU)
void set_ICaL_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm, double Cai)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	set_ICaL_mCRN_va_rates(p, &var->ICaL_va_ss, &var->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->ICaL_va_tau        *= p.ICaL_va_tau_scale;

	set_ICaL_mCRN_vi_rates(p, &var->ICaL_vi_ss, &var->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->ICaL_vi_tau        *= p.ICaL_vi_tau_scale;

	// calcium inactivation
	var->ICaL_ci_ss		= 0.29+0.8/(1.0+exp((Cai -1.2e-4)/0.00006));
	var->ICaL_ci_tau	= 2.0;
}

void set_ICaL_mCRN_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss        = sigmoid(Vm_ss, -2.0, -5.0*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	if (fabs(Vm_tau - 10) < 1e-10)
	{
		*va_tau                     = 0.763;
	}
	else *va_tau                    = (1.0/(1.0+exp((Vm_tau +10.0)/-6.24)))*(1.0-exp((Vm_tau +10.0)/-6.24))/(0.035*(Vm_tau +10.0));
}

void set_ICaL_mCRN_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss      = sigmoid(Vm_ss, -34, 6.3*kscale);  // V, V1/2, k 1/(1+exp((V-V1/2)/k)
	*vi_tau  	= 400.0/(1+4.5*exp(-0.0007*pow(Vm_tau -9,2.0)));
}

// Everything else follows the same format
void update_gates_ICaL_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->ICaL_va       		= rush_larsen(s->ICaL_va, var->ICaL_va_ss, var->ICaL_va_tau, dt);
	s->ICaL_vi    			= rush_larsen(s->ICaL_vi, var->ICaL_vi_ss, var->ICaL_vi_tau, dt);
	s->ICaL_ci    			= rush_larsen(s->ICaL_ci, var->ICaL_ci_ss, var->ICaL_ci_tau, dt);
}

void compute_ICaL_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaL   			= p.gCaL * s->ICaL_va * s->ICaL_vi * s->ICaL_ci * (Vm - 60);
	var->ICaL 				*= p.GCaL;
}
// End ICaL =================================================================//|

// IKur =====================================================================\\|
void set_IKur_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;   // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;  // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;  // Voltage modified by shift applied to inactivation time constant

	var->IKur_va_ss			= pow((1+exp((Vm_ac_ss  -(-2.81))/(-9.49*p.IKur_va_ss_kscale))),(-1.0/3.0));
	var->IKur_va_al			= 1.47/(exp((Vm_ac_tau  - (-33.20))/-30.63)+exp((Vm_ac_tau  -27.6)/-30.65));	
	var->IKur_va_bet		= 0.42/(exp((Vm_ac_tau  -(-26.64))/2.49)+exp((Vm_ac_tau  +44.41)/20.36));
	var->IKur_va_tau        = 1.0/(var->IKur_va_al +  var->IKur_va_bet);
	var->IKur_va_tau		*= p.IKur_va_tau_scale;

	var->IKur_vi_ss         = sigmoid(Vm_inac_ss, 99.45, 27.48*p.IKur_vi_ss_kscale);
	var->IKur_vi_al         = 1.0/(21.0+exp((Vm_inac_tau  -185.0)/-28.0)); 
	var->IKur_vi_bet        = exp((Vm_inac_tau  -158.0)/16.0);
	var->IKur_vi_tau        = 1.0/(var->IKur_vi_al +  var->IKur_vi_bet);
	var->IKur_vi_tau        *= p.IKur_vi_tau_scale;

	// Dynamic conductance factor
	var->IKur_dynamic_g		= 1.0+3.0/(1.0+exp((Vm -14.0)/-6.0));
}

void update_gates_IKur_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKur_va                      = rush_larsen(s->IKur_va, var->IKur_va_ss, var->IKur_va_tau, dt);
	s->IKur_vi                      = rush_larsen(s->IKur_vi, var->IKur_vi_ss, var->IKur_vi_tau, dt);
}

void compute_IKur_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKur 				= p.gKur * var->IKur_dynamic_g * pow(s->IKur_va, 3) * s->IKur_vi * (Vm - var->EK);
	var->IKur				*= p.GKur;
}
// End IKur =================================================================//|

// IKr ======================================================================\\|
void set_IKr_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKr_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKr_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKr_va_ss			= sigmoid(Vm_ac_ss, -7.645, -5.377*p.IKr_va_ss_kscale);
	var->IKr_va_al			= 0.04*(Vm_ac_tau  -248.0)/(1.0-exp((Vm_ac_tau  -248.0)/-28.0));
	var->IKr_va_bet			= 0.028*(Vm_ac_tau  +163.0)/(exp((Vm_ac_tau  +163.0)/21.0)-1.0);
	var->IKr_va_tau         = 1.0/(var->IKr_va_al +  var->IKr_va_bet);
	var->IKr_va_tau			*= p.IKr_va_tau_scale;

	// Time-independant gate
	// controlled by het set parameters
	var->IKr_vi_ti        = 0.6+1.0/(0.5 + p.IKr_rect_exp_factor*exp((Vm - (26.0 + p.IKr_rect_vshift))/(20.4*p.IKr_rect_grad_scale)));

	// explict het
	/*if (strcmp(p.Celltype, "LA") == 0 || strcmp(p.Celltype, "PV") == 0)
	  {
	  var->IKr_vi_ti		= 0.6+1.0/(0.5+0.5*exp((Vm+8.0)/24.4));
	  }
	  else // RA/BB
	  {
	  var->IKr_vi_ti 		=  0.6+1.0/(0.5+exp((Vm-26.0)/20.4));
	  }*/
}

void update_gates_IKr_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKr_va               = rush_larsen(s->IKr_va, var->IKr_va_ss, var->IKr_va_tau, dt);
}

void compute_IKr_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKr 				= p.gKr * s->IKr_va * var->IKr_vi_ti * (Vm - var->EK);
	var->IKr 				*= p.GKr;
}
// End IKr ==================================================================//|

// IKs ======================================================================\\|
void set_IKs_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKs_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKs_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKs_va_ss  		= sqrt(1.0/(1.0+exp((Vm_ac_ss - (13.0))/(-12.0*p.IKs_va_ss_kscale))));
	if (fabs(Vm_ac_tau+28.5) < 1e-10) /* denominator = 0 */
	{
		var->IKs_va_al              = 0.00115;
		var->IKs_va_bet             = 0.000759;
	}
	else
	{
		var->IKs_va_al           = 0.00001*(Vm_ac_tau +28.5)/(1.0-exp((Vm_ac_tau +28.5)/-115.0));
		var->IKs_va_bet          = 0.00023*(Vm_ac_tau +28.5)/(exp((Vm_ac_tau +28.5)/3.3)-1.0);
	}
	var->IKs_va_tau              = 1.0/(var->IKs_va_al + var->IKs_va_bet);
	var->IKs_va_tau				 *= p.IKs_va_tau_scale;
}

void update_gates_IKs_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKs_va               = rush_larsen(s->IKs_va, var->IKs_va_ss, var->IKs_va_tau, dt);
}

void compute_IKs_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKs		= p.gKs * s->IKs_va*s->IKs_va * (Vm - var->EK);
	var->IKs		*= p.GKs;
}
// End IKs ==================================================================//|

// IKACh =====================================================================\\|
void set_IKACh_mCRN_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKACh_va_ss_shift;   // Voltage modified by shift applied to activation steady state
	double Vm_ac_tau        = Vm - p.IKACh_va_tau_shift;  // Voltage modified by shift applied to activation time constant

	var->IKACh_va_ss        = 1.0/(1.0+exp((-93.0- Vm_ac_ss)/(-15.2*p.IKACh_va_ss_kscale)));
	var->IKACh_va_tau       = 360.0+130.0*(1.0-exp(-(Vm_ac_tau +130.0)/50.0));
	var->IKACh_va_tau       *= p.IKACh_va_tau_scale;

	var->IKACh_v_ti			= 1.0/(0.1+exp(0.078*(Vm - var->EK -65.0)));
}

void update_gates_IKACh_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->IKACh_va               = rush_larsen(s->IKACh_va, var->IKACh_va_ss, var->IKACh_va_tau, dt);
}

void compute_IKACh_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IKACh        = p.gKACh * s->IKACh_va * var->IKACh_v_ti * (Vm - var->EK);
	var->IKACh        *= p.GKACh;
}
// End IKACh ==================================================================//|

// IK ======================================================================\\|
void set_IK1_mCRN_variables(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_in        = Vm - p.IK1_va_shift;
	var->IK1_va_ti      = (p.IK1_den_add_factor*1.0 + exp(p.IK1_VEexp_scale*0.072*(Vm_in-(var->EK))));
}

void compute_IK1_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IK1            = p.gK1 *  (Vm - (var->EK + p.IK1_Erev_shift))/var->IK1_va_ti; // NOTE!!!
	var->IK1            *= p.GK1;
}
// End IK1 ==================================================================//|

// Ca2+ handling, background and pump currents ==============================\\|
// INCX
void compute_INCX_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double Cai          = s->Cai;
	double Cao          = s->Cao;
	double Nai          = s->Nai;
	double Nao          = s->Nao;

	var->INCX           = p.INCX_bar *( exp(p.INCX_gamma* p.FoRT*Vm) * Nai*Nai*Nai * Cao - exp((p.INCX_gamma-1)* p.FoRT*Vm) * Nao*Nao*Nao*Cai)      \
						  / ((pow(p.INCX_kNao,3) + Nao*Nao*Nao)*(p.INCX_kCao + Cao) * (1 + p.INCX_k * exp((p.INCX_gamma-1)* p.FoRT*Vm)));
	var->INCX           *= p.GNCX;

}

// INaK
void compute_INaK_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	double sigma        = (exp(s->Nao/67.3)-1.0)/7.0;
	double FNaK         = pow(1.0+0.1245*exp(-0.1*p.F*Vm/(p.R*p.T))+0.0365*sigma*exp(-Vm*p.FoRT ), -1.0);
	var->INaK			= p.INaK_bar * FNaK * p.Ko / (1.0 + pow(p.INaK_kNa / s->Nai, 4.0)) / (p.Ko + p.INaK_kK);
	var->INaK           *= p.GNaK;
}

// ICaP
void compute_ICaP_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICaP           = (p.ICaP_bar * s->Cai) / (s->Cai + p.ICaP_kCa);
	var->ICaP           *= p.GCaP;
}

// INab
void compute_INab_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->INab           = p.gNab * (Vm - var->ENa);
	var->INab           *= p.GNab;
}

// ICab
void compute_ICab_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->ICab           = p.gCab * (Vm - var->ECa);
	var->ICab           *= p.GCab;
}

/*void compute_IClb_Varela_dAM(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->IClb                         = 8.0e-4 * (Vm - ((p.R * p.T)/p.F)*log(29.26/132)); // static || cli not updated
	var->IClb	*= p.GClb;
}*/
// End Ca2+ handling currents ===============================================//|

// Homeostasis ==============================================================\\|
void comp_homeostasis_mCRN(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	// All calculated currents can be accessed in here through var e.g. var->ICaL, var->INCX or var->INCX_sl/j etc

	// Differentials
	var->dNai 		= p.Cm*(-3*var->INaK - var->INCX - var->INab - var->INa)/(p.F*p.Vcyto);
	var->dKi		= p.Cm*(2*var->INaK- var->IK1 - var->Ito - var->IKur - var->IKr - var->IKs - var->IKACh)/(p.F*p.Vcyto);

	// Buffering
	double cmdndot              = 200.0 *s->Cai   *(1.0-s->cmdn /p.cmdnbar) - 0.476*s->cmdn /p.cmdnbar;  // ms
	double trpndot              = 78.4  *s->Cai   *(1.0-s->trpn /p.trpnbar) - 0.392*s->trpn /p.trpnbar;
	double csqndot              = 0.48  *s->CajSR *(1.0-s->csqn /p.csqnbar) - 0.400*s->csqn /p.csqnbar;

	s->cmdn                 += dt * cmdndot * p.cmdnbar; // mM 
	s->trpn                 += dt * trpndot * p.trpnbar;
	s->csqn                 += dt * csqndot * p.csqnbar;

	// Jrel / RyRs
	var->J_rel 				= p.J_rel_max * s->RyRo*s->RyRo * s->RyRr * s->RyRi * (s->CajSR - s->Cai);

	var->Cai_Fn				= 1e3*(1e-15*p.VjSR * var->J_rel - 1e-15*p.Cm*(0.5*var->ICaL - 0.2*var->INCX)/(2.0*p.F));

	s->RyRo                 = rush_larsen(s->RyRo, 1.0/(1.0+exp(-(var->Cai_Fn-3.4175e-13 )/13.67e-16)) , 11.2, dt);
	s->RyRr                 = rush_larsen(s->RyRr, 1.0-1.0/(1.0+exp(-(var->Cai_Fn-6.835e-14)/13.67e-16)) , 1.91+2.09/(1.0+exp(-(var->Cai_Fn-3.4175e-13)/13.67e-16)), dt);
	if (fabs(Vm -7.9) < 1e-10)
	{
		s->RyRi             = rush_larsen(s->RyRi, 1.0-1.0/(1.0+exp(-(Vm -40.0)/17.0)) , 0.9231, dt);
	}
	else s->RyRi            = rush_larsen(s->RyRi, 1.0-1.0/(1.0+exp(-(Vm -40.0)/17.0)) , 6.0*(1.0-exp(-(Vm -7.9)/5.0))/(1.0+0.3*exp(-(Vm -7.9)/5.0))/(Vm -7.9), dt);

	// SERCA and Jleak
	var->J_SERCA 	= p.J_SERCA_max / (1.0+p.J_SERCA_kCa/s->Cai);
	var->J_leak		= p.J_leak_max * s->CanSR / p.J_leak_kCaSR;
	var->J_SERCA	*= p.Gup;
	var->J_leak		*= p.Gleak;

	// SR compart transfer
	var->J_jsr_nsr          = (s->CanSR - s->CajSR)/p.J_jsr_nsr_tau;

	// Update concentrations
	var->dCai 			= p.Cm*(2.0*var->INCX - var->ICaP - var->ICaL - var->ICab)/(2.0*p.F*p.Vcyto) + (var->J_leak - var->J_SERCA)*p.VnSR/p.Vcyto + var->J_rel*p.VjSR/p.Vcyto - p.trpnbar*trpndot - p.cmdnbar*cmdndot;
	s->CanSR   			+= dt* (var->J_SERCA - var->J_jsr_nsr*p.VjSR/p.VnSR - var->J_leak);
	s->CajSR  			+= dt* (var->J_jsr_nsr - var->J_rel - 31.0*csqndot);
	s->Cai   			+= dt* var->dCai;
	s->Nai  			+= dt* var->dNai;
	s->Ki  				+= dt* var->dKi;
}
// End Homeostasis ==========================================================//|
// End Current formulations =====================================================================//|
