// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Implementation of a novel hybrid-minimal ====  //
// model for integration with biophysically detailed ======  //
// spatial Ca2+ handling, originally presented with this ==  //
// framework. No additional references required. ==========  //
// In-code identifier: minimal ============================  //
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
void set_parameters_native_minimal(Cell_parameters *p)
{
	// Ion current conductances
	p->gIp0d			= 16;		// s/mF
	p->gIp1r			= 0.1652;	// s/mF
	p->gIp2d			= 0.75;		// s/mF
	p->gIp2r			= 0.017622;	// s/mF
	p->gIp3r			= 0.05;		// s/mF
	p->gIp4r			= 0.3; 		// s/mF

	// (and biphys versions, just for outputs)
	p->gNa				= p->gIp0d;
	p->gto				= p->gIp1r;
	p->gKur				= p->gIp2r;
	p->gCaL				= p->gIp2d;
	p->gKr				= p->gIp3r;
	p->gKur				= p->gIp4r;

	// Default celltype
	p->Celltype			= "EPI";

	// Default ISO model
	p->ISO_model        = "Toy_methods_demonstration";
    p->ACh_model        = "test";
}

// Initial conditions
void initial_conditions_native_minimal(State_variables *s, Cell_parameters p)
{
	s->Vm				= -85;
	s->Ip0d_va			= 0.001231;
	s->Ip0d_vi_1		= 0.992842;
	s->Ip0d_vi_2		= 0.988864;
	s->Ip1r_va			= 0.0;
	s->Ip1r_vi			= 0.5;
	s->Ip2d_va			= 0;
	s->Ip2d_vi			= 1;
	s->Ip2r_va			= 0;
	s->Ip2r_vi			= 1;
	s->Ip3r_va			= 0.000175;

	// Assign state from param (even if constant)
	s->Nai              = p.Nai;
	s->Nao              = p.Nao;
	s->Ki               = p.Ki;
	s->Ko               = p.Ko;
	s->Cao              = p.Cao;
}


// Het and modulation functions ===================================\\|
// Parent function
void set_het_mod_minimal(Cell_parameters *p)
{
	if (									   p->Het_set_ref == 0)			set_celltype_native_minimal(p);
	if (p->ISO > 0.0 						&& p->ISO_set_ref == 0) 		set_modulation_ISO_native_minimal(p);
	if (strcmp(p->Agent, "none") != 0 		&& p->Agent_set_ref == 0) 		set_modulation_Agent_native_minimal(p); 
	if (strcmp(p->Remodelling, "none") != 0 && p->Remodelling_set_ref == 0) set_modulation_Remodelling_native_minimal(p);
	if (strcmp(p->Mutation, "none") != 0	&& p->Mutation_set_ref == 0)	set_modulation_Mutation_native_minimal(p);
	if (p->ACh > 0.0						&& p->ACh_set_ref == 0)			set_modulation_ACh_minimal(p);
}

// Heterogeneity
void set_celltype_native_minimal(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired
	if (strcmp(p->Celltype, "EPI") == 0)
	{
		p->Gto			*= 0.5;
		p->GKur			*= 0.2;
	}
	else if (strcmp(p->Celltype, "M") == 0)
	{
		p->GKur   		*= 0.2;
		p->Gto   		*= 0.5;
		p->GKr  		*= 0.6;
		p->GCaL 		*= 1.1;
	}
	else if (strcmp(p->Celltype, "ENDO") == 0)
	{
		p->GKur 		*= 0.2;
		p->Gto  		*= 0.5;
		p->GKr  		*= 0.9;
		p->Gup  		*= 1.1;
		p->Gleak  		*= 1.1;
	}
	else if (strcmp(p->Celltype, "PK") == 0)
	{
		p->GKur 		*= 0.2;
		p->Gto  		*= 0.5;
		p->GKr  		*= 0.9;
		p->Gup  		*= 1.1;
		p->Gleak  		*= 1.1;
        p->GCaL         *= 1.5;
        p->GNa          *= 2;
	}
	else if (strcmp(p->Celltype, "RA") == 0)
	{
		p->Gto                  *= 1.75;
		p->GKur                 *= 1.1;//1.2;
		p->Ito_va_ss_shift     	*= -5;
		p->GK1                  *= 0.5;
		p->GKr                  *= 0.25;
		p->GKs                  *= 0.25;
		p->GCaL      			*= 0.75;
		p->ICaL_va_ss_shift     += 2.5;
		p->Gup                 	*= 1.7;
		p->Gleak               	*= 1.3;
        p->IK1_Erev_shift       += 5;
	}
    else if (strcmp(p->Celltype, "RA_small_mammal") == 0)
    {
        p->Gto                  *= 1.75;
        p->GKur                 *= 1.5;
        p->Ito_va_ss_shift      *= -5;
        p->GK1                  *= 0.75;
        p->GKs                  *= 0.5;
        p->GCaL                 *= 0.75;
        p->GLTCC_kva1_va2       *= 0.7;
        p->ICaL_va_ss_shift     += 2.5;
        p->Gup                  *= 3.75;
        p->Gleak                *= 1.3;
    }
	else if (strcmp(p->Celltype, "RA_reduced_INa_NCX") == 0)
	{
		p->GNa                  *= 0.3;
		p->Gto                  *= 1.75;
		p->GKur                 *= 1.1;//1.2;
		p->Ito_va_ss_shift     	*= -5;
		p->GK1                  *= 0.5;
		p->GKr                  *= 0.25;
		p->GKs                  *= 0.25;
		p->GCaL      			*= 0.75;
		p->ICaL_va_ss_shift     += 2.5;
		p->Gup                 	*= 1.7;
		p->Gleak               	*= 1.3;
        p->IK1_Erev_shift       += 5;
        p->GNCX                 *= 0.5;
	}
    else if (strcmp(p->Celltype, "EPI_CCS") == 0)
    {
        p->Gto          *= 0.5;
        p->GKur         *= 0.2;
        p->GNa          *= 0.4;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
    }
    else if (strcmp(p->Celltype, "M_CCS") == 0)
    {
        p->GKur         *= 0.2;
        p->Gto          *= 0.5;
        p->GKr          *= 0.6;
        p->GCaL         *= 1.1;
        p->GNa          *= 0.4;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
    }
    else if (strcmp(p->Celltype, "ENDO_CCS") == 0)
    {
        p->GKur         *= 0.2;
        p->Gto          *= 0.5;
        p->GKr          *= 0.9;
        p->Gup          *= 1.1;
        p->Gleak        *= 1.1;
        p->GNa          *= 0.4;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
    }
    else if (strcmp(p->Celltype, "PK_CCS") == 0)
    {
        p->GKur         *= 0.2;
        p->Gto          *= 0.5;
        p->GKr          *= 0.9;
        p->Gup          *= 1.1;
        p->Gleak        *= 1.1;
        p->GCaL         *= 1.5;
        p->GNa          *= 2;
        p->GNa          *= 0.4;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
    }
    else if (strcmp(p->Celltype, "RA_CCS") == 0)
    {
        p->Gto                  *= 1.75;
        p->GKur                 *= 1.1;//1.2;
        p->Ito_va_ss_shift      *= -5;
        p->GK1                  *= 0.5;
        p->GKr                  *= 0.25;
        p->GKs                  *= 0.25;
        p->GCaL                 *= 0.75;
        p->ICaL_va_ss_shift     += 2.5;
        p->Gup                  *= 1.7;
        p->Gleak                *= 1.3;
        p->IK1_Erev_shift       += 5;
        p->GNa          *= 0.4;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
    }
    else if (strcmp(p->Celltype, "Pace_CCS") == 0)
    {
        p->GNa                  *= 0.25;
        p->GNa          *= 0.4*1.75*1.5*1.35;
        p->GCaL         *= 3;
        p->GKr          *= 0.3;
        p->GKs          *= 0.3;
        p->GK1          *= 0.5;
        p->IK1_Erev_shift       += 20;
    }
	else
	{
		printf("ERROR: \"%s\" is not a valid Celltype for the minimal model. Please check Model_minimal.cpp for options\n\n", p->Celltype);
		exit(1);
	}
}

void set_modulation_ISO_native_minimal(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired

	if (strcmp(p->ISO_model, "Toy_methods_demonstration") == 0)
	{
		//p->gIp2d		*= (1.0 + p->ISO*2); // etc, so ISO = 0 gives control, ISO = 1 gives full difference
		if (strcmp(p->Celltype, "EPI") == 0 || strcmp(p->Celltype, "M") == 0 || strcmp(p->Celltype, "ENDO") == 0)
		{
			p->GKur				*= (1.0 + p->ISO*(2.0-1));
			p->GKr				*= (1.0 + p->ISO*(2.0-1));
			p->Gup				*= (1.0 + p->ISO*(2.0-1));
			p->GLTCC_kva1_va2	*= (1.0 + p->ISO*(1.5-1)); 	// this will be passed to ICaL_scale for native models
		}
		else if (strcmp(p->Celltype, "RA") == 0)
		{
			p->GKur				*= (1.0 + p->ISO*(1.5-1));
			p->Gto				*= (1.0 + p->ISO*(1.5-1));
			p->Gup				*= (1.0 + p->ISO*(1.2-1));
			p->GLTCC_kva1_va2	*= (1.0 + p->ISO*(2.0-1));
		}
	}
	else 
	{
		printf("ERROR: \"%s\" is not a valid ISO model for \"%s\". Please check Model_minimal.cpp for options\n\n", p->ISO_model, p->Model);
		exit(1);
	}
}


void set_modulation_Agent_native_minimal(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired
	if (strcmp(p->Agent, "none") == 0); // do nothing
	else if (strcmp(p->Agent, "test") == 0) // model specific option
	{
		p->gIp2d		*= (1.0 + p->Agent_prop*(2 - 1));		
	} 
	else
	{
		printf("ERROR: \"%s\" is not a valid Pharmacological Agent for the minimal model. Please check Model_minimal.cpp for options\n\n", p->Agent);
		exit(1);
	}
}

void set_modulation_Remodelling_native_minimal(Cell_parameters *p)
{
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired
	if (strcmp(p->Remodelling, "none") == 0); // do nothing
	else if (strcmp(p->Remodelling, "test") == 0) // model specific option
	{
		// All mods
		p->gIp2d        *= (1.0 + p->Remodelling_prop*(2 - 1));
	}
    else if (strcmp(p->Remodelling, "AF") == 0)
    {
        // Toy AF remodelling model for minimal model (gives main features of AF)
        p->GNa                  *= (1.0 + p->Remodelling_prop*(0.5  -1));
        p->GCaL                 *= (1.0 + p->Remodelling_prop*(0.35  -1));
        p->Gto                  *= (1.0 + p->Remodelling_prop*(0.2  -1));
        p->GKur                 *= (1.0 + p->Remodelling_prop*(0.75  -1));
        p->GKr                  *= (1.0 + p->Remodelling_prop*(3  -1));
        p->GKs                  *= (1.0 + p->Remodelling_prop*(3  -1));

        p->Gup                  *= (1.0 + p->Remodelling_prop*(2.5  -1));
        p->GRyR_kCO             *= (1.0 + p->Remodelling_prop*(2  -1));
    }
	else
	{	
		printf("ERROR: \"%s\" is not a valid Remodelling for the minimal model. Please check Model_minimal.cpp for options\n\n", p->Remodelling);
		exit(1);
	}
}

void set_modulation_Mutation_native_minimal(Cell_parameters *p)
{   
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G")
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code - if a parameter is set (rather than modified) it will overwrite previous
	// settings; this may or may not be desired
	if (strcmp(p->Mutation, "none") == 0); // do nothing
	else if (strcmp(p->Mutation, "test") == 0) // model specific option
	{
		// All mods
		p->gIp2d        *= 2;
	}
	else 
	{
		printf("ERROR: \"%s\" is not a valid Mutation for the minimal model. Please check Model_minimal.cpp for options\n\n", p->Mutation);
		exit(1);
	}
}

// ACh ===================================\\|
void set_modulation_ACh_minimal(Cell_parameters *p)
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
		printf("ERROR: \"%s\" is not a valid ACh for the minimal model. Please check Model_minimal.cpp for options\n\n", p->ACh_model);
		exit(1);
	}

}
// End ACh ===============================//|
// End Het and modulation functions ===============================//|
// end Parameters and specific settings =========================================================//|

// Compute model functions ======================================================================\\|
void compute_model_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	set_gate_rates_minimal_native(p, var, Vm);
	update_gating_variables_minimal_native(p, var, s, Vm, dt);
	compute_Itot_minimal_native(p, var, s, Vm);
}

void compute_model_minimal_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	set_gate_rates_minimal_native(p, var, Vm);
	update_gating_variables_minimal_native(p, var, s, Vm, dt);
	compute_Itot_minimal_integrated(p, var, s, Vm);
}

void set_gate_rates_minimal_native(Cell_parameters p, Model_variables *var, double Vm)
{
	set_Ip0d_rates(p, var, Vm);
	set_Ip1r_rates(p, var, Vm);
	set_Ip2d_rates(p, var, Vm);
	set_Ip2r_rates(p, var, Vm);
	set_Ip3r_rates(p, var, Vm);
	set_Ip4r_variables(p, var, Vm);

	// set to biophys variable equivilents for outputs
	var->INa_va_al			= var->Ip0d_va_al;
	var->INa_va_bet			= var->Ip0d_va_bet;
	var->INa_va_ss			= var->Ip0d_va_ss;
	var->INa_va_tau			= var->Ip0d_va_tau;
	var->INa_vi_1_al		= var->Ip0d_vi_1_al;
	var->INa_vi_1_bet		= var->Ip0d_vi_1_bet;
	var->INa_vi_1_ss		= var->Ip0d_vi_1_ss;
	var->INa_vi_1_tau		= var->Ip0d_vi_1_tau;
	var->INa_vi_2_al		= var->Ip0d_vi_2_al;
	var->INa_vi_2_bet		= var->Ip0d_vi_2_bet;
	var->INa_vi_2_ss		= var->Ip0d_vi_2_ss;
	var->INa_vi_2_tau		= var->Ip0d_vi_2_tau;
	var->Ito_va_ss			= var->Ip1r_va_ss;
	var->Ito_va_tau			= var->Ip1r_va_tau;
	var->Ito_vi_ss			= var->Ip1r_vi_ss;
	var->Ito_vi_tau			= var->Ip1r_vi_tau;
	var->ICaL_va_ss			= var->Ip2d_va_ss;
	var->ICaL_va_tau		= var->Ip2d_va_tau;
	var->ICaL_vi_ss			= var->Ip2d_vi_ss;
	var->ICaL_vi_tau		= var->Ip2d_vi_tau;
	var->IKur_va_ss			= var->Ip2r_va_ss;
	var->IKur_va_tau		= var->Ip2r_va_tau;
	var->IKur_vi_ss			= var->Ip2r_vi_ss;
	var->IKur_vi_tau		= var->Ip2r_vi_tau;
	var->IKr_va_ss			= var->Ip3r_va_ss;
	var->IKr_va_tau			= var->Ip3r_va_tau;
	var->IKr_vi_ti			= var->Ip3r_vi_ti;
	var->IK1_va_ti			= var->Ip4r_va_ti;
}

void update_gating_variables_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	update_gates_Ip0d(p, var, s, Vm, dt);
	update_gates_Ip1r(p, var, s, Vm, dt);
	update_gates_Ip2d(p, var, s, Vm, dt);
	update_gates_Ip2r(p, var, s, Vm, dt);
	update_gates_Ip3r(p, var, s, Vm, dt);
	// No update for Ip4r
}

void compute_Itot_minimal_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot	= 0;

	compute_Ip0d(p, var, s, Vm);
	compute_Ip1r(p, var, s, Vm);
	compute_Ip2d(p, var, s, Vm);
	compute_Ip2r(p, var, s, Vm);
	compute_Ip3r(p, var, s, Vm);
	compute_Ip4r(p, var, s, Vm);

	var->Itot 	= var->Ip0d + var->Ip1r + var->Ip2d + var->Ip2r + var->Ip3r + var->Ip4r;
}

void compute_Itot_minimal_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Itot   = 0;

	compute_Ip0d(p, var, s, Vm);
	compute_Ip1r(p, var, s, Vm);
	compute_Ip2r(p, var, s, Vm);
	compute_Ip3r(p, var, s, Vm);
	compute_Ip4r(p, var, s, Vm);

	var->Itot   = var->Ip0d + var->Ip1r + var->Ip2r + var->Ip3r + var->Ip4r;
	var->Itot	+= var->INCX + var->ICaP + var->ICab + var->ICaL;	// Adds Ca currents computed in CRU
}
// End Compute model functions ==================================================================//|

// Current formulations =========================================================================\\|
// Ip0d =====================================================================\\|
void set_Ip0d_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac 		= Vm;
	double Vm_inac 		= Vm;

	// Set Activation gate alpha and beta
	var->Ip0d_va_al  				= 0.32*(Vm_ac+47.13)/(1-exp(-0.1*(Vm_ac+47.13)));
	var->Ip0d_va_bet				= 0.08*exp(-Vm_ac/11);

	// Set inactivation gates alphas and betas
	if (Vm_inac < -40.0)
	{
		var->Ip0d_vi_1_al			= 0.135*exp((80+Vm_inac)/-6.8);
		var->Ip0d_vi_1_bet  		= 3.56*exp(0.079*Vm_inac)+310000*exp(0.35*Vm_inac);
		var->Ip0d_vi_2_al			= (-127140*exp(0.2444*Vm_inac)-0.00003474*exp(-0.04391*Vm_inac))*((Vm_inac+37.78)/(1+exp(0.311*(Vm_inac+79.23))));
		var->Ip0d_vi_2_bet 			= (0.1212*exp(-0.01052*Vm_inac))/(1+exp(-0.1378*(Vm_inac+40.14)));
	}
	else
	{
		var->Ip0d_vi_1_al			= 0;
		var->Ip0d_vi_1_bet			= 1/(0.13*(1+exp((Vm_inac+10.66)/-11.1)));
		var->Ip0d_vi_2_al			= 0;
		var->Ip0d_vi_2_bet			= (0.3*exp(-0.0000002535*Vm_inac))/(1+exp(-0.1*(Vm_inac+32)));
	}

	// Set tau and SS from alpha and beta
	var->Ip0d_va_tau 				= 1/(var->Ip0d_va_al + var->Ip0d_va_bet); // 1/(a+b)
	var->Ip0d_vi_1_tau 				= 1/(var->Ip0d_vi_1_al + var->Ip0d_vi_1_bet);
	var->Ip0d_vi_2_tau				= 1/(var->Ip0d_vi_2_al + var->Ip0d_vi_2_bet);
	var->Ip0d_va_ss 				= var->Ip0d_va_al * var->Ip0d_va_tau; // a*tau
	var->Ip0d_vi_1_ss				= var->Ip0d_vi_1_al * var->Ip0d_vi_1_tau;
	var->Ip0d_vi_2_ss				= var->Ip0d_vi_2_al * var->Ip0d_vi_2_tau;
}

void update_gates_Ip0d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ip0d_va     					= rush_larsen(s->Ip0d_va, var->Ip0d_va_ss, var->Ip0d_va_tau, dt); // lib/Membrane.c
	s->Ip0d_vi_1 					= rush_larsen(s->Ip0d_vi_1, var->Ip0d_vi_1_ss, var->Ip0d_vi_1_tau, dt);
	s->Ip0d_vi_2  					= rush_larsen(s->Ip0d_vi_2, var->Ip0d_vi_2_ss, var->Ip0d_vi_2_tau, dt);

	// For outputs
	s->INa_va		= s->Ip0d_va;
	s->INa_vi_1		= s->Ip0d_vi_1;
	s->INa_vi_2		= s->Ip0d_vi_2;
}

void compute_Ip0d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip0d  		= p.gIp0d * pow(s->Ip0d_va, 3) * s->Ip0d_vi_1 * s->Ip0d_vi_2 * (Vm - 76);
	var->Ip0d		*= p.GNa;

	var->INa		= var->Ip0d;
}
// End Ip0d =================================================================//|

// Ip1r =====================================================================\\|
void set_Ip1r_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.Ito_va_ss_shift; 	// Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.Ito_vi_ss_shift;	// Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.Ito_va_tau_shift;	// Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.Ito_vi_tau_shift;	// Voltage modified by shift applied to inactivation time constant

	// Voltage activation
	var->Ip1r_va_ss			= sigmoid(Vm_ac_ss, 1.0, -11.0*p.Ito_va_ss_kscale); // Vm, V1/2, k
	var->Ip1r_va_tau		= 1000*(0.0035 * exp(-(Vm_ac_tau/30.0)*2) + 0.0015);  // ms
	var->Ip1r_va_tau		*= p.Ito_va_tau_scale;

	// Voltage inactivation
	var->Ip1r_vi_ss			= sigmoid(Vm_inac_ss, -40.5, 11.5*p.Ito_vi_ss_kscale); // Vm, V1/2, k
	var->Ip1r_vi_tau		= 1000*(0.025635 * exp (-((Vm_inac_tau - (-52.45))/15.8827)*2) + 0.01414); // ms
	var->Ip1r_vi_tau		*= p.Ito_vi_tau_scale;
}

void update_gates_Ip1r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ip1r_va				= rush_larsen(s->Ip1r_va, var->Ip1r_va_ss, var->Ip1r_va_tau, dt);
	s->Ip1r_vi				= rush_larsen(s->Ip1r_vi, var->Ip1r_vi_ss, var->Ip1r_vi_tau, dt);

	// For outputs
	s->Ito_va				= s->Ip1r_va;
	s->Ito_vi				= s->Ip1r_vi;
}

void compute_Ip1r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip1r				= p.gIp1r * s->Ip1r_va * s->Ip1r_vi * (Vm - (-88));
	var->Ip1r				*= p.Gto;

	var->Ito				= var->Ip1r; // This is for outputs, as all other models use real current names
}
// End Ip1r =================================================================//|

// Ip2d =====================================================================\\|
void set_Ip2d_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	set_Ip2d_va_rates(p, &var->Ip2d_va_ss, &var->Ip2d_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
	var->Ip2d_va_tau		*= p.ICaL_va_tau_scale;

	set_Ip2d_vi_rates(p, &var->Ip2d_vi_ss, &var->Ip2d_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
	var->Ip2d_vi_tau		*= p.ICaL_vi_tau_scale;
}

void set_Ip2d_va_rates(Cell_parameters p, double *va_ss, double *va_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*va_ss					= sigmoid(Vm_ss, 5.0, -6.24*kscale); // V, V1/2, k
	if (Vm_tau != 5.00)
	{
		*va_tau 			= *va_ss *(1-exp(-(Vm_tau-5.0)/6.24))/(0.035*(Vm_tau-5.0));  // Time constant
	}
	else *va_tau			= *va_ss *(1-exp(-(1e-10)/6.24))/(0.035*(1e-10));
}

void set_Ip2d_vi_rates(Cell_parameters p, double *vi_ss, double *vi_tau, double Vm_ss, double Vm_tau, double kscale)
{
	*vi_ss					= sigmoid(Vm_ss, -32.06, 8.6*kscale);  // V, V1/2, k
	*vi_tau  				= 1.0/(0.0197*exp(- (0.0337*(Vm_tau - (-7)))*(0.0337*(Vm_tau - (-7))) ) + 0.02 );
}

void update_gates_Ip2d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ip2d_va          	= rush_larsen(s->Ip2d_va, var->Ip2d_va_ss, var->Ip2d_va_tau, dt);
	s->Ip2d_vi          	= rush_larsen(s->Ip2d_vi, var->Ip2d_vi_ss, var->Ip2d_vi_tau, dt);

	// For outputs
	s->ICaL_va           	= s->Ip2d_va;
	s->ICaL_vi           	= s->Ip2d_vi;
}

void compute_Ip2d(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip2d           	= p.gIp2d * s->Ip2d_va * s->Ip2d_vi * (Vm - (22));
	var->Ip2d           	*= p.GCaL;

	var->ICaL           	= var->Ip2d; // This is for outputs, as all other models use real current names
}
// End Ip2d =================================================================//|

// Ip2r =====================================================================\\|
void set_Ip2r_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac_ss         = Vm - p.IKur_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac_ss       = Vm - p.IKur_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
	double Vm_ac_tau        = Vm - p.IKur_va_tau_shift;   // Voltage modified by shift applied to activation time constant
	double Vm_inac_tau      = Vm - p.IKur_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

	// Voltage activation
	var->Ip2r_va_ss     	= sigmoid(Vm_ac_ss, -6, -8.6*p.IKur_va_ss_kscale); // Vm, V1/2, k
	var->Ip2r_va_tau    	= 1000*(0.009/(1.0 + exp((Vm_ac_tau - (-5))/12)) + 0.0005);  // ms
	var->Ip2r_va_tau    	*= p.IKur_va_tau_scale;

	// Voltage inactivation
	var->Ip2r_vi_ss     	= sigmoid(Vm_inac_ss, -7.5, 10*p.IKur_vi_ss_kscale); // Vm, V1/2, k
	var->Ip2r_vi_tau    	= 1000*(0.59/(1 + exp((Vm_inac_tau - (-60.0))/10)) + 3.05); // ms
	var->Ip2r_vi_tau    	*= p.IKur_vi_tau_scale;
}

void update_gates_Ip2r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ip2r_va          	= rush_larsen(s->Ip2r_va, var->Ip2r_va_ss, var->Ip2r_va_tau, dt);
	s->Ip2r_vi          	= rush_larsen(s->Ip2r_vi, var->Ip2r_vi_ss, var->Ip2r_vi_tau, dt);

	// For outputs
	s->IKur_va          	= s->Ip2r_va;
	s->IKur_vi          	= s->Ip2r_vi;
}

void compute_Ip2r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip2r           	= p.gIp2r * s->Ip2r_va * s->Ip2r_vi * (Vm - (-88));
	var->Ip2r           	*= p.GKur;

	var->IKur            	= var->Ip2r; // This is for outputs, as all other models use real current names
}
// End Ip2r =================================================================//|

// Ip3r =====================================================================\\|
void set_Ip3r_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac         	= Vm - p.IKr_va_ss_shift;    // Voltage modified by shift applied to activation steady state
	double Vm_inac         	= Vm - p.IKr_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state

	// Alpha
	if (fabs(Vm_ac+14.1)< 1e-10) 		var->Ip3r_va_al  = 0.0015; // Denominator = 0 clause
	else 								var->Ip3r_va_al  = 0.0003*(Vm_ac+14.1)/(1-exp((Vm_ac+14.1)/-5));

	// Beta
	if (fabs(Vm_ac-3.3328) < 1e-10) 	var->Ip3r_va_bet = 3.7836118e-4;
	else 								var->Ip3r_va_bet = 0.000073898*(Vm_ac-3.3328)/(exp((Vm_ac-3.3328)/5.1237)-1);

	// time constant and steady state
	var->Ip3r_va_tau		=  1/(var->Ip3r_va_al+var->Ip3r_va_bet);
	var->Ip3r_va_tau    	*= p.IKr_va_tau_scale;

	var->Ip3r_va_ss			= sigmoid(Vm_ac, -14.10, -6.5*p.IKr_va_ss_kscale); // V, V1/2, k

	// Time-independent inactivation gate
	var->Ip3r_vi_ti			= sigmoid(Vm_inac, -15, 22.4*p.IKr_vi_ss_kscale);
	var->IKr_vi_ti			= var->Ip3r_vi_ti; // for outputs

}

void update_gates_Ip3r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->Ip3r_va          	= rush_larsen(s->Ip3r_va, var->Ip3r_va_ss, var->Ip3r_va_tau, dt);

	// For outputs
	s->IKr_va          		= s->Ip3r_va;
}

void compute_Ip3r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip3r           	= p.gIp3r * s->Ip3r_va * var->Ip3r_vi_ti * (Vm - (-88));
	var->Ip3r           	*= p.GKr;

	var->IKr            	= var->Ip3r; // This is for outputs, as all other models use real current names
}
// End Ip3r =================================================================//|

// Ip4r =====================================================================\\|
void set_Ip4r_variables(Cell_parameters p, Model_variables *var, double Vm)
{
	var->Ip4r_va_ti			= (1+exp(0.07*(Vm - (-80 + p.IK1_va_shift))));
	var->IK1_va_ti			= var->Ip4r_va_ti; // for outputs
}
void compute_Ip4r(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->Ip4r				= p.gIp4r * (Vm - (-88 + p.IK1_Erev_shift))/var->Ip4r_va_ti;
	var->Ip4r				*= p.GK1;

	var->IK1				= var->Ip4r;
}
// End Ip4r =================================================================//|
// End Current formulations =====================================================================//|
