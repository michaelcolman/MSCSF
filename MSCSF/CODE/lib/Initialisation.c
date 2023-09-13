// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Initialisation of parameters and variables ==  //
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

#include "Initialisation.h"
#include "Structs.h"
#include "Arguments.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Function list ================================================================================\\|
//	set_simulation_defaults()
//	set_simulation_settings()
//	set_model_conditions()
//	set_model_group_variables()
//	set_local_model_conditions()
//	set_default_parameters()
//	set_parameters_spatial_Ca_defaults()
//	update_parameters_spatial_Ca_0D()
//	initialise_measurement_variables()
//	set_modification_defaults_native()
//	assign_modification_from_arguments()
//	assign_concentrations_from_arguments()
// End Function list ============================================================================//|

// Simulation settings ==========================================================================\\|
// Defaults
void set_simulation_defaults(Simulation_parameters *sim, double dt)
{
	sim->BCL				= 1000; 						// ms
	sim->NBeats				= 10;
	sim->Paced_time			= (sim->NBeats-1) * sim->BCL + 5; // up to end of last applied stimulus (which can be up to 5 ms long)
	sim->Total_time			= sim->Paced_time + 2000 - 5;  	// 2000 ms quiescent period
	sim->dt					= dt;							// ms

	sim->reference			= "";
	sim->results_reference	= "";
	sim->state_reference_read  = "none";
	sim->state_reference_write = "none";

	// S2
	sim->S2_CL				= 0; 		// No S2
	sim->NS2				= 1;		// Default to 1 S2 stimulus if applied
	sim->S2_time			= 0;

	// Spatial output
	sim->Spatial_output_interval_vtk    = 5;	// 5 ms is default
	sim->Spatial_output_interval_data   = 1;	// 1 ms is default
    sim->Spatial_output_start_time      = 0;    
    sim->Spatial_output_end_time        = sim->Total_time;

	sim->Delayed_CaSR_IC    = "Off";
	sim->CaSR_IC_delay      = 1000; // ms
	sim->CaSR_set           = false;
}

// Sets stim variables, model type etc dependant on input arguments
void set_simulation_settings(Simulation_parameters *sim, Argument_parameters A, const char * Model_type)
{
	// NOTE: throughout this function, native models are automatically set with no quiescnet period,
	// whereas integrated models are set to 2000 ms quiescent period. Pass "Total_time" and "Paced_time"
	// as arguments to control these directly

	// Sim time and pacing=========================================\\|
	if (A.reference_arg == true)                sim->reference			= A.reference;
	if (A.results_reference_arg == true)        sim->results_reference	= A.results_reference;
	if (A.state_reference_read_arg == true)     sim->state_reference_read	= A.state_reference_read;
	if (A.state_reference_write_arg == true)    sim->state_reference_write	= A.state_reference_write;

	if (A.BCL_arg == true) 
	{	
		sim->BCL 						= A.BCL;
		sim->Paced_time         		= (sim->NBeats-1) * sim->BCL + 5;
		if (strcmp(Model_type, "native") == 0) 		sim->Total_time		= sim->NBeats * sim->BCL;
		if (strcmp(Model_type, "integrated") == 0) 	sim->Total_time 	= sim->Paced_time + 2000 - 5;
	}
	if (A.NBeats_arg == true)
	{
		sim->NBeats						= A.NBeats;
		sim->Paced_time         		= (sim->NBeats-1) * sim->BCL + 5;
		if (A.Total_time_arg == true)
			sim->Total_time				= A.Total_time; 	
		else
		{	
			if (strcmp(Model_type, "native") == 0) 		sim->Total_time = sim->NBeats * sim->BCL;
			if (strcmp(Model_type, "integrated") == 0) 	sim->Total_time = sim->Paced_time + 2000 - 5;
		}
	}
	else if (A.Total_time_arg == true)
	{
		sim->Total_time             	= A.Total_time;
		if (A.Paced_time_arg == true)
			sim->Paced_time				= A.Paced_time;
		else
		{ 
			if (strcmp(Model_type, "native") == 0)			sim->Paced_time	= sim->Total_time;
			if (strcmp(Model_type, "integrated") == 0)		sim->Paced_time	= sim->Total_time - 2000;	
		}
		sim->NBeats						= (sim->Paced_time + (sim->BCL-2))/sim->BCL;
		sim->Paced_time					= (sim->NBeats-1) * (sim->BCL) + 5;
	}
	if (sim->NBeats < 0) sim->NBeats 	= 0;
	if (A.S2_arg == true)
	{
		sim->S2_CL						= A.S2_CL;
		if (A.NS2_arg == true) sim->NS2	= A.NS2;
		sim->S2_time					= (sim->NS2) * sim->S2_CL + sim->Paced_time;
		if (A.Total_time_arg == false)	sim->Total_time = sim->S2_time + 2000 - 5;
	}

    // Update output time range based on total time
    sim->Spatial_output_end_time        = sim->Total_time;
	// End Sim time and pacing=====================================//|

	// dt
	if (A.dt_arg	== true)	sim->dt = A.dt;

	// Vclamp on/off
	sim->Vclamp		= A.Vclamp;

	// Read/Write state
	sim->Write_state	= A.Write_state;
	sim->Read_state		= A.Read_state;

	// Spatial output interval
	if (A.SOI_arg 	== true) 	sim->Spatial_output_interval_vtk 	= A.SOI;
	if (A.SOId_arg 	== true) 	sim->Spatial_output_interval_data 	= A.SOId;
    if (A.SORs_arg  == true)    sim->Spatial_output_start_time      = A.SORs;
    if (A.SORe_arg  == true)    sim->Spatial_output_end_time        = A.SORe;

	// Delayed CaSR IC functionality
	if (A.Delayed_CaSR_IC_arg == true) 	sim->Delayed_CaSR_IC 	= A.Delayed_CaSR_IC;
	if (A.CaSR_IC_delay_arg == true)	sim->CaSR_IC_delay		= A.CaSR_IC_delay;
}
// End simulation settings ======================================================================//|

// Model conditions (Model type, remodelling etc) ===============================================\\|
void set_model_conditions(Cell_parameters *p, Argument_parameters A)
{
	// Defaults
	p->Model                = "minimal";  // minimal model is default model type
	p->Agent                = "none";
	p->Agent_prop		    = 1.0; 	// full agent effect  (only does anything if Agent is not "none")
	p->Remodelling          = "none";
	p->Remodelling_prop	    = 1.0;	// full remodelling effect (same as above)
	p->ISO_model		    = "default";
	p->ISO                  = 0.0;	// uM / 0 - saturating
	p->ACh_model		    = "default";
	p->ACh				    = 0.0;	// uM / 0 - saturating
	p->spatial_gradient     = "none";
	p->spatial_gradient_prop = 0.0; // one end of gradient (e.g. apex if apico-basal het, where 1 = base)
	p->Mutation			    = "none";

	// (hAM only) 
	p->hAM                = false;
	p->environment        = "intact"; // Needs a default for state files anyway

    // integrated model specific
    p->Ca_cellular_het    = "On";

	// Set by arguments
	if (A.Model_arg == true)   				p->Model          			= A.Model;
	if (A.ISO_arg == true)              	p->ISO            			= A.ISO;
	if (A.ACh_arg == true)              	p->ACh            			= A.ACh;
	if (A.Remodelling_arg == true)      	p->Remodelling    			= A.Remodelling;
	if (A.Agent_arg == true)            	p->Agent          			= A.Agent;
	if (A.Mutation_arg == true)      		p->Mutation   	  			= A.Mutation;
	if (A.Remodelling_prop_arg == true) 	p->Remodelling_prop 		= A.Remodelling_prop;
	if (A.Agent_prop_arg == true) 			p->Agent_prop 				= A.Agent_prop;
	if (A.spatial_gradient_arg == true) 	p->spatial_gradient			= A.spatial_gradient;
	if (A.spatial_gradient_prop_arg == true)p->spatial_gradient_prop	= A.spatial_gradient_prop;
	//if (A.ISO_model_arg == true)		p->ISO_model	  = A.ISO_model;  // Set in Main after parameters are set to allow model-dependant default and argument overwriting

	// Spatial cell model specific
	p->tau_ss_type	= "slow";		// slowest/control -> argument overwriting done later as may also be celltype/remodelling dependent

    // Calcium model heterogeneity
    if(A.Ca_cellular_het_arg == true)       p->Ca_cellular_het = A.Ca_cellular_het;
}

void set_model_group_variables(Cell_parameters *p, Argument_parameters A)
{
	// hAM specific settings
	if (strcmp(p->Model, "hAM_CRN") == 0 || strcmp(p->Model, "hAM_GB") == 0 || strcmp(p->Model, "hAM_NG") == 0 || strcmp(p->Model, "hAM_MT") == 0 \
			|| strcmp(p->Model, "hAM_WL_CRN") == 0 || strcmp(p->Model, "hAM_WL_GB") == 0 || strcmp(p->Model, "hAM_CRN_mWL") == 0 || strcmp(p->Model, "hAM_GB_mWL") == 0 || strcmp(p->Model, "hAM_NG_mWL") == 0)
	{
		p->hAM    = true;
	}

	if (p->hAM    == true)
	{
		if (A.environment_arg 	== true) 	p->environment  = A.environment;
	}
}

void set_local_model_conditions(Cell_parameters p_in, Cell_parameters *p)
{
	// p_in is Global_params; p is local (array) of params
	p->Model			= p_in.Model;
	p->Agent			= p_in.Agent;
	p->Remodelling		= p_in.Remodelling;
	p->Agent_prop		= p_in.Agent_prop;
	p->Remodelling_prop	= p_in.Remodelling_prop;
	p->ISO				= p_in.ISO;
	p->ISO_model		= p_in.ISO_model;
	p->Mutation			= p_in.Mutation;
	p->hAM				= p_in.hAM;
	p->environment		= p_in.environment;
	p->ACh				= p_in.ACh;
	p->ACh_model		= p_in.ACh_model;
}
// End Model conditions (Model type, remodelling etc) ===========================================//|

// Default parameters (parameters which apply to all models; model specific can overwrite =======\\|
// Native Ca handling =======================================================\\|
void set_default_parameters(Cell_parameters *p)
{
	// Constants ==================================================\\|
	p->R			= 8.314;    // J mol^-1 K^-1
	p->F			= 96.485;   // C mmol-1
	p->T			= 310;      // K
	p->FoRT			= p->F/(p->R*p->T);
	// End Constants ==============================================//|

	// Stimulus parameters ========================================\\|
	p->stimduration	= 5.0; 		// ms
	p->stimmag		= -12.5;	// pA/pF
	// End Stimulus parameters ====================================//|

	// Concentrations / environment ===============================\\|
	// Concentrations (constant OR initial condition)
	p->Nai              = 7.95;     // mM
	p->Nao              = 140;      // mM
	p->Ki               = 143.103;  // mM
	p->Ko               = 4.5;      // mM
	p->Cai              = 0.000102; // mM 
	p->Cao              = 1.8;      // mM
	p->Cm				= 100;		// pF
	// End Concentrations =========================================//|

	// Reference for multi files for het and mod ==================\\| 
	p->Het_set_ref			= 0; // 0 = not set (replace with booleans??)
	p->ISO_set_ref			= 0;
	p->Agent_set_ref		= 0;
	p->Remodelling_set_ref 	= 0;
	p->Mutation_set_ref		= 0;
	p->ACh_set_ref			= 0;
	p->Celltype             = "default";	// This is overwritten in sepcific params if default is a specific cell type
	// End reference for multi files for het and mod ==============//| 
}
// End native Ca handling ===================================================//|

// Spatial Ca handling ======================================================\\|
void set_parameters_spatial_Ca_defaults(Cell_parameters *p)
{
	// Cell-structure (volumes) ===================================\\|
	p->Cm_CRU		= 0.005;					// pF
	p->vcyt_CRU     = 1.0;       				// micro m^3
	p->vss_CRU      = 0.0175;					// micro m^3	
	p->vnsr_CRU     = 0.1;						// micro m^3
	p->vjsr_CRU     = 0.015;					// micro m^3
	p->vds_CRU_mean = 1.712e-3;					// micro m^3
	p->v_ss_cyt     = p->vss_CRU/p->vcyt_CRU;   // ratio
	p->v_cyt_nsr    = p->vcyt_CRU/p->vnsr_CRU;
	p->v_jsr_nsr    = p->vjsr_CRU/p->vnsr_CRU;
	// End cell structure =========================================//|

	// Initial conditions for Cai =================================\\|
	p->Cai			= 0.076;        // uM
	p->CaSR			= 750;			// uM
	// End Initial conditions for Cai =============================//|

	// Dyad =======================================================\\|
	// Structure
	p->NRyR_mean	    = 100;			// Number of RyRs per dyad
    p->NRyR_propvar     = 0.4;          // 100 +/- 100*0.4 RyRs per dyad
	p->NLTCC_mean	    = 13;			// Number of LTCCs per dyad	
    p->NLTCC_propvar    = 0.28;          // 13 +/- 13*0.28 LTCCs per dyad

	// RyR
	p->J_rel_max	= 0.000205;			// (um^3/ms)
	p->RyR_kCO_A	= 1.58125e-4;		// (uM^H/ms)
	p->RyR_kOC		= 1.0;				// (ms^-1)
	p->RyR_Cads_H	= 2.5;
	p->RyR_kOC_det	= 1.0;				// multiplier for deterministic model

	p->RyR_mon_tau 		= 25;           // ms
	p->RyR_mi_tau    	= 30;           // ms
	p->RyR_mon_beta_tau	= 200;          // ms
	p->RyR_mi_beta_tau  = 75;           // ms
	p->RyR_mon_grad 	= 2.0;

	// LTCC
	p->J_LTCC_max		= 0.0204;		// micro mol C-1 ms^-1 * (micro m^3) || 11.9 micro mol C-1 ms^-1 * 1.712e-3 micro m^3 
	p->LTCC_kva1_va2	= 0.3;			// ms^-1
	p->LTCC_kva2_va1	= 6.0;			// ms^-1
	p->LTCC_Ca_bar		= 6.0;			// uM
	// End dyad ===================================================//|

	// Inter compatrment transfer =================================\\|
	p->tau_ds       = 0.022;    // ms
	p->tau_ss_cyt	= 0.1;		// ms
	p->tau_nsr_jsr	= 5;		// ms
	// End Inter compatrment transfer =============================//|

	// Spatial coupling ===========================================\\|
	p->tau_ss_trans_f	= 1.0;		// ms  // fast
	p->tau_ss_long_f	= 1.5;		// ms
	p->tau_ss_trans_mf	= 1.15;		// ms  // medium-fast
	p->tau_ss_long_mf	= 1.75;		// ms
	p->tau_ss_trans_m	= 1.25;		// ms  // medium
	p->tau_ss_long_m	= 1.95;		// ms
	p->tau_ss_trans_ms	= 1.3;		// ms  // medim-slow
	p->tau_ss_long_ms	= 2.1;		// ms
	p->tau_ss_trans_s	= 1.35;		// ms  // slow
	p->tau_ss_long_s	= 2.2;		// ms
	p->tau_cyto_trans	= 2.3;		// ms
	p->tau_cyto_long	= 2.9;		// ms
	p->tau_nsr_trans	= 7;		// ms
	p->tau_nsr_long		= 12;		// ms
	//p->tau_ss_trans		= p->tau_ss_trans_s;// ms  assigned in function
	//p->tau_ss_long		= p->tau_ss_long_s;	// ms  from above params
	// End Spatial coupling =======================================//|

	// Buffering ==================================================\\|
	p->Kcam         = 7.0;      // uM
	p->Bcam         = 24.0;     // uM
	p->Kbsr         = 0.6;      // uM
	p->Bbsr         = 47.0;     // uM
	p->Kmca         = 0.033;    // uM
	p->Bmca         = 140.0;    // uM
	p->Kmmg         = 3.64;     // uM
	p->Bmmg         = 140.0;    // uM

	p->Kcsqn        = 0.8;      // mM
	p->Bcsqn        = 10;       // mM
	// End Buffering ==============================================//|

	// SERCA ======================================================\\|
	p->J_SERCA_max		= 0.3729;   	// uM/ms
	p->J_SERCA_kCa		= 0.1605;		// uM
	p->J_SERCA_kCaSR	= 1700;			// uM
	p->J_leak_max		= 1.41265*1e-5; // ms^-1
	p->J_leak_kCaSR		= 450;			// uM
	// End SERCA ==================================================//|

	// Membrane fluxes ============================================\\|
	p->Fjunc            = 0.5;			// propotion of fluxes in subspace

	// NCX
	p->INCX_bar			= 0.3726;		// (um^3.uM.ms^-1)
	p->INCX_kCai   		= 3.59;        	// uM
	p->INCX_kCao		= 1.3;			// mM
	p->INCX_kNai		= 12.3;			// mM
	p->INCX_kNao   		= 87.5;         // mM
	p->INCX_kda			= 0.11;			// uM
	p->INCX_ksat		= 0.27;
	p->INCX_gamma		= 0.35;

	// Cab
	p->ICab_bar			= 1.8237*1e-5;	// (um^3 . uM . ms^-1)

	// CaP
	p->ICaP_bar			= 0.137*1e-2;	// (um^3 . uM . ms^-1)
	p->ICaP_kCa         = 0.5;       	// uM
	// End membrane fluxes ========================================//|
}

void update_parameters_spatial_Ca_0D(Cell_parameters *p)
{
	// This whole thing I really do not like, and will not be necessary in model V2!
	// Move all of this to model specific functions???!!!??!!?

	p->RyR_kCO_A	= 2.371875e-3; // (uM^(1/H)/ms)

	// Multiplier for deterministic model 
	// These are due to not amazing correlation between stochastic and determinitic RyR models
	// this will be improved in a future version of the model
	if (strcmp(p->Model, "minimal") == 0)
	{
		if (strcmp(p->Celltype, "EPI") == 0)						p->RyR_kOC_det  = 0.9; 
		else if (strcmp(p->Celltype, "ENDO") == 0 || strcmp(p->Celltype, "M") == 0)  
		{
			p->GRyR_kCO			*= 1.2;
			if (strcmp(p->Celltype, "M") == 0 && strcmp(p->Remodelling, "none") == 0)
			{
				p->RyR_kOC_det  = 0.65;	
			}
		}
		else if (strcmp(p->Celltype, "RA") == 0)
		{
			if (p->ISO > 0.0 && strcmp(p->Remodelling, "none") == 0) 	p->RyR_kOC_det  = 1.25;
			else if (strcmp(p->Remodelling, "RSERCA_NCX") == 0) 		p->RyR_kOC_det  = 0.9;
		}
	}	
	else if (strcmp(p->Model, "hVM_ORD_s") == 0)
	{
		if (strcmp(p->Remodelling, "RSERCA_NCX") == 0) 	p->RyR_kOC_det  = 0.85;	
	}
	else if (strcmp(p->Model, "hAM_CAZ_s") == 0)
	{
		if (strcmp(p->Remodelling, "none") == 0)
		{
			if (p->ISO > 0.0)	p->RyR_kOC_det  = 2.25;
			else 				p->RyR_kOC_det  = 1.65; // NEW TEST
		}
	}
	else if (strcmp(p->Model, "dAM_VA") == 0)
	{

		p->RyR_kOC_det  = 0.7;
		if (strcmp(p->Remodelling, "short_AP_CICR") == 0)  p->RyR_kCO_A *= 1.75;
		//if (p->ISO > 0.0)   p->RyR_kOC_det  = 3.5;
	}
}
// End Spatial Ca handling ==================================================//|

// Initialise measreument varaibles =========================================\\|
void initialise_measurement_variables(Model_variables *var)
{
	for (int i = 0; i < 9; i++) 
		var->APD_p_switch[i]= -1;       
	var->ex_switch          = 0;             
	var->APD_t_switch       = -1;           
	var->t_ex               = -100;            // Time at which cell was excited   (ms)
	var->dvdt               = 0;               // Rate of change of voltage        (mV/ms)
	var->dvdt_max           = 0;               // Maximum rate of change of voltage(mV/ms)
	var->dvdt_max_prev      = 0;               // Maximum rate of change of voltage(mV/ms)
	var->Vmax               = -80;             // Maximum voltage                  (mV)
	var->Vmax_prev          = -80;             // Maximum voltage, previous beat   (mV)
	var->Vmin               = 50;              // Minimum voltage                  (mV)
	var->Vmin_prev          = 50;              // Minimum voltage, previous beat   (mV) *well, actually at point of excitation for current beat
	var->Vmin_prev_prev     = 50;              // Minimum voltage, previous beat   (mV) *point of excitation for previous beat
	var->Vamp               = 0;               // Amplitude voltage                (mV)
	var->Vamp_prev          = 0;               // Amplitude voltage, previous beat (mV)
	var->APD_t              = 0;               // APD at threshold voltage         (ms)
	for (int i = 0; i < 9; i++)
		var->APD_p[i]       = 0;               // APD at multiple % repolarisation (ms)
	var->APD_t_prev         = 0;               // previous beat APD at threshold   (ms)
	for (int i = 0; i < 9; i++)
		var->APD_p_prev[i]  = 0;               // previous beat APD at %           (ms)
	var->CaT_min            = 1;               // Minimum intracellular Ca2+       (mM)
	var->CaT_max            = 0;               // Maximum intracellular Ca2+       (mM)
	var->CaSR_min           = 5000;            // Minimum SR Ca2+                  (mM)
	var->CaSR_max           = 0;               // Maximum SR Ca2+                  (mM)
	var->CaT_min_prev       = 1;   		       // Minimum intracellular Ca2+ prev  (mM)
	var->CaT_max_prev       = 0;           	   // Maximum intracellular Ca2+ prev  (mM)
	var->CaSR_min_prev      = 5000;            // Minimum SR Ca2+ prev             (mM)
	var->CaSR_max_prev      = 0;               // Maximum SR Ca2+ prev             (mM)
}
// End Initialise measreument varaibles =====================================//|
// End Default parameters =======================================================================//|

// Current modification variables ===============================================================\\|
void set_modification_defaults_native(Cell_parameters *p)
{
	// Scale factors all defaulted to 1
	p->GNa						= 1.0;
	p->GNaL						= 1.0;
	p->Gto						= 1.0;
	p->GCaL						= 1.0;
	p->GKur						= 1.0;
	p->GKr						= 1.0;
	p->GKs						= 1.0;
	p->GK1						= 1.0;
	p->GNCX						= 1.0;
	p->GCaP						= 1.0;
	p->GNab						= 1.0;
	p->GCab						= 1.0;
	p->GKb						= 1.0;
	p->GNaK						= 1.0;
	p->GClCa					= 1.0;
	p->GClb						= 1.0;
	p->GRyR_kCO					= 1.0;
	p->GLTCC_kva1_va2			= 1.0;
	p->GKACh					= 1.0;
	p->Gf						= 1.0;

	p->INa_va_tau_scale			= 1.0;
	p->INa_vi_1_tau_scale		= 1.0;
	p->INa_vi_2_tau_scale		= 1.0;
	p->INaL_va_tau_scale		= 1.0;
	p->INaL_vi_tau_scale		= 1.0;
	p->Ito_va_tau_scale			= 1.0;
	p->Ito_vi_tau_scale			= 1.0;
	p->ICaL_va_tau_scale		= 1.0;
	p->ICaL_vi_tau_scale		= 1.0;
	p->IKur_va_tau_scale		= 1.0;
	p->IKur_vi_tau_scale		= 1.0;
	p->IKr_va_tau_scale			= 1.0;
	p->IKr_vi_tau_scale			= 1.0;
	p->IKs_va_tau_scale			= 1.0;
	p->IKACh_va_tau_scale		= 1.0;
	p->IKACh_vi_tau_scale		= 1.0;
	p->If_va_tau_scale			= 1.0;

	p->INaL_va_ss_kscale		= 1.0;
	p->INaL_vi_ss_kscale		= 1.0;
	p->Ito_va_ss_kscale			= 1.0;
	p->Ito_vi_ss_kscale			= 1.0;
	p->ICaL_va_ss_kscale		= 1.0;
	p->ICaL_vi_ss_kscale		= 1.0;
	p->IKur_va_ss_kscale		= 1.0;
	p->IKur_vi_ss_kscale		= 1.0;
	p->IKr_va_ss_kscale			= 1.0;
	p->IKr_vi_ss_kscale			= 1.0;
	p->IKs_va_ss_kscale			= 1.0;
	p->IKACh_va_ss_kscale		= 1.0;
	p->IKACh_vi_ss_kscale		= 1.0;
	p->If_va_ss_kscale			= 1.0;

	// Voltage shifts all defaulted to 0
	p->INa_va_shift				= 0.0;
	p->INa_vi_shift				= 0.0;
	p->INaL_va_shift			= 0.0;
	p->INaL_vi_shift			= 0.0;

	p->INaL_va_ss_shift			= 0.0;
	p->INaL_vi_ss_shift			= 0.0;
	p->INaL_va_tau_shift		= 0.0;
	p->INaL_vi_tau_shift		= 0.0;
	p->Ito_va_ss_shift			= 0.0;
	p->Ito_vi_ss_shift			= 0.0;
	p->Ito_va_tau_shift			= 0.0;
	p->Ito_vi_tau_shift			= 0.0;
	p->ICaL_va_ss_shift			= 0.0;
	p->ICaL_vi_ss_shift			= 0.0;
	p->ICaL_va_tau_shift		= 0.0;
	p->ICaL_vi_tau_shift		= 0.0;
	p->IKur_va_ss_shift			= 0.0;
	p->IKur_vi_ss_shift			= 0.0;
	p->IKur_va_tau_shift		= 0.0;
	p->IKur_vi_tau_shift		= 0.0;
	p->IKr_va_ss_shift        	= 0.0;
	p->IKr_va_tau_shift       	= 0.0;
	p->IKr_vi_ss_shift        	= 0.0;
	p->IKr_vi_tau_shift       	= 0.0;
	p->IKs_va_ss_shift        	= 0.0;
	p->IKs_va_tau_shift       	= 0.0;
	p->IK1_va_shift				= 0.0;
	p->IKACh_va_ss_shift		= 0.0;
	p->IKACh_va_tau_shift		= 0.0;
	p->IKACh_vi_ss_shift		= 0.0;
	p->IKACh_vi_tau_shift		= 0.0;
	p->If_va_ss_shift			= 0.0;
	p->If_va_tau_shift			= 0.0;

	// Ca handling
	p->Gup						= 1.0;
	p->Gleak					= 1.0;
	p->Grel						= 1.0;

	// Current equation modifiers
	p->IKr_rect_vshift			= 0.0;
	p->IKr_rect_grad_scale		= 1.0;
	p->IKr_rect_exp_factor		= 1.0;
	p->IK1_VEexp_scale			= 1.0;
	p->IK1_den_add_factor		= 1.0;
	p->IK1_Erev_shift			= 0.0;
}

void assign_concentrations_from_arguments(Cell_parameters *p, Argument_parameters A)
{
    // If the argument has been passed, set the parameters to the defined argument, else leave as defaults
    if (A.Nai_arg       == true)    p->Nai      = A.Nai;
    if (A.Nao_arg       == true)    p->Nao      = A.Nao;
    if (A.Ki_arg        == true)    p->Ki       = A.Ki;
    if (A.Ko_arg        == true)    p->Ko       = A.Ko;
    if (A.Cao_arg       == true)    p->Cao      = A.Cao;
}

void assign_modification_from_arguments(Cell_parameters *p, Argument_parameters A)
{
	// If the argument has been passed, set the parameters to the defined argument, else leave as defaults
	if (A.GNa_arg		== true)	p->GNa		*= A.GNa;
	if (A.GNaL_arg		== true)	p->GNaL		*= A.GNaL;
	if (A.Gto_arg		== true)	p->Gto		*= A.Gto;
	if (A.GCaL_arg		== true)	p->GCaL		*= A.GCaL;
	if (A.GKur_arg		== true)	p->GKur		*= A.GKur;
	if (A.GKr_arg		== true)	p->GKr		*= A.GKr;
	if (A.GKs_arg		== true)	p->GKs		*= A.GKs;
	if (A.GK1_arg		== true)	p->GK1		*= A.GK1;
	if (A.GNCX_arg		== true)	p->GNCX		*= A.GNCX;
	if (A.GCaP_arg		== true)	p->GCaP		*= A.GCaP;
	if (A.GNab_arg		== true)	p->GNab		*= A.GNab;
	if (A.GCab_arg		== true)	p->GCab		*= A.GCab;
	if (A.GKb_arg		== true)	p->GKb		*= A.GKb;
	if (A.GNaK_arg		== true)	p->GNaK		*= A.GNaK;
	if (A.GClCa_arg		== true)	p->GClCa	*= A.GClCa;
	if (A.GClb_arg		== true)	p->GClb		*= A.GClb;
	if (A.GKACh_arg		== true)	p->GKACh	*= A.GKACh;
	if (A.Gup_arg		== true)	p->Gup		*= A.Gup;
	if (A.Gleak_arg		== true)	p->Gleak	*= A.Gleak;
	if (A.Grel_arg		== true)	p->Grel		*= A.Grel;

	if (A.RyR_Po_arg    == true)    p->GRyR_kCO 		*= A.RyR_Po;
	if (A.LTCC_Po_arg   == true)    p->GLTCC_kva1_va2 	*= A.LTCC_Po;

	if (A.INa_va_tau_scale_arg		== true)	p->INa_va_tau_scale		*= A.INa_va_tau_scale;
	if (A.INa_vi_1_tau_scale_arg	== true)	p->INa_vi_1_tau_scale	*= A.INa_vi_1_tau_scale;
	if (A.INa_vi_2_tau_scale_arg	== true)	p->INa_vi_2_tau_scale	*= A.INa_vi_2_tau_scale;
	if (A.INaL_va_tau_scale_arg		== true)	p->INaL_va_tau_scale	*= A.INaL_va_tau_scale;
	if (A.INaL_vi_tau_scale_arg		== true)	p->INaL_vi_tau_scale	*= A.INaL_vi_tau_scale;
	if (A.Ito_va_tau_scale_arg		== true)	p->Ito_va_tau_scale		*= A.Ito_va_tau_scale;
	if (A.Ito_vi_tau_scale_arg		== true)	p->Ito_vi_tau_scale		*= A.Ito_vi_tau_scale;
	if (A.ICaL_va_tau_scale_arg		== true)	p->ICaL_va_tau_scale	*= A.ICaL_va_tau_scale;
	if (A.ICaL_vi_tau_scale_arg		== true)	p->ICaL_vi_tau_scale	*= A.ICaL_vi_tau_scale;
	if (A.IKur_va_tau_scale_arg		== true)	p->IKur_va_tau_scale	*= A.IKur_va_tau_scale;
	if (A.IKur_vi_tau_scale_arg		== true)	p->IKur_vi_tau_scale	*= A.IKur_vi_tau_scale;
	if (A.IKr_va_tau_scale_arg		== true)	p->IKr_va_tau_scale		*= A.IKr_va_tau_scale;
	if (A.IKs_va_tau_scale_arg		== true)	p->IKs_va_tau_scale		*= A.IKs_va_tau_scale;
	if (A.IKACh_va_tau_scale_arg	== true)	p->IKACh_va_tau_scale	*= A.IKACh_va_tau_scale;

	if (A.INa_va_shift_arg			== true)	p->INa_va_shift			+= A.INa_va_shift;
	if (A.INa_vi_shift_arg			== true)	p->INa_vi_shift			+= A.INa_vi_shift;
	if (A.INaL_va_shift_arg			== true)	p->INaL_va_shift		+= A.INaL_va_shift;
	if (A.INaL_vi_shift_arg			== true)	p->INaL_vi_shift		+= A.INaL_vi_shift;
	if (A.Ito_va_ss_shift_arg		== true)	p->Ito_va_ss_shift		+= A.Ito_va_ss_shift;
	if (A.Ito_va_tau_shift_arg		== true)	p->Ito_va_tau_shift		+= A.Ito_va_tau_shift;
	if (A.Ito_vi_ss_shift_arg		== true)	p->Ito_vi_ss_shift		+= A.Ito_vi_ss_shift;
	if (A.Ito_vi_tau_shift_arg		== true)	p->Ito_vi_tau_shift		+= A.Ito_vi_tau_shift;
	if (A.ICaL_va_ss_shift_arg		== true)	p->ICaL_va_ss_shift		+= A.ICaL_va_ss_shift;
	if (A.ICaL_va_tau_shift_arg		== true)	p->ICaL_va_tau_shift	+= A.ICaL_va_tau_shift;
	if (A.ICaL_vi_ss_shift_arg		== true)	p->ICaL_vi_ss_shift		+= A.ICaL_vi_ss_shift;
	if (A.ICaL_vi_tau_shift_arg		== true)	p->ICaL_vi_tau_shift	+= A.ICaL_vi_tau_shift;
	if (A.IKur_va_ss_shift_arg		== true)	p->IKur_va_ss_shift		+= A.IKur_va_ss_shift;
	if (A.IKur_va_tau_shift_arg		== true)	p->IKur_va_tau_shift	+= A.IKur_va_tau_shift;
	if (A.IKur_vi_ss_shift_arg		== true)	p->IKur_vi_ss_shift		+= A.IKur_vi_ss_shift;
	if (A.IKur_vi_tau_shift_arg		== true)	p->IKur_vi_tau_shift	+= A.IKur_vi_tau_shift;
	if (A.IKr_va_ss_shift_arg		== true)	p->IKr_va_ss_shift		+= A.IKr_va_ss_shift;
	if (A.IKr_va_tau_shift_arg		== true)	p->IKr_va_tau_shift		+= A.IKr_va_tau_shift;
	if (A.IKr_vi_ss_shift_arg		== true)	p->IKr_vi_ss_shift		+= A.IKr_vi_ss_shift;
	if (A.IKs_va_ss_shift_arg		== true)	p->IKs_va_ss_shift		+= A.IKs_va_ss_shift;
	if (A.IKs_va_tau_shift_arg		== true)	p->IKs_va_tau_shift		+= A.IKs_va_tau_shift;
	if (A.IKACh_va_ss_shift_arg		== true)	p->IKACh_va_ss_shift	+= A.IKACh_va_ss_shift;
	if (A.IKACh_va_tau_shift_arg	== true)	p->IKACh_va_tau_shift	+= A.IKACh_va_tau_shift;
	if (A.IK1_va_shift_arg			== true)	p->IK1_va_shift			+= A.IK1_va_shift;
	if (A.IK1_Erev_shift_arg		== true)	p->IK1_Erev_shift		+= A.IK1_Erev_shift;
	if (A.Ito_shift_arg				== true)
	{
		p->Ito_va_ss_shift			+= A.Ito_shift;
		p->Ito_vi_ss_shift			+= A.Ito_shift;
		p->Ito_va_tau_shift			+= A.Ito_shift;
		p->Ito_vi_tau_shift			+= A.Ito_shift;
	}
	if (A.ICaL_shift_arg			== true)
	{
		p->ICaL_va_ss_shift			+= A.ICaL_shift;
		p->ICaL_vi_ss_shift			+= A.ICaL_shift;
		p->ICaL_va_tau_shift		+= A.ICaL_shift;
		p->ICaL_vi_tau_shift		+= A.ICaL_shift;
	}
	if (A.IKur_shift_arg			== true)
	{
		p->IKur_va_ss_shift			+= A.IKur_shift;
		p->IKur_vi_ss_shift			+= A.IKur_shift;
		p->IKur_va_tau_shift		+= A.IKur_shift;
		p->IKur_vi_tau_shift		+= A.IKur_shift;
	}

	if (A.Ito_va_ss_kscale_arg		== true)	p->Ito_va_ss_kscale		*= A.Ito_va_ss_kscale;
	if (A.Ito_vi_ss_kscale_arg		== true)	p->Ito_vi_ss_kscale		*= A.Ito_vi_ss_kscale;
	if (A.ICaL_va_ss_kscale_arg		== true)	p->ICaL_va_ss_kscale	*= A.ICaL_va_ss_kscale;
	if (A.ICaL_vi_ss_kscale_arg		== true)	p->ICaL_vi_ss_kscale	*= A.ICaL_vi_ss_kscale;
	if (A.IKur_va_ss_kscale_arg		== true)	p->IKur_va_ss_kscale	*= A.IKur_va_ss_kscale;
	if (A.IKur_vi_ss_kscale_arg		== true)	p->IKur_vi_ss_kscale	*= A.IKur_vi_ss_kscale;
	if (A.IKr_va_ss_kscale_arg		== true)	p->IKr_va_ss_kscale		*= A.IKr_va_ss_kscale;
	if (A.IKr_vi_ss_kscale_arg		== true)	p->IKr_vi_ss_kscale		*= A.IKr_vi_ss_kscale;
	if (A.IKs_va_ss_kscale_arg		== true)	p->IKs_va_ss_kscale		*= A.IKs_va_ss_kscale;
	if (A.IKACh_va_ss_kscale_arg	== true)	p->IKACh_va_ss_kscale	*= A.IKACh_va_ss_kscale;
}
// End current modification variables ===========================================================//|

