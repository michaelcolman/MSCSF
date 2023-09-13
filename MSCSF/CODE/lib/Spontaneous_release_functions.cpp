// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Everything for spontaneous release functions=  //
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
// Note: notation: "3D_cell" refers to dynamic models =====  //
// which are fit to the 3D_cell model; "General" is all ===  //
// other controllable dynamic models ======================  //

#include "Spontaneous_release_functions.h"
#include "Structs.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>

// Function list ================================================================================\\|
//	Defaults and setup
//	    set_SRF_defaults()
//	    SRF_setup()
//	    SRF_tissue_heterogeneity()
//	
//	Parameter sets and settings   ** These are where to add new SRF models **
//	    set_SRF_Direct_Control_parameters()
//	    set_SRF_dynamic_parameters_3D_cell()
//	    set_SRF_dynamic_parameters_controllable()
//	
//	Dynamic model functions which determine the parameters for distributions
//	    3D cell fit functions
//	        determine_SRF_3D_cell_PSCR()
//	        determine_SRF_3D_cell_ti_sep()
//	        determine_SRF_3D_cell_ti_dist_widths()
//	        determine_SRF_3D_cell_MD()
//	
//	    Dynamic, controllable
//	        determine_SRF_dynamic_general_PSCR()
//	        determine_SRF_dynamic_general_ti_sep()
//	        determine_SRF_dynamic_general_ti_dist_widths()
//	        determine_dynamic_general_cell_MD()
//	        determine_SRF_dynamic_general_duration_dist_widths()
//	
//	    determine_SRF_params_from_CaSR()
//	    need_to_recalculate()
//	
//	Run SRF functions
//	    set_and_run_SRF()
//	    run_SRF()
//	    calc_SRF_mults()
//
//	    read_SRF_settings_from_file()
//	
//	SRF waveform functions
//	    Determine_waveform_parameters()
//	    Determine_waveform_parameters_long()
//	    Waveform()
//	    Waveform_plateau()
//	
//	Distributions and inverses
//	    Determine_ti()
//	    Determine_duration()
//	    Determine_duration_ks_from_MD()
//	    Determine_NRyRo_peak()
//	    Determine_NRyRo_plateau()
//	
//	Test functions
//	    test_and_produce_distributions()
//	    test_and_produce_CaSR_dependency()
// End function list ============================================================================//|

// Defaults and setup========================================================\\|
void set_SRF_defaults(Spontaneous_release_functions *srf, Argument_parameters A)
{
	srf->Mode		= "Off";		//"Dynamic" or "Direct_Control" or "Read"
	srf->Model		= "3D_cell"; 	// Dynamic only, "3D_cell" or "General"
	srf->Pset		= "Control";
	srf->SRF_het	= "Off";		// "Off", "General_1_small", "General_1_large"
    srf->write_SRF_settings = "Off";

	// Default switch settings
	srf->srf_set				= -3;	// initial value, doesn't mean anything other than ensures cell has undergone an AP before SRF may be set
	srf->srf_calc				= -1;
	srf->waveform_init			= 0;
	srf->CaSR_t_calc			= 0; // mM
    srf->init_write_flag        = false; // does not need to be written to file

	// Dist width scaling
	srf->t_to_peak_dist_scale		= 1; // 0-1 only is valid
	srf->NRyRo_peak_dist_scale		= 1; // 0-2 is valid (2 being twice the observed variance)
	srf->NRyRo_plateau_dist_scale	= 1; // same as above

    if (A.t_to_peak_dist_scale_arg  == true)        srf->t_to_peak_dist_scale       = A.t_to_peak_dist_scale;
    if (A.NRyRo_peak_dist_scale_arg == true)        srf->NRyRo_peak_dist_scale      = A.NRyRo_peak_dist_scale;
    if (A.NRyRo_plateau_dist_scale_arg  == true)    srf->NRyRo_plateau_dist_scale   = A.NRyRo_plateau_dist_scale;

    // Dynamic recalculation
    srf->recalc_SR_diff             = 10; // microM

	// and to overwrite with args
	if (A.SRF_mode_arg == true)
	{
		srf->Mode 	= A.SRF_mode;
		// Set defaults for Model and Pset based on mode
		if (strcmp(srf->Mode, "Dynamic") == 0)
		{
			srf->Model = "3D_cell";
			srf->Pset = "Control";
		}
		else srf->Pset = "General_1";
	}
	if (A.SRF_model_arg == true)	
	{
		srf->Model 	= A.SRF_model;
		// Set default Pset based on which type (if dynamic; will be overwritten on next line if specific Pset argument passed)
		if (strcmp(srf->Mode, "Dynamic") == 0)
		{
			if (strcmp(srf->Model, "3D_cell") == 0) 		srf->Pset = "Control";
			else if (strcmp(srf->Model, "General") == 0) 	srf->Pset = "General_1p10";
		}
		else srf->Pset = "General_1"; // for "Direct_Control" model
	}
	if (A.SRF_Pset_arg == true)	srf->Pset 	= A.SRF_Pset;

	if (strcmp(srf->Mode, "Dynamic") == 0 && strcmp(srf->Model, "General") == 0 && A.SRF_het_arg == true) srf->SRF_het = A.SRF_het;

    if (A.write_SRF_arg == true) srf->write_SRF_settings = A.write_SRF_settings;
}

// Call apprpriate set parameters function
void SRF_setup(Spontaneous_release_functions *srf, Argument_parameters A)
{
	if (strcmp(srf->Mode, "Off") != 0)
	{
		if (strcmp(srf->Mode, "Direct_Control") == 0) set_SRF_Direct_Control_parameters(srf, A);
		else if (strcmp(srf->Mode, "Dynamic") == 0)
		{
			if (strcmp(srf->Model, "3D_cell") == 0) 		set_SRF_dynamic_parameters_3D_cell(srf);
			else if (strcmp(srf->Model, "General") == 0) 	set_SRF_dynamic_parameters_controllable(srf, A);
			else
			{
				printf("ERROR: %s is not a valid Spontaneous Release Function model for \"Dynamic\". Please select from \"3D_cell\" or \"General\"\n", srf->Model);
				exit(1);
			}
		}
		else
		{
			printf("ERROR: %s is not a valid Spontaneous Release Function mode. Please select from \"Direct_Control\" or \"Dynamic\"\n", srf->Mode);
			exit(1); 
		}
	}
}

// Randomly assigned tissue heterogeneity
void SRF_tissue_heterogeneity(Spontaneous_release_functions *srf, double rand)
{
	if (strcmp(srf->SRF_het, "Off") == 0); // do nothing
	else if (strcmp(srf->SRF_het, "General_1_small") == 0)
	{
		if (rand < 0.6)        srf->Pset = "General_1p10";
		else if (rand < 0.7)   srf->Pset = "General_1p05";
		else if (rand < 0.8)   srf->Pset = "General_1p00";
		else if (rand < 0.9)   srf->Pset = "General_1p15";
		else                   srf->Pset = "General_1p20";
	}
	else if (strcmp(srf->SRF_het, "General_1_large") == 0)
	{
		if (rand < 0.2)        srf->Pset = "General_1p10";
		else if (rand < 0.4)   srf->Pset = "General_1p05";
		else if (rand < 0.6)   srf->Pset = "General_1p00";
		else if (rand < 0.8)   srf->Pset = "General_1p15";
		else                   srf->Pset = "General_1p20";
	}
	else 
	{
		printf("ERROR: %s is not a valid SRF_het version\n", srf->SRF_het);
		exit(1);
	}
}
// End Defaults =============================================================//|

// Parameter sets ===========================================================\\|
// Static/Direct_Control model ===========================================\\|
void set_SRF_Direct_Control_parameters(Spontaneous_release_functions *srf, Argument_parameters A)
{
	// Default duration width from median function parameters
	srf->DW_k1_A		= 276.72; // ms
	srf->DW_k1_a		= 324.11; // ms
	srf->DW_k1_k		= 70.67;
	srf->DW_k1_min		= 40;	  // ms
	srf->DW_k2_A		= 233.26; // ms
	srf->DW_k2_a		= 160.74; // ms
	srf->DW_k2_k		= 14.6;
	srf->DW_k2_min		= 53;	  // ms
	srf->duration_width = -10;	  // dummy value meaning it will be set from MD

	// RyRo vs duration params
	srf->NRyRo_peak_median_A	= 692.99;
	srf->NRyRo_peak_median_min	= 0.059;
	srf->NRyRo_peak_median_H	= -1.6;

	if (strcmp(srf->Pset, "none") == 0) srf->PSCRE	= 0;
	else if (strcmp(srf->Pset, "User_control") == 0) // Assign ALL parameters from user specified arguments
	{
		if (A.SRF_N_DC_set >= 7)
		{
			srf->PSCRE          = A.SRF_PSCRE;
			srf->CF_ti_sep      = A.SRF_CF_ti_sep;
			srf->ti_sep         = A.SRF_ti_sep;      		// relative to time of last excitation
			srf->k_ti_F1_ms     = A.SRF_ti_dist_width_1;    // ms, ti-this = earliest SCRE 
			srf->k_ti_F2_ms     = A.SRF_ti_dist_width_2;    // ms, ti+this = latest SCRE
			srf->MD             = A.SRF_MD;
			srf->duration_width	= A.SRF_duration_width;
			//srf->duration_width = 150;      

			// convert ms to gradient (approx)
			srf->k_ti_F1    = 0.145*srf->k_ti_F1_ms;
			srf->k_ti_F2    = 0.145*srf->k_ti_F2_ms;
		}
		else 
		{
			printf("ERROR: User_control Pset has been set for Direct_control SRF model, but the wrong number of parameters have bene set;  %d vs the required 7\n", A.SRF_N_DC_set);
			printf("\tFor this model, you must pass in all of: SRF_DC_{PSCRE/CF/ti_sep/ti_W1/ti_W2/MD/duration_W}\n");
			exit(1);
		}
	}
	else if (strcmp(srf->Pset, "General_1") == 0) // Example implementation
	{
		srf->PSCRE			= 1.0;
		srf->CF_ti_sep		= 0.1;
		srf->ti_sep			= 500;		// relative to time of last excitation
		srf->k_ti_F1_ms		= 200;		// ms, ti-this = earliest SCRE 
		srf->k_ti_F2_ms		= 800;		// ms, ti+this = latest SCRE
		srf->MD				= 250;
		//srf->duration_width = 150;	  

		// convert ms to gradient (approx)
		srf->k_ti_F1	= 0.145*srf->k_ti_F1_ms;
		srf->k_ti_F2	= 0.145*srf->k_ti_F2_ms;
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid SRF Pset for the Direct_Control model.\n", srf->Pset);
		exit(1);
	}	
}
// End Static/Direct_Control model =======================================//|

// Dynamic model - fit to 3D cell==================================\\|
void set_SRF_dynamic_parameters_3D_cell(Spontaneous_release_functions *srf)
{
	srf->duration_width = -10;    // dummy value meaning it will be set from MD

	// Defaults for duration distribution from median (also defined for Direct_Control model)
	srf->DW_k1_A        = 276.72; // ms
	srf->DW_k1_a        = 324.11; // ms
	srf->DW_k1_k        = 70.67;
	srf->DW_k1_min      = 40;     // ms
	srf->DW_k2_A        = 233.26; // ms
	srf->DW_k2_a        = 160.74; // ms
	srf->DW_k2_k        = 14.6;
	srf->DW_k2_min      = 53;     // ms

	// RyRo vs duration params
	srf->NRyRo_peak_median_A	= 692.99;
	srf->NRyRo_peak_median_min	= 0.059;
	srf->NRyRo_peak_median_H	= -1.6;

	// Defaults / control 3D cell (	tau_ss = slow; no NCX/RyR remodelling )
	// Probability
	srf->PSCR_threshold		= 1.091;		// mM
	srf->PSCR_k				= 0.00577;

	// ti distribution
	srf->ti_sep_a			= 1.157;		// mM
	srf->ti_sep_A			= 37.78;
	srf->ti_sep_k			= 0.041;
	srf->ti_sep_min			= 31.83;		// ms
	srf->CF_ti_sep_a		= 1.192;		// mM
	srf->CF_ti_sep_k		= 0.0311;
	srf->CF_ti_sep_A		= 0.25;

	srf->KF1_a				= 1.136;		// mM
	srf->KF1_A				= 8.619;
	srf->KF1_k				= 0.0449;
	srf->KF1_min			= 2.967;

	srf->KF2_A				= 583.603;
	srf->KF2_B				= 295149.038;
	srf->KF2_H1				= -18.979;
	srf->KF2_H2				= -68.306;
	srf->KF2_min			= 5;

	// duration distribution
	srf->MD_a				= 1.136;		// mM
	srf->MD_A				= 208.57;
	srf->MD_k				= 0.05255;
	srf->MD_min				= 90;			// ms 

	if (strcmp(srf->Pset, "Control") == 0); 		// Do nothing for control, as that is to what above parameters correspond
	else if (strcmp(srf->Pset, "RSERCA_NCX") == 0)	// Remodelling SERCA/NCX 
	{
		// Probability
		srf->PSCR_threshold     = 0.993;        // mM
		srf->PSCR_k             = 0.007;

		// ti distribution
		srf->ti_sep_a           = 1.096;        // mM
		srf->ti_sep_A           = 27.6264;
		srf->ti_sep_k           = 0.0581;
		srf->ti_sep_min         = 18;        	// ms
		srf->CF_ti_sep_a        = 1.124;        // mM
		srf->CF_ti_sep_k        = 0.0481;
		srf->CF_ti_sep_A        = 0.25;

		srf->KF1_a              = 1.100;        // mM
		srf->KF1_A              = 3.955;
		srf->KF1_k              = 0.0704;
		srf->KF1_min            = 1.532;

		srf->KF2_A              = 173.361;
		srf->KF2_B              = 338.7416;
		srf->KF2_H1             = -21.099;
		srf->KF2_H2             = -118.7711;
		srf->KF2_min            = 3.9;

		// duration distribution
		srf->MD_a               = 0.953;        // mM
		srf->MD_A               = 801.86;
		srf->MD_k               = 0.063;
		srf->MD_min             = 90;           // ms 	

		// Duration distribution from median
		srf->DW_k1_A        	= 273.53; 		// ms
		srf->DW_k1_a        	= 268.012; 		// ms
		srf->DW_k1_k        	= 40.24;
		srf->DW_k1_min      	= 48.26;     	// ms
		srf->DW_k2_A        	= 369.44; 		// ms
		srf->DW_k2_a        	= 188.14; 		// ms
		srf->DW_k2_k        	= 18.716;
		srf->DW_k2_min      	= 60.56;
	}
	else if (strcmp(srf->Pset, "RCRU") == 0)  // Remodelling CRU coupling
	{
		// Probability
		srf->PSCR_threshold     = 0.963;        // mM
		srf->PSCR_k             = 0.00749;

		// ti distribution
		srf->ti_sep_a           = 1.017;        // mM
		srf->ti_sep_A           = 65.047;
		srf->ti_sep_k           = 0.0674;
		srf->ti_sep_min         = 27;           // ms
		srf->CF_ti_sep_a        = 0.99;        // mM
		srf->CF_ti_sep_k        = 0.0111;
		srf->CF_ti_sep_A        = 0.15;

		srf->KF1_a              = 1.04;        // mM
		srf->KF1_A              = 10.93;
		srf->KF1_k              = 0.0659;
		srf->KF1_min            = 1.529;

		srf->KF2_A              = 181.12;
		srf->KF2_B              = 0;
		srf->KF2_H1             = -31.02;
		srf->KF2_H2             = 1;
		srf->KF2_min            = 6.7;

		// duration distribution
		srf->MD_a               = 1.0312;        // mM
		srf->MD_A               = 136.86;
		srf->MD_k               = 0.08544;
		srf->MD_min             = 64.12;         // ms   

		// Duration distribution from median
		srf->DW_k1_A            = 155.61;       // ms
		srf->DW_k1_a            = 223.98;      	// ms
		srf->DW_k1_k            = 57.27;
		srf->DW_k1_min          = 40;        	// ms
		srf->DW_k2_A            = 109.18;       // ms
		srf->DW_k2_a            = 140.46;       // ms
		srf->DW_k2_k            = 10.7;
		srf->DW_k2_min          = 70.422;
	}
	else
	{
		printf("ERROR: %s is not a valid Spontaneous Release Function Pset for the 3D_cell model fit implementation\n", srf->Pset);
		exit(1);
	}
}
// End Dynamic model - fit to 3D cell==============================//|

// Dynamic model - general/controllable ============================\\|
void set_SRF_dynamic_parameters_controllable(Spontaneous_release_functions *srf, Argument_parameters A)
{
	// Global
	// RyRo vs duration params
	srf->NRyRo_peak_median_A    = 692.99;
	srf->NRyRo_peak_median_min  = 0.059;
	srf->NRyRo_peak_median_H    = -1.6;

	if (strcmp(srf->Pset, "User_control") == 0) // Assign all parameters from user settings
	{
		if (A.SRF_N_Dyn_set == 12)
		{
			srf->PSCR_threshold         = A.SRF_PSCR_threshold;
			srf->CaSR_max               = A.SRF_CaSR_max;
			srf->CaSR_Prange            = A.SRF_CaSR_Prange;
			srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
			srf->ti_sep_max             = A.SRF_ti_sep_max;
			srf->ti_sep_min             = A.SRF_ti_sep_min;
			srf->ti_width_max           = A.SRF_ti_width_max;
			srf->ti_width_min           = A.SRF_ti_width_min;
			srf->MD_max                 = A.SRF_MD_max;
			srf->MD_min                 = A.SRF_MD_min;
			srf->duration_width_max     = A.SRF_duration_width_max;
			srf->duration_width_min     = A.SRF_duration_width_min;
			srf->CaSR_width_H           = A.SRF_CaSR_width_H;
		}
		else
		{
			printf("ERROR: User_control Pset has been set for General Dynamic SRF model, but the wrong number of parameters have bene set;  %d vs the required 12\n", A.SRF_N_Dyn_set);
			printf("\tFor this model, you must pass in all of: SRF_Dyn_{PSCR_threshold/CaSR_max/CaSR_Prange/ti_sep_max/ti_sep_min/ti_width_max/ti_width_min/MD_max/MD_min/duration_width_max/duration_width_min/H}\n");
			exit(1);
		}
	}

	else if (strcmp(srf->Pset, "General_1p10") == 0) // similar to dynamic fit model
	{
		srf->PSCR_threshold 		= 1.1;
		srf->CaSR_max  				= 1.7;
		srf->CaSR_Prange			= 0.05;
		srf->CaSR_min  				= srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max  			= 870;
		srf->ti_sep_min   			= 30;
		srf->ti_width_max			= 1000;
		srf->ti_width_min  			= 20;
		srf->MD_max  				= 800;
		srf->MD_min  				= 50;
		srf->duration_width_max 	= 300;
		srf->duration_width_min 	= 20;
		srf->CaSR_width_H			= 2.5;
	}
	else if (strcmp(srf->Pset, "General_1p05") == 0) // General 1, lower CaSR threshold (-0.05 mM)
	{
		srf->PSCR_threshold         = 1.05;
		srf->CaSR_max               = 1.65;
		srf->CaSR_Prange            = 0.05;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 870;
		srf->ti_sep_min             = 30;
		srf->ti_width_max           = 1000;
		srf->ti_width_min           = 20;
		srf->MD_max                 = 800;
		srf->MD_min                 = 50;
		srf->duration_width_max     = 300;
		srf->duration_width_min     = 20;
		srf->CaSR_width_H           = 2.5;
	}
	else if (strcmp(srf->Pset, "General_1p00") == 0) // General 1, lower CaSR threshold (-0.1 mM)
	{
		srf->PSCR_threshold         = 1.00;
		srf->CaSR_max               = 1.60;
		srf->CaSR_Prange            = 0.05;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 870;
		srf->ti_sep_min             = 30;
		srf->ti_width_max           = 1000;
		srf->ti_width_min           = 20;
		srf->MD_max                 = 800;
		srf->MD_min                 = 50;
		srf->duration_width_max     = 300;
		srf->duration_width_min     = 20;
		srf->CaSR_width_H           = 2.5;
	}
	else if (strcmp(srf->Pset, "General_1p15") == 0) // General 1, higher CaSR threshold (+0.05 mM)
	{
		srf->PSCR_threshold         = 1.15;
		srf->CaSR_max               = 1.75;
		srf->CaSR_Prange            = 0.05;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 870;
		srf->ti_sep_min             = 30;
		srf->ti_width_max           = 1000;
		srf->ti_width_min           = 20;
		srf->MD_max                 = 800;
		srf->MD_min                 = 50;
		srf->duration_width_max     = 300;
		srf->duration_width_min     = 20;
		srf->CaSR_width_H           = 2.5;
	}
	else if (strcmp(srf->Pset, "General_1p20") == 0) // General 1, higher CaSR threshold (+0.1 mM)
	{
		srf->PSCR_threshold         = 1.2;
		srf->CaSR_max               = 1.8;
		srf->CaSR_Prange            = 0.05;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 870;
		srf->ti_sep_min             = 30;
		srf->ti_width_max           = 1000;
		srf->ti_width_min           = 20;
		srf->MD_max                 = 800;
		srf->MD_min                 = 50;
		srf->duration_width_max     = 300;
		srf->duration_width_min     = 20;
		srf->CaSR_width_H           = 2.5;
	}
	else if (strcmp(srf->Pset, "hAM_GB_native") == 0) // Example set to work with hAM_GB native model
	{
		srf->PSCR_threshold         = 0.6;
		srf->CaSR_max               = 0.7;
		srf->CaSR_Prange            = 0.01;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 1000;
		srf->ti_sep_min             = 200;
		srf->ti_width_max           = 800;
		srf->ti_width_min           = 150;
		srf->MD_max                 = 800;
		srf->MD_min                 = 100;
		srf->duration_width_max     = 300;
		srf->duration_width_min     = 50;
		srf->CaSR_width_H           = 2.5;
	}
	else if (strcmp(srf->Pset, "rbAM_ITO_DC") == 0) // Example parameterised to data of Workman et al. 2012 Ito dynamic clamp
	{
		srf->PSCR_threshold         = 0.72;
		srf->CaSR_max               = 1.07;
		srf->CaSR_Prange            = 0.01;
		srf->CaSR_min               = srf->PSCR_threshold - srf->CaSR_Prange;
		srf->ti_sep_max             = 1340-400;
		srf->ti_sep_min             = 445-400;
		srf->ti_width_max           = 262;
		srf->ti_width_min           = 8;
		srf->MD_max                 = 250;
		srf->MD_min                 = 75;
		srf->duration_width_max     = 150;
		srf->duration_width_min     = 15;
		srf->CaSR_width_H           = 2.5;
	}
	else
	{
		printf("ERROR: %s is not a valid Spontaneous Release Function Pset for the General, dynamic controllable implementation\n", srf->Pset);
		exit(1);
	}		
}
// End Dynamic model - general/controllable ========================//|
// End Parameter sets =======================================================//|

// Dynamic model functions which determine the parameters for distributions =\\|
// 3D cell fit functions ==========================================\\|
void determine_SRF_3D_cell_PSCR(Spontaneous_release_functions *srf, double CaSR)
{
	srf->PSCRE		= 1.0/(1 + exp(-(CaSR - srf->PSCR_threshold)/srf->PSCR_k));
}

void determine_SRF_3D_cell_ti_sep(Spontaneous_release_functions *srf, double CaSR)
{
	// ti sep
	srf->ti_sep			= srf->ti_sep_A * exp(-(CaSR - srf->ti_sep_a)/srf->ti_sep_k) + srf->ti_sep_min;
	// impose maximal constraint  -> investigate this??? why did I put this here?
	if (srf->ti_sep > 877) srf->ti_sep = 877;

	// CF at that point
	srf->CF_ti_sep      = srf->CF_ti_sep_A/(1 + exp(-(CaSR - srf->CF_ti_sep_a)/srf->CF_ti_sep_k)) + 0.05;
}

void determine_SRF_3D_cell_ti_dist_widths(Spontaneous_release_functions *srf, double CaSR)
{
	srf->k_ti_F1        = srf->KF1_A * exp(-(CaSR - srf->KF1_a)/srf->KF1_k) + srf->KF1_min; // min ADDED since data
	srf->k_ti_F2        = srf->KF2_A * pow(CaSR, srf->KF2_H1) + srf->KF2_B * pow(CaSR, srf->KF2_H2) + srf->KF2_min;
	if (srf->k_ti_F1 > 80) 	srf->k_ti_F1 = 80;
	if (srf->k_ti_F2 > 450) srf->k_ti_F2 = 450;
}

void determine_SRF_3D_cell_MD(Spontaneous_release_functions *srf, double CaSR)
{
	srf->MD   = srf->MD_A * exp(-(CaSR - srf->MD_a)/srf->MD_k) + srf->MD_min;
	if (srf->MD > 800) srf->MD = 800;
}
// End 3D cell fit functions ======================================//|

// Dynamic, controllable ==========================================\\|
void determine_SRF_dynamic_general_PSCR(Spontaneous_release_functions *srf, double CaSR)
{
	srf->PSCRE      = 1.0/(1 + exp(-(CaSR - srf->PSCR_threshold)/(srf->CaSR_Prange*0.1)));
}

void determine_SRF_dynamic_general_ti_sep(Spontaneous_release_functions *srf, double CaSR)
{
	srf->ti_sep         = (srf->ti_sep_max - srf->ti_sep_min)*exp(- (5*(CaSR - srf->CaSR_min)/(srf->CaSR_max - srf->CaSR_min)))  + srf->ti_sep_min;
	srf->CF_ti_sep		= 0.4; // reasonable for most disitrbutions actually observed
}

void determine_SRF_dynamic_general_ti_dist_widths(Spontaneous_release_functions *srf, double CaSR)
{
	double width_F1;		// width of function 1, ms -> gradients derived from this
	width_F1	=	(srf->ti_width_max - srf->ti_width_min)*pow((srf->ti_sep - srf->ti_sep_min)/(srf->ti_sep_max - srf->ti_sep_min), srf->CaSR_width_H) + srf->ti_width_min;

	// ms -> gradient parameter
	srf->k_ti_F1	= 0.145*width_F1;
	srf->k_ti_F2	= 1.5 * srf->k_ti_F1;	// Preserved ti_sep ~= peak of hist, with CF_ti_sep = 0.4
}

void determine_dynamic_general_cell_MD(Spontaneous_release_functions *srf, double CaSR)
{
	srf->MD = (srf->MD_max - srf->MD_min)*exp(- (5*(CaSR - srf->CaSR_min)/(srf->CaSR_max - srf->CaSR_min)))  + srf->MD_min;
	if (srf->MD > srf->MD_max) srf->MD = srf->MD_max;
}

void determine_SRF_dynamic_general_duration_dist_widths(Spontaneous_release_functions *srf, double CaSR)
{
	// this is in ms -> converted to ks in run function
	srf->duration_width = (srf->duration_width_max - srf->duration_width_min)*pow( (srf->MD - srf->MD_min)/(srf->MD_max - srf->MD_min)  , srf->CaSR_width_H) + srf->duration_width_min;
}
// End Dynamic, controllable ======================================//|

// Function select function
void determine_SRF_params_from_CaSR(Spontaneous_release_functions *srf, double CaSR)
{
	if (strcmp(srf->Model, "3D_cell") == 0)
	{
		determine_SRF_3D_cell_ti_sep(srf, CaSR-0.05);		// shift of 0.05 gives better fit for paced (vs Ca clamp on which derived)
		determine_SRF_3D_cell_ti_dist_widths(srf, CaSR-0.05);
		determine_SRF_3D_cell_MD(srf, CaSR-0.05);
	}
	else // Dynamic, controllable 
	{
		determine_SRF_dynamic_general_ti_sep(srf, CaSR);
		determine_SRF_dynamic_general_ti_dist_widths(srf, CaSR);
		determine_dynamic_general_cell_MD(srf, CaSR);
		determine_SRF_dynamic_general_duration_dist_widths(srf, CaSR);
	}
}

// Finally, recalculate conditions
void need_to_recalculate(Spontaneous_release_functions *srf, double CaSR, double threshold, double NRyR)
{
	if (NRyR < 3e-5)
	{
		if ( (CaSR - srf->CaSR_t_calc) > threshold || (srf->CaSR_t_calc - CaSR) > threshold)
		{
			srf->srf_set = 0;
		}
	}	
}
// End Dynamic model functions which determine disttribution parameters  ====//|

// run SRF functions ========================================================\\|
// Set params and switches ====================\\|
void set_and_run_SRF(Spontaneous_release_functions *srf, Dyad_variables *d, const char * Mode, RAND *rand, int ex_switch, double sim_time, double CaJSR)
{
	// Direct_Control model ================\\|
	if (strcmp(Mode, "Direct_Control") == 0)
	{
		if (srf->srf_set == 0) // if ready to be set
		{
			srf->rand[0] = rand->mtrand1();
			if (srf->rand[0] < srf->PSCRE)	// Only need to produce other rands if event will happen
			{
				srf->srf_calc	=	1;
				for (int j = 0; j < 5; j++) srf->rand[j] = rand->mtrand1();	
			}
			else
			{
				srf->srf_calc 	= 0;
				srf->srf_set	= -1; // Not ready to be set until next time (or it will check against probabilty every dt)
			}
		}
		run_SRF(srf, d, sim_time);
	} 
	// end Direct_Control model ============//|
	// Dynamic model ================\\|
	else if (strcmp(Mode, "Dynamic") == 0)
	{
		if (srf->srf_set == 1 && ex_switch == 1)	// if it has been set, but cell is currently in excitation state (could be spont AP)
		{
			srf->NRyRo			= 0;				// now in excitation mode so don't want SRF to contribute to NRyRo
			srf->ti				= 100000;			// Long so it doesn't do anything until actually been set!
			srf->srf_set		= 0;				// ready to be set again for next time 
		}
		if (srf->srf_set == 0 && ex_switch == 0) 	// not been set and not in excitation state
		{
			srf->CaSR_t_calc	= CaJSR;			// set the CaSR_t_calc to current CaJSR
			srf->rand[0] = rand->mtrand1();
			// Set the probability of SCRE from CaSR
			if (strcmp(srf->Model, "3D_cell") == 0)			determine_SRF_3D_cell_PSCR(srf, 1e-3*CaJSR); 	// CaSR in mM
			else if (strcmp(srf->Model, "General") == 0)	determine_SRF_dynamic_general_PSCR(srf, 1e-3*CaJSR);
			if (srf->rand[0] < srf->PSCRE)	// Only need to produce other rands if event will happen
			{
				srf->srf_calc   =   1;
				for (int j = 0; j < 5; j++) srf->rand[j] = rand->mtrand1();
				determine_SRF_params_from_CaSR(srf, 1e-3*CaJSR);	// Set dist params from CaSR	
			}
			else
			{
				srf->srf_calc   = 0;
				srf->srf_set    = -2;	// now need -2, or on next timestep will attempt to set again -> this tells it to wait for large enough change in CaSR
			}
		}
		run_SRF(srf, d, sim_time);
		if (srf->srf_set == -1 && d->Mi < 0.2 && ex_switch == 0) srf->srf_set = 0;  // if recovered and not in excitation mode, prepare to be set again
		if (srf->srf_set == -2) need_to_recalculate(srf, CaJSR, srf->recalc_SR_diff /*micro M*/,  srf->NRyRo);   // if CaSR has changed by more than threshold since set, re-set
	}
    else if (strcmp(Mode, "Read") == 0)
    {
        run_SRF(srf, d, sim_time); // if read, only need to run, as parameters have already been set directly from the file read in
    }
	// end Dynamic model ============//|
}
// End set params and switches ================//|

// Actual run SRF =============================\\|
void run_SRF(Spontaneous_release_functions *srf, Dyad_variables *d, double sim_time)
{
	double t_to_peak_rand_factor;
	srf->NRyRo	= 0;

	// if function ready to be set but not already been set
	if (srf->srf_set == 0 && srf->srf_calc == 1)
	{
		srf->tset	= sim_time; 	// store sim time at which function is set

		// Set waveform parameters
		// ti
		Determine_ti(srf, srf->rand[0]);

		// duration
		if (srf->duration_width == -10) Determine_duration_ks_from_MD(srf);                           // if default value, determine duration dist from average
		else                            srf->k_D_F1 = srf->k_D_F2 = 1.8*0.145*srf->duration_width;    // else set k from width in ms
		Determine_duration(srf, srf->rand[1]);

		// t to peak and decay time
		t_to_peak_rand_factor = srf->t_to_peak_dist_scale*(srf->rand[2]) + (1-srf->t_to_peak_dist_scale)*0.5; // dist_scale = 0 returns 0.5, = 1 retrurns rand (0-1), intermediary returns rand within smaller range (0.5+/- C; C < 0.5)
		srf->t_to_peak	= 24.0 + t_to_peak_rand_factor*(srf->duration - 52);
		srf->decay_time = srf->duration - srf->t_to_peak;

		// NRyRo_peak
		Determine_NRyRo_peak(srf, srf->duration, srf->rand[3]);

		// Long waveform adjustments
		if (srf->duration >= 300)
		{
			Determine_NRyRo_plateau(srf, srf->duration, srf->rand[4]);
			srf->NRyRo_peak = srf->NRyRo_peak - srf->NRyRo_plateau;	// this is now the spike on top of plateau
			srf->ti_plateau = srf->ti;								// start of plateau is whole waveform so want original ti
			srf->ti         = srf->t_to_peak + srf->ti_plateau;		// this is now ti_spike
			srf->t_to_peak  = 50;									// of spike
			srf->decay_time = 35;									// of spike
			srf->ti         = srf->ti - srf->t_to_peak;				// ti of spike (modded from prev calculated ti)
		}
		else srf->NRyRo_plateau = 0; // for outputs to be able to retain original NRyRo_peak

		// Set actual waveform parameters
		Determine_waveform_parameters(srf);
		if (srf->duration >= 300) Determine_waveform_parameters_long(srf);

		srf->srf_set = 1; // indicates that waveform paramteers have been set
	} // end function ready to be set IF

	// set NRyRo from waveform
	if (srf->srf_set == 1)
	{
		Waveform(srf, sim_time);
		if (srf->duration >= 300) Waveform_plateau(srf, sim_time);
	}
	else srf->NRyRo	= 0.0;

	// Check if waveform has actually intiated a spontaneous release
	if (srf->srf_set == 1 && srf->waveform_init == 0 && srf->NRyRo > 0.0002)
	{
		srf->waveform_init 	= 1;
		srf->tinit			= sim_time;
        srf->init_write_flag    = true;     // it has been set, so needs to be written
	}

	// Calculate approximate proportion of active CRUs
    // ??? Add a maxmimum time, beyond which made re-available? 
	if (srf->waveform_init == 1)
	{
		srf->SRF_prop_active 	= (sim_time - srf->tinit)/(srf->duration);
		if (srf->SRF_prop_active > 1.0) srf->SRF_prop_active = 1.0;
	}
	else srf->SRF_prop_active	= 0;

	// Check to see if wave has finished and reset if so
	if (srf->NRyRo < 0.0002 && srf->waveform_init == 1)
	{
		srf->srf_set                = -1; // not set but not ready to be set again
		srf->waveform_init    	 	= 0;
		srf->srf_calc               = 0;
	}

	// set the open RyR seen by dyad functions appropriately
	if (srf->NRyRo > 0) d->NRyR_O_SRF = srf->NRyRo;
	else                d->NRyR_O_SRF = 0;
}

// Print SRF properties to file
void print_SRF_properties_to_file(Spontaneous_release_functions *srf, std::ostream& out, int n)
{
    if (srf->init_write_flag == true)
    {
        //if (srf->duration < 300) out<<srf->tset + srf->ti<<"   "<<n<<"    "<<srf->NRyRo_peak<<"  "<<srf->duration<<" "<<srf->t_to_peak+srf->tinit<<std::endl;
        //else out<<srf->tset + srf->ti_plateau<<"   "<<n<<"    "<<srf->NRyRo_peak+srf->NRyRo_plateau<<"  "<<srf->duration<<" "<<srf->ti+srf->tinit<<std::endl;

        out << n << " " << srf->duration << " " << srf->tset + srf->ti << " " << srf->tset + srf->ti_plateau << " " << srf->NRyRo_peak \
            << " " << srf->NRyRo_plateau << " " << srf->t_to_peak+srf->tinit << " " << srf->ti+srf->tinit \
            << " " << srf->thalf_1 << " " <<  srf->thalf_2 << " " << srf->k1_waveform << " " << srf->k2_waveform \
            << " " << srf->thalf_plateau_1 << " " << srf->thalf_plateau_2 << " " << srf->k1_plateau << " " << srf->k2_plateau << std::endl;

        srf->init_write_flag = false; // has been written so no longer needs to be
    }
}

// Read SRF settings from a file for reproducing given simulations
void read_SRF_settings_from_file(Spontaneous_release_functions *srf, const char *input_file)
{
    FILE *in;

    in = fopen(input_file, "r");

    if (!in){
        printf("File %s not found. Exiting\n", input_file);
        exit(1);
    }
    else printf("File %s has been read\n", input_file);
    
    int n; // cell identifier
    float temp_float;
    float duration, NRyRo_peak, NRyRo_plateau, thalf_1, thalf_2, k1_waveform, k2_waveform;
    float thalf_plateau_1, thalf_plateau_2, k1_plateau, k2_plateau;

    while (!feof(in))
    {
        fscanf(in, "%d ", &n);
        
        //if srf[n] set is set then read into temps only as only want first SCRE for a given n. Otherwise, read in proper
        if (srf[n].srf_set == 1) for (int i = 0; i < 15; i++) fscanf(in, "%f ", &temp_float);
        else
        {
            fscanf(in, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &duration, &temp_float, &temp_float, &NRyRo_peak, \
                    &NRyRo_plateau, &temp_float, &temp_float, &thalf_1, &thalf_2, &k1_waveform, \
                    &k2_waveform, &thalf_plateau_1, &thalf_plateau_2, &k1_plateau, &k2_plateau);
            
            srf[n].duration         = duration;
            srf[n].NRyRo_peak       = NRyRo_peak;
            srf[n].NRyRo_plateau    = NRyRo_plateau;
            srf[n].thalf_1          = thalf_1;
            srf[n].thalf_2          = thalf_2;
            srf[n].k1_waveform      = k1_waveform;
            srf[n].k2_waveform      = k2_waveform;
            srf[n].thalf_plateau_1  = thalf_plateau_1;
            srf[n].thalf_plateau_2  = thalf_plateau_2;
            srf[n].k1_plateau       = k1_plateau;
            srf[n].k2_plateau       = k2_plateau;

            srf[n].Mode             = "Read"; // and set mode to read for node n
        }

        // note that this node has now had its parameters set
        srf[n].srf_set = 1;

        //printf("%d %f %f %f %f %f %f %f %f %f %f %f\n", n, srf[n].duration, srf[n].NRyRo_peak, srf[n].NRyRo_plateau, srf[n].thalf_1, srf[n].thalf_2, srf[n].k1_waveform, srf[n].k2_waveform, srf[n].thalf_plateau_1, srf[n].thalf_plateau_2, srf[n].k1_plateau, srf[n].k2_plateau);
    }

    fclose(in);
}

// End actual run SRF =========================//|

// NCX and Krel multipliers ===================\\|
void calc_SRF_mults(Spontaneous_release_functions *srf, Membrane_fluxes *mem, Dyad_variables *dyad)
{
    // This is due to a global RyR waveform not inducing same behaviour
    // as the local-propagation which produces the waveform; too much
    // release and too large an INaCa -> where longer duration SRF
    // have the largest difference. Hence, reduce both based on
    // SRF duration (derived by using exact 3D model waveforms in 0D
    // model and adjusting until Ca_ds, Ca_Ss and INaCa_SS were the same
    // then approximating values vs duration
    if (srf->NRyRo > 0.002)
    {
        mem->NCX_SRF_mult     =  0.4/(1 + exp((srf->duration - 176)/46)) + 0.25;
        dyad->Krel_SRF_mult   =  0.35/(1 + exp((srf->duration - 113)/18)) + 0.4;
    }
    else
    {
        mem->NCX_SRF_mult  	= 1.0;
        dyad->Krel_SRF_mult = 1.0;
    }
}
// End NCX and Krel multipliers ===============//|
// End run SRF functions ====================================================//|

// Actual waveform functions ================================================\\|
void Determine_waveform_parameters(Spontaneous_release_functions *srf)
{
    srf->k1_waveform    = 0.16980607*srf->t_to_peak  +  0.00254852;
    srf->k2_waveform    = 0.16980607*srf->decay_time +  0.00254852;

    srf->thalf_1 = srf->tset + srf->ti + (srf->t_to_peak/2);
    srf->thalf_2 = srf->tset + srf->ti + srf->t_to_peak + (srf->decay_time/2);
}

void Determine_waveform_parameters_long(Spontaneous_release_functions *srf)
{
    double 	tau_plateau = 50; 	// ms -> for rise and decay of the plateau

    srf->k1_plateau		= 0.16980607*tau_plateau +  0.00254852; 
    srf->k2_plateau		= 0.16980607*tau_plateau +  0.00254852;

    srf->thalf_plateau_1	=  srf->tset + srf->ti_plateau + (tau_plateau/2);
    srf->thalf_plateau_2	=  srf->tset + srf->ti_plateau + tau_plateau + (srf->duration-2*tau_plateau)+ (tau_plateau/2);
}

void Waveform(Spontaneous_release_functions *srf, double t)
{
    srf->Fn1    	= 1/(1 + exp(-(t - srf->thalf_1)/srf->k1_waveform) );
    srf->Fn2    	= 1/(1 + exp( (t - srf->thalf_2)/srf->k2_waveform) );

    srf->NRyRo      = srf->NRyRo_peak * srf->Fn1 * srf->Fn2;
}

void Waveform_plateau(Spontaneous_release_functions *srf, double t)
{
    srf->Fn3		= 1/(1 + exp(-(t - srf->thalf_plateau_1)/srf->k1_plateau) );
    srf->Fn4	 	= 1/(1 + exp((t - srf->thalf_plateau_2)/srf->k2_plateau) );

    // If plateau, the normal waveform params (above) sets the spike; this sets the plateau on which
    // the spike occurs. Thus, need to add to the already set NRyRo = f(time) set above
    srf->NRyRo 		+= (srf->NRyRo_plateau *  srf->Fn3 * srf->Fn4);
}
// End Actual waveform functions ============================================//|

// Distribution functions and their inverses ================================\\|
// Initiation time, ti ========================\\|
void Determine_ti(Spontaneous_release_functions *srf, double rand)
{
    if (rand < 0.0002) 				srf->ti = -srf->k_ti_F1*log(( 2 * (  srf->CF_ti_sep) /  0.0002                         ) - 1) + srf->ti_sep;
    else if (rand > 0.9998)			srf->ti = -srf->k_ti_F2*log(((2 * (1-srf->CF_ti_sep))/((0.9998 + 1) - 2*srf->CF_ti_sep)) - 1) + srf->ti_sep;
    else
    {
        if (rand < srf->CF_ti_sep)	srf->ti = -srf->k_ti_F1*log(( 2 * (  srf->CF_ti_sep) /  rand                           ) - 1) + srf->ti_sep;
        else 						srf->ti	= -srf->k_ti_F2*log(((2 * (1-srf->CF_ti_sep))/((rand   + 1) - 2*srf->CF_ti_sep)) - 1) + srf->ti_sep;
    }
    if (srf->ti < 0.0) srf->ti = 0.0; // maybe print error here?
}
// End Initiation time, ti ====================//|

// Duration ===================================\\|
void Determine_duration(Spontaneous_release_functions *srf, double rand)
{
    if (rand < 0.0002)      srf->duration =  -srf->k_D_F1*log(1.0/0.0002 - 1)	+ srf->MD; // extreme clauses
    else if (rand > 0.9998) srf->duration =  -srf->k_D_F2*log(1.0/0.9998 - 1) 	+ srf->MD;
    else if (rand < 0.5)    srf->duration =  -srf->k_D_F1*log(1.0/rand   - 1) 	+ srf->MD; // general population
    else                    srf->duration =  -srf->k_D_F2*log(1.0/rand   - 1) 	+ srf->MD;

    // Ensuring not outside expected ranges
    if (srf->duration < 20) 	srf->duration = 1500;//printf("negative duration %f\n", srf->duration);  // just to ensure we don't get too many short waves; should be a rare case, this ensures it has no real effect
    if (srf->duration > 1500) 	srf->duration = 1500; // maximum physiological duration

    //if (srf->duration < srf->MD - (srf->k_D_F1/(1.8*0.145)) - 10 )
    //{
	//	srf->duration = srf->MD -(srf->k_D_F1/(1.8*0.145) - 10);
	//printf("duration outside range catch %f %f\n", (srf->k1_duration/(1.8*0.145)), (srf->k1_duration/(1.8*0.145)-10));
	//}
}

// for when duration width as a function of median, not pre-defined
void Determine_duration_ks_from_MD(Spontaneous_release_functions *srf)
{
	srf->k_D_F1		= srf->DW_k1_A/(1 + exp(-((srf->MD - srf->DW_k1_a))/srf->DW_k1_k)) + srf->DW_k1_min;
	srf->k_D_F2		= srf->DW_k2_A/(1 + exp(-((srf->MD - srf->DW_k2_a))/srf->DW_k2_k)) + srf->DW_k2_min;

	// ms -> gradient param conversion
	srf->k_D_F1		*= 1.8*0.145;
	srf->k_D_F2		*= 1.8*0.145;
}
// End Duration ===============================//|

// NRyR peak ==================================\\|
void Determine_NRyRo_peak(Spontaneous_release_functions *srf, double duration, double rand)
{
	srf->NRyRo_peak_median	= srf->NRyRo_peak_median_A * pow(duration, srf->NRyRo_peak_median_H) + srf->NRyRo_peak_median_min;
	srf->NRyRo_peak			= srf->NRyRo_peak_median + (rand - 0.5)*0.05*srf->NRyRo_peak_dist_scale;	// uniform variance approximation; reasonable ; 0.05*dist_scale is maximum difference from median
	if (srf->NRyRo_peak	> 0.9) srf->NRyRo_peak = 0.9; // maximal value
}

void Determine_NRyRo_plateau(Spontaneous_release_functions *srf, double duration, double rand)
{
	srf->NRyRo_plateau		= 31.09*pow(0.01*duration, -7.39) + 0.034;
	srf->NRyRo_plateau		+= srf->NRyRo_plateau_dist_scale*(rand-0.5)*(0.25*srf->NRyRo_plateau); // dist_scale *0.25*NRyR_plateau is maximum variance
}
// End NRyR peak ==============================//|
// End Distribution functions and their inverses ============================//|

// Test and output distributions ============================================\\|
void test_and_produce_distributions(Spontaneous_release_functions *srf, const char * directory, RAND *r, double CaSR)
{
	FILE * out;
	char str[1000];

	double rand;

	// Initiation time distribution ===========\\|
	int v_hist[500];
	int v_int;
	// First, produce inverse function
	if (CaSR == 0.0) sprintf(str, "%s/SRF_distributions/ti_inverse_function.dat", directory);
	else sprintf(str, "%s/SRF_distributions/ti_inverse_function_%0.2fmM.dat", directory, CaSR);
	out = fopen(str, "wt");
	for (double x = 0; x <= 1.0; x+= 0.01)
	{
		Determine_ti(srf, x);	
		fprintf(out, "%f %f\n", srf->ti, x);
	}
	fclose(out);
	// Now distribtion from random sampling
	for (int i = 0; i < 500; i++) v_hist[i] = 0;
	for (int i = 0; i < 10000; i++)
	{
		rand = r->mtrand1();
		Determine_ti(srf, rand);
		v_int = int(srf->ti)/10;
		v_hist[v_int] ++; // ti intervals of 10 ms	
	}
	if (CaSR == 0.0) sprintf(str, "%s/SRF_distributions/ti_distribution.dat", directory);
	else sprintf(str, "%s/SRF_distributions/ti_distribution_%0.2fmM.dat", directory, CaSR);
	out = fopen(str, "wt");
	for (int i = 0; i < 500; i++) fprintf(out, "%d %d\n", (i*10)+5, v_hist[i]);
	fclose(out);
	// End Initiation time distribution =======//|

	// Duration distribtion ===================\\|
	for (int i = 0; i < 500; i++) v_hist[i] = 0;
	if (srf->duration_width == -10) Determine_duration_ks_from_MD(srf);    									// if default value, determine duration dist from average
	else                            srf->k_D_F1 = srf->k_D_F2 = 1.8*0.145*srf->duration_width;    // else set k from width in ms
	// First, produce inverse function
	if (CaSR == 0.0) sprintf(str, "%s/SRF_distributions/duration_inverse_function.dat", directory);
	else sprintf(str, "%s/SRF_distributions/duration_inverse_function_%0.2fmM.dat", directory, CaSR);
	out = fopen(str, "wt");
	for (double x = 0; x <= 1.0; x+= 0.01)
	{
		Determine_duration(srf, x);
		fprintf(out, "%f %f\n", srf->duration, x);
	}
	fclose(out);
	// Now distribtion from random sampling
	for (int i = 0; i < 10000; i++)
	{
		rand = r->mtrand1();
		Determine_duration(srf, rand);
		v_int = int(srf->duration)/10;
		v_hist[v_int] ++; // ti intervals of 10 ms
	}
	if (CaSR == 0.0) sprintf(str, "%s/SRF_distributions/duration_distribution.dat", directory);
	else sprintf(str, "%s/SRF_distributions/duration_distribution_%0.2fmM.dat", directory, CaSR);
	out = fopen(str, "wt");
	for (int i = 0; i < 500; i++) fprintf(out, "%d %d\n", (i*10)+5, v_hist[i]);
	fclose(out);	
	// End Duration distribtion ===============//|

	// NRyRo_peak =============================\\|
	// Produce the NRyRo_peak vs duration function
	sprintf(str, "%s/SRF_distributions/NRyRo_peak_median_vs_duration_function.dat", directory);
	out = fopen(str, "wt");	
	for (double x = 0; x <= 1000; x+= 1)
	{
		Determine_NRyRo_peak(srf, x, 0.5);
		fprintf(out, "%f %f\n", x, srf->NRyRo_peak); // NRyRo_peak = median as rand = 0.5
	}
	fclose(out);
	for (int i = 0; i < 500; i++) v_hist[i] = 0;
	// Now distribtion from random sampling (this dist depends on two random samples - duration and NRyR
	for (int i = 0; i < 10000; i++)
	{
		rand = r->mtrand1();
		Determine_duration(srf, rand); // First, sample duration
		rand = r->mtrand1();
		Determine_NRyRo_peak(srf, srf->duration, rand); // and now sample NRyRo
		v_int = int(10*srf->NRyRo_peak);
		v_hist[v_int] ++; // ti intervals of 10 ms
	}
	if (CaSR == 0.0) sprintf(str, "%s/SRF_distributions/NRyRo_peak_distribution.dat", directory);
	else sprintf(str, "%s/SRF_distributions/NRyRo_peak_distribution_%0.2fmM.dat", directory, CaSR);
	out = fopen(str, "wt");
	for (int i = 0; i < 10; i++) fprintf(out, "%f %d\n", float(i)/10, v_hist[i]);
	fclose(out);
	// End NRyRo_peak =========================//|
}

// Below is for Dynamic models, to produce dependency on CaSR and dists at various CaSR
void test_and_produce_CaSR_dependency(Spontaneous_release_functions *srf, const char * directory, RAND *r)
{
	FILE * out;
    FILE * out2;
	char str[1000];
    char str2[1000];

	double rand, CaJSR;
    int CaSR_int;

	// Plot how aves and widths vary with CaSR  and at 0.1mM intervals, calculate inverses and actual distributions 
	sprintf(str, "%s/SRF_distributions/CaSR_dependence.dat", directory);
	out = fopen(str, "wt");
	//SRF_setup(srf, A);		// Define the dist parameters based on 3D_cell, General and specific Pset

	for (CaJSR = 500; CaJSR < 2000; CaJSR += 1)
	{
		if (strcmp(srf->Model, "3D_cell") == 0)         determine_SRF_3D_cell_PSCR(srf, 1e-3*CaJSR);    // CaSR in mM
		else if (strcmp(srf->Model, "General") == 0)    determine_SRF_dynamic_general_PSCR(srf, 1e-3*CaJSR);
		determine_SRF_params_from_CaSR(srf, 1e-3*CaJSR);

		if (srf->duration_width == -10) Determine_duration_ks_from_MD(srf);                           // if default value, determine duration dist from average
		else                            srf->k_D_F1 = srf->k_D_F2 = 1.8*0.145*srf->duration_width;    // else set k from width in ms
		if (strcmp(srf->Model, "3D_cell") == 0) srf->CaSR_Prange = 10*srf->PSCR_k; // for General Prange is what is already defined and k = 0.1*Prange

		if (1e-3*CaJSR > srf->PSCR_threshold - srf->CaSR_Prange) fprintf(out, "%f %f %f %f %f %f %f %f\n", 1e-3*CaJSR, srf->PSCRE, srf->ti_sep, srf->k_ti_F1/0.145, srf->k_ti_F2/0.145, srf->MD, srf->k_D_F1/(1.8*0.145), srf->k_D_F2/(1.8*0.145));
		else fprintf(out, "%f %f 0 0 0 0 0 0\n", 1e-3*CaJSR, srf->PSCRE);	

        // write settings files if option is selected
        if (strcmp(srf->write_SRF_settings, "On") == 0)
        {
            CaSR_int = CaJSR;
            sprintf(str2, "%s/SRF_distributions/SRF_Settings_file_%04d.txt", directory, CaSR_int);
            out2 = fopen(str2, "wt");
            fprintf(out2, "9\nSRF_mode Direct_Control\nSRF_Pset User_control\nSRF_DC_PSCRE %f\nSRF_DC_CF 0.4\nSRF_DC_ti_sep %f\nSRF_DC_ti_W1 %f\nSRF_DC_ti_W2 %f\nSRF_DC_MD %f\nSRF_DC_duration_W %f\n", \
                    srf->PSCRE, srf->ti_sep, srf->k_ti_F1/0.145, srf->k_ti_F2/0.145, srf->MD, srf->k_D_F1/(1.8*0.145) );
            fclose(out2);
        }

		// calculate inverses and distributions at CaSR intervals
		if (int(CaJSR)%100 == 0 && 1e-3*CaJSR 	> srf->PSCR_threshold - srf->CaSR_Prange) test_and_produce_distributions(srf, directory, r, 1e-3*CaJSR);
	}
	fclose(out);
}
// End test and output distributions ========================================//|
