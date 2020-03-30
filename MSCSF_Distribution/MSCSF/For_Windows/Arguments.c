// ****REVIEW/BETA version of the code. NOT FINAL. NOT TO BE REDISTRIBUTED*****
// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Command-line and settings-file argument =====  //
// handling. ==============================================  //
// ========================================================  //
// COPYRIGHT MICHAEL A. COLMAN 2015-2019.==================  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND =  //
// CONTRIBUTORS "AS-IS" AND ANY EXPRESS OR IMPLIED ========  //
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED =  //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A ========  //
// PARTICULAR PURPOSE ARE DISCLAIMED. =====================  //
// ========================================================  //
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

#include "Arguments.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h> 

// Functions  =========================================================================\\|
//	set_argument_defaults()
//	call_argument_functions()
//	set_arguments()
// End Functions  =====================================================================//|

// Set Argument flag defaults =========================================================\\|
void set_argument_defaults(Argument_parameters *A)
{
	// Simulation setttings =========\\|
	A->reference_arg	        	= false;
	A->results_reference_arg    	= false;
	A->state_reference_read_arg  	= false;
	A->state_reference_write_arg 	= false;
	A->BCL_arg                  	= false;
	A->Total_time_arg           	= false;
	A->Paced_time_arg           	= false;
	A->NBeats_arg               	= false;
	A->dt_arg                   	= false;
	A->S2_arg			        	= false;
	A->NS2_arg			        	= false;
	A->SOI_arg			        	= false;
	A->SOId_arg			        	= false;
	A->Multi_stim_arg	        	= false;
	A->settings_file            	= false;
	// End sim settings =============//|

	// Model and cell conditions=====\\|
	A->Model_arg        			= false;
	A->Celltype_arg     			= false;
	A->Agent_arg					= false;
	A->Remodelling_arg				= false;
	A->ISO_arg						= false;
	A->ISO_model_arg				= false;
	A->ACh_arg						= false;
	A->ACh_model_arg				= false;
	A->Mutation_arg					= false;
	A->Ihyp_arg						= false;
	A->environment_arg				= false;
	A->Remodelling_prop_arg 		= false;
	A->Agent_prop_arg				= false;
	A->spatial_gradient_arg 		= false;
	A->spatial_gradient_prop_arg 	= false;

	A->Vclamp						= "Off";
	A->Write_state					= "Off";
	A->Read_state					= "Off";
	// End model and cell ===========//|   

	// Tissue settings ==============\\|
	A->Tissue_order_arg			= false;
	A->Tissue_model_arg			= false;
	A->Tissue_type_arg			= false;
	A->Orientation_type_arg		= false;
	A->D_uniformity_arg			= false;
	A->Stimulus_type_arg		= false;
	A->S2_Stimulus_type_arg		= false;
	A->Dscale_arg				= false;
	A->Dscale_mod_map_arg   	= false;
	A->D_AR_scale_arg       	= false;
	A->D_AR_scale_mod_map_arg   = false;

	A->S1_x_loc_arg				= false;
	A->S1_x_size_arg			= false;
	A->S1_y_loc_arg				= false;
	A->S1_y_size_arg			= false;
	A->S1_z_loc_arg				= false;
	A->S1_z_size_arg			= false;
	A->S2_x_loc_arg				= false;
	A->S2_x_size_arg			= false;
	A->S2_y_loc_arg				= false;
	A->S2_y_size_arg			= false;
	A->S2_z_loc_arg				= false;
	A->S2_z_size_arg			= false;

	A->S1_shape_arg				= false;
	A->S2_shape_arg				= false;

	A->D1_arg							= false;
	A->D_AR_arg							= false;
	A->dx_arg							= false;
	A->OX_arg							= false;
	A->OY_arg							= false;
	A->OZ_arg							= false;
	A->Global_orientation_direction_arg = false;


	A->Remodelling_map_arg				= false;
	A->ISO_map_arg						= false;
	A->ACh_map_arg						= false;
	A->SRF_map_arg						= false;
	A->Direct_modulation_map_arg		= false;

	A->map_shape_arg    			    = false;
	A->map_in_type_arg  			    = false;
	A->map_x_loc_arg    			    = false;
	A->map_y_loc_arg    			    = false;
	A->map_z_loc_arg    			    = false;
	A->map_x_size_arg   			    = false;
	A->map_y_size_arg   			    = false;
	A->map_z_size_arg   			    = false;

	A->stim_file_arg					= false;
	A->S2_stim_file_arg					= false;
	A->phase_file_arg					= false;
	A->Dscale_base_map_file_arg			= false;
	A->Dscale_mod_map_file_arg			= false;
	A->D_AR_scale_base_map_file_arg		= false;
	A->D_AR_scale_mod_map_file_arg		= false;
	A->ISO_map_file_arg					= false;
	A->remod_map_file_arg				= false;
	A->ACh_map_file_arg					= false;	
	A->SRF_map_file_arg					= false;
	A->Direct_modulation_map_file_arg	= false;
	A->spatial_gradient_map_file_arg	= false;
	A->Tissue_model_2_arg   			= false;
	A->Multiple_models_arg  			= false;
	// End Tissue settings ==========//|

	// Spatial cell models ==========\\|
	A->Cell_size_arg		= false;
	A->Sim_cell_size_arg	= false;
	A->Cai_IC_arg			= false;
	A->CaSR_IC_arg			= false;
	A->RyR_Po_arg			= false;
	A->LTCC_Po_arg			= false;
	A->Detub_arg			= false;
	A->tau_ss_arg 			= false;
	A->SERCA_het_arg		= false;
	A->NCX_het_arg		    = false;
	A->RyR_het_arg		    = false;
	A->LTCC_het_arg		    = false;
	A->TT_map_file_arg		= false;
	A->SERCA_map_file_arg	= false;
	A->NCX_map_file_arg		= false;
	A->RyR_het_map_file_arg	= false;
	A->LTCC_map_file_arg	= false;
    A->volds_het_arg        = false;

	A->Delayed_CaSR_IC_arg	= false;
	A->CaSR_IC_delay_arg	= false;
	// Spatial cell models ==========//|

	// Spontaneous release functions \\|
	A->SRF_mode_arg			= false;
	A->SRF_model_arg		= false;
	A->SRF_Pset_arg			= false;
	A->SRF_het_arg			= false;
	A->SRF_N_DC_set			= 0;
	A->SRF_N_Dyn_set		= 0;
	// Spontaneous release functions //|

	// Direct control current mod ===\\|
	A->DC_current_mod_arg		= false;	   	
	A->GNa_arg					= false;                 
	A->GNaL_arg					= false;                
	A->Gto_arg					= false;               
	A->GCaL_arg					= false;                
	A->GKur_arg					= false;                
	A->GKr_arg					= false;                 
	A->GKs_arg					= false;                 
	A->GK1_arg					= false;                 
	A->GNCX_arg					= false;                
	A->GCaP_arg					= false;                
	A->GNab_arg					= false;                
	A->GCab_arg					= false;                
	A->GKb_arg					= false;                 
	A->GNaK_arg					= false;                
	A->GClCa_arg				= false;               
	A->INa_va_tau_scale_arg		= false;    
	A->INa_vi_1_tau_scale_arg	= false;  
	A->INa_vi_2_tau_scale_arg	= false;  
	A->INaL_va_tau_scale_arg	= false;   
	A->INaL_vi_tau_scale_arg	= false;   
	A->Ito_va_tau_scale_arg		= false;    
	A->Ito_vi_tau_scale_arg		= false;    
	A->ICaL_va_tau_scale_arg	= false;   
	A->ICaL_vi_tau_scale_arg	= false;   
	A->IKur_va_tau_scale_arg	= false;   
	A->IKur_vi_tau_scale_arg	= false;   
	A->IKr_va_tau_scale_arg		= false;    
	A->IKs_va_tau_scale_arg		= false;    
	A->INa_va_shift_arg			= false;        
	A->INa_vi_shift_arg			= false;        
	A->INaL_va_shift_arg		= false;       
	A->INaL_vi_shift_arg		= false;       
	A->Ito_va_ss_shift_arg		= false;     
	A->Ito_vi_ss_shift_arg		= false;     
	A->Ito_va_tau_shift_arg		= false;    
	A->Ito_vi_tau_shift_arg		= false;    
	A->Ito_va_ss_kscale_arg		= false;    
	A->Ito_vi_ss_kscale_arg		= false;    
	A->ICaL_va_ss_shift_arg		= false;    
	A->ICaL_vi_ss_shift_arg		= false;    
	A->ICaL_va_tau_shift_arg	= false;   
	A->ICaL_vi_tau_shift_arg	= false;   
	A->ICaL_va_ss_kscale_arg	= false;   
	A->ICaL_vi_ss_kscale_arg	= false;   
	A->IKur_va_ss_shift_arg		= false;    
	A->IKur_vi_ss_shift_arg		= false;    
	A->IKur_va_tau_shift_arg	= false;   
	A->IKur_vi_tau_shift_arg	= false;   
	A->IKur_va_ss_kscale_arg	= false;   
	A->IKur_vi_ss_kscale_arg	= false;   
	A->IKr_va_ss_shift_arg		= false;     
	A->IKr_va_tau_shift_arg		= false;    
	A->IKr_va_ss_kscale_arg		= false;    
	A->IKr_vi_ss_shift_arg		= false;     
	A->IKr_vi_ss_kscale_arg		= false;    
	A->IKs_va_ss_shift_arg		= false;     
	A->IKs_va_tau_shift_arg		= false;    
	A->IKs_va_ss_kscale_arg		= false;    
	A->IK1_va_shift_arg			= false;        
	A->IK1_Erev_shift_arg		= false;        
	A->Gup_arg					= false;                 
	A->Gleak_arg				= false;                
	A->Grel_arg					= false;                
	A->Ito_shift_arg			= false;
	A->ICaL_shift_arg			= false;
	A->IKur_shift_arg			= false;
	A->GClb_arg					= false;
	A->GKACh_arg				= false;
	A->IKACh_va_tau_scale_arg	= false;
	A->IKACh_va_ss_shift_arg	= false;
	A->IKACh_va_ss_kscale		= false;
	A->IKACh_va_tau_shift_arg	= false;
	// End Direct control current mod //|
}
// End Set Argument defaults ==========================================================//|

// Call argumebt functions (settings file and command line) ===========================\\|
void call_argument_functions(int argc, char *argin[], Argument_parameters *A, char const *Version)
{
	// Print time and model type to Log file
	FILE *out;
	out = fopen("Log.dat", "a");
	time_t rawtime;
	time (&rawtime);
	fprintf(out, "Sim time: %sModel version ran: %s\nSettings:\t", ctime (&rawtime), Version);

	// Settings file ====================================\\|
	// Settings_file variables
	const char *settings_filename;
	A->Narg_settings = 1;
	FILE *settings_in;
	settings_in = NULL; // todo Jakub - to make sure this is initialised.

	// First, check to see if we want to read in a settings file
	if (argc > 1) if (strcmp(argin[1], "Settings_file") == 0) // check for if first argument is "Settings_file"
	{ 
		settings_filename = argin[2];
		settings_in = fopen(settings_filename, "r");
        if (settings_in == NULL)
        {
            printf("ERROR: Settings file \"%s\" could not be found. Is it in this directory?\n", settings_filename);
            exit(1);
        }
		A->settings_file = true;
		fscanf(settings_in, "%d\n", &A->Narg_settings);
		printf("Settings file %s read in; N_arguments = %d\n", settings_filename, A->Narg_settings);
		A->Narg_settings *= 2; // 2 options for each argument (arg and value)
		A->Narg_settings += 1; // count from 1 for consistency between command line and file-read versions
	}

	int counter;
	if (A->settings_file == true) // if from file, read from file to array to pass into argument function
	{
		for(counter = 1; counter < A->Narg_settings; counter++)
		{
			// I'm sure there is a way to do this without using 2 arrays
			// However, this seemed easiest not to lose information and preserve the array type as same as command-line for set_arguments function
			fscanf(settings_in, "%s ", A->argin_sf_in[counter]);
			A->argin_sf[counter] = A->argin_sf_in[counter];
		}

		fclose(settings_in);
	}
	// End Settings file ================================//|

	// Set arguments ====================================\\|
	// Now pass command line and/or settings files arrays into set arguments function to set Argument_parameters struct
	if (A->settings_file == false) set_arguments(argc, argin, A, Version, out);    	// Sets the arguments struct from  command-line arguments
	else
	{
		set_arguments(A->Narg_settings, A->argin_sf, A, Version, out);    			// Sets the arguments struct from file-arguments
		if (argc > 3) set_arguments(argc, argin, A, Version, out);  				// read in additional arguments from command-line
	}
	// End Set arguments ================================//|

	// Finalise log file
	fprintf(out, "\n\n");
	fclose(out);
}
// End call argumebt functions (settings file and command line) =======================//|

// Set arguments ======================================================================\\|
void set_arguments(int Narg, char *argin[], Argument_parameters *A, char const *Version, FILE *out)
{

	// Structure of argument types in this loop:
	//		1 - Simulation settings (references, BCL, S2, Vclamp etc)
	//		2 - Mdoel and cell conditions (Model, Celltype, ISO etc)
	//		3 - Tissue model arguments (Tissue model, diffusion aprameters, modualtion maps, map files etc)
	//		4 - Spatial cell model options (cellsize, simcellszie, tau_ss_type etc)
	//		5 - Spontaneous release function options
	//		6 - Direct control of current and flux modulation

	// Read arguments (command line or from file) into Argument_parameters struct
	int counter = 1;
	while (counter < Narg) // loop until Number of arguments has been reached
	{
		bool isFound = false;
		// Simulation setttings =============================================\\|
		if (strcmp(argin[counter], "Reference") == 0 || strcmp(argin[counter], "reference") == 0 || strcmp(argin[counter], "Ref") == 0 || strcmp(argin[counter], "ref") == 0) // if there is no difference between argument and "BCL"
		{
			A->reference              	= argin[counter+1];       	// Reads argument into reference variable
			A->reference_arg			= true;						// Argument has been passed
			fprintf(out, "reference   %s ", argin[counter+1]);    	// Print option to log file
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Results_Reference") == 0 || strcmp(argin[counter], "results_reference") == 0 || strcmp(argin[counter], "Results_reference") == 0 || strcmp(argin[counter], "results_ref") == 0 || strcmp(argin[counter], "Results_ref") == 0) // if there is no difference between argument and "BCL"
		{
			A->results_reference                = argin[counter+1];         // Reads argument into reference variable
			A->results_reference_arg            = true;                     // Argument has been passed
			fprintf(out, "results_reference   %s ", argin[counter+1]);      // Print option to log file
			counter++; isFound = true;
		}       
		if (strcmp(argin[counter], "State_Reference_read") == 0 || strcmp(argin[counter], "state_reference_read") == 0 || strcmp(argin[counter], "state_ref_read") == 0 || strcmp(argin[counter], "State_ref_read") == 0) // if there is no difference between argument and "BCL"
		{
			A->state_reference_read                = argin[counter+1];         // Reads argument into reference variable
			A->state_reference_read_arg            = true;                     // Argument has been passed
			fprintf(out, "state_reference_read   %s ", argin[counter+1]);      // Print option to log file
			counter++; isFound = true;
		} 
		if (strcmp(argin[counter], "State_Reference_write") == 0 || strcmp(argin[counter], "state_reference_write") == 0 || strcmp(argin[counter], "state_ref_write") == 0 || strcmp(argin[counter], "State_ref_write") == 0) // if there is no difference between argument and "BCL"
		{
			A->state_reference_write                = argin[counter+1];         // Reads argument into reference variable
			A->state_reference_write_arg            = true;                     // Argument has been passed
			fprintf(out, "state_reference_write   %s ", argin[counter+1]);      // Print option to log file
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Vclamp") == 0)
		{
			A->Vclamp           = argin[counter+1];
			fprintf(out, "Vclamp   %s ", argin[counter+1]);
			if (strcmp(A->Vclamp, "On") != 0 && strcmp(A->Vclamp, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Vclamp argument. Please pass only \"Off\" or \"On\"\n\n", A->Vclamp);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "BCL") == 0)
		{
			A->BCL				= atoi(argin[counter+1]);	
			A->BCL_arg			= true;						
			fprintf(out, "BCL %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Total_time") == 0 || strcmp(argin[counter], "total_time") == 0 || strcmp(argin[counter], "Total_Time") == 0)
		{
			A->Total_time		= atoi(argin[counter+1]); 
			A->Total_time_arg	= true;
			fprintf(out, "Total_time %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Paced_time") == 0 || strcmp(argin[counter], "paced_time") == 0 || strcmp(argin[counter], "Paced_Time") == 0)
		{
			A->Paced_time       = atoi(argin[counter+1]);
			A->Paced_time_arg   = true;
			fprintf(out, "Paced_time %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "NBeats") == 0 || strcmp(argin[counter], "Nbeats") == 0 || strcmp(argin[counter], "Beats") == 0 || strcmp(argin[counter], "beats") == 0)
		{
			A->NBeats       	= atoi(argin[counter+1]);
			A->NBeats_arg   	= true;
			fprintf(out, "NBeats %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "dt") == 0)         
		{
			A->dt              	= atof(argin[counter+1]);       
			A->dt_arg          = true;                         
			fprintf(out, "dt %s ", argin[counter+1]);                
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2") == 0)
		{
			A->S2_CL            = atoi(argin[counter+1]);
			A->S2_arg          = true;
			fprintf(out, "S2   %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "NS2") == 0)
		{
			A->NS2            = atoi(argin[counter+1]);
			A->NS2_arg        = true;
			fprintf(out, "NS2   %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Write_state") == 0)
		{
			A->Write_state           	= argin[counter+1];
			fprintf(out, "Write_state   %s ", argin[counter+1]);
			if (strcmp(A->Write_state, "On") != 0 && strcmp(A->Write_state, "Off") != 0 && strcmp(A->Write_state, "phase") != 0 && strcmp(A->Write_state, "ave") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Read/Write state argument. Please pass only \"Off\" or \"On\" or (tissue only): \"phase\" or \"ave\"\n\n", A->Write_state);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Read_state") == 0)
		{
			A->Read_state      	     	= argin[counter+1];
			fprintf(out, "Read_state   %s ", argin[counter+1]);
			if (strcmp(A->Read_state, "On") != 0 && strcmp(A->Read_state, "Off") != 0 && strcmp(A->Read_state, "phase") != 0 && strcmp(A->Read_state, "ave") != 0 && strcmp(A->Read_state, "single_cell") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Read/Write state argument. Please pass only \"Off\" or \"On\" or (tissue only): \"phase\" , \"ave\" or \"single_cell\"\n\n", A->Read_state);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Spatial_output_interval_vtk") == 0)
		{
			A->SOI            = atoi(argin[counter+1]);
			A->SOI_arg        = true;
			fprintf(out, "Spatial_output_interval_vtk   %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Spatial_output_interval_data") == 0)
		{
			A->SOId            = atoi(argin[counter+1]);
			A->SOId_arg        = true;
			fprintf(out, "Spatial_output_interval_data   %s ", argin[counter+1]);
			counter++; isFound = true;
		}	
		if (strcmp(argin[counter], "Multi_stim") == 0)
		{
			A->Multi_stim               = argin[counter+1];
			A->Multi_stim_arg			= true;
			fprintf(out, "Multi_stim   %s ", argin[counter+1]);
			if (strcmp(A->Multi_stim, "On") != 0 && strcmp(A->Multi_stim, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Multi_stim argument. Please pass only \"Off\" or \"On\"\n\n", A->Multi_stim);
				exit(1);
			}
			counter++; isFound = true;
		}
		// End Simulation setttings =========================================//|

		// Model and cell conditions ========================================\\|
		if (strcmp(argin[counter], "Model") == 0)
		{
			//printf("%s %s\n", argin[counter], argin[counter+1]);
			A->Model			= argin[counter+1];
			A->Model_arg       	= true;
			fprintf(out, "Model %s ", argin[counter+1]);
			//printf("%s\n", A->Model);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Celltype") == 0)
		{
			A->Celltype			= argin[counter+1];
			A->Celltype_arg    	= true;
			fprintf(out, "Celltype %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Agent") == 0)
		{
			A->Agent			= argin[counter+1];
			A->Agent_arg       	= true;
			fprintf(out, "Agent %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Remodelling") == 0)
		{
			A->Remodelling		= argin[counter+1];
			A->Remodelling_arg 	= true;
			fprintf(out, "Remodelling %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Mutation") == 0)
		{
			A->Mutation      	= argin[counter+1];
			A->Mutation_arg  	= true;
			fprintf(out, "Mutation %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ISO") == 0)
		{
			A->ISO				= atof(argin[counter+1]);
			A->ISO_arg       	= true;
			fprintf(out, "ISO %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ISO_model") == 0)
		{
			A->ISO_model        = argin[counter+1];
			A->ISO_model_arg    = true;
			fprintf(out, "ISO_model %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ACh") == 0)
		{
			A->ACh              = atof(argin[counter+1]);
			A->ACh_arg          = true;
			fprintf(out, "ACh %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ACh_model") == 0)
		{
			A->ACh_model        = argin[counter+1];
			A->ACh_model_arg    = true;
			fprintf(out, "ACh_model %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Remodelling_proportion") == 0)
		{
			A->Remodelling_prop              = atof(argin[counter+1]);
			A->Remodelling_prop_arg          = true;
			fprintf(out, "Remodelling_proportion %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Agent_proportion") == 0)
		{
			A->Agent_prop              = atof(argin[counter+1]);
			A->Agent_prop_arg          = true;
			fprintf(out, "Agent_proportion %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Spatial_gradient") == 0)
		{
			A->spatial_gradient              = (argin[counter+1]);
			A->spatial_gradient_arg          = true;
			fprintf(out, "Spatial_gradient %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Spatial_gradient_proportion") == 0)
		{
			A->spatial_gradient_prop              = atof(argin[counter+1]);
			A->spatial_gradient_prop_arg          = true;
			fprintf(out, "Spatial_gradient_proportion %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		// hAM single cell specific settings
		if (strcmp(argin[counter], "Environment") == 0 || strcmp(argin[counter], "environment") == 0)
		{
			A->environment          = argin[counter+1];
			A->environment_arg    	= true;
			fprintf(out, "environment %s ", argin[counter+1]);
			if (strcmp(A->environment, "isolated") != 0 && strcmp(A->environment, "intact") != 0)
			{
				printf("ERROR: \"%s\" is not a valid environment argument. Please pass only \"intact\" or \"isolated\"\n\n", A->environment);
				exit(1);
			}
			counter++; isFound = true;
		}
		// End model and cell conditions ====================================//|

		// Applied hyperpolarising current ==================================\\|
		if (strcmp(argin[counter], "Ihyp") == 0)
		{
			A->AIhyp            = atof(argin[counter+1]);
			A->Ihyp_arg			= true;
			fprintf(out, "Ihyp %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		// End Applied hyperpolarising current ==============================//|

		// Tissue model arguments ===========================================\\|
		if (strcmp(argin[counter], "Tissue_order") == 0)
		{
			A->Tissue_order        = argin[counter+1];
			A->Tissue_order_arg    = true;
			fprintf(out, "Tissue_order %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Tissue_model") == 0)
		{
			A->Tissue_model        = argin[counter+1];
			A->Tissue_model_arg    = true;
			fprintf(out, "Tissue_model %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Tissue_type") == 0)
		{
			A->Tissue_type        = argin[counter+1];
			A->Tissue_type_arg    = true;
			fprintf(out, "Tissue_type %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Orientation_type") == 0)
		{
			A->Orientation_type        = argin[counter+1];
			A->Orientation_type_arg    = true;
			fprintf(out, "Orientation_type %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "D_uniformity") == 0)
		{
			A->D_uniformity        = argin[counter+1];
			A->D_uniformity_arg    = true;
			fprintf(out, "D_uniformity %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Stimulus_location_type") == 0)
		{
			A->Stimulus_loc_type   = argin[counter+1];
			A->Stimulus_type_arg    = true;
			fprintf(out, "Stimulus_location_type %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_Stimulus_location_type") == 0)
		{
			A->S2_Stimulus_loc_type   = argin[counter+1];
			A->S2_Stimulus_type_arg   = true;
			fprintf(out, "S2_Stimulus_location_type %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Dscale") == 0)
		{
			A->Dscale        = atof(argin[counter+1]);
			A->Dscale_arg    = true;
			fprintf(out, "Dscale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "D_AR_scale") == 0)
		{
			A->D_AR_scale        = atof(argin[counter+1]);
			A->D_AR_scale_arg    = true;
			fprintf(out, "D_AR_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Multiple_models") == 0)
		{
			A->Multiple_models            = argin[counter+1];
			A->Multiple_models_arg        = true;
			fprintf(out, "Multiple_models %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Tissue_model_2") == 0)
		{
			A->Tissue_model_2            = argin[counter+1];
			A->Tissue_model_2_arg        = true;
			fprintf(out, "Tissue_model_2 %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		// Stimulus shape loc and size
		if (strcmp(argin[counter], "S1_shape") == 0)
		{
			A->S1_shape   	= argin[counter+1];
			A->S1_shape_arg   = true;
			fprintf(out, "S1_shape %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_x_loc") == 0)
		{
			A->S1_x_loc      = atoi(argin[counter+1]);
			A->S1_x_loc_arg  = true;
			fprintf(out, "S1_x_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_x_size") == 0)
		{
			A->S1_x_size      = atoi(argin[counter+1]);
			A->S1_x_size_arg  = true;
			fprintf(out, "S1_x_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_y_loc") == 0)
		{
			A->S1_y_loc      = atoi(argin[counter+1]);
			A->S1_y_loc_arg  = true;
			fprintf(out, "S1_y_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_y_size") == 0)
		{
			A->S1_y_size      = atoi(argin[counter+1]);
			A->S1_y_size_arg  = true;
			fprintf(out, "S1_y_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_z_loc") == 0)
		{
			A->S1_z_loc      = atoi(argin[counter+1]);
			A->S1_z_loc_arg  = true;
			fprintf(out, "S1_z_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S1_z_size") == 0)
		{
			A->S1_z_size      = atoi(argin[counter+1]);
			A->S1_z_size_arg  = true;
			fprintf(out, "S1_z_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		if (strcmp(argin[counter], "S2_shape") == 0)
		{
			A->S2_shape   	= argin[counter+1];
			A->S2_shape_arg   = true;
			fprintf(out, "S2_shape %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_x_loc") == 0)
		{
			A->S2_x_loc      = atoi(argin[counter+1]);
			A->S2_x_loc_arg  = true;
			fprintf(out, "S2_x_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_x_size") == 0)
		{
			A->S2_x_size      = atoi(argin[counter+1]);
			A->S2_x_size_arg  = true;
			fprintf(out, "S2_x_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}	
		if (strcmp(argin[counter], "S2_y_loc") == 0)
		{
			A->S2_y_loc      = atoi(argin[counter+1]);
			A->S2_y_loc_arg  = true;
			fprintf(out, "S2_y_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_y_size") == 0)
		{
			A->S2_y_size      = atoi(argin[counter+1]);
			A->S2_y_size_arg  = true;
			fprintf(out, "S2_y_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_z_loc") == 0)
		{
			A->S2_z_loc      = atoi(argin[counter+1]);
			A->S2_z_loc_arg  = true;
			fprintf(out, "S2_z_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "S2_z_size") == 0)
		{
			A->S2_z_size      = atoi(argin[counter+1]);
			A->S2_z_size_arg  = true;
			fprintf(out, "S2_z_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		// Diffusion parameters
		if (strcmp(argin[counter], "D1") == 0)
		{
			A->D1      = atof(argin[counter+1]);
			A->D1_arg  = true;
			fprintf(out, "D1 %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "D_AR") == 0)
		{
			A->D_AR      = atof(argin[counter+1]);
			A->D_AR_arg  = true;
			fprintf(out, "D_AR %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "dx") == 0)
		{
			A->dx      = atof(argin[counter+1]);
			A->dx_arg  = true;
			fprintf(out, "dx %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "OX") == 0)
		{
			A->OX      = atof(argin[counter+1]);
			A->OX_arg  = true;
			fprintf(out, "OX %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "OY") == 0)
		{
			A->OY      = atof(argin[counter+1]);
			A->OY_arg  = true;
			fprintf(out, "OY %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "OZ") == 0)
		{
			A->OZ      = atof(argin[counter+1]);
			A->OZ_arg  = true;
			fprintf(out, "OZ %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Global_orientation_direction") == 0)
		{
			A->Global_orientation_direction      = (argin[counter+1]);
			A->Global_orientation_direction_arg  = true;
			fprintf(out, "Global_orientation_direction %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		// remodelling/iso etc maps
		if (strcmp(argin[counter], "Dscale_mod_map") == 0)
		{
			A->Dscale_mod_map_on       = argin[counter+1];
			A->Dscale_mod_map_arg      = true;
			fprintf(out, "Dscale_mod_map   %s ", argin[counter+1]);
			if (strcmp(A->Dscale_mod_map_on, "On") != 0 && strcmp(A->Dscale_mod_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Dscale_mod_map argument. Please pass only \"Off\" or \"On\"\n\n", A->Dscale_mod_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "D_AR_scale_mod_map") == 0)
		{
			A->D_AR_scale_mod_map_on       = argin[counter+1];
			A->D_AR_scale_mod_map_arg      = true;
			fprintf(out, "D_AR_scale_mod_map   %s ", argin[counter+1]);
			if (strcmp(A->D_AR_scale_mod_map_on, "On") != 0 && strcmp(A->D_AR_scale_mod_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid D_AR_scale_mod_map argument. Please pass only \"Off\" or \"On\"\n\n", A->D_AR_scale_mod_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Remodelling_map") == 0)
		{
			A->Remodelling_map_on     	= argin[counter+1];
			A->Remodelling_map_arg		= true;
			fprintf(out, "Remodelling_map   %s ", argin[counter+1]);
			if (strcmp(A->Remodelling_map_on, "On") != 0 && strcmp(A->Remodelling_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Remodelling_map argument. Please pass only \"Off\" or \"On\"\n\n", A->Remodelling_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ISO_map") == 0)
		{
			A->ISO_map_on       = argin[counter+1];
			A->ISO_map_arg      = true;
			fprintf(out, "ISO_map   %s ", argin[counter+1]);
			if (strcmp(A->ISO_map_on, "On") != 0 && strcmp(A->ISO_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid ISO_map argument. Please pass only \"Off\" or \"On\"\n\n", A->ISO_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ACh_map") == 0)
		{
			A->ACh_map_on       = argin[counter+1];
			A->ACh_map_arg      = true;
			fprintf(out, "ACh_map   %s ", argin[counter+1]);
			if (strcmp(A->ACh_map_on, "On") != 0 && strcmp(A->ACh_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid ACh_map argument. Please pass only \"Off\" or \"On\"\n\n", A->ACh_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_map") == 0)
		{
			A->SRF_map_on       = argin[counter+1];
			A->SRF_map_arg      = true;
			fprintf(out, "SRF_map   %s ", argin[counter+1]);
			if (strcmp(A->SRF_map_on, "On") != 0 && strcmp(A->SRF_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid SRF_map argument. Please pass only \"Off\" or \"On\"\n\n", A->SRF_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Direct_modulation_map") == 0)
		{
			A->Direct_modulation_map_on       = argin[counter+1];
			A->Direct_modulation_map_arg      = true;
			fprintf(out, "Direct_modulation_map   %s ", argin[counter+1]);
			if (strcmp(A->Direct_modulation_map_on, "On") != 0 && strcmp(A->Direct_modulation_map_on, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Direct_modulation_map argument. Please pass only \"Off\" or \"On\"\n\n", A->Direct_modulation_map_on);
				exit(1);
			}
			counter++; isFound = true;
		}

		// Idealised map shape and coords
		if (strcmp(argin[counter], "map_shape") == 0)
		{
			A->map_shape     = argin[counter+1];
			A->map_shape_arg   = true;
			fprintf(out, "map_shape %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_in_type") == 0)
		{
			A->map_in_type     = argin[counter+1];
			A->map_in_type_arg   = true;
			fprintf(out, "map_in_type %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_x_loc") == 0)
		{
			A->map_x_loc      = atoi(argin[counter+1]);
			A->map_x_loc_arg  = true;
			fprintf(out, "map_x_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_x_size") == 0)
		{
			A->map_x_size      = atoi(argin[counter+1]);
			A->map_x_size_arg  = true;
			fprintf(out, "map_x_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_y_loc") == 0)
		{
			A->map_y_loc      = atoi(argin[counter+1]);
			A->map_y_loc_arg  = true;
			fprintf(out, "map_y_loc %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_y_size") == 0)
		{
			A->map_y_size      = atoi(argin[counter+1]);
			A->map_y_size_arg  = true;
			fprintf(out, "map_y_size %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "map_z_loc") == 0)
		{
			A->map_z_loc      = atoi(argin[counter+1]);
            A->map_z_loc_arg  = true;
            fprintf(out, "map_z_loc %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "map_z_size") == 0)
        {
            A->map_z_size      = atoi(argin[counter+1]);
            A->map_z_size_arg  = true;
            fprintf(out, "map_z_size %s ", argin[counter+1]);
            counter++; isFound = true;
        }

		// Direct control of input files (all except geo and orientation)
		if (strcmp(argin[counter], "stim_file") == 0)
        {
            A->stim_file      = argin[counter+1];
            A->stim_file_arg  = true;
            fprintf(out, "stim_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }	
		if (strcmp(argin[counter], "S2_stim_file") == 0)
        {
            A->S2_stim_file      = argin[counter+1];
            A->S2_stim_file_arg  = true;
            fprintf(out, "S2_stim_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "phase_file") == 0)
        {
            A->phase_file      = argin[counter+1];
            A->phase_file_arg  = true;
            fprintf(out, "phase_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "Dscale_base_map_file") == 0)
        {
            A->Dscale_base_map_file      = argin[counter+1];
            A->Dscale_base_map_file_arg  = true;
            fprintf(out, "Dscale_base_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Dscale_mod_map_file") == 0)
        {
            A->Dscale_mod_map_file      = argin[counter+1];
            A->Dscale_mod_map_file_arg  = true;
            fprintf(out, "Dscale_mod_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "D_AR_scale_base_map_file") == 0)
        {
            A->D_AR_scale_base_map_file      = argin[counter+1];
            A->D_AR_scale_base_map_file_arg  = true;
            fprintf(out, "D_AR_scale_base_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "D_AR_scale_mod_map_file") == 0)
        {
            A->D_AR_scale_mod_map_file      = argin[counter+1];
            A->D_AR_scale_mod_map_file_arg  = true;
            fprintf(out, "D_AR_scale_mod_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ISO_map_file") == 0)
        {
            A->ISO_map_file      = argin[counter+1];
            A->ISO_map_file_arg  = true;
            fprintf(out, "ISO_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "remod_map_file") == 0)
        {
            A->remod_map_file      = argin[counter+1];
            A->remod_map_file_arg  = true;
            fprintf(out, "remod_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ACh_map_file") == 0)
        {
            A->ACh_map_file      = argin[counter+1];
            A->ACh_map_file_arg  = true;
            fprintf(out, "ACh_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_map_file") == 0)
        {
            A->SRF_map_file      = argin[counter+1];
            A->SRF_map_file_arg  = true;
            fprintf(out, "SRF_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "Direct_modulation_map_file") == 0)
        {
            A->Direct_modulation_map_file      = argin[counter+1];
            A->Direct_modulation_map_file_arg  = true;
            fprintf(out, "Direct_modulation_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "Spatial_gradient_map_file") == 0)
        {
            A->spatial_gradient_map_file      = argin[counter+1];
            A->spatial_gradient_map_file_arg  = true;
            fprintf(out, "Spatial_gradient_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		// End Tissue model arguments =======================================//|

		// Spatial cell models ==============================================\\|
		if (strcmp(argin[counter], "Cell_size") == 0)
        {
            A->Cell_size        = argin[counter+1];
            A->Cell_size_arg    = true;
            fprintf(out, "Cell_size %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Sim_cell_size") == 0)
        {
            A->Sim_cell_size        = argin[counter+1];
            A->Sim_cell_size_arg    = true;
            fprintf(out, "Sim_cell_size %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Cai") == 0)
        {
            A->Cai_IC        	= atof(argin[counter+1]);
            A->Cai_IC_arg    	= true;
            fprintf(out, "Cai %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "CaSR") == 0)
        {
            A->CaSR_IC           = atof(argin[counter+1]);
            A->CaSR_IC_arg       = true;
            fprintf(out, "CaSR %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "RyR_Po") == 0)
        {
            A->RyR_Po           = atof(argin[counter+1]);
            A->RyR_Po_arg       = true;
            fprintf(out, "RyR_Po %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "LTCC_Po") == 0)
        {
            A->LTCC_Po           = atof(argin[counter+1]);
            A->LTCC_Po_arg       = true;
            fprintf(out, "LTCC_Po %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Detub") == 0)
        {
            A->Detub        = argin[counter+1];
            A->Detub_arg    = true;
            fprintf(out, "Detub %s ", argin[counter+1]);
			if (strcmp(A->Detub, "On") != 0 && strcmp(A->Detub, "Off") != 0)
            {
                printf("ERROR: \"%s\" is not a valid Detub argument. Please pass only \"Off\" or \"On\"\n\n", A->Detub);
                exit(1);
            }
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "TT_map_file") == 0)
        {
            A->TT_map_file        = argin[counter+1];
            A->TT_map_file_arg    = true;
            fprintf(out, "TT_map %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "SERCA_het") == 0)
        {
            A->SERCA_het        = argin[counter+1];
            A->SERCA_het_arg    = true;
            fprintf(out, "SERCA_het %s ", argin[counter+1]);
            if (strcmp(A->SERCA_het, "On") != 0 && strcmp(A->SERCA_het, "Off") != 0)
            {
                printf("ERROR: \"%s\" is not a valid SERCA_het argument. Please pass only \"Off\" or \"On\"\n\n", A->SERCA_het);
                exit(1);
            }
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "SERCA_map_file") == 0)
        {
            A->SERCA_map_file        = argin[counter+1];
            A->SERCA_map_file_arg    = true;
            fprintf(out, "SERCA_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "NCX_het") == 0)
        {
            A->NCX_het        = argin[counter+1];
            A->NCX_het_arg    = true;
            fprintf(out, "NCX_het %s ", argin[counter+1]);
            if (strcmp(A->NCX_het, "On") != 0 && strcmp(A->NCX_het, "Off") != 0)
            {
                printf("ERROR: \"%s\" is not a valid NCX_het argument. Please pass only \"Off\" or \"On\"\n\n", A->NCX_het);
                exit(1);
            }
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "NCX_map_file") == 0)
        {
            A->NCX_map_file        = argin[counter+1];
            A->NCX_map_file_arg    = true;
            fprintf(out, "NCX_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "RyR_het") == 0)
        {
            A->RyR_het        = argin[counter+1];
            A->RyR_het_arg    = true;
            fprintf(out, "RyR_het %s ", argin[counter+1]);
            if (strcmp(A->RyR_het, "Off") != 0 && strcmp(A->RyR_het, "random") != 0 && strcmp(A->RyR_het, "map") != 0)
            {
                printf("ERROR: \"%s\" is not a valid RyR_het argument. Please pass only \"Off\", \"map\" or \"random\"\n\n", A->RyR_het);
                exit(1);
            }
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "RyR_het_map_file") == 0)
        {
            A->RyR_het_map_file        = argin[counter+1];
            A->RyR_het_map_file_arg    = true;
            fprintf(out, "RyR_het_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "LTCC_het") == 0)
        {
            A->LTCC_het        = argin[counter+1];
            A->LTCC_het_arg    = true;
            fprintf(out, "LTCC_het %s ", argin[counter+1]);
            if (strcmp(A->LTCC_het, "Off") != 0 && strcmp(A->LTCC_het, "random") != 0 && strcmp(A->LTCC_het, "map") != 0)
            {
                printf("ERROR: \"%s\" is not a valid LTCC_het argument. Please pass only \"Off\", \"map\" or \"random\"\n\n", A->LTCC_het);
                exit(1);
            }
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "LTCC_map_file") == 0)
        {
            A->LTCC_map_file        = argin[counter+1];
            A->LTCC_map_file_arg    = true;
            fprintf(out, "LTCC_map_file %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "volds_het") == 0)
        {
            A->volds_het        = argin[counter+1];
            A->volds_het_arg    = true;
            fprintf(out, "volds_het %s ", argin[counter+1]);
            if (strcmp(A->volds_het, "On") != 0 && strcmp(A->volds_het, "Off") != 0)
            {
                printf("ERROR: \"%s\" is not a valid volds_het argument. Please pass only \"Off\" or \"On\"\n\n", A->volds_het);
                exit(1);
            }
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "tau_ss_type") == 0)
        {
            A->tau_ss_type        = argin[counter+1];
            A->tau_ss_arg    = true;
            fprintf(out, "tau_ss_type %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Delayed_CaSR_IC") == 0)
		{
			A->Delayed_CaSR_IC		= argin[counter+1];
			A->Delayed_CaSR_IC_arg	= true;
			fprintf(out, "Delayed_CaSR_IC %s ", argin[counter+1]);
            if (strcmp(A->Delayed_CaSR_IC, "On") != 0 && strcmp(A->Delayed_CaSR_IC, "Off") != 0)
			{
				printf("ERROR: \"%s\" is not a valid Delayed_CaSR_IC argument. Please pass only \"Off\" or \"On\"\n\n", A->Detub);
				exit(1);
			}
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "CaSR_IC_delay") == 0)
		{
			A->CaSR_IC_delay      = atof(argin[counter+1]);
			A->CaSR_IC_delay_arg  = true;
			fprintf(out, "CaSR_IC_delay %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		// End Spatial cell models ==========================================//|

		// Spontaneous Release Function options =============================\\i|
		if (strcmp(argin[counter], "SRF_mode") == 0)
		{
			A->SRF_mode        = argin[counter+1];
			A->SRF_mode_arg    = true;
			fprintf(out, "SRF_mode %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_model") == 0)
		{
			A->SRF_model       = argin[counter+1];
			A->SRF_model_arg   = true;
			fprintf(out, "SRF_model %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_Pset") == 0)
		{
			A->SRF_Pset        = argin[counter+1];
			A->SRF_Pset_arg    = true;
			fprintf(out, "SRF_Pset %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_het") == 0)
		{
			A->SRF_het        = argin[counter+1];
			A->SRF_het_arg    = true;
			fprintf(out, "SRF_het %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		// User control, DC model
		if (strcmp(argin[counter], "SRF_DC_PSCRE") == 0)
		{
			A->SRF_PSCRE       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_PSCRE %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_DC_CF") == 0)
		{
			A->SRF_CF_ti_sep       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_CF %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_DC_ti_sep") == 0)
		{
			A->SRF_ti_sep       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_ti_sep %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_DC_ti_W1") == 0)
		{
			A->SRF_ti_dist_width_1       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_ti_W1 %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		// todo Jakub - to prevent excessive nesting.
		if (strcmp(argin[counter], "SRF_DC_ti_W2") == 0)
		{
			A->SRF_ti_dist_width_2       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_ti_W2 %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_DC_MD") == 0)
		{
			A->SRF_MD       = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_MD %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "SRF_DC_duration_W") == 0)
		{
			A->SRF_duration_width   = atof(argin[counter+1]);
			fprintf(out, "SRF_DC_duration_W %s ", argin[counter+1]);
			A->SRF_N_DC_set++;
			counter++; isFound = true;
        }

		// User control, dynamic model
		if (strcmp(argin[counter], "SRF_Dyn_PSCR_threshold") == 0)
        {
            A->SRF_PSCR_threshold       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_PSCR_threshold %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_CaSR_max") == 0)
        {
            A->SRF_CaSR_max       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_CaSR_max %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_CaSR_Prange") == 0)
        {
            A->SRF_CaSR_Prange       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_CaSR_Prange %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_ti_sep_max") == 0)
        {
            A->SRF_ti_sep_max       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_ti_sep_max %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_ti_sep_min") == 0)
        {
            A->SRF_ti_sep_min       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_ti_sep_min %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_ti_width_max") == 0)
        {
            A->SRF_ti_width_max       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_ti_width_max %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_ti_width_min") == 0)
        {
            A->SRF_ti_width_min       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_ti_width_min %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_MD_max") == 0)
        {
            A->SRF_MD_max       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_MD_max %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_MD_min") == 0)
        {
            A->SRF_MD_min       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_MD_min %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_duration_width_max") == 0)
        {
            A->SRF_duration_width_max       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_duration_width_max %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_duration_width_min") == 0)
        {
            A->SRF_duration_width_min       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_duration_width_min %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "SRF_Dyn_H") == 0)
        {
            A->SRF_CaSR_width_H       = atof(argin[counter+1]);
            fprintf(out, "SRF_Dyn_H %s ", argin[counter+1]);
            A->SRF_N_Dyn_set++;
            counter++; isFound = true;
        }
		// End Spontaneous Release Function options =========================//|

		// Direct control current modidification ============================\\|
		if (strcmp(argin[counter], "INa_scale") == 0)
        {
            A->GNa              = atof(argin[counter+1]);
			A->GNa_arg				= true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "INa_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INaL_scale") == 0)
        {
            A->GNaL  			= atof(argin[counter+1]);
            A->GNaL_arg  		= true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "INaL_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Ito_scale") == 0)
        {
            A->Gto             = atof(argin[counter+1]);
            A->Gto_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "Ito_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaL_scale") == 0)
        {
            A->GCaL             = atof(argin[counter+1]);
            A->GCaL_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "ICaL_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_scale") == 0)
        {
            A->GKur             = atof(argin[counter+1]);
            A->GKur_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IKur_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }	
		if (strcmp(argin[counter], "IKr_scale") == 0)
        {
            A->GKr             = atof(argin[counter+1]);
            A->GKr_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IKr_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKs_scale") == 0)
        {
            A->GKs             = atof(argin[counter+1]);
            A->GKs_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IKs_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IK1_scale") == 0)
        {
            A->GK1             = atof(argin[counter+1]);
            A->GK1_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IK1_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INCX_scale") == 0)
        {
            A->GNCX             = atof(argin[counter+1]);
            A->GNCX_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "INCX_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaP_scale") == 0)
        {
            A->GCaP             = atof(argin[counter+1]);
            A->GCaP_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "ICaP_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INab_scale") == 0)
        {
            A->GNab             = atof(argin[counter+1]);
            A->GNab_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "INab_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICab_scale") == 0)
        {
            A->GCab             = atof(argin[counter+1]);
            A->GCab_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "ICab_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKb_scale") == 0)
        {
            A->GKb             = atof(argin[counter+1]);
            A->GKb_arg         = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IKb_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INaK_scale") == 0)
        {
            A->GNaK            = atof(argin[counter+1]);
            A->GNaK_arg        = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "INaK_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IClCa_scale") == 0)
        {
            A->GClCa           = atof(argin[counter+1]);
            A->GClCa_arg       = true;
			A->DC_current_mod_arg 	= true;
            fprintf(out, "IClCa_scale %s ", argin[counter+1]);
            counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IClb_scale") == 0)
		{
			A->GClb           = atof(argin[counter+1]);
			A->GClb_arg       = true;
			A->DC_current_mod_arg   = true;
			fprintf(out, "IClb_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKACh_scale") == 0)
		{
			A->GKACh           = atof(argin[counter+1]);
			A->GKACh_arg       = true;
			A->DC_current_mod_arg   = true;
			fprintf(out, "IKACh_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Jup_scale") == 0)
		{
			A->Gup            = atof(argin[counter+1]);
			A->Gup_arg        = true;
			A->DC_current_mod_arg   = true;
			fprintf(out, "Jup_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Jleak_scale") == 0)
		{
			A->Gleak            = atof(argin[counter+1]);
			A->Gleak_arg        = true;
			A->DC_current_mod_arg   = true;
			fprintf(out, "Jleak_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Jrel_scale") == 0)
		{
			A->Grel            = atof(argin[counter+1]);
			A->Grel_arg        = true;
			A->DC_current_mod_arg   = true;
			fprintf(out, "Jrel_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}

		if (strcmp(argin[counter], "INa_va_tau_scale") == 0)
		{
			A->INa_va_tau_scale			= atof(argin[counter+1]);
			A->INa_va_tau_scale_arg		= true;
			A->DC_current_mod_arg   	= true;
			fprintf(out, "INa_va_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}	
		if (strcmp(argin[counter], "INa_vi_1_tau_scale") == 0)
		{
			A->INa_vi_1_tau_scale			= atof(argin[counter+1]);
			A->INa_vi_1_tau_scale_arg		= true;
			A->DC_current_mod_arg   	= true;
			fprintf(out, "INa_vi_1_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}	
		if (strcmp(argin[counter], "INa_vi_2_tau_scale") == 0)
		{
			A->INa_vi_2_tau_scale			= atof(argin[counter+1]);
			A->INa_vi_2_tau_scale_arg		= true;
			A->DC_current_mod_arg   	= true;
			fprintf(out, "INa_vi_2_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "INaL_va_tau_scale") == 0)
		{
			A->INaL_va_tau_scale        = atof(argin[counter+1]);
			A->INaL_va_tau_scale_arg    = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "INaL_va_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "INaL_vi_tau_scale") == 0)
		{
			A->INaL_vi_tau_scale        = atof(argin[counter+1]);
			A->INaL_vi_tau_scale_arg    = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "INaL_vi_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Ito_va_tau_scale") == 0)
		{
			A->Ito_va_tau_scale        = atof(argin[counter+1]);
			A->Ito_va_tau_scale_arg    = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "Ito_va_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Ito_vi_tau_scale") == 0)
		{
			A->Ito_vi_tau_scale        = atof(argin[counter+1]);
			A->Ito_vi_tau_scale_arg    = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "Ito_vi_tau_scale %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ICaL_va_tau_scale") == 0)
		{
            A->ICaL_va_tau_scale        = atof(argin[counter+1]);
            A->ICaL_va_tau_scale_arg    = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "ICaL_va_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaL_vi_tau_scale") == 0)
        {
            A->ICaL_vi_tau_scale        = atof(argin[counter+1]);
            A->ICaL_vi_tau_scale_arg    = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "ICaL_vi_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_va_tau_scale") == 0)
        {
            A->IKur_va_tau_scale        = atof(argin[counter+1]);
            A->IKur_va_tau_scale_arg    = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKur_va_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_vi_tau_scale") == 0)
        {
            A->IKur_vi_tau_scale        = atof(argin[counter+1]);
            A->IKur_vi_tau_scale_arg    = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKur_vi_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKr_va_tau_scale") == 0)
        {
            A->IKr_va_tau_scale        = atof(argin[counter+1]);
            A->IKr_va_tau_scale_arg    = true;
            A->DC_current_mod_arg      = true;
            fprintf(out, "IKr_va_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "IKs_va_tau_scale") == 0)
        {
            A->IKs_va_tau_scale        = atof(argin[counter+1]);
            A->IKs_va_tau_scale_arg    = true;
            A->DC_current_mod_arg      = true;
            fprintf(out, "IKs_va_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKACh_va_tau_scale") == 0)
        {
            A->IKACh_va_tau_scale        = atof(argin[counter+1]);
            A->IKACh_va_tau_scale_arg    = true;
            A->DC_current_mod_arg      = true;
            fprintf(out, "IKACh_va_tau_scale %s ", argin[counter+1]);
            counter++; isFound = true;
        }

		if (strcmp(argin[counter], "INa_va_shift") == 0)
        {
            A->INa_va_shift        		= atof(argin[counter+1]);
            A->INa_va_shift_arg   		= true;
            A->DC_current_mod_arg      	= true;
            fprintf(out, "INa_va_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INa_vi_shift") == 0)
        {
            A->INa_vi_shift        		= atof(argin[counter+1]);
            A->INa_vi_shift_arg   		= true;
            A->DC_current_mod_arg      	= true;
            fprintf(out, "INa_vi_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "INaL_va_shift") == 0)
        {
            A->INaL_va_shift            = atof(argin[counter+1]);
            A->INaL_va_shift_arg        = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "INaL_va_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "INaL_vi_shift") == 0)
        {
            A->INaL_vi_shift            = atof(argin[counter+1]);
            A->INaL_vi_shift_arg        = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "INaL_vi_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }	
		if (strcmp(argin[counter], "Ito_va_ss_shift") == 0)
        {
            A->Ito_va_ss_shift          = atof(argin[counter+1]);
            A->Ito_va_ss_shift_arg      = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "Ito_va_ss_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Ito_va_tau_shift") == 0)
        {
            A->Ito_va_tau_shift         = atof(argin[counter+1]);
            A->Ito_va_tau_shift_arg     = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "Ito_va_tau_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Ito_vi_ss_shift") == 0)
        {
            A->Ito_vi_ss_shift          = atof(argin[counter+1]);
            A->Ito_vi_ss_shift_arg      = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "Ito_vi_ss_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Ito_vi_tau_shift") == 0)
        {
            A->Ito_vi_tau_shift         = atof(argin[counter+1]);
            A->Ito_vi_tau_shift_arg     = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "Ito_vi_tau_shift %s ", argin[counter+1]);
            counter++; isFound = true;
		}
		if (strcmp(argin[counter], "Ito_shift") == 0)
        {
            A->Ito_shift          		= atof(argin[counter+1]);
            A->Ito_shift_arg      		= true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "Ito_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaL_va_ss_shift") == 0)
		{
			A->ICaL_va_ss_shift          = atof(argin[counter+1]);
			A->ICaL_va_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "ICaL_va_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ICaL_va_tau_shift") == 0)
		{
			A->ICaL_va_tau_shift         = atof(argin[counter+1]);
			A->ICaL_va_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "ICaL_va_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ICaL_vi_ss_shift") == 0)
		{
			A->ICaL_vi_ss_shift          = atof(argin[counter+1]);
			A->ICaL_vi_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "ICaL_vi_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ICaL_vi_tau_shift") == 0)
		{
			A->ICaL_vi_tau_shift         = atof(argin[counter+1]);
			A->ICaL_vi_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "ICaL_vi_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "ICaL_shift") == 0)
        {
            A->ICaL_shift          		= atof(argin[counter+1]);
            A->ICaL_shift_arg      		= true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "ICaL_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_va_ss_shift") == 0)
		{
			A->IKur_va_ss_shift          = atof(argin[counter+1]);
			A->IKur_va_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKur_va_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKur_va_tau_shift") == 0)
		{
			A->IKur_va_tau_shift         = atof(argin[counter+1]);
			A->IKur_va_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKur_va_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKur_vi_ss_shift") == 0)
		{
			A->IKur_vi_ss_shift          = atof(argin[counter+1]);
			A->IKur_vi_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKur_vi_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKur_vi_tau_shift") == 0)
		{
			A->IKur_vi_tau_shift         = atof(argin[counter+1]);
			A->IKur_vi_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKur_vi_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKur_shift") == 0)
        {
            A->IKur_shift          		= atof(argin[counter+1]);
            A->IKur_shift_arg      		= true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKur_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKr_va_ss_shift") == 0)
		{
			A->IKr_va_ss_shift          = atof(argin[counter+1]);
			A->IKr_va_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKr_va_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKr_va_tau_shift") == 0)
		{
			A->IKr_va_tau_shift         = atof(argin[counter+1]);
			A->IKr_va_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKr_va_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKr_vi_ss_shift") == 0)
        {
            A->IKr_vi_ss_shift         = atof(argin[counter+1]);
            A->IKr_vi_ss_shift_arg     = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKr_vi_ss_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKs_va_ss_shift") == 0)
		{
			A->IKs_va_ss_shift          = atof(argin[counter+1]);
			A->IKs_va_ss_shift_arg      = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKs_va_ss_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IKs_va_tau_shift") == 0)
		{
			A->IKs_va_tau_shift         = atof(argin[counter+1]);
			A->IKs_va_tau_shift_arg     = true;
			A->DC_current_mod_arg       = true;
			fprintf(out, "IKs_va_tau_shift %s ", argin[counter+1]);
			counter++; isFound = true;
		}
		if (strcmp(argin[counter], "IK1_va_shift") == 0)
        {
            A->IK1_va_shift         = atof(argin[counter+1]);
            A->IK1_va_shift_arg     = true;
            A->DC_current_mod_arg       = true; 
            fprintf(out, "IK1_va_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IK1_Erev_shift") == 0)
        {
            A->IK1_Erev_shift         = atof(argin[counter+1]);
            A->IK1_Erev_shift_arg     = true;
            A->DC_current_mod_arg       = true; 
            fprintf(out, "IK1_Erev_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKACh_va_ss_shift") == 0)
        {
            A->IKACh_va_ss_shift          = atof(argin[counter+1]);
            A->IKACh_va_ss_shift_arg      = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKACh_va_ss_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }
        if (strcmp(argin[counter], "IKACh_va_tau_shift") == 0)
        {
            A->IKACh_va_tau_shift         = atof(argin[counter+1]);
            A->IKACh_va_tau_shift_arg     = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKACh_va_tau_shift %s ", argin[counter+1]);
            counter++; isFound = true;
        }

		if (strcmp(argin[counter], "Ito_va_ss_kscale") == 0)
        {
            A->Ito_va_ss_kscale         = atof(argin[counter+1]);
            A->Ito_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "Ito_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "Ito_vi_ss_kscale") == 0)
        {
            A->Ito_vi_ss_kscale         = atof(argin[counter+1]);
            A->Ito_vi_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "Ito_vi_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaL_va_ss_kscale") == 0)
        {
            A->ICaL_va_ss_kscale         = atof(argin[counter+1]);
            A->ICaL_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "ICaL_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "ICaL_vi_ss_kscale") == 0)
        {
            A->ICaL_vi_ss_kscale         = atof(argin[counter+1]);
            A->ICaL_vi_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "ICaL_vi_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_va_ss_kscale") == 0)
        {
            A->IKur_va_ss_kscale         = atof(argin[counter+1]);
            A->IKur_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "IKur_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKur_vi_ss_kscale") == 0)
        {
            A->IKur_vi_ss_kscale         = atof(argin[counter+1]);
            A->IKur_vi_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "IKur_vi_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKr_va_ss_kscale") == 0)
        {
            A->IKr_va_ss_kscale         = atof(argin[counter+1]);
            A->IKr_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "IKr_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKr_vi_ss_kscale") == 0)
        {
            A->IKr_vi_ss_kscale         = atof(argin[counter+1]);
            A->IKr_vi_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "IKr_vi_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKs_va_ss_kscale") == 0)
        {
            A->IKs_va_ss_kscale         = atof(argin[counter+1]);
            A->IKs_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;  
            fprintf(out, "IKs_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		if (strcmp(argin[counter], "IKACh_va_ss_kscale") == 0)
        {
            A->IKACh_va_ss_kscale         = atof(argin[counter+1]);
            A->IKACh_va_ss_kscale_arg     = true;
            A->DC_current_mod_arg       = true;
            fprintf(out, "IKACh_va_ss_kscale %s ", argin[counter+1]);
            counter++; isFound = true;
        }
		// End Direct control current modidification ========================//|

        if (strcmp(argin[counter], "Settings_file") == 0)
        {
            if (counter > 1) 
            {
                printf("ERROR: settings file must be passed first\n");
                exit(1);
            }
            printf("Additional settings to settings_file being read in\n");
            //fprintf(out, "Settings from settings file \"%s\" and additional from command line: ", argin[counter+1]);
            counter++; isFound = true;
        }

		// Valid argument check =============================================\\|
		if (!isFound)
		{
			printf("ERROR: \"%s\" is not a valid argument.\n", argin[counter]);
			printf("Please select from the following options:\n\n");
            printf("\tSettings_file [filename]\n");
			printf("\tReference [text]\tResults_Reference [text]\tState_Reference_read [text]\tState_Reference_write [text]\tVclamp [On/Off]\t{Read/Write}_state [On/Off/phase/single_cell/ave] (phase for tissue 2D+ only; single_cell/ave for tissue models only)\n\n");
			printf("[Simulation settings]:\n");
			printf("\tBCL [x (ms)]\tTotal_time [x (ms)]\tPaced_time [x (ms)]\tNBeats [n]\tdt [x (ms)]\n");
			printf("\tS2  [x (ms)]\tNS2 [n]\n\n");
			printf("[Model and cell conditions]:\n");
			printf("\tModel [text]\tCelltype [text]\tAgent [text]\tRemodelling [text]\tISO [x (0-1uM)]\tISO_model [text]\n");
			printf("\tACh [0-1]\tACh_model [text]\n");
			printf("\tIhyp [x (pA/pF)]\n");
			printf("\tenvironment [intact/isolated]\n\n");
			printf("\tRemodelling_proportion [0-1]\tAgent_proportion [0-1]\n");
            printf("\tSpatial_gradient [string]\tSpatial_gradient_proportion [0-1]\n");
			printf("[Direct control current modification]:\n");
			printf("\t{INa/INaL/Ito/ICaL/IKur/IKr/IKs/IK1/INCX/ICaP/INab/ICab/IKb/INaK/IClCa/IKACh/IClb}_scale [x]\n\n");
			printf("\t{Jup/Jleak/Jrel}_scale [x]\n\n");
			printf("\t{INa_va/INa_vi_1/INa_vi_2/INaL_va/INaL_vi/Ito_va/Ito_vi/ICaL_va/ICaL_vi/IKur_va/IKur_vi/IKr_va/IKs_va/IKACh_va}_tau_scale [x]\n\n");
			printf("\t{INa_va/INa_vi/INaL_va/INaL_vi/Ito_va_ss/Ito_va_tau/Ito_vi_ss/Ito_vi_tau}_shift [x (mV)]\n\n");
			printf("\t{ICaL_va_ss/ICaL_va_tau/ICaL_vi_ss/ICaL_vi_tau}_shift [x (mV)]\n\n");
			printf("\t{IKur_va_ss/IKur_va_tau/IKur_vi_ss/IKur_vi_tau}_shift [x (mV)]\n\n");
			printf("\t{IKr_va_ss/IKr_va_tau/IKr_vi_ss/IKs_va_ss/IKs_va_tau/IK1_va/IK1_Erev/IKACh_va_ss}_shift [x (mV)]\n\n");
			printf("\t{Ito_va/Ito_vi/ICaL_va/ICaL_vi/IKur_va/IKur_vi/IKr_va/IKr_vi/IKs_va/IKACh_va}_ss_kscale [x]\n\n");

			if (strcmp(Version, "Tissue_native") == 0 || strcmp(Version, "Tissue_integrated") == 0)
			{
				printf("Additional tissue model options:\n");
				printf("\tSpatial_output_interval_{vtk/data} [int ms]\n");
				printf("\tTissue_order	[1D/2D/3D/geo]\t Tissue_model [basic, ...]\t Tissue_type [homogeneous/heterogeneous]\n");
				printf("\tOrientation_type [isotropic/anisotropic]\t D_uniformity [uniform/regional/map]\n");
                printf("\tSpatial_output_interval_{vtk/data} [int ms]\n");
				printf("\tStimulus_location_type [edge/centre/cross_field/{other specific string}]\n");
				printf("\tS2_Stimulus_location_type [S1/cross_field/edge/{other specific string}]\n");
				printf("\t{S1/S2}_shape [cuboid/sphere]\n");
				printf("\t{S1/S2}_{x/y/z}_loc [n] {S1/S2}_{x/y/z}_size [n]\n");
				printf("\tMulti_stim [On/Off]\n");
                printf("\tMultiple_models [On/Off] Tissue_model_2 [model string]\n");
				printf("\tDscale [double]\tD1 [double]\tD_AR [double]\tD_AR_scale [double]\t dx [double]\n");
				printf("\t{OX/OY/OZ} [double; 0-1]\tGlobal_orientation_direction [string: X/Y/Z/{XY/XZ/YZ}_plus/{XY/XZ/YZ}_minus/XYZ_{ppp/ppm/pmp/mpp}]\n");
				printf("\t{ISO/ACh/Remodelling/Dscale_mod/D_AR_scale_mod/Direct_modulation}_map [On/Off]\n");
                printf("\tmap_{x/y/z}_shape [cuboid/sphere] map_in_type [file/coords] map_{x/y/z}_loc [n] map_{x/y/z}_size [n]\n");
				printf("\t{stim/S2_stim/phase/Dscale_base_map/D_AR_scale_base_map/Dscale_mod_map/D_AR_scale_mod_map/ISO_map/ACh_map/remod_map/Direct_modulation_map/Spatial_gradient_map}_file [string]\n\n");
			}

			if (strcmp(Version, "Single_cell_3D") == 0)
            {
                printf("Additional spatial cell model model options:\n");
				printf("\tSpatial_output_interval_{vtk/data} [int ms]\n");
				printf("\tCell_size [string]\tSim_cell_size [string]\tCai [uM]\tCaSR [uM]\n");
				printf("\tDetub [On/Off]\tTT_map_file [string]\n");
                printf("\tSERCA_het [On/Off]\tNCX_het [On/Off]\tRyR_het [Off/random/map]\tLTCC_het [Off/random/map]\tvolds_het [On/Off]\n");
                printf("\tSERCA_map_file [string]\tNCX_map_file [string]\tRyR_het_map_file [string]\tLTCC_map_file [string]\n");
				printf("\ttau_ss_type [slow/medium_slow/medium/medium_fast/fast]\n\n");
				printf("\tDelayed_CaSR_IC [Off/On]\t CaSR_IC_delay [x (ms)]\n\n");
			}

			if (strcmp(Version, "Single_cell_0D") == 0 || strcmp(Version, "Tissue_integrated") == 0)
            {
                printf("Additional 0D integrated model options (single cell or tissue):\n");
                printf("Cai [uM]\tCaSR [uM]\n");
				printf("\tSRF_mode [Off/Direct_Control/Dynamic]\tSRF_model [3D_cell/General] (Dynamic mode only)\tSRF_Pset [User_control/string]\n");
				printf("\tSRF_het [Off/General_1_small/General_1_large]\n");
				printf("\tIF \"User_control\" Pset is passed (DC mode), must also define all of: SRF_DC_{PSCRE/CF/ti_sep/ti_W1/ti_W2/MD/duration_W} [double]\n");
				printf("\tIF \"User_control\" Pset is passed (Dynamic+General mode/model), must also define all of: SRF_Dyn_{PSCR_threshold/CaSR_max/CaSR_Prange/ti_sep_max/ti_sep_min/ti_width_max/ti_width_min/MD_max/MD_min/duration_width_max/duration_width_min/H} [double]\n\n");
				printf("\tDelayed_CaSR_IC [Off/On]\t CaSR_IC_delay [x (ms)]\n\n");

                if (strcmp(Version, "Tissue_integrated") == 0) 
                {
                    printf("\tSRF_map [On/Off]\n");
                    printf("\tSRF_map_file [string]\n\n");
                }
            }

			exit(1);
		}
		// End Valid argument check =========================================//|

		counter++; 
	}
}
// End set arguments ==================================================================//|

