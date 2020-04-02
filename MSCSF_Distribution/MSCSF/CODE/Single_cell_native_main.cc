// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Main file for single cell simulations, ======  //
// using native Ca2+ handling systems. ====================  //
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

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>

#include "lib/Arguments.h"
#include "lib/Initialisation.h"
#include "lib/Structs.h"
#include "lib/Model.h"
#include "lib/Read_write_state.h"
#include "lib/Outputs.h"

using namespace std;

// Main *****************************************************************************************\\|
int main(int argc, char *argv[])
{
	// Read in path for geometry and state files ========\\|
	char PATH[1000];
	FILE *path_in;
	path_in = fopen("PATH.txt", "r");
	if (path_in == NULL) // Check if file can be opened
	{
		printf("ERROR: Cannot open \"PATH.txt\"\n");
		printf("Ensure it is present in the current directory\n");
		exit(1);
	}
	fscanf(path_in, "%[^\n]", PATH); // read up to new line into PATH char
	printf("\n*PATH location is %s*\n", PATH);
	fclose(path_in);
	// End read path ====================================//|

	// Output model version to screen ===================\\|
	printf("\n");
	printf("|============================================================|\n");
	printf("|Multi-scale simulation of cardiac electrophysiology ========|\n");
	printf("|Model version: Single-cell - 0D - native Ca handling =======|\n");
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
	// These functions populate Argin struct from command-line arguments 
	// All in lib/Arguments.c
	Argument_parameters Argin;									// Initialise struct 				
	set_argument_defaults(&Argin);								// Sets all boolean switches to false 
    call_argument_functions(argc, argv, &Argin, "Single_cell_native"); // Reads settings file and calls set args
	// End Argument handling ============================//|

	// Set and assign simulation settings ===============\\|
	// Assigns Sim.{Model, BCL, reference, Beats, Total_time, read/write_state, dt, ..} from Argin struct
	Simulation_parameters Sim;									// Initialise sim parameters struct || lib/Structs.h
	set_simulation_defaults(&Sim, 0.02);						// (sim struct, dt)					|| lib/Initialisation.c
	set_simulation_settings(&Sim, Argin, "native");				// Set from arguments				|| lib/Initialisation.c
	// End set and assign settings ======================//|

	// Output files =====================================\\|
	// Create folders by default or references passed in
	char * directory 	= (char*)malloc(500);
    char * results_dir  = (char*)malloc(500);
    char * res_dir_full = (char*)malloc(500); // full path including directory name
    char * params_dir   = (char*)malloc(500);
	char * mkdirectory	= (char*)malloc(500);
	char * mkfile		= (char*)malloc(500);
	if (Argin.reference_arg == true) sprintf(directory, "Outputs_single_native_%s", Sim.reference);
	else sprintf(directory, "Outputs_single_native");
	sprintf(mkdirectory, "mkdir -p %s", directory);
	system(mkdirectory);
    if (Argin.results_reference_arg == true) sprintf(results_dir, "Results_%s", Sim.results_reference);
    else sprintf(results_dir, "Results");
	sprintf(mkdirectory, "mkdir -p %s/%s", directory, results_dir);
	system(mkdirectory);
    sprintf(res_dir_full, "%s/%s", directory, results_dir);
    if (Argin.results_reference_arg == true) sprintf(params_dir, "%s/Parameters_%s", directory, Sim.results_reference);
    else sprintf(params_dir, "%s/Parameters", directory);
	sprintf(mkdirectory, "mkdir -p %s", params_dir);
	system(mkdirectory);
	
	// Now create actual output files
	printf(">Creating output files...\n");

    sprintf(mkfile, "%s/%s/Currents.dat", directory, results_dir);
	ofstream out_cu(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

    sprintf(mkfile, "%s/%s/Properties.dat", directory, results_dir);
	ofstream out_ex(mkfile); 	 // Contains all AP properties related outputs
	printf("\t %s\n", mkfile);

	printf("\n");
	free(mkfile);
	free(mkdirectory);
	// End Output files =================================//|

	// Initialise simulation structs and variables ======\\|
	double		sim_time				= 0;
	int			iteration_counter		= 0; 	// Number of iterations over dt
	int			outcount				= 0; 	// Output number reference

	// All below defined in lib/Structs.h
	Cell_parameters					Params;		// Parameters/constants	
	State_variables					State;		// Time-dependent state variables
	Model_variables					Variables;	// Calculated variables
	double 							Vm;			// Global copy of voltage
	printf(">Variables and structs declared\n");
	// End Initialise simulation structs and variables ==//|

	// Set model condition parameters from arguments | lib/Initialisation.c
	// Sets Params.{Model, modulation(ISO,remodelling, mutation etc), model group} to defaults
	// then overwrites with Argin values, if arguments have been passed 
	set_model_conditions(&Params, Argin);

	// Currently only sets human atrial flag to true if a human atrial model is set 
    set_model_group_variables(&Params, Argin); // lib/Initialisation.c 

	// Check Model is appropriate for this code version =\\|
	// Add your new model clause here so it doesn't return an error
	if      (strcmp(Params.Model, "minimal") == 0);			// hybrid minimal model; Colman 2019
	else if (strcmp(Params.Model, "hAM_GB") == 0);			// human atrial cell; Grandi et al 2011 Circ Res 2011; 109(9):1055-66
	else if (strcmp(Params.Model, "hAM_CRN") == 0);			// human atrial cell; Courtemanche et al 1998 Am J Physiol 1998; 275(1 Pt 2):H301-21.
	else if (strcmp(Params.Model, "hAM_NG") == 0);			// human atrial cell; Nygren et al 1998 Circ Res 1998;82(1):63-81
	else if (strcmp(Params.Model, "hAM_MT") == 0);			// human atrial cell; Maleckar et al 2009  Am. J. Physiol Heart Circ Physiol 2009;297(4):H1398-410
	else if (strcmp(Params.Model, "hAM_WL_CRN") == 0);		// human atrial cell; Colman et al 2018 Front Physiol 9:1211; CRN Ca2+ handling
	else if (strcmp(Params.Model, "hAM_WL_GB") == 0);		// human atrial cell; Colman et al 2018 Front Physiol 9:1211; GB Ca2+ handling
	else if (strcmp(Params.Model, "hAM_GB_mWL") == 0);		// human atrial cell; Colman et al 2018 Front Physiol 9:1211; GB modified
	else if (strcmp(Params.Model, "hAM_CRN_mWL") == 0);		// human atrial cell; Colman et al 2018 Front Physiol 9:1211; CRN modified
	else if (strcmp(Params.Model, "hAM_NG_mWL") == 0);		// human atrial cell; Colman et al 2018 Front Physiol 9:1211; NG modified
	else if (strcmp(Params.Model, "dAM_VA") == 0);			// dog atrial cell; Varela et al 2016 PLOS Computational Biology 12(12): e1005245
	else if (strcmp(Params.Model, "hAM_GB_modded") == 0);	// human atrial cell; Grandi et al 2011 Circ Res 2011; 109(9):1055-66 MODDED
	else
	{
		printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params.Model);
		exit(1);
	}
	// End appropriate check ============================//|

	// Set parameters (defaults and model specific) =====\\|
	// Set default parameters (constants etc); can be overwritten by model-specific later
	set_default_parameters(&Params);				// lib/Initialisation.c
	printf(">Default parameters set\n");

	// Set model specific full parameters
	Params.dt = Sim.dt; 	// Set before "set_params" called, which may explicitly set dt, for checking if dt has changed
	set_parameters_native(&Params, Params.Model);	// lib/Model.c and dependants (may set dt for model-specific)

	// Set model condition params from arguments where passed; only for those which are set in set_parameters_native() to overwrite with argument value
	if (Argin.Celltype_arg 	== true)	Params.Celltype 	= Argin.Celltype; 
	if (Argin.ISO_model_arg == true)	Params.ISO_model 	= Argin.ISO_model; 
	if (Argin.ACh_model_arg == true)    Params.ACh_model    = Argin.ACh_model; 
	if (Argin.Ihyp_arg      == true)   	Params.AIhyp        = Argin.AIhyp;

	// Update Sim.dt if Params.dt has been explicitly set in "set_parameters" (thus Sim.dt != Params.dt), and dt has NOT been passed as a command-line argument.
	if (Argin.dt_arg    	== false && Params.dt != Sim.dt) 	Sim.dt = Params.dt;
	printf(">Model and version specific parameters set\n");
	// End set parameters (defaults and model specific) =//|

	// Set current modification =========================\\|
	// Default modifiers || sets all scale factors to 1 and shifts to 0 so they can be multiplicatively applied by various modifications
	set_modification_defaults_native(&Params);		// lib/Initialisation.c

	// lib/Initialisation.c || sets the mod variables (Gx, x_shift/tau_scale etc) from arguments 
	assign_modification_from_arguments(&Params, Argin);
	
	// lib/Model.c  -> lib/Model_X.cp; calls functions which set modification variables for het and modulation
	// Updates the modifier variables (scales, shifts etc) using the defined settings for any/all het and modulation
	set_heterogeneity_and_modulation_native(&Params);

	// set expression scale by open rate scale, as both are equivilent in native models
	Params.GCaL *= Params.GLTCC_kva1_va2;
	Params.Grel	*= Params.GRyR_kCO;
	printf(">Heterogeneity and modulation parameters set\n");
	// end set current modification =====================\\| 

	// Initialise stimulus ==============================\\|
	stimulus_setup(Params, &Variables, Sim.dt, Sim.BCL, Sim.S2_CL, Sim.Paced_time); // lib/Model.c
	printf(">Stimulus settings set\n");
	// End initialise stimulus ==========================//|

	// Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params, argc, argv);  // lib/Outputs.c

	// Setup complete, simulation running ======\\|
	time_t rawtime;
	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Setup complete:\n");
	printf("|Code is now running. Started at %s", ctime (&rawtime));
	printf("|============================================================|\n\n");
	// End Setup complete, simulation running ==//|

	// Initial conditions and model function outs ===============================================\\|
	// Set initial conditions of state variables
	// Function in lib/Model.c calls specific functions in lib/Model_X.cpp
	initial_conditions_native(&State, Params, Params.Model);
	Vm = State.Vm;					// Global Vm copy to state variable
	printf("Initial conditions set\n");
    
	// Initialise measurement variables and flags
	initialise_measurement_variables(&Variables);  // lib/Initialisation.c
	
	// Output final parameters and calculate and output 
	// voltage-dependant functions, as used in the simulation
	compute_and_output_current_functions(Params, &Variables, params_dir); // lib/Model.c
	
	// run voltage clamp if option is set
	if(strcmp(Sim.Vclamp, "On") == 0)
	{
		run_voltage_clamp(Params, &Variables, &State, params_dir, Sim.dt); 	// lib/Model.c
		initial_conditions_native(&State, Params, Params.Model);    		// lib/Model.c || Re-set initial conditions
		printf("Voltage clamp completed\n");
	}

	// Read state from file  || this must be after Vclamp as that resets ICs for each step
	if (strcmp(Sim.Read_state, "On") == 0)
	{
		Read_state_single_cell_native(&State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_read); //lib/Read_write_state.c
		printf("Initial conditions / state read in from file\n");
		Vm = State.Vm; // set global voltage to state again
	}
	// End Initial conditions and model function outs ===========================================//|

	// Time loop ================================================================================\\|
	printf("Time loop started:\nTime = %.0fms\n",sim_time);
	for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
	{
		// Compute stimulus current || lib/Model.c || sets Istims to 0 or stimmag dependant on time
		compute_Istim(Params, &Variables, Sim.Paced_time, Sim.S2_time, sim_time, iteration_counter);

		// Solve the model || lib/Model.c -> lib/Model_X.cpp
		// This sets and updates all gates, and calculates Itot
		compute_model_native(Params, &Variables, &State, Vm, Sim.dt);

		// Update Voltage
		State.Vm	= State.Vm + Sim.dt*(-(Variables.Itot + Variables.Istim + Variables.Istim_S2));

		// Excitation state and measurements | lib/Model.c  | "State.Vm" is voltage at t, "Vm" is voltage at t-dt
		determine_excitation_state(&Variables, Vm, sim_time);							
		calculate_measurement_properties(&Variables, Vm, State.Vm, sim_time, Sim.dt, -70, State.Cai, State.CanSR); // -70 is APD V threshold	

		// Assign global voltage to state voltage (now both = V at t)
		Vm			= State.Vm;

		// Output data to files
		if (iteration_counter%(int)(1/Sim.dt) == 0) // if sim_time is an integer (i.e. per ms)
		{
			output_currents(out_cu, sim_time, Variables, State, Vm);		// lib/Outputs.cpp || V, currents, gating variables, concs etc
			output_excitation_properties(out_ex, sim_time, Variables, Vm);	// lib/Outputs.cpp || APD, excitation state, dv/dt etc
		}

		iteration_counter ++; // number of steps in dt
		if (iteration_counter%(2000*((int)(1/Sim.dt))) == 0) printf("Time = %.0fms\n",sim_time); // output every 2000 ms
	}
	// End Time loop ============================================================================//|

	// Print final time in simulation land
	printf("Final Time = %.0fms\n\n",sim_time);

	// Write state 
	if (strcmp(Sim.Write_state, "On") == 0)
	{
		Write_state_single_cell_native(State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_write); //lib/Read_write_state.c
		printf("State written to file\n");
	}
	// End Write state

	// Output final 2 beat properties to file and screen || APD, dvdt_max etc
	char * log_reference 	= (char*)malloc(500);
	sprintf(log_reference, "%s/Properties_log.dat", directory);
	output_properties_to_screen(log_reference, Variables, Sim); 	// lib/Outputs.cpp
	free(log_reference);

	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Code has now finished. Finished at %s", ctime (&rawtime));
	printf("|============================================================|\n");
    
    // Output disclaimer
    printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    output_disclaimer_citations(Params, Sim);   // lib/Outputs.c
    printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

	// free memory
	free(directory);
	free(results_dir);
    free(res_dir_full);
	free(params_dir);
} 
// End Main *************************************************************************************//|

