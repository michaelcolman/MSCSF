// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Main file for single cell simulations, ======  //
// using the Colman-lab spatial Ca2+ handling system, with=  //
// non-spatial model reduction and spontaneous release ====  //
// function implementation. ===============================  //
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
#include <omp.h>

#include "lib/Arguments.h"
#include "lib/Initialisation.h"
#include "lib/Structs.h"
#include "lib/Model.h"
#include "lib/Read_write_state.h"
#include "lib/Outputs.h"
#include "lib/CRU.h"
#include "lib/MersenneTwister.h"
#include "lib/Spontaneous_release_functions.h"
#include "lib/myofilament.hpp"

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
	printf("|Model version: Single-cell - 0D - integrated Ca handling====|\n");
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
 	// These functions populate Argin struct from command-line arguments
    // All in lib/Arguments.c
	Argument_parameters Argin;                                  // Initialise struct 
	set_argument_defaults(&Argin);                              // Sets all boolean switches to false 
	call_argument_functions(argc, argv, &Argin, "Single_cell_0D"); // Reads settings file and calls set args
	// End Argument handling ============================//|

	// Set and assign simulation settings ===============\\|
	// Assigns Sim.{Model, BCL, reference, Beats, Total_time, read/write_state, dt, ..} from Argin struct
	Simulation_parameters Sim;									// Initialise sim parameters struct || lib/Structs.h
	set_simulation_defaults(&Sim, 0.01);						// (sim struct, dt)					|| lib/Initialisation.c
	set_simulation_settings(&Sim, Argin, "integrated");			// Set from arguments				|| lib/Initialisation.c
	// End set and assign settings ======================//|

    // Determine and set OS parameters ==================\\|
    Sim.Windows  = false;
    Sim.Mac      = false;
    Sim.Linux    = false;
    Sim.OS_set   = false;

    #ifdef _WIN32
    //printf("You have Windows Operating System");
    Sim.Windows = true;
    Sim.OS_set = true;
    #endif

    #ifdef __APPLE__
    //printf("You have Mac Operating System");
    Sim.Mac = true;
    Sim.OS_set = true;
    #endif

    #ifdef __linux__
    //printf("You have Linux Operating System");
    Sim.Linux = true;
    Sim.OS_set = true;
    #endif

    printf("*****\n");
    if (Sim.Windows == true) printf("OS is WINDOWS\n");
    if (Sim.Mac     == true) printf("OS is MAC\n");
    if (Sim.Linux   == true) printf("OS is LINUX\n");
    if (Sim.OS_set  == false)
    {
        printf("ERROR: OS is not set / valid. Please see main file to add other OS\n");
        exit(1);
    }
    printf("*****\n");
    printf("\n");
    // End determine and set OS parameters ==============//|

	// Output files =====================================\\|
	// Create folders by default or references passed in
	char * directory 	= (char*)malloc(500);
	char * results_dir  = (char*)malloc(500);
    char * res_dir_full = (char*)malloc(500); // full path including directory name
	char * params_dir   = (char*)malloc(500);
	char * mkdirectory	= (char*)malloc(500);
	char * mkfile		= (char*)malloc(500);
	if (Argin.reference_arg == true) sprintf(directory, "Outputs_0Dcell_%s", Sim.reference);
	else sprintf(directory, "Outputs_0Dcell");

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkdirectory, "mkdir -p %s", directory);
    else if (Sim.Windows == true)               sprintf(mkdirectory, "mkdir %s", directory);
    system(mkdirectory);

    if (Argin.results_reference_arg == true) sprintf(results_dir, "Results_%s", Sim.results_reference);
    else sprintf(results_dir, "Results");

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkdirectory, "mkdir -p %s/%s", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkdirectory, "mkdir %s\\%s", directory, results_dir);
    system(mkdirectory);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(res_dir_full, "%s/%s", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(res_dir_full, "%s\\%s", directory, results_dir);

    if (Sim.Mac == true || Sim.Linux == true)
    {
        if (Argin.results_reference_arg == true) sprintf(params_dir, "%s/Parameters_%s", directory, Sim.results_reference);
        else sprintf(params_dir, "%s/Parameters", directory);
        sprintf(mkdirectory, "mkdir -p %s", params_dir);
    }
    else if (Sim.Windows == true)
    {
        if (Argin.results_reference_arg == true) sprintf(params_dir, "%s\\Parameters_%s", directory, Sim.results_reference);
        else sprintf(params_dir, "%s\\Parameters", directory);
        sprintf(mkdirectory, "mkdir %s", params_dir);
    }
    system(mkdirectory);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkdirectory, "mkdir -p %s/SRF_distributions", params_dir);
    else if (Sim.Windows == true)               sprintf(mkdirectory, "mkdir %s\\SRF_distributions", params_dir);
	system(mkdirectory);
	
	// Now create actual output files
	printf(">Creating output files...\n");
    
    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/Currents.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\Currents.txt", directory, results_dir);
    ofstream out_cu(mkfile);    // Contains all current related outputs
    printf("\t %s\n", mkfile);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/Properties.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\Properties.txt", directory, results_dir);
    ofstream out_ex(mkfile);     // Contains all AP properties related outputs
    printf("\t %s\n", mkfile);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/CRU.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\CRU.txt", directory, results_dir);
    ofstream out_cru(mkfile);    // Contains all current related outputs
    printf("\t %s\n", mkfile);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/SRF_properties.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\SRF_properties.txt", directory, results_dir);
    ofstream out_srf_prop(mkfile);    // Contains a list of ti and NRyR of actually occuring SRF
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
	Ca_variables					Ca;			// Ca for spatial Ca-handling model
	CRU_variables					CRU;		// sizes and averages for 3D cell model
	Dyad_variables					Dyad;		// RyR and LTCC flux variables, local dyad params
	SR_fluxes						SR;			// SERCA and Jleak fluxes variables
	Membrane_fluxes					MEM;		// NCX, ICab, ICap
	RAND							Rand;		// alt radnom number 
	Spontaneous_release_functions	SRF;		// Spontaneous release function variables and parameters
	double 							Vm;			// Global voltage

    // Myofilament and force model
    Myofilament                     myofil; // lib/myofilament.cpp
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
	if		(strcmp(Params.Model, "minimal") == 0);
	else if	(strcmp(Params.Model, "hVM_ORD_s") == 0);
	else if	(strcmp(Params.Model, "hAM_CAZ_s") == 0);
	else if	(strcmp(Params.Model, "dAM_VA") == 0);
	else if	(strcmp(Params.Model, "mCRN") == 0);
	else
	{
		printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params.Model);
		exit(1);
	}
	// End appropriate check ============================//|

	// Set parameters (defaults and model specific) =====\\|
	// Set cellsize and spontaneous release function defaults
	spatial_cell_settings(&CRU, Argin);                 // lib/CRU.cpp
	set_SRF_defaults(&SRF, Argin);						// lib/Spontaneous_release_functions.cpp

	// Default modifiers || sets all scale factors to 1 and shifts to 0 so they can be multiplicatively applied by various modifications
	set_modification_defaults_native(&Params);		// lib/Initialisation.c

	// Set default parameters (constants etc); can be overwritten by model-specific later
    set_default_parameters(&Params);                // lib/Initialisation.c
    printf(">Default parameters set\n");

	// Set model specific parameters
	Params.dt = Sim.dt;     // Set before "set_params" called, which may explicitly set dt, for checking if dt has changed
	set_parameters_native(&Params, Params.Model);   // lib/Model.c and dependants (may set dt for model-specific)

	// Set model condition params from arguments where passed; only for those which are set in set_parameters_native() to overwrite with argument value
	if (Argin.Celltype_arg  == true)    Params.Celltype     = Argin.Celltype;
	if (Argin.ISO_model_arg == true)    Params.ISO_model    = Argin.ISO_model;
	if (Argin.ACh_model_arg == true)    Params.ACh_model    = Argin.ACh_model;
	if (Argin.Ihyp_arg      == true)    Params.AIhyp        = Argin.AIhyp;

    assign_concentrations_from_arguments(&Params, Argin);

	// Update Sim.dt if Params.dt has been explicitly set in "set_parameters" (thus Sim.dt != Params.dt), and dt has NOT been passed as a command-line argument.
	if (Argin.dt_arg    	== false && Params.dt != Sim.dt) 	Sim.dt = Params.dt;

	myofil.LSODA_set(); 		// lib/myofilament.cpp
	myofil.dt_myof = Sim.dt; 	// dt for LSODA solver same as for whole model
	printf(">Model and version specific parameters set\n");

	// Now set the default and specific integrated Ca2+ handling parameters - overwrites similar parameters set in native
    set_parameters_spatial_Ca_defaults(&Params);        // lib/Initialisation.c
    set_parameters_spatial_Ca(&Params, Params.Model);   // lib/Model.c and dependants
	update_parameters_spatial_Ca_0D(&Params);			// lib/Initialisation.c

	// Overwrite initial conditions of Cai and CaSR if argument passed
	if (Argin.Cai_IC_arg	== true)	Params.Cai			= Argin.Cai_IC;
	if (Argin.CaSR_IC_arg	== true)	Params.CaSR			= Argin.CaSR_IC;
	printf(">Default parameters set - integrated Ca2+ handling\n");
	// End set parameters (defaults and model specific) =//|

	// Set modifications ================================\\|
	// lib/Initialisation.c || sets the mod variables (Gx, x_shift/tau_scale etc) from arguments 
	assign_modification_from_arguments(&Params, Argin);

	// lib/Model.c -> lib/Model_X.cpp; calls functions which set modification variables for het and modulation
	// Updates the modifier variables (scales, shifts etc) using the defined settings for any/all het and modulation
	set_heterogeneity_and_modulation_native(&Params);
	if (strcmp(Params.Ca_cellular_het, "On") == 0) update_heterogeneity_and_modulation_integrated(&Params); // lib/Model.c 

	// scale channel numbers by expression scale	
	// (Grel and GCaL refer to expression i.e. NRyR/NLTCC)
	Params.NRyR_mean 	*= Params.Grel;
	Params.NLTCC_mean 	*= Params.GCaL;

	printf(">Heterogeneity and modulation parameters set\n");
	// end set current modification =====================\\|

	// Membrane capacitance as a function of cell size ==\\|
	Params.Cm           = Params.Cm_CRU * CRU.NTOT_CRUs;
	printf(">Cm total for whole cell = %.2f pF\n", Params.Cm);
	// End Membrane capacitance / cell size =============//|

	// Local variables from global parameters
	// Remnant of 3D model; here just assigns Dyad, Mem and SR variables from Params
	set_sub_cellular_local_scale(Params, &Dyad, &MEM, &SR); // lib/CRU.cpp

	// Initialise stimulus ==============================\\|
	stimulus_setup(Params, &Variables, Sim.dt, Sim.BCL, Sim.S2_CL, Sim.Paced_time); // lib/Model.c
	printf(">Stimulus settings set\n");
	// End initialise stimulus ==========================//|

	// Spontaneous release functions ==========\\|
	SRF_setup(&SRF, Argin); // lib/Spontaneous_release_functions.cpp -> calls appropriate set SRF parameters function 
	
	// Output SRF probabiliy distributions as used in code; static or vs CaSR
	if (strcmp(SRF.Mode, "Direct_Control") == 0) 
	{
		test_and_produce_distributions(&SRF, params_dir, &Rand, 0.0);	// lib/Spontaneous_release_functions.cpp
	}
	else if (strcmp(SRF.Mode, "Dynamic") == 0)
	{
		test_and_produce_CaSR_dependency(&SRF, params_dir, &Rand);		// lib/Spontaneous_release_functions.cpp
	}
	printf("Spontaneous release function parameters set; distribtutions written to file\n");
	// End spontaneous release functions ======//|

	// Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params, argc, argv);  // lib/Outputs.c
	output_settings_0D_cell(Params, Sim, CRU, directory, SRF);			               // lib/Outputs.c

	// Setup complete, simulation running ======\\|
	time_t rawtime;
	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Setup complete:\n");
	printf("|Code is now running. Started at %s", ctime (&rawtime));
	printf("|============================================================|\n\n");
	// End Setup complete, simulation running ==//|

	// Setup dyad params to global params as just one CRU (here is where dyad het set in 3D)
	Dyad.vol_ds		= Params.vds_CRU_mean;
	Dyad.NRyR		= Params.NRyR_mean;
	Dyad.NLTCC		= Params.NLTCC_mean;
	printf("NRyR and LLTCC set\n");

	// Set initial conditions of state variables
	// Function in lib/Model.c calls specific functions in lib/Model_X.cpp
    initial_conditions_native(&State, Params, Params.Model);                        // lib/Model.c

    // Integrated calcium handling conditions (0D)
	initial_conditions_calcium_0D(&Ca, Params);										// lib/CRU.cpp
	initial_conditions_dyad_det(&Dyad);												// lib/CRU.cpp
 
	Vm = State.Vm;
	printf("Initial conditions set\n");

	// Initialise measurement variables and flags
	initialise_measurement_variables(&Variables);                                   // lib/Initialisation.c
	MEM.NCX_SRF_mult		= 1.0;
	Dyad.Krel_SRF_mult		= 1.0;

	// Output final parameters and calculate and output 
	// voltage-dependant functions, as used in the simulation
	compute_and_output_current_functions(Params, &Variables, params_dir); // lib/Model.c

	// run voltage clamp if option is set
	if(strcmp(Sim.Vclamp, "On") == 0) 
	{
		printf("Voltage clamp being performed ....\n");
		run_voltage_clamp_0Dcell(Params, &Variables, &State, &Dyad, &Ca, &CRU, &SR, &MEM, params_dir, Sim.dt); // lib/CRU.cpp
		initial_conditions_native(&State, Params, Params.Model);                        // lib/Model.c
		initial_conditions_calcium_0D(&Ca, Params);                                     // lib/CRU.cpp
		initial_conditions_dyad_det(&Dyad);
		printf("Voltage clamp finished\n");
	}

	// Read state  ============================\\|
	if (strcmp(Sim.Read_state, "On") == 0)
	{
		Read_state_single_cell_integrated_0D(&State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_read); //lib/Read_write_state.c
		assign_CRU_variables_from_state_read(&Dyad, &Ca, State); // lib/CRU.cpp
		printf("Initial conditions / state read in from file\n");
	}
	// End Read state  ========================//|
	Vm = State.Vm;

	// Time loop ================================================================================\\|
	printf("Time loop started:\nTime = %.0fms\n",sim_time);
	for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
	{
		// Impose CaSR at specified time if argument passed (allows precise setting of CaSR during simulation)
		if (Sim.CaSR_set == false && strcmp(Sim.Delayed_CaSR_IC, "On") == 0 && sim_time >= Sim.CaSR_IC_delay) 
		{ Ca.NSR = Ca.JSR = Argin.CaSR_IC; Sim.CaSR_set = true; }

		// Assign Ca state variables (seen by ionic model) from integrated whole-cell ave variables
		State.Cai		= 1e-3*Ca.CYTO;		// Ca dependent currents, Cai (in mM not uM)
		State.CanSR		= 1e-3*Ca.NSR;		// Ca dependent currents, Cansr (in mM not uM)
		State.CajSR		= 1e-3*Ca.JSR;		// Ca dependent currents, Cajsr (in mM not uM)

		// Excitation state (necessary for SRF) | lib/Model.c 
		determine_excitation_state_integrated_0D(&Variables, Vm, sim_time, &Dyad.Ca_JSR_t_ex, Ca.JSR, &Dyad.SRF_prop_active,  SRF.SRF_prop_active, &SRF.waveform_init, &SRF.srf_set, SRF.Mode);
		Dyad.ex_switch	= Variables.ex_switch;

		// Spontaneous release functions || lib/Spontaneous_release_functions.cpp
		set_and_run_SRF(&SRF, &Dyad, SRF.Mode, &Rand, Variables.ex_switch, sim_time, Ca.JSR);
		calc_SRF_mults(&SRF, &MEM, &Dyad);

		// Compute stimulus current || lib/Model.c || sets Istims to 0 or stimmag dependant on time
		compute_Istim(Params, &Variables, Sim.Paced_time, Sim.S2_time, sim_time, iteration_counter);  	// lib/Model.c

		// Zero reaction terms so they can be sequentially modified
		Ca.SS_reac = Ca.CYTO_reac = Ca.NSR_reac = Ca.JSR_reac = 0;

		// Inter-compartment transfer || lib/CRU.cpp
		comp_J_ds_ss(Params, Ca.DS, Ca.SS, Dyad.vol_ds, &Ca.SS_reac);
		comp_J_ss_cyto(Params, Ca.SS, Ca.CYTO, &Ca.SS_reac, &Ca.CYTO_reac);	
		comp_J_nsr_jsr(Params, Ca.NSR, Ca.JSR, &Ca.NSR_reac, &Ca.JSR_reac);	

		// Comp dyad || lib/CRU.cpp -> computes and solves JCaL, Jrel and Cads/jsr fluxes
		comp_dyad_0D(Params, &Dyad, Ca.DS, Ca.JSR, Ca.SS /*to which ds is coupled*/, &Ca.JSR_reac, Vm, Sim.dt, Params.Model);

		// Buffering || lib/CRU.cpp
		comp_buffering(Params, &Ca.Bcyto, &Ca.Bss, &Ca.Bjsr, Ca.CYTO, Ca.SS, Ca.JSR);

		// Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
		comp_SR_fluxes(Params, &SR, Ca.CYTO, Ca.NSR, &Ca.CYTO_reac, &Ca.NSR_reac);

		// Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
		comp_membrane_fluxes(Params, &MEM, State, Ca.CYTO, Ca.SS, &Ca.CYTO_reac, &Ca.SS_reac, Vm, MEM.NCX_SRF_mult);

		// trpn || lib/myofilament.cpp || this is general needs to be looked at
		myofil.run_step_myofilament(1e-3*Ca.CYTO, 8, 0.015);
		Ca.CYTO_reac += -myofil.Jtrpn;

		// Update concentrations
		Ca.DS   	= (Ca.SS + Params.tau_ds*(Dyad.K_rel*Ca.JSR + Dyad.J_CaL))/(1 + Params.tau_ds*Dyad.K_rel); // quasi-steady-state approx
		Ca.SS 		= Ca.SS 		+ Ca.Bss		* 	Sim.dt*(Ca.SS_reac);	
		Ca.CYTO 	= Ca.CYTO 		+ Ca.Bcyto 		*	Sim.dt*(Ca.CYTO_reac);
		Ca.NSR 		= Ca.NSR 		+ 					Sim.dt*(Ca.NSR_reac);	
		Ca.JSR 		= Ca.JSR		+ Ca.Bjsr		*	Sim.dt*(Ca.JSR_reac);

		// Whole-cell averages || including computing currents from Ca fluxes
		calc_whole_cell_values_including_currents_from_flux_0D(Params, Ca, &CRU, Dyad, SR, MEM, CRU.NTOT_CRUs);	// lib/CRU.cpp

		// Assign currents for use in AP model
		Variables.ICaL 	= CRU.I_CAL;
		Variables.INCX	= CRU.I_NCX_bulk + CRU.I_NCX_ss;
		Variables.ICaP	= CRU.I_CaP_bulk + CRU.I_CaP_ss;
		Variables.ICab	= CRU.I_Cab_bulk + CRU.I_Cab_ss;
		// End Spatial Ca handling ====================================================//|

		// Solve the AP model || lib/Model.c -> lib/Model_X.cpp
		// This sets and updates all gates, and calculates Itot
		compute_model_integrated(Params, &Variables, &State, Vm, Sim.dt);				// lib/Model.c

		// Update Voltage
		State.Vm	= State.Vm + Sim.dt*(-(Variables.Itot + Variables.Istim + Variables.Istim_S2));

		// measurement properties | lib/Model.c  | "State.Vm" is voltage at t, "Vm" is voltage at t-dt
		calculate_measurement_properties(&Variables, Vm, State.Vm, sim_time, Sim.dt, -70, State.Cai, State.CanSR);		// -70 is APD V threshold

        // Assign global voltage to state voltage
        Vm			= State.Vm;

        // Output data to files
        if (iteration_counter%(int)(1/Sim.dt) == 0) // if sim_time is an integer (i.e. per ms)
        {
            output_currents(out_cu, sim_time, Variables, State, Vm);        // lib/Outputs.cpp || V, currents, gating variables, concs etc
            output_excitation_properties(out_ex, sim_time, Variables, Vm);  // lib/Outputs.cpp || APD, excitation state, dv/dt etc
            output_CRU(out_cru, sim_time, Ca, CRU, Vm);                     // lib/Outputs.cpp || Ca concentrations, Jrel, JCaL, membrane and SR fluxes
		}

        // Print SRF ti and NRyRopeak to file for every actually induced SCRE
        if (iteration_counter%(50*((int)(1/Sim.dt))) == 0) // as 50 is less than time between successive SRF, we only need to sample at 50 ms intervals
        print_SRF_properties_to_file(&SRF, out_srf_prop, 1);

		iteration_counter ++; // number of steps in dt
		if (iteration_counter%(2000*((int)(1/Sim.dt))) == 0) printf("Time = %.0fms\n",sim_time); // output every 2000 ms
	}
	// End Time loop ============================================================================//|

	// Print final time in simulation land
	printf("Final Time = %.0fms\n\n",sim_time);

	// Write state 
	if (strcmp(Sim.Write_state, "On") == 0)
	{
		// First, assign averages from CRU model to state for outputs
		assign_state_variables_from_CRU_write(Dyad, Ca, &State);	// lib/CRU.cpp
		Write_state_single_cell_integrated_0D(State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_write); //lib/Read_write_state.c
		printf("State written to file\n");
	}
	// End Write state

	// Output final beat properties to file and screen || APD, dvdt_max etc
	char * log_reference    = (char*)malloc(500);
	if (Sim.Mac == true || Sim.Linux == true)   sprintf(log_reference, "%s/Properties_log.txt", directory);
    else if (Sim.Windows == true)               sprintf(log_reference, "%s\\Properties_log.txt", directory);
    output_properties_to_screen(log_reference, Variables, Sim);     // lib/Outputs.cpp
	free(log_reference);

	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Code has now finished. Finished at %s", ctime (&rawtime));
	printf("|============================================================|\n");
	
    // Output disclaimer to screen
	printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
	output_disclaimer_citations(Params, Sim);   // lib/Outputs.c
	output_disclaimer_citations_spatial_cell(Params, Sim);  // lib/Outputs.c
	printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

	// free memory
	free(directory);	
	free(results_dir);	
    free(res_dir_full);
	free(params_dir);	
} 
// End Main *************************************************************************************//|

