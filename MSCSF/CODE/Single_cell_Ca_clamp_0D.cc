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
// Ca2+ clamp protocol ====================================  //
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
	printf("|Ca2+ clamp mode=============================================|\n");
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
    Argument_parameters Argin;                                  // Initialise struct                            || lib/Arguments.c
    set_argument_defaults(&Argin);                              // Sets all boolean switches to false           || lib/Arguments.c
    call_argument_functions(argc, argv, &Argin, "Single_cell_0D"); // Reads settings file and calls set args || lib/Arguments.c 
	if (Argin.Cai_IC_arg    != 1)     { printf("***ERROR Cai IC must be specified for Ca clamp\n"); exit(EXIT_FAILURE); }
    if (Argin.CaSR_IC_arg  != 1)     { printf("***ERROR CaSR IC must be specified for Ca clamp\n"); exit(EXIT_FAILURE); }
	// End Argument handling ============================//|
	
	// Set and assign simulation settings ===============\\|
	Simulation_parameters Sim;									// Initialise sim parameters struct || lib/Structs.h
	set_simulation_defaults(&Sim, 0.025);						// (sim struct, dt)					|| lib/Initialisation.c
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
    char * directory    = (char*)malloc(500);
    char * results_dir  = (char*)malloc(500);
    char * res_dir_full = (char*)malloc(500); // full path including directory name
    char * params_dir   = (char*)malloc(500);
    char * mkdirectory  = (char*)malloc(500);
    char * mkfile       = (char*)malloc(500);
    if (Argin.reference_arg == true) sprintf(directory, "Outputs_Ca_clamp_0Dcell_%s", Sim.reference);
    else sprintf(directory, "Outputs_Ca_clamp_0Dcell");

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
	Cell_parameters					Params;	
	State_variables					State;
	Model_variables					Variables;
	Ca_variables					Ca;		// Ca for spatial Ca-handling model
	CRU_variables					CRU;	// sizes and averages for 3D cell model
	Dyad_variables					Dyad;	// RyR and LTCC flux variables, local dyad params
	SR_fluxes						SR;		// SERCA and Jleak fluxes variables
	Membrane_fluxes					MEM;	// NCX, ICab, ICap
	RAND							Rand;	// alt radnom number 
	Spontaneous_release_functions	SRF;	// Spontaneous release function variables and parameters
	double 							Vm;		// Global voltage

    // Myofilament and force model
    Myofilament                     myofil; // lib/myofilament.cpp

	// Ca clamp specific
	int 							activated_switch = 0;
	printf(">Variables and structs declared\n");
	// End Initialise simulation structs and variables ==//|

	// Set model condition parameters from arguments | lib/Initialisation.c
	set_model_conditions(&Params, Argin); // Set Model, ISO, remodelling etc to default or argument controlled

    // Currently only sets human atrial flag to true if a human atrial model is set
    set_model_group_variables(&Params, Argin); // lib/Initialisation.c

	// Check Model is appropriate for this code version =\\|
	if		(strcmp(Params.Model, "minimal") == 0);
	else if	(strcmp(Params.Model, "hVM_ORD_s") == 0);
	else if	(strcmp(Params.Model, "hAM_CAZ_s") == 0);
	else
	{
		printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params.Model);
		exit(1);
	}
	// End appropriate check ============================//|

	// Set parameters (defaults and model specific) =====\\|
	spatial_cell_settings(&CRU, Argin);                 // lib/CRU.cpp
	set_SRF_defaults(&SRF, Argin);						// lib/Spontaneous_release_functions.cpp
	SRF.srf_set    = -1;								// so it doesn't need excitation to initiate
	if (strcmp(SRF.Mode, "Dynamic") != 0)  { printf("***ERROR SRF mode must be dynamic for Ca clamp!\n"); exit(EXIT_FAILURE); }

	// Default modifiers
	set_modification_defaults_native(&Params);		// lib/Initialisation.c

	// Set default global parameters (may want to overwrite a modifier here, hence defaulted above)
	set_default_parameters(&Params);				// lib/Initialisation.c
	printf(">Default parameters set\n");

	// Set model specific parameters
	Params.dt = Sim.dt; 	// Set before "set_params" called, which may explicitly set dt, for checking if dt has changed
	set_parameters_native(&Params, Params.Model);	// lib/Model.c and dependants (may set dt for model-specific)
	
     // Set model condition params from arguments where passed; only for those which are set in set_parameters_native() to overwrite with argument value
    if (Argin.Celltype_arg  == true)    Params.Celltype     = Argin.Celltype;
    if (Argin.ISO_model_arg == true)    Params.ISO_model    = Argin.ISO_model;
    if (Argin.ACh_model_arg == true)    Params.ACh_model    = Argin.ACh_model;
    if (Argin.Ihyp_arg      == true)    Params.AIhyp        = Argin.AIhyp;
    
    // Update Sim.dt if Params.dt has been explicitly set in "set_parameters" (thus Sim.dt != Params.dt), and dt has NOT been passed as a command-line argument.
	if (Argin.dt_arg    	== false && Params.dt != Sim.dt) 	Sim.dt = Params.dt;

    myofil.LSODA_set();         // lib/myofilament.cpp
    myofil.dt_myof = Sim.dt;    // dt for LSODA solver same as for whole model
	printf(">Model and version specific parameters set\n");

	// Below is after set specific native so that native Ca handling params are overwritten	
	set_parameters_spatial_Ca_defaults(&Params);	// lib/Initialisation.c
	set_parameters_spatial_Ca(&Params, Params.Model); // lib/Model.c and dependants
	update_parameters_spatial_Ca_0D(&Params);		// lib/Initialisation.c

	if (Argin.Cai_IC_arg	== true)	Params.Cai			= Argin.Cai_IC;
	if (Argin.CaSR_IC_arg	== true)	Params.CaSR			= Argin.CaSR_IC;
	printf(">Default parameters set - integrated Ca2+ handling\n");
	// End set parameters (defaults and model specific) =//|

	// Set current modification =========================\\| 
	assign_modification_from_arguments(&Params, Argin);	// lib/Initialisation.c

	set_heterogeneity_and_modulation_native(&Params);	// lib/Model.c
    update_heterogeneity_and_modulation_integrated(&Params); // lib/Model.c

	Params.NRyR_mean 	*= Params.Grel;
	//Params.NLTCC_mean 	*= Params.GCaL;
	Params.NLTCC_mean	= 0;
	Params.GNa			= 0;
	printf(">Heterogeneity and modulation parameters set\n");
	// end set current modification =====================\\|

	// Membrane capacitance as a function of cell size ==\\|
	Params.Cm           = Params.Cm_CRU * CRU.NTOT_CRUs;
	printf(">Cm total for whole cell = %.2f pF\n", Params.Cm);
	// End Membrane capacitance / cell size =============//|

	// Local variables from global parameters
	set_sub_cellular_local_scale(Params, &Dyad, &MEM, &SR); // Scale all local factors by 1 compared to global (so same functions can be called)
    printf("local factors set bu params\n");

	// Spontaneous release functions ==========\\|
	SRF_setup(&SRF, Argin); // Dynamic only so no IF below
    test_and_produce_CaSR_dependency(&SRF, params_dir, &Rand);		// lib/Spontaneous_release_functions.cpp
    printf("Spontaneous release function parameters set; distribtutions written to file\n");
    // End spontaneous release functions ======//|

    // Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params, argc, argv);  // lib/Outputs.c
    output_settings_0D_cell(Params, Sim, CRU, res_dir_full, SRF);			           // lib/Outputs.c

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
	initial_conditions_native(&State, Params, Params.Model); 						// lib/Model.c
	initial_conditions_calcium_0D(&Ca, Params);										// lib/CRU.cpp
	initial_conditions_dyad_det(&Dyad);												// lib/CRU.cpp
    initialise_measurement_variables(&Variables); 
	Vm = State.Vm;
	MEM.NCX_SRF_mult		= 1.0;
	Dyad.Krel_SRF_mult		= 1.0;
	printf("Initial conditions set\n");


	// Time loop ================================================================================\\|
	printf("Time loop started:\nTime = %.0fms\n",sim_time);
	for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
	{
		// Excitation state and measurements | lib/Model.c  | "State.Vm" is voltage at t, "Vm" is voltage at t-dt
		//determine_excitation_state(&Variables, Vm, sim_time);							
		State.Cai		= 1e-3*Ca.CYTO;		// Ca dependent currents, Cai (in mM not uM)
		State.CanSR		= 1e-3*Ca.NSR;		// Ca dependent currents, Cansr (in mM not uM)
		State.CajSR		= 1e-3*Ca.JSR;		// Ca dependent currents, Cajsr (in mM not uM)

		// Spontaneous release functions
		set_and_run_SRF(&SRF, &Dyad, SRF.Mode, &Rand, Variables.ex_switch, sim_time, Ca.JSR);					// lib/Spontaneous_release_functions.cpp
		calc_SRF_mults(&SRF, &MEM, &Dyad);		// lib/Spontaneous_release_functions.cpp

		// Zero reaction terms
		Ca.SS_reac = Ca.CYTO_reac = Ca.NSR_reac = Ca.JSR_reac = 0;

		Dyad.ex_switch = 0;//currents.ex_switch;
		// Ca clamp settings =======
		if (Dyad.active == 1) activated_switch = 1; // if dyad i is active, set switch to 1 and no longer clamp
		if (activated_switch == 0) // if not in activate state, clamp Ca
		{
			if (Ca.CYTO  <= Argin.Cai_IC) Ca.CYTO     = Argin.Cai_IC; // IC arg is the clamp value
			if (Ca.SS    <= Argin.Cai_IC) Ca.SS      = Argin.Cai_IC; // no Ca.ds as it is set to ca.ss + fluxes
			if (Dyad.NRyR_O == 0) Ca.JSR = Ca.NSR  = Argin.CaSR_IC; // if no ryrs open at all, clamp Ca SR
		}
		if (Ca.JSR > Argin.CaSR_IC) Ca.JSR = Argin.CaSR_IC; // ensure CaSR does not load above clamp value (e.g. due to leak + clamp Ca)
		if (Ca.NSR > Argin.CaSR_IC) Ca.NSR = Argin.CaSR_IC;
		// End Ca clamp ==========//

		// Inter-compartment transfer || lib/CRU.cpp
		comp_J_ds_ss(Params, Ca.DS, Ca.SS, Dyad.vol_ds, &Ca.SS_reac);
		comp_J_ss_cyto(Params, Ca.SS, Ca.CYTO, &Ca.SS_reac, &Ca.CYTO_reac);	
		comp_J_nsr_jsr(Params, Ca.NSR, Ca.JSR, &Ca.NSR_reac, &Ca.JSR_reac);	

		// Comp dyad || lib/CRU.cpp
		comp_dyad_0D(Params, &Dyad, Ca.DS, Ca.JSR, Ca.SS /*to which ds is coupled*/, &Ca.JSR_reac, Vm, Sim.dt, Params.Model);

		// Buffering || lib/CRU.cpp
		comp_buffering(Params, &Ca.Bcyto, &Ca.Bss, &Ca.Bjsr, Ca.CYTO, Ca.SS, Ca.JSR);

		// Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
		comp_SR_fluxes(Params, &SR, Ca.CYTO, Ca.NSR, &Ca.CYTO_reac, &Ca.NSR_reac);

		// Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
		comp_membrane_fluxes(Params, &MEM, State, Ca.CYTO, Ca.SS, &Ca.CYTO_reac, &Ca.SS_reac, Vm, MEM.NCX_SRF_mult);

		// Update local concentrations
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

		// Solve the AP model 
		compute_model_integrated(Params, &Variables, &State, Vm, Sim.dt);				// lib/Model.c

		// Update Voltage
		State.Vm  = State.Vm + Sim.dt*(-(Variables.Itot));

		// Assign global voltage to state voltage
		Vm			= State.Vm;

		// Output data to files
		if (iteration_counter%(int)(1/Sim.dt) == 0) // if sim_time is an integer (i.e. per ms)
		{
			output_CRU(out_cru, sim_time, Ca, CRU, Vm);										// lib/Outputs.cpp
			outcount++;	
		}

		iteration_counter ++;
		if (iteration_counter%(500*((int)(1/Sim.dt))) == 0) printf("Time = %.0fms\n",sim_time); // output every 500 ms
	}
	// End Time loop ============================================================================//|

	printf("Final Time = %.0fms\n\n",sim_time);

	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Code has now finished. Finished at %s", ctime (&rawtime));
	printf("|============================================================|\n");

    // Output disclaimer to screen
    printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    output_disclaimer_citations(Params, Sim);   // lib/Outputs.c
    output_disclaimer_citations_spatial_cell(Params, Sim);  // lib/Outputs.c
    printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

	free(directory);	
    free(res_dir_full);
} 
// End Main *************************************************************************************//|

