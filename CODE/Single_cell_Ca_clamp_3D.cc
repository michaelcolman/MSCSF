// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Main file for single cell simulations, ======  //
// using the Colman-lab spatial Ca2+ handling system.======  //
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
#include "lib/Spatial_coupling.h"
#include "lib/CRU.h"
#include "lib/MersenneTwister.h"
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
	printf("|Model version: Single-cell - 3D - integrated Ca handling====|\n");
	printf("|Ca2+ clamp mode=============================================|\n");
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
    Argument_parameters Argin;                                  // Initialise struct                            || lib/Arguments.c
    set_argument_defaults(&Argin);                              // Sets all boolean switches to false           || lib/Arguments.c
    call_argument_functions(argc, argv, &Argin, "Single_cell_3D"); // Reads settings file and calls set args || lib/Arguments.c
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
    char * sr_dir       = (char*)malloc(500);
    char * mkdirectory  = (char*)malloc(500);
    char * mkfile       = (char*)malloc(500);
    if (Argin.reference_arg == true) sprintf(directory, "Outputs_Ca_clamp_3Dcell_%s", Sim.reference);
    else sprintf(directory, "Outputs_Ca_clamp_3Dcell");

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

    sprintf(sr_dir, "Spatial_%s", results_dir);
    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkdirectory, "mkdir -p %s/%s", directory, sr_dir);
    else if (Sim.Windows == true)               sprintf(mkdirectory, "mkdir %s\\%s", directory, sr_dir);
    system(mkdirectory);

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

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/Ca_linescan_Z.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\Ca_linescan_Z.txt", directory, results_dir);
    ofstream out_lsZ(mkfile);     // Contains linescan of Ca, longitudinal
    printf("\t %s\n", mkfile);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/Ca_linescan_Y.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\Ca_linescan_Y.txt", directory, results_dir);
    ofstream out_lsY(mkfile);     // Contains linescan of Ca, transverse
    printf("\t %s\n", mkfile);

    if (Sim.Mac == true || Sim.Linux == true)   sprintf(mkfile, "%s/%s/Ca_linescan_X.txt", directory, results_dir);
    else if (Sim.Windows == true)               sprintf(mkfile, "%s\\%s\\Ca_linescan_X.txt", directory, results_dir);
    ofstream out_lsX(mkfile);     // Contains linescan of Ca, transverse
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
    SC_variables                    SC;     // Spatial coupling
    CRU_variables					CRU;	// sizes and averages for 3D cell model
    Dyad_variables					*Dyad;	// RyR and LTCC flux variables, local dyad params
    SR_fluxes						*SR;	// SERCA and Jleak fluxes variables
    Membrane_fluxes					*MEM;	// NCX, ICab, ICap
    RAND							*Rand;	// radnom number 
    double 							Vm;		// Global voltage

    // Ca clamp only
    int                     *activated_switch;

    // Myofilament and force model
    Myofilament             		*myofil; // lib/myofilament.cpp
    double                  		Force;
    printf(">Variables and structs declared\n");
    // End Initialise simulation structs and variables ==//|

    // Set model condition parameters from arguments | lib/Initialisation.c
    set_model_conditions(&Params, Argin); // Set Model, ISO, remodelling etc to default or argument controlled

    // Check Model is appropriate for this code version =\\|
    if      (strcmp(Params.Model, "minimal") == 0);
    else if	(strcmp(Params.Model, "hVM_ORD_s") == 0);
    else if	(strcmp(Params.Model, "hAM_CAZ_s") == 0);
    else
    {
        printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params.Model);
        exit(1);
    }
    // End appropriate check ============================//|

    // Spatial cell setup - settings and allocation =====\\|
    // Set cellsize etc
    spatial_cell_settings(&CRU, Argin);					// lib/CRU.cpp	
    SC_set_array_sizes(&SC, CRU.NX, CRU.NY, CRU.NZ);	// lib/Spatial_coupling.cpp
    SC_array_allocation_N3(&SC, SC.NX, SC.NY, SC.NZ);   // Allocates arrays of size NX*NY*NZ
    SC.N	= SC.NX * SC.NY * SC.NZ;				
    for (int n = 0; n < SC.N; n++) SC.geo[n] = 1;		// Idealised cell; no empty space

    // Array allocation
    Ca_array_allocation(SC.N, &Ca);				// lib/CRU.cpp
    SC_array_allocation_Ncell(&SC, SC.N);       // lib/Spatial_coupling.cpp
    Dyad	= new Dyad_variables[SC.N];			
    SR		= new SR_fluxes[SC.N];			
    MEM		= new Membrane_fluxes[SC.N];
    activated_switch = new int [SC.N];

    // myofilament and force
    myofil          = new Myofilament[SC.N];
    for (int n = 0; n < SC.N; n++) myofil[n].LSODA_set(); // lib/myofilament.cpp

    printf(">Spatial arrays allocated | NCRU = %d NX = %d NY = %d NZ = %d\n", SC.N, SC.NX, SC.NY, SC.NZ);

    // Cell index and neighbours
    SC_set_index_and_geo_linear(&SC);               // lib/Spatial_coupling.cpp
    SC_set_neighbours(&SC);                         // lib/Spatial_coupling.cpp
    printf(">Linear index and neighbours set\n");
    // End spatial cell setup ===========================//|

    // Set parameters (defaults and model specific) =====\\|
    // Default modifiers
    set_modification_defaults_native(&Params);		// lib/Initialisation.c

    // Set default global parameters (may want to overwrite a modifier here, hence defaulted above)
    set_default_parameters(&Params);				// lib/Initialisation.c
    printf(">Default parameters set - native\n");

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
    for (int n =0; n < SC.N; n++) myofil[n].dt_myof = Sim.dt; // dt for LSODA solver same as for whole model
    printf(">Model and version specific parameters set\n");

    set_parameters_spatial_Ca_defaults(&Params);	// lib/Initialisation.c
    set_parameters_spatial_Ca(&Params, Params.Model);	// lib/Model.c and dependants

    if (Argin.Cai_IC_arg	== true)	Params.Cai			= Argin.Cai_IC;
    if (Argin.CaSR_IC_arg	== true)	Params.CaSR			= Argin.CaSR_IC;

    printf(">Default parameters set - integrated Ca2+ handling\n");
    // End set parameters (defaults and model specific) =//|

    // Set modifications ================================\\| 
    assign_modification_from_arguments(&Params, Argin);

    // lib/Model.c -> lib/Model_X.cpp; calls functions which set modification variables for het and modulation
    // Updates the modifier variables (scales, shifts etc) using the defined settings for any/all het and modulation
    set_heterogeneity_and_modulation_native(&Params);
    update_heterogeneity_and_modulation_integrated(&Params); // lib/Model.c -> lib/Model_X.cpp

    // scale channel numbers by expression scale
    // (Grel and GCaL refer to expression i.e. NRyR/NLTCC)
    Params.NRyR_mean        *= Params.Grel;
    Params.NLTCC_mean       *= Params.GCaL;
    printf(">Heterogeneity and modulation parameters set\n");

    // Assign tau ss, which may be defined by model, modification or argumnet
    if (Argin.tau_ss_arg == true) Params.tau_ss_type = Argin.tau_ss_type; //otherwise default or set in model/modification function
    set_tau_ss(&Params);
    printf(">Subspace coupling time constants set\n");
    // end set current modification =====================//|

    // Membrane capacitance as a function of cell size ==\\|
    Params.Cm			= Params.Cm_CRU * CRU.NTOT_CRUs;
    printf(">Cm total for whole cell (pre detubulation) = %.2f pF\n", Params.Cm);
    // End Membrane capacitance / cell size =============//|

    // Set local scale params from global params ========\\|
    // Allocate arrays, apply global params to local params
    // set local scale factors from global
    CRU_map_array_allocation(SC.N, &CRU);

    // set local scale factors from global
    for (int n = 0; n < SC.N; n++) set_sub_cellular_local_scale(Params, &Dyad[n], &MEM[n], &SR[n]);  // lib/CRU.cpp

    // scale by local het map || lib/CRU.cpp
    initialise_sub_cellular_het_maps(SC.N, &CRU); // set all maps to 1
    read_sub_cellular_het_maps(SC, &CRU, PATH, directory);  // read and assign maps, where On
    set_sub_cellular_local_het_scale(Params, SC.N, Dyad, MEM, SR, CRU); // scale local flux rates by map
    for (int n = 0; n < SC.N; n++) Dyad[n].NRyR     = Params.NRyR_mean  * CRU.RyR_het_map[n];   // no need for IF, as map[n] = 1 if RyR het is not "map"
    for (int n = 0; n < SC.N; n++) Dyad[n].NLTCC    = Params.NLTCC_mean * CRU.LTCC_map[n];      // no need for IF, as map[n] = 1 if LTCC het is not "map"
    printf(">Local sub-cellular scaling parameters set\n");
    // End Set local scale params from global params ====//|

    // Allocate array for random numbers
    printf("Allocating rand array; this may take a while (but significantly improves parallelisation performance)\n");
    Rand    = new RAND [SC.N];

    // Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params, argc, argv);  // lib/Outputs.c
    output_settings_3D_cell(Params, Sim, CRU, res_dir_full);                           // lib/Outputs.c

    // Setup complete, simulation running ======\\|
    time_t rawtime;
    time (&rawtime);
    printf("|============================================================|\n");
    printf("|Setup complete:\n");
    printf("|Code is now running. Started at %s", ctime (&rawtime));
    printf("|============================================================|\n\n");
    // End Setup complete, simulation running ==//|

    // Setup dyad heterogeneity and allocate local dyad arrays
    for (int n = 0; n < SC.N; n++)
    {
        // This must be first, as sets NRyR and NLTCC per dyad (if het = random)
        // sets vol_ds based on ave and random  (if het = On)
        // rand 0.5 so no random het
        set_dyad_heterogeneity(Params, &Dyad[n], 0.5, n, params_dir, &CRU.dyad_het_map[n], CRU);    // lib/CRU.cpp
        Dyad_array_allocation(&Dyad[n]);
    }
    printf("Local NRyR and LLTCC set and dyad arrays allocated\n");

    // Set initial conditions of state variables
    // Function in lib/Model.c calls specific functions in lib/Model_X.cpp
    initial_conditions_native(&State, Params, Params.Model);                        // lib/Model.c

    // Integrated calcium handling conditions
    initial_conditions_calcium(&Ca, Params, SC.N);                                  // lib/CRU.cpp
    for (int n = 0; n < SC.N; n++) initial_conditions_dyad_stochastic(&Dyad[n]);    // lib/CRU.cpp

    Vm = State.Vm;
    printf("Initial conditions set\n");

    // Time loop ================================================================================\\|
    printf("Time loop started:\nTime = %.0fms\n",sim_time);
    for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
    {
        State.Cai		= 1e-3*Ca.CYTO;		// Ca dependent currents, Cai (in mM not uM)
        State.CanSR		= 1e-3*Ca.NSR;		// Ca dependent currents, Cansr (in mM not uM)
        State.CajSR		= 1e-3*Ca.JSR;		// Ca dependent currents, Cajsr (in mM not uM)

        // Spatial Ca handling ========================================================\\|
        // First, populate random number array; can be done in ====\\|
        // series or parallel independent of parallelisation of whole model
#pragma omp parallel for default(none) shared(SC, Dyad, Rand) //private(mtrand1)
        for (int n = 0; n < SC.N; n++)
        {
            for (int j = 0; j < Dyad[n].NRyR; j++)  Dyad[n].rand_RyR[j]   	= Rand[n].mtrand1(); // allows faster parallelisation
            for (int j = 0; j < Dyad[n].NLTCC; j++) Dyad[n].rand_LTCC[j]   	= Rand[n].mtrand1();
        }
        // End populate random number array =======================//|

        // Spatial loop - 1 =======================================\\|
#pragma omp parallel for default(none) shared(SC, Vm, Params, Variables, State, Sim, Ca, sim_time, Dyad, SR, MEM, myofil, activated_switch, Argin)
        for (int n = 0; n < SC.N; n++)
        {
            // Zero reaction terms
            Ca.ss_reac[n] = Ca.cyto_reac[n] = Ca.nsr_reac[n] = Ca.jsr_reac[n] = 0;

            if (Dyad[n].active == 1) activated_switch[n] = 1; // if dyad i is active, set switch to 1 and no longer clamp
            if (activated_switch[n] == 0) // if not in activate state, clamp Ca
            {
                if (Ca.cyto[n]  <= Argin.Cai_IC) Ca.cyto[n]     = Argin.Cai_IC; // IC arg is the clamp value
                if (Ca.ss[n]    <= Argin.Cai_IC) Ca.ss[n]      = Argin.Cai_IC; // no Ca.ds as it is set to ca.ss + fluxes
                if (Dyad[n].NRyR_O == 0) Ca.JSR = Ca.nsr[n]  = Argin.CaSR_IC; // if no ryrs open at all, clamp Ca SR
            }
            if (Ca.jsr[n] > Argin.CaSR_IC) Ca.jsr[n] = Argin.CaSR_IC; // ensure CaSR does not load above clamp value (e.g. due to leak + clamp Ca)
            if (Ca.nsr[n] > Argin.CaSR_IC) Ca.nsr[n] = Argin.CaSR_IC;

            // Inter-compartment transfer || lib/CRU.cpp
            comp_J_ds_ss(Params, Ca.ds[n], Ca.ss[n], Dyad[n].vol_ds, &Ca.ss_reac[n]);
            comp_J_ss_cyto(Params, Ca.ss[n], Ca.cyto[n], &Ca.ss_reac[n], &Ca.cyto_reac[n]);	
            comp_J_nsr_jsr(Params, Ca.nsr[n], Ca.jsr[n], &Ca.nsr_reac[n], &Ca.jsr_reac[n]);	

            // Comp dyad || lib/CRU.cpp
            comp_dyad_3D(Params, &Dyad[n], Ca.ds[n], Ca.jsr[n], Ca.ss[n] /*to which ds is coupled*/, &Ca.jsr_reac[n], Vm, Sim.dt, Params.Model);

            // Buffering || lib/CRU.cpp
            comp_buffering(Params, &Ca.bcyto[n], &Ca.bss[n], &Ca.bjsr[n], Ca.cyto[n], Ca.ss[n], Ca.jsr[n]);

            // Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
            comp_SR_fluxes(Params, &SR[n], Ca.cyto[n], Ca.nsr[n], &Ca.cyto_reac[n], &Ca.nsr_reac[n]);

            // Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
            comp_membrane_fluxes(Params, &MEM[n], State, Ca.cyto[n], Ca.ss[n], &Ca.cyto_reac[n], &Ca.ss_reac[n], Vm, 1.0);  // 1.0 is SRF mult as not relevant here

            // trpn and force || lib/myofilament.cpp
            myofil[n].run_step_myofilament(1e-3*Ca.cyto[n], 8, 0.015);
            Ca.cyto_reac[n] += -myofil[n].Jtrpn;

            // Spatial coupling ===============\\|
            // lib/Spatial_coupling.cpp
            calc_diff_FDM_tau(&SC, Ca.ss, n, Params.tau_ss_trans, Params.tau_ss_long);
            Ca.ss_reac[n] += SC.diff[n];

            calc_diff_FDM_tau(&SC, Ca.cyto, n, Params.tau_cyto_trans, Params.tau_cyto_long);
            Ca.cyto_reac[n] += SC.diff[n];

            calc_diff_FDM_tau(&SC, Ca.nsr, n, Params.tau_nsr_trans, Params.tau_nsr_long);
            Ca.nsr_reac[n] += SC.diff[n];
            // End spatial coupling ===========//|
        }
        // End spatial loop - 1 ===================================//|

#pragma omp parallel for default(none) shared(SC, Ca, Params, Dyad, Sim)
        // Spatial loop - 2 =======================================\\|
        for (int n = 0; n < SC.N; n++)
        {
            // Update local concentrations
            Ca.ds[n]   	= (Ca.ss[n] + Params.tau_ds*(Dyad[n].K_rel*Ca.jsr[n] + Dyad[n].J_CaL))/(1 + Params.tau_ds*Dyad[n].K_rel); // quasi-steady-state approx
            Ca.ss[n] 	= Ca.ss[n] 		+ Ca.bss[n]		* 	Sim.dt*(Ca.ss_reac[n]);	
            Ca.cyto[n] 	= Ca.cyto[n] 	+ Ca.bcyto[n] 	*	Sim.dt*(Ca.cyto_reac[n]);
            Ca.nsr[n] 	= Ca.nsr[n] 	+ 					Sim.dt*(Ca.nsr_reac[n]);	
            Ca.jsr[n] 	= Ca.jsr[n]		+ Ca.bjsr[n]	*	Sim.dt*(Ca.jsr_reac[n]);
        }
        // End spatial loop - 2 ===================================//|

        // Whole-cell averages || including computing currents from Ca fluxes
        calc_whole_cell_values_including_currents_from_flux(SC.N, Params, &Ca, &CRU, Dyad, SR, MEM, CRU.NTOT_CRUs);	// lib/CRU.cpp

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
        if (iteration_counter % Variables.dtinv == 0) // if sim_time is an integer (i.e. per ms)
        {
            output_CRU(out_cru, sim_time, Ca, CRU, Vm);										// lib/Outputs.cpp

            // Linescan
            linescan_out_X(out_lsX, SC, Ca.cyto, int(float(SC.NY/2)), int(float(SC.NZ/2))); 		// lib/Outputs.cpp
            linescan_out_Y(out_lsY, SC, Ca.cyto, int(float(SC.NX/2)), int(float(SC.NZ/2))); 		// lib/Outputs.cpp
            linescan_out_Z(out_lsZ, SC, Ca.cyto, int(float(SC.NX/2)), int(float(SC.NY/2))); 		// lib/Outputs.cpp

            if (Sim.Spatial_output_interval_vtk  > 0) // such that setting to zero means no spatial outputs
            {
                if (outcount %Sim.Spatial_output_interval_vtk == 0)
                {
                    vtk_3D_output("Ca", directory, "Spatial_Results", Ca.cyto, SC, outcount);   // every x ms, output vtk file
                    vtk_3D_output("CaSR", directory, "Spatial_Results", Ca.nsr, SC, outcount);  // every x ms, output vtk file
                }
            }
            if (Sim.Spatial_output_interval_data  > 0) // such that setting to zero means no spatial outputs
            {
                if (outcount %Sim.Spatial_output_interval_data == 0)
                {
                    array_1D_output("Ca", directory, "Spatial_Results", Ca.cyto, SC, outcount);   // every x ms, output vtk file
                    array_1D_output("CaSR", directory, "Spatial_Results", Ca.nsr, SC, outcount);  // every x ms, output vtk file
                }
            }
            outcount++;	
        }

        iteration_counter ++;
        if (iteration_counter%(100*(Variables.dtinv)) == 0) printf("Time = %.0fms\n",sim_time); // output every 500 ms
    }
    // End Time loop ============================================================================//|

    printf("Final Time = %.0fms\n\n",sim_time);

    time (&rawtime);
    printf("|============================================================|\n");
    printf("|Code has now finished. Finished at %s", ctime (&rawtime));
    printf("|============================================================|\n");

    // Output disclaimer to screen
    printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    output_disclaimer_citations(Params, Sim);               // lib/Outputs.c
    output_disclaimer_citations_spatial_cell(Params, Sim);  // lib/Outputs.c
    printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

    free(directory);	
    free(res_dir_full);
    SC_array_deallocation(&SC);         // lib/Spatial_coupling.cpp
    Ca_array_deallocation(&Ca);			// lib/CRU.cpp
    CRU_map_array_deallocation(&CRU);   // lib/CRU.cpp
    for (int n = 0; n < SC.N; n++)  Dyad_array_deallocation(&Dyad[n]);	// lib/CRU.cpp
    delete[] Dyad;
    delete[] SR;
    delete[] MEM;
    delete[] Rand;
    delete[] CRU.TT_map;
    delete[] myofil;
    delete[] activated_switch;
} 
// End Main *************************************************************************************//|

