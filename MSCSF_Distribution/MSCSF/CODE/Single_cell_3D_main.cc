// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Main file for single cell simulations, ======  //
// using the Colman-lab spatial Ca2+ handling system.======  //
// ========================================================  //
// GNU 3 LICENSE TEXT =====================================  //
// COPYRIGHT (C) 2016-2020 MICHAEL A. COLMAN ==============  //
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
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
	// These functions populate Argin struct from command-line arguments
	// All in lib/Arguments.c 
	Argument_parameters Argin;									// Initialise struct 		
	set_argument_defaults(&Argin);								// Sets all boolean switches to false 
	call_argument_functions(argc, argv, &Argin, "Single_cell_3D"); // Reads settings file and calls set args
	// End Argument handling ============================//|

	// Set and assign simulation settings ===============\\|
	// Assigns Sim.{Model, BCL, reference, Beats, Total_time, read/write_state, dt, ..} from Argin struct
	Simulation_parameters Sim;									// Initialise sim parameters struct || lib/Structs.h
	set_simulation_defaults(&Sim, 0.025);						// (sim struct, dt)					|| lib/Initialisation.c
	set_simulation_settings(&Sim, Argin, "integrated");			// Set from arguments				|| lib/Initialisation.c
	// End set and assign settings ======================//|

	// Output files =====================================\\|
	// Create folders by default or references passed in
	char * directory 	= (char*)malloc(500);
	char * results_dir  = (char*)malloc(500);
    char * res_dir_full = (char*)malloc(500); // full path including directory name
	char * params_dir   = (char*)malloc(500);
	char * sr_dir       = (char*)malloc(500);
	char * mkdirectory	= (char*)malloc(500);
	char * mkfile		= (char*)malloc(500);
	if (Argin.reference_arg == true) sprintf(directory, "Outputs_3Dcell_%s", Sim.reference);
	else sprintf(directory, "Outputs_3Dcell");
	sprintf(mkdirectory, "mkdir -p %s", directory);
	system(mkdirectory);
	if (Argin.results_reference_arg == true) sprintf(results_dir, "Results_%s", Sim.results_reference);
	else sprintf(results_dir, "Results");
	sprintf(mkdirectory, "mkdir -p %s/%s", directory, results_dir);
	system(mkdirectory);
    sprintf(res_dir_full, "%s/%s", directory, results_dir);
	sprintf(sr_dir, "Spatial_%s", results_dir);
	sprintf(mkdirectory, "mkdir -p %s/%s", directory, sr_dir);
	system(mkdirectory);
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
	ofstream out_ex(mkfile);     // Contains all AP properties related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/CRU.dat", directory, results_dir);
	ofstream out_cru(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Ca_linescan_Z.dat", directory, results_dir);
	ofstream out_lsZ(mkfile);     // Contains linescan of Ca, longitudinal
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Ca_linescan_Y.dat", directory, results_dir);
	ofstream out_lsY(mkfile);     // Contains linescan of Ca, transverse
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Ca_linescan_X.dat", directory, results_dir);
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
	Cell_parameters					Params;		// Parameters/constants	
	State_variables					State;		// Time-dependent state variables
	Model_variables					Variables;	// Calculated variables
	Ca_variables					Ca;			// Ca for spatial Ca-handling model
	SC_variables                    SC;     	// Spatial coupling
	CRU_variables					CRU;		// sizes and averages for 3D cell model
	Dyad_variables					*Dyad;		// RyR and LTCC flux variables, local dyad params
	SR_fluxes						*SR;		// SERCA and Jleak fluxes variables
	Membrane_fluxes					*MEM;		// NCX, ICab, ICap
	RAND							*Rand;		// radnom number 
	double 							Vm;			// Global voltage

	// Myofilament and force model
	Myofilament             		*myofil; // lib/myofilament.cpp
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
	else if	(strcmp(Params.Model, "hVM_ORD_s") == 0);		// Simplified version of O'Hara et al. 2011 PLOS Comp. Biol. 7, e1002061
	else if	(strcmp(Params.Model, "hAM_CAZ_s") == 0);		// Simplified version of Colman et al. 2013 J. Physiol. 591, 4249-4272
	else
	{
		printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params.Model);
		exit(1);
	}
	// End appropriate check ============================//|

	// Spatial cell setup - settings and allocation =====\\|
	// Set the cellsize and sim_cellsize type, which sets cell dimensions
	spatial_cell_settings(&CRU, Argin);					// lib/CRU.cpp	
	SC_set_array_sizes(&SC, CRU.NX, CRU.NY, CRU.NZ);	// lib/Spatial_coupling.cpp ; sets SC.N from CRU.N
	SC_array_allocation_N3(&SC, SC.NX, SC.NY, SC.NZ);   // Allocates arrays of size NX*NY*NZ
	SC.N	= SC.NX * SC.NY * SC.NZ;					// Idealised cell; no empty space	
	for (int n = 0; n < SC.N; n++) SC.geo[n] = 1;		// Idealised cell; no empty space

	// Array allocation =============\\|
	// Allocates arrays for local values of Ca in each compartment
	Ca_array_allocation(SC.N, &Ca);				// lib/CRU.cpp

	// Allocates neighbour and diff arrays
	SC_array_allocation_Ncell(&SC, SC.N);       // lib/Spatial_coupling.cpp

	// Allocate arrays of structs which are spatially local
	Dyad	= new Dyad_variables[SC.N];		// dyad parameters per CRU
	SR		= new SR_fluxes[SC.N];			// SR fluxes per CRU
	MEM		= new Membrane_fluxes[SC.N];	// membrane fluxes per CRU

	// myofilament and force
	myofil          = new Myofilament[SC.N];
	for (int n = 0; n < SC.N; n++) myofil[n].LSODA_set(); // lib/myofilament.cpp

	printf(">Spatial arrays allocated | NCRU = %d NX = %d NY = %d NZ = %d\n", SC.N, SC.NX, SC.NY, SC.NZ);
	// End Array allocation =========//|

	// Cell index and neighbours
	SC_set_index_and_geo_linear(&SC);               // lib/Spatial_coupling.cpp
	SC_set_neighbours(&SC);                         // lib/Spatial_coupling.cpp
	printf(">Linear index and neighbours set\n");
	// End spatial cell setup ===========================//|

	// Set parameters (defaults and model specific) =====\\|
	// Default modifiers || sets all scale factors to 1 and shifts to 0 so they can be multiplicatively applied by various modifications
	set_modification_defaults_native(&Params);		// lib/Initialisation.c

	// Set default parameters (constants etc); can be overwritten by model-specific later	
	set_default_parameters(&Params);				// lib/Initialisation.c
	printf(">Default parameters set - native\n");

	// Set model specific parameters
	Params.dt = Sim.dt; 	// Set before "set_params" called, which may explicitly set dt, for checking if dt has changed
	set_parameters_native(&Params, Params.Model);	// lib/Model.c and dependants (may set dt for model-specific)

	// Set model condition params from arguments where passed; only for those which are set in set_parameters_native() to overwrite with argument value
	if (Argin.Celltype_arg 	== true)	Params.Celltype 	= Argin.Celltype; 
	if (Argin.ISO_model_arg == true)	Params.ISO_model 	= Argin.ISO_model;
	if (Argin.ACh_model_arg == true)    Params.ACh_model    = Argin.ACh_model; 
	if (Argin.Ihyp_arg      == true)   	Params.AIhyp        = Argin.AIhyp;

    assign_concentrations_from_arguments(&Params, Argin);

	// Update Sim.dt if Params.dt has been explicitly set in "set_parameters" (thus Sim.dt != Params.dt), and dt has NOT been passed as a command-line argument.
	if (Argin.dt_arg    	== false && Params.dt != Sim.dt) 	Sim.dt = Params.dt;
	for (int n =0; n < SC.N; n++) myofil[n].dt_myof = Sim.dt; // dt for LSODA solver same as for whole model
	printf(">Model and version specific parameters set\n");

	// Now set the default and specific integrated Ca2+ handling parameters - overwrites similar parameters set in native
	set_parameters_spatial_Ca_defaults(&Params);		// lib/Initialisation.c
	set_parameters_spatial_Ca(&Params, Params.Model);	// lib/Model.c and dependants

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
	update_heterogeneity_and_modulation_integrated(&Params); // lib/Model.c -> lib/Model_X.cpp	

	// scale channel numbers by expression scale
	// (Grel and GCaL refer to expression i.e. NRyR/NLTCC)
	Params.NRyR_mean 	    *= Params.Grel;
	Params.NLTCC_mean 	    *= Params.GCaL;
	printf(">Heterogeneity and modulation parameters set\n");

	// Assign tau ss, which may be defined by model, modification or argumnet
    if (Argin.tau_ss_arg == true) Params.tau_ss_type = Argin.tau_ss_type; //otherwise default or set in model/modification function
    // Sets actual time constants from type reference ("slow" to "fast")
    set_tau_ss(&Params); // lib/CRU.cpp:
    printf(">Subspace coupling time constants set\n");
    // end set modification =============================//|

    // Membrane capacitance as a function of cell size ==\\|
    Params.Cm			= Params.Cm_CRU * CRU.NTOT_CRUs; // Cm CRU not used in calcs other than this
    CRU.NCRUs_MEM       = CRU.NTOT_CRUs;                 // by default, all CRUs have membrane component
    printf(">Cm total for whole cell (pre detubulation) = %.2f pF\n", Params.Cm);
    // End Membrane capacitance / cell size =============//|

    // Set local scale params from global params ========\\|
    // Allocate arrays, apply global params to local params
    // for sub-cellular heterogeneity.
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
    Rand	= new RAND [SC.N];

    // Initialise stimulus ==============================\\|
    stimulus_setup(Params, &Variables, Sim.dt, Sim.BCL, Sim.S2_CL, Sim.Paced_time); // lib/Model.c
    printf(">Stimulus settings set\n");
    // End initialise stimulus ==========================//|

    // Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params, argc, argv);  // lib/Outputs.c
    output_settings_3D_cell(Params, Sim, CRU, res_dir_full);			               // lib/Outputs.c

    // Setup complete, simulation running ======\\|
    time_t rawtime;
    time (&rawtime);
    printf("|============================================================|\n");
    printf("|Setup complete:\n");
    printf("|Code is now running. Started at %s", ctime (&rawtime));
    printf("|============================================================|\n\n");
    // End Setup complete, simulation running ==//|

    // Setup dyad heterogeneity and allocate local dyad arrays
    if (strcmp(CRU.volds_het, "On") == 0 || strcmp(CRU.RyR_het, "random") == 0 || strcmp(CRU.LTCC_het, "random") == 0)
    { 
        for (int n = 0; n < SC.N; n++)
        {
            // Set NRyR and LTCC from random variation if het = randomm
            // sets vol_ds based on ave and random  (if het = On)
            set_dyad_heterogeneity(Params, &Dyad[n], Rand[n].mtrand1(), n, params_dir, &CRU.dyad_het_map[n], CRU);	// lib/CRU.cpp
        }
        write_random_dyad_het_vtk(SC, CRU.dyad_het_map, directory); // lib/CRU.cpp
    }
    for (int n = 0; n < SC.N; n++)
    {
        Dyad_array_allocation(&Dyad[n]);
    }
    printf("Local NRyR and LLTCC set and dyad arrays allocated\n");

    // Set initial conditions of state variables
    // Function in lib/Model.c calls specific functions in lib/Model_X.cpp
    initial_conditions_native(&State, Params, Params.Model); 						// lib/Model.c

    // Integrated calcium handling conditions
    initial_conditions_calcium(&Ca, Params, SC.N);									// lib/CRU.cpp
    for (int n = 0; n < SC.N; n++) initial_conditions_dyad_stochastic(&Dyad[n]); 	// lib/CRU.cpp

    Vm = State.Vm;
    printf("Initial conditions set\n");

    // Initialise measurement variables and flags
    initialise_measurement_variables(&Variables);                                   // lib/Initialisation.c

    // Output final parameters and calculate and output 
    // voltage-dependant functions, as used in the simulation
    compute_and_output_current_functions(Params, &Variables, params_dir); // lib/Model.c

    // run voltage clamp if option is set
    if(strcmp(Sim.Vclamp, "On") == 0)
    {
        printf("Voltage clamp being performed ....\n");
        run_voltage_clamp_3Dcell(Params, &Variables, &State, Dyad, &Ca, &CRU, &SC, SR, MEM, Rand, params_dir, Sim.dt); // lib/CRU.cpp
        initial_conditions_native(&State, Params, Params.Model);                        // lib/Model.c
        initial_conditions_calcium(&Ca, Params, SC.N);                                  // lib/CRU.cpp
        for (int n = 0; n < SC.N; n++) initial_conditions_dyad_stochastic(&Dyad[n]);
        printf("Voltage clamp finished\n");
    }

    // Read state from file  || this must be after Vclamp as that resets ICs for each step
    if (strcmp(Sim.Read_state, "On") == 0)
    {
        // lib/Read_write_state.c
        Read_state_single_cell_integrated_spatial(&State, Params, Dyad, &Ca, Sim.BCL, PATH, Params.Model, Sim.state_reference_read, SC.NX, SC.NY, SC.NZ, SC.N);
        Ca.CYTO     = State.Cai;
        Ca.SS       = State.Cai_sl;
        Ca.DS       = State.Cai_j;
        Ca.NSR      = State.CanSR;
        Ca.JSR      = State.CajSR;
        printf("Initial conditions / state read in from file\n");
    }
    // Reads in ave values and assigns all local values to the ave    
    else if (strcmp(Sim.Read_state, "ave") == 0)
    {
        Read_state_single_cell_integrated(&State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_read); //lib/Read_write_state.c
        // Now copy global state variables into local CRU and Ca
        for (int n = 0; n < SC.N; n++)
        {
            Ca.cyto[n]		= State.Cai;
            Ca.ss[n]		= State.Cai_sl;
            Ca.ds[n]		= State.Cai_j;
            Ca.nsr[n]		= State.CanSR;
            Ca.jsr[n]		= State.CajSR;
            Dyad[n].Monomer = State.Myo_m;
            Dyad[n].Mi		= State.Myo_c;
        }
        Ca.CYTO 	= State.Cai;
        Ca.SS		= State.Cai_sl;
        Ca.DS		= State.Cai_j;
        Ca.NSR		= State.CanSR;
        Ca.JSR		= State.CajSR;
        printf("Initial conditions / state read in from file\n");
    }
    // End Read state =========================//|
    Vm          = State.Vm;

    // Time loop ================================================================================\\|
    printf("Time loop started:\nTime = %.0fms\n",sim_time);
    for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
    {
        // Assign Ca state variables (seen by ionic model) from integrated whole-cell ave variables
        State.Cai		= 1e-3*Ca.CYTO;		// Ca dependent currents, Cai (in mM not uM)
        State.CanSR		= 1e-3*Ca.NSR;		// Ca dependent currents, Cansr (in mM not uM)
        State.CajSR		= 1e-3*Ca.JSR;		// Ca dependent currents, Cajsr (in mM not uM)

        // Compute stimulus current || lib/Model.c || sets Istims to 0 or stimmag dependant on time
        compute_Istim(Params, &Variables, Sim.Paced_time, Sim.S2_time, sim_time, iteration_counter);  	// lib/Model.c

        // Spatial Ca handling ========================================================\\|
        // First, populate random number array; can be done in ====\\|
        // series or parallel independent of parallelisation of whole model
#pragma omp parallel for default(none) shared(SC, Dyad, Rand) //private(mtrand1)
        for (int n = 0; n < SC.N; n++)
        {
            for (int j = 0; j < Dyad[n].NRyR; j++)  Dyad[n].rand_RyR[j]   	= Rand[n].mtrand1(); // allows faster parallelisation
            for (int j = 0; j < Dyad[n].NLTCC; j++) Dyad[n].rand_LTCC[j]   	= Rand[n].mtrand1(); // as calling mtrand within functions seems slower
        }
        // End populate random number array =======================//|

        // Spatial loop - 1 =======================================\\|
#pragma omp parallel for default(none) shared(SC, Vm, Params, Variables, State, Sim, Ca, sim_time, Dyad, SR, MEM, myofil, Argin)
        for (int n = 0; n < SC.N; n++)
        {
            // Impose CaSR at specified time if argument passed (allows precise setting of CaSR during simulation)
            if (Sim.CaSR_set == false && strcmp(Sim.Delayed_CaSR_IC, "On") == 0 && sim_time >= Sim.CaSR_IC_delay)
            { Ca.nsr[n] = Ca.jsr[n] = Argin.CaSR_IC; Sim.CaSR_set = true; }

            // Zero reaction terms so they can be sequentially modified
            Ca.ss_reac[n] = Ca.cyto_reac[n] = Ca.nsr_reac[n] = Ca.jsr_reac[n] = 0;

            // Inter-compartment transfer || lib/CRU.cpp
            comp_J_ds_ss(Params, Ca.ds[n], Ca.ss[n], Dyad[n].vol_ds, &Ca.ss_reac[n]);
            comp_J_ss_cyto(Params, Ca.ss[n], Ca.cyto[n], &Ca.ss_reac[n], &Ca.cyto_reac[n]);	
            comp_J_nsr_jsr(Params, Ca.nsr[n], Ca.jsr[n], &Ca.nsr_reac[n], &Ca.jsr_reac[n]);	

            // Comp dyad || lib/CRU.cpp -> computes and solves JCaL, Jrel and Cads/jsr fluxes
            comp_dyad_3D(Params, &Dyad[n], Ca.ds[n], Ca.jsr[n], Ca.ss[n] /*to which ds is coupled*/, &Ca.jsr_reac[n], Vm, Sim.dt, Params.Model);

            // Buffering || lib/CRU.cpp
            comp_buffering(Params, &Ca.bcyto[n], &Ca.bss[n], &Ca.bjsr[n], Ca.cyto[n], Ca.ss[n], Ca.jsr[n]);

            // Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
            comp_SR_fluxes(Params, &SR[n], Ca.cyto[n], Ca.nsr[n], &Ca.cyto_reac[n], &Ca.nsr_reac[n]);

            // Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
            comp_membrane_fluxes(Params, &MEM[n], State, Ca.cyto[n], Ca.ss[n], &Ca.cyto_reac[n], &Ca.ss_reac[n], Vm, 1.0);  // 1.0 is SRF mult as not relevant here

            // trpn  || lib/myofilament.cpp || this is general needs to be looked at
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

		// Assign currents for use in AP model from whole-cell averages
		Variables.ICaL 	= CRU.I_CAL;
		Variables.INCX	= CRU.I_NCX_bulk + CRU.I_NCX_ss;
		Variables.ICaP	= CRU.I_CaP_bulk + CRU.I_CaP_ss;
		Variables.ICab	= CRU.I_Cab_bulk + CRU.I_Cab_ss;
		// End Spatial Ca handling ====================================================//|

		// Solve the AP model || lib/Model.c -> lib/Model_X.cpp
		// This sets and updates all gates, and calculates Itot 
		compute_model_integrated(Params, &Variables, &State, Vm, Sim.dt);

		// Update Voltage
		State.Vm	= State.Vm + Sim.dt*(-(Variables.Itot + Variables.Istim + Variables.Istim_S2));

		// Excitation state and measurements | lib/Model.c  | "State.Vm" is voltage at t, "Vm" is voltage at t-dt
		determine_excitation_state(&Variables, Vm, sim_time);							
		calculate_measurement_properties(&Variables, Vm, State.Vm, sim_time, Sim.dt, -70, State.Cai, State.CanSR);		// -70 is APD V threshold	

		// Assign global voltage to state voltage
		Vm			= State.Vm;

		// Output data to files
		if (iteration_counter%(int)(1/Sim.dt) == 0) // if sim_time is an integer (i.e. per ms)
		{
			output_currents(out_cu, sim_time, Variables, State, Vm);		// lib/Outputs.cpp || V, currents, gating variables, concs etc
			output_excitation_properties(out_ex, sim_time, Variables, Vm);	// lib/Outputs.cpp || APD, excitation state, dv/dt etc
			output_CRU(out_cru, sim_time, Ca, CRU, Vm);						// lib/Outputs.cpp || Ca concentrations, Jrel, JCaL, membrane and SR fluxes

			// Spatial data out ===============\\|
			// Linescan || linescan through centre of cell in each direction
			linescan_out_X(out_lsX, SC, Ca.cyto, int(float(SC.NY/2)), int(float(SC.NZ/2))); 		// lib/Outputs.cpp
			linescan_out_Y(out_lsY, SC, Ca.cyto, int(float(SC.NX/2)), int(float(SC.NZ/2))); 		// lib/Outputs.cpp
			linescan_out_Z(out_lsZ, SC, Ca.cyto, int(float(SC.NX/2)), int(float(SC.NY/2))); 		// lib/Outputs.cpp

            // Full 3D spatial data
            if (sim_time >= Sim.Spatial_output_start_time && sim_time <= Sim.Spatial_output_end_time)
            {
                if (Sim.Spatial_output_interval_vtk  > 0) // such that setting to zero means no spatial outputs
                {
                    if (outcount %Sim.Spatial_output_interval_vtk == 0)
                    {
                        vtk_3D_output("Ca", directory, sr_dir, Ca.cyto, SC, outcount);   // every x ms, output vtk file
                        // Feel free to add any new spatial variables here (J_SERCA, CaJSR etc)
                    }
                }
                if (Sim.Spatial_output_interval_data  > 0) // such that setting to zero means no spatial outputs
                {
                    if (outcount %Sim.Spatial_output_interval_data == 0)
                    {
                        array_1D_output("Ca", directory, sr_dir, Ca.cyto, SC, outcount);   // every x ms, output binary data file
                        array_1D_output("CaSR", directory, sr_dir, Ca.nsr, SC, outcount);  // every x ms, output binary data file
                        array_1D_output("CaDS", directory, sr_dir, Ca.ds, SC, outcount);   // every x ms, output binary data file
                        // Feel free to add any new spatial variables here (J_SERCA, CaJSR etc)
                    }
                }
            }
            // End Spatial data out ===========//|
            outcount++; // ms couter	
        }

        iteration_counter ++; // number of steps in dt
        if (iteration_counter%(500*((int)(1/Sim.dt))) == 0) printf("Time = %.0fms\n",sim_time); // output every 500 ms
    }
    // End Time loop ============================================================================//|

    // Print final time in simulation land
    printf("Final Time = %.0fms\n\n",sim_time);

    // Write state 
    if (strcmp(Sim.Write_state, "On") == 0)
    {
        // First, assign averages from CRU model to state for outputs
        State.Cai     = Ca.CYTO;
        State.Cai_sl  = Ca.SS;
        State.Cai_j       = Ca.DS;
        State.CanSR       = Ca.NSR;
        State.CajSR       = Ca.JSR;
        State.Myo_m       = CRU.Monomer;
        State.Myo_c       = CRU.Mi;
        Write_state_single_cell_integrated_spatial(State, Params, Dyad, Ca, Sim.BCL, PATH, Params.Model, Sim.state_reference_write, SC.NX, SC.NY, SC.NZ, SC.N);
        printf("State written to file - full spatial\n");
    }
    else if (strcmp(Sim.Write_state, "ave") == 0)
    {
        // First, assign averages from CRU model to state for outputs
        State.Cai		= Ca.CYTO;
        State.Cai_sl	= Ca.SS;
        State.Cai_j		= Ca.DS;
        State.CanSR		= Ca.NSR;
        State.CajSR		= Ca.JSR;
        State.Myo_m		= CRU.Monomer;
        State.Myo_c		= CRU.Mi;
        Write_state_single_cell_integrated(State, Params, Sim.BCL, PATH, Params.Model, Sim.state_reference_write); //lib/Read_write_state.c
        printf("State written to file - cell aves\n");
    }
    // End Write state

    // Output final beat properties to file and screen || APD, dvdt_max etc
    char * log_reference    = (char*)malloc(500);
    sprintf(log_reference, "%s/Properties_log.dat", directory);
    output_properties_to_screen(log_reference, Variables, Sim);     // lib/Outputs.cpp
    free(log_reference);

    time (&rawtime);
    printf("|============================================================|\n");
    printf("|Code has now finished. Finished at %s", ctime (&rawtime));
    printf("|============================================================|\n");

    // Output disclaimer to screen
    printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    output_disclaimer_citations(Params, Sim);   			// lib/Outputs.c
    output_disclaimer_citations_spatial_cell(Params, Sim);  // lib/Outputs.c
    printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

    // free memory
    free(directory);	
    free(results_dir);	
    free(res_dir_full);
    free(params_dir);	
    free(sr_dir);	
    SC_array_deallocation(&SC);         // lib/Spatial_coupling.cpp
    Ca_array_deallocation(&Ca);			// lib/CRU.cpp
    CRU_map_array_deallocation(&CRU);   // lib/CRU.cpp
    for (int n = 0; n < SC.N; n++)  Dyad_array_deallocation(&Dyad[n]);	// lib/CRU.cpp
    delete[] Dyad;
    delete[] SR;
    delete[] MEM;
    delete[] Rand;
    delete[] myofil;
} 
// End Main *************************************************************************************//|

