// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Main file for tissue-scale simulations, =====  //
// using the Colman-lab spatial Ca2+ handling system, with=  //
// non-spatial model reduction and spontaneous release ====  //
// function implementation. ==============================   //
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
#include "lib/Spatial_coupling.h"
#include "lib/Tissue.h"
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
	printf("|Model version: Tissue - 0D - integrated Ca handling ========|\n");
	printf("|============================================================|\n");
	printf("\n");
	// End output to screen =============================//|

	// Command line argument initialise and read in =====\\|
	// These functions populate Argin struct from command-line arguments
	// All in lib/Arguments.c
	Argument_parameters Argin;                                  // Initialise struct                            || lib/Arguments.c
	set_argument_defaults(&Argin);                              // Sets all boolean switches to false           || lib/Arguments.c
	call_argument_functions(argc, argv, &Argin, "Tissue_integrated"); // Reads settings file and calls set args || lib/Arguments.c
	// End Argument handling ============================//|

	// Set and assign simulation settings ===============\\|
	// Assigns Sim.{BCL, reference, Beats, Total_time, read/write_state, dt, ..} from Argin struct
	Simulation_parameters Sim;									// Initialise sim parameters struct || lib/Structs.h
	set_simulation_defaults(&Sim, 0.01);						// (sim struct, dt)					|| lib/Initialisation.c
	set_simulation_settings(&Sim, Argin, "integrated");			// Set from arguments				|| lib/Initialisation.c
	// End set and assign settings ======================//|

	// Output files =====================================\\|
	// Create folders by default or references passed in
	char * directory 	= (char*)malloc(500);
	char * results_dir  = (char*)malloc(500);
    char * res_dir_full = (char*)malloc(500); // full path including directory name
	char * sr_dir       = (char*)malloc(500);
	char * mkdirectory	= (char*)malloc(500);
	char * mkfile		= (char*)malloc(500);
	if (Argin.reference_arg == true) sprintf(directory, "Outputs_0Dtissue_%s", Sim.reference);
	else sprintf(directory, "Outputs_0Dtissue");
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

	// Now create actual output files
	printf(">Creating output files...\n");

	sprintf(mkfile, "%s/%s/Currents_cell1.dat", directory, results_dir);
	ofstream out_cu(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/CRU_cell1.dat", directory, results_dir);
	ofstream out_cru1(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Properties_cell1.dat", directory, results_dir);
	ofstream out_ex(mkfile); 	 // Contains all AP properties related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Currents_cell2.dat", directory, results_dir);
	ofstream out_cu2(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/CRU_cell2.dat", directory, results_dir);
	ofstream out_cru2(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Properties_cell2.dat", directory, results_dir);
	ofstream out_ex2(mkfile);     // Contains all AP properties related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Currents_cell3.dat", directory, results_dir);
	ofstream out_cu3(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/CRU_cell3.dat", directory, results_dir);
	ofstream out_cru3(mkfile);    // Contains all current related outputs
	printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Properties_cell3.dat", directory, results_dir);
	ofstream out_ex3(mkfile);     // Contains all AP properties related outputs
	printf("\t %s\n", mkfile);

    sprintf(mkfile, "%s/%s/SRF_properties.dat", directory, results_dir);
    ofstream out_srf_prop(mkfile);    // Contains a list of ti and NRyR of actually occuring SRF
    printf("\t %s\n", mkfile);

	sprintf(mkfile, "%s/%s/Vm_linescan_x.dat", directory, results_dir);
	ofstream out_ls(mkfile);     // Contains linescan of Vm
	printf("\t %s\n", mkfile);

	printf("\n");
	free(mkfile);
	free(mkdirectory);
	// End Output files =================================//|

	// Initialise simulation structs and variables ======\\|
	double		sim_time				= 0;
	int			iteration_counter		= 0; 	// Number of iterations over dt
	int			outcount				= 0; 	// Output number reference (iterates per ms)
	int 		phase_counter			= 0;	// counter for phase outputs if run

	// All below defined in lib/Structs.h
	Cell_parameters					Params_global;	// Parameters/constants || Global settings
	Cell_parameters					*Params;		// Parameters/constants || Cell-specific settings
	State_variables					*State;			// Time-dependent state variables
	Model_variables					*Variables;		// Calculated variables
	Ca_variables                    *Ca;        	// Ca for spatial Ca-handling model (Ca_i/ss/ds/nsr/jsr; reactions; buffering)
	CRU_variables                   *CRU;   		// cell dimensins; cell geo specifiers; averages of local values
	Dyad_variables                  *Dyad;  		// RyR and LTCC flux variables, local dyad params
	SR_fluxes                       *SR;        	// SERCA and Jleak fluxes variables
	Membrane_fluxes                 *MEM;   		// NCX, ICab, ICap
	SC_variables					SC;				// Spatial coupling (neighbour maps, D arrays, coupling functions)
	Tissue_parameters				Tissue;			// Tissue settings (tissue model and dimension, array sizes, diffusion params, anisotropy etc)
	Spontaneous_release_functions   *SRF;    		// Spontaneous release function variables and parameters
	RAND                            *Rand;   		// random number array
	double 							*Vm;			// Global copy of voltage

	// Myofilament and force model
	Myofilament                     *myofil;        // lib/myofilament.cpp

	// For Ca spatial outputs
	double 							*Cai;			// Global copy of Cai
	double 							*CaSR;			// Global copy of CaSR
	printf(">Variables and structs declared\n");
	// End Initialise simulation structs and variables ==//|

	// Set model condition parameters from arguments | lib/Initialisation.c
	// Sets Params.{Model, modulation(ISO,remodelling, mutation etc), model group} to defaults
	// then overwrites with Argin values, if arguments have been passed
	set_model_conditions(&Params_global, Argin);	// Set Model, modulation (ISO, remodelling etc), environment to default or direct argument controlled
	// set_model_group_variables now below, after default model check

	// Check Model is appropriate for this code version =\\|
	if      (strcmp(Params_global.Model, "minimal") == 0);		// hybrid minimal model; Colman 2018
	else if (strcmp(Params_global.Model, "hVM_ORD_s") == 0);	// human atrial cell; Courtemanche et al 1998 Am J Physiol 1998; 275(1 Pt 2):H301-21.
	else if (strcmp(Params_global.Model, "hAM_CAZ_s") == 0);	// human atrial cell; Colman et al 2018 Front Physiol 9:1211; CRN Ca2+ handling
	else
	{
		printf("ERROR: \"%s\" is not a valid Model selection for this code version\n", Params_global.Model);
		exit(1);
	}
	// End appropriate check ============================//|

	// Tissue setup - geometry, arrays and memory =======\\|
	// Sets tissue conditions (Tissue_order, Tissue_model, Tissue_type, Orientation_type, maps on/off etc)
	// to default values, then assigns any from args (where passed in).
	// Some of these define which conditions are selected in tissue setup; some are related to implementation
	// at run time. Most can be overwritten by tissue model defaults and/or arguments
	set_tissue_model_conditions(&Tissue, Argin);	// lib/Tissue.cpp

	// Set specific tissue settings (array sizes, diffusion and fibre parameters, stimulus coords and type etc)
	// lib/Tissue.cpp
	if (strcmp(Tissue.Tissue_order, "geo") == 0) set_tissue_settings_anatomical(Params_global, &Tissue); // for models with realistic geometries (filenames etc)
	else set_tissue_settings_idealised(Params_global, &Tissue);	// for idealised models

	// Set model to tissue default unless Model argument is passed
	// If Tissue.Default_model is "none", single cell default will be used
	if (strcmp(Tissue.Default_model, "none") != 0 && Argin.Model_arg == false) Params_global.Model = Tissue.Default_model;
	set_model_group_variables(&Params_global, Argin); // Model dependent so can't be called earlier || // lib/Initialisation.c

	// set stimulus from predefined type if not already defined by default, or if argument overwrites default
	set_coord_stim_and_map_from_defined_type(Params_global, &Tissue, Argin);  // lib/Tissue.cpp

	// Set global orientation if ideal and set by directional argument (rather than explicitly and locally from file)
	if (strcmp(Tissue.Tissue_order, "geo") != 0) set_global_orientation_direction_from_arg(Params_global, &Tissue, Argin); // lib/Tissue.cpp

	// Set to arguments if expliticlty passed to overwite defaults (inc model specifc) || assuming you'd only pass an argument if you want it to be used!
	overwrite_tissue_properties_from_args(Params_global, &Tissue, Argin); // lib/Tissue.cpp

	// set CV cell indexes -> location of cells for conduction velocity calculation, if CV tissue model is selected
	if (strcmp(Tissue.Tissue_model, "conduction_velocity") == 0) set_CV_cells(Params_global, &Tissue);  // lib/Tissue.cpp

	// IF isotropic:
	// Set D1 to Diso and AR to 1 such that D2, used for isotropic diffusion, will = Diso after being set as D1/AR
	if (strcmp(Tissue.Orientation_type,"isotropic") == 0) { Tissue.D1 = Tissue.Diso; Tissue.D_AR = 1.0; }

	// Allocate arrays || lib/Spatial_coupling.cpp
	// Tissue.NX set by tissue model; passed into Spatial Coupling now finalised
	SC_set_array_sizes(&SC, Tissue.NX, Tissue.NY, Tissue.NZ); 	// Sets array sizes in SC struct from tissue settings
	SC_array_allocation_N3(&SC, SC.NX, SC.NY, SC.NZ); 			// Allocates arrays of size NX*NY*NZ || geo and 3D->1D geo index
	printf(">Spatial coupling NX*NY*NZ arrays allocated\n");

	// Read/Create geometry and calculate Ncell
	select_tissue_geometry_function(Tissue, &SC, PATH, directory); 	// lib/Tissue.cpp
	printf("\tGeometry size (X*Y*Z, %d * %d * %d) || Ncells = %d\n\n", SC.NX, SC.NY, SC.NZ, SC.N);

	// Allocate arrays size Ncell
	SC_array_allocation_Ncell(&SC, SC.N);		// lib/Spatial_coupling.cpp || geo_linear, D arrays, neighbour map, orientation, laplacian components
	tissue_array_allocation(&Tissue, SC.N);		// lib/Tissue.cpp || stim/ISO/remodelling etc map arrays
	printf(">Spatial coupling Ncell arrays allocated\n");

	// Allocate model structs and variable arrays || these are size N as ony require entries for real tissue, not all space
	Params		= new Cell_parameters[SC.N];
	State		= new State_variables[SC.N];
	Variables	= new Model_variables[SC.N];
	Ca          = new Ca_variables[SC.N];
	CRU         = new CRU_variables[SC.N];
	Dyad        = new Dyad_variables[SC.N];
	SR          = new SR_fluxes[SC.N];
	MEM         = new Membrane_fluxes[SC.N];
	SRF			= new Spontaneous_release_functions[SC.N];
	Vm			= new double[SC.N];
	Cai			= new double[SC.N];
	CaSR		= new double[SC.N];
	myofil      = new Myofilament[SC.N];
	for (int n = 0; n < SC.N; n++) myofil[n].LSODA_set(); // lib/myofilament.cpp
	printf(">Ncell struct arrays allocated\n");

	// Cell index and neighbours (geo_index[3D_ref] returns 1D ref; geo_3D_index[1D_ref] returns 3D_ref; geo_linear[1D_ref] = geo[3D_ref]
	SC_set_index_and_geo_linear(&SC);				// lib/Spatial_coupling.cpp
	SC_set_neighbours(&SC);							// lib/Spatial_coupling.cpp
	printf(">Linear index and neighbours set\n");

	// Create stimulus area
	// Check for multi site/times stimulus setting; otherwise call regular function
	if (strcmp(Tissue.Multi_stim, "On") == 0) select_stimulus_area_function_multi_stim(&Tissue, SC, PATH, directory, Sim.S2_CL, Tissue.Nstims);
	else select_stimulus_area_function(&Tissue, SC, PATH, directory, Sim.S2_CL);	// lib/Tissue.cpp
	printf("\tStimulus area created | Nstim = %d\n", Tissue.Nstim); 
	if (Sim.S2_CL != 0) printf("\tS2 Stimulus area created | Nstim = %d\n", Tissue.Nstim_S2);

	// Assign refs for single cell outputs || these are the 1D cell refs for the 3 cells chosen to output detailed data
	int cell1ref, cell2ref, cell3ref;
	cell1ref = 5;
	cell2ref = int(float(SC.N/3)); // in idealised model, not likley to be x-edge (/2, 4 or 5 is)
	cell3ref = SC.N - 5;

	// setup diffusion coefficient arrays
	set_D_dx_global(&SC, Tissue.dx, Tissue.dy, Tissue.dz, Tissue.D1, Tissue.D_AR); // sets D and dx from tissue mdoel settings || lib/Spatial_coupling.cpp

	// Fibre orientation
	set_orientation(&SC, Tissue, PATH, Tissue.Tissue_order);    // lib/Tissue.cpp || This sets orientation to 0, then sets/reads in IF set to anisotropic
	if (strcmp(Tissue.Orientation_type, "anisotropic") == 0) output_fibre_orientation(SC, Tissue, directory); // lib/Tissue.cpp

	// Baseline, non-uniform Dscale (celltype or map; map for geo only)
	// This is for regional or continuous/complex gradient in D1 and/or DAR (inherehent to tissue model)
	if (strcmp(Tissue.D_uniformity, "uniform") != 0) update_D_arrays_Dscale_baseline(&SC, &Tissue,PATH, directory); // lib/Tissue.cpp

	// Modification Dscale (homogeneous or map; ideal or geo)
	// This is for scaling D1 and/or DAR locally associated with modulation (e.g. remodelling)
	update_D_arrays_Dscale_mod(&SC, &Tissue, PATH, directory); // lib/Tissue.cpp

	// now set the D components spatial array from D1 and D2 arrays and orientation
	set_D_array_anisotropic(&SC);                       // lib/Spatial_coupling.cpp
	output_D1_and_D2(SC, directory);                    // lib/Spatial_coupling.cpp  || Outputs vtk for inspection

	// Phase re-entry map
	if (strcmp(Sim.Read_state, "phase") == 0) // needs to create phase map if reading phase ICs
		create_read_phasemap(Tissue.phasemap, Tissue, SC, PATH, directory);
	// End tissue setup - geometry, arrays and memory ===//|

	// Global and local settings from maps etc ==========\\|
	// Default modifiers || sets all scale factors to 1 and shifts to 0 so they can be multiplicatively applied by various modifications
	set_modification_defaults_native(&Params_global);       			// lib/Initialisation.c

	// Set default parameters (constants etc); can be overwritten by model-specific later
	// Set for Global params (which may determine other options); later actually set for each cell
	set_default_parameters(&Params_global);								// lib/Initialisation.c	
	set_parameters_spatial_Ca_defaults(&Params_global);    				// lib/Initialisation.c
	set_parameters_spatial_Ca(&Params_global, Params_global.Model); 	// lib/Model.c and dependants
	update_parameters_spatial_Ca_0D(&Params_global);       				// lib/Initialisation.c

	// Set model condition parameters local defaults from global | lib/Initialisation.c
	for (int n = 0; n < SC.N; n++) set_local_model_conditions(Params_global, &Params[n]); // Set Model, ISO, remodelling etc locally

	// Update conditions locally for heterogeneous conditions (e.g. ISO/remodelling maps etc)
	// ISO
	if (strcmp(Tissue.ISO_map_on, "On") == 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.ISO_map, Tissue.ISO_map_file, "ISO"); // lib/Tissue.cpp
		for (int n = 0; n < SC.N; n++) Params[n].ISO *= Tissue.ISO_map[n]; // NOTE: multiplies global ISO value (which local is defaulted to) by local ISO SCALING
		// i.e. if local map is 1, then Params.ISO = global ISO value; if local map is 0, then Params.ISO = 0
		printf(">ISO map read and local ISO conc set\n");
	}
	// Remodelling
	if (strcmp(Tissue.remod_map_on, "On") == 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.remod_map, Tissue.remod_map_file, "remodelling"); // lib/Tissue.cpp
		for (int n = 0; n < SC.N; n++) Params[n].Remodelling_prop *= Tissue.remod_map[n]; // Again, multiplies global value by local map value, thus scaling between 0 and global value
		printf(">Remodelling map read and local remodelling condition set\n");
	}
	// ACh
	if (strcmp(Tissue.ACh_map_on, "On") == 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.ACh_map, Tissue.ACh_map_file, "ACh"); // lib/Tissue.cpp
		for (int n = 0; n < SC.N; n++) Params[n].ACh *= Tissue.ACh_map[n]; // NOTE: multiplies global ACh value (which local is defaulted to) by local ACh SCALING
		printf(">ACh map read and local ACh conc set\n");
	}
	//SRF
	if (strcmp(Tissue.SRF_map_on, "On") == 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.SRF_map, Tissue.SRF_map_file, "SRF"); // lib/Tissue.cpp
		printf(">SRF map read\n");
	}
	// Direct_modulation || this refers to command line DC mods such as "ICaL_scale" etc (applies to all DC mods passed)
	if (strcmp(Tissue.Direct_modulation_map_on, "On") == 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.Direct_modulation_map, Tissue.Direct_modulation_map_file, "Direct_modulation"); // lib/Tissue.cpp
		printf(">Direct_modulation map read\n");
	}
	// Spatial gradient
	if (Argin.spatial_gradient_arg == true) Tissue.spatial_gradient_map_on = Argin.spatial_gradient; // else use default set by tissue model
	if (strcmp(Tissue.spatial_gradient_map_on, "none") == 0) Tissue.spatial_gradient_map_on = "Off"; //will only = none if arg is set to none to turn it off
	if (strcmp(Tissue.spatial_gradient_map_on, "Off") != 0)
	{
		create_or_read_map_double(&Tissue, SC, PATH, directory, Tissue.spatial_gradient_map, Tissue.spatial_gradient_map_file, "Spatial_gradient"); // lib/Tissue.cpp
		for (int n = 0; n < SC.N; n++)
		{
			Params[n].spatial_gradient = Tissue.spatial_gradient_map_on;
			Params[n].spatial_gradient_prop = Tissue.spatial_gradient_map[n]; // note: sets directly from map value
		}
		printf(">Spatial_gradient map read\n");
	}
	else // spatial gradient is off, so local param must be set to none
	{
		for (int n = 0; n < SC.N; n++)
		{
			Params[n].spatial_gradient = "none";
			Params[n].spatial_gradient_prop = 0;
		}
	}
	// End Global and local settings from maps etc ======//|

	// Loop of tissue for cell-by-cell setup ======================\\|
	for (int n = 0; n < SC.N; n++)
	{
		// Set parameters (defaults and model specific) =====\\|
		// Set cellsize and spontaneous release function defaults
		spatial_cell_settings(&CRU[n], Argin);                 // lib/CRU.cpp
		set_SRF_defaults(&SRF[n], Argin);                      // lib/Spontaneous_release_functions.cpp
		if (strcmp(Tissue.SRF_map_on, "On") == 0) if (Tissue.SRF_map[n] == 0.0) SRF[n].Mode = "Off"; // set non SRF regions to be Off || others will be global SRF parameters

		// Default modifiers || sets all scale factors to 1 and shifts to 0 so they can be multiplicatively applied by various modifications
		set_modification_defaults_native(&Params[n]);		// lib/Initialisation.c

		// Set default global parameters (may want to overwrite a modifier here, hence defaulted above)
		if (n == 0) printf(">Setting default parameters...\n");
		set_default_parameters(&Params[n]);					// lib/Initialisation.c
		if (n == SC.N -1) printf(">Default parameters set\n");

		// Set model specific parameters
		Params[n].dt = Sim.dt; 	// Set before "set_params" called, which may explicitly set dt, for checking if dt has changed

		// Select local baseline model if multiple models is on
		if (strcmp(Tissue.Multiple_models, "On") == 0) 
		{
			if (strcmp(Tissue.Tissue_type, "homogeneous") == 0)
			{
				printf("ERROR: Multiple Models cannot be run with homogeneous Tissue_type; heterogeneity must exist to assign regions to two Models!\n");
				exit(1);
			}
			// For regions assigned "Model_2", set local Model to the entry held in "Tissue_model_2" (set in Tissue model settings or by argument)
			// No need to do anything for regions assigned "Model_1" as this is what Params[n].Model already contains
			if (strcmp(Tissue.Modeltype_number[SC.geo_linear[n]], "Model_2") == 0)  Params[n].Model = Tissue.Tissue_model_2; 
		}

		set_model_group_variables(&Params[n], Argin);  // Model dependent so needs to be called here (as Params[n].Model hmay have changed

		// Set default parameters (constants etc); can be overwritten by model-specific later
		// Can be set now that local conditions (Model, ISO, Remodelling etc) have be set
		set_parameters_native(&Params[n], Params[n].Model);	// lib/Model.c and dependants (may set dt for model-specific)

		// Set model condition params from arguments where passed; only for those which are set in set_parameters_native() to overwrite with argument value
		if (Argin.Celltype_arg 	== true)	Params[n].Celltype 	= Argin.Celltype; 	
		if (Argin.ISO_model_arg == true)	Params[n].ISO_model	= Argin.ISO_model; 	
		if (Argin.ACh_model_arg == true)    Params[n].ACh_model = Argin.ACh_model; 

		// Update Sim.dt if Params.dt has been explicitly set in "set_parameters" (thus Sim.dt != Params.dt), and dt has NOT been passed as a command-line argument.
		if (Argin.dt_arg    	== false && Params[n].dt != Sim.dt) 	Sim.dt = Params[n].dt;
		myofil[n].dt_myof       = Sim.dt; // dt for LSODA solver same as for whole model
		if (n == SC.N -1) printf(">Model and version specific parameters set\n");

		// Now set the default and specific integrated Ca2+ handling parameters - overwrites similar parameters set in native
		set_parameters_spatial_Ca_defaults(&Params[n]);    		// lib/Initialisation.c
		set_parameters_spatial_Ca(&Params[n], Params[n].Model); // lib/Model.c and dependants
		update_parameters_spatial_Ca_0D(&Params[n]);       		// lib/Initialisation.c

		// Overwrite initial conditions of Cai and CaSR if argument passed
		if (Argin.Cai_IC_arg    == true)    Params[n].Cai          = Argin.Cai_IC;
		if (Argin.CaSR_IC_arg   == true)    Params[n].CaSR         = Argin.CaSR_IC;
		if (n == SC.N -1) printf(">Default parameters set - integrated Ca2+ handling\n");
		// End set parameters (defaults and model specific) =//|

		// Set current modification ==========================\\|
		// lib/Initialisation.c || sets the mod variables (Gx, x_shift/tau_scale etc) from arguments
		// These are the "DC" modulation variables, controlled by Tissue.Direct_modulation_map_on if set
		// If map is off, then simply update modifiers with argument values for all tissue; if map is on, do it only within map region
		if (strcmp(Tissue.Direct_modulation_map_on, "On") == 0)  // if map is on, only apply to map regions
		{
			if (Tissue.Direct_modulation_map[n] > 0.0) assign_modification_from_arguments(&Params[n], Argin); 
		}
		else assign_modification_from_arguments(&Params[n], Argin); // lib/Initialisation.c || assign to all nodes 

		// Now set local celltype if heterogeneity is on (else all cells will be set to default or argument celltype)
		// SC.geo_linear[n] = cellnumber at n; celltype_number[cellnumber] = string of celltype defined in Tissue model settings
		if (strcmp(Tissue.Tissue_type, "heterogeneous") == 0) Params[n].Celltype = Tissue.celltype_number[SC.geo_linear[n]];


		// Now call set heterogneiety and modulation (local celltype and modulation are already set by maps if relevant)
		// lib/Model.c  -> lib/Model_X.cp; calls functions which set modification variables for het and modulation
		// Updates the modifier variables (scales, shifts etc) using the defined settings for any/all het and modulation	
		set_heterogeneity_and_modulation_native(&Params[n]);		// lib/Model.c
		update_heterogeneity_and_modulation_integrated(&Params[n]); // lib/Model.c

		// scale channel numbers by expression scale    
		// (Grel and GCaL refer to expression i.e. NRyR/NLTCC)
		Params[n].NRyR_mean    *= Params[n].Grel; // "Grel/CaL" refers to scaling expression and thus corresponds to Nx
		Params[n].NLTCC_mean   *= Params[n].GCaL;
		if (n == SC.N -1) printf(">Heterogeneity and modulation parameters set\n");
		// end set current modification ======================//| 

		// Membrane capacitance as a function of cell size ==\\|
		Params[n].Cm           = Params[n].Cm_CRU * CRU[n].NTOT_CRUs;
		if (n == SC.N -1) printf(">Cm total for whole cell = %.2f pF\n", Params[n].Cm);
		// End Membrane capacitance / cell size =============//|

		// Local variables from global parameters
		// Remnant of 3D model; here just assigns Dyad, Mem and SR variables from Params
		set_sub_cellular_local_scale(Params[n], &Dyad[n], &MEM[n], &SR[n]); 
	}
	// End loop of tissue for cell-by-cell setup ==================//|

	// Initialise stimulus ==============================\\|
	// Stimulus settings use Params and Variables[0], but do not correspond to cell at element 0
	// Cells to apply stimulus is determined by stimulus map
    if (strcmp(Tissue.Multi_stim, "On") == 0) Sim.Paced_time += Tissue.stim_delay[Tissue.Nstims-1]; // add final stimulus delay to paced time
	stimulus_setup(Params[0], &Variables[0], Sim.dt, Sim.BCL, Sim.S2_CL, Sim.Paced_time); // lib/Model.c

	// If multi timed stim is on, we use Params and Variables for all of the stims
	// Again, the stimulus settings in Params[x] and Variables[x] do not correspond to those cells
	// Cells to apply stimulus is determined by stimulus map
	if (strcmp(Tissue.Multi_stim, "On") == 0) for (int n = 1; n < Tissue.Nstims; n++) stimulus_setup(Params[n], &Variables[n], Sim.dt, Sim.BCL, Sim.S2_CL, Sim.Paced_time); // lib/Model.c
	printf(">Stimulus settings set\n");
	// End initialise stimulus ==========================//|

	// Spontaneous release functions ==========\\|
	for (int n = 0; n < SC.N; n++) 
	{
		SRF_setup(&SRF[n], Argin); // lib/Spontaneous_release_functions.cpp -> calls appropriate set SRF parameters function
	}
	printf("Spontaneous release function parameters set.\n");
	// End spontaneous release functions ======//|

	// Output settings to screen and file || done here so can output actual settings (rather than inputs) for confidence
	Params_global.Celltype 	= Params[0].Celltype; // Just to output baseline model, independent of which model is in cell 0 if multi models used
	Params_global.CaSR 		= Params[0].CaSR;
	Params_global.Cai		= Params[0].Cai;
    Params_global.ISO_model = Params[0].ISO_model;
    Params_global.ACh_model = Params[0].ACh_model;
	assign_modification_from_arguments(&Params_global, Argin);					                // lib/Outputs.c
    output_settings(Sim, res_dir_full, Argin.DC_current_mod_arg, Params_global, argc, argv);    // lib/Outputs.c
	output_settings_tissue(Sim, Tissue, res_dir_full);								                // lib/Outputs.c
	output_settings_0D_cell(Params_global, Sim, CRU[0], res_dir_full, SRF[0]);		                // lib/Outputs.c

	// Setup complete, simulation running ======\\|
	time_t rawtime;
	time (&rawtime);
	printf("|============================================================|\n");
	printf("|Setup complete:\n");
	printf("|Code is now running. Started at %s", ctime (&rawtime));
	printf("|============================================================|\n\n");
	// End Setup complete, simulation running ==//|

	printf("Setting initial conditions....\n");
	// Set initial conditions of state variables
	for (int n = 0; n < SC.N; n++)
	{
		// Setup dyad params to global params as just one CRU per cell (here is where dyad het set in 3D)
		Dyad[n].vol_ds      = Params[n].vds_CRU_mean;
		Dyad[n].NRyR        = Params[n].NRyR_mean;
		Dyad[n].NLTCC       = Params[n].NLTCC_mean;
		if (n == SC.N -1) printf("NRyR and LLTCC set\n");

		// Set initial conditions of state variables
		// Function in lib/Model.c calls specific functions in lib/Model_X.cpp
		initial_conditions_native(&State[n], Params[n], Params[n].Model); 	

		// Integrated calcium handling conditions (0D)
		initial_conditions_calcium_0D(&Ca[n], Params[n]);                   // lib/CRU.cpp 
		initial_conditions_dyad_det(&Dyad[n]);                              // lib/CRU.cpp

		Vm[n] = State[n].Vm;

		// Initialise measurement variables and flags
		initialise_measurement_variables(&Variables[n]);                    // lib/Initialisation.c
		MEM[n].NCX_SRF_mult        = 1.0;	// Normalised except for when needed
		Dyad[n].Krel_SRF_mult      = 1.0;
	}
	printf("Initial conditions set\n");

	// Read state
	if (strcmp(Sim.Read_state, "On") == 0)        // Reads whole tissue
	{
		// Reads whole tissue -> state file must have been written using same tissue model!!
		Read_state_tissue_integrated_whole_tissue(State, Params, Sim.BCL, PATH, Params_global.Model, SC.N, Tissue.Tissue_order, Tissue.Tissue_model, Tissue.Tissue_type, Tissue.Orientation_type, Sim.state_reference_read); //lib/Read_write_state.c
		for (int n = 0; n < SC.N; n++)
		{
			assign_CRU_variables_from_state_read(&Dyad[n], &Ca[n], State[n]); // lib/CRU.cpp
		}
		printf("Initial conditions / state read in from file - whole tissue\n");
	}
	else if (strcmp(Sim.Read_state, "Off") != 0) // other versions of read state
	{
		for (int n = 0; n < SC.N; n++)
		{
			// Reads in file written by single cell model to all tissue (needs file for each celltype and condition present)
			if (strcmp(Sim.Read_state, "single_cell") == 0) 	Read_state_single_cell_integrated_0D(&State[n], Params[n], Sim.BCL, PATH, Params[n].Model, Sim.state_reference_read); 		//lib/Read_write_state.c

			// Reads state from just one coupled cell to whole tissue (same as single cell except written by coupled)
			else if (strcmp(Sim.Read_state, "ave") == 0)		Read_state_tissue_integrated_ave_tissue(&State[n], Params[n], Sim.BCL, PATH, Params[n].Model, Sim.state_reference_read); 	//lib/Read_write_state.c 

			// Reads in single cell phase file into tissue for phase re-entry
			else if (strcmp(Sim.Read_state, "phase") == 0)		Read_state_phase(&State[n], Params[n], Sim.BCL, PATH, Params[n].Model, Tissue.phasemap[n], Sim.state_reference_read); 		//lib/Read_write_state.c

			assign_CRU_variables_from_state_read(&Dyad[n], &Ca[n], State[n]); // lib/CRU.cpp
		}
		if (strcmp(Sim.Read_state, "single_cell") == 0) 	printf("Initial conditions / state read in from file - single cell to tissue\n");
		if (strcmp(Sim.Read_state, "ave") == 0) 			printf("Initial conditions / state read in from file - coupled cell to tissue\n");
		if (strcmp(Sim.Read_state, "phase") == 0)			printf("NOTE:: As phase re-entry, have you set Beats = 0 (and Total_time = x) to ensure no applied stimuli??\n");
		if (strcmp(Sim.Read_state, "phase") == 0) 			printf("Initial conditions / state read in from file - phase distribution\n");
	}
	// Error check for phase write
	if (strcmp(Sim.Write_state, "phase") == 0 && strcmp(Tissue.Tissue_order, "1D") != 0)
	{
		printf("ERROR: Writing phasemap MUST be performed on a 1D strand; pass in \"Tissue_order 1D\"\n");
		exit(1);
	}
	for (int n = 0; n < SC.N; n++) Vm[n] = State[n].Vm;

	// Calculate diffusion tensor differentials and laplacian =====\\|
	printf("Calculating d differential and laplacian\n");
	for (int n = 0; n < SC.N; n++)
	{
		calc_dD_anisotropic_3D(&SC, n); // lib/Spatial_coupling.cpp
		calc_laplacian_and_BCs(&SC, n); // lib/Spatial_coupling.cpp
	}
	// End Calculate diffusion tensor differentials and laplacian =//|

	// Rand array allocation (here so all files have already been read in and checked, rather than spending time here only to throw out an error later)
	bool SRF_check = false;
	for (int n = 0; n < SC.N; n++) if (strcmp(SRF[n].Mode, "Off") != 0) SRF_check = true; // any nodes have non-Off SRF?
	if (SRF_check == true) 
	{
		printf("Allocating rand array; this may take a while (but significantly improves parallelisation performance)\n");
		Rand    = new RAND [SC.N];
	}

	// Spontaneous release functions || het ===\\|
	if (strcmp(SRF[0].SRF_het, "Off") != 0)
	{
		for (int n = 0; n < SC.N; n++) 
		{	
			// Set heterogneity of SRF model
			SRF_tissue_heterogeneity(&SRF[n], Rand[n].mtrand1()); 	// lib/Spontaneous_release_functions.cpp

			// now set SRF params local from local SRF model
			SRF_setup(&SRF[n], Argin);								// lib/Spontaneous_release_functions.cpp
		}
		printf("Spontaneous release function parameter tissue heterogeneity set.\n");
	}
	// End spontaneous release functions ======//|

	// Time loop ================================================================================\\|
	printf("Time loop started:\nTime = %.0fms\n",sim_time);
	for (sim_time = 0.0; sim_time <= (float)Sim.Total_time; sim_time += Sim.dt)
	{
		// Compute stimulus current || lib/Model.c || sets Istims to 0 or stimmag dependant on time
		// Note: outside of tissue loop as indexes do not correspond with cell indexes
		compute_Istim(Params[0], &Variables[0], Sim.Paced_time, Sim.S2_time, sim_time, iteration_counter);
		if (strcmp(Tissue.Multi_stim, "On") == 0) for (int m = 1; m < Tissue.Nstims; m++) compute_Istim(Params[m], &Variables[m], Sim.Paced_time, Sim.S2_time, sim_time, iteration_counter - Tissue.stim_delay[m]*(int)(1.0/Sim.dt));

		// Impose CaSR at specified time if argument passed (allows precise setting of CaSR during simulation)
		if (Sim.CaSR_set == false && strcmp(Sim.Delayed_CaSR_IC, "On") == 0 && sim_time >= Sim.CaSR_IC_delay)
		{ for (int n = 0; n < SC.N; n++) { Ca[n].NSR = Ca[n].JSR = Argin.CaSR_IC; Ca[n].CYTO = Ca[n].SS = Ca[n].DS = Argin.Cai_IC; } Sim.CaSR_set = true; }

		// Loop over all tissue - 1 ===============================\\|
#pragma omp parallel for default(none) shared(SC, Vm, Params, Variables, State, Sim, Tissue, sim_time, Dyad, MEM, SR, CRU, Ca, SRF, Rand, Argin, myofil)
		for (int n = 0; n < SC.N; n++)
		{
			// Compute spatial differential || lib/Spatial_coupling.cpp
			// calculates "SC.diff[n]" 
			calc_diff_from_lap(&SC, Vm, n);

			// Assign Ca state variables (seen by ionic model) from integrated whole-cell ave variables
			State[n].Cai       = 1e-3*Ca[n].CYTO;     // Ca dependent currents, Cai (in mM not uM)
			State[n].CanSR     = 1e-3*Ca[n].NSR;      // Ca dependent currents, Cansr (in mM not uM)
			State[n].CajSR     = 1e-3*Ca[n].JSR;      // Ca dependent currents, Cajsr (in mM not uM)

			// Excitation state (necessary for SRF) | lib/Model.c 
			determine_excitation_state_integrated_0D(&Variables[n], Vm[n], sim_time, &Dyad[n].Ca_JSR_t_ex, Ca[n].JSR, &Dyad[n].SRF_prop_active,  SRF[n].SRF_prop_active, &SRF[n].waveform_init, &SRF[n].srf_set, SRF[n].Mode);
			Dyad[n].ex_switch  = Variables[n].ex_switch;

			// Spontaneous release functions || lib/Spontaneous_release_functions.cpp
			set_and_run_SRF(&SRF[n], &Dyad[n], SRF[n].Mode, &Rand[n], Variables[n].ex_switch, sim_time, Ca[n].JSR);
			calc_SRF_mults(&SRF[n], &MEM[n], &Dyad[n]);      

			// Spatial Ca handling ========================================================\\|
			// Zero reaction terms
			Ca[n].SS_reac = Ca[n].CYTO_reac = Ca[n].NSR_reac = Ca[n].JSR_reac = 0;

			// Inter-compartment transfer || lib/CRU.cpp
			comp_J_ds_ss(Params[n], Ca[n].DS, Ca[n].SS, Dyad[n].vol_ds, &Ca[n].SS_reac);
			comp_J_ss_cyto(Params[n], Ca[n].SS, Ca[n].CYTO, &Ca[n].SS_reac, &Ca[n].CYTO_reac);
			comp_J_nsr_jsr(Params[n], Ca[n].NSR, Ca[n].JSR, &Ca[n].NSR_reac, &Ca[n].JSR_reac);

			// Comp dyad || lib/CRU.cpp -> computes and solves JCaL, Jrel and Cads/jsr fluxes
			comp_dyad_0D(Params[n], &Dyad[n], Ca[n].DS, Ca[n].JSR, Ca[n].SS /*to which ds is coupled*/, &Ca[n].JSR_reac, Vm[n], Sim.dt, Params[n].Model);

			// Buffering || lib/CRU.cpp
			comp_buffering(Params[n], &Ca[n].Bcyto, &Ca[n].Bss, &Ca[n].Bjsr, Ca[n].CYTO, Ca[n].SS, Ca[n].JSR);

			// Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
			comp_SR_fluxes(Params[n], &SR[n], Ca[n].CYTO, Ca[n].NSR, &Ca[n].CYTO_reac, &Ca[n].NSR_reac);

			// Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
			comp_membrane_fluxes(Params[n], &MEM[n], State[n], Ca[n].CYTO, Ca[n].SS, &Ca[n].CYTO_reac, &Ca[n].SS_reac, Vm[n], MEM[n].NCX_SRF_mult);

			// trpn  || lib/myofilament.cpp || this is general needs to be looked at
			myofil[n].run_step_myofilament(1e-3*Ca[n].CYTO, 8, 0.015);
			Ca[n].CYTO_reac += -myofil[n].Jtrpn;

			// Update local concentrations
			Ca[n].DS        = (Ca[n].SS + Params[n].tau_ds*(Dyad[n].K_rel*Ca[n].JSR + Dyad[n].J_CaL))/(1 + Params[n].tau_ds*Dyad[n].K_rel); // quasi-steady-state approx
			Ca[n].SS        = Ca[n].SS      + Ca[n].Bss     *   Sim.dt*(Ca[n].SS_reac);
			Ca[n].CYTO      = Ca[n].CYTO    + Ca[n].Bcyto   *   Sim.dt*(Ca[n].CYTO_reac);
			Ca[n].NSR       = Ca[n].NSR     +                   Sim.dt*(Ca[n].NSR_reac);
			Ca[n].JSR       = Ca[n].JSR     + Ca[n].Bjsr    *   Sim.dt*(Ca[n].JSR_reac);

			// Whole-cell averages || including computing currents from Ca fluxes
			calc_whole_cell_values_including_currents_from_flux_0D(Params[n], Ca[n], &CRU[n], Dyad[n], SR[n], MEM[n], CRU[n].NTOT_CRUs);    // lib/CRU.cpp

			// Assign currents for use in AP model
			Variables[n].ICaL   = CRU[n].I_CAL;
			Variables[n].INCX   = CRU[n].I_NCX_bulk + CRU[n].I_NCX_ss;
			Variables[n].ICaP   = CRU[n].I_CaP_bulk + CRU[n].I_CaP_ss;
			Variables[n].ICab   = CRU[n].I_Cab_bulk + CRU[n].I_Cab_ss;
			// End Spatial Ca handling ====================================================//|

			// Solve the model || lib/Model.c -> lib/Model_X.cpp
			// This sets and updates all gates, and calculates Itot
			compute_model_integrated(Params[n], &Variables[n], &State[n], Vm[n], Sim.dt);   // lib/Model.c

			// Update local Voltage from Itot and stimulus current
			// Note [0].Istim is correct, as only calculated once; stim_area determines whether to actually apply stimulus to cell n
			State[n].Vm	= State[n].Vm + Sim.dt*(-(Variables[n].Itot + Variables[0].Istim*Tissue.stim_area[n] + Variables[0].Istim_S2*Tissue.S2_stim_area[n]));

			// Add multi_stim if set
			// If stim map is on, then now Istim[x] corresponds to stim_map = x, so region x will be stimulated when Istim[x] is non-zero
			if (strcmp(Tissue.Multi_stim, "On") == 0) for (int m = 1; m < Tissue.Nstims; m++) State[n].Vm += -(Sim.dt * Variables[m].Istim * Tissue.multi_stim_area[m][n]); 

			// Update local voltage due to spatial coupling
			State[n].Vm = State[n].Vm + Sim.dt*SC.diff[n];

			// Excitation state and measurements | lib/Model.c | "State.Vm" is voltage at t, "Vm" is voltage at t-dt
			calculate_measurement_properties(&Variables[n], Vm[n], State[n].Vm, sim_time, Sim.dt, -70, State[n].Cai, State[n].CanSR);		// -70 is APD V threshold	
		} 
		// End tissue loop - 1 ====================================//|

		// Loop over all tissue - 2 ===============================\\|
#pragma omp parallel for default(none) shared(SC, Vm, State, Cai, CaSR)
		for (int n = 0; n < SC.N; n++)
		{
			// Assign global voltage to state voltage
			// MUST be outside above loop
			Vm[n]			= State[n].Vm;

			// These are for spatial outputs
			Cai[n]			= 1e3*State[n].Cai; // uM
			CaSR[n]			= State[n].CanSR;
		}
		// End tissue loop - 2 ====================================//|

		// Output data to files - average and linescan ============\\|
		if (iteration_counter%(int)(1/Sim.dt) == 0) // if sim_time is an integer (i.e. per ms)
		{
			// Whole cell averages
			output_currents(out_cu, sim_time, Variables[cell1ref], State[cell1ref], Vm[cell1ref]);		// lib/Outputs.cpp || V, currents, gating variables, concs etc
			output_currents(out_cu2, sim_time, Variables[cell2ref], State[cell2ref], Vm[cell2ref]);		// lib/Outputs.cpp
			output_currents(out_cu3, sim_time, Variables[cell3ref], State[cell3ref], Vm[cell3ref]);		// lib/Outputs.cpp
			output_excitation_properties(out_ex, sim_time, Variables[cell1ref], Vm[cell1ref]);			// lib/Outputs.cpp || APD, excitation state, dv/dt etc
			output_excitation_properties(out_ex2, sim_time, Variables[cell2ref], Vm[cell2ref]);			// lib/Outputs.cpp	
			output_excitation_properties(out_ex3, sim_time, Variables[cell3ref], Vm[cell3ref]);			// lib/Outputs.cpp	
			output_CRU(out_cru1, sim_time, Ca[cell1ref], CRU[cell1ref], Vm[cell1ref]);                  // lib/Outputs.cpp || Ca concentrations, Jrel, JCaL, membrane and SR fluxes
			output_CRU(out_cru2, sim_time, Ca[cell2ref], CRU[cell2ref], Vm[cell2ref]);                  // lib/Outputs.cpp
			output_CRU(out_cru3, sim_time, Ca[cell3ref], CRU[cell3ref], Vm[cell3ref]);                  // lib/Outputs.cpp

			// Spatial data out ===============\\|
			// Linescan (idealised models only)
			if (strcmp(Tissue.Tissue_order, "geo") != 0) linescan_out_X(out_ls, SC, Vm, int(float(SC.NY/2)), int(float(SC.NZ/2)));// lib/Outputs.cpp

            // Full 3D spatial data (per unit output time)
            if (sim_time >= Sim.Spatial_output_start_time && sim_time <= Sim.Spatial_output_end_time)
            {
                if (Sim.Spatial_output_interval_vtk  > 0) // such that setting to zero means no spatial outputs
                {
                    if (outcount %Sim.Spatial_output_interval_vtk == 0) 
                    {
                        vtk_3D_output("Vm", directory, sr_dir, Vm, SC, outcount); 		// every x ms, output vtk file
                        vtk_3D_output("Ca", directory, sr_dir, Cai, SC, outcount); 		// every x ms, output vtk file
                        vtk_3D_output("CaSR", directory, sr_dir, CaSR, SC, outcount); 	// every x ms, output vtk file
                    }
                }
                if (Sim.Spatial_output_interval_data  > 0) // such that setting to zero means no spatial outputs
                {
                    if (outcount %Sim.Spatial_output_interval_data == 0) 
                    {
                        array_1D_output("Vm", directory, sr_dir, Vm, SC, outcount); 		// every x ms, output bin data array
                        array_1D_output("Ca", directory, sr_dir, Cai, SC, outcount); 		// every x ms, output bin data array
                        array_1D_output("CaSR", directory, sr_dir, CaSR, SC, outcount); 	// every x ms, output bin data array
                    }
                }
            }
            // End Spatial data out ===========//|

            // If phase output is set, and times are appropriate, output state to phase files || numbered 0-200
            if (strcmp(Sim.Write_state, "phase") == 0)	
            {
                if (sim_time > (Sim.NBeats-1)*Sim.BCL && sim_time < (Sim.NBeats -1)*Sim.BCL + 402)
                {
                    if (phase_counter%2 == 0) 
                    {
                        assign_state_variables_from_CRU_write(Dyad[5], Ca[5], &State[5]);
                        Write_state_phase(State[5], Params[5], Sim.BCL, PATH, Params[5].Model, 200-(phase_counter/2), Sim.state_reference_write); //lib/Read_write_state.c
                    }
                    printf("Written phase file %d\n", 200-(phase_counter/2));
                    phase_counter++;		
                }
            }

            outcount++;
        }
        // End Output data to files - average and linescan ========//|

        // Print SRF ti and NRyRopeak to file for every actually induced SCRE
        if (iteration_counter%(50*((int)(1/Sim.dt))) == 0) // as 50 is less than time between successive SRF, we only need to sample at 50 ms intervals
            for (int n = 0; n < SC.N; n++) print_SRF_properties_to_file(&SRF[n], out_srf_prop, n);

        iteration_counter ++;	// number of steps in dt
        if (iteration_counter%(500 *((int)(1/Sim.dt))) == 0) printf("Time = %.0fms\n",sim_time); // output every 500 ms
    }
    // End Time loop ============================================================================//|

    // Print final time in simulation land
    printf("Final Time = %.0fms\n\n",sim_time);

    // Write state
    if (strcmp(Sim.Write_state, "On") == 0) // whole tissue dump
    {
        for (int n = 0; n < SC.N; n++)
        {
            assign_state_variables_from_CRU_write(Dyad[n], Ca[n], &State[n]);	// lib/CRU.cpp
        }
        Write_state_tissue_integrated_whole_tissue(State, Params, Sim.BCL, PATH, Params_global.Model, SC.N, Tissue.Tissue_order, Tissue.Tissue_model, Tissue.Tissue_type, Tissue.Orientation_type, Sim.state_reference_write); //lib/Read_write_state.c
        printf("State written to file\n");
    }
    else if (strcmp(Sim.Write_state, "ave") == 0) // writes state for just one cell in the tissue (for region x)
    {
        if (strcmp(Tissue.Tissue_type, "homogeneous") != 0)
        {
            printf("ERROR: average tissue state write must be performed on homogeneous tissue - if wanting to apply to heterogeneous, run 1D homogeneous model for each celltype\n");
            exit(1);
        }
        else 
        {
            assign_state_variables_from_CRU_write(Dyad[10], Ca[10], &State[10]); // lib/CRU.cpp
            Write_state_tissue_integrated_ave_tissue(State[10], Params[10], Sim.BCL, PATH, Params[10].Model, Sim.state_reference_write); //lib/Read_write_state.c
        }
        printf("State written to file - one coupled cell\n");
    }
    // End Write state

    // Output final beat properties to file and screen  || APD, dvdt_max etc
    char * log_reference    = (char*)malloc(500);
    sprintf(log_reference, "%s/Properties_log.dat", directory);
    output_properties_to_screen(log_reference, Variables[cell2ref], Sim);     // lib/Outputs.cpp
    free(log_reference);

    // Output ativation map, final beat, vtk and datafile || lib/Outputs.cpp
    Output_activation(directory, sr_dir, Variables, SC);

    // Calculate and output conduction velocity || lib/Tissue.cpp
    if (strcmp(Tissue.Tissue_model, "conduction_velocity") == 0) calculate_CV(Tissue, Variables, directory);

    // Calculate and output conduction success for VW (1D only - any tissue model - only makes sense for S1-S2 pacing, and so only calculates after S2)
    if (strcmp(Tissue.Tissue_order, "1D") == 0 && Sim.S2_CL != 0) compute_conduction_success(Tissue, Variables , SC.N, Sim.S2_time, Sim.S2_CL, directory); // lib/Tissue.cpp

    time (&rawtime);
    printf("|============================================================|\n");
    printf("|Code has now finished. Finished at %s", ctime (&rawtime));
    printf("|============================================================|\n");

	// Output disclaimer to screen
	printf("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
	output_disclaimer_citations(Params_global, Sim);   // lib/Outputs.c
	output_disclaimer_citations_tissue(Params_global, Tissue);
	output_disclaimer_citations_spatial_cell(Params_global, Sim); 
	printf("--\n/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n");

	// Free memory
	free(directory);
	free(results_dir);
    free(res_dir_full);
	free(sr_dir);
	SC_array_deallocation(&SC);			// lib/Spatial_coupling.cpp
	tissue_array_deallocation(&Tissue);	// lib/Tissue.cpp
	delete [] Params;
	delete [] State;
	delete [] Variables;
	delete [] Ca;
	delete [] CRU;
	delete [] Dyad;
	delete [] SR;
	delete [] MEM;
	delete [] Rand;
	delete [] SRF;
	delete [] Vm;
	delete [] Cai;
	delete [] CaSR;
	delete [] myofil;
} 
// End Main *************************************************************************************//|

