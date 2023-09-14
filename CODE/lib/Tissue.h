// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Tissue model setttings and options ==========  //
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

#ifndef TISSUE_H
#define TISSUE_H

#include "Structs.h"
#include <math.h>

// Set model and order
void set_tissue_model_conditions(Tissue_parameters *t, Argument_parameters A);

// Set tissue setttings
void set_tissue_settings_idealised(Cell_parameters p, Tissue_parameters *t);
void set_tissue_settings_anatomical(Cell_parameters p, Tissue_parameters *t);
void set_coord_stim_and_map_from_defined_type(Cell_parameters p, Tissue_parameters *t, Argument_parameters A);
void overwrite_tissue_properties_from_args(Cell_parameters p, Tissue_parameters *t, Argument_parameters A);
void set_global_orientation_direction_from_arg(Cell_parameters p, Tissue_parameters *t, Argument_parameters A);

// Array allocation and deallocation 
void tissue_array_allocation(Tissue_parameters *t, int Ncell);
void tissue_array_deallocation(Tissue_parameters *t);

// Create or read geometries
void select_tissue_geometry_function(Tissue_parameters t, SC_variables *sc, const char *PATH, const char* Output_dir);
void create_idealised_geometry_homogeneous(SC_variables *sc);
void create_idealised_geometry_heterogeneous(SC_variables *sc, Tissue_parameters t);

// Create or read stimulus area
void select_stimulus_area_function(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, int S2_CL);
void select_stimulus_area_function_multi_stim(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, int S2_CL, int Nstims);
void create_stimulus_area(SC_variables sc, const char * Tissue_order, int * stim_area, int x, int xs, int y, int yz, int z, int zs, int *Nstim, const char* Output_dir, const char *ref);
void create_stimulus_area_sphere(SC_variables sc, const char *Tissue_order, int * stim_area, int x, int xs, int y, int z, int *Nstim, const char* Output_dir, const char *ref);

// Create or read orientation
void set_orientation(SC_variables *sc, Tissue_parameters t, const char *PATH, const char *Tissue_order);
void create_orientation_ideal(SC_variables *sc, Tissue_parameters t);
void read_orientation_anatomical(SC_variables *sc, Tissue_parameters t, const char *PATH);
void read_orientation_anatomical_transverse(SC_variables *sc, Tissue_parameters t, const char *PATH);
void output_fibre_orientation(SC_variables sc, Tissue_parameters t, const char* Output_dir);
void output_fibre_orientation_o2(SC_variables sc, Tissue_parameters t, const char* Output_dir);
void output_fibre_orientation_o3(SC_variables sc, Tissue_parameters t, const char* Output_dir);

// General maps
void create_or_read_map_double(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, double *map, const char *map_file, const char *ref);
void create_map_patch(SC_variables sc, const char *Tissue_order, double * map_patch, int x, int xs, int y, int ys, int z, int zs, const char* Output_dir, const char *ref, const char * shape);

// Non-uniform, celltype Diffusion coefficients
void update_D_arrays_Dscale_baseline(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory);
void update_D_arrays_Dscale_mod(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory);

// Create phase map for phase re-entry
void create_read_phasemap(int * phase, Tissue_parameters t, SC_variables sc, const char* PATH, const char* Output_dir);
void create_phasemap_2D(int * phase, SC_variables sc, const char* Output_dir);
void create_phasemap_3D(int * phase, SC_variables sc, const char* Output_dir);

// CV calculation
void set_CV_cells(Cell_parameters p, Tissue_parameters *t);
void set_CV_cells_const_distance(Cell_parameters p, Tissue_parameters *t);
void calculate_CV(Tissue_parameters t, Model_variables *var, const char* directory);

// Conduction success calculation
void compute_conduction_success(Tissue_parameters t, Model_variables *var, int N, double S2_time, double S2_CL, const char* directory);

// Disconnect regions
void Modify_neighbours_region_disconnect(SC_variables *sc, Tissue_parameters *t);

//NETWORK model
void update_G_arrays_Dscale_baseline(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory);
void update_G_arrays_Dscale_mod(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory);

#endif

