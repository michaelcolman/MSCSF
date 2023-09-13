// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Spatial coupling. Includes functions for: ===  //
//      array allocation ==================================  //
//      reading geometries and maps =======================  //
//      setup of diffusion coefficient arrays =============  //
//      neighbour and array type mapps ====================  //
//      implementation of the anisotropic FDM method ======  //
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

#ifndef SC_H
#define SC_H

#include "Structs.h"
#include <math.h>


// Array sizes and allocation/deallocation
void SC_set_array_sizes(SC_variables *sc, int NX, int NY, int NZ);
void SC_array_allocation_N3(SC_variables *sc, int NX, int NY, int NZ);
void SC_array_allocation_Ncell(SC_variables *sc, int Ncell);
void SC_array_deallocation(SC_variables *sc);

// Read files | returns Ncells or NMap
int read_geo_file(SC_variables *sc, int *geo, const char * filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref, int Ncelltypes);
int read_map_file(SC_variables sc, int *map, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref);
int read_map_file_double(SC_variables sc, double *map, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref);

// D and dx
void set_D_dx_global(SC_variables *sc, double dx, double dy, double dz, double D1, double D_AR);
void set_D_array_anisotropic(SC_variables *sc);
void output_D1_and_D2(SC_variables sc, const char* Output_dir);

// Cell index and neighbouhood arrays
void SC_set_index_and_geo_linear(SC_variables *sc);
void SC_set_neighbours(SC_variables *sc);

// Finite difference method functions
// uniform, isotropic
void calc_diff_FDM_iso(SC_variables *sc, double *v, int n);
void calc_diff_FDM_tau(SC_variables *sc, double *v, int n, double tau_trans, double tau_long);

// anisotropic/non-unfirm
void calc_dD_anisotropic_3D(SC_variables *sc, int n);
void calc_diff_FDM_anisotropic(SC_variables *sc, double *v, int n);

// BCs and solver
void calc_laplacian_and_BCs(SC_variables *sc, int n);
void calc_diff_from_lap(SC_variables *sc, double *v, int n);

// All network model functions
void zero_orientation_ideal(SC_variables *sc);
void SC_array_allocation_Njunc(SC_variables *sc, int N);
void SC_array_deallocation_Njunc(SC_variables *sc);
void set_G_dx(SC_variables *sc, double dx, double dy, double dz, double Gt, double Gl, double Gscale); // now unused?
//void set_G_dx_global(SC_variables *sc, double dx, double dy, double dz, double Gl, double Gt, double Gt2);
void set_G_dx_global(SC_variables *sc, double dx, double dy, double dz, double Gl, double Gt, double Gt2, double symm_fac_diag_long, double symm_fac_diag_trans, double symm_fac_corner);
void set_gGgap_array(SC_variables *sc, char const* Orientation_type);
void update_junctions(SC_variables *sc);
void calc_N_junctions(SC_variables *sc);
void set_gjunc_and_junc_maps(SC_variables *sc, const char* Output_dir, Tissue_parameters *t);
void calc_IGap(SC_variables *sc, double *v, int n);
void output_junc_maps(SC_variables *sc, const char* Output_dir, Tissue_parameters *t);
void default_junction_maps(SC_variables *sc);
void read_map_file_Njunc(SC_variables *sc, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, double *map);

#endif
