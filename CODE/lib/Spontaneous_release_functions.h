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

#ifndef SRF_H
#define SRF_H

#include "Structs.h"
#include <math.h>

// Defaults and setup
void set_SRF_defaults(Spontaneous_release_functions *srf, Argument_parameters A);
void SRF_tissue_heterogeneity(Spontaneous_release_functions *srf, double rand);
void SRF_setup(Spontaneous_release_functions *srf, Argument_parameters A);

// Parameter sets - Defined
void set_SRF_Direct_Control_parameters(Spontaneous_release_functions *srf, Argument_parameters A);

// Parameter sets - Dynamic
void set_SRF_dynamic_parameters_3D_cell(Spontaneous_release_functions *srf);
void set_SRF_dynamic_parameters_controllable(Spontaneous_release_functions *srf, Argument_parameters A);

// Distribution and probablity parameters - Dynamic
void determine_SRF_3D_cell_PSCR(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_3D_cell_ti_sep(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_3D_cell_ti_dist_widths(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_3D_cell_MD(Spontaneous_release_functions *srf, double CaSR);

void determine_SRF_dynamic_general_PSCR(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_dynamic_general_ti_sep(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_dynamic_general_ti_dist_widths(Spontaneous_release_functions *srf, double CaSR);
void determine_dynamic_general_cell_MD(Spontaneous_release_functions *srf, double CaSR);
void determine_SRF_dynamic_general_duration_dist_widths(Spontaneous_release_functions *srf, double CaSR);

void determine_SRF_params_from_CaSR(Spontaneous_release_functions *srf, double CaSR);
void need_to_recalculate(Spontaneous_release_functions *srf, double CaSR, double threshold, double NRyR);

// Set and run
void set_and_run_SRF(Spontaneous_release_functions *srf, Dyad_variables *d, const char * Mode, RAND *rand, int ex_switch, double sim_time, double CaJSR);
void run_SRF(Spontaneous_release_functions *srf, Dyad_variables *d, double sim_time);
void print_SRF_properties_to_file(Spontaneous_release_functions *srf, std::ostream& out, int n);
void calc_SRF_mults(Spontaneous_release_functions *srf, Membrane_fluxes *mem, Dyad_variables *dyad);

// Read
void read_SRF_settings_from_file(Spontaneous_release_functions *srf, const char *input_file);

// Waveform
void Determine_waveform_parameters(Spontaneous_release_functions *srf);
void Determine_waveform_parameters_long(Spontaneous_release_functions *srf);
void Waveform(Spontaneous_release_functions *srf, double t);
void Waveform_plateau(Spontaneous_release_functions *srf, double t);

// Distribution functions
void Determine_ti(Spontaneous_release_functions *srf, double rand);

void Determine_duration(Spontaneous_release_functions *srf, double rand);
void Determine_duration_ks_from_MD(Spontaneous_release_functions *srf);

void Determine_NRyRo_peak(Spontaneous_release_functions *srf, double duration, double rand);
void Determine_NRyRo_plateau(Spontaneous_release_functions *srf, double duration, double rand);

// Test and produce
void test_and_produce_distributions(Spontaneous_release_functions *srf, const char * directory, RAND *r, double CaSR);
void test_and_produce_CaSR_dependency(Spontaneous_release_functions *srf, const char * directory, RAND *r);

#endif

