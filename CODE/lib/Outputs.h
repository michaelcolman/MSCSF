// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Outputs. For outputting settings, data and ==  //
// disclaimer references to file and screen ===============  //
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

#ifndef OUTPUTS_H
#define OUTPUTS_H

#include "Structs.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>


// Global output functions
void output_properties_to_screen(const char * log_reference, Model_variables var, Simulation_parameters Sim);
void output_currents(std::ostream& out, double sim_time, Model_variables var, State_variables s, double Vm);
void output_excitation_properties(std::ostream& out, double sim_time, Model_variables var, double Vm);

void output_currents_csv(std::ostream& out, double sim_time, Model_variables var, State_variables s, double Vm);
void output_excitation_properties_csv(std::ostream& out, double sim_time, Model_variables var, double Vm);

// Integrated Ca handling models only
void output_CRU(std::ostream& out, double sim_time, Ca_variables Ca, CRU_variables cru, double Vm);
void output_CRU_csv(std::ostream& out, double sim_time, Ca_variables Ca, CRU_variables cru, double Vm);

// Spatial outputs (3D cell and tissue models)
void linescan_out_X(std::ostream& out, SC_variables sc, double * variable, int y, int z);
void linescan_out_Y(std::ostream& out, SC_variables sc, double * variable, int x, int z);
void linescan_out_Z(std::ostream& out, SC_variables sc, double * variable, int x, int y);
void data_2D_XYslice_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count, int z);
void data_2D_XZslice_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count, int y);
void data_2D_YZslice_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count, int x);
void vtk_3D_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count);
void vtk_3D_output_region_specific(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count, int region);
void array_1D_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count);
void Output_activation(const char * dir, const char * dir2, Model_variables *v, SC_variables sc);
void data_3D_output(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count);
void array_1D_binary_read(const char *string,  const char * dir, const char * dir2, double *variable, SC_variables sc, int count);

// Settings
void output_settings(Simulation_parameters sim, char const * directory, bool DC_current_mod_arg, Cell_parameters p, int argc, char *argin[]);
void output_settings_tissue(Simulation_parameters sim, Tissue_parameters t, char const * directory);
void output_settings_3D_cell(Cell_parameters p, Simulation_parameters sim, CRU_variables cru, char const * directory);
void output_settings_0D_cell(Cell_parameters p, Simulation_parameters sim, CRU_variables cru, char const * directory, Spontaneous_release_functions srf);

// Disclaimer and citations
void output_disclaimer_citations(Cell_parameters p, Simulation_parameters sim);
void output_disclaimer_citations_tissue(Cell_parameters p, Tissue_parameters t);
void output_disclaimer_citations_tissue_network(Cell_parameters p, Tissue_parameters t);
void output_disclaimer_citations_spatial_cell(Cell_parameters p, Simulation_parameters sim);

#endif


