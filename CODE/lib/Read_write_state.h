// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Read and write state. Outputs state variables  //
// for single cell and tissue in mulitple ways ============  //
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

#ifndef RW_STATE_H
#define RW_STATE_H

#include "Structs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void Write_state_variables_native(State_variables s, FILE *out, const char * Model);
void Read_state_variables_native(State_variables *s, FILE *in, const char * Model);
void Read_spatial_Ca_system_state(Dyad_variables *d, Ca_variables *Ca, int n, FILE *in);

void Write_state_single_cell_native(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_single_cell_native(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Write_state_single_cell_integrated(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_single_cell_integrated(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Write_state_single_cell_integrated_0D(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_single_cell_integrated_0D(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);

void Write_state_single_cell_integrated_spatial(State_variables s, Cell_parameters p, Dyad_variables *d, Ca_variables Ca, int BCL, const char * PATH, const char *Model, const char * State_ref, int NX, int NY, int NZ, int N);
void Read_state_single_cell_integrated_spatial(State_variables *s, Cell_parameters p, Dyad_variables *d, Ca_variables *Ca, int BCL, const char * PATH, const char *Model, const char * State_ref, int NX, int NY, int NZ, int N);

void Write_state_tissue_native_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Read_state_tissue_native_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Write_state_tissue_integrated_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Read_state_tissue_integrated_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Write_state_tissue_native_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_tissue_native_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Write_state_tissue_integrated_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_tissue_integrated_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);

void Write_state_tissue_native_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Read_state_tissue_native_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Write_state_tissue_integrated_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Read_state_tissue_integrated_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char* Orientation_type, const char * State_ref);
void Write_state_tissue_native_net_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_tissue_native_net_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Write_state_tissue_integrated_net_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);
void Read_state_tissue_integrated_net_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref);

void Write_state_phase(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref);
void Read_state_phase(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref);
void Write_state_phase_0D(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref);
void Read_state_phase_0D(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref);

#endif
