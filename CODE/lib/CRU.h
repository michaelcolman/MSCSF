// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Spatial and 0D intracellular Ca2+ handling ==  //
// settings and dynamics, header.==========================  //
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

#ifndef CRU_H
#define CRU_H

#include <math.h>
#include "Structs.h"
#include "MersenneTwister.h"

// Whole CRU functions ============================================\\|
// Array allocation and deallocation
void Ca_array_allocation(int NCRU, Ca_variables *Ca);
void Ca_array_deallocation(Ca_variables *Ca);
void CRU_map_array_allocation(int NCRU, CRU_variables *cru);
void CRU_map_array_deallocation(CRU_variables *cru);
void Dyad_array_allocation(Dyad_variables *d);
void Dyad_array_deallocation(Dyad_variables *d);

// 3D cell settings
void spatial_cell_settings(CRU_variables *cru, Argument_parameters A);

// tau ss type
void set_tau_ss(Cell_parameters *p);

// Local scaling including TT map scales
void initialise_sub_cellular_het_maps(int N, CRU_variables *cru);
void set_sub_cellular_local_scale(Cell_parameters p, Dyad_variables *d, Membrane_fluxes *m, SR_fluxes *sr);
void read_sub_cellular_het_maps(SC_variables sc, CRU_variables *cru, const char *PATH, const char *Outputs_dir);
void set_sub_cellular_local_het_scale(Cell_parameters p, int N, Dyad_variables *d, Membrane_fluxes *m, SR_fluxes *sr, CRU_variables cru);

// Set dyad heterogeneity
void set_dyad_heterogeneity(Cell_parameters p, Dyad_variables *d, double rand, int n, char const *directory, double *map, CRU_variables cru);
void write_random_dyad_het_vtk(SC_variables sc, double *map, const char* Output_dir);

// Initial conditions
void initial_conditions_calcium(Ca_variables *Ca, Cell_parameters p, int NCRU);
void initial_conditions_calcium_0D(Ca_variables *Ca, Cell_parameters p);
void initial_conditions_dyad_stochastic(Dyad_variables *d);
void initial_conditions_dyad_det(Dyad_variables *d);
void assign_CRU_variables_from_state_read(Dyad_variables *d, Ca_variables *Ca, State_variables s);
void assign_state_variables_from_CRU_write(Dyad_variables d, Ca_variables Ca, State_variables *s);

// Inter-compartment transfer
void comp_J_ds_ss(Cell_parameters p, double Ca_ds, double Ca_ss, double vol_ds, double *reac_ss);
void comp_J_ss_cyto(Cell_parameters p, double Ca_ss, double Ca_cyto, double *reac_ss, double *reac_cyto);
void comp_J_nsr_jsr(Cell_parameters p, double Ca_nsr, double Ca_jsr, double *reac_nsr, double *reac_jsr);

// Buffering
void comp_buffering(Cell_parameters p, double *Bcyto, double *Bss, double *Bjsr, double Ca_cyto, double Ca_ss, double Ca_jsr);
void buffering_cyto(Cell_parameters p, double *Bcyto, double Ca);
void buffering_subspace(Cell_parameters p, double *Bss, double Ca);
void buffering_JSR(Cell_parameters p, double *Bjsr, double Ca);

// Whole cell averages and currents
void calc_whole_cell_values_including_currents_from_flux(int N, Cell_parameters p, Ca_variables *Ca, CRU_variables *cru, Dyad_variables *d, SR_fluxes *sr, Membrane_fluxes *m, int NTOT);
void calc_whole_cell_values_including_currents_from_flux_0D(Cell_parameters p, Ca_variables Ca, CRU_variables *cru, Dyad_variables d, SR_fluxes sr, Membrane_fluxes m, int NTOT);

// Current from flux
double compute_current_from_flux(Cell_parameters p, double flux, int valence, double vol, int NCRUs);
// End whole CRU functions ========================================//|

// Dyad fluxes functions ==========================================\\|
void comp_dyad_3D(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Ca_jsr, double Ca_cyto /*to which ds is coupled*/, double *reac_jsr, double Vm, double dt, const char *Model);
void comp_dyad_0D(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Ca_jsr, double Ca_cyto /*to which ds is coupled*/, double *reac_jsr, double Vm, double dt, const char *Model);

// RyR
void set_and_update_monomer_state(Cell_parameters p, Dyad_variables *d, double Ca_jsr, double dt);
void set_RyR_rates(Cell_parameters p, Dyad_variables *d, double Ca_ds);
void update_RyR_stochastic(Dyad_variables *d, double dt);

// LTCC
void set_LTCC_rates(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Vm, const char * Model);
void comp_LTCC_bar(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Vm);
void update_gates_LTCC_det(Cell_parameters p, Dyad_variables *d, double dt);
void update_LTCC_stochastic(Dyad_variables *d, double dt);
// End dyad fluxes functions ======================================//|

// SR fluxes functions ============================================\\|
void comp_SR_fluxes(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr, double *reac_cyto, double *reac_nsr);
void comp_Jup(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr);
void comp_Jleak(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr);
// End SR fluxes functions ========================================//|

// MEM fluxes functions ===========================================\\|
void comp_membrane_fluxes(Cell_parameters p, Membrane_fluxes *m, State_variables s, double Ca_cyto, double Ca_ss, double *reac_cyto, double *reac_ss, double Vm, double SRF_mult);
void comp_JMEM(Cell_parameters p, Membrane_fluxes *m, State_variables s, double Ca_cyto, double Ca_ss, double Vm, double SRF_mult);
double comp_JNCX(Cell_parameters p, Membrane_fluxes *m, double Cai, State_variables s, double Vm);
double comp_JCaP(Cell_parameters p, double Cai);
double comp_JCab(Cell_parameters p, double Cai, double Vm);
// End MEM fluxes functions =======================================//|

// Voltage clamp 3D and 0D cell
void run_voltage_clamp_3Dcell(Cell_parameters p, Model_variables *var, State_variables *s, Dyad_variables *d, Ca_variables *Ca, CRU_variables *cru, SC_variables *sc, SR_fluxes *sr, Membrane_fluxes *m, RAND *rand, char const *directory, double dt);
void run_voltage_clamp_0Dcell(Cell_parameters p, Model_variables *var, State_variables *s, Dyad_variables *d, Ca_variables *Ca, CRU_variables *cru, SR_fluxes *sr, Membrane_fluxes *m, char const *directory, double dt);

#endif

