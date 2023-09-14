// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Spatial and 0D intracellular Ca2+ handling ==  //
// settings and dynamics. =================================  //
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

#include "CRU.h"
#include "Structs.h"
#include "Model.h"
#include "Spatial_coupling.h"
#include <omp.h>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>

// Function list  =====================================================================\\|
//	Array allocation and deallocation
//     Ca_array_allocation()
//     Ca_array_deallocation()
//     CRU_map_array_allocation()
//     CRU_map_array_deallocation()
//     Dyad_array_allocation()
//     Dyad_array_deallocation()
//	
//	spatial_cell_settings()
//	set_tau_ss()
//	set_sub_cellular_local_scale()
//  read_sub_cellular_het_maps()
//  set_sub_cellular_local_het_scale()
//	set_dyad_heterogeneity()
//	write_random_dyad_het_vtk()
//	
//	initial conditions
//	    initial_conditions_calcium()
//	    initial_conditions_calcium_0D()
//	    initial_conditions_dyad_stochastic()
//	    initial_conditions_dyad_det()
//	    assign_CRU_variables_from_state_read()
//      assign_state_variables_from_CRU_write()
//	
//	Inter-compartment transfer functions
//	    comp_J_ds_ss()
//	    comp_J_ss_cyto()
//	    comp_J_nsr_jsr()
//	
//	Buffering
//	    comp_buffering()
//	    buffering_cyto()
//	    buffering_subspace()
//	    buffering_JSR()
//	
//	calc_whole_cell_values_including_currents_from_flux()
//	calc_whole_cell_values_including_currents_from_flux_0D()
//	compute_current_from_flux()
//	
//	Fluxes functions
//	    comp_dyad_3D()
//	    comp_dyad_0D()
//	    set_and_update_monomer_state()
//	    set_RyR_rates()
//	    set_LTCC_rates()
//	    comp_LTCC_bar()
//	    update_gates_LTCC_det()
//	    comp_SR_fluxes()
//	    comp_Jup()
//	    comp_Jleak()
//		comp_membrane_fluxes()
//	    comp_JMEM()
//	    comp_JNCX()
//	    comp_JCab()
//	    comp_JCaP()
//	
//	Stochastic integration / state update
//	    update_RyR_stochastic()
//	    update_LTCC_stochastic()
//	
//	Voltagae clamp
//	    run_voltage_clamp_3Dcell()
//	    run_voltage_clamp_0Dcell()
// End Function list ==================================================================//|

// Whole CRU / global functions =================================================================\\|
// Array allocation and deallocation ==================================================\\|
void Ca_array_allocation(int NCRU, Ca_variables *Ca)
{
	Ca->ds          = new double [NCRU];
	Ca->ss          = new double [NCRU];
	Ca->cyto        = new double [NCRU];
	Ca->nsr         = new double [NCRU];
	Ca->jsr         = new double [NCRU];

	Ca->ss_reac     = new double [NCRU];
	Ca->cyto_reac   = new double [NCRU];
	Ca->nsr_reac    = new double [NCRU];
	Ca->jsr_reac    = new double [NCRU];

	Ca->bcyto       = new double [NCRU];
	Ca->bss         = new double [NCRU];
	Ca->bjsr        = new double [NCRU];
}

void Ca_array_deallocation(Ca_variables *Ca)
{
	delete  [] Ca->ds;
	delete  [] Ca->ss;
	delete  [] Ca->cyto;
	delete  [] Ca->jsr;
	delete  [] Ca->nsr;

	delete  [] Ca->ss_reac;
	delete  [] Ca->cyto_reac;
	delete  [] Ca->nsr_reac;
	delete  [] Ca->jsr_reac;

	delete  [] Ca->bcyto;
	delete  [] Ca->bss;
	delete  [] Ca->bjsr;
}

void CRU_map_array_allocation(int NCRU, CRU_variables *cru)
{
    cru->TT_map         = new int    [NCRU];
    cru->SERCA_map      = new double [NCRU];
    cru->NCX_map        = new double [NCRU];
    cru->RyR_het_map    = new double [NCRU];
    cru->LTCC_map       = new double [NCRU];
    cru->dyad_het_map  = new double [NCRU];
}

void CRU_map_array_deallocation(CRU_variables *cru)
{
    delete [] cru->TT_map;
    delete [] cru->SERCA_map;
    delete [] cru->NCX_map;
    delete [] cru->RyR_het_map;
    delete [] cru->LTCC_map;
    delete [] cru->dyad_het_map;
}

// Arrays with multiple elements per dyad
void Dyad_array_allocation(Dyad_variables *d)
{
	d->RyR_state		= new int 		[d->NRyR]; // this is local NRyR as per dyad
	d->rand_RyR			= new double 	[d->NRyR];

	d->LTCC_va_state	= new int 		[d->NLTCC];
	d->LTCC_vi_state	= new int 		[d->NLTCC];
	d->LTCC_ci_state	= new int 		[d->NLTCC];
	d->rand_LTCC		= new double	[d->NLTCC];
}

void Dyad_array_deallocation(Dyad_variables *d)
{
    delete [] d->RyR_state;
    delete [] d->rand_RyR;
    delete [] d->LTCC_va_state;
    delete [] d->LTCC_vi_state;
    delete [] d->LTCC_ci_state;
    delete [] d->rand_LTCC;
}
// End Array allocation and deallocation ==============================================//|

// 3D cell settings ===================================================================\\|
void spatial_cell_settings(CRU_variables *cru, Argument_parameters A)
{
	// Defaults
	cru->Cell_size	 		= "standard"; 	// currently only "standard" and "thin" as options
	cru->Sim_Cell_size	 	= "full";	    // "testing", "portion" or "full" (or other option if added by you!)
	cru->Detub				= "Off";
    cru->LTCC_redist        = "On";
    cru->SERCA_het          = "Off";
    cru->NCX_het            = "Off";
    cru->RyR_het            = "Off";
    cru->LTCC_het           = "Off";
    cru->volds_het          = "Off";
    cru->TT_map_file        = "none";
    cru->SERCA_map_file     = "none";
    cru->NCX_map_file       = "none";
    cru->RyR_het_map_file   = "none";
    cru->LTCC_map_file      = "none";

	// Arguments
	if (A.Cell_size_arg == true)     	cru->Cell_size     	    = A.Cell_size;
	if (A.Sim_cell_size_arg == true)    cru->Sim_Cell_size 	    = A.Sim_cell_size;
	if (A.Detub_arg == true)			cru->Detub			    = A.Detub;
	if (A.LTCC_redist_arg == true)		cru->LTCC_redist	    = A.LTCC_redist;
	if (A.TT_map_file_arg == true)		cru->TT_map_file	    = A.TT_map_file;
	if (A.SERCA_het_arg == true)		cru->SERCA_het		    = A.SERCA_het;
	if (A.SERCA_map_file_arg == true)	cru->SERCA_map_file	    = A.SERCA_map_file;
	if (A.NCX_het_arg == true)		    cru->NCX_het		    = A.NCX_het;
	if (A.NCX_map_file_arg == true)	    cru->NCX_map_file	    = A.NCX_map_file;
	if (A.RyR_het_arg == true)		    cru->RyR_het		    = A.RyR_het;
	if (A.RyR_het_map_file_arg == true)	cru->RyR_het_map_file	= A.RyR_het_map_file;
	if (A.LTCC_het_arg == true)		    cru->LTCC_het		    = A.LTCC_het;
	if (A.LTCC_map_file_arg == true)	cru->LTCC_map_file	    = A.LTCC_map_file;
	if (A.volds_het_arg == true)		cru->volds_het		    = A.volds_het;

	if (strcmp(cru->Detub, "On") == 0 && strcmp(cru->Sim_Cell_size, "full") != 0)
	{
		printf("ERROR: Simulation cell size must be \"full\" to use a TT map. Currently set to %s\n", cru->Sim_Cell_size);
		exit(1);	
	}
    if (strcmp(cru->SERCA_het, "On") == 0 && strcmp(cru->Sim_Cell_size, "full") != 0)
    {
        printf("ERROR: Simulation cell size must be \"full\" to use a SERCA map. Currently set to %s\n", cru->Sim_Cell_size);
        exit(1);
    }
    if (strcmp(cru->NCX_het, "On") == 0 && strcmp(cru->Sim_Cell_size, "full") != 0)
    {
        printf("ERROR: Simulation cell size must be \"full\" to use a NCX map. Currently set to %s\n", cru->Sim_Cell_size);
        exit(1);
    }
    if (strcmp(cru->LTCC_het, "map") == 0 && strcmp(cru->Sim_Cell_size, "full") != 0)
    {
        printf("ERROR: Simulation cell size must be \"full\" to use a LTCC map. Currently set to %s\n", cru->Sim_Cell_size);
        exit(1);
    }
    if (strcmp(cru->RyR_het, "map") == 0 && strcmp(cru->Sim_Cell_size, "full") != 0)
    {
        printf("ERROR: Simulation cell size must be \"full\" to use a RyR map. Currently set to %s\n", cru->Sim_Cell_size);
        exit(1);
    }

	// Actual cell dimensions based on cell size
	if (strcmp(cru->Cell_size, "standard") == 0)
	{
		cru->NX2        = 15;
		cru->NY2        = 20;
		cru->NZ2        = 65; // longitudinal   
		cru->NTOT_CRUs	= cru->NX2*cru->NY2*cru->NZ2;
	}
	else if (strcmp(cru->Cell_size, "thin") == 0)
	{
		cru->NX2        = 10;
		cru->NY2        = 15;
		cru->NZ2        = 95; // longitudinal   
		cru->NTOT_CRUs	= cru->NX2*cru->NY2*cru->NZ2;
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid cell size reference. Please check CRU.cpp for options\n\n", cru->Cell_size);
		exit(1);
	}

	// Cm scaling performed after this function in Main

	// Set cell sizes for simulation (if using smaller portion)
	if (strcmp(cru->Sim_Cell_size, "full") == 0)
	{	
		cru->NX         = cru->NX2;
		cru->NY         = cru->NY2;
		cru->NZ         = cru->NZ2; // longitudinal
	}
	else if (strcmp(cru->Sim_Cell_size, "portion") == 0) // cross-section; longitudinal portion
	{
		cru->NX         = cru->NX2;
		cru->NY         = cru->NY2;
		cru->NZ         = 20; // longitudinal
	}
	else if (strcmp(cru->Sim_Cell_size, "testing") == 0)
	{
		cru->NX         = 5;
		cru->NY         = 5;
		cru->NZ         = 10; // longitudinal   
	}
	else if (strcmp(cru->Sim_Cell_size, "debug") == 0)
	{
		cru->NX         = 1;
		cru->NY         = 1;
		cru->NZ         = 2; // longitudinal   
	}
	else
	{
		printf("ERROR: \"%s\" is not a valid sim cell size reference. Please check CRU.cpp for options\n\n", cru->Sim_Cell_size);
		exit(1);
	}

}
// End 3D cell settings ===============================================================//|

// tau SS type assignment =============================================================\\|
// Simpy sets the actual parameters for tau ss from the different options in Parameters
void set_tau_ss(Cell_parameters *p)
{
	if (strcmp(p->tau_ss_type, "slow") == 0) 
	{
		p->tau_ss_trans = p->tau_ss_trans_s;
		p->tau_ss_long 	= p->tau_ss_long_s;
	}
	else if (strcmp(p->tau_ss_type, "medium_slow") == 0) 
	{
		p->tau_ss_trans = p->tau_ss_trans_ms;
		p->tau_ss_long 	= p->tau_ss_long_ms;
	}
	else if (strcmp(p->tau_ss_type, "medium") == 0) 
	{
		p->tau_ss_trans = p->tau_ss_trans_m;
		p->tau_ss_long 	= p->tau_ss_long_m;
	}
	else if (strcmp(p->tau_ss_type, "medium_fast") == 0) 
	{
		p->tau_ss_trans = p->tau_ss_trans_mf;
		p->tau_ss_long 	= p->tau_ss_long_mf;
	}
	else if (strcmp(p->tau_ss_type, "fast") == 0) 
	{
		p->tau_ss_trans = p->tau_ss_trans_f;
		p->tau_ss_long 	= p->tau_ss_long_f;
	}
	else
	{
		printf("ERROR: tau_ss_type \"%s\" is not valid\n", p->tau_ss_type);
		exit(1);
	}
}
// End tau SS type assignment =========================================================//|

// Sub-cellular local scale factors ===================================================\\|
void initialise_sub_cellular_het_maps(int N, CRU_variables *cru)
{
    for (int n = 0; n < N; n++)
    {
        cru->TT_map[n]          = 1;
        cru->SERCA_map[n]       = 1.0;
        cru->NCX_map[n]         = 1.0;
        cru->RyR_het_map[n]     = 1.0;
        cru->LTCC_map[n]        = 1.0;
        cru->dyad_het_map[n]    = 0.5; // this will contain rand, where 0.5 = normal
    }
}

void set_sub_cellular_local_scale(Cell_parameters p, Dyad_variables *d, Membrane_fluxes *m, SR_fluxes *sr)
{
    // Default local values to global values
    d->GRyR_kCO         = p.GRyR_kCO;
    d->GLTCC_kva1_va2   = p.GLTCC_kva1_va2;
    sr->Gup             = p.Gup;
    sr->Gleak           = p.Gleak;
    m->GNCX             = p.GNCX;
    m->GCab             = p.GCab;
    m->GCaP             = p.GCaP;
    d->vol_ds           = p.vds_CRU_mean;
}

void read_sub_cellular_het_maps(SC_variables sc, CRU_variables *cru, const char *PATH, const char *Outputs_dir)
{
    double Nmap;
    if (strcmp(cru->SERCA_het, "On") == 0)
    {
        if (strcmp(cru->SERCA_map_file, "none") == 0)
        {
            printf("ERROR: SERCA heterogeneity is On but no map file has been specified\n");
            exit(1);
        }
        else Nmap    =  read_map_file_double(sc, cru->SERCA_map, cru->SERCA_map_file, "Sub_cellular_het_geometries", PATH, Outputs_dir, "SERCA"); // lib/Spatial_coupling.cpp
    }
    if (strcmp(cru->NCX_het, "On") == 0)
    {
        if (strcmp(cru->NCX_map_file, "none") == 0)
        {
            printf("ERROR: NCX heterogeneity is On but no map file has been specified\n");
            exit(1);
        }
        else Nmap    =  read_map_file_double(sc, cru->NCX_map, cru->NCX_map_file, "Sub_cellular_het_geometries", PATH, Outputs_dir, "NCX"); // lib/Spatial_coupling.cpp
    }
    if (strcmp(cru->RyR_het, "map") == 0)
    {
        if (strcmp(cru->RyR_het_map_file, "none") == 0)
        {
            printf("ERROR: RyR heterogeneity is On but no map file has been specified\n");
            exit(1);
        }
        else Nmap    =  read_map_file_double(sc, cru->RyR_het_map, cru->RyR_het_map_file, "Sub_cellular_het_geometries", PATH, Outputs_dir, "RyR"); // lib/Spatial_coupling.cpp
    }
    if (strcmp(cru->LTCC_het, "map") == 0)
    {
        if (strcmp(cru->LTCC_map_file, "none") == 0)
        {
            printf("ERROR: LTCC heterogeneity is On but no map file has been specified\n");
            exit(1);
        }
        else Nmap    =  read_map_file_double(sc, cru->LTCC_map, cru->LTCC_map_file, "Sub_cellular_het_geometries", PATH, Outputs_dir, "LTCC"); // lib/Spatial_coupling.cpp
    }
    /*if (strcmp(cru->Detub, "On") == 0)
    {
        if (strcmp(cru->TT_map_file, "none") == 0)
        {
            printf("ERROR: TT map is On but no map file has been specified\n");
            exit(1);
        }
        else int temp = read_map_file(sc, cru->TT_map, cru->TT_map_file, "Sub_cellular_het_geometries" , PATH, Outputs_dir, "TT_map");   // lib/Spatial_coupling.cpp
    }*/
}

void set_sub_cellular_local_het_scale(Cell_parameters p, int N, Dyad_variables *d, Membrane_fluxes *m, SR_fluxes *sr, CRU_variables cru)
{
    // scale local values (which already have global scaling applied by previous function) by heterogeneity maps
    for (int n = 0; n < N; n++)
    {
        sr[n].Gup       *= cru.SERCA_map[n];
        //sr[n].Gleak     *= cru.SERCA_map[n];
        m[n].GNCX       *= cru.NCX_map[n];
        m[n].GNCX       *= cru.TT_map[n];
        m[n].GCab       *= cru.TT_map[n];
        m[n].GCaP       *= cru.TT_map[n];
    }
}
// End Sub-cellular local scale factors ===============================================//|

// Dyad heterogeneity =================================================================\\|
void set_dyad_heterogeneity(Cell_parameters p, Dyad_variables *d, double rand, int n, char const *directory, double *map, CRU_variables cru)
{
    FILE * dyad_het_out;
    char * filename       = (char*)malloc(500);
    sprintf(filename, "%s/Random_dyad_het_distributions.dat", directory);

    double num, N;
    double k, mean, propvar;

    // First, produce disttributions from which actual values will be calculated (for outputs and checking)
    if (n == 0)
    {
        dyad_het_out = fopen(filename, "wt");
        for (num = 0.000001; num <= 1.000001; num += 0.005)
        {
            // vol_ds
            mean    = p.vds_CRU_mean;
            propvar = 0.4; // up to +/- 0.4*vol
            k       = -0.09*(propvar/0.5); // -0.09 scales by 0.5*mean, so scale by propvar relative to this
            N       = (log((1.0/num) - 1)*k + 1)*mean;      // inverse
            fprintf(dyad_het_out, "%f %f ", num, N);

            // RyR
            mean    = p.NRyR_mean;
            propvar = p.NRyR_propvar; 
            k       = -0.09*(propvar/0.5); // -0.09 scales by 0.5*mean, so scale by propvar relative to this
            N       = (log((1.0/num) - 1)*k + 1)*mean;
            fprintf(dyad_het_out, "%f ", N);

            // LTCC
            mean    = p.NLTCC_mean;
            propvar = p.NLTCC_propvar; 
            k       = -0.09*(propvar/0.5); // -0.09 scales by 0.5*mean, so scale by propvar relative to this
            N       = (log((1.0/num) - 1)*k + 1)*mean;
            fprintf(dyad_het_out, "%f\n", N);
        }
        fclose(dyad_het_out);
    }

    sprintf(filename, "%s/Dyad_vol_NRyR_NLTCC_list.dat", directory);
    dyad_het_out = fopen(filename, "a");

    // Now set actual values according to dists, if random is on (want same random number in so that dyad_vol, NRyR and NLTCC all vary together if at all)
    // Vol ds
    if (strcmp(cru.volds_het, "On") == 0)
    {
        mean        = p.vds_CRU_mean;
        propvar     = 0.4;
        k           = -0.09*(propvar/0.5);
        d->vol_ds   = (log((1.0/rand) - 1)*k + 1)*mean;
    }
    //else d->vol_ds  = p.vds_CRU_mean; // already set to mean by set_local_scale

    // RyR
    if (strcmp(cru.RyR_het, "random") == 0)
    {
        mean        = p.NRyR_mean;
        propvar     = p.NRyR_propvar;
        k           = -0.09*(propvar/0.5);
        N           = (log((1.0/rand) - 1)*k + 1)*mean;
        d->NRyR     = (int)N;
    }
    // don't want else here, as don't want to overwrite value from map if read in.
    // Has been set to mean by default already

    // NLTCC
    if (strcmp(cru.LTCC_het, "random") == 0)
    {
        mean        = p.NLTCC_mean;
        propvar     = p.NLTCC_propvar;
        k           = -0.09*(propvar/0.5);
        N           = (log((1.0/rand) - 1)*k + 1)*mean;
        d->NLTCC    = (int)N;
    }

    if (strcmp(cru.volds_het, "On") == 0 || strcmp(cru.RyR_het, "random") == 0 || strcmp(cru.LTCC_het, "random") == 0)
        *map = rand;

    fprintf(dyad_het_out, "%f %d %f %d %d\n", rand, n, d->vol_ds, d->NRyR, d->NLTCC);

    fclose(dyad_het_out);
    free(filename);
}

void write_random_dyad_het_vtk(SC_variables sc, double *map, const char* Output_dir)
{
    char *string = (char*)malloc(500);

    // Write map vtk to Outputs, with map type in filename (this is "ref" here)
    sprintf(string,"%s/Dyad_random_het_map.vtk", Output_dir);

    FILE *out;
    out = fopen(string, "wt");
    
    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    fprintf(out, "SCALARS %s float 1\n", string);
    fprintf(out, "LOOKUP_TABLE default\n");


    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                int idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
                fprintf(out, "%f ", map[idx]);
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }
    fclose(out);
    free(string);
}
// End Dyad heterogeneity =============================================================//|

// Initial conditions =================================================================\\|
void initial_conditions_calcium(Ca_variables *Ca, Cell_parameters p, int NCRU)
{
    // p.Cai/CaSR defined in lib/Initialisation.c
    Ca->DS			= p.Cai;		// uM
    Ca->SS			= p.Cai;		// uM
    Ca->CYTO		= p.Cai;		// uM
    Ca->NSR			= p.CaSR;		// uM
    Ca->JSR			= p.CaSR;		// uM

    for (int n = 0; n < NCRU; n++)
    {
        Ca->ds[n]	= Ca->DS;
        Ca->ss[n]	= Ca->SS;
        Ca->cyto[n]	= Ca->CYTO;
        Ca->nsr[n]	= Ca->NSR;
        Ca->jsr[n]	= Ca->JSR;
    }
}

void initial_conditions_calcium_0D(Ca_variables *Ca, Cell_parameters p)
{   
    // p.Cai/CaSR defined in lib/Initialisation.c
    Ca->DS          = p.Cai;        // uM
    Ca->SS          = p.Cai;        // uM
    Ca->CYTO        = p.Cai;        // uM
    Ca->NSR         = p.CaSR;       // uM
    Ca->JSR         = p.CaSR;       // uM
}

void initial_conditions_dyad_stochastic(Dyad_variables *d)
{
    for (int i = 0; i < d->NRyR; i++) 
    {
        d->RyR_state[i] = 0; // corresponds to CA state
    }
    d->Monomer = d->Mi 	= 0;
    d->active			= 0;

    for (int i = 0; i < d->NLTCC; i++) 
    {
        d->LTCC_va_state[i] = 0; // initially in closed state of voltage inactivation
        d->LTCC_vi_state[i] = 1; // initially in not inactivated state
        d->LTCC_ci_state[i] = 1;
    }
}

void initial_conditions_dyad_det(Dyad_variables *d)
{
    d->ICaL_va_0    = 1.0;
    d->ICaL_va_1    = 0.0;
    d->ICaL_va_2    = 0.0;
    d->ICaL_vi      = 1.0;
    d->ICaL_ci      = 1.0;
    d->NRyR_O1_det  = 0.0;
    d->Monomer      = 0.0;
    d->Mi           = 0.0;
    d->SRF_prop_active = 0.0;
}

void assign_CRU_variables_from_state_read(Dyad_variables *d, Ca_variables *Ca, State_variables s)
{
    // State struct is passed into read/write - so we need to relate our actual variables to the state equivilents
    d->Monomer    = s.Myo_m;
    d->Mi         = s.Myo_c;
    Ca->CYTO       = s.Cai;
    Ca->SS         = s.Cai_sl;
    Ca->DS         = s.Cai_j;
    Ca->NSR        = s.CanSR;
    Ca->JSR        = s.CajSR;
    d->ICaL_va_1  = s.ICaL_va;
    d->ICaL_va_0  = s.ICaL_vi_s;
    d->ICaL_va_2  = 1 - d->ICaL_va_1 - d->ICaL_va_0;
    d->ICaL_vi    = s.ICaL_vi;
    d->ICaL_ci    = s.ICaL_ci;
}

void assign_state_variables_from_CRU_write(Dyad_variables d, Ca_variables Ca, State_variables *s)
{
    // State struct is passed into read/write - so we need to relate our actual variables to the state equivilents
    s->Cai       = Ca.CYTO;
    s->Cai_sl    = Ca.SS;
    s->Cai_j     = Ca.DS;
    s->CanSR     = Ca.NSR;
    s->CajSR     = Ca.JSR;
    s->Myo_m     = d.Monomer;
    s->Myo_c     = d.Mi;
    s->ICaL_va   = d.ICaL_va_1;
    s->ICaL_vi_s = d.ICaL_va_0; // NOTE: vi_s set for va_0
    s->ICaL_vi   = d.ICaL_vi;
    s->ICaL_ci   = d.ICaL_ci;
}
// End Initial conditions =============================================================//|

// Inter-compartment transfer functions ===============================================\\|
void comp_J_ds_ss(Cell_parameters p, double Ca_ds, double Ca_ss, double vol_ds, double *reac_ss)
{
    *reac_ss		+= ((Ca_ds - Ca_ss)/p.tau_ds) * (vol_ds/p.vss_CRU);
}

void comp_J_ss_cyto(Cell_parameters p, double Ca_ss, double Ca_cyto, double *reac_ss, double *reac_cyto)
{
    *reac_ss      += -((Ca_ss - Ca_cyto)/p.tau_ss_cyt);          	// J_ss_cyto, SS reaction
    *reac_cyto    +=  ((Ca_ss - Ca_cyto)/p.tau_ss_cyt)*p.v_ss_cyt;  // J_ss_cyto, CYTO reaction
}

void comp_J_nsr_jsr(Cell_parameters p, double Ca_nsr, double Ca_jsr, double *reac_nsr, double *reac_jsr)
{
    *reac_jsr     +=  ((Ca_nsr - Ca_jsr)/p.tau_nsr_jsr);            	// J_nsr_jsr, JSR reaction
    *reac_nsr     += -((Ca_nsr - Ca_jsr)/p.tau_nsr_jsr)*p.v_jsr_nsr;    // J_nsr_jsr, NSR reaction
}
// End Inter-compartment transfer functions ===========================================//|

// Buffering ==========================================================================\\|
void comp_buffering(Cell_parameters p, double *Bcyto, double *Bss, double *Bjsr, double Ca_cyto, double Ca_ss, double Ca_jsr)
{
    // Based on work by:
    // Wagner and Keizer, Biophys J. 1994 Jul;67(1):447-56
    // Nivala et al. Front. Physiol. 2012;3:114.

    buffering_cyto(p, Bcyto, Ca_cyto);
    buffering_subspace(p, Bss, Ca_ss);
    buffering_JSR(p, Bjsr, Ca_jsr);
}

void buffering_cyto(Cell_parameters p, double *Bcyto, double Ca)
{
    *Bcyto   =  (p.Kcam*p.Bcam)/(pow(Ca + p.Kcam, 2));
    *Bcyto   += (p.Kbsr*p.Bbsr)/(pow(Ca + p.Kbsr, 2));
    *Bcyto   += (p.Kmca*p.Bmca)/(pow(Ca + p.Kmca, 2));
    *Bcyto   += (p.Kmmg*p.Bmmg)/(pow(Ca + p.Kmmg, 2));
    *Bcyto   += 1;
    *Bcyto   = 1/(*Bcyto);	
}

void buffering_subspace(Cell_parameters p, double *Bss, double Ca)
{
    *Bss   =  (p.Kcam*p.Bcam)/(pow(Ca + p.Kcam, 2));
    *Bss   += (1.5*p.Kbsr*p.Bbsr)/(pow(Ca + 1.5*p.Kbsr, 2));
    *Bss   += (0.5*p.Kmca*0.5*p.Bmca)/(pow(Ca + 0.5*p.Kmca, 2));
    *Bss   += (1.5*p.Kmmg*p.Bmmg)/(pow(Ca + 1.5*p.Kmmg, 2));
    *Bss   *= 0.1;
    *Bss   += 1;
    *Bss   = 1/(*Bss);
}

void buffering_JSR(Cell_parameters p, double *Bjsr, double Ca)
{
    // 1e-3 is to convert Ca_jsr to mM (buffering is dimensionless; params are in mM)
    *Bjsr    = 1/( 1 + (p.Kcsqn*p.Bcsqn)/pow((p.Kcsqn + 1e-3*Ca),2));
}
// End buffering ======================================================================//|

// Whole cell averages and currents ===================================================\\|
void calc_whole_cell_values_including_currents_from_flux(int N, Cell_parameters p, Ca_variables *Ca, CRU_variables *cru, Dyad_variables *d, SR_fluxes *sr, Membrane_fluxes *m, int NTOT)
{
    Ca->CYTO = Ca->SS = Ca->DS = Ca->NSR = Ca->JSR = 0;
    cru->J_SERCA = cru->J_LEAK = 0;
    cru->J_NCX_bulk	= cru->J_CaP_bulk = cru->J_Cab_bulk  = cru->J_NCX_ss = cru->J_CaP_ss = cru->J_Cab_ss = 0;
    cru->J_REL = 0;
    cru->PRyR_OA = cru->PRyR_OI = cru->PRyR_CA = cru->PRyR_CI = 0;
    cru->Monomer = cru->Mi = 0;
    cru->Nactive = 0;
    cru->Pactive = 0;
    cru->J_CAL = cru->PLTCC_O = cru->I_CAL = 0;

    for (int n = 0; n < N; n++)
    {
        Ca->CYTO		+= Ca->cyto[n]/N;
        Ca->SS			+= Ca->ss[n]/N;
        Ca->DS			+= Ca->ds[n]/N;
        Ca->NSR			+= Ca->nsr[n]/N;
        Ca->JSR			+= Ca->jsr[n]/N;

        cru->J_REL		+= d[n].J_rel/N;

        cru->PRyR_OA	+= float(float(d[n].NRyR_OA)/(N*d[n].NRyR));
        cru->PRyR_OI	+= float(float(d[n].NRyR_OI)/(N*d[n].NRyR));
        cru->PRyR_CA	+= float(float(d[n].NRyR_CA)/(N*d[n].NRyR));
        cru->PRyR_CI	+= float(float(d[n].NRyR_CI)/(N*d[n].NRyR));
        cru->Monomer	+= d[n].Monomer/N;
        cru->Mi			+= d[n].Mi/N;

        cru->Nactive	+= d[n].active; 	// == 0 if not active, 1 if active, sum = total active

        cru->J_CAL		+= d[n].J_CaL/N;
        cru->PLTCC_O	+= float(float(d[n].NLTCC_O)/(N*d[n].NLTCC));

        cru->J_SERCA 	+= sr[n].J_SERCA/N;
        cru->J_LEAK		+= sr[n].J_leak/N;

        // Membrane curents
        cru->J_NCX_bulk += m[n].J_NCX_bulk/N;
        cru->J_CaP_bulk += m[n].J_CaP_bulk/N;
        cru->J_Cab_bulk += m[n].J_Cab_bulk/N;
        cru->J_NCX_ss 	+= m[n].J_NCX_ss/N;
        cru->J_CaP_ss 	+= m[n].J_CaP_ss/N;
        cru->J_Cab_ss 	+= m[n].J_Cab_ss/N;

        // ICaL     || N = total CRUs in sim cell; NTOT = total in full cell
        // so for full cell, N = NTOT, but for a portion, N < NTOT and NTOT scales summed
        // current to whole cell
        cru->I_CAL		+= (-compute_current_from_flux(p, d[n].J_CaL, 2, d[n].vol_ds, NTOT))/N;

        // Note: For detub, we would only want to calc ave over NMEM (as any element without a TT
        // will have a current/flux of 0). E.g. if Nmem = Ntot/2, then the ave we would get summing as above
        // (over Ntot) would be the real ave/2. However, we would then scale up to full cell by multiplying by NMEM.
        // NMEM * correct ave is the same as Ntot * incorrect ave (= 2NMEM * correct_ave/2 = NMEM*correct ave)
        // Thus, we can just leave it as it is - our ave will be too small, but multiplied by an NTOT
        // which is too large by the same factor, cancelling out.
        // This is because we can only use a detub map with full sim cell size
    }

    cru->Pactive		= float(cru->Nactive)/N;

    // Other currents from fluxes
    cru->I_NCX_bulk		= compute_current_from_flux(p, cru->J_NCX_bulk, 1, p.vcyt_CRU, NTOT);
    cru->I_CaP_bulk		= compute_current_from_flux(p, cru->J_CaP_bulk, 2, p.vcyt_CRU, NTOT);
    cru->I_Cab_bulk		= compute_current_from_flux(p, cru->J_Cab_bulk, 2, p.vcyt_CRU, NTOT);

    cru->I_NCX_ss		= compute_current_from_flux(p, cru->J_NCX_ss, 1, p.vss_CRU, NTOT);
    cru->I_CaP_ss		= compute_current_from_flux(p, cru->J_CaP_ss, 2, p.vss_CRU, NTOT);
    cru->I_Cab_ss		= compute_current_from_flux(p, cru->J_Cab_ss, 2, p.vss_CRU, NTOT);
}

void calc_whole_cell_values_including_currents_from_flux_0D(Cell_parameters p, Ca_variables Ca, CRU_variables *cru, Dyad_variables d, SR_fluxes sr, Membrane_fluxes m, int NTOT)
{
    // Assign CRU value (for outputs) from Dyad, SR or MEM value
    cru->J_REL      = d.J_rel;

    cru->PRyR_OA    =  d.NRyR_O1_det*(1 - 0.75*d.Mi); 
    cru->PRyR_OI    = d.NRyR_O1_det * (0.25 + 0.75*d.Mi);
    cru->PRyR_CA    = (1 -  d.NRyR_O1_det)*(1 - 0.75*d.Mi);
    cru->PRyR_CI    = (1 -  d.NRyR_O1_det)*(0.25 + 0.75*d.Mi);
    cru->Monomer    = d.Monomer;
    cru->Mi         = d.Mi;

    cru->Nactive    = d.active;     // == 0 if not active, 1 if active, sum = total active
    cru->Pactive	= cru->Nactive; 

    cru->J_CAL      = d.J_CaL;
    cru->PLTCC_O    = d.ICaL_va_2;

    cru->J_SERCA    = sr.J_SERCA;
    cru->J_LEAK     = sr.J_leak;

    cru->J_NCX_bulk = m.J_NCX_bulk;
    cru->J_CaP_bulk = m.J_CaP_bulk;
    cru->J_Cab_bulk = m.J_Cab_bulk;
    cru->J_NCX_ss   = m.J_NCX_ss;
    cru->J_CaP_ss   = m.J_CaP_ss;
    cru->J_Cab_ss   = m.J_Cab_ss;

    // Currents from flux
    cru->I_CAL      	= (-compute_current_from_flux(p, cru->J_CAL, 2, d.vol_ds, NTOT));

    cru->I_NCX_bulk     = compute_current_from_flux(p, cru->J_NCX_bulk, 1, p.vcyt_CRU, NTOT);
    cru->I_CaP_bulk     = compute_current_from_flux(p, cru->J_CaP_bulk, 2, p.vcyt_CRU, NTOT);
    cru->I_Cab_bulk     = compute_current_from_flux(p, cru->J_Cab_bulk, 2, p.vcyt_CRU, NTOT);

    cru->I_NCX_ss       = compute_current_from_flux(p, cru->J_NCX_ss, 1, p.vss_CRU, NTOT);
    cru->I_CaP_ss       = compute_current_from_flux(p, cru->J_CaP_ss, 2, p.vss_CRU, NTOT);
    cru->I_Cab_ss       = compute_current_from_flux(p, cru->J_Cab_ss, 2, p.vss_CRU, NTOT);
}
// End Whole cell averages and currents ===============================================//|

// Current from flux ==================================================================\\|
double compute_current_from_flux(Cell_parameters p, double flux, int valence, double vol, int NCRUs)
{
    double I; // current, A/F
    I = flux * (1e-3/p.Cm) * valence * p.F * (vol) * NCRUs; //Cm in pF; F in C mmol-1; vol in micro m^3
    return I;
}
// End Current from flux ==============================================================//|
// End Whole CRU / global functions =============================================================//|

// Dyad fluxes functions ========================================================================\\|
// comp dyad ======================================================\\|
void comp_dyad_3D(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Ca_jsr, double Ca_cyto /*to which ds is coupled*/, double *reac_jsr, double Vm, double dt, const char *Model)
{
    // RyR model (stochastic) =======\\|
    set_and_update_monomer_state(p, d, 1e-3*Ca_jsr, dt);        	// Ca_jsr in mM || updates and sets monomer rates	
    set_RyR_rates(p, d, Ca_ds);										// sets transition rates
    update_RyR_stochastic(d, dt);                               	// Update state, monte-carlo
    d->K_rel		= d->NRyR_OA * (p.J_rel_max/(d->vol_ds));	// K_rel term
    //d->K_rel		*= d->Grel;										// Scale according to J_rel scaling  NO! - Grel now scales NRyR
    d->J_rel		= d->K_rel * (Ca_jsr - Ca_ds);					// Compute J_rel flux (uM/ms)
    *reac_jsr		+= -d->J_rel * (d->vol_ds/p.vjsr_CRU);			// jsr reaction term
    if (d->NRyR_OA > int(float(0.1*d->NRyR))) 	d->active = 1;		// in active state
    else 										d->active = 0;		// not in active state
    // End RyR model (stochastic) ===//|

    // LTCC model (stochastic) ======\\|
    set_LTCC_rates(p, d, Ca_ds, Vm, Model);							// Sets transition rates
    comp_LTCC_bar(p, d, Ca_ds, Vm);									// Sets dynamic flux rate
    update_LTCC_stochastic(d, dt);									// Update states, monte-carlo
    d->J_CaL		= d->LTCC_bar * -d->NLTCC_O;					// Flux through LTCC
    //d->J_CaL		*= p.GCaL;										// Scales flux  || NO! Have GCaL scale NLTCC, as more accurate
    // End TCC model (stochastic) ===//|
}

void comp_dyad_0D(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Ca_jsr, double Ca_cyto /*to which ds is coupled*/, double *reac_jsr, double Vm, double dt, const char *Model)
{
    // RyR model (deterministic) ====\\|
    // This definitely needs to be improved in a future model
    for (int i = 0 ; i < 5 ; i++)
    {
        set_RyR_rates(p, d, Ca_cyto + p.tau_ds*(d->J_CaL));              // sets transition rates -> passes in Ca_ds no flux approx with J_rel =0 to make it ICaL dependant only
        d->RyR_kOC_A	= p.RyR_kOC_det*(1 + 0.4/(1 + exp(-(d->Ca_JSR_t_ex - 1100)/32))); // CaSR open rate scale
        d->RyR_kOC_A 	*= 1 - 0.1*p.ISO;
        d->NRyR_O1_det += 0.2*dt*(d->RyR_kCO*(1-d->NRyR_O1_det)  - d->RyR_kOC_A*d->RyR_kOC*d->NRyR_O1_det);
    }

    // monomer inactivation
    set_and_update_monomer_state(p, d, 1e-3*Ca_jsr, dt);            // Ca_jsr in mM || updates and sets monomer rates
    d->Mi		+= dt* (d->mi_al*(1-d->Mi) - d->mi_b*d->Mi); 		// update Mi deterministically

    // Actual NRyR open
    d->NRyR_O_det	= d->NRyR * d->NRyR_O1_det * ( 1 - 0.75*d->SRF_prop_active) * (1 - (0.75*d->Mi)); // Mi contributes up to 75% total inac 

    // Select whether open RyR is det or SRF
    if (d->ex_switch == 0)
    {
        d->NRyR_O		= d->NRyR_O_SRF * d->NRyR;
        d->NRyR_O1_det	= d->NRyR_O_SRF;	
    }
    else d->NRyR_O		=  d->NRyR_O_det;	

    d->K_rel        = d->NRyR_O * (p.J_rel_max/(d->vol_ds));   	// K_rel term
    d->K_rel		*= d->Krel_SRF_mult;							// 1 for normal, < 1 if SRF in progress when excited
    d->J_rel        = d->K_rel * (Ca_jsr - Ca_ds);           		// Compute J_rel flux (uM/ms)
    *reac_jsr       += -d->J_rel * (d->vol_ds/p.vjsr_CRU);          // jsr reaction term

    if (d->NRyR_O > float(0.05*d->NRyR)) 	d->active = 1;
    else 									d->active = 0;
    // End RyR model ================//|

    // LTCC model (det; HH) =========\\|
    set_LTCC_rates(p, d, Ca_ds, Vm, Model);                         // Sets transition rates
    comp_LTCC_bar(p, d, Ca_ds, Vm);                                 // Sets dynamic flux rate
    update_gates_LTCC_det(p, d, dt);								// Updates gates
    //printf("%d %f %f %f %f\n", d->NLTCC, d->LTCC_bar, d->ICaL_va_2, d->ICaL_va_1 , d->ICaL_va_0);
    d->J_CaL		=  - (d->NLTCC * d->LTCC_bar * d->ICaL_va_2 * d->ICaL_vi * d->ICaL_ci ); // number channels * flux through single channel		
    // End LTCC model (det; HH) =====//|
}
// End comp dyad ==================================================//|

// RyR functions ==============================\\|
void set_and_update_monomer_state(Cell_parameters p, Dyad_variables *d, double Ca_jsr, double dt)
{
    d->csqn         = (p.Bcsqn*p.Kcsqn)/(Ca_jsr + p.Kcsqn);
    d->csqn_ca      = p.Bcsqn - d->csqn;

    d->mon_ss       = 1/(1.0 + exp(-6.5*(d->csqn-6.37)));
    d->mon_al       = d->mon_ss/p.RyR_mon_tau;
    d->mon_b        = (1-d->mon_ss)/(p.RyR_mon_beta_tau);

    d->Monomer      = d->Monomer + dt * ( d->mon_al*(1-d->Monomer) - d->mon_b * d->Monomer);

    d->mi_ss        = 1/(1.0 + exp(p.RyR_mon_grad*12.0*(d->Monomer-0.5)));
    d->mi_al        = (1-d->mi_ss)/p.RyR_mi_tau; // alpha is transition to mi (not away from it)
    d->mi_b         = d->mi_ss/(p.RyR_mi_beta_tau);
}

void set_RyR_rates(Cell_parameters p, Dyad_variables *d, double Ca_ds)
{
    // Model based on previous work:
    // Stern et al. J Gen Physiol. 1999 Mar;113(3):469-89.
    // Restrepo et al. Biophys J. 2008 Oct;95(8):3767-89.

    // Ca2+ activation
    d->RyR_kCO		= p.RyR_kCO_A * pow(Ca_ds, p.RyR_Cads_H); 	// transition rate closed to open
    d->RyR_kCO		*= d->GRyR_kCO;	// Scale rate	
    d->RyR_kOC		= p.RyR_kOC;								// transition rate open to closed			

    // Ca2+ SR (CSQN) inactivation
    d->RyR_kAI 		= d->mi_al;   // inactivation
    d->RyR_kIA 		= d->mi_b;    // de-inactivation
}

// Update RyR stochastic at end of file |||
// End RyR functions ==========================//|

// LTCC functions =============================\\|
void set_LTCC_rates(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Vm, const char * Model)
{
    // Model similar to the Markovian version of the HH model, presnted by
    // Song et al. Biophys J. 2015 Apr 21;108(8)1908-21.	

    double Vm_ac_ss         = Vm - p.ICaL_va_ss_shift;    // Voltage modified by shift applied to activation steady state
    double Vm_inac_ss       = Vm - p.ICaL_vi_ss_shift;    // Voltage modified by shift applied to inactivation steady state
    double Vm_ac_tau        = Vm - p.ICaL_va_tau_shift;   // Voltage modified by shift applied to activation time constant
    double Vm_inac_tau      = Vm - p.ICaL_vi_tau_shift;   // Voltage modified by shift applied to inactivation time constant

    // Voltage activation ===========\\|
    if (strcmp(Model, "hVM_ORD_s") == 0) 		set_ICaL_hVM_ORD_simple_va_rates(p, &d->ICaL_va_ss, &d->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);		
    else if (strcmp(Model, "hAM_CAZ_s") == 0) 	set_ICaL_hAM_CAZ_simple_va_rates(p, &d->ICaL_va_ss, &d->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);
    else 								 		set_Ip2d_va_rates(p, &d->ICaL_va_ss, &d->ICaL_va_tau, Vm_ac_ss, Vm_ac_tau, p.ICaL_va_ss_kscale);

    d->ICaL_va_al_01                = d->ICaL_va_ss/d->ICaL_va_tau;         // Rate from state va0 to va1 (V-dependent)
    d->ICaL_va_b_01                 = (1-d->ICaL_va_ss)/d->ICaL_va_tau;     // Rate from state va1 to va0
    d->ICaL_va_al_12                = p.LTCC_kva1_va2; 			            // Rate from va1 to va2 (V-independent)
    d->ICaL_va_al_12				*= d->GLTCC_kva1_va2;					// Scale rate
    d->ICaL_va_b_12                 = p.LTCC_kva2_va1;                      // Rate from va2 to va1 (V-independent)		
    // End Voltage activation =======//|

    // Voltage inactivation =========\\|
    if (strcmp(Model, "hAM_CAZ_s") == 0) 	set_ICaL_hAM_CAZ_simple_vi_rates(p, &d->ICaL_vi_ss, &d->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);
    else                                 	set_Ip2d_vi_rates(p, &d->ICaL_vi_ss, &d->ICaL_vi_tau, Vm_inac_ss, Vm_inac_tau, p.ICaL_vi_ss_kscale);

    d->ICaL_vi_al                   = d->ICaL_vi_ss/d->ICaL_vi_tau;
    d->ICaL_vi_b                    = (1-d->ICaL_vi_ss)/d->ICaL_vi_tau;
    // End Voltage inactivation =====//|

    // Ca activation ================\\|
    d->Ca_Ca_bar                    = Ca_ds/p.LTCC_Ca_bar;
    d->ICaL_ci_tau                  = 15; // ms
    if (strcmp(Model, "hAM_CAZ_s") == 0) 	d->ICaL_ci_tau = 30; // ms

    d->ICaL_ci_ss                   = 1/(1 + d->Ca_Ca_bar*d->Ca_Ca_bar);
    d->ICaL_ci_al                   = d->ICaL_ci_ss/d->ICaL_ci_tau;
    d->ICaL_ci_b                    = (1-d->ICaL_ci_ss)/d->ICaL_ci_tau;
    // End Ca activation ============//|
}

void comp_LTCC_bar(Cell_parameters p, Dyad_variables *d, double Ca_ds, double Vm)
{
    // Model based on Luo and Rudy Circ Res 1994 Jun;74(6):1071-96
    double Vm_in;
    Vm_in = Vm;
    if (Vm == 0.0) Vm_in = 1e-10; // otherwise ICaL_bar = 0!
    double z        =  Vm_in*p.FoRT;  // (Vm*F)(R*T)
    double Ca_eff;
    double gamma 	= 0.341;

    if (Ca_ds > 40) Ca_eff = 40; // maximum Ca "seen" by LTCC || for cleanliness
    else Ca_eff = Ca_ds;

    d->LTCC_bar = (1.0/(d->vol_ds)) * 4* p.J_LTCC_max * z * p.F * ( (0.5 * gamma * Ca_eff*1e-3 * exp(2*z) - gamma*p.Cao)/(exp(2*z) - 1));
}

// Update LTCC stochastic at end of file |||

// Update LTCC deterministic
void update_gates_LTCC_det(Cell_parameters p, Dyad_variables *d, double dt)
{
    // 2 state gates updated via rush larsen
    d->ICaL_vi          = rush_larsen(d->ICaL_vi, d->ICaL_vi_ss, d->ICaL_vi_tau, dt);
    d->ICaL_ci          = rush_larsen(d->ICaL_ci, d->ICaL_ci_ss, d->ICaL_ci_tau, dt);

    // 3 state gate updated Forward Euler
    //        f0 = f0           + dt*  (beta*f0                         - alpha*fi)
    d->ICaL_va_0 = d->ICaL_va_0 + dt * (d->ICaL_va_b_01 * d->ICaL_va_1  - d->ICaL_va_al_01 * d->ICaL_va_0);
    d->ICaL_va_2 = d->ICaL_va_2 + dt * (d->ICaL_va_al_12 * d->ICaL_va_1  - d->ICaL_va_b_12 * d->ICaL_va_2); // this acts from gate 1, hence alpha and beta other way round (alpha describes transition from 1 to 2)
    d->ICaL_va_1 = 1.0 - d->ICaL_va_0 - d->ICaL_va_2; // using that gates must sum to 1 to avoid the longer computation of gate va1
}
// End LTCC functions =========================//|
// End Dyad fluxes functions ====================================================================//|

// SR fluxes functions ==========================================================================\\|
void comp_SR_fluxes(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr, double *reac_cyto, double *reac_nsr)
{
    // Based on work by:
    // Shiferaw et al. Biophys J, 2003 Dec;85(6):3666-86.
    // Shannon et al. Biophys J. 2004 Nov;87(5):3351-71.

    comp_Jup(p, sr, Ca_cyto, Ca_nsr);                   		// compute Jup
    comp_Jleak(p, sr, Ca_cyto, Ca_nsr);                   		// compute Jleak
    *reac_cyto       += -(sr->J_SERCA - sr->J_leak);            // Update the reaction of bulk cyto with Jup-Jleak
    *reac_nsr        += (sr->J_SERCA - sr->J_leak)*p.v_cyt_nsr; // Update the reaction of nsr with Jup-Jleak || recall v_cyt_nsr is the ratio of volumes
}

void comp_Jup(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr)
{
    sr->cai_term	= pow(Ca_cyto/p.J_SERCA_kCa, 2.0);
    sr->casr_term	= pow(Ca_nsr/p.J_SERCA_kCaSR, 2.0);
    sr->J_SERCA		= p.J_SERCA_max * ((sr->cai_term - sr->casr_term)/(1 + sr->cai_term + sr->casr_term));
    sr->J_SERCA		*= sr->Gup; 
}

void comp_Jleak(Cell_parameters p, SR_fluxes *sr, double Ca_cyto, double Ca_nsr)
{
    sr->J_leak		= p.J_leak_max * ( (Ca_nsr*Ca_nsr)/( (Ca_nsr*Ca_nsr) + (p.J_leak_kCaSR*p.J_leak_kCaSR) ) )*(Ca_nsr - Ca_cyto);
    sr->J_leak		*= sr->Gleak;
}
// End SR fluxes functions ======================================================================//|

// Membrane fluxes functions ====================================================================\\|
void comp_membrane_fluxes(Cell_parameters p, Membrane_fluxes *m, State_variables s, double Ca_cyto, double Ca_ss, double *reac_cyto, double *reac_ss, double Vm, double SRF_mult)
{
    // Based on work by:
    // Shiferaw et al. Biophys J, 2003 Dec;85(6):3666-86.
    // Shannon et al. Biophys J. 2004 Nov;87(5):3351-71.

    comp_JMEM(p, m, s, Ca_cyto, Ca_ss,  Vm, SRF_mult);
    *reac_cyto	+= m->J_MEM_bulk;
    *reac_ss	+= m->J_MEM_ss;
}

void comp_JMEM(Cell_parameters p, Membrane_fluxes *m, State_variables s, double Ca_cyto, double Ca_ss, double Vm, double SRF_mult)
{
    double factor 		= (1.0/p.v_ss_cyt);

    m->J_NCX_bulk		= (1 - p.Fjunc) * 		   comp_JNCX(p, m, Ca_cyto, s, Vm); 	
    m->J_NCX_ss			= (	   p.Fjunc) * factor * comp_JNCX(p, m, Ca_ss, s, Vm); 
    m->J_NCX_bulk		*= m->GNCX;	
    m->J_NCX_ss			*= m->GNCX*SRF_mult;	
    m->J_NCX			= m->J_NCX_bulk + m->J_NCX_ss;

    m->J_Cab_bulk       = (1 - p.Fjunc) *          comp_JCab(p, Ca_cyto, Vm);
    m->J_Cab_ss         = (    p.Fjunc) * factor * comp_JCab(p, Ca_ss, Vm);
    m->J_Cab_bulk		*= m->GCab;	
    m->J_Cab_ss			*= m->GCab;	
    m->J_Cab            = m->J_Cab_bulk + m->J_Cab_ss;

    m->J_CaP_bulk		= (1 - p.Fjunc) * 		   comp_JCaP(p, Ca_cyto);
    m->J_CaP_ss			= (    p.Fjunc) * factor * comp_JCaP(p, Ca_ss);
    m->J_CaP_bulk		*= m->GCaP;	
    m->J_CaP_ss			*= m->GCaP;	
    m->J_CaP			= m->J_CaP_bulk + m->J_CaP_ss;

    m->J_MEM_bulk		= m->J_NCX_bulk - 	m->J_Cab_bulk 	- m->J_CaP_bulk;
    m->J_MEM_ss			= m->J_NCX_ss 	-	m->J_Cab_ss		- m->J_CaP_ss;
}

double comp_JNCX(Cell_parameters p, Membrane_fluxes *m, double Cai, State_variables s, double Vm)
{
    // Formulation from Shannon et al 2004, Biophys J 87:3351-3371
    double JNCX;

    double Ca_eff;
    Ca_eff = Cai;

    double Nai = s.Nai;
    double Nao = s.Nao;

    m->z = Vm*p.FoRT;  // (Vm*F)(R*T)
    m->Ka = 1/(1+pow((p.INCX_kda/Ca_eff),3));

    m->t1 = p.INCX_kCai * pow(Nao,3) * (1 + pow((Nai/p.INCX_kNai),3));
    m->t2 = pow(p.INCX_kNao,3) * Ca_eff * (1 + (Ca_eff/p.INCX_kCai));
    m->t3 = (p.INCX_kCao*pow(Nai,3)) + (pow(Nai,3) * p.Cao) + (pow(Nao,3) * Ca_eff);

    m->denomenator = (m->t1 + m->t2 + m->t3) * (1 + p.INCX_ksat * exp((p.INCX_gamma-1)*m->z));
    m->numerator = p.INCX_bar * m->Ka * ( (exp(p.INCX_gamma*m->z) * pow(Nai,3) * p.Cao) - (exp((p.INCX_gamma-1)*m->z) * pow(Nao,3) * Ca_eff )  );

    //printf("%.16f %.16f\n", (exp((p.INCX_gamma-1)*m->z) * pow(Nao,3) * Ca_eff ), m->z);

    JNCX = m->numerator/m->denomenator;

    return JNCX;
}

double comp_JCab(Cell_parameters p, double Cai, double Vm)
{
    double JCab;

    double ecab = ((p.R*p.T)/(2*p.F))*log(p.Cao/(1e-3*Cai)); // convert Cai to mM 
    JCab        = p.ICab_bar*(Vm-ecab);
    return JCab;
}

double comp_JCaP(Cell_parameters p, double Cai)
{
    double JCap;
    JCap        = (p.ICaP_bar*Cai)/(p.ICaP_kCa + Cai);
    return JCap;
}
// End Membrane fluxes functions ================================================================//|

// Stochastic integration =======================================================================\\|
// RyR ================================================================================\\|
void update_RyR_stochastic(Dyad_variables *d, double dt)
{
    d->NRyR_OA = d->NRyR_CA = d->NRyR_OI = d->NRyR_CI = 0; // zero total state occupancy

    for (int i = 0; i < d->NRyR; i++)
    {
        //d->rand = d->mtrand1(); // assign rand to  randomly generated number
        d->rand     = d->rand_RyR[i];
        if (d->rand > 1.0) printf("RAND GREATER THAN 1");
        else if (d->rand < 0.0) printf("RAND LESS THAN 0");

        switch (d->RyR_state[i])
        {
            case 0: // if currently in closed state, not inactivated (CA)
                if (d->rand <= d->RyR_kCO*dt) // if random number is less than probability of transition within dt, do the transition
                {
                    d->RyR_state[i] = 1; // do transition
                    d->NRyR_OA ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= (d->RyR_kCO + d->RyR_kAI)*dt) // else, if the random number is between the prob of transition 1 and 1+2, do transition 2
                {
                    d->RyR_state[i] = 2; // do transition -> state "2" is CI
                    d->NRyR_CI ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= 1) // otherwise, same state
                {
                    d->RyR_state[i] = 0; // no transition
                    d->NRyR_CA ++;          // count it as in closed state. No need to set state, as already there (but will anyway!)
                    break;
                }

            case 1: // if currently in open state, not inactivated (OA)
                if (d->rand <= d->RyR_kOC*dt) // if random number is less than probability of transition within dt, do the transition
                {
                    d->RyR_state[i] = 0; // do transition
                    d->NRyR_CA ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= (d->RyR_kOC + d->RyR_kAI)*dt) // else, if the random number is between the prob of transition 1 and 1+2, do transition 2
                {
                    d->RyR_state[i] = 3; // do transition -> state "3" is OI
                    d->NRyR_OI ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= 1) // 
                {
                    d->RyR_state[i] = 1; // no transition
                    d->NRyR_OA ++;          // count it as in closed state. No need to set state, as already there (but will anyway!)
                    break;
                }

            case 2: // if currently in closed state, inactivated (CI)
                if (d->rand <= d->RyR_kCO*dt) // if random number is less than probability of transition within dt, do the transition
                {
                    d->RyR_state[i] = 3; // do transition
                    d->NRyR_OI ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= (d->RyR_kCO + d->RyR_kIA)*dt) // else, if the random number is between the prob of transition 1 and 1+2, do transition 2
                {
                    d->RyR_state[i] = 0; // do transition -> state "2" is CI
                    d->NRyR_CA ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= 1) // otherwise, same state
                {
                    d->RyR_state[i] = 2; // no transition
                    d->NRyR_CI ++;          // count it as in closed state. No need to set state, as already there (but will anyway!)
                    break;
                }

            case 3: // if currently in open state, inactivated (OI)
                if (d->rand <= d->RyR_kOC*dt) // if random number is less than probability of transition within dt, do the transition
                {
                    d->RyR_state[i] = 2; // do transition
                    d->NRyR_CI ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= (d->RyR_kOC + d->RyR_kIA)*dt) // else, if the random number is between the prob of transition 1 and 1+2, do transition 2
                {
                    d->RyR_state[i] = 1; // do transition -> state "3" is OI
                    d->NRyR_OA ++;          // count transition
                    break;                  // break case statement
                }
                else if (d->rand <= 1) // 
                {
                    d->RyR_state[i] = 3; // no transition
                    d->NRyR_OI ++;          // count it as in closed state. No need to set state, as already there (but will anyway!)
                    break;
                }
        } // end switch 
    } // end for
} // end RyR stochastic
// End RyR ============================================================================//|

// LTCC ===============================================================================\\|
void update_LTCC_stochastic(Dyad_variables *d, double dt)
{
    d->NLTCC_O = 0;

    for (int i = 0; i < d->NLTCC; i++)
    {
        //d->rand = d->mtrand1(); // assign rand to  randomly generated number
        d->rand     = d->rand_LTCC[i];
        if (d->rand > 1.0) printf("RAND GREATER THAN 1");
        else if (d->rand < 0.0) printf("RAND LESS THAN 0");

        // This works by cyling through states, if no transition occurs in first state, check for transition in second state, then third
        if (d->LTCC_va_state[i] == 0) // closed state 0 voltage activation
        {
            if (d->rand <= d->ICaL_va_al_01 * dt) d->LTCC_va_state[i] = 1; // if rand if less than alpha rate va state 0-1, transition to state 1
            else
            {
                if (d->LTCC_vi_state[i] == 0) // if voltage inactivation state is in 0
                {
                    if (d->rand <= (d->ICaL_vi_al + d->ICaL_va_al_01) * dt) d->LTCC_vi_state[i] = 1;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_al + d->ICaL_va_al_01) * dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_al + d->ICaL_va_al_01) * dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
                else if (d->LTCC_vi_state[i] == 1)
                {
                    if (d->rand <= (d->ICaL_vi_b + d->ICaL_va_al_01) * dt) d->LTCC_vi_state[i] = 0;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_b + d->ICaL_va_al_01) * dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_b + d->ICaL_va_al_01) * dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
            }
        }
        else if (d->LTCC_va_state[i] == 1)
        {
            if (d->rand <= d->ICaL_va_b_01 * dt) d->LTCC_va_state[i] = 0;
            else if (d->rand <= (d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_va_state[i] = 2;
            else
            {
                if (d->LTCC_vi_state[i] == 0)
                {
                    if (d->rand <= (d->ICaL_vi_al + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_vi_state[i] = 1;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_al + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_al + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
                else if (d->LTCC_vi_state[i] == 1)
                {
                    if (d->rand <= (d->ICaL_vi_b + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_vi_state[i] = 0;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_b + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_b + d->ICaL_va_al_12 + d->ICaL_va_b_01) * dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
            }
        }
        else if (d->LTCC_va_state[i] == 2)
        {
            if (d->rand <= d->ICaL_va_b_12*dt) d->LTCC_va_state[i] = 1;
            else
            {
                if (d->LTCC_vi_state[i] == 0)
                {
                    if (d->rand <= (d->ICaL_vi_al + d->ICaL_va_b_12)*dt) d->LTCC_vi_state[i] = 1;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_al + d->ICaL_va_b_12)*dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_al + d->ICaL_va_b_12)*dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
                else if (d->LTCC_vi_state[i] == 1)
                {
                    if (d->rand <= (d->ICaL_vi_b + d->ICaL_va_b_12)*dt) d->LTCC_vi_state[i] = 0;
                    else
                    {
                        if (d->LTCC_ci_state[i] == 0)
                        {
                            if (d->rand <= (d->ICaL_ci_al + d->ICaL_vi_b + d->ICaL_va_b_12) * dt) d->LTCC_ci_state[i] = 1;
                        }
                        else if (d->LTCC_ci_state[i] == 1)
                        {
                            if (d->rand <= (d->ICaL_ci_b + d->ICaL_vi_b + d->ICaL_va_b_12 ) * dt) d->LTCC_ci_state[i] = 0;
                        }
                    }
                }
            }
        }

        if (d->LTCC_va_state[i] == 2 && d->LTCC_vi_state[i] == 1 && d->LTCC_ci_state[i] == 1) d->NLTCC_O ++; // if open state 2, inactivation state 1 for v and Ca, then channel is open
    } // end for
} // End LTCC stochastic
// End LTCC============================================================================//|
// End stochastic integration ===================================================================//|

// 3D cell voltage clamp ========================================================================\\|
void run_voltage_clamp_3Dcell(Cell_parameters p, Model_variables *var, State_variables *s, Dyad_variables *d, Ca_variables *Ca, CRU_variables *cru, SC_variables *sc, SR_fluxes *sr, Membrane_fluxes *m, RAND *rand, char const *directory, double dt)
{
    double Vm, Vclamp, Vhold, time, clamp_time, Vstart, Vend;
    double Ipeak, Ipeak2, Ipeak3, Ipeak4;
    char * filename       = (char*)malloc(500);

    FILE * I_out;
    FILE * I_out_individual;
    FILE * IV_out;

    // ICaL V clamp ===========================\\|
    printf("Applying ICaL voltage clamp\n");
    sprintf(filename, "%s/ICaL_Vclamp_traces.dat", directory);
    I_out = fopen(filename, "wt");
    sprintf(filename, "%s/ICaL_IV.dat", directory);
    IV_out = fopen(filename, "wt");

    // Vclamp settings
    Vstart      = -50;
    Vend        = 50;
    Vhold       = -50;
    clamp_time  = 500;

    for (Vclamp = -50; Vclamp < Vend+1; Vclamp +=5)
    {
        Ipeak = 0;
        initial_conditions_native(s, p, p.Model);    // lib/Model.c
        initial_conditions_calcium(Ca, p, sc->N);                                  
        for (int n = 0; n < sc->N; n++) initial_conditions_dyad_stochastic(&d[n]);    

        if (int(Vclamp)%10 == 0)
        {
            sprintf(filename, "%s/ICaL_Vclamp_trace_%.0f.dat", directory, Vclamp);
            I_out_individual = fopen(filename, "wt");
        }
        for (time = 0.0; time < 3*clamp_time; time += dt)
        {
            if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
            else Vm = Vclamp;

            // Assign Ca state variables (seen by ionic model) from integrated whole-cell ave variables
            s->Cai       = 1e-3*Ca->CYTO;     // Ca dependent currents, Cai (in mM not uM)
            s->CanSR     = 1e-3*Ca->NSR;      // Ca dependent currents, Cansr (in mM not uM)
            s->CajSR     = 1e-3*Ca->JSR;      // Ca dependent currents, Cajsr (in mM not uM)

#pragma omp parallel for default(none) shared(sc, d, rand)
            for (int n = 0; n < sc->N; n++)
            {
                for (int j = 0; j < d[n].NRyR; j++)  d[n].rand_RyR[j]     = rand[n].mtrand1(); // allows faster parallelisation
                for (int j = 0; j < d[n].NLTCC; j++) d[n].rand_LTCC[j]    = rand[n].mtrand1(); // as calling mtrand within functions seems slower
            }

#pragma omp parallel for default(none) shared(sc, Vm, p, var, s, Ca, time, d, sr, m, dt)
            // Spatial loop 1
            for (int n = 0; n < sc->N; n++)
            {
                // Zero reaction terms so they can be sequentially modified
                Ca->ss_reac[n] = Ca->cyto_reac[n] = Ca->nsr_reac[n] = Ca->jsr_reac[n] = 0;

                // Inter-compartment transfer || lib/CRU.cpp
                comp_J_ds_ss(p, Ca->ds[n], Ca->ss[n], d[n].vol_ds, &Ca->ss_reac[n]);
                comp_J_ss_cyto(p, Ca->ss[n], Ca->cyto[n], &Ca->ss_reac[n], &Ca->cyto_reac[n]);
                comp_J_nsr_jsr(p, Ca->nsr[n], Ca->jsr[n], &Ca->nsr_reac[n], &Ca->jsr_reac[n]);

                // Comp dyad || lib/CRU.cpp -> computes and solves JCaL, Jrel and Cads/jsr fluxes
                comp_dyad_3D(p, &d[n], Ca->ds[n], Ca->jsr[n], Ca->ss[n] /*to which ds is coupled*/, &Ca->jsr_reac[n], Vm, dt, p.Model);

                // Buffering || lib/CRU.cpp
                comp_buffering(p, &Ca->bcyto[n], &Ca->bss[n], &Ca->bjsr[n], Ca->cyto[n], Ca->ss[n], Ca->jsr[n]);

                // Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
                comp_SR_fluxes(p, &sr[n], Ca->cyto[n], Ca->nsr[n], &Ca->cyto_reac[n], &Ca->nsr_reac[n]);

                // Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
                comp_membrane_fluxes(p, &m[n], *s, Ca->cyto[n], Ca->ss[n], &Ca->cyto_reac[n], &Ca->ss_reac[n], Vm, 1.0);  // 1.0 is SRF mult as not relevant here

                // Spatial coupling ===============\\|
                // lib/Spatial_coupling.cpp
                calc_diff_FDM_tau(sc, Ca->ss, n, p.tau_ss_trans, p.tau_ss_long);
                Ca->ss_reac[n] += sc->diff[n];

                calc_diff_FDM_tau(sc, Ca->cyto, n, p.tau_cyto_trans, p.tau_cyto_long);
                Ca->cyto_reac[n] += sc->diff[n];

                calc_diff_FDM_tau(sc, Ca->nsr, n, p.tau_nsr_trans, p.tau_nsr_long);
                Ca->nsr_reac[n] += sc->diff[n];
                // End spatial coupling ===========//|
            }

            // Spatial loop 2
            for (int n = 0; n < sc->N; n++)
            {
                // Update local concentrations
                Ca->ds[n]    = (Ca->ss[n] + p.tau_ds*(d[n].K_rel*Ca->jsr[n] + d[n].J_CaL))/(1 + p.tau_ds*d[n].K_rel); // quasi-steady-state approx
                Ca->ss[n]    = Ca->ss[n]      + Ca->bss[n]     *   dt*(Ca->ss_reac[n]);
                Ca->cyto[n]  = Ca->cyto[n]    + Ca->bcyto[n]   *   dt*(Ca->cyto_reac[n]);
                Ca->nsr[n]   = Ca->nsr[n]     +                    dt*(Ca->nsr_reac[n]);
                Ca->jsr[n]   = Ca->jsr[n]     + Ca->bjsr[n]    *   dt*(Ca->jsr_reac[n]);
            }

            // Whole-cell averages || including computing currents from Ca fluxes
            calc_whole_cell_values_including_currents_from_flux(sc->N, p, Ca, cru, d, sr, m, cru->NTOT_CRUs); // lib/CRU.cpp

            // Assign currents for use in AP model from whole-cell averages
            var->ICaL  = cru->I_CAL;
            var->INCX  = cru->I_NCX_bulk + cru->I_NCX_ss;
            var->ICaP  = cru->I_CaP_bulk + cru->I_CaP_ss;
            var->ICab  = cru->I_Cab_bulk + cru->I_Cab_ss;
            // End Spatial Ca handling ====================================================//|

            // Solve the AP model || lib/Model.c -> lib/Model_X.cpp
            // This sets and updates all gates, and calculates Itot
            compute_model_integrated(p, var, s, Vm, dt);

            if (Vm == Vclamp && var->ICaL < Ipeak) Ipeak = var->ICaL;

            fprintf(I_out, "%f %f %f %f %f\n", time, Vm, var->ICaL, cru->PRyR_OA, 1000*s->Cai);
            if (int(Vclamp)%10 == 0) fprintf(I_out_individual, "%f %f %f %f %f\n", time, Vm, var->ICaL, cru->PRyR_OA, 1000*s->Cai);
        }
        fprintf(IV_out, "%f %f\n", Vclamp, Ipeak);
        if (int(Vclamp)%10 == 0) fclose(I_out_individual);
    }

    fclose(I_out);
    fclose(IV_out);
    // End ICaL V clamp =======================//|

    // Ito/IKur/IKr/IKs V clamp================\\|
    printf("Applying IK voltage clamp\n");
    sprintf(filename, "%s/IK_Vclamp_traces.dat", directory);
    I_out = fopen(filename, "wt");
    sprintf(filename, "%s/IK_IV.dat", directory);
    IV_out = fopen(filename, "wt");

    // Vclamp settings
    Vstart      = -50;
    Vend        = 70;
    Vhold       = -50;
    clamp_time  = 100;

    for (Vclamp = -50; Vclamp < Vend+1; Vclamp +=10)
    {
        Ipeak = Ipeak2 = Ipeak3 = Ipeak4 = 0;
        initial_conditions_native(s, p, p.Model);    // lib/Model.c
        sprintf(filename, "%s/Ito_Vclamp_trace_%.0f.dat", directory, Vclamp);
        I_out_individual = fopen(filename, "wt");
        for (time = 0.0; time < 3*clamp_time; time += dt)
        {
            if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
            else Vm = Vclamp;

            compute_model_integrated(p, var, s, Vm, dt);

            if (Vm == Vclamp && var->Ito  > Ipeak)  Ipeak  = var->Ito;
            if (Vm == Vclamp && var->IKur > Ipeak2) Ipeak2 = var->IKur;
            if (Vm == Vclamp && var->IKr  > Ipeak3) Ipeak3 = var->IKr;
            if (Vm == Vclamp && var->IKs  > Ipeak4) Ipeak4 = var->IKs;

            fprintf(I_out, "%f %f %f %f %f %f %f\n", time, Vm, var->Ito, var->IKur, var->Ito+var->IKur, var->IKr, var->IKs);
            fprintf(I_out_individual, "%f %f %f %f %f %f %f\n", time, Vm, var->Ito, var->IKur, var->Ito+var->IKur, var->IKr, var->IKs);
        }
        fprintf(IV_out, "%f %f %f %f %f\n", Vclamp, Ipeak, Ipeak2, Ipeak3, Ipeak4);
        fclose(I_out_individual);
    }

    fclose(I_out);
    fclose(IV_out);
    // End Ito/Kur/Kr/Ks V clamp===============//|

    // IK1 ====================================\\|
    printf("Applying IK1 voltage clamp\n");
    sprintf(filename, "%s/IK1_IV.dat", directory);
    IV_out = fopen(filename, "wt");
    for (Vclamp = -150; Vclamp < 100; Vclamp +=0.1)
    {
        Vm = Vclamp;
        compute_model_integrated(p, var, s, Vm, dt);
        fprintf(IV_out, "%f %f\n", Vm, var->IK1);
    }
    fclose(IV_out);
    // End IK1 ================================//|

    free(filename);
}
// End 3D cell voltage clamp ====================================================================//|

// 0D cell voltage clamp ========================================================================\\|
void run_voltage_clamp_0Dcell(Cell_parameters p, Model_variables *var, State_variables *s, Dyad_variables *d, Ca_variables *Ca, CRU_variables *cru, SR_fluxes *sr, Membrane_fluxes *m, char const *directory, double dt)
{
    double Vm, Vclamp, Vhold, time, clamp_time, Vstart, Vend;
    double Ipeak, Ipeak2;
    char * filename       = (char*)malloc(500);

    FILE * I_out;
    FILE * I_out_individual;
    FILE * IV_out;

    // ICaL V clamp ===========================\\|
    printf("Applying ICaL voltage clamp\n");
    sprintf(filename, "%s/ICaL_Vclamp_traces.dat", directory);
    I_out = fopen(filename, "wt");
    sprintf(filename, "%s/ICaL_IV.dat", directory);
    IV_out = fopen(filename, "wt");

    // Vclamp settings
    Vstart      = -50;
    Vend        = 50;
    Vhold       = -50;
    clamp_time  = 500;

    for (Vclamp = -50; Vclamp < Vend+1; Vclamp +=5)
    {
        Ipeak = 0;
        initial_conditions_native(s, p, p.Model);    // lib/Model.c
        initial_conditions_calcium_0D(Ca, p);       // lib/CRU.cpp
        initial_conditions_dyad_det(d);

        if (int(Vclamp)%10 == 0)
        {
            sprintf(filename, "%s/ICaL_Vclamp_trace_%.0f.dat", directory, Vclamp);
            I_out_individual = fopen(filename, "wt");
        }
        for (time = 0.0; time < 3*clamp_time; time += dt)
        {
            if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
            else Vm = Vclamp;

            // Assign Ca state variables (seen by ionic model) from integrated whole-cell ave variables
            s->Cai       = 1e-3*Ca->CYTO;     // Ca dependent currents, Cai (in mM not uM)
            s->CanSR     = 1e-3*Ca->NSR;      // Ca dependent currents, Cansr (in mM not uM)
            s->CajSR     = 1e-3*Ca->JSR;      // Ca dependent currents, Cajsr (in mM not uM)

            // Zero reaction terms so they can be sequentially modified
            Ca->SS_reac = Ca->CYTO_reac = Ca->NSR_reac = Ca->JSR_reac = 0;

            // Inter-compartment transfer || lib/CRU.cpp
            comp_J_ds_ss(p, Ca->DS, Ca->SS, d->vol_ds, &Ca->SS_reac);
            comp_J_ss_cyto(p, Ca->SS, Ca->CYTO, &Ca->SS_reac, &Ca->CYTO_reac);
            comp_J_nsr_jsr(p, Ca->NSR, Ca->JSR, &Ca->NSR_reac, &Ca->JSR_reac);

            // Comp dyad || lib/CRU.cpp -> computes and solves JCaL, Jrel and Cads/jsr fluxes
            comp_dyad_0D(p, d, Ca->DS, Ca->JSR, Ca->SS /*to which ds is coupled*/, &Ca->JSR_reac, Vm, dt, p.Model);

            // Buffering || lib/CRU.cpp
            comp_buffering(p, &Ca->Bcyto, &Ca->Bss, &Ca->Bjsr, Ca->CYTO, Ca->SS, Ca->JSR);

            // Comp SR fluxes || Jup, Jleak (SERCA) || lib/CRU.cpp
            comp_SR_fluxes(p, sr, Ca->CYTO, Ca->NSR, &Ca->CYTO_reac, &Ca->NSR_reac);

            // Comp Membrane fluxes || JNCX, JCaP, JCab || lib/CRU.cpp
            comp_membrane_fluxes(p, m, *s, Ca->CYTO, Ca->SS, &Ca->CYTO_reac, &Ca->SS_reac, Vm, 1.0);  // 1.0 is SRF mult as not relevant here

            // Update local concentrations
            Ca->DS    = (Ca->SS + p.tau_ds*(d->K_rel*Ca->JSR + d->J_CaL))/(1 + p.tau_ds*d->K_rel); // quasi-steady-state approx
            Ca->SS    = Ca->SS      + Ca->Bss     *   dt*(Ca->SS_reac);
            Ca->CYTO  = Ca->CYTO    + Ca->Bcyto   *   dt*(Ca->CYTO_reac);
            Ca->NSR   = Ca->NSR     +                 dt*(Ca->NSR_reac);
            Ca->JSR   = Ca->JSR     + Ca->Bjsr    *   dt*(Ca->JSR_reac);

            // Whole-cell averages || including computing currents from Ca fluxes
            calc_whole_cell_values_including_currents_from_flux_0D(p, *Ca, cru, *d, *sr, *m, cru->NTOT_CRUs); // lib/CRU.cpp

            // Assign currents for use in AP model from whole-cell averages
            var->ICaL  = cru->I_CAL;
            var->INCX  = cru->I_NCX_bulk + cru->I_NCX_ss;
            var->ICaP  = cru->I_CaP_bulk + cru->I_CaP_ss;
            var->ICab  = cru->I_Cab_bulk + cru->I_Cab_ss;
            // End Spatial Ca handling ====================================================//|

            // Solve the AP model || lib/Model.c -> lib/Model_X.cpp
            // This sets and updates all gates, and calculates Itot
            compute_model_integrated(p, var, s, Vm, dt);

            if (Vm == Vclamp && var->ICaL < Ipeak) Ipeak = var->ICaL;

            fprintf(I_out, "%f %f %f %f %f\n", time, Vm, var->ICaL, cru->PRyR_OA, 1000*s->Cai);
            if (int(Vclamp)%10 == 0) fprintf(I_out_individual, "%f %f %f %f %f\n", time, Vm, var->ICaL, cru->PRyR_OA, 1000*s->Cai);
        }
        fprintf(IV_out, "%f %f\n", Vclamp, Ipeak);
        if (int(Vclamp)%10 == 0) fclose(I_out_individual);
    }

    fclose(I_out);
    fclose(IV_out);
    // End ICaL V clamp =======================//|

    // Ito/IKur V clamp =======================\\|
    printf("Applying Ito voltage clamp\n");
    sprintf(filename, "%s/Ito_Vclamp_traces.dat", directory);
    I_out = fopen(filename, "wt");
    sprintf(filename, "%s/Ito_IV.dat", directory);
    IV_out = fopen(filename, "wt");

    // Vclamp settings
    Vstart      = -50;
    Vend        = 70;
    Vhold       = -50;
    clamp_time  = 100;

    for (Vclamp = -50; Vclamp < Vend+1; Vclamp +=10)
    {
        Ipeak = Ipeak2 = 0;
        initial_conditions_native(s, p, p.Model);    // lib/Model.c
        sprintf(filename, "%s/Ito_Vclamp_trace_%.0f.dat", directory, Vclamp);
        I_out_individual = fopen(filename, "wt");
        for (time = 0.0; time < 3*clamp_time; time += dt)
        {
            if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
            else Vm = Vclamp;

            compute_model_integrated(p, var, s, Vm, dt);

            if (Vm == Vclamp && var->Ito > Ipeak) Ipeak = var->Ito;
            if (Vm == Vclamp && var->IKur > Ipeak2) Ipeak2 = var->IKur;

            fprintf(I_out, "%f %f %f %f %f\n", time, Vm, var->Ito, var->IKur, var->Ito+var->IKur);
            fprintf(I_out_individual, "%f %f %f %f %f\n", time, Vm, var->Ito, var->IKur, var->Ito+var->IKur);
        }
        fprintf(IV_out, "%f %f %f\n", Vclamp, Ipeak, Ipeak2);
        fclose(I_out_individual);
    }

    fclose(I_out);
    fclose(IV_out);
    // End Ito/Kur V clamp ====================//|

    // IK1 ====================================\\|
    printf("Applying IK1 voltage clamp\n");
    sprintf(filename, "%s/IK1_IV.dat", directory);
    IV_out = fopen(filename, "wt");
    for (Vclamp = -150; Vclamp < 100; Vclamp +=0.1)
    {
        Vm = Vclamp;
        compute_model_integrated(p, var, s, Vm, dt);
        fprintf(IV_out, "%f %f\n", Vm, var->IK1);
    }
    fclose(IV_out);
    // End IK1 ================================//|

    free(filename);
}
// End 0D cell voltage clamp ====================================================================//|
