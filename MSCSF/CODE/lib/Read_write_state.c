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

#include "Structs.h"
#include "Read_write_state.h"

// Function list ================================================================================\\|
//	Actual read/write functions
//	    Write_state_variables_native()
//	    Read_state_variables_native()
//	    Read_spatial_Ca_system_state()
//	
//	Filename and function call
//	    Write_state_single_cell_native()
//	    Read_state_single_cell_native()
//	    Write_state_single_cell_integrated()
//	    Read_state_single_cell_integrated()
//	    Write_state_single_cell_integrated_0D()
//	    Read_state_single_cell_integrated_0D()
//	
//	    Write_state_tissue_native_whole_tissue()
//	    Read_state_tissue_native_whole_tissue()
//	    Write_state_tissue_integrated_whole_tissue()
//	    Read_state_tissue_integrated_whole_tissue()
//	    Write_state_tissue_native_ave_tissue()
//	    Read_state_tissue_native_ave_tissue()
//	    Write_state_tissue_integrated_ave_tissue()
//	    Read_state_tissue_integrated_ave_tissue()
//	
//	    Write_state_phase()
//	    Read_state_phase()
//	    Write_state_phase_0D()
//	    Read_state_phase_0D()
// End Function list ============================================================================//|

// Functions which write the state ====================================================\\|
void Write_state_variables_native(State_variables s, FILE *out, const char * Model)
{
	fprintf(out, "%lf\t", s.Vm);
	fprintf(out, "%lf\t", s.INa_va);
	fprintf(out, "%lf\t", s.INa_vi_1);
	fprintf(out, "%lf\t", s.INa_vi_2);
	fprintf(out, "%lf\t", s.INaL_va);
	fprintf(out, "%lf\t", s.INaL_vi);
	fprintf(out, "%lf\t", s.Ito_va);
	fprintf(out, "%lf\t", s.Ito_vi);
	fprintf(out, "%lf\t", s.Ito_vi_s);
	fprintf(out, "%lf\t", s.Ito_vi_3);
	fprintf(out, "%lf\t", s.ICaL_va);
	fprintf(out, "%lf\t", s.ICaL_vi);
	fprintf(out, "%lf\t", s.ICaL_vi_s);
	fprintf(out, "%lf\t", s.ICaL_ci);
	fprintf(out, "%lf\t", s.ICaL_ci_j);
	fprintf(out, "%lf\t", s.IKur_va);
	fprintf(out, "%lf\t", s.IKur_vi);
	fprintf(out, "%lf\t", s.IKr_va);
	fprintf(out, "%lf\t", s.IKr_vi);
	fprintf(out, "%lf\t", s.IKs_va);		
	fprintf(out, "%lf\t", s.IKs_va_2);		
	fprintf(out, "%lf\t", s.IK1_va);		
	fprintf(out, "%lf\t", s.IKACh_va);		
	fprintf(out, "%lf\t", s.IKACh_vi);		
	fprintf(out, "%lf\t", s.If_va);		// 25

	fprintf(out, "%lf\t", s.Cai);
	fprintf(out, "%lf\t", s.Cai_j);
	fprintf(out, "%lf\t", s.Cai_sl);
	fprintf(out, "%lf\t", s.CajSR);
	fprintf(out, "%lf\t", s.CanSR);			// 30

	fprintf(out, "%lf\t", s.Nai);
	fprintf(out, "%lf\t", s.Nai_j);
	fprintf(out, "%lf\t", s.Nai_sl);
	fprintf(out, "%lf\t", s.Ki);			// 34

	fprintf(out, "%lf\t", s.Cao);
	fprintf(out, "%lf\t", s.Nao);
	fprintf(out, "%lf\t", s.Ko);			// 37

	fprintf(out, "%lf\t", s.RyRo);			
	fprintf(out, "%lf\t", s.RyRi);
	fprintf(out, "%lf\t", s.RyRr);
	fprintf(out, "%lf\t", s.Myo_c);
	fprintf(out, "%lf\t", s.Myo_m);
	fprintf(out, "%lf\t", s.Tn_CHc);
	fprintf(out, "%lf\t", s.Tn_CHm);
	fprintf(out, "%lf\t", s.Tn_CL);			// 45

	fprintf(out, "%lf\t", s.cmdn);
	fprintf(out, "%lf\t", s.trpn);
	fprintf(out, "%lf\t", s.csqn);			// 48

	fprintf(out, "%lf\t", s.CaCal);
	fprintf(out, "%lf\t", s.Catrop);
	fprintf(out, "%lf\t", s.Camg);
	fprintf(out, "%lf\t", s.Mgmg);
	fprintf(out, "%lf\t", s.CaCalse);		// 53

	// If you are adding a new model which has a new state variable, add it here as an If statement.
	// This is so that previous state files are still valid.
	// It was decided this would be still cleaner than having a specific function for every individual model
	// Ensure you add it to both write and read!!
	// if (strcmp(Model, "model ref") == 0)
	//{
	//	fprintf(out, "%lf\t", new gate);
	//	fprintf(out, "%lf\t", new gate 2;
	//	fprintf(out, "%lf\t", new gate 3);
	//}
    // else if (strcmp(Model, "model_ref") == 0) printf(....)
}
// End functions which write the state ================================================//|

// Functions which read the state =====================================================\\|
void Read_state_variables_native(State_variables *s, FILE *in, const char * Model)
{
	double temp;  // variable to read empty states into

	fscanf(in, "%lf\t", &s->Vm);
	fscanf(in, "%lf\t", &s->INa_va);
	fscanf(in, "%lf\t", &s->INa_vi_1);
	fscanf(in, "%lf\t", &s->INa_vi_2);
	fscanf(in, "%lf\t", &s->INaL_va);
	fscanf(in, "%lf\t", &s->INaL_vi);
	fscanf(in, "%lf\t", &s->Ito_va);
	fscanf(in, "%lf\t", &s->Ito_vi);
	fscanf(in, "%lf\t", &s->Ito_vi_s);
	fscanf(in, "%lf\t", &s->Ito_vi_3);
	fscanf(in, "%lf\t", &s->ICaL_va);
	fscanf(in, "%lf\t", &s->ICaL_vi);
	fscanf(in, "%lf\t", &s->ICaL_vi_s);
	fscanf(in, "%lf\t", &s->ICaL_ci);
	fscanf(in, "%lf\t", &s->ICaL_ci_j);
	fscanf(in, "%lf\t", &s->IKur_va);
	fscanf(in, "%lf\t", &s->IKur_vi);
	fscanf(in, "%lf\t", &s->IKr_va);
	fscanf(in, "%lf\t", &s->IKr_vi);
	fscanf(in, "%lf\t", &s->IKs_va);		
	fscanf(in, "%lf\t", &s->IKs_va_2);		
	fscanf(in, "%lf\t", &s->IK1_va);		
	fscanf(in, "%lf\t", &s->IKACh_va);
    fscanf(in, "%lf\t", &s->IKACh_vi);
    fscanf(in, "%lf\t", &s->If_va);		// 25

	fscanf(in, "%lf\t", &s->Cai);
	fscanf(in, "%lf\t", &s->Cai_j);
	fscanf(in, "%lf\t", &s->Cai_sl);
	fscanf(in, "%lf\t", &s->CajSR);
	fscanf(in, "%lf\t", &s->CanSR);			// 30

	fscanf(in, "%lf\t", &s->Nai);
	fscanf(in, "%lf\t", &s->Nai_j);
	fscanf(in, "%lf\t", &s->Nai_sl);
	fscanf(in, "%lf\t", &s->Ki);			// 34

	fscanf(in, "%lf\t", &s->Cao);
	fscanf(in, "%lf\t", &s->Nao);
	fscanf(in, "%lf\t", &s->Ko);			// 37

	fscanf(in, "%lf\t", &s->RyRo);
	fscanf(in, "%lf\t", &s->RyRi);
	fscanf(in, "%lf\t", &s->RyRr);
	fscanf(in, "%lf\t", &s->Myo_c);
	fscanf(in, "%lf\t", &s->Myo_m);
	fscanf(in, "%lf\t", &s->Tn_CHc);
	fscanf(in, "%lf\t", &s->Tn_CHm);
	fscanf(in, "%lf\t", &s->Tn_CL);			// 45

	fscanf(in, "%lf\t", &s->cmdn);
	fscanf(in, "%lf\t", &s->trpn);
	fscanf(in, "%lf\t", &s->csqn);			// 48

	fscanf(in, "%lf\t", &s->CaCal);
	fscanf(in, "%lf\t", &s->Catrop);
	fscanf(in, "%lf\t", &s->Camg);
	fscanf(in, "%lf\t", &s->Mgmg);
	fscanf(in, "%lf\t", &s->CaCalse);		// 53

    // If you are adding a new model which has a new state variable, add it here as an If statement
    // This is so that previous state files are still valid
    // It was decided this would be still cleaner than having a specific function for every individual model
    //if (strcmp(Model, "model_ref") == 0)
    //{
    //    fscanf(in, "%lf\t", &new gate);
    //    fscanf(in, "%lf\t", &new gate 2);
	//	fscanf(in, "%lf\t", &new gate 3);
    //}
    // else if (strcmp(Model, "model_ref") == 0) fscanf(....)


	// And for minimal model, assign to alternatively named state
    s->Ip0d_va          = s->INa_va;
    s->Ip0d_vi_1        = s->INa_vi_1;
    s->Ip0d_vi_2        = s->INa_vi_2;
    s->Ip1r_va          = s->Ito_va;
    s->Ip1r_vi          = s->Ito_vi;
    s->Ip2d_va          = s->ICaL_va;
    s->Ip2d_vi          = s->ICaL_vi;
    s->Ip2r_va          = s->IKur_va;
    s->Ip2r_vi          = s->IKur_vi;
    s->Ip3r_va          = s->IKr_va;
}

void Read_spatial_Ca_system_state(Dyad_variables *d, Ca_variables *Ca, int n, FILE *in)
{
        int m;
        fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %d\n", &Ca->cyto[n], &Ca->ss[n], &Ca->ds[n], &Ca->nsr[n], &Ca->jsr[n], &d->Monomer, &d->Mi, &d->active);
        for (m = 0; m < d->NRyR; m++) fscanf(in, "%d ", &d->RyR_state[m]);
        for (m = 0; m < d->NLTCC; m++) fscanf(in, "%d %d %d ", &d->LTCC_va_state[m], &d->LTCC_vi_state[m], &d->LTCC_ci_state[m]);
}
// End functions which read the state =================================================//|

// Functions which set filename and call write/read state =============================\\|
// Single cell ================================\\|
void Write_state_single_cell_native(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
	FILE *out;
	char *string = (char*)malloc(500);

	sprintf(string, "%s/State_files/Single_cell/Native_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

	out = fopen(string, "wt");

	if (out == NULL)
	{
		printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
		exit(1);
	}

	Write_state_variables_native(s, out, Model);

	fclose(out);
}

void Read_state_single_cell_native(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Single_cell/Native_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}

void Write_state_single_cell_integrated(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Single_cell/Integrated_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_single_cell_integrated(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Single_cell/Integrated_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}

void Write_state_single_cell_integrated_spatial(State_variables s, Cell_parameters p, Dyad_variables *d, Ca_variables Ca, int BCL, const char * PATH, const char *Model, const char * State_ref, int NX, int NY, int NZ, int N)
{
    FILE *out;
    char *string = (char*)malloc(500);
    int n, m;

    sprintf(string, "%s/State_files/Single_cell/Integrated_model_spatial_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_dimen_%d_%d_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, NX, NY, NZ, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    // First, write whole-cell variables
    Write_state_variables_native(s, out, Model);
    fprintf(out, "\n");

    // Then spatial
    for(n = 0; n < N; n++)
    {
        fprintf(out, "%lf %lf %lf %lf %lf %lf %lf %d\n", Ca.cyto[n], Ca.ss[n], Ca.ds[n], Ca.nsr[n], Ca.jsr[n], d[n].Monomer, d[n].Mi, d[n].active);
        for (m = 0; m < d[n].NRyR; m++) fprintf(out, "%d ", d[n].RyR_state[m]);
        fprintf(out, "\n");
        for (m = 0; m < d[n].NLTCC; m++) fprintf(out, "%d %d %d ", d[n].LTCC_va_state[m], d[n].LTCC_vi_state[m], d[n].LTCC_ci_state[m]);
        fprintf(out, "\n");
    }

    fclose(out);
}

void Read_state_single_cell_integrated_spatial(State_variables *s, Cell_parameters p, Dyad_variables *d, Ca_variables *Ca, int BCL, const char * PATH, const char *Model, const char * State_ref, int NX, int NY, int NZ, int N)
{
    FILE *in;
    char *string = (char*)malloc(500);
    int n, m;

    sprintf(string, "%s/State_files/Single_cell/Integrated_model_spatial_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_dimen_%d_%d_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, NX, NY, NZ, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    // First, read whole-cell variables
    Read_state_variables_native(s, in, Model);

    // Then spatial
    for(n = 0; n < N; n++)
    {
        Read_spatial_Ca_system_state(&d[n], Ca, n, in);
    }

    fclose(in);
}

void Write_state_single_cell_integrated_0D(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Single_cell/Integrated_0D_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_single_cell_integrated_0D(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Single_cell/Integrated_0D_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}
// End Single cell ============================//|

// tissue =====================================\\|
void Write_state_tissue_native_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Native_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Write_state_variables_native(s[n], out, p[n].Model);
        fprintf(out, "\n");
    }
    fclose(out);
}

void Read_state_tissue_native_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Native_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Read_state_variables_native(&s[n], in, p[n].Model);
    }
    fclose(in);
}

void Write_state_tissue_integrated_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Integrated_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Write_state_variables_native(s[n], out, p[n].Model);
        fprintf(out, "\n");
    }
    fclose(out);
}

void Read_state_tissue_integrated_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Integrated_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Read_state_variables_native(&s[n], in, p[n].Model);
    }
    fclose(in);
}

void Write_state_tissue_native_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Native_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_tissue_native_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Native_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}

void Write_state_tissue_integrated_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Integrated_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_tissue_integrated_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Integrated_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}
// End tissue =================================//|

// Tissue network model =======================\\|
void Write_state_tissue_native_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Native_net_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Write_state_variables_native(s[n], out, p[n].Model);
        fprintf(out, "\n");
    }
    fclose(out);
}

void Read_state_tissue_native_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Native_net_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Read_state_variables_native(&s[n], in, p[n].Model);
    }
    fclose(in);
}

void Write_state_tissue_integrated_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Integrated_net_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Write_state_variables_native(s[n], out, p[n].Model);
        fprintf(out, "\n");
    }
    fclose(out);
}

void Read_state_tissue_integrated_net_whole_tissue(State_variables *s, Cell_parameters *p, int BCL, const char * PATH, const char *Model, int N, const char* Tissue_order, const char* Tissue_model, const char* Tissue_type, const char *Orientation_type, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);
    int n;

    sprintf(string, "%s/State_files/Tissue/Integrated_net_model_%s_BCL_%d_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_%s_%s_%s_%s_ref_%s_state.dat", PATH, Model, BCL, p[0].ISO, p[0].ACh, p[0].Remodelling, p[0].Agent, p[0].Mutation, Tissue_order, Tissue_model, Tissue_type, Orientation_type, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    for(n = 0; n < N; n++)
    {
        Read_state_variables_native(&s[n], in, p[n].Model);
    }
    fclose(in);
}

void Write_state_tissue_native_net_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Native_net_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_tissue_native_net_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Native_net_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}

void Write_state_tissue_integrated_net_ave_tissue(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Integrated_net_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_tissue_integrated_net_ave_tissue(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/Tissue/Integrated_net_ave_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_env_%s_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, p.environment, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}
// End tissue network model ===================//|

// phase ======================================\\|
void Write_state_phase(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/phase_files/Native_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_phase_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, phase, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_phase(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/phase_files/Native_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_phase_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, phase, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}

void Write_state_phase_0D(State_variables s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref)
{
    FILE *out;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/phase_files/Integrated_0D_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_phase_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, phase, State_ref);

    out = fopen(string, "wt");

    if (out == NULL)
    {
        printf("Cannot create state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Write_state_variables_native(s, out, Model);
    fclose(out);
}

void Read_state_phase_0D(State_variables *s, Cell_parameters p, int BCL, const char * PATH, const char *Model, int phase, const char * State_ref)
{
    FILE *in;
    char *string = (char*)malloc(500);

    sprintf(string, "%s/State_files/phase_files/Integrated_0D_model_%s_BCL_%d_region_%s_ISO_%.2f_ACh_%.2f_remodelling_%s_drug_%s_mut_%s_phase_%d_ref_%s_state.dat", PATH, Model, BCL, p.Celltype, p.ISO, p.ACh, p.Remodelling, p.Agent, p.Mutation, phase, State_ref);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot open state file %s\t does the folder exist?? Is your path %s correct?\n", string, PATH);
        exit(1);
    }

    Read_state_variables_native(s, in, Model);
    fclose(in);
}
// End phase ==================================//|
// End Functions which set filename and call write/read state =========================//|

