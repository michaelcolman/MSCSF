// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Spatial coupling. Includes functions for: ===  //
//		array allocation ==================================  //
// 		reading geometries and maps =======================  //
//		setup of diffusion coefficient arrays =============  //
//		neighbour and array type mapps ====================  //
//		implementation of the anisotropic FDM method ======  //
// ========================================================  //
// GNU 3 LICENSE TEXT =====================================  //
// COPYRIGHT (C) 2016-2022 MICHAEL A. COLMAN ==============  //
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

#include "Spatial_coupling.h"
#include "Structs.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>

// Function list ================================================================================\\|
//	Array allocation
//	    SC_set_array_sizes()
//	    SC_array_allocation_N3()
//	    SC_array_allocation_Ncell()
//	    SC_array_deallocation()
//	
//	Read geometry/maps
//	    read_geo_file()
//	    read_map_file()
//	    read_map_file_double()
//	
//	Set diffusion tensor global and local
//	    set_D_dx_global()
//	    set_D_array_anisotropic()
//	    output_D1_and_D2()
//	
//	Set cell index and neigbours
//	    SC_set_index_and_geo_linear()
//	    SC_set_neighbours()
//	
//	FDM methods
//	    calc_diff_FDM_iso()
//	    calc_diff_FDM_tau()
//	    calc_dD_anisotropic_3D()
//	    calc_laplacian_and_BCs()
//
//  NETWORK model functions
//      list them all here
// End function list ============================================================================//|

// Allocate and deallocate spatial arrays =======================================================\\|
// Set array sizes ================================================\\|
void SC_set_array_sizes(SC_variables *sc, int NX, int NY, int NZ)
{
	sc->NX = NX;
	sc->NY = NY;
	sc->NZ = NZ;
}
// End set array sizes ============================================//|

// Arrays of size NX*NY*NZ ========================================\\|
void SC_array_allocation_N3(SC_variables *sc, int NX, int NY, int NZ)
{
	int N = NX * NY * NZ;

	// Geometry arrays
	sc->geo			= new int[N];
	sc->geo_index	= new int[N];
}
// End arrays of size NX*NY*NZ ====================================//|

// Arrays of size Ncell ===========================================\\|
void SC_array_allocation_Ncell(SC_variables *sc, int N)
{
	// Geometry arrays	
	sc->geo_linear 		= new int [N];  // linear array of celltypes
	sc->geo_3D_index	= new int [N];	// returns 3D index at 1D cell element
	sc->x_index			= new int [N];	// returns x coordinate of cell n
	sc->y_index			= new int [N];
	sc->z_index			= new int [N];

	// Diffusion coefficient array and D derivates etc
	sc->D     			= new double [N];
	sc->D1     			= new double [N];
	sc->D2     			= new double [N];
	sc->Dxx            	= new double [N];
	sc->Dyy            	= new double [N];
	sc->Dzz            	= new double [N];
	sc->Dxy            	= new double [N];
	sc->Dxz            	= new double [N];
	sc->Dyz            	= new double [N];

	sc->dDxx_dx        	= new double [N];
	sc->dDxy_dx        	= new double [N];
	sc->dDxz_dx        	= new double [N];
	sc->dDyy_dy        	= new double [N];
	sc->dDxy_dy        	= new double [N];
	sc->dDyz_dy        	= new double [N];
	sc->dDzz_dz        	= new double [N];
	sc->dDxz_dz        	= new double [N];
	sc->dDyz_dz        	= new double [N];

	// Differential
	sc->diff			= new double [N];

	// Neighbours
	sc->xp				= new int [N];
	sc->xm				= new int [N];
	sc->yp				= new int [N];
	sc->ym				= new int [N];
	sc->zp				= new int [N];
	sc->zm				= new int [N];

	sc->xp_yp          	= new int[N];
	sc->xp_ym          	= new int[N];
	sc->xp_zp          	= new int[N];
	sc->xp_zm          	= new int[N];

	sc->xm_yp          	= new int[N];
	sc->xm_ym          	= new int[N];
	sc->xm_zp          	= new int[N];
	sc->xm_zm          	= new int[N];

	sc->yp_zp          	= new int[N];
	sc->yp_zm          	= new int[N];
	sc->ym_zp          	= new int[N];
	sc->ym_zm          	= new int[N];

	sc->xm_ym_zm       	= new int[N];
	sc->xm_ym_zp       	= new int[N];
	sc->xm_yp_zm       	= new int[N];
	sc->xm_yp_zp       	= new int[N];
	sc->xp_ym_zm       	= new int[N];
	sc->xp_ym_zp       	= new int[N];
	sc->xp_yp_zm       	= new int[N];
	sc->xp_yp_zp       	= new int[N];

	// Orientation
	sc->ox				= new double [N];
	sc->oy				= new double [N];
	sc->oz				= new double [N];

	sc->ox2				= new double [N];
	sc->oy2				= new double [N];
	sc->oz2				= new double [N];
	sc->ox3				= new double [N];
	sc->oy3				= new double [N];
	sc->oz3				= new double [N];

	// Laplacian
	sc->lap_self 		= new double[N];
	sc->lap_xm 			= new double[N];
	sc->lap_xp 			= new double[N];
	sc->lap_ym 			= new double[N];
	sc->lap_yp 			= new double[N];
	sc->lap_zm 			= new double[N];
	sc->lap_zp 			= new double[N];
	sc->lap_xm_ym 		= new double[N];
	sc->lap_xm_yp 		= new double[N];
	sc->lap_xp_ym 		= new double[N];
	sc->lap_xp_yp 		= new double[N];
	sc->lap_xm_zm 		= new double[N];
	sc->lap_xm_zp 		= new double[N];
	sc->lap_xp_zm 		= new double[N];
	sc->lap_xp_zp 		= new double[N];
	sc->lap_ym_zm 		= new double[N];
	sc->lap_ym_zp 		= new double[N];
	sc->lap_yp_zm 		= new double[N];
	sc->lap_yp_zp 		= new double[N];

    // NETWORK
    sc->Gl                  = new double [N];
    sc->Gt                  = new double [N];
    sc->Gt2                 = new double [N];
    sc->gGap_node_xx       = new double [N];
    sc->gGap_node_yy       = new double [N];
    sc->gGap_node_zz       = new double [N];
    sc->gGap_node_xypp     = new double [N];
    sc->gGap_node_xypm     = new double [N];
    sc->gGap_node_xzpp     = new double [N];
    sc->gGap_node_xzpm     = new double [N];
    sc->gGap_node_yzpp     = new double [N];
    sc->gGap_node_yzpm     = new double [N];
    sc->gGap_node_xyzppp   = new double [N];
    sc->gGap_node_xyzppm   = new double [N];
    sc->gGap_node_xyzpmp   = new double [N];
    sc->gGap_node_xyzmpp   = new double [N];

    sc->connection_type_node_xx       = new int [N];
    sc->connection_type_node_yy       = new int [N];
    sc->connection_type_node_zz       = new int [N];
    sc->connection_type_node_xypp     = new int [N];
    sc->connection_type_node_xypm     = new int [N];
    sc->connection_type_node_xzpp     = new int [N];
    sc->connection_type_node_xzpm     = new int [N];
    sc->connection_type_node_yzpp     = new int [N];
    sc->connection_type_node_yzpm     = new int [N];
    sc->connection_type_node_xyzppp   = new int [N];
    sc->connection_type_node_xyzppm   = new int [N];
    sc->connection_type_node_xyzpmp   = new int [N];
    sc->connection_type_node_xyzmpp   = new int [N];
}
// End arrays of size Ncell =======================================//|

// Deallocate all arrays ==========================================\\|
void SC_array_deallocation(SC_variables *sc)
{
    // Geometry arrays, NX*NY*NZ
    delete []	sc->geo;
    delete []	sc->geo_index;

    // Geometry, Ncell
    delete [] 	sc->geo_linear;
    delete [] 	sc->geo_3D_index;
    delete [] 	sc->x_index;
    delete [] 	sc->y_index;
    delete [] 	sc->z_index;

    // Diffusion coefficient and spatial D derivatives
    delete []	sc->D;
    delete []	sc->D1;
    delete []	sc->D2;
    delete [] 	sc->Dxx; 
    delete [] 	sc->Dyy;
    delete [] 	sc->Dzz;
    delete [] 	sc->Dxy;
    delete [] 	sc->Dxz;
    delete [] 	sc->Dyz;

    delete [] 	sc->dDxx_dx;   
    delete [] 	sc->dDxy_dx;
    delete [] 	sc->dDxz_dx;
    delete [] 	sc->dDyy_dy;
    delete []	sc->dDxy_dy;
    delete [] 	sc->dDyz_dy;
    delete [] 	sc->dDzz_dz;
    delete [] 	sc->dDxz_dz;
    delete [] 	sc->dDyz_dz;

    // Differential
    delete []	sc->diff;

    // Neighbours
    delete []	sc->xp;
    delete []	sc->xm;
    delete []	sc->yp;
    delete []	sc->ym;
    delete []	sc->zp;
    delete []	sc->zm;

    delete [] 	sc->xp_yp;
    delete [] 	sc->xp_ym;
    delete [] 	sc->xp_zp;
    delete [] 	sc->xp_zm;
    delete [] 	sc->xm_yp;
    delete [] 	sc->xm_ym;
    delete [] 	sc->xm_zp;
    delete [] 	sc->xm_zm;
    delete [] 	sc->yp_zp;
    delete [] 	sc->yp_zm;
    delete [] 	sc->ym_zp;
    delete [] 	sc->ym_zm;

    delete [] 	sc->xm_ym_zm;
    delete [] 	sc->xm_ym_zp;
    delete [] 	sc->xm_yp_zm;
    delete [] 	sc->xm_yp_zp;
    delete [] 	sc->xp_ym_zm;
    delete [] 	sc->xp_ym_zp;
    delete [] 	sc->xp_yp_zm;
    delete [] 	sc->xp_yp_zp;

	// Orientation
	delete [] 	sc->ox;
	delete [] 	sc->oy;
	delete [] 	sc->oz;
	delete [] 	sc->ox2;
	delete [] 	sc->oy2;
	delete [] 	sc->oz2;
	delete [] 	sc->ox3;
	delete [] 	sc->oy3;
	delete [] 	sc->oz3;

	// Laplacian
	delete [] 	sc->lap_self;
	delete [] 	sc->lap_xm;
	delete [] 	sc->lap_xp;
	delete [] 	sc->lap_ym;
	delete [] 	sc->lap_yp;
	delete [] 	sc->lap_zm;
	delete [] 	sc->lap_zp;
	delete [] 	sc->lap_xm_ym;
	delete [] 	sc->lap_xm_yp;
	delete [] 	sc->lap_xp_ym;
	delete [] 	sc->lap_xp_yp;
	delete [] 	sc->lap_xm_zm;
	delete [] 	sc->lap_xm_zp;
	delete [] 	sc->lap_xp_zm;
	delete [] 	sc->lap_xp_zp;
	delete [] 	sc->lap_ym_zm;
	delete [] 	sc->lap_ym_zp;
	delete [] 	sc->lap_yp_zm;
	delete [] 	sc->lap_yp_zp;

    // NETWORK
    delete []   sc->Gl;
    delete []   sc->Gt;
    delete []   sc->Gt2;
    delete []   sc->gGap_node_xx;
    delete []   sc->gGap_node_yy;
    delete []   sc->gGap_node_zz;
    delete []   sc->gGap_node_xypp;
    delete []   sc->gGap_node_xypm;
    delete []   sc->gGap_node_xzpp;
    delete []   sc->gGap_node_xzpm;
    delete []   sc->gGap_node_yzpp;
    delete []   sc->gGap_node_yzpm;
    delete []   sc->gGap_node_xyzppp;
    delete []   sc->gGap_node_xyzppm;
    delete []   sc->gGap_node_xyzpmp;
    delete []   sc->gGap_node_xyzmpp;

    delete []   sc->connection_type_node_xx;
    delete []   sc->connection_type_node_yy;
    delete []   sc->connection_type_node_zz;
    delete []   sc->connection_type_node_xypp;
    delete []   sc->connection_type_node_xypm;
    delete []   sc->connection_type_node_xzpp;
    delete []   sc->connection_type_node_xzpm;
    delete []   sc->connection_type_node_yzpp;
    delete []   sc->connection_type_node_yzpm;
    delete []   sc->connection_type_node_xyzppp;
    delete []   sc->connection_type_node_xyzppm;
    delete []   sc->connection_type_node_xyzpmp;
    delete []   sc->connection_type_node_xyzmpp;
}
// End deallocate all arrays ======================================//|

// NETWORK
void SC_array_allocation_Njunc(SC_variables *sc, int N)
{
    // Geometry arrays
    sc->gGap_jn            = new double [N];
    sc->IGap               = new double [N];
    sc->jn_map_minus        = new int [N];
    sc->jn_map_plus         = new int [N];
    sc->connection_type_jn  = new int [N];
    sc->gGgap_mod_map       = new double [N];
    sc->gGgap_base_map       = new double [N];
}
void SC_array_deallocation_Njunc(SC_variables *sc)
{
    delete []   sc->gGap_jn;
    delete []   sc->IGap;
    delete []   sc->jn_map_minus;
    delete []   sc->jn_map_plus;
    delete []   sc->connection_type_jn;
    delete []   sc->gGgap_mod_map;
    delete []   sc->gGgap_base_map;
}
// End Allocate and deallocate spatial arrays ===================================================//|

// Read geometry/maps from file into arrays =====================================================\\|
int read_geo_file(SC_variables *sc, int *geo, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref, int Ncelltypes)
{
    FILE *in;
    char *string = (char*)malloc(500);

    // Assign filename to string
    sprintf(string, "%s/%s/%s", PATH, fileroot, filein);

    // Read in file
    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot load geometry file %s\t :: is the path correct? Does the file exist in that path?\n", string);
        exit(1);
    }
    else printf("File loaded %s\n", string);

    int temp; // temporary int for reading from file
    int idx, cell_count;
    cell_count = 0;

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {	
                idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);
                geo[idx] = 0;
                fscanf(in, "%d ", &temp);
                geo[idx] = temp;
                if (geo[idx] > 0) cell_count++; // how many real cells
            }
        }
    }

    printf("Geometry file %s read || Ncells = %d\n", filein, cell_count);

    fclose(in);

    // Write vtk of geoemtry to Outputs directory
    sprintf(string, "%s/Geometry_%s.vtk", Output_dir, ref);

    FILE *out;
    out = fopen(string, "wt");

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY, sc->NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc->NX*sc->NY*sc->NZ);
    fprintf(out, "SCALARS geometry float 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
				idx = i + (sc->NX*j) + (sc->NX*sc->NY*k);
				fprintf(out, "%d ", geo[idx]); 
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);

    // Now, write geo file for each celltype individually
    // Write vtk of geoemtry to Outputs directory
    for (int ncell = 1; ncell < Ncelltypes+1; ncell++) // from 1 to N, not 0 to N-1
    {
        sprintf(string, "%s/Geometry_%s_celltype_%d.vtk", Output_dir, ref, ncell);

        FILE *out;
        out = fopen(string, "wt");

        fprintf(out, "# vtk DataFile Version 3.0\n");
        fprintf(out, "vtk output\n");
        fprintf(out, "ASCII\n");
        fprintf(out, "DATASET STRUCTURED_POINTS\n");
        fprintf(out, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY, sc->NZ);
        fprintf(out, "SPACING 1 1 1\n");
        fprintf(out, "ORIGIN 0 0 0\n");
        fprintf(out, "POINT_DATA %d\n", sc->NX*sc->NY*sc->NZ);
        fprintf(out, "SCALARS geometry float 1\n");
        fprintf(out, "LOOKUP_TABLE default\n");

        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX*sc->NY*k);
                    if (geo[idx] == ncell) fprintf(out, "%d ", geo[idx]);
                    else fprintf(out, "0 ");
                }
                fprintf(out, "\n");
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }

	return cell_count;
}

// The difference here is that the array being read into is already size N  (not X*Y*Z - which the file input is)
// Reads integer map
int read_map_file(SC_variables sc, int *map, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref)
{
	FILE *in;
	char *string = (char*)malloc(500);

	// Assign filename to string
	sprintf(string, "%s/%s/%s", PATH, fileroot, filein);

	in = fopen(string, "r");

	if (in == NULL)
	{
		printf("Cannot load geometry map file %s\t :: is the path correct? Does the file exist in that path?\n", string);
		exit(1);
	}
	else printf("File loaded %s\n", string);

	int temp; // temporary int for reading from file
	int idx, cell_count, map_count;
	cell_count = 0;
	map_count = 0;

	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
				fscanf(in, "%d ", &temp);			

				if (sc.geo[idx] > 0)
				{
					map[cell_count] = temp; 	// reading into 1D array of size N
					if (map[cell_count] > 0) map_count++;
					cell_count++;
				}
			}
		}
	}

	printf("Map file %s read || Nmap = %d, recheck of Ncells = %d\n", filein, map_count, cell_count);

	fclose(in);

	// Write map vtk to Outputs, with map type in filename (this is "ref" here)
	sprintf(string,"%s/Map_%s.vtk", Output_dir, ref);

	FILE *out;
	out = fopen(string, "wt");

	sprintf(string,"%s", ref);

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

	cell_count = 0;

	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
				if (sc.geo[idx] > 0)
				{
					fprintf(out, "%d ", map[cell_count]);
					cell_count++;
				}
				else fprintf(out, "-100 ");
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);

	return map_count;
} // end read map

// And now, the difference is we are reading in a non-integer
int read_map_file_double(SC_variables sc, double *map, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref)
{
	FILE *in;
	char *string = (char*)malloc(500);

	// Assign filename to string
	sprintf(string, "%s/%s/%s", PATH, fileroot, filein);

	in = fopen(string, "r");

	if (in == NULL)
	{
		printf("Cannot load geometry map file %s\t :: is the path correct? Does the file exist in that path?\n", string);
		exit(1);
	}
	else printf("File loaded %s\n", string);

	double temp; // temporary int for reading from file
	int idx, cell_count, map_count;
	cell_count = 0;
	map_count = 0;

	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
				fscanf(in, "%lf ", &temp);

				if (sc.geo[idx] > 0)
				{
					map[cell_count] = temp;		// reading into 1D array of size N
					if (map[cell_count] >= 0.0) map_count++;
					cell_count++;
				}
			}
		}
	}

	printf("Map file %s read || Nmap = %d, recheck of Ncells = %d\n", filein, map_count, cell_count);

	fclose(in);

	// Write map vtk to Outputs, with map type in filename (this is "ref" here)
	sprintf(string,"%s/Map_%s.vtk", Output_dir, ref);

	FILE *out;
	out = fopen(string, "wt");
	sprintf(string,"%s", ref);

	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "vtk output\n");
	fprintf(out, "ASCII\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
	fprintf(out, "SPACING 1 1 1\n");
	fprintf(out, "ORIGIN 0 0 0\n");
	fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
	//fprintf(out, "SCALARS ImageFile float 1\n");
	fprintf(out, "SCALARS %s float 1\n", string);
	fprintf(out, "LOOKUP_TABLE default\n");

	cell_count = 0;

	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
				if (sc.geo[idx] > 0)
				{
					fprintf(out, "%f ", map[cell_count]);
					cell_count++;
				}
				else fprintf(out, "-100 ");
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);

	return map_count;
}

// NETWORK
void read_map_file_Njunc(SC_variables *sc, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, double *map)
{
    FILE *in;
    char *string = (char*)malloc(500);
    double temp;

    // Assign filename to string
    sprintf(string, "%s/%s/%s", PATH, fileroot, filein);

    in = fopen(string, "r");

    if (in == NULL)
    {
        printf("Cannot load geometry map file %s\t :: is the path correct? Does the file exist in that path?\n", string);
        exit(1);
    }
    else printf("File loaded %s\n", string);

    for (int i = 0; i < sc->Njunc; i++)
    {
        fscanf(in, "%lf ", &temp);
        map[i] = temp;
    }

    fclose(in);

    // Write map vtk to Outputs, with map type in filename (this is "ref" here)
    sprintf(string,"%s/Junction_Map.dat", Output_dir);

    FILE *out;
    out = fopen(string, "wt");
    for (int i = 0; i < sc->Njunc; i++)
    {
        fprintf(out, "%f\n", map[i]);
    }
    fclose(out);

}

// End Read geometry/maps from file into arrays =================================================//|

// Set D1, D2, dx and local values of these, from Tissue or Cell defined parameters =============\\|
void set_D_dx_global(SC_variables *sc, double dx, double dy, double dz, double D1, double D_AR)
{
	// Assigns dx, dy, dz and D1, D2 in SC_variables struct (using which caluclations are performed)
	// from values set in Tissue or CRU structs (or else) as set by tissue or cell model
	sc->dx = dx;
	sc->dy = dy;
	sc->dz = dz;

	for (int n = 0; n < sc->N; n++)
	{
		sc->D1[n] = D1;
		sc->D2[n] = D1/D_AR;
	}
}

// Sets components of the diffusion tensor from D1 and D2 arrays and fibre orientation
void set_D_array_anisotropic(SC_variables *sc)
{
	double D1, D2;
	for (int n = 0; n < sc->N; n++)
	{
		D1 = sc->D1[n];
		D2 = sc->D2[n];

		// Principal directions
		sc->Dxx[n]         = D2 +  (D1 - D2) * sc->ox[n]*sc->ox[n];
		sc->Dyy[n]         = D2 +  (D1 - D2) * sc->oy[n]*sc->oy[n];
		sc->Dzz[n]         = D2 +  (D1 - D2) * sc->oz[n]*sc->oz[n];

		// Diagonals
		sc->Dxy[n]         = (D1 - D2) * sc->ox[n]*sc->oy[n];
		sc->Dxz[n]         = (D1 - D2) * sc->ox[n]*sc->oz[n];
		sc->Dyz[n]         = (D1 - D2) * sc->oy[n]*sc->oz[n];
        sc->Dxy[n]         = fabs(sc->Dxy[n]);
        sc->Dxz[n]         = fabs(sc->Dxz[n]);
        sc->Dyz[n]         = fabs(sc->Dyz[n]);
	}
}

// print vtk files of D1 and D2 arrays (with all spatial het incorporated)
void output_D1_and_D2(SC_variables sc, const char* Output_dir)
{
	char *string = (char*)malloc(500);
	FILE *out;
	int idx;
	int cell_count;

	sprintf(string,"%s/Map_D1.vtk", Output_dir);
	out = fopen(string, "wt");

	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "vtk output\n");
	fprintf(out, "ASCII\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
	fprintf(out, "SPACING 1 1 1\n");
	fprintf(out, "ORIGIN 0 0 0\n");
	fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
	fprintf(out, "SCALARS D1 float 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");

	cell_count = 0;
	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
				if (sc.geo[idx] > 0)
				{
					fprintf(out, "%f ", sc.D1[cell_count]);
					cell_count++;
				}
				else fprintf(out, "-100 ");
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);

	sprintf(string,"%s/Map_D2.vtk", Output_dir);
	out = fopen(string, "wt");

	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "vtk output\n");
	fprintf(out, "ASCII\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
	fprintf(out, "SPACING 1 1 1\n");
	fprintf(out, "ORIGIN 0 0 0\n");
	fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
	fprintf(out, "SCALARS D2 float 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");

	cell_count = 0;
	for (int k = 0; k < sc.NZ; k++)
	{
		for (int j = 0; j < sc.NY; j++)
		{
			for (int i = 0; i < sc.NX; i++)
			{
				idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
				if (sc.geo[idx] > 0)
				{
					fprintf(out, "%f ", sc.D2[cell_count]);
					cell_count++;
				}
				else fprintf(out, "-100 ");
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);
}
// End Set D1, D2, dx and local values ==========================================================//|

// Set cell index and neigbours==================================================================\\|
void SC_set_index_and_geo_linear(SC_variables *sc)
{
	int count = 0; 	// counter of number of cells
	int idx;		// 3D identifier for each cell

	for (int k = 0; k < sc->NZ; k++)
	{
		for (int j = 0; j < sc->NY; j++)
		{
			for (int i = 0; i < sc->NX; i++)
			{
				idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);	// 3D identifier for each cell
				if (sc->geo[idx] > 0) // if it is an actual cell/node
				{
					sc->geo_index[idx] 		= count;		// gives the 1D>Ncell cell count at 3D idx
					sc->geo_3D_index[count] = idx;			// gives 3D index for 2D index(count)
					sc->x_index[count]		= i;
					sc->y_index[count]		= j;
					sc->z_index[count]		= k;
					sc->geo_linear[count] 	= sc->geo[idx]; 	// Set linear cell to 3D cell
					count++;
				}
			}
		}
	}
}

// Sets the neighbours for each cell; implements BCs by setting empty space neighbours to itself
void SC_set_neighbours(SC_variables *sc)
{
	int count = 0;  // counter of number of cells
	int idx;		// 3D identifier for each cell
	int idx_xm;		// 3D identifier for x-1 cell
	int idx_xp;		// 3D identifier for x+1 cell
	int idx_ym;		// 3D identifier for y-1 cell
	int idx_yp;		// 3D identifier for y+1 cell
	int idx_zm;		// 3D identifier for z-1 cell
	int idx_zp;		// 3D identifier for z+1 cell
	int idx_xm_ym;	// 3D identifier for x-1, y-1 cell
	int idx_xm_yp;
	int idx_xm_zm;
	int idx_xm_zp;
	int idx_xp_ym;	
	int idx_xp_yp;
	int idx_xp_zm;
	int idx_xp_zp;
	int idx_ym_zm;
	int idx_ym_zp;
	int idx_yp_zm;
	int idx_yp_zp;
	int idx_xm_ym_zm;
	int idx_xm_ym_zp;
	int idx_xm_yp_zm;
	int idx_xm_yp_zp;
	int idx_xp_ym_zm;
	int idx_xp_ym_zp;
	int idx_xp_yp_zm;
	int idx_xp_yp_zp;

	for (int k = 0; k < sc->NZ; k++)
	{
		for (int j = 0; j < sc->NY; j++)
		{
			for (int i = 0; i < sc->NX; i++)
			{
				idx 		=  i 	+ (sc->NX *  j)		+ (sc->NX * sc->NY *  k   ); 
				idx_xm		= (i-1) + (sc->NX *  j)		+ (sc->NX * sc->NY *  k   );
				idx_xp		= (i+1) + (sc->NX *  j)		+ (sc->NX * sc->NY *  k   );
				idx_ym		=  i	+ (sc->NX * (j-1))	+ (sc->NX * sc->NY *  k   );
				idx_yp		=  i	+ (sc->NX * (j+1))	+ (sc->NX * sc->NY *  k   );
				idx_zm		=  i	+ (sc->NX *  j)		+ (sc->NX * sc->NY * (k-1));
				idx_zp		=  i	+ (sc->NX *  j)		+ (sc->NX * sc->NY * (k+1));

				idx_xm_ym	= (i-1) + (sc->NX * (j-1)) 	+ (sc->NX * sc->NY *  k   );
				idx_xm_yp	= (i-1) + (sc->NX * (j+1)) 	+ (sc->NX * sc->NY *  k   );
				idx_xm_zm	= (i-1) + (sc->NX *  j) 	+ (sc->NX * sc->NY * (k-1));
				idx_xm_zp	= (i-1) + (sc->NX *  j) 	+ (sc->NX * sc->NY * (k+1));

				idx_xp_ym	= (i+1) + (sc->NX * (j-1)) 	+ (sc->NX * sc->NY *  k   );
				idx_xp_yp	= (i+1) + (sc->NX * (j+1)) 	+ (sc->NX * sc->NY *  k   );
				idx_xp_zm	= (i+1) + (sc->NX *  j) 	+ (sc->NX * sc->NY * (k-1));
				idx_xp_zp	= (i+1) + (sc->NX *  j) 	+ (sc->NX * sc->NY * (k+1));

				idx_ym_zm	=  i 	+ (sc->NX * (j-1)) 	+ (sc->NX * sc->NY * (k-1));
				idx_ym_zp	=  i	+ (sc->NX * (j-1)) 	+ (sc->NX * sc->NY * (k+1));
				idx_yp_zm	=  i	+ (sc->NX * (j+1)) 	+ (sc->NX * sc->NY * (k-1));
				idx_yp_zp	=  i	+ (sc->NX * (j+1)) 	+ (sc->NX * sc->NY * (k+1));

				idx_xm_ym_zm = (i-1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k-1));
				idx_xm_ym_zp = (i-1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
				idx_xm_yp_zm = (i-1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k-1));
				idx_xm_yp_zp = (i-1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));
				idx_xp_ym_zm = (i+1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k-1));
				idx_xp_ym_zp = (i+1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
				idx_xp_yp_zm = (i+1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k-1));
				idx_xp_yp_zp = (i+1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));

				if (sc->geo[idx] > 0) // if it is an actual cell/node
				{
					// Principal directions =======================\\|
					// x direction ======================\\|
					// x-1 (xm) 
					if (i == 0)						sc->xm[count] = count;					// if small x edge, xminus returns itself
					else if (sc->geo[idx_xm] < 1)	sc->xm[count] = count;					// if xminus neighbour is empty space, also return itself
					else							sc->xm[count] = sc->geo_index[idx_xm];	// xminus is a cell; return linear index of cell at idx x-1 

					// x+1 (xp) 
					if (i == sc->NX-1)              sc->xp[count] = count;                  // if large x edge, xplus returns itself
					else if (sc->geo[idx_xp] < 1)   sc->xp[count] = count;                  // if xplus neighbour is empty space, also return itself
					else                            sc->xp[count] = sc->geo_index[idx_xp];  // xplus is a cell; return linear index of cell at idx x+1
					// End x direction ==================//|

					// y direction ======================\\|
					// y-1 (ym) 
					if (j == 0)                     sc->ym[count] = count;                  // if small y edge, yminus returns itself
					else if (sc->geo[idx_ym] < 1)   sc->ym[count] = count;                  // if yminus neighbour is empty space, also return itself
					else                            sc->ym[count] = sc->geo_index[idx_ym];  // yminus is a cell; return linear index of cell at idx y-1 

					// y+1 (yp) 
					if (j == sc->NY-1)              sc->yp[count] = count;                  // if large y edge, yplus returns itself
					else if (sc->geo[idx_yp] < 1)   sc->yp[count] = count;                  // if yplus neighbour is empty space, also return itself
					else                            sc->yp[count] = sc->geo_index[idx_yp];  // yplus is a cell; return linear index of cell at idx y+1
					// End y direction ==================//|

					// z direction ======================\\|
					// z-1 (zm) 
					if (k == 0)                     sc->zm[count] = count;                  // if small z edge, zminus returns itself
					else if (sc->geo[idx_zm] < 1)   sc->zm[count] = count;                  // if zminus neighbour is empty space, also return itself
					else                            sc->zm[count] = sc->geo_index[idx_zm];  // zminus is a cell; return linear index of cell at idx z-1 

					// z+1 (zp) 
					if (k == sc->NZ-1)              sc->zp[count] = count;                  // if large z edge, zplus returns itself
					else if (sc->geo[idx_zp] < 1)   sc->zp[count] = count;                  // if zplus neighbour is empty space, also return itself
					else                            sc->zp[count] = sc->geo_index[idx_zp];  // zplus is a cell; return linear index of cell at idx z+1
					// End z direction ==================//|
					// end Principal directions ===================//|

					// Diagonals ==================================\\|
					// xm directions ====================\\|
					// x-1, y-1
					if (i == 0 || j == 0)            	sc->xm_ym[count] = count;           // if either x or y is lowest bound, x-1 y-1 cannot exist
					else if (sc->geo[idx_xm_ym] < 1)   	sc->xm_ym[count] = count;                  
					else                            	sc->xm_ym[count] = sc->geo_index[idx_xm_ym];  

					// x-1, y+1
					if (i == 0 || j == sc->NY-1)       	sc->xm_yp[count] = count;           // if either x is lowest or y is highest, cannot exist  
					else if (sc->geo[idx_xm_yp] < 1)   	sc->xm_yp[count] = count;                  
					else                            	sc->xm_yp[count] = sc->geo_index[idx_xm_yp]; 

					// x-1, z-1
					if (i == 0 || k == 0)               sc->xm_zm[count] = count;           // if either x or z is lowest bound, x-1 z-1 cannot exist
					else if (sc->geo[idx_xm_zm] < 1)    sc->xm_zm[count] = count;
					else                                sc->xm_zm[count] = sc->geo_index[idx_xm_zm];

					// x-1, z+1
					if (i == 0 || k == sc->NZ-1)        sc->xm_zp[count] = count;           // if either x is lowest or z is highest, cannot exist  
					else if (sc->geo[idx_xm_zp] < 1)    sc->xm_zp[count] = count;
					else                                sc->xm_zp[count] = sc->geo_index[idx_xm_zp]; 
					// End xm directions ================//|

					// xp directions ====================\\|
					// x+1, y-1
					if (i == sc->NX-1 || j == 0)       	sc->xp_ym[count] = count;           
					else if (sc->geo[idx_xp_ym] < 1)   	sc->xp_ym[count] = count;                  
					else                            	sc->xp_ym[count] = sc->geo_index[idx_xp_ym];  

					// x+1, y+1
					if (i == sc->NX-1 || j == sc->NY-1)	sc->xp_yp[count] = count;            
					else if (sc->geo[idx_xp_yp] < 1)   	sc->xp_yp[count] = count;                  
					else                            	sc->xp_yp[count] = sc->geo_index[idx_xp_yp]; 

					// x+1, z-1
					if (i == sc->NX-1 || k == 0)        sc->xp_zm[count] = count;           
					else if (sc->geo[idx_xp_zm] < 1)    sc->xp_zm[count] = count;
					else                                sc->xp_zm[count] = sc->geo_index[idx_xp_zm];

					// x+1, z+1
					if (i == sc->NX-1 || k == sc->NZ-1) sc->xp_zp[count] = count;          
					else if (sc->geo[idx_xp_zp] < 1)    sc->xp_zp[count] = count;
					else                                sc->xp_zp[count] = sc->geo_index[idx_xp_zp]; 
					// End xp directions ================//|

					// ym, z directions =================\\|
					// y-1, z-1
					if (j == 0 || k == 0)               sc->ym_zm[count] = count;           
					else if (sc->geo[idx_ym_zm] < 1)    sc->ym_zm[count] = count;
					else                                sc->ym_zm[count] = sc->geo_index[idx_ym_zm];

					// y-1, z+1
					if (j == 0 || k == sc->NZ-1)        sc->ym_zp[count] = count;          
					else if (sc->geo[idx_ym_zp] < 1)    sc->ym_zp[count] = count;
					else                                sc->ym_zp[count] = sc->geo_index[idx_ym_zp];
					// End ym, z directions =============//|

					// yp, z directions =================\\|
					// y+1, z-1
					if (j == sc->NY-1 || k == 0)        sc->yp_zm[count] = count;           
					else if (sc->geo[idx_yp_zm] < 1)    sc->yp_zm[count] = count;
					else                                sc->yp_zm[count] = sc->geo_index[idx_yp_zm];

					// y+1, z+1
					if (j == sc->NY-1 || k == sc->NZ-1) sc->yp_zp[count] = count;           
					else if (sc->geo[idx_yp_zp] < 1)    sc->yp_zp[count] = count;
					else                                sc->yp_zp[count] = sc->geo_index[idx_yp_zp];
					// End ym, z directions =============//|
					// End Diagonals ==============================//|

					// Corners ====================================\\|
					// xm directions ====================\\|
					// xm ym ==================\\|
					// x-1, y-1, z-1
					if (i == 0 || j == 0 || k == 0)     				sc->xm_ym_zm[count] = count;      
					else if (sc->geo[idx_xm_ym_zm] < 1) 				sc->xm_ym_zm[count] = count;
					else                                				sc->xm_ym_zm[count] = sc->geo_index[idx_xm_ym_zm];

					// x-1, y-1, z+1
					if (i == 0 || j == 0 || k == sc->NZ-1) 				sc->xm_ym_zp[count] = count;      
					else if (sc->geo[idx_xm_ym_zp] < 1) 				sc->xm_ym_zp[count] = count;
					else                                				sc->xm_ym_zp[count] = sc->geo_index[idx_xm_ym_zp];
					// end xm ym ==============//|

					// xm yp ==================\\|
					// x-1, y+1, z-1
					if (i == 0 || j == sc->NY-1 || k == 0)  			sc->xm_yp_zm[count] = count;      
					else if (sc->geo[idx_xm_yp_zm] < 1) 				sc->xm_yp_zm[count] = count;
					else                                				sc->xm_yp_zm[count] = sc->geo_index[idx_xm_yp_zm];

					// x-1, y+1, z+1
					if (i == 0 || j == sc->NY-1 || k == sc->NZ-1) 		sc->xm_yp_zp[count] = count;      
					else if (sc->geo[idx_xm_yp_zp] < 1) 				sc->xm_yp_zp[count] = count;
					else                                				sc->xm_yp_zp[count] = sc->geo_index[idx_xm_yp_zp];
					// end xm yp ==============//|

					// xp directions ====================\\|
					// xp ym ==================\\|
					// x+1, y-1, z-1
					if (i == sc->NX-1 || j == 0 || k == 0)  			sc->xp_ym_zm[count] = count;      
					else if (sc->geo[idx_xp_ym_zm] < 1) 				sc->xp_ym_zm[count] = count;
					else                                				sc->xp_ym_zm[count] = sc->geo_index[idx_xp_ym_zm];

					// x+1, y-1, z+1
					if (i == sc->NX-1 || j == 0 || k == sc->NZ-1) 		sc->xp_ym_zp[count] = count;      
					else if (sc->geo[idx_xp_ym_zp] < 1) 				sc->xp_ym_zp[count] = count;
					else                                				sc->xp_ym_zp[count] = sc->geo_index[idx_xp_ym_zp];
					// end xp ym ==============//|

					// xp yp ==================\\|
					// x+1, y+1, z-1
					if (i == sc->NX-1 || j == sc->NY-1 || k == 0)  		sc->xp_yp_zm[count] = count;      
					else if (sc->geo[idx_xp_yp_zm] < 1) 				sc->xp_yp_zm[count] = count;
					else                                				sc->xp_yp_zm[count] = sc->geo_index[idx_xp_yp_zm];

					// x+1, y+1, z+1
					if (i == sc->NX-1 || j == sc->NY-1 || k == sc->NZ-1)sc->xp_yp_zp[count] = count;      
					else if (sc->geo[idx_xp_yp_zp] < 1) 				sc->xp_yp_zp[count] = count;
					else                                				sc->xp_yp_zp[count] = sc->geo_index[idx_xp_yp_zp];
					// end xp yp ==============//|
					// End Corners ================================//|

					count++;
				}	
			}
		}
	}
}
// End set cell index and neigbours==============================================================//|

// Finite difference method =====================================================================\\|
// simple, isotropic
void calc_diff_FDM_iso(SC_variables *sc, double *v, int n)
{
	//Note: in 1D, ym = yp = zm = zp = return self; in 2D, zm and zp return self
	// (v_x-1 + v_x+1 - 2*v_x) / (dx^2)
	// Recall that xm[n] = index of cell x -1 to cell ; v[xm[n]] = the variable value in that cell
	sc->diff[n]    =   sc->D[n] * ( (v[sc->xm[n]] + v[sc->xp[n]] - 2*v[n])/(sc->dx*sc->dx) );
	sc->diff[n]    +=  sc->D[n] * ( (v[sc->ym[n]] + v[sc->yp[n]] - 2*v[n])/(sc->dy*sc->dy) );
	sc->diff[n]    +=  sc->D[n] * ( (v[sc->zm[n]] + v[sc->zp[n]] - 2*v[n])/(sc->dz*sc->dz) );
}

// simple, isotropic, tau method
void calc_diff_FDM_tau(SC_variables *sc, double *v, int n, double tau_trans, double tau_long)
{
	sc->diff[n]    =   (v[sc->xm[n]] + v[sc->xp[n]] - 2*v[n])/tau_trans;
	sc->diff[n]    +=  (v[sc->ym[n]] + v[sc->yp[n]] - 2*v[n])/tau_trans;
	sc->diff[n]    +=  (v[sc->zm[n]] + v[sc->zp[n]] - 2*v[n])/tau_long;
}

// Anisotropic, non-uniform
// D differential
void calc_dD_anisotropic_3D(SC_variables *sc, int n)
{   
	// Recall that xm[n] = index of cell x -1 to cell ; D[xm[n]] = the variable value in that cell

	// First, dDxx/dx, dDyy/dy and dDzz/dz (exaclty same as iso + het D above)
	// dDxx/dx
	if (sc->xp[n] == n || sc->xm[n] == n)  // don't need the "and" clause, as if both are non-tissue, then both return itself and this naturally goes to zero
		sc->dDxx_dx[n] 	= (sc->Dxx[sc->xp[n]] - sc->Dxx[sc->xm[n]])/(sc->dx);
	else    sc->dDxx_dx[n]	= (sc->Dxx[sc->xp[n]] - sc->Dxx[sc->xm[n]])/(2*sc->dx);

	// dDyy/dy
	if (sc->yp[n] == n || sc->ym[n] == n)
		sc->dDyy_dy[n] 	= (sc->Dyy[sc->yp[n]] - sc->Dyy[sc->ym[n]])/(sc->dy);  
	else    sc->dDyy_dy[n]  = (sc->Dyy[sc->yp[n]] - sc->Dyy[sc->ym[n]])/(2*sc->dy);

	// dDzz/dz
	if (sc->zp[n] == n || sc->zm[n] == n)
			sc->dDzz_dz[n] 	= (sc->Dzz[sc->zp[n]] - sc->Dzz[sc->zm[n]])/(sc->dz);  
	else    sc->dDzz_dz[n]	= (sc->Dzz[sc->zp[n]] - sc->Dzz[sc->zm[n]])/(2*sc->dz);

	// Now the mixed terms:
	// dDxy/dx
	if (sc->xp[n] == n || sc->xm[n] == n) // still just x, as gradient is in x direction
			sc->dDxy_dx[n]	= (sc->Dxy[sc->xp[n]] - sc->Dxy[sc->xm[n]])/(sc->dx);
	else 	sc->dDxy_dx[n]  = (sc->Dxy[sc->xp[n]] - sc->Dxy[sc->xm[n]])/(2*sc->dx);

	// dDxz/dx
	if (sc->xp[n] == n || sc->xm[n] == n) // still just x, as gradient is in x direction
			sc->dDxz_dx[n]  = (sc->Dxz[sc->xp[n]] - sc->Dxz[sc->xm[n]])/(sc->dx);
	else    sc->dDxz_dx[n]  = (sc->Dxz[sc->xp[n]] - sc->Dxz[sc->xm[n]])/(2*sc->dx);

	// dDyx/dy (= dDxy/dy)
	if (sc->yp[n] == n || sc->ym[n] == n) // still just y, as gradient is in y direction
			sc->dDxy_dy[n]  = (sc->Dxy[sc->yp[n]] - sc->Dxy[sc->ym[n]])/(sc->dy);
	else    sc->dDxy_dy[n]  = (sc->Dxy[sc->yp[n]] - sc->Dxy[sc->ym[n]])/(2*sc->dy);

	// dDyz/dy
	if (sc->yp[n] == n || sc->ym[n] == n) // still just y, as gradient is in y direction
			sc->dDyz_dy[n]  = (sc->Dyz[sc->yp[n]] - sc->Dyz[sc->ym[n]])/(sc->dy);
	else    sc->dDyz_dy[n]  = (sc->Dyz[sc->yp[n]] - sc->Dyz[sc->ym[n]])/(2*sc->dy);

	// dDzx/dz (=dDxz/dz)
	if (sc->zp[n] == n || sc->zm[n] == n) // still just z, as gradient is in z direction
			sc->dDxz_dz[n]  = (sc->Dxz[sc->zp[n]] - sc->Dxz[sc->zm[n]])/(sc->dz);
	else    sc->dDxz_dz[n]  = (sc->Dxz[sc->zp[n]] - sc->Dxz[sc->zm[n]])/(2*sc->dz);

	// dDzy/dz (=dDyz/dz)
	if (sc->zp[n] == n || sc->zm[n] == n) // still just z, as gradient is in z direction
			sc->dDyz_dz[n]  = (sc->Dyz[sc->zp[n]] - sc->Dyz[sc->zm[n]])/(sc->dz);
	else    sc->dDyz_dz[n]  = (sc->Dyz[sc->zp[n]] - sc->Dyz[sc->zm[n]])/(2*sc->dz);
}  

// update dvdt / calc diff
void calc_diff_FDM_anisotropic(SC_variables *sc, double *v, int n)
{
    sc->diff[n] =   (sc->Dxx[n] * ((v[sc->xm[n]] + v[sc->xp[n]] - 2*v[n])/(sc->dx*sc->dx)) ); // d2udx 
    sc->diff[n] +=  (sc->Dyy[n] * ((v[sc->ym[n]] + v[sc->yp[n]] - 2*v[n])/(sc->dy*sc->dy)) ); // y
    sc->diff[n] +=  (sc->Dzz[n] * ((v[sc->zm[n]] + v[sc->zp[n]] - 2*v[n])/(sc->dz*sc->dz)) ); // z

    sc->diff[n] +=  (((v[sc->xp[n]] - v[sc->xm[n]])/(2*sc->dx)) * (sc->dDxx_dx[n] + sc->dDxy_dy[n] + sc->dDxz_dz[n]));   // dudx
    sc->diff[n] +=  (((v[sc->yp[n]] - v[sc->ym[n]])/(2*sc->dy)) * (sc->dDyy_dy[n] + sc->dDxy_dx[n] + sc->dDyz_dz[n]));   // y
    sc->diff[n] +=  (((v[sc->zp[n]] - v[sc->zm[n]])/(2*sc->dz)) * (sc->dDzz_dz[n] + sc->dDxz_dx[n] + sc->dDyz_dy[n]));   // z

    sc->diff[n] += (2*sc->Dxy[n]*(v[sc->xp_yp[n]] + v[sc->xm_ym[n]] - v[sc->xp_ym[n]] - v[sc->xm_yp[n]])/(4*sc->dx*sc->dy) ); //d2udxy
    sc->diff[n] += (2*sc->Dxz[n]*(v[sc->xp_zp[n]] + v[sc->xm_zm[n]] - v[sc->xp_zm[n]] - v[sc->xm_zp[n]])/(4*sc->dx*sc->dz) );
    sc->diff[n] += (2*sc->Dyz[n]*(v[sc->yp_zp[n]] + v[sc->ym_zm[n]] - v[sc->yp_zm[n]] - v[sc->ym_zp[n]])/(4*sc->dy*sc->dz) );
}

// Create laplacian and implement BCs
void calc_laplacian_and_BCs(SC_variables *sc, int n)
{
	// This function computes the contribution of each neighbour to the voltage differential
	// In implementation, the factor will multiply the neighbour's voltage to determine effect on present cell	
	double factor;

	// First, default laplacian to 0
	sc->lap_self[n] 	= 0;
	sc->lap_xm[n] 		= 0;
	sc->lap_xp[n] 		= 0;
	sc->lap_ym[n] 		= 0;
	sc->lap_yp[n] 		= 0;
	sc->lap_zm[n] 		= 0;
	sc->lap_zp[n] 		= 0;
	sc->lap_xm_ym[n] 	= 0;
	sc->lap_xm_yp[n] 	= 0;
	sc->lap_xp_ym[n]	= 0;
	sc->lap_xp_yp[n]	= 0;
	sc->lap_xm_zm[n]	= 0;
	sc->lap_xm_zp[n]	= 0;
	sc->lap_xp_zm[n]	= 0;
	sc->lap_xp_zp[n]	= 0;
	sc->lap_ym_zm[n]	= 0;
	sc->lap_ym_zp[n]	= 0;
	sc->lap_yp_zm[n]	= 0;
	sc->lap_yp_zp[n]	= 0;

	// dudx
	factor = (1.0 / 2.0) * (1 / sc->dx) * (sc->dDxx_dx[n]+sc->dDxy_dy[n]+sc->dDxz_dz[n]); // this is what multiplies (V_xp - V_xm)

	if (sc->xp[n] != n && sc->xm[n] != n) // this means if both neighbours are real tissue, as returns self if not real tissue
	{
		// implementing upwind scheme
		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_xp[n] 		+= 2.0 * factor;
			sc->lap_xm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_xp[n] 		+= 0.0 * factor;
			sc->lap_xm[n] 		+= -2.0 * factor;
		}

	} else if(sc->xp[n] != n && sc->xm[n] == n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_xp[n] 		+= 2.0 * factor;
			sc->lap_xm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_xp[n] 		+= 0.0 * factor;
			sc->lap_xm[n] 		+= 0.0 * factor;
		}

	} else if(sc->xp[n] == n && sc->xm[n] != n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_xp[n] 		+= 0.0 * factor;
			sc->lap_xm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_xp[n] 		+= 0.0 * factor;
			sc->lap_xm[n] 		+= -2.0 * factor;
		}

	} else { // is this else necessary? already set to 0!
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_xp[n] 			+= 0.0 * factor;
		sc->lap_xm[n] 			+= 0.0 * factor;
	}

	// dudy
	factor = (1.0 / 2.0) * (1 / sc->dy) * (sc->dDxy_dx[n]+sc->dDyy_dy[n]+sc->dDyz_dz[n]);

	if(sc->yp[n] != n && sc->ym[n] != n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_yp[n] 		+= 2.0 * factor;
			sc->lap_ym[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_yp[n] 		+= 0.0 * factor;
			sc->lap_ym[n] 		+= -2.0 * factor;
		}

	} else if(sc->yp[n] != n && sc->ym[n] == n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_yp[n] 		+= 2.0 * factor;
			sc->lap_ym[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_yp[n] 		+= 0.0 * factor;
			sc->lap_ym[n] 		+= 0.0 * factor;
		}

	} else if(sc->yp[n] == n && sc->ym[n] != n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_yp[n] 		+= 0.0 * factor;
			sc->lap_ym[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_yp[n] 		+= 0.0 * factor;
			sc->lap_ym[n] 		+= -2.0 * factor;
		}

	} else {
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_yp[n] 			+= 0.0 * factor;
		sc->lap_ym[n] 			+= 0.0 * factor;
	}


	// dudz
	factor = (1.0 / 2.0) * (1 / sc->dz) * (sc->dDxz_dx[n]+sc->dDyz_dy[n]+sc->dDzz_dz[n]);

	if(sc->zp[n] != n && sc->zm[n] != n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_zp[n] 		+= 2.0 * factor;
			sc->lap_zm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_zp[n] 		+= 0.0 * factor;
			sc->lap_zm[n] 		+= -2.0 * factor;
		}

	} else if(sc->zp[n] != n && sc->zm[n] == n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= -2.0 * factor;
			sc->lap_zp[n] 		+= 2.0 * factor;
			sc->lap_zm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_zp[n] 		+= 0.0 * factor;
			sc->lap_zm[n] 		+= 0.0 * factor;
		}

	} else if(sc->zp[n] == n && sc->zm[n] != n){

		if(factor>=0.0){
			sc->lap_self[n] 	+= 0.0 * factor;
			sc->lap_zp[n] 		+= 0.0 * factor;
			sc->lap_zm[n] 		+= 0.0 * factor;
		} else {
			sc->lap_self[n] 	+= 2.0 * factor;
			sc->lap_zp[n] 		+= 0.0 * factor;
			sc->lap_zm[n] 		+= -2.0 * factor;
		}

	} else {
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_zp[n] 			+= 0.0 * factor;
		sc->lap_zm[n] 			+= 0.0 * factor;
	}

	// dudx2
	factor = (1 / sc->dx) * (1 / sc->dx) * (sc->Dxx[n]);

	if(sc->xp[n] != n && sc->xm[n] != n){
		sc->lap_self[n] 		+= -2.0 * factor;
		sc->lap_xp[n] 			+= 1.0 * factor;
		sc->lap_xm[n] 			+= 1.0 * factor;
	} else if(sc->xp[n] != n && sc->xm[n] == n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_xp[n] 			+= 1.0 * factor;
		sc->lap_xm[n] 			+= 0.0 * factor;
	} else if(sc->xp[n] == n && sc->xm[n] != n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_xp[n] 			+= 0.0 * factor;
		sc->lap_xm[n] 			+= 1.0 * factor;
	} else {
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_xp[n] 			+= 0.0 * factor;
		sc->lap_xm[n] 			+= 0.0 * factor;
	}

	// dudy2
	factor = (1 / sc->dy) * (1 / sc->dy) * (sc->Dyy[n]);

	if(sc->yp[n] != n && sc->ym[n] != n){
		sc->lap_self[n] 		+= -2.0 * factor;
		sc->lap_yp[n] 			+= 1.0 * factor;
		sc->lap_ym[n] 			+= 1.0 * factor;
	} else if(sc->yp[n] != n && sc->ym[n] == n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_yp[n] 			+= 1.0 * factor;
		sc->lap_ym[n] 			+= 0.0 * factor;
	} else if(sc->yp[n] == n && sc->ym[n] != n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_yp[n] 			+= 0.0 * factor;
		sc->lap_ym[n] 			+= 1.0 * factor;
	} else {
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_yp[n] 			+= 0.0 * factor;
		sc->lap_ym[n] 			+= 0.0 * factor;
	}

	// dudz2
	factor = (1 / sc->dz) * (1 / sc->dz) * (sc->Dzz[n]);

	if(sc->zp[n] != n && sc->zm[n] != n){
		sc->lap_self[n] 		+= -2.0 * factor;
		sc->lap_zp[n] 			+= 1.0 * factor;
		sc->lap_zm[n] 			+= 1.0 * factor;
	} else if(sc->zp[n] != n && sc->zm[n] == n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_zp[n] 			+= 1.0 * factor;
		sc->lap_zm[n] 			+= 0.0 * factor;
	} else if(sc->zp[n] == n && sc->zm[n] != n){
		sc->lap_self[n] 		+= -1.0 * factor;
		sc->lap_zp[n] 			+= 0.0 * factor;
		sc->lap_zm[n] 			+= 1.0 * factor;
	} else {
		sc->lap_self[n] 		+= 0.0 * factor;
		sc->lap_zp[n] 			+= 0.0 * factor;
		sc->lap_zm[n] 			+= 0.0 * factor;
	}

	// dudxdy
	factor = (1.0 / 4.0) * (1 / sc->dx) * (1 / sc->dy) * 2.0 * (sc->Dxy[n]);

	if(sc->xp_yp[n] != n  && sc->xp_ym[n] != n && sc->xm_ym[n] != n && sc->xm_ym[n] != n) {
		sc->lap_xp_yp[n] 		+= 1.0 * factor;
		sc->lap_xm_ym[n] 		+= 1.0 * factor;
		sc->lap_xp_ym[n] 		+= -1.0 * factor;
		sc->lap_xm_yp[n] 		+= -1.0 * factor;
	} else {
		sc->lap_xp_yp[n] 		+= 0.0 * factor;
		sc->lap_xm_ym[n] 		+= 0.0 * factor;
		sc->lap_xp_ym[n] 		+= 0.0 * factor;
		sc->lap_xm_yp[n] 		+= 0.0 * factor;
	}

	// dudxdz
	factor = (1.0 / 4.0) * (1 / sc->dx) * (1 / sc->dz) * 2.0 * (sc->Dxz[n]);

	if(sc->xp_zp[n] != n && sc->xp_zm[n] != n && sc->xm_zp[n] != n && sc->xm_zm[n] != n) {
		sc->lap_xp_zp[n] 		+= 1.0 * factor;
		sc->lap_xm_zm[n] 		+= 1.0 * factor;
		sc->lap_xp_zm[n] 		+= -1.0 * factor;
		sc->lap_xm_zp[n] 		+= -1.0 * factor;
	} else {
		sc->lap_xp_zp[n] 		+= 0.0 * factor;
		sc->lap_xm_zm[n] 		+= 0.0 * factor;
		sc->lap_xp_zm[n] 		+= 0.0 * factor;
		sc->lap_xm_zp[n] 		+= 0.0 * factor;
	}

	// dudydz
	factor = (1.0 / 4.0) * (1 / sc->dy) * (1 / sc->dz) * 2.0 * (sc->Dyz[n]);

	if(sc->yp_zp[n] != n && sc->yp_zm[n] != n && sc->ym_zp[n] != n && sc->ym_zm[n] != n) {
		sc->lap_yp_zp[n] 		+= 1.0 * factor;
		sc->lap_ym_zm[n] 		+= 1.0 * factor;
		sc->lap_yp_zm[n] 		+= -1.0 * factor;
		sc->lap_ym_zp[n] 		+= -1.0 * factor;
	} else {
		sc->lap_yp_zp[n] 		+= 0.0 * factor;
		sc->lap_ym_zm[n] 		+= 0.0 * factor;
		sc->lap_yp_zm[n] 		+= 0.0 * factor;
		sc->lap_ym_zp[n] 		+= 0.0 * factor;
	}
}

void calc_diff_from_lap(SC_variables *sc, double *v, int n)
{
	sc->diff[n]		= v[n]*sc->lap_self[n];

	sc->diff[n]		+= v[sc->xm[n]]*sc->lap_xm[n];
	sc->diff[n]		+= v[sc->xp[n]]*sc->lap_xp[n];
	sc->diff[n]		+= v[sc->ym[n]]*sc->lap_ym[n];
	sc->diff[n]		+= v[sc->yp[n]]*sc->lap_yp[n];
	sc->diff[n]		+= v[sc->zm[n]]*sc->lap_zm[n];
	sc->diff[n]		+= v[sc->zp[n]]*sc->lap_zp[n];

	sc->diff[n]		+= v[sc->xm_ym[n]]*sc->lap_xm_ym[n];
	sc->diff[n]		+= v[sc->xm_yp[n]]*sc->lap_xm_yp[n];
	sc->diff[n]		+= v[sc->xp_ym[n]]*sc->lap_xp_ym[n];
	sc->diff[n]		+= v[sc->xp_yp[n]]*sc->lap_xp_yp[n];
	
	sc->diff[n]		+= v[sc->xm_zm[n]]*sc->lap_xm_zm[n];
	sc->diff[n]		+= v[sc->xm_zp[n]]*sc->lap_xm_zp[n];
	sc->diff[n]		+= v[sc->xp_zm[n]]*sc->lap_xp_zm[n];
	sc->diff[n]		+= v[sc->xp_zp[n]]*sc->lap_xp_zp[n];
	
	sc->diff[n]		+= v[sc->ym_zm[n]]*sc->lap_ym_zm[n];
	sc->diff[n]		+= v[sc->ym_zp[n]]*sc->lap_ym_zp[n];
	sc->diff[n]		+= v[sc->yp_zm[n]]*sc->lap_yp_zm[n];
	sc->diff[n]		+= v[sc->yp_zp[n]]*sc->lap_yp_zp[n];
}
// End alternative implementation
// End finite difference method =================================================================//|

// NETWORK model implementation =================================================================\\|
// set parameters of conduction and resolution
void set_G_dx_global(SC_variables *sc, double dx, double dy, double dz, double Gl, double Gt, double Gt2, double symm_fac_diag_long, double symm_fac_diag_trans, double symm_fac_corner)
{
    // Assigns dx, dy, dz and D1, D2 in SC_variables struct (using which caluclations are performed)
    // from values set in Tissue or CRU structs (or else) as set by tissue or cell model
    sc->dx = dx;
    sc->dy = dy;
    sc->dz = dz;

    for (int n = 0; n < sc->N; n++)
    {
        sc->Gl[n]   = Gl;
        sc->Gt[n]   = Gt;
        sc->Gt2[n]  = Gt2; 
    }

    // Factor in dx here??

    sc->symm_fac_diag_long      =  symm_fac_diag_long;
    sc->symm_fac_diag_trans     =  symm_fac_diag_trans;
    sc->symm_fac_corner         =  symm_fac_corner;
}

// Set gGgap node array (apply orientation weights) ==========================\\|
void set_gGgap_array(SC_variables *sc, char const* Orientation_type)
{
    double Gt, Gl, Gt2; // local copies of Gtransverse and Glong
    double theta_x_y;   // angle in x-y plane
    double theta_z_x;
    double theta_z_y;
    double W_axis_x_y; // Weight towards principal (axis; away from diagonal) in x-y plane
    double W_axis_x_z; // Weight towards principal (axis; away from diagonal) in x-z plane
    double W_axis_y_z; // Weight towards principal (axis; away from diagonal) in y-z plane
    double W_axis;     // Weight towards principal direction
    double W_diag_plane; // Diagonal, in plane
    double W_diag_elevation; // Diagonal, elevated
    double W_corner; // corner
    double W_xx, W_yy, W_zz;
    double W_xypp, W_xypm, W_xzpp, W_xzpm, W_yzpp, W_yzpm;
    double W_xyzppp, W_xyzppm, W_xyzpmp, W_xyzmpp;
    double x, y, z;
    double x_t, y_t, z_t;
    double x_t2, y_t2, z_t2;
    double dx, dy, dz;
    double theta, phi;
    double theta_t, phi_t;
    double theta_t2, phi_t2;

    // "Pointing" parameters (1 if pointing in said direction, 0 otherwise)
    double Px, Py, Pz;
    double Pxyp, Pxzp, Pyzp;
    double Px_t, Py_t, Pz_t;
    double Pxyp_t, Pxzp_t, Pyzp_t;
    double Px_t2, Py_t2, Pz_t2;
    double Pxyp_t2, Pxzp_t2, Pyzp_t2;

    // Weights for transverse connections
    double W_t_xx, W_t_yy, W_t_zz;
    double W_t_xypp, W_t_xypm, W_t_xzpp, W_t_xzpm, W_t_yzpp, W_t_yzpm;
    double W_t_xyzppp, W_t_xyzppm, W_t_xyzpmp, W_t_xyzmpp;
    double W_t2_xx, W_t2_yy, W_t2_zz;
    double W_t2_xypp, W_t2_xypm, W_t2_xzpp, W_t2_xzpm, W_t2_yzpp, W_t2_yzpm;
    double W_t2_xyzppp, W_t2_xyzppm, W_t2_xyzpmp, W_t2_xyzmpp;

    printf("Assigning gGgap according to network model\n");

    //FILE * debug_outputs;
    //debug_outputs = fopen("Debug_network.txt", "a");
    //FILE * debug_outputs2;
    //debug_outputs2 = fopen("t1.vtk", "wt");
    //FILE * debug_outputs3;
    //debug_outputs3 = fopen("t2.vtk", "wt");

    dx = sc->dx;
    dy = sc->dy;
    dz = sc->dz;
        
    // symmetry factors
    double symm_fac_diag_long           = sc->symm_fac_diag_long;
    double symm_fac_diag_trans          = sc->symm_fac_diag_trans;
    double symm_fac_corner              = sc->symm_fac_corner;
    printf("Symmetry factors are: long = %f trans = %f corner = %f\n", symm_fac_diag_long, symm_fac_diag_trans, symm_fac_corner);

    bool nofibreerrorflag = false;

    // Loop over all tissue ===================================================================================================\\|
    for (int n = 0; n < sc->N; n++)
    {
        // Assign local Gl and Gt
        Gl      = sc->Gl[n];
        Gt      = sc->Gt[n];
        Gt2     = sc->Gt[n];

        // Assign local fibre variables
        x = sc->ox[n];
        y = sc->oy[n];
        z = sc->oz[n];

        // Transverse orientation components 
        if (strcmp(Orientation_type, "three_eigenvectors") == 0)
        {
            // for this version gt and gt2 are different
            Gt2     = sc->Gt2[n]; 

            x_t = sc->ox2[n];
            y_t = sc->oy2[n];
            z_t = sc->oz2[n];

            x_t2 = sc->ox3[n];
            y_t2 = sc->oy3[n];
            z_t2 = sc->oz3[n];
        }
        else if (strcmp(Orientation_type, "anisotropic") == 0)  // Only primary is given so need to calculate the transverse 
        {
            if (z == 1) // special case as theta undefined. But easy as set to other axes
            {
                x_t     = 1;
                y_t     = 0;
                z_t     = 0;

                x_t2    = 0;
                y_t2    = 1;
                z_t2    = 0;
            }
            else 
            {
                // First, transform by 90 degree rotation in phi
                phi = asin(z);
                phi_t2 = phi + M_PI/2;
                x_t2 = (x/cos(phi))*cos(phi_t2);
                y_t2 = (y/cos(phi))*cos(phi_t2);
                z_t2 = sin(phi_t2);

                // Now, find normal to that vector and the original (from cross product of the 2)
                x_t = y*z_t2 - z*y_t2;
                y_t = z*x_t2 - x*z_t2;
                z_t = x*y_t2 - y*x_t2;
                //printf("x %f y %f z %f xt %f yt %f zt %f phi %f phi 2 %f x2 %f y2 %f z2 %f\n", x, y, z, x_t, y_t, z_t, phi, phi_t2, x_t2, y_t2, z_t2);
            }
        }
        else
        {
            printf("ERROR: Orientation type must be \"anisotropic\" or \"three_eigenvectors\" for the network model to be valid\n");
            exit(1);
        }

        // Default all weights to 0, so only need to add to relevant ones
        W_axis = W_diag_plane = W_diag_elevation = W_corner = 0;
        W_xx = W_yy = W_zz = 0;
        W_xypp = W_xypm = W_xzpp = W_xzpm = W_yzpp = W_yzpm = 0;
        W_xyzppp = W_xyzppm = W_xyzpmp = W_xyzmpp = 0;

        W_t_xx = W_t_yy = W_t_zz = 0;
        W_t_xypp = W_t_xypm = W_t_xzpp = W_t_xzpm = W_t_yzpp = W_t_yzpm = 0;
        W_t_xyzppp = W_t_xyzppm = W_t_xyzpmp = W_t_xyzmpp = 0;
        W_t2_xx = W_t2_yy = W_t2_zz = 0;
        W_t2_xypp = W_t2_xypm = W_t2_xzpp = W_t2_xzpm = W_t2_yzpp = W_t2_yzpm = 0;
        W_t2_xyzppp = W_t2_xyzppm = W_t2_xyzpmp = W_t2_xyzmpp = 0;

        // Default gGgap to zero so also can just be added where relevant
        sc->gGap_node_xx[n] = sc->gGap_node_yy[n] = sc->gGap_node_zz[n] = 0;
        sc->gGap_node_xypp[n] = sc->gGap_node_xypm[n] = sc->gGap_node_xzpp[n] = 0;
        sc->gGap_node_xzpm[n] = sc->gGap_node_yzpp[n] = sc->gGap_node_yzpm[n] = 0;
        sc->gGap_node_xyzppp[n] = sc->gGap_node_xyzppm[n] = sc->gGap_node_xyzpmp[n] = sc->gGap_node_xyzmpp[n] = 0;

        // Default all pointing parameters to 0
        Px = Py = Pz = Pxyp = Pxzp = Pyzp = 0;
        Px_t = Py_t = Pz_t = Pxyp_t = Pxzp_t = Pyzp_t = 0;
        Px_t2 = Py_t2 = Pz_t2 = Pxyp_t2 = Pxzp_t2 = Pyzp_t2 = 0;

        if (x*x + y*y + z*z > 0) // do fibres exist?
        {
            // Primary eigenvector ===========================================================\\|
            // First, calculate theta and weights in each plane
            theta_x_y       = asin(fabs(y)/(sqrt(x*x + y*y))); 
            theta_z_x       = asin(fabs(z)/(sqrt(z*z + x*x)));
            theta_z_y       = asin(fabs(z)/(sqrt(z*z + y*y)));
            if (x == 0 && y == 0)  theta_x_y = 0; // override any nans from a /0
            if (x == 0 && z == 0)  theta_z_x = 0;
            if (y == 0 && z == 0)  theta_z_y = 0;
            W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
            W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
            W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

            //printf("%f %f %f %f %f %f\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

            // Set pointing parameters
            // Main axes (only one can be non-zero at a time)
            if (x*x >= y*y && x*x >= z*z)   Px = 1;
            else if (y*y >= z*z)            Py = 1;
            else                            Pz = 1;

            // Diagonals (more than one can be non-zero)
            if (x*y >= 0)                   Pxyp = 1;
            if (x*z >= 0)                   Pxzp = 1;
            if (y*z >= 0)                   Pyzp = 1;

            // Assign quadrant weights
            W_axis              = Px*W_axis_x_y*W_axis_x_z          + Py*W_axis_x_y*W_axis_y_z          + Pz*W_axis_x_z*W_axis_y_z;
            W_diag_plane        = Px*(1-W_axis_x_y)*(W_axis_x_z)    + Py*(1-W_axis_x_y)*(W_axis_y_z)    + Pz*(1-W_axis_x_z)*W_axis_y_z;
            W_diag_elevation    = Px*W_axis_x_y*(1-W_axis_x_z)      + Py*W_axis_x_y*(1-W_axis_y_z)      + Pz*(W_axis_x_z)*(1-W_axis_y_z);
            W_corner            = Px*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz*(1-W_axis_x_z)*(1-W_axis_y_z);

            //printf("W %f %f %f %f\n", W_axis, W_diag_plane, W_diag_elevation, W_corner);

            // Assign weights to directions
            W_xx                = Px * W_axis;    // only non-zero if pointing primarily towards x
            W_yy                = Py * W_axis;
            W_zz                = Pz * W_axis;

            W_xypp              = (Px + Py) * Pxyp      * W_diag_plane;
            W_xypm              = (Px + Py) * (1-Pxyp)  * W_diag_plane;

            W_xzpp              = Px * Pxzp     * W_diag_elevation      + Pz * Pxzp     * W_diag_plane;
            W_xzpm              = Px * (1-Pxzp) * W_diag_elevation      + Pz * (1-Pxzp) * W_diag_plane;

            W_yzpp              = (Py + Pz) * Pyzp      * W_diag_elevation;
            W_yzpm              = (Py + Pz) * (1-Pyzp)  * W_diag_elevation;

            W_xyzppp            = Pxzp      *   Pyzp    * W_corner;
            W_xyzpmp            = (1-Pxyp)  *   Pxzp    * W_corner;
            W_xyzmpp            = (1-Pxyp)  *  (1-Pxzp) * W_corner;
            W_xyzppm            = Pxyp      *  (1-Pxzp) * W_corner;
            // End Primary eigenvector =======================================================//|

            // Transverse directions =========================================================\\|
            // transverse 1 =================================================\\|
            theta_x_y       = asin(fabs(y_t)/(sqrt(x_t*x_t + y_t*y_t)));    // again, put y on top
            theta_z_x       = asin(fabs(z_t)/(sqrt(z_t*z_t + x_t*x_t)));
            theta_z_y       = asin(fabs(z_t)/(sqrt(z_t*z_t + y_t*y_t)));
            if (x_t == 0 && y_t == 0)  theta_x_y = 0; // override any nans from a /0
            if (x_t == 0 && z_t == 0)  theta_z_x = 0;
            if (y_t == 0 && z_t == 0)  theta_z_y = 0;
            W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
            W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
            W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

            //printf("%f %f %f %f %f %f\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

            if (x_t*x_t >= y_t*y_t && x_t*x_t >= z_t*z_t)   Px_t = 1;
            else if (y_t*y_t >= z_t*z_t)                    Py_t = 1;
            else                                            Pz_t = 1;

            // Diagonals (more than one can be non-zero)
            if (x_t*y_t >= 0)                   Pxyp_t = 1;
            if (x_t*z_t >= 0)                   Pxzp_t = 1;
            if (y_t*z_t >= 0)                   Pyzp_t = 1;

            // Assign quadrant weights
            W_axis              = Px_t*W_axis_x_y*W_axis_x_z          + Py_t*W_axis_x_y*W_axis_y_z          + Pz_t*W_axis_x_z*W_axis_y_z;
            W_diag_plane        = Px_t*(1-W_axis_x_y)*(W_axis_x_z)    + Py_t*(1-W_axis_x_y)*(W_axis_y_z)    + Pz_t*(1-W_axis_x_z)*W_axis_y_z;
            W_diag_elevation    = Px_t*W_axis_x_y*(1-W_axis_x_z)      + Py_t*W_axis_x_y*(1-W_axis_y_z)      + Pz_t*(W_axis_x_z)*(1-W_axis_y_z);
            W_corner            = Px_t*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py_t*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz_t*(1-W_axis_x_z)*(1-W_axis_y_z);

            // Assign weights to directions
            W_t_xx                = Px_t * W_axis;    // only non-zero if pointing primarily towards x
            W_t_yy                = Py_t * W_axis;
            W_t_zz                = Pz_t * W_axis;

            W_t_xypp              = (Px_t + Py_t) * Pxyp_t      * W_diag_plane;
            W_t_xypm              = (Px_t + Py_t) * (1-Pxyp_t)  * W_diag_plane;

            W_t_xzpp              = Px_t * Pxzp_t     * W_diag_elevation      + Pz_t * Pxzp_t     * W_diag_plane;
            W_t_xzpm              = Px_t * (1-Pxzp_t) * W_diag_elevation      + Pz_t * (1-Pxzp_t) * W_diag_plane;

            W_t_yzpp              = (Py_t + Pz_t) * Pyzp_t      * W_diag_elevation;
            W_t_yzpm              = (Py_t + Pz_t) * (1-Pyzp_t)  * W_diag_elevation;

            W_t_xyzppp            = Pxzp_t      *   Pyzp_t    * W_corner;
            W_t_xyzpmp            = (1-Pxyp_t)  *   Pxzp_t    * W_corner;
            W_t_xyzmpp            = (1-Pxyp_t)  *  (1-Pxzp_t) * W_corner;
            W_t_xyzppm            = Pxyp_t      *  (1-Pxzp_t) * W_corner;
            // End transverse 1 =============================================//|

            // Transverse 2 =================================================\\|
            theta_x_y       = asin(fabs(y_t2)/(sqrt(x_t2*x_t2 + y_t2*y_t2)));
            theta_z_x       = asin(fabs(z_t2)/(sqrt(z_t2*z_t2 + x_t2*x_t2)));
            theta_z_y       = asin(fabs(z_t2)/(sqrt(z_t2*z_t2 + y_t2*y_t2)));
            if (x_t2 == 0 && y_t2 == 0)  theta_x_y = 0; // override any nans from a /0
            if (x_t2 == 0 && z_t2 == 0)  theta_z_x = 0;
            if (y_t2 == 0 && z_t2 == 0)  theta_z_y = 0;
            W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
            W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
            W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

            //printf("%f %f %f %f %f %f\n\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

            if (x_t2*x_t2 >= y_t2*y_t2 && x_t2*x_t2 >= z_t2*z_t2)   Px_t2 = 1;
            else if (y_t2*y_t2 >= z_t2*z_t2)                        Py_t2 = 1;
            else                                                    Pz_t2 = 1;

            // Diagonals (more than one can be non-zero)
            if (x_t2*y_t2 >= 0)                   Pxyp_t2 = 1;
            if (x_t2*z_t2 >= 0)                   Pxzp_t2 = 1;
            if (y_t2*z_t2 >= 0)                   Pyzp_t2 = 1;

            //printf("Points t %f %f %f\n", Px_t, Py_t, Pz_t);
            //printf("Points t2 %f %f %f\n", Px_t2, Py_t2, Pz_t2);

            // Assign quadrant weights
            W_axis              = Px_t2*W_axis_x_y*W_axis_x_z          + Py_t2*W_axis_x_y*W_axis_y_z          + Pz_t2*W_axis_x_z*W_axis_y_z;
            W_diag_plane        = Px_t2*(1-W_axis_x_y)*(W_axis_x_z)    + Py_t2*(1-W_axis_x_y)*(W_axis_y_z)    + Pz_t2*(1-W_axis_x_z)*W_axis_y_z;
            W_diag_elevation    = Px_t2*W_axis_x_y*(1-W_axis_x_z)      + Py_t2*W_axis_x_y*(1-W_axis_y_z)      + Pz_t2*(W_axis_x_z)*(1-W_axis_y_z);
            W_corner            = Px_t2*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py_t2*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz_t2*(1-W_axis_x_z)*(1-W_axis_y_z);

            //printf("W %f %f %f %f\n", W_axis, W_diag_plane, W_diag_elevation, W_corner);

            // Assign weights to directions
            W_t2_xx                = Px_t2 * W_axis;    // only non-zero if pointing primarily towards x
            W_t2_yy                = Py_t2 * W_axis;
            W_t2_zz                = Pz_t2 * W_axis;

            W_t2_xypp              = (Px_t2 + Py_t2) * Pxyp_t2      * W_diag_plane;
            W_t2_xypm              = (Px_t2 + Py_t2) * (1-Pxyp_t2)  * W_diag_plane;

            W_t2_xzpp              = Px_t2 * Pxzp_t2     * W_diag_elevation      + Pz_t2 * Pxzp_t2     * W_diag_plane;
            W_t2_xzpm              = Px_t2 * (1-Pxzp_t2) * W_diag_elevation      + Pz_t2 * (1-Pxzp_t2) * W_diag_plane;

            W_t2_yzpp              = (Py_t2 + Pz_t2) * Pyzp_t2      * W_diag_elevation;
            W_t2_yzpm              = (Py_t2 + Pz_t2) * (1-Pyzp_t2)  * W_diag_elevation;

            W_t2_xyzppp            = Pxzp_t2      *   Pyzp_t2    * W_corner;
            W_t2_xyzpmp            = (1-Pxyp_t2)  *   Pxzp_t2    * W_corner;
            W_t2_xyzmpp            = (1-Pxyp_t2)  *  (1-Pxzp_t2) * W_corner;
            W_t2_xyzppm            = Pxyp_t2      *  (1-Pxzp_t2) * W_corner;
            // End transverse 2 =============================================//|
            // End Transverse directions =====================================================//|
        } // end checking if fibre exists
        else // if no fibre
        {
            //printf("ERROR: No fibre found for node %d || model will work but is desgined for anisotropic simulations\n", n);
            // Assign fibre as being in all direcrions to ensure no reduced coupling
            W_xx = W_yy = W_zz = 1;
            W_xypp = W_xypm = W_xzpp = W_xzpm = W_yzpp = W_yzpm = 1;
            W_xyzppp = W_xyzppm = W_xyzpmp = W_xyzmpp = 1;
            nofibreerrorflag = true;
        }

        // Assign gGgap based on weights =================================================\\|
        sc->gGap_node_xx[n]        = Gt*W_t_xx     + Gt2*W_t2_xx     + Gl*W_xx; // simply = weight in transverse * gt + weight along fibre * Gl
        sc->gGap_node_yy[n]        = Gt*W_t_yy     + Gt2*W_t2_yy     + Gl*W_yy; // note that both weights can be 0, giving no contribution to said direction
        sc->gGap_node_zz[n]        = Gt*W_t_zz     + Gt2*W_t2_zz     + Gl*W_zz;

        sc->gGap_node_xypp[n]      = (Gt*W_t_xypp   + Gt2*W_t2_xypp)*symm_fac_diag_trans   + Gl*W_xypp * symm_fac_diag_long;
        sc->gGap_node_xypm[n]      = (Gt*W_t_xypm   + Gt2*W_t2_xypm)*symm_fac_diag_trans   + Gl*W_xypm * symm_fac_diag_long;
        sc->gGap_node_xzpp[n]      = (Gt*W_t_xzpp   + Gt2*W_t2_xzpp)*symm_fac_diag_trans   + Gl*W_xzpp * symm_fac_diag_long;
        sc->gGap_node_xzpm[n]      = (Gt*W_t_xzpm   + Gt2*W_t2_xzpm)*symm_fac_diag_trans   + Gl*W_xzpm * symm_fac_diag_long;
        sc->gGap_node_yzpp[n]      = (Gt*W_t_yzpp   + Gt2*W_t2_yzpp)*symm_fac_diag_trans   + Gl*W_yzpp * symm_fac_diag_long;
        sc->gGap_node_yzpm[n]      = (Gt*W_t_yzpm   + Gt2*W_t2_yzpm)*symm_fac_diag_trans   + Gl*W_yzpm * symm_fac_diag_long;

        sc->gGap_node_xyzppp[n]    = Gt*W_t_xyzppp + Gt2*W_t2_xyzppp + Gl*W_xyzppp * symm_fac_corner;
        sc->gGap_node_xyzppm[n]    = Gt*W_t_xyzppm + Gt2*W_t2_xyzppm + Gl*W_xyzppm * symm_fac_corner;
        sc->gGap_node_xyzpmp[n]    = Gt*W_t_xyzpmp + Gt2*W_t2_xyzpmp + Gl*W_xyzpmp * symm_fac_corner;
        sc->gGap_node_xyzmpp[n]    = Gt*W_t_xyzmpp + Gt2*W_t2_xyzmpp + Gl*W_xyzmpp * symm_fac_corner;
        // End Assign gGgap based on weights ============================================//|

        // Determine type of connection 0 = none, 1 = transverse, 2 = axial/long
        sc->connection_type_node_xx[n]                                                  = 0;
        if (W_xx > 0) sc->connection_type_node_xx[n]                                    = 2;
        else if (W_t_xx > 0 || W_t2_xx > 0) sc->connection_type_node_xx[n]              = 1;
        sc->connection_type_node_yy[n]                                                  = 0;
        if (W_yy > 0) sc->connection_type_node_yy[n]                                    = 2;
        else if (W_t_yy > 0 || W_t2_yy > 0) sc->connection_type_node_yy[n]              = 1;
        sc->connection_type_node_zz[n]                                                  = 0;
        if (W_zz > 0) sc->connection_type_node_zz[n]                                    = 2;
        else if (W_t_zz > 0 || W_t2_zz > 0) sc->connection_type_node_zz[n]              = 1;
        sc->connection_type_node_xypp[n]                                                = 0;
        if (W_xypp > 0) sc->connection_type_node_xypp[n]                                = 2;
        else if (W_t_xypp > 0 || W_t2_xypp > 0) sc->connection_type_node_xypp[n]        = 1;
        sc->connection_type_node_xypm[n]                                                = 0;
        if (W_xypm > 0) sc->connection_type_node_xypm[n]                                = 2;
        else if (W_t_xypm > 0 || W_t2_xypm > 0) sc->connection_type_node_xypm[n]        = 1;
        sc->connection_type_node_xzpp[n]                                                = 0;
        if (W_xzpp > 0) sc->connection_type_node_xzpp[n]                                = 2;
        else if (W_t_xzpp > 0 || W_t2_xzpp > 0) sc->connection_type_node_xzpp[n]        = 1;
        sc->connection_type_node_xzpm[n]                                                = 0;
        if (W_xzpm > 0) sc->connection_type_node_xzpm[n]                                = 2;
        else if (W_t_xzpm > 0 || W_t2_xzpm > 0) sc->connection_type_node_xzpm[n]        = 1;
        sc->connection_type_node_yzpp[n]                                                = 0;
        if (W_yzpp > 0) sc->connection_type_node_yzpp[n]                                = 2;
        else if (W_t_yzpp > 0 || W_t2_yzpp > 0) sc->connection_type_node_yzpp[n]        = 1;
        sc->connection_type_node_yzpm[n]                                                = 0;
        if (W_yzpm > 0) sc->connection_type_node_yzpm[n]                                = 2;
        else if (W_t_yzpm > 0 || W_t2_yzpm > 0) sc->connection_type_node_yzpm[n]        = 1;
        sc->connection_type_node_xyzppp[n]                                              = 0;
        if (W_xyzppp > 0) sc->connection_type_node_xyzppp[n]                            = 2;
        else if (W_t_xyzppp > 0 || W_t2_xyzppp > 0) sc->connection_type_node_xyzppp[n]  = 1;
        sc->connection_type_node_xyzppm[n]                                              = 0;
        if (W_xyzppm > 0) sc->connection_type_node_xyzppm[n]                            = 2;
        else if (W_t_xyzppm > 0 || W_t2_xyzppm > 0) sc->connection_type_node_xyzppm[n]  = 1;
        sc->connection_type_node_xyzpmp[n]                                              = 0;
        if (W_xyzpmp > 0) sc->connection_type_node_xyzpmp[n]                            = 2;
        else if (W_t_xyzpmp > 0 || W_t2_xyzpmp > 0) sc->connection_type_node_xyzpmp[n]  = 1;
        sc->connection_type_node_xyzmpp[n]                                              = 0;
        if (W_xyzmpp > 0) sc->connection_type_node_xyzmpp[n]                            = 2;
        else if (W_t_xyzmpp > 0 || W_t2_xyzmpp > 0) sc->connection_type_node_xyzmpp[n]  = 1;

        /*fprintf(debug_outputs, "V3 ps:\t");
          fprintf(debug_outputs, "%f %f %f ", x, y, z);
          fprintf(debug_outputs, "Wx %f Wy %f Wz %f ", sc->gGap_node_xx[n], sc->gGap_node_yy[n], sc->gGap_node_zz[n]);
          fprintf(debug_outputs, "Wxypp %f Wxypm %f Wxzpp %f Wxzpm %f Wyzpp %f Wyzpm %f ", sc->gGap_node_xypp[n], sc->gGap_node_xypm[n], sc->gGap_node_xzpp[n], sc->gGap_node_xzpm[n], sc->gGap_node_yzpp[n], sc->gGap_node_yzpm[n]);
          fprintf(debug_outputs, "xyzppp %f xzyppm %f xyzpmp %f xzy mpp %f ", sc->gGap_node_xyzppp[n], sc->gGap_node_xyzppm[n], sc->gGap_node_xyzpmp[n], sc->gGap_node_xyzmpp[n]);

          float sum2D = sc->gGap_node_xx[n] + sc->gGap_node_yy[n] + sc->gGap_node_xypp[n] + sc->gGap_node_xypm[n];

          float sum =  sc->gGap_node_xx[n] + sc->gGap_node_yy[n] + sc->gGap_node_zz[n] +
          sc->gGap_node_xypp[n] + sc->gGap_node_xypm[n] + sc->gGap_node_xzpp[n] +  sc->gGap_node_xzpm[n] + sc->gGap_node_yzpp[n] +  sc->gGap_node_yzpm[n] +
          sc->gGap_node_xyzppp[n] + sc->gGap_node_xyzppm[n] + sc->gGap_node_xyzpmp[n] + sc->gGap_node_xyzmpp[n];

          fprintf(debug_outputs, "sum 2D %f 3D %f \n\n", sum2D, sum);*/

        sc->gGap_node_xx[n]     *= 1.0/dx;
        sc->gGap_node_yy[n]     *= 1.0/dy;
        sc->gGap_node_zz[n]     *= 1.0/dz;

        sc->gGap_node_xypp[n]   *= 1.0/sqrt(dx*dx + dy*dy);
        sc->gGap_node_xypm[n]   *= 1.0/sqrt(dx*dx + dy*dy);
        sc->gGap_node_xzpp[n]   *= 1.0/sqrt(dx*dx + dz*dz);
        sc->gGap_node_xzpm[n]   *= 1.0/sqrt(dx*dx + dz*dz);
        sc->gGap_node_yzpp[n]   *= 1.0/sqrt(dy*dy + dz*dz);
        sc->gGap_node_yzpm[n]   *= 1.0/sqrt(dy*dy + dz*dz);

        sc->gGap_node_xyzppp[n] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
        sc->gGap_node_xyzppm[n] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
        sc->gGap_node_xyzpmp[n] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
        sc->gGap_node_xyzmpp[n] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);

    } // end loop over all tissue =============================================================================================//|
    
    if (nofibreerrorflag == true) printf("ERROR: at least one node did not have an associated myocyte orientation; defaulted to axial connections in all directions to avoid reduced coupling implications\n");

    /*fprintf(debug_outputs2, "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS 3 3 3\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA 27\nSCALARS ImageFile float 3\nLOOKUP_TABLE default\n");
    for (int newcount = 0; newcount < sc->N; newcount++) fprintf(debug_outputs2, "%f %f %f ", x_t, y_t, z_t);
    fprintf(debug_outputs3, "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS 3 3 3\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA 27\nSCALARS ImageFile float 3\nLOOKUP_TABLE default\n");
    for (int newcount = 0; newcount < sc->N; newcount++) fprintf(debug_outputs3, "%f %f %f ", x_t2, y_t2, z_t2);
    fclose(debug_outputs);
    fclose(debug_outputs2);
    fclose(debug_outputs3);*/
}
// End set gGgap node array (apply orientation weights) ======================//|

// Calculate nunber of junctions (cellular connections) ======================\\|
void calc_N_junctions(SC_variables *sc)
{
    // simply returns Njunc for dynamic array allocation
    // Only need to search in positive direction, as negative direction neighbours will have already been
    // set from positive direction of its negative neighbour
    int idx;                // 3D identifier for each cell
    int idx_xp;             // 3D identifier for x+1 cell
    int idx_yp;             // 3D identifier for y+1 cell
    int idx_zp;             // 3D identifier for z+1 cell
    int idx_xm_yp;
    int idx_xm_zp;
    int idx_xp_yp;
    int idx_xp_zp;
    int idx_ym_zp;
    int idx_yp_zp;

    int idx_xm_ym_zp;
    int idx_xm_yp_zp;
    int idx_xp_ym_zp;
    int idx_xp_yp_zp;

    sc->Njunc = 0;

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
                idx             =  i    + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_xp          = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_yp          =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_zp          =  i    + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));

                idx_xm_yp       = (i-1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xm_zp       = (i-1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_xp_yp       = (i+1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xp_zp       = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_ym_zp       =  i    + (sc->NX * (j-1))      + (sc->NX * sc->NY * (k+1));
                idx_yp_zp       =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY * (k+1));

                idx_xm_ym_zp = (i-1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xm_yp_zp = (i-1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_ym_zp = (i+1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_yp_zp = (i+1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));

                if (sc->geo[idx] > 0) // if it is an actual cell/node
                {
                    if (i < sc->NX-1 &&                                 sc->geo[idx_xp] > 0)        sc->Njunc++;
                    if (j < sc->NY-1 &&                                 sc->geo[idx_yp] > 0)        sc->Njunc++;
                    if (k < sc->NZ-1 &&                                 sc->geo[idx_zp] > 0)        sc->Njunc++;
                    if (i < sc->NX-1 && j < sc->NY-1 &&                 sc->geo[idx_xp_yp] > 0)     sc->Njunc++;
                    if (i > 0        && j < sc->NY-1 &&                 sc->geo[idx_xm_yp] > 0)     sc->Njunc++;
                    if (i < sc->NX-1 && k < sc->NZ-1 &&                 sc->geo[idx_xp_zp] > 0)     sc->Njunc++;
                    if (i > 0        && k < sc->NZ-1 &&                 sc->geo[idx_xm_zp] > 0)     sc->Njunc++;
                    if (j < sc->NY-1 && k < sc->NZ-1 &&                 sc->geo[idx_yp_zp] > 0)     sc->Njunc++;
                    if (j > 0        && k < sc->NZ-1 &&                 sc->geo[idx_ym_zp] > 0)     sc->Njunc++;
                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xp_yp_zp] > 0)  sc->Njunc++;
                    if (i > 0        && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xm_yp_zp] > 0)  sc->Njunc++;
                    if (i < sc->NX-1 && j > 0        && k < sc->NZ-1 && sc->geo[idx_xp_ym_zp] > 0)  sc->Njunc++;
                    if (i > 0        && j > 0        && k < sc->NZ-1 && sc->geo[idx_xm_ym_zp] > 0)  sc->Njunc++;
                } // end if
            } // end x for
        } // end y for
    } // end z for
} // end function
// End Calculate nunber of junctions (cellular connections) ====================//|

// Set junctional conduction and connection maps ===============================\\|
void set_gjunc_and_junc_maps(SC_variables *sc, const char* Output_dir, Tissue_parameters *t)
{
    int count = 0;      // counter of number of cells
    int count_junc = 0; // junction counter

    int idx;                // 3D identifier for each cell
    int idx_xp;             // 3D identifier for x+1 cell
    int idx_yp;             // 3D identifier for y+1 cell
    int idx_zp;             // 3D identifier for z+1 cell
    int idx_xm_yp;
    int idx_xm_zp;
    int idx_xp_yp;
    int idx_xp_zp;
    int idx_ym_zp;
    int idx_yp_zp;
    int idx_xm_ym_zp;
    int idx_xm_yp_zp;
    int idx_xp_ym_zp;
    int idx_xp_yp_zp;

    bool xx;
    bool yy;
    bool zz;
    bool xypp;
    bool xypm;
    bool xzpp;
    bool xzpm;
    bool yzpp;
    bool yzpm;
    bool xyzppp;
    bool xyzppm;
    bool xyzpmp;
    bool xyzmpp;

    printf("Creating network connections and assigning gap conductance\n");

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
                idx             =  i    + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_xp          = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_yp          =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_zp          =  i    + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));

                idx_xm_yp       = (i-1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xm_zp       = (i-1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_xp_yp       = (i+1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xp_zp       = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_ym_zp       =  i    + (sc->NX * (j-1))      + (sc->NX * sc->NY * (k+1));
                idx_yp_zp       =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY * (k+1));

                idx_xm_ym_zp = (i-1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xm_yp_zp = (i-1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_ym_zp = (i+1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_yp_zp = (i+1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));

                if (sc->geo[idx] > 0) // if it is an actual cell/node
                {
                    xx      = false;
                    yy      = false;
                    zz      = false;
                    xypp    = false;
                    xypm    = false;
                    xzpp    = false;
                    xzpm    = false;
                    yzpp    = false;
                    yzpm    = false;
                    xyzppp  = false;
                    xyzppm  = false;
                    xyzpmp  = false;
                    xyzmpp  = false;

                    // if this positive neighbour exists, set junction to true
                    if (i < sc->NX-1 &&                     sc->geo[idx_xp] > 0)    xx      = true;
                    if (j < sc->NY-1 &&                     sc->geo[idx_yp] > 0)    yy      = true;
                    if (k < sc->NZ-1 &&                     sc->geo[idx_zp] > 0)    zz      = true;
                    if (i < sc->NX-1 && j < sc->NY-1 &&     sc->geo[idx_xp_yp] > 0) xypp    = true;
                    if (i > 0        && j < sc->NY-1 &&     sc->geo[idx_xm_yp] > 0) xypm    = true;
                    if (i < sc->NX-1 && k < sc->NZ-1 &&     sc->geo[idx_xp_zp] > 0) xzpp    = true;
                    if (i > 0        && k < sc->NZ-1 &&     sc->geo[idx_xm_zp] > 0) xzpm    = true;
                    if (j < sc->NY-1 && k < sc->NZ-1 &&     sc->geo[idx_yp_zp] > 0) yzpp    = true;
                    if (j > 0        && k < sc->NZ-1 &&     sc->geo[idx_ym_zp] > 0) yzpm    = true;

                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xp_yp_zp] > 0)  xyzppp = true;
                    if (i > 0        && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xm_yp_zp] > 0)  xyzmpp = true;
                    if (i < sc->NX-1 && j > 0        && k < sc->NZ-1 && sc->geo[idx_xp_ym_zp] > 0)  xyzpmp = true;
                    if (i > 0        && j > 0        && k < sc->NZ-1 && sc->geo[idx_xm_ym_zp] > 0)  xyzppm = true;

                    // Unset junctions if regions should not be electrically coupled
                    if (t->disconnect_regions_flag == true)
                    {
                        for (int DN = 0; DN < t->Ndisconnected_regions; DN++) // loop over number of disconnected region pairs
                        {
                            // if current cell is region 1 and neighbour region 2, or current cell is reg 2 and nei reg 1, uncouple by setting that junction flag back to false
                            if (i < sc->NX-1)
                                if( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp] == t->disconnect_regions[DN][0]) ) xx = false;
                            if (j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp] == t->disconnect_regions[DN][0]) ) yy = false;
                            if (k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_zp] == t->disconnect_regions[DN][0]) ) zz = false;
                            if (i < sc->NX-1 && j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][0]) ) xypp = false;
                            if (i > 0 && j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][0]) ) xypm = false;
                            if (i < sc->NX-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][0]) ) xzpp = false;
                            if (i > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][0]) ) xzpm = false;
                            if (j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][0]) ) yzpp = false;
                            if (j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][0]) ) yzpm = false;
                            if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][0]) ) xyzppp = false;
                            if (i > 0 && j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][0]) ) xyzmpp = false;
                            if (i < sc->NX-1 && j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][0]) ) xyzpmp = false;
                            if (i > 0 && j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][0]) ) xyzppm = false;
                        } // end DN for
                    } // end disconnect uncoupled regions

                    // calculate condutance at each junction
                    // Principal directions
                    if (xx == true)
                    {
                        // find 1D,Ncell index of plus and minus neighbours
                        // we want jn_map, which is of Njunc and incremented with each junc, to return the 1DNcell index (which the voltage array is of)
                        // : minus is itself (as the junction is between itself and plus neighour) and plus is plus neighbour
                        sc->jn_map_minus[count_junc]    = count; // or = sc->geo_index[idx], same thing
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xp]; // sets the 1D index of xplus to the i+1 node
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xx[count] + sc->gGap_node_xx[sc->geo_index[idx_xp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xx[count] == 1 && sc->connection_type_node_xx[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xx[count] == 2 && sc->connection_type_node_xx[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xx[count]*sc->connection_type_node_xx[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    if (yy == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_yp]; // sets the 1D index of xplus to the i+1 node
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_yy[count] + sc->gGap_node_yy[sc->geo_index[idx_yp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_yy[count] == 1 && sc->connection_type_node_yy[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_yy[count] == 2 && sc->connection_type_node_yy[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_yy[count]*sc->connection_type_node_yy[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    if (zz == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_zz[count] + sc->gGap_node_zz[sc->geo_index[idx_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_zz[count] == 1 && sc->connection_type_node_zz[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_zz[count] == 2 && sc->connection_type_node_zz[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_zz[count]*sc->connection_type_node_zz[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    // Diagonals
                    if (xypp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xp_yp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xypp[count] + sc->gGap_node_xypp[sc->geo_index[idx_xp_yp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xypp[count] == 1 && sc->connection_type_node_xypp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xypp[count] == 2 && sc->connection_type_node_xypp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xypp[count]*sc->connection_type_node_xypp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    if (xypm == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xm_yp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xypm[count] + sc->gGap_node_xypm[sc->geo_index[idx_xm_yp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xypm[count] == 1 && sc->connection_type_node_xypm[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xypm[count] == 2 && sc->connection_type_node_xypm[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xypm[count]*sc->connection_type_node_xypm[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    if (xzpp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xp_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xzpp[count] + sc->gGap_node_xzpp[sc->geo_index[idx_xp_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xzpp[count] == 1 && sc->connection_type_node_xzpp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xzpp[count] == 2 && sc->connection_type_node_xzpp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xzpp[count]*sc->connection_type_node_xzpp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    if (xzpm == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xm_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xzpm[count] + sc->gGap_node_xzpm[sc->geo_index[idx_xm_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xzpm[count] == 1 && sc->connection_type_node_xzpm[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xzpm[count] == 2 && sc->connection_type_node_xzpm[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xzpm[count]*sc->connection_type_node_xzpm[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    if (yzpp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_yp_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_yzpp[count] + sc->gGap_node_yzpp[sc->geo_index[idx_yp_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_yzpp[count] == 1 && sc->connection_type_node_yzpp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_yzpp[count] == 2 && sc->connection_type_node_yzpp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_yzpp[count]*sc->connection_type_node_yzpp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    if (yzpm == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_ym_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_yzpm[count] + sc->gGap_node_yzpm[sc->geo_index[idx_ym_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_yzpm[count] == 1 && sc->connection_type_node_yzpm[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_yzpm[count] == 2 && sc->connection_type_node_yzpm[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_yzpm[count]*sc->connection_type_node_yzpm[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    // Corners
                    if (xyzppp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xp_yp_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xyzppp[count] + sc->gGap_node_xyzppp[sc->geo_index[idx_xp_yp_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xyzppp[count] == 1 && sc->connection_type_node_xyzppp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xyzppp[count] == 2 && sc->connection_type_node_xyzppp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xyzppp[count]*sc->connection_type_node_xyzppp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }

                    if (xyzppm == true) // ppm = mmp
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xm_ym_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xyzppm[count] + sc->gGap_node_xyzppm[sc->geo_index[idx_xm_ym_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xyzppm[count] == 1 && sc->connection_type_node_xyzppm[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xyzppm[count] == 2 && sc->connection_type_node_xyzppm[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xyzppm[count]*sc->connection_type_node_xyzppm[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)


                        count_junc++;
                    }
                    if (xyzpmp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xp_ym_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xyzpmp[count] + sc->gGap_node_xyzpmp[sc->geo_index[idx_xp_ym_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xyzpmp[count] == 1 && sc->connection_type_node_xyzpmp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xyzpmp[count] == 2 && sc->connection_type_node_xyzpmp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xyzpmp[count]*sc->connection_type_node_xyzpmp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    if (xyzmpp == true)
                    {
                        sc->jn_map_minus[count_junc]    = count;
                        sc->jn_map_plus[count_junc]     = sc->geo_index[idx_xm_yp_zp];
                        sc->gGap_jn[count_junc]        = (sc->gGap_node_xyzmpp[count] + sc->gGap_node_xyzmpp[sc->geo_index[idx_xm_yp_zp]])/2; // junction g is average of the two coupled cells

                        sc->connection_type_jn[count_junc] = 0;
                        if (sc->connection_type_node_xyzmpp[count] == 1 && sc->connection_type_node_xyzmpp[sc->geo_index[idx_xp]] == 1)     sc->connection_type_jn[count_junc] = 1;  // connection is transverse
                        if (sc->connection_type_node_xyzmpp[count] == 2 && sc->connection_type_node_xyzmpp[sc->geo_index[idx_xp]] == 2)     sc->connection_type_jn[count_junc] = 2;  // connection is axial
                        if (sc->connection_type_node_xyzmpp[count]*sc->connection_type_node_xyzmpp[sc->geo_index[idx_xp]] == 2)             sc->connection_type_jn[count_junc] = 3;  // connection is mixed (1*2 or 2*1)

                        count_junc++;
                    }
                    count++;
                } // end if actua node
            } // end x loop
        } // end y loop
    } // end z loop
    printf("Junction count check = %d\n", count_junc);
}
// End Set junctional conduction and connection maps ===========================//

// Calculate gap jucntion flux and assign to each neighbour ====================\\|
void calc_IGap(SC_variables *sc, double *v, int n)
{
    sc->IGap[n]    = sc->gGap_jn[n] * (v[sc->jn_map_plus[n]] - v[sc->jn_map_minus[n]]);
    sc->diff[sc->jn_map_plus[n]]    += -sc->IGap[n];
    sc->diff[sc->jn_map_minus[n]]   +=  sc->IGap[n];
}
// Calculate gao jucntion flux and assign to each neighbour ====================//|

// Remove/modify network connections ===========================================\\|
void default_junction_maps(SC_variables *sc)
{
    for (int n = 0; n < sc->Njunc; n++) sc->gGgap_mod_map[n] = 1.0;
    for (int n = 0; n < sc->Njunc; n++) sc->gGgap_base_map[n] = 1.0;
}

// maps then read in from file if applied to scale the above

void update_junctions(SC_variables *sc) 
{
    for (int n = 0; n < sc->Njunc; n++)
    {
        sc->gGap_jn[n] *= (sc->gGgap_mod_map[n]*sc->gGgap_base_map[n]); // Diff_mod_map should = 1 unless specifically set
    }
}
// End Remove network connections ==============================================//|

// Write network connections visualisation files ===============================\\|
void output_junc_maps(SC_variables *sc, const char* Output_dir, Tissue_parameters *t)
{
    int count = 0;      // counter of number of cells
    int count_junc = 0; // junction counter

    int idx;                // 3D identifier for each cell
    int idx_xp;             // 3D identifier for x+1 cell
    int idx_yp;             // 3D identifier for y+1 cell
    int idx_zp;             // 3D identifier for z+1 cell
    int idx_xm_yp;
    int idx_xm_zp;
    int idx_xp_yp;
    int idx_xp_zp;
    int idx_ym_zp;
    int idx_yp_zp;
    int idx_xm_ym_zp;
    int idx_xm_yp_zp;
    int idx_xp_ym_zp;
    int idx_xp_yp_zp;

    bool xx;
    bool yy;
    bool zz;
    bool xypp;
    bool xypm;
    bool xzpp;
    bool xzpm;
    bool yzpp;
    bool yzpm;
    bool xyzppp;
    bool xyzppm;
    bool xyzpmp;
    bool xyzmpp;

    // files for writing for vis
    char *string = (char*)malloc(500);

    sprintf(string, "%s/Connections_xx.txt", Output_dir);
    FILE *outxx;
    outxx = fopen(string, "wt");
    sprintf(string, "%s/Connections_xx.vtk", Output_dir);
    FILE *outxx_vtk;
    outxx_vtk = fopen(string, "wt");

    fprintf(outxx_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxx_vtk, "vtk outxx_vtkput\n");
    fprintf(outxx_vtk, "ASCII\n");
    fprintf(outxx_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxx_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY, sc->NZ);
    fprintf(outxx_vtk, "SPACING 1 1 1\n");
    fprintf(outxx_vtk, "ORIGIN 0.5 0 0\n");
    fprintf(outxx_vtk, "POINT_DATA %d\n", (sc->NX-1)*sc->NY*sc->NZ);
    fprintf(outxx_vtk, "SCALARS connections float 1\n");
    fprintf(outxx_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_yy.txt", Output_dir);
    FILE *outyy;
    outyy = fopen(string, "wt");
    sprintf(string, "%s/Connections_yy.vtk", Output_dir);
    FILE *outyy_vtk;
    outyy_vtk = fopen(string, "wt");

    fprintf(outyy_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outyy_vtk, "vtk outyy_vtkput\n");
    fprintf(outyy_vtk, "ASCII\n");
    fprintf(outyy_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outyy_vtk, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY-1, sc->NZ);
    fprintf(outyy_vtk, "SPACING 1 1 1\n");
    fprintf(outyy_vtk, "ORIGIN 0 0.5 0\n");
    fprintf(outyy_vtk, "POINT_DATA %d\n", sc->NX*(sc->NY-1)*sc->NZ);
    fprintf(outyy_vtk, "SCALARS connections float 1\n");
    fprintf(outyy_vtk, "LOOKUP_TABLE default\n");


    sprintf(string, "%s/Connections_xypp.txt", Output_dir);
    FILE *outxypp;
    outxypp = fopen(string, "wt");
    sprintf(string, "%s/Connections_xypp.vtk", Output_dir);
    FILE *outxypp_vtk;
    outxypp_vtk = fopen(string, "wt");

    fprintf(outxypp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxypp_vtk, "vtk outxypp_vtkput\n");
    fprintf(outxypp_vtk, "ASCII\n");
    fprintf(outxypp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxypp_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ);
    fprintf(outxypp_vtk, "SPACING 1 1 1\n");
    fprintf(outxypp_vtk, "ORIGIN 0.5 0.5 0\n");
    fprintf(outxypp_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*sc->NZ);
    fprintf(outxypp_vtk, "SCALARS connections float 1\n");
    fprintf(outxypp_vtk, "LOOKUP_TABLE default\n");


    sprintf(string, "%s/Connections_xypm.txt", Output_dir);
    FILE *outxypm;
    outxypm = fopen(string, "wt");
    sprintf(string, "%s/Connections_xypm.vtk", Output_dir);
    FILE *outxypm_vtk;
    outxypm_vtk = fopen(string, "wt");

    fprintf(outxypm_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxypm_vtk, "vtk outxypm_vtkput\n");
    fprintf(outxypm_vtk, "ASCII\n");
    fprintf(outxypm_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxypm_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ);
    fprintf(outxypm_vtk, "SPACING 1 1 1\n");
    fprintf(outxypm_vtk, "ORIGIN 0.5 0.5 0\n");
    fprintf(outxypm_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*sc->NZ);
    fprintf(outxypm_vtk, "SCALARS connections float 1\n");
    fprintf(outxypm_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_zz.vtk", Output_dir);
    FILE *outzz_vtk;
    outzz_vtk = fopen(string, "wt");

    fprintf(outzz_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outzz_vtk, "vtk outzz_vtkput\n");
    fprintf(outzz_vtk, "ASCII\n");
    fprintf(outzz_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outzz_vtk, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY, sc->NZ-1);
    fprintf(outzz_vtk, "SPACING 1 1 1\n");
    fprintf(outzz_vtk, "ORIGIN 0 0.5 0\n");
    fprintf(outzz_vtk, "POINT_DATA %d\n", sc->NX*sc->NY*(sc->NZ-1));
    fprintf(outzz_vtk, "SCALARS connections float 1\n");
    fprintf(outzz_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_xzpp.vtk", Output_dir);
    FILE *outxzpp_vtk;
    outxzpp_vtk = fopen(string, "wt");

    fprintf(outxzpp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxzpp_vtk, "vtk outxzpp_vtkput\n");
    fprintf(outxzpp_vtk, "ASCII\n");
    fprintf(outxzpp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxzpp_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY, sc->NZ-1);
    fprintf(outxzpp_vtk, "SPACING 1 1 1\n");
    fprintf(outxzpp_vtk, "ORIGIN 0.5 0 0.5\n");
    fprintf(outxzpp_vtk, "POINT_DATA %d\n", (sc->NX-1)*sc->NY*(sc->NZ-1));
    fprintf(outxzpp_vtk, "SCALARS connections float 1\n");
    fprintf(outxzpp_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_xzpm.vtk", Output_dir);
    FILE *outxzpm_vtk;
    outxzpm_vtk = fopen(string, "wt");

    fprintf(outxzpm_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxzpm_vtk, "vtk outxzpm_vtkput\n");
    fprintf(outxzpm_vtk, "ASCII\n");
    fprintf(outxzpm_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxzpm_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY, sc->NZ-1);
    fprintf(outxzpm_vtk, "SPACING 1 1 1\n");
    fprintf(outxzpm_vtk, "ORIGIN 0.5 0 0.5\n");
    fprintf(outxzpm_vtk, "POINT_DATA %d\n", (sc->NX-1)*sc->NY*(sc->NZ-1));
    fprintf(outxzpm_vtk, "SCALARS connections float 1\n");
    fprintf(outxzpm_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_yzpp.vtk", Output_dir);
    FILE *outyzpp_vtk;
    outyzpp_vtk = fopen(string, "wt");

    fprintf(outyzpp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outyzpp_vtk, "vtk outyzpp_vtkput\n");
    fprintf(outyzpp_vtk, "ASCII\n");
    fprintf(outyzpp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outyzpp_vtk, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY-1, sc->NZ-1);
    fprintf(outyzpp_vtk, "SPACING 1 1 1\n");
    fprintf(outyzpp_vtk, "ORIGIN 0 0.5 0.5\n");
    fprintf(outyzpp_vtk, "POINT_DATA %d\n", sc->NX*(sc->NY-1)*(sc->NZ-1));
    fprintf(outyzpp_vtk, "SCALARS connections float 1\n");
    fprintf(outyzpp_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_yzpm.vtk", Output_dir);
    FILE *outyzpm_vtk;
    outyzpm_vtk = fopen(string, "wt");

    fprintf(outyzpm_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outyzpm_vtk, "vtk outyzpm_vtkput\n");
    fprintf(outyzpm_vtk, "ASCII\n");
    fprintf(outyzpm_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outyzpm_vtk, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY-1, sc->NZ-1);
    fprintf(outyzpm_vtk, "SPACING 1 1 1\n");
    fprintf(outyzpm_vtk, "ORIGIN 0 0.5 0.5\n");
    fprintf(outyzpm_vtk, "POINT_DATA %d\n", sc->NX*(sc->NY-1)*(sc->NZ-1));
    fprintf(outyzpm_vtk, "SCALARS connections float 1\n");
    fprintf(outyzpm_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_xyzppp.vtk", Output_dir);
    FILE *outxyzppp_vtk;
    outxyzppp_vtk = fopen(string, "wt");

    fprintf(outxyzppp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxyzppp_vtk, "vtk outxyzppp_vtkput\n");
    fprintf(outxyzppp_vtk, "ASCII\n");
    fprintf(outxyzppp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxyzppp_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ-1);
    fprintf(outxyzppp_vtk, "SPACING 1 1 1\n");
    fprintf(outxyzppp_vtk, "ORIGIN 0.5 0.5 0.5\n");
    fprintf(outxyzppp_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*(sc->NZ-1));
    fprintf(outxyzppp_vtk, "SCALARS connections float 1\n");
    fprintf(outxyzppp_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_xyzppm.vtk", Output_dir);
    FILE *outxyzppm_vtk;
    outxyzppm_vtk = fopen(string, "wt");

    fprintf(outxyzppm_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxyzppm_vtk, "vtk outxyzppm_vtkput\n");
    fprintf(outxyzppm_vtk, "ASCII\n");
    fprintf(outxyzppm_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxyzppm_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ-1);
    fprintf(outxyzppm_vtk, "SPACING 1 1 1\n");
    fprintf(outxyzppm_vtk, "ORIGIN 0.5 0.5 0.5\n");
    fprintf(outxyzppm_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*(sc->NZ-1));
    fprintf(outxyzppm_vtk, "SCALARS connections float 1\n");
    fprintf(outxyzppm_vtk, "LOOKUP_TABLE default\n");

    sprintf(string, "%s/Connections_xyzpmp.vtk", Output_dir);
    FILE *outxyzpmp_vtk;
    outxyzpmp_vtk = fopen(string, "wt");

    fprintf(outxyzpmp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxyzpmp_vtk, "vtk outxyzpmp_vtkput\n");
    fprintf(outxyzpmp_vtk, "ASCII\n");
    fprintf(outxyzpmp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxyzpmp_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ-1);
    fprintf(outxyzpmp_vtk, "SPACING 1 1 1\n");
    fprintf(outxyzpmp_vtk, "ORIGIN 0.5 0.5 0.5\n");
    fprintf(outxyzpmp_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*(sc->NZ-1));
    fprintf(outxyzpmp_vtk, "SCALARS connections float 1\n");
    fprintf(outxyzpmp_vtk, "LOOKUP_TABLE default\n");


    sprintf(string, "%s/Connections_xyzmpp.vtk", Output_dir);
    FILE *outxyzmpp_vtk;
    outxyzmpp_vtk = fopen(string, "wt");

    fprintf(outxyzmpp_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(outxyzmpp_vtk, "vtk outxyzmpp_vtkput\n");
    fprintf(outxyzmpp_vtk, "ASCII\n");
    fprintf(outxyzmpp_vtk, "DATASET STRUCTURED_POINTS\n");
    fprintf(outxyzmpp_vtk, "DIMENSIONS %d %d %d\n", sc->NX-1, sc->NY-1, sc->NZ-1);
    fprintf(outxyzmpp_vtk, "SPACING 1 1 1\n");
    fprintf(outxyzmpp_vtk, "ORIGIN 0.5 0.5 0.5\n");
    fprintf(outxyzmpp_vtk, "POINT_DATA %d\n", (sc->NX-1)*(sc->NY-1)*(sc->NZ-1));
    fprintf(outxyzmpp_vtk, "SCALARS connections float 1\n");
    fprintf(outxyzmpp_vtk, "LOOKUP_TABLE default\n");

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
                idx             =  i    + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_xp          = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY *  k   );
                idx_yp          =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_zp          =  i    + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));

                idx_xm_yp       = (i-1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xm_zp       = (i-1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_xp_yp       = (i+1) + (sc->NX * (j+1))      + (sc->NX * sc->NY *  k   );
                idx_xp_zp       = (i+1) + (sc->NX *  j)         + (sc->NX * sc->NY * (k+1));
                idx_ym_zp       =  i    + (sc->NX * (j-1))      + (sc->NX * sc->NY * (k+1));
                idx_yp_zp       =  i    + (sc->NX * (j+1))      + (sc->NX * sc->NY * (k+1));

                idx_xm_ym_zp = (i-1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xm_yp_zp = (i-1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_ym_zp = (i+1) + (sc->NX * (j-1)) + (sc->NX * sc->NY * (k+1));
                idx_xp_yp_zp = (i+1) + (sc->NX * (j+1)) + (sc->NX * sc->NY * (k+1));

                if (sc->geo[idx] > 0) // if it is an actual cell/node
                {
                    xx      = false;
                    yy      = false;
                    zz      = false;
                    xypp    = false;
                    xypm    = false;
                    xzpp    = false;
                    xzpm    = false;
                    yzpp    = false;
                    yzpm    = false;
                    xyzppp  = false;
                    xyzppm  = false;
                    xyzpmp  = false;
                    xyzmpp  = false;

                    // if this positive neighbour exists, set junction to true
                    if (i < sc->NX-1 &&                     sc->geo[idx_xp] > 0)    xx      = true;
                    if (j < sc->NY-1 &&                     sc->geo[idx_yp] > 0)    yy      = true;
                    if (k < sc->NZ-1 &&                     sc->geo[idx_zp] > 0)    zz      = true;
                    if (i < sc->NX-1 && j < sc->NY-1 &&     sc->geo[idx_xp_yp] > 0) xypp    = true;
                    if (i > 0        && j < sc->NY-1 &&     sc->geo[idx_xm_yp] > 0) xypm    = true;
                    if (i < sc->NX-1 && k < sc->NZ-1 &&     sc->geo[idx_xp_zp] > 0) xzpp    = true;
                    if (i > 0        && k < sc->NZ-1 &&     sc->geo[idx_xm_zp] > 0) xzpm    = true;
                    if (j < sc->NY-1 && k < sc->NZ-1 &&     sc->geo[idx_yp_zp] > 0) yzpp    = true;
                    if (j > 0        && k < sc->NZ-1 &&     sc->geo[idx_ym_zp] > 0) yzpm    = true;

                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xp_yp_zp] > 0)  xyzppp = true;
                    if (i > 0        && j < sc->NY-1 && k < sc->NZ-1 && sc->geo[idx_xm_yp_zp] > 0)  xyzmpp = true;
                    if (i < sc->NX-1 && j > 0        && k < sc->NZ-1 && sc->geo[idx_xp_ym_zp] > 0)  xyzpmp = true;
                    if (i > 0        && j > 0        && k < sc->NZ-1 && sc->geo[idx_xm_ym_zp] > 0)  xyzppm = true;

                    // Unset junctions if regions should not be electrically coupled
                    if (t->disconnect_regions_flag == true)
                    {
                        for (int DN = 0; DN < t->Ndisconnected_regions; DN++) // loop over number of disconnected region pairs
                        {
                            // if current cell is region 1 and neighbour region 2, or current cell is reg 2 and nei reg 1, uncouple by setting that junction flag back to false
                            if (i < sc->NX-1)
                                if( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp] == t->disconnect_regions[DN][0]) ) xx = false;
                            if (j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp] == t->disconnect_regions[DN][0]) ) yy = false;
                            if (k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_zp] == t->disconnect_regions[DN][0]) ) zz = false;
                            if (i < sc->NX-1 && j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][0]) ) xypp = false;
                            if (i > 0 && j < sc->NY-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][0]) ) xypm = false;
                            if (i < sc->NX-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][0]) ) xzpp = false;
                            if (i > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][0]) ) xzpm = false;
                            if (j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][0]) ) yzpp = false;
                            if (j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][0]) ) yzpm = false;
                            if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][0]) ) xyzppp = false;
                            if (i > 0 && j < sc->NY-1 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][0]) ) xyzmpp = false;
                            if (i < sc->NX-1 && j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][0]) ) xyzpmp = false;
                            if (i > 0 && j > 0 && k < sc->NZ-1)
                                if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][0]) ) xyzppm = false;
                        } // end DN for
                    } // end disconnect uncoupled regions

                    // calculate condutance at each junction
                    // Principal directions
                    if (xx == true)
                    {
                        fprintf(outxx, "%d %d %d %f\n", i, j, k, sc->gGap_jn[count_junc]);
                        if (i < sc->NX-1) fprintf(outxx_vtk, "%f ", sc->gGap_jn[count_junc]); // only print vtk for NX-1 as connections are between nodes
                        count_junc++;
                    }
                    else if (i < sc->NX-1) fprintf(outxx_vtk, "0 ");

                    if (yy == true)
                    {
                        fprintf(outyy, "%d %d %d %f\n", i, j, k, sc->gGap_jn[count_junc]);
                        if (j < sc->NY-1) fprintf(outyy_vtk, "%f ", sc->gGap_jn[count_junc]);
                        count_junc++;
                    }
                    else if (j < sc->NY-1) fprintf(outyy_vtk, "0 ");

                    if (zz == true)
                    {
                        if (k < sc->NZ-1) fprintf(outzz_vtk, "%f ", sc->gGap_jn[count_junc]);
                        count_junc++;
                    }
                    else if (k < sc->NZ-1) fprintf(outzz_vtk, "0 ");

                    // Diagonals
                    if (xypp == true)
                    {
                        fprintf(outxypp, "%d %d %d %f\n", i, j, k, sc->gGap_jn[count_junc]*sqrt(2));
                        if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypp_vtk, "0 ");

                    if (xypm == true)
                    {
                        fprintf(outxypm, "%d %d %d %f\n", i-1, j, k, sc->gGap_jn[count_junc]*sqrt(2));
                        if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypm_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypm_vtk, "0 ");

                    if (xzpp == true)
                    {
                        if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpp_vtk, "0 ");

                    if (xzpm == true)
                    {
                        if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpm_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpm_vtk, "0 ");

                    if (yzpp == true)
                    {
                        if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpp_vtk, "0 ");

                    if (yzpm == true)
                    {
                        if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpm_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(2));
                        count_junc++;
                    }
                    else if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpm_vtk, "0 ");

                    // Corners
                    if (xyzppp == true)
                    {
                        if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(3));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppp_vtk, "0 ");

                    if (xyzppm == true) // ppm = mmp
                    {
                        if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppm_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(3));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppm_vtk, "0 ");

                    if (xyzpmp == true)
                    {
                        if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzpmp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(3));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzpmp_vtk, "0 ");

                    if (xyzmpp == true)
                    {
                        if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzmpp_vtk, "%f ", sc->gGap_jn[count_junc]*sqrt(3));
                        count_junc++;
                    }
                    else if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzmpp_vtk, "0 ");

                    count++;
                } // end if actual node
                else
                {
                    if (i < sc->NX-1) fprintf(outxx_vtk, "-100 ");
                    if (j < sc->NY-1) fprintf(outyy_vtk, "-100 ");
                    if (k < sc->NZ-1) fprintf(outzz_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypp_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1) fprintf(outxypm_vtk, "-100 ");
                    if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpp_vtk, "-100 ");
                    if (i < sc->NX-1 && k < sc->NZ-1) fprintf(outxzpm_vtk, "-100 ");
                    if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpp_vtk, "-100 ");
                    if (j < sc->NY-1 && k < sc->NZ-1) fprintf(outyzpm_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppp_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzppm_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzpmp_vtk, "-100 ");
                    if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1) fprintf(outxyzmpp_vtk, "-100 ");
                }
            } // end x loop
            fprintf(outxx_vtk, "\n");
            fprintf(outyy_vtk, "\n");
            fprintf(outzz_vtk, "\n");
            fprintf(outxypp_vtk, "\n");
            fprintf(outxypm_vtk, "\n");
            fprintf(outxzpp_vtk, "\n");
            fprintf(outxzpm_vtk, "\n");
            fprintf(outyzpp_vtk, "\n");
            fprintf(outyzpm_vtk, "\n");
            fprintf(outxyzppp_vtk, "\n");
            fprintf(outxyzppm_vtk, "\n");
            fprintf(outxyzpmp_vtk, "\n");
            fprintf(outxyzmpp_vtk, "\n");
        } // end y loop
        fprintf(outxx_vtk, "\n");
        fprintf(outyy_vtk, "\n");
        fprintf(outzz_vtk, "\n");
        fprintf(outxypp_vtk, "\n");
        fprintf(outxypm_vtk, "\n");
        fprintf(outxzpp_vtk, "\n");
        fprintf(outxzpm_vtk, "\n");
        fprintf(outyzpp_vtk, "\n");
        fprintf(outyzpm_vtk, "\n");
        fprintf(outxyzppp_vtk, "\n");
        fprintf(outxyzppm_vtk, "\n");
        fprintf(outxyzpmp_vtk, "\n");
        fprintf(outxyzmpp_vtk, "\n");
    } // end z loop

    printf("Junction visualisation files written. Junction count check = %d\n", count_junc);

    fclose(outxx);
    fclose(outxx_vtk);
    fclose(outyy_vtk);
    fclose(outzz_vtk);
    fclose(outxypp_vtk);
    fclose(outxypm_vtk);
    fclose(outyy);
    fclose(outxypp);
    fclose(outxypm);
    fclose(outxzpp_vtk);
    fclose(outxzpm_vtk);
    fclose(outyzpp_vtk);
    fclose(outyzpm_vtk);
    fclose(outxyzppp_vtk);
    fclose(outxyzppm_vtk);
    fclose(outxyzpmp_vtk);
    fclose(outxyzmpp_vtk);
}
// End write network connections visualisation files ===========================//|
// End NETWORK model implementation ==============================================================//|

