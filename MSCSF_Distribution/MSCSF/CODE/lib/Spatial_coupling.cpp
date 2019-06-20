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
}
// End deallocate all arrays ======================================//|
// End Allocate and deallocate spatial arrays ===================================================//|

// Read geometry/maps from file into arrays =====================================================\\|
int read_geo_file(SC_variables *sc, int *geo, const char *filein, const char * fileroot, const char *PATH, const char* Output_dir, const char * ref)
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

