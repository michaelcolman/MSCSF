// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Tool for converting spatial binary data to ==  //
// plain text and vtk data - tissue models. ===============  //
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
//      PERTAINS TO THE CITATION OF COLMAN 2022 PLOS COMP =  //
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

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include <omp.h>

#include "lib/Arguments.h"
#include "lib/Initialisation.h"
#include "lib/Structs.h"
#include "lib/Model.h"
#include "lib/Read_write_state.h"
#include "lib/Outputs.h"
#include "lib/Spatial_coupling.h"
#include "lib/Tissue.h"

using namespace std;

int main(int argc, char *argv[])
{
    // Read in path for geometry and state files ========\\|
    char PATH[1000];
    FILE *path_in;
    path_in = fopen("PATH.txt", "r");
    if (path_in == NULL) // Check if file can be opened
    {
        printf("ERROR: Cannot open \"PATH.txt\"\n");
        printf("Ensure it is present in the current directory\n");
        exit(1);
    }
    fscanf(path_in, "%[^\n]", PATH); // read up to new line into PATH char
    printf("\n*PATH location is %s*\n", PATH);
    fclose(path_in);
    // End read path ====================================//|

    Argument_parameters Argin;                                  // Initialise struct                            || lib/Arguments.h
    set_argument_defaults(&Argin);
    int counter = 1;

    double remove_axial_perc;
    double remove_trans_perc;

    remove_axial_perc = 0.2;
    remove_trans_perc = 0.8;

    double dist_width;
    dist_width = 0.3;

    double dist_mean;
    dist_mean = 1.0;

    double threshold_upper_axial = 1;
    double threshold_lower_axial = 0;
    double threshold_upper_trans = 1;
    double threshold_lower_trans = 0;

    // MAP
    char * map_file = (char*)malloc(500);
    bool apply_by_map = false;

    // Relevant arguments for this tool
    while (counter < argc)
    {
        if (strcmp(argv[counter], "Tissue_order") == 0)
        {
            Argin.Tissue_order        = argv[counter+1];
            Argin.Tissue_order_arg    = true;
            counter++;
        }
        else if (strcmp(argv[counter], "Tissue_model") == 0)
        {
            Argin.Tissue_model        = argv[counter+1];
            Argin.Tissue_model_arg    = true;
            counter++;
        }
        else if (strcmp(argv[counter], "axial_perc") == 0)
        {
            remove_axial_perc         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "trans_perc") == 0)
        {
            remove_trans_perc         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "distribution_width") == 0)
        {
            dist_width                    = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "distribution_mean") == 0)
        {
            dist_mean                    = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "threshold_upper_axial") == 0)
        {
            threshold_upper_axial         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "threshold_lower_axial") == 0)
        {
            threshold_lower_axial         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "threshold_upper_trans") == 0)
        {
            threshold_upper_trans         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "threshold_lower_trans") == 0)
        {
            threshold_lower_trans         = atof(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "Apply_map_file") == 0)
        {
            map_file            = argv[counter+1];
            apply_by_map        = true;
            counter++;
        }
        else
        {
            printf("ERROR: \"%s\" is not a valid argument for this tool\n", argv[counter]);
            printf("Please use ONLY:\n");
            printf("\tTissue_order  [1D/2D/3D/geo]\t Tissue_model [basic, ...]\n");
            printf("\taxial_perc [0-1]\ttrans_perc [0-1]\n");
            printf("\tdistribution_width [x]\tdistribution_mean [x]\n");
            printf("\tthreshold_upper_axial [x]\tthreshold_lower_axial [x]\tthreshold_upper_trans [x]\tthreshold_lower_trans [x]\n");
            printf("\tApply_map_file [filename]\n");
            exit(1);
        }
        counter++;
    }

    printf("\n\n*********** TOOL SETTINGS **************\n");
    printf("This tool will always create three maps (obvioulsy you don't need to use all of them)\n");
    printf("\tDiscrete map (connections are on or off)\n");
    printf("\t\tAxial proportion to be removed = %f\tTransverse = %f\n\n", remove_axial_perc, remove_trans_perc);
    printf("\tRand map (connections are scaled by normal disttribution\n");
    printf("\t\tdistribution_width = %f\tdistribution_mean = %f\n\n", dist_width, dist_mean);
    printf("\tContinuous map (connections are scaled by defined distributions; map be different for axial or transverse\n");
    printf("\t\tthreshold_upper_axial = %f\tthreshold_lower_axial = %f\tthreshold_upper_trans = %f\tthreshold_lower_trans = %f\n", threshold_upper_axial, threshold_lower_axial, threshold_upper_trans, threshold_lower_trans);
    if (apply_by_map == false) printf("\t\tModifications will be applied to all tissue\n\n");
    if (apply_by_map == true) printf("\t\tModifications will be applied only to regions specified by remodelling map %s\n\n", map_file);

    // Simulation settings, including refs which may be relevant for identifying the folders
    Simulation_parameters Sim;                                  // Initialise sim parameters struct || lib/Structs.h
    set_simulation_defaults(&Sim, 0.02);                        // (sim struct, dt)                 || lib/Initialisation.c
    set_simulation_settings(&Sim, Argin, "native");             // Set from arguments               || lib/Initialisation.c

    Cell_parameters                 Params_global;
    Tissue_parameters               Tissue;
    SC_variables                    SC;
    RAND                            *Rand;      // radnom number
    int                             *Scale_map; // map of size Njunc which determines whether to apply modification or not                                            

    set_model_conditions(&Params_global, Argin);
    set_tissue_model_conditions(&Tissue, Argin);    // lib/Tissue.cpp
    if (strcmp(Tissue.Tissue_order, "geo") == 0) set_tissue_settings_anatomical(Params_global, &Tissue);
    else set_tissue_settings_idealised(Params_global, &Tissue);

    printf("Orientation type is %s\n", Tissue.Orientation_type);
    
    // Create directory to hold maps
    char * mkdirectory  = (char*)malloc(500);
    sprintf(mkdirectory, "mkdir -p %s", "Heterogeneous_connection_maps");
    system(mkdirectory);

    // Allocate arrays || lib/Spatial_coupling.cpp
    SC_set_array_sizes(&SC, Tissue.NX, Tissue.NY, Tissue.NZ);   // Sets array sizes in SC struct from tissue settings
    SC_array_allocation_N3(&SC, SC.NX, SC.NY, SC.NZ);           // Allocates arrays of size NX*NY*NZ || geo and 3D->1D geo index
    printf(">Spatial coupling NX*NY*NZ arrays allocated\n");
    select_tissue_geometry_function(Tissue, &SC, PATH, "Heterogeneous_connection_maps");  // lib/Tissue.cpp
    printf("\tGeometry size (X*Y*Z, %d * %d * %d) || Ncells = %d\n\n", Tissue.NX, SC.NY, SC.NZ, SC.N);
    
    // NETWORK calc Njunc based on geo
    calc_N_junctions(&SC);
    printf("\tN junctions calculated  = %d\n", SC.Njunc);

    // Allocate arrays size Ncell
    SC_array_allocation_Ncell(&SC, SC.N);       // lib/Spatial_coupling.cpp || geo_linear, D arrays, neighbour map, orientation, laplacian components
    tissue_array_allocation(&Tissue, SC.N);     // lib/Tissue.cpp || stim/ISO/remodelling etc map arrays
    printf(">Spatial coupling Ncell arrays allocated\n");

    // NETWORK Allocate Njunc arrays
    SC_array_allocation_Njunc(&SC, SC.Njunc);
    printf(">Spatial coupling Njunc arrays allocated\n");


    SC_set_index_and_geo_linear(&SC);               // lib/Spatial_coupling.cpp
    SC_set_neighbours(&SC);                         // lib/Spatial_coupling.cpp
    printf(">Linear index and neighbours set\n");

    set_G_dx_global(&SC, Tissue.dx, Tissue.dy, Tissue.dz, Tissue.Gl, Tissue.Gt, Tissue.Gt2, 1.0, 1.0, 1.0);

    set_orientation(&SC, Tissue, PATH, Tissue.Tissue_order);    // lib/Tissue.cpp || This sets orientation to 0, then sets/reads in IF set to anisotropic
    output_fibre_orientation(SC, Tissue, "Heterogeneous_connection_maps");

    set_gGgap_array(&SC, Tissue.Orientation_type);
    set_gjunc_and_junc_maps(&SC, "Heterogeneous_connection_maps", &Tissue);

    default_junction_maps(&SC);
    
    // Read map if selected
    if (apply_by_map == true)
    {
        Tissue.map_in_type = "file";
        create_or_read_map_double(&Tissue, SC, PATH, "Heterogeneous_connection_maps", Tissue.Dscale_mod_map, map_file, "Dmod"); // lib/Tissue.cpp
    }
    else for (int i = 0; i < SC.N; i++) Tissue.Dscale_mod_map[i] = 1; // set to mod everywhere 

    // Assign nodal map to junction map
    Scale_map = new int [SC.Njunc];
    for (int i = 0; i < SC.Njunc; i++)
    {
        if (Tissue.Dscale_mod_map[SC.jn_map_plus[i]] > 0.0 || Tissue.Dscale_mod_map[SC.jn_map_minus[i]] > 0.0) Scale_map[i] = 1;
        else Scale_map[i] = 0;
        //Scale_map[i] = 1;
    }

    // Random numbers
    Rand    = new RAND [SC.Njunc];
    double rand;

    // Now for the actual maps - this tool will always make multiple types of map, you can use whichever you like
    char * discrete_map_name = (char*)malloc(500);
    sprintf(discrete_map_name, "Heterogeneous_connection_maps/%s_Junction_discrete_map_%f_%f.dat", Tissue.Tissue_model, remove_axial_perc, remove_trans_perc);
    char * rand_map_name = (char*)malloc(500);
    sprintf(rand_map_name, "Heterogeneous_connection_maps/%s_Junction_rand_map_%f_%f.dat", Tissue.Tissue_model, dist_mean, dist_width);
    char * cont_map_name = (char*)malloc(500);
    sprintf(cont_map_name, "Heterogeneous_connection_maps/%s_Junction_cont_map_%f_%f_%f_%f.dat", Tissue.Tissue_model, threshold_upper_axial, threshold_lower_axial, threshold_upper_trans, threshold_lower_trans);
    
    // Discrete map (connection on or off)
    FILE *discrete_map; 
    discrete_map = fopen(discrete_map_name, "wt");

    for (int n = 0; n < SC.Njunc; n++)
    {
        if (Scale_map[n] == 1)
        {
            rand = Rand[n].mtrand1();
            //cout << "rand " << rand << endl;
            if (SC.connection_type_jn[n] == 1) // connection is transverse
            {
                if (rand < remove_trans_perc) SC.gGgap_mod_map[n] = 0;
            }
            else if (SC.connection_type_jn[n] == 2) // connection is axial
            {
                if (rand < remove_axial_perc) SC.gGgap_mod_map[n] = 0;
            }
            else if (SC.connection_type_jn[n] == 3) // connection is mixed -> base threshold on random
            {
                if (rand < ((remove_trans_perc+remove_axial_perc)/2)) SC.gGgap_mod_map[n] = 0;
            }
        }
        fprintf(discrete_map, "%f\n", SC.gGgap_mod_map[n]);
    }

    // set maps back to default
    default_junction_maps(&SC);

    // Random sample from simple normal distribution (width and mean defined; axial and transverse treated identically)
    FILE *rand_map;
    rand_map = fopen(rand_map_name, "wt"); 

    for (int n = 0; n < SC.Njunc; n++)
    {
        if (Scale_map[n] == 1)
        {
            rand = Rand[n].mtrand1();
            //cout << "rand 2 " << rand << endl;
            SC.gGgap_base_map[n] = -dist_width*log(1.0/rand   - 1) + dist_mean;
            if (SC.gGgap_base_map[n] < 0) SC.gGgap_base_map[n] = 0;
        }
        fprintf(rand_map, "%f\n", SC.gGgap_base_map[n]);
    }

    // set maps back to default
    default_junction_maps(&SC);

    // From defined distribution
    FILE *cont_map;
    cont_map = fopen(cont_map_name, "wt");

    double width_a, mid_a; // axial map
    double width_t, mid_t; // transverse map

    width_a = (threshold_upper_axial - threshold_lower_axial)/2;
    mid_a = threshold_upper_axial - width_a;

    width_t = (threshold_upper_trans - threshold_lower_trans)/2;
    mid_t = threshold_upper_trans - width_t;

    for (int n = 0; n < SC.Njunc; n++)
    {
        if (Scale_map[n] == 1)
        {
            rand = Rand[n].mtrand1();

            if (SC.connection_type_jn[n] == 1) // connection is transverse
            {
                SC.gGgap_mod_map[n] = -width_t*0.2*log((1.0/rand) -1) + mid_t;
                //SC.gGgap_mod_map[n] = 1.0/(exp(-(rand-mid_t)/(width_t*0.2))+1);
                //printf("trans: %f\n", SC.gGgap_mod_map[n]);
            }
            else if (SC.connection_type_jn[n] == 2) // connection is axial
            {
                SC.gGgap_mod_map[n] = -width_a*0.2*log((1.0/rand) -1) + mid_a;
                //SC.gGgap_mod_map[n] = 1.0/(exp(-(rand-mid_a)/(width_a*0.2))+1);
                //printf("axial: %f\n", SC.gGgap_mod_map[n]);
            }
            else if (SC.connection_type_jn[n] == 3) // connection is mixed -> base threshold on random
            {
                SC.gGgap_mod_map[n] = -((width_a+width_t)/2)*0.2*log((1/rand) -1) + (mid_a+mid_t)/2;
                //SC.gGgap_mod_map[n] = 1.0/(exp(-(rand-((mid_a+mid_t)/2))/(((width_a+width_t)/2)*0.2))+1);
                //printf("mixed: %f\n", SC.gGgap_mod_map[n]);
            }
            if (SC.gGgap_mod_map[n] < 0) SC.gGgap_mod_map[n] = 0; // important as dists can easily go below zero
            if (SC.gGgap_mod_map[n] > 1.0) SC.gGgap_mod_map[n] = 1.0; // important as dists can easily go above one
        }
        fprintf(cont_map, "%f\n", SC.gGgap_mod_map[n]);

    }
    // delete
    delete [] SC.geo;
    delete [] SC.geo_index;

    delete [] Rand;
    delete [] Scale_map;

    fclose(discrete_map);
    fclose(rand_map);
    fclose(cont_map);

} // end main
