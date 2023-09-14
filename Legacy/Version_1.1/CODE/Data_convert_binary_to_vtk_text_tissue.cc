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

    int start_time, end_time, interval; // over which to convert data
    const char * variable;              // variable to read/write (Vm/Cai/CaSR)
    bool write_vtk = true;              // Whether to write vtk file
    bool write_data = false;            // Whether to write plain text data file
    //bool model_type_native = true;      // for native or integrated tissue models
    const char * model_type = "native";
    variable = "Vm";  // default

	start_time 	= 0;
	end_time 	= 10;
	interval	= 10;

    Argument_parameters Argin;                                  // Initialise struct                            || lib/Arguments.h
    set_argument_defaults(&Argin);
    int counter = 1;
    while (counter < argc)
    {
        if (strcmp(argv[counter], "Reference") == 0 || strcmp(argv[counter], "reference") == 0) // if there is no difference between argument and "BCL"
        {
            Argin.reference                = argv[counter+1];         // Reads argument into reference variable
            Argin.reference_arg            = true;                     // Argument has been passed
            counter++;
        }
        else if (strcmp(argv[counter], "Results_Reference") == 0 || strcmp(argv[counter], "results_reference") == 0) // if there is no difference between argument and "BCL"
        {
            Argin.results_reference                = argv[counter+1];         // Reads argument into reference variable
            Argin.results_reference_arg            = true;                     // Argument has been passed
            counter++;
        }
        else if (strcmp(argv[counter], "Model") == 0)
        {
            Argin.Model            = argv[counter+1];
            Argin.Model_arg        = true;
            counter++;
        }
        else if (strcmp(argv[counter], "Model_type") == 0)
        {
            model_type            = argv[counter+1];
            counter++;
        }
        else if (strcmp(argv[counter], "Tissue_order") == 0)
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

        else if (strcmp(argv[counter], "Variable") == 0)
        {
            variable        = argv[counter+1];
            counter++;
        }

        else if (strcmp(argv[counter], "start_time") == 0)
        {
            start_time        = atoi(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "end_time") == 0)
        {
            end_time        = atoi(argv[counter+1]);
            counter++;
        }
        else if (strcmp(argv[counter], "interval") == 0)
        {
            interval        = atoi(argv[counter+1]);
            counter++;
        }

        else if (strcmp(argv[counter], "Write_vtk") == 0)
        {
            if (strcmp(argv[counter+1], "On") == 0) write_vtk = true;
            else if (strcmp(argv[counter+1], "Off") == 0) write_vtk = false;
            else
            {
                printf("ERROR: Write_vtk can only be On or Off\n");
            }
            counter++;
        }
        else if (strcmp(argv[counter], "Write_data") == 0)
        {
            if (strcmp(argv[counter+1], "On") == 0) write_data = true;
            else if (strcmp(argv[counter+1], "Off") == 0) write_data = false;
            else
            {
                printf("ERROR: Write_data can only be On or Off\n");
            }
            counter++;
        }

        else
        {
            printf("ERROR: \"%s\" is not a valid argument for this post processing\n", argv[counter]);
            printf("Please use ONLY:\n");
            printf("\tReference [text]\tResults_Reference [text]\tModel [text]\tTissue_order  [1D/2D/3D/geo]\t Tissue_model [basic, ...]\tModel_type [native/integrated]\n");
            printf("\tVarianble [Vm/Cai/CaSR]\tstart_time [int]\tend_time [int]\tinterval [n ms]\n");
            printf("\tWrite_vtk [On/Off]\tWrite_data [On/Off]");
            exit(1);
        }
        counter++;
    }


    // Simulation settings, including refs which may be relevant for identifying the folders
    Simulation_parameters Sim;                                  // Initialise sim parameters struct || lib/Structs.h
    set_simulation_defaults(&Sim, 0.02);                        // (sim struct, dt)                 || lib/Initialisation.c
    set_simulation_settings(&Sim, Argin, "native");             // Set from arguments               || lib/Initialisation.c

    // Set the folders to read data from
    // Output files =====================================\\|
    char * directory    = (char*)malloc(500);
    char * results_dir  = (char*)malloc(500);
    char * sr_dir       = (char*)malloc(500);
    char * filename_out = (char*)malloc(500);

    if (strcmp(model_type, "native") == 0) 
    {
        if (Argin.reference_arg == true) sprintf(directory, "Outputs_tissue_native_%s", Sim.reference);
        else sprintf(directory, "Outputs_tissue_native");
    }
    else if (strcmp(model_type, "integrated") == 0)
    {
        if (Argin.reference_arg == true) sprintf(directory, "Outputs_0Dtissue_%s", Sim.reference);
        else sprintf(directory, "Outputs_0Dtissue");
    }
    else
    {
        printf("ERROR: model_type \"%s\" is not valid: chose only \"native\" or \"integated\"\n");
        exit(1);
    }

    if (Argin.results_reference_arg == true) sprintf(results_dir, "Results_%s", Sim.results_reference);
    else sprintf(results_dir, "Results");
    sprintf(sr_dir, "Spatial_%s", results_dir);

    Cell_parameters                 Params_global;
    Tissue_parameters               Tissue;
    double                          *V;            // Whichever variable is to be written

    set_model_conditions(&Params_global, Argin);
    set_tissue_model_conditions(&Tissue, Argin);    // lib/Tissue.cpp
    if (strcmp(Tissue.Tissue_order, "geo") == 0) set_tissue_settings_anatomical(Params_global, &Tissue);
    else set_tissue_settings_idealised(Params_global, &Tissue);

    printf("\nbinary data to vtk being written for the following settings:\n");
    printf("\tTissue order: %s; Tissue model: %s ", Tissue.Tissue_order, Tissue.Tissue_model);
    if (strcmp(Tissue.Tissue_order, "geo") == 0) printf("and tissue geometry file is: %s\n", Tissue.geo_file);
    else printf("\n");
    printf("\tstart time: %d\n\tend time: %d\n\tinterval: %d\n", start_time, end_time, interval);
    printf("\tOutputs directory: %s\n\tSpatial results directory: %s\n", directory, sr_dir);
    printf("\tWriting for variable: %s ", variable);
    if (write_vtk == true) printf("\tWriting to vtk");
    if (write_data == true) printf("\tWriting to text data file");
    printf("\n\n\n");

    SC_variables                    SC;
    // Allocate arrays || lib/Spatial_coupling.cpp
    SC_set_array_sizes(&SC, Tissue.NX, Tissue.NY, Tissue.NZ);   // Sets array sizes in SC struct from tissue settings
    SC_array_allocation_N3(&SC, SC.NX, SC.NY, SC.NZ);           // Allocates arrays of size NX*NY*NZ || geo and 3D->1D geo index
    printf(">Spatial coupling NX*NY*NZ arrays allocated\n");
    select_tissue_geometry_function(Tissue, &SC, PATH, directory);  // lib/Tissue.cpp
    printf("\tGeometry size (X*Y*Z, %d * %d * %d) || Ncells = %d\n\n", Tissue.NX, SC.NY, SC.NZ, SC.N);

    // Allocate variable
    V = new double [SC.N];

    int iteration;

    for (iteration = start_time; iteration <= end_time; iteration += interval)
    {
        // Read in binary data || lib/Outputs.cpp
        array_1D_binary_read(variable, directory, sr_dir, V, SC, iteration);

        // Output as vtk
        if (write_vtk == true) 
        {
            sprintf(filename_out, "%s/%s/%s_output_%04d.vtk", directory, sr_dir, variable, iteration);
            printf("Creating visualisation file %s\n", filename_out);
            vtk_3D_output(variable, directory, sr_dir, V, SC, iteration);
        }

        // Output as data
        if (write_data == true) 
        {
            sprintf(filename_out, "%s/%s/%s_output_%04d.dat", directory, sr_dir, variable, iteration);
            printf("Creating visualisation file %s\n", filename_out);
            data_3D_output(variable, directory, sr_dir, V, SC, iteration);
        }
    }

    // delete
    delete [] V;
    delete [] SC.geo;
    delete [] SC.geo_index;
    free(directory);
    free(results_dir);
    free(sr_dir);

} // end main
