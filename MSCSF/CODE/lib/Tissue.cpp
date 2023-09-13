// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Tissue model setttings and options ==========  //
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

#include "Tissue.h"
#include "Structs.h"
#include "Spatial_coupling.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>

// Function list ================================================================================\\|
//	Setup and tissue model
//	    set_tissue_model_conditions()
//	    tissue_array_allocation()
//	    tissue_array_deallocation()
//	    set_tissue_settings_idealised() 	** This is where to add a new model **
//	    set_tissue_settings_anatomical()	** This is where to add a new model **
//	    set_coord_stim_and_map_from_defined_type()
//	    set_global_orientation_direction_from_arg()
//	    overwrite_tissue_properties_from_args()
//	
//	Create and read geometry, stimulus and map functions
//	    select_tissue_geometry_function()
//	    create_idealised_geometry_homogeneous()
//	    create_idealised_geometry_heterogeneous()
//	    select_stimulus_area_function()
//	    select_stimulus_area_function_multi_stim()
//	    create_stimulus_area()
//	    create_stimulus_area_sphere()
//	    set_orientation()
//	    create_orientation_ideal()
//	    read_orientation_anatomical()
//	    output_fibre_orientation()
//	    create_or_read_map_double()
//	    create_map_patch()
//	
//	Diffusion tensor non-uniformity
//	    update_D_arrays_Dscale_baseline()
//	    update_D_arrays_Dscale_mod()
//	
//	phase map for phase re-entry intiation
//	    create_read_phasemap()
//	    create_phasemap_2D()
//	    create_phasemap_3D()
//	
//	conduction velocity model parameters and calculation
//	    set_CV_cells()
//	    calculate_CV()
//	
//	compute_conduction_success()
// End Function list ============================================================================//|

// Set tissue model and type ====================================================================\\|
void set_tissue_model_conditions(Tissue_parameters *t, Argument_parameters A)
{
	// Defaults
	t->Tissue_order		= "2D";				// Dimension of model; 1D-3D + geo
	t->Tissue_model		= "basic"; 			// Model specifier
	t->Tissue_type		= "homogeneous";	// Eletrophysiology homogeneous or heterogeneous
	t->Orientation_type	= "isotropic";		// Orientation type: isotropic, ansitropic, orthotropic
	t->D_uniformity		= "uniform";		// Whether D magnitude varies in space , uniform, regional, map
	t->S1_loc_type		= "edge";			// Reference for basic stimulus settings
	t->S2_loc_type		= "S1";				// Reference for basic stimulus settings
	t->Ncelltypes		= 1;				// homogeneous
	t->stim_set			= false;			// stimulus has not yet been set
	t->S2_stim_set		= false;			// stimulus has not yet been set
	t->Dscale			= 1.0;				// i.e. not scaled
	t->Dscale_map_on    = "Off";            // Apply Dscale homogeneously, not according to map
	t->D_AR_scale       = 1.0;
	t->D_AR_scale_map_on    = "Off";        // Apply Dscale homogeneously, not according to map
	t->ISO_map_on		= "Off";			// Apply ISO homogeneously, not according to map
	t->remod_map_on		= "Off";			// Apply remodelling homogeneously, not according to map
	t->ACh_map_on		= "Off";			// Apply ACh homogeneously, not according to map
	t->SRF_map_on		= "Off";			// Apply SRF homogeneously, not according to map
	t->Direct_modulation_map_on		= "Off";// Apply Direct_modulation homogeneously, not according to map
	t->ideal_map_set    = false;
	t->S1_shape			= "cuboid";
	t->S2_shape			= "cuboid";
	t->ideal_map_shape  = "sphere";
	t->map_in_type      = "coords";         // this is auto set to file for geo models
	t->Multi_stim		= "Off";
	t->Nstims			= 1;
	t->Global_orientation_direction = "Off";
	t->spatial_gradient_map_on = "Off"; 
	t->Default_model    = "none";       // do not set from tissue defaults
	t->Tissue_model_2  = "none";
	t->Multiple_models  = "Off";

    t->Tissue_settings_set = false;

    // NETWORK
    t->junction_mod_map_on          = "Off";
    t->junction_base_map_on         = "Off";

    // Disconnect regions
    t->disconnect_regions_flag         = false;
    t->Ndisconnected_regions           = 0;

	// Overwrite from arguments || may need to do again if tissue settings set some of these
	if (A.Tissue_order_arg == true) 	t->Tissue_order		= A.Tissue_order;
	if (A.Tissue_model_arg == true) 	t->Tissue_model		= A.Tissue_model;
	if (A.Tissue_type_arg == true) 		t->Tissue_type 		= A.Tissue_type;
	if (A.Orientation_type_arg == true) t->Orientation_type	= A.Orientation_type;
	if (A.D_uniformity_arg == true) 	t->D_uniformity		= A.D_uniformity;
	if (A.Stimulus_type_arg == true) 	t->S1_loc_type		= A.Stimulus_loc_type;
	if (A.S2_Stimulus_type_arg == true)	t->S2_loc_type		= A.S2_Stimulus_loc_type;
	if (A.Dscale_arg == true)			t->Dscale			= A.Dscale;
	if (A.D_AR_scale_arg == true)   	t->D_AR_scale		= A.D_AR_scale;
	if (A.S1_shape_arg == true)			t->S1_shape			= A.S1_shape;
	if (A.S2_shape_arg == true)			t->S2_shape			= A.S2_shape;
	if (A.Tissue_model_2_arg == true)	t->Tissue_model_2	= A.Tissue_model_2;
	if (A.Multiple_models_arg == true)	t->Multiple_models	= A.Multiple_models;
}
// End set tissue model and type ================================================================//|

// Array allocation and deallocation ============================================================\\|
void tissue_array_allocation(Tissue_parameters *t, int N)
{
	t->stim_area	= new int [N];
	t->S2_stim_area	= new int [N];
	t->phasemap		= new int [N];
	t->ISO_map		= new double [N];
	t->Dscale_base_map	= new double [N];
	t->Dscale_mod_map	= new double [N];
	t->D_AR_scale_base_map	= new double [N];
	t->D_AR_scale_mod_map	= new double [N];
	t->remod_map	= new double [N];
	t->ACh_map		= new double [N];
	t->SRF_map		= new double [N];
	t->Direct_modulation_map		= new double [N];
	t->spatial_gradient_map		= new double [N];

	t->multi_stim_area	= new int *[20];
	for (int n = 0; n < 20; n++) t->multi_stim_area[n] = new int [N];
}
void tissue_array_deallocation(Tissue_parameters *t)
{
	delete [] t->stim_area;
	delete [] t->S2_stim_area;
	delete [] t->phasemap;
	delete [] t->ISO_map;
	delete [] t->Dscale_base_map;
	delete [] t->Dscale_mod_map;
	delete [] t->D_AR_scale_base_map;
	delete [] t->D_AR_scale_mod_map;
	delete [] t->remod_map;
	delete [] t->ACh_map;
	delete [] t->SRF_map;
	delete [] t->Direct_modulation_map;
	delete [] t->spatial_gradient_map;

	for (int n = 0; n < 20; n++) delete t->multi_stim_area[n];
	delete t->multi_stim_area;
}
// End array allocation and deallocation ========================================================//|

// Set tissue settings from model and type - IDEALISED ==========================================\\|
void set_tissue_settings_idealised(Cell_parameters p, Tissue_parameters *t)
{
	// In general, homogeneous tissue models can apply globally, but heterogeneous are model 
	// (or model-group)-specific, as the celltypes must be defined
	// NOTE: settings overwritten at end of function if arguments passed; so these define defaults

	// Default global orientations
	t->OX = t->OY = t->OZ = 0.0;

	// Defaulat sizes so that NX=0 does not happen
	t->NX = t->NY = t->NZ = 1;

    // Network model, default symmetry factors
    t->symm_fac_diag_long       = 1.0;
    t->symm_fac_diag_trans      = 1.0;
    t->symm_fac_corner          = 1.0;

    t->symm_fac_diag_long       = 0.77;//0.79;
    t->symm_fac_diag_trans      = 0.85;//0.91;
    t->symm_fac_corner          = 0.72;//0.72;
    
    t->apply_symmetry_factor    = "Off";
    
    t->Tissue_settings_set = false;

    // Fully functional example - use as a template for new mdoels ====================\\|
    // Remember: these settings are defaults only and can be overwritten by arguments
	if (strcmp(t->Tissue_model, "functional_model_test") == 0)  
	{
		// Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient 
		t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
		t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
		t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
		t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"
        t->Gl   = t->D1/t->dx;          // nS/pF; conversion from D to g; NETWORK model
        t->Gt   = t->Gl/t->D_AR;        // nS/pF

		// Geometry size
		if (strcmp(t->Tissue_order, "1D") == 0)
		{
			t->NX = 100;
			t->NY = t->NZ = 1;
		}
		else if (strcmp(t->Tissue_order, "2D") == 0)
		{
			t->NX = 300;
			t->NY = 300;
			t->NZ = 1;
		}
		else if (strcmp(t->Tissue_order, "3D") == 0)
		{
			t->NX = 100;
			t->NY = 100;
			t->NZ = 10;
		}
		// End Parameters which MUST be set =======================//|

		// Optional settings to differ from default ===============\\|
		// Default cell model
		t->Default_model    = "minimal";    // Unless argument is passed, will set cell model to this || default is "none" to set using single-cell defaults

		// Stimulus parameters ==========================\\|
		// shape based on coord parameters
		t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
		t->S2_shape         = "sphere";

		// Set coord parameters by:        
		// 1 - Select from pre defined (set in  set_coord_stim_and_map_from_defined_type)
		t->S1_loc_type      = "centre";     // edge, centre or whole_tissue || default is edge
		t->S2_loc_type      = "cross_field";     // idential to S1, edge, centre or cross_field || default is S1

		// OR 2 - define explicitly
		// S1
		//t->S1_x_size    = 5;                        // Extent either side of location: Total size = 5+5+1 = 11
		//t->S1_x_loc     = t->S1_x_size;             // Centre of stimulus in x; example is S1_size+1 th cell (from 0)
		//t->S1_y_loc     = int(float(t->NY/2.0));    // Centre of stimulus in x; example centre in y
		//t->S1_y_size    = int(float(t->NY/2.0));    // Extent either side of location: example is NY/2 + NY/2 + 1 >= NY
		//t->S1_z_loc     = int(float(t->NZ/2.0));    // entre of stimulus in x; example is centre in z
		//t->S1_z_size    = int(float(t->NZ/2.0));    // Extent either side of location: example is NZ/2 + NZ/2 + 1 >= NZ
		//t->stim_set     = true; // NEED this, otherwise it will overwrite with default type at end; you can still overwrite with argument

		// S2
		t->S2_x_size    = 5;                        // Extent either side of location: Total size = 5+5+1 = 11
		t->S2_x_loc     = int(float(t->NX/2.0));;             // Centre of stimulus in x; example is S2_size+1 th cell (from 0)
		t->S2_y_loc     = int(float(t->NY/3.0));    // Centre of stimulus in x; example centre in y
		t->S2_y_size    = 5;    // Extent either side of location: example is NY/2 + NY/2 + 1 >= NY
	    t->S2_z_loc     = int(float(t->NZ/2.0));    // entre of stimulus in x; example is centre in z
		t->S2_z_size    = int(float(t->NZ/2.0));    // Extent either side of location: example is NZ/2 + NZ/2 + 1 >= NZ
		t->S2_stim_set  = true; // NEED this, otherwise it will overwrite with default type at end; you can still overwrite with argument
		// End Stimulus parameters ======================//|

		// Electrophysuological heterogeneity ===========\\|
		//t->Tissue_type          = "heterogeneous"; // homogeneous or heterogeneous || default is homogeneous. Controlled by "Tissue_type [homogeneous/heterogeneous]"
        // uncomment above to make your tissue model heterogeneous by default

        // Below must be specified if you want to include heterogeneity in your model, regardless of its default being heterogeneous or homogeneous
		// Heterogeneity settings || Heterogeneity in x-direction only (use stim_location to change this direction relative to wave-propagation)
		t->Ncelltypes                   = 3;        // Number of cell types
		t->celltype_number[1]           = "ENDO";   // NOTE: starts from 1 as 0 is empty space
		t->celltype_number[2]           = "M";      // Must ensure you have a heterogeneity type for your cell model for each one listed here 
		t->celltype_number[3]           = "RA";
		t->het_junction_X_location[1]   = t->NX/3;      // End of celltype 1 is NX/3
		t->het_junction_X_location[2]   = 2*t->NX/3;    // End of celltype 2 is 2NX/3

		// Additional het option:
		// Multiple cell models || must be for each celltype, and celltype (het) must exist for that model!
		// This could be for atria and ventricles orto use two different models for different regions (e.g SAN vs atria)
		//t->Multiple_models              = "On"; // uncomment to default to two models for your tissue model. Controlled by "Multiple_models [Off/On]"
		t->Tissue_model_2  = "hAM_CRN";    // This is what will be assigned to areas given "Model_2" below
		//Like Default model, this is default and can be overwritten by arguments - controlled by "Tissue_model_2 [model_identifier]" 

        // Below must be specified to run two models in your tissue model, regardless of whether Multiple_models is On or Off by default.
		t->Modeltype_number[1]          = "Model_1"; // Apply the main cell model  // example is vent
		t->Modeltype_number[2]          = "Model_1"; // vent
		t->Modeltype_number[3]          = "Model_2"; // apply cell model in Tissue_model_2 (which could be differnet to that set above if
		// overwritten by arguments); this example is is atria
		// End Electrophysuological heterogeneity =======//|

		// Anisotropy ===================================\\|
		//t->Orientation_type     = "anisotropic"; // isotropic or anisotropic. Controlled by "Orientation_type [isotropic/anisotropic]"
        // uncomment above to make your tissue model anistropic by default

		// Global fibres (if anisptopic; has no effect if isotropic); must be set if anisotropic is set or passed as an argument
		t->OX = sqrt(0.5); // x-component of fibre orientation. Controlled by "OX [value 0-1]"
		t->OY = -sqrt(0.5);
		t->OZ = sqrt(1 - t->OX*t->OX - t->OY*t->OY); // normalised;
		// End Anisotropy ===============================//|

		// Non-uniform diffusion coefficient ============\\|
		// Heterogeneity settings (Ncelltypes, celltype_number[x], het_junction_X_location[x]) must be set to apply this functionality
		// Note: you can still apply homogeneous electrophysiology by setting Tissue_type to homogeneous here
		//t->D_uniformity         = "regional";     // uniform or regional (map is only valid for non-idealised, geo models); default is uniform
        // Controlled by "D_uniformity [uniform/regional]"
        // uncomment above to make your tissue model regional D_uniformity by default

        // Below must be specified to include regional D_uniformity in your tissue model, regardless of the default D_uniformtiy settings
		// Celltype dependant D and anisotropy ratio - must have settings for Ncelltypes
		// These will apply the scale factor to D1 and D_AR as set globally above
		t->non_uniform_D1_scale[1]          = 0.5;//1.0;  // ENDO
		t->non_uniform_AR_scale[1]          = 1.0; // ENDO
		t->non_uniform_D1_scale[2]          = 1.0;//1.0;  // M
		t->non_uniform_AR_scale[2]          = 2.0; // M
		t->non_uniform_D1_scale[3]          = 1.5;//0.5;  // EPI
		t->non_uniform_AR_scale[3]          = 0.5; // EPI
		// End Non-uniform diffusion coefficient ========//|

		// Spatially heterogeneous modulation map =======\\|
		// Modulations: ISO, ACh, Remodelling, Drug, SRF, Dscale, D_AR_scale, Direct_modulation (direct argument control of parameters e.g. Jup_scale, Ito_va_ss_vshift etc)
		// The same map is applied to all modulations for which map is on; 
		// for multiple map shapes, or maps with gradients, you must use a non-idealised geo model with map files
		// Note that the map will determine where the full magntiude of modulation applies; 
		// this magntiude is determined by the relevant condition variables (ISO, ACh, {Remodelling/Drug}_prop, Dscale, D_AR_Scale)
		// SRF and Direct_modulation maps are either on or off and do not have a scale amount

		// Set coordinates and shape of map; this one applies to right half of tissue
		t->ideal_map_x_loc      = t->NX - int(float(t->NX/4.0));
		t->ideal_map_y_loc      = int(float(t->NY/2.0));
		t->ideal_map_z_loc      = int(float(t->NZ/2.0));
		t->ideal_map_x_size     = int(float(t->NX/4.0));
		t->ideal_map_y_size     = int(float(t->NY/2.0));
		t->ideal_map_z_size     = int(float(t->NZ/2.0));
		t->ideal_map_shape      = "cuboid"; // default is sphere
		t->ideal_map_set        = true; // so it isn't overwritten by default

		// In general, this will only be applied if you pass arguments to apply modulation
		// However, you could also default these maps to be on for a specific reused model
		// You can't default the magntiude conditon variable (ISO, ACh, {Remodelling/Drug}_prop) here 
		// as this is set in other parts of the model
        // Uncomment below to default any map as being on in your tissue model
		//t->ISO_map_on		            = "On";		"ISO_map [On/Off]"	
		//t->remod_map_on		            = "On";		"Remodelling_map [On/Off]"
		//t->ACh_map_on		            = "On";		"ACh_map [On/Off]"
		//t->SRF_map_on		            = "On";		"SRF_map [On/Off]"
		//t->Direct_modulation_map_on		= "On";	 "Direct_modulation_map [On/Off]"

		// For Dscale and D_AR_scale, you can set a default value of the magntiude here (which can still be overwritten by arguments)
		//t->Dscale_map_on        = "On"; "Dscale_mod_map [On/Off]" 
		//t->D_AR_scale_map_on    = "On"; "D_AR_scale_mod_map [On/Off]"
		//t->Dscale			    = 0.5;  // multiplies local D1 (set by global D1 or regional D_uniformity) in map region
		//t->D_AR_scale           = 2.0;  // multiples local D_AR in map region
		// End Spatially heterogeneous modulation map ===//|
		
        // Disconnect regions ===========================\\|
        //disconnect_regions_flag              = true;
        //Ndisconnected_regions           = 2;
        
        // region pair one (disconnect 2 and 3)
        //disconnect_regions[0][0]        = 2;
        //disconnect_regions[0][1]        = 3;
        
        // region pair two (disconnect 2 and 3)
        //disconnect_regions[0][0]        = 2;
        //disconnect_regions[0][1]        = 3;
        
        // End disconnect regions =======================//|

        // End Optional settings to differ from default ===========//|
	}
	// End Template example || copy and paste where needed ============================//|

	// Global, homogeneous settings ===================================================\\|
	// Basic tissue model settings ==========================================\\|
	if (strcmp(t->Tissue_model, "basic") == 0 || strcmp(t->Tissue_model, "basic_small") == 0) // basic_small is for rapid testing)
	{
        t->Tissue_settings_set = true; // now the code knows the settings have been found

		// Settings
		t->Tissue_type          = "homogeneous";
		t->Orientation_type     = "isotropic";
		t->D_uniformity         = "uniform";

		// Global params
		t->dx = t->dy = t->dz = 0.2; 	// mm
		t->D1	= 0.2;					// mm^2/ms
		t->D_AR = 4;                    // anisotropy ratio (D1/D2)
		t->Diso	= 0.1;					// mm^2/ms
        t->Gl   = t->D1/t->dx;          // nS/pF; conversion from D to g 
        t->Gt   = t->Gl/t->D_AR;

		t->OX = 0.70711;
		t->OY = -0.70711;
		t->OZ = 0;
    
        if (strcmp(t->Tissue_model, "basic") == 0)
        {
            // Geometry size
            if (strcmp(t->Tissue_order, "1D") == 0)
            {
                t->NX = 100;
                t->NY = t->NZ = 1;
            }
            else if (strcmp(t->Tissue_order, "2D") == 0)
            {
                t->NX = 100;
                t->NY = 100;
                t->NZ = 1;
            }
            else if (strcmp(t->Tissue_order, "3D") == 0)
            {
                t->NX = 100;
                t->NY = 100;
                t->NZ = 10;
            }
        }
        // Overwrite sizes for basic small
        else if (strcmp(t->Tissue_model, "basic_small") == 0)
        {
            if (strcmp(t->Tissue_order, "1D") == 0)
            {
                t->NX = 10;
                t->NY = t->NZ = 1;
            }
            else if (strcmp(t->Tissue_order, "2D") == 0)
            {
                t->NX = 5;
                t->NY = 5;
                t->NZ = 1;
            }
            else if (strcmp(t->Tissue_order, "3D") == 0)
            {
                t->NX = 5;
                t->NY = 5;
                t->NZ = 2;
            }
        }
    }
    // End basic tissue model settings ======================================//|

    // Conduction velocity ==================================================\\|
    if (strcmp(t->Tissue_model, "conduction_velocity") == 0)
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Global params
        t->dx = t->dy = t->dz = 0.25;   // mm
        t->D1   = 0.4;                  // mm^2/ms
        t->D_AR = 4;                    // anisotropy ratio (D1/D2)
        t->Diso = 0.1;                  // mm^2/ms

        // NETWORK
        t->Gl   = t->D1/t->dx;          // nS/pF; conversion from D to g 
        t->Gt   = t->Gl/t->D_AR;

        t->symm_fac_diag_long       = 0.77;//0.79;
        t->symm_fac_diag_trans      = 0.85;//0.91;
        t->symm_fac_corner          = 0.72;//0.72;
        t->apply_symmetry_factor    = "On";

        t->OX = 1;
        t->OY = 0;
        t->OZ = 0;

        t->S1_loc_type      = "centre"; 	// specifics set at end of this function
        t->S2_loc_type      = "centre"; 	// specifics set at end of this function
        t->S1_shape			= "sphere";
        t->S2_shape			= "sphere";

        t->Tissue_type          = "homogeneous";
        t->Orientation_type     = "anisotropic";
        t->D_uniformity         = "uniform";

        if (strcmp(t->Tissue_order, "2D") == 0)
        {
            t->NX = 301;//201;
            t->NY = 301;//201;
            t->NZ = 1;
        }
        else if (strcmp(t->Tissue_order, "3D") == 0)
        {
            t->NX = 126;
            t->NY = 126;
            t->NZ = 126;	
        }
        else 
        {
            printf("ERROR: Conduction velocity tissue model only suitable for 2D or 3D\n");
            exit(1);
        }
    }
    // End Conduction velocity ==============================================//|

    // Fibre orientation tests ==============================================\\|
    if (strcmp(t->Tissue_model, "fibre_orientation") == 0)
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Global params
        t->dx = t->dy = t->dz = 0.25;    // mm
        t->D1   = 0.4;                  // mm/ms
                                        //t->D2   = 0.05;                 // mm/ms
        t->D_AR = 4;
        t->Diso = 0.25;

        // NETWORK base below on AR like above
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        t->OX = 1/sqrt(2);//-0.23;//1/(sqrt(2));
        t->OY = 1/(sqrt(2));//sqrt(1 - t->OZ*t->OZ);//0.7071;
        t->OZ = 0;//sqrt(1 - t->OX*t->OX - t->OY*t->OY);

        t->S1_loc_type      = "centre";     // specifics set at end of this function
        t->S2_loc_type      = "centre";     // specifics set at end of this function
        t->S1_shape         = "sphere";
        t->S2_shape         = "sphere";

        t->Tissue_type          = "homogeneous";
        t->Orientation_type     = "anisotropic";
        t->D_uniformity         = "uniform";

        if (strcmp(t->Tissue_order, "2D") == 0)
        {
            t->NX = 201;
            t->NY = 201;
            t->NZ = 1;
        }
        else if (strcmp(t->Tissue_order, "3D") == 0)
        {
            t->NX = 101;
            t->NY = 101;
            t->NZ = 101;
        }
        else
        {
            printf("ERROR: fibre_orientation tissue model only suitable for 2D or 3D\n");
            exit(1);
        }

    }
    // End Fibre orientation tests ==========================================//|

    // Idealised re-entry ===================================================\\|
    if (strcmp(t->Tissue_model, "re-entry") == 0 || strcmp(t->Tissue_model, "re-entry_large") == 0)
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Global params
        t->dx = t->dy = t->dz = 0.25;    // mm
        t->D1   = 0.2;                  // mm^2/ms
        t->D_AR = 4;                    // anisotropy ratio (D1/D2)
        t->Diso = 0.1;                  // mm^2/ms
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        t->S1_loc_type      = "edge";    
        t->S2_loc_type      = "cross_field";     // specifics set at end of this function
        t->S1_shape         = "cuboid";
        t->S2_shape         = "cuboid";

        t->Tissue_type          = "homogeneous";
        t->Orientation_type     = "isotropic";
        t->D_uniformity         = "uniform";

        t->OX = 0;//0.70711;
        t->OY = 0;//-0.70711;
        t->OZ = 0;

        // Geometry size
        if (strcmp(t->Tissue_order, "1D") == 0)
        {
            t->NX = 100;
            t->NY = 1;
            t->NZ = 1;
            //printf("ERROR: re-entry must be 2D or 3D\n");
            //exit(1);
        }
        else if (strcmp(t->Tissue_order, "2D") == 0)
        {
            t->NX = 200;
            t->NY = 200;
            t->NZ = 1;
        }
        else if (strcmp(t->Tissue_order, "3D") == 0)
        {
            t->NX = 200;
            t->NY = 200;
            t->NZ = 10;
        }

        if (strcmp(t->Tissue_model, "re-entry_large") == 0)
        {
            if (strcmp(t->Tissue_order, "2D") == 0)
            {
                t->NX = 400;
                t->NY = 400;
                t->NZ = 1;
            }
            else if (strcmp(t->Tissue_order, "3D") == 0)
            {
                t->NX = 400;
                t->NY = 400;
                t->NZ = 10;
            }
        }

        // Map for SRF; to remove edge effects
        t->ideal_map_x_loc      = int(float(t->NX/2.0));
        t->ideal_map_y_loc      = int(float(t->NY/2.0));
        t->ideal_map_z_loc      = int(float(t->NZ/2.0));
        t->ideal_map_x_size     = int(float(t->NX/2.0))-4;
        t->ideal_map_y_size     = int(float(t->NY/2.0))-4;
        t->ideal_map_z_size     = int(float(t->NZ/2.0))-4;
        t->ideal_map_shape      = "cuboid"; // default is sphere
        t->ideal_map_set        = true; // so it isn't overwritten by default
    }
    // end Idealised re-entry ===============================================//|
    // End Global, homogeneous ========================================================//|

    // More specific models============================================================\\|
    // Transmural strand/slice/slab through ventricular wall ================\\|
    if (strcmp(t->Tissue_model, "Vent_transmural") == 0) // Transmural strand/slice/slab through ventricular wall
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Global params
        t->Default_model        = "minimal";
        t->dx = t->dy = t->dz = 0.2;    // mm
        t->D1   = 0.2;                  // mm^2/ms
        t->D_AR = 4;                    // anisotropy ratio (D1/D2)
        t->Diso = 0.1;                  // mm^2/ms		
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        t->Tissue_type			= "heterogeneous";	
        t->Orientation_type		= "isotropic";
        t->D_uniformity			= "uniform";
        t->OX = t->OY = t->OZ = 0.0;
        t->OX = 1;

        // Geometry size
        if (strcmp(t->Tissue_order, "1D") == 0)
        {
            t->NX = 200;
            t->NY = t->NZ = 1;
        }
        else if (strcmp(t->Tissue_order, "2D") == 0)
        {
            t->NX = 200;
            t->NY = 100;
            t->NZ = 1;
        }
        else if (strcmp(t->Tissue_order, "3D") == 0)
        {
            t->NX = 100;
            t->NY = 100;
            t->NZ = 25;
        }			

        // Heterogeneity || Heterogeneity in x-direction only (use stim_location to change this direction relative to wave-propagation)
        // Note: no need to set if homogeneous, as celltype will then be set by Model default or command line-argument
        t->Ncelltypes					= 3;
        t->celltype_number[1]			= "ENDO";	// NOTE: starts from 1 as 0 is empty space
        t->celltype_number[2]			= "M";
        t->celltype_number[3]			= "EPI";
        t->het_junction_X_location[1]	= t->NX/3;		// End of celltype 1 is NX/3
        t->het_junction_X_location[2]	= 2*t->NX/3;	// End of celltype 2 is 2NX/3
        // junction location should always be celltypes-1

        // Celltype dependant D and anisotropy ratio
        //t->D_uniformity			= "regional";
        t->non_uniform_D1_scale[1]          = 0.5;//1.0;  // ENDO
        t->non_uniform_AR_scale[1]          = 1.0; // ENDO == global AR set
        t->non_uniform_D1_scale[2]          = 1.0;//1.0;  // M
        t->non_uniform_AR_scale[2]          = 2.0; // M
        t->non_uniform_D1_scale[3]          = 1.5;//0.5;  // EPI
        t->non_uniform_AR_scale[3]          = 0.5; // EPI
    }
    // End Transmural strand/slice/slab through ventricular wall ============//|

    // Heterogeneous uman atria  ============================================\\|
    // Currently, this has not been used for anything official and serves purpose as a template for full idealised controls - adjust as you please!
    if (strcmp(t->Tissue_model, "Heterogeneous_human_atria") == 0) // Strand/tissue with all human atrial celltypes
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

		// Global params
		t->Default_model        = "hAM_CRN";
		t->dx = t->dy = t->dz = 0.2;    // mm
		t->D1   = 0.2;                  // mm^2/ms
		t->D_AR = 10;                    // anisotropy ratio (D1/D2)
		t->Diso = 0.1;                  // mm^2/ms

        // NETWORK
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

		t->OX = t->OY = t->OZ = 0.0;
		t->OX = 1.0;

		t->Tissue_type          = "heterogeneous";
		t->Orientation_type     = "anisotropic";
		t->D_uniformity         = "regional";

		// Geometry size
		if (strcmp(t->Tissue_order, "1D") == 0)
		{
			t->NX = 200;
			t->NY = t->NZ = 1;
		}
		else if (strcmp(t->Tissue_order, "2D") == 0)
		{
			t->NX = 200;
			t->NY = 100;
			t->NZ = 1;
		}
		else if (strcmp(t->Tissue_order, "3D") == 0)
		{
			t->NX = 100;
			t->NY = 100;
			t->NZ = 10;		// Thin atrial wall!
		}

        // Stimulus
        t->S1_loc_type      = "centre";

		// Heterogeneity || Heterogeneity in x-direction only (use stim_location to change this direction relative to wave-propagation)
		t->Ncelltypes                   = 10;
		t->celltype_number[1]           = "RA";  // NOTE: starts from 1 as 0 is empty space
		t->celltype_number[2]           = "PM";
		t->celltype_number[3]           = "CT";
		t->celltype_number[4]           = "RAA";
		t->celltype_number[5]           = "AVR";
		t->celltype_number[6]           = "BB";
		t->celltype_number[7]           = "LA";
		t->celltype_number[8]           = "AS";
		t->celltype_number[9]           = "LAA";
		t->celltype_number[10]          = "PV";

		for (int i = 1; i < t->Ncelltypes; i++)
		{
			t->het_junction_X_location[i] = i*t->NX/t->Ncelltypes; // evenly distribute celltypes within geometry
		}
		// max junction location ref should always be celltypes-1

		// Non_uniform D scaling || Ensure numbers match celltypes defined ^
		// will only matter if D_uniformity is overwritten to regional
		t->non_uniform_D1_scale[1]			= 1.0;	// RA
		t->non_uniform_AR_scale[1]	        = 1.0;	// RA
		t->non_uniform_D1_scale[2]			= 1.2;	// PM
		t->non_uniform_AR_scale[2]			= 1.0;	// PM
		t->non_uniform_D1_scale[3]			= 1.0;	// CT
		t->non_uniform_AR_scale[3]			= 1.5;	// CT
		t->non_uniform_D1_scale[4]			= 0.8;	// RAA
		t->non_uniform_AR_scale[4]			= 0.5;	// RAA
		t->non_uniform_D1_scale[5]			= 0.5;	// AVR
		t->non_uniform_AR_scale[5]			= 0.3;	// AVR
		t->non_uniform_D1_scale[6]			= 1.0;	// BB
		t->non_uniform_AR_scale[6]			= 1.5;	// BB
		t->non_uniform_D1_scale[7]			= 1.1;	// LA
		t->non_uniform_AR_scale[7]			= 0.8;	// LA
		t->non_uniform_D1_scale[8]			= 0.5;	// AS
		t->non_uniform_AR_scale[8]			= 0.3;	// AS
		t->non_uniform_D1_scale[9]			= 0.5;	// LAA
		t->non_uniform_AR_scale[9]			= 0.3;	// LAA
		t->non_uniform_D1_scale[10]			= 0.8;	// PV
		t->non_uniform_AR_scale[10]			= 1.8;	// PV


        // Disconnect regions ===========================\\|
        //t->disconnect_regions_flag         = true;
        //t->Ndisconnected_regions           = 2; // region junctions really, as it is N pairs of regions to disconnect

        // region pair one (disconnect 3 and 4 CT and RAA)
        //t->disconnect_regions[0][0]        = 2;
        //t->disconnect_regions[0][1]        = 3;

        // region pair two (disconnect 9 and 10, LAA and PV)
        //t->disconnect_regions[1][0]        = 9;
        //t->disconnect_regions[1][1]        = 10;
        // End disconnect regions =======================//|

	}
	// End Heterogeneous Human atria  ========================================//|

	// End Model-specific =============================================================//|
	if (t->Tissue_settings_set == false) // if not found, then tissue model has not been found
	{
		printf("ERROR: \"%s\" is not a valid Idealised Tissue model. Please see lib/Tissue.cpp for options\n\n", t->Tissue_model);
		exit(1);
	}
}
// End Set tissue settings from model and type - IDEALISED ======================================//|

// Set tissue settings from model and type - ANATOMICAL =========================================\\|
void set_tissue_settings_anatomical(Cell_parameters p, Tissue_parameters *t)
{
	t->S1_loc_type      = "file"; // default for anatomical is file. 
	t->S2_loc_type      = "file"; // This will be overwritten if argument passed, or set differently below
	t->map_in_type      = "file";
	t->stim_set			= true; // indicating that it doesn't need to be set by predefined loc type as with ideal
	t->S2_stim_set		= true;

    // Network model, default symmetry factors
    t->symm_fac_diag_long       = 1.0;
    t->symm_fac_diag_trans      = 1.0;
    t->symm_fac_corner          = 1.0;

    t->symm_fac_diag_long       = 0.77;//0.79;
    t->symm_fac_diag_trans      = 0.85;//0.91;
    t->symm_fac_corner          = 0.72;//0.72;

    t->apply_symmetry_factor    = "Off";

    t->Tissue_settings_set = false; 

	// Fully functional example - use as a template for new mdoels ====================\\|
	// Remember: all Optional settings are defaults only and can be overwritten by arguments
	//
	// NOTE: All geometries/fibres/stim/map files etc must be in your path: PATH/MSCSF_state_and_geometry_files/Tissue_geometries
	// If you have multiple files (e.g. different ectopic stim sites or remodelling maps) then it would be very helpful to list
	// all of them in a comment in the relevant section below.
	//
	if (strcmp(t->Tissue_model, "functional_model_test") == 0) 
	{
		// Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient 
		t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
		t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
		t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
		t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

		// Dimensions | cannot be overwritten
		t->NX       = 101;
		t->NY       = 101;
		t->NZ       = 30;

		// Files
		t->geo_file     = "functional_model_test_geo.dat";  // cannot be overwritten 
		t->stim_file    = "functional_model_test_stim.dat"; // you can set default stim to be coords and not set a file here; file can also be overwritten by argument "stim_file [filename]"
		t->S2_stim_file = "functional_model_test_stim.dat"; // can default to S1 stim file; can be overwritten by argument "S2_stim_file [filename]"
		// Geometry files should be 3D grids with 0 representing empty space and 1-N representing tissue of N celltypes
		// Stimulus files can contain anything for empty space (e.g. negative number), 0 for non-stimulus regions, and 1 for stimulus regions
		// More complex stim file for multiple sites with different timings is below in the Multi-timed stimulus section
		// End Parameters which MUST be set here ==================//|

		// Optional settings to differ from default ===============\\|
		// Default cell model
		t->Default_model    = "minimal";

		// Basic Stimulus parameters ====================\\|
		// Different S2 stimulus file to S1 file | can be overwritten by arguments
		t->S2_stim_file = "functional_model_test_S2.dat"; // overwritten by argument "S2_stim_file [filename]" 
        // Other available S2 files:
        // "functional_model_test_S2_2.dat"

		// Define S1 and/or S2 by coordiantes rather than file
		//t->S1_loc_type      = "coords";  // we want to preserve file as default here
		//t->S2_loc_type      = "coords"; 
		//t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
		//t->S2_shape         = "sphere";       // overwritten by argument "{S1/S2}_shape [sphere/cuboid]"
		// uncomment above to default to coords
		// OR if Argument "Stimulus_location_type coords" / 
		// "S2_Stimulus_location_type coords" is passed
		// the settings below will apply

		// S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
		t->S1_x_size    = 25;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
		t->S1_x_loc     = 60;                       // Location of the centre of the stimulus, x-coordiante
		t->S1_y_size    = 20;                       
		t->S1_y_loc     = 75;             
		t->S1_z_size    = 25;                        
		t->S1_z_loc     = 20;             

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ
        // End Basic Stimulus parameters ================//|

        // Heterogeneity ================================\\|
        t->Tissue_type          = "heterogeneous";
        // Celltypes 
        // These must correspond with your geoemtry file;
        // non-tissue should be denoted with a 0 
        // tissue should be denoted sequentially from 1 for the number of celltypes required
        // then set which celltype stirng reference corresponds to that number in the geoemtry file
        t->Ncelltypes         = 3;
        t->celltype_number[1] = "ENDO"; // cell 1 is vent endo
        t->celltype_number[2] = "M";
        t->celltype_number[3] = "RA";

        // Additional het option:
        // Multiple cell models || must be for each celltype, and celltype (het) must exist for that model!
        // This could be for atria and ventricles orto use two different models for different regions (e.g SAN vs atria)
        //t->Multiple_models              = "On";
        t->Tissue_model_2  = "hAM_CRN";    // This is what will be assigned to areas given "Model_2" below
        //Like Default model, this is default and can be overwritten by arguments

        t->Modeltype_number[1]          = "Model_1"; // Apply the main cell model  // example is vent
        t->Modeltype_number[2]          = "Model_1"; // vent
        t->Modeltype_number[3]          = "Model_2"; // apply cell model in Tissue_model_2 (which could be differnet to that set above if
        // overwritten by arguments); this example is is atria
        // End Heterogeneity ============================//|

        // Spatial gradient heterogeneity ===============\\|
        // Raher than just On or Off, as with other maps, 
        // here we can aslo define the default spatial
        // het model which applies to this tissue model
        // e.g. continuous apico-basal for vent, 
        // or distance from SAN for atria
        // set your het map model in lib/Model.c -> set_spatial_gradient()
        t->spatial_gradient_map_on = "Off"; // just default here; set as below if you want to apply it
        //t->spatial_gradient_map_on = "apico_basal_example"; // set your het map model in lib/Model.c -> set_spatial_gradient()
        t->spatial_gradient_map_file = "functional_model_test_spatial_gradient.dat"; // File 
        // Maps set the gradient proportion (0-1) directly
        // Maps should be 3D arrays the same dimensions as the geometry file; 
        // non-tissue can be denoted with anything, suggested a negative number
        // tissue must contain doubles 0-1 
        // End Spatial gradient heterogeneity ===========//|

        // Anisotropy ===================================\\|
        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "functional_model_test"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by arguments
        // orientation means a single file with 3 numbers (x, y, z) for each geometry location; file suffix orientation (i.e. file = {root}_orientation.dat
        // xyz means separate files for each component; suffixes OX OY OZ
        // Values are normalised components, total length = 1
        // angles means theta and phi files; suffixes theta phi; phi = angle from z; ox  = sin(theta)*cos(phi); oy  = cos(theta)*cos(phi);
        // angles_short_axis suffixes theta and phi; phi = angle from z; sc-oy = sin(theta)*cos(phi); ox  = cos(theta)*cos(phi);
        // angles in radians
        // All files must be 3D arrays of the same size as the geometry, with empty space represented by anything (e.g. negative numbers)
        // End Anisotropy ===============================//|

        // Non-uniform diffusion coefficient ============\\|
        // 1 - Non-uniformity applied by celltpye
        // Heterogeneity settings (Ncelltypes, celltype_number[x]) must be set to apply this functionality
        // Note: you can still apply homogeneous electrophysiology by setting Tissue_type to homogeneous here
        t->D_uniformity         = "regional";    

        // Celltype dependant D and anisotropy ratio - must have settings for Ncelltypes
        // These will apply the scale factor to D1 and D_AR as set globally above
        t->non_uniform_D1_scale[1]          = 0.5;//1.0;  // ENDO
        t->non_uniform_AR_scale[1]          = 1.0; // ENDO
        t->non_uniform_D1_scale[2]          = 1.0;//1.0;  // M
        t->non_uniform_AR_scale[2]          = 2.5; // M
        t->non_uniform_D1_scale[3]          = 1.5;//0.5;  // EPI
        t->non_uniform_AR_scale[3]          = 0.5; // EPI

        // OR 2 - Non-uniformity by map | all maps can be overwritten by arguments
        // Maps scale D1 and AR directly, so should contain 1 for areas to apply baseline parameters already set
        // Maps should be 3D arrays the same dimensions as the geometry file; 
        // non-tissue can be denoted with anything, suggested a negative number
        // tissue can contain any double, where that value will scale the local D or D_AR (i.e. 1.0 = baseline)
        //t->D_uniformity             = "map";     // set default here or set/overwrite with arguments
        t->Dscale_base_map_file     = "functional_model_test_Dscale_base_map.dat";      // must be set | file can contain all 1s if AR non-uniformity only is wanted 
        t->D_AR_scale_base_map_file = "functional_model_test_D_AR_scale_base_map.dat";  // must be set | file can contain all 1s if D1 non-uniformity only is wanted

        // OR 3 -Non-uniformity by map and region
        // for example, to apply an apico-basal gradient in D on top of heterogeneous values in ENDO, M and EPI 
        // Need map files AND regional settings (if not already set above, as they have been here)
        //t->D_uniformity             = "regional_map"; // set default here or set/overwrite with arguments
        //t->Dscale_base_map_file     = "functional_model_test_Dscale_base_map.dat";      // must be set | file can contain all 1s if AR non-uniformity only is wanted 
        //t->D_AR_scale_base_map_file = "functional_model_test_D_AR_scale_base_map.dat";  // must be set | file can contain all 1s if D1 non-uniformity only is wanted
        //t->non_uniform_D1_scale[1]          = 0.5;//1.0;  // ENDO
        //t->non_uniform_AR_scale[1]          = 1.0; // ENDO
        //t->non_uniform_D1_scale[2]          = 1.0;//1.0;  // M
        //t->non_uniform_AR_scale[2]          = 2.0; // M
        //t->non_uniform_D1_scale[3]          = 1.5;//0.5;  // EPI
        //t->non_uniform_AR_scale[3]          = 0.5; // EPI
        // End Non-uniform diffusion coefficient ========//|

        // Phase re-entry map(s) ========================\\|
        t->phase_file   = "functional_model_test_phasemap.dat";
        // Your phase file should be of the format with negative number for empty space, and
        // phase 0 - 2 Pi scaled to integers 0-200. 
        // Other options list here to keep track
        // End Phase re-entry map(s) ====================//|

        // Modulation maps ==============================\\|
        // Only need map files for components you want to apply by a map! Set as needed
        // All maps can be overwritten by arguments - these just provide a handy default so you don't have to pass the filename argument
        // All maps are 3D arrays of the same size as the geometry, with empty space represented by anything (e.g. negative numbers)
        // Maps contain doubles from 0 to 1 for all tissue, determining where the full magntiude of modulation applies; 
        // this magntiude is determined by the relevant condition variables (ISO, ACh, {Remodelling/Drug}_prop, Dscale, D_AR_Scale)
        // *The SRF and Direct_modulation maps, while variable type double, are essentially integeter with 1s and 0s only as these parameters cannot scale continuously
        t->ISO_map_file                 = "functional_model_test_ISO_map.dat";
        t->ACh_map_file                 = "functional_model_test_ACh_map.dat";
        t->remod_map_file               = "functional_model_test_remod_map.dat";
        t->SRF_map_file                 = "functional_model_test_SRF_map.dat";
        t->Dscale_mod_map_file          = "functional_model_test_Dscale_mod_map.dat";
        t->D_AR_scale_mod_map_file      = "functional_model_test_D_AR_scale_mod_map.dat";
        t->Direct_modulation_map_file   = "functional_model_test_Direct_modulation_map.dat";

        // In general, this will only be applied if you pass arguments to apply modulation
        // However, you could also default these maps to be on for a specific reused model
        // You can't default the magntiude conditon variable (ISO, ACh, {Remodelling/Drug}_prop) here 
        // as this is set in other parts of the model
        //t->ISO_map_on		= "On";			
        //t->remod_map_on		= "On";		
        //t->ACh_map_on		= "On";			
        //t->SRF_map_on		= "On";			
        //t->Direct_modulation_map_on		= "On";			

        // For Dscale and D_AR_scale, you can set a default value of the magntiude here (which can still be overwritten by arguments)
        //t->Dscale_map_on        = "On";  
        //t->D_AR_scale_map_on    = "On";   
        //t->Dscale			    = 0.5;  // multiplies local D1 (set by global D1 or regional D_uniformity) in map region
        //t->D_AR_scale           = 2.0;  // multiples local D_AR in map region

        // or you can apply modulation maps by coordinates
        // All maps are either file or coordinates; combinations are not possible
        // If using coordinates, the same map applies to all modulations
        //t->map_in_type          = "coords"; // can set default to coords, or leave default as map and have coords set if argument passed
        // Set coordinates and shape of map; this one applies to right half of tissue
        t->ideal_map_x_loc      = t->NX - int(float(t->NX/4.0));
        t->ideal_map_y_loc      = int(float(t->NY/2.0));
        t->ideal_map_z_loc      = int(float(t->NZ/2.0));
        t->ideal_map_x_size     = int(float(t->NX/4.0));
        t->ideal_map_y_size     = int(float(t->NY/2.0));
        t->ideal_map_z_size     = int(float(t->NZ/2.0));
        t->ideal_map_shape      = "cuboid"; // default is sphere
        t->ideal_map_set        = true; // so it isn't overwritten by default
        // End modulation maps ===========================//|

        // Multi-timed stimulus ==========================\\|
        // Stim  map should be anything empty space;
        // 1 -> Nstims where we want stimuli to apply, where the number denotes the sequence (1 = first, Nstim = last)
        //t->Multi_stim   = "On"; // can default to Multi_timed sim being on, don't have to
        // Availabe file for multi stim: t->stim_file    = "functional_model_test_multi_stim.dat";	
        t->Nstims       = 5; // number of differently timed stimuli
        // delay settings || Need Nstims-1 entries here
        // stim file element = 1 has no stim delay as this is intial site
        // stim_delay[1] applies delay to second stim area (stim file element = 2)
        // delay is in ms relative to first stimulus;
        t->stim_delay[1] = 6;
        t->stim_delay[2] = 10;
        t->stim_delay[3] = 12;
        t->stim_delay[4] = 25;
        // End Multi-timed stimulus ======================//|
        // End Optional settings to differ from default ===========//|
    }
    // End Fully functional example - use as a template for new mdoels ================//|

    // Human ventricular wedge reconstruction =========================================\\|
    // Benson et al. 2007 Chaos 17(1):015105 
    if ((strcmp(t->Tissue_model, "Human_vent_wedge") == 0))
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Settings
        t->Tissue_type          = "heterogeneous";
        t->Orientation_type     = "isotropic";
        t->D_uniformity         = "uniform";

        // Default cell model
        t->Default_model    = "minimal";

        // Files
        t->geo_file     = "Human_vent_wedge_geo.dat";
        t->stim_file    = "Human_vent_wedge_stim.dat";

        // Dimensions
        t->NX       = 102;
        t->NY       = 102;
        t->NZ       = 102;

        // Diffusion properties
        t->dx = t->dy   = 0.2125;
        t->dz           = 0.25;    // mm
        t->D1           = 0.1171;
        t->D_AR         = 9;               // anisotropy ratio (D1/D2)
        t->Diso         = 0.04;            // for isotropic
        t->Gl           = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt           = t->Gl/t->D_AR;

        // Celltypes
        t->Ncelltypes         = 3;
        t->celltype_number[1] = "ENDO"; // cell 1 is vent endo
        t->celltype_number[2] = "M";
        t->celltype_number[3] = "EPI";

        // Anisotropy
        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Human_vent_wedge"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "xyz";
    }
    // End Human ventricular wedge reconstruction =====================================//|

    // Whole canine ventricle reconstruction ==========================================\\|
    // Benson et al. 2008 Progress in Biophysics and Molecular Biology, 
    // Cardiovascular Physiome, 96(1): 187-208
    if ((strcmp(t->Tissue_model, "Canine_vent") == 0))
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // Settings
        t->Tissue_type          = "homogeneous";
        t->Orientation_type     = "isotropic";
        t->D_uniformity         = "uniform";

        // Default cell model
        t->Default_model    = "minimal";

        // Files
        t->geo_file     = "Canine_vent_geo.dat";
        t->stim_file    = "Canine_vent_multi_stim.dat";

        // Dimensions
        t->NX       = 128;
        t->NY       = 128;
        t->NZ       = 115;

        // Diffusion properties
        t->dx = t->dy = t->dz = 0.78;    // mm
        t->Diso     = 0.1;  // for isotropic
        t->D1       = 0.576;
        t->D_AR     = 9;
        t->Gl       = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt       = t->Gl/t->D_AR;

        // Multi-stim settings
        t->Multi_stim           = "On";
        t->Nstims               = 3;
        t->stim_delay[1]        = 5;//10;
        t->stim_delay[2]        = 10;//30;

        // Anisotropy
        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Canine_vent"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "xyz";
    }
    // End Whole canine ventricle reconstruction ======================================//|

    // Vector fibre and fibrosis network geos==========================================\\|
    if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_vis") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        // Dimensions | cannot be overwritten
        t->NX       = 20;
        t->NY       = 20;
        t->NZ       = 1;

        // Files
        t->geo_file                 = "Fibrosis_geometry_20_20_1.dat";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 5;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 5;
        t->S1_y_loc     = int(float(t->NY/2.0));
        t->S1_z_size    = 5;
        t->S1_z_loc     = int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Fibrosis_geometry_20_20_1"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi
    }

    if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_vis_2") == 0 || strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_vis_2_X") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        // Dimensions | cannot be overwritten
        t->NX       = 50;
        t->NY       = 50;
        t->NZ       = 1;

        // Files
        t->geo_file                 = "Fibrosis_geometry_50_50_1.dat";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 5;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 5;
        t->S1_y_loc     = int(float(t->NY/2.0));
        t->S1_z_size    = 5;
        t->S1_z_loc     = int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_vis_2") == 0) t->orientation_file_root            = "Fibrosis_geometry_50_50_1"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        else if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_vis_2_X") == 0) t->orientation_file_root     = "Fibrosis_geometry_50_50_1_X";
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi

    }

    if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_2D_sims") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        // Dimensions | cannot be overwritten
        t->NX       = 200;
        t->NY       = 200;
        t->NZ       = 1;

        // Files
        t->geo_file                 = "Fibrosis_geometry_200_200_1.dat";  // cannot be overwritten
        t->junction_mod_map_file    = "Fibrosis_geometry_200_200_1_junction_map_0.2_0.8.dat";       // can be overwritten and only get used if map on arg is passed
        t->junction_base_map_file   = "Fibrosis_geometry_200_200_1_junction_baseline_het_map.dat";

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 10;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 10;
        t->S1_y_loc     = int(float(t->NY/2.0));
        t->S1_z_size    = 10;
        t->S1_z_loc     = int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Fibrosis_geometry_200_200_1"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi
    }

    if (strcmp(t->Tissue_model, "Vector_field_and_fibrosis_3D") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        // Dimensions | cannot be overwritten
        t->NX       = 100;
        t->NY       = 100;
        t->NZ       = 100;

        // Files
        t->geo_file     = "Fibrosis_geometry_100_100_100.dat";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 5;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 5;
        t->S1_y_loc     = int(float(t->NY/2.0));
        t->S1_z_size    = 5;
        t->S1_z_loc     = int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Fibrosis_geometry_100_100_100"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi
    }
    // End Vector fibre and fibrosis network geos======================================//|
    
    // 300x300 2D sheet for fibrosis paper ============================================\\|
    // This is a homogeneous sheet designed for atria (diffusion params set for that AR)
    if (strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_OY") == 0 || strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_field_control") == 0 || strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_field_remodelled") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.3125;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz  -> to match 3D atria from Auckland for directly comparable simulations
        t->D1   = 0.4;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 7;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;

        // Dimensions | cannot be overwritten
        t->NX       = 300;
        t->NY       = 300;
        t->NZ       = 1;

        // Files
        t->geo_file                 = "2D_human_atria_fibrosis_300x300.txt";  // cannot be overwritten
        //t->junction_mod_map_file    = "Fibrosis_geometry_200_200_1_junction_map_0.2_0.8.dat";       // can be overwritten and only get used if map on arg is passed
        //t->junction_base_map_file   = "Fibrosis_geometry_200_200_1_junction_baseline_het_map.dat";

        // map file options (fibrosis and SRF maps)
        //2D_human_atria_fibrosis_map_ln_10_perc_25.txt
        //2D_human_atria_fibrosis_map_ln_10_perc_50.txt
        //2D_human_atria_fibrosis_map_ln_25_perc_25.txt
        //2D_human_atria_fibrosis_map_ln_25_perc_50.txt
        //2D_human_atria_fibrosis_map_ln_37_perc_25.txt
        //2D_human_atria_fibrosis_map_ln_37_perc_50.txt
        //2D_human_atria_fibrosis_map_ln_50_perc_25.txt
        //2D_human_atria_fibrosis_map_ln_50_perc_50.txt

        //2D_human_atria_fibrosis_300x300_OY_JNmap_10_25.txt etc for all of the above and the three orientations

        // Default cell model
        t->Default_model    = "minimal";

        // Settings
        t->Tissue_type          = "heterogeneous";
        t->Orientation_type     = "anisotropic";
        t->D_uniformity         = "uniform";

        // Celltypes
        t->Ncelltypes         = 1; // only 1 region, but set to heterogeneous to default this region to RA for simpicity
        t->celltype_number[1] = "RA"; // cell 1 is vent endo

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape     = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 10;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 10;
        t->S1_y_loc     = int(float(t->NY/2.0));
        t->S1_z_size    = 10;
        t->S1_z_loc     = int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        // Fibre files
        if (strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_OY") == 0)
            t->orientation_file_root    = "2D_human_atria_fibrosis_300x300_OY"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        if (strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_field_control") == 0)
            t->orientation_file_root    = "2D_human_atria_fibrosis_300x300_field_control"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        if (strcmp(t->Tissue_model, "2D_human_atria_fibrosis_300x300_field_remodelled") == 0)
            t->orientation_file_root    = "2D_human_atria_fibrosis_300x300_field_remodelled"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi
    }
    // End 300x300 2D sheet for fibrosis paper ========================================//|

    // Rat three eigenvector gro test =================================================\\|
    // Model from: Whittaker DG, Benson AP, Teh I, Schneider JE, Colman MA 
    // "Investigation of the Role of Myocyte Orientations in Cardiac Arrhythmia Using Image-Based Models" Biophysical Journal 2019, 117 (12), 2396-2408 
    if (strcmp(t->Tissue_model, "Rat_vent_DTI_2") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.1;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;
        t->Gt2  = t->Gt/2.0;

        // Dimensions | cannot be overwritten
        t->NX       = 105;
        t->NY       = 115;
        t->NZ       = 115;

        // Files
        t->geo_file                 = "Rat_vent_DTI_2.dat";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 20;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = 55;//int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 20;
        t->S1_y_loc     = 55;//int(float(t->NY/2.0));
        t->S1_z_size    = 20;
        t->S1_z_loc     = 10;//int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Rat_vent_DTI_2"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "orientation"; // or xyz or angles | cannot be overwritten by argumentsi
    }
    // End Rat three eigenvector gro test =============================================//|

    // Rat DTI Benson 2011 ============================================================\\|
    // ADD citation
    if (strcmp(t->Tissue_model, "Rat_vent_DTI_Benson_2011") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;//0.4;//0.17;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;
        t->Gt2  = t->Gt/2.0;

        // Dimensions | cannot be overwritten
        t->NX       = 70;
        t->NY       = 85;
        t->NZ       = 75;

        // Files
        t->geo_file                 = "Rat_vent_DTI_Benson_2011.txt";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 10;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = 48;//int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_loc     = 62;//int(float(t->NY/2.0));
        t->S1_z_loc     = 31;//int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "anisotropic";
        t->orientation_file_root    = "Rat_vent_DTI_Benson_2011"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "xyz"; // or xyz or angles | cannot be overwritten by arguments

        t->D_uniformity             = "map";     // set default here or set/overwrite with arguments
        t->Dscale_base_map_file     = "Rat_vent_DTI_Benson_Dbaseline.txt";      // must be set | file can contain all 1s if AR non-uniformity only is wanted
        t->D_AR_scale_base_map_file = "Rat_vent_DTI_Benson_DAR_baseline.txt";
    }
    // End Rat DTI Benson 2011 ========================================================//|

    // Rat bi-ventricular geo 2  ======================================================\\|
    // Model from: Whittaker DG, Benson AP, Teh I, Schneider JE, Colman MA
    // "Investigation of the Role of Myocyte Orientations in Cardiac Arrhythmia Using Image-Based Models" Biophysical Journal 2019, 117 (12), 2396-2408
    if (strcmp(t->Tissue_model, "Rat_vent_three_eigenvectors_control_base") == 0)
    {
        // Parameters which MUST be set here ======================\\|
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        t->dx = t->dy = t->dz = 0.2;    // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 0.2;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 4;                    // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"

        // gl and gt
        t->Gl   = t->D1/t->dx;  // nS/pF; conversion from D to g
        t->Gt   = t->Gl/t->D_AR;
        t->Gt2  = t->Gt/2.0;

        // Dimensions | cannot be overwritten
        t->NX       = 63;
        t->NY       = 54;
        t->NZ       = 60;

        // Files
        t->geo_file                 = "Rat_vent_three_eigenvectors_control_base.dat";  // cannot be overwritten

        // Default cell model
        t->Default_model    = "minimal";

        // S1 || each setting can be individually overwritten by arguments "{S1/S2}_{x/y/x}_{loc/size} [n]"
        t->S1_shape         = "sphere";     // sphere or cuboid || default is cuboid
        t->S1_loc_type  = "coords";
        t->S1_x_size    = 5;                        // Extent of stimulus around location centre; example Total size = 5+5+1 = 11
        t->S1_x_loc     = 20;//int(float(t->NX/2.0));    // Location of the centre of the stimulus, x-coordiante
        t->S1_y_size    = 5;
        t->S1_y_loc     = 20;//int(float(t->NY/2.0));
        t->S1_z_size    = 5;
        t->S1_z_loc     = 20;//int(float(t->NZ/2.0));

        // S2 - setup for cross-field
        t->S2_x_loc     = int(float(t->NX/2.0));    // half of y, all of x and z
        t->S2_x_size    = int(float(t->NX/2.0));
        t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
        t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
        t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
        t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ

        t->Orientation_type         = "three_eigenvectors";
        t->orientation_file_root    = "Rat_vent_three_eigenvectors_control_base"; // cannot be overwritten by arguments. Filnames are {root}_suffix.dat
        t->orientation_file_type    = "xyz"; // or xyz or angles | cannot be overwritten by argumentsi
    }
    // End Rat bi-ventricular geo 2 ===================================================//|

    // Human whole heart Halina Auckland ==============================================\\|
    if (strcmp(t->Tissue_model, "Human_whole_heart_1_test") == 0 || strcmp(t->Tissue_model, "Human_whole_heart_1_geo_vis") == 0 || strcmp(t->Tissue_model, "Human_whole_heart_MI") == 0)
    {
        t->Tissue_settings_set = true; // now the code knows the settings have been found

        // space step and global/baseline diffusion coefficient
        // Control values:
        //t->dx = t->dy = t->dz = 0.219; // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        //t->dx = t->dy = t->dz = 0.365; // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->dx = t->dy = t->dz = 0.292; // mm | can be overwritten by argument "dx [x]", but only for dx=dy=dz
        t->D1   = 1;                  // mm^2/ms | can be overwritten by argument "D1 [x]" and modulated by argument "Dscale [x]"
        t->D_AR = 10;                   // anisotropy ratio (D1/D2) | can be overwrriten by argument "D_AR [x]" and modulated by "D_AR_scale [x]"
        t->Diso = 0.3;                  // mm^2/ms | can be overwritten by argument "D1 [x]" (if anisotropy = isotropic) and modulated by "Dscale [x]"
        t->Gl   = t->D1/t->dx;          // nS/pF; conversion from D to g; NETWORK model
        t->Gt   = t->Gl/t->D_AR;        // nS/p

        // um 219 dimensions
        /*t->NX       = 558
          t->NY       = 471;
          t->NZ       = 350;*/

        // um 365 dimensions
        //t->NX       = 335;
        //t/->NY       = 283;
        //t->NZ       = 210;

        // um 292 dimensions
        t->NX       = 419;
        t->NY       = 353;
        t->NZ       = 263;

        if (strcmp(t->Tissue_model, "Human_whole_heart_MI") == 0)
        {
            t->NX       = 251;
            t->NY       = 311;
            t->NZ       = 374;
        }

        // Default cell model
        t->Default_model    = "minimal";

        // Files
        t->geo_file     = "Human_whole_heart_control_292um.txt";//"AVN_Bundle_Branches_only.txt";
        t->stim_file    = "Human_whole_heart_control_292um_SAN_stim.txt";
        if (strcmp(t->Tissue_model, "Human_whole_heart_1_geo_vis") == 0) t->geo_file = "Human_whole_heart_control_292um_VIS_geo.txt";
        if (strcmp(t->Tissue_model, "Human_whole_heart_MI") == 0) t->geo_file = "Human_whole_heart_MI.txt";
        if (strcmp(t->Tissue_model, "Human_whole_heart_MI") == 0) t->stim_file = "Human_whole_heart_MI_SAN_stim.txt";

        t->Orientation_type         = "isotropic";
        //t->Tissue_type              = "homogeneous";
        // Heterogeneity ================================\\|
        t->Tissue_type          = "heterogeneous";
        t->Ncelltypes           = 10;
        t->celltype_number[1]   = "Pace_CCS"; // SAN, non-coupled
        t->celltype_number[2]   = "Pace_CCS"; // SAN, coupled to atria
        t->celltype_number[3]   = "RA_CCS"; // AM
        t->celltype_number[4]   = "Pace_CCS"; // AVN, non coupled
        t->celltype_number[5]   = "Pace_CCS"; // AVN, coupled
        t->celltype_number[6]   = "PK_CCS"; // RPK, non coupled
        t->celltype_number[7]   = "PK_CCS"; // RPK, coupled
        t->celltype_number[8]   = "PK_CCS"; // LPK, non coupled
        t->celltype_number[9]   = "PK_CCS"; // LPK, coupled
        t->celltype_number[10]  = "EPI_CCS"; // VM
        // End Heterogeneity ============================//|

        // Disconnect regions ===========================\\|
        t->disconnect_regions_flag         = true;
        t->Ndisconnected_regions           = 9; // region junctions really, as it is N pairs of regions to disconnect

        // region pair one (disconnect 3 and 10, atria and ventricles)
        t->disconnect_regions[0][0]        = 3;
        t->disconnect_regions[0][1]        = 10;
        
        // region pair two (disconnect 5 and 10, AVN and ventricles)
        t->disconnect_regions[1][0]        = 5;
        t->disconnect_regions[1][1]        = 10;
        
        // region pair three (disconnect 1 and 3, SAN non couple and AM)
        t->disconnect_regions[2][0]        = 1;
        t->disconnect_regions[2][1]        = 3;

        // region pair four (disconnect 6 and 10, RPK non couple and VM)
        t->disconnect_regions[3][0]        = 6;
        t->disconnect_regions[3][1]        = 10;
        
        // region pair five (disconnect 7 and 10, LPK non couple and VM)
        t->disconnect_regions[4][0]        = 7;
        t->disconnect_regions[4][1]        = 10;
        
        // region pair four (disconnect 6 and 10, RPK non couple and AM)
        t->disconnect_regions[5][0]        = 6;
        t->disconnect_regions[5][1]        = 3;
        
        // region pair five (disconnect 7 and 10, LPK non couple and AM)
        t->disconnect_regions[6][0]        = 7;
        t->disconnect_regions[6][1]        = 3;
        
        // region pair four (disconnect 6 and 10, RPK couple and AM)
        t->disconnect_regions[7][0]        = 8;
        t->disconnect_regions[7][1]        = 3;
        
        // region pair five (disconnect 7 and 10, LPK couple and AM)
        t->disconnect_regions[8][0]        = 9;
        t->disconnect_regions[8][1]        = 3;
        // End disconnect regions =======================//|

        // 1 - Non-uniformity applied by celltpye
        // Heterogeneity settings (Ncelltypes, celltype_number[x]) must be set to apply this functionality
        // Note: you can still apply homogeneous electrophysiology by setting Tissue_type to homogeneous here
        t->D_uniformity         = "regional";
        //t->D_uniformity         = "uniform";


        // Celltype dependant D and anisotropy ratio - must have settings for Ncelltypes
        // These will apply the scale factor to D1 and D_AR as set globally above
        t->non_uniform_D1_scale[1]          = 0.33;     // SAN 
        t->non_uniform_AR_scale[1]          = 1.0;      // SAN
        t->non_uniform_D1_scale[2]          = 0.33;     // SAN 
        t->non_uniform_AR_scale[2]          = 1.0;      // SAN
        t->non_uniform_D1_scale[3]          = 1.0;      // AM
        t->non_uniform_AR_scale[3]          = 1.0;      // AM
        t->non_uniform_D1_scale[4]          = 0.33;     // AVN
        t->non_uniform_AR_scale[4]          = 1.0;      // AVN
        t->non_uniform_D1_scale[5]          = 0.33;     // AVN
        t->non_uniform_AR_scale[5]          = 1.0;      // AVN
        t->non_uniform_D1_scale[6]          = 1.0;      // PK
        t->non_uniform_AR_scale[6]          = 1.0;      // PK
        t->non_uniform_D1_scale[7]          = 1.0;      // PK
        t->non_uniform_AR_scale[7]          = 1.0;      // PK
        t->non_uniform_D1_scale[8]          = 1.0;      // PK
        t->non_uniform_AR_scale[8]          = 1.0;      // PK
        t->non_uniform_D1_scale[9]          = 1.0;      // PK
        t->non_uniform_AR_scale[9]          = 1.0;      // PK
        t->non_uniform_D1_scale[10]         = 1.0;      // VM
        t->non_uniform_AR_scale[10]         = 1.0;      // VM

        // CT, BB and PM regions want D1 scale 2 (0.3 -> 0.6)
    }
    // End human whole heart Halina Auckland ==========================================//|

    if (t->Tissue_settings_set == false)
    {
        printf("ERROR: \"%s\" is not a valid Tissue model for geo. Please see lib/Tissue.cpp for options\n\n", t->Tissue_model);
        exit(1);
    }

}
// End Set tissue settings from model and type - ANATOMICAL =====================================//|

// Set idealised coordiante stimulus parameters =================================================\\|
void set_coord_stim_and_map_from_defined_type(Cell_parameters p, Tissue_parameters *t, Argument_parameters A)
{
    // These need to be here to unset "stim_set" if it has been set to true by the tissue model defaults but we want to overwrite whole shape with argument
    if (A.Stimulus_type_arg 	== true)    
    {
        t->S1_loc_type      = A.Stimulus_loc_type;
        t->stim_set         = false;
    }
    if (A.S2_Stimulus_type_arg 	== true)    
    {
        t->S2_loc_type      = A.S2_Stimulus_loc_type;
        t->S2_stim_set      = false;
    }

    // General/global stimulus settings =====================================\\|
    if (t->stim_set == false)
    {
        if (strcmp(t->S1_loc_type, "edge") == 0)    // Stimulus from x -> 0 edge
        {
            t->S1_x_size    = 5;                        // Total size = 5+5+1 = 11
            t->S1_x_loc     = t->S1_x_size;             // S1_size+1 th cell (from 0)
            t->S1_y_loc     = int(float(t->NY/2.0));    // centre in y
            t->S1_y_size    = int(float(t->NY/2.0));    // NY/2 + NY/2 + 1 >= NY
            t->S1_z_loc     = int(float(t->NZ/2.0));    // centre in z
            t->S1_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ
            t->stim_set 	= true;
        }
        else if (strcmp(t->S1_loc_type, "centre") == 0)    // Stimulus from centre
        {
            // Size 5+5+1 = 11, centre each dimension
            t->S1_x_loc     = int(float(t->NX/2.0));
            t->S1_x_size    = 7;//12;
            t->S1_y_loc     = int(float(t->NY/2.0));
            t->S1_y_size    = 7;//12;
            t->S1_z_loc     = int(float(t->NZ/2.0));
            t->S1_z_size    = 7;//12;
            t->stim_set 	= true;
        }
        else if (strcmp(t->S1_loc_type, "whole_tissue") == 0)    // Stimulus from x -> 0 edge
        {
            t->S1_x_size    = int(float(t->NX/2.0));    
            t->S1_x_loc     = int(float(t->NX/2.0));             
            t->S1_y_loc     = int(float(t->NY/2.0));    // centre in y
            t->S1_y_size    = int(float(t->NY/2.0));    // NY/2 + NY/2 + 1 >= NY
            t->S1_z_loc     = int(float(t->NZ/2.0));    // centre in z
            t->S1_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ
            t->stim_set     = true;
        }
        else if (strcmp(t->S1_loc_type, "coords") == 0 || strcmp(t->S1_loc_type, "file") == 0) // This is for anatomically detailed models
        {
            t->stim_set     = true;
        }
    }

    // Default S2 = S1, centre or edge-> can be set in model specific OR arguments -> will only set these if not set elsewhere
    // Arguments will only overwrite the aspect defined by the arg
    if (t->S2_stim_set == false)
    {
        if (strcmp(t->S2_loc_type, "S1") == 0)    // same as S1
        {
            t->S2_x_loc		= t->S1_x_loc;
            t->S2_x_size	= t->S1_x_size;
            t->S2_y_loc		= t->S1_y_loc;
            t->S2_y_size	= t->S1_y_size;
            t->S2_z_loc		= t->S1_z_loc;
            t->S2_z_size	= t->S1_z_size;
            t->S2_stim_set 	= true;	
        }
        else if (strcmp(t->S2_loc_type, "cross_field") == 0)    // cross-field
        {
            t->S2_x_loc     = int(float(t->NX/2.0));	// half of y, all of x and z
            t->S2_x_size    = int(float(t->NX/2.0));
            t->S2_y_loc     = int(float(t->NY/4.0));    // centre in y
            t->S2_y_size    = int(float(t->NY/4.0));    // NY/2 + NY/2 + 1 >= NY
            t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
            t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ
            t->S2_stim_set  = true;
        }
        else if (strcmp(t->S2_loc_type, "edge") == 0)    // Stimulus from x -> 0 edge
        {
            t->S2_x_size    = 5;                        // Total size = 5+5+1 = 11
            t->S2_x_loc     = t->S2_x_size;             // S2_size+1 th cell (from 0)
            t->S2_y_loc     = int(float(t->NY/2.0));    // centre in y
            t->S2_y_size    = int(float(t->NY/2.0));    // NY/2 + NY/2 + 1 >= NY
            t->S2_z_loc     = int(float(t->NZ/2.0));    // centre in z
            t->S2_z_size    = int(float(t->NZ/2.0));    // NZ/2 + NZ/2 + 1 >= NZ
            t->S2_stim_set  = true;
        }
        else if (strcmp(t->S2_loc_type, "centre") == 0)    // Stimulus from centre
        {
            // Size 5+5+1 = 11, centre each dimension
            t->S2_x_loc     = int(float(t->NX/2.0));
            t->S2_x_size    = 7;
            t->S2_y_loc     = int(float(t->NY/2.0));
            t->S2_y_size    = 7;
            t->S2_z_loc     = int(float(t->NZ/2.0));
            t->S2_z_size    = 7;
            t->S2_stim_set  = true;
        }
        else if (strcmp(t->S2_loc_type, "coords") == 0) // This is for anatomically detailed models
        {
            t->S2_stim_set  = true;
        }
    }
    // End General/global stimulus settings =================================//|

    // Set default map (ISO, Dscale, remod, ACh, SRF) coords ================\\|
    // This creates a patch in the centre of the tissue
    if (t->ideal_map_set == false)
    {
        t->ideal_map_x_loc      = int(float(t->NX/2.0));
        t->ideal_map_y_loc      = int(float(t->NY/2.0));
        t->ideal_map_z_loc      = int(float(t->NZ/2.0));
        t->ideal_map_x_size     = int(float(t->NX/4.0));
        t->ideal_map_y_size     = int(float(t->NY/4.0));
        t->ideal_map_z_size     = int(float(t->NZ/4.0));
        t->ideal_map_set        = true;
    }
    // End Set default map (ISO, Dscale, remod, ACh, SRF) coords ============//|

    if (t->stim_set == false)
    {
        printf("ERROR: \"%s\" is not a valid Stimulus model for tissue model %s. Please see lib/Tisseue.cpp for options\n", t->S1_loc_type, t->Tissue_model);
        printf("(If this is unexpected, and you have set model specific stim type, have you set \"stim_set = true\" where settings are defined?)\n\n");
        exit(1);
    }
} // end set tissue settings
  // End Set idealised coordiante stimulus parameters =============================================//|

  // Set global fibre orientation from arguments  =================================================\\|
void set_global_orientation_direction_from_arg(Cell_parameters p, Tissue_parameters *t, Argument_parameters A)
{
    if (A.Global_orientation_direction_arg == true) 
    {
        t->Global_orientation_direction = A.Global_orientation_direction;

        if (strcmp(t->Tissue_order, "2D") == 0)
        {
            if (strcmp(t->Global_orientation_direction, "X") == 0)
            {
                t->OX = 1.0;
                t->OY = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "Y") == 0)
            {
                t->OX = 0.0;
                t->OY = 1.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XY_plus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OY = sqrt(0.5);
            }
            else if (strcmp(t->Global_orientation_direction, "XY_minus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OY = -sqrt(0.5);
            }
            else 
            {
                printf("ERROR: Global_orientation_direction %s is invalid for 2D models; must be X / Y / XY_plus / XY_minus\n", t->Global_orientation_direction);
                exit(1);
            }
        }

        else if (strcmp(t->Tissue_order, "3D") == 0)
        {
            if (strcmp(t->Global_orientation_direction, "X") == 0)
            {
                t->OX = 1.0;
                t->OY = 0.0;
                t->OZ = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "Y") == 0)
            {
                t->OX = 0.0;
                t->OY = 1.0;
                t->OZ = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "Z") == 0)
            {
                t->OX = 0.0;
                t->OY = 0.0;
                t->OZ = 1.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XY_plus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OY = sqrt(0.5);
                t->OZ = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XY_minus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OY = -sqrt(0.5);
                t->OZ = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XZ_plus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OZ = sqrt(0.5);
                t->OY = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XZ_minus") == 0)
            {
                t->OX = sqrt(0.5);
                t->OZ = -sqrt(0.5);
                t->OY = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "YZ_plus") == 0)
            {
                t->OY = sqrt(0.5);
                t->OZ = sqrt(0.5);
                t->OX = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "YZ_minus") == 0)
            {
                t->OY = sqrt(0.5);
                t->OZ = -sqrt(0.5);
                t->OX = 0.0;
            }
            else if (strcmp(t->Global_orientation_direction, "XYZ_ppp") == 0)
            {
                t->OX = sqrt(1.0/3.0);
                t->OY = sqrt(1.0/3.0);
                t->OZ = sqrt(1.0/3.0);
            }
            else if (strcmp(t->Global_orientation_direction, "XYZ_ppm") == 0)
            {
                t->OX = sqrt(1.0/3.0);
                t->OY = sqrt(1.0/3.0);
                t->OZ = -sqrt(1.0/3.0);
            }
            else if (strcmp(t->Global_orientation_direction, "XYZ_pmp") == 0)
            {
                t->OX = sqrt(1.0/3.0);
                t->OY = -sqrt(1.0/3.0);
                t->OZ = sqrt(1.0/3.0);
            }
            else if (strcmp(t->Global_orientation_direction, "XYZ_mpp") == 0)
            {
                t->OX = -sqrt(1.0/3.0);
                t->OY = sqrt(1.0/3.0);
                t->OZ = sqrt(1.0/3.0);
            }
            else 
            {
                printf("ERROR: Global_orientation_direction %s is invalid for 3D models; must be X / Y / Z / {XY/XZ/YZ}_plus / {XY/XZ/YZ}_minus / XYZ_{ppp/ppm/pmp/mpp}\n", t->Global_orientation_direction);
                exit(1);
            }
        }
    }
}
// End Set global fibre orientation from arguments  =============================================//|

// overwrite properties from arguments ==========================================================\\|
void overwrite_tissue_properties_from_args(Cell_parameters p, Tissue_parameters *t, Argument_parameters A)
{
    // Diffusion parameters
    if (A.D1_arg == true)       t->D1   = A.D1;
    if (A.D_AR_arg == true)     t->D_AR = A.D_AR;
    if (A.D1_arg == true)       t->Diso = A.D1; // as D1 set to Diso if isotropic
    if (A.dx_arg == true)       t->dx   = t->dy = t->dz = A.dx;
    if (A.Dscale_arg == true)			t->Dscale			= A.Dscale;
    if (A.D_AR_scale_arg == true)   	t->D_AR_scale		= A.D_AR_scale;
    //if (A.D1_arg == true)       t->Gl   = A.D1;
    //if (A.D2_arg == true)       t->Gt   = A.D2;

    // NETWORK MODEL
    if (A.D1_arg == true)
    {
        t->Gl = t->D1/t->dx;
        t->Gt = t->Gl/t->D_AR;
    }
    if (A.D_AR_arg == true)
    {
        t->Gt = t->Gl/t->D_AR;
    }

    // Idealised model fibre orientation
    if (A.OX_arg == true)       t->OX   = A.OX;
    if (A.OY_arg == true)       t->OY   = A.OY;
    if (A.OZ_arg == true)       t->OZ   = A.OZ;

    // Ensure normalised fibres / set final orientaiton from others
    if (strcmp(t->Tissue_order, "2D") == 0)
    {
        t->OZ = 0;
        if (A.OX_arg == true && A.OY_arg == false) t->OY = sqrt(1 - t->OX*t->OX);
        if (A.OX_arg == false && A.OY_arg == true) t->OX = sqrt(1 - t->OY*t->OY);
    }
    else
    {
        if (A.OX_arg == true && A.OY_arg == true && A.OZ_arg == false) t->OZ = sqrt(1 - t->OX*t->OX - t->OY*t->OY);
        if (A.OX_arg == true && A.OY_arg == false && A.OZ_arg == true) t->OY = sqrt(1 - t->OX*t->OX - t->OZ*t->OZ);
        if (A.OX_arg == false && A.OY_arg == true && A.OZ_arg == true) t->OX = sqrt(1 - t->OZ*t->OZ - t->OY*t->OY);
    }

    // Stimulus properties
    if (A.S1_shape_arg  == true)    t->S1_shape     = A.S1_shape;
    if (A.S2_shape_arg  == true)    t->S2_shape     = A.S2_shape;

    if (A.S1_x_loc_arg  == true)    t->S1_x_loc     = A.S1_x_loc;
    if (A.S1_x_size_arg == true)    t->S1_x_size    = A.S1_x_size;
    if (A.S1_y_loc_arg  == true)    t->S1_y_loc     = A.S1_y_loc;
    if (A.S1_y_size_arg == true)    t->S1_y_size    = A.S1_y_size;
    if (A.S1_z_loc_arg  == true)    t->S1_z_loc     = A.S1_z_loc;
    if (A.S1_z_size_arg == true)    t->S1_z_size    = A.S1_z_size;

    if (A.S2_x_loc_arg  == true)    t->S2_x_loc     = A.S2_x_loc;
    if (A.S2_x_size_arg == true)    t->S2_x_size    = A.S2_x_size;
    if (A.S2_y_loc_arg  == true)    t->S2_y_loc     = A.S2_y_loc;
    if (A.S2_y_size_arg == true)    t->S2_y_size    = A.S2_y_size;
    if (A.S2_z_loc_arg  == true)    t->S2_z_loc     = A.S2_z_loc;
    if (A.S2_z_size_arg == true)    t->S2_z_size    = A.S2_z_size;

    if (A.Multi_stim_arg == true)   t->Multi_stim   = A.Multi_stim;

    // Tissue model conditions
    if (A.Tissue_type_arg 		== true)    t->Tissue_type      = A.Tissue_type;
    if (A.Orientation_type_arg	== true)    t->Orientation_type = A.Orientation_type;
    if (A.D_uniformity_arg 		== true)    t->D_uniformity     = A.D_uniformity;
    if (A.Dscale_mod_map_arg    == true)    t->Dscale_map_on    = A.Dscale_mod_map_on;
    if (A.D_AR_scale_mod_map_arg    == true)    t->D_AR_scale_map_on    = A.D_AR_scale_mod_map_on;
    if (A.ISO_map_arg 			== true)    t->ISO_map_on       = A.ISO_map_on;
    if (A.Remodelling_map_arg 	== true)    t->remod_map_on     = A.Remodelling_map_on;
    if (A.ACh_map_arg 			== true)    t->ACh_map_on       = A.ACh_map_on;
    if (A.SRF_map_arg 			== true)    t->SRF_map_on       = A.SRF_map_on;
    if (A.Direct_modulation_map_arg 			== true)    t->Direct_modulation_map_on       = A.Direct_modulation_map_on;
    if (A.stim_file_arg 		== true)    t->stim_file        = A.stim_file;
    if (A.S2_stim_file_arg 		== true)    t->S2_stim_file     = A.S2_stim_file;
    if (A.phase_file_arg 		== true)    t->phase_file       = A.phase_file;
    if (A.ISO_map_file_arg 		== true)    t->ISO_map_file     = A.ISO_map_file;
    if (A.remod_map_file_arg 	== true)    t->remod_map_file   = A.remod_map_file;
    if (A.ACh_map_file_arg 		== true)    t->ACh_map_file     = A.ACh_map_file;
    if (A.SRF_map_file_arg 		== true)    t->SRF_map_file     = A.SRF_map_file;
    if (A.Direct_modulation_map_file_arg 		== true)    t->Direct_modulation_map_file     = A.Direct_modulation_map_file;
    //if (A.Stimulus_type_arg 	== true)    t->S1_loc_type      = A.Stimulus_loc_type; // these two are now in the previous function set_coord_stim_and_map_from_defined_type
    //if (A.S2_Stimulus_type_arg 	== true)    t->S2_loc_type      = A.S2_Stimulus_loc_type;
    if (A.Dscale_base_map_file_arg 	== true)t->Dscale_base_map_file = A.Dscale_base_map_file;
    if (A.Dscale_mod_map_file_arg 	== true)t->Dscale_mod_map_file  = A.Dscale_mod_map_file;
    if (A.D_AR_scale_base_map_file_arg 	== true)t->D_AR_scale_base_map_file = A.D_AR_scale_base_map_file;
    if (A.D_AR_scale_mod_map_file_arg 	== true)t->D_AR_scale_mod_map_file  = A.D_AR_scale_mod_map_file;
    if (A.spatial_gradient_map_file_arg 		== true)    t->spatial_gradient_map_file     = A.spatial_gradient_map_file;
    if (A.Multiple_models_arg == true)	t->Multiple_models	= A.Multiple_models;
    if (A.Tissue_model_2_arg == true)	t->Tissue_model_2	= A.Tissue_model_2;

    // Idealised map settings
    if (A.map_shape_arg  == true)    t->ideal_map_shape    = A.map_shape;
    if (A.map_in_type_arg == true)   t->map_in_type        = A.map_in_type;
    if (A.map_x_loc_arg  == true)    t->ideal_map_x_loc    = A.map_x_loc;
    if (A.map_x_size_arg == true)    t->ideal_map_x_size   = A.map_x_size;
    if (A.map_y_loc_arg  == true)    t->ideal_map_y_loc    = A.map_y_loc;
    if (A.map_y_size_arg == true)    t->ideal_map_y_size   = A.map_y_size;
    if (A.map_z_loc_arg  == true)    t->ideal_map_z_loc    = A.map_z_loc;
    if (A.map_z_size_arg == true)    t->ideal_map_z_size   = A.map_z_size;

    //NETWORK
    if (A.junction_mod_map_on_arg       == true)  t->junction_mod_map_on    = A.junction_mod_map_on;
    if (A.junction_mod_map_file_arg     == true)  t->junction_mod_map_file  = A.junction_mod_map_file;
    if (A.junction_base_map_on_arg      == true)  t->junction_base_map_on   = A.junction_base_map_on;
    if (A.junction_base_map_file_arg    == true)  t->junction_base_map_file = A.junction_base_map_file;
    if (A.apply_symmetry_factor_arg     == true)  t->apply_symmetry_factor  = A.apply_symmetry_factor;

    // Now that Tissue type and orientation type have their final values, we can error check:
    if (strcmp(t->Orientation_type, "isotropic") != 0 && strcmp(t->Orientation_type, "anisotropic") != 0 && strcmp(t->Orientation_type, "three_eigenvectors") != 0) 
    {
        printf("ERROR! Orientation_type %s is invalid; must be \"isotropic\", \"anisotropic\" or \"three_eigenvectors\"\n", t->Orientation_type);
    }
    if (strcmp(t->Tissue_type, "homogeneous") != 0 && strcmp(t->Tissue_type, "heterogeneous") != 0) 
    {
        printf("ERROR! Tissue_type %s is invalid; must be \"homogeneous\" or \"heterogeneous\"\n", t->Tissue_type);
    }

}
// End overwrite properties from arguments ======================================================//|

// Create or read geometries ====================================================================\\|
// Select appropriate function
void select_tissue_geometry_function(Tissue_parameters t, SC_variables *sc, const char * PATH, const char* Output_dir)
{
    if (strcmp(t.Tissue_order, "geo") == 0) sc->N = read_geo_file(sc, sc->geo, t.geo_file, "Tissue_geometries", PATH, Output_dir, "anatomy", t.Ncelltypes);  // read geo || lib/Spatial_coupling.cc
    else if (strcmp(t.Tissue_order, "1D") == 0 || strcmp(t.Tissue_order, "2D") == 0 || strcmp(t.Tissue_order, "3D") == 0)
    {
        printf(">Creating idealised geometry....\n");
        if (strcmp(t.Tissue_type, "heterogeneous") == 0 || strcmp(t.D_uniformity, "regional") == 0 || strcmp(t.D_uniformity, "regional_map") == 0) 	create_idealised_geometry_heterogeneous(sc, t);
        else if (strcmp(t.Tissue_type, "homogeneous") == 0) create_idealised_geometry_homogeneous(sc);
        else
        {
            printf("ERROR: \"%s\" is not a valid Tissue type. Please select \"homogeneous\" or \"heterogeneous\"\n\n", t.Tissue_type);
            exit(1);
        }

        // Write vtk of geo file
        char *string = (char*)malloc(500);
        FILE *out;
        int idx;
        sprintf(string,"%s/Geometry_idealised.vtk", Output_dir);
        out = fopen(string, "wt");

        fprintf(out, "# vtk DataFile Version 3.0\n");
        fprintf(out, "vtk output\n");
        fprintf(out, "ASCII\n");
        fprintf(out, "DATASET STRUCTURED_POINTS\n");
        fprintf(out, "DIMENSIONS %d %d %d\n", sc->NX, sc->NY, sc->NZ);
        fprintf(out, "SPACING 1 1 1\n");
        fprintf(out, "ORIGIN 0 0 0\n");
        fprintf(out, "POINT_DATA %d\n", sc->NX*sc->NY*sc->NZ);
        fprintf(out, "SCALARS geo float 1\n");
        fprintf(out, "LOOKUP_TABLE default\n");

        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX*sc->NY*k);
                    if (sc->geo[idx] > 0)
                    {
                        fprintf(out, "%d ", sc->geo[idx]);
                    }
                    else fprintf(out, "-100 ");
                }
                fprintf(out, "\n");
            }
            fprintf(out, "\n");
        }
        fclose(out);
        free(string);
    }
    else
    {
        printf("ERROR: \"%s\" is not a valid Tissue order. Please select from \"1D\", \"2D\", \"3D\", \"geo\"\n\n", t.Tissue_order);
        exit(1);
    }
}

// Idealised create functions
void create_idealised_geometry_homogeneous(SC_variables *sc)
{
    int cell_count, idx;
    cell_count = 0;

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
                idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);
                sc->geo[idx] = 0;       // Default to empty space
                sc->geo[idx] = 1;		// In this case, all tissue is one celltype	
                cell_count++;
            }
        }
    }
    sc->N = cell_count;
    //printf("\tHomogeneous geometry size (X*Y*Z, %d * %d * %d) created || Ncells = %d\n", sc->NX, sc->NY, sc->NZ, sc->Ncell);	
}

// will create an idealised geometry segmented into Ncelltypes with junction boundaries as set
void create_idealised_geometry_heterogeneous(SC_variables *sc, Tissue_parameters t)
{
    int cell_count, idx;
    cell_count = 0;

    for (int k = 0; k < sc->NZ; k++)
    {
        for (int j = 0; j < sc->NY; j++)
        {
            for (int i = 0; i < sc->NX; i++)
            {
                idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);
                sc->geo[idx]	= 0;		// Default to empty space

                // Set first celltype which starts from x=0 bounds
                if (i < t.het_junction_X_location[1]) 	sc->geo[idx] = 1;
                else
                {
                    for (int c = 2; c < t.Ncelltypes; c++) // start from celltypes 2 up to second-to-last celltype
                    {
                        if (i >= t.het_junction_X_location[c-1] && i < t.het_junction_X_location[c]) sc->geo[idx] = c;
                    }
                    // And final celltype, x->NX
                    if (i >= t.het_junction_X_location[t.Ncelltypes-1]) sc->geo[idx] = t.Ncelltypes;
                }
                if (sc->geo[idx] > 0) cell_count++;	
            }
        }
    }
    sc->N = cell_count;
}

// geo read function
// End create or read geometries ================================================================//|

// Create or read stimulus ======================================================================\\|
void select_stimulus_area_function(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, int S2_CL)
{
    if (strcmp(t->Tissue_order, "geo") == 0) 
    {
        if (strcmp(t->S1_loc_type, "file") == 0)
        {
            t->Nstim 	= read_map_file(sc, t->stim_area, t->stim_file, "Tissue_geometries", PATH, Output_dir, "S1"); // lib/Spatial_coupling.cpp
        }
        else if (strcmp(t->S1_loc_type, "coords") == 0)
        {
            if (strcmp(t->S1_shape, "sphere") == 0) create_stimulus_area_sphere(sc,t->Tissue_order, t->stim_area, t->S1_x_loc, t->S1_x_size, t->S1_y_loc, t->S1_z_loc, &t->Nstim, Output_dir, "S1");
            else create_stimulus_area(sc,t->Tissue_order, t->stim_area, t->S1_x_loc, t->S1_x_size, t->S1_y_loc, t->S1_y_size, t->S1_z_loc, t->S1_z_size, &t->Nstim, Output_dir, "S1");
        }
        else
        {
            printf("ERROR: S1 location type %s is not valid for anatomical geoemtry models. Please select either \"file\" or \"coords\"\n", t->S1_loc_type);
            exit(1);
        }
        if (S2_CL != 0)
        {
            if (strcmp(t->S2_loc_type, "file") == 0)
            {
                t->Nstim_S2 = read_map_file(sc, t->S2_stim_area, t->S2_stim_file, "Tissue_geometries", PATH, Output_dir, "S2"); // lib/Spatial_coupling.cpp
            }
            else if (strcmp(t->S2_loc_type, "coords") == 0)
            {
                if (strcmp(t->S2_shape, "sphere") == 0) create_stimulus_area_sphere(sc,t->Tissue_order, t->S2_stim_area, t->S2_x_loc, t->S2_x_size, t->S2_y_loc, t->S2_z_loc, &t->Nstim_S2, Output_dir, "S2");
                else create_stimulus_area(sc,t->Tissue_order, t->S2_stim_area, t->S2_x_loc, t->S2_x_size, t->S2_y_loc, t->S2_y_size, t->S2_z_loc, t->S2_z_size, &t->Nstim_S2, Output_dir, "S2");
            }
            else
            {
                printf("ERROR: S2 location type %s is not valid for anatomical geoemtry models. Please select either \"file\" or \"coords\"\n", t->S2_loc_type);
                exit(1);
            }
        }
    }
    else 	//idealised; no need to repeat error checking performed when geomtry was read in
    {
        printf(">Creating idealised stimulus area....\n");
        if (strcmp(t->S1_shape, "sphere") == 0) create_stimulus_area_sphere(sc,t->Tissue_order, t->stim_area, t->S1_x_loc, t->S1_x_size, t->S1_y_loc, t->S1_z_loc, &t->Nstim, Output_dir, "S1");
        else create_stimulus_area(sc,t->Tissue_order, t->stim_area, t->S1_x_loc, t->S1_x_size, t->S1_y_loc, t->S1_y_size, t->S1_z_loc, t->S1_z_size, &t->Nstim, Output_dir, "S1");
        if (strcmp(t->S2_shape, "sphere") == 0) create_stimulus_area_sphere(sc,t->Tissue_order, t->S2_stim_area, t->S2_x_loc, t->S2_x_size, t->S2_y_loc, t->S2_z_loc, &t->Nstim_S2, Output_dir, "S2");
        else create_stimulus_area(sc,t->Tissue_order, t->S2_stim_area, t->S2_x_loc, t->S2_x_size, t->S2_y_loc, t->S2_y_size, t->S2_z_loc, t->S2_z_size, &t->Nstim_S2, Output_dir, "S2");
    }
}

// Multi timed and site stimulus functions
void select_stimulus_area_function_multi_stim(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, int S2_CL, int Nstims)
{
    // Cimplete this function - ideal and read in; read in = reads in number 1-N then creates map 0 or 1 for stim n
    ///*if (strcmp(t->S1_shape, "sphere") == 0) */create_stimulus_area_sphere(sc,t->Tissue_order, t->stim_area, t->S1_x_loc, t->S1_x_size, t->S1_y_loc, t->S1_z_loc, &t->Nstim);
    int idx;
    int cell_count = 0;

    if (strcmp(t->Tissue_order, "geo") == 0)
    {
        if (strcmp(t->S1_loc_type, "file") == 0)
        {
            // First, read in stim file to stim area as normally
            t->Nstim    = read_map_file(sc, t->stim_area, t->stim_file, "Tissue_geometries", PATH, Output_dir, "S1"); // lib/Spatial_coupling.cpp	

            // Now, need to loop over and assign stim_area[n] = 1 when stim_area = n
            // 1: assign stimmap[1->N] from stim map
            // 2: reassign stim map to just be for 1
            cell_count = 0;
            for (int k = 0; k < sc.NZ; k++)
            {
                for (int j = 0; j < sc.NY; j++)
                {
                    for (int i = 0; i < sc.NX; i++)
                    {
                        idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                        if (sc.geo[idx] > 0) // if it is an actual cell/node
                        {
                            // default all maps to zero
                            for (int n = 1; n < Nstims; n++) t->multi_stim_area[n][cell_count] = 0;	

                            // Assign correct map
                            if (t->stim_area[cell_count] > Nstims) { printf("ERROR: Multi-stim map has an entry greater than max number of different sites\n"); exit(1); }
                            if (t->stim_area[cell_count] > 1) t->multi_stim_area[t->stim_area[cell_count]-1][cell_count] = 1; // stim_area[c]-1 because we want multi_stim_area[1][n] to be if stim map = 2 (stim map = 1 is normal stim area)

                            cell_count++;
                        }	
                    }
                }
            }
            cell_count = 0;
            for (int k = 0; k < sc.NZ; k++)
            {
                for (int j = 0; j < sc.NY; j++)
                {
                    for (int i = 0; i < sc.NX; i++)
                    {
                        idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                        if (sc.geo[idx] > 0) // if it is an actual cell/node
                        {
                            if (t->stim_area[cell_count] != 1) t->stim_area[cell_count] = 0; // only keep this map IF equal to 1
                            cell_count++;
                        }
                    }
                }
            }
        }
        else
        {
            printf("ERROR: S1 location type %s is not valid for anatomical geoemtry models with multi-time stimuli. Please select \"file\" only\n", t->S1_loc_type);
            exit(1);
        }
    }
    else    //idealised; no need to repeat error checking performed when geomtry was read in
    {
        printf("ERROR: Multi-site stimulus is setup for 3D geo models only\n");
        exit(1);
        //create_stimulus_area(sc,t->Tissue_order, t->multi_stim_area[n], t->S1_x_loc+(n*20), t->S1_x_size, t->S1_y_loc, t->S1_y_size, t->S1_z_loc, t->S1_z_size, &t->Nstim);
    }
}

void create_stimulus_area(SC_variables sc, const char *Tissue_order, int * stim_area, int x, int xs, int y, int ys, int z, int zs, int *Nstim, const char* Output_dir, const char *ref)
{
    int cell_count	= 0;
    int stim_count	= 0;
    int idx;

    double r;

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                if (sc.geo[idx] > 0) // if it is an actual cell/node
                {
                    if (strcmp(Tissue_order, "1D") == 0) 
                    {
                        if (i >= x-xs && i <= x+xs)
                        {
                            stim_area[cell_count] = 1;		// Apply stimulus
                            stim_count ++;
                        }			
                        else stim_area[cell_count] = 0;		// Don't apply stimulus
                    }
                    else if (strcmp(Tissue_order, "2D") == 0)
                    {
                        if (i >= x-xs && i <= x+xs && j >= y-ys && j <= y+ys)
                        {
                            stim_area[cell_count] = 1;      // Apply stimulus
                            stim_count ++;
                        }
                        else stim_area[cell_count] = 0;     // Don't apply stimulus
                    }	
                    else if (strcmp(Tissue_order, "3D") == 0 || strcmp(Tissue_order, "geo") == 0)
                    {
                        if (i >= x-xs && i <= x+xs && j >= y-ys && j <= y+ys && k >= z-zs && k <= z+zs)
                        {
                            stim_area[cell_count] = 1;      // Apply stimulus
                            stim_count ++;
                        }
                        else stim_area[cell_count] = 0;     // Don't apply stimulus
                    }
                    cell_count ++;
                }
            }
        }
    }
    *Nstim 	= stim_count;

    char *string = (char*)malloc(500);
    FILE *out;

    sprintf(string,"%s/Map_%s.vtk", Output_dir, ref);
    out = fopen(string, "wt");

    printf("%s\n", string);

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    //fprintf(out, "SCALARS ImageFile float 1\n");
    fprintf(out, "SCALARS stim_%s float 1\n", ref);
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
                    fprintf(out, "%d ", stim_area[cell_count]);
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

void create_stimulus_area_sphere(SC_variables sc, const char *Tissue_order, int * stim_area, int x, int xs, int y, int z, int *Nstim, const char* Output_dir, const char *ref)
{
    int cell_count  = 0;
    int stim_count  = 0;
    int idx;

    double r;

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                if (sc.geo[idx] > 0) // if it is an actual cell/node
                {
                    if (strcmp(Tissue_order, "2D") == 0)
                    {
                        r = sqrt((i-x)*(i-x) + (j-y)*(j-y));
                        if (r < xs)
                        {
                            stim_area[cell_count] = 1;      // Apply stimulus
                            stim_count ++;
                        }
                        else stim_area[cell_count] = 0;
                    }
                    else if (strcmp(Tissue_order, "3D") == 0 || strcmp(Tissue_order, "geo") == 0)
                    {
                        r = sqrt((i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z));
                        if (r < xs)
                        {
                            stim_area[cell_count] = 1;      // Apply stimulus
                            stim_count ++;
                        }
                        else stim_area[cell_count] = 0;
                    }
                    cell_count ++;
                }
            }
        }
    }
    *Nstim  = stim_count;

    char *string = (char*)malloc(500);
    FILE *out;

    sprintf(string,"%s/Map_%s.vtk", Output_dir, ref);
    out = fopen(string, "wt");

    printf("%s\n", string);

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    //fprintf(out, "SCALARS ImageFile float 1\n");
    fprintf(out, "SCALARS stim_%s float 1\n", ref);
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
                    fprintf(out, "%d ", stim_area[cell_count]);
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
// End create or read stimulus ==================================================================//|

// Create or read orientation ===================================================================\\|
void set_orientation(SC_variables *sc, Tissue_parameters t, const char *PATH, const char *Tissue_order)
{
    // set orientation as a default to 0 (no fibre)
    for (int n = 0; n < sc->N; n++)
    {
        sc->ox[n]   = 0;
        sc->oy[n]   = 0;
        sc->oz[n]   = 0;
        
        sc->ox2[n]   = 0;
        sc->oy2[n]   = 0;
        sc->oz2[n]   = 0;
        sc->ox3[n]   = 0;
        sc->oy3[n]   = 0;
        sc->oz3[n]   = 0;
    }

    // Only write/read fibres if it will be anisotropic
    if (strcmp(t.Orientation_type, "anisotropic") == 0)
    {
        // now, overwrite with actual fibres if anisotropic
        if (strcmp(Tissue_order, "geo") == 0)	read_orientation_anatomical(sc, t, PATH);
        else                                   	create_orientation_ideal(sc, t);
    }
    else if(strcmp(t.Orientation_type, "three_eigenvectors") == 0)
    {
        // now, overwrite with actual fibres if anisotropic
        if (strcmp(Tissue_order, "geo") == 0)  
        {
            read_orientation_anatomical(sc, t, PATH);
            read_orientation_anatomical_transverse(sc, t, PATH);
        }
        else 
        {
            printf("three eigenvectors anisotropy type must be tissue type geo\n");
            exit(1);
        }
    }
}

void create_orientation_ideal(SC_variables *sc, Tissue_parameters t)
{
    for (int n = 0; n < sc->N; n++)
    {
        sc->ox[n]	= t.OX;		// local orientation = global orientation
        sc->oy[n]	= t.OY;		// local orientation = global orientation
        sc->oz[n]	= t.OZ;		// local orientation = global orientation
    }
}

void read_orientation_anatomical(SC_variables *sc, Tissue_parameters t, const char *PATH)
{
    FILE *in1, *in2, *in3;
    char *string = (char*)malloc(500);

    int idx;
    int count = 0;

    if (strcmp(t.orientation_file_type, "xyz") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_OX.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OY.dat", PATH, t.orientation_file_root);
        in2 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OZ.dat", PATH, t.orientation_file_root);
        in3 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL || in2 == NULL || in3 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f ", &X);
                    fscanf(in2, "%f ", &Y);
                    fscanf(in3, "%f ", &Z);

                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error  - fibre present where there is no geometry. Continuing..\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0) printf("Fibre error\n");
                        sc->ox[count]  = X;
                        sc->oy[count]  = Y;
                        sc->oz[count]  = Z;

                        //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox[count] > 1.0) sc->ox[count] = 1.0;
                        if (sc->oy[count] > 1.0) sc->oy[count] = 1.0;
                        if (sc->oz[count] > 1.0) sc->oz[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
            }
        }

        fclose(in1);
        fclose(in2);
        fclose(in3);
    } // end "xyz" style

    else if (strcmp(t.orientation_file_type, "orientation") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_orientation.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f %f %f ", &X, &Y, &Z);
                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0 && Y == 0 && Z == 0) printf("Fibre error\n");
                        sc->ox[count]  = X;
                        sc->oy[count]  = Y;
                        sc->oz[count]  = Z;

                        //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox[count] > 1.0) sc->ox[count] = 1.0;
                        if (sc->oy[count] > 1.0) sc->oy[count] = 1.0;
                        if (sc->oz[count] > 1.0) sc->oz[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
                //fscanf(in1, "\n");
            }
            //fscanf(in1, "\n");
        }
        fclose(in1);
    } // end "orientation" style

    else if (strcmp(t.orientation_file_type, "angles") == 0)
    {
        float theta, phi;

        sprintf(string, "%s/Tissue_geometries/%s_theta.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");

        sprintf(string, "%s/Tissue_geometries/%s_phi.dat", PATH, t.orientation_file_root);
        in2 = fopen(string, "r");

        if (in1 == NULL || in2 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read and convert components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f ", &theta);
                    fscanf(in2, "%f ", &phi);
                    // radians

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        sc->ox[count]  = sin(theta)*cos(phi);
                        sc->oy[count]  = cos(theta)*cos(phi);
                        sc->oz[count]  = sin(phi);

                        if (theta < -10 || phi < -10) printf("ERROR -10\n");

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
            }
        }
        fclose(in1);
        fclose(in2);
    }  // end "angles" style

    // angles style, where angles are defined from the short axis
    else if (strcmp(t.orientation_file_type, "angles_short_axis") == 0)
    {
        float theta, phi;

        sprintf(string, "%s/Tissue_geometries/%s_theta.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");

        sprintf(string, "%s/Tissue_geometries/%s_phi.dat", PATH, t.orientation_file_root);
        in2 = fopen(string, "r");

        if (in1 == NULL || in2 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read and convert components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f ", &theta);
                    fscanf(in2, "%f ", &phi);

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        sc->ox[count]  = cos(theta)*cos(phi);
                        sc->oy[count]  = sin(theta)*cos(phi);
                        sc->oz[count]  = sin(phi);

                        if (theta < -10 || phi < -10) printf("ERROR -10\n");

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
            }
        }
        fclose(in1);
        fclose(in2);
    }  // end "angles_short_axis" style
}

void read_orientation_anatomical_transverse(SC_variables *sc, Tissue_parameters t, const char *PATH)
{
    FILE *in1, *in2, *in3;
    char *string = (char*)malloc(500);

    int idx;
    int count = 0;

    // Transverse 1 (in sheet) ================================================================\\|

    if (strcmp(t.orientation_file_type, "xyz") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_OX2.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OY2.dat", PATH, t.orientation_file_root);
        in2 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OZ2.dat", PATH, t.orientation_file_root);
        in3 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL || in2 == NULL || in3 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f ", &X);
                    fscanf(in2, "%f ", &Y);
                    fscanf(in3, "%f ", &Z);

                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error  - fibre present where there is no geometry. Continuing..\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0) printf("Fibre error\n");
                        sc->ox2[count]  = X;
                        sc->oy2[count]  = Y;
                        sc->oz2[count]  = Z;

                        //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox2[count] > 1.0) sc->ox2[count] = 1.0;
                        if (sc->oy2[count] > 1.0) sc->oy2[count] = 1.0;
                        if (sc->oz2[count] > 1.0) sc->oz2[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
            }
        }

        fclose(in1);
        fclose(in2);
        fclose(in3);
    } // end "xyz" style

    else if (strcmp(t.orientation_file_type, "orientation") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_orientation2.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f %f %f ", &X, &Y, &Z);
                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0 && Y == 0 && Z == 0) printf("Fibre error\n");
                        sc->ox2[count]  = X;
                        sc->oy2[count]  = Y;
                        sc->oz2[count]  = Z;

                        //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox2[count] > 1.0) sc->ox2[count] = 1.0;
                        if (sc->oy2[count] > 1.0) sc->oy2[count] = 1.0;
                        if (sc->oz2[count] > 1.0) sc->oz2[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                         }
                }
                //fscanf(in1, "\n");
            }
            //fscanf(in1, "\n");
        }
        fclose(in1);
    } // end "orientation" style
    // End Transverse 1 (in sheet) =============================================================//|

    count = 0;
    // Transverse 2 ============================================================================\\|
        if (strcmp(t.orientation_file_type, "xyz") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_OX3.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OY3.dat", PATH, t.orientation_file_root);
        in2 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        sprintf(string, "%s/Tissue_geometries/%s_OZ3.dat", PATH, t.orientation_file_root);
        in3 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL || in2 == NULL || in3 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f ", &X);
                    fscanf(in2, "%f ", &Y);
                    fscanf(in3, "%f ", &Z);

                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error  - fibre present where there is no geometry. Continuing..\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0) printf("Fibre error\n");
                        sc->ox3[count]  = X;
                        sc->oy3[count]  = Y;
                        sc->oz3[count]  = Z;

                           //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox3[count] > 1.0) sc->ox3[count] = 1.0;
                        if (sc->oy3[count] > 1.0) sc->oy3[count] = 1.0;
                        if (sc->oz3[count] > 1.0) sc->oz3[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;
                    }
                }
            }
        }

        fclose(in1);
        fclose(in2);
        fclose(in3);
    } // end "xyz" style

    else if (strcmp(t.orientation_file_type, "orientation") == 0)
    {
        printf("Reading fibres, xyz style, from file...\n");
        float X, Y, Z;

        // open fibre files
        sprintf(string, "%s/Tissue_geometries/%s_orientation3.dat", PATH, t.orientation_file_root);
        in1 = fopen(string, "r");
        printf("Fibre file %s open\n", string);

        if (in1 == NULL)
        {
            printf("Cannot load fibre files %s\t :: is the path correct? Does the file exist in that path?\n", string);
            exit(1);
        }

        // read in components to sc array
        for (int k = 0; k < sc->NZ; k++)
        {
            for (int j = 0; j < sc->NY; j++)
            {
                for (int i = 0; i < sc->NX; i++)
                {
                    idx = i + (sc->NX*j) + (sc->NX * sc->NY * k);

                    fscanf(in1, "%f %f %f ", &X, &Y, &Z);
                    //if (X != 0 && sc->geo[idx] == 0)  printf("Fibre error\n");

                    if (sc->geo[idx] > 0) // note that ox, oy, oz arrays are only size Ncell, not NX*NY*NZ
                    {
                        //if (X == 0 && Y == 0 && Z == 0) printf("Fibre error\n");
                        sc->ox3[count]  = X;
                        sc->oy3[count]  = Y;
                        sc->oz3[count]  = Z;

                        //ensure "normalised" (this is a soft check of each component; not a full check of the vector)
                        if (sc->ox3[count] > 1.0) sc->ox3[count] = 1.0;
                        if (sc->oy3[count] > 1.0) sc->oy3[count] = 1.0;
                        if (sc->oz3[count] > 1.0) sc->oz3[count] = 1.0;

                        //printf("%f %f %f\n", sc->ox[count], sc->oy[count], sc->oz[count]);

                        count ++;

                                         }
                }
                //fscanf(in1, "\n");
            }
            //fscanf(in1, "\n");
        }
        fclose(in1);
    } // end "orientation" style
    // End Transverse 2 ========================================================================//|
}// end read orientation anatomical transverse

void output_fibre_orientation(SC_variables sc, Tissue_parameters t, const char* Output_dir)
{
    char *string = (char*)malloc(500);
    sprintf(string, "%s/Orientation.vtk", Output_dir);

    FILE *out;
    out = fopen(string, "wt");

    int count = 0;
    int idx;

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    fprintf(out, "SCALARS ImageFile float 3\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
                if (sc.geo[idx] > 0)
                { 
                    fprintf(out, "%f %f %f ", sc.ox[count], sc.oy[count], sc.oz[count]);
                    count++;
                }
                else fprintf(out, "0 0 0 ");
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }
    fclose(out);
}

void output_fibre_orientation_o2(SC_variables sc, Tissue_parameters t, const char* Output_dir)
{
    char *string = (char*)malloc(500);
    sprintf(string, "%s/Orientation2.vtk", Output_dir);

    FILE *out;
    out = fopen(string, "wt");

    int count = 0;
    int idx;

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    fprintf(out, "SCALARS ImageFile float 3\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
                if (sc.geo[idx] > 0)
                {
                    fprintf(out, "%f %f %f ", sc.ox2[count], sc.oy2[count], sc.oz2[count]);
                    count++;
                }
                else fprintf(out, "0 0 0 ");
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }
    fclose(out);
}

void output_fibre_orientation_o3(SC_variables sc, Tissue_parameters t, const char* Output_dir)
{
    char *string = (char*)malloc(500);
    sprintf(string, "%s/Orientation3.vtk", Output_dir);

    FILE *out;
    out = fopen(string, "wt");

    int count = 0;
    int idx;

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    fprintf(out, "SCALARS ImageFile float 3\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX*sc.NY*k);
                if (sc.geo[idx] > 0)
                {
                    fprintf(out, "%f %f %f ", sc.ox3[count], sc.oy3[count], sc.oz3[count]);
                    count++;
                }
                else fprintf(out, "0 0 0 ");
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }
    fclose(out);
}
// End Create or read orientation ===============================================================//|

// Read other maps ==============================================================================\\|
void create_or_read_map_double(Tissue_parameters *t, SC_variables sc, const char *PATH, const char* Output_dir, double *map, const char *map_file, const char *ref)
{
    int Nmap;
    if (strcmp(t->Tissue_order, "geo") == 0)	
    {
        if (strcmp(t->map_in_type, "file") == 0) Nmap    =  read_map_file_double(sc, map, map_file, "Tissue_geometries", PATH, Output_dir, ref); // lib/Spatial_coupling.cpp
        else  if (strcmp(t->map_in_type, "coords") == 0)    
        {
            create_map_patch(sc, t->Tissue_order, map, t->ideal_map_x_loc, t->ideal_map_x_size, t->ideal_map_y_loc, t->ideal_map_y_size, t->ideal_map_z_loc, t->ideal_map_z_size, Output_dir, ref, t->ideal_map_shape);

            if (strcmp(ref, "ISO") == 0);
            else if (strcmp(ref, "remodelling") == 0);
            else if (strcmp(ref, "ACh") == 0);
            else if (strcmp(ref, "SRF") == 0);
            else if (strcmp(ref, "Dmod") == 0);
            else if (strcmp(ref, "D_AR_mod") == 0);
            else if (strcmp(ref, "Direct_modulation") == 0);
            else
            {
                printf("ERROR: The functionality for map type (%s) is currently NOT implemented for coordinates and can only be used with geo and a map file\n", ref);
                exit(1);
            }
        }
        else
        {
            printf("ERROR: map in type %s is not valid for anatomical geoemtry models. Please select either \"file\" or \"coords\"\n", t->S1_loc_type);
            exit(1);
        }
    }
    else // idealised geo models
    {
        if (strcmp(t->map_in_type, "coords") == 0)
            create_map_patch(sc, t->Tissue_order, map, t->ideal_map_x_loc, t->ideal_map_x_size, t->ideal_map_y_loc, t->ideal_map_y_size, t->ideal_map_z_loc, t->ideal_map_z_size, Output_dir, ref, t->ideal_map_shape);
        else
        {
            printf("ERROR: map in type %s is not valid for idealised geoemtry models. Please select \"coords\" only\n", t->S1_loc_type);
            exit(1);
        }

        if (strcmp(ref, "ISO") == 0);
        else if (strcmp(ref, "remodelling") == 0);
        else if (strcmp(ref, "ACh") == 0);
        else if (strcmp(ref, "SRF") == 0);
        else if (strcmp(ref, "Dmod") == 0);
        else if (strcmp(ref, "D_AR_mod") == 0);
        else if (strcmp(ref, "Direct_modulation") == 0);
        else 
        {
            printf("ERROR: The functionality for map type (%s) is currently NOT in idealised models and can only be used with geo and a map file\n", ref);
            exit(1);
        }
    }
}

void create_map_patch(SC_variables sc, const char *Tissue_order, double * map_patch, int x, int xs, int y, int ys, int z, int zs, const char* Output_dir, const char *ref, const char *shape)
{
    int cell_count  = 0;
    int idx;

    double r;

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                if (sc.geo[idx] > 0) // if it is an actual cell/node
                {
                    if (strcmp(Tissue_order, "1D") == 0)
                    {
                        if (i >= x-xs && i <= x+xs)
                        {
                            map_patch[cell_count] = 1;      
                        }
                        else map_patch[cell_count] = 0;     
                    }
                    else if (strcmp(Tissue_order, "2D") == 0)
                    {
                        r = sqrt((i-x)*(i-x) + (j-y)*(j-y));

                        if (strcmp(shape, "sphere") == 0)
                        {
                            if (r < xs) map_patch[cell_count] = 1;
                            else map_patch[cell_count] = 0;
                        }
                        else 
                        {
                            if (i >= x-xs && i <= x+xs && j >= y-ys && j <= y+ys)
                            {
                                map_patch[cell_count] = 1;     
                            }
                            else map_patch[cell_count] = 0;     
                        }
                    }
                    else if (strcmp(Tissue_order, "3D") == 0 || strcmp(Tissue_order, "geo") == 0)
                    {
                        r = sqrt((i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z));

                        if (strcmp(shape, "sphere") == 0)
                        {
                            if (r < xs) map_patch[cell_count] = 1;
                            else map_patch[cell_count] = 0;
                        }
                        else
                        {
                            if (i >= x-xs && i <= x+xs && j >= y-ys && j <= y+ys && k >= z-zs && k <= z+zs)
                            {
                                map_patch[cell_count] = 1;      
                            }
                            else map_patch[cell_count] = 0;    
                        }
                    }
                    cell_count ++;
                }
            }
        }
    }

    // output vtk of map
    FILE *out;
    char *string = (char*)malloc(500);
    sprintf(string,"%s/Map_%s.vtk", Output_dir, ref);
    out = fopen(string, "wt");

    sprintf(string, "%s", ref);

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
                    fprintf(out, "%f ", map_patch[cell_count]);
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
// End Create or read other maps ================================================================//|

// Diffusion tensor non-uniformity scale factor setup ===========================================\\|
// Baseline non-uniformity
void update_D_arrays_Dscale_baseline(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory)
{
    if (strcmp(t->D_uniformity, "map") == 0 || strcmp(t->D_uniformity, "regional_map") == 0)
    {
        create_or_read_map_double(t, *sc, PATH, directory, t->Dscale_base_map, t->Dscale_base_map_file, "Duniformity"); // lib/Tissue.cpp
        create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_base_map, t->D_AR_scale_base_map_file, "D_AR_uniformity"); // lib/Tissue.cpp
    }

    for (int n = 0; n < sc->N; n++)
    {
        if (strcmp(t->D_uniformity, "regional") == 0)    // Scale D1 and D2 by regional settings
        {
            //sc->geo_linear[n] is celltype number of cell n; non_uniform_D1_scale[celltype] = scale factor for D for celltype
            sc->D1[n]       *= t->non_uniform_D1_scale[sc->geo_linear[n]];
            sc->D2[n]       = sc->D1[n] / (t->D_AR * t->non_uniform_AR_scale[sc->geo_linear[n]]); // global AR x local AR scale
        }
        else if (strcmp(t->D_uniformity, "map") == 0) // scale by map
        {
            sc->D1[n]       *= t->Dscale_base_map[n];
            sc->D2[n]       = sc->D1[n] / (t->D_AR * t->D_AR_scale_base_map[n]); 
        }
        else if (strcmp(t->D_uniformity, "regional_map") == 0)    // Scale by region then by map
        {
            //sc->geo_linear[n] is celltype number of cell n; non_uniform_D1_scale[celltype] = scale factor for D for celltype
            sc->D1[n]       *= t->non_uniform_D1_scale[sc->geo_linear[n]] * t->Dscale_base_map[n];
            sc->D2[n]       = sc->D1[n] / (t->D_AR * t->non_uniform_AR_scale[sc->geo_linear[n]] * t->D_AR_scale_base_map[n]); // global AR x local AR scale
        }
        else
        {
            printf("ERROR: D_uniformity \"%s\" is invalid: must be \"uniform\", \"regional\", \"map\" or \"regional_map\"\n", t->D_uniformity);
            exit(1);
        }
    }
}

// Scale factor due to modification
void update_D_arrays_Dscale_mod(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory)
{
    double Dscale;
    if (strcmp(t->Dscale_map_on, "Off") == 0) 
    {
        if (strcmp(t->D_AR_scale_map_on, "On") == 0) create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_mod_map, t->D_AR_scale_mod_map_file, "D_AR_mod"); // lib/Tissue.cpp
        for (int n = 0; n < sc->N; n++)
        { 
            Dscale = t->Dscale;
            sc->D1[n] *= Dscale;
            if (strcmp(t->D_AR_scale_map_on, "Off") == 0) sc->D2[n] *= Dscale/t->D_AR_scale; // D2 must first be rescaled by D1scale, then D_AR_scale homogeneously applied
            else if (strcmp(t->D_AR_scale_map_on, "On") == 0)
            {
                sc->D2[n] *= Dscale/( (1 - t->D_AR_scale_mod_map[n]) + (t->D_AR_scale_mod_map[n] * t->D_AR_scale) ); // map = 0 (control) -> scale = 1; map = 1 -> scale = Dscale; 0<map<1, linear combination
            }
        }
    }
    else if (strcmp(t->Dscale_map_on, "On") == 0)
    {
        create_or_read_map_double(t, *sc, PATH, directory, t->Dscale_mod_map, t->Dscale_mod_map_file, "Dmod"); // lib/Tissue.cpp
        if (strcmp(t->D_AR_scale_map_on, "On") == 0) create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_mod_map, t->D_AR_scale_mod_map_file, "D_AR_mod"); // lib/Tissue.cpp
        for (int n = 0; n < sc->N; n++)
        {
            Dscale = ( (1 - t->Dscale_mod_map[n]) + (t->Dscale_mod_map[n] * t->Dscale) ); // map = 0 (control) -> scale = 1; map = 1 -> scale = Dscale; 0<map<1, linear combination
            sc->D1[n] *= Dscale;
            if (strcmp(t->D_AR_scale_map_on, "Off") == 0) sc->D2[n] *= Dscale/t->D_AR_scale; // D2 must first be rescaled by D1scale, then D_AR_scale homogeneously applied
            else if (strcmp(t->D_AR_scale_map_on, "On") == 0)
            {
                sc->D2[n] *= Dscale/( (1 - t->D_AR_scale_mod_map[n]) + (t->D_AR_scale_mod_map[n] * t->D_AR_scale) ); 
            }

            //sc->D1[n] *= ( (1 - t->Dscale_mod_map[n]) + (t->Dscale_mod_map[n] * t->Dscale) ); // map = 0 (control) -> D1; map = 1 -> D1.Dscale; 0<map<1, linear combination
            //sc->D2[n] = sc->D1[n]/(t->D_AR * t->D_AR_scale); // preserve AR after modification D1
        }
    }
}
// End Diffusion tensor setup ===================================================================//|

// Modify neighbours for region disconnect ======================================================\\|
void Modify_neighbours_region_disconnect(SC_variables *sc, Tissue_parameters *t)
{
    int count = 0;  // counter of number of cells
    int idx;        // 3D identifier for each cell
    int idx_xm;     // 3D identifier for x-1 cell
    int idx_xp;     // 3D identifier for x+1 cell
    int idx_ym;     // 3D identifier for y-1 cell
    int idx_yp;     // 3D identifier for y+1 cell
    int idx_zm;     // 3D identifier for z-1 cell
    int idx_zp;     // 3D identifier for z+1 cell
    int idx_xm_ym;  // 3D identifier for x-1, y-1 cell
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
                idx         =  i    + (sc->NX *  j)     + (sc->NX * sc->NY *  k   );
                idx_xm      = (i-1) + (sc->NX *  j)     + (sc->NX * sc->NY *  k   );
                idx_xp      = (i+1) + (sc->NX *  j)     + (sc->NX * sc->NY *  k   );
                idx_ym      =  i    + (sc->NX * (j-1))  + (sc->NX * sc->NY *  k   );
                idx_yp      =  i    + (sc->NX * (j+1))  + (sc->NX * sc->NY *  k   );
                idx_zm      =  i    + (sc->NX *  j)     + (sc->NX * sc->NY * (k-1));
                idx_zp      =  i    + (sc->NX *  j)     + (sc->NX * sc->NY * (k+1));

                idx_xm_ym   = (i-1) + (sc->NX * (j-1))  + (sc->NX * sc->NY *  k   );
                idx_xm_yp   = (i-1) + (sc->NX * (j+1))  + (sc->NX * sc->NY *  k   );
                idx_xm_zm   = (i-1) + (sc->NX *  j)     + (sc->NX * sc->NY * (k-1));
                idx_xm_zp   = (i-1) + (sc->NX *  j)     + (sc->NX * sc->NY * (k+1));

                idx_xp_ym   = (i+1) + (sc->NX * (j-1))  + (sc->NX * sc->NY *  k   );
                idx_xp_yp   = (i+1) + (sc->NX * (j+1))  + (sc->NX * sc->NY *  k   );
                idx_xp_zm   = (i+1) + (sc->NX *  j)     + (sc->NX * sc->NY * (k-1));
                idx_xp_zp   = (i+1) + (sc->NX *  j)     + (sc->NX * sc->NY * (k+1));

                idx_ym_zm   =  i    + (sc->NX * (j-1))  + (sc->NX * sc->NY * (k-1));
                idx_ym_zp   =  i    + (sc->NX * (j-1))  + (sc->NX * sc->NY * (k+1));
                idx_yp_zm   =  i    + (sc->NX * (j+1))  + (sc->NX * sc->NY * (k-1));
                idx_yp_zp   =  i    + (sc->NX * (j+1))  + (sc->NX * sc->NY * (k+1));

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
                    for (int DN = 0; DN < t->Ndisconnected_regions; DN++) // loop over number of disconnected region pairs
                    {
                        // x-1  -> if current node is region 0 for disconnect pair DN and xminus is region 1, or current node is region 1 and xminus is region 0, then disconnect (return self to neighbour map)
                        if (i > 0) 
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm] == t->disconnect_regions[DN][0]) ) sc->xm[count] = count;

                        // x+1  -> if current node is region 0 for disconnect pair DN and xplus is region 1, or current node is region 1 and xplus is region 0, then disconnect (return self to neighbour map)
                        if (i < sc->NX-1) 
                            if( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp] == t->disconnect_regions[DN][0]) ) sc->xp[count] = count;

                        // y-1
                        if (j > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_ym] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_ym] == t->disconnect_regions[DN][0]) ) sc->ym[count] = count;

                        // y+1
                        if (j < sc->NY-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp] == t->disconnect_regions[DN][0]) ) sc->yp[count] = count;

                        // z-1
                        if (k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_zm] == t->disconnect_regions[DN][0]) ) sc->zm[count] = count;

                        // z+1
                        if (k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_zp] == t->disconnect_regions[DN][0]) ) sc->zp[count] = count;

                        // x-1 y-1
                        if (i > 0 && j > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_ym] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_ym] == t->disconnect_regions[DN][0]) ) sc->xm_ym[count] = count;

                        // x-1 y+1
                        if (i > 0 && j < sc->NY-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp] == t->disconnect_regions[DN][0]) ) sc->xm_yp[count] = count;

                        // x-1 z-1
                        if (i > 0 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_zm] == t->disconnect_regions[DN][0]) ) sc->xm_zm[count] = count;

                        // x-1 z+1
                        if (i > 0 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_zp] == t->disconnect_regions[DN][0]) ) sc->xm_zp[count] = count;

                        // x+1 y-1
                        if (i < sc->NX-1 && j > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_ym] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_ym] == t->disconnect_regions[DN][0]) ) sc->xp_ym[count] = count;

                        // x+1 y+1
                        if (i < sc->NX-1 && j < sc->NY-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp] == t->disconnect_regions[DN][0]) ) sc->xp_yp[count] = count;

                        // x+1 z-1
                        if (i < sc->NX-1 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_zm] == t->disconnect_regions[DN][0]) ) sc->xp_zm[count] = count;

                        // x+1 z+1
                        if (i < sc->NX-1 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_zp] == t->disconnect_regions[DN][0]) ) sc->xp_zp[count] = count;

                        // y-1 z-1
                        if (j > 0 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_ym_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_ym_zm] == t->disconnect_regions[DN][0]) ) sc->ym_zm[count] = count;

                        // y-1 z+1
                        if (j > 0 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_ym_zp] == t->disconnect_regions[DN][0]) ) sc->ym_zp[count] = count;

                        // y+1 z-1
                        if (j < sc->NY-1 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp_zm] == t->disconnect_regions[DN][0]) ) sc->yp_zm[count] = count;

                        // y+1 z+1
                        if (j < sc->NY-1 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_yp_zp] == t->disconnect_regions[DN][0]) ) sc->yp_zp[count] = count;

                        // x, y,z
                        if (i > 0 && j > 0 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_ym_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_ym_zm] == t->disconnect_regions[DN][0]) ) sc->xm_ym_zm[count] = count;
                        if (i > 0 && j > 0 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_ym_zp] == t->disconnect_regions[DN][0]) ) sc->xm_ym_zp[count] = count;
                        if (i > 0 && j < sc->NY-1 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp_zm] == t->disconnect_regions[DN][0]) ) sc->xm_yp_zm[count] = count;
                        if (i > 0 && j < sc->NY-1 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xm_yp_zp] == t->disconnect_regions[DN][0]) ) sc->xm_yp_zp[count] = count;
                        if (i < sc->NX-1 && j > 0 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_ym_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_ym_zm] == t->disconnect_regions[DN][0]) ) sc->xp_ym_zm[count] = count;
                        if (i < sc->NX-1 && j > 0 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_ym_zp] == t->disconnect_regions[DN][0]) ) sc->xp_ym_zp[count] = count;
                        if (i < sc->NX-1 && j < sc->NY-1 && k > 0)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp_zm] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp_zm] == t->disconnect_regions[DN][0]) ) sc->xp_yp_zm[count] = count;
                        if (i < sc->NX-1 && j < sc->NY-1 && k < sc->NZ-1)
                            if ( (sc->geo[idx] == t->disconnect_regions[DN][0] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][1]) || (sc->geo[idx] == t->disconnect_regions[DN][1] && sc->geo[idx_xp_yp_zp] == t->disconnect_regions[DN][0]) ) sc->xp_yp_zp[count] = count;

                    } // end DN for
                    count++; // ensure count is updated for every node, independent of connection type
                } // end geo if
            } // end NX
        } // end NY
    } // end NZ
}
// End Modify neighbours for region disconnect ==================================================//|


// NETWORK version of D setup ===================================================================\\|
// These apply the same as above using "D base/mod" -> specific ones for the network map aplied elsewhere
// Baseline non-uniformity
void update_G_arrays_Dscale_baseline(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory)
{
    if (strcmp(t->D_uniformity, "map") == 0 || strcmp(t->D_uniformity, "regional_map") == 0)
    {
        create_or_read_map_double(t, *sc, PATH, directory, t->Dscale_base_map, t->Dscale_base_map_file, "Duniformity"); // lib/Tissue.cpp
        create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_base_map, t->D_AR_scale_base_map_file, "D_AR_uniformity"); // lib/Tissue.cpp
    }

    for (int n = 0; n < sc->N; n++)
    {
        if (strcmp(t->D_uniformity, "regional") == 0)    // Scale D1 and D2 by regional settings
        {
            //sc->geo_linear[n] is celltype number of cell n; non_uniform_D1_scale[celltype] = scale factor for D for celltype
            sc->Gl[n]       *= t->non_uniform_D1_scale[sc->geo_linear[n]];
            sc->Gt[n]       = sc->Gt[n] / (t->D_AR * t->non_uniform_AR_scale[sc->geo_linear[n]]); // global AR x local AR scale
        }
        else if (strcmp(t->D_uniformity, "map") == 0) // scale by map
        {
            sc->Gl[n]       *= t->Dscale_base_map[n];
            sc->Gt[n]       = sc->Gl[n] / (t->D_AR * t->D_AR_scale_base_map[n]);
        }
        else if (strcmp(t->D_uniformity, "regional_map") == 0)    // Scale by region then by map
        {
            //sc->geo_linear[n] is celltype number of cell n; non_uniform_D1_scale[celltype] = scale factor for D for celltype
            sc->Gl[n]       *= t->non_uniform_D1_scale[sc->geo_linear[n]] * t->Dscale_base_map[n];
            sc->Gt[n]       = sc->Gl[n] / (t->D_AR * t->non_uniform_AR_scale[sc->geo_linear[n]] * t->D_AR_scale_base_map[n]); // global AR x local AR scale
        }
        else
        {
            printf("ERROR: D_uniformity \"%s\" is invalid: must be \"uniform\", \"regional\", \"map\" or \"regional_map\"\n", t->D_uniformity);
            exit(1);
        }
    }
}

// Scale factor due to modification
void update_G_arrays_Dscale_mod(SC_variables *sc, Tissue_parameters *t, const char *PATH, const char* directory)
{
    double Dscale;
    if (strcmp(t->Dscale_map_on, "Off") == 0)
    {
        if (strcmp(t->D_AR_scale_map_on, "On") == 0) create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_mod_map, t->D_AR_scale_mod_map_file, "D_AR_mod"); // lib/Tissue.cpp
        for (int n = 0; n < sc->N; n++)
        {
            Dscale = t->Dscale;
            sc->Gl[n] *= Dscale;
            if (strcmp(t->D_AR_scale_map_on, "Off") == 0) sc->Gt[n] *= Dscale/t->D_AR_scale; // D2 must first be rescaled by D1scale, then D_AR_scale homogeneously applied
            else if (strcmp(t->D_AR_scale_map_on, "On") == 0)
            {
                sc->Gt[n] *= Dscale/( (1 - t->D_AR_scale_mod_map[n]) + (t->D_AR_scale_mod_map[n] * t->D_AR_scale) ); // map = 0 (control) -> scale = 1; map = 1 -> scale = Dscale; 0<map<1, linear combination
            }
        }
    }
    else if (strcmp(t->Dscale_map_on, "On") == 0)
    {
        create_or_read_map_double(t, *sc, PATH, directory, t->Dscale_mod_map, t->Dscale_mod_map_file, "Dmod"); // lib/Tissue.cpp
        if (strcmp(t->D_AR_scale_map_on, "On") == 0) create_or_read_map_double(t, *sc, PATH, directory, t->D_AR_scale_mod_map, t->D_AR_scale_mod_map_file, "D_AR_mod"); // lib/Tissue.cpp
        for (int n = 0; n < sc->N; n++)
        {
            Dscale = ( (1 - t->Dscale_mod_map[n]) + (t->Dscale_mod_map[n] * t->Dscale) ); // map = 0 (control) -> scale = 1; map = 1 -> scale = Dscale; 0<map<1, linear combination
            sc->Gl[n] *= Dscale;
            if (strcmp(t->D_AR_scale_map_on, "Off") == 0) sc->Gt[n] *= Dscale/t->D_AR_scale; // D2 must first be rescaled by D1scale, then D_AR_scale homogeneously applied
            else if (strcmp(t->D_AR_scale_map_on, "On") == 0)
            {
                sc->Gt[n] *= Dscale/( (1 - t->D_AR_scale_mod_map[n]) + (t->D_AR_scale_mod_map[n] * t->D_AR_scale) );
            }

            //sc->D1[n] *= ( (1 - t->Dscale_mod_map[n]) + (t->Dscale_mod_map[n] * t->Dscale) ); // map = 0 (control) -> D1; map = 1 -> D1.Dscale; 0<map<1, linear combination
            //sc->D2[n] = sc->D1[n]/(t->D_AR * t->D_AR_scale); // preserve AR after modification D1
        }
    }
}
// End NETWORK version of D setup ===============================================================//|

// Create/read phase map for phase re-entry =====================================================\\|
void create_read_phasemap(int * phase, Tissue_parameters t, SC_variables sc, const char* PATH, const char* Output_dir)
{
    if (strcmp(t.Tissue_order, "1D") == 0)
    {
        printf("ERROR!: 1D is not an appropriate tissue type for phase re-entry map\n");
        exit(1);
    }
    else if (strcmp(t.Tissue_order, "2D") == 0)
    {
        create_phasemap_2D(phase, sc, Output_dir);
    }
    else if (strcmp(t.Tissue_order, "3D") == 0)
    {
        create_phasemap_3D(phase, sc, Output_dir);
    }
    else if (strcmp(t.Tissue_order, "geo") == 0)
    {
        read_map_file(sc, phase, t.phase_file, "3D_phasemaps", PATH, Output_dir, "phase"); // read -> lib/Spatial_coupling.cpp
    }
}

void create_phasemap_2D(int * phase, SC_variables sc, const char* Output_dir)
{
    int Cx, Cy;
    Cx = sc.NX/2;
    Cy = sc.NY/2;

    printf("Centre of phase map is %d %d\n", Cx, Cy);

    int cell_count = 0;

    double x, y;
    double theta;

    for (int j = 0; j < sc.NY; j++)
    {
        for (int i = 0; i < sc.NX; i++)
        {
            int idx = i + (sc.NX*j);
            if (sc.geo[idx] > 0)
            {
                x = i - Cx;
                y = j - Cy;

                if (x == 0 && y >= 0) theta = 1.5*M_PI;//-M_PI/2;
                else if (x == 0) theta = M_PI/2;
                else theta = atan(y/x);
                if (x <= 0) theta -= M_PI;
                theta += 1.5*M_PI;
                theta *= (100/3.15);///M_PI); // scales 0-> 2PI to 0->200
                phase[cell_count] = (int)theta;//*(90/M_PI);   
                if (phase[cell_count] >= 200) phase[cell_count] = 0;
                if (phase[cell_count] < 0) phase[cell_count] = 0;
                //printf("phase %d = %f\n", cell_count, phase[cell_count]); 
                cell_count++;
            }
        }
    }

    // And output
    char *string = (char*)malloc(500);
    FILE *out;

    sprintf(string,"%s/Phase_Map_2D.vtk", Output_dir);
    out = fopen(string, "wt");

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, 1);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*1);
    fprintf(out, "SCALARS phase_map float 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    cell_count = 0;
    for (int k = 0; k < 1; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                int idx = i + (sc.NX*j);
                if (sc.geo[idx] > 0)
                {
                    fprintf(out, "%d ", phase[cell_count]);
                    cell_count++;
                }
                else fprintf(out, "-100 ");
            }
        }
    }
    fclose(out);
}

void create_phasemap_3D(int * phase, SC_variables sc, const char* Output_dir)
{
    int Cx, Cy;
    Cx = sc.NX/2;
    Cy = sc.NY/2;

    printf("Centre of phase map is %d %d\n", Cx, Cy);

    int cell_count = 0;

    double x, y;
    double theta;

    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                int idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                if (sc.geo[idx] > 0)
                {
                    x = i - Cx;
                    y = j - Cy;

                    if (x == 0 && y >= 0) theta = 1.5*M_PI;//-M_PI/2;
                    else if (x == 0) theta = M_PI/2;
                    else theta = atan(y/x);
                    if (x <= 0) theta -= M_PI;
                    theta += 1.5*M_PI;
                    theta *= (100/3.15);///M_PI); // scales 0-> 2PI to 0->200
                    phase[cell_count] = (int)theta;//*(90/M_PI);   
                    if (phase[cell_count] >= 200) phase[cell_count] = 0;
                    if (phase[cell_count] < 0) phase[cell_count] = 0;
                    //printf("phase %d = %f\n", cell_count, phase[cell_count]); 
                    cell_count++;
                }
            }
        }
    }

    // And output
    char *string = (char*)malloc(500);
    FILE *out;

    sprintf(string,"%s/Phase_Map_3D.vtk", Output_dir);
    out = fopen(string, "wt");

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", sc.NX, sc.NY, sc.NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", sc.NX*sc.NY*sc.NZ);
    fprintf(out, "SCALARS phase_map float 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    cell_count = 0;
    for (int k = 0; k < sc.NZ; k++)
    {
        for (int j = 0; j < sc.NY; j++)
        {
            for (int i = 0; i < sc.NX; i++)
            {
                int idx = i + (sc.NX*j) + (sc.NX * sc.NY * k);
                if (sc.geo[idx] > 0)
                {
                    fprintf(out, "%d ", phase[cell_count]);
                    cell_count++;
                }
                else fprintf(out, "-100 ");
            }
        }
    }
    fclose(out);
}
// End Create phase map for phase re-entry ======================================================//|

// conduction velocity calculation ==============================================================\\|
void set_CV_cells(Cell_parameters p, Tissue_parameters *t)
{
    // Needs to be here such that correct dx is used

    int dist_S_1 = 10;  // distance in x or y from centre (where stimulus is) to first cell (independent of NX/Y/Z)
    int dist_2_E_X = int(float((t->NX-1)/2)) - dist_S_1; // same distance from edge of tissue as cell 1 is from stimulus site
    int dist_2_E_Y = int(float((t->NY-1)/2)) - dist_S_1;
    int dist_2_E_Z = int(float((t->NZ-1)/2)) - dist_S_1;

    int dist_X = dist_2_E_X - dist_S_1;
    int dist_Y = dist_2_E_Y - dist_S_1;
    int dist_Z = dist_2_E_Z - dist_S_1;
    t->dist_x	= dist_X*t->dx;
    t->dist_y	= dist_Y*t->dy;
    t->dist_z	= dist_Z*t->dz;
    t->dist_xy	= sqrt( (t->dist_x*t->dist_x) + (t->dist_y*t->dist_y) );
    t->dist_xz	= sqrt( (t->dist_x*t->dist_x) + (t->dist_z*t->dist_z) );
    t->dist_yz	= sqrt( (t->dist_y*t->dist_y) + (t->dist_z*t->dist_z) );
    t->dist_xyz = sqrt( (t->dist_x*t->dist_x) + (t->dist_y*t->dist_y) + (t->dist_z*t->dist_z) );

    // cell index for cells to evaluate CV
    // as idealised, location in 1D = linearised cell index; don't need the Ncell index
    // Basically here: CV_x_1 is the first cell to calc CV in x direction, which is offset centre in x, and centre in y and z; second cell is further from stimulus, before edge
    t->CV_x_1	= ((int(float((t->NX-1)/2.0))) + dist_S_1) 	    + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)))); // coords = (NX/2+dist_1, NY/2, NZ/2)  i.e. dist_S1 right of centre stimulus)
    t->CV_x_2	= ((int(float((t->NX-1)/2.0))) + dist_2_E_X) 	+ (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)))); // coords = (NX/2+dist_2, NY/2, NZ/2)  i.e. NX-dist_S1 right of centre, dist_S1 left of edge
    t->CV_y_1	= ((int(float((t->NX-1)/2.0))))		 	        + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));
    t->CV_y_2	= ((int(float((t->NX-1)/2.0))))		 	        + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));
    t->CV_z_1	= ((int(float((t->NX-1)/2.0)))) 	            + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 )); 
    t->CV_z_2	= ((int(float((t->NX-1)/2.0)))) 	            + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z )); 

    // Diagonal directions; in plane
    t->CV_xyp_1  = ((int(float((t->NX-1)/2.0))) + dist_S_1)     + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));  // coords = (NX/2+dist_1, NY/2+dist_1, NZ/2)
    t->CV_xyp_2  = ((int(float((t->NX-1)/2.0))) + dist_2_E_X)   + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));  // coords = (NX/2+dist_2, NY/2+dist_2, NZ/2)

    t->CV_xym_1  = ((int(float((t->NX-1)/2.0))) - dist_S_1)     + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));  // coords = (NX/2-dist_1, NY/2+dist_1, NZ/2)
    t->CV_xym_2  = ((int(float((t->NX-1)/2.0))) - dist_2_E_X)   + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0))));  // coords = (NX/2-dist_2, NY/2+dist_2, NZ/2)

    t->CV_xzp_1  = ((int(float((t->NX-1)/2.0))) + dist_S_1)     + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));
    t->CV_xzp_2  = ((int(float((t->NX-1)/2.0))) + dist_2_E_X)   + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z ));

    t->CV_xzm_1  = ((int(float((t->NX-1)/2.0))) - dist_S_1)     + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));
    t->CV_xzm_2  = ((int(float((t->NX-1)/2.0))) - dist_2_E_X)   + (t->NX* (int(float((t->NY-1)/2.0))))                  + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z ));

    t->CV_yzp_1  = ((int(float((t->NX-1)/2.0))))                + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1))       + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));
    t->CV_yzp_2  = ((int(float((t->NX-1)/2.0))))                + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y))     + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z ));

    t->CV_yzm_1  = ((int(float((t->NX-1)/2.0))))                + (t->NX* (int(float((t->NY-1)/2.0)) - dist_S_1))       + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));
    t->CV_yzm_2  = ((int(float((t->NX-1)/2.0))))                + (t->NX* (int(float((t->NY-1)/2.0)) - dist_2_E_Y))     + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z ));

    // Diagonals; 3D
    t->CV_xyzppp_1 = ((int(float((t->NX-1)/2.0))) + dist_S_1)   + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));   // coords = (NX/2+dist_1, NY/2+dist_1, NZ/2+dist_1)
    t->CV_xyzppp_2 = ((int(float((t->NX-1)/2.0))) + dist_2_E_X) + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z )); // coords = (NX/2+dist_2, NY/2+dist_2, NZ/2+dist_2)

    t->CV_xyzppm_1 = ((int(float((t->NX-1)/2.0))) + dist_S_1)   + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) - dist_S_1 ));   // coords = (NX/2+dist_1, NY/2+dist_1, NZ/2-dist_1)
    t->CV_xyzppm_2 = ((int(float((t->NX-1)/2.0))) + dist_2_E_X) + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) - dist_2_E_Z )); // coords = (NX/2+dist_2, NY/2+dist_2, NZ/2-dist_2)

    t->CV_xyzpmp_1 = ((int(float((t->NX-1)/2.0))) + dist_S_1)   + (t->NX* (int(float((t->NY-1)/2.0)) - dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));   // coords = (NX/2+dist_1, NY/2-dist_1, NZ/2+dist_1)
    t->CV_xyzpmp_2 = ((int(float((t->NX-1)/2.0))) + dist_2_E_X) + (t->NX* (int(float((t->NY-1)/2.0)) - dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z )); // coords = (NX/2+dist_2, NY/2-dist_2, NZ/2+dist_2)

    t->CV_xyzmpp_1 = ((int(float((t->NX-1)/2.0))) - dist_S_1)   + (t->NX* (int(float((t->NY-1)/2.0)) + dist_S_1 ))      + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_S_1 ));   // coords = (NX/2-dist_1, NY/2+dist_1, NZ/2+dist_1)
    t->CV_xyzmpp_2 = ((int(float((t->NX-1)/2.0))) - dist_2_E_X) + (t->NX* (int(float((t->NY-1)/2.0)) + dist_2_E_Y ))    + (t->NX*t->NY* (int(float((t->NZ-1)/2.0)) + dist_2_E_Z )); // coords = (NX/2-dist_2, NY/2+dist_2, NZ/2+dist_2)
}

void set_CV_cells_const_distance(Cell_parameters p, Tissue_parameters *t)
{
    // set radius for starting and end cell
    double radius_S = 4; // mm
    double radius_E = 10;//

    double distance = radius_E - radius_S;
    t->dist_x   = distance;
    t->dist_y   = distance;
    t->dist_z   = distance;
    t->dist_xy  = distance;
    t->dist_xz  = distance;
    t->dist_yz  = distance;
    t->dist_xyz = distance;

    // Set cells to measure CV relative to (starting point, outside of stimulus region)
    // offset for in axes measurment
    double offset_S_axes        = radius_S/t->dx;
    int offset_S_axes_int       = (int)(offset_S_axes+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_S, offset_S_axes, offset_S_axes_int);

    // offset for diagonal measurement
    double offset_S_diag        = sqrt((radius_S*radius_S)/2)/t->dx;
    int offset_S_diag_int       = (int)(offset_S_diag+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_S, offset_S_diag, offset_S_diag_int);

    // offset for double diagonal measurement
    double offset_S_dd    = sqrt((radius_S*radius_S)/3)/t->dx;
    int offset_S_dd_int   = (int)(offset_S_dd+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_S, offset_S_dd, offset_S_dd_int);

    // starting points for cells (offset from stimulus site)
    //                  x                                         y                                      z
    t->CV_x_1       = (t->S1_x_loc + offset_S_axes_int) + (t->NX * t->S1_y_loc)                         + (t->NX*t->NY * t->S1_z_loc); // 1D index is x + NX*y + NZ*NY*z
    t->CV_y_1       = t->S1_x_loc                       + (t->NX * (t->S1_y_loc + offset_S_axes_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_z_1       = t->S1_x_loc                       + (t->NX * t->S1_y_loc)                         + (t->NX*t->NY * (t->S1_z_loc + offset_S_axes_int));

    t->CV_xyp_1     = (t->S1_x_loc + offset_S_diag_int) + (t->NX * (t->S1_y_loc + offset_S_diag_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_xym_1     = (t->S1_x_loc - offset_S_diag_int) + (t->NX * (t->S1_y_loc + offset_S_diag_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_xzp_1     = (t->S1_x_loc + offset_S_diag_int) + (t->NX * (t->S1_y_loc))                       + (t->NX*t->NY * (t->S1_z_loc + offset_S_diag_int));
    t->CV_xzm_1     = (t->S1_x_loc - offset_S_diag_int) + (t->NX * (t->S1_y_loc))                       + (t->NX*t->NY * (t->S1_z_loc + offset_S_diag_int));
    t->CV_yzp_1     = t->S1_x_loc                       + (t->NX * (t->S1_y_loc + offset_S_diag_int))   + (t->NX*t->NY * (t->S1_z_loc + offset_S_diag_int));
    t->CV_yzm_1     = t->S1_x_loc                       + (t->NX * (t->S1_y_loc - offset_S_diag_int))   + (t->NX*t->NY * (t->S1_z_loc + offset_S_diag_int));

    t->CV_xyzppp_1  = (t->S1_x_loc + offset_S_dd_int)   + (t->NX * (t->S1_y_loc + offset_S_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_S_dd_int));
    t->CV_xyzppm_1  = (t->S1_x_loc + offset_S_dd_int)   + (t->NX * (t->S1_y_loc + offset_S_dd_int))     + (t->NX*t->NY * (t->S1_z_loc - offset_S_dd_int));
    t->CV_xyzpmp_1  = (t->S1_x_loc + offset_S_dd_int)   + (t->NX * (t->S1_y_loc - offset_S_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_S_dd_int));
    t->CV_xyzmpp_1  = (t->S1_x_loc - offset_S_dd_int)   + (t->NX * (t->S1_y_loc + offset_S_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_S_dd_int));

    // Set cells to measure CV from (distance is differece between these cells and the starting points)
    // offset for in axes measurment
    double offset_E_axes        = radius_E/t->dx;
    int offset_E_axes_int       = (int)(offset_E_axes+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_E, offset_E_axes, offset_E_axes_int);

    // offset for diagonal measurement
    double offset_E_diag        = sqrt((radius_E*radius_E)/2)/t->dx;
    int offset_E_diag_int       = (int)(offset_E_diag+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_E, offset_E_diag, offset_E_diag_int);

    // offset for double diagonal measurement
    double offset_E_dd    = sqrt((radius_E*radius_E)/3)/t->dx;
    int offset_E_dd_int   = (int)(offset_E_dd+0.5); // so that 0.8 rounds to 1 not 0
    printf("rad %f ncells %f %d\n", radius_E, offset_E_dd, offset_E_dd_int);

    // ending points for cells (offset from stimulus site)
    //                  x                                         y                                      z
    t->CV_x_2       = (t->S1_x_loc + offset_E_axes_int) + (t->NX * t->S1_y_loc)                         + (t->NX*t->NY * t->S1_z_loc); // 1D index is x + NX*y + NZ*NY*z
    t->CV_y_2       = t->S1_x_loc                       + (t->NX * (t->S1_y_loc + offset_E_axes_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_z_2       = t->S1_x_loc                       + (t->NX * t->S1_y_loc)                         + (t->NX*t->NY * (t->S1_z_loc + offset_E_axes_int));

    t->CV_xyp_2     = (t->S1_x_loc + offset_E_diag_int) + (t->NX * (t->S1_y_loc + offset_E_diag_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_xym_2     = (t->S1_x_loc - offset_E_diag_int) + (t->NX * (t->S1_y_loc + offset_E_diag_int))   + (t->NX*t->NY * t->S1_z_loc);
    t->CV_xzp_2     = (t->S1_x_loc + offset_E_diag_int) + (t->NX * (t->S1_y_loc))                       + (t->NX*t->NY * (t->S1_z_loc + offset_E_diag_int));
    t->CV_xzm_2     = (t->S1_x_loc - offset_E_diag_int) + (t->NX * (t->S1_y_loc))                       + (t->NX*t->NY * (t->S1_z_loc + offset_E_diag_int));
    t->CV_yzp_2     = t->S1_x_loc                       + (t->NX * (t->S1_y_loc + offset_E_diag_int))   + (t->NX*t->NY * (t->S1_z_loc + offset_E_diag_int));
    t->CV_yzm_2     = t->S1_x_loc                       + (t->NX * (t->S1_y_loc - offset_E_diag_int))   + (t->NX*t->NY * (t->S1_z_loc + offset_E_diag_int));

    t->CV_xyzppp_2  = (t->S1_x_loc + offset_E_dd_int)   + (t->NX * (t->S1_y_loc + offset_E_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_E_dd_int));
    t->CV_xyzppm_2  = (t->S1_x_loc + offset_E_dd_int)   + (t->NX * (t->S1_y_loc + offset_E_dd_int))     + (t->NX*t->NY * (t->S1_z_loc - offset_E_dd_int));
    t->CV_xyzpmp_2  = (t->S1_x_loc + offset_E_dd_int)   + (t->NX * (t->S1_y_loc - offset_E_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_E_dd_int));
    t->CV_xyzmpp_2  = (t->S1_x_loc - offset_E_dd_int)   + (t->NX * (t->S1_y_loc + offset_E_dd_int))     + (t->NX*t->NY * (t->S1_z_loc + offset_E_dd_int));

}

void calculate_CV(Tissue_parameters t, Model_variables *var, const char* directory)
{
    double CV_x, CV_y, CV_z, CV_xyp, CV_xym, CV_xzp, CV_xzm, CV_yzp, CV_yzm;    // recall that Tissue.CV_x_1/2 returns index of a specified cell
    double CV_xyzppp, CV_xyzppm, CV_xyzpmp, CV_xyzmpp;

    // 2D and 3D
    CV_x = t.dist_x/(var[t.CV_x_2].t_ex - var[t.CV_x_1].t_ex);
    CV_y = t.dist_y/(var[t.CV_y_2].t_ex - var[t.CV_y_1].t_ex);

    CV_xyp = t.dist_xy/(var[t.CV_xyp_2].t_ex - var[t.CV_xyp_1].t_ex);
    CV_xym = t.dist_xy/(var[t.CV_xym_2].t_ex - var[t.CV_xym_1].t_ex);

    // 3D only
    if (strcmp(t.Tissue_order, "3D") == 0)
    {
        CV_z = t.dist_z/(var[t.CV_z_2].t_ex - var[t.CV_z_1].t_ex);

        CV_xzp = t.dist_xz/(var[t.CV_xzp_2].t_ex - var[t.CV_xzp_1].t_ex);
        CV_xzm = t.dist_xz/(var[t.CV_xzm_2].t_ex - var[t.CV_xzm_1].t_ex);

        CV_yzp = t.dist_yz/(var[t.CV_yzp_2].t_ex - var[t.CV_yzp_1].t_ex);
        CV_yzm = t.dist_yz/(var[t.CV_yzm_2].t_ex - var[t.CV_yzm_1].t_ex);

        CV_xyzppp = t.dist_xyz/(var[t.CV_xyzppp_2].t_ex - var[t.CV_xyzppp_1].t_ex);
        CV_xyzppm = t.dist_xyz/(var[t.CV_xyzppm_2].t_ex - var[t.CV_xyzppm_1].t_ex);
        CV_xyzpmp = t.dist_xyz/(var[t.CV_xyzpmp_2].t_ex - var[t.CV_xyzpmp_1].t_ex);
        CV_xyzmpp = t.dist_xyz/(var[t.CV_xyzmpp_2].t_ex - var[t.CV_xyzmpp_1].t_ex);
    }

    // Print to screen
    printf("Conduction velocities:\n\tX = %.3f m/s Y = %.3f m/s Z = %.3f m/s\n", CV_x, CV_y, CV_z);
    printf("\tXYp = %.3f m/s XYm = %.3f m/s XZp = %.3f m/s XZm = %.3f m/s YZp = %.3f m/s YZm  = %.3f m/s\n", CV_xyp, CV_xym, CV_xzp, CV_xzm, CV_yzp, CV_yzm);
    printf("\tXYZppp = %.3f m/s XYZppm = %.3f m/s XYZpmp = %.3f m/s XYZmpp = %.3f m/s\n", CV_xyzppp, CV_xyzppm, CV_xyzpmp, CV_xyzmpp);

    // Save to file
    char * log_reference    = (char*)malloc(500);
    sprintf(log_reference, "%s/Conduction_velocity_log.dat", directory);
    FILE *CV_out;
    CV_out = fopen(log_reference, "a");
    fprintf(CV_out, "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t.Tissue_order, t.D1, t.D1/t.D_AR, t.dx, t.OX, t.OY, t.OZ, CV_x, CV_y, CV_z, CV_xyp, CV_xym, CV_xzp, CV_xzm, CV_yzp, CV_yzm, CV_xyzppp, CV_xyzppm, CV_xyzpmp, CV_xyzmpp);
    fclose(CV_out);
    free(log_reference);
}
// End conduction velocity calculation ==========================================================//|

// Conduction success calculation ===============================================================\\|
void compute_conduction_success(Tissue_parameters t, Model_variables *var, int N, double S2_time, double S2_CL, const char* directory)
{
    int left_ex, right_ex, ex_type;

    // if its most recent excitation is after the time of S2 stimulus, it was successful
    // This is of course not infallable, if your situation is such that at time of S2, S1 is still propagating
    // But works for almost all cases and doesn't seem worth the additional functionality

    // Left success
    if (var[3].t_ex > S2_time-5) left_ex = 1; 
    else left_ex = 0;

    // Right success
    if (var[N - 3].t_ex > S2_time-5) right_ex = 1; 
    else right_ex = 0;

    // Determine excitation type
    ex_type = left_ex + right_ex; // 0 for no conduction, 1 for uni block, 2 for full conduction

    // Print to screen
    printf("S2 loc = %d\t S2 = %f\tleft excitation = %d\tright excitation = %d\tconduction success type = %d\n", t.S2_x_loc, S2_CL, left_ex, right_ex, ex_type);

    // Save to file
    char * log_reference    = (char*)malloc(500);
    sprintf(log_reference, "%s/1D_conduction_success_log.dat", directory);

    FILE *CB_out;
    CB_out = fopen(log_reference, "a");
    fprintf(CB_out, "%d %f %d %d %d\n", t.S2_x_loc, S2_CL, left_ex, right_ex, ex_type);
    fclose(CB_out);

    free(log_reference);
}
// End Conduction success calculation ===========================================================//|

