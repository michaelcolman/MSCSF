// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Global Model file, which contains functions =  //
// which: 1) selects appropriate set and solve functions; =  //
// 2) sets global and selects appropriate model-specific ==  //
// heterogeneity and modulation; 3) determines stimulus and  //
// measurement variables; 4) contains the Luo-Rudy ========  //
// implementation of INa: =================================  //
// Luo-Rudy Circulation Research. 1991;68:1501-1526 =======  //
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

#include "Model.h"
#include "Structs.h"

// Function list ================================================================================\\|
//	Model specific function call functions:
//	    set_parameters_native()
//	    set_parameters_spatial_Ca()
//	    initial_conditions_native()
//	    compute_model_native()
//	    compute_model_integrated()
//	    compute_and_output_current_functions()
//	
//	stimulus:
//	    stimulus_setup()
//	    compute_Istim()
//	
//	Heterogeneity and modulation - global and function call  ** these are where to add new het and mod models **
//	    set_heterogeneity_and_modulation_native()
//	    update_heterogeneity_and_modulation_integrated()
//	    set_MODIFIER_X_Y()  - template for global/common modifier
//	    set_global_Agents()
//	    set_celltype_hAM()
//	    set_ISO_hAM()
//	    set_global_remodelling()
//	    set_remodelling_hAM()
//	    set_mutation_hAM()
//	    set_ACh_global()
//	
//	compute_reversal_potentials()
//	
//	determine_excitation_state()
//	determine_excitation_state_integrated_0D()
//	calculate_flux_integrals()
//	calculate_measurement_properties()
//	
//	run_voltage_clamp()
//	
//	rush_larsen()
//	sigmoid()
//	
//	Luo-rudy INa functions
//		set_INa_LR_rates()
//		update_gates_INa_LR()
//		compute_INa_LR()
// End Function list ============================================================================//|

// Functions to select appropriate specific functions ===========================================\\|
void set_parameters_native(Cell_parameters *p, char const *Model)
{
	// 1 - models with own specific parameters
	if (strcmp(Model, "minimal") == 0)					set_parameters_native_minimal(p);			// lib/Model_minimal.cpp
	else if (strcmp(Model, "hAM_GB") == 0)				set_parameters_native_hAM_GB(p);			// lib/Model_hAM_GB.cpp
	else if (strcmp(Model, "hAM_CRN") == 0)				set_parameters_native_hAM_CRN(p);			// lib/Model_hAM_CRN.cpp
	else if (strcmp(Model, "hAM_NG") == 0)				set_parameters_native_hAM_NG(p);			// lib/Model_hAM_NG.cpp
	else if (strcmp(Model, "hVM_ORD_s") == 0)			set_parameters_native_hVM_ORD_simple(p);	// lib/Model_hVM_ORD_simple.cpp
	else if (strcmp(Model, "hAM_CAZ_s") == 0)			set_parameters_native_hAM_CAZ_simple(p);	// lib/Model_hAM_CAZ_simple.cpp
	else if (strcmp(Model, "dAM_VA") == 0)				set_parameters_native_dAM_VA(p);	        // lib/Model_dAM_VA.cpp
	//else if (strcmp(Model, "speciesCELL_MODEL") == 0)	set_parameters_native_speciesCELL_MODEL(p);	// lib/Model_speciesCELL_MODEL.cpp // NEW MODEL

	// 2 - models which inherit some or all parameters from other models
	else if (strcmp(Model, "hAM_MT") == 0)					
	{
		set_parameters_native_hAM_NG(p);															// lib/Model_hAM_NG.cpp
		update_parameters_native_hAM_MT(p);															// lib/Model_hAM_MT.cpp
	}
	else if (strcmp(Model, "hAM_WL_CRN") == 0 || strcmp(Model, "hAM_CRN_mWL") == 0)					
	{
		set_parameters_native_hAM_CRN(p);          													// lib/Model_hAM_CRN.cpp
		update_parameters_native_hAM_WL(p);        													// lib/Model_hAM_WL.cpp
	}
	else if (strcmp(Model, "hAM_WL_GB") == 0 || strcmp(Model, "hAM_GB_mWL") == 0)
	{
		set_parameters_native_hAM_GB(p);	 														// lib/Model_hAM_GB.cpp
		update_parameters_native_hAM_WL(p);															// lib/Model_hAM_WL.cpp
	}
	else if (strcmp(Model, "hAM_NG_mWL") == 0)
	{
		set_parameters_native_hAM_NG(p);															// lib/Model_hAM_NG.cpp
		update_parameters_native_hAM_WL(p);															// lib/Model_hAM_WL.cpp
	}
	//else if (strcmp(Model, "speciesCELL_MODEL") == 0)
	//{
	//    set_parameters_native_X(p);                                                               // lib/Model_X.cpp
	//    update_parameters_native_speciesCELL_MODEL(p);                                            // lib/Model_speciesCELL_MODEL.cpp
	//}
	else
	{
		printf("ERROR: \"%s\" is not a valid model type, parameters cannot be set\n\n", Model);
		exit(1);
	}
}

void set_parameters_spatial_Ca(Cell_parameters *p, char const *Model)
{
	// Default parameters for spatial handling model are already set; this is obnly required for speicific models which have different parameters
	if (strcmp(Model, "minimal") == 0);			// do nothing -> same as default
	else if (strcmp(Model, "hVM_ORD_s") == 0); 	// do nothing -> same as default
	else if (strcmp(Model, "hAM_CAZ_s") == 0)	update_parameters_integrated_hAM_CAZ_simple(p);		// lib/Model_hAM_CAZ_simple.cpp
	//else if (strcmp(Model, "speciesCELL_MODEL") == 0)		update_parameters_integrated_speciesCELL_model(p);	// lib/Model_speciesCELL_MODEL.cpp
	else
	{
		printf("ERROR: \"%s\" is not a valid spatial Ca model type. See \"set_parameters_spatial_Ca()\" in \"lib/Model.c\" for options\n\n", Model);
		exit(1);
	}
}

void initial_conditions_native(State_variables *s, Cell_parameters p, char const *Model)
{
	if (strcmp(Model, "minimal") == 0)      				initial_conditions_native_minimal(s, p);    // lib/Model_minimal.cpp
	else if (strcmp(Model, "hAM_CRN") == 0)  				initial_conditions_native_hAM_CRN(s, p);    // lib/Model_hAM_CRN.cpp
	else if (strcmp(Model, "hAM_GB") == 0)  				initial_conditions_native_hAM_GB(s, p);     // lib/Model_hAM_GB.cpp
	else if (strcmp(Model, "hAM_NG") == 0)  				initial_conditions_native_hAM_NG(s, p);    	// lib/Model_hAM_NG.cpp
	else if (strcmp(Model, "hAM_MT") == 0)  				initial_conditions_native_hAM_MT(s, p);    	// lib/Model_hAM_MT.cpp
	else if (strcmp(Model, "hAM_WL_CRN") == 0) 				initial_conditions_native_hAM_WL(s, p);   	// lib/Model_hAM_WL.cpp
	else if (strcmp(Model, "hAM_CRN_mWL") == 0)				initial_conditions_native_hAM_WL(s, p);    	// lib/Model_hAM_WL.cpp
	else if (strcmp(Model, "hAM_WL_GB") == 0) 				initial_conditions_native_hAM_WL(s, p);     // lib/Model_hAM_WL.cpp
	else if (strcmp(Model, "hAM_GB_mWL") == 0) 				initial_conditions_native_hAM_WL(s, p);     // lib/Model_hAM_WL.cpp
	else if (strcmp(Model, "hAM_NG_mWL") == 0) 				initial_conditions_native_hAM_WL(s, p);    	// lib/Model_hAM_WL.cpp
	else if (strcmp(Model, "hVM_ORD_s") == 0) 				initial_conditions_native_hVM_ORD_simple(s, p);   	// lib/Model_hVM_ORD_simple.cpp
	else if (strcmp(Model, "hAM_CAZ_s") == 0) 				initial_conditions_native_hAM_CAZ_simple(s, p);   	// lib/Model_hAM_CAZ_simple.cpp
	else if (strcmp(Model, "dAM_VA") == 0) 					initial_conditions_native_dAM_VA(s, p);     		// lib/Model_dAM_VA.cpp // NEW MODEL
	//else if (strcmp(Model, "speciesCELL_MODEL") == 0) 	initial_conditions_native_speciesCELL_MODEL(s, p);   // lib/Model_speciesCELL_MODEL.cpp // NEW MODEL
	else
	{
		printf("ERROR: \"%s\" is not a valid model type, initial conditions cannot be set\n\n", Model);
		exit(1);
	}
}

void compute_model_native(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	if (strcmp(p.Model, "minimal") == 0)      					compute_model_minimal_native(p, var, s, Vm, dt); 		// lib/Model_minimal.cpp
	else if (strcmp(p.Model, "hAM_CRN") == 0)   				compute_model_hAM_CRN_native(p, var, s, Vm, dt); 		// lib/Model_hAM_CRN.cpp
	else if (strcmp(p.Model, "hAM_GB") == 0)    				compute_model_hAM_GB_native(p, var, s, Vm, dt); 		// lib/Model_hAM_GB.cpp
	else if (strcmp(p.Model, "hAM_NG") == 0)   					compute_model_hAM_NG_native(p, var, s, Vm, dt); 		// lib/Model_hAM_NG.cpp
	else if (strcmp(p.Model, "hAM_MT") == 0)   					compute_model_hAM_MT_native(p, var, s, Vm, dt); 		// lib/Model_hAM_MT.cpp
	else if (strcmp(p.Model, "hAM_WL_CRN") == 0)   				compute_model_hAM_WL_native(p, var, s, Vm, dt); 		// lib/Model_hAM_WL.cpp
	else if (strcmp(p.Model, "hAM_CRN_mWL") == 0)   			compute_model_hAM_WL_native(p, var, s, Vm, dt); 		// lib/Model_hAM_WL.cpp
	else if (strcmp(p.Model, "hAM_WL_GB") == 0)    				compute_model_hAM_WL_native(p, var, s, Vm, dt); 		// lib/Model_hAM_WL.cpp
	else if (strcmp(p.Model, "hAM_GB_mWL") == 0)    			compute_model_hAM_WL_native(p, var, s, Vm, dt); 		// lib/Model_hAM_WL.cpp
	else if (strcmp(p.Model, "hAM_NG_mWL") == 0) 				compute_model_hAM_WL_native(p, var, s, Vm, dt); 		// lib/Model_hAM_WL.cpp
	else if (strcmp(p.Model, "dAM_VA") == 0)   					compute_model_dAM_VA_native(p, var, s, Vm, dt);	// lib/Model_dAM_VA.cpp // NEW MODEL
	//else if (strcmp(p.Model, "speciesCELL_MODEL") == 0)   	compute_model_speciesCELL_MODEL_native(p, var, s, Vm, dt);	// lib/Model_speciesCELL_MODEL.cpp // NEW MODEL
	else
	{
		printf("ERROR: \"%s\" is not a valid model type, model cannot be computed. See \"compute_model_native()\" in \"lib/Model.c\" for options\n\n", p.Model);
		exit(1);
	}
}

void compute_model_integrated(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	if (strcmp(p.Model, "minimal") == 0)                        compute_model_minimal_integrated(p, var, s, Vm, dt);        // lib/Model_minimal.cpp
	else if (strcmp(p.Model, "hVM_ORD_s") == 0)                 compute_model_hVM_ORD_simple_integrated(p, var, s, Vm, dt); // lib/Model_hVM_ORD_simple.cpp
	else if (strcmp(p.Model, "hAM_CAZ_s") == 0)                 compute_model_hAM_CAZ_simple_integrated(p, var, s, Vm, dt); // lib/Model_hAM_CAZ_simple.cpp
    //else if (strcmp(p.Model, "speciesCELL_MODEL") == 0)         compute_model_speciesCELL_MODEL_integrated(p, var, s, Vm, dt); // lib/Model_speciesCELL_MODEL.cpp // NEW MODEL
	else
	{
		printf("ERROR: \"%s\" is not a valid model type, model cannot be computed. See \"compute_model_integrated()\" in \"lib/Model.c\" for options\n\n", p.Model);
		exit(1);
	}
}

void compute_and_output_current_functions(Cell_parameters p, Model_variables *var, char const *directory)
{
	// This function outputs the conductance/flux rate and the voltage dependence
	// variables for the ion currents, as used in the simulation under all celltype
	// and modualtions employed, as well as the final values of scaling a shifts
	// (with all modifications applied sequentially)
	// This function is used to output all the voltage dependent variables as 
	// used in the simulation (i.e. with all modifications applied).
	// These can be used to check correct implementations and for figures.
	// Create files for each current for which you want to output variables 
	char * filename       = (char*)malloc(500);

	// Calculate and output all voltage-dependent functions for ion current gating
	// and calculation variables
	FILE * INa_out;
	FILE * INaL_out;
	FILE * Ito_out;
	FILE * IKur_out;
	FILE * ICaL_out;
	FILE * IKr_out;
	FILE * IKs_out;
	FILE * IKACh_out;
	FILE * If_out;
	sprintf(filename, "%s/INa.txt", directory);
	INa_out = fopen(filename, "wt");
	sprintf(filename, "%s/INaL.txt", directory);
	INaL_out = fopen(filename, "wt");
	sprintf(filename, "%s/Ito.txt", directory);
	Ito_out = fopen(filename, "wt");
	sprintf(filename, "%s/IKur.txt", directory);
	IKur_out = fopen(filename, "wt");
	sprintf(filename, "%s/ICaL.txt", directory);
	ICaL_out = fopen(filename, "wt");
	sprintf(filename, "%s/IKr.txt", directory);
	IKr_out = fopen(filename, "wt");
	sprintf(filename, "%s/IKs.txt", directory);
	IKs_out = fopen(filename, "wt");
	sprintf(filename, "%s/IKACh.txt", directory);
	IKACh_out = fopen(filename, "wt");
	sprintf(filename, "%s/If.txt", directory);
	If_out = fopen(filename, "wt");

	// Loop over voltages and calculate functions
	double Vm;
	for (Vm = -100; Vm < 100; Vm+=1.0)
	{

		// Call set gate rates function to calculate variables at Voltage
		if (strcmp(p.Model, "minimal") == 0)                      set_gate_rates_minimal_native(p, var, Vm);   	 			// lib/Model_minimal.cpp
		else if (strcmp(p.Model, "hAM_CRN") == 0)                 set_gate_rates_hAM_CRN_native(p, var, Vm, p.Cai);    		// lib/Model_hAM_CRN.cpp
		else if (strcmp(p.Model, "hAM_GB") == 0)                  set_gate_rates_hAM_GB_native(p, var, Vm, p.Cai);     		// lib/Model_hAM_GB.cpp
		else if (strcmp(p.Model, "hAM_NG") == 0)                  set_gate_rates_hAM_NG_native(p, var, Vm, p.Cai, p.Ko); 	// lib/Model_hAM_NG.cpp
		else if (strcmp(p.Model, "hAM_MT") == 0)                  set_gate_rates_hAM_MT_native(p, var, Vm, p.Cai, p.Ko);	// lib/Model_hAM_MT.cpp
		else if (strcmp(p.Model, "hAM_WL_CRN") == 0)              set_gate_rates_hAM_WL_native(p, var, Vm, p.Cai);    		// lib/Model_hAM_WL.cpp
		else if (strcmp(p.Model, "hAM_CRN_mWL") == 0)             set_gate_rates_hAM_WL_native(p, var, Vm, p.Cai);    		// lib/Model_hAM_WL.cpp
		else if (strcmp(p.Model, "hAM_WL_GB") == 0)               set_gate_rates_hAM_WL_native(p, var, Vm, p.Cai);     		// lib/Model_hAM_WL.cpp
		else if (strcmp(p.Model, "hAM_GB_mWL") == 0)              set_gate_rates_hAM_WL_native(p, var, Vm, p.Cai);    		// lib/Model_hAM_WL.cpp
		else if (strcmp(p.Model, "hAM_NG_mWL") == 0)              set_gate_rates_hAM_WL_native(p, var, Vm, p.Cai);     		// lib/Model_hAM_WL.cpp
		else if (strcmp(p.Model, "hVM_ORD_s") == 0)               set_gate_rates_hVM_ORD_simple_native(p, var, Vm, p.Cai);  // lib/Model_hVM_ORD_simple.cpp
		else if (strcmp(p.Model, "hAM_CAZ_s") == 0)               set_gate_rates_hAM_CAZ_simple_native(p, var, Vm, p.Cai);  // lib/Model_hAM_CAZ_simple.cpp
		else if (strcmp(p.Model, "dAM_VA") == 0)                  set_gate_rates_dAM_VA_native(p, var, Vm, p.Cai);     		// lib/Model_dAM_VA.cpp // NEW MODEL
		//else if (strcmp(Model, "speciesCELL_MODEL") == 0)     set_gate_rates_speciesCELL_MODEL_native(p, var, Vm, p.Cai); // lib/Model_speciesCELL_MODEL.cpp // NEW MODEL

		// INa
		fprintf(INa_out,    "%f %f %f %f %f ", 	Vm, var->INa_va_al, var->INa_va_bet, var->INa_va_ss, var->INa_va_tau);
		fprintf(INa_out, 	"%f %f %f %f ", 	var->INa_vi_1_al, var->INa_vi_1_bet, var->INa_vi_1_ss, var->INa_vi_1_tau);
		fprintf(INa_out, 	"%f %f %f %f\n",   	var->INa_vi_2_al, var->INa_vi_2_bet, var->INa_vi_2_ss, var->INa_vi_2_tau);

		// INaL
		fprintf(INaL_out,    "%f %f %f %f %f ", Vm, var->INaL_va_al, var->INaL_va_bet, var->INaL_va_ss, var->INaL_va_tau);
		fprintf(INaL_out,    "%f %f %f %f ",    var->INaL_vi_al, var->INaL_vi_bet, var->INaL_vi_ss, var->INaL_vi_tau);

		// Ito
		fprintf(Ito_out,    "%f %f %f %f %f ",  Vm, var->Ito_va_ss, var->Ito_va_tau, var->Ito_vi_ss, var->Ito_vi_tau);
		fprintf(Ito_out, 	"%f %f %f %f\n", 	var->Ito_vi_3_ss, var->Ito_vi_s_tau, var->Ito_vi_3_tau, var->Ito_vi_Fs);

		// ICaL
		fprintf(ICaL_out,   "%f %f %f %f %f %f\n", Vm, var->ICaL_va_ss, var->ICaL_va_tau, var->ICaL_vi_ss, var->ICaL_vi_tau, var->ICaL_vi_s_tau);

		// IKur
		fprintf(IKur_out,   "%f %f %f %f %f %f\n",  Vm, var->IKur_va_ss, var->IKur_va_tau, var->IKur_vi_ss, var->IKur_vi_tau, var->IKur_dynamic_g);

		// IKr
		fprintf(IKr_out,    "%f %f %f %f\n",     Vm, var->IKr_va_ss, var->IKr_va_tau, var->IKr_vi_ti);

		// IKs
		fprintf(IKs_out,    "%f %f %f %f\n",	Vm, var->IKs_va_ss, var->IKs_va_tau, var->IKs_va_2_tau);

		// IKACh
		fprintf(IKACh_out,    "%f %f %f %f\n",	Vm, var->IKACh_va_ss, var->IKACh_va_tau, var->IKACh_v_ti);

		// If
		fprintf(If_out,    "%f %f %f\n",	Vm, var->If_va_ss, var->If_va_tau);

	}

	fclose(INa_out);
	fclose(INaL_out);
	fclose(Ito_out);
	fclose(IKur_out);
	fclose(ICaL_out);
	fclose(IKr_out);
	fclose(IKs_out);
	fclose(IKACh_out);
	fclose(If_out);
	// End compute voltage dependent variables ==============================//|

	// Output the final parameters for conductance/maximal flux rates =======\\|
	FILE * max_out;
	sprintf(filename, "%s/Magnitude_parameters.txt", directory);
	max_out = fopen(filename, "wt");
	fprintf(max_out, "gNa %f\n", 		p.gNa 	* p.GNa);	
	fprintf(max_out, "gNaL %f\n", 		p.gNaL 	* p.GNaL);	
	fprintf(max_out, "gto %f\n", 		p.gto 	* p.Gto);	
	fprintf(max_out, "gCaL %f\n", 		p.gCaL 	* p.GCaL);	
	fprintf(max_out, "pCaL %f\n", 		p.pCaL 	* p.GCaL);	
	fprintf(max_out, "gKur %f\n", 		p.gKur 	* p.GKur);	
	fprintf(max_out, "gKr %f\n", 		p.gKr 	* p.GKr);	
	fprintf(max_out, "gKs %f\n", 		p.gKs 	* p.GKs);	
	fprintf(max_out, "gK1 %f\n", 		p.gK1 	* p.GK1);	
	fprintf(max_out, "gKACh %f\n", 		p.gKACh	* p.GKACh);	
	fprintf(max_out, "INCX_max %f\n", 	p.INCX_bar		* p.GNCX);	
	fprintf(max_out, "INaK_max %f\n", 	p.INaK_bar		* p.GNaK);	
	fprintf(max_out, "ICaP_max %f\n", 	p.ICaP_bar		* p.GCaP);	
	fprintf(max_out, "ICab_max %f\n", 	p.ICab_bar		* p.GCab);	
	fprintf(max_out, "J_rel_max %f\n", 	p.J_rel_max 	* p.Grel);	
	fprintf(max_out, "J_up_max %f\n", 	p.J_SERCA_max 	* p.Gup);	
	fprintf(max_out, "J_leak_max %f\n", p.J_leak_max 	* p.Gleak);	
	fclose(max_out);
	// End Output the final parameters for conductance/maximal flux rates ===//|

	// Shifts and scale factors =============================================\\|
	FILE *modifiers_out;
	sprintf(filename, "%s/Modifier_parameters.txt", directory);
	modifiers_out = fopen(filename, "wt");
	fprintf(modifiers_out,"INa_scale %.8f\nINaL_scale %.8f\nIto_scale %.8f\nICaL_scale %.8f\nIKur_scale %.8f\n", p.GNa, p.GNaL, p.Gto, p.GCaL, p.GKur);
	fprintf(modifiers_out,"IKr_scale %.8f\nIKs_scale  %.8f\nIK1_scale %.8f\nINCX_scale %.8f\nICaP_scale %.8f\n", p.GKr, p.GKs, p.GK1, p.GNCX, p.GCaP);
	fprintf(modifiers_out,"INab_scale %.8f\nICab_scale %.8f\nIKb_scale %.8f\nINaK_scale %.8f\nIClCa_scale %.8f\n", p.GNab, p.GCab, p.GKb, p.GNaK, p.GClCa);
	fprintf(modifiers_out,"Jup_scale %.8f\nJleak_scale %.8f\nJrel_scale %0.8f\n", p.Gleak, p.Gup, p.Grel);

	fprintf(modifiers_out,"INa_va_tau_scale %.8f\nINa_vi_1_tau_scale %.8f\nINa_vi_2_tau_scale %.8f\nINaL_va_tau_scale %.8f\nINaL_vi_tau_scale %.8f\n", p.INa_va_tau_scale, p.INa_vi_1_tau_scale, p.INa_vi_2_tau_scale, p.INaL_va_tau_scale, p.INaL_vi_tau_scale);
	fprintf(modifiers_out,"Ito_va_tau_scale %.8f\nIto_vi_tau_scale  %.8f\nICaL_va_tau_scale %.8f\nICaL_vi_tau_scale %.8f\n", p.Ito_va_tau_scale, p.Ito_vi_tau_scale, p.ICaL_va_tau_scale, p.ICaL_vi_tau_scale);
	fprintf(modifiers_out,"IKur_va_tau_scale %.8f\nIKur_vi_tau_scale %.8f\nIKr_va_tau_scale  %0.8f\nIKs_va_tau_scale %0.8f\n", p.IKur_va_tau_scale, p.IKur_vi_tau_scale, p.IKr_va_tau_scale, p.IKs_va_tau_scale);

	fprintf(modifiers_out,"INa_va_shift    %.8f\nINa_vi_shift     %.8f\nINaL_va_shift   %.8f\nINaL_vi_shift    %.8f\n", p.INa_va_shift, p.INa_vi_shift, p.INaL_va_shift, p.INaL_vi_shift);
	fprintf(modifiers_out,"Ito_va_ss_shift %.8f\nIto_va_tau_shift %.8f\nIto_vi_ss_shift %.8f\nIto_vi_tau_shift %.8f\n", p.Ito_va_ss_shift, p.Ito_va_tau_shift, p.Ito_vi_ss_shift, p.Ito_vi_tau_shift);
	fprintf(modifiers_out,"ICaL_va_ss_shift %.8f\nICaL_va_tau_shift %.8f\nICaL_vi_ss_shift %.8f\nICaL_vi_tau_shift %.8f\n", p.ICaL_va_ss_shift, p.ICaL_va_tau_shift, p.ICaL_vi_ss_shift, p.ICaL_vi_tau_shift);
	fprintf(modifiers_out,"IKur_va_ss_shift %.8f\nIKur_va_tau_shift %.8f\nIKur_vi_ss_shift %.8f\nIKur_vi_tau_shift %.8f\n", p.IKur_va_ss_shift, p.IKur_va_tau_shift, p.IKur_vi_ss_shift, p.IKur_vi_tau_shift);
	fprintf(modifiers_out,"IKr_va_ss_shift %.8f\nIKr_va_tau_shift %.8f\nIKr_vi_ss_shift %.8f\nIKs_va_ss_shift %.8f\nIKs_va_tau_shift %.8f\nIK1_va_shift    %.8f\n", p.IKr_va_ss_shift, p.IKr_va_tau_shift, p.IKr_vi_ss_shift, p.IKs_va_ss_shift, p.IKs_va_tau_shift, p.IK1_va_shift);

	fprintf(modifiers_out,"Ito_va_ss_kscale %.8f\nIto_vi_ss_kscale %.8f\nICaL_va_ss_kscale %.8f\nICaL_vi_ss_kscale %.8f\nIKur_va_ss_kscale %.8f\nIKur_vi_ss_kscale %.8f\n", p.Ito_va_ss_kscale, p.Ito_vi_ss_kscale, p.ICaL_va_ss_kscale, p.ICaL_vi_ss_kscale, p.IKur_va_ss_kscale, p.IKur_vi_ss_kscale);
	fprintf(modifiers_out,"IKr_va_ss_kscale %.8f\nIKr_vi_ss_kscale %.8f\nIKs_va_ss_kscale %.8f\n", p.IKr_va_ss_kscale, p.IKr_vi_ss_kscale, p.IKs_va_ss_kscale);
	fclose(modifiers_out);
	// End Shifts and scale factors =========================================//|
}
// End Functions to select appropriate specific functions =======================================//|

// Stimulus current =============================================================================\\|
void stimulus_setup(Cell_parameters p, Model_variables *var, double dt, int BCL, int S2, int Paced_time)
{
	var->stimduration_int		= p.stimduration*(1/dt);
	var->BCL_int				= BCL * (int)(1/dt);
	var->stimflag				= false;
	var->S2_stimflag			= false;
	var->stimcount				= 0;
	var->stimcount_S2			= 0;
	var->Istim					= 0.0;
	var->Istim_S2				= 0.0;
	var->S2_int					= S2 * (int)(1/dt);
	var->Paced_time_int			= Paced_time * (int)(1/dt);
}

void compute_Istim(Cell_parameters p, Model_variables *var, double Paced_time, double S2_time, double time, int time_int)
{
	// S1
	if( (time_int == 0 || time_int % var->BCL_int == 0) && time < Paced_time)
	{
		var->stimflag			= true;
		var->ex_switch     		= 0;	// To reset measurement values, even if AP not < rep threshold at time of stimulus
	}

	if (var->stimflag == true)
	{
		var->Istim				= p.stimmag;
		var->stimcount++;
		if (var->stimcount >= var->stimduration_int)
		{
			var->stimflag		= false;
			var->stimcount		= 0;
		} // end stimcount > duration_int IF
	} // end stimflag = true IF
	else var->Istim				= 0.0;

	// S2
	if (time > Paced_time && var->S2_int > 0)
	{
		if( ((time_int - (var->Paced_time_int-5) )% var->S2_int == 0) && time < S2_time)
		{
			var->S2_stimflag           = true;
			var->ex_switch     			= 0;	// To reset measurement values, even if AP not < rep threshold at time of stimulus
		}

		if (var->S2_stimflag == true)
		{
			var->Istim_S2           = p.stimmag;
			var->stimcount_S2++;
			if (var->stimcount_S2 >= var->stimduration_int)
			{
				var->S2_stimflag       = false;
				var->stimcount_S2      = 0;
			} // end stimcount > duration_int IF
		} // end stimflag = true IF
		else var->Istim_S2            = 0.0;
	}
}
// End Stimulus current =========================================================================//|

// Current modification variables | Het and modulation ==========================================\\|
// Global selection function ==========================================================\\|
void set_heterogeneity_and_modulation_native(Cell_parameters *p)
{
	// Sets global or common het and/or modulation (functions below in this file)
	// then calls model specific functions for single-model het and/or modulation if global functions did not actually set the parameter
	// (All "X_set_ref" start off = 0; as global/common function is called, set to 1 to indicate function call. If the specific option
	// is not caught by one IF statement in global/common function (i.e. not been set) the "X_set_ref" set back to 0
	// and specific functions are called. If not specific function, error is returned; if option not caught by IF within specific function
	// error also returned (within specific function found in Model_X.cpp))
	// Only one value may be contained in each type (celltype, ISO/ISO_model, Agent, Remodelling, Mutation) at a time; multiple different
	// types may be specified simulataneously

	// First, set all set_refs to 0
	p->Het_set_ref = p->ISO_set_ref = p->Agent_set_ref = p->Remodelling_set_ref = p->Mutation_set_ref = p->ACh_set_ref = 0; 

	// Global
	set_global_Agents(p);				// Sets any pharma agents which apply to all models
	set_global_remodelling(p);			// Sets any remodelling which applies to all models
	if (p->ACh > 0.0) set_ACh_global(p);// Sets ACh which applies to all models
	//set_MODIFIER_X_Y(Cell_parameters *p);

	// common: human atrial models
	if (p->hAM == true && p->Het_set_ref == 0)    				set_celltype_hAM(p);
	if (p->ISO > 0.0 && p->hAM == true &&  p->ISO_set_ref == 0)	set_ISO_hAM(p);
	if (p->hAM == true && p->Remodelling_set_ref == 0) 			set_remodelling_hAM(p); 
	if (p->hAM == true && p->Mutation_set_ref == 0) 			set_mutation_hAM(p);
	// could have list of specific model strings if similar p->hAM variable does not exist for your model group.

	// Model-specific || note: each of these functions checks the "set_ref" variales before calling model-specific function
	if (strcmp(p->Model, "minimal") == 0)			set_het_mod_minimal(p);			// lib/Model_minimal.cpp
	else if (strcmp(p->Model, "hAM_CRN") == 0)		set_het_mod_hAM_CRN(p);			// lib/Model_hAM_CRN.cpp
	else if (strcmp(p->Model, "hAM_GB") == 0)		set_het_mod_hAM_GB(p);			// lib/Model_hAM_GB/cpp
	else if (strcmp(p->Model, "hAM_NG") == 0)		set_het_mod_hAM_NG(p);			// lib/Model_hAM_NG.cpp
	else if (strcmp(p->Model, "hAM_MT") == 0)		set_het_mod_hAM_MT(p);			// lib/Model_hAM_MT.cpp
	else if (strcmp(p->Model, "hAM_WL_CRN") == 0) 	set_het_mod_hAM_WL(p);  		// lib/Model_hAM_WL.cpp
	else if (strcmp(p->Model, "hAM_CRN_mWL") == 0)	set_het_mod_hAM_WL(p); 			// lib/Model_hAM_WL.cpp
	else if (strcmp(p->Model, "hAM_WL_GB") == 0) 	set_het_mod_hAM_WL(p);			// lib/Model_hAM_WL.cpp
	else if (strcmp(p->Model, "hAM_GB_mWL") == 0)	set_het_mod_hAM_WL(p); 			// lib/Model_hAM_WL.cpp
	else if (strcmp(p->Model, "hAM_NG_mWL") == 0)	set_het_mod_hAM_WL(p); 			// lib/Model_hAM_WL.cpp
	else if (strcmp(p->Model, "hVM_ORD_s") == 0) 	set_het_mod_hVM_ORD_simple(p); 	// lib/Model_hVM_ORD_simple.cpp
	else if (strcmp(p->Model, "hAM_CAZ_s") == 0)	set_het_mod_hAM_CAZ_simple(p); 	// lib/Model_hAM_CAZ_simple.cpp
	else if (strcmp(p->Model, "dAM_VA") == 0) 		set_het_mod_dAM_VA(p); 	// lib/Model_dAM_VA.cpp
	//else if (strcmp(p->Model, "speciesCELL_MODEL") == 0) set_het_mod_speciesCELL_MODEL(p); 	// lib/Model_speciesCELL_MODEL.cpp NEW

	// Spatial gradient (global as will be geo and type, rather than model, dependent)
	// For any non-celltype depedent, global continuous heterogeneity 
	// e.g. apico-basal, distance from SAN; want to apply on top of
	// all previous celltype/modulation settings
	set_spatial_gradient(p);
}

void update_heterogeneity_and_modulation_integrated(Cell_parameters *p)
{
	// Only further modifications to above which are different for integrated compared to native
    //if (strcmp(p->Model, "speciesCELL_MODEL") == 0)	update_het_and_mod_speciesCELL_MODEL_integrated(p);	// lib/Model_speciesCELL_MODEL.cpp
}
// End Global selection function ======================================================//|

// Global or common functions =========================================================\\|
// Template ==============================\\|
// Copy and rename this for any modifier function
void set_MODIFIER_X_Y(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->{ACh/ISO}" is concentration of {ACh/ISO}; "p->{ACh/ISO}_model" is the {ACh/ISO} model to be applied
	// "p->{Remodelling/Agent}_prop" is the proportion to maximal value of Remodellibng/Agent; 
	// "p->{Remodelling/Agent}" is the model to apply
	// ===========================================
	// This approach assumes linear transition between X = 0 and full effect at X = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double X_scaled = 1/ ( 1 + exp((-4.362*(log(p->X) + 3.27)) ) )  and use "X_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	// Start =====
	// Note: setting X_set_ref to 1 and then 0 if not set is important for error checking
	// and for being able to jump into global/common AND specific functions properly
	// p->X_set_ref = 1; // Indicates X(=ISO/Agent/Celltype/Remodeling/Mutation) has been set
	// if (strcmp(p->X, "SPECIFIER") == 0) { do stuff }
	// else if (strcmp(p->X, "SPECIFIER2") == 0) { do stuff }

	// Replace below with correct scale (ISO/ACh or Remodelling/Agent prop) or no prop for mutation, celltype
	//else if (strcmp(p->modifier, "test") == 0)
	//{
	//  p->Gto              *= (1.0 + p->mod_prop*(3.0-1.0));    // At maximal effect, x3, 
	//  p->IKr_va_ss_kscale *= (1.0 + p->mod_prop*(0.5-1.0));    // At maximal effect, x0.5
	//  p->ICaL_vi_ss_shift += p->mod_prop*5;                    // At maximal effect, + 5
	//  p->Ito_va_tau_scale *= (1.0 + p->mod_prop*(2.3-1.0));    // At maximal effect x2.3
	//  p->Gup              *= (1.0 + p->mod_prop*(2.5-1.0));    // At maximal effect, x2.5
	//	p->Ko				= 5.4; // explicitly set extracellular potassium

	// RyR and LTCC slightly more complex:
	//    p->GCaL, p->Grel *= x ; // use GCaL and Grel for expression not activity
	//	  p->GLTCC_kva1_va2, p->GRyR_kCO *= x // use for channel activity	  
	//	  Both the above multiply maximum flux rate/channel conductance for native models
	//	  For integrated models, they do different things.
	//}

	// etc
	// else p->X_set_ref = 0; // Indicates X=(Mutation/Agent/Celltype/Remodeling/Mutation) has not actually been set
	// End =======
}
// End template ==========================//|

// Pharma agents =========================\\|
void set_global_Agents(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->Agent_prop" is proportion of maximal Agent effect; "p->Agent" is the Agent to be applied
	// ===========================================
	// This approach assumes linear transition between Agent_prop = 0 and full effect at Agent_prop = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double Agent_scaled = 1/ ( 1 + exp((-4.362*(log(p->Agent_prop) + 3.27)) ) )  and use "Agent_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	p->Agent_set_ref			= 1; // Indicates Agent has been set (is set back to zero if not caught by any below IF statement)

	if (strcmp(p->Agent, "MC-II-157c") == 0)
	{
		// References:  Colman et al. 2017 Front. Physiol 8; Guo et al. 2014 PLOS One 9, e105553
		p->IKr_va_ss_shift      += 14 * p->Agent_prop;  			// 14 is maximum effect
		p->GKr                  *= (1.0 + p->Agent_prop*(0.88-1));	// x 0.88 is maximum effect
	}
	// testing exmaple illustration of implementation
	else if (strcmp(p->Agent, "test_global") == 0)
	{
		p->Gto              *= (1.0 + p->Agent_prop*(2.0  -1));    	// Scale factor = Multiplies previous settings  | 2 is scale factor, 2-1 is additive scale factor
		p->IKr_va_ss_kscale *= (1.0 + p->Agent_prop*(1.25 -1)); 	// Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->IKs_va_tau_scale *= (1.0 + p->Agent_prop*(0.75 -1));    	// Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= (1.0 + p->Agent_prop*(2.0  -1));     // Scales intracellular upatke rate. Multiplies previous settings.
		p->ICaL_vi_ss_shift += 5*p->Agent_prop;       				// Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->gKur             = 0.003;    							// Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods.
		p->Ko               = 5.4; 									// explicitly set extracellular potassium 
	}
	// Add new pharmacological agent here: else if (strcmp(p->Agent, "X") == 0) {   }

	else p->Agent_set_ref		= 0;	// Sets back to zero to indicate Agent has not been set despite function call
}
// End Pharma agents =====================//|

// Heterogeneity =========================\\|
void set_celltype_hAM(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ================================================================================//|

	p->Het_set_ref		= 1; // Indicates heterogeneity has been set (is set back to zero if not caught by any below IF statement)

	// Heterogeneity according to Colman et al. 2013 J Physiol 591(17):4249-72
	if (strcmp(p->Celltype, "RA") == 0); // do nothing as default (here for error checking)
	else if (strcmp(p->Celltype, "PM") == 0)
	{
		p->GCaL     *= 0.94;
	}
	else if (strcmp(p->Celltype, "CT") == 0)
	{
		p->GCaL     *= 1.68;
		p->Gto      *= 1.35;
	}
	else if (strcmp(p->Celltype, "RAA") == 0)
	{
		p->GCaL     			*= 1.68;
		p->Gto      			*= 0.4;
		p->IKur_va_ss_shift 	+= 14;
		p->IKur_vi_ss_shift 	+= -25;
		p->IKur_vi_ss_kscale 	*= 1.769;
	}
	else if (strcmp(p->Celltype, "AVR") == 0)
	{
		p->GCaL     *= 0.67;
		p->Gto      *= 0.6;
		p->GKr		*= 1.63;
	}
	else if (strcmp(p->Celltype, "BB") == 0)
	{
		p->GCaL     *= 2.32;
		p->Gto      *= 1.17;
		p->GK1		*= 0.85;
	}
	else if (strcmp(p->Celltype, "LA") == 0)
	{
		p->GKr     				*= 1.6;
		p->Gto     				*= 0.53;
		p->GKs					*= 1.8;
		p->IKs_va_ss_shift 		+= -13.43;
		p->IKs_va_ss_kscale 	*= 0.5;
	}
	else if (strcmp(p->Celltype, "AS") == 0)
	{
		p->GCaL             	*= 0.4;
		p->Gto              	*= 0.4*0.4;
		p->GKur					*= 0.677;
		p->IKur_va_ss_shift 	+= 14;
		p->IKur_vi_ss_shift 	+= -25;
		p->IKur_vi_ss_kscale 	*= 1.769;
	}
	else if (strcmp(p->Celltype, "LAA") == 0)
	{
		p->GCaL             	*= 1.68;
		p->Gto              	*= 0.4;
		p->GKur					*= 0.8;
		p->GKr      			*= 1.6;
		p->Gto      			*= 0.53;
		p->GKs      			*= 1.8;
		p->IKur_va_ss_shift 	+= 14;
		p->IKur_vi_ss_shift 	+= -25;
		p->IKur_va_ss_kscale 	*= 2.13;
		p->IKur_vi_ss_kscale 	*= 1.769;
	}
	else if (strcmp(p->Celltype, "PV") == 0)
	{
		p->Gto 	*= 0.75;
		p->GCaL *= 0.75;
		p->GKr 	*= 2.4;
		p->GKs 	*= 1.87;
		p->GK1 	*= 0.67;
	}
	// End Col 2013 het

	// testing exmaple illustration of implementation
	else if (strcmp(p->Celltype, "test_global") == 0)
	{
		p->Gto              *= 2;       // Scale factor = Multiplies previous settings
		p->gKur             = 0.003;    // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
		p->IKr_va_ss_kscale *= 1.25;    // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->ICaL_vi_ss_shift += 5;       // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->IKs_va_tau_scale *= 0.75;    // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= 2;       // Scales intracellular upatke rate. Multiplies previous settings.
	}
	// Add new celltypes here: else if (strcmp(p->Celltype, "X") == 0) {   }

	else p->Het_set_ref = 0;  // Sets back to zero to indicate celltype has not been set despite function call
}
// End heterogeneity =====================//|

// ISO ===================================\\|
void set_ISO_hAM(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->ISO" is concentration of ISO; "p->ISO_model" is the ISO model to be applied
	// ===========================================
	// This approach assumes linear transition between ISO = 0 and full effect at ISO = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double ISO_scaled = 1/ ( 1 + exp((-4.362*(log(p->ISO) + 3.27)) ) )  and use "ISO_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	p->ISO_set_ref            = 1; // Indicates ISO has been set (is set back to zero if not caught by any below IF statement)

	if (strcmp(p->ISO_model, "Col") == 0)
	{
		// Implementation from Colman et al. 2013 J. Physiol. 591, 4249-4272
		//p->GCaL             *= (1.0 + p->ISO*(2.0	-1));		// x 2 at maximum ISO
		p->GLTCC_kva1_va2   *= (1.0 + p->ISO*(2.0	-1));		// x 2 at maximum ISO -> this is GCaL in native and open rate scale in integrated
		p->GKur             *= (1.0 + p->ISO*(1.6	-1));
		p->GKs              *= (1.0 + p->ISO*(2.5	-1));
		p->Gup              *= (1.0 + p->ISO*(2.5	-1));
		p->GNa              *= (1.0 + p->ISO*(1.25	-1));

		p->ICaL_va_ss_shift += (p->ISO*-5);
	}
	else if (strcmp(p->ISO_model, "GB") == 0)
	{
		// Implementation from Grandi et al. 2011 Circ. Res. 109, 1055-1066.
		p->GLTCC_kva1_va2       *= (1.0 + p->ISO*(1.5-1));
		p->GKur                 *= (1.0 + p->ISO*(3.0-1));
		p->GKs                  *= (1.0 + p->ISO*(3.0-1));

		p->ICaL_va_ss_shift     += (p->ISO*-3);
		p->ICaL_vi_ss_shift     += (p->ISO*-3);
		p->ICaL_va_tau_shift    += (p->ISO*-3);
		p->ICaL_vi_tau_shift    += (p->ISO*-3);
		p->IKs_va_ss_shift      += (p->ISO*-40);
		p->IKs_va_tau_shift     += (p->ISO*-40);
		p->koCa                 += (p->ISO*10); // GB Ca handling only

		p->INaK_kNa         	*= (1.0 + p->ISO*(0.75	-1));
		p->Kmf					*= (1.0 + p->ISO*(0.5 	-1));
		p->koff_tncl			*= (1.0 + p->ISO*(1.5 	-1));
	}
	else p->ISO_set_ref       = 0; // Sets back to zero to indicate ISO has not been set despite function call
}
// End ISO ===============================//|

// Remodelling ===========================\\|
// ALL models
void set_global_remodelling(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->Remodelling_prop" is proportion of maximal Remodelling effect; "p->Remodelling" is the Remodelling to be applied
	// ===========================================
	// This approach assumes linear transition between Remodelling_prop = 0 and full effect at Remodelling_prop = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double Remodelling_scaled = 1/ ( 1 + exp((-4.362*(log(p->Remodelling_prop) + 3.27)) ) )  and use "Remodelling_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	// Start =====
	p->Remodelling_set_ref = 1; // Indicates X(=ISO/Agent/Celltype/Remodeling/Mutation) has been set
	if (strcmp(p->Remodelling, "RSERCA_NCX") == 0)
	{
		// From main paper Colman 2019 PLOS Comp Biol
		p->GK1				*= (1.0 + p->Remodelling_prop*(0.5 - 1));
		p->IK1_va_shift		+= 10*p->Remodelling_prop;
		p->IK1_Erev_shift	+= 10*p->Remodelling_prop; // shift Erev and V function	
		p->Gup				*= (1.0 + p->Remodelling_prop*(1.5  -1));
		p->GCaL				*= (1.0 + p->Remodelling_prop*(0.5  -1));
		p->GNCX				*= (1.0 + p->Remodelling_prop*(0.5  -1));
		p->GKr				*= (1.0 + p->Remodelling_prop*(1.5  -1));
		p->GNa				*= 0.5;

		// Model/cell specific mods
		if (strcmp(p->Model, "minimal") == 0 && strcmp(p->Celltype, "RA") == 0)
		{
			p->GKr          *= (1.0 + p->Remodelling_prop*(0.5  -1));
			p->GK1          *= (1.0 + p->Remodelling_prop*(0.5  -1));
			p->GCaL         *= (1.0 + p->Remodelling_prop*(2	-1));
		}
		else if (strcmp(p->Model, "hVM_ORD_s") == 0) 
		{
			p->GKr			*= (1.0 + p->Remodelling_prop*((1.0/1.5) - 1));
			p->Gup			*= (1.0 + p->Remodelling_prop*(1.25 	 - 1));
		}
	}
	else if (strcmp(p->Remodelling, "RCRU") == 0)
	{
		// From main paper Colman 2019 PLOS Comp Biol
		p->tau_ss_type		= "fast";
	}
	// else if (strcmp(p->X, "SPECIFIER2") == 0) { do stuff }
	// etc
	else p->Remodelling_set_ref = 0; // Indicates X=(Mutation/Agent/Celltype/Remodeling/Mutation) has not actually been set
	// End =======
}

// Human atrial models
void set_remodelling_hAM(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->Remodelling_prop" is proportion of maximal Remodelling effect; "p->Remodelling" is the Remodelling to be applied
	// ===========================================
	// This approach assumes linear transition between Remodelling_prop = 0 and full effect at Remodelling_prop = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double Remodelling_scaled = 1/ ( 1 + exp((-4.362*(log(p->Remodelling_prop) + 3.27)) ) )  and use "Remodelling_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	p->Remodelling_set_ref            = 1; // Indicates Remodelling has been set (is set back to zero if not caught by any below IF statement)

	// AF remodelling models ========\\|
	// Grandi et al. 2011 Circ. Res. 109, 1055-1066
	if (strcmp(p->Remodelling, "AF_GB") == 0)
	{	
		p->GNa					*= (1.0 + p->Remodelling_prop*(0.9  -1));
		p->gNaL					= 0.0025; // (s/mF) NOTE: this is zero if no remodelling, so must be set rather than multiplied!
		p->GKs					*= (1.0 + p->Remodelling_prop*(2.0  -1));
		p->Gto					*= (1.0 + p->Remodelling_prop*(0.3  -1));
		p->GK1					*= (1.0 + p->Remodelling_prop*(2.0  -1));
		p->GCaL					*= (1.0 + p->Remodelling_prop*(0.3  -1));
		p->GNCX					*= (1.0 + p->Remodelling_prop*(1.4  -1)); 
		p->GKur					*= (1.0 + p->Remodelling_prop*(0.5  -1));
		p->Grel					*= (1.0 + p->Remodelling_prop*(3.0  -1));	// (equiv of +20 to baseline of 10 in GB model)
		p->Gleak				*= (1.0 + p->Remodelling_prop*(1.25 -1));
	}
	//Colman et al. 2013 J Physiol 591(17):4249-72 models	
	else if (strcmp(p->Remodelling, "AF_Col_1") == 0) // 
	{
		p->GCaL					*= (1.0 + p->Remodelling_prop*(0.3  -1));
		p->Gto					*= (1.0 + p->Remodelling_prop*(0.3  -1));
		p->GK1  				*= (1.0 + p->Remodelling_prop*(2.0  -1));
		p->Ito_va_ss_shift 		+= 16*p->Remodelling_prop;
		p->ICaL_vi_tau_scale 	*= (1.0 + p->Remodelling_prop*(1.62 -1));	
	}
	else if (strcmp(p->Remodelling, "AF_Col_2") == 0) // 
	{
		p->GCaL					*= (1.0 + p->Remodelling_prop*(0.37  -1));
		p->Gto					*= (1.0 + p->Remodelling_prop*(0.34  -1));
		p->GKur					*= (1.0 + p->Remodelling_prop*(0.51  -1));
		p->GK1					*= (1.0 + p->Remodelling_prop*(2.06  -1));
	}
	else if (strcmp(p->Remodelling, "AF_Col_3") == 0) // 
	{
		p->GCaL                 *= (1.0 + p->Remodelling_prop*(0.35  -1));
		p->Gto                  *= (1.0 + p->Remodelling_prop*(0.35  -1));
		p->GK1                  *= (1.0 + p->Remodelling_prop*(1.75  -1));
	}
	else if (strcmp(p->Remodelling, "AF_Col_4") == 0) // 
	{
		p->GCaL                 *= (1.0 + p->Remodelling_prop*(0.3  -1));
		p->Gto 	                *= (1.0 + p->Remodelling_prop*(0.34 -1));
		p->GKur                 *= (1.0 + p->Remodelling_prop*(0.50 -1));
		p->GK1                  *= (1.0 + p->Remodelling_prop*(2.0  -1));
		p->GKs					*= (1.0 + p->Remodelling_prop*(2.0  -1));
		p->GNCX					*= (1.0 + p->Remodelling_prop*(1.55 -1));
		p->Gleak				*= (1.0 + p->Remodelling_prop*(1.25 -1));
		p->Grel					*= (1.0 + p->Remodelling_prop*(3.0  -1));
	}
	// End AF remodelling models =====//|

	// testing exmaple illustration of implementation
	else if (strcmp(p->Mutation, "test_global") == 0)
	{
		p->Gto              *= (1.0 + p->Remodelling_prop*(2.0  -1));     // Scale factor = Multiplies previous settings  | 2 is scale factor, 2-1 is additive scale factor
		p->IKr_va_ss_kscale *= (1.0 + p->Remodelling_prop*(1.25 -1));     // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->IKs_va_tau_scale *= (1.0 + p->Remodelling_prop*(0.75 -1));     // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= (1.0 + p->Remodelling_prop*(2.0  -1));     // Scales intracellular upatke rate. Multiplies previous settings.
		p->ICaL_vi_ss_shift += 5*p->Remodelling_prop;                     // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->gKur             = 0.003;                                // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods.
	}
	// Add new mutation here: else if (strcmp(p->mutation, "X") == 0) {   }


	else p->Remodelling_set_ref       = 0; // Sets back to zero to indicate Remodelling has not been set despite function call
}
// End Remodelling =======================//|

// Mutations =============================\\|
void set_mutation_hAM(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ================================================================================//|

	p->Mutation_set_ref            = 1; // Indicates Mutation has been set (is set back to zero if not caught by any below IF statement)

	// IKur mutations as in Colman, Ni et al. 2017 PLOS Comp Biol 13(6):e1005587
	if (strcmp(p->Mutation, "D322H") == 0)
	{
		p->GKur					*= 1.78;
		p->IKur_va_ss_shift		+= -3.26;
		p->IKur_va_ss_kscale	*= 0.809;
		p->IKur_vi_ss_shift		+= 9.61;
		p->IKur_vi_ss_kscale	*= 0.801;
	}
	else if (strcmp(p->Mutation, "Y155C") == 0)
	{
		p->GKur					*= 0.4745;
		p->IKur_va_ss_shift		+= 0.89;
		p->IKur_va_ss_kscale	*= 0.75;
		p->IKur_vi_ss_shift		+= 5.01;
		p->IKur_vi_ss_kscale	*= 0.78;
	}
	// End IKur mutations

	// testing exmaple illustration of implementation
	else if (strcmp(p->Mutation, "test") == 0)
	{
		p->Gto              *= 2;       // Scale factor = Multiplies previous settings
		p->gKur             = 0.003;    // Actual conductance explicitly set. Will overwrite any previous settings of g, but not scale factor mods. 
		p->IKr_va_ss_kscale *= 1.25;    // Multiplies gradient parameter for voltage-activation steady-state. Multiplies previous settings
		p->ICaL_vi_ss_shift += 5;       // Shifts the voltage dependence of the steady state of inactivation gate. Summed to previous settings
		p->IKs_va_tau_scale *= 0.75;    // Multiplies time constant of voltage activation. Multiplies previous settings.
		p->Gup              *= 2;       // Scales intracellular upatke rate. Multiplies previous settings.
	}
	// Add new mutation here: else if (strcmp(p->mutation, "X") == 0) {   }

	else p->Mutation_set_ref       = 0;    // Sets back to zero to indicate Mutation has not been set despite function call
}
// End Mutations =========================//|

// ACh ===================================\\|
void set_ACh_global(Cell_parameters *p)
{
	// ================================================================================\\|
	// Can modify parameters (e.g. conductance) directly, or scaling factors.
	// (conductance denoted "g", scale factor "G"; current = f(g*G))
	// In general, use the scale factors as these are multiplicative/additive throughout
	// the code. 
	// Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
	// ===========================================
	// "p->ACh" is concentration of ACh; "p->ACh_model" is the ACh model to be applied
	// ===========================================
	// This approach assumes linear transition between ACh = 0 and full effect at ACh = 1. 
	// If non-linear effect is reqired, multiply by another function:
	// e.g. double ACh_scaled = 1/ ( 1 + exp((-4.362*(log(p->ACh) + 3.27)) ) )  and use "ACh_scaled" in all below functions
	// This could be done per target for different dose dependency
	// ================================================================================//|

	// NOTE: if the model does not have IKACh (or If) in, then modifying them won't make a difference!

	p->ACh_set_ref = 1; // Indicates X(=ISO/Agent/Celltype/Remodeling/Mutation) has been set
	if (strcmp(p->ACh_model, "test_global") == 0)  // just an example of how to implement: NOT an actual implementation
	{
		p->GCaL		*= (1.0 + (0.5 - 1)*p->ACh);
		p->GK1		*= (1.0 + (2.0 - 1)*p->ACh);
		p->gKACh	= p->gKACh_max * ( pow(p->ACh, 1.5) / ( pow(2.8e-1, 1.5) + pow(p->ACh, 1.5) ) ); 
		p->If_va_ss_shift += -7.2*p->ACh;
	}
	else p->ACh_set_ref = 0; // not been set
}
// End ACh ===============================//|

// Spatial gradient ======================\\|
void set_spatial_gradient(Cell_parameters *p)
{
    // ================================================================================\\|
    // Can modify parameters (e.g. conductance) directly, or scaling factors.
    // (conductance denoted "g", scale factor "G"; current = f(g*G))
    // In general, use the scale factors as these are multiplicative/additive throughout
    // the code.
    // Only use the parameter itself when wanting to set its baseline value (which will still be modded by the modifier)
    // ================================================================================//|

    // NOTE: these provided examples are not based on any data and are purely to
    // illustrate implementation

    // in single cells, it is "prop" - value of prop will be determined in tissue by
    // a map, so "prop" and "map" are used interchangably in these comments
    // Implement as: factor = value_at_map=0 * (1-map) + value_at_map=1 * map
    // If either end of your map corresponds to the baseline cell model, simply
    // set these values as 1. 

    if (strcmp(p->spatial_gradient, "none") == 0);    // Do nothing
    else if (strcmp(p->spatial_gradient, "apico_basal_example") == 0)
    {
		// Just example of implementation - not real!
        //         value at apex (map = 0)            value at base (map = 1)
        p->GKr      *= 0.5*(1-p->spatial_gradient_prop)   + 1.5*p->spatial_gradient_prop;
        p->GCaL     *= 1.25*(1-p->spatial_gradient_prop)  + 0.75*p->spatial_gradient_prop;
    }
    else if (strcmp(p->spatial_gradient, "SAN_distance_example") == 0)
    {
		// Just example of implementation - not real!
        //         value adjacent to SAN = baseline (map = 0)   value away from SAN (PV, LAA etc) (map = 1)
        p->GKr      *= 1*(1-p->spatial_gradient_prop)   + 1.5*p->spatial_gradient_prop;
        p->GCaL     *= 1*(1-p->spatial_gradient_prop)   + 0.75*p->spatial_gradient_prop;
    }
	else 
	{
		printf("ERROR: \"%s\" is not a valid spatial_gradient Please check Model.c for options\n\n", p->spatial_gradient);
		exit(1);
	}
}
// End Spatial gradient ==================//|
// End Global or common functions =====================================================//|
// End Current modification variables | Het and modulation ======================================//|

// Reveral potentials ===========================================================================\\|
void compute_reversal_potentials(Cell_parameters p, Model_variables *var, State_variables *s)
{
	var->ENa			= 		((p.R * p.T)/p.F)*log(s->Nao/s->Nai);
	var->EK				= 		((p.R * p.T)/p.F)*log(s->Ko/s->Ki);
	var->EKs            =       ((p.R * p.T)/p.F)*log((s->Ko + 0.01833*s->Nao)/(s->Ki + 0.01833*s->Nai));
	var->EKs_ORD        =       ((p.R * p.T)/p.F)*log((s->Ko +         s->Nao)/(s->Ki + 0.01833*s->Nai));
	var->ECa			= 0.5*	((p.R * p.T)/p.F)*log(s->Cao/s->Cai);
	var->ECl			= 		((p.R * p.T)/p.F)*log(15.0/150.0);

	var->ENa_j			= 		((p.R * p.T)/p.F)*log(s->Nao/s->Nai_j);
	var->ENa_sl			= 		((p.R * p.T)/p.F)*log(s->Nao/s->Nai_sl);
	var->ECa_j			= 0.5*	((p.R * p.T)/p.F)*log(s->Cao/s->Cai_j);
	var->ECa_sl			= 0.5*	((p.R * p.T)/p.F)*log(s->Cao/s->Cai_sl);
}
// Compute Reveral potentials ===================================================================//|

// Excitation properties / measurements =========================================================\\|
void determine_excitation_state(Model_variables *var, double Vm, double time)
{
	if (var->ex_switch == 0)	// If currently not excited
	{
		if (Vm > -30)			// threshold to determine excited state || set defaults
		{
			var->ex_switch		= 1; 	// In excitation state
			var->APD_t_switch	= 0; 	// Ready to be calculated
			var->t_ex			= time;

			var->dvdt_max_prev	= var->dvdt_max;
			var->dvdt_max		= 0;

			var->Vmax_prev		= var->Vmax;
			var->Vmax			= -80; 	// Ensure below threshold
			var->Vmin_prev_prev = var->Vmin_prev;
			var->Vmin_prev		= var->Vmin;
			var->Vmin			= 50;
			var->Vamp_prev		= var->Vamp;

			for (int i = 0; i < 9; i++) var->APD_p_switch[i] = 0;

			var->CaT_min_prev	= var->CaT_min;
			var->CaT_min		= 1;	// Ensure above values
			var->CaT_max_prev   = var->CaT_max;
			var->CaT_max		= 0;	// Ensure below values
			var->CaSR_min_prev	= var->CaSR_min;	
			var->CaSR_min		= 5000;	// Ensure above values
			var->CaSR_max_prev	= var->CaSR_max;
			var->CaSR_max		= 0;	// Ensure below values

			// Set previous APD
			var->APD_t_prev		= var->APD_t;
			for (int i = 0; i < 9; i++) var->APD_p_prev[i] = var->APD_p[i];

            var->J_SERCA_integral   = 0;
            var->J_rel_integral     = 0;
		}
	}
	else if (var->ex_switch == 1)    // If currently excited
	{
		if (Vm < -45)			// threshold to determine repolarised
		{
			var->ex_switch 		= 0; // no longer in excitation state
		}
	}
}

void determine_excitation_state_integrated_0D(Model_variables *var, double Vm, double time, double *Ca_JSR_t_ex, double Ca_JSR, double *dyad_SRF_prop_active, double srf_SRF_prop_active, int *srf_init, int *srf_set, const char *SRF_Mode)
{   
	if (var->ex_switch == 0)    // If currently not excited
	{
		if (Vm > -30)           // threshold to determine excited state || set defaults
		{
			var->ex_switch      = 1;    // In excitation state
			var->APD_t_switch   = 0;    // Ready to be calculated
			var->t_ex           = time;

			var->dvdt_max_prev  = var->dvdt_max;
			var->dvdt_max       = 0;

			var->Vmax_prev      = var->Vmax;
			var->Vmax           = -80;  // Ensure below threshold
			var->Vmin_prev_prev = var->Vmin_prev;
			var->Vmin_prev      = var->Vmin;
			var->Vmin           = 50;
			var->Vamp_prev      = var->Vamp;

			for (int i = 0; i < 9; i++) var->APD_p_switch[i] = 0;

			var->CaT_min_prev   = var->CaT_min;
			var->CaT_min        = 1;    // Ensure above values
			var->CaT_max_prev   = var->CaT_max;
			var->CaT_max        = 0;    // Ensure below values
			var->CaSR_min_prev  = var->CaSR_min;    
			var->CaSR_min       = 5000; // Ensure above values
			var->CaSR_max_prev  = var->CaSR_max;
			var->CaSR_max       = 0;    // Ensure below values

			// Set previous APD 
			var->APD_t_prev     = var->APD_t;
			for (int i = 0; i < 9; i++) var->APD_p_prev[i] = var->APD_p[i];

			// Set CaJSR at time of exciataion
			*Ca_JSR_t_ex    = Ca_JSR;

			// SRF stuff
			*dyad_SRF_prop_active	= srf_SRF_prop_active;	
			if (strcmp(SRF_Mode, "Direct_Control") == 0)
			{
				*srf_init      = 0; // comment out if want just one static beat
				*srf_set       = 0;
			}
			else if (strcmp(SRF_Mode, "Dynamic") == 0)
			{
				*srf_set    = -1;
			}	
		}
	}
	else if (var->ex_switch == 1)    // If currently excited
	{
		if (Vm < -65)           // threshold to determine repolarised
		{
			var->ex_switch      	= 0; // no longer in excitation state
			*dyad_SRF_prop_active	= 0;		
		}
	}
}

// Calculate integrals
void calculate_flux_integrals(Cell_parameters p, Model_variables *var, double J_SERCA, double J_NCX, double J_rel, double J_LTCC)
{
    var->J_SERCA_integral   += p.dt*J_SERCA;
    var->J_NCX_integral     += p.dt*J_NCX;
    var->J_LTCC_integral    += p.dt*J_LTCC;
    var->J_rel_integral     += p.dt*J_rel;
}

void calculate_measurement_properties(Model_variables *var, double Vm1, double Vm2, double time, double dt, double APD_threshold, double CaT, double CaSR)
{
	// dv/dt and dv/dt_max
	var->dvdt		= (Vm2 - Vm1)/dt;
	if (var->dvdt > var->dvdt_max) var->dvdt_max = var->dvdt;

	// Max, min and amplitude
	if (Vm2 > var->Vmax)	var->Vmax = Vm2;
	if (var->dvdt <= 0.0 )if (Vm2 < var->Vmin)	var->Vmin = Vm2;  // calculate Vmin if AP is repolarising only
	var->Vamp				= var->Vmax - var->Vmin_prev; // such that amplitude is V just before stimulus to peak, not on repolarisation

	// Maxmimum and minimum Ca2+ properties
	if (CaT > var->CaT_max)		var->CaT_max = CaT;
	if (CaT < var->CaT_min)		var->CaT_min = CaT; 
	if (CaSR > var->CaSR_max)	var->CaSR_max = CaSR;
	if (CaSR < var->CaSR_min)	var->CaSR_min = CaSR; 

	// APD to set threshold
	if (var->APD_t_switch == 0)	// if in state where APD needs to be calculated (i.e. has been excited but not calculated)
	{
		if (Vm2 < APD_threshold) // if the voltage is now below the threshold
		{
			var->APD_t			= time - var->t_ex;
			var->APD_t_switch	= 1;	// Been calculated
		}
	}

	// APD to different percentages
	double perc, threshold;
	for (int i = 0; i < 9; i++)
	{
		if (var->APD_p_switch[i] == 0)	// if in state where APD needs to be calculated (i.e. has been excited but not calculated)
		{
			perc = (i+1)*10; 	// converts 0 to 10%, 8 to 90%
			perc *= 0.01;		// convers % to proportion
			threshold = var->Vmax - perc*(var->Vmax - var->Vmin_prev); // Vmin_prev is V just before excitation; current Vmax
			if (Vm2 < threshold) // if the voltage is now below the threshold
			{
				var->APD_p[i]			= time - var->t_ex;
				var->APD_p_switch[i]	= 1;	// Been calculated
			}
		}
	}
}


// End Excitation properties / measurements =====================================================//|

// Voltage clamp ================================================================================\\|
void run_voltage_clamp(Cell_parameters p, Model_variables *var, State_variables *s, char const *directory, double dt)
{
	double Vm, Vclamp, Vhold, time, clamp_time, Vstart, Vend;
	double Ipeak, Ipeak2, Ipeak3, Ipeak4;
	char * filename       = (char*)malloc(500);

	FILE * I_out;
	FILE * I_out_individual;
	FILE * IV_out;

	// ICaL V clamp ===========================\\|
	printf("Applying ICaL voltage clamp\n");
	sprintf(filename, "%s/ICaL_Vclamp_traces.txt", directory);
	I_out = fopen(filename, "wt");
	sprintf(filename, "%s/ICaL_IV.txt", directory);
	IV_out = fopen(filename, "wt");

	// Vclamp settings
	Vstart 		= -50;
	Vend		= 50;
	Vhold		= -80;//-50;
	clamp_time 	= 500;

	for (Vclamp = -50; Vclamp < Vend+1; Vclamp +=5)
	{
		Ipeak = 0;
		initial_conditions_native(s, p, p.Model);    // lib/Model.c
		if (int(Vclamp)%10 == 0)
		{	
			sprintf(filename, "%s/ICaL_Vclamp_trace_%.0f.txt", directory, Vclamp);
			I_out_individual = fopen(filename, "wt");
		}
		for (time = 0.0; time < 3*clamp_time; time += dt)
		{
			if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
			else Vm = Vclamp;

			compute_model_native(p, var, s, Vm, dt);

			if (Vm == Vclamp && var->ICaL < Ipeak) Ipeak = var->ICaL;

			//fprintf(I_out, "%f %f %f %f\n", time, Vm, var->ICaL, s->Cai);
            fprintf(I_out, "%f %f %f %f %f %f %f ", time, Vm, var->ICaL, s->Cai, s->ICaL_va, s->ICaL_vi, s->ICaL_ci);
            fprintf(I_out, "%f %f %f %f %f\n", s->ICaL_5sm_O, s->ICaL_5sm_C1, s->ICaL_5sm_C2, s->ICaL_5sm_I1, s->ICaL_5sm_I2);
			if (int(Vclamp)%10 == 0) 
            {
                fprintf(I_out_individual, "%f %f %f %f %f %f %f ", time, Vm, var->ICaL, s->Cai, s->ICaL_va, s->ICaL_vi, s->ICaL_ci);
                fprintf(I_out_individual, "%f %f %f %f %f\n", s->ICaL_5sm_O, s->ICaL_5sm_C1, s->ICaL_5sm_C2, s->ICaL_5sm_I1, s->ICaL_5sm_I2);
            }

		}
		fprintf(IV_out, "%f %f\n", Vclamp, Ipeak);
		if (int(Vclamp)%10 == 0) fclose(I_out_individual);
	}

	fclose(I_out);
	fclose(IV_out);
	// End ICaL V clamp =======================//|

	// Ito/IKur/IKr/IKs V clamp ===============\\|
	printf("Applying IK voltage clamp\n");
	sprintf(filename, "%s/IK_Vclamp_traces.txt", directory);
	I_out = fopen(filename, "wt");
	sprintf(filename, "%s/IK_IV.txt", directory);
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
		sprintf(filename, "%s/IK_Vclamp_trace_%.0f.txt", directory, Vclamp);
		I_out_individual = fopen(filename, "wt");
		for (time = 0.0; time < 3*clamp_time; time += dt)
		{
			if (time < clamp_time || time > 2*clamp_time) Vm = Vhold;
			else Vm = Vclamp;

			compute_model_native(p, var, s, Vm, dt);

			if (Vm == Vclamp && var->Ito  > Ipeak)  Ipeak  = var->Ito;
			if (Vm == Vclamp && var->IKur > Ipeak2) Ipeak2 = var->IKur;
			if (Vm == Vclamp && var->IKr  > Ipeak3) Ipeak3 = var->IKr;
			if (Vm == Vclamp && var->IKs > Ipeak4)  Ipeak4 = var->IKs;

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
	sprintf(filename, "%s/IK1_IV.txt", directory);
	IV_out = fopen(filename, "wt");
	for (Vclamp = -150; Vclamp < 100; Vclamp +=0.1)
	{
		Vm = Vclamp;
		compute_model_native(p, var, s, Vm, dt);
		fprintf(IV_out, "%f %f\n", Vm, var->IK1);
	}
	fclose(IV_out);
	// End IK1 ================================//|

	free(filename);
}
// End Voltage clamp ============================================================================//|

// Frequently used functions ====================================================================\\|
double rush_larsen(double y, double ss, double tau, double dt)
{
	double gate;
	gate = ss - (ss-y)*exp(-dt/tau);
	return gate;
}

double sigmoid(double V, double V_half, double k)
{
	double value;
	value = 1/(1 + exp((V - V_half)/k) );
	return value;
}
// End Frequently used functions ================================================================//|

// Luo-Rudy 1991 INa formulation ================================================================\\|
// Reference: ********************************************||
// "A model of the ventricular cardiac action potential. 
// Depolarization, repolarization, and their interaction"
// Circulation Research. 1991;68:1501-1526
// *******************************************************||
void set_INa_LR_rates(Cell_parameters p, Model_variables *var, double Vm)
{
	double Vm_ac        = Vm - p.INa_va_shift; 	// Shift of the voltage used to calculate alpha and beta, activation
	double Vm_inac      = Vm - p.INa_vi_shift;	// Shift of the voltage used to calculate alpha and beta, inactivation

	// Set Activation gate alpha and beta
	var->INa_va_al                 = 0.32*(Vm_ac+47.13)/(1-exp(-0.1*(Vm_ac+47.13)));
	if (fabs(Vm_ac + 47.13) < 1e-10) var->INa_va_al = 3.2;
	var->INa_va_bet                = 0.08*exp(-Vm_ac/11.0);

	// Set inactivation gates alphas and betas
	if (Vm_inac < -40.0)
	{
		var->INa_vi_1_al           = 0.135*exp((80+Vm_inac)/-6.8);
		var->INa_vi_1_bet          = 3.56*exp(0.079*Vm_inac)+310000*exp(0.35*Vm_inac);
		var->INa_vi_2_al           = (-127140*exp(0.2444*Vm_inac)-0.00003474*exp(-0.04391*Vm_inac))*((Vm_inac+37.78)/(1+exp(0.311*(Vm_inac+79.23))));
		var->INa_vi_2_bet          = (0.1212*exp(-0.01052*Vm_inac))/(1+exp(-0.1378*(Vm_inac+40.14)));
	}
	else
	{
		var->INa_vi_1_al           = 0;
		var->INa_vi_1_bet          = 1.0/(0.13*(1+exp((Vm_inac+10.66)/-11.1)));
		var->INa_vi_2_al           = 0;
		var->INa_vi_2_bet          = (0.3*exp(-0.0000002535*Vm_inac))/(1+exp(-0.1*(Vm_inac+32)));
	}

	// Set tau and SS from alpha and beta
	var->INa_va_tau                = 1.0/(var->INa_va_al + var->INa_va_bet); // 1/(a+b)
	var->INa_vi_1_tau              = 1.0/(var->INa_vi_1_al + var->INa_vi_1_bet);
	var->INa_vi_2_tau              = 1.0/(var->INa_vi_2_al + var->INa_vi_2_bet);
	var->INa_va_ss                 = var->INa_va_al * var->INa_va_tau; // a*tau
	var->INa_vi_1_ss               = var->INa_vi_1_al * var->INa_vi_1_tau;
	var->INa_vi_2_ss               = var->INa_vi_2_al * var->INa_vi_2_tau;

	var->INa_va_tau 				*= p.INa_va_tau_scale;
	var->INa_vi_1_tau 				*= p.INa_vi_1_tau_scale;
	var->INa_vi_2_tau 				*= p.INa_vi_2_tau_scale;
}

void update_gates_INa_LR(Cell_parameters p, Model_variables *var, State_variables *s, double Vm, double dt)
{
	s->INa_va                      = rush_larsen(s->INa_va, var->INa_va_ss, var->INa_va_tau, dt); // lib/Membrane.c
	s->INa_vi_1                    = rush_larsen(s->INa_vi_1, var->INa_vi_1_ss, var->INa_vi_1_tau, dt);
	s->INa_vi_2                    = rush_larsen(s->INa_vi_2, var->INa_vi_2_ss, var->INa_vi_2_tau, dt);
}

void compute_INa_LR(Cell_parameters p, Model_variables *var, State_variables *s, double Vm)
{
	var->INa    	= p.gNa * pow(s->INa_va, 3) * s->INa_vi_1 * s->INa_vi_2 * (Vm - var->ENa);
	var->INa		*= p.GNa;
}
// End Luo-Rudy 1991 INa formulation ============================================================//|
