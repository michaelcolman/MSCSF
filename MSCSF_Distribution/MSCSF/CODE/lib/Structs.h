// ****REVIEW/BETA version of the code. NOT FINAL. NOT TO BE REDISTRIBUTED*****
// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Declaration of all variable structs =========  //
// ========================================================  //
// COPYRIGHT MICHAEL A. COLMAN 2015-2019.==================  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND =  //
// CONTRIBUTORS "AS-IS" AND ANY EXPRESS OR IMPLIED ========  //
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED =  //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A ========  //
// PARTICULAR PURPOSE ARE DISCLAIMED. =====================  //
// ========================================================  //
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

#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdbool.h>
#include "MersenneTwister.h"

// Struct list:
// struct{}Smulation_parameters;
// struct{}Cell_parameters;
// struct{}State_variables;
// struct{}Model_variables;
// struct{}Ca_variables;
// struct{}CRU_variables;
// struct{}Dyad_variables;
// struct{}SR_fluxes;
// struct{}Membrane_fluxes;
// struct{}RAND;
// struct{}SC_variables;
// struct{}Spontaneous_release_functions;
// struct{}Tissue_parameters;
// struct{}Argument_parameters;

// Define the simulation parameters struct ======================================================\\|
typedef struct{

	// Reference
	char const	*reference;			// Any reference to identify simulation
	char const  *results_reference; // To name a specific results directory within the Outputs directory
	char const  *state_reference_read;   // To name state files for particular cases within conditions
	char const  *state_reference_write;   // To name state files for particular cases within conditions

	// Simulation pacing and time
	int BCL;				// ms
	int Total_time;			// ms
	int Paced_time;			// ms
	int NBeats;				// N
	double dt;				// ms

	int S2_CL;				// ms
	int NS2;				// N
	int S2_time;			// ms 

	// Simulation conditions
	char const	*Vclamp;		// "On" or "Off"
	char const	*Write_state;	// "On" or "Off"
	char const 	*Read_state; 	// "On" or "Off"

	// Interval to output spatial files
	int Spatial_output_interval_vtk;	// ms // for VTK
	int Spatial_output_interval_data;	// ms // for array
    int Spatial_output_start_time;      // lower bound of time to output spatial data
    int Spatial_output_end_time;        // upper bound of time to output spatial data

	// Delayed impose CaSR functionality
	const char *Delayed_CaSR_IC; 	// "On" or "Off"
	double		CaSR_IC_delay;		// ms
	bool		CaSR_set;			// true or false if already been set

}Simulation_parameters;
// End Define the simulation parameters struct ==================================================//|

// Define the Cell_parameters struct (set once) =================================================\\|
// Contains model parameters and constants and scaling/shift variables (as not dynamically determined)
typedef struct{

	// Global control variables ===================================\\|
	double 		dt;						// Integration time-step	
	char const* Model;					// The baseline model
	char const* Celltype;				// Region or other celltype
	char const* Agent;					// Pharmacological agent
	double		Agent_prop;				// Proportion of pharma agent to set (linear 0-1)
	char const* Remodelling;			// Disease remodelling
	double 		Remodelling_prop;		// Proportion of maximum remodelling to set (linear 0-1)
	char const* ISO_model;				// Defines which ISO model to use (if multiple models are available for a specific cell model)
	double 		ISO;					// Isoprenaline concentration
	char const* ACh_model;				// ACh model to be implemented
	double		ACh;					// Acetlycholine concentration
	char const* spatial_gradient;       // Model for a spatial gradient, such as apico-basal or distance from SAN
	double      spatial_gradient_prop;  // Proportion of spatial gradient
	char const* Mutation;				// Mutation
	int 		Het_set_ref;			// Reference for tracking if celltype has been set
	int			ISO_set_ref;			// Reference for tracking if ISO has been set
	int			Agent_set_ref;			// Reference for tracking if Agent has been set
	int			Remodelling_set_ref;	// Reference for tracking if Remodelling has been set
	int 		Mutation_set_ref;		// Reference for tracking if Mutation has been set
	int 		ACh_set_ref;			// Reference for tracking if ACh has been set

	// Global control, spatial cell model only
	char const*	tau_ss_type;			// Identifier for time constant of coupling of sub-space

	// Multi-model, single species trackers
	bool	hAM;				// True if a hAM model; false for all others

	// hAM specific settings
	char const* Ca_handling; 	// Ca handling model, hAM_WL (chose between GB and CRN)
	char const*	environment;	// isolated vs intact
	// End global control variables ===============================//|	

	// Constants ==================================================\\|
	double R;					// J mol^-1 K^-1
	double F;					// C mmol-1
	double T;					// K
	double FoRT;	
	// End Constants ==============================================//|

	// Cell structure==============================================\\|
	double Cm;					// Membrane capacitance 						(pF)
	double Cm_F;				// Membrane capacitance							(F)
	double Cm_CRU;				// Membrane capacitance per CRU					(pF)

	double Vcell;				// volume of whole cell (um^3)
	double Vcyto;				// volume of intracellular Ca2+ space
	double VjSR;				// volume of jSR / release Ca2+ space
	double VnSR;				// volume of nSR / uptake Ca2+ space
	double Vjunc;				// Volume of intracellular junctional release space
	double Vsl;					// Volume of intracellular sub-sarcolemmal space
	double Vc;					// Volume of extracellular cleft space

	double Fjunc;				// Proportion of currents in junctional vs non-junctional compartments
	double Fjunc_ICaL;			// Proportion of currents in junctional vs non-junctional compartments, LTCCs/ICaL

	// Spatial Ca handling cell-structure parameters
	double 	vds_CRU_mean;		// Volume of single dyadic cleft space 		(micro m^3)
	double 	vss_CRU;			// Volume of single intracellular subspace 	(micro m^3)
	double 	vcyt_CRU;			// Volume of single bulk intracellular space(micro m^3)
	double 	vnsr_CRU;			// Volume of single network SR 				(micro m^3)
	double 	vjsr_CRU;			// Volume of single junctional SR			(micro m^3)
	double 	v_ss_cyt;			// Ratio of ss to cyto volumes
	double 	v_cyt_nsr;			// Ratio of cyto to nsr volumes
	double 	v_jsr_nsr;			// Ratio of jsr to nsr volumes
	int 	NRyR_mean;			// mean number of RyRs per dyad
	int 	NLTCC_mean;			// mean number of LTCCs per dyad
	double 	NRyR_propvar;		// maximum differences for NRyR as proportion of NRyR
	double 	NLTCC_propvar;		// maximum difference for NLTCC as proportion of NLTCC
	// End cell structure  ========================================//|

	// Stimulus parameters ========================================\\|
	double stimduration;		// duration of applied stimulus (ms)
	double stimmag;				// magnitude of applied stimulus (A/F || pA/pF)
	// End Stimulus parameters ====================================//|

	// Current parameters =========================================\\|
	// Hybrid minimal model ===================\\|
	double gIp0d;				// Conductance of phase 0 depolarising current (s/mF | ns/pF)
	double gIp1r;				// Conductance of phase 1 repolarising current (s/mF | ns/pF)
	double gIp2d;				// Conductance of phase 2 depolarising current (s/mF | ns/pF)
	double gIp2r;				// Conductance of phase 2 repolarising current (s/mF | ns/pF)
	double gIp3r;				// Conductance of phase 3 repolarising current (s/mF | ns/pF)
	double gIp4r;				// Conductance of phase 4 repolarising current (s/mF | ns/pF) 
	// End Hybrid minimal model ===============//|

	// Biophysically detailed models ==========\\|
	// Conductances
	double gNa;					// Conductance of fast-sodium current 			(s/mF)
	double gNaL;				// Conductance of fast-sodium current 			(s/mF)
	double gto;					// Conductance of transient-outward K+ current 	(s/mF)
	double gCaL;				// Conductance of L-type ca current 			(s/mF)
	double pCaL;				// Permeability constant fot L-type ca current	(cm/s)
	double pCaL_K;				// Permeability constant fot ICaL, K+ component (cm/s)
	double pCaL_Na;				// Permeability constant fot ICaL, Na+ component(cm/s)
	double gKur;				// Conductance of ultrarpid K+ current 			(s/mF)
	double gKr;					// Conductance of rapid delayed K+ current 		(s/mF)
	double gKs;					// Conductance of slow delayed K+ current 		(s/mF)
	double gK1;					// Conductance of time-independent K+ current 	(s/mF)
	double gClCa;				// Conductance of Ca2+-activated chloride		(s/mF)
	double gNab;				// Conductance of background sodium current 	(s/mF)
	double gCab;				// Conductance of background calcium current 	(s/mF)
	double gKb;					// Conductance of background potassium current 	(s/mF)
	double gKACh;				// Conductance of ACh-activated K+ current 		(s/mF)
	double gKACh_max;			// Conductance of ACh-activated K+ current 		(s/mF) (maximal; conc independent)
	double gf;					// Conductance of "funny" current				(s/mF)
	double gClb;                // Conductance of background cloride current    (s/mF)

	double AIhyp;				// Magnitude of applied hyperpolarising current (pA/pF)

	// ICaL
	double ICaL_ci_tau;			// Calcium induced activation, time constant	(ms)
	double ICaL_ci_al;			// Calcium induced activation, alpha rate 		(ms^-1)
	double ICaL_ci_bet;			// Calcium induced activation, beta rate		(ms^-1)
    double ICaL_5sm_Or;         // constant opening rate, 5 state LTCC markov chain model
    double ICaL_5sm_Cr;         // constant opening rate, 5 state LTCC markov chain model

	// NCX params
	double INCX_bar;			// Maximal current scale factor for NCX 		(pA/pF) or (um^3.uM.ms^-1)
	double INCX_kNao;			// [Na+]o saturation constant for INCX 			(mM)
	double INCX_kNai;			// [Na+]i saturation constant for INCX          (mM)
	double INCX_kCao;			// [Ca2+]o saturation constant for INCX 		(mM)
	double INCX_kCai;			// [Ca2+]i saturation constant for INCX 		(mM)	or (uM)
	double INCX_k;				// Saturation factor for INCX 					(dimensionless)
	double INCX_gamma;			// Voltage dependence factor for INCX 			(dimensionless)
	double INCX_kda;			// Magnitude calcium saturation constant 		(mM)	or (uM)
	double INCX_ksat;			// Saturation constant							(mM)

	// NaK params
	double INaK_bar;			// Maximal current scale factor for INaK		(pA/pF)
	double INaK_kK;				// Half-Saturation constant for Ko INaK 		(mM)
	double INaK_kNa;			// Half-saturation constant for Nai INaK		(mM)

	// Membrane calcium pump (PMCA) params
	double ICaP_bar;			// Conductance factor for calcium pump (pA/pF)	or (um^3 . uM . ms^-1)
	double ICaP_kCa;			// Saturation constant for calcium pump (mM)	or (uM)

	// Background calcium current, flux version
	double ICab_bar;			// (um^3.uM.ms^-1)  CHECK

	//IClCa params
	double IClCa_kd;			// Saturation constant for Ca2+-activated Cl current 

	// ICaL kinetics
	double ICaL_vi_Fs;			// Fraction of slow voltage inactivation
	// End Biophysically detailed model =======//|
	// End Current parameters =====================================//|

	// Concentrations =============================================\\|
	// Many may also be updated state variables; the param in that case
	// becomes the IC for the state variable
	double Nai;					// Intracellular sodium concentration (mM)
	double Nao;					// Extracellular sodium concentration (mM)
	double Ki;					// Intracellular potassium concentration (mM)
	double Ko;					// Extracellular potassium concentration (mM)
	double Cai;					// Intracellular calcium concentration (mM) (or uM for spatial cell models)
	double Cao;					// Extracellular calcium concentration (mM)
	double CaSR;				// SR calcium concenration			   (mM) (or uM for spatial cell models)

	double Nai_sl;					// Intracellular sodium concentration (mM)
	double Nai_j;					// Intracellular sodium concentration (mM)
	double Cai_sl;					// Intracellular calcium concentration (mM) (or uM for spatial cell models)
	double Cai_j;					// Intracellular calcium concentration (mM) (or uM for spatial cell models)
	// End concentrations =========================================//|

	// Ca2+ handling ==============================================\\|
	// General/common =========================\\|
	double J_rel_max;			// Maximal flux rate of intracellular Ca2+ release 	(mM/ms) or (um^3/ms)
	double J_SERCA_max;			// Maximal flux rate of intracellular Ca2+ uptake 	(mM/ms) or (uM/ms)
	double J_SERCA_kCa;			// SERCA Cai constant								(mM)	or (uM)
	double J_SERCA_kCaSR;		// SERCA CaSR constant								(uM)
	double J_leak_max;			// Maximal flux rate of intracellular Ca2+ leak		(mM/ms) or ms^-1
	double J_leak_kCaSR;		// Jleak CaSR2+ constant							(-)		or uM
	double J_jsr_nsr_tau;		// Time constant of transfer between jsr and nsr	(ms)
	// End General/common =====================//|

	// Spatial Ca handling specific ===========\\|
	// RyR Ca-dependent
	double RyR_kCO_A;			// RyR closed->open rate multiplier		(uM^H/ms)
	double RyR_kOC;				// RyR open->closed rate				(ms^-1)
	double GRyR_kCO;			// Scales open transition rate
	double RyR_Cads_H;			// Power to raise Ca_ds by in CO
	double RyR_kOC_det;			// multiplier for deterministic model

	// RyR CSQN dependent
	double RyR_mon_tau;			// (ms)
	double RyR_mi_tau;			// (ms)
	double RyR_mon_beta_tau;	// (ms)
	double RyR_mi_beta_tau;		// (ms)
	double RyR_mon_grad;		

	// LTCC
	double J_LTCC_max;			// Maximal flux rate (NOT conductance!)	(micro mol C-1 ms^-1 * (micro m^3))
	double LTCC_kva1_va2;		// transition rate to activated state (voltage-independent) 	(ms^-1)
	double LTCC_kva2_va1;		// transition rate to not activated state (voltage-independent)	(ms^-1)
	double LTCC_Ca_bar;			// Calcium colnstant for Ca-induced inactivation				(uM)
	double GLTCC_kva1_va2;		// Scales tansition rate

	// inter-compartment transfer
	double tau_ds;				// dyadic cleft transfer time constant	(ms)
	double tau_ss_cyt;			// transfer between ss and cyto			(ms)
	double tau_nsr_jsr;			// transfer between nsr and jsr			(ms)

	// Spatial coupling transfer
	double tau_ss_trans_f;		// transverse time constant of ss-ss(fast) (ms)
	double tau_ss_long_f;		// transverse time constant of ss-ss(fast) (ms)
	double tau_ss_trans_mf;		// transverse time constant of ss-ss(mfast)(ms)
	double tau_ss_long_mf;		// transverse time constant of ss-ss(mfast)(ms)
	double tau_ss_trans_m;		// transverse time constant of ss-ss(med) (ms)
	double tau_ss_long_m;		// transverse time constant of ss-ss(med) (ms)
	double tau_ss_trans_ms;		// transverse time constant of ss-ss(med) (ms)
	double tau_ss_long_ms;		// transverse time constant of ss-ss(med) (ms)
	double tau_ss_trans_s;		// transverse time constant of ss-ss(mslow)(ms)
	double tau_ss_long_s;		// transverse time constant of ss-ss(mslow)(ms)
	double tau_ss_trans;		// transverse time constant of ss-ss	(ms)
	double tau_ss_long;			// longitudinal time cnstant ss-ss		(ms)
	double tau_cyto_trans;		// transverse time constant of cyt-cyt	(ms)
	double tau_cyto_long;		// longitudinal time cnstant cyt-cyt	(ms)
	double tau_nsr_trans;		// transverse time constant of nsr-nsr	(ms)
	double tau_nsr_long;		// longitudinal time cnstant nsr-nsr	(ms)

	// Buffering (similar to GB lab, but not identical)
	double Kcam;				// uM
	double Bcam;				// uM
	double Kbsr;				// uM
	double Bbsr;				// uM
	double Kmca;				// uM
	double Bmca;				// uM
	double Kmmg;				// uM
	double Bmmg;				// uM

	double Kcsqn;				// mM
	double Bcsqn;				// mM
	// End Spatial Ca handling specific =======//|

	// Nattel-lab parameters ==================\\|
	double cmdnbar;				// Total calmodulin concentration
	double trpnbar;				// Total troponin concentration
	double csqnbar;				// Total csqn concentration
	double cmdn_k;				// half-saturation constant for calmodulin 	(mM)
	double trpn_k;				// half-saturation constant for troponin 	(mM)
	double csqn_k;				// half-saturation constant for calsequestrin (mM)
	// End Nattel-lab parameters ===============//|

	// Grandi-Bers lab parameters ==============\\|
	// Cell structure, sophisticated and environment
	double junctionLength;		// um
	double junctionRadius;		// um
	double distSLcyto;			// um
	double distJuncSL;			// um

	double DcaJuncSL;			// cm^2/s
	double DcaSLcyto;			// cm^2/s
	double DnaJuncSL;			// cm^2/s
	double DnaSLcyto;			// cm^2/s

	double SAjunc;				// um^2
	double SAsl;				// um^2

	double J_ca_juncsl;			// L/ms
	double J_ca_slmyo;			// L/ms
	double J_na_juncsl;			// L/ms
	double J_na_slmyo;			// L/ms

	// Transport - SR (
	double ks;					// ms^-1
	double Kmf;					// mM
	double Kmr;					// mM
	double hillSRCaP;
	double koCa;				// ms^-1 mM^-1
	double kiCa;				// ms^-1 mM^-1
	double kim;					// ms^-1
	double kom;					// ms^-1
	double ec50SR;				// mM
	double MaxSR;				// mM (?)
	double MinSR;				// mM (?)
	// J_SERCA_max mM/ms already defined for all models above

	// Buffering
	double Bmax_SLlowj;			// mM
	double koff_sll;			// ms^-1
	double kon_sll;				// ms^-1 mM^-1
	double Bmax_SLhighj;		// mM
	double kon_slh;				// ms^-1 mM^-1
	double koff_slh;			// ms^-1
	double Bmax_SLlowsl;		// mM
	double Bmax_SLhighsl;		// mM
	double koff_myoca;			// ms^-1
	double kon_myoca;			// ms^-1 mM^-1
	double Bmax_myosin;			// mM
	double Mgi;
	double kon_myomg;			// ms^-1 mM^-1
	double koff_myomg;			// ms^-1
	double kon_tnchca;			// ms^-1 mM^-1
	double koff_tnchca;			// ms^-1
	double Bmax_TnChigh;		// mM
	double kon_tnchmg;			// ms^-1 mM^-1
	double koff_tnchmg;			// ms^-1
	double kon_tncl;			// ms^-1 mM^-1
	double koff_tncl;			// ms^-1
	double Bmax_TnClow;			// mM
	double Bmax_CaM;			// mM
	double kon_cam;				// ms^-1 mM^-1
	double koff_cam;			// ms^-1
	double Bmax_SR;				// mM
	double kon_sr;				// ms^-1 mM^-1
	double koff_sr;				// ms^-1
	double Bmax_Csqn;			// mM
	double kon_csqn;			// ms^-1 mM^-1
	double koff_csqn;			// ms^-1
	double Bmax_Naj;			// mM
	double Bmax_Nasl;			// mM
	double koff_na;				// ms^-1
	double kon_na;				// ms^-1 mM^-1
	// End Grandi-Bers lab parameters ==========//|

	// colman 2013 / reduced Koivumaki =========\\|
	// Need to tidy with unit etc!!!!
	double BCa;
	double SLlow;
	double SLhigh;

	double KdBCa;
	double KdSLlow;
	double KdSLhigh;

	double DCa;
	double DCaSR;
	double Aj_nj;
	double xj_nj;
	// End colman 2013 / reduced Koivumaki =====//|
	// End Ca2+ handling ==========================================//|	

	// MODIFIERS ==================================================\\|
	// Simple scaling 
	double GNa;					// Scale factor for fast sodium current (and Ip0d)
	double GNaL;				// Scale factor for late sodium current
	double Gto;					// Scale factor for Ito (and Ip1r)
	double GCaL;				// Scale factor for ICaL (and Ip2d)
	double GKur;				// Scale factor for IKur (and Ip2r)
	double GKr;					// Scale factor for IKr (and Ip3r)
	double GKs;					// Scale factor for IKs
	double GK1;					// Scale factor for IK1
	double GNCX;				// Scale factor for INCX
	double GCaP;				// Scale factor for ICaP
	double GNab;				// Scale factor for INab
	double GCab;				// Scale factor for ICab
	double GKb;					// Scale factor for IKb
	double GNaK;				// Scale factor for INaK
	double GClCa;				// Scale factor for IClCa
	double GClb;				// Scale factor for IClb
	double GKACh;				// Scale factor for IKACh
	double Gf;					// Scale factor for If

	// Time constant scaling
	double INa_va_tau_scale;	// Scales time constant for INa voltage activation
	double INa_vi_1_tau_scale;	// Scales time constant for INa voltage inactivation 1
	double INa_vi_2_tau_scale;	// Scales time constant for INa voltage inactivation 2
	double INaL_va_tau_scale;	// Scales time constant for INaL voltage activation
	double INaL_vi_tau_scale;	// Scales time constant for INaL voltage inactivation
	double Ito_va_tau_scale;	// Scales time constant for Ito voltage activation
	double Ito_vi_tau_scale;	// Scales time constant for Ito voltage inactivation
	double ICaL_va_tau_scale;	// Scales time constant for ICaL voltage activation
	double ICaL_vi_tau_scale;	// Scales time constant for ICaL voltage inactivation
	double IKur_va_tau_scale;	// Scales time constant for IKur voltage activation
	double IKur_vi_tau_scale;	// Scales time constant for IKur voltage inactivation
	double IKr_va_tau_scale;	// Scales time constant for IKr voltage activation
	double IKr_vi_tau_scale;	// Scales time constant for IKr voltage inactivation
	double IKs_va_tau_scale;	// Scales time constant for IKr voltage activation
	double IKACh_va_tau_scale;	// Scales time constant for IKACh voltage activation
	double IKACh_vi_tau_scale;	// Scales time constant for IKACh voltage inactivation
	double If_va_tau_scale;		// Scales time constant for If volatge activation

	// Voltage dependence
	// (note: approach is to shift the V1/2 in the equations; equivalent is to minus shift from Voltage input into equations)
	double INa_va_shift; 		// Shift of the V of activation || alpha and beta  	mV
	double INa_vi_shift;		// Shift of the V of inactivation || alpha and beta mV

	double INaL_va_shift; 		// Shift of the V of activation || alpha and beta  	mV
	double INaL_vi_shift;		// Shift of the V of inactivation || alpha and beta mV
	double INaL_va_ss_shift;	// Steady state ac
	double INaL_vi_ss_shift;	// Steady state inac
	double INaL_va_tau_shift;	// time constant ac
	double INaL_vi_tau_shift;	// time constant inac
	double INaL_va_ss_kscale;	// gradient of ac ss 
	double INaL_vi_ss_kscale;	// gradient of inac ss

	double Ito_va_ss_shift;		// Shift of the V1/2 of acitvation steady state     mV
	double Ito_vi_ss_shift;		// Shift of the V1/2 of inactivation steady state	mV
	double Ito_va_tau_shift;	// Shift of the voltage dependence of time constant mV
	double Ito_vi_tau_shift;	// Shift of the voltage dependence of time constant mV
	double Ito_va_ss_kscale;	// Scales the gradient parameter of activation steady state
	double Ito_vi_ss_kscale;	// Scales the gradient parameter of inactivation steady state

	double ICaL_va_ss_shift;    // Shift of the V1/2 of acitvation steady state     mV
	double ICaL_vi_ss_shift;    // Shift of the V1/2 of inactivation steady state   mV
	double ICaL_va_tau_shift;   // Shift of the voltage dependence of time constant mV
	double ICaL_vi_tau_shift;   // Shift of the voltage dependence of time constant mV
	double ICaL_va_ss_kscale;   // Scales the gradient parameter of activation steady state
	double ICaL_vi_ss_kscale;   // Scales the gradient parameter of inactivation steady state

	double IKur_va_ss_shift;    // Shift of the V1/2 of acitvation steady state     mV
	double IKur_vi_ss_shift;    // Shift of the V1/2 of inactivation steady state   mV
	double IKur_va_tau_shift;   // Shift of the voltage dependence of time constant mV
	double IKur_vi_tau_shift;   // Shift of the voltage dependence of time constant mV
	double IKur_va_ss_kscale;   // Scales the gradient parameter of activation steady state
	double IKur_vi_ss_kscale;   // Scales the gradient parameter of inactivation steady state

	double IKr_va_ss_shift;    	// Shift of the V1/2 of acitvation steady state     mV
	double IKr_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double IKr_va_ss_kscale;    // Scales the gradient parameter of activation steady state
	double IKr_vi_ss_shift;    	// Shift of the V1/2 of acitvation steady state     mV
	double IKr_vi_ss_kscale;    // Scales the gradient parameter of activation steady state
	double IKr_vi_tau_shift;    // Shift of the voltage dependence of time constant mV

	double IKACh_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double IKACh_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double IKACh_va_ss_kscale;    // Scales the gradient parameter of activation steady state
	double IKACh_vi_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double IKACh_vi_ss_kscale;    // Scales the gradient parameter of activation steady state
	double IKACh_vi_tau_shift;    // Shift of the voltage dependence of time constant mV

	double IKs_va_ss_shift;    	// Shift of the V1/2 of acitvation steady state     mV
	double IKs_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double IKs_va_ss_kscale;    // Scales the gradient parameter of activation steady state

	double If_va_ss_shift;    	// Shift of the V1/2 of acitvation steady state     mV
	double If_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double If_va_ss_kscale;    // Scales the gradient parameter of activation steady state

	double IK1_va_shift;    	// Shift of the V1/2 of acitvation steady state     mV

	// Other variables which essentially change form of the equation
	double IKr_rect_vshift;
	double IKr_rect_grad_scale;
	double IKr_rect_exp_factor;
	double IK1_VEexp_scale;
	double IK1_den_add_factor;
	double IK1_Erev_shift;

	// Ca handling modification
	double Gup;					// Scale factor for intracellular, SERCA Ca2+ uptake
	double Gleak;				// Scale factor for intracellular Ca2+ leak
	double Grel;				// Scale factor for intracellualr Ca2+ release 

	double BC;					// bulk const from colman 2013 hAM model
	// End MODIFIERS ==============================================//|
}Cell_parameters;
// End Define the parameters struct (set once) ==================================================//|

// Define the State variables (time-dependent) ==================================================\\|
// Gating variables, concentrations, and othe time-dependent variables
typedef struct{

	double Vm;					// Whole-cell membrane potential (mV)

	// Current gating variables ===================================\\|
	// Biophysically detailed currents ========\\|
	// INa
	double INa_va;				// voltage activation
	double INa_vi_1;			// voltage inactivation
	double INa_vi_2;			// voltage inactivation 2

	//INaL
	double INaL_va;				// voltage activation
	double INaL_vi;				// voltage inactivation

	// Ito
	double Ito_va;				// voltage activation
	double Ito_vi;				// voltage inactivation (fast if slow also used)
	double Ito_vi_s;			// voltage inactivation, slow
	double Ito_vi_3;			// third Ito voltage inactivation gate

	// ICaL
	double ICaL_va;				// voltage activation
	double ICaL_vi;				// voltage inactivation (fast if slow also used)
	double ICaL_vi_s;			// voltage inactivation, slow
	double ICaL_ci;				// calcium inactivation
	double ICaL_ci_j;			// calcium inactivation, junctional

	// IKur
	double IKur_va;				// voltage activation
	double IKur_vi;				// voltage inactivation

	// IKr
	double IKr_va;				// voltage activation
	double IKr_vi;				// used where IKr is used for IKf in Rabbit atria || also used for va_slow

	// IKs
	double IKs_va;				// voltage activation
	double IKs_va_2;			// voltage activation 2

	// IK1 (ORD has time-dependent component)
	double IK1_va;

	// IKACh
	double IKACh_va;			// voltage activation
	double IKACh_vi;			// voltage inactivation

	// If
	double If_va;				// voltage activation
	// End Biophysically detailed currents ====//|

	// Hybrid minimal model currents ==========\\|
	// Phase 0 depolarising current
	double Ip0d_va;				// voltage activation
	double Ip0d_vi_1;			// voltage inactivation 1
	double Ip0d_vi_2;			// voltage inactivation 2

	// Phase 1 repolarising current
	double Ip1r_va;				// voltage activation
	double Ip1r_vi;				// voltage inactivation

	// Phase 2 depolarising current
	double Ip2d_va;				// voltage activation
	double Ip2d_vi;				// voltage inactivation

	// Phase 2 repolarising current
	double Ip2r_va;				// voltage activation
	double Ip2r_vi;				// voltage activation

	// Phase 3 repolarising current
	double Ip3r_va;				// voltage activation
	// End Hybrid minimal model currents ======//|
    
    // ICaL 5-state Markov model ==============\\|
    double ICaL_5sm_C1; // closed 1
    double ICaL_5sm_C2; // closed 2 
    double ICaL_5sm_I1; // inactivated 1
    double ICaL_5sm_I2; // inactivated 2
    double ICaL_5sm_O;  // open
    // End ICaL 5-state Markov model ==========//|
	// End Current gating variables ===============================//|

	// Concentrations =============================================\\|
	// If not updated, then state variable will be set once from the
	// param. Note that the models use the state variable, so it must
	// exist even if it is a constant
	double Nai;                 // Intracellular sodium concentration (mM)
	double Nao;                 // Extracellular sodium concentration (mM)
	double Ki;                  // Intracellular potassium concentration (mM)
	double Ko;                  // Extracellular potassium concentration (mM)
	double Cai;                 // Intracellular calcium concentration (mM)
	double Cao;                 // Extracellular calcium concentration (mM)

	double Nai_j;				// Intracellular sodium concentration, junctional compartment (mM)
	double Nai_sl;				// Intracellular sodium concentration, junctional compartment (mM)
	double Cai_j;				// Intracellular calcium concentration, junctional compart 	(mM)
	double Cai_sl;				// Intracellular calcium concentration, sub-sarcolemmal (mM)
	// End concentrations =========================================//| 

	// Ca2+ handling ==============================================\\|
	double cmdn;				// calmodulin concentration 	(mM)
	double trpn;				// troponin concentration		(mM)
	double csqn;				// calsequestrin concenrration 	(mM)

	double CajSR;				// junctional SR / release compartment Ca concentration (mM)
	double CanSR;				// network SR / uptake compartment Ca concentration  	(mM)

	double RyRo;				// proportion open RyRs / activation gate
	double RyRr;				// proportion refactory RyRs / inactivation gate 1
	double RyRi;				// proportion inactivated RyRs / inactvation gate 2

	// Grandi-Bers Ca handling
	double	Tn_CHm;
	double	Tn_CHc;
	double	Myo_m;
	double	Myo_c;
	double	Tn_CL;

	// Nygren-Giles Ca handling
	double CaCalse;
	double CaCal;
	double Catrop;
	double Camg;
	double Mgmg;
	// End Ca2+ handling ==========================================//|

}State_variables;
// End Define the state variables (time-dependent) ==============================================//|

// Define the model variables struct  ===========================================================\\|
// All variables which are calculated within the model, and may need to be seen from multiple functions
// Includes the currents themselves, the current variables (steady states, taus), stimulation variables
typedef struct{

	// Stimulus variables =========================================\\|
	double 	Istim;				// Stimulus current (pA/pF)
	bool	stimflag;			// true IF stimulus is being applied
	int 	BCL_int;			// integer value of BCL	
	int 	stimduration_int;	// integer value of stim duration
	int 	stimcount;			// number of times stimulus has been applied
	int 	stimcount_S2;		// number of times stimulus has been applied
	int 	S2_int;				// integer value of S2 coupling interval
	double	Istim_S2;			// S2 stimulus current (pA/pF)
	bool	S2_stimflag;		// true IF S2 stimulus is being applied
	int 	Paced_time_int;		// integer value of paced time (s1)
	// end Stimulus variables =====================================//|

	// Reversal potentials ========================================\\|
	double ENa;					// reversal potential, sodium ions 				(mV)
	double EK;					// reversal potential, potassium ions 			(mV)
	double EKs;					// reversal potential, potassium ions, IKs		(mV)
	double EKs_ORD;				// reversal potential, potassium ions, IKs, ORD	(mV)
	double ECa;					// reversal potential, calcium ions 			(mV)
	double ENa_j;				// reversal potential, junctional sodium 		(mV)
	double ENa_sl;				// reversal potential, sub-sarcolemmal sodium	(mV)
	double ECa_j;				// reversal potential, junctional calcium 		(mV)
	double ECa_sl;				// reversal potential, sub-sarcolemmal calcium 	(mV)
	double ECl;					// reversal potential, chloride 				(mV)
	// End Reversal potentials ====================================//|

	// Currents ===================================================\\|
	double Itot;				// Total ionic current (pA/pF)

	// Real currents (all pA/pF) 
	double INa;					// Fast sodium current
	double INaL;				// Late sodium current
	double Ito;					// Transient outward potassium current
	double ICaL;				// L-type calcium current
	double IKur;				// Ultrarpid potassium current
	double IKr;					// Rapid delayed rectifier potassium current
	double IKs;					// Slow delayed rectifier potassium current
	double IK1;					// Time-independent rectifier potassium current
	double INCX;				// Sodium-calcium exchanger
	double INaK;				// Sodim-potassium pump
	double ICaP;				// membrane calcium pump (PMCA)
	double INab;				// Background sodium current
	double ICab;				// Background calcium current
	double IKb;					// Background potassium current
	double IClCa;				// Ca2+-activated chloride current
	double IClb;				// Background chloride current
	double IKACh;				// ACh-activated K+ current
	double If;					// "Funny" current
	//double Ihyp;				// Magnitude of applied hyperpolarising current

	// Currents, more sophisticated structure models
	double INa_sl;				// Fast sodium current, sub-sarcolemmal
	double INaL_sl;				// Late sodium current, sub-sarcolemmal
	double INab_sl;				// Background sodium current, sub-sarcolemmal
	double ICab_sl;				// Background calcium current, sub-sarcolemmal
	double ICaP_sl;				// membrane calcium pump (PMCA), sub-sarcolemmal
	double INCX_sl;				// Sodium-calcium exchanger, sub-sarcolemmal
	double ICaL_sl;				// L-type calcium current, sub-sarcolemmal
	double INaK_sl;				// Sodim-potassium pump, sub-sarcolemmal
	double IClCa_sl;			// Ca2+-activated chloride current, sub-sarcolemmal
	double IKs_sl;				// Slow delayed rectifier potassium current, sub-sarcolemmal
	double INa_j;				// Fast sodium current, junctional
	double INaL_j;				// Late sodium current, junctional
	double INab_j;				// Background sodium current, junctional
	double ICab_j;				// Background calcium current, junctional
	double ICaP_j;				// membrane calcium pump (PMCA), junctional
	double INCX_j;				// Sodium-calcium exchanger, junctional
	double ICaL_j;				// L-type calcium current, junctional
	double INaK_j;				// Sodim-potassium pump, junctional
	double IClCa_j;				// Ca2+-activated chloride current, junctional
	double IKs_j;				// Slow delayed rectifier potassium current, junctional

	double ICaL_Ca_j;			// Calcium specific component, ICaL, junctional
	double ICaL_Ca_sl;			// Calcium specific component, ICaL, sub-sarcolemmal
	double ICaL_K;				// Potassium component, ICaL
	double ICaL_Na_j;			// Sodium specific component, ICaL, junctional
	double ICaL_Na_sl;			// Sodium specific component, ICaL, sub-sarcolemmal

	// Hybrid minimal model current names
	double Ip0d;				// Phase-0 depolarising current
	double Ip1r;				// Phase-1 repolarising current
	double Ip2d;				// Phase-2 depolarising current
	double Ip2r;				// Phase-2 repolarising current
	double Ip3r;				// Phase-3 repolarising current
	double Ip4r;				// Phase-4 repolarising current
	// End currents ===============================================//|

	// Current calculation variables ==============================\\|
	// Dynamic conductance
	double ICaL_bar_sl;			// Maximal current through the L-type calcium channels, sub-sarcolemma 	(Ca-specific component if appropriate)
	double ICaL_bar_j;			// Maximal current through the L-type calcium channels, junctional		(Ca-specific component if appropriate)
	double ICaL_bar_K;			// Maximal current through the L-type calcium channels, Potassium component
	double ICaL_bar_Na_sl;		// Maximal current through the L-type calcium channels, Sodium, sub-sarcolemma
	double ICaL_bar_Na_j;		// Maximal current through the L-type calcium channels, Sodium, junctional

	// Steady states and time constants -  biophysically detailed models
	double INa_va_ss;          	// voltage activation, steady state
	double INa_va_tau;         	// voltage activation, time constant
	double INa_vi_1_ss;        	// voltage inactivation, steady state
	double INa_vi_1_tau;       	// voltage inactivation, time constant
	double INa_vi_2_ss;        	// voltage inactivation, steady state
	double INa_vi_2_tau;       	// voltage inactivation, time constant
	double INa_va_al;          	// voltage activation, alpha transition rate (1-y-> y)
	double INa_va_bet;         	// voltage activation, beta transition rate  (y -> 1-y)
	double INa_vi_1_al;        	// voltage activation, alpha transition rate (1-y -> y)
	double INa_vi_1_bet;       	// voltage activation, beta transition rate  (y -> 1-y)
	double INa_vi_2_al;        	// voltage activation, alpha transition rate (1-y -> y)
	double INa_vi_2_bet;       	// voltage activation, beta transition rate  (y -> 1-y)

	double INaL_va_ss;         	// voltage activation, steady state
	double INaL_va_tau;        	// voltage activation, time constant
	double INaL_vi_ss;         	// voltage inactivation, steady state
	double INaL_vi_tau;        	// voltage inactivation, time constant
	double INaL_va_al;          // voltage activation, alpha transition rate (1-y-> y)
	double INaL_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
	double INaL_vi_al;          // voltage inactivation, alpha transition rate (1-y-> y)
	double INaL_vi_bet;         // voltage inactivation, beta transition rate  (y -> 1-y)

	double Ito_va_ss;          	// voltage activation, steady state
	double Ito_va_tau;         	// voltage activation, time constant
	double Ito_vi_ss;          	// voltage inactivation, steady state
	double Ito_vi_3_ss;        	// voltage inactivation, steady state
	double Ito_vi_tau;         	// voltage inactivation, time constant
	double Ito_vi_s_tau;       	// voltage inactivation, time constant, slow
	double Ito_vi_3_tau;       	// voltage inactivation, time constant, slow
	double Ito_va_al;			// voltage activation, alpha transition rate (1-y-> y)
	double Ito_va_bet;			// voltage activation, beta transition rate  (y -> 1-y)
	double Ito_vi_al;			// voltage activation, alpha transition rate (1-y-> y)
	double Ito_vi_bet;			// voltage activation, beta transition rate  (y -> 1-y)
	double Ito_vi_Fs;			// Fraction of slow voltage inactivation

	double ICaL_va_ss;          // voltage activation, steady state
	double ICaL_va_tau;         // voltage activation, time constant
	double ICaL_vi_ss;          // voltage inactivation, steady state
	double ICaL_vi_tau;         // voltage inactivation, time constant
	double ICaL_vi_s_tau;       // voltage inactivation, time constant, slow
	double ICaL_ci_ss;			// calcium inactivation, steady state
	double ICaL_ci_j_ss;		// calcium inactivation, steady state, junctional
	double ICaL_ci_tau;			// calcium inactivation, time constant		
	double ICaL_ci_al;          // Calcium induced activation, alpha rate       (ms^-1)
	double ICaL_ci_bet;         // Calcium induced activation, beta rate        (ms^-1)	

	double IKur_va_ss;          // voltage activation, steady state
	double IKur_va_tau;         // voltage activation, time constant
	double IKur_vi_ss;          // voltage inactivation, steady state
	double IKur_vi_tau;         // voltage inactivation, time constant
	double IKur_va_al;          // voltage activation, alpha transition rate (1-y-> y)
	double IKur_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
	double IKur_vi_al;          // voltage activation, alpha transition rate (1-y-> y)
	double IKur_vi_bet;         // voltage activation, beta transition rate  (y -> 1-y)
	double IKur_dynamic_g;		// Dynamic conductance factor (s/mF)

	double IKr_va_ss;			// voltage activation, steady state
	double IKr_va_tau;			// voltage activation, time constant
	double IKr_va_al;			// alpha transition rate
	double IKr_va_bet;			// beta transition rate

	double IKr_vi_ss;			// voltage inactivation, steady state
	double IKr_vi_tau;			// voltage inactivation, time constant
	double IKr_vi_al;			// alpha transition rate
	double IKr_vi_bet;			// beta transition rate
	double IKr_va_Fs;			// Fraction of gates if using dynamic and static

	double IKs_va_ss;			// voltage activation, steady state
	double IKs_va_tau;			// voltage activation, time constant
	double IKs_va_2_tau;		// voltage activation 2, time constant
	double IKs_va_al;			// voltage activation, alpha transition rate
	double IKs_va_bet;			// voltage activation, beta transition rate

	double IKr_vi_ti;			// Time-independent inactivation gate
	double IK1_va_ti;			// Time-independent activation
	double IK1_va_ss;			// Steady state for time-dependent component
	double IK1_va_tau;			// time constant for time-dependent component

	double IKACh_va_ss;			// voltage activation, steady state
	double IKACh_va_tau;		// voltage activation, time constant
	double IKACh_vi_ss;			// voltage activation, steady state
	double IKACh_vi_tau;		// voltage activation, time constant
	double IKACh_v_ti;			// time-independent factor

	double If_va_ss;			// voltage activation, steady state
	double If_va_tau;			// voltage activation, time constant

	// Steady states and time constants - hybrid minimal model
	double Ip0d_va_ss;			// voltage activation, steady state
	double Ip0d_va_tau;			// voltage activation, time constant
	double Ip0d_vi_1_ss;		// voltage inactivation, steady state
	double Ip0d_vi_1_tau;		// voltage inactivation, time constant
	double Ip0d_vi_2_ss;		// voltage inactivation, steady state
	double Ip0d_vi_2_tau;		// voltage inactivation, time constant
	double Ip0d_va_al;			// voltage activation, alpha transition rate (1-y-> y)
	double Ip0d_va_bet;			// voltage activation, beta transition rate  (y -> 1-y)
	double Ip0d_vi_1_al;		// voltage activation, alpha transition rate (1-y -> y)
	double Ip0d_vi_1_bet;		// voltage activation, beta transition rate  (y -> 1-y)
	double Ip0d_vi_2_al;		// voltage activation, alpha transition rate (1-y -> y)
	double Ip0d_vi_2_bet;		// voltage activation, beta transition rate  (y -> 1-y)

	double Ip1r_va_ss;			// voltage activation, steady state
	double Ip1r_va_tau;			// voltage activation, time constant
	double Ip1r_vi_ss;			// voltage inactivation, steady state
	double Ip1r_vi_tau;			// voltage inactivation, time constant

	double Ip2d_va_ss;			// voltage activation, steady state
	double Ip2d_va_tau;			// voltage activation, time constant
	double Ip2d_vi_ss;			// voltage inactivation, steady state
	double Ip2d_vi_tau;			// voltage inactivation, time constant

	double Ip2r_va_ss;			// voltage activation, steady state
	double Ip2r_va_tau;			// voltage activation, time constant
	double Ip2r_vi_ss;			// voltage inactivation, steady state
	double Ip2r_vi_tau;			// voltage inactivation, time constant

	double Ip3r_va_ss;			// voltage activation, steady state
	double Ip3r_va_tau;			// voltage activation, time constant
	double Ip3r_va_al;			// voltage activation, alpha transition rate (1-y -> y)
	double Ip3r_va_bet;			// voltage activation, beta transition rate  (y -> 1-y)
	double Ip3r_vi_ti;			// Time-independent inactivation gate

	double Ip4r_va_ti;			// Time-independent activation gate

    // 5 state ICaL markov chain rates
    double ICaL_5sm_rate_C1_C2;
    double ICaL_5sm_rate_C2_C1;
    double ICaL_5sm_rate_I1_I2;
    double ICaL_5sm_rate_I2_I1;
    double ICaL_5sm_rate_C1_I1;
    double ICaL_5sm_rate_I1_C1;
    double ICaL_5sm_rate_C2_I2;
    double ICaL_5sm_rate_I2_C2;
    double ICaL_5sm_rate_C2_O;
    double ICaL_5sm_rate_O_C2;
    double ICaL_5sm_rate_I2_O;
    double ICaL_5sm_rate_O_I2;
    double ICaL_O;      // open, all models
	// End current calculation variables ==========================//|

	// Homeostasis ================================================\\|
	// Concentrations
	double dNai;				// Differential of intracellular sodium concentration 		(mM/ms)
	double dKi;					// Differential of intracellular potassium concentration 	(mM/ms)
	double dCai;				// Differential of intracellular calcium concentration 		(mM/ms)

	// Ca2+ handling
	double J_rel;				// Intracellular calcium release flux 	(mM/ms)
	double J_SERCA;				// Intracellular Ca2+ uptake flux		(mM/ms)
	double J_leak;				// Intracellular Ca2+ leak flux			(mM/ms)
	double J_jsr_nsr;			// jSR - nSR transfer flux				(mM/ms)
	double Cai_Fn;				// intracellular Ca2+ flux term			(mM/ms)

	double dMyo_c;	// Buffering differentials
	double dMyo_m;
	double dTn_CHc;
	double dTn_CHm;
	double dTn_CL;
	double J_CaB_cytosol;
	double kCaSR;
	double koSRCa;
	double kiSRCa;
	double RI;
	// End Homeostasis ============================================//|

	// Properties and measurement variables =======================\\|
	int 	ex_switch;				// Tracks whether cell is currently excited
	int 	APD_t_switch;			// Tracks whether APD at threshold has been calculated
	int		APD_p_switch[9];		// Tracks whether APD to % repolarisation has been calculated
	double 	t_ex;					// Time at which cell was excited	(ms)
	double 	dvdt;					// Rate of change of voltage		(mV/ms)
	double 	dvdt_max;				// Maximum rate of change of voltage(mV/ms)
	double 	dvdt_max_prev;			// Maximum rate of change of voltage(mV/ms)
	double	Vmax;					// Maximum voltage					(mV)
	double 	Vmax_prev;				// Maximum voltage, previous beat	(mV)
	double  Vmin;					// Minimum voltage					(mV)
	double 	Vmin_prev;				// Minimum voltage, previous beat 	(mV) *well, actually at point of excitation for current beat
	double 	Vmin_prev_prev;			// Minimum voltage, previous beat 	(mV) *point of excitation for previous beat
	double	Vamp;					// Amplitude voltage				(mV)
	double	Vamp_prev;				// Amplitude voltage, previous beat (mV)
	double 	APD_t;					// APD at threshold voltage			(ms)
	double 	APD_p[9];				// APD at multiple % repolarisation (ms)
	double 	APD_t_prev;				// previous beat APD at threshold 	(ms)
	double 	APD_p_prev[9];			// previous beat APD at % 			(ms)
	double 	CaT_min;				// Minimum intracellular Ca2+		(mM)
	double 	CaT_max;				// Maximum intracellular Ca2+		(mM)
	double 	CaSR_min;				// Minimum SR Ca2+					(mM)
	double 	CaSR_max;				// Maximum SR Ca2+					(mM)
	double 	CaT_min_prev;			// Minimum intracellular Ca2+ prev	(mM)
	double 	CaT_max_prev;			// Maximum intracellular Ca2+ prev	(mM)
	double 	CaSR_min_prev;			// Minimum SR Ca2+ prev				(mM)
	double 	CaSR_max_prev;			// Maximum SR Ca2+ prev				(mM)

    // Flux integrals
    double  J_SERCA_integral;
    double  J_NCX_integral;
    double  J_LTCC_integral;
    double  J_rel_integral;
	// End Properties and measurement variables ===================//|

}Model_variables;
// End Define the model variables struct ========================================================//|

// Structs used in novel (3D + integrated) Ca handling - unused in "native" models ========================\\|
// Define the Ca struct (3D and integrated cell models ==========================================\\|
typedef struct{

	// Convention: CAPS = whole cell average; lowercase = local, spatial
	// Whole cell averages	|| uM
	double DS; 		// concentration in the dyadic cleft
    double SS; 		// subspace
    double CYTO; 	// bulk-intracellular space
    double JSR; 	// Junctional SR 
    double NSR; 	// Network SR

	// Local concentrations
	double *ds;
    double *ss;
    double *cyto;
	double *jsr;
	double *nsr;

	// Reaction term variables, global and local
	double SS_reac;
	double CYTO_reac;
	double JSR_reac;
	double NSR_reac;
	double *ss_reac;
	double *cyto_reac;
	double *jsr_reac;
	double *nsr_reac;

	// Buffering
	double Bcyto;       // Buffering in the cytoplasm
	double Bss;         // Buffering in subspace
	double Bjsr;        // Buffering in JSR

	// Local buffering
	double *bcyto;
	double *bss;
	double *bjsr;

} Ca_variables;
// End Define the Ca struct (3D and integrated cell models ======================================//|

// CRU variables (whole-cell settings and averages) =============================================\\|
typedef struct{

	char const *Cell_size;		// Actual size of cell
	char const *Sim_Cell_size;	// "full", "portion" + other specifics > proportion of cell to be simulated

	int NX, NY, NZ;			// NCRU(x,y,z) dimensions - modelled
	int NX2, NY2, NZ2;		// NCRU(x,y,z) dimensions - full
	int NTOT_CRUs;			// For scaling to whole cell when solving cell portion
    int NCRUs_surface;      // Number of CRUs assocaited with surface sarcolemma
    int NCRUs_interior;     // Number of CRUs assocaited with cell interior
    int NCRUs_MEM;          // Number of CRUs which have a membrane component (ICaL, JMEM)
    int NCRUs_interior_MEM; // Number of interior CRUs with membrane compoenent
    double NLTCC_redist_factor;  // Factor to redistribute LTCCs due to detubulation

	double J_REL;			// Whole-cell average of J_Rel
	double PRyR_OA;			// Proportion of RyRs in O1; whole cell average
	double PRyR_OI;			// Proportion of RyRs in O1; whole cell average
	double PRyR_CA;			// Proportion of RyRs in O1; whole cell average
	double PRyR_CI;			// Proportion of RyRs in O1; whole cell average
	double Monomer;			// Whole-cell average state of monomer
	double Mi;				// Whole-cell average state of Mi
	int	   Nactive; 		// Total number of CRUs in active state
	double Pactive;			// Proportion of CRUs in active state

	double J_CAL;			// Whole-cell average of J_CaL
	double I_CAL;			// Whole-cell average of I_CaL
	double PLTCC_O;			// Proportion of open LTCCs; whole cell average
	double PLTCC_va1;		// Proportion of LTCCs in votage activation, but not active, state; whole cell average
	double PLTCC_va2;		// Proportion of LTCCs in voltage activation AND active state; whole cell average
	double PLTCC_vi;		// Proportion of LTCCs in voltage inactivation state; whole cell average
	double PLTCC_ci;		// Proportion of LTCCs in calcium inactivation state; whole cell average

	double J_SERCA;			// Whole-cell average of J_SERCA
	double J_LEAK;			// Whole-cell average of J_leak
	double J_NCX_bulk;		// Whole-cell average of J_NCX
	double J_CaP_bulk;		// Whole-cell average of J_CaP
	double J_Cab_bulk;		// Whole-cell average of J_Cab
	double J_NCX_ss;		// Whole-cell average of J_NCX
	double J_CaP_ss;		// Whole-cell average of J_CaP
	double J_Cab_ss;		// Whole-cell average of J_Cab

	double I_NCX_bulk;		// Whole cell current (pA/pF)
	double I_CaP_bulk;		// Whole cell current (pA/pF)
	double I_Cab_bulk;		// Whole cell current (pA/pF)
	double I_NCX_ss;		// Whole cell current (pA/pF)
	double I_CaP_ss;		// Whole cell current (pA/pF)
	double I_Cab_ss;		// Whole cell current (pA/pF)
	
	// sub-cellular maps =============\\|
	const char	*Detub;				// "On" or "Off" to read in a detubulation map
    const char  *LTCC_redist;       // "On" or "Off" to redistribute LTCCs when detub is on
	char const 	*TT_map_file;		// Filename for reading in a TT map
	int 		*TT_map;			// TT_map
	char const 	*SERCA_map_file;	// Filename for reading in a SERCA map
	const char	*SERCA_het;			// "On" or "Off" to read in a SERCA map
	double 		*SERCA_map;			// SERCA_map
	char const 	*NCX_map_file;	    // Filename for reading in a NCX map
	const char	*NCX_het;			// "On" or "Off" to read in a NCX map
	double 		*NCX_map;			// NCX_map
	char const 	*RyR_het_map_file;	// Filename for reading in a RyR map
	const char	*RyR_het;			// "Off", "random" or "map"
	double 		*RyR_het_map;		// RyR_map
	char const 	*LTCC_map_file;	    // Filename for reading in a LTCC map
	const char	*LTCC_het;			// "On" or "random" or "map"
	double 		*LTCC_map;			// LTCC_map
    const char  *volds_het;         // "Off" or "random" to apply vold_ds het randomly
    double      *dyad_het_map;      // For outputting the random dyad heterogeneity
	// End sub-cellular maps==========//|
	
	// Single-dyad model ============\\|
	const char * Dyad_geo_file;
	const char * RyR_map_file;
	const char * RyR_phos_map_file;
	double 	dx;	// spatial step in um
	double  D;	// Diffusion coefficient; uM/um
	// End Single-dyad model ========//|
}CRU_variables;
// End CRU variables ============================================================================//|

// Dyad variables (local vol, LTCC and RyR fluxes ===============================================\\|
typedef struct{

	int     ex_switch;   // Tracks whether cell is currently excited

	// Dyad structure, size and state
	double 	vol_ds;		// volume of local dyadic space
	int 	NRyR;		// Number of RyRs in local dyad
	int 	NLTCC;		// Number of L-type calcium channels in local dyad
	int		active;		// 0 for not active, 1 for active
	double rand;		// local random number

	// RyR model ====================\\|
	double J_rel;		// Intracellular Ca2+ release 	(uM/ms)
	double K_rel;		// Intracellular Ca2+ release parameter for Ca_ds no flux approx
	double *rand_RyR;	// Array of random numbers for RyR model
	double Grel;		// Local copy of Grel ->J_rel scale
	double GRyR_kCO;	// Local copy of GRyR_kCO ->Open rate scale
	double Krel_SRF_mult; // multiplier for SRF

	// States and rates (stochastic model)
	double RyR_kCO;  	// Rate from closed to open (activated or inactivated) (ms^-1)
    double RyR_kOC;     // Rate from open to closed                            (ms^-1)
    double RyR_kAI;     // Rate from activated to inactivated                  (ms^-1)
    double RyR_kIA;     // Rate from inactivated to activated                  (ms^-1)

	int *RyR_state;		// Current state (0-3) of single RyR
	int	NRyR_CA;		// Number of channels (per dyad) in state Closed, Active
	int	NRyR_CI;		// Number of channels (per dyad) in state Closed, Inactivated
	int	NRyR_OA;		// Number of channels (per dyad) in state Open, Active (*flux state)
	int	NRyR_OI;		// Number of channels (per dyad) in state Openm Inactivated

	// 0D specific varaibles
	double RyR_kOC_A;	// Open-closed rate multiplier
	double NRyR_O1_det;	// Proportion of open channels (det model, 2state, not including monomer inac)
	double NRyR_O_det;	// open channels (actual including monomer inac)
	double NRyR_O_SRF;	// open channels from spontaneous release functions	
	double NRyR_O;		// actual open proportion (SRF OR det)
	double Ca_JSR_t_ex;	// Ca_JSR concentration at time of excitation
	double SRF_prop_active;	// The proportion of active channels assocaited with SRF release
	// End RyR model ================//|

	// LTCC model ===================\\|
	double J_CaL;		// Ca2+ flux through LTCCs				(uM/ms)
	double *rand_LTCC; 	// Array of random numbers for LTCC model 

	// State occupancy (stochastic)
	int NLTCC_O;			// Number of LTCCs in open state
	int *LTCC_va_state;     // 0-2, voltage activation state (state 2 is open)
	int *LTCC_vi_state;     // 0-1, voltage inactivation state (0 is inactivated; 1 is NOT inactivated)
	int *LTCC_ci_state;     // 0-1, Ca inactivation state (same as vi)
	double GLTCC_kva1_va2;	// Local copy of GLTCC_kva1_va2 -> open rate scale

	// State variables (deterministic)
	double ICaL_va_0;   	// Voltage-dependent activation, state 0 (3 state gate)
    double ICaL_va_1;   	// state 1
    double ICaL_va_2;   	// State 2 (active | open state)
    double ICaL_vi;     	// Voltage-dependent inactivation
    double ICaL_ci;     	// Ca2+ -dependent inactivation

	// steady-state, rates, tau
	double ICaL_va_al_01;   // alpha rate from state va 0 to 1						(ms^-1)
	double ICaL_va_b_01;    // beta rate state va 0 to 1 (i.e., rate 1 -> 0) 		(ms^-1)
	double ICaL_va_al_12;   // alpha, 1->2											(ms^-1)
	double ICaL_va_b_12;    // beta, 1->2											(ms^-1)
	double ICaL_va_ss;      // steady-state va (voltage dependent part, states 0-1)	(ms^-1)
	double ICaL_va_tau;     // time constant ms (voltage dependent part, states 0-1)(ms^-1)

	double ICaL_vi_al;      // alpha vi												(ms^-1)
	double ICaL_vi_b;		//														(ms^-1)
	double ICaL_vi_ss;		//														(ms^-1)
	double ICaL_vi_tau;		//														(ms^-1)
	
	double ICaL_ci_al;      // Ca dependent inactivation, alpha						(ms^-1)
	double ICaL_ci_b;		//														(ms^-1)
	double ICaL_ci_ss;		//														(ms^-1)
	double ICaL_ci_tau;		//														(ms^-1)

	double Ca_Ca_bar;       // Ca divided by Ca_bar
	double LTCC_bar;		// dynamic flux rate									(uM/ms)
	// End LTCC model ===============//|

	// CSQN/monomer inactivation ====\\|
	double csqn;        // concenration of csqn 				(mM)
	double csqn_ca;		// concentration of csqn bound to Ca 	(mM) 
	double Monomer;     // proportion of monomer form
	double Mi;          // proportion of monomer bound to inactivate RyR
	double mon_al;		// alpha transition rate				(ms^-1)
	double mon_b;		// beta									(ms^-1)
	double mon_ss;		// steady state
	double mi_ss;		// state state
	double mi_al;		// alpha								(ms^-1)
	double mi_b;		// beta									(ms^-1)
	// End CSQN/monomer inactivation=//|

}Dyad_variables;
// End dyad variables (local vol, LTCC and RyR fluxes ===========================================//|

// SR fluxes variables (SERCA and Jleak =========================================================\\|
typedef struct{

	double J_SERCA;		// SERCA uptake flux 	(uM/ms)
	double J_leak;		// SR Ca2+ leak			(uM/ms)
	double cai_term;	// Cai term of Jup
	double casr_term;	// CaSR term of Jup
	double Gup;			// Local copy of J_up_scale
	double Gleak;		// Local copy of J_leak scale

}SR_fluxes;
// End SR fluxes ================================================================================//|

// Membrane Ca2+ fluxes (NCX, Cap, Cab) =========================================================\\|
typedef struct{

	double J_NCX;    	// Flux through NCX				(uM/ms)
	double J_NCX_bulk;	// Flux through NCX, bulk space	(uM/ms)
	double J_NCX_ss;	// Flux through NCX, subspace	(uM/ms)
	double GNCX;		// Local copy of scale factor
	double NCX_SRF_mult;// multiplier specific to SRF

	double J_Cab;    	// Flux through Cab				(uM/ms)
	double J_Cab_bulk;	// Flux through Cab, bulk space	(uM/ms)
	double J_Cab_ss;	// Flux through Cab, subspace	(uM/ms)
	double GCab;		// Local copy of scale factor

	double J_CaP;    	// Flux through Cap				(uM/ms)
	double J_CaP_bulk;	// Flux through Cap, bulk space	(uM/ms)
	double J_CaP_ss;	// Flux through Cap, subspace	(uM/ms)
	double GCaP;		// Local copy of scale factor

	double J_MEM_bulk;	// Total flux through membrane	(uM/ms)
	double J_MEM_ss;	// Total flux through membrane	(uM/ms)

	double z, Ka, t1, t2, t3, numerator, denomenator; // NCX calculated variables

}Membrane_fluxes;
// End MEM_fluxes ===============================================================================//|

typedef struct{
    MTRand mtrand1;     // Random number
    double rand;
} RAND;
// End Structs used in novel (3D + integrated) Ca handling - unused in "native" models ====================//|

// Define the Spatial_coupling struct ===========================================================\\|
typedef struct{

	// Geometry and array sizes
	int NX, NY, NZ;
	int N;	// Total number of cells/nodes
	int Njunc;	// total number of cell-cell connections

	// Space step
	double dx, dy, dz;

	// Arrays =====================================================\\|
	// Geometry
	int *geo;         // contains the geometry (NX*NY*NZ)
	int *geo_linear;  // contains linearised geometry (Ncell)
	int *geo_index;   // returns the index of the 1D array (i; geo_linear) when passed 3D array (x, y, z; geo) value
	int *geo_3D_index;	// returns 3D index from Ncell index
	int *x_index;		// returns the x value at each Ncell
	int *y_index;		// returns the y value at each Ncell
	int *z_index;		// returns the z value at each Ncell

	// Diffusion coefficient
	double *D;		// Baseline D, isotropic
	double *D1;		// Primary diffusion coefficient (longitudinal to fibre)
	double *D2;		// Restricted diffusion coefficinet (transverse)
	double *Dxx; // x-component of diffusion tensor array, anisotropic
	double *Dyy;
	double *Dzz;
	double *Dxy;
	double *Dxz;
	double *Dyz;

	// anisotropic D differentials
	double *dDxx_dx;    // differential in x direction of Dxx
	double *dDxy_dx;
	double *dDxz_dx;
	double *dDyy_dy;
	double *dDxy_dy;
	double *dDyz_dy;
	double *dDzz_dz;
	double *dDxz_dz;
	double *dDyz_dz;

	//double Dscale;	    // Scaling of D (global)
    //double D_AR_scale;  // Scaling of AR (global)

	// For anisotropy and non-uniformity
	//double *D1_factor;			// Multiplication factor specifically for D1
	//double *D2_factor;          // Multiplication factor specifically for D2

	// Differential
	double *diff;	// Contains the voltage differential for diffusion (mV/ms)

	// Neighbours (1D cell indexes)
	int *xp;		// 1D index of cell at ((x+1), y, z)
	int *xm;		// 1D index of cell at ((x-1), y, z)
	int *yp;		// 1D index of cell at (x, (y+1), z)
	int *ym;		// 1D index of cell at (x, (y-1), z)
	int *zp;		// 1D index of cell at (x, y, (z+1))
	int *zm;		// 1D index of cell at (x, y, (z-1))

	int *xp_yp;		// ((x+1), (y+1), z)
    int *xp_ym;		// ((x+1), (y-1), z)
    int *xp_zp;		// ((x+1), y, (z+1))
    int *xp_zm;		// ((x+1), y, (z-1))
				
    int *xm_yp;		// etc
    int *xm_ym;
    int *xm_zp;
    int *xm_zm;

    int *yp_zp;
    int *yp_zm;
    int *ym_zp;
    int *ym_zm;

	int *xm_ym_zm;
	int *xm_ym_zp;
	int *xm_yp_zm;
	int *xm_yp_zp;
	int *xp_ym_zm;
	int *xp_ym_zp;
	int *xp_yp_zm;
	int *xp_yp_zp;

	// local orientations
	double *ox;		// Component of orientation vector in x direction
	double *oy;		// Component of orientation vector in y direction
	double *oz;		// Component of orientation vector in z direction

	// laplacian arrays
	double *lap_self;
	double *lap_xm;
	double *lap_xp;
	double *lap_ym;
	double *lap_yp;
	double *lap_zm;
	double *lap_zp;
	double *lap_xm_ym;
	double *lap_xm_yp;
	double *lap_xp_ym;
	double *lap_xp_yp;
	double *lap_xm_zm;
	double *lap_xm_zp;
	double *lap_xp_zm;
	double *lap_xp_zp;
	double *lap_ym_zm;
	double *lap_ym_zp;
	double *lap_yp_zm;
	double *lap_yp_zp;
	// can also add the corners here if needed
	// End arrays =================================================//|


	// Network model arrays || Ncell ================\\|
	double  Gt; // transverse conductance
	double  Gl; // longitudinal conductance

	double  *gDiff_node_xx;
	double  *gDiff_node_yy;
	double  *gDiff_node_zz;
	double  *gDiff_node_xypp;
	double  *gDiff_node_xypm;
	double  *gDiff_node_xzpp;
	double  *gDiff_node_xzpm;
	double  *gDiff_node_yzpp;
	double  *gDiff_node_yzpm;
	double  *gDiff_node_xyzppp;
	double  *gDiff_node_xyzppm;
	double  *gDiff_node_xyzpmp;
	double  *gDiff_node_xyzmpp;

	// Njunc arrays
	double  *gDiff_jn;  // junction conductance
	double  *IDiff;     // junction current
	int     *jn_map_minus;  // returns Ncell array index of the minus cell for IDiff[n]
	int     *jn_map_plus;
	// End network model arrays =====================//|
}SC_variables;
// End define the spatial coupling struct =======================================================//|

// Define the Tissue_parameters struct ==========================================================\\|
typedef struct{

	// Global settings
	char const *Tissue_order;		// 1D, 2D, 3D, geometry
	char const *Tissue_model;		// specific tissue model to be run
	char const *Tissue_type;		// homogeneous or heterogeneous
	char const *Orientation_type;	// isotropic, aniostropic, orthotropic
	char const *D_uniformity;		// uniform, non-uniform (regional or map)	
	char const *S1_loc_type;		// "edge", "centre" or "specified"
	char const *S2_loc_type;		// "S1", "edge", "centre" or "specified"
	char const *S1_shape;			// cube or sphere
	char const *S2_shape;			// cube or sphere
	double		Dscale;				// scaling of D (homogeneous or by map)
    double      D_AR_scale;         // scaling of anisotropi ratio
    const char *Dscale_map_on;      // "On" or "Off" to apply Dscale homogeneosuly or by map
    const char *D_AR_scale_map_on;  // "On" or "Off" to apply Dscale homogeneosuly or by map
	char const *ISO_map_on;			// "On" or "Off" to assign ISO from map or homogeneous
	char const *remod_map_on;		// "On" or "Off" to assign remodelling from map or homogeneous
	char const *ACh_map_on;			// "On" or "Off" to assign ACh from map or homogeneous
	char const *SRF_map_on;			// "On" or "Off" to assign SRF from map or homogeneous
    char const *Direct_modulation_map_on; // "On" or "Off" to assign command line direct mod from map or homogeneous
    char const *Global_orientation_direction;
    char const *spatial_gradient_map_on; // "Off", "On" or a map specifier 
    char const *Default_model;      // Whether we want a tissue model to have a defaulted cell model
    char const *Tissue_model_2;    // Whether we want a tissue model to have a defaulted cell model for second region
    char const *Multiple_models;    // To use two different cell models for different regions in tissue

	// Model-specific settings
	int NX, NY, NZ;			    // Dimension sizes
	double D1, D2, D_AR, Diso;	// Such that can be set model-specifically
	double dx, dy, dz;		    // Spatial steps
	double Gt, Gl, Giso;	    // network model conductances

	// Stimulus location
	int S1_x_loc;		// x-location of stimulus centre
	int S1_x_size;		// extent EITHER SIDE of stimulus centre
	int S1_y_loc;		// y-location of stimulus centre
	int S1_y_size;		// extent EITHER SIDE of stimulus centre
	int S1_z_loc;		// z-location of stimulus centre
	int S1_z_size;		// extent EITHER SIDE of stimulus centre

	int S2_x_loc;		// x-location of stimulus centre
	int S2_x_size;		// extent EITHER SIDE of stimulus centre
	int S2_y_loc;		// y-location of stimulus centre
	int S2_y_size;		// extent EITHER SIDE of stimulus centre
	int S2_z_loc;		// z-location of stimulus centre
	int S2_z_size;		// extent EITHER SIDE of stimulus centre

	// Stimulus and other location specific arrays
	int Nstim;			// Number of nodes to stimulate
	int Nstim_S2;		// Number of nodes to stimulate
	int *stim_area;		// area to apply stimulus (= 1 for yes, = 0 for no)
	int *S2_stim_area;	// area to apply stimulus (= 1 for yes, = 0 for no)
	bool stim_set;		// tracker of whether stimulus loc/szie has been set
	bool S2_stim_set;	// tracker of whether S2 stimulus loc/size has been set

	int **multi_stim_area;	// Stim area for multiple stimulus sites and timings
	char const 	*Multi_stim;	// "On" or "Off"
	int Nstims;				// Number of different stim sites/timings
	int stim_delay[20];		// Delay (relative to first stim) for each stim

	// ISO/remodelling/D map
	double 	*ISO_map;		// ISO conc
	double 	*Dscale_base_map;	// D scaling (inherrent e.g. scaling gradient)
    double  *Dscale_mod_map;    // D_scaling modification (e.g. fibrosis patch)
	double 	*D_AR_scale_base_map;	// D_AR scaling (inherrent e.g. scaling gradient)
    double  *D_AR_scale_mod_map;    // D_AR scaling modification (e.g. fibrosis patch)
	double 	*remod_map;		// remodelling
	double 	*ACh_map;		// ACh conc
	double 	*SRF_map;		// SRF
    double  *Direct_modulation_map; // Where to apply command line direct mod arguments (e.g. Jup_Scale, INa_va_ss_vshfit, Ito_vi_tau_scale, if passed in directly)
    double  *spatial_gradient_map;
    char const *map_in_type;    // file or coords (geo only)

    // Variables for creating idealisd maps (only one map can be created at once; can apply to multiple varables
    int ideal_map_x_loc;    // x-location of stimulus centre
    int ideal_map_x_size;   // extent EITHER SIDE of stimulus centre
    int ideal_map_y_loc;    // y-location of stimulus centre
    int ideal_map_y_size;   // extent EITHER SIDE of stimulus centre
    int ideal_map_z_loc;    // z-location of stimulus centre
    int ideal_map_z_size;   // extent EITHER SIDE of stimulus centre
    char const *ideal_map_shape; // cube or sphere
    bool ideal_map_set;     // track whether the above have been set, or to apply default

	// Phase map for re-entry
	int *phasemap;		// map for phase re-entry - 0-200 or so

	// Heterogeneity and non-uniformity parameters
	int Ncelltypes;						// Number of different celltypes present in the tissue model
	const char * celltype_number[20];	// Contains the string referencing each celltype at cell i
	int	het_junction_X_location[20];	// Location (in x) to split celltypes
	double non_uniform_D1_scale[20];	// Scales D2 (or Diso) for non-uniform diffusion coefficients by celltype
	double non_uniform_AR_scale[20];	// Scales the anisotropy ratio (D1 relative to D2)
	const char * Modeltype_number[20];	// Contains the string referencing each Model at cell i

	// Global orientation directions
	double OX;		// Component of orientation in x-direction, etc
	double OY;
	double OZ;

	// Conduction velocity measurement variables
	double 	dist_x;     // distance between x_1 and x_2 in x direction
	double 	dist_y;
    double  dist_z;
    double  dist_xy;
    double  dist_xz;
    double  dist_yz;
    double  dist_xyz;
	int 	CV_x_1;		// 1D index reference for cell x 1
	int	 	CV_x_2;
	int 	CV_y_1;
	int 	CV_y_2;
    int     CV_z_1;
    int     CV_z_2;
	int     CV_xyp_1;
    int     CV_xyp_2;
	int     CV_xym_1;
    int     CV_xym_2;
    int     CV_xzp_1;
    int     CV_xzp_2;
    int     CV_xzm_1;
    int     CV_xzm_2;
    int     CV_yzp_1;
    int     CV_yzp_2;
    int     CV_yzm_1;
    int     CV_yzm_2;
    int     CV_xyzppp_1;
    int     CV_xyzppp_2;
    int     CV_xyzppm_1;
    int     CV_xyzppm_2;
    int     CV_xyzpmp_1;
    int     CV_xyzpmp_2;
    int     CV_xyzmpp_1;
    int     CV_xyzmpp_2;

	// Anatomically detailed model file references
	const char * geo_file;
	const char * stim_file;
	const char * S2_stim_file;
	const char * orientation_file_root;
	const char * orientation_file_type;
	const char * phase_file;
	const char * Dscale_base_map_file;
	const char * Dscale_mod_map_file;
	const char * D_AR_scale_base_map_file;
	const char * D_AR_scale_mod_map_file;
	const char * ISO_map_file;
	const char * remod_map_file;
	const char * ACh_map_file;
	const char * SRF_map_file;
	const char * Direct_modulation_map_file;
	const char * spatial_gradient_map_file;

}Tissue_parameters;
// End define the tissue parameters struct ======================================================//|

// Define the Spontaneous Release Functions =====================================================\\|
typedef struct{

	// Modes, models, parameter sets
	const char * Mode;			// "Defined" or "Dynamic"
	const char * Model;			// "3D_cell" or "General"  ; Dynamic modes only
	const char * Pset;			// Specific parameter set identifier (applies to defined and dynamic)
	const char * SRF_het;		// tissue heterogeneity of SRF

	// Switches and tracking variables
	int 	srf_set;			// Whether or not waveform parameters have been SET (not initiated); 0 = ready to be set, 1 = set, -1 means not ready, -2 means been set as not initiating so don't recalc
	int 	srf_calc;			// = 1 to indicate need to calculate SRF
	int		waveform_init;		// Whether the waveform has actually been initiated
    bool    init_write_flag;    // Whether the properties of the initiated SRF need to be written to file
	double 	tset;				// time at which SRF function was set (and thus ti etc is relative to)
	double 	tinit;				// time at which SRF function is actually started (NRyR > 0)
	double 	SRF_prop_active;	// Proportion of active CRUs associated with SRF waveform (approximate)
	double 	CaSR_t_calc;		// CaSR at the time SRF is calculated (required as gobal CaSR != local CaSR during wave)

	// Random number array
	double rand[10];

	// SRF variables
	// Open RyR associated with SRF
	double NRyRo;				// Proportion of open RyRs associated with SRF
	double Fn1;					// Sigmoidal function 1 describing SRF waveform
	double Fn2;					// Sigmoidal function 2 describing SRF waveform
	double Fn3;					// Sigmoidal function 3 for SRF waveform - plateau function 1
	double Fn4;					// Sigmoidal function 4 for SRF waveform - plateau function 2

	// Probability of event
	double PSCRE;				// Probability of spontaneous calcium release, 0-1

	// Waveform parameters
	double ti;					// Initiation time of the waveform
	double duration;			// Duration of the waveform
	double t_to_peak;			// Time from ti to peak of waveform
	double decay_time;			// Time from peak to end of waveform
	double NRyRo_peak;			// peak of open RyR
	double k1_waveform;			// gradient parameter function 1 of waveform
	double k2_waveform;			// gradient parameter function 2 of waveform
	double thalf_1;				// half maximal activation time for function 1 of waveform
	double thalf_2;				// half maximal activation time for function 2 of waveform

	// Long waveform parameters
	double NRyRo_plateau;		// Proportion of RyR open at plateau
	double ti_plateau;			// Initiation time of plateau (ti is used for the spike)
	double thalf_plateau_1;		// half maximal activation time for function 1 of waveform plateau
	double thalf_plateau_2;		// half maximal activation time for function 2 of waveform plateau
	double k1_plateau;			// gradient parameter function 1 of plateau waveform
	double k2_plateau;			// gradient parameter function 2 of plateau waveform

	// Distribution parameters
	// ti
	double ti_sep;				// initiation time at separation point of two functions
	double CF_ti_sep;			// Cumulative frequency at the separation point of the two functions
	double k_ti_F1;				// Initiation time gradient parameter, F1 (left of ti_sep)
	double k_ti_F2;				// Initiation time gradient parameter, F2 (right of ti_sep)
	double k_ti_F1_ms;			// width in ms to define gradient parameter
	double k_ti_F2_ms;			// width in ms to define gradient parameter

	// duration
	double MD;					// Median duartion (ms)
	double k_D_F1;				// Width of duration dist function 1
	double k_D_F2;				// Width of duration dist function 2
	double k_D_F1_ms;			// Width of duration dist function 1 in ms
	double k_D_F2_ms;			// Width of duration dist function 2 in ms
	double DW_k1_A;				// Coefficient for duration width 1 as a function of MD
	double DW_k1_a;				// half-activation constant for duration width 1 as a function of MD
	double DW_k1_min;			// minimal value for duration width 1 as a function of MD
	double DW_k1_k;				// gradient parameter for duration width 1 as a function of MD
	double DW_k2_A;				// Coefficient for duration width 2 as a function of MD
	double DW_k2_a;				// half-activation constant for duration width 2 as a function of MD
	double DW_k2_min;			// minimal value for duration width 2 as a function of MD
	double DW_k2_k;				// gradient parameter for duration width 2 as a function of MD
	double duration_width;		// width of dist if not defining from MS, ms

	// RyR
	double NRyRo_peak_median;		// median value for the peak open RyR
	double NRyRo_peak_median_A;		// Coefficient for determining NRyRo_peak_median from duration
	double NRyRo_peak_median_min;	// minimal value of NRyRo_peak
	double NRyRo_peak_median_H;		// power coefficient

	// Dynamic model parameters which define distribution parameters from CaSR
	// Probabiliry
	double  PSCR_threshold;			// CaSR value for probability of release = 0.5 || mM -> losely defines the threshold for release (which is actually when this fn > 0)
	double  PSCR_k;					// gradient parameter for PSCR sigmoid

	// ti distribution
	double  ti_sep_a;				// CaSR reference value for ti = f(CaSR)
	double  ti_sep_A;				// Ceofficinet for ti = f(CaSR) (= maximum value - minimum value)
	double  ti_sep_k;				// gradient parameter for this function
	double  ti_sep_min;				// minimum value of ti
	double  CF_ti_sep_a;			// CaSR refernece value for CF at ti = f(CaSR)
	double  CF_ti_sep_A;			// Coefficient for CF at ti = f(CaSR) function (= maximum value - 0.05)
	double  CF_ti_sep_k;			// gradient parameter for this function
	double  KF1_a;					// CaSR reference value for the function defining the width of ti distribution to left of ti
	double  KF1_A;					// Coefficient for width of ti disitrbution function
	double  KF1_k;					// Gradient parameter for this function
	double  KF1_min;				// Minimum value of this width
	double  KF2_A;					// Coefficient 1 for function defining width of ti distribution to the right of ti	
	double  KF2_B;					// Coefficient 2 for function defining width of ti distribution to the right of ti
	double  KF2_H1;					// Power of first term in this function
	double  KF2_H2;					// Power of second term in this function
	double  KF2_min;				// Minimum valule of this width

	// duration distribution
	double  MD_a;					// CaSR reference for the function defining median duration (MD) as a funtion of CaSR
	double  MD_A;					// Coefficient of this function		
	double  MD_k;					// Gradient parameter for this function
	double  MD_min;					// Miminum value of MD 

	// General, controllable dynamic model (some are same as prior; only new here)
	double 	CaSR_min;				// minimum value of CaSR for which PSCR > 0 (specifically = 0.02)
	double 	CaSR_max;				// maximum value of CaSR, above which distributions converge
	double 	CaSR_Prange;			// range in CaSR over which probability goes from 0 (0.02; PSCR_min) to 1 (0.98; PSCR_threshold + PSCR_min)
	double 	ti_sep_max;				// maximum value of ti_sep (ms)	(i.e. at CaSR_min)
	double 	ti_width_min;			// minimum widthn of ti distribution (i.e. at CaSR_max)
	double 	ti_width_max;			// maximum widthn of ti distribution (i.e. at CaSR_min)
	double 	MD_max;					// maximum value of median duration (ms) (i.e. at CaSR_min)
	double 	duration_width_min;		// minimum widthn of duration distribution (i.e. at CaSR_max)
	double 	duration_width_max;		// maximum widthn of duration distribution (i.e. at CaSR_min)
	double 	CaSR_width_H;			// Power of function defiining how rapidly, between min and max CaSR, dist widths change

	// Threshold change for resetting dynamic
	double 	recalc_SR_diff;			// uM (!)

	// Variability scaling parameters (to fix or maximal scale constant dists)
	double t_to_peak_dist_scale;	// Scales the timing of peak from randomly between bounds (1) to no variance, in middle (0)
	double NRyRo_peak_dist_scale;	// Scales the amount NRyRo_peak can vary for a given duration
	double NRyRo_plateau_dist_scale;// Scales the amount NRyRo_plateau can vary for a given duration

}Spontaneous_release_functions;
// End the Spontaneous Release Functions ========================================================//|

// Define the Arguments struct ==================================================================\\|
typedef struct {

	// Reference ==================================================\\|
	char const	*reference;			    // Any reference to identify simulation global directory
	bool		reference_arg;		    // true IF argument is passed
	char const	*results_reference;		// Any reference to identify simulation results directory
	bool		results_reference_arg;	// true IF argument is passed
	char const	*state_reference_read;		// Any reference to identify specific state files
	bool		state_reference_read_arg;	// true IF argument is passed
	char const	*state_reference_write;		// Any reference to identify specific state files
	bool		state_reference_write_arg;	// true IF argument is passed
	// End reference ==============================================//|

	// Simulation settings ========================================\\|
    char        argin_sf_in[500][500]; // Array to hold values read from settings file
    char        *argin_sf[500];      // Array to pass settings file values into set args function
    bool        settings_file;      // whether settings file has been passed
    int         Narg_settings;      // Number of arguments, settings file
	int         BCL;                // Basic cycle length, ms
	bool        BCL_arg;            // True IF BCL argument has been passed
	int         Total_time;         // Total simulation time, ms
	bool        Total_time_arg;     // True IF Total time argument has been passed
	int         Paced_time;         // Time during which stimulus is applied
	bool        Paced_time_arg;     // True IF paced time argument has been passed
	int         NBeats;             // Number of applied stimuli
	bool        NBeats_arg;         // True IF NBeats argument has been passed
	int			S2_CL;				// CL of S2 stimulus
	bool 		S2_arg;				// True IF S2 CL argument passed
	int 		NS2;				// Number S2 stimuli
	bool		NS2_arg;			// True IF NS2 argument passed
	double      dt;                 // Simulation integration time-step
	bool        dt_arg;             // True IF dt argument has been passed
	char const	*Vclamp;			// "On" or "Off"
	char const  *Write_state;		// "On" or "Off"	
	char const  *Read_state;		// "On" or "Off"
	int			SOI;				// Spatial Output Interval (VTK)
	bool		SOI_arg;			// True IF argument passed (VTK)
	int			SOId;				// Spatial Output Interval (data)
	bool		SOId_arg;			// True IF argument passed (data)
    int         SORs;               // Spatial output range start time
    bool        SORs_arg;           // True IF argument passed
    int         SORe;               // Spatial output range start time
    bool        SORe_arg;           // True IF argument passed
	char const 	*Multi_stim;		// "On" or "Off" for multiple stim sites
	bool		Multi_stim_arg;		//	True IF argument passed 
	// End simulation settings ====================================//|

	// Model and cell conditions ==================================\\|
	char const  *Model;             // String containing selected model
	bool        Model_arg;          // True IF Model argument has been passed
	char const  *Celltype;          // String containing selected model
	bool        Celltype_arg;       // True IF argument has been passed
	char const	*Agent;				// String containing pharma agent
	bool		Agent_arg;			// True IF argument passed
	char const	*Remodelling;		// String containing remodelling
	bool		Remodelling_arg;	// True IF argument passed
	double      ISO;				// Concentration of ISO
	bool		ISO_arg;			// True IF argument passed
	char const  *ISO_model;			// Which ISO model is to be used
	bool		ISO_model_arg;		// True IF argument passed
	char const	*Mutation;			// Mutation
	bool 		Mutation_arg;		// True IF argument has been passed
	double      ACh;                // Concentration of ACh
	bool        ACh_arg;            // True IF argument passed
	char const  *ACh_model;         // Which ACh model is to be used
	bool        ACh_model_arg;      // True IF argument passed

	char const 	*environment;		// For human single cell only; isolated vs intact cell environment 
	bool		environment_arg;	// True IF environment argument is passed

	double		AIhyp;				// Magnitude of applied hyperpolarising current (pA/pF)
	bool		Ihyp_arg;			// True IF argument has been passed

	double		Remodelling_prop;	// Proportion of remodelling to implement (linear scale control - remod)
	bool		Remodelling_prop_arg;// True IF argument has been passed 
	double		Agent_prop;			// Proportion of agent to implement (linear scale control - saturating)
	bool		Agent_prop_arg;		// True IF argument has been passed 
    
    const char  *spatial_gradient;
    bool        spatial_gradient_arg;
    double      spatial_gradient_prop;
    bool        spatial_gradient_prop_arg;

	// Global control, spatial cell model only
	char const* tau_ss_type;        // Identifier for time constant of coupling of sub-space
	bool		tau_ss_arg;			// True IF argument has been passed
	// End model and cell conditions ==============================//|

	// Tissue model settings ======================================\\|
	char const 	*Tissue_order;			// 1D, 2D, 3D or geo
	bool		Tissue_order_arg;		// True IF argument has been passed
	char const 	*Tissue_model;			// Model type specifier
	bool		Tissue_model_arg;		// True IF argument has been passed
	char const	*Tissue_type;			// Homogeneous or heterogeneous
	bool		Tissue_type_arg;		// True IF argument has been passed
	char const 	*Orientation_type;		// isotropic, anisotropic or orthotropic
	bool		Orientation_type_arg;	// True IF argument has been passed
	char const	*Stimulus_loc_type;		// Specifier for type of stimulus (idealised models only)
	bool		Stimulus_type_arg;		// True IF argument has been passed
	char const	*S2_Stimulus_loc_type;	// Specifier for type of stimulus (idealised models only)
	bool		S2_Stimulus_type_arg;	// True IF argument has been passed
	char const 	*S1_shape;				// shape of S1 stimulus
	bool		S1_shape_arg;			// True IF arugmnet passed
	char const 	*S2_shape;				// shape of S2 stimulus
	bool		S2_shape_arg;			// True IF arugmnet passed
	char const 	*D_uniformity;			// uniform, non-uniform	
	bool 		D_uniformity_arg;		// True IF argument has been passed
	double 		Dscale;					// Scales D homogeneously
	bool		Dscale_arg;				// True IF argument has been passed
	double 		D_AR_scale;				// Scales AR homogeneously
	bool		D_AR_scale_arg;		    // True IF argument has been passed
    const char  *Dscale_mod_map_on;     // On or Off to apply Dscale by map
    bool        Dscale_mod_map_arg;     // True IF argument has been passed
    const char  *D_AR_scale_mod_map_on;     // On or Off to apply D_AR_scale by map
    bool        D_AR_scale_mod_map_arg;     // True IF argument has been passed
	char const	*Remodelling_map_on;	// Whether ot not to employ a remodelling map
	bool		Remodelling_map_arg;	// True IF argument has been passed
	char const	*ISO_map_on;			// Whether ot not to employ a remodelling map
	bool		ISO_map_arg;			// True IF argument has been passed
	char const	*ACh_map_on;			// Whether ot not to employ a remodelling map
	bool		ACh_map_arg;			// True IF argument has been passed
	char const	*SRF_map_on;			// Whether ot not to employ a remodelling map
	bool		SRF_map_arg;			// True IF argument has been passed
	char const	*Direct_modulation_map_on;			// Whether ot not to employ a remodelling map
	bool		Direct_modulation_map_arg;			// True IF argument has been passed
    char const  *Multiple_models;    // To use two different cell models for different regions in tissue
    bool        Multiple_models_arg;
    char const  *Tissue_model_2;       // Second tissue model
    bool        Tissue_model_2_arg;

	// Controlling S1 and S2 via arguments
	int 	S1_x_loc;       // x-location of stimulus centre
	bool 	S1_x_loc_arg;	// True IF argument has been passed
	int 	S1_x_size;      // extent EITHER SIDE of stimulus centre
	bool 	S1_x_size_arg;	// True IF argument has been passed
	int 	S1_y_loc;       // y-location of stimulus centre
	bool 	S1_y_loc_arg;	// True IF argument has been passed
	int 	S1_y_size;      // extent EITHER SIDE of stimulus centre
	bool 	S1_y_size_arg;	// True IF argument has been passed
	int 	S1_z_loc;       // z-location of stimulus centre
	bool 	S1_z_loc_arg;	// True IF argument has been passed
	int 	S1_z_size;      // extent EITHER SIDE of stimulus centre
	bool 	S1_z_size_arg;	// True IF argument has been passed
	int 	S2_x_loc;       // x-location of stimulus centre
	bool 	S2_x_loc_arg;	// True IF argument has been passed
	int 	S2_x_size;      // extent EITHER SIDE of stimulus centre
	bool 	S2_x_size_arg;	// True IF argument has been passed
	int 	S2_y_loc;       // y-location of stimulus centre
	bool 	S2_y_loc_arg;	// True IF argument has been passed
	int 	S2_y_size;      // extent EITHER SIDE of stimulus centre
	bool 	S2_y_size_arg;	// True IF argument has been passed
	int 	S2_z_loc;       // z-location of stimulus centre
	bool 	S2_z_loc_arg;	// True IF argument has been passed
	int 	S2_z_size;      // extent EITHER SIDE of stimulus centre
	bool 	S2_z_size_arg;	// True IF argument has been passed

	// Controlling D and dx directly (+ orientation)
	double	D1;				// primary component of diffusion
	bool	D1_arg;
	double 	D_AR;
	bool	D_AR_arg;
	double 	dx;
	bool	dx_arg;
	double  OX;
	bool 	OX_arg;
	double  OY;
	bool 	OY_arg;
	double  OZ;
	bool 	OZ_arg;
    char const  *Global_orientation_direction;
    bool        Global_orientation_direction_arg; 

    // ISO/ACh/remod/SRF/D_map idealised coordinates
	char const 	*map_shape;				
	bool		map_shape_arg;			
    char const  *map_in_type;           
    bool        map_in_type_arg;
    double      map_x_loc;
    bool        map_x_loc_arg;
    double      map_y_loc;
    bool        map_y_loc_arg;
    double      map_z_loc;
    bool        map_z_loc_arg;
    double      map_x_size;
    bool        map_x_size_arg;
    double      map_y_size;
    bool        map_y_size_arg;
    double      map_z_size;
    bool        map_z_size_arg;

	// Map files
	const char * stim_file;
	const char * S2_stim_file;
	const char * phase_file;
	//const char * Dscale_map_file;
    const char * Dscale_base_map_file;
    const char * Dscale_mod_map_file;
    const char * D_AR_scale_base_map_file;
    const char * D_AR_scale_mod_map_file;
	const char * ISO_map_file;
	const char * remod_map_file;
	const char * ACh_map_file;
	const char * SRF_map_file;
	const char * Direct_modulation_map_file;
	const char * spatial_gradient_map_file;
	bool 	stim_file_arg;
	bool 	S2_stim_file_arg;
	bool 	phase_file_arg;
	//bool	Dscale_map_file_arg;
	bool	Dscale_base_map_file_arg;
	bool	Dscale_mod_map_file_arg;
	bool	D_AR_scale_base_map_file_arg;
	bool	D_AR_scale_mod_map_file_arg;
	bool 	ISO_map_file_arg;
	bool 	remod_map_file_arg;
	bool 	ACh_map_file_arg;
	bool 	SRF_map_file_arg;
	bool 	Direct_modulation_map_file_arg;
	bool 	spatial_gradient_map_file_arg;
	// Tissue model settings ======================================//|

	// Spatial single cell model settings =========================\\|
	char const 	*Cell_size;				// Reference for cell size
	bool		Cell_size_arg;			// True IF argument has been passed
	char const 	*Sim_cell_size;			// Reference for simulation cell size
	bool        Sim_cell_size_arg;      // True IF argument has been passed
	double		Cai_IC;					// Initial condition for Cai
	bool		Cai_IC_arg;				// True IF argument has been passed
	double		CaSR_IC;				// Initial condition for CaSR
	bool		CaSR_IC_arg;			// True IF argument has been passed
	const char  *Detub;             	// "On" or "Off" to read in a detubulation map
	bool		Detub_arg;				// True IF argument has been passed
    const char  *LTCC_redist;           // "On" or "Off" to redistribute LTCCs when detub is on
    bool        LTCC_redist_arg;        // True IF argument has been passed
	const char  *SERCA_het;            	// "On" or "Off" to read in a detubulation map
	bool		SERCA_het_arg;			// True IF argument has been passed
	const char  *NCX_het;            	// "On" or "Off" to read in a detubulation map
	bool		NCX_het_arg;			// True IF argument has been passed
	const char  *RyR_het;            	// "Off", "random" or "map" 
	bool		RyR_het_arg;			// True IF argument has been passed
	const char  *LTCC_het;            	// "Off", "random" or "map"
	bool		LTCC_het_arg;			// True IF argument has been passed
	char const  *TT_map_file;       	// Filename for reading in a TT map
	bool		TT_map_file_arg;		// True IF argument has been passed
	char const  *SERCA_map_file;       	// Filename for reading in a TT map
	bool		SERCA_map_file_arg;		// True IF argument has been passed
	char const  *NCX_map_file;       	// Filename for reading in a TT map
	bool		NCX_map_file_arg;		// True IF argument has been passed
	char const  *RyR_het_map_file;      // Filename for reading in a TT map
	bool		RyR_het_map_file_arg;	// True IF argument has been passed
	char const  *LTCC_map_file;       	// Filename for reading in a TT map
	bool		LTCC_map_file_arg;		// True IF argument has been passed
    const char  *volds_het;             // "On" or "Off" to apply vold_ds het randomly
    bool        volds_het_arg;          // True IF argument has been passed
	// End Spatial single cell model settings =====================//|

	// Spontaneous Release Function control =======================\\|
	char const 	*SRF_mode;				// Off, Defined or Dynamic
	bool		SRF_mode_arg;			// True IF argument has been passed
	char const 	*SRF_model;				// 3D_cell or General - dynamic only
	bool		SRF_model_arg;			// True IF argument has been passed
	char const 	*SRF_Pset;				// Specific SRF parameter set
	bool		SRF_Pset_arg;			// True IF argument has been passed
	char const 	*SRF_het;				// SRF heterogeneity 	
	bool 		SRF_het_arg;
	
	// User control =================\\|
	int 		SRF_N_DC_set;
	int 		SRF_N_Dyn_set;

	// DC model variables
	double 		SRF_PSCRE;
	double 		SRF_CF_ti_sep;
	double 		SRF_ti_sep;
	double 		SRF_ti_dist_width_1;		// left of ti sep	
	double 		SRF_ti_dist_width_2;		// right of ti sep
	double 		SRF_MD;
	double 		SRF_duration_width;

	// Dynamic model variables
	double 		SRF_PSCR_threshold;
	double 		SRF_CaSR_max;
	double		SRF_CaSR_Prange;
	double 		SRF_ti_sep_max;
	double		SRF_ti_sep_min;
	double 		SRF_ti_width_max;
	double 		SRF_ti_width_min;
	double		SRF_MD_max;
	double		SRF_MD_min;
	double 		SRF_duration_width_max;
	double		SRF_duration_width_min;
	double 		SRF_CaSR_width_H;
	// End user control =============//|
	// End Spontaneous Release Function control ===================//|

	// Direct current modification ================================\\|
	// Simple scaling 
	double GNa;                 // Scale factor for fast sodium current (and Ip0d)
	double GNaL;                // Scale factor for late sodium current
	double Gto;                 // Scale factor for Ito (and Ip1r)
	double GCaL;                // Scale factor for ICaL (and Ip2d)
	double GKur;                // Scale factor for IKur (and Ip2r)
	double GKr;                 // Scale factor for IKr (and Ip3r)
	double GKs;                 // Scale factor for IKs
	double GK1;                 // Scale factor for IK1
	double GNCX;                // Scale factor for INCX
	double GCaP;                // Scale factor for ICaP
	double GNab;                // Scale factor for INab
	double GCab;                // Scale factor for ICab
	double GKb;                 // Scale factor for IKb
	double GNaK;                // Scale factor for INaK
	double GClCa;               // Scale factor for IClCa
	double GKACh;              // Scale factor for IKACh
	double GClb;               	// Scale factor for IClb

	// Time constant scaling
	double INa_va_tau_scale;   	// Scales time constant for INa voltage activation
	double INa_vi_1_tau_scale;  // Scales time constant for INa voltage inactivation 1
	double INa_vi_2_tau_scale;  // Scales time constant for INa voltage inactivation 2
	double INaL_va_tau_scale;   // Scales time constant for INaL voltage activation
	double INaL_vi_tau_scale;   // Scales time constant for INaL voltage inactivation
	double Ito_va_tau_scale;    // Scales time constant for Ito voltage activation
	double Ito_vi_tau_scale;    // Scales time constant for Ito voltage inactivation
	double ICaL_va_tau_scale;   // Scales time constant for ICaL voltage activation
	double ICaL_vi_tau_scale;   // Scales time constant for ICaL voltage inactivation
	double IKur_va_tau_scale;   // Scales time constant for IKur voltage activation
	double IKur_vi_tau_scale;   // Scales time constant for IKur voltage inactivation
	double IKr_va_tau_scale;    // Scales time constant for IKr voltage activation
	double IKs_va_tau_scale;    // Scales time constant for IKr voltage activation
	double IKACh_va_tau_scale;	// Scales time constant for IKACh voltage activation

	// Voltage dependence
	double INa_va_shift;        // Shift of the V of activation || alpha and beta   mV
	double INa_vi_shift;        // Shift of the V of inactivation || alpha and beta mV

	double INaL_va_shift;       // Shift of the V of activation || alpha and beta   mV
	double INaL_vi_shift;       // Shift of the V of inactivation || alpha and beta mV

	double Ito_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double Ito_vi_ss_shift;     // Shift of the V1/2 of inactivation steady state   mV
	double Ito_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double Ito_vi_tau_shift;    // Shift of the voltage dependence of time constant mV
	double Ito_va_ss_kscale;    // Scales the gradient parameter of activation steady state
	double Ito_vi_ss_kscale;    // Scales the gradient parameter of inactivation steady state
	double Ito_shift;			// Applies to all V shifts 							mV

	double ICaL_va_ss_shift;    // Shift of the V1/2 of acitvation steady state     mV
	double ICaL_vi_ss_shift;    // Shift of the V1/2 of inactivation steady state   mV
	double ICaL_va_tau_shift;   // Shift of the voltage dependence of time constant mV
	double ICaL_vi_tau_shift;   // Shift of the voltage dependence of time constant mV
	double ICaL_va_ss_kscale;   // Scales the gradient parameter of activation steady state
	double ICaL_vi_ss_kscale;   // Scales the gradient parameter of inactivation steady state
	double ICaL_shift;			// Applies to all V shifts 							mV

	double IKur_va_ss_shift;    // Shift of the V1/2 of acitvation steady state     mV
	double IKur_vi_ss_shift;    // Shift of the V1/2 of inactivation steady state   mV
	double IKur_va_tau_shift;   // Shift of the voltage dependence of time constant mV
	double IKur_vi_tau_shift;   // Shift of the voltage dependence of time constant mV
	double IKur_va_ss_kscale;   // Scales the gradient parameter of activation steady state
	double IKur_vi_ss_kscale;   // Scales the gradient parameter of inactivation steady state
	double IKur_shift;			// Applies to all V shifts 							mV

	double IKr_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double IKr_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double IKr_va_ss_kscale;    // Scales the gradient parameter of activation steady state
	double IKr_vi_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double IKr_vi_ss_kscale;    // Scales the gradient parameter of activation steady state

	double IKs_va_ss_shift;     // Shift of the V1/2 of acitvation steady state     mV
	double IKs_va_tau_shift;    // Shift of the voltage dependence of time constant mV
	double IKs_va_ss_kscale;    // Scales the gradient parameter of activation steady state

	double IKACh_va_ss_shift;	// Shift of the V1/2 of acitvation steady state     mV
	double IKACh_va_ss_kscale;	// Scales the gradient parameter of activation steady state
	double IKACh_va_tau_shift;	// Shift of the tau of acitvation steady state     mV

	double IK1_va_shift;        // Shift of the V1/2 of acitvation steady state     mV
	double IK1_Erev_shift;      // Shift of the V term in V - Erev
	// End direct current modification ============================//|

	// Ca handling modification ===================================\\|
	double Gup;                 // Scale factor of intracellular, SERCA Ca2+ uptake
	double Gleak;               // Scalr factor of intracellular Ca2+ leak 
	double Grel;				// Scale factor for intracellualr Ca2+ release 

	// Spatial model only
	double 	RyR_Po;				// Scales open rate of RyRs
	bool	RyR_Po_arg;			// True IF argument has been passed
	double 	LTCC_Po;			// SCales open rate of LTCCs
	bool	LTCC_Po_arg;		// True IF argument has been passed

	// Delayed impose CaSR functionality
	const char *Delayed_CaSR_IC;    // "On" or "Off"
	bool 		Delayed_CaSR_IC_arg;
	double      CaSR_IC_delay;      // ms
	bool 		CaSR_IC_delay_arg;
	// End Ca handling modification ===============================//|

	// Boolean switches if modulation arguments have been passed ==\\|
	bool DC_current_mod_arg;		// True IF ANY direct current mod argument has been passed
	bool GNa_arg;                	// Scale factor for fast sodium current (and Ip0d)
	bool GNaL_arg;                	// Scale factor for late sodium current
	bool Gto_arg;                 	// Scale factor for Ito (and Ip1r)
	bool GCaL_arg;                	// Scale factor for ICaL (and Ip2d)
	bool GKur_arg;                	// Scale factor for IKur (and Ip2r)
	bool GKr_arg;                 	// Scale factor for IKr (and Ip3r)
	bool GKs_arg;                 	// Scale factor for IKs
	bool GK1_arg;                 	// Scale factor for IK1
	bool GNCX_arg;                	// Scale factor for INCX
	bool GCaP_arg;                	// Scale factor for ICaP
	bool GNab_arg;                	// Scale factor for INab
	bool GCab_arg;                	// Scale factor for ICab
	bool GKb_arg;                 	// Scale factor for IKb
	bool GNaK_arg;                	// Scale factor for INaK
	bool GClCa_arg;               	// Scale factor for IClCa
	bool INa_va_tau_scale_arg;    	// Scales time constant for INa voltage activation
	bool INa_vi_1_tau_scale_arg;  	// Scales time constant for INa voltage inactivation 1
	bool INa_vi_2_tau_scale_arg;  	// Scales time constant for INa voltage inactivation 2
	bool INaL_va_tau_scale_arg;   	// Scales time constant for INaL voltage activation
	bool INaL_vi_tau_scale_arg;   	// Scales time constant for INaL voltage inactivation
	bool Ito_va_tau_scale_arg;    	// Scales time constant for Ito voltage activation
	bool Ito_vi_tau_scale_arg;    	// Scales time constant for Ito voltage inactivation
	bool ICaL_va_tau_scale_arg;   	// Scales time constant for ICaL voltage activation
	bool ICaL_vi_tau_scale_arg;   	// Scales time constant for ICaL voltage inactivation
	bool IKur_va_tau_scale_arg;   	// Scales time constant for IKur voltage activation
	bool IKur_vi_tau_scale_arg;   	// Scales time constant for IKur voltage inactivation
	bool IKr_va_tau_scale_arg;    	// Scales time constant for IKr voltage activation
	bool IKs_va_tau_scale_arg;    	// Scales time constant for IKr voltage activation
	bool INa_va_shift_arg;        	// Shift of the V of activation || alpha and beta   mV
	bool INa_vi_shift_arg;        	// Shift of the V of inactivation || alpha and beta mV
	bool INaL_va_shift_arg;       	// Shift of the V of activation || alpha and beta   mV
	bool INaL_vi_shift_arg;       	// Shift of the V of inactivation || alpha and beta mV
	bool Ito_va_ss_shift_arg;     	// Shift of the V1/2 of acitvation steady state     mV
	bool Ito_vi_ss_shift_arg;     	// Shift of the V1/2 of inactivation steady state   mV
	bool Ito_va_tau_shift_arg;    	// Shift of the voltage dependence of time constant mV
	bool Ito_vi_tau_shift_arg;    	// Shift of the voltage dependence of time constant mV
	bool Ito_va_ss_kscale_arg;    	// Scales the gradient parameter of activation steady state
	bool Ito_vi_ss_kscale_arg;    	// Scales the gradient parameter of inactivation steady state
	bool Ito_shift_arg;    
	bool ICaL_va_ss_shift_arg;    	// Shift of the V1/2 of acitvation steady state     mV
	bool ICaL_vi_ss_shift_arg;    	// Shift of the V1/2 of inactivation steady state   mV
	bool ICaL_va_tau_shift_arg;   	// Shift of the voltage dependence of time constant mV
	bool ICaL_vi_tau_shift_arg;   	// Shift of the voltage dependence of time constant mV
	bool ICaL_va_ss_kscale_arg;   	// Scales the gradient parameter of activation steady state
	bool ICaL_vi_ss_kscale_arg;   	// Scales the gradient parameter of inactivation steady state
	bool ICaL_shift_arg;    
	bool IKur_va_ss_shift_arg;    	// Shift of the V1/2 of acitvation steady state     mV
	bool IKur_vi_ss_shift_arg;    	// Shift of the V1/2 of inactivation steady state   mV
	bool IKur_va_tau_shift_arg;   	// Shift of the voltage dependence of time constant mV
	bool IKur_vi_tau_shift_arg;   	// Shift of the voltage dependence of time constant mV
	bool IKur_va_ss_kscale_arg;   	// Scales the gradient parameter of activation steady state
	bool IKur_vi_ss_kscale_arg;   	// Scales the gradient parameter of inactivation steady state
	bool IKur_shift_arg;    
	bool IKr_va_ss_shift_arg;     	// Shift of the V1/2 of acitvation steady state     mV
	bool IKr_va_tau_shift_arg;    	// Shift of the voltage dependence of time constant mV
	bool IKr_va_ss_kscale_arg;    	// Scales the gradient parameter of activation steady state
	bool IKr_vi_ss_shift_arg;     	// Shift of the V1/2 of acitvation steady state     mV
	bool IKr_vi_ss_kscale_arg;    	// Scales the gradient parameter of activation steady state
	bool IKs_va_ss_shift_arg;     	// Shift of the V1/2 of acitvation steady state     mV
	bool IKs_va_tau_shift_arg;    	// Shift of the voltage dependence of time constant mV
	bool IKs_va_ss_kscale_arg;    	// Scales the gradient parameter of activation steady state
	bool IK1_va_shift_arg;        	// Shift of the V1/2 of acitvation steady state     mV
	bool IK1_Erev_shift_arg;      	// Shift of the V term in V - Erev                  mV
	bool Gup_arg;                 	// Scale factor of intracellular, SERCA Ca2+ uptake
	bool Gleak_arg;               	// Scale factor of intracellular Ca2+ leak 
	bool Grel_arg;               	// Scale factor of intracellular Ca2+ release
	bool GClb_arg;
	bool GKACh_arg;
	bool IKACh_va_tau_scale_arg;
	bool IKACh_va_ss_shift_arg;
	bool IKACh_va_ss_kscale_arg;
	bool IKACh_va_tau_shift_arg;
	// end boolean switches =======================================\\|

}Argument_parameters;
// End Define the arguments struct ==============================================================//|


#endif

