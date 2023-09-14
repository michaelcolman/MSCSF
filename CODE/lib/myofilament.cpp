// Source code associated with  ===========================  //
// "Multi-scale cardiac simulation framework" =============  //
// For simulation of cardiac cellular and tissue dynamics =  //
// from the spatial cellular to full organ scales. ========  //
// With implementation of multiple, published cell models =  //
// as well as novel models developed in my lab. ===========  //
// ========================================================  //
// This file: Myofilament, trpn and force =================  //
// CODE BY GARETH JONES from the original paper ===========  //
// Gauthier, Greenstein, Winslow Front Physiol 2012;3:244 =  //
// Incorporating Rice, Winslow, Hunter Am J Physiol 1999 ==  //
// May;276(5 pt 2):H1734-1754. ============================  //
// NOTE: I am not convinced/satisfied by this specific ====  //
// implementation and incorporation into the spatial model=  //
// As buffering is not a key, imrpotant component of this =  //
// initial version of the model, it is left as-is for now.=  //
// However, to apply this model to something specific about  //
// buffering, this needs to be looked at. =================  //
// I also intend to update it for "version 2" of the model=  //
// whenever that is finished..... =========================  //
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

#include "myofilament.hpp"

void Myofilament::run_step_myofilament(const double cai, const double ATP, const double ADP)
{
	myo_cai = cai;
	myo_ATPi = ATP;
	myo_ADP = ADP;
	LSODA_solve();
	state_update();
}

void Myofilament::myofilament_ODE(double t, double *y, double *ydot, void *data)
{
	// Implementation of the force generation model
	double khtrpn_pos = 100; // mM-1 ms-1 - Ca2+ on-rate for troponin high affinity sites
	double khtrpn_neg = 3.3e-4; // ms-1 - Ca2+ off-rate for troponin high affinity sites
	double kltrpn_pos = 100; // mM-1 ms-1 - Ca2+ on-rate for troponin low affinity sites
	double kltrpn_neg = 4e-2; // ms-1 - Ca2+ off-rate for troponin low affinity sites
	double HTRPN_tot = 0.14; // mM - Total troponin high-affinity sites
	double LTRPN_tot = 0.7; // mM - Total troponin low-affinity sitess
	double kpn_trpn = 0.04; // ms-1 - Transition rate from tropomyosin permissive to non-permissive
	double SL_Cort = 2.15; // uM - Sarcomere length
	double fXB = 0.05; // ms-1 - Transition rate from weak to strong cross bridge
	double gXB_min = 0.1; // ms-1 - Minimum transition rate from strong to weak cross bridge
	double zeta = 0.1; // N mm-2 - Conversion factor normalizing to physiological force - fudge factor
	double VAM_max = 7.2e-3; // mM ms-1 - Maximal rate of ATP hydrolysis by myofibrils (AM ATPase)
	double KMAM_ATP = 0.03; // mM - ATP half saturation constant of AM ATPase
	double KiAM = 0.26; // mM - ADP inhibition constant of AM ATPase

	double phi_myo = 1 + ((2.3 - SL_Cort)/pow((2.3-1.7),1.6));

	double f01 = 3*fXB;
	double f12 = 10*fXB; 
	double f23 = 7*fXB;
	double g01 = 1*gXB_min;
	double g12 = 2*gXB_min;
	double g23 = 3*gXB_min;
	double g01_SL = 1*phi_myo*gXB_min;
	double g12_SL = 2*phi_myo*gXB_min;
	double g23_SL = 3*phi_myo*gXB_min;

	double KCa_trpn = kltrpn_neg/kltrpn_pos;
	double Ntrpn = 3.5*SL_Cort - 2.0;
	double Khalf_trpn = 1/(1+(KCa_trpn/(1.4e-3 - 0.8e-3*((SL_Cort-1.7)/0.6))));
	double knp_trpn = kpn_trpn * pow(((LTRPNCa)/(Khalf_trpn * LTRPN_tot)),Ntrpn);

	double path_sum = g01*g12*g23 + f01*g12*g23 + f01*f12*g23 + f01*f12*f23;
	double P1_max = f01*g12*g23/path_sum;
	double P2_max = f01*f12*g23/path_sum;
	double P3_max = f01*f12*f23/path_sum;

	// This is what's in their supplement - noooope  
	// N0 = 1 - (N1 + P0 + P1 + P2 + P3);

	// dP0/dt
	ydot[0] = -(kpn_trpn + f01)*P0 + knp_trpn * N0 + g01_SL * P1;

	// std::cout << knp_trpn << "  " << P0 << "  " << N0 << "  " << P1 << std::endl;

	// dP1/dt
	ydot[1] = -(kpn_trpn + f12 + g01_SL)* P1 + knp_trpn * N1 + f01 * P0 + g12_SL* P2;

	// dP2/dt
	ydot[2] = -(f23 + g12_SL)*P2 + f12*P1 + g23_SL*P3;

	// dP3/dt
	ydot[3] = -g23_SL*P3 + f23*P2;

	// This is in their supplement as well ????
	// N1 += dt*(kpn_trpn*P1 + (knp_trpn + g01_SL)*N1);

	// dN1/dt
	ydot[4] =  kpn_trpn*P1 - (knp_trpn + g01_SL)*N1;

	// dN0/dt
	ydot[5] = -ydot[0] - ydot[1] - ydot[2] - ydot[3] - ydot[4];    

	// dHTRPNCa/dt
	ydot[6] = ((khtrpn_pos* myo_cai * (HTRPN_tot - HTRPNCa)) - (khtrpn_neg*HTRPNCa));
	// dLTRPNCa/dt
	ydot[7] = (kltrpn_pos* myo_cai * (LTRPN_tot - LTRPNCa)) - (kltrpn_neg*(1-(2/3*Norm_force))*LTRPNCa);

	// Troponin flux
	Jtrpn = ydot[6] + ydot[7];

	Force = zeta * ((P1 + N1 + 2*P2 + 3*P3)/(P1_max + 2*P2_max + 3*P3_max));

	Norm_force = (P1 + N1 + 2*P2 + 3*P3)/(P1_max + 2*P2_max + 3*P3_max);

	// VAM
	V_AM = VAM_max * ((f01*P0 + f12*P1 + f23*P2)/(f01 + f12 + f23)) / ( 1 + KMAM_ATP/myo_ATPi*(1 + myo_ADP/KiAM));

}

void Myofilament::state_update()
{
	P0 = y_myofilament[1];
	P1 = y_myofilament[2];
	P2 = y_myofilament[3];
	P3 = y_myofilament[4];
	N1 = y_myofilament[5];
	N0 = y_myofilament[6];
	HTRPNCa = y_myofilament[7];
	LTRPNCa = y_myofilament[8];
}

// Calls the solver and passes the function below containing ODEs to be solved
void Myofilament::LSODA_solve()
{   
	lsoda_integrator_myofilament.lsoda(&Myofilament::myofilament_ODE, y_myofilament, &t_myofilament, tout,0);

	tout = tout + dt_myof;
}

void Myofilament::LSODA_set()
{
	// Set initial conditions for solver
	y_myofilament[0] = 0.0; //Dummy
	y_myofilament[1] = P0;
	y_myofilament[2] = P1;
	y_myofilament[3] = P2;
	y_myofilament[4] = P3;
	y_myofilament[5] = N1;
	y_myofilament[6] = N0;
	y_myofilament[7] = HTRPNCa;
	y_myofilament[8] = LTRPNCa;

	t_myofilament = 0;
	tout =  t_myofilament + dt_myof;

}

