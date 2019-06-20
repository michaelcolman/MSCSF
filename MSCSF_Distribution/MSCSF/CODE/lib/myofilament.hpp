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

#include <cmath>
#include "lsoda.hpp"

class Myofilament
{
	public:
		Myofilament() : lsoda_integrator_myofilament(this, NEQ_myofilament)
	{
		LSODA_set();
		//std::cout << "Set initial conditions for myofilament solver." << std::endl;

		P0 = 1.6672012132596455e-03;// 1.04e-5;
		P1 = 1.4416741873793939e-03;// 9.02e-6;
		P2 = 2.6920336520753884e-03;// 1.68e-5;
		P3 = 2.3449372176401837e-03;// 1.46e-5;
		N0 = 9.9041768766118021e-01; // 0.99;
		N1 = 1.4356068652011372e-03; // 8.99e-6;
		HTRPNCa = 1.3055735570840798e-01;
		LTRPNCa = 1.8066410206828074e-02;
		myo_Pim = 2.0;
		myo_cai = 0;
		myo_ATPi = 0;
		myo_ADP = 0;

		dt_myof = 0.025;
		t = 0;
	}

		~Myofilament()
		{
			lsoda_integrator_myofilament.n_lsoda_terminate();
		}

		void myo_out(std::ofstream &out);

		void run_step_myofilament(const double cai, const double ATP, const double ADP);
		void myofilament_ODE(double t, double *y, double *ydot, void *data);
		void state_update();
		void LSODA_solve();
		void LSODA_set();

		double dt_myof;
		double t;
		double V_AM;
		double Jtrpn;
		double Force;
		double Norm_force;

	private:
		// These long ones are from some code they sent me
		// The others are from their Qince Li 2014 paper
		double P0;
		double P1;
		double P2;
		double P3;
		double N0;
		double N1;
		double HTRPNCa;
		double LTRPNCa;
		double myo_Pim;
		double myo_cai;
		double myo_ATPi;
		double myo_ADP;

		static const int NEQ_myofilament = 8;
		// static const int NEQ_myofilament_N = 5;
		// double atol[NEQ_myofilament+1], rtol[NEQ_myofilament+1], 
		double t_myofilament;
		double tout;
		double y_myofilament[NEQ_myofilament+1];
		// int    iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
		// int    itol, itask, istate, iopt, jt, iout;
		LSODA<Myofilament> lsoda_integrator_myofilament;

};
