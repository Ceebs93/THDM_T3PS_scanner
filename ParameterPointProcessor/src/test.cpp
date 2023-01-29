#include "../../packages/HiggsSignals-2.6.2/include/HiggsSignals.h"
#include "2HDMC/THDM.h"                                                                            
#include "2HDMC/SM.h"                                                                              
#include "2HDMC/HBHS.h"                                                                            
#include "2HDMC/Constraints.h"                                                                     
#include "2HDMC/DecayTable.h"                                                                      
#include "2HDMC/HBHS.h"
#include <iostream>                                                                                
#include <string>                                                                                  
#include <fstream>                                                                                 
#include <cmath>
#include "EWPO.h"                                                                                  
                                                                                                     
using namespace std;                                                                               
                                                                                                     
int main(int argc, char* argv[])                                                                   
{
	int       yt_in;  // - Type
	double    Z7_in;  // - Z7
	double    mH_in;  // - mH
	double   mHc_in;  // - mHc
	double    mA_in;  // - mA
	double   cba_in;  // - cos(b-a)
	double    tb_in;  // - tan(b)
	
	if     ( argc == 2 )   // - filename as input
	{
		ifstream file( argv[1]);
		file >> yt_in;
		file >> Z7_in;
		file >> mH_in;
		file >> mHc_in;
		file >> mA_in;
		file >> cba_in;
		file >> tb_in;
		file.close();
	}

	else if ( argc == 8 )   // - parameters as input
	{
	  	yt_in      = (int)   atoi(argv[1]);
		Z7_in      = (double)atof(argv[2]);
		mH_in      = (double)atof(argv[3]);
		mHc_in     = (double)atof(argv[4]);
	  	mA_in      = (double)atof(argv[5]);
	   	cba_in     = (double)atof(argv[6]);
		tb_in      = (double)atof(argv[7]);
	}
	
	else
	{
	  printf("ParameterScan_T3PS usage:\n");
	  printf("ParameterScan_T3PS <yt_in> <Z7> <mH> <mHc> <mA> <cba> <tanb>");
	}

	int nHiggsneut = 3;
	int nHiggsplus = 1;

	double vev    = 246.2206;
	double mh_ref = 125.09;
	double Z5_c   = ( mh_ref*mh_ref*cba_in*cba_in + mH_in*mH_in*(1.0-cba_in*cba_in) - mA_in*mA_in)/vev/vev;
	double Z4_c   = (2.0*( (mA_in*mA_in - mHc_in*mHc_in))/vev/vev) + Z5_c;

	const HBHSResult *hbhsres_ptr = nullptr;

	  HBHS hbhs{};

	  // -- Create SM and set parameters
	  SM sm;

	  // - Matching HiggsBounds
  	  sm.set_qmass_pole(6, 172.6);		
  	  sm.set_qmass_pole(5, 4.60);		
  	  sm.set_qmass_pole(4, 1.40);	
  	  sm.set_lmass_pole(3, 1.7771);	
  	  sm.set_alpha(1./127.934);
  	  sm.set_alpha0(1./137.0359895);
  	  sm.set_alpha_s(0.118);
  	  sm.set_MZ(91.187);
  	  sm.set_MW(80.41);
  	  sm.set_gamma_Z(2.490);
  	  sm.set_gamma_W(2.080);
  	  sm.set_GF(1.16639E-5);

	  //Create 2HDM and set SM parameters
	  THDM model;
	  model.set_SM(sm);

	  bool pset = model.set_param_hybrid(mh_ref,mH_in,cba_in,Z4_c,Z5_c,Z7_in,tb_in);

 	  if (!pset) {
   	    cerr << "The specified parameters are not valid\n";
   	    return -1;
 	  }

 	  // Set Yukawa couplings
	  model.set_yukawas_type(yt_in);

	  // -- Prepare to calculate observables
	  Constraints constr(model);

	  // -- EWPO
	  double S,T,U,V,W,X;
          constr.oblique_param(mh_ref,S,T,U,V,W,X);

	  // -- Theory
	  bool sta      = constr.check_stability();
          bool uni      = constr.check_unitarity();
          bool per_4pi = constr.check_perturbativity(4.0*M_PI);
          bool per_8pi = constr.check_perturbativity(8.0*M_PI);

	  const HBHSResult hbhs_result = hbhs.check(model);
	  double chi2_HS = hbhs_result.hs.chisq;

	  // - S,T,U fit
	  double chi2_ST_hepfit, chi2_ST_gfitter;
	  get_chi2_ST(S,T, chi2_ST_hepfit, chi2_ST_gfitter);

//	initialize_HiggsSignals(nHiggsneut, nHiggsplus, "latestresults");

//	void run_HiggsSignals(double *Chisq_mu, double *Chisq_mh, double *Chisq,
 //                          int *nobs, double *Pvalue);

	return 0;

}
