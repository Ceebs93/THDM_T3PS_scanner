#include "../../packages/HiggsSignals-2.6.2/include/HiggsSignals.h"
#include "2HDMC/HBHS.h"
#include <iostream>                                                                                
#include <string>                                                                                  
#include <fstream>                                                                                 
#include <cmath>                                                                                   
// include "EWPO.h"                                                                                  
                                                                                                      
#define VERBOSE                                                                                    
                                                                                                     
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

	const HBHSResult *hbhsres_ptr = nullptr;

	  HBHS hbhs{};

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



}
