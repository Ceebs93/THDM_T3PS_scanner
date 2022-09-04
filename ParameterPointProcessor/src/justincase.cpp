#include "2HDMC/THDM.h"                                                                            
#include "2HDMC/SM.h"                                                                              
#include "2HDMC/HBHS.h"                                                                            
#include "2HDMC/Constraints.h"                                                                     
#include "2HDMC/DecayTable.h"                                                                      
#include <iostream>                                                                                
#include <string>                                                                                  
#include <fstream>                                                                                 
#include <cmath>                                                                                   
#include "EWPO.h"                                                                                  
                                                                                                      
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
	  
	
	// -- Conversion -- //
	// -- Z4, Z5 calculated
	
	double vev    = 246.2206;
	double mh_ref = 125.09;
	double Z5_c   = ( mh_ref*mh_ref*cba_in*cba_in + mH_in*mH_in*(1.0-cba_in*cba_in) - mA_in*mA_in)/vev/vev;
	double Z4_c   = (2.0*( (mA_in*mA_in - mHc_in*mHc_in))/vev/vev) + Z5_c;

	/////////////////
	// -- 2HDMC -- //
	/////////////////
	


	const HBHSResult *hbhsres_ptr = nullptr;

	#if defined HiggsBounds
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

	  Constraints check(model);

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

	  # ifdef FAST
          if ( !(chi2_ST_hepfit < 40.0 || chi2_ST_gfitter < 40.0) )
          {
                  exit( -1 );
          }
// Error is triggered now, but the likelihood value is still changing each turn indicating that the error isn't actually the issue. Values are just being generated randomly in the python processor so it must be that the values that are being passed back from this file aren't formatted correctly or are  not varying each loop for some reason. The former makes more sense than the latter as input is changing so why would the output not be changing? Plus it seems to be impossible to avoid the output to the console from 2HDMC so this is inevitably going through the stream. Propose adding a marker variable that python can search for and parse from. This should allow for filtering out this extra output consistently and finding the needed values in consistent positions. If this is unsuccessful then perhaps a return to the idea of having an output file and reading that in, but unlike in previous attempts we use the same method to run this file from the python code but we do nothing with the output python reads back from it. This is not ideal, however as adding a unique identifier is problematic due to the way data is currently handed over to this file. We would perhaps need to add a new variable much higher up the flow diagram of this software, which might mean large changes for the structure - problematic.	
	  printf("Chi2 has a value of     ");
	  std::cout << chi2_HS << endl;

	  double chi2_Tot_hepfit, chi2_Tot_gfitter;
          chi2_Tot_gfitter = chi2_HS + chi2_ST_gfitter;
          chi2_Tot_hepfit = chi2_HS + chi2_ST_hepfit;
        
          # ifdef FAST
           if ( chi2_HS > 300.0 )
           {
                std::cerr << "chi_HS over 300!" << std::endl;
                exit( -1 );
           }
          # endif

	  hbhsres_ptr = &hbhs_result;

	  int    hbres[6];
	  int    hbchan[6];
	  double hbobs[6];
	  int    hbcomb[6];  

	
	  double tot_hbobs = hbobs[0];
	  double sens_ch  = hbchan[0];

	  for(int i=0; i<6; i++)
          {
		printf("sens_ch_%d: %f \n", i, hbchan[i]), "\n";
          }

	  double delta_rho = constr.delta_rho(mh_ref);
	  double delta_amu = constr.delta_amu();	

	  /////////////////////
	  // -- COUPLINGS -- //
	  /////////////////////

	  std::complex<double> g_HpHmh_;
	  model.get_coupling_hhh(1,4,4, g_HpHmh_);
 	  double g_HpHmh = std::abs( g_HpHmh_ );

	  /////////////////
	  // -- DECAY -- //
	  /////////////////
			
	  // -- Prepare to calculate decay widths

	  DecayTable table(model);
	  
	  // --Branching fractions
	  struct BR BRfrac_h, BRfrac_H, BRfrac_A, BRfrac_Hp;

	  table.geth_BR( 1, BRfrac_h);
	  table.geth_BR( 2, BRfrac_H );
	  table.geth_BR( 3, BRfrac_A );
	  table.geth_BR( 4, BRfrac_Hp);

	  double br_h_bb     = BRfrac_h.brdd[3][3];
	  double br_h_tautau = BRfrac_h.brll[3][3];
	  double br_h_gg     = BRfrac_h.brhgg;
	  table.geth_BR( 4, BRfrac_Hp);

	  double br_h_bb     = BRfrac_h.brdd[3][3];
	  double br_h_tautau = BRfrac_h.brll[3][3];
	  double br_h_gg     = BRfrac_h.brhgg;
	  double br_h_WW     = BRfrac_h.brvv[3];
	  double br_h_ZZ     = BRfrac_h.brvv[2];
	  double br_h_gaga   = BRfrac_h.brvv[1];

	  double br_H_tt     = BRfrac_H.bruu[3][3];
	  double br_H_bb     = BRfrac_H.brdd[3][3];
	  double br_H_gg     = BRfrac_H.brhgg;
	  double br_H_mumu   = BRfrac_H.brll[3][3];
	  double br_H_tautau = BRfrac_H.brll[3][3];
	  double br_H_Zga    = BRfrac_H.brhZga;
	  double br_H_Zh     = BRfrac_H.brvh[2][1];
	  double br_H_WW     = BRfrac_H.brvv[3];
	  double br_H_ZZ     = BRfrac_H.brvv[2];
	  double br_H_ZA     = BRfrac_H.brvh[2][3];
	  double br_H_AA     = BRfrac_H.brhh[3];
	  double br_H_hh     = BRfrac_H.brhh[1];
	  double br_H_gaga   = BRfrac_H.brvv[1];

	  double br_A_tt     = BRfrac_A.bruu[3][3];
	  double br_A_bb     = BRfrac_A.brdd[3][3];
	  double br_A_gg     = BRfrac_A.brhgg;
	  double br_A_mumu   = BRfrac_A.brll[2][2];
	  double br_A_tautau = BRfrac_A.brll[3][3];
	  double br_A_Zga    = BRfrac_A.brhZga;
	  double br_A_Zh     = BRfrac_A.brvh[2][1];
	  double br_A_ZH     = BRfrac_A.brvh[2][2];
	  double br_A_gaga   = BRfrac_A.brvv[1];

	  double br_Hp_Wh    = BRfrac_Hp.brvh[3][1];
	  double br_Hp_WH    = BRfrac_Hp.brvh[3][2];
	  double br_Hp_WA    = BRfrac_Hp.brvh[3][3];
	  double br_Hp_tb    = BRfrac_Hp.brdu[3][3];
	  double br_Hp_taunu = BRfrac_Hp.brln[3][3];

	  // -- From 2HDMC:
	  const char *hnames[6] = {" ","h ", "H ", "A ", "H+", "H-"};
	  double Gamma_h  = table.get_gammatot_h(1);
	  double Gamma_H  = table.get_gammatot_h(2);
	  double Gamma_A  = table.get_gammatot_h(3);
	  double Gamma_Hc = table.get_gammatot_h(4);

	  double mh,mH,mA,mHc,sinba,m12_2,tb;
	  double l1,l2,l3,l4,l5,l6,l7;

	  model.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tb);
	  model.get_param_phys(mh,mH,mA,mHc,sinba,l6,l7,m12_2,tb);

	  //printf("m12_2 has a value of     ");
          //std::cout << m12_2 << endl;


	  // Please note that here sin(b-a) comes from the physical basis 
	  // so it can be negative too!

	  double k_huu, k_hdd;
	  k_huu = abs(sinba) + cba_in/tb;
	  k_hdd = abs(sinba) - cba_in*tb;

	# endif

	//printf("                             ");
	//printf("I'm printing out the stuff in the DEBUG section now!");
	//printf("                              ");
	double mh_hybrid,mH_hybrid,cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid;
	double mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys;
	double l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen;
	
	model.get_param_phys(mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys);
	model.get_param_hybrid(mh_hybrid, mH_hybrid, cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid);
	model.get_param_gen(l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen);
	
	//printf("Inside ParameterScan_T3PS.cpp\n");
	//printf("yt_in: %d\n", yt_in );
	
	printf("\nComparison of variables\n");
//	printf("| cba   | %+8.3e |            |            | %+8.3e |            |\n", cba_in, cba_hybrid );
//	printf("| sba   |            |            | %+8.3e |            |            |\n", sba_phys);
//	printf("| tb    | %+8.3e | %+8.3e | %+8.3e | %+8.3e |            |\n", tb_in, tb_gen, tb_phys, tb_hybrid );
//	printf("| m12_2 |            | %+8.3e | %+8.3e |            |            |\n", m12_2_gen, m12_2_phys );
//	printf("|-------------------------------------------------------------------------\n" );
	
	printf("S:           %8.4f\n", S);
	printf("T:           %8.4f\n", T);
	printf("U:           %8.4f\n", U);
	printf("delta_rho:   %8.4f\n", delta_rho);
	printf("\nConstraints:\n");
	printf("  Potential stability: %s\n", 
	  (constr.check_stability() ? "OK" : "Not OK"));
	printf(" Tree-level unitarity: %s\n", 
	  (constr.check_unitarity() ? "OK" : "Not OK"));
	printf("       Perturbativity: %s\n", 
	  (constr.check_perturbativity() ? "OK" : "Not OK")); 

     // #if defined HiggsBounds
//	HBHS hbhs{};
      #endif

	std:cout
	
	// -- Input
	<< Z7_in << " "			//1
	<< mH_in << " "			//2
	<< mHc_in << " "		//3
        << mA_in  << " "                // 4
	<< cba_in << " "                // 5
	<< tb_in  << " "                // 6
	     	  
	<< std::endl;

	return 0;

}
