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

//#define VERBOSE

using namespace std;

int main(int argc, char* argv[])
{

	/////////////////
	//             //
	// -- INPUT -- //
	//             //
	/////////////////

	int       yt_in=2;  // - Type
	double    Z7_in=0;  // - Z7
	double    mH_in=500;  // - mH
	double   mHc_in=650;  // - mHc
	double    mA_in=400;  // - mA
	double   cba_in=0.15;  // - cos(b-a)
	double    tb_in=7;  // - tan(b)
//	const char file = "tester_ciara.txt"
	
	ifstream file( argv[1]);
	//  file.open( fname.c_str(), ios::in);
	file >> yt_in;
	file >> Z7_in;
	file >> mH_in;
	file >> mHc_in;
	file >> mA_in;
	file >> cba_in;
	file >> tb_in;
	file.close();

	if ( argc == 8 )   // - parameters as input
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

	// -- Create 2HDM and set SM parameters
	THDM model;
	model.set_SM(sm);
	
	// -- Parameter set validation check
	bool pset = model.set_param_hybrid(mh_ref,mH_in,cba_in,Z4_c,Z5_c,Z7_in,tb_in);
	
	if (!pset)
	{
	  std::cerr << "The specified parameters are not valid\n";
	  return -1;
	}
	
	// -- Set Yukawa couplings
	model.set_yukawas_type( yt_in );
	
	///////////////////////////
	/// 			///
	/// --- CONSTRAINTS --- ///
	/// 			///
	///////////////////////////
	
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
	
	// -- Call HB & HS initialization
	const HBHSResult *hbhsres_ptr = nullptr;

 	 HBHS hbhs{};

  	 const HBHSResult hbhs_result = hbhs.check(model);
  	 hbhs_result.hb.print();
 	 hbhs_result.hs.print();
 	 hbhsres_ptr = &hbhs_result;
      
	//Write LesHouches-style output
 	model.write_LesHouches("test.lha", true, true, true, hbhsres_ptr);

	// - S,T,U fit
 	double chi2_ST_hepfit, chi2_ST_gfitter;
	get_chi2_ST(S,T, chi2_ST_hepfit, chi2_ST_gfitter);

  	# ifdef FAST
  	if ( !(chi2_ST_hepfit < 40.0 || chi2_ST_gfitter < 40.0) )
  	{
  		exit( -1 );
  	}
  	# endif

	////////////////////////////////////
	// -- HiggsBounds/HiggsSignals -- //
	////////////////////////////////////
	
	// - !!! Important note !!!  -  //
	// It seems that if you call constr.delta_amu() before
	// HB/HS then HB/HS won't work properly!
	
	// -- See HiggsSignals manual for more information
	
	int mass_pdf = 2;
	
	int    hbres[6];
	int    hbchan[6];
	double hbobs[6];
	int    hbcomb[6];  
	
	double tot_hbobs = hbobs[0];
	double sens_ch  = hbchan[0];

   printf("sens_ch: %f \n", sens_ch);


	double csqmu;
	double csqmh_ref;
	double chi2_HS;
	int nobs;
	double pval;
	
	double dMh[3]={0., 0., 0.,};
	//HS_set_mass_uncertainties(dMh);
	
  	# ifdef FAST
  	if ( chi2_HS > 300.0 )
  	{
  	   std::cerr << "chi_HS over 300!" << std::endl;
  		exit( -1 );
 	}
  	# endif

	double chi2_Tot_hepfit, chi2_Tot_gfitter;
	chi2_Tot_gfitter = chi2_HS + chi2_ST_gfitter;
	chi2_Tot_hepfit = chi2_HS + chi2_ST_hepfit;

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
	
	// -- From 2HDMC:
	double Gamma_h  = table.get_gammatot_h(1);
	double Gamma_H  = table.get_gammatot_h(2);
	double Gamma_A  = table.get_gammatot_h(3);
	double Gamma_Hc = table.get_gammatot_h(4);
	
	double mh,mH,mA,mHc,sinba,m12_2,tb;
	double l1,l2,l3,l4,l5,l6,l7;
	
	model.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tb);
	model.get_param_phys(mh,mH,mA,mHc,sinba,l6,l7,m12_2,tb);

	// Please note that here sin(b-a) comes from the physical basis so it can be
	// negative too!
	
	double k_huu, k_hdd;
	k_huu = abs(sinba) + cba_in/tb;
	k_hdd = abs(sinba) - cba_in*tb;


	/////////////////
	//             //
	// -- DEBUG -- //
	//             //
	/////////////////
	
	# ifdef VERBOSE
	double mh_hybrid,mH_hybrid,cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid;
	double mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys;
	double l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen;
	
	model.get_param_phys(mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys);
	model.get_param_hybrid(mh_hybrid, mH_hybrid, cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid);
	model.get_param_gen(l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen);
	
	printf("Inside ParameterScan_T3PS.cpp\n");
	printf("yt_in: %d\n", yt_in );
	
	printf("\nComparison of variables\n");
	printf("|-------------------------------------------------------------------------\n" );
	printf("| Var   | Input      | Gen        |  Phys      | Hybrid     | Calc       |\n" );
	printf("|-------------------------------------------------------------------------\n" );
	printf("| Z4    |            |            |            | %+8.3e | %+8.3e |\n", Z4_hybrid, Z4_c );
	printf("| Z5    |            |            |            | %+8.3e | %+8.3e |\n", Z5_hybrid, Z5_c );
	printf("| Z7    | %+8.3e |            |            | %+8.3e |            |\n", Z7_in, Z7_hybrid );
	printf("| mH    | %+8.3e |            | %+8.3e | %+8.3e |            |\n", mH_in, mH_phys, mH_hybrid );
	printf("| mHc   | %+8.3e |            | %+8.3e |            |            |\n", mHc_in, mHc_phys);
	printf("| mA    | %+8.3e |            | %+8.3e |            |            |\n", mA_in, mA_phys );
	printf("| cba   | %+8.3e |            |            | %+8.3e |            |\n", cba_in, cba_hybrid );
	printf("| sba   |            |            | %+8.3e |            |            |\n", sba_phys);
	printf("| tb    | %+8.3e | %+8.3e | %+8.3e | %+8.3e |            |\n", tb_in, tb_gen, tb_phys, tb_hybrid );
	printf("| m12_2 |            | %+8.3e | %+8.3e |            |            |\n", m12_2_gen, m12_2_phys );
	printf("|-------------------------------------------------------------------------\n" );
	
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

	constr.print_all(mh_ref);

	table.print_decays(1);
	table.print_decays(2);
	table.print_decays(3);
	table.print_decays(4);

	printf("Z7_in:       %.3e\n", Z7_in);
	printf("mH_in:       %.3e\n", mH_in);                 
	printf("mHc_in:      %.3e\n", mHc_in);                 
	printf("mA_in:       %.3e\n", mA_in);                 
	printf("cba_in:      %.3e\n", cba_in);                 
	printf("tb_in:       %.3e\n", tb_in);                 
	
	// -- Auxiliary
	printf("sinba:       %.3e\n", sinba);                    
	printf("Z4_c:        %.3e\n", Z4_c);                    
	printf("Z5_c:        %.3e\n", Z5_c);                    
	printf("m12_2:       %.3e\n", m12_2);                    
	
	// -- lambdas
	printf("l1:          %.3e\n", l1);                         
	printf("l2:          %.3e\n", l2);                         
	printf("l3:          %.3e\n", l3);                         
	printf("l4:          %.3e\n", l4);                         
	printf("l5:          %.3e\n", l5);                         
	printf("l6:          %.3e\n", l6);                         
	printf("l7:          %.3e\n", l7);                         

	// -- Theory
	printf("sta:         %d\n", sta);                       
	printf("uni:         %d\n", uni);                       
	printf("per_4pi:     %d\n", per_4pi);                       
	printf("per_8pi:     %d\n", per_8pi);                       

	// -- EWPO
	printf("S:           %.3e\n", S);                           
	printf("T:           %.3e\n", T);                           
	printf("U:           %.3e\n", U);                           
	printf("V:           %.3e\n", V);                           
	printf("W:           %.3e\n", W);                           
	printf("X:           %.3e\n", X);                           
	printf("delta_rho:   %.3e\n", delta_rho);           
	
	// -- (g-2)
	printf("delta_amu:   %.3e\n", delta_amu);           
	
	// -- HiggsBounds/HiggsSignals
	printf("tot_hbobs:   %.3e\n", tot_hbobs);           
	printf("sens_ch:     %f\n", sens_ch);           
	printf("chi2_HS:     %.3e\n", chi2_HS);           

	# endif
	
	//////////////////
	//              //
	// -- OUTPUT -- //
	//              //
	//////////////////
	
	std::cout                                    
	
	// -- Input
	<< Z7_in  << " "                    // 1
	<< mH_in  << " "                    // 2
	<< mHc_in << " "                    // 3
	<< mA_in  << " "                    // 4
	<< cba_in << " "                    // 5
	<< tb_in  << " "                    // 6
	
	// -- Auxiliary
	<< sinba << " "                     // 7
	<< Z4_c << " "                      // 8
	<< Z5_c << " "                      // 9
	<< m12_2 << " "                     // 10
	
	// -- lambdas
	<< l1 << " "                        // 11
	<< l2 << " "                        // 12
	<< l3 << " "                        // 13
	<< l4 << " "                        // 14
	<< l5 << " "                        // 15
	<< l6 << " "                        // 16
	<< l7 << " "                        // 17

	// -- Coupling
	<< g_HpHmh << " "                   // 18
	
	// -- Widths
	<< Gamma_h  << " "                  // 19
	<< Gamma_H  << " "                  // 20
	<< Gamma_Hc << " "                  // 21
	<< Gamma_A  << " "                  // 22

//   << br_h_bb     << " "               // 23
//   << br_h_tautau << " "               // 24
//   << br_h_gg     << " "               // 25
//   << br_h_WW     << " "               // 26
//   << br_h_ZZ     << " "               // 27
//   << br_h_gaga   << " "               // 28
	
   // -- BR(A -> XX)
//	<< br_A_tt     << " "               // 29
//	<< br_A_bb     << " "               // 30
//	<< br_A_gg     << " "               // 31
//	<< br_A_mumu   << " "               // 32
//	<< br_A_tautau << " "               // 33
//	<< br_A_Zga    << " "               // 34
//	<< br_A_Zh     << " "               // 35
//	<< br_A_ZH     << " "               // 36
//	<< br_A_gaga   << " "               // 37
           
   // -- BR(H -> XX)
//	<< br_H_tt       << " "             // 38
//	<< br_H_bb       << " "             // 39
//	<< br_H_gg       << " "             // 40
//	<< br_H_mumu     << " "             // 41
//	<< br_H_tautau   << " "             // 42
//	<< br_H_Zga      << " "             // 43
//	<< br_H_Zh       << " "             // 44
//	<< br_H_WW       << " "             // 45
//	<< br_H_ZZ       << " "             // 46
//	<< br_H_ZA       << " "             // 47
//	<< br_H_AA       << " "             // 48
//	<< br_H_hh       << " "             // 49
//	<< br_H_gaga     << " "             // 50
           
   // -- BR(H+ -> XX)
//	<< br_Hp_tb     << " "              // 51
//	<< br_Hp_taunu  << " "              // 52
//	<< br_Hp_Wh     << " "              // 53
//	<< br_Hp_WH     << " "              // 54
//	<< br_Hp_WA     << " "              // 55
	
	// -- Theory
	<< sta     << " "                   // 56
	<< uni     << " "                   // 57
	<< per_4pi << " "                   // 58
	<< per_8pi << " "                   // 59

	// -- EWPO
	<< S << " "                         // 60
	<< T << " "                         // 61
	<< U << " "                         // 62
	<< V << " "                         // 63
	<< W << " "                         // 64
	<< X << " "                         // 65
	<< delta_rho << " "                 // 66
	
	// -- (g-2)
	<< delta_amu << " "                 // 67
	
	// -- HiggsBounds
	<< tot_hbobs << " "                 // 68
	<< sens_ch   << " "                 // 69

	// -- HiggsSignals
	<< chi2_HS << " "                   // 70

	<< chi2_ST_hepfit  << " "           // 71
	<< chi2_ST_gfitter << " "           // 72

	<< chi2_Tot_hepfit  << " "          // 73
	<< chi2_Tot_gfitter << " "          // 74

	<< k_huu  << " "                    // 75
	<< k_hdd  <<                        // 76

	std::endl;
	
	return 0;

}
