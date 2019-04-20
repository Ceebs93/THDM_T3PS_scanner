#include "2HDMC/THDM.h"
#include "2HDMC/SM.h"
#include "2HDMC/HBHS.h"
#include "2HDMC/Constraints.h"
#include "2HDMC/DecayTable.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

//#define VERBOSE

using namespace std;

int main(int argc, char* argv[]) {

  int yt;
  double Z7_in;
  double mH_in;
  double mHc_in;
  double mA_in;
  double cba_in;
  double tb_in;

  if ( argc == 2)
  {

	ifstream file( argv[1]);
	//  file.open( fname.c_str(), ios::in);
  	file >> yt;
  	file >> Z7_in;
  	file >> mH_in;
  	file >> mHc_in;
  	file >> mA_in;
  	file >> cba_in;
  	file >> tb_in;
  	file.close();

  }
  else if ( argc == 8)
  {
 		 yt         = (int)   atoi(argv[1]);
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
	 printf("ParameterScan_T3PS <yt> <Z7> <mH> <mHc> <mA> <cba> <tanb>");
  }
  

  ////////////////////////////////////

  double vev = 246.2206;
  double mh_ref  = 125.;
  double Z5_c = ( mh_ref*mh_ref*cba_in*cba_in + mH_in*mH_in*(1.0-cba_in*cba_in) - mA_in*mA_in)/vev/vev;
  double Z4_c = (2.0*( (mA_in*mA_in - mHc_in*mHc_in))/vev/vev) + Z5_c;
  
  // Create SM and set parameters
  SM sm;
  sm.set_qmass_pole(6, 172.5);		
  sm.set_qmass_pole(5, 4.75);		
  sm.set_qmass_pole(4, 1.42);	
  sm.set_lmass_pole(3, 1.77684);	
  sm.set_alpha(1./127.934);
  sm.set_alpha0(1./137.0359997);
  sm.set_alpha_s(0.119);
  sm.set_MZ(91.15349);
  sm.set_MW(80.36951);
  sm.set_gamma_Z(2.49581);
  sm.set_gamma_W(2.08856);
  sm.set_GF(1.16637E-5);

  // Create 2HDM and set SM parameters
  THDM model;
  model.set_SM(sm);

  // Parameter set validation check
  bool pset = model.set_param_hybrid(mh_ref,mH_in,cba_in,Z4_c,Z5_c,Z7_in,tb_in);
  
  if (!pset) {
    cerr << "The specified parameters are not valid\n";
    return -1;
  }

  // Set Yukawa couplings
  model.set_yukawas_type(yt);

  // Print the parameters in different parametrizations to stdout
  //model.print_param_phys();
  //model.print_param_gen();
  //model.print_param_higgs();
  //model.print_param_hybrid();

  
  ////////////////////////////////////
  /// 									  ///
  /// --- Checking constraints --- ///
  /// 									  ///
  ////////////////////////////////////
  // -   !!! Important note !!!   - //
  // You need to call HB&HS initialization before constr(mode)!

  // Call HB & HS initialization
  HB_init();
  // Prepare to calculate observables
  Constraints constr(model);

  double S,T,U,V,W,X;   

  constr.oblique_param(mh_ref,S,T,U,V,W,X);

  bool BitAllowedStability      = constr.check_stability();
  bool BitAllowedUnitarity      = constr.check_unitarity();
//bool BitAllowedPerturbativity = constr.check_perturbativity(); // the default is 4.0*M_PI
  bool BitAllowedPerturbativity = constr.check_perturbativity(8.0*M_PI);

 
  /////////////////////
  // - HiggsBounds - //
  /////////////////////
  // - !!! Important note !!!  -  //
  // It seems that if you call constr.delta_amu() before
  // HB/HS then HB/HS won't work properly!

// See HiggsSignals manual for more information
  
// Share couplings of 2HDM model with HiggsBounds/HiggsSignals
  HB_set_input_effC(model);
  
  int hbres[6];
  int hbchan[6];
  double hbobs[6];
  int hbcomb[6];  

// Run HiggsBounds 'full', i.e. with each Higgs result separately  
  HB_run_full(hbres, hbchan, hbobs, hbcomb);
//  printf("\nHiggsBounds results (full):\n");
//  printf("  Higgs  res  chan       ratio        ncomb\n");
//  for (int i=1;i<=4;i++) {
//    printf("%5d %5d %6d %16.8E %5d   %s\n", i, hbres[i],hbchan[i],hbobs[i],hbcomb[i],hbobs[i]<1 ? "Allowed" : "Excluded");
//  }
//  printf("------------------------------------------------------------\n");
//  printf("  TOT %5d %6d %16.8E %5d   %s\n", hbres[0],hbchan[0],hbobs[0],hbcomb[0],hbobs[0]<1 ? "ALLOWED" : "EXCLUDED");

  double tot_hbobs = hbobs[0];
  double mostsensitivech = hbchan[0];

  double csqmu;
  double csqmh_ref;
  double csqtot;
  int nobs;
  double pval;

  //printf("\nHiggsSignals results:\n");
  //printf(" Chi^2 from rates: %16.8E\n", csqmu);
  //printf("  Chi^2 from mass: %16.8E\n", csqmh_ref);
  //printf("      Total chi^2: %16.8E\n", csqtot);
  //printf("    # observables: %16d\n\n", nobs);


//	std::complex <double> cs;
//	std::complex <double> cp;

//	model.get_coupling_vvh(2,2,1,coupling);
//  	model.get_coupling_hdd(1,3,3,cs,cp);

//  printf("Hvev:       %8.4f\n", sqrt(Hvev_2));
//  printf("mA:       %8.4f\n", mA);


  // Warning: Beware with these ones as they seem to destroy some HB/HS
  // functionality!!!
  double delta_rho = constr.delta_rho(mh_ref);
  double delta_amu = constr.delta_amu();

// constr.print_all(mh_ref);

  //////////////////////
  // - Gamma widths - //
  //////////////////////

  // Prepare to calculate decay widths
  DecayTable table(model);

  // Print to stdout total widths of Higgs bosons
  // table.print_width(1);
  // table.print_width(2);
  // table.print_width(3);
  // table.print_width(4);	

  struct BR brf;
  table.geth_BR(3,brf); // 3 - stands for the pseudoscalar A

  double br_A_tt     = brf.bruu[3][3];
  double br_A_bb     = brf.brdd[3][3];
  double br_A_gg     = brf.brhgg;
  double br_A_tautau = brf.brll[3][3];
  double br_A_zh     = brf.brvh[2][1];

//  printf("br tt:      %.3e\n", br_A_tt    );
//  printf("br bb:      %.3e\n", br_A_bb    );
//  printf("br gg:      %.3e\n", br_A_gg    );
//  printf("br tau tau: %.3e\n", br_A_tautau);
//  printf("br zh:      %.3e\n", br_A_zh    );

  // You can cross-check the branching fractions here with the full table
  // printf("br tt:      %.3e\n", brf.bruu[3][3]);
  // printf("br bb:      %.3e\n", brf.brdd[3][3]);
  // printf("br gg:      %.3e\n", brf.brhgg);
  // printf("br tau tau: %.3e\n", brf.brll[3][3]);
  // printf("br zh:      %.3e\n", brf.brvh[2][1]);
  // table.print_decays(3);

  // table.print_decays(1);
  // table.print_decays(2);
  // table.print_decays(3);
  // table.print_decays(4);



  double Hvev_2 = sm.get_v2();

  // From 2HDMC:
  //  const char *hnames[6] = {" ","h ", "H ", "A ", "H+", "H-"};
  double Gamma_h  = table.get_gammatot_h(1);
  double Gamma_H  = table.get_gammatot_h(2);
  double Gamma_A  = table.get_gammatot_h(3);
  double Gamma_Hc = table.get_gammatot_h(4);

   double mh, mH,mA,mHc,sinba,l6,l7,m12_2,tb;
   double l1,l2,l3,l4,l5;

	model.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tb);
	model.get_param_phys(mh,mH,mA,mHc,sinba,l6,l7,m12_2,tb);
	// Please note that here sin(b-a) comes from the physical basis so it can be
	// negative too!

  # ifdef VERBOSE
  double mh_hybrid,mH_hybrid,cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid;
  double mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys;
  double l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen;

  model.get_param_phys(mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys);
  model.get_param_hybrid(mh_hybrid, mH_hybrid, cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid);
  model.get_param_gen(l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen);

  printf("Inside ParameterScan_T3PS.cpp\n");
  printf("yt: %d\n", yt );

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
//  printf("sba - cba*tb:%8.4f\n", (sinba-cba*tanb));
//  printf("cs/(sba - cba*tb):%8.4f\n", std::abs(cs/(sinba-cba*tanb))/(sqrt(2.0)*4.75/246.0) );
//  std::cout << "vvh:     " << coupling << std::endl;
//  std::cout << "hdd cs:     " << cs << std::endl;
//  std::cout << "hdd cp:     " << cp << std::endl;
  # endif


  // 
  
  std::cout                                    

	// - Input
	<< Z7_in << " "                         			//
	<< mH_in << " "                         			//
	<< mHc_in   << " "                         		//
	<< mA_in    << " "                         		//
	<< cba_in << " "                        			//
	<< tb_in << " "                           		//

	// -
	<< sinba << " "                        			//
	<< Z4_c << " "                         			//
	<< Z5_c << " "                         			//
	<< m12_2 << " "                         			//

	// -- Lambdas
	<< l1 << " "                            			//
	<< l2 << " "                            			//
	<< l3 << " "                            			//
	<< l4 << " "                            			//
	<< l5 << " "                            			//
	<< l6 << " "                            			//
	<< l7 << " "                            			//
	
	// -- Widths
	<< Gamma_h  << " "                      			//
	<< Gamma_H  << " "                      			//
	<< Gamma_Hc << " "                      			//
	<< Gamma_A  << " "                      			//

   // -- BR(A->)
	<<  br_A_tt     << " "                      		//
	<<  br_A_bb     << " "                      		//
	<<  br_A_gg     << " "                      		//
	<<  br_A_tautau << " "                      		//
	<<  br_A_zh     << " "                      		//

	// Theoretical test - allowed/excluded 
	<< BitAllowedStability << " "                   //
	<< BitAllowedUnitarity << " "                   //
	<< BitAllowedPerturbativity << " "              //
	
	// Oblique
	<< S << " "                             			//
	<< T << " "                             			//
	<< U << " "                             			//
	<< V << " "                             			//
	<< W << " "                             			//
	<< X << " "                             			//
	<< delta_rho << " "                             //

	// -- g-2
	<< delta_amu << " "                             //

	// Experimental test - statistics
	<< tot_hbobs << " "                             //
   << mostsensitivech << " "                       //
	
	<< std::endl;

  HB_finish();

	return 0;

}
