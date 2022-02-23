#include "THDM.h"
#include "SM.h"
#include "HBHS.h"
#include "Constraints.h"
#include "DecayTable.h"
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

	int       yt_in;  // - Type
	double    Z7_in;  // - Z7
	double    mH_in;  // - mH
	double   mHc_in;  // - mHc
	double    mA_in;  // - mA
	double   cba_in;  // - cos(b-a)
	double    tb_in;  // - tan(b)
//	const char file = "tester_ciara.txt"
	
	if     ( argc == 2 )   // - filename as input
	{
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
	
	// -- Write model information to the screen
//	model.print_param_phys();
//	model.print_param_gen();
	Constraints check(model);
	check.print_all(mh_ref);

	///////////////////////////
	/// 			///
	/// --- CONSTRAINTS --- ///
	/// 			///
	///////////////////////////
	
	// -   !!! Important note !!!   - //
	// You need to call HB&HS initialization before constr(mode)!

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

       #if defined HiggsBounds	
 	 HBHS hbhs{};

  	 const HBHSResult hbhs_result = hbhs.check(model);
	
	 double chi2_HS = hbhs_result.hs.chisq;
	 cout << chi2_HS;
 
         // - S,T,U fit

         double chi2_ST_hepfit, chi2_ST_gfitter;
         get_chi2_ST(S,T, chi2_ST_hepfit, chi2_ST_gfitter);

         cout << chi2_ST_hepfit;
	 cout << chi2_ST_gfitter;

         double chi2_Tot_hepfit, chi2_Tot_gfitter;
         chi2_Tot_gfitter = chi2_HS + chi2_ST_gfitter;
         chi2_Tot_hepfit = chi2_HS + chi2_ST_hepfit;

         # ifdef FAST
         if ( !(chi2_ST_hepfit < 40.0 || chi2_ST_gfitter < 40.0) )
         {
                 exit( -1 );
         }
         # endif

         # ifdef FAST
         if ( chi2_HS > 300.0 )
         {
         std::cerr << "chi_HS over 300!" << std::endl;
              exit( -1 );
         }
         # endif

        // double chi2_Tot_hepfit, chi2_Tot_gfitter;
        // chi2_Tot_gfitter = chi2_HS + chi2_ST_gfitter;
        // chi2_Tot_hepfit = chi2_HS + chi2_ST_hepfit;
              
  //	 hbhs_result.hb.print();
 //	 hbhs_result.hs.print();
 	 hbhsres_ptr = &hbhs_result;
         
       #endif
      
	//Write LesHouches-style output
 	model.write_LesHouches("/scratch/cb27g11/Ciara_tst.lha", true, true, true, hbhsres_ptr);

	// - S,T,U fit
	//double chi2_ST_hepfit, chi2_ST_gfitter;
	//get_chi2_ST(S,T, chi2_ST_hepfit, chi2_ST_gfitter);

  	//# ifdef FAST
  	//if ( !(chi2_ST_hepfit < 40.0 || chi2_ST_gfitter < 40.0) )
  	//{
  	//	exit( -1 );
  	//}
  	//# endif

	////////////////////////////////////
	// -- HiggsBounds/HiggsSignals -- //
	////////////////////////////////////
	
	// - !!! Important note !!!  -  //
	// It seems that if you call constr.delta_amu() before
	// HB/HS then HB/HS won't work properly!
	
	// -- See HiggsSignals manual for more information
	
	//int mass_pdf = 2;
//	HS_set_pdf(mass_pdf);
//	HS_setup_assignment_range_massobservables(2.);
//	HS_set_output_level(0);


	// -- Share couplings of 2HDM model with HiggsBounds/HiggsSignals
//	HB_set_input_effC(model);

	
	//int    hbres[6];
	//int    hbchan[6];
	//double hbobs[6];
	//int    hbcomb[6];  
	
	// -- Run HiggsBounds 'full', i.e. with each Higgs result separately  
//	HB_run_full(hbres, hbchan, hbobs, hbcomb);

	
	//double tot_hbobs = hbobs[0];
	//double sens_ch  = hbchan[0];

   //printf("sens_ch: %f \n", sens_ch);


//	for(int i=0; i<6; i++)
//   {
//		  printf("sens_ch_%d: %f \n", i, hbchan[i]);
//   }
	
//	double csqmu;
//	double csqmh_ref;
//	double chi2_HS;
//	int nobs;
//	double pval;
	
//	double dMh[3]={0., 0., 0.,};
	//HS_set_mass_uncertainties(dMh);
	
	
	//run_HiggsSignals_full(&csqmu, &csqmh_ref, &chi2_HS, &nobs, &pval);




//  	# ifdef FAST
//  	if ( chi2_HS > 300.0 )
//  	{
//  	   std::cerr << "chi_HS over 300!" << std::endl;
//  		exit( -1 );
// 	}
//  	# endif

//	double chi2_Tot_hepfit, chi2_Tot_gfitter;
//	chi2_Tot_gfitter = chi2_HS + chi2_ST_gfitter;
//	chi2_Tot_hepfit = chi2_HS + chi2_ST_hepfit;

	
	// std::complex <double> cs;
	// std::complex <double> cp;
	// model.get_coupling_vvh(2,2,1,coupling);
	// model.get_coupling_hdd(1,3,3,cs,cp);
	
	// Warning: Beware with these ones as they seem to destroy some HB/HS
	// functionality!!!
//	double delta_rho = constr.delta_rho(mh_ref);
//	double delta_amu = constr.delta_amu();

	/////////////////////
	// -- COUPLINGS -- //
	/////////////////////
	
//	std::complex<double> g_HpHmh_;
//	model.get_coupling_hhh(1,4,4, g_HpHmh_);
//	double g_HpHmh = std::abs( g_HpHmh_ );


	/////////////////
	// -- DECAY -- //
	/////////////////
	
	// -- Prepare to calculate decay widths
//	DecayTable table(model);
	
	// const char *dnames[4] = {" ","d ", "s ", "b "};
	// const char *unames[4] = {" ","u ", "c ", "t "};
	// const char *lnames[4] = {" ","e ", "mu", "ta"};
   // const char *nnames[4] = {" ","ve", "vm", "vt"};
	// const char *hnames[6] = {" ","h ", "H ", "A ", "H+", "H-"};
	// const char *vnames[5] = {" ","ga", "Z ", "W+", "W-"};
	// param h Index of Higgs boson (1,2,3,4 = h,H,A,H+)
	
	// -- Branching fractions
//	struct BR BRfrac_h, BRfrac_H, BRfrac_A, BRfrac_Hp;

//	table.geth_BR( 1, BRfrac_h);
//	table.geth_BR( 2, BRfrac_H );
//	table.geth_BR( 3, BRfrac_A );
//	table.geth_BR( 4, BRfrac_Hp);
	
//	double br_h_bb     = BRfrac_h.brdd[3][3];
//	double br_h_tautau = BRfrac_h.brll[3][3];
//	double br_h_gg     = BRfrac_h.brhgg;
//	double br_h_WW     = BRfrac_h.brvv[3];
//	double br_h_ZZ     = BRfrac_h.brvv[2];
//	double br_h_gaga   = BRfrac_h.brvv[1];

//	double br_H_tt     = BRfrac_H.bruu[3][3];
//	double br_H_bb     = BRfrac_H.brdd[3][3];
//	double br_H_gg     = BRfrac_H.brhgg;
//	double br_H_mumu   = BRfrac_H.brll[3][3];
//	double br_H_tautau = BRfrac_H.brll[3][3];
//	double br_H_Zga    = BRfrac_H.brhZga;
//	double br_H_Zh     = BRfrac_H.brvh[2][1];
//	double br_H_WW     = BRfrac_H.brvv[3];
//	double br_H_ZZ     = BRfrac_H.brvv[2];
//	double br_H_ZA     = BRfrac_H.brvh[2][3];
//	double br_H_AA     = BRfrac_H.brhh[3];
//	double br_H_hh     = BRfrac_H.brhh[1];
//	double br_H_gaga   = BRfrac_H.brvv[1];

//	double br_A_tt     = BRfrac_A.bruu[3][3];
//	double br_A_bb     = BRfrac_A.brdd[3][3];
//	double br_A_gg     = BRfrac_A.brhgg;
//	double br_A_mumu   = BRfrac_A.brll[2][2];
//	double br_A_tautau = BRfrac_A.brll[3][3];
//	double br_A_Zga    = BRfrac_A.brhZga;
//	double br_A_Zh     = BRfrac_A.brvh[2][1];
//	double br_A_ZH     = BRfrac_A.brvh[2][2];
//	double br_A_gaga   = BRfrac_A.brvv[1];

//	double br_Hp_Wh    = BRfrac_Hp.brvh[3][1];
//	double br_Hp_WH    = BRfrac_Hp.brvh[3][2];
//	double br_Hp_WA    = BRfrac_Hp.brvh[3][3];
//	double br_Hp_tb    = BRfrac_Hp.brdu[3][3];
//	double br_Hp_taunu = BRfrac_Hp.brln[3][3];
	
	// -- From 2HDMC:
	// const char *hnames[6] = {" ","h ", "H ", "A ", "H+", "H-"};
//	double Gamma_h  = table.get_gammatot_h(1);
//	double Gamma_H  = table.get_gammatot_h(2);
//	double Gamma_A  = table.get_gammatot_h(3);
//	double Gamma_Hc = table.get_gammatot_h(4);
	
//	double mh,mH,mA,mHc,sinba,m12_2,tb;
//	double l1,l2,l3,l4,l5,l6,l7;
	
//	model.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tb);
//	model.get_param_phys(mh,mH,mA,mHc,sinba,l6,l7,m12_2,tb);

	// Please note that here sin(b-a) comes from the physical basis so it can be
	// negative too!
	
//	double k_huu, k_hdd;
//	k_huu = abs(sinba) + cba_in/tb;
//	k_hdd = abs(sinba) - cba_in*tb;


	/////////////////
	//             //
	// -- DEBUG -- //
	//             //
	/////////////////
	
	# ifdef VERBOSE
//	double mh_hybrid,mH_hybrid,cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid;
//	double mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys;
//	double l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen;
	
//	model.get_param_phys(mh_phys,mH_phys,mA_phys,mHc_phys,sba_phys,l6_phys,l7_phys,m12_2_phys,tb_phys);
//	model.get_param_hybrid(mh_hybrid, mH_hybrid, cba_hybrid, Z4_hybrid, Z5_hybrid, Z7_hybrid, tb_hybrid);
//	model.get_param_gen(l1_gen,l2_gen,l3_gen,l4_gen,l5_gen,l6_gen,l7_gen,m12_2_gen,tb_gen);
	
//	printf("Inside ParameterScan_T3PS.cpp\n");
//	printf("yt_in: %d\n", yt_in );
	
//	printf("\nComparison of variables\n");
//	printf("|-------------------------------------------------------------------------\n" );
//	printf("| Var   | Input      | Gen        |  Phys      | Hybrid     | Calc       |\n" );
//	printf("|-------------------------------------------------------------------------\n" );
//	printf("| Z4    |            |            |            | %+8.3e | %+8.3e |\n", Z4_hybrid, Z4_c );
//	printf("| Z5    |            |            |            | %+8.3e | %+8.3e |\n", Z5_hybrid, Z5_c );
//	printf("| Z7    | %+8.3e |            |            | %+8.3e |            |\n", Z7_in, Z7_hybrid );
//	printf("| mH    | %+8.3e |            | %+8.3e | %+8.3e |            |\n", mH_in, mH_phys, mH_hybrid );
//	printf("| mHc   | %+8.3e |            | %+8.3e |            |            |\n", mHc_in, mHc_phys);
//	printf("| mA    | %+8.3e |            | %+8.3e |            |            |\n", mA_in, mA_phys );
//	printf("| cba   | %+8.3e |            |            | %+8.3e |            |\n", cba_in, cba_hybrid );
//	printf("| sba   |            |            | %+8.3e |            |            |\n", sba_phys);
//	printf("| tb    | %+8.3e | %+8.3e | %+8.3e | %+8.3e |            |\n", tb_in, tb_gen, tb_phys, tb_hybrid );
//	printf("| m12_2 |            | %+8.3e | %+8.3e |            |            |\n", m12_2_gen, m12_2_phys );
//	printf("|-------------------------------------------------------------------------\n" );
//	
//	printf("S:           %8.4f\n", S);
//	printf("T:           %8.4f\n", T);
//	printf("U:           %8.4f\n", U);
//	printf("delta_rho:   %8.4f\n", delta_rho);
//	printf("\nConstraints:\n");
//	printf("  Potential stability: %s\n", 
//	  (constr.check_stability() ? "OK" : "Not OK"));
//	printf(" Tree-level unitarity: %s\n", 
//	  (constr.check_unitarity() ? "OK" : "Not OK"));
//	printf("       Perturbativity: %s\n", 
//	  (constr.check_perturbativity() ? "OK" : "Not OK"));

	// Print the parameters in different parametrizations to stdout
	//model.print_param_phys();
	//model.print_param_gen();
	//model.print_param_higgs();
	//model.print_param_hybrid();
	
	//printf("\nHiggsSignals results:\n");
	//printf(" Chi^2 from rates: %16.8E\n", csqmu);
	//printf("  Chi^2 from mass: %16.8E\n", csqmh_ref);
	//printf("      Total chi^2: %16.8E\n", csqtot);
	//printf("    # observables: %16d\n\n", nobs);
	
	//printf("\nHiggsBounds results (full):\n");
	//printf("  Higgs  res  chan       ratio        ncomb\n");
	//for (int i=1;i<=4;i++)
	//{
	//  printf("%5d %5d %6d %16.8E %5d   %s\n", i, hbres[i],hbchan[i],hbobs[i],hbcomb[i],hbobs[i]<1 ? "Allowed" : "Excluded");
	//}
	//printf("------------------------------------------------------------\n");
	//printf("  TOT %5d %6d %16.8E %5d   %s\n", hbres[0],hbchan[0],hbobs[0],hbcomb[0],hbobs[0]<1 ? "ALLOWED" : "EXCLUDED");

//	constr.print_all(mh_ref);

	//printf("br tt:      %.3e\n", br_A_tt    );
	//printf("br bb:      %.3e\n", br_A_bb    );
	//printf("br gg:      %.3e\n", br_A_gg    );
	//printf("br tau tau: %.3e\n", br_A_tautau);
	//printf("br zh:      %.3e\n", br_A_Zh    );
	//
	//You can cross-check the branching fractions here with the full table
//	printf("br tt:      %.3e\n", BRfrac_A.bruu[3][3]);
//	printf("br bb:      %.3e\n", BRfrac_A.brdd[3][3]);
//	printf("br gg:      %.3e\n", BRfrac_A.brhgg);
//	printf("br tau tau: %.3e\n", BRfrac_A.brll[3][3]);
//	printf("br zh:      %.3e\n", BRfrac_A.brvh[2][1]);
	//table.print_decays(3);
	//
//	table.print_decays(1);
//	table.print_decays(2);
//	table.print_decays(3);
//	table.print_decays(4);

//	printf("Z7_in:       %.3e\n", Z7_in);
//	printf("mH_in:       %.3e\n", mH_in);                 
//	printf("mHc_in:      %.3e\n", mHc_in);                 
//	printf("mA_in:       %.3e\n", mA_in);                 
//	printf("cba_in:      %.3e\n", cba_in);                 
//	printf("tb_in:       %.3e\n", tb_in);                 
	
	// -- Auxiliary
//	printf("sinba:       %.3e\n", sinba);                    
//	printf("Z4_c:        %.3e\n", Z4_c);                    
//	printf("Z5_c:        %.3e\n", Z5_c);                    
//	printf("m12_2:       %.3e\n", m12_2);                    
	
	// -- lambdas
//	printf("l1:          %.3e\n", l1);                         
//	printf("l2:          %.3e\n", l2);                         
//	printf("l3:          %.3e\n", l3);                         
//	printf("l4:          %.3e\n", l4);                         
//	printf("l5:          %.3e\n", l5);                         
//	printf("l6:          %.3e\n", l6);                         
//	printf("l7:          %.3e\n", l7);                         

	// -- Couplings
	//printf("g_HpHmh.r:   %.3e\n", g_HpHmh.real());             
	//printf("g_HpHmh.i:   %.3e\n", g_HpHmh.imag());             

	// -- Widths
//	printf("Gamma_h:     %.3e\n", Gamma_h);             
//	printf("Gamma_H:     %.3e\n", Gamma_H);             
//	printf("Gamma_Hc:    %.3e\n", Gamma_Hc);             
//	printf("Gamma_A:     %.3e\n", Gamma_A);             

	// -- BR(h->XX)
//	printf("br_h_bb:     %.3e\n", br_h_bb);      
//	printf("br_h_gg:     %.3e\n", br_h_gg);      
//	printf("br_h_gaga:   %.3e\n", br_h_gaga);      
//	printf("br_h_tautau: %.3e\n", br_h_tautau);      
//      printf("br_h_WW:     %.3e\n", br_h_WW);      
//	printf("br_h_ZZ:     %.3e\n", br_h_ZZ);

	
   // -- BR(A->XX)
//	printf("br_A_tt:     %.3e\n", br_A_tt);      
//	printf("br_A_bb:     %.3e\n", br_A_bb);      
//	printf("br_A_gg:     %.3e\n", br_A_gg);      
//	printf("br_A_gaga:   %.3e\n", br_A_gaga);      
//	printf("br_A_tautau: %.3e\n", br_A_tautau);      
//	printf("br_A_Zh:     %.3e\n", br_A_Zh);      

	// -- BR(H->XX)
//	printf("br_H_tt:     %.3e\n", br_H_tt);      
//	printf("br_H_bb:     %.3e\n", br_H_bb);      
//	printf("br_H_gg:     %.3e\n", br_H_gg);      
//	printf("br_H_gaga:   %.3e\n", br_H_gaga);      
//	printf("br_H_tautau: %.3e\n", br_H_tautau);      
//	printf("br_H_Zh:     %.3e\n", br_H_Zh);      
//	printf("br_H_ZA:     %.3e\n", br_H_ZA);      
//	printf("br_H_AA:     %.3e\n", br_H_AA);      
//	printf("br_H_WW:     %.3e\n", br_H_WW);      
//	printf("br_H_ZZ:     %.3e\n", br_H_ZZ);

	// -- BR(H+->XX)
//	printf("br_Hp_Wh:     %.3e\n", br_Hp_Wh);
//	printf("br_Hp_WA:     %.3e\n", br_Hp_WA);      
//	printf("br_Hp_tb:     %.3e\n", br_Hp_tb);
//	printf("br_Hp_taunu:  %.3e\n", br_Hp_taunu);


	// -- Theory
//	printf("sta:         %d\n", sta);                       
//	printf("uni:         %d\n", uni);                       
//	printf("per_4pi:     %d\n", per_4pi);                       
//	printf("per_8pi:     %d\n", per_8pi);                       

	// -- EWPO
//	printf("S:           %.3e\n", S);                           
//	printf("T:           %.3e\n", T);                           
//	printf("U:           %.3e\n", U);                           
//	printf("V:           %.3e\n", V);                           
//	printf("W:           %.3e\n", W);                           
//	printf("X:           %.3e\n", X);                           
//	printf("delta_rho:   %.3e\n", delta_rho);           
	
	// -- (g-2)
//	printf("delta_amu:   %.3e\n", delta_amu);           
	
	// -- HiggsBounds/HiggsSignals
//	printf("tot_hbobs:   %.3e\n", tot_hbobs);           
//	printf("sens_ch:     %f\n", sens_ch);           
//	printf("chi2_HS:     %.3e\n", chi2_HS);           



	# endif
	
	//////////////////
	//              //
	// -- OUTPUT -- //
	//              //
	//////////////////
	
//	double br_H_sum =	br_H_tt    
//	+ br_H_bb    
//	+ br_H_gg    
//	+ br_H_mumu  
//	+ br_H_tautau
//	+ br_H_Zga   
//	+ br_H_Zh    
//	+ br_H_WW    
//	+ br_H_ZZ    
//	+ br_H_ZA    
//	+ br_H_AA    
//	+ br_H_hh
//	+ br_H_gaga;

//	double br_A_sum = 	br_A_tt     
//	+ br_A_bb     
//	+ br_A_gg     
//	+ br_A_mumu   
//	+ br_A_tautau 
//	+ br_A_Zga    
//	+ br_A_Zh     
//	+ br_A_ZH     
//	+ br_A_gaga   ;

//	std::cout                                    
	
	// -- Input
//	<< Z7_in  << " "                    // 1
//	<< mH_in  << " "                    // 2
//	<< mHc_in << " "                    // 3
//	<< mA_in  << " "                    // 4
//	<< cba_in << " "                    // 5
//	<< tb_in  << " "                    // 6
	
	// -- Auxiliary
//	<< sinba << " "                     // 7
//	<< Z4_c << " "                      // 8
//	<< Z5_c << " "                      // 9
//	<< m12_2 << " "                     // 10
	
	// -- lambdas
//	<< l1 << " "                        // 11
//	<< l2 << " "                        // 12
//	<< l3 << " "                        // 13
//	<< l4 << " "                        // 14
//	<< l5 << " "                        // 15
//	<< l6 << " "                        // 16
//	<< l7 << " "                        // 17

	// -- Coupling
//	<< g_HpHmh << " "                   // 18
	
	// -- Widths
//	<< Gamma_h  << " "                  // 19
//	<< Gamma_H  << " "                  // 20
//	<< Gamma_Hc << " "                  // 21
//	<< Gamma_A  << " "                  // 22

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
//	<< sta     << " "                   // 56
//	<< uni     << " "                   // 57
//	<< per_4pi << " "                   // 58
//	<< per_8pi << " "                   // 59

	// -- EWPO
//	<< S << " "                         // 60
//	<< T << " "                         // 61
//	<< U << " "                         // 62
//	<< V << " "                         // 63
//	<< W << " "                         // 64
//	<< X << " "                         // 65
//	<< delta_rho << " "                 // 66
	
	// -- (g-2)
//	<< delta_amu << " "                 // 67
	
	// -- HiggsBounds
//	<< tot_hbobs << " "                 // 68
//	<< sens_ch   << " "                 // 69

	// -- HiggsSignals
//	<< chi2_HS << " "                   // 70

//	<< chi2_ST_hepfit  << " "           // 71
//	<< chi2_ST_gfitter << " "           // 72

//	<< chi2_Tot_hepfit  << " "          // 73
//	<< chi2_Tot_gfitter << " "          // 74

//	<< k_huu  << " "                    // 75
//	<< k_hdd  <<                        // 76

//	<< br_H_sum  << " "                 // tmp
//	<< br_A_sum  <<                     // tmp

//	std::endl;
	
	//~HBHS();
	
	return 0;

}
