/*******************************************************************************
 2HDMC - two-Higgs-doublet model calculator
 Demo program
 
 http://2hdmc.hepforge.org
*******************************************************************************/
#include "THDM.h"
#include "SM.h"
#include "HBHS.h"
#include "Constraints.h"
#include "DecayTable.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;

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

  // Set parameters of the 2HDM in the 'physical' basis
  double mh       = 125.;
  double mH       = 400.;
  double mA       = 500.;
  double mC       = 550.;
  double sba      = 0.999;
  double lambda_6 = 0.;
  double lambda_7 = 0.;
  double m12_2    = 15800.; 
  double tb       = 10.;
  
  bool pset = model.set_param_phys(mh,mH,mA,mC,sba,lambda_6,lambda_7,m12_2,tb);
  
  if (!pset) {
    cerr << "The specified parameters are not valid\n";
    return -1;
  }

#if defined HiggsBounds
  HB_init();
  HS_init();
#endif

  // Set Yukawa couplings to type II
  model.set_yukawas_type(2);
  
  // Print the parameters in different parametrizations
  model.print_param_phys();
  model.print_param_gen();
  model.print_param_higgs();
  model.print_param_hybrid();

  // Prepare to calculate observables
  Constraints constr(model);

  double S,T,U,V,W,X;   

  constr.oblique_param(mh_ref,S,T,U,V,W,X);
  
  printf("\nConstraints:\n");
  printf("  Potential stability: %s\n", 
	 (constr.check_stability() ? "OK" : "Not OK"));
  printf(" Tree-level unitarity: %s\n", 
	 (constr.check_unitarity() ? "OK" : "Not OK"));
  printf("       Perturbativity: %s\n", 
	 (constr.check_perturbativity() ? "OK" : "Not OK"));

  printf("\n");
  printf(" Oblique S: %12.5e\n", S);  
  printf(" Oblique T: %12.5e\n", T);
  printf(" Oblique U: %12.5e\n", U);  
  printf(" Oblique V: %12.5e\n", V);
  printf(" Oblique W: %12.5e\n", W);  
  printf(" Oblique X: %12.5e\n", X);
  printf(" Delta_rho: %12.5e\n", constr.delta_rho(mh_ref));
  printf("\n");
  printf(" Delta_amu: %12.5e\n\n", constr.delta_amu());

#if defined HiggsBounds

// See HiggsSignals manual for more information
  int mass_pdf = 2;
  HS_set_pdf(mass_pdf);
  HS_setup_assignment_range_massobservables(2.);
  HS_set_output_level(0);

// Share couplings of 2HDM model with HiggsBounds/HiggsSignals
  HB_set_input_effC(model);
  
  // Arrays hold the HiggsBounds results 
  int hbres[6];
  double hbobs[6];
  int hbchan[6];
  int hbcomb[6];  

// Run HiggsBounds 'full', i.e. with each Higgs result separately  
  HB_run_full(hbres, hbchan, hbobs, hbcomb);
  printf("\nHiggsBounds results (full):\n");
  printf("  Higgs  res  chan       ratio        ncomb\n");
  for (int i=1;i<=4;i++) {
    printf("%5d %5d %6d %16.8E %5d   %s\n", i, hbres[i],hbchan[i],hbobs[i],hbcomb[i],hbobs[i]<1 ? "Allowed" : "Excluded");
  }
  printf("------------------------------------------------------------\n");
  printf("  TOT %5d %6d %16.8E %5d   %s\n", hbres[0],hbchan[0],hbobs[0],hbcomb[0],hbobs[0]<1 ? "ALLOWED" : "EXCLUDED");
  
  double csqmu;
  double csqmh;
  double csqtot;
  int nobs;
  double pval;
  
  double dMh[3]={0., 0., 0.,};
  HS_set_mass_uncertainties(dMh);
 
  HS_run(&csqmu, &csqmh, &csqtot, &nobs, &pval);

  printf("\nHiggsSignals results:\n");
  printf(" Chi^2 from rates: %16.8E\n", csqmu);
  printf("  Chi^2 from mass: %16.8E\n", csqmh);
  printf("      Total chi^2: %16.8E\n", csqtot);
  printf("    # observables: %16d\n\n", nobs);

#endif  

  // Prepare to calculate decay widths
  DecayTable table(model);

  // Print total widths of Higgs bosons
  table.print_width(1);
  table.print_width(2);
  table.print_width(3);
  table.print_width(4);	

  table.print_decays(1);
  table.print_decays(2);
  table.print_decays(3);
  table.print_decays(4);

  // Write output to LesHouches file
  model.write_LesHouches("Demo_out.lha", 1, 0, 1, 1);

#if defined HiggsBounds
  HB_finish();
  HS_finish();
#endif
}

