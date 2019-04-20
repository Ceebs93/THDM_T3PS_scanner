#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "HBHS.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 11) {
    cout << "Usage: ./CalcPhys mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 m_12^2 tan_beta yukawas_type output_filename\n";
    return -1;
  }

  double mh_in = (double)atof(argv[1]);
  double mH_in = (double)atof(argv[2]);
  double mA_in = (double)atof(argv[3]);
  double mHp_in = (double)atof(argv[4]);
  double sba_in = (double)atof(argv[5]);
  double l6_in = (double)atof(argv[6]);
  double l7_in = (double)atof(argv[7]);
  double m12_2_in = (double)atof(argv[8]); 
  double tb_in = (double)atof(argv[9]);
  int yt_in = (int)atof(argv[10]);
  char* file = argv[11];
  

  THDM model;

  SM sm; 
  model.set_SM(sm);

  bool pset = model.set_param_phys(mh_in,mH_in,mA_in,mHp_in,sba_in,l6_in,l7_in,m12_2_in,tb_in);

  if (!pset) {
	cerr << "The parameters you have specified were not valid\n";
	return -1;
  }

#if defined HiggsBounds
  HB_init();
  HS_init();
#endif


  model.set_yukawas_type(yt_in);

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;  

  // Write model information to the screen
  model.print_param_phys();
  model.print_param_gen();	  

  Constraints check(model);
  check.print_all(mh_ref);


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

  // Write LesHouches-style output
  model.write_LesHouches(file,true,true,true,true);

  // Print Higgs decays to the screen
  DecayTable table(model);
  table.print_decays(1);
  table.print_decays(2);
  table.print_decays(3);
  table.print_decays(4);

// Print parameters than can be used as input for HDECAY
// model.print_hdecay();

#if defined HiggsBounds
  HB_finish();
  HS_finish();
#endif


}
