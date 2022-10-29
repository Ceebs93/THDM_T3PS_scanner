/*******************************************************************************
 2HDMC - two-Higgs-doublet model calculator
 Demo program

 http://2hdmc.hepforge.org
*******************************************************************************/
#include "Constraints.h"
#include "DecayTable.h"
#include "HBHS.h"
#include "SM.h"
#include "THDM.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;

  // Create SM and set parameters
  SM sm;
  sm.set_qmass_pole(6, 172.5);
  sm.set_qmass_pole(5, 4.75);
  sm.set_qmass_pole(4, 1.42);
  sm.set_lmass_pole(3, 1.77684);
  sm.set_alpha(1. / 127.934);
  sm.set_alpha0(1. / 137.0359997);
  sm.set_alpha_s(0.119);
  sm.set_MZ(91.15349);
  sm.set_MW(80.36951);
  sm.set_gamma_Z(2.49581);
  sm.set_gamma_W(2.08856);
  sm.set_GF(1.16637E-5);

    if (argc < 11) {
    cout << "Usage: ./check_point_userinput mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 "
            "m_12^2 tan_beta yukawas_type output_filename\n";
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
  char *file = argv[11];

  THDM model;

  bool pset = model.set_param_phys(mh_in, mH_in, mA_in, mHp_in, sba_in, l6_in,
                                   l7_in, m12_2_in, tb_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  model.set_yukawas_type(yt_in);

  // Prepare to calculate observables
  Constraints constr(model);

  double S, T, U, V, W, X;

  constr.oblique_param(mh_ref, S, T, U, V, W, X);

  const HBHSResult *hbhsres_ptr = nullptr;
  #if defined HiggsBounds
    HBHS hbhs{};
    const HBHSResult hbhs_result = hbhs.check(model);
    hbhsres_ptr = &hbhs_result;
  #endif

  // Prepare to calculate decay widths
  DecayTable table(model);

  if((constr.check_stability() == 1) & (constr.check_unitarity() == 1) & (constr.check_perturbativity() == 1)
  & (hbhs_result.hb.result[0] == 1))
  {
    cout << "POINT FOUND!!!";
    model.write_LesHouches(file, true, true, true, hbhsres_ptr);
  }
}
