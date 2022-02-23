#include "Constraints.h"
#include "DecayTable.h"
#include "HBHS.h"
#include "THDM.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 11) {
    cout << "Usage: ./CalcPhys mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 "
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

  SM sm;
  model.set_SM(sm);

  bool pset = model.set_param_phys(mh_in, mH_in, mA_in, mHp_in, sba_in, l6_in,
                                   l7_in, m12_2_in, tb_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  model.set_yukawas_type(yt_in);

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;

  // Write model information to the screen
  model.print_param_phys();
  model.print_param_gen();

  Constraints check(model);
  check.print_all(mh_ref);

  const HBHSResult *hbhsres_ptr = nullptr;
#if defined HiggsBounds
  HBHS hbhs{};

  const HBHSResult hbhs_result = hbhs.check(model);
  hbhs_result.hb.print();
  hbhs_result.hs.print();
  hbhsres_ptr = &hbhs_result;
#endif

  // Write LesHouches-style output
  model.write_LesHouches(file, true, true, true, hbhsres_ptr);

  // Print Higgs decays to the screen
  DecayTable table(model);
  table.print_decays(1);
  table.print_decays(2);
  table.print_decays(3);
  table.print_decays(4);

  // Print parameters than can be used as input for HDECAY
  // model.print_hdecay();
}
