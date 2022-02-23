#include "Constraints.h"
#include "HBHS.h"
#include "THDM.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 8) {
    cout << "Usage: ./CalcInert mh(SM-like) mH(inert) mA(inert) mH+(inert) "
            "Lambda_2 Lambda_3 output_filename\n";
    return -1;
  }

  double mh_in = (double)atof(argv[1]);
  double mH_in = (double)atof(argv[2]);
  double mA_in = (double)atof(argv[3]);
  double mHp_in = (double)atof(argv[4]);
  double l2_in = (double)atof(argv[5]);
  double l3_in = (double)atof(argv[6]);
  char *file = argv[7];

  if ((mh_in <= 0) || (mH_in <= 0) || (mA_in <= 0) || (mHp_in <= 0)) {
    cout << "ERROR: All mass parameters must be positive\n";
    return -1;
  }

  THDM model;
  SM sm;

  bool pset = model.set_inert(mh_in, mH_in, mA_in, mHp_in, l2_in, l3_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  Constraints check(model);

  model.print_param_phys();
  model.print_param_gen();
  model.print_param_higgs();

  check.print_all(mh_ref);

  const HBHSResult *hbhsres_ptr = nullptr;
#if defined HiggsBounds
  HBHS hbhs{};

  const HBHSResult hbhs_result = hbhs.check(model);
  hbhs_result.hb.print();
  hbhs_result.hs.print();
  hbhsres_ptr = &hbhs_result;
#endif

  model.write_LesHouches(file, true, true, true, hbhsres_ptr);
}
