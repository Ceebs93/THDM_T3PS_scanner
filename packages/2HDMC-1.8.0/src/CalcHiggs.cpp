#include "THDM.h"
#include "Constraints.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 11) {
    cout << "Usage: ./CalcHiggs Lambda_1 Lambda_2 Lambda_3 Lambda_4 Lambda_5 Lambda_6 Lambda_7 m_Hp yukawas_type output_filename\n";
    return -1;
  }

  double l1_in = (double)atof(argv[1]);
  double l2_in = (double)atof(argv[2]);
  double l3_in = (double)atof(argv[3]);
  double l4_in = (double)atof(argv[4]);
  double l5_in = (double)atof(argv[5]);
  double l6_in = (double)atof(argv[6]);
  double l7_in = (double)atof(argv[7]);
  double mHp_in = (double)atof(argv[8]);
  int yt_in = (int)atof(argv[9]);
  char* file = argv[10];

  THDM model;
  SM sm;

  bool pset = model.set_param_higgs(l1_in,l2_in,l3_in,l4_in,l5_in,l6_in,l7_in,mHp_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  model.set_yukawas_type(yt_in);

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  Constraints check(model);

  model.print_param_phys();
  model.print_param_higgs();

  check.print_all(mh_ref);

  model.write_LesHouches(file,true,true,true);


}
