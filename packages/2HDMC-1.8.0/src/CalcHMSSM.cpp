#include "THDM.h"
#include "Constraints.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 4) {
    cout << "Usage: ./CalcHMSSM Mh MA tanb output_filename\n";
    return -1;
  }

  double mh_in = (double)atof(argv[1]);
  double mA_in = (double)atof(argv[2]);
  double tanb_in = (double)atof(argv[3]);
  char* file = argv[4];

  THDM model;
  SM sm;

  bool pset = model.set_hMSSM(mh_in,mA_in,tanb_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  Constraints check(model);

  model.print_param_phys();
  model.print_param_gen();
  model.print_param_hybrid();

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  check.print_all(mh_ref);

  model.write_LesHouches(file,true,true,true);

}
