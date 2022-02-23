#include "THDM.h"
#include "Constraints.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 8) {
    cout << "Usage: ./CalcHybrid Mh MH cos(beta-alpha) Z4 Z5 Z7 tanb yukawas_type output_filename\n";
    return -1;
  }

  double mh_in = (double)atof(argv[1]);
  double mH_in = (double)atof(argv[2]);
  double cba_in = (double)atof(argv[3]);
  double Z4_in = (double)atof(argv[4]);
  double Z5_in = (double)atof(argv[5]);
  double Z7_in = (double)atof(argv[6]);
  double tanb_in = (double)atof(argv[7]);
  int yt_in = (int)atof(argv[8]);
  char* file = argv[9];

  THDM model;
  SM sm;

  bool pset = model.set_param_hybrid(mh_in,mH_in,cba_in,Z4_in,Z5_in,Z7_in,tanb_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  model.set_yukawas_type(yt_in);

  model.print_param_phys();
  model.print_param_gen();
  model.print_param_hybrid();

  Constraints check(model);
  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  check.print_all(mh_ref);

  model.write_LesHouches(file,true,true,true);

}
