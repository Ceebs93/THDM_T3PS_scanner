#include "THDM.h"
#include "Constraints.h"
#include <iostream>
#include "DecayTable.h"

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 4) {
    cout << "Usage: CalcMSSM.x mA tan_beta output_filename\n";
    return -1;
  }

  double mA_in = (double)atof(argv[1]);
  double tb_in = (double)atof(argv[2]);
  char* file = argv[3];

  THDM model;
  SM sm;

  bool pset = model.set_MSSM(mA_in,tb_in);

  if (!pset) {
    cerr << "The values given for one or more parameters were invalid\n";
    return -1;
  }

  model.print_param_phys();
  model.print_param_gen();
  model.print_param_hybrid();

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  Constraints check(model);
  check.print_all(mh_ref);

  model.write_LesHouches(file,true,true,true);


}
