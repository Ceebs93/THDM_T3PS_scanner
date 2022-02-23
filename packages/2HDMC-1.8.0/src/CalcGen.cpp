#include "THDM.h"
#include "Constraints.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 12) {
    cout << "Usage: ./CalcGen lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6 lambda_7 m_12^2 tan_beta yukawas_type output_filename\n";
    return -1;
  }

  double l1_in = (double)atof(argv[1]);
  double l2_in = (double)atof(argv[2]);
  double l3_in = (double)atof(argv[3]);
  double l4_in = (double)atof(argv[4]);
  double l5_in = (double)atof(argv[5]);
  double l6_in = (double)atof(argv[6]);
  double l7_in = (double)atof(argv[7]);
  double m12_2_in = (double)atof(argv[8]);
  double tb_in = (double)atof(argv[9]);
  int yt_in = (int)atof(argv[10]);
  char* file = argv[11];

  THDM model;
  SM sm;

  bool pset = model.set_param_gen(l1_in,l2_in,l3_in,l4_in,l5_in,l6_in,l7_in,m12_2_in,tb_in);

  if (!pset) {
    cerr << "The parameters you have specified were not valid\n";
    return -1;
  }

  model.set_yukawas_type(yt_in);

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  Constraints check(model);

  model.print_param_phys();
  model.print_param_gen();

  check.print_all(mh_ref);


   model.write_LesHouches(file,true,true,true);

}
