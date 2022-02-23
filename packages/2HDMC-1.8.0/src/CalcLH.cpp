#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include <iostream>

using namespace std;

void pc(const char* text, complex <double> c1, complex <double> c2);


int main(int argc, char* argv[]) {

  if (argc < 3) {
    cout << "Usage: ./CalcLH input_file output_file\n";
    return -1;
  }

  char* in_file = argv[1];
  char* out_file = argv[2];

  THDM model;

  bool pset = model.read_LesHouches(in_file);

  if (!pset) {
    cerr << "Parameters have not been set from the LesHouches file\n";
    return -1;
  }

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;

  Constraints check(model);

  model.print_param_gen();
  model.print_param_higgs();
  model.print_param_phys();

  check.print_all(mh_ref);

  complex <double> cs,cp;
  model.get_coupling_hll(1,3,3,cs,cp);
  pc("h-tautau",cs,cp);
  model.get_coupling_hll(2,3,3,cs,cp);
  pc("H-tautau",cs,cp);
  model.get_coupling_hll(3,3,3,cs,cp);
  pc("A-tautau",cs,cp);
  model.write_LesHouches(out_file,true,true,true);

  DecayTable table(model);
  table.print_decays(1);
  table.print_decays(2);
  table.print_decays(3);
  table.print_decays(4);


}


void pc(const char* text, complex <double> c1, complex <double> c2) {
  const char* s = "S: ";
  const char* p = "P: ";

  printf("%5s%10s(%7.4f)+i(%7.4f)%10s(%7.4f)+i(%7.4f)\n", text, s,real(c1), imag(c1), p,real(c2), imag(c2));
}
