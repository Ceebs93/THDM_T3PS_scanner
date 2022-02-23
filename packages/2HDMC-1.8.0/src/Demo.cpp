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

  // Create 2HDM and set SM parameters
  THDM model;
  model.set_SM(sm);

  // Set parameters of the 2HDM in the 'physical' basis
  double mh = 125.;
  double mH = 400.;
  double mA = 420.;
  double mC = 440.;
  double sba = 0.999;
  double lambda_6 = 0.;
  double lambda_7 = 0.;
  double m12_2 = 40000.;
  double tb = 3.;

  bool pset =
      model.set_param_phys(mh, mH, mA, mC, sba, lambda_6, lambda_7, m12_2, tb);

  if (!pset) {
    cerr << "The specified parameters are not valid\n";
    return -1;
  }

  // Set Yukawa couplings to type II
  model.set_yukawas_type(2);

  // Print the parameters in different parametrizations
  model.print_param_phys();
  model.print_param_gen();
  model.print_param_higgs();
  model.print_param_hybrid();

  // Prepare to calculate observables
  Constraints constr(model);

  double S, T, U, V, W, X;

  constr.oblique_param(mh_ref, S, T, U, V, W, X);

  printf("\nConstraints:\n");
  printf("  Potential stability: %s\n",
         (constr.check_stability() ? "OK" : "Not OK"));
  printf(" Tree-level unitarity: %s\n",
         (constr.check_unitarity() ? "OK" : "Not OK"));
  printf("       Perturbativity: %s\n",
         (constr.check_perturbativity() ? "OK" : "Not OK"));

  printf("\n");
  printf(" Oblique S: %12.5e\n", S);
  printf(" Oblique T: %12.5e\n", T);
  printf(" Oblique U: %12.5e\n", U);
  printf(" Oblique V: %12.5e\n", V);
  printf(" Oblique W: %12.5e\n", W);
  printf(" Oblique X: %12.5e\n", X);
  printf(" Delta_rho: %12.5e\n", constr.delta_rho(mh_ref));
  printf("\n");
  printf(" Delta_amu: %12.5e\n\n", constr.delta_amu());

  const HBHSResult *hbhsres_ptr = nullptr;
#if defined HiggsBounds
  HBHS hbhs{};

  const HBHSResult hbhs_result = hbhs.check(model);
  hbhs_result.hb.print();
  hbhs_result.hs.print();
  hbhsres_ptr = &hbhs_result;
#endif

  // Prepare to calculate decay widths
  DecayTable table(model);

  // Print total widths of Higgs bosons
  table.print_width(1);
  table.print_width(2);
  table.print_width(3);
  table.print_width(4);

  table.print_decays(1);
  table.print_decays(2);
  table.print_decays(3);
  table.print_decays(4);

  // Write output to LesHouches file
  model.write_LesHouches("Demo_out.lha", 1, 0, 1, hbhsres_ptr);
}
