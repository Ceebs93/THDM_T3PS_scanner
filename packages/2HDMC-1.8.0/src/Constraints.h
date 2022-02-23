#if !defined(CONSTRAINTS_H)
#define CONSTRAINTS_H

#include "SM.h"
#include "THDM.h"
#include "DecayTable.h"

using namespace std;

/** @brief Calculates observables useful for constraining the 2HDM.
*
*   This class implements various physical observables that can be used to
*   constrain the underlying Two-Higgs doublet model. The 2HDM is specified by a
*   THDM object. Methods are provided to access the NMSSMTools library of search
*   channels for Higgs bosons at LEP (option).
*/
class Constraints {

 public:

  /**
  * @brief Default constructor without model specification
  *
  * The default constructor initializes the Constraints object, but the model
  * on which to operate has to be specified at a later stage through the
  * set_model method.
  *
  */
  Constraints();


  /**
  * @brief Constructor with model specification
  *
  * This constructor initializes the Constraints object to calculate constraints
  * on a general 2HDM. The model on which to operate is specified as an
  * argument to the constructor.
  *
  * @param  mod THDM object for which the constraints should be calculated
  */
  Constraints(THDM mod);


  /**
  * @brief Sets the model for which to perform the calculations
  *
  * This method sets the underlying THDM object on which the other methods
  * defined in this class operate when performing calculations.
  *
  * @param  mod THDM object to use for the underlying model
  */
  void set_THDM(THDM mod);


  /**
  * @brief Checks tree-level unitarity constraints
  *
  * This method checks whether the parameters of the Higgs potential results
  * in an S-matrix for Higgs-Higgs scattering that fullfills tree-level unitarity,
  * as discussed in <a href="http://xxx.lanl.gov/abs/hep-ph/0508020">hep-ph/0508020</a>.
  *
  * @param unitarity_limit optional value to be used for unitarity limit (default = \f$ 16\pi \f$)
  *
  * @returns Boolean \a true if the unitarity constraint is satisfied, \a false otherwise
  */
  bool check_unitarity(double unitarity_limit = 16*M_PI);

  /**
  * @brief Checks perturbativity
  *
  * This method checks wether the couplings of the Higgs and Gauge bosons fullfills perturbativity
  *
  * @param  perturbativity_limit optional value to be used for perturbativity limit (default = \f$ 4\pi \f$)
  *
  * @returns Boolean \a true if the perturbativity constraint is satisfied, \a false otherwise
  */
  bool check_perturbativity(double perturbativity_limit = 4*M_PI);

  /**
  * @brief Checks the stability of the Higgs potential
  *
  * This method checks whether the Higgs potential is stable, i.e. if no
  * directions exist in field space for which \f$ V\to -\infty \f$ for large
  * values of the fields. This is done through a combination of analytical and numerical
  * algorithms, depending on the generality of the specified potential.
  *
  * @returns Boolean \a true if the potential is deemed stable, \a false otherwise
  */
  bool check_stability();

  /**
  * @brief Checks the stability of the Higgs potential
  *
  * This method does the exact same this as check_stability()
  *
	* @see check_stability
	*
  * @returns Boolean \a true if the potential is deemed stable, \a false otherwise
  */
  bool check_positivity();

  /**
  * @brief Checks the 2HDM against the combined mass constraints from LEP and Tevatron
  *
  * This method checks the charged Higgs mass of the 2HDM against the built in
  * implementation routine and also checks the neutral Higgs masses of the
  * 2HDM against the existing constraints from LEP run II using the NMSSMTools
  * library and/or against the existing constraints from LEP run II and
  * Tevatron using the HiggsBounds library.
  * To use one of these library requires linking with NMSSMTools and/or
  * HiggsBounds, see the full <a href="http://arxiv.org/abs/0902.0851">manual</a> for details
  * on how to do this.
  *
  * @returns  Boolean \a true if the 2HDM is compliant with ALL CHECKED mass limits from
  *           LEP and Tevatron, \a false otherwise. In case neither NMSSMTools
  *           nor HiggsBounds is linked to 2HMDC, the method only checks the
  *           internal charged Higgs limits.
  *
  * @see      check_charged, check_NMSSMTools, check_HiggsBounds
  */
  bool check_masses();

  /**
  * @brief Checks the 2HDM against the combined mass constraints from LEP and Tevatron
  *
  * @deprecated   Use check_masses().
  */
  bool check_lep();


  /**
  * @brief Checks the 2HDM against the charged Higgs mass constraints from LEP
  *
  * This method checks the charged Higgs mass of the 2HDM against the existing
  * constraints from DELPHI.
  *
  * The result is given as a series of boolean variables, one for each search
  * channel, with value \a true if the model is consistent with the DELPHI
  * results, and \a false otherwise. The logical AND operation applied to
  * these variables is returned as output from the method.
  *
  * @param    HpHp     Channel \f$ Z \to H^+ H^- \f$, model independent
  * @param    HpHptau  Channel \f$ H^+H^- \to \tau^+ \nu_\tau \tau^- \bar{\nu}_\tau \f$
  * @param    HpHpcs   Channel \f$ H^+H^- \to c \bar{s} \bar{c} s \f$
  *
  * @returns  Boolean \a true if the model is compliant with the constraints
  *
  * @see      check_masses
  *           from every channel
  */
  bool check_charged(bool &HpHp, bool &HpHptau, bool &HpHpcs);

  /**
  * @brief Checks the 2HDM against mass constraints from LEP using NMSSMTools
  *
  * This method checks the neutral Higgs masses of the 2HDM against the
  * existing constraints from LEP run II, channel by channel, using the
  * NMSSMTools library.
  * To use this library requires linking with NMSSMTools, see the full
  * <a href="http://arxiv.org/abs/0902.0851">manual</a> for details on how to do this.
  *
  * The result is given as a series of boolean variables, one for each search
  * channel, with value \a true if the model is consistent with the LEP results,
  *  and \a false otherwise. The logical AND operation applied to these variables
  * is returned as output from the method.
  *
  * @param    hZ  Channel \f$ Z \to Zh \f$, independent of decay mode for h
  *
  * @returns  Boolean \a true if the model is compliant with the constraints
  *           from every channel
  *
  * @see      check_masses
  */
  bool check_NMSSMTools(bool &hZ, bool &hZ2b, bool &hZ2tau, bool &hZinv, bool &hZ2j,
		 bool &hZ2gamma, bool &hZ4b, bool &hZ4tau, bool &hZ2b2tau,
		 bool &hA, bool &hA4b, bool &hA4tau, bool &hA2b2tau,
		 bool &hA6b, bool &hA6tau, bool &ZhZjj);

  /**
  * @brief Checks the 2HDM against mass constraints from LEP using HiggsBounds
  *
  * This method checks the neutral Higgs masses of the 2HDM against the
  * existing constraints from LEP run II and Tevatron using the
  * HiggsBounds library.
  * To use this library requires linking with HiggsBounds, see the full
  * <a href="http://arxiv.org/abs/0902.0851">manual</a> for details on how to do this.
  *
  *
  * @param    HBresult   \a 0 if the model is excluded, \a 1 if the model is
  *                      not excluded, \a -1 invalid model
  * @param    chan       Channel with highest statistical sensitivity, see Key.dat
  * @param    obsratio   Ratio of the model rate to the observed limit for
  *                      this process
  * @param    ncombined  Number of Higgs bosons which contributed to the model
  *                      rate
  *
  * @returns  Boolean \a true if the model is compliant with the constraints
  *
  * @see      check_masses
  */
  bool check_HiggsBounds(int &HBresult, int &chan, double &obsratio,
			 int &ncombined);


  /**
  * @brief Calculates the anomalous magnetic moment of the muon
  *
  * This method calculates the anomalous magnetic moment of the muon,
  * \f$ \delta a_\mu=(g-2)_\mu \f$, in the 2HDM. The calculation includes both
  * the one-loop contributions and the (often dominant) two-loop contributions
  * from Barr-Zee type diagrams with a neutral Higgs and a photon.
  *
  * @returns Value for \f$ \delta a_\mu \f$
  */
  double delta_amu();


  /**
  * @brief Calculates \f$ \Delta\rho \f$
  *
  * This method calculates the \f$ \Delta\rho \f$ relation in the 2HDM. It is
  * related through the oblique parameter T through \f$ \Delta\rho = \alpha T
  * \f$, where \f$ \alpha \f$ is the EM coupling.
  *
  * @param    mh Mass of the %SM Higgs, for which the contribution to the
  *              oblique parameters is subtracted.
  * @returns  Value for \f$ \Delta\rho \f$
  *
  * @see   oblique_param
  */
  double delta_rho(double mh);


  /**
  * @brief Calculates oblique EW parameters S,T,U,V,W,X
  *
  * This method calculates the oblique EW parameters (S,T,U,V,W,X) in the 2HDM.
  * The %SM contribution is subtracted. Conventions from <a
  * href="http://xxx.lanl.gov/abs/0802.4353">arXiv:0802.4353</a>.
  *
  * @param mh Mass of the %SM Higgs, for which the contribution to the oblique
  *           parameters is subtracted.
  *
  * @param S  Result for S parameter (set by the method)
  * @param T  Result for T parameter (set by the method)
  * @param U  Result for U parameter (set by the method)
  * @param V  Result for V parameter (set by the method)
  * @param W  Result for W parameter (set by the method)
  * @param X  Result for X parameter (set by the method)
  */
  void oblique_param(double mh,double &S, double &T, double &U, double &V, double &W, double &X);

  /**
  * @brief Prints all constraints and observables
  *
  * This method prints the results from checking all constraints and the value
  * of all observables defined in the class.
  *
  * @param mh_ref Mass of the %SM Higgs, for which the contribution to the
  *               oblique parameters is subtracted.
  */
  void print_all(double mh_ref);

 private:
  THDM model;
  SM sm;
  DecayTable table;

  double dmu_f(double z);
  double dmu_g(double z);
  double dmu_L(double z, int h);

  double Fdrho(double x, double y);
  double G_fcn(double x, double y, double q);
  double Gtilde_fcn(double x, double y, double q);
  double Ghat_fcn(double x, double q);
  double f_fcn(double t, double r);
  double H_fcn(double x, double y, double q);
  double Htilde_fcn(double x, double y, double q);
  double Hhat_fcn(double x, double q);

  void init_Hp();
	void init_externals();


  double *mHp1, *mHp2;
  int 		nHp1, nHp2;
  double *valHp1, *valHp2;

  double delta;

  constexpr static double Z_LIMIT_MCH = 39.6;

};

// Structs used by NMSSMTools Fortran subroutines and common blocks
extern "C"
{
  extern void initialize_();
  extern void subexp_(double par[24], double prob[37]);

  extern struct {
    double alsmz;
    double alemz;
    double gf;
    double g1;
    double g2;
    double s2tw;
  } gauge_;


  extern struct {
    double ms;
    double mc;
    double mb;
    double mbp;
    double mt;
    double mtau;
    double mmuon;
    double mz;
    double mw;
  } smspec_;

  extern struct {
    double mgl;
    double mch[2];
    double u[2][2];
    double v[2][2];
    double mneu[5];
    double neu[5][5];
  } susyspec_;

  extern struct {
    double smass[3];
    double scomp[3][3];
    double pmass[2];
    double pcomp[2][2];
    double cmass;
  } higgspec_;

  extern struct {
    double mur;
    double mul;
    double mdr;
    double mdl;
    double mlr;
    double mll;
    double mnl;
    double mst1;
    double mst2;
    double msb1;
    double msb2;
    double msl1;
    double msl2;
    double msnt;
    double cst;
    double csb;
    double csl;
    double msmu1;
    double msmu2;
    double msmunt;
    double csmu;
  } sfspec_;

  extern struct {
    double brjj[5];
    double brmm[5];
    double brll[5];
    double brss[5];
    double brcc[5];
    double brbb[5];
    double brtt[5];
    double brww[3];
    double brzz[3];
    double brgg[5];
    double brzg[5];
    double brhhh[4];
    double brhaa[3][3];
    double brhchc[3];
    double brhaz[2][3];
    double braha[3];
    double brahz[3][2];
    double brhcw[5];
    double brhiggs[5];
    double brneu[5][5][5];
    double brcha[3][5];
    double brhsq[10][3];
    double brhsl[7][3];
    double brasq[6][2];
    double brasl[3][2];
    double brsusy[5];
    double width[5];
  } brn_;

}

#endif
