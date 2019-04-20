#if !defined(DECAY_H)
#define DECAY_H

#include "SM.h"
#include "THDM.h"
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>

struct BR
{
  double brdd[4][4];
  double bruu[4][4];
	double brdu[4][4];
	double brll[4][4];
	double brln[4][4];
	double brvv[4];
	double brvh[4][5];
	double brhh[5];
	double brhZga;
	double brhgg;
};


using namespace std;

/**
* @brief Calculates the decay modes of 2HDM Higgs bosons
* 
* Given a THDM object, a DecayTable can be generated. From this table, the
* Higgs boson decay widths and branching ratios are obtained. For the complete
* list of available decay modes, we refer to the complete documentation, or
* the list of member methods below.
*/
class DecayTable {

 public: 
    
  /**
  * @brief Default constructor
  * 
  * This default constructor takes a THDM object as argument for which 
  * the decays are to be calculated. %SM properties are taken from the SM
  * object in the THDM.
  * 
  * @param mod Two-Higgs doublet model for which to calculate decay modes
  */
  DecayTable(THDM mod);


  /**
  * @brief Sets underlying 2HDM
  * 
  * This method sets the THDM underlying the DecayTable.
  * 
  * @param model Two-Higgs doublet model for which to calculate decay modes
  */
	void set_model(THDM model);


  /**
  * @brief Underlying 2HDM
  * 
  * Use to obtain the underlying THDM object
  * 
  * @returns The THDM object on which this DecayTable operates
  */
	THDM get_model();  


  /**
  * @brief Decay width \f$\Gamma(h\to u_1\overline{u}_2) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to a pair of up-type quarks. QCD corrections have been included.
  * 
  * @param h  Index of Higgs boson (1,2,3 = h,H,A)
  * @param u1 Index of up-type quark (1,2,3 = \f$ u,c,t \f$)
  * @param u2 Index of up-type antiquark (1,2,3 = \f$ \bar{u},\bar{c},\bar{t} \f$)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_huu(int h, int u1, int u2);


  /**
  * @brief Decay width \f$ \Gamma(h\to d_1\overline{d}_2) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to a pair of down-type quarks. QCD corrections have been included.
  * 
  * @param h  Index of Higgs boson (1,2,3 = h,H,A)
  * @param d1 Index of down-type quark (1,2,3 = \f$ d,s,b \f$)
  * @param d2 Index of down-type antiquark (1,2,3 = \f$ \bar{d},\bar{s},\bar{b} \f$)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hdd(int h, int d1, int d2);

  /**
  * @brief Decay width \f$ \Gamma(h\to l_1\overline{l}_2) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to a pair of charged leptons.
  * 
  * @param h  Index of Higgs boson (1,2,3 = h,H,A)
  * @param l1 Index of lepton (1,2,3 = \f$ e,\mu,\tau \f$)
  * @param l2 Index of antilepton (1,2,3 = \f$ \bar{e},\bar{\mu},\bar{\tau} \f$)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hll(int h, int l1, int l2);


  /**
  * @brief Decay width \f$ \Gamma(h\to u\overline{d}) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to a pair of quarks. QCD corrections have been included.
  * 
  * @param h  Index of Higgs boson (4 = H+)
  * @param d  Index of down-type antiquark (1,2,3 = \f$ \bar{d},\bar{s},\bar{b} \f$)
  * @param u  Index of up-type quark (1,2,3 = \f$ u,c,t \f$)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hdu(int h, int d, int u);


  /**
  * @brief Decay width \f$ \Gamma(h\to l\overline{\nu}_{l'}) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to lepton-neutrino.
  * 
  * @param h  Index of Higgs boson (4 = H+)
  * @param l  Index of charged lepton (1,2,3 = \f$ e^+,\mu^+,\tau^+ \f$)
  * @param n  Index of neutrino (1,2,3 = \f$ \nu_e,\nu_\mu,\nu_\tau \f$)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hln(int h, int l, int n);  
  

  /**
  * @brief Decay width \f$ \Gamma(h\to H_1 H_2) \f$
  * 
  * This method calculates the on-shell decay width for the decay of Higgs
  * boson \a h to a pair of Higgs bosons
  * 
  * @param h  Index of decaying Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h1 Index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h2 Index of second Higgs boson (1,2,3,4 = h,H,A,H+)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hhh(int h, int h1, int h2);


  /**
  * @brief Decay width \f$ \Gamma(h\to VV) \f$
  * 
  * This method calculates the decay width for the Higgs boson \a h 
  * to a pair of vector bosons. The decay mode with one vector
  * boson off-shell is included.
  * 
  * @param h  Index of decaying Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param v  Index of vector bosons (1,2,3 = \f$\gamma \f$,Z,W)
  * 
  * @returns The decay width in GeV
  */
  double  get_gamma_hvv(int h, int v);

  /**
  * @brief Decay width \f$ \Gamma(h\to VH) \f$
  * 
  * This method calculates the decay width for the Higgs boson \a h 
  * to one massive vector and one Higgs boson. Decay with the vector
  * boson off-shell is included.
  * 
  * @param h  Index of decaying Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param V  Index of vector boson (2,3 = \f$\gamma \f$,Z,W)
  * @param H  Index of final-state Higgs boson (1,2,3,4 = h,H,A,H+)
	*
  * @returns The decay width in GeV
  */
  double  get_gamma_hvh(int H, int V, int h);

  /**
  * @brief Decay width \f$ \Gamma(h\to gg) \f$
  * 
  * This method calculates the decay width for the Higgs boson \a h 
  * to a pair of gluons. LO QCD corrections are included.
  * 
  * @param h  Index of decaying Higgs boson (1,2,3,4 = h,H,A,H+)
	*
  * @returns The decay width in GeV
  */
  double  get_gamma_hgg(int h);


  /**
  * @brief Decay width \f$ \Gamma(h\to \gamma\gamma) \f$
  * 
  * This method calculates the decay width for the neutral
  * Higgs boson \a h to a pair of photons.
  * 
  * @param h  Index of decaying Higgs boson (1,2,3 = h,H,A)
	*
  * @returns The decay width in GeV
  */
	double  get_gamma_hgaga(int h);

  /**
  * @brief Decay width \f$ \Gamma(h\to Z\gamma) \f$
  * 
  * This method calculates the decay width for the neutral
  * Higgs boson \a h to a Z and a photon.
  * 
  * @param h  Index of decaying Higgs boson (1,2,3 = h,H,A)
	*
  * @returns The decay width in GeV
  */
	double  get_gamma_hZga(int h);

  /**
  * @brief Total width \f$ \Gamma_h \f$
  * 
  * Calculates the total decay width of Higgs boson \a h 
  * 
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
	*
  * @returns Total width in GeV
  */
  double  get_gammatot_h(int h);

  /**
  * @brief Total width \f$ \Gamma_V \f$
  * 
  * Returns the total decay width of vector boson \a v 
  * 
  * @param v  Index of vector boson (2,3 = Z,W)
	*
  * @returns Total width in GeV
  */
  double  get_gammatot_v(int v);

  /**
  * @brief Total width \f$ \Gamma_t \f$
  * 
  * Returns the total decay width of the top quark 
  * 
  * @returns Total width in GeV
  */
  double  get_gammatot_top();


  /**
  * @brief Decay width for \f$ t \to H^+X \f$
  * 
  * Returns the decay width of the top quark in the charged Higgs mode
  * 
  * @returns Decay width in GeV
  * 
  * @param u Index of decaying quark (1,2,3 = \f$ u,c,t \f$)
  * @param h Index of Higgs boson (4 = H+)
	* @param d Index of down-type quark (1,2,3 = \f$ d,s,b \f$)
  */
  double  get_gamma_uhd(int u, int h, int d);


  double  get_gamma_uhu(int u1, int h, int u2);

	void geth_BR(int h, struct BR &br);
  
  /**
  *	@brief Prints the decay modes of a Higgs boson
  *
  * The decay modes of Higgs boson \a h are printed to stdout
  * 
  * @param h Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  */
  void    print_decays(int h);

  /**
  *	@brief Prints the total width of a Higgs boson
  *
  * The total decay width of Higgs boson \a h are printed to stdout
  * 
  * @param h Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  */
  void    print_width(int h);

  /**
  *	@brief Prints the decay modes of the top quark
  *
  * The decay modes of the top quark are printed to stdout
  */
  void    print_top_decays();

  /**
  *	@brief Prints decay information for a Higgs boson in LesHouches format
  *
  * The decay information for Higgs boson \a h are printed to a file in
  * LesHouches format
  *
  * @param output The name of the output file to write
  * @param h 			Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param full 	If \a true, all decay modes are printed, if \a false only the total width
  */
  void    print_decays_LesHouches(FILE* output, int h, bool full);

  /**
  *	@brief Prints decay information for the top quark in LesHouches format
  *
  * The decay information for the top quark are printed to a file in
  * LesHouches format
  *
  * @param output The name of the output file to write
  * @param full 	If \a true, all decay modes are printed, if \a false only the total width
  */
  void    print_top_decays_LesHouches(FILE* output, bool full);

  /**
  *     @brief Turns QCD corrections on or off
  *
  * This method is used to turn QCD corrections on or off. If the output is 
  * meant to be used with the MadGraph/MadEvent 2HDMC model 
  * QCD corrections should be turned off to get a consistent result. 
  *
  * @param set  If \a true QCD corrections are turned on, 
  *             if \a false QCD corrections are turned off
  */
  void    set_qcd(bool set);

  static double  DHp(double ui, double uj, double xi, double xj, double sqL);
  static double  DHm(double ui, double uj, double xi, double xj, double sqL);
  static double  BHp(double ui, double uj, double xi, double xj, double sqL);


 private:
  THDM model;
  SM sm;
  
  complex <double> F_sf(double t);
  complex <double> F_pf(double t);
  complex <double> F_0(double t);
  complex <double> F_1(double t);
  complex <double> ftau(double t);
  complex <double> gtau(double t);
  complex <double> I_1(double tau, double lambda);
  complex <double> I_2(double tau, double lambda);
  complex <double> FF_s(double tau, double lambda);
  complex <double> FF_p(double tau, double lambda);
  complex <double> FW(double tau, double lambda);
  complex <double> FHp(double tau, double lambda);

  double gammatot_h[5];
  double gamma_uhd[5][5][5];
  double gamma_uhu[5][5][5];
  double gamma_hdd[5][5][5];
  double gamma_huu[5][5][5];
  double gamma_hdu[5][5][5];
  double gamma_hll[5][5][5];
  double gamma_hln[5][5][5];
  double gamma_hgg[5];
  double gamma_hgaga[5];
  double gamma_hZga[5];
  double gamma_hvv[5][5];
  double gamma_hvh[5][5][5];
  double gamma_hhh[5][5][5];

  double  br(double dG, double G);

  double  hvv_onshell(int h, int V, double M);
  double  hvv_offshell(int h, int V, double M);
  double  hvv_all(int h, int V, double M);
  double  hff_onshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc, int h, bool tt);
  double  hpff_onshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc, int h);
  double  htt_onshell(double M, double m1, double m2, int Nc, int h);
  double  htt_offshell(double M, double m1, double m2, int Nc, int h);
  double  htb_offshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc);

  double  hvh_onshell(int H, int V, int h, double M);
  double  hvh_offshell(int H, int V, int h, double M);
  double  hdu_offshell(int h, int d, int u, double M);
  double  hgaga(int h);
  double  hZga(int h);
  double  hgg(int h);
  double  PS2(double M, double m1, double m2);
  
  double  interp(double R, double x, double y, double c);

  bool    qcd_on;
  int     err_code;

  void    print_decay_LesHouches(FILE* output, double br, int id1, int id2);
  void    print_decay(const char *h, const char *id1, const char *id2, double g, double br);
  void    print_decays(FILE* output, int h, bool full, bool lh);



};

#endif
