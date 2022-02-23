#if !defined(THDM_H)
#define THDM_H

#include "SM.h"
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix.h>

using namespace std;

struct HBHSResult;

/**
* \mainpage 2HDMC Class documentation
* 2HDMC is a general-purpose calculator for the two-Higgs doublet model.
* It allows parametrization of the Higgs potential in many different ways,
* convenient specification of generic Yukawa sectors, the evaluation of decay
* widths, theoretical constraints and much more.
*/

/**
* @brief Specifies a general two-Higgs doublet model
*
* This class implements a description of a general two-Higgs doublet model
* in terms of the potential parameters and Yukawa couplings necessary for
* a full specification of the model. Several methods are available to set
* the model parameters to various special cases.
*
* From a THDM object, all the couplings between physical states can be
* accessed, including Higgs-fermion, Higgs-Vector, Higgs-Higgs, 3H and 4H
* couplings.
*/
//static int thdmc_set_param(int key, int npara, double *para, double *res);

int thdmc_set_param(int key, double smpara[], double para[], double res[], int slha);

class THDM {

 public:



  /**
  * @brief Default constructor.
  *
  * Empty default constructor which does nothing except initialization of this
  * THDM object. The model is not specified in any way by using this constructor.
  */
  THDM();


  /**
  * @brief Sets the underlying SM
  *
  * This method can be used to specify the SM parameters which this 2HDM
  * is using.
  *
  * @param sm_in SM object specifying the parameters to use in this 2HDM
  */
  void set_SM(SM sm_in);

  /**
  * @brief Returns the underlying SM
  *
  * @returns The SM object specifying the parameters used in this 2HDM
  */
  SM get_SM();

  /**
  * @brief Specifies 2HDM using generic potential
  *
  * This method lets the user specify the 2HDM using a generic basis.
  * For details on the basis choices available for the Higgs potential, we
  * refer to the complete <a href="http://arxiv.org/abs/0902.0851">manual</a>.
  *
  * @param lambda1 Value of \f$ \lambda_1 \f$
  * @param lambda2 Value of \f$ \lambda_2 \f$
  * @param lambda3 Value of \f$ \lambda_3 \f$
  * @param lambda4 Value of \f$ \lambda_4 \f$
  * @param lambda5 Value of \f$ \lambda_5 \f$
  * @param lambda6 Value of \f$ \lambda_6 \f$
  * @param lambda7 Value of \f$ \lambda_7 \f$
  * @param m12_2 Soft \f$ Z_2 \f$-breaking parameter
  * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_gen(double lambda1, double lambda2, double lambda3,
                     double lambda4, double lambda5, double lambda6,
                     double lambda7, double m12_2, double tan_beta);

  /**
  * @brief Specifies 2HDM in the Higgs basis
  *
  * This method lets the user specify the 2HDM in the Higgs basis \f$ (v_2\equiv 0) \f$.
  * For details on the basis choices available for the Higgs potential, we
  * refer to the complete <a href="http://arxiv.org/abs/0902.0851">manual</a>.
  *
  * @param Lambda1 Value of \f$ \Lambda_1 \f$
  * @param Lambda2 Value of \f$ \Lambda_2 \f$
  * @param Lambda3 Value of \f$ \Lambda_3 \f$
  * @param Lambda4 Value of \f$ \Lambda_4 \f$
  * @param Lambda5 Value of \f$ \Lambda_5 \f$
  * @param Lambda6 Value of \f$ \Lambda_6 \f$
  * @param Lambda7 Value of \f$ \Lambda_7 \f$
  * @param m_Hp    Mass of the charged Higgs
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_higgs(double Lambda1, double Lambda2, double Lambda3,
                       double Lambda4, double Lambda5, double Lambda6,
                       double Lambda7, double m_Hp);

  /**
  * @brief Specifies 2HDM in the hybrid basis of 1507.04281
  *
  * @param mh Value of light CP-even Higgs mass \f$ m_h \f$
  * @param mH Value of heavy CP-even Higgs mass \f$ m_H \f$
  * @param cba Mixing parameter \f$ \cos(\beta-\alpha) \f$. Valid range: -1 < cba <= 1.
  * @param Z_4 Value of \f$ Z_4 \f$
  * @param Z_5 Value of \f$ Z_5 \f$
  * @param Z_7 Value of \f$ Z_7 \f$
  * @param tanb Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_hybrid(double mh, double mH, double cba, double Z4,double Z5, double Z7, double tanb);

  /**
  * @brief Specifies 2HDM in the hybrid basis of 1507.04281, but with sba instead of cba
  *
  * @param mh Value of light CP-even Higgs mass \f$ m_h \f$
  * @param mH Value of heavy CP-even Higgs mass \f$ m_H \f$
  * @param sba Mixing parameter \f$ \sin(\beta-\alpha) \f$. Valid range: -1 < cba <= 1.
  * @param Z_4 Value of \f$ Z_4 \f$
  * @param Z_5 Value of \f$ Z_5 \f$
  * @param Z_7 Value of \f$ Z_7 \f$
  * @param tanb Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_hybrid_sba(double mh, double mH, double sba, double Z4,double Z5, double Z7, double tanb);

  /**
  * @brief Specifies 2HDM using potential from the Higgs Hunter's Guide
  *
  * This method lets the user specify the 2HDM using the CP-conserving form
  * of the potential given by Eq. (4.8) in "The Higgs Hunter's Guide".
  * For details on the basis choices available for the Higgs potential, we
  * refer to the complete <a href="http://arxiv.org/abs/0902.0851">manual</a>.
  *
  * @param lambda1 Value of \f$ \lambda_1 \f$
  * @param lambda2 Value of \f$ \lambda_2 \f$
  * @param lambda3 Value of \f$ \lambda_3 \f$
  * @param lambda4 Value of \f$ \lambda_4 \f$
  * @param lambda5 Value of \f$ \lambda_5 \f$
  * @param lambda6 Value of \f$ \lambda_6 \f$
  * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_HHG(double lambda1, double lambda2, double lambda3,
                     double lambda4, double lambda5, double lambda6,
                     double tan_beta);

  /**
  * @brief Specifies 2HDM in the physical basis
  *
  * This method lets the user specify the 2HDM using the basis of physical
  * Higgs masses. For details on the basis choices available for the Higgs
  * potential, we refer to the complete <a href="http://arxiv.org/abs/0902.0851">manual</a>.
  *
  * @param m_h  Mass of lightest CP-even Higgs \f$ h \f$
  * @param m_H  Mass of heavier CP-even Higgs \f$ H \f$
  * @param m_A  Mass of CP-odd Higgs \f$ A \f$
  * @param m_Hp Mass of charged Higgs
  * @param sba  Mixing parameter \f$ \sin(\beta-\alpha) \f$. NB: Correct sign on \f$ \sin(\beta-\alpha) \f$ must be determined from the condition \f$ \cos(\beta-\alpha) \geq 0 \f$
  * @param lambda6  Value of \f$ \lambda_6 \f$ in generic potential
  * @param lambda7  Value of \f$ \lambda_7 \f$ in generic potential
  * @param m12_2    Soft \f$ Z_2 \f$-breaking parameter
  * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_param_phys(double m_h,double m_H, double m_A, double m_Hp,
                      double sba, double lambda6, double lambda7,
                      double m12_2, double tan_beta);

  bool set_param_sm(double mh);


  /**
  * @brief Specifies the 2HDM of the tree-level MSSM
  *
  * This method lets the user specify a 2HDM with the properties of the tree-level
  * MSSM in terms of masses and coupling relations. The Yukawa sector is also
  * automatically selected to be of type II.
  *
  * @param m_A      Mass of CP-odd Higgs \f$ A \f$
  * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_MSSM(double m_A, double tan_beta);

  /**
  * @brief Specifies the 2HDM of the tree-level MSSM + mass corrections to h ("hMSSM")
  *
  * This method lets the user specify a 2HDM with the properties of the tree-level
  * MSSM, with additional mass corrections from the lambda_2 contribution to the (2,2)
  * element of the mass matrix. These can be identified as the leading (mt^4) MSSM
  * corrections. The Yukawa sector is automatically selected to be of type II.
  *
  * @param m_h      Mass of light CP-even Higgs \f$ h \f$
  * @param m_A      Mass of CP-odd Higgs \f$ A \f$
  * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_hMSSM(double mh, double mA, double tanb);


  /**
  * @brief Specifies the 2HDM with one "inert Higgs" doublet
  *
  * This method lets the user specify a 2HDM with special properties as follows:
  * only one of the doublets acquires a vev, there is an exact \f$ Z_2 \f$-symmetry,
  * preventing interdoublet mixing, and the doublet without a vev has no Yukawa
  * couplings. This results in one doublet being "inert" (or dark). The lightest
  * Higgs originating from this doublet is then stable, thus a dark matter candidate.
  * The conventions used for specifying the parameters for this model are the same as
  * in <a href="http://arxiv.org/abs/0810.3924">arXiv:0810.3924</a>.
  *
  * @param m_h  Mass of "SM-like" Higgs - NB: Changed meaning compared to "physical" basis
  * @param m_H  Mass of CP-even "inert" Higgs - NB: Changed meaning compared to "physical" basis
  * @param m_A  Mass of CP-odd "inert" Higgs
  * @param m_Hp Mass of charged "inert" Higgs
  * @param lambda2 Value of \f$ \lambda_2 \f$ parameter
  * @param lambda3 Value of \f$ \lambda_2 \f$ parameter
  *
  * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
  */
  bool set_inert(double m_h,double m_H, double m_A, double m_Hp, double lambda2, double lambda3);


  /**
  * @brief Returns parameter set in the generic basis
  *
  * This method returns a consistent set of parameter values describing the
  * current model in the generic basis.
  *
  * @param lambda1 Returned value of \f$ \lambda_1 \f$
  * @param lambda2 Returned value of \f$ \lambda_2 \f$
  * @param lambda3 Returned value of \f$ \lambda_3 \f$
  * @param lambda4 Returned value of \f$ \lambda_4 \f$
  * @param lambda5 Returned value of \f$ \lambda_5 \f$
  * @param lambda6 Returned value of \f$ \lambda_6 \f$
  * @param lambda7 Returned value of \f$ \lambda_7 \f$
  * @param m12_2   Returned value of \f$ m_{12}^2 \f$
  * @param tan_beta Returned value of \f$ \tan\beta \f$
  */
  void get_param_gen(double &lambda1, double &lambda2, double &lambda3,
                     double &lambda4, double &lambda5, double &lambda6,
                     double &lambda7, double &m12_2,   double &tan_beta);


  /**
  * @brief Returns parameter set in the Higgs basis
  *
  * This method returns a consistent set of parameter values describing the
  * current model in the Higgs basis.
  *
  * @param Lambda1 Returned value of \f$ \Lambda_1 \f$
  * @param Lambda2 Returned value of \f$ \Lambda_2 \f$
  * @param Lambda3 Returned value of \f$ \Lambda_3 \f$
  * @param Lambda4 Returned value of \f$ \Lambda_4 \f$
  * @param Lambda5 Returned value of \f$ \Lambda_5 \f$
  * @param Lambda6 Returned value of \f$ \Lambda_6 \f$
  * @param Lambda7 Returned value of \f$ \Lambda_7 \f$
  * @param m_Hp    Returned value of the charged Higgs mass
  */
  void get_param_higgs(double &Lambda1, double &Lambda2, double &Lambda3,
                       double &Lambda4, double &Lambda5, double &Lambda6,
                       double &Lambda7, double &m_Hp);

  void get_param_hybrid(double &m_h, double &m_H, double &sba,
                        double &Z4, double &Z5, double &Z7, double &tan_beta);


  /**
  * @brief Returns parameter set in Higgs Hunter's Guide basis
  *
  * This method returns a consistent set of parameter values describing the
  * current model in the basis used in the Higgs Hunter's Guide. NB: The use
  * of this method assumes \f$ \lambda_6=\lambda_7=0 \f$ in the generic basis.
  *
  * @param lambda1 Returned value of \f$ \lambda_1 \f$
  * @param lambda2 Returned value of \f$ \lambda_2 \f$
  * @param lambda3 Returned value of \f$ \lambda_3 \f$
  * @param lambda4 Returned value of \f$ \lambda_4 \f$
  * @param lambda5 Returned value of \f$ \lambda_5 \f$
  * @param lambda6 Returned value of \f$ \lambda_6 \f$
  * @param tan_beta Returned value of \f$ \tan\beta \f$
  */
  void get_param_HHG(double &lambda1, double &lambda2, double &lambda3,
                     double &lambda4, double &lambda5, double &lambda6,
                     double &tan_beta);


  /**
  * @brief Returns parameter set in physical basis
  *
  * This method returns a consistent set of parameter values describing the
  * current model in the "physical" basis.
  *
  * @param m_h      Returned value of lightest CP-even Higgs mass
  * @param m_H      Returned value of heaviest CP-even Higgs mass
  * @param m_A      Returned value of CP-odd Higgs mass
  * @param m_Hp     Returned value of charged Higgs mass
  * @param sba      Returned value of \f$ \sin(\beta-\alpha) \f$
  * @param lambda6  Returned value of \f$ \lambda_6 \f$
  * @param lambda7  Returned value of \f$ \lambda_7 \f$
  * @param m12_2    Returned value of \f$ m_{12}^2 \f$
  * @param tan_beta Returned value of \f$ \tan\beta \f$
  */
  void get_param_phys(double &m_h,double &m_H, double &m_A, double &m_Hp,
                      double &sba, double &lambda6, double &lambda7,
                      double &m12_2, double &tan_beta);


  /**
  * @brief Changes basis for the Higgs doublets
  *
  * This method performs a change of basis, changing the value of
  * \f$ \tan\beta \f$. This results in a recalculation of all the other
  * potential parameters, but the physical couplings remain unchanged.
  * Yukawa couplings are not modified.
  *
  * @param tan_beta Value of \f$ \tan\beta \f$ to use in the new basis
  */
  void recalc_tan_beta(double tan_beta);

  /**
  * @brief Returns the mass of physical Higgs bosons
  *
  * This method returns the mass of a physical Higgs Boson. The numbering
  * convention corresponds to (1,2,3,4 = h,H,A,H+) where always \f$ m_h < m_H \f$.
  *
  * @param h Index of Higgs boson
  *
  * @returns The mass of Higgs boson \a h
  */
  double get_hmass(int h);

  /**
  * @brief Returns invariant \f$ \sin(\beta-\alpha) \f$
  *
  * @returns Value of \f$ \sin(\beta-\alpha) \f$ in model
  */
  double get_sba();

  /**
  * @brief Returns invariant \f$ \cos(\beta-\alpha) \f$
  *
  * @returns Value of \f$ \cos(\beta-\alpha) \f$ in model
  */
  double get_cba();


  /**
  * @brief Higgs coupling \f$ q_{ki} \f$ factors
  *
  * This method returns the invariants \f$ q_{ki} \f$ (<a href="
  * http://arxiv.org/abs/hep-ph/0602242">hep-ph/0602242</a>) which are used
  * for the triple and quartic Higgs couplings.
  *
  * @param k Higgs index (1--4)
  * @param i Coupling index (1--2)
  *
  * @returns The value of the coefficient \f$ q_{ki} \f$
  */
  complex <double> get_qki(int k, int i);


  /**
  * @brief Initializes complete Yukawa sector to a specific type
  *
  * This method is used to specify a type of Yukawa sector for the 2HDM.
  * The types (1-4) implemented follow the convention of <a href=""></a>.
  *
  * @param type Type of Yukawa sector (1--4)
  */
  void set_yukawas_type(int type);


  /**
   * @brief Gets the Yukawa type.
   *
   * @return int 1--4 for the standard Yukawa types, 0 for the IDM, -1 for general Yukawas
   */
  int get_yukawas_type();

  /**
  * @brief Initializes and sets diagonal Yukawa couplings for down-type quarks
  *
  * This method initializes and sets the Yukawa couplings for the down-type quarks.
  * The diagonal matrix \f$ \kappa^D \f$ is automatically specified to comply with
  * the quark pole masses, whereas the diagonal elements of \f$ \rho^D \f$ can be
  * specified as input.
  *
  * @param rhod Diagonal element of \f$ \rho^D \f$ for the \f$ d \f$ quark
  * @param rhos Diagonal element of \f$ \rho^D \f$ for the \f$ s \f$ quark
  * @param rhob Diagonal element of \f$ \rho^D \f$ for the \f$ b \f$ quark
  *
  * @see set_yukawas_up, set_yukawas_lepton
  */
  void set_yukawas_down(double rhod, double rhos, double rhob);


  /**
  * @brief Initializes and sets diagonal Yukawa couplings for up-type quarks
  *
  * This method initializes and sets the Yukawa couplings for the up-type quarks.
  * The diagonal matrix \f$ \kappa^U \f$ is automatically specified to comply with
  * the quark pole masses, whereas the diagonal elements of \f$ \rho^U \f$ can be
  * specified as input.
  *
  * @param rhou Diagonal element of \f$ \rho^U \f$ for the \f$ u \f$ quark
  * @param rhoc Diagonal element of \f$ \rho^U \f$ for the \f$ c \f$ quark
  * @param rhot Diagonal element of \f$ \rho^U \f$ for the \f$ t \f$ quark
  *
  * @see set_yukawas_down, set_yukawas_lepton
  */
  void set_yukawas_up(double rhou, double rhoc, double rhot);

  /**
  * @brief Initializes and sets diagonal Yukawa couplings for the charged leptons
  *
  * This method initializes and sets the Yukawa couplings for the leptons.
  * The diagonal matrix \f$ \kappa^L \f$ is automatically specified to comply with
  * the lepton pole masses, whereas the diagonal element of \f$ \rho^L \f$ can be
  * specified as input.
  *
  * @param rhoe   Diagonal element of \f$ \rho^L \f$ for \f$ e \f$
  * @param rhomu  Diagonal element of \f$ \rho^L \f$ for \f$ \mu \f$
  * @param rhotau Diagonal element of \f$ \rho^L \f$ for \f$ \tau \f$
  *
  * @see set_yukawas_down, set_yukawas_up
  */
  void set_yukawas_lepton(double rhoe, double rhomu, double rhotau);


  /**
  * @brief Initializes and sets Yukawa couplings for down-type quarks
  *
  * This method initializes and sets the Yukawa couplings for the down-type quarks.
  * The diagonal matrix \f$ \kappa^D \f$ is automatically specified to comply with
  * the quark pole masses, whereas the symmetric matrix \f$ \rho^D \f$ is given by
  * the six elements specified as input.
  *
  * @param rho11 Yukawa coupling \f$ \rho^D_{dd} \f$
  * @param rho22 Yukawa coupling \f$ \rho^D_{ss} \f$
  * @param rho33 Yukawa coupling \f$ \rho^D_{bb} \f$
  * @param rho12 Yukawa coupling \f$ \rho^D_{ds} \f$
  * @param rho13 Yukawa coupling \f$ \rho^D_{db} \f$
  * @param rho23 Yukawa coupling \f$ \rho^D_{sb} \f$
  */
  void set_yukawas_down(double rho11,double rho22,double rho33,double rho12,double rho13,double rho23);


  /**
  * @brief Initializes and sets Yukawa couplings for up-type quarks
  *
  * This method initializes and sets the Yukawa couplings for the up-type quarks.
  * The diagonal matrix \f$ \kappa^U \f$ is automatically specified to comply with
  * the quark pole masses, whereas the symmetric matrix \f$ \rho^U \f$ is given by
  * the six elements specified as input.
  *
  * @param rho11 Yukawa coupling \f$ \rho^U_{uu} \f$
  * @param rho22 Yukawa coupling \f$ \rho^U_{cc} \f$
  * @param rho33 Yukawa coupling \f$ \rho^U_{tt} \f$
  * @param rho12 Yukawa coupling \f$ \rho^U_{uc} \f$
  * @param rho13 Yukawa coupling \f$ \rho^U_{ut} \f$
  * @param rho23 Yukawa coupling \f$ \rho^U_{ct} \f$
  */
  void set_yukawas_up(double rho11,double rho22,double rho33,double rho12,double rho13,double rho23);


  /**
  * @brief Initializes and sets Yukawa couplings for leptons
  *
  * This method initializes and sets the Yukawa couplings for the leptons.
  * The diagonal matrix \f$ \kappa^L \f$ is automatically specified to comply with
  * the lepton masses, whereas the symmetric matrix \f$ \rho^L \f$ is given by
  * the six elements specified as input.
  *
  * @param rho11 Yukawa coupling \f$ \rho^L_{ee} \f$
  * @param rho22 Yukawa coupling \f$ \rho^L_{\mu\mu} \f$
  * @param rho33 Yukawa coupling \f$ \rho^L_{\tau\tau} \f$
  * @param rho12 Yukawa coupling \f$ \rho^L_{e\mu} \f$
  * @param rho13 Yukawa coupling \f$ \rho^L_{e\tau} \f$
  * @param rho23 Yukawa coupling \f$ \rho^L_{\mu\tau} \f$
  */
  void set_yukawas_lepton(double rho11,double rho22,double rho33,double rho12,double rho13,double rho23);

  /**
  * @brief Initializes and sets Yukawa couplings for the inert 2HDM
  *
  * This method initializes and sets the Yukawa couplings for the "inert"
  * model. That means there are no Yukawa couplings of any type, except what
  * is necessary for the mass terms, i.e. all \f$ \rho = 0 \f$.
  */
  void set_yukawas_inert();


  /**
  * @brief Returns diagonal elements of Yukawa matrix for down-type quarks
  *
  * @param rhod Returned value of \f$ \rho^D_{dd} \f$
  * @param rhos Returned value of \f$ \rho^D_{ss} \f$
  * @param rhob Returned value of \f$ \rho^D_{bb} \f$
  */
  void get_yukawas_down(double &rhod, double &rhos, double &rhob);


  /**
  * @brief Returns diagonal elements of Yukawa matrix for up-type quarks
  *
  * @param rhou Returned value of \f$ \rho^U_{uu} \f$
  * @param rhoc Returned value of \f$ \rho^U_{cc} \f$
  * @param rhot Returned value of \f$ \rho^U_{tt} \f$
  */
  void get_yukawas_up(double &rhou, double &rhoc, double &rhot);


  /**
  * @brief Returns diagonal elements of Yukawa matrix for leptons
  *
  * @param rhoe Returned value of \f$ \rho^L_{ee} \f$
  * @param rhomu Returned value of \f$ \rho^L_{\mu\mu} \f$
  * @param rhotau Returned value of \f$ \rho^L_{\tau\tau} \f$
  */
  void get_yukawas_lepton(double &rhoe, double &rhomu, double &rhotau);


  void get_kappa_down(double &kd, double &ks, double &kb);
  void get_kappa_up(double &ku, double &kc, double &kt);
  void get_kappa_lepton(double &ke, double &kmu, double &ktau);

  void get_kappa_down(double mu, double &kd, double &ks, double &kb);
  void get_kappa_up(double mu, double &ku, double &kc, double &kt);
  void get_kappa_lepton(double mu, double &ke, double &kmu, double &ktau);
  void get_rho_down(double mu, double &rd, double &rs, double &rb);
  void get_rho_up(double mu, double &ru, double &rc, double &rt);
  void get_rho_lepton(double mu, double &re, double &rmu, double &rtau);

  /**
  * @brief Returns Yukawa matrix for down-type quarks
  *
  * @param rho_D_out Returned Yukawa matrix \f$ \rho^D \f$
  */
  void get_yukawas_down(gsl_matrix *rho_D_out);


  /**
  * @brief Returns Yukawa matrix for up-type quarks
  *
  * @param rho_U_out Returned Yukawa matrix \f$ \rho^U \f$
  */
  void get_yukawas_up(gsl_matrix *rho_U_out);


  /**
  * @brief Returns Yukawa matrix for leptons
  *
  * @param rho_L_out Returned Yukawa matrix \f$ \rho^L \f$
  */
  void get_yukawas_lepton(gsl_matrix *rho_L_out);


  /**
  * @brief Couplings of Higgses to down-type fermions
  *
  * Calculates the coupling \f$ hdd \f$ between one physical Higgs state, specified
  * by \a h, and two down-type quarks \a f1 and \a f2.
  *
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param f1 First fermion (1,2,3 = d,s,b)
  * @param f2 Second fermion (1,2,3 = d,s,b)
  * @param cs Returned (complex) value for scalar coupling
  * @param cp Returned (complex) value for pseudoscalar coupling
  */
  void get_coupling_hdd(int h,int f1,int f2,complex <double> &cs, complex <double> &cp);


  /**
  * @brief Couplings of Higgses to up-type fermions
  *
  * Calculates the coupling \f$ huu \f$ between one physical Higgs state, specified
  * by \a h, and two up-type quarks \a f1 and \a f2.
  *
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param f1 First fermion (1,2,3 = u,c,t)
  * @param f2 Second fermion (1,2,3 = u,c,t)
  * @param cs Returned (complex) value for scalar coupling
  * @param cp Returned (complex) value for pseudoscalar coupling
  */
  void get_coupling_huu(int h,int f1,int f2,complex <double> &cs, complex <double> &cp);


  /**
  * @brief Couplings of Higgses to mixed type fermions
  *
  * Calculates the coupling \f$ hdu \f$ between one physical Higgs state, specified
  * by \a h (only relevant one is charged Higgs), and two quarks \a f1 and \a f2.
  *
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param d Down-type fermion (1,2,3 = d,s,b)
  * @param u Up-type fermion (1,2,3 = u,c,t)
  * @param cs Returned (complex) value for scalar coupling
  * @param cp Returned (complex) value for pseudoscalar coupling
  */
  void get_coupling_hdu(int h,int d,int u,complex <double> &cs, complex <double> &cp);

  /**
  * @brief Couplings of Higgses to charged leptons
  *
  * Calculates the coupling \f$ hll \f$ between one physical Higgs state, specified
  * by \a h, and two leptons \a f1 and \a f2.
  *
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param f1 First fermion (1,2,3 = \f$ e,\mu,\tau \f$)
  * @param f2 Second fermion (1,2,3 = \f$ e,\mu,\tau \f$)
  * @param cs Returned (complex) value for scalar coupling
  * @param cp Returned (complex) value for pseudoscalar coupling
  */
  void get_coupling_hll(int h,int f1,int f2,complex <double> &cs, complex <double> &cp);

  /**
  * @brief Couplings of Higgses to mixed leptons
  *
  * Calculates the coupling \f$ hl\nu_l \f$ between one physical Higgs state, specified
  * by \a h (only relevant one is charged Higgs), and two leptons \a f1 and \a f2.
  *
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param l  Charged lepton (1,2,3 = \f$ e,\mu,\tau \f$)
  * @param n  Neutrino (1,2,3 = \f$ \nu_e,\nu_\mu,\nu_\tau \f$)
  * @param cs Returned (complex) value for scalar coupling
  * @param cp Returned (complex) value for pseudoscalar coupling
  */
  void get_coupling_hln(int h,int l,int n,complex <double> &cs, complex <double> &cp);


  /**
  * @brief Couplings of Higgses to pairs of vector bosons
  *
  * Calculates the coupling \f$ hV_1 V_2 \f$ between one physical Higgs state, specified
  * by \a h , and two vector bosons \a v1 and \a v2. NB: Neutral Higgses have no coupling
  * to photons implemented, but the loop mediated decays \f$ h\to \gamma \gamma \f$ can
  * nevertheless be calculated using the DecayTable class. Conventions are according to hep-ph/0602242.
  *
  * @param v1 Index of first vector boson (1,2,3 = \f$ \gamma,Z,W^+ \f$)
  * @param v2 Index of second vector boson (1,2,3 = \f$ \gamma,Z,W^+ \f$)
  * @param h  Index of Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param c  Returned (complex) value for coupling
  */
  void get_coupling_vvh(int v1,int v2,int h,complex <double> &c);

  /**
  * @brief Couplings of vector bosons to pairs of Higgses
  *
  * Calculates the coupling \f$ Vh_1 h_2 \f$ between two physical Higgs states, specified
  * by \a h1 and \a h2, and one vector boson \a v. Conventions are according to hep-ph/0602242.
  *
  * @param h1 Index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h2 Index of second Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param v  Index of vector boson (1,2,3 = \f$ \gamma,Z,W^+ \f$)
  * @param c  Returned (complex) value for coupling
  */
  void get_coupling_vhh(int v,int h1,int h2,complex <double> &c);


  /**
  * @brief Triple Higgs couplings
  *
  * Calculates the coupling \f$ h_1 h_2 h_3 \f$ between three physical Higgs states.
  * Conventions are according to hep-ph/0602242
  *
  * @param h1 Index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h2 Index of second Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h3 Index of third Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param c  Returned (complex) value for coupling
  */
  void get_coupling_hhh(int h1,int h2,int h3,complex <double> &c);

  /**
  * @brief Couplings of two vector bosons and two Higgses
  *
  * Calculates the coupling \f$ V_1 V_2 h_1 h_2 \f$ between two vector bosons and
  * two physical Higgs states. Conventions are according to hep-ph/0602242
  *
  * @param v1 Index of first vector boson (1,2,3 = \f$ \gamma,Z,W^+ \f$)
  * @param v2 Index of second vector boson (1,2,3 = \f$ \gamma,Z,W^+ \f$)
  * @param h1 Index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h2 Index of second Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param c  Returned (complex) value for coupling
  */
  void get_coupling_vvhh(int v1,int v2,int h1,int h2,complex <double> &c);

  /**
  * @brief Quartic Higgs couplings
  *
  * Calculates the coupling \f$ h_1 h_2 h_3 h_4 \f$ between four physical Higgs states.
  * Conventions are according to hep-ph/0602242
  *
  * @param h1 Index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h2 Index of second Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h3 Index of third Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param h4 Index of fourth Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param c  Returned (complex) value for coupling
  */
  void get_coupling_hhhh(int h1,int h2,int h3,int h4,complex <double> &c);


  /**
  * @brief Calculates tree-level unitarity constraints
  *
  * This method calculates the eigenvalues of the S-matrix for Higgs-Higgs scattering
  * as defined in <a href="http://xxx.lanl.gov/abs/hep-ph/0508020">hep-ph/0508020</a>.
  *
  * @returns Value of largest eigenvalue
  */
  double calc_unitarity();

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
  * @brief Calculates largest quartic Higgs coupling
  *
  * This method calculates all the quartic Higgs boson couplings and returns the largest one as well as the indices of the four Higgs bosons
  *
  * @param gmax  Returned (complex) value for coupling
  * @param imax  Returned index of first Higgs boson (1,2,3,4 = h,H,A,H+)
  * @param jmax  Returned index of second boson (1,2,3,4 = h,H,A,H+)
  * @param kmax  Returned index of third boson (1,2,3,4 = h,H,A,H+)
  * @param lmax  Returned index of fourth boson (1,2,3,4 = h,H,A,H+)
  */
  void calc_perturbativity(complex <double> &gmax,int &imax,int &jmax,int &kmax,int &lmax);

  /**
  * @brief Checks perturbativity
  *
  * This method checks whether the couplings of the Higgs and Gauge bosons fullfills perturbativity
  *
  * @param perturbativity_limit optional value to be used for perturbativity limit (default = \f$ 4\pi \f$)
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
  * @brief Prints the potential parameters in the generic basis to stdout
  */
  void print_param_gen();

  /**
  * @brief Prints the potential parameters in the Higgs basis to stdout
  */
  void print_param_higgs();


  /**
  * @brief Prints the potential parameters in the hybrid basis to stdout
  */
  void print_param_hybrid();


  /**
  * @brief Prints the potential parameters in the physical basis to stdout
  */
  void print_param_phys();

  /**
  * @brief Prints the potential parameters in the Higgs Hunter's Guide basis to stdout
  */
  void print_param_HHG();

  /**
  * @brief Prints the Yukawa matrices to stdout
  */
  void print_yukawas();

  void print_hdecay();

  /**
  * @brief Reads the 2HDM parameters from a LesHouches compliant file
  *
  * This method reads the parameters necessary to specify the 2HDM (and the %SM)
  * from a file complying with the LesHouches standard. There are several options
  * for how to prepare this file (see <a href="http://arxiv.org/abs/0902.0851">manual</a> for details,
  * or the code distribution for examples).
  *
  * @param file The name of the file to read
  *
  * @returns Boolean \a true if the file could be opened and read correctly AND
  *          a 2HDM was correctly specified in the file, \a false otherwise
  */
  bool read_LesHouches(const char* file);

  /**
  * @brief Writes output in a LesHouches compliant file format
  *
  * This method writes the potential parameters and other selected output to
  * a specified file in LesHouches-compliant file format which can then be
  * further processed by other codes.
  *
  * @param file The name of the file which is to be written
  * @param fulldecay If \a true, decay modes of the Higgs bosons in LesHouches
  *                  format are written to the file
  * @param couplings If \a true, all couplings (> 200 values) for the Higgs bosons
  *                  are written to the file. Should be used when supplying input
  *                  for MadGraph/MadEvent 2HDMC model
  * @param qcd_on    Turns QCD corrections on or off, default is on. Should be
  *                  turned off when supplying input for MadGraph/MadEvent 2HDMC model
  */
  void write_LesHouches(const char* file, bool fulldecay, bool couplings, bool qcd_on=true, const HBHSResult *hbhs=nullptr);

  void write_model(const char* file);

  double get_alpha();

  /**
	* @brief Small value
	*
	*  Minimum value used for widths, branching ratios etc. to determine when
	*  something should be considered zero.
	*/
  constexpr static double EPS = 1E-12;

 private:
  double      lambda[8];
  double      beta;
  double      m22_2;
  double      sinba;
  bool        params_set;
  double      v2;
  gsl_matrix *kappa_D;
  gsl_matrix *kappa_U;
  gsl_matrix *kappa_L;
  gsl_matrix *rho_D;
  gsl_matrix *rho_U;
  gsl_matrix *rho_L;
  gsl_matrix *rho_N;
  int         yukawas_type;

  static bool first_run;
  const static char *version;

  SM sm;

  void init();
  double get_m12_2();
  void set_kappa();
  void set_kappa_D();
  void set_kappa_U();
  void set_kappa_L();

  void print_info();

};

#endif
