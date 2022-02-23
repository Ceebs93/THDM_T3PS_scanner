#if !defined(SM_H)
#define SM_H

#include <gsl/gsl_matrix.h>

/**
* @brief Class describing the Standard Model
*
* This class holds all the parameters specifying the Standard Model, and
* methods to access them in convenient ways.
*
*/
class SM {

 public:

  /// Default value for EM coupling \f$ \alpha_{\rm{EM}}(M_Z) \f$
  constexpr static double alpha    = 1./127.934;

  /// Default value for EM coupling at low energy
  constexpr static double alpha0   = 1./137.0359997;

  /// Default value for the Fermi constant \f$ G_F \f$ (\f$ \rm{GeV}^{-2} \f$)
  constexpr static double GF      = 1.16637E-5;

  /// Default value for \f$ M_Z \f$ (GeV)
  constexpr static double MZ      = 91.15349;

  /// Default value for \f$ M_W \f$ (GeV)
  constexpr static double MW      = 80.36951;

  /// Default value for the strong coupling \f$ \alpha_s(M_Z) \f$
  constexpr static double alpha_s = 0.119;

  /// Default value for total width \f$ \Gamma_Z \f$ (GeV)
  constexpr static double GammaZ  = 2.49581;

  /// Default value for total width \f$ \Gamma_W \f$ (GeV)
  constexpr static double GammaW  = 2.08856;

  /// Default value for quark pole mass \f$ m_d \f$ (GeV)
  constexpr static double md_p    = 0.0;

  /// Default value for quark pole mass \f$ m_u \f$ (GeV)
  constexpr static double mu_p    = 0.0;

  /// Default value for quark pole mass \f$ m_s \f$ (GeV)
  constexpr static double ms_p    = 0.1;

  // Scale for strange quark mass (GeV)
  constexpr static double Q_ms    = 2.0;

  /// Default value for quark pole mass \f$ m_c \f$ (GeV)
  constexpr static double mc_p    = 1.42;

  /// Default value for quark pole mass \f$ m_b \f$ (GeV)
  constexpr static double mb_p    = 4.75;

  /// Default value for quark pole mass \f$ m_t \f$ (GeV)
  constexpr static double mt_p    = 172.5;

  /// Default value for electron pole mass \f$ m_e \f$ (GeV)
  constexpr static double me_p    = 0.510998918E-3;

  /// Default value for muon pole mass \f$ m_\mu \f$ (GeV)
  constexpr static double mmu_p   = 0.105658367;

  /// Default value for tau pole mass \f$ m_\tau \f$ (GeV)
  constexpr static double mtau_p  = 1.77684;

  // PDG booklet 2012
  /// Default value for CKM matrix element \f$ V_{ud} \f$
  constexpr static double Vud     = 0.97427;
  /// Default value for CKM matrix element \f$ V_{us} \f$
  constexpr static double Vus     = 0.22534;
  /// Default value for CKM matrix element \f$ V_{ub} \f$
  constexpr static double Vub     = 0.00351;
  /// Default value for CKM matrix element \f$ V_{cd} \f$
  constexpr static double Vcd     = 0.22520;
  /// Default value for CKM matrix element \f$ V_{cs} \f$
  constexpr static double Vcs     = 0.97334;
  /// Default value for CKM matrix element \f$ V_{cb} \f$
  constexpr static double Vcb     = 0.0412;
  /// Default value for CKM matrix element \f$ V_{td} \f$
  constexpr static double Vtd     = 0.00867;
  /// Default value for CKM matrix element \f$ V_{ts} \f$
  constexpr static double Vts     = 0.0404;
  /// Default value for CKM matrix element \f$ V_{tb} \f$
  constexpr static double Vtb     = 0.999146;


  /**
  * @brief Default constructor initializing the Standard Model with default parameter
  *        values
  *
  * This constructor initializes a new SM object with default parameters. The
  * values used as defaults are specified by public constants in this class.
  */
  SM();


  /**
  * @brief Sets the EM coupling \f$ \alpha \f$
  *
  * @param alpha_in New value for \f$ \alpha(M_Z) \f$
  */
  void set_alpha(double alpha_in);

  /**
  * @brief Sets the EM coupling \f$ \alpha_0 \f at low energy$
  *
  * @param alpha_in New value for \f$ \alpha(Q^2=0) \f$
  */
  void set_alpha0(double alpha_in);


  /**
  * @brief Sets the Fermi constant \f$ G_F \f$
  *
  * @param GF_in New value for \f$ G_F \f$ (\f$ GeV^{-2} \f$)
  */
  void set_GF(double GF_in);

  /**
  * @brief Sets \f$ M_Z \f$
  *
  * @param MZ_in New value for \f$ M_Z \f$
  */
  void set_MZ(double MZ_in);

  /**
  * @brief Sets \f$ M_W \f$
  *
  * @param MW_in New value for \f$ M_W \f$
  */
  void set_MW(double MW_in);

  /**
  * @brief Sets \f$ \gamma_Z \f$
  *
  * @param GammaZ_in New value for \f$ \Gamma_Z \f$
  */
  void set_gamma_Z(double GammaZ_in);

  /**
  * @brief Sets \f$ \gamma_W \f$
  *
  * @param GammaW_in New value for \f$ \Gamma_W \f$
  */
  void set_gamma_W(double GammaW_in);

  /**
  * @brief Sets the strong coupling \f$ \alpha_s \f$
  *
  * @param alpha_s_in New value for \f$ \alpha_s(M_Z) \f$
  */
  void set_alpha_s(double alpha_s_in);

  /**
  * @brief Current value of \f$ \alpha \f$
  *
  * @returns Value of \f$ \alpha \f$
  */
  double get_alpha();
  double get_alpha0();

  /**
  * @brief Current value of \f$ G_F \f$
  *
  * @returns Value of \f$ G_F \f$ (\f$ \rm{GeV}^{-2} \f$)
  */
  double get_GF();

  /**
  * @brief Current value of \f$ \sin\theta_W \f$
  *
  * @returns Value of \f$ \sin\theta_W \f$
  */
  double get_sintw();

  /**
  * @brief Current value of \f$ \cos\theta_W \f$
  *
  * @returns Value of \f$ \cos\theta_W \f$
  */
  double get_costw();

  /**
  * @brief Current value of weak SU(2) coupling \f$ g \f$
  *
  * @returns Value of \f$ g \f$
  */
  double get_g();

  /**
  * @brief Current value of hypercharge U(1) coupling \f$ g' \f$
  *
  * @returns Value of \f$ g' \f$
  */
  double get_gprime();

  /**
  * @brief Current value of EM coupling \f$ e \f$
  *
  * @returns Value of \f$ g' \f$
  */
  double get_e();

  /**
  * @brief Current value of \f$ \alpha_s \f$
  *
  * @returns Value of \f$ \alpha_s \f$
  */
  double get_alpha_s();

  /**
  * @brief Current value of \f$ M_W \f$
  *
  * @returns Value of \f$ M_W \f$ (GeV)
  */
  double get_MW();

  /**
  * @brief Current value of \f$ \Gamma_W \f$
  *
  * @returns Value of \f$ \Gamma_W \f$ (GeV)
  */
  double get_gamma_W();

  /**
  * @brief Current value of \f$ M_Z \f$
  *
  * @returns Value of \f$ M_Z \f$ (GeV)
  */
  double get_MZ();

  /**
  * @brief Current value of \f$ \Gamma_Z \f$
  *
  * @returns Value of \f$ \Gamma_Z \f$ (GeV)
  */
  double get_gamma_Z();

  /**
  * @brief Current value of the vev \f$ v=1/\sqrt{\sqrt{2}G_F} \f$
  *
  * @returns Value of \f$ v \f$ (GeV)
  */
  double get_v();

  /**
  * @brief Current value of the vev squared \f$ v^2=1/(\sqrt{2}G_F) \f$
  *
  * @returns Value of \f$ v^2 \f$ (\f$ \rm{GeV}^2 \f$)
  */
  double get_v2();

  /**
  * @brief Vector boson mass
  *
  * This method returns the mass of a given vector boson. The numbering
  * convention used throughout 2HDMC is such that (1 = photon, 2 = Z0, 3 = W+).
  *
  * @param v Index (1-3) for the vector boson
  *
  * @returns Value of \f$ m_V \f$ (GeV)
  */
  double get_vmass(int v);

  /**
  * @brief Vector boson width
  *
  * This method returns the decay width of a given vector boson. The numbering
  * convention is such that (1 = photon, 2 = Z0, 3 = W+).
  *
  * @param v Index (1-3) for the vector boson
  *
  * @returns Value of \f$ \Gamma_V \f$ (GeV)
  */
  double get_gamma_V(int v);

 /**
  * @brief Sets the CKM matrix to a 3x3 unit matrix
  *
  * This method sets the CKM matrix to a 3x3 unit matrix, i.e. eliminating
  * all charged current interactions between different generations.
  */
  void set_diagonal_CKM();

  /**
  * @brief Full CKM matrix
  *
  * This method can be used to retrieve the full CKM matrix as a gsl_matrix.
  *
  * @returns CKM matrix
  *
  * @see get_MD, get_MU, get_ML, get_CKM_element
  */
  gsl_matrix* get_CKM_matrix();

  /**
  * @brief Element of the CKM matrix
  *
  * This method can be used to retrieve a certain element of the CKM matrix.
  *
  * @param i Row index (1,2,3 = u,c,t) of CKM matrix element
  * @param j Column index (1,2,3 = d,s,b) of CKM matrix element
  * @returns The requested CKM element
  *
  * @see get_CKM
  */
  double get_CKM_element(int i, int j);

  /**
  * @brief Sets the pole mass of a lepton
  *
  * This method sets the pole mass of lepton \a l to the supplied value.
  *
  * @param l Index of lepton (1,2,3 = \f$ e,\mu,\tau \f$)
  * @param lmass_in New mass of lepton
  */
  void set_lmass_pole(int l, double lmass_in);

  /**
  * @brief Sets the pole mass of a quark
  *
  * This method sets the pole mass of quark \a q to the supplied value.
  *
  * @param q Index of quark (1,2,3,4,5,6 = \f$ d,u,s,c,b,t \f$)
  * @param qmass_in New mass of quark
  */
  void set_qmass_pole(int q, double qmass_in);

  void set_qmass_msbar(int flav, double qmass_in);

  /**
  * @brief Sets the pole mass of a down-type quark
  *
  * This method sets the pole mass of quark \a d to the supplied value.
  *
  * @param d Index of quark (1,2,3 = \f$ d,s,b \f$)
  * @param dmass_in New mass of quark
  */
  void set_dmass_pole(int d, double dmass_in);

  /**
  * @brief Sets the pole mass of an up-type quark
  *
  * This method sets the pole mass of quark \a u to the supplied value.
  *
  * @param u Index of quark (1,2,3 = \f$ u,c,t \f$)
  * @param umass_in New mass of quark
  */
  void set_umass_pole(int u, double umass_in);

  /**
  * @brief Quark pole mass
  *
  * @param q Index (1,2,3,4,5,6 = \f$ d,u,s,c,b,t \f$) of quark
  *
  * @returns The pole mass \f$ m_q \f$
  */
  double get_qmass_pole(int q);

  /**
  * @brief Down-type quark pole mass
  *
  * @param d Index (1,2,3 = \f$ d,s,b \f$) of down-type quark
  *
  * @returns Quark pole mass \f$ m_d \f$
  */
  double get_dmass_pole(int d);

  /**
  * @brief Up-type quark pole mass
  *
  * @param u Index (1,2,3 = \f$ u,c,t \f$) of up-type quark
  *
  * @returns Quark pole mass \f$ m_u \f$
  */
  double get_umass_pole(int u);

  /**
  * @brief Lepton pole mass
  *
  * @param l Index (1,2,3 = \f$ e,\mu,\tau \f$) of lepton
  *
  * @returns The pole mass \f$ m_l \f$
  */
  double get_lmass_pole(int l);


  /**
  * @brief Mass matrix for down-type quarks
  *
  * This method can be used to retrieve, in matrix form, the current pole masses
  * used for the down-type quarks. The masses are returned as the diagonal
  * elements of a 3x3 gsl_matrix that must be allocated by the user.
  *
  * @returns Diagonal mass matrix for the down-type quarks
  *
  * @see get_MU, get_ML, get_CKM
  */
  gsl_matrix* get_MD();

  /**
  * @brief Mass matrix for up-type quarks
  *
  * This method can be used to retrieve, in matrix form, the current pole masses
  * used for the up-type quarks. The masses are returned as the diagonal
  * elements of a 3x3 gsl_matrix that must be allocated by the user.
  *
  * @returns Diagonal mass matrix for the up-type quarks
  *
  * @see get_MD, get_ML, get_CKM
  */
  gsl_matrix* get_MU();

  /**
  * @brief Mass matrix for leptons
  *
  * This method can be used to retrieve, in matrix form, the current pole masses
  * used for the leptons. The masses are returned as the diagonal
  * elements of a 3x3 gsl_matrix that must be allocated by the user.
  *
  * @returns Diagonal mass matrix for the leptons
  *
  * @see get_MD, get_MU, get_CKM
  */
  gsl_matrix* get_ML();


  /**
  * @brief Quark \f$ \overline{\rm{MS}}\f$ mass
  *
  * This method determines the \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_q(\overline{m}_q) \f$
  *
  * @param q Index (1,2,3,4,5,6 = \f$ d,u,s,c,b,t \f$) of quark
  *
  * @returns Quark \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_q(\overline{m}_q) \f$
  *
  * @see run_qmass_MSbar
  */
  double get_qmass_MSbar(int q);

  /**
  * @brief Down-type quark \f$ \overline{\rm{MS}}\f$ mass
  *
  * This method determines the \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_d(\overline{m}_d) \f$
  *
  * @param d Index (1,2,3 = \f$ d,s,b \f$) of down-type quark
  *
  * @returns Quark \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_d(\overline{m}_d) \f$
  *
  * @see get_umass_MSbar, run_qmass_MSbar
  */
  double get_dmass_MSbar(int d);

  /**
  * @brief Up-type quark \f$ \overline{\rm{MS}}\f$ mass
  *
  * This method determines the \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_u(\overline{m}_u) \f$
  *
  * @param u Index (1,2,3 = \f$ u,c,t \f$) of up-type quark
  *
  * @returns Quark \f$ \overline{\rm{MS}}\f$ mass \f$ \overline{m}_u(\overline{m}_u) \f$
  *
  * @see get_dmass_MSbar, run_qmass_MSbar
  */
  double get_umass_MSbar(int u);

  /**
  * @brief Evaluates running quark masses
  *
  * This method evaluates the running quark_mass \a quark_mass, initially specified
  * at the input scale \a Qinit, at the new scale \a Qfin. The calculation is performed
  * in the \f$ \overline{\rm{MS}} \f$ scheme with variable number of active
  * flavours and threshold matching.
  *
  * @param quark_mass The input value for the running quark mass
  * @param Qinit      Starting scale
  * @param Qfin       Final scale at which the mass is to be evaluated
  * @param mtop       \f$ \overline{\rm{MS}} \f$ top mass (for thresholds)
  * @param mbot       \f$ \overline{\rm{MS}} \f$ bottom mass (for thresholds)
  *
  * @returns The running mass at the scale \a Qfin
  */
  double run_qmass_MSbar(double quark_mass, double Qinit, double Qfin, double mtop, double mbot);


  /**
  * @brief Evaluates the running strong coupling
  *
  * This method evaluates the running strong coupling \f$ \alpha_s \f$ in the
  * \f$ \overline{\rm{MS}} \f$ scheme at the scale specified by \a Q.
  *
  * @param Q          Scale at which to evaluate the running coupling
  * @param mtop       \f$ \overline{\rm{MS}} \f$ top mass (for thresholds)
  * @param mbot       \f$ \overline{\rm{MS}} \f$ bottom mass (for thresholds)
  *
  * @returns The strong coupling \f$ \alpha_s \f$ at the scale \a Q
  */
  double run_alphas_MSbar(double Q, double mtop, double mbot);

  /**
  * @brief Number of active flavours
  *
  * This method returns the number of active quark flavours to be used at
  * a certain mass scale when decoupling the heavier quarks. The thresholds
  * are determined from the quark \f$ \overline{\rm{MS}} \f$ masses.
  *
  * @param M The scale for which the number of flavours is desired
  *
  * @returns The number of quark flavours active at this scale
  */
  int get_Nactivef(double M);


  /**
  * @brief %SM width of the top quark
  *
  * Calculates the width of the top quark in the %SM, i.e. no contributions from
  * \f$ t\to H^+ b \f$ decays etc. are included.
  *
  * @returns The top width \f$ \Gamma_t \f$
  */
  double get_gamma_top();
  double get_gamma_tWd(int d);




  // --- These masses are for compatibility test with HDECAY ---
  // Scale for HD quark masses (HD)
  constexpr static double Q_HD    = 5.0;

  // switch for HD quark masses (HD)
  constexpr static bool b_HD    = false;

  /// Default value for running quark mass at 5 GeV \f$ m_s(5) \f$ (GeV)
  constexpr static double ms_5    = 0.08160951656;

  /// Default value for running quark mass at 5 GeV \f$ m_c(5) \f$ (GeV)
  constexpr static double mc_5    = 0.8960809463;

  /// Default value for running quark mass at 5 GeV \f$ m_b(5) \f$ (GeV)
  constexpr static double mb_5    = 4.310679848;

  /// Default value for running quark mass at 5 GeV \f$ m_t(5) \f$ (GeV)
  constexpr static double mt_5    = 246.5552058;
  // -------------------------------------------------------------


 private:

  // Internal variables which store the current values of the SM parameters
  double m_alpha;
  double m_alpha0;
  double m_GF;
  double m_s2tW;
  double m_alpha_s;
  double m_md_p[3];
  double m_mu_p[3];
  double m_ml_p[3];
  double m_MW;
  double m_MZ;
  double m_GammaW;
  double m_GammaZ;
  double m_CKM[3][3];

  double m_qmass_ms[7];

	// Support functions for the running mass calculation
  double m_threshold(double as);
  double RQ(double as, int nf, int loops);
  double QCD_beta(int c, int nf);

  void clear_lookup();


};

#endif
