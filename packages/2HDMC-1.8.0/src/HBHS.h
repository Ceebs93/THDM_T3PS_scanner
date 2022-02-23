#ifndef HBHS_H
#define HBHS_H

#include "DecayTable.h"
#include "THDM.h"
#include <array>

struct HBHSResult;

//! Class wrapping HiggsBounds and HiggsSignals.
class HBHS {
public:
  static constexpr size_t nHzero = 3; //! number of neutral Higgs bosons
  static constexpr size_t nHplus = 1; //! number of charged Higgs bosons

  /**
   * @brief Construct a new HBHS object.
   *
   * Initializes HiggsBounds and HiggsSignals and optionally sets the
   * HiggsSignals PDF.
   *
   * @param higgssignals_pdf the kind of pdf to use (see set_pdf() in the
   * HiggsSignals documentation)
   */
  HBHS(int higgssignals_pdf = 2);

  //! Finishes HiggsBounds and HiggsSignals
  ~HBHS();

  /**
   * @brief set the mass uncertainties for use in HiggsBounds and
   * HiggsSignals.
   *
   * @param dMh neutral Higgs mass uncertainties
   * @param dMHp charged Higgs mass uncertainties
   */
  void set_mass_uncertainties(const std::array<double, nHzero> &dMh,
                              const std::array<double, nHplus> &dMHp);

  /**
   * @brief Checks the specified model with HiggsBounds and HiggsSignals
   *
   * @param model the model data
   * @return HBHSResult contains all HiggsBounds and HiggsSignals results
   */
  HBHSResult check(THDM &model) const;

private:
  using ArrayHzero = std::array<double, nHzero>;
  using MatrixHzero = std::array<std::array<double, nHzero>, nHzero>;

  struct EffC {
    ArrayHzero ghjss_s;
    ArrayHzero ghjss_p;
    ArrayHzero ghjcc_s;
    ArrayHzero ghjcc_p;
    ArrayHzero ghjbb_s;
    ArrayHzero ghjbb_p;
    ArrayHzero ghjtt_s;
    ArrayHzero ghjtt_p;
    ArrayHzero ghjmumu_s;
    ArrayHzero ghjmumu_p;
    ArrayHzero ghjtautau_s;
    ArrayHzero ghjtautau_p;
    ArrayHzero ghjWW;
    ArrayHzero ghjZZ;
    ArrayHzero ghjZga;
    ArrayHzero ghjgaga;
    ArrayHzero ghjgg;
    MatrixHzero ghjhiZ;
  };

  /**
   * @brief Calculates the effective couplings.
   *
   * Estimates the loop induced effective couplings through the SM-normalized
   * partial widths.
   *
   * @param model the model data
   * @return HB_effC the obtained effective couplings
   */
  EffC effective_couplings(THDM &model) const;

  struct NonSMBR {
    ArrayHzero BR_hjinvisible;
    std::array<std::array<std::array<double, nHzero>, nHzero>, nHzero>
        BR_hkhjhi;
    MatrixHzero BR_hjhiZ;
    std::array<std::array<double, nHplus>, nHzero> BR_hjHpW;
    ArrayHzero BR_hjemu;
    ArrayHzero BR_hjetau;
    ArrayHzero BR_hjmutau;
  };

  /**
   * @brief Obtains the branching ratios into new physics for HiggsBounds.
   *
   * @param model the model data
   * @return NonSMBR obtained branching ratios
   */
  static NonSMBR nonSM_branching_ratios(THDM &model);

  /**
   * @brief Handles the neutral Higgs input for HiggsBounds.
   *
   * Uses the effective coupling input scheme.
   *
   * @param model the model data
   */
  void neutral_input(THDM &model) const;

  /**
   * @brief Handles the charged Higgs input for HiggsBounds
   *
   * Uses the tbH+ cross section tabulated in HiggsBounds.
   *
   * @param model the model data
   */
  void charged_input(THDM &model) const;
};

/**
 * @brief HiggsBounds result.
 *
 * Contains all data returned by run_HiggsBounds_full(). See the HiggsBounds
 * documentation.
 */
struct HBResult {
  std::array<int, HBHS::nHzero + HBHS::nHplus + 1> result;
  std::array<double, HBHS::nHzero + HBHS::nHplus + 1> obsratio;
  std::array<int, HBHS::nHzero + HBHS::nHplus + 1> channel;
  std::array<int, HBHS::nHzero + HBHS::nHplus + 1> ncombined;

  void print() const;
};

/**
 * @brief HiggsSignals result.
 *
 * Contains all data returned by run_HiggsSignals_full() as well as the 13TeV
 * LHC signal rates for the most important channels as returned by
 * get_Rvalues(). See the HiggsSignals documentation.
 */
struct HSResult {
  double chisq_mu;
  double chisq_mh;
  double chisq;
  int nobs;
  double pvalue;

  std::array<double, HBHS::nHzero> r_H_WW;
  std::array<double, HBHS::nHzero> r_H_ZZ;
  std::array<double, HBHS::nHzero> r_H_gaga;
  std::array<double, HBHS::nHzero> r_H_tautau;
  std::array<double, HBHS::nHzero> r_H_bb;
  std::array<double, HBHS::nHzero> r_VH_bb;

  void print() const;
};

//! Combined HiggsBounds and HiggsSignals results.
struct HBHSResult {
  HBResult hb;
  HSResult hs;
};
#endif
