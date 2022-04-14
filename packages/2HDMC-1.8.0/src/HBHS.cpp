#include "HBHS.h"
#include "DecayTable.h"
#include "/scratch/cb27g11/THDM_T3PS_scanner/packages/HiggsBounds-5.10.1/include/HiggsBounds.h"
#include "/scratch/cb27g11/THDM_T3PS_scanner/packages/HiggsSignals-2.6.2/include/HiggsSignals.h"
#include "SM.h"
#include "THDM.h"
#include <array>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#if defined HiggsBounds

namespace {
double zeroIfNaN(double x) { return x == x ? x : 0.; }

constexpr bool debug = false;
} // namespace

HBHS::HBHS(int hs_pdf) {
  initialize_HiggsBounds(nHzero, nHplus, 3 /* LandH */);
  initialize_HiggsSignals_latestresults(nHzero, nHplus);
  HiggsSignals_setup_pdf(hs_pdf);
}

HBHS::~HBHS() {
  finish_HiggsBounds();
  finish_HiggsSignals();
}

HBHSResult HBHS::check(THDM &model) const {
  neutral_input(model);
  charged_input(model);

  HBHSResult result{};
  run_HiggsBounds_full(result.hb.result.data(), result.hb.channel.data(),
                       result.hb.obsratio.data(), result.hb.ncombined.data());

  if (debug) {
    for (size_t h = 1; h <= nHzero; ++h) {
      double BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, BR_hjmumu, BR_hjtautau,
          BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg;
      HiggsBounds_get_neutral_BR(h, &BR_hjss, &BR_hjcc, &BR_hjbb, &BR_hjtt,
                                 &BR_hjmumu, &BR_hjtautau, &BR_hjWW, &BR_hjZZ,
                                 &BR_hjZga, &BR_hjgaga, &BR_hjgg);
      printf(
          "HiggsBounds BRs:\n %2ld %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E "
          "%16.8E %16.8E %16.8E %16.8E %16.8E\n",
          h, BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, BR_hjmumu, BR_hjtautau,
          BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg);
      double singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH,
          ggZH;
      HiggsBounds_get_neutral_hadr_CS(h, LHC13, &singleH, &ggH, &bbH, &VBF, &WH,
                                      &ZH, &ttH, &tH_tchan, &tH_schan, &qqZH,
                                      &ggZH);
      printf(
          "HiggsBounds LHC13 cxns:\n %2ld %16.8E %16.8E %16.8E %16.8E %16.8E "
          "%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
          h, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH,
          ggZH);
    }
  }

  run_HiggsSignals_full(&result.hs.chisq_mu, &result.hs.chisq_mh,
                        &result.hs.chisq, &result.hs.nobs, &result.hs.pvalue);

  for (int i = 0; i != nHzero; ++i)
    get_HiggsSignals_Rvalues(i + 1, 4 /* LHC13 */, &result.hs.r_H_WW[i],
                             &result.hs.r_H_ZZ[i], &result.hs.r_H_gaga[i],
                             &result.hs.r_H_tautau[i], &result.hs.r_H_bb[i],
                             &result.hs.r_VH_bb[i]);
  return result;
}

HBHS::EffC HBHS::effective_couplings(THDM &model) const {
  EffC effC{};
  std::array<double, 3> Mh = {model.get_hmass(1), model.get_hmass(2),
                              model.get_hmass(3)};

 // DecayTable table(model);
  THDM sm_like;
  DecayTable table(model);

  for (size_t h = 1; h <= nHzero; ++h) {
    sm_like.set_param_sm(Mh[h - 1]);
    sm_like.set_yukawas_type(1);

    std::complex<double> cs, cp, cs_sm, cp_sm;
    // h -> ss
    model.get_coupling_hdd(h, 2, 2, cs, cp);
    sm_like.get_coupling_hdd(1, 2, 2, cs_sm, cp_sm);
    effC.ghjss_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjss_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "ss", effC.ghjss_s[h - 1],
             effC.ghjss_p[h - 1]);

    model.get_coupling_hdd(h, 3, 3, cs, cp);
    sm_like.get_coupling_hdd(1, 3, 3, cs_sm, cp_sm);
    effC.ghjbb_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjbb_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "bb", effC.ghjbb_s[h - 1],
             effC.ghjbb_p[h - 1]);

    model.get_coupling_huu(h, 2, 2, cs, cp);
    sm_like.get_coupling_huu(1, 2, 2, cs_sm, cp_sm);
    effC.ghjcc_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjcc_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "cc", effC.ghjcc_s[h - 1],
             effC.ghjcc_p[h - 1]);

    model.get_coupling_huu(h, 3, 3, cs, cp);
    sm_like.get_coupling_huu(1, 3, 3, cs_sm, cp_sm);
    effC.ghjtt_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjtt_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "tt", effC.ghjtt_s[h - 1],
             effC.ghjtt_p[h - 1]);

    model.get_coupling_hll(h, 2, 2, cs, cp);
    sm_like.get_coupling_hll(1, 2, 2, cs_sm, cp_sm);
    effC.ghjmumu_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjmumu_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "mumu", effC.ghjmumu_s[h - 1],
             effC.ghjmumu_p[h - 1]);

    model.get_coupling_hll(h, 3, 3, cs, cp);
    sm_like.get_coupling_hll(1, 3, 3, cs_sm, cp_sm);
    effC.ghjtautau_s[h - 1] = cs.imag() / cs_sm.imag();
    effC.ghjtautau_p[h - 1] = -cp.real() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E %16.8E\n", h, "tata", effC.ghjtautau_s[h - 1],
             effC.ghjtautau_p[h - 1]);

    model.get_coupling_vvh(2, 2, h, cs);
    sm_like.get_coupling_vvh(2, 2, 1, cs_sm);
    effC.ghjZZ[h - 1] = cs.imag() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "ZZ", effC.ghjZZ[h - 1]);

    model.get_coupling_vvh(3, 3, h, cs);
    sm_like.get_coupling_vvh(3, 3, 1, cs_sm);
    effC.ghjWW[h - 1] = cs.imag() / cs_sm.imag();
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "WW", effC.ghjWW[h - 1]);

    DecayTable sm_table(sm_like);
    double hgaga = table.get_gamma_hgaga(h);
    double hgaga_sm = sm_table.get_gamma_hgaga(1);
    effC.ghjgaga[h - 1] = sqrt(hgaga / hgaga_sm);
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "gaga", effC.ghjgaga[h - 1]);

    double hZga = table.get_gamma_hZga(h);
    double hZga_sm = sm_table.get_gamma_hZga(1);
    effC.ghjZga[h - 1] = zeroIfNaN(sqrt(hZga / hZga_sm));
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "Zga", effC.ghjZga[h - 1]);

    double hgg = table.get_gamma_hgg(h);
    double hgg_sm = sm_table.get_gamma_hgg(1);
    effC.ghjgg[h - 1] = sqrt(hgg / hgg_sm);
    if (debug)
      printf("%2ld %5s %16.8E\n", h, "gg", effC.ghjgg[h - 1]);

    for (size_t j = 1; j <= nHzero; ++j) {
      std::complex<double> c;
      model.get_coupling_vhh(2, h, j, c);
      effC.ghjhiZ[h - 1][j - 1] =
          c.real() / (model.get_SM().get_g() / 2. / model.get_SM().get_costw());
    }
    if (debug)
      printf("%2ld %5s %16.8E %16.8E %16.8E\n", h, "Zhj", effC.ghjhiZ[h - 1][0],
             effC.ghjhiZ[h - 1][1], effC.ghjhiZ[h - 1][2]);
  }
  return effC;
}

HBHS::NonSMBR HBHS::nonSM_branching_ratios(THDM &model) {
  DecayTable table(model);
  NonSMBR nonSMBR{};
  for (size_t h = 1; h <= nHzero; ++h) {
    nonSMBR.BR_hjinvisible[h - 1] = 0.0;
    nonSMBR.BR_hjHpW[h - 1][0] =
        zeroIfNaN(table.get_gamma_hvh(h, 3, 4) / table.get_gammatot_h(h));
    if (debug)
      printf("BR h%2ld -> H+- W-+ = %16.8E\n", h, nonSMBR.BR_hjHpW[h - 1][0]);
    nonSMBR.BR_hjemu[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 1, 2) / table.get_gammatot_h(h));
    nonSMBR.BR_hjetau[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 1, 3) / table.get_gammatot_h(h));
    nonSMBR.BR_hjmutau[h - 1] =
        zeroIfNaN(table.get_gamma_hll(h, 2, 3) / table.get_gammatot_h(h));
  }
  for (size_t j = 1; j <= nHzero; j++) {
    for (size_t i = 1; i <= nHzero; i++) {
      for (size_t k = 1; k <= nHzero; k++) {
        nonSMBR.BR_hkhjhi[k - 1][j - 1][i - 1] =
            zeroIfNaN(table.get_gamma_hhh(k, j, i) / table.get_gammatot_h(k));
        if (debug)
          printf("BR h%2ld -> h%2ld h%2ld = %16.8E\n", k, j, i,
                 nonSMBR.BR_hkhjhi[k - 1][j - 1][i - 1]);
      }
      nonSMBR.BR_hjhiZ[j - 1][i - 1] =
          zeroIfNaN(table.get_gamma_hvh(j, 2, i) / table.get_gammatot_h(j));
      if (debug)
        printf("BR h%2ld -> h%2ld Z = %16.8E\n", j, i,
               nonSMBR.BR_hjhiZ[j - 1][i - 1]);
    }
  }
  if (model.get_yukawas_type() == 0) {
    if (model.get_hmass(1) < model.get_hmass(3)) {
      nonSMBR.BR_hjinvisible[1] = nonSMBR.BR_hkhjhi[1][0][0];
      nonSMBR.BR_hkhjhi[1][0][0] = 0;
    } else if (model.get_hmass(1) > model.get_hmass(3)) {
      nonSMBR.BR_hjinvisible[1] = nonSMBR.BR_hkhjhi[1][2][2];
      nonSMBR.BR_hkhjhi[1][2][2] = 0;
    } else {
      nonSMBR.BR_hjinvisible[1] =
          nonSMBR.BR_hkhjhi[1][0][0] + nonSMBR.BR_hkhjhi[1][2][2];
      nonSMBR.BR_hkhjhi[1][0][0] = 0;
      nonSMBR.BR_hkhjhi[1][2][2] = 0;
    }
    if (debug)
      printf("BR_hSM->invisible = %16.8E\n", nonSMBR.BR_hjinvisible[1]);
  }
  return nonSMBR;
}

void HBHS::neutral_input(THDM &model) const {
  // properties
  const std::array<double, nHzero> Mh{model.get_hmass(1), model.get_hmass(2),
                                      model.get_hmass(3)};

  // let HiggsBounds calculate the total widths internally
  std::array<double, nHzero> GammaTotal{-1., -1., -1.};
  // unless we are in the intert model
  if (model.get_yukawas_type() == 0) {
    DecayTable table{model};
    GammaTotal[0] = table.get_gammatot_h(1);
    GammaTotal[2] = table.get_gammatot_h(3);
  }
  const std::array<int, nHzero> CP_value{+1, +1, -1};
  HiggsBounds_neutral_input_properties(Mh.data(), GammaTotal.data(),
                                       CP_value.data());
  // effective couplings
  const EffC effC = effective_couplings(model);
  HiggsBounds_neutral_input_effC(
      effC.ghjss_s.data(), effC.ghjss_p.data(), effC.ghjcc_s.data(),
      effC.ghjcc_p.data(), effC.ghjbb_s.data(), effC.ghjbb_p.data(),
      effC.ghjtt_s.data(), effC.ghjtt_p.data(), effC.ghjmumu_s.data(),
      effC.ghjmumu_p.data(), effC.ghjtautau_s.data(), effC.ghjtautau_p.data(),
      effC.ghjWW.data(), effC.ghjZZ.data(), effC.ghjZga.data(),
      effC.ghjgaga.data(), effC.ghjgg.data(), effC.ghjhiZ[0].data());

  // non-SM decay branching ratios
  const auto nonSMBR = nonSM_branching_ratios(model);
  HiggsBounds_neutral_input_nonSMBR(
      nonSMBR.BR_hjinvisible.data(), nonSMBR.BR_hkhjhi[0][0].data(),
      nonSMBR.BR_hjhiZ[0].data(), nonSMBR.BR_hjemu.data(),
      nonSMBR.BR_hjetau.data(), nonSMBR.BR_hjmutau.data(),
      nonSMBR.BR_hjHpW[0].data());
}

void HBHS::charged_input(THDM &model) const {
  DecayTable table{model};
  SM sm = model.get_SM();

  const std::array<double, nHplus> MHp{model.get_hmass(4)};
  const std::array<double, nHplus> HPlusGammaTot{table.get_gammatot_h(4)};
  const std::array<double, nHplus> CS_lep_HpjHmi_ratio{1.};
  const std::array<double, nHplus> BR_tHpjb{table.get_gamma_uhd(3, 4, 3) /
                                            table.get_gammatot_top()};
  const std::array<double, nHplus> BR_Hpjcs{table.get_gamma_hdu(4, 2, 2) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hpjcb{table.get_gamma_hdu(4, 3, 2) /
                                            table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hptaunu{table.get_gamma_hln(4, 3, 3) /
                                              table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_Hptb{table.get_gamma_hdu(4, 3, 3) /
                                           table.get_gammatot_h(4)};
  const std::array<double, nHplus> BR_HpWZ{0.};
  std::array<std::array<double, nHzero>, nHplus> BR_HphiW{};
  for (size_t i = 1; i <= nHzero; i++) {
    BR_HphiW[0][i - 1] = table.get_gamma_hvh(4, 3, i) / table.get_gammatot_h(4);
  }
  HiggsBounds_charged_input(
      MHp.data(), HPlusGammaTot.data(), CS_lep_HpjHmi_ratio.data(),
      sm.get_gamma_top() / table.get_gammatot_top(), BR_tHpjb.data(),
      BR_Hpjcs.data(), BR_Hpjcb.data(), BR_Hptaunu.data(), BR_Hptb.data(),
      BR_HpWZ.data(), BR_HphiW[0].data());

  // get charged Higgs couplings
  std::complex<double> cs, cp;
  model.get_coupling_hdu(4, 3, 3, cs, cp);
  const double mb = sm.get_dmass_MSbar(3);
  double mbrun = 0.;
  if (mb > 0) {
    mbrun = sm.run_qmass_MSbar(mb, mb, model.get_hmass(4), sm.get_qmass_pole(6),
                               sm.get_qmass_pole(5));
  }
  const double kappa_b =
      mb > 0 ? (cs + cp).imag() * sm.get_v() / (sqrt(2) * mbrun) : 0;

  const double mt = sm.get_umass_MSbar(3);
  const double mtrun = sm.run_qmass_MSbar(
      mt, mt, model.get_hmass(4), sm.get_qmass_pole(6), sm.get_qmass_pole(5));
  const double kappa_t = -(-cs + cp).imag() * sm.get_v() / (sqrt(2) * mtrun);
  if (debug)
    printf("Hccoups %16.8E %16.8E\n", kappa_t, kappa_b);
  // pass cross sections (mostly zero)
  std::array<double, nHplus> CS_Hpmjtb{
      HCCS_tHc(model.get_hmass(4), kappa_t, kappa_b,
               table.get_gamma_uhd(3, 4, 3) / table.get_gammatot_top())};
  const std::array<double, nHzero> zero{};
  HiggsBounds_charged_input_hadr(LHC13, CS_Hpmjtb.data(), zero.data(),
                                 zero.data(), zero.data(), zero.data(),
                                 zero.data(), zero.data(), zero.data(),
                                 zero.data(), zero.data());
}

void HBHS::set_mass_uncertainties(
    const std::array<double, HBHS::nHzero> &dMh,
    const std::array<double, HBHS::nHplus> &dMHp) {
  HiggsBounds_set_mass_uncertainties(dMh.data(), dMHp.data());
  HiggsSignals_neutral_input_MassUncertainty(dMh.data());
}

void HBResult::print() const {
  printf("\nHiggsBounds results:\n");
  printf("  Higgs  res  chan       ratio        ncomb\n");
  for (int i = 1; i <= 4; i++) {
    printf("%5d %5d %6d %16.8E %5d   %s\n", i, result[i], channel[i],
           obsratio[i], ncombined[i], result[i] == 1 ? "Allowed" : "Excluded");
  }
  printf("------------------------------------------------------------\n");
  printf("  TOT %5d %6d %16.8E %5d   %s\n", result[0], channel[0], obsratio[0],
         ncombined[0], result[0] == 1 ? "ALLOWED" : "EXCLUDED");
}

void HSResult::print() const {
  printf("\nHiggsSignals results:\n");
  printf(" Chi^2 from rates: %16.8E\n", chisq_mu);
  printf("  Chi^2 from mass: %16.8E\n", chisq_mh);
  printf("      Total chi^2: %16.8E\n", chisq);
  printf("    # observables: %16d\n\n", nobs);
}

#endif
