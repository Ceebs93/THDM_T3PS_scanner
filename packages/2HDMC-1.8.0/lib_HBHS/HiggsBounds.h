/**
 * @file
 * @brief C interface to the HiggsBounds subroutines.
 *
 * This API handles conversions in array conventions, string arguments and
 * defines enums for some flags.
 *
 *  - **Argument types:** All `intent(in)` scalar values (i.e. non-arrays) are
 *    passed by value, all `intent(out)` scalar values are passed by pointer
 *    (e.g. `int *HBresult`). For arrays, `intent(in)` arguments are passed as
 *    const arrays (e.g. `const double Mh[]`) and `intent(out)` arguments are
 *    passed as arrays (e.g. `int HBresult[]` in run_HiggsBounds_full()).
 *
 *  - **Array conventions:** Multidimensional arrays are stored in row major
 *    order in C vs column major order in Fortran. All array arguments have to
 *    be contiguous blocks of memory representing an array in row major order.
 *    This is true for native C arrays as well as the C accessors returned by
 *    most C++ array/tensor libraries (e.g. `data()` in Eigen3, Xtensor). All
 *    multidimensional function arguments are documented with detailed
 *    information. A `C` example:
 *
 *        // 3 neutral and 1 charged Higgs bosons, e.g. 2HDM
 *        initialize_HiggsBounds(3, 1, LandH);
 *        ...
 *        double BR_HpjhiW[3][1];
 *        BR_HpjhiW[1][0] = 0.3;
 *        ...
 *        HiggsBounds_neutral_input_nonSMBR( ... , BR_HpjhiW);
 *
 *  - **Any index values passed follow the Fortran convention and start at 1**
 *    (see e.g. HiggsBounds_get_neutral_hadr_CS()).
 *
 *  - **String Arguments:** We avoid passing strings from C to Fortran since it
 *    is very messy. T
 *     - There is no wrapper for initialize_higgsbounds(). Instead we provide a
 *       wrapper around initialize_higgsbounds_int() called
 *       initialize_HiggsBounds().
 *     - We replaced some string arguments by enums (e.g.
 *       #HiggsBounds_likelihood_type).
 *     - We do not provide wrappers for the `_single` input subroutines, use the
 *       corresponding full input subroutines instead. (The exception is
 *       HiggsBounds_neutral_input_hadr_channelrates_single(), which takes no
 *       string arguments. Here only a wrapper around the `_single` verion is
 *       provided as it usually more useful than the full verion.)
 *
 *  - **Enums:** We define enums for some of the integer flags taken as
 *    arguments by the Fortran subroutines. This includes the dataset selection
 *    in initialize_HiggsBounds() (#HiggsBounds_analyses_flag), the collider
 *    selection for the hadronic cross sections (#HiggsBounds_collider_id), and
 *    the likelihood type (#HiggsBounds_likelihood_type).
 *
 *
 * This also includes interfaces to the functions in access_SM.f90 and
 * access_effC.f90. All cross section values are in pb.
 *
 */

#pragma once

#ifdef __cplusplus
#include <cstddef>
#else
#include <stddef.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Select which experimental analyses to use with
 * initialize_HiggsBounds()
 */
enum HiggsBounds_analyses_flag {
    onlyL, //!< only Lep data
    onlyH, //!< only Tevatron and LHC data
    LandH, //!< LEP, Tevatron and LHC data
    onlyP, //!< only published results (i.e. those with an arXiv-id)
    list   //!< only data defined in usefulbits::analysislist
};

/**
 * @brief Select which collider the data refers to.
 *
 * Used with HiggsBounds_neutral_input_hadr() and
 * HiggsBounds_charged_input_hadr(). This value corresponds to the center of
 * mass energy in TeV.
 */
enum HiggsBounds_collider_id {
    Tevatron = 2, //!< 2 TeV \f$p \bar{p}\f$
    LHC7 = 7,     //!< 7 TeV \f$pp\f$
    LHC8 = 8,     //!< 8 TeV \f$pp\f$
    LHC13 = 13    //!< 13 TeV \f$pp\f$
};

/**
 * @brief Select the kind of likelihood to return.
 * Used with the get_likelihood functions (e.g. HiggsBounds_get_likelihood())
 *
 */
enum HiggsBounds_likelihood_type {
    pred = 0, //!< predicted likelihood
    obs = 1   //!< observed likelihood
};

/**
 * @brief Interface to initialize_higgsbounds_int()
 *
 * @param nHiggsneut Number of neutral Higgs bosons in the model
 * @param nHiggsplus Number of charged Higgs bosons in the model
 * @param analyses_flag dataset from #HiggsBounds_analyses_flag
 */
void initialize_HiggsBounds(size_t nHiggsneut, size_t nHiggsplus,
                            int analyses_flag);

/**
 * @brief Interface to higgsbounds_neutral_input_properties()
 * @copydoc higgsbounds_neutral_input_properties()
 */
void HiggsBounds_neutral_input_properties(const double Mh[],
                                          const double GammaTotal_hj[],
                                          const int CP_value[]);

/**
 * @brief Interface to higgsbounds_neutral_input_effc().
 *
 * **The parameter ghjhiZ is a 2-dimensional `nHiggsneut x nHiggsneut` symmetric
 * array.**
 *
 * @copydoc higgsbounds_neutral_input_effc()
 */
void HiggsBounds_neutral_input_effC(
    const double ghjss_s[], const double ghjss_p[], const double ghjcc_s[],
    const double ghjcc_p[], const double ghjbb_s[], const double ghjbb_p[],
    const double ghjtt_s[], const double ghjtt_p[], const double ghjmumu_s[],
    const double ghjmumu_p[], const double ghjtautau_s[],
    const double ghjtautau_p[], const double ghjWW[], const double ghjZZ[],
    const double ghjZga[], const double ghjgaga[], const double ghjgg[],
    const double ghjhiZ[]);

/**
 * @brief Interface to higgsbounds_neutral_input_smbr()
 */
void HiggsBounds_neutral_input_SMBR(
    const double BR_hjss[], const double BR_hjcc[], const double BR_hjbb[],
    const double BR_hjtt[], const double BR_hjmumu[],
    const double BR_hjtautau[], const double BR_hjWW[], const double BR_hjZZ[],
    const double BR_hjZga[], const double BR_hjgaga[], const double BR_hjgg[]);

/**
 * @brief Interface to higgsbounds_neutral_input_nonsmbr()
 *
 *  - **The parameter `BR_hkhjhi` is a 3-dimensional `nHiggsneut x nHiggsneut x
 *    nHiggsneut` array. The element `BR_hkhjhi[k][j][i]` will be accessible in
 *    Fortran as `BR_hkhjhi(k,j,i)`.**
 *  - **The parameter `BR_hjhiZ` is a 2-dimensional `nHiggsneut x nHiggsneut`
 *    array. The element at `BR_hjhiZ[j][i]` will be accessible in Fortran as
 *    `BR_hjhiZ(j,i)`.**
 *  - **The parameter `BR_hjHpiW` is a  2-dimensional `nHiggsneut x nHiggsplus`
 *    array. The element at `BR_hjHpiW[j][i]` will be accessible in Fortran as
 *    `BR_hjHpiW(j,i)`.**
 *
 * @copydoc higgsbounds_neutral_input_nonsmbr()
 *
 */
void HiggsBounds_neutral_input_nonSMBR(
    const double BR_hjinvisible[], const double BR_hkhjhi[],
    const double BR_hjhiZ[], const double BR_hjemu[], const double BR_hjetau[],
    const double BR_hjmutau[], const double BR_hjHpiW[]);

/**
 * @brief Interface to higgsbounds_neutral_input_effc_firstgen()
 * @copydoc higgsbounds_neutral_input_effc_firstgen()
 */
void HiggsBounds_neutral_input_effC_firstgen(
    const double ghjuu_s[], const double ghjuu_p[], const double ghjdd_s[],
    const double ghjdd_p[], const double ghjee_s[], const double ghjee_p[]);

/**
 * @brief Interface to higgsbounds_neutral_input_effc_fv()
 * @copydoc higgsbounds_neutral_input_effc_fv()
 */
void HiggsBounds_neutral_input_effC_FV(
    const double ghjuc_s[], const double ghjuc_p[], const double ghjut_s[],
    const double ghjut_p[], const double ghjct_s[], const double ghjct_p[],
    const double ghjds_s[], const double ghjds_p[], const double ghjdb_s[],
    const double ghjdb_p[], const double ghjsb_s[], const double ghjsb_p[]);

/**
 * @brief Interface to higgsbounds_neutral_input_firstgenbr()
 * @copydoc higgsbounds_neutral_input_firstgenbr()
 */
void HiggsBounds_neutral_input_firstgenBR(const double BR_hjuu[],
                                          const double BR_hjdd[],
                                          const double BR_hjee[]);

/**
 * @brief Interface to higgsbounds_neutral_input_fvbr()
 * @copydoc higgsbounds_neutral_input_fvbr()
 */
void HiggsBounds_neutral_input_FVBR(
    const double BR_hjuc[], const double BR_hjds[], const double BR_hjut[],
    const double BR_hjdb[], const double BR_hjct[], const double BR_hjsb[]);

/**
 * @brief Interface to higgsbounds_neutral_input_lep()
 *
 * **The parameter `XS_ee_hjhi_ratio` is a 2-dimensional `nHiggsneut x
 * nHiggsneut` symmetric array.**
 *
 * @copydoc higgsbounds_neutral_input_lep()
 */
void HiggsBounds_neutral_input_LEP(const double XS_ee_hjZ_ratio[],
                                   const double XS_ee_bbhj_ratio[],
                                   const double XS_ee_tautauhj_ratio[],
                                   const double XS_ee_hjhi_ratio[]);

/**
 * @brief Interface to higgsbounds_neutral_input_hadr()
 *
 *  - **The parameter collider should be specified through the corresponding
 *    #HiggsBounds_collider_id.**
 *  - **The parameter `CS_hjhi` is a 2-dimensional `nHiggsneut x nHiggsneut`
 *    symmetric array.**
 *
 * @copydoc higgsbounds_neutral_input_hadr()
 */
void HiggsBounds_neutral_input_hadr(
    int collider, const double CS_hj_ratio[], const double CS_gg_hj_ratio[],
    const double CS_bb_hj_ratio[], const double CS_hjW_ratio[],
    const double CS_hjZ_ratio[], const double CS_vbf_ratio[],
    const double CS_tthj_ratio[], const double CS_thj_tchan_ratio[],
    const double CS_thj_schan_ratio[], const double CS_qq_hjZ_ratio[],
    const double CS_gg_hjZ_ratio[], const double CS_tWhj_ratio[],
    const double CS_hjhi[]);

/**
 * @brief Interface to higgsbounds_charged_input()
 *
 * **The parameter `BR_HpjhiW` is a 2-dimensional `nHiggsplus x nHiggsneut`
 * array. The element at `BR_HpjhiW[j][i]` will be accessible in Fortran as
 * `BR_HpjhiW(j,i)`.**
 *
 * @copydoc higgsbounds_charged_input()
 */
void HiggsBounds_charged_input(const double Mhplus[],
                               const double GammaTotal_Hpj[],
                               const double CS_ee_HpjHmj_ratio[],
                               double BR_tWpb, const double BR_tHpjb[],
                               const double BR_Hpjcs[], const double BR_Hpjcb[],
                               const double BR_Hpjtaunu[],
                               const double BR_Hpjtb[], const double BR_HpjWZ[],
                               const double BR_HpjhiW[]);

/**
 * @brief Interface to higgsbounds_charged_input_firstgenbr()
 * @copydoc higgsbounds_charged_input_firstgenbr()
 */
void HiggsBounds_charged_input_firstgenBR(const double BR_Hpjud[],
                                          const double BR_Hpjus[],
                                          const double BR_Hpjcd[],
                                          const double BR_Hpjub[],
                                          const double BR_Hpjenu[],
                                          const double BR_Hpjmunu[]);

/**
 * @brief Interface to higgsbounds_charged_input_hadr()
 *
 *  - **The parameter collider should be specified through the corresponding
 *    #HiggsBounds_collider_id.**
 *  - **The parameter `CS_Hpmjhi` is a  2-dimensional `nHiggsplus x nHiggsneut`
 *    array. The element at `CS_Hpmjhi[j][i]` will be accessible in Fortran as
 *    `CS_Hpmjhi(j,i)`.**
 *
 * @copydoc higgsbounds_charged_input_hadr()
 */
void HiggsBounds_charged_input_hadr(
    int collider, const double CS_Hpmjtb[], const double CS_Hpmjcb[],
    const double CS_Hpmjbjet[], const double CS_Hpmjcjet[],
    const double CS_Hpmjjetjet[], const double CS_HpmjW[],
    const double CS_HpmjZ[], const double CS_vbf_Hpmj[],
    const double CS_HpjHmj[], const double CS_Hpmjhi[]);

/**
 * @brief Interface to higgsbounds_charged_input_effc_fermions()
 * @copydoc higgsbounds_charged_input_effc_fermions()
 */
void HiggsBounds_charged_input_effC_fermions(
    const double hcjud_L[], const double hcjud_R[], const double hcjcs_L[],
    const double hcjcs_R[], const double hcjtb_L[], const double hcjtb_R[],
    const double hcjus_L[], const double hcjus_R[], const double hcjub_L[],
    const double hcjub_R[], const double hcjcd_L[], const double hcjcd_R[],
    const double hcjcb_L[], const double hcjcb_R[], const double hcjtd_L[],
    const double hcjtd_R[], const double hcjts_L[], const double hcjts_R[]);

/**
 * @brief Interface to higgsbounds_charged_input_exoticbr()
 * @copydoc higgsbounds_charged_input_exoticbr()
 */
void HiggsBounds_charged_input_exoticBR(const double BR_Hpjud[],
                                        const double BR_Hpjus[],
                                        const double BR_Hpjcd[],
                                        const double BR_Hpjub[],
                                        const double BR_Hpjenu[],
                                        const double BR_Hpjmunu[]);

/**
 * @brief Interface to higgsbounds_get_neutral_hadr_cs()
 *
 *  - **The index `i` of the requested Higgs boson starts at 1.**
 *  - **The parameter collider should be specified through the corresponding
 *    #HiggsBounds_collider_id.**
 *
 * @copydoc higgsbounds_get_neutral_hadr_cs()
 */
void HiggsBounds_get_neutral_hadr_CS(int i, int collider, double *singleH,
                                     double *ggH, double *bbH, double *VBF,
                                     double *WH, double *ZH, double *ttH,
                                     double *tH_tchan, double *tH_schan,
                                     double *qqZH, double *ggZH);

/**
 * @brief Interface to higgsbounds_get_neutral_br()
 *
 * **The index `i` of the requested Higgs boson starts at 1.**
 *
 * @copydoc higgsbounds_get_neutral_br()
 */
void HiggsBounds_get_neutral_BR(int i, double *BR_hjss, double *BR_hjcc,
                                double *BR_hjbb, double *BR_hjtt,
                                double *BR_hjmumu, double *BR_hjtautau,
                                double *BR_hjWW, double *BR_hjZZ,
                                double *BR_hjZga, double *BR_hjgaga,
                                double *BR_hjgg);

/**
 * @brief Interface to higgsbounds_set_mass_uncertainties()
 * @copydoc higgsbounds_set_mass_uncertainties()
 */
void HiggsBounds_set_mass_uncertainties(const double dMhneut[],
                                        const double dMhch[]);

/**
 * @brief Interface to run_higgsbounds()
 * @copydoc run_higgsbounds()
 */
void run_HiggsBounds(int *HBresult, int *chan, double *obsratio,
                     int *ncombined);

/**
 * @brief Interface to run_higgsbounds_single()
 * @copydoc run_higgsbounds_single()
 */
void run_HiggsBounds_single(int h, int *HBresult, int *chan, double *obsratio,
                            int *ncombined);

/**
 * @brief Interface to run_higgsbounds_full()
 * @copydoc run_higgsbounds_full()
 */
void run_HiggsBounds_full(int HBresult[], int chan[], double obsratio[],
                          int ncombined[]);

/**
 * @brief Interface to run_higgsbounds_classic()
 * @copydoc run_higgsbounds_classic()
 */
void run_HiggsBounds_classic(int *HBresult, int *chan, double *obsratio,
                             int *ncombined);

/**
 * @brief Interface to higgsbounds_get_most_sensitive_channels_per_higgs()
 * @copydoc higgsbounds_get_most_sensitive_channels_per_higgs()
 */
void HiggsBounds_get_most_sensitive_channels_per_Higgs(int nH, int pos,
                                                       int *HBresult, int *chan,
                                                       double *obsratio,
                                                       double *predratio,
                                                       int *ncombined);

/**
 * @brief Interface to higgsbounds_get_most_sensitive_channels()
 * @copydoc higgsbounds_get_most_sensitive_channels()
 */
void HiggsBounds_get_most_sensitive_channels(int pos, int *HBresult, int *chan,
                                             double *obsratio,
                                             double *predratio, int *ncombined);

/**
 * @brief Interface to higgsbounds_get_likelihood()
 * - **The parameter obspred should be specified through the
 * #HiggsBounds_likelihood_type.**
 *
 * @copydoc higgsbounds_get_likelihood()
 */
void HiggsBounds_get_likelihood(int analysisID, int *Hindex, int *nc, int *cbin,
                                double *M, double *llh, int obspred);

/**
 * @brief Interface to higgsbounds_get_likelihood_for_higgs()
 * - **The parameter obspred should be specified through the
 * #HiggsBounds_likelihood_type.**
 *
 * @copydoc higgsbounds_get_likelihood_for_higgs()
 */
void HiggsBounds_get_likelihood_for_Higgs(int analysisID, int cbin_in,
                                          int Hindex, int *nc, int *cbin,
                                          double *M, double *llh, int obspred);

/**
 * @brief Interface to higgsbounds_get_likelihood_for_comb()
 * - **The parameter obspred should be specified through the
 * #HiggsBounds_likelihood_type.**
 *
 * @copydoc higgsbounds_get_likelihood_for_comb()
 */
void HiggsBounds_get_likelihood_for_comb(int analysisID, int cbin_in,
                                         int Hindex, int *nc, int *cbin,
                                         double *M, double *llh, int obspred);

#ifdef enableCHISQ
/**
 * @brief Interface to initialize_higgsbounds_chisqtables()
 *
 * @copydoc initialize_higgsbounds_chisqtables()
 */
void initialize_HiggsBounds_chisqtables();

/**
 * @brief Interface to higgsbounds_get_lepchisq()
 *
 * @copydoc higgsbounds_get_lepchisq()
 */
void HiggsBounds_get_LEPChisq(double theory_uncertainty_1s,
                              double *chisq_withouttheory,
                              double *chisq_withtheory, int *channel);

/**
 * @brief Interface to finish_higgsbounds_chisqtables()
 *
 * @copydoc finish_higgsbounds_chisqtables()
 */
void finish_HiggsBounds_chisqtables();
#endif

/**
 * @brief Interface to finish_higgsbounds()
 *
 * @copydoc finish_higgsbounds()
 */
void finish_HiggsBounds();

/**
 * @brief Interface to higgsbounds_neutral_input_hadr_channelrates_single()
 *
 *  - **We do not provide a C interface to
 *    higgsbounds_neutral_input_hadr_channelrates(), use this function
 * instead.**
 *  - **The parameter collider should be specified through the corresponding
 *    #HiggsBounds_collider_id.**
 *
 * @copydoc higgsbounds_neutral_input_hadr_channelrates_single()
 */
void HiggsBounds_neutral_input_hadr_channelrates_single(int collider,
                                                        int nHiggs, int p,
                                                        int d, double val);

/**
 * @brief Interface to higgsbounds_neutral_input_hadr_channelrates_clean()
 *
 * @copydoc higgsbounds_neutral_input_hadr_channelrates_clean()
 */
void HiggsBounds_neutral_input_hadr_channelrates_clean();

/**
 * @brief Interface to smbr_hww()
 * @copydoc smbr_hww()
 */
double SMBR_HWW(double Mh);

/**
 * @brief Interface to smbr_hzz()
 * @copydoc smbr_hzz()
 */
double SMBR_HZZ(double Mh);

/**
 * @brief Interface to smbr_hbb()
 * @copydoc smbr_hbb()
 */
double SMBR_Hbb(double Mh);

/**
 * @brief Interface to smbr_htautau()
 * @copydoc smbr_htautau()
 */
double SMBR_Htautau(double Mh);

/**
 * @brief Interface to smbr_hgamgam()
 * @copydoc smbr_hgamgam()
 */
double SMBR_Hgamgam(double Mh);

/**
 * @brief Interface to smbr_hgg()
 * @copydoc smbr_hgg()
 */
double SMBR_Hgg(double Mh);

/**
 * @brief Interface to smbr_htoptop()
 * @copydoc smbr_htoptop()
 */
double SMBR_Htoptop(double Mh);

/**
 * @brief Interface to smbr_hcc()
 * @copydoc smbr_hcc()
 */
double SMBR_Hcc(double Mh);

/**
 * @brief Interface to smbr_hss()
 * @copydoc smbr_hss()
 */
double SMBR_Hss(double Mh);

/**
 * @brief Interface to smbr_hmumu()
 * @copydoc smbr_hmumu()
 */
double SMBR_Hmumu(double Mh);

/**
 * @brief Interface to smbr_hzgam()
 * @copydoc smbr_hzgam()
 */
double SMBR_HZgam(double Mh);

/**
 * @brief Interface to smgamma_h()
 * @copydoc smgamma_h()
 */
double SMGamma_H(double Mh);

/**
 * @brief Interface to smgamma_twpb()
 * @copydoc smgamma_twpb()
 */
double SMGamma_tWpb(double mt);

/**
 * @brief Interface to smcs_tev_hw()
 * @copydoc smcs_tev_hw()
 */
double SMCS_tev_HW(double Mh);

/**
 * @brief Interface to smcs_tev_hz()
 * @copydoc smcs_tev_hz()
 */
double SMCS_tev_HZ(double Mh);

/**
 * @brief Interface to smcs_tev_gg_h()
 * @copydoc smcs_tev_gg_h()
 */
double SMCS_tev_gg_H(double Mh);

/**
 * @brief Interface to smcs_tev_bb_h()
 * @copydoc smcs_tev_bb_h()
 */
double SMCS_tev_bb_H(double Mh);

/**
 * @brief Interface to smcs_tev_vbf_h()
 * @copydoc smcs_tev_vbf_h()
 */
double SMCS_tev_vbf_H(double Mh);

/**
 * @brief Interface to smcs_tev_bg_hb()
 * @copydoc smcs_tev_bg_hb()
 */
double SMCS_tev_bg_Hb(double Mh);

/**
 * @brief Interface to smcs_tev_tth()
 * @copydoc smcs_tev_tth()
 */
double SMCS_tev_ttH(double Mh);

/**
 * @brief Interface to smcs_lhc7_hw()
 * @copydoc smcs_lhc7_hw()
 */
double SMCS_lhc7_HW(double Mh);

/**
 * @brief Interface to smcs_lhc7_hz()
 * @copydoc smcs_lhc7_hz()
 */
double SMCS_lhc7_HZ(double Mh);

/**
 * @brief Interface to smcs_lhc7_gg_h()
 * @copydoc smcs_lhc7_gg_h()
 */
double SMCS_lhc7_gg_H(double Mh);

/**
 * @brief Interface to smcs_lhc7_bb_h()
 * @copydoc smcs_lhc7_bb_h()
 */
double SMCS_lhc7_bb_H(double Mh);

/**
 * @brief Interface to smcs_lhc7_vbf_h()
 * @copydoc smcs_lhc7_vbf_h()
 */
double SMCS_lhc7_vbf_H(double Mh);

/**
 * @brief Interface to smcs_lhc7_tth()
 * @copydoc smcs_lhc7_tth()
 */
double SMCS_lhc7_ttH(double Mh);

/**
 * @brief Interface to smcs_lhc8_hw()
 * @copydoc smcs_lhc8_hw()
 */
double SMCS_lhc8_HW(double Mh);

/**
 * @brief Interface to smcs_lhc8_hz()
 * @copydoc smcs_lhc8_hz()
 */
double SMCS_lhc8_HZ(double Mh);

/**
 * @brief Interface to smcs_lhc8_gg_h()
 * @copydoc smcs_lhc8_gg_h()
 */
double SMCS_lhc8_gg_H(double Mh);

/**
 * @brief Interface to smcs_lhc8_bb_h()
 * @copydoc smcs_lhc8_bb_h()
 */
double SMCS_lhc8_bb_H(double Mh);

/**
 * @brief Interface to smcs_lhc8_vbf_h()
 * @copydoc smcs_lhc8_vbf_h()
 */
double SMCS_lhc8_vbf_H(double Mh);

/**
 * @brief Interface to smcs_lhc8_tth()
 * @copydoc smcs_lhc8_tth()
 */
double SMCS_lhc8_ttH(double Mh);

/**
 * @brief Interface to smcs_lhc13_hw()
 * @copydoc smcs_lhc13_hw()
 */
double SMCS_lhc13_HW(double Mh);

/**
 * @brief Interface to smcs_lhc13_hz()
 * @copydoc smcs_lhc13_hz()
 */
double SMCS_lhc13_HZ(double Mh);

/**
 * @brief Interface to smcs_lhc13_gg_h()
 * @copydoc smcs_lhc13_gg_h()
 */
double SMCS_lhc13_gg_H(double Mh);

/**
 * @brief Interface to smcs_lhc13_bb_h()
 * @copydoc smcs_lhc13_bb_h()
 */
double SMCS_lhc13_bb_H(double Mh);

/**
 * @brief Interface to smcs_lhc13_vbf_h()
 * @copydoc smcs_lhc13_vbf_h()
 */
double SMCS_lhc13_vbf_H(double Mh);

/**
 * @brief Interface to smcs_lhc13_tth()
 * @copydoc smcs_lhc13_tth()
 */
double SMCS_lhc13_ttH(double Mh);

/**
 * @brief Interface to smcs_effc_hz()
 *
 *  **The parameter collider should be specified through the corresponding
 *  #HiggsBounds_collider_id.**
 *
 * @copydoc smcs_effc_hz()
 */
double SMCS_effC_HZ(double Mh, int collider, double ghZZ, double ghtt_s,
                    double ghbb_s, double ghtt_p, double ghbb_p);

/**
 * @brief Interface to smcs_effc_gg_hz()
 *
 *  **The parameter collider should be specified through the corresponding
 *  #HiggsBounds_collider_id.**
 *
 * @copydoc smcs_effc_gg_hz()
 */
double SMCS_effC_gg_HZ(double Mh, int collider, double ghZZ, double ghtt_s,
                       double ghbb_s, double ghtt_p, double ghbb_p);

/**
 * @brief Interface to smcs_effc_qq_hz()
 *
 *  **The parameter collider should be specified through the corresponding
 *  #HiggsBounds_collider_id.**
 *
 * @copydoc smcs_effc_qq_hz()
 */
double SMCS_effC_qq_HZ(double Mh, int collider, double ghZZ, double ghtt_s,
                       double ghbb_s, double ghtt_p, double ghbb_p);

/**
 * @brief Interface to smcs_effc_bb_hz()
 *
 *  **The parameter collider should be specified through the corresponding
 *  #HiggsBounds_collider_id.**
 *
 * @copydoc smcs_effc_bb_hz()
 */
double SMCS_effC_bb_HZ(double Mh, int collider, double ghbb_s, double ghbb_p);

/**
 * @brief Interface to smcs_effc_hw()
 *
 *  **The parameter collider should be specified through the corresponding
 *  #HiggsBounds_collider_id.**
 *
 * @copydoc smcs_effc_hw()
 */
double SMCS_effC_HW(double Mh, int collider, double ghWW, double ghtt_s,
                    double ghbb_s);

/**
 * @brief Interface to hccs_thc()
 * @copydoc hccs_thc()
 */
double HCCS_tHc(double MHc, double gHcjt, double gHcjb, double BR_tHpjb);

#ifdef __cplusplus
} // extern "C"
#endif
