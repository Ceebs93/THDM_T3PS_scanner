/**
 * @file
 * @brief C interface to the HiggsSignals subroutines.
 *
 * Provides access to the most important subroutines from
 * HiggsSignals_subroutines.F90.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Interface to initialize_higgssignals_latestresults()
 * @copydoc initialize_higgssignals_latestresults
 */
void initialize_HiggsSignals_latestresults(int nHiggsneut, int nHiggsplus);

/**
 * @brief Interface to higgssignals_neutral_input_massuncertainty()
 * @copydoc higgssignals_neutral_input_massuncertainty
 */
void HiggsSignals_neutral_input_MassUncertainty(const double dMh[]);

/**
 * @brief Interface to setup_pdf()
 * @copydoc setup_pdf
 */
void HiggsSignals_setup_pdf(int pdf_in);

/**
 * @brief Interface to run_higgssignals_lhc_run1_combination()
 * @copydoc run_higgssignals_lhc_run1_combination
 */
void run_HiggsSignals_LHC_Run1_combination(double *Chisq_mu, double *Chisq_mh,
                                           double *Chisq, int *nobs,
                                           double *Pvalue);

/**
 * @brief Interface to run_higgssignals_stxs()
 * @copydoc run_higgssignals_stxs
 */
void run_HiggsSignals_STXS(double *Chisq_STXS_rates, double *Chisq_STXS_mh,
                           double *Chisq_STXS, int *nobs_STXS,
                           double *Pvalue_STXS);

/**
 * @brief Interface to run_higgssignals()
 * @copydoc run_higgssignals
 */
void run_HiggsSignals(double *Chisq_mu, double *Chisq_mh, double *Chisq,
                      int *nobs, double *Pvalue);

/**
 * @brief Interface to run_higgssignals_full()
 * @copydoc run_higgssignals_full
 */
void run_HiggsSignals_full(double *Chisq_mu, double *Chisq_mh, double *Chisq,
                           int *nobs, double *Pvalue);

/**
 * @brief Interface to get_rvalues()
 * @copydoc get_rvalues
 */
void get_HiggsSignals_Rvalues(int ii, int collider, double *R_H_WW,
                              double *R_H_ZZ, double *R_H_gaga,
                              double *R_H_tautau, double *R_H_bb,
                              double *R_VH_bb);

/**
 * @brief Interface to complete_hs_results()
 * @copydoc complete_hs_results
 */
void complete_HS_results();

/**
 * @brief Interface to finish_higgssignals()
 * @copydoc finish_higgssignals
 */
void finish_HiggsSignals();

#ifdef __cplusplus
}
#endif
