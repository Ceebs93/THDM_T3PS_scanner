subroutine initialize_HiggsSignals_latestresults_C(nHiggsneut, nHiggsplus) &
    bind(c, name='initialize_HiggsSignals_latestresults')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(in), value :: nHiggsneut, nHiggsplus
    call initialize_HiggsSignals(nHiggsneut, nHiggsplus, "latestresults")
end subroutine

subroutine HiggsSignals_neutral_input_MassUncertainty_C(dMh) &
    bind(c, name='HiggsSignals_neutral_input_MassUncertainty')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: dMh(np(Hneut))
    call HiggsSignals_neutral_input_MassUncertainty(dMh)
end subroutine

subroutine setup_pdf_C(pdf_in) &
    bind(c, name='HiggsSignals_setup_pdf')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(in), value :: pdf_in
    call setup_pdf(pdf_in)
end subroutine

subroutine run_HiggsSignals_LHC_Run1_combination_C(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue) &
    bind(c, name='run_HiggsSignals_LHC_Run1_combination')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: nobs
    real(kind=c_double), intent(out) :: Pvalue, Chisq, Chisq_mu, Chisq_mh
    call run_HiggsSignals_LHC_Run1_combination(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
end subroutine

subroutine run_HiggsSignals_STXS_C(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS) &
    bind(c, name='run_HiggsSignals_STXS')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: nobs_STXS
    real(kind=c_double), intent(out) :: Pvalue_STXS, Chisq_STXS, Chisq_STXS_rates, Chisq_STXS_mh
    call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)
end subroutine

subroutine run_HiggsSignals_C(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue) &
    bind(c, name='run_HiggsSignals')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: nobs
    real(kind=c_double), intent(out) :: Pvalue, Chisq, Chisq_mu, Chisq_mh
    call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
end subroutine

subroutine run_HiggsSignals_full_C(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue) &
    bind(c, name='run_HiggsSignals_full')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: nobs
    real(kind=c_double), intent(out) :: Pvalue, Chisq, Chisq_mu, Chisq_mh
    call run_HiggsSignals_full(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
end subroutine

subroutine get_Rvalues_C(ii, collider, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb) &
    bind(c, name='get_HiggsSignals_Rvalues')
    use, intrinsic :: iso_c_binding

    integer(kind=c_int), intent(in), value :: ii, collider
    real(kind=c_double), intent(out) :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb
    call get_Rvalues(ii, collider, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
end subroutine

subroutine complete_HS_results_C() bind(c, name='complete_HS_results')
    use, intrinsic :: iso_c_binding

    implicit none
    call complete_HS_results()
end subroutine

subroutine finish_HiggsSignals_C() bind(c, name = 'finish_HiggsSignals')
    use, intrinsic :: iso_c_binding

    implicit none
    call finish_HiggsSignals()
end subroutine
