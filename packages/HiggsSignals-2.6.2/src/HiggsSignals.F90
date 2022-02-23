!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 31/01/2013)
!--------------------------------------------------------------------
program HiggsSignals
!--------------------------------------------------------------------
! Creates the command-line executable of HiggsSignals
!
! HiggsSignals performs a chi^2 test of Higgs signal rate and mass
! predictions of an arbitrary Higgs sector with measurements from
! the Tevatron and LHC experiments.
!--------------------------------------------------------------------
    use usefulbits, only: inputmethod, np, Hneut, Hplus, whichanalyses, vsmall
    use input, only: do_input
    use usefulbits_hs, only: output_level, HiggsSignals_info, Exptdir, Nparam, HSres
    use io, only: setup_input_for_hs, do_output_for_hs
    use STXS, only: calculate_model_predictions_for_STXS, get_chisq_from_STXS, &
                    get_number_of_STXS_observables, print_STXS
    use numerics, only: gammp

    implicit none

    double precision :: Pvalue, Chisq, Chisq_mu, Chisq_mh
    double precision :: Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, Chisq_LHCRun1, Pvalue_LHCRun1
    double precision :: Chisq_STXS, Chisq_STXS_rates, Chisq_STXS_mh, Pvalue_STXS
    double precision :: Pvalue_tot
    integer :: nobs, nobs_LHCRun1, nobs_STXS, mode

    inputmethod = 'datfile'
    output_level = 0
    whichanalyses = 'onlyH'

    call setup_input_for_hs

    call initialize_HiggsSignals(np(Hneut), np(Hplus), Exptdir)

    call do_input

!  write(*,*) "Run HS..."
    call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)

!  write(*,*) "Run HS Run1 ..."
    call run_HiggsSignals_LHC_Run1_combination(Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1)

!  write(*,*) "Run HS STXS ..."
    call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)

!  call calculate_model_predictions_for_STXS

!  call get_chisq_from_STXS(Chisq_STXS, Pvalue_STXS)
!
!  call get_number_of_STXS_observables(nobs_STXS)

!  call print_STXS()

    if ((Chisq + Chisq_LHCRun1 + Chisq_STXS) .gt. vsmall .and. (nobs + nobs_LHCRun1 + nobs_STXS - Nparam) .gt. 0) then
        Pvalue_tot = 1 - gammp(dble(nobs + nobs_LHCRun1 + nobs_STXS - Nparam)/2.0D0, (Chisq + Chisq_LHCRun1 + Chisq_STXS)/2.0D0)
    endif

    call complete_HS_results()

!  write(*,*) "Chisq_mu, Chisq_mu_LHCRun1, Chisq_STXS = ",Chisq_mu,Chisq_mu_LHCRun1,Chisq_STXS
!  write(*,*) "Chisq_mh, Chisq_mh_LHCRun1 = ",Chisq_mh,Chisq_mh_LHCRun1
!  write(*,*) "nobs, nobs_LHCRun1, nobs_STXS, Nparam = ", nobs, nobs_LHCRun1, nobs_STXS, Nparam
!  write(*,*) "Pvalue_tot = ", Pvalue_tot

    call do_output_for_hs

    call finish_HiggsSignals

end program
