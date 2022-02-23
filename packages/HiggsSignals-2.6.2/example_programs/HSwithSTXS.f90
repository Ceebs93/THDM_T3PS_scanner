!******************************************************
! This example program is part of HiggsSignals-2 (TS 09/08/2018).
!******************************************************
program HSwithSTXS
! This example program runs on a model with inclusive XS equal to the SM.
! We then modify the relative signal efficiency (r-value) of one STXS bin
! in the range [0,2].
!******************************************************

    use STXS
    use theory_BRfunctions
    implicit none

    double precision :: BR_hjss, BR_hjcc, BR_hjbb, &
        BR_hjtt, BR_hjmumu, BR_hjtautau, BR_hjWW, &
        BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg, GammaTotal, &
        chisq_tot, pval
    double precision :: scale

    integer :: i

    call initialize_HiggsSignals_empty(1, 0)

    call load_STXS("LHC13_CMS_H-ZZ_36.1fb01_Higgs-pT")

! Take SM predictions at 125 GeV
    GammaTotal = BRSM_GammaTot(125.0D0)
    BR_hjss = BRSM_Hss(125.0D0)
    BR_hjcc = BRSM_Hcc(125.0D0)
    BR_hjbb = BRSM_Hbb(125.0D0)
    BR_hjtt = BRSM_Htoptop(125.0D0)
    BR_hjmumu = BRSM_Hmumu(125.0D0)
    BR_hjtautau = BRSM_Htautau(125.0D0)
    BR_hjWW = BRSM_HWW(125.0D0)
    BR_hjZZ = BRSM_HZZ(125.0D0)
    BR_hjZga = BRSM_HZga(125.0D0)
    BR_hjgaga = BRSM_Hgaga(125.0D0)
    BR_hjgg = BRSM_Hgg(125.0D0)

    call HiggsBounds_neutral_input_properties(125.0D0, GammaTotal, 1)

!do i=11,11
!scale = 0.1*(i-1)

    scale = 1.0D0

    call HiggsBounds_neutral_input_hadr_single(13, "XS_hj_ratio", 1.0D0)
    call HiggsBounds_neutral_input_hadr_single(13, "XS_vbf_ratio", scale)
    call HiggsBounds_neutral_input_hadr_single(13, "XS_hjZ_ratio", scale)
    call HiggsBounds_neutral_input_hadr_single(13, "XS_hjW_ratio", scale)
    call HiggsBounds_neutral_input_hadr_single(13, "XS_tthj_ratio", scale)

    call HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, &
                                        BR_hjtt, BR_hjmumu, &
                                        BR_hjtautau, BR_hjWW, &
                                        BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                        BR_hjgg)

    write (*, *) "Modifying r-value of STXS observable with ID = 99365 (bin 5):"
    do i = 1, 21
        call assign_modification_factor_to_STXS(99365, 5, (/0.0D0 + 0.1D0*(i - 1), 1.0D0, 1.0D0, 1.0D0, 1.0D0/))
        call calculate_model_predictions_for_STXS
        call get_chisq_from_STXS(chisq_tot, pval)
        write (*, *) "r-value, chisq, pval = ", 0.0D0 + 0.1D0*(i - 1), chisq_tot, pval
    enddo

! call load_STXS("STXS_ATL_H-ZZ_13TeV")
! call load_STXS("LHC13")

!  call calculate_model_predictions_for_STXS
!  call get_chisq_from_STXS(chisq_tot, pval)
!  write(*,*) "chisq, pval = ", chisq_tot, pval

!enddo

! call print_STXS()

end program HSwithSTXS
