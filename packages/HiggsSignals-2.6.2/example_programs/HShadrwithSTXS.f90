!******************************************************
! This example program is part of HiggsSignals-2 (TS 29/03/2017).
!******************************************************
program HShadrwithSTXS
! This example program uses the hadronic cross section input format
! to scan over the two scale factors:
!   scale_ggf scales the SM single Higgs production cross section
!   scale_VBF scales the SM VBF production cross section
! The remaining cross sections (ttH, VH) are fixed to their SM values.
! We use a set of Simplified Template Cross Section (STXS) measurements
! for the chi-2 test.
!
! The output is written into /results/HShadrwithSTXS.dat, which can be plotted with the
! python script plot_HShadr.py in the results folder.
!******************************************************
    use STXS
    implicit none

    integer :: nHzero, nHplus, ndf, i, j, k, ii, jj, CP_value
    double precision :: Pvalue, Chisq, Chisq_mu, Chisq_mh, scale_ggf, scale_VBF
    double precision :: SMGammaTotal, SMGamma_h
    double precision :: SMBR_Htoptop, SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Htt, &
        SMBR_Hmumu, SMBR_Htautau, SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg
! Entries of CS arrays: TEV, LHC7, LHC8, LHC13
    double precision :: Mh, GammaTotal, CS_hj_ratio(4), &
        CS_gg_hj_ratio(4), CS_bb_hj_ratio(4), &
        CS_hjW_ratio(4), CS_hjZ_ratio(4), &
        CS_vbf_ratio(4), CS_tthj_ratio(4), &
        CS_hjhi(4), CS_thj_schan_ratio(4), &
        CS_thj_tchan_ratio(4), &
        CS_qq_hjZ_ratio(4), CS_gg_hjZ_ratio(4), &
        BR_hjss, BR_hjcc, &
        BR_hjbb, BR_hjtt, &
        BR_hjmumu, &
        BR_hjtautau, &
        BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, &
        BR_hjgg
    character(len=100)::filename
    double precision :: dm
    integer :: collider, collider_s

    nHzero = 1
    nHplus = 0

    dm = 0.0D0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
    call initialize_HiggsSignals_empty(nHzero, nHplus)
    call load_STXS("LHC13")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)    ----!
!  call setup_output_level(0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)        ----!
!  call setup_pdf(2)
    !---- Set the assignment range for the peak-centered method (optional)             ----!
!  call setup_assignmentrange_massobservables(4.0D0)
!---- Pass the Higgs mass uncertainty to HiggsSignals                                ----!
!  call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set number of free model parameters                                            ----!
    call setup_Nparam(2)
!---- Open output text file                                                          ----!
    filename = 'results/HShadrwithSTXS.dat'
    open (21, file=filename)
    write (21, *) '# Mh        scale_ggf        scale_VBF        Chisq        Pvalue'
    write (21, *) '#--------------------------------------------------------------------------'

    do i = 1, 81
        do j = 1, 201
            scale_ggf = 0.0D0 + (i - 1)*0.025D0
            scale_VBF = 0.0D0 + (j - 1)*0.05D0

!  do i=1,2
!   do j=1,2
!    scale_ggf = 1.0D0 +(i-1)*1.D0
!    scale_VBF = 1.0D0 +(j-1)*1.D0

            write (*, *) "kappa_ggf , kappa_VBF = ", scale_ggf, scale_VBF

            Mh = 125.09D0
            SMGammaTotal = SMGamma_h(Mh)

            if (.not. (SMGammaTotal .lt. 0)) then
                GammaTotal = SMGammaTotal
! CP even
                CP_value = 1
! This applies to all 4 elements:
                CS_hj_ratio = 1.0D0*scale_ggf
                CS_gg_hj_ratio = 1.0D0*scale_ggf
                CS_bb_hj_ratio = 1.0D0*scale_ggf
                CS_hjW_ratio = 1.0D0
                CS_hjZ_ratio = 1.0D0
                CS_vbf_ratio = 1.0D0*scale_VBF
                CS_tthj_ratio = 1.0D0
                CS_hjhi = 0.0D0
                CS_thj_tchan_ratio = 1.0D0
                CS_thj_schan_ratio = 1.0D0
                CS_qq_hjZ_ratio = 1.0D0
                CS_gg_hjZ_ratio = 1.0D0

                BR_hjss = SMBR_Hss(Mh)
                BR_hjcc = SMBR_Hcc(Mh)
                BR_hjbb = SMBR_Hbb(Mh)
                BR_hjmumu = SMBR_Hmumu(Mh)
                BR_hjtautau = SMBR_Htautau(Mh)
                BR_hjWW = SMBR_HWW(Mh)
                BR_hjZZ = SMBR_HZZ(Mh)
                BR_hjZga = SMBR_HZgam(Mh)
                BR_hjgaga = SMBR_Hgamgam(Mh)
                BR_hjgg = SMBR_Hgg(Mh)

                call HiggsBounds_neutral_input_properties(Mh, GammaTotal, CP_value)

                do collider = 1, 4
                    select case (collider)
                    case (1)
                        collider_s = 2
                    case (2)
                        collider_s = 7
                    case (3)
                        collider_s = 8
                    case (4)
                        collider_s = 13
                    end select

                    call HiggsBounds_neutral_input_hadr(collider_s, CS_hj_ratio(collider), &
                                                        CS_gg_hj_ratio(collider), CS_bb_hj_ratio(collider), &
                                                        CS_hjW_ratio(collider), CS_hjZ_ratio(collider), &
                                                        CS_vbf_ratio(collider), CS_tthj_ratio(collider), &
                                                        CS_thj_tchan_ratio(collider), CS_thj_schan_ratio(collider), &
                                                        CS_qq_hjZ_ratio(collider), CS_gg_hjZ_ratio(collider), &
                                                        CS_hjhi(collider))

                enddo

                call HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, &
                                                    BR_hjmumu, BR_hjtautau, BR_hjWW, BR_hjZZ, &
                                                    BR_hjZga, BR_hjgaga, BR_hjgg)

                call calculate_model_predictions_for_STXS()

                call get_chisq_from_STXS(Chisq, Pvalue)

!     call print_STXS()

!      call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

                write (21, *) Mh, scale_ggf, scale_VBF, Chisq, Pvalue

            endif
        enddo
    enddo

    close (21)

    call finish_HiggsSignals

end program HShadrwithSTXS
