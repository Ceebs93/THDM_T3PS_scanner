!******************************************************
program HBwithchannelrates
!
!******************************************************
    use theory_XS_SM_functions
    use theory_BRfunctions
    use channels, only: HiggsBounds_deactivate_analyses, HiggsBounds_activate_all_analyses
    use output, only: createKey
    implicit none

    integer :: nHzero, nHplus, i, CP_value(3)
    double precision :: Mh(3), GammaTotal(3), CS_hjhi(3, 3), &
        CS_hj_ratio(3), CS_gg_hj_ratio(3), CS_bb_hj_ratio(3), &
        CS_hjW_ratio(3), CS_hjZ_ratio(3), CS_thj_schan_ratio(3), &
        CS_vbf_ratio(3), CS_tthj_ratio(3), CS_thj_tchan_ratio(3), &
        CS_qq_hjZ_ratio(3), CS_gg_hjZ_ratio(3), &
        CS_tWhj_ratio(3), &
        BR_hjss(3), BR_hjcc(3), &
        BR_hjbb(3), BR_hjtt(3), &
        BR_hjmumu(3), &
        BR_hjtautau(3), &
        BR_hjWW(3), BR_hjZZ(3), BR_hjZga(3), BR_hjgaga(3), &
        BR_hjgg(3)
    double precision :: ggH(3), bbH(3), tb
    character(len=100)::filename_in, filename_out
    integer :: HBresult0, chan0, ncombined0, HBresult1, chan1, ncombined1
    double precision :: obsratio0, obsratio1, llh0, llh1
    integer :: error, status, n, cbin
    double precision :: M_av
    integer ::  Hindex, nc
    double precision :: interfere_mod, totcs

    nHzero = 3
    nHplus = 0

    interfere_mod = -0.05D0

    call initialize_HiggsBounds(nHzero, nHplus, 'onlyH')

! Optionally, deactivate the 95%CL limit extraction for CMS MSSM h/H/A->tautau, since
! we want to use the likelihood information instead, as well as all relevant previous
! results from non-standard Higgs to tau tau searches from ATLAS/CMS:

! call HiggsBounds_deactivate_analyses((/3316, 2014049, 20140492/))

! Note: Deactivation of analyses can be changed before every HiggsBounds run
! (as currently done, see below). If all analyses need to be activated again, just call
! call HiggsBounds_activate_all_analyses

    filename_in = "../example_data/Mh125/mh125_8.tsv"
    filename_out = "Mh125_HBwithchannelrates.dat"

    call system('rm -f Mh125_HBwithchannelrates.dat')

    open (432, file=trim(adjustl(filename_in)), action='read', status='old', iostat=status)
    open (433, file=trim(adjustl(filename_out)), action='write', status='new')

    if (status .ne. 0) then
        write (*, *) 'Bad status', status, 'with the following file:'
        write (*, *) trim(adjustl(filename_in))
        stop
    endif
    read(432, *) ! skip header

    n = 0
    do
! Read in the relevant cross section / BR predictions from the data grid for 8 TeV.
        read (432, *, iostat=error) Mh, tb, ggH, bbH, BR_hjtautau
        if (error == -1) exit
        n = n + 1
!----
! QUICK STOP FOR TESTING
!    if (n.le.0) cycle
!    if (n.ge.5) exit

        if (n .le. 1000) cycle
        if (n .ge. 2001) exit
!----
        if (mod(n, 100) .eq. 0) write (*, *) "number of processed points: ", n, " MA,TB = ", Mh(3), tb

        CP_value(1) = 1
        CP_value(2) = 1
        CP_value(3) = -1

        do i = 1, 3
            GammaTotal(i) = BRSM_GammaTot(Mh(i))
! Normalize the cross sections to the SM predictions using the 8 TeV predictions
! (using HiggsBounds SM cross section routines)
            CS_gg_hj_ratio(i) = ggH(i)/XS_lhc8_gg_H_SM(Mh(i))
            CS_bb_hj_ratio(i) = bbH(i)/XS_lhc8_bb_H_SM(Mh(i))
        enddo
! Set the other predictions to zero here.
        CS_hjW_ratio = 0.0D0
        CS_hjZ_ratio = 0.0D0
        CS_vbf_ratio = 0.0D0
        CS_tthj_ratio = 0.0D0
        CS_thj_tchan_ratio = 0.0D0
        CS_thj_schan_ratio = 0.0D0
        CS_qq_hjZ_ratio = 0.0D0
        CS_gg_hjZ_ratio = 0.0D0
        CS_tWhj_ratio = 0.D0
        CS_hjhi = 0.0D0
        BR_hjss = 0.0D0
        BR_hjcc = 0.0D0
        BR_hjbb = 0.0D0
        BR_hjtt = 0.0D0
        BR_hjmumu = 0.0D0
        BR_hjWW = 0.0D0
        BR_hjZZ = 0.0D0
        BR_hjZga = 0.0D0
        BR_hjgaga = 0.0D0
        BR_hjgg = 0.0D0
!           BR_hjinvisible=0.0D0
!       BR_hjhihi_nHbynH=0.0D0

        call HiggsBounds_neutral_input_properties(Mh, GammaTotal, CP_value)

! Set 8 TeV cross section input
        call HiggsBounds_neutral_input_hadr(8, CS_hj_ratio, &
                                            CS_gg_hj_ratio, CS_bb_hj_ratio, &
                                            CS_hjW_ratio, CS_hjZ_ratio, &
                                            CS_vbf_ratio, CS_tthj_ratio, &
                                            CS_thj_tchan_ratio, CS_thj_schan_ratio, &
                                            CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, &
                                            CS_tWhj_ratio, &
                                            CS_hjhi)

! Set 13 TeV cross section input (assume SM-normalized XS are same as for 8 TeV)
        call HiggsBounds_neutral_input_hadr(13, CS_hj_ratio, &
                                            CS_gg_hj_ratio, CS_bb_hj_ratio, &
                                            CS_hjW_ratio, CS_hjZ_ratio, &
                                            CS_vbf_ratio, CS_tthj_ratio, &
                                            CS_thj_tchan_ratio, CS_thj_schan_ratio, &
                                            CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, &
                                            CS_tWhj_ratio, &
                                            CS_hjhi)

! Set BR input for SM final states
        call HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, &
                                            BR_hjtt, BR_hjmumu, &
                                            BR_hjtautau, BR_hjWW, &
                                            BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                            BR_hjgg)

! Run the standard HiggsBounds routine considering all analyses
        call run_HiggsBounds(HBresult0, chan0, obsratio0, ncombined0)

! Get observed likelihood
        call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh0, 'obs')

! Calculate total cross section (normalized) for H/A include modification from interference.
! Then give HB zero cross section for H, and total H/A (including modification) for A.

        write (*, *) CS_gg_hj_ratio(2), CS_bb_hj_ratio(2), BR_hjtautau(2)
        write (*, *) CS_gg_hj_ratio(3), CS_bb_hj_ratio(3), BR_hjtautau(3)

! First for gg->phi:
        totcs = (CS_gg_hj_ratio(2)*BR_hjtautau(2) + &
                 CS_gg_hj_ratio(3)*BR_hjtautau(3))*(1.0D0 + interfere_mod)

        call HiggsBounds_neutral_input_hadr_channelrates_single(8, 2, 6, 4, 0.0D0)
        call HiggsBounds_neutral_input_hadr_channelrates_single(8, 3, 6, 4, totcs)

        call HiggsBounds_neutral_input_hadr_channelrates_single(13, 2, 6, 4, 0.0D0)
        call HiggsBounds_neutral_input_hadr_channelrates_single(13, 3, 6, 4, totcs)

        write (*, *) "sigma(gg->H/A->tautau)/sigma_SM(gg->H): ", totcs
! Then for bb->phi:
        totcs = (CS_bb_hj_ratio(2)*BR_hjtautau(2) + &
                 CS_bb_hj_ratio(3)*BR_hjtautau(3))*(1.0D0 + interfere_mod)

! Enter zero cross section for H, and total H/A (including modification) for A:
        call HiggsBounds_neutral_input_hadr_channelrates_single(8, 2, 7, 4, 0.0D0)
        call HiggsBounds_neutral_input_hadr_channelrates_single(8, 3, 7, 4, totcs)

        call HiggsBounds_neutral_input_hadr_channelrates_single(13, 2, 7, 4, 0.0D0)
        call HiggsBounds_neutral_input_hadr_channelrates_single(13, 3, 7, 4, totcs)

        write (*, *) "sigma(bb->H/A->tautau)/sigma_SM(bb->H): ", totcs

! Run the standard HiggsBounds routine considering all analyses
        call run_HiggsBounds(HBresult1, chan1, obsratio1, ncombined1)

! Get observed likelihood
        call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh1, 'obs')

        write (433, *) n, Mh, tb, HBresult0, chan0, obsratio0, ncombined0, llh0, &
            HBresult1, chan1, obsratio1, ncombined1, llh1
    enddo

!--------------------------------------------------------
! Example for the auxiliary chi^2 functions (run on the last point of the scan)
!--------------------------------------------------------
! Get observed likelihood for lightest Higgs

!       write(*,*) "#------------------------------------------------------------------#"
!       call HiggsBounds_get_likelihood_for_Higgs(14029, 0, 1, nc, cbin, M_av, llh, 'obs')
!       write(*,*) "The observed likelihood value for the light Higgs h is ",llh
!       write(*,*) "The binary code and number of Higgses of the formed combination is ", cbin, nc
!       write(*,*) "The likelihood has been evaluated at an average mass value of ", M_av
!       write(*,*) "#------------------------------------------------------------------#"

! Get observed likelihood for a subset of available Higgs bosons
! (here, e.g., exclude h and H from possible combination)

!       call HiggsBounds_get_likelihood_for_comb(14029, 3, Hindex, nc, cbin, M_av, llh, 'obs')
!       write(*,*) "The observed likelihood value is ",llh
!       write(*,*) "The binary code, number of Higgses and Higgs index of the formed combination is ", cbin, nc, Hindex
!       write(*,*) "The likelihood has been evaluated at an average mass value of ", M_av
!       write(*,*) "#------------------------------------------------------------------#"
!--------------------------------------------------------

! Write out the key with the used analyses (indicating possibly deactivated analyses)
!    call createKey("HB_with_deactivated_analyses_")
    close (432)
    close (433)

    call finish_HiggsBounds
end program HBwithchannelrates
