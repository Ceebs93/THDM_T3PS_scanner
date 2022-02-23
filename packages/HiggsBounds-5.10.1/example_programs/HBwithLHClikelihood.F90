!******************************************************
program HBwithLHClikelihood
!
! In this example we evaluate the exclusion likelihoods ATLAS and CMS Higgs searches
! with tautau final states for the mhmod+ scenario. The input is provided from two
! datafiles (8 TeV and 13 Tev) created with SusHi, obtained here:
! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMNeutral
!
! We also demonstrate how certain analyses can be deactivated in the standard
! HiggsBounds run. In particular, we deactivate here the latest ATLAS and CMS tau tau
! 95% CL limits from the standard HiggsBounds procedure, as we want to use the
! likelihood instead.
!
! After a successful run the output can be plotted with the python script
! plot_mhmodp_llh.py (needs matplotlib package!)
!
! (TS 24/11/2017)
!******************************************************

    use theory_XS_SM_functions
    use theory_BRfunctions
    use channels, only: HiggsBounds_deactivate_analyses, HiggsBounds_activate_all_analyses
    use output, only: createKey
    implicit none

    integer :: nHzero, nHplus, i, CP_value(3)
    double precision :: Mh(3), GammaTotal(3), CS_gg_hj_ratio(3), CS_bb_hj_ratio(3), &
        BR_hjss(3), BR_hjcc(3), BR_hjbb(3), BR_hjtt(3), BR_hjmumu(3), BR_hjtautau(3), &
        BR_hjgg(3), BR_hjWW(3), BR_hjZZ(3), BR_hjZga(3), BR_hjgaga(3)
    double precision :: ggH8TeV(3), bbH8TeV(3), ggH13TeV(3), bbH13TeV(3), tb
    character(len=100)::filename_in_8TeV, filename_in_13TeV, filename_out
    integer :: HBresult, chan, ncombined, HBresult_all, chan_all, ncombined_all
    double precision :: obsratio, obsratio_all
    integer :: error, status, n, cbin
    double precision :: M_av, llh_CMS8, llh_CMS13, llh_ATLAS13, llh_ATLAS20, &
        llh_exp_CMS8, llh_exp_CMS13, llh_exp_ATLAS13, llh_exp_ATLAS20
    integer ::  Hindex, nc

    nHzero = 3
    nHplus = 0

    call initialize_HiggsBounds(nHzero, nHplus, 'onlyH')

    filename_in_8TeV = "../example_data/Mh125/mh125_8.tsv"
    filename_in_13TeV = "../example_data/Mh125/mh125_13.tsv"
    filename_out = "Mh125_HBwithLHClikelihood.dat"
    call system('rm -f Mh125_HBwithLHClikelihood.dat')



    open (432, file=trim(adjustl(filename_in_8TeV)), action='read', status='old', iostat=status)
    if (status .ne. 0) then
        write (*, *) 'Bad status', status, 'with the following file:'
        write (*, *) trim(adjustl(filename_in_8TeV))
        stop
    endif
    read(432, *) ! skip header

    open (433, file=trim(adjustl(filename_in_13TeV)), action='read', status='old', iostat=status)
    if (status .ne. 0) then
        write (*, *) 'Bad status', status, 'with the following file:'
        write (*, *) trim(adjustl(filename_in_13TeV))
        stop
    endif
    read(433, *) ! skip header

    open (434, file=trim(adjustl(filename_out)), action='write', status='new')

    n = 0
    do
! Read in the relevant cross section and BR predictions from the data grid for 8 TeV.
! The grid contains values for MA <= 1 TeV. Set values to zero for MA values beyond.
        read (432, *, iostat=error)  Mh, tb,  ggH8TeV,  bbH8TeV, BR_hjtautau
        ! if (error == -1) then
        !     ggH8TeV = (/0.0D0, 0.0D0, 0.0D0/)
        !     bbH8TeV = (/0.0D0, 0.0D0, 0.0D0/)
        ! endif
! Read in the relevant cross section and BR predictions from the data grid for 13 TeV.
        read (433, *, iostat=error)  Mh, tb,  ggH13TeV,  bbH13TeV, BR_hjtautau
        if (error == -1) exit

        n = n + 1
!----
! QUICK STOP FOR TESTING
!    if (n.le.141) cycle
        ! if (n .ge. 501) exit
!----
        if (mod(n, 100) .eq. 0) write (*, *) "number of processed points: ", n, " MA,TB = ", Mh(3), tb

! This is technically not correct, but not needed here anyways:
        do i = 1, 3
            GammaTotal(i) = BRSM_GammaTot(Mh(i))
        enddo
! Also, not really needed here:
        CP_value(1) = 1
        CP_value(2) = 1
        CP_value(3) = -1

! Normalize the cross sections to the SM predictions using the 8 TeV predictions
! (using HiggsBounds SM cross section routines)
        do i = 1, 3
            CS_gg_hj_ratio(i) = ggH8TeV(i)/XS_lhc8_gg_H_SM(Mh(i))
            CS_bb_hj_ratio(i) = bbH8TeV(i)/XS_lhc8_bb_H_SM(Mh(i))
        enddo

        ! Set 8 TeV cross section input for gg->phi and bb->phi
        call HiggsBounds_neutral_input_hadr_single(8, "XS_gg_hj_ratio", CS_gg_hj_ratio)
        call HiggsBounds_neutral_input_hadr_single(8, "XS_bb_hj_ratio", CS_bb_hj_ratio)

! Normalize the cross sections to the SM predictions using the 13 TeV predictions
! (using HiggsBounds SM cross section routines)
        do i = 1, 3
            CS_gg_hj_ratio(i) = ggH13TeV(i)/XS_lhc13_gg_H_SM(Mh(i))
            CS_bb_hj_ratio(i) = bbH13TeV(i)/XS_lhc13_bb_H_SM(Mh(i))
        enddo

        ! Set 8 TeV cross section input
        call HiggsBounds_neutral_input_hadr_single(13, "XS_gg_hj_ratio", CS_gg_hj_ratio)
        call HiggsBounds_neutral_input_hadr_single(13, "XS_bb_hj_ratio", CS_bb_hj_ratio)

! Set the other predictions to zero here.
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

        call HiggsBounds_neutral_input_properties(Mh, GammaTotal, CP_value)

! Set BR input for SM final states
        call HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, &
                                            BR_hjtt, BR_hjmumu, &
                                            BR_hjtautau, BR_hjWW, &
                                            BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                            BR_hjgg)

! Activate all analyses (in case some of them have been deactivated before)
        call HiggsBounds_activate_all_analyses

! Run the standard HiggsBounds routine considering all analyses
        call run_HiggsBounds(HBresult_all, chan_all, obsratio_all, ncombined_all)

! Deactivate current CMS and ATLAS searches for non-standard Higgs to tautau
        call HiggsBounds_deactivate_analyses((/14029, 2014049, 20140492, 170907242, 200212223/))
! Standard HiggsBounds run (gives 95% CL limit):
        call run_HiggsBounds(HBresult, chan, obsratio, ncombined)

! Obtain exclusion-likelihood from H->tautau searches:
!
!  The arguments of the following subroutines mean the following:
!        Analysis-ID --> 14029 for CMS 8 TeV results,
!                   17020 for CMS 13 TeV results,
!                   170907242 for ATLAS 13 TeV results.
!                   200212223 for new ATLAS 13 TeV results
!   Hindex  -> Index of the Higgs boson that was selected as most sensitive [int output]
!   M_av    -> mass position where limit is extracted (a signal strength-weighted mass
!              average in case of combined Higgs bosons) [dbl output]
!   nc      -> number of combined Higgs bosons [int output]
!   cbin    -> binary code of the combined Higgs bosons (see manual) [int output]
!   llh     -> -2ln L value [dbl output]
!   obspred -> 'obs' or 'pred' to chose whether the observed or expected likelihood should be
!              extracted. [char input]
!
! Get expected/predicted likelihood
        call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh_exp_CMS8, 'pred')
! ! Get observed likelihood
        call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh_CMS8, 'obs')

! Get expected/predicted likelihood
        call HiggsBounds_get_likelihood(17020, Hindex, nc, cbin, M_av, llh_exp_CMS13, 'pred')
! Get observed likelihood
        call HiggsBounds_get_likelihood(17020, Hindex, nc, cbin, M_av, llh_CMS13, 'obs')

! Get expected/predicted likelihood
        call HiggsBounds_get_likelihood(170907242, Hindex, nc, cbin, M_av, llh_exp_ATLAS13, 'pred')
! Get observed likelihood
        call HiggsBounds_get_likelihood(170907242, Hindex, nc, cbin, M_av, llh_ATLAS13, 'obs')

! Get expected/predicted likelihood
        call HiggsBounds_get_likelihood(200212223, Hindex, nc, cbin, M_av, llh_exp_ATLAS20, 'pred')
! Get observed likelihood
        call HiggsBounds_get_likelihood(200212223, Hindex, nc, cbin, M_av, llh_ATLAS20, 'obs')

        write (434, *) n, Mh, tb, HBresult, chan, obsratio, ncombined, &
            HBresult_all, chan_all, obsratio_all, ncombined_all, &
            Hindex, M_av, nc, cbin, llh_CMS8, llh_exp_CMS8, llh_CMS13, &
            llh_exp_CMS13, llh_ATLAS13, llh_exp_ATLAS13, llh_ATLAS20, llh_exp_ATLAS20

    enddo

! Write out the key with the used analyses (indicating possibly deactivated analyses)
    call createKey("HB_with_deactivated_analyses_")
    close (432)
    close (433)
    close (434)

    call finish_HiggsBounds
end program HBwithLHClikelihood
