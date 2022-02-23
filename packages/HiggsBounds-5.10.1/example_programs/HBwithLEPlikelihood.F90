!--------------------------------------------------------------------------------------
! This example program is part of HiggsBounds-5 (TS 23/02/2016).
!--------------------------------------------------------------------------------------
program HBwithLEPlikelihood
! This example program scans the SM Higgs mass from 100 GeV to 120 GeV and evaluates
! the usual HiggsBounds 95% C.L. limit as well as the chi^2 value from the LEP exclusions.
!--------------------------------------------------------------------------------------
    use theory_colliderSfunctions
    implicit none

    integer :: nHzero, nHplus
    integer :: HBresult, chan, ncombined, chan2
    integer, parameter :: fileid = 78
    double precision :: obsratio, theory_uncertainty_1s, chisq_withouttheory, chisq_withtheory
    double precision :: SMGamma_h, SMGammaTotal
    double precision :: Mh, ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
        ghjmumu_s, ghjmumu_p, ghjtautau_s, ghjtautau_p, &
        ghjWW, ghjZZ, ghjZga, ghjgaga, ghjgg, ghjhiZ
    character(len=100)::filename
    integer :: i

    nHzero = 1
    nHplus = 0
    theory_uncertainty_1s = 1.5D0

!--Setting up the output
    filename = 'HBchisq-output.dat'
    open (fileid, file=filename)
    write (fileid, *) '# mh   HBres   channel   obsratio   ncombined   chisq_wo_uncertainty ', &
        '  chisq_w_uncertainty channel(chisq)'
    write (fileid, *) '#--------------------------------------------------------------------', &
        '----------------------#'

!---- Initialize HiggsBounds with the LEP chisq tables ----!
    call initialize_HiggsBounds_chisqtables
!---- Initialize HiggsBounds with only LEP results ----!
    call initialize_HiggsBounds(nHzero, nHplus, "onlyL")

    do i = 1, 101
        Mh = 100.0D0 + (i - 1)*0.2D0

        SMGammaTotal = SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.
        if (.not. (SMGammaTotal .lt. 0)) then
            ghjss_s = 1.0d0
            ghjss_p = 0.0d0
            ghjcc_s = 1.0d0
            ghjcc_p = 0.0d0
            ghjbb_s = 1.0d0
            ghjbb_p = 0.0d0
            ghjtt_s = 1.0d0
            ghjtt_p = 0.0d0
            ghjmumu_s = 1.0d0
            ghjmumu_p = 0.0d0
            ghjtautau_s = 1.0d0
            ghjtautau_p = 0.0d0
            ghjWW = 1.0d0
            ghjZZ = 1.0d0
            ghjZga = 1d0
            ghjgg = 1.0d0
            ghjhiZ = 0d0
            ghjgaga = 1.0d0

            call HiggsBounds_neutral_input_properties(Mh, SMGammaTotal)

            call HiggsBounds_neutral_input_effC(ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                                                ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                                                ghjmumu_s, ghjmumu_p, &
                                                ghjtautau_s, ghjtautau_p, &
                                                ghjWW, ghjZZ, ghjZga, &
                                                ghjgaga, ghjgg, ghjhiZ)

!! n.b.: The LEP chi^2 extension requires a preceding HiggsBounds.
            call run_HiggsBounds_full(HBresult, chan, obsratio, ncombined)

            call HiggsBounds_get_LEPChisq(theory_uncertainty_1s, chisq_withouttheory, chisq_withtheory, chan2)

            write (fileid, *) Mh, HBresult, chan, obsratio, ncombined, &
                chisq_withouttheory, chisq_withtheory, chan2

        endif

    enddo

    close (fileid)

    call finish_HiggsBounds_chisqtables
    call finish_HiggsBounds

end program HBwithLEPlikelihood
