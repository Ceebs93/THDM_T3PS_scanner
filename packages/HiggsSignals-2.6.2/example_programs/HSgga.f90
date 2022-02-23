!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals-2 (TS 29/03/2017).
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
program HSgga
!--------------------------------------------------------------------------------------
! In this example we scan over two kappa scale factors, kappag and kappaga, of the
! 125 GeV Higgs boson, using the effective couplings input.
!
! The output is written into /results/HSgga.dat, which can be plotted with the
! python script plot_HSgga.py in the results folder.
!--------------------------------------------------------------------------------------
    use theory_colliderSfunctions
    use usefulbits, only: vsmall
    implicit none

    integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
    double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
    double precision :: SMGammaTotal
    double precision :: kappag, kappaga
    double precision :: Mh, GammaTotal, ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
        ghjmumu_s, ghjmumu_p, ghjtautau_s, ghjtautau_p, &
        ghjWW, ghjZZ, ghjZga, ghjgaga, ghjgg, &
        ghjhiZ
    character(len=100)::filename
    double precision :: dm
    integer                  :: pdf
!-HiggsBounds internal functions to obtain SM branching ratios
    double precision :: SMBR_Htoptop, SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau, &
        SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg, SMGamma_h
    double precision :: Htogaga_rate, HtoVV_rate, HtoFF_rate
    nHzero = 1
    nHplus = 0

!--Setting up the output
    filename = 'results/HSgga.dat'
    open (21, file=filename)
    write (21, *) '# mh kappag kappaga Chisq_mu Chisq ndf Htogaga_rate HtoVV_rate HtoFF_rate'
    write (21, *) '#--------------------------------------------------------------------#'

!--Enter the Higgs mass and its theory uncertainty here:
    Mh = 125.09D0
    dm = 0.0D0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
    call initialize_HiggsSignals(nHzero, nHplus, "LHC13_ATL_Hgaga_STXS")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)    ----!
    call setup_output_level(0)
    !---- Set the assignment range for the peak-centered method (optional)              ----!
    call setup_assignmentrange_massobservables(4.0D0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)        ----!
    pdf = 2
    call setup_pdf(pdf)
!---- Pass the Higgs mass uncertainty to HiggsSignals (if relevant)                  ----!
! call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set number of free model parameters ----!
    call setup_Nparam(2)

    do i = 1, 81
        do j = 1, 81
            kappag = 0.0D0 + (i - 1)*0.025D0
            kappaga = 0.0D0 + (j - 1)*0.025D0

            SMGammaTotal = SMGamma_h(Mh)

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
                ghjZga = 1.0d0
                ghjgg = kappag
                ghjhiZ = 0d0
                ghjgaga = kappaga

!----Calculate the new total decay width:
                GammaTotal = SMGammaTotal*(1 + &
                                           (ghjWW**2.0 - 1)*SMBR_HWW(Mh) + (ghjZZ**2.0 - 1)*SMBR_HZZ(Mh) + &
                                           (ghjgg**2.0 - 1)*SMBR_Hgg(Mh) + (ghjtt_s**2.0 - 1)*SMBR_Htoptop(Mh) + &
                                           (ghjbb_s**2.0 - 1)*SMBR_Hbb(Mh) + (ghjtautau_s**2.0 - 1)*SMBR_Htautau(Mh) + &
                                           (ghjss_s**2.0 - 1)*SMBR_Hss(Mh) + (ghjcc_s**2.0 - 1)*SMBR_Hcc(Mh) + &
                                           (ghjZga**2.0 - 1)*SMBR_HZgam(Mh) + (ghjmumu_s**2.0 - 1)*SMBR_Hmumu(Mh) + &
                                           (ghjgaga**2.0 - 1)*SMBR_Hgamgam(Mh))

                call HiggsBounds_neutral_input_properties(Mh, GammaTotal)

                call HiggsBounds_neutral_input_effC( &
                    ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                    ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                    ghjmumu_s, ghjmumu_p, &
                    ghjtautau_s, ghjtautau_p, &
                    ghjWW, ghjZZ, ghjZga, &
                    ghjgaga, ghjgg, ghjhiZ)

                call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

! This will get the SM normalized rates for inclusive Higgs production,
! with H-> gamma gamma, VV and FF decays:
                call get_rates(1, 4, 5, (/"1.1", "2.1", "3.1", "4.1", "5.1"/), Htogaga_rate)
                call get_rates(1, 4, 5, (/"1.2", "2.2", "3.2", "4.2", "5.2"/), HtoVV_rate)
                call get_rates(1, 4, 5, (/"1.4", "2.4", "3.4", "4.4", "5.4"/), HtoFF_rate)

! This will collect the main HiggsSignals results together into one file
                write (21, *) mh, kappag, kappaga, Chisq_mu, Chisq, ndf, Htogaga_rate, HtoVV_rate, HtoFF_rate

            endif
        enddo
    enddo

    close (21)

    write (*, *) "Finishing HiggsSignals..."
    call finish_HiggsSignals

contains

!**************************************************************
    function get_g2hgaga(ghbb, ghtt, ghtautau, ghWW, ghZZ)
! Evaluates g2hgaga from other effective couplings, using partial widths informations
! at a Higgs mass of 126 GeV (calculated with HDECAY and taken from
! http://people.web.psi.ch/spira/higgscoup/ ).
!**************************************************************
        double precision, intent(in) :: ghbb, ghtt, ghtautau, ghWW, ghZZ
        double precision :: get_g2hgaga

        get_g2hgaga = (ghtt**2)*0.70904D-01 + (ghbb**2)*0.18760D-04 + (ghWW**2)*1.5863 + &
                      ghtt*ghbb*(-0.17319D-02) + ghtt*ghWW*(-0.67074) + &
                      ghbb*ghWW*0.82093D-02 + (ghtautau**2)*0.22663E-04 + &
                      ghtt*ghtautau*(-0.18696E-02) + ghbb*ghtautau*0.41239E-04 + &
                      ghtautau*ghWW*0.88634E-02

    end function get_g2hgaga
!**************************************************************
end program HSgga
