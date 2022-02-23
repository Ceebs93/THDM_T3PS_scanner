!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 30/10/2014).
!
! In this example we scan over the masses of two Higgs bosons, for fixed values of
! universal signal strength scale factors (scale1, scale2) and theoretical mass un-
! certainties (dm(1), dm(2)). We loop over the three different pdf choices.
! The output (total chi^2) can be plotted by using the gnuplot script plot_2Higgses.gnu.
!--------------------------------------------------------------------------------------
program HS_2Higgses
!--------------------------------------------------------------------------------------
    use theory_colliderSfunctions
    use usefulbits, only: vsmall
    use usefulbits_hs, only: HSres
    use pc_chisq, only: print_inverse_cov_mh_to_file, get_masschisq_from_separation
    implicit none

    integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
    double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
    double precision :: Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, Chisq_LHCRun1, Pvalue_LHCRun1
    double precision :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, Pvalue_STXS
    integer ::  nobs_LHCRun1, nobs_STXS
    double precision :: SMGammaTotal(2)
    double precision :: scale_bbh, scale_ggh, dggh, dbbh
    double precision :: Mh(2), GammaTotal(2), ghjss_s(2), ghjss_p(2), ghjcc_s(2), ghjcc_p(2), &
        ghjbb_s(2), ghjbb_p(2), ghjtt_s(2), ghjtt_p(2), &
        ghjmumu_s(2), ghjmumu_p(2), ghjtautau_s(2), ghjtautau_p(2), &
        ghjWW(2), ghjZZ(2), ghjZga(2), ghjgaga(2), ghjgg(2), &
        ghjhiZ(2, 2)
    character(len=100)::filename
    character(len=1)::pdfchar
    double precision :: dm(2)
    integer :: pdf
    double precision :: scale1, scale2, csq_sep
    nHzero = 2
    nHplus = 0

!--Enter the Higgs mass and its theory uncertainty here:
    Mh = (/125.0D0, 125.0D0/)
    dm = (/0.5D0, 0.5D0/)

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
    call initialize_HiggsSignals(nHzero, nHplus, "latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
    call setup_output_level(0)
! call setup_assignmentrange_massobservables(2.0D0)
! Always normalize rate prediction w.r.t. to predicted Higgs mass
!  call setup_rate_normalization(.False.,.False.)
!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
    call HiggsSignals_neutral_input_MassUncertainty(dm)
! Set number of free parameters (for p-value evaluation)
    call setup_nparam(2)

!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
    call setup_pdf(2)
    open (21, file='results/2Higgses_pdf2.dat')
    write (21, *) ' Mh1 Mh2 dm1 dm2 scale1 scale2 chisq_mu chisq_mh chisq ndf Pvalue chisq_sep'

    scale1 = sqrt(5D-1)
    scale2 = sqrt(1-scale1**2)

    ghjss_s(1) = scale1
    ghjss_s(2) = scale2
    ghjss_p = 0.0d0
    ghjcc_s(1) = scale1
    ghjcc_s(2) = scale2
    ghjcc_p = 0.0d0
    ghjbb_s(1) = scale1
    ghjbb_s(2) = scale2
    ghjbb_p = 0.0d0
    ghjtt_s(1) = scale1
    ghjtt_s(2) = scale2
    ghjtt_p = 0.0d0
    ghjmumu_s(1) = scale1
    ghjmumu_s(2) = scale2
    ghjmumu_p = 0.0d0
    ghjtautau_s(1) = scale1
    ghjtautau_s(2) = scale2
    ghjtautau_p = 0.0d0
    ghjWW(1) = scale1
    ghjWW(2) = scale2
    ghjZZ(1) = scale1
    ghjZZ(2) = scale2
    ghjZga(1) = scale1
    ghjZga(2) = scale2
    ghjgg(1) = scale1
    ghjgg(2) = scale2
!     ghjggZ(1) = scale1
!     ghjggZ(2) = scale2
    ghjhiZ = 0d0
    ghjgaga(1) = scale1
    ghjgaga(2) = scale2

    GammaTotal = -1D0

    do i = 1, 17
        do j = 1, 17
            Mh(1) = 123D0 + (i - 1)*0.25D0
            Mh(2) = 123D0 + (j - 1)*0.25D0

            call HiggsBounds_neutral_input_properties(Mh, GammaTotal, (/1, 1/))

            call HiggsBounds_neutral_input_effC( &
                ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                ghjmumu_s, ghjmumu_p, &
                ghjtautau_s, ghjtautau_p, &
                ghjWW, ghjZZ, ghjZga, &
                ghjgaga, ghjgg, ghjhiZ)

            call run_HiggsSignals_full(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

            call get_masschisq_from_separation(csq_sep) ! gets only the mass separation-chi^2 from peak obs.
            csq_sep = csq_sep + HSres(1)%Chisq_LHCRun1_mhsep ! adds the LHC Run-1 mass separation-chi^2

            write (21, *) Mh, dm, scale1, scale2, Chisq_mu, Chisq_mh, Chisq, &
                ndf - 2, Pvalue, csq_sep
        end do
    end do
    close (21)

    write (*, *) "Finishing HiggsSignals..."
    call finish_HiggsSignals

end program HS_2Higgses
