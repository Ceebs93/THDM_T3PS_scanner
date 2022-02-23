!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 26/09/2013).
!--------------------------------------------------------------------------------------
program HS_mass
! In this example the peak-centered chi^2 method is applied to a SM-like Higgs boson
! (overall signal strength scale factor mu) within the mass range 110 - 140 GeV.
! All three mass pdf choices are considered. Theoretical mass uncertainties and
! assignment range can be changed.
!--------------------------------------------------------------------------------------
    use theory_colliderSfunctions
    use usefulbits, only: vsmall
    use pc_chisq, only: print_cov_mh_to_file, print_cov_mu_to_file, print_inverse_cov_mh_to_file, &
                        get_peakchisq, print_corr_mu_to_file
    use io, only: get_number_of_observables, get_ID_of_peakobservable, get_peakinfo_from_HSresults
    implicit none

    integer :: nHzero, nHplus, i, j, k, ii, jj
    double precision :: obsratio, mass, Lambda
    double precision :: SMGammaTotal
    double precision :: Mh, GammaTotal, ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
        ghjmumu_s, ghjmumu_p, ghjtautau_s, ghjtautau_p, &
        ghjWW, ghjZZ, ghjZga, ghjgaga, ghjgg, &
        ghjhiZ

    character(len=100)::filename
    double precision :: dm
    integer                  :: pdf
    integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID, domH, nHcomb, Nassigned

    double precision :: Chisq_peak_mu, Chisq_peak_mh, Chisq_peak, Pvalue_peak
    double precision :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, Pvalue_STXS
    double precision :: Chisq_LHCRun1_mu, Chisq_LHCRun1_mh, Chisq_LHCRun1, Pvalue_LHCRun1
    integer :: nobs_STXS, nobs_peak, nobs_LHCRun1

    double precision :: mupred
    double precision, allocatable :: csqmu(:), csqmh(:), csqmax(:), csqtot(:)
    integer, allocatable :: ncomb(:)
!-HiggsBounds internal functions to obtain SM branching ratios
    double precision :: SMBR_Htoptop, SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau, &
        SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg, SMGamma_h
    logical :: use_assignmentrange

    nHzero = 1
    nHplus = 0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
    call initialize_HiggsSignals(nHzero, nHplus, "latestresults")
! call initialize_HiggsSignals(nHzero,nHplus,"latestresults-1.4.0-LHCinclusive")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
    call setup_output_level(0)

!--Enter the Higgs mass and its theory uncertainty here:
    dm = 0.0D0
    Lambda = 2.0D0
    use_assignmentrange = .False.

!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
    call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set the assignment range for various observables                                 ----!
! Could have different values, too...
    if (use_assignmentrange) then
        call setup_assignmentrange(Lambda)
        call setup_assignmentrange_massobservables(Lambda)
    else
        Lambda = 0.0D0
    endif

! Normalize rate prediction w.r.t. predicted Higgs mass (for any mass):
    ! call setup_rate_normalization(.False., .False.)
! Instead, normalize rate prediction w.r.t. observed Higgs mass (for any predicted mass):
!   call setup_rate_normalization(.True.,.True.)

    do pdf = 1, 3
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
        call setup_pdf(pdf)
        select case (pdf)
        case (1)
            filename = 'results/HS_mass_pdf1.dat'
        case (2)
            filename = 'results/HS_mass_pdf2.dat'
        case (3)
            filename = 'results/HS_mass_pdf3.dat'
        case default
        end select

        open (21, file=filename)
        write (21, *) '# mh   dmh     Chisq_mu    Chisq_mh    Chisq    Nassigned    ndf   Lambda'
        write (21, *) '#----------------------------------------------------'

        do j = 1, 101 !181,181!
!    mh = 120.0D0 +(j-1)*0.1D0
            mh = 122.5D0 + (j - 1)*0.05D0
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
                ghjZga = 1.0d0
                ghjgg = 1.0d0
                ghjhiZ = 0d0
                ghjgaga = 1.0d0

                call HiggsBounds_neutral_input_properties(Mh, SMGammaTotal)

                call HiggsBounds_neutral_input_effC( &
                    ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                    ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                    ghjmumu_s, ghjmumu_p, &
                    ghjtautau_s, ghjtautau_p, &
                    ghjWW, ghjZZ, ghjZga, &
                    ghjgaga, ghjgg, ghjhiZ)

                call run_HiggsSignals(Chisq_peak_mu, Chisq_peak_mh, Chisq_peak, nobs_peak, &
                                      Pvalue_peak)

                call run_HiggsSignals_LHC_Run1_combination(Chisq_LHCRun1_mu, Chisq_LHCRun1_mh, &
                                                           Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1)

                call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, &
                                           nobs_STXS, Pvalue_STXS)

                call complete_HS_results()

! PEAK OBSERVABLES ONLY!
                call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

                allocate (csqmu(npeakmu), csqmh(npeakmu), csqmax(npeakmu), csqtot(npeakmu), ncomb(npeakmu))

                Nassigned = 0
                do ii = 1, npeakmu
                    call get_ID_of_peakobservable(ii, ID)
                    call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
                    ncomb(ii) = nHcomb
                    call get_peakchisq(ID, csqmu(ii), csqmh(ii), csqmax(ii), csqtot(ii))
                    Nassigned = Nassigned + nHcomb
                enddo

                deallocate (csqmu, csqmh, csqmax, csqtot, ncomb)

                write (21, *) mh, dm, Chisq_peak_mu + Chisq_STXS_rates + Chisq_LHCRun1_mu, &
                    Chisq_peak_mh + Chisq_STXS_mh + Chisq_LHCRun1_mh, &
                    Chisq_peak + Chisq_STXS + Chisq_LHCRun1, Nassigned, &
                    nobs_peak + nobs_STXS + nobs_LHCRun1, Lambda

            endif
        enddo
        close (21)
    enddo

    write (*, *) "Finishing HiggsSignals..."
    call finish_HiggsSignals

end program HS_mass
