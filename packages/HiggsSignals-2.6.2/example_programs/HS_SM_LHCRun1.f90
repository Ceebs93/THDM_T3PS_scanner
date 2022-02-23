!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 26/09/2013).
!--------------------------------------------------------------------------------------
program HS_SM_LHCRun1
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

    integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
    double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh, Lambda
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
    double precision :: mupred
    double precision, allocatable :: csqmu(:), csqmh(:), csqmax(:), csqtot(:)
    integer, allocatable :: ncomb(:)
!-HiggsBounds internal functions to obtain SM branching ratios
    double precision :: SMBR_Htoptop, SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau, &
        SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg, SMGamma_h

    nHzero = 1
    nHplus = 0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
     call initialize_HiggsSignals(nHzero,nHplus,"latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
    call setup_output_level(0)

!--Enter the Higgs mass and its theory uncertainty here:
    dm = 1.0D0
    Lambda = 10.0D0

!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
    call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set the assignment range for the peak-centered method (optional)                                 ----!
!     This can be done either to all observables or only to the
!     mass-sensitive observables, which contribute to the Higgs mass chi^2
! call setup_assignmentrange(Lambda)

!  call setup_correlations(.true.,.true.) ! mu, mass

    do pdf = 1, 3
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
        call setup_pdf(pdf)
        select case (pdf)
        case (1)
            filename = 'results/HS_SM_LHCrun1_mass_pdf1.dat'
        case (2)
            filename = 'results/HS_SM_LHCrun1_mass_pdf2.dat'
        case (3)
            filename = 'results/HS_SM_LHCrun1_mass_pdf3.dat'
        case default
        end select

        open (21, file=filename)
        write (21, *) '# mh   dmh     Chisq_mu    Chisq_mh    Chisq    Nassigned    ndf   Lambda'
        write (21, *) '#----------------------------------------------------'

        do j = 1, 101 !181,181!
            mh = 120.D0 + (j - 1)*0.1D0

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

                call run_HiggsSignals_LHC_Run1_combination(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!     call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!     call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

!     allocate(csqmu(npeakmu),csqmh(npeakmu),csqmax(npeakmu),csqtot(npeakmu),ncomb(npeakmu))

                !   Nassigned=0
!     do ii=1,npeakmu
!      call get_ID_of_peakobservable(ii, ID)
!      call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
!      ncomb(ii)=nHcomb
!      call get_peakchisq(ID, csqmu(ii), csqmh(ii), csqmax(ii), csqtot(ii))
!      Nassigned=Nassigned+nHcomb
!     enddo
!
!     deallocate(csqmu,csqmh,csqmax,csqtot,ncomb)

                write (21, *) mh, dm, Chisq_mu, Chisq_mh, Chisq, ndf, Lambda

            endif
        enddo
        close (21)
    enddo

    write (*, *) "Finishing HiggsSignals..."
    call finish_HiggsSignals

end program HS_SM_LHCRun1
