!--------------------------------------------------------------------------------------
program HSwithSLHA
!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 29/01/2013).
!
! In this example we demonstrate how HiggsSignals can be run on SLHA files. The SLHA
! file has to contain the two HiggsBounds SLHA input blocks
!        HiggsBoundsInputHiggsCouplingsBosons
!        HiggsBoundsInputHiggsCouplingsFermions
! (see HiggsBounds (version 3 or more) manual for more details)
!
! Run with
!   ./HSwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output:
! The HiggsSignals SLHA output blocks will be added to each SLHA file.
! The results are summarized in an additional file <stem>-fromHS
!
! We furthermore demonstrate how to get the signal rates directly from HiggsSignals
! after a successful run.
!--------------------------------------------------------------------------------------
    use io, only: HiggsSignals_SLHA_output, get_peakinfo_from_HSresults
    use pc_chisq, only: print_cov_mu_to_file, print_peaks_to_file
    implicit none
    integer :: nH, nHplus, nobs_peak, nobs_STXS, nobs_LHCRun1
    double precision :: Pvalue_peak, Chisq_peak, Chisq_peak_mu, Chisq_peak_mh
    double precision :: Pvalue_STXS, Chisq_STXS, Chisq_STXS_rates, Chisq_STXS_mh
    double precision :: Pvalue_LHCRun1, Chisq_LHCRun1, Chisq_LHCRun1_mu, Chisq_LHCRun1_mh
    double precision :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb, &
        totalrate, ggf_rate, hgaga_rate, gghgg_rate
    double precision, allocatable :: dMh(:), masses(:), dmasses(:)
    integer :: i, npoints
    integer, parameter :: fileid = 78, fileid2 = 79
    character(len=8) :: istring
    character(len=300) :: inputfilename, outputfilename
    character(len=300) :: stem
    character(LEN=300) :: temp, tmpstring
    integer :: number_args, ios
!--------------------------------------------------------------------------------------
    nH = 3
    nHplus = 1

    allocate (dMh(nH), masses(nH), dmasses(nH))
! These are default values. In SLHA mode, these values are read in from the Block "DMASS".
    dMh = (/2.0D0, 2.0D0, 0.0D0/)
!----------------------------------- preprocessing --------------------------------------!
    number_args = IARGC()
    if (number_args .ne. 2) then
        stop "Incorrect number of arguments given to HSwithSLHA"
    endif
    ! Read arguments into text strings.
    i = 1
    temp = ""
    call GETARG(i, temp)
    read (temp, *) npoints
    i = i + 1
    temp = ""
    call GETARG(i, temp)
    stem = ""
    stem = trim(temp)
!------------------------------ HiggsSignals options ------------------------------------!
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
    call initialize_HiggsSignals(nH, nHplus, "latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)          ----!
    call setup_output_level(1)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)                  ----!
    call setup_pdf(2)
    !---- Set number of free model parameters (used for p-value calculation only)        ----!
    call setup_nparam(0)
!---- Set the assignment range for the peak-centered method (optional)                                 ----!
!  call setup_assignmentrange_massobservables(2.0D0)
!---- Use symmetric rate errors? (0: original(default), 1: averaged-symmetrical)     ----!
! call setup_symmetricerrors(.false.)
!----

    outputfilename = trim(adjustl(stem))//'-fromHS'
    open (fileid, file=trim(outputfilename))

    do i = 1, npoints
        if (i .gt. 99999999) stop 'need to increase the size of istring in HSwithSLHA'
        write (istring, '(I8)') i
        inputfilename = trim(adjustl(stem))//'.'//trim(adjustl(istring))
!--Test if input file exists and is non-empty
        open (fileid2, file=inputfilename, form='formatted')
        read (fileid2, '(A)', iostat=ios) tmpstring
        print *, ios
        if (ios .eq. 0) then
            close (fileid2)

!-------------------------------- HiggsSignals run --------------------------------------!
!---- Feed HiggsSignals with the the model input using HiggsBounds subroutine                  ----!
            call HiggsBounds_input_SLHA(inputfilename)
!---- Checking the Higgs mass uncertainty                                                                             ----!
! The theoretical Higgs mass uncertainties are read in from the SLHA Block "DMASS"
! in the call of the subroutine HiggsBounds_input_SLHA. If the block "DMASS" is absent,
! they are set to zero. If the user wants to change the values obtained from the SLHA
! file, he/she can call the subroutine
!         HiggsSignals_neutral_input_MassUncertainty
! AFTER reading in the SLHA file.
! Here, we only print out the values which have been stored already in
! HiggsBounds/HiggsSignals. If their sum is <= zero, the Block DMASS was probably absent
! and we set the default values as specified above.
!
            call get_neutral_Higgs_masses(masses, dmasses)
            write (*, *) "Neutral Higgs boson mass spectrum (from SLHA): "
            write (*, '(2X,A,F8.2,A,F4.2)') "mass(h0) = ", masses(1), " +- ", dmasses(1)
            write (*, '(2X,A,F8.2,A,F4.2)') "mass(H0) = ", masses(2), " +- ", dmasses(2)
            write (*, '(2X,A,F8.2,A,F4.2)') "mass(A0) = ", masses(3), " +- ", dmasses(3)
            if (sum(dmasses) .le. 0.0D0) then
                write (*, *) "BLOCK DMASS not found, changing mass uncertainies to ", dMh
                call HiggsSignals_neutral_input_MassUncertainty(dMh)
            endif
!---- Run HiggsSignals on peak observables (13 TeV)
            call run_HiggsSignals(Chisq_peak_mu, Chisq_peak_mh, Chisq_peak, nobs_peak, Pvalue_peak)
!---- Run HiggsSignals on STXS observables (13 TeV)
            call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)
!---- Run HiggsSignals on LHC Run 1 observables        (7/8 TeV)
            call run_HiggsSignals_LHC_Run1_combination(Chisq_LHCRun1_mu, Chisq_LHCRun1_mh, &
                                                       Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1)
            call complete_HS_results()
!----------------------------- HiggsSignals output --------------------------------------!
!---- Attach HiggsSignals SLHA blocks to SLHA file                                                                         ----!
!   integer argument gives level of details:
!   0 : writes only HiggsSignalsResults block
!   else : writes all blocks
            call HiggsSignals_SLHA_output(0)
!---- Now, some examples of how to read out the signal rates. Note, that these subroutines
!     have to be called after run_HiggsSignals
!
!---- Get signal-rate ratios (without efficiencies) for lightest Higgs boson and LHC13----!
            call get_Rvalues(1, 4, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
!---- Get the total signal rate (without efficiencies)                                                                  ----!
            call get_rates(1, 4, 25, (/"1.1", "1.2", "1.3", "1.4", "1.5", &
                                       "2.1", "2.2", "2.3", "2.4", "2.5", &
                                       "3.1", "3.2", "3.3", "3.4", "3.5", &
                                       "4.1", "4.2", "4.3", "4.4", "4.5", &
                                       "5.1", "5.2", "5.3", "5.4", &
                                       "5.5"/), totalrate)
!---- Get the gluon gluon fusion rate (without efficiencies)                                                  ----!
            call get_rates(1, 4, 1, (/"1.0"/), ggf_rate)
!---- Get the H -> gamma gamma rate (without efficiencies)                                                          ----!
            call get_rates(1, 4, 1, (/"0.1"/), hgaga_rate)
!        NEW SINCE HiggsSignals-1.1.0: more decay modes accessible via get_rates:
!          Decay mode ID (Final state):  6 (Zgamma), 7 (cc), 8 (mumu), 9 (gg)
!---- Get the gg->H->gg (without efficiencies)                                                                      ----!
            call get_rates(1, 4, 1, (/"1.9"/), gghgg_rate)

            write (*, '(A,F10.4)') "R_H_WW     = ", R_H_WW
            write (*, '(A,F10.4)') "R_H_ZZ     = ", R_H_ZZ
            write (*, '(A,F10.4)') "R_H_gaga   = ", R_H_gaga
            write (*, '(A,F10.4)') "R_H_tautau = ", R_H_tautau
            write (*, '(A,F10.4)') "R_H_bb     = ", R_H_bb
            write (*, '(A,F10.4)') "R_VH_bb    = ", R_VH_bb
            write (*, '(A,F10.4)') "totalrate  = ", totalrate
            write (*, '(A,F10.4)') "ggf_rate   = ", ggf_rate
            write (*, '(A,F10.4)') "h->gaga    = ", hgaga_rate
            write (*, '(A,F10.4)') "gg->h->gg  = ", gghgg_rate

!--This will collect the main HiggsSignals results together into one file
            write (fileid, *) i, Chisq_peak_mu + Chisq_STXS_rates + Chisq_LHCRun1_mu, &
                Chisq_peak_mh + Chisq_STXS_mh + Chisq_LHCRun1_mh, &
                nobs_peak + nobs_STXS + nobs_LHCRun1
        else
            close (fileid2)
            call system("rm -f "//inputfilename)
        endif

    enddo

    close (fileid)
! call print_peaks_to_file
! call print_cov_mu_to_file

    call finish_HiggsSignals

end program HSwithSLHA
