!******************************************************
! This example program is part of HiggsSignals (TS 05/03/2013).
!******************************************************
program HBandHSwithSLHA
!
! In this example we run both HiggsBounds and HiggsSignals simultaneously
! on SLHA file(s). NOTE: The feature of selecting different experimental
! data in HiggsBounds is not fully supported here. In particular, the onlyL
! option of HiggsBounds does not work (it will be turned into LandH).
! If you want to take into account only LEP limits in HiggsBounds, you have
! to run both programs separately, using e.g. the example programs
! HBwithSLHA and HSwithSLHA.
!
! The SLHA file(s) has to contain the two HiggsBounds SLHA
! input blocks
!        HiggsBoundsInputHiggsCouplingsBosons
!        HiggsBoundsInputHiggsCouplingsFermions
! (see HiggsBounds (version 3 or more) manual for more details)
!
! Run with
!   ./HBandHSwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output
! The block HiggsSignalsResults will be added to each SLHA file.
!
!******************************************************
    use io, only: HiggsSignals_SLHA_output, get_peakinfo_from_HSresults

    implicit none
    integer :: nH, nHplus, nobs_peak, nobs_STXS, nobs_LHCRun1, HBresult, chan, ncombined
    double precision :: Pvalue_peak, Chisq_peak, Chisq_peak_mu, Chisq_peak_mh
    double precision :: Pvalue_STXS, Chisq_STXS, Chisq_STXS_rates, Chisq_STXS_mh
    double precision :: Pvalue_LHCRun1, Chisq_LHCRun1, Chisq_LHCRun1_mu, Chisq_LHCRun1_mh
    double precision :: obsratio
    double precision, allocatable :: dMh(:), dCS(:), dBR(:)
    integer :: i, npoints
    integer, parameter :: fileid = 78, fileid2 = 79
    character(len=8) :: istring
    character(len=300) :: inputfilename, outputfilename
    character(len=300) :: stem
    character(LEN=300) :: temp, tmpstring
    integer :: number_args, ios

    nH = 3
    nHplus = 1

    allocate (dMh(nH))
!--n.b. have to set theoretical uncertainties on Higgs masses dMh (in GeV) for h,H,A:
    dMh = (/2.0D0, 2.0D0, 0.0D0/)
!-------------------------- preprocessing ------------------------------!
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
!---------------------------- HiggsBounds and HiggsSignals ------------------------------!
!---- Initialize HiggsBounds and specify the dataset it should use                                      ----!
    call initialize_HiggsBounds(nH, nHplus, 'LandH')
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
    call initialize_HiggsSignals(nH, nHplus, "latestresults")
!------------------------------ HiggsSignals options ------------------------------------!
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)          ----!
    call setup_output_level(0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)                  ----!
    call setup_pdf(2)

    outputfilename = trim(adjustl(stem))//'-fromHBandHS'

    open (fileid, file=trim(outputfilename))

    do i = 1, npoints

        if (i .gt. 99999999) stop 'need to increase the size of istring in HSwithSLHA'
        write (istring, '(I8)') i

        inputfilename = trim(adjustl(stem))//'.'//trim(adjustl(istring))

!--Test if input file exists and is non-empty
        open (fileid2, file=inputfilename, form='formatted')
        read (fileid2, '(A)', iostat=ios) tmpstring

        if (ios .eq. 0) then
            close (fileid2)

!---------------------- HiggsBounds and HiggsSignals run --------------------------------!
!---- Feed HiggsBounds/Signals with the the model input using HiggsBounds subroutine ----!
            call HiggsBounds_input_SLHA(inputfilename)
!---- We want to use the mass variation treatment for the theoretical uncertainty    ----!
!          HiggsBounds. Thus we set the mass uncertainties here (neutral Higgses mass errors
!     set to dMh, charged Higgs mass error set to 0.0D0)
            call HiggsBounds_set_mass_uncertainties(dMh, 0.0D0)
!---- First, run HiggsBounds                                                         ----!
            call run_HiggsBounds(HBresult, chan, obsratio, ncombined)
!---- Now, we have to fill again the input for the HiggsSignals run                  ----!
            call HiggsBounds_input_SLHA(inputfilename)
!---- Pass the Higgs mass uncertainty to HiggsSignals                                ----!
            call HiggsSignals_neutral_input_MassUncertainty(dMh)
!---- Run HiggsSignals on peak observables (13 TeV)                                  ----!
            call run_HiggsSignals(Chisq_peak_mu, Chisq_peak_mh, Chisq_peak, nobs_peak, Pvalue_peak)
!---- Run HiggsSignals on STXS observables (13 TeV)                                  ----!
            call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)
!---- Run HiggsSignals on LHC Run 1 observables        (7/8 TeV)                            ----!
            call run_HiggsSignals_LHC_Run1_combination(Chisq_LHCRun1_mu, Chisq_LHCRun1_mh, &
                                                       Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1)
            call complete_HS_results()
!----------------------------- HiggsSignals output --------------------------------------!
!---- Attach HiggsBounds SLHA output block to SLHA file                                                           ----!
            call HiggsBounds_SLHA_output
!---- Attach HiggsSignals SLHA output blocks to SLHA file                                                         ----!
!   integer argument gives level of details:
!   0 : writes only HiggsSignalsResults block
!   else : writes all blocks
            call HiggsSignals_SLHA_output(0)
!---- This will collect the main HiggsSignals results together into one file                 ----!
            write (fileid, *) i, Chisq_peak_mu + Chisq_STXS_rates + Chisq_LHCRun1_mu, &
                Chisq_peak_mh + Chisq_STXS_mh + Chisq_LHCRun1_mh, &
                nobs_peak + nobs_STXS + nobs_LHCRun1, HBresult, chan, obsratio, ncombined
        else
            close (fileid2)
            call system("rm -f "//inputfilename)
        endif
    enddo

    close (fileid)

    call finish_HiggsBounds
    call finish_HiggsSignals

end program HBandHSwithSLHA
