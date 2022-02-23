!--------------------------------------------------------------------------------------
program HS_efficiencies
!
! This example program is part of HiggsSignals (TS 27/08/2013).
!--------------------------------------------------------------------------------------
! This example program demonstrates how models with different channel efficiencies
! of the Higgs searches can be treated in HiggsSignals (HS).
!--------------------------------------------------------------------------------------
! use theory_colliderSfunctions
! use usefulbits, only : vsmall
! use usefulbits_hs,only : analyses
! use pc_chisq, only : print_peaks_to_LaTeX
    use io, only: get_peakinfo_from_HSresults, get_number_of_observables, &
                  get_ID_of_peakobservable, get_peak_channels

    implicit none
    integer :: nH, nHplus, ndf, jj, kk, ll, mm
    double precision :: Chisq, Chisq_mu, Chisq_mh, Pvalue
!  double precision, allocatable :: dMh(:)
!  double precision :: dCS(5),dBR(5),dggh, dbbh
    double precision :: SMGamma_h
    double precision :: Mh, GammaTotal, ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
        ghjmumu_s, ghjmumu_p, ghjtautau_s, ghjtautau_p, &
        ghjWW, ghjZZ, ghjZga, ghjgaga, ghjgg, &
        ghjhiZ
    double precision :: mupred, etaZ, sf, mu_average
    integer :: domH, nHcomb
    integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID, pID
    integer :: i, npoints
    character(len=8) :: istring
    character(len=300) :: inputfilename, outputfilename
    character(len=300) :: stem
    character(LEN=300) :: temp
    integer :: number_args, stat
    integer :: Nc
    integer, allocatable :: channel_p_IDs(:), channel_d_IDs(:)
    double precision, allocatable :: efficiencies(:), eff_ratio(:)

!---Note here: We only run HiggsSignals on the lightest Higgs boson. This can be easily
!---extended to all 3 MSSM neutral Higgs bosons. In that case, the effective couplings
!---and mass uncertainties have to be given as arrays of size=nH (Cf. the Higgsbounds
!---manual for HB-3.x.x for how to call HiggsBounds_neutral_input_effC correctly!)
    nH = 1
    nHplus = 0

!   allocate(dMh(nH))
!--n.b. have to set theoretical uncertainties on Higgs mass dMh (in GeV):
!   dMh = (/ 0.0D0 /)
!-------------------------- HiggsSignals ------------------------------!

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
    call initialize_HiggsSignals(nH, nHplus, "latestresults")
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)                  ----!
    call setup_pdf(2)
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)          ----!
    call setup_output_level(0)
!---- Pass the Higgs mass uncertainty to HiggsSignals                                                                  ----!
!   call HiggsSignals_neutral_input_MassUncertainty(dMh)
!---- Use symmetric rate errors? (0: original(default), 1: averaged-symmetrical)     ----!
!  call setup_symmetricerrors(.false.)
!---- Setup a wider assignment range                                                 ----!
!  call setup_assignmentrange(10.0D0)

!----HiggsBounds/Signals effective couplings input.
!           These have to be inserted for the model which we want to test, i.e. we would have
!    to write an interface to set via arguments in the executables call, or reading
!    in a text file, etc.
!----For now, we set them by hand to the SM values except for the vector boson couplings,
!    where we set them to a scale factor sf (for demonstration):

    do mm = 1, 3
        select case (mm)
        case (1)
            outputfilename = "results/HS_efficiencies_kv0p9.dat"
        case (2)
            outputfilename = "results/HS_efficiencies_kv1p0.dat"
        case (3)
            outputfilename = "results/HS_efficiencies_kv1p1.dat"
!    case(4)
!     outputfilename="results/HS_efficiencies_kv2p0.dat"
        case default
        end select
        open (21, file=outputfilename)

        sf = 0.9D0 + (mm - 1)*0.1D0
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
        ghjWW = sf
        ghjZZ = sf
        ghjZga = 1.0d0
        ghjgg = 1.0d0
        ghjhiZ = 0d0
        ghjgaga = 1.0d0

        Mh = 125.09D0
        GammaTotal = SMGamma_h(Mh)

!-Set the HiggsSignals input
        call HiggsBounds_neutral_input_properties(Mh, GammaTotal)

        call HiggsBounds_neutral_input_effC( &
            ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
            ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
            ghjmumu_s, ghjmumu_p, &
            ghjtautau_s, ghjtautau_p, &
            ghjWW, ghjZZ, ghjZga, &
            ghjgaga, ghjgg, ghjhiZ)

!-Run HS on the original experimental data in order to evaluate the model predictions
        call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!-Get the number of the peak-observables (Don't care about ntotal, npeakmh, nmpred, nanalyses)
        call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

! We now want to set different efficiencies as evaluated e.g. from a MC analysis.
! Note that these are given as relative changes with respect to the SM efficiencies, i.e.
! in the MC you should both test your model and the SM with the analysis cuts and
! estimate the ratio eta = efficiency(model)/efficiency(SM) for each Higgs channel.
!
! As an example, we want to give all VBF,WH,ZH channels a factor of etaZ in all implemented
! searches where these channels are present.

        do ll = 1, 51
            etaZ = 0.0D0 + (ll - 1)*0.05D0
!-Loop over the number of peak observables
            do kk = 1, npeakmu
!--Get for each peak observable its unique ID:
                call get_ID_of_peakobservable(kk, ID)
!--Get the predicted signal strength modifier (mupred) for this peak observable:
!   call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
!--Get channel information
!    (Nc),dIDs(Nc),efficiencies(Nc))
                call get_peak_channels(ID, Nc, channel_p_IDs, channel_d_IDs, efficiencies)
!   write(*,*) ID, Nc
!   write(*,*) channel_p_IDs
!   write(*,*) channel_d_IDs
!   write(*,*) efficiencies
                allocate (eff_ratio(Nc))
                eff_ratio = 1.0D0
                do jj = lbound(channel_p_IDs, dim=1), ubound(channel_p_IDs, dim=1)
!---Get the ID for the production mode (2:VBF, 3:WH, 4:ZH):
                    pID = channel_p_IDs(jj)
                    if (pID .eq. 2 .or. pID .eq. 3 .or. pID .eq. 4) then
                        eff_ratio(jj) = etaZ
                    endif
                enddo
!---Hand over the efficiency ratios:
                call assign_modelefficiencies_to_peak(ID, Nc, eff_ratio)
                deallocate (channel_p_IDs, channel_d_IDs, efficiencies, eff_ratio)
            enddo

!-Set the HiggsSignals input (again!)
            call HiggsBounds_neutral_input_properties(Mh, GammaTotal)

            call HiggsBounds_neutral_input_effC( &
                ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                ghjmumu_s, ghjmumu_p, &
                ghjtautau_s, ghjtautau_p, &
                ghjWW, ghjZZ, ghjZga, &
                ghjgaga, ghjgg, ghjhiZ)

            call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

            mu_average = 0.0D0
            do kk = 1, npeakmu
!--Get for each peak observable its unique ID:
                call get_ID_of_peakobservable(kk, ID)
!--Get the predicted signal strength modifier (mupred) for this peak observable:
                call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
                mu_average = mu_average + mupred
            enddo
            mu_average = mu_average/npeakmu

            write (21, *) sf, etaZ, Chisq_mu, mu_average

        enddo
        close (21)
    enddo

    call finish_HiggsSignals

!--------------------------------------------------------------------------------------
end program HS_efficiencies
!--------------------------------------------------------------------------------------
