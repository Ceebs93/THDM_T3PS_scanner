!> @file
!! @brief HiggsSignals user subroutines
!!
!! These subroutines are the recommended way of interfacing your code with HiggsSignals.
!! A typical use looks as follows:
!!   1. At the beginning, call an initialization routine (usually initialize_higgssignals_latestresults()),
!!      only do this once in your code. Don't forget to initialize HiggsBounds as well.
!!   2. For every parameter point, use the HiggsBounds input routines to provide the model data.
!!      Then run the analysis using run_higgssignals(), run_higgssignals_stxs(),
!!      and/or run_higgssignals_lhc_run1_combination(). We recommend to call all three of them as
!!      they are complementary. You can then obtain the total \f$\chi^2\f$ by summing the three contributions
!!      `chisq = chisq_peak + chisq_stxs + chisq_run1` and use this in your fit.
!!   3. When you are done, call finish_higgssignals() and finish_higgsbounds().
!!
!! Many of the subroutines in this file (mainly those named `setup_...`) modify the behaviour
!! of these core routines. We strive to provide a good default behaviour, so you should only need to use
!! these if you want to study very specific effects.
!!
!! Several subroutines (those named `get_...`) provide additional information that may be of interest to the user.
!! In particular get_rvalues() gives access to the values of the most important signal rates used in HiggsSignals.
!! Saving these can be useful to get a better understanding of how your model deviates from the SM expectation.

!> Initialize HiggsSignals
!! The first HiggsSignals subroutine that should be called
!! by the user. It calls subroutines to read in the tables of Standard
!! Model decay and production rates from HiggsBounds, sets up the
!! experimental data from Tevatron and LHC, allocate arrays, etc.
!! @param nHiggsneut Number of neutral Higgs bosons in the model
!! @param nHiggsplus Number of charged Higgs bosons in the model
!! @param Expt_string name of experimental dataset to be used
subroutine initialize_HiggsSignals(nHiggsneut, nHiggsplus, Expt_string)
    use usefulbits, only: np, Hneut, Hplus, Chineut, Chiplus, debug, inputmethod, &
                          theo, whichanalyses, just_after_run
    use usefulbits_HS, only: HiggsSignals_info, nanalys, eps, Exptdir, obs
    use datatables, only: setup_observables, setup_LHC_Run1_combination
    use STXS, only: load_STXS
    use input, only: check_number_of_particles, check_whichanalyses
    use io, only: setup_input_for_hs, setup_output_for_hs
    use theory_BRfunctions, only: setup_BRSM, BRSM
    use theory_XS_SM_functions, only: setup_XSSM, XSSM
    use install_data_hs, only: pathname_HS

    implicit none
    !--------------------------------------input
    integer, intent(in) :: nHiggsneut
    integer, intent(in) :: nHiggsplus
    character(LEN=*), intent(in) :: Expt_string
    !-----------------------------------internal
    integer :: status
    logical :: exptdirpresent = .False.
    !----------------------------------parameter
    eps = 5.0D0
    np(Hneut) = nHiggsneut

    np(Hplus) = nHiggsplus

    if (Expt_string .ne. 'none') then
        Exptdir = Expt_string
        exptdirpresent = .True.
    endif

    call system("ls "//pathname_HS//"/Expt_tables/"//Expt_string//"&>/dev/null", status)
    if (status .ne. 0 .and. exptdirpresent) then
        write (*, *) "Error: Experimental dataset '"//trim(adjustl(Expt_string))// &
            "' does not exist in "//pathname_HS//"Expt_tables/"
        stop
    endif

    np(Chineut) = 0! not considering bounds on neutralinos here
    np(Chiplus) = 0! not considering bounds on charginos here

    debug = .False.

    select case (whichanalyses)
    case ('onlyL')
        whichanalyses = 'LandH'
    case ('onlyH', 'onlyP', 'list ', 'LandH')
    case default
        whichanalyses = 'onlyH'
    end select

    call HiggsSignals_info
    if (inputmethod == 'subrout') then
        if (allocated(theo)) then
            if (debug) write (*, *) "HiggsBounds/HiggsSignals internal structure already initialized!"
        else
            if (debug) write (*, *) 'doing other preliminary tasks...'; call flush (6)
            call setup_input_for_hs
        endif
    endif

    if (debug) write (*, *) 'reading in Standard Model tables...'; call flush (6)
    if (.not. allocated(BRSM)) call setup_BRSM
    if (.not. allocated(XSSM)) call setup_XSSM
    call setup_uncertainties

    if (debug) write (*, *) 'reading in experimental data...'; call flush (6)
    if (exptdirpresent) call setup_observables
    if (exptdirpresent) call load_STXS(Expt_string)

    call setup_LHC_Run1_combination

    if (debug) write (*, *) 'sorting out processes to be checked...'; call flush (6)
    nanalys = size(obs)

    if (debug) write (*, *) 'preparing output arrays...'; call flush (6)
    call setup_output_for_hs

    if (debug) write (*, *) 'HiggsSignals has been initialized...'; call flush (6)

    just_after_run = .False.

contains
    subroutine setup_uncertainties
        use install_data_hs, only: pathname_HS
        use io, only: read_matrix_from_file
        use usefulbits_hs, only: delta_rate

        logical :: BRmodel, BRSM, XSmodel, XSSM, XScorr, XScorrSM

        call read_matrix_from_file(9, pathname_HS//"BRcov.in", delta_rate%BRcov, BRmodel)
        call read_matrix_from_file(9, pathname_HS//"BRcovSM.in", delta_rate%BRcovSM, BRSM)
        call read_matrix_from_file(12, pathname_HS//"XScov.in", delta_rate%CScov, XSmodel)
        call read_matrix_from_file(12, pathname_HS//"XScovSM.in", delta_rate%CScovSM, XSSM)
        call read_matrix_from_file(12, pathname_HS//"XScov_13TeV.in", delta_rate%CS13cov, XSmodel)
        call read_matrix_from_file(12, pathname_HS//"XScovSM_13TeV.in", delta_rate%CS13covSM, XSSM)
        call read_matrix_from_file(12, pathname_HS//"XScorr.in", delta_rate%CScorr, XScorr)
        call read_matrix_from_file(12, pathname_HS//"XScorrSM.in", delta_rate%CScorrSM, XScorrSM)
        call read_matrix_from_file(12, pathname_HS//"XScorr_13TeV.in", delta_rate%CS13corr, XScorr)
        call read_matrix_from_file(12, pathname_HS//"XScorrSM_13TeV.in", delta_rate%CS13corrSM, XScorrSM)

        if (BRmodel .and. BRSM) then
            delta_rate%BRcov_ok = .True.
            write (*, *) "Covariance matrix for relative branching ratio uncertainties read in successfully."
        else
            write (*, *) "Covariance matrix for relative branching ratio uncertainties not provided. Using default values."
        endif
        if (XSmodel .and. XSSM) then
            delta_rate%CScov_ok = .True.
            write (*, *) "Covariance matrix for relative cross section uncertainties read in successfully."
        else
            write (*, *) "Covariance matrix for relative cross section uncertainties not provided. Using default values."
        endif
        if (XScorr .and. XScorrSM) then
            delta_rate%CScorr_ok = .True.
            write (*, *) "Correlation matrix for relative cross section uncertainties read in successfully."
        else
            write (*, *) "Correlation matrix for relative cross section uncertainties not provided. Using default values."
        endif
    end subroutine setup_uncertainties
end subroutine initialize_HiggsSignals

!> Initializes HiggsSignals using the latest results.
!! Wrapper subroutine to intitialize HiggsSignals with the experimental
!! dataset "latestresults", avoiding to specify this via a string argument.
!! @param nHiggsneut Number of neutral Higgs bosons in the model
!! @param nHiggsplus Number of charged Higgs bosons in the model
subroutine initialize_HiggsSignals_latestresults(nHiggsneut, nHiggsplus)
    implicit none
    integer, intent(in) :: nHiggsneut
    integer, intent(in) :: nHiggsplus
    call initialize_HiggsSignals(nHiggsneut, nHiggsplus, "latestresults")
end subroutine initialize_HiggsSignals_latestresults

!> Initializes HiggsSignals without a dataset.
!! @param nHiggsneut Number of neutral Higgs bosons in the model
!! @param nHiggsplus Number of charged Higgs bosons in the model
subroutine initialize_HiggsSignals_empty(nHiggsneut, nHiggsplus)
    integer, intent(in) :: nHiggsneut
    integer, intent(in) :: nHiggsplus
    call initialize_HiggsSignals(nHiggsneut, nHiggsplus, "none")
end subroutine initialize_HiggsSignals_empty

!> Sets the theoretical mass uncertainty of the Higgs bosons.
!! @param dMh theoretical mass uncertainties
subroutine HiggsSignals_neutral_input_MassUncertainty(dMh)
    use usefulbits, only: theo, np, Hneut

    implicit none
    double precision, intent(in) :: dMh(np(Hneut))

    if (.not. allocated(theo)) then
        stop 'subroutine HiggsSignals_initialize must be called first'
    endif
    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsSignal_neutral_input_MassUncertainty should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsSignal_neutral_input_MassUncertainty'
    endif

    theo(1)%particle(Hneut)%dM = dMh

end subroutine HiggsSignals_neutral_input_MassUncertainty

!> Switch rate normalization between observed and predicted mass.
!! @param normalize_to_refmass normalize to observed (true) or predicted (false) mass
!! @param normalize_to_refmass_outside_dmtheo normalize to observed (true) or predicted (false) mass if masses do not agree within uncertainty
subroutine setup_rate_normalization(normalize_to_refmass, normalize_to_refmass_outside_dmtheo)
    use usefulbits_hs, only: normalize_rates_to_reference_position, &
                             normalize_rates_to_reference_position_outside_dmtheo
    implicit none
    logical, intent(in) :: normalize_to_refmass
    logical, intent(in) ::  normalize_to_refmass_outside_dmtheo

    if (normalize_to_refmass) then
        write (*, *) "Using SM rate prediction at observed mass for signal strength calculation."
    else
        write (*, *) "Using SM rate prediction at predicted mass for signal strength calculation."
    endif
    if (normalize_to_refmass_outside_dmtheo) then
        write (*, *) "If predicted mass and observed mass do not agree within theory uncertainty:", &
            " SM rate prediction at observed mass is used for signal strength calculation."
    else
        write (*, *) "If predicted mass and observed mass do not agree within theory uncertainty:", &
            " SM rate prediction at predicted mass is used for signal strength calculation."
    endif
    normalize_rates_to_reference_position = normalize_to_refmass
    normalize_rates_to_reference_position_outside_dmtheo = normalize_to_refmass_outside_dmtheo

end subroutine setup_rate_normalization

!> Switch correlations on/off.
!! @param corr_mu use correlations in signal strength measurements
!! @param corr_mh use correlations in mass measurements
subroutine setup_correlations(corr_mu, corr_mh)
    use usefulbits_hs, only: correlations_mu, correlations_mh
    implicit none
    logical, intent(in) :: corr_mu, corr_mh

    if (.not. corr_mu) then
        write (*, *) 'Correlations in signal strength observables are switched off.'
    endif
    correlations_mu = corr_mu

    if (.not. corr_mh) then
        write (*, *) 'Correlations in Higgs mass observables are switched off.'
    endif
    correlations_mh = corr_mh
end subroutine setup_correlations

!> Symmetrize measured rate uncertainties.
!! @param symm switch, default: false
subroutine setup_symmetricerrors(symm)
    use usefulbits_hs, only: symmetricerrors
    implicit none
    logical, intent(in) :: symm
    if (symm) then
        write (*, *) "Using averaged (symmetrical) experimental rate uncertainties."
    else
        write (*, *) "Using original (asymmetrical) experimental rate uncertainties."
    endif
    symmetricerrors = symm
end subroutine setup_symmetricerrors

!> Treat measured rate uncertainties as absolute uncertainties.
!! Treated as relative uncertainties otherwise.
!! @param absol switch, default: false
subroutine setup_absolute_errors(absol)
    use usefulbits_hs, only: absolute_errors
    implicit none
    logical, intent(in) :: absol
    if (absol) then
        write (*, *) "Using absolute experimental rate uncertainties."
    else
        write (*, *) "Using relative experimental rate uncertainties."
    endif
    absolute_errors = absol
end subroutine setup_absolute_errors

!> Sets the mass range (in standard deviations) for assignment of a Higgs to peak observables.
!! Sets value of usefulbits_hs::assignmentrange.
!! @param range number of standard deviations
subroutine setup_assignmentrange(range)
    use usefulbits_hs, only: assignmentrange, assignmentrange_massobs, pdf
    implicit none

    double precision, intent(in) :: range

    if (range .le. 0.0D0) then
        write (*, *) "Error: Bad assignment range ", range
        write (*, *) "Keeping the value ", assignmentrange
    else
        assignmentrange = range
        assignmentrange_massobs = range
    endif
    if (assignmentrange .gt. 1.0D0 .and. pdf .eq. 1) then
        write (*, *) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
    endif
end subroutine setup_assignmentrange

!> Sets the mass range (in standard deviations) for assignment of a Higgs to mass observables.
!! Sets value of usefulbits_hs::assignmentrange_massobs.
!! @param range number of standard deviations
subroutine setup_assignmentrange_massobservables(range)
    use usefulbits_hs, only: assignmentrange_massobs, pdf
    implicit none

    double precision, intent(in) :: range

    if (range .le. 0.0D0) then
        write (*, *) "Error: Bad assignment range ", range
        write (*, *) "Keeping the value ", assignmentrange_massobs
    else
        assignmentrange_massobs = range
    endif

    if (assignmentrange_massobs .gt. 1.0D0 .and. pdf .eq. 1) then
        write (*, *) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
    endif
end subroutine setup_assignmentrange_massobservables

!> Sets the number of free parameters for the P-value calculation.
!! @param Np number of free parameters
subroutine setup_nparam(Np)
    use usefulbits_hs, only: Nparam
    implicit none
    integer, intent(in) :: Np
    Nparam = Np
end subroutine setup_nparam

!> Sets the number of iterations for the Higgs-to-peak-assignment.
!! @param iter number of iterations
subroutine setup_Higgs_to_peaks_assignment_iterations(iter)
    use usefulbits_hs, only: iterations
    implicit none
    integer, intent(in) :: iter
    iterations = iter
end subroutine setup_Higgs_to_peaks_assignment_iterations

!> Passes a logical switch to whether the central observed values are set to the SM expectation.
!! @set_sm integer value to be set to 0 (use observed values) or 1 (use SM expected values)
subroutine setup_SM_expected(set_sm)
    use usefulbits_hs, only : set_SM_expected
    implicit none
    integer, intent(in) :: set_sm
    if(set_sm.eq.0) then
      set_SM_expected = .False.
      write(*,*) "Central values of measurements are as observed."
    else if (set_sm.eq.1) then
      set_SM_expected = .True.
      write(*,*) "Central values of measurements are set to the SM prediction!"
    else
      write(*,*) "WARNING: subroutine setup_SM_expected requires input of 0 (.False.) or 1. (.True.)!"
    endif
end subroutine setup_SM_expected
!> Controls the level of output
!! @param level
!! value | effect
!! ------|-------
!!   0   | silent mode
!!   1   | output for each analysis with observables and their predictions
!!   2   | + detailed analysis information
!!   3   | + creates files peak_information.txt and peak_massesandrates.txt
subroutine setup_output_level(level)
    use usefulbits_hs, only: output_level, additional_output
    implicit none
    integer, intent(in) :: level

    if (level .eq. 0 .or. level .eq. 1 .or. level .eq. 2 .or. level .eq. 3) then
        output_level = level
    else
        stop 'Error in subroutine setup_output_level: level not equal to 0,1,2 or 3.'
    endif
    if (level .eq. 3) additional_output = .True.

end subroutine setup_output_level

!> Sets the probability density function for the Higgs mass uncertainty
!! @param pdf_in
!! value | effect
!! ------|------
!!   1   | box-shaped pdf
!!   2   | Gaussian pdf
!!   3   | box-shaped theory error + Gaussian experimental pdf
subroutine setup_pdf(pdf_in)
    use usefulbits_hs, only: pdf
    implicit none
    integer, intent(in) :: pdf_in
    character(LEN=13) :: pdf_desc(3) = (/'box         ', 'Gaussian    ', 'box+Gaussian'/)

    pdf = pdf_in
    if ((pdf .eq. 1) .or. (pdf .eq. 2) .or. (pdf .eq. 3)) then
        write (*, '(1X,A,A,1I1,A)') 'Use a '//trim(pdf_desc(pdf))//' probability density function ', &
            'for the Higgs mass(es) (pdf=', pdf, ')'
    endif
end subroutine setup_pdf

!> Assigns toy values for the mass and mu value to a peak observable.
!! @param ID observable ID
!! @param mu_obs toy value for mu
!! @param mh_obs toy valye for mh
subroutine assign_toyvalues_to_peak(ID, mu_obs, mh_obs)
    use usefulbits_hs, only: obs, usetoys
    implicit none

    integer, intent(in) :: ID
    double precision, intent(in) :: mh_obs, mu_obs
    integer :: pos, ii

    pos = -1
    do ii = lbound(obs, dim=1), ubound(obs, dim=1)
        if (obs(ii)%id .eq. ID) then
            pos = ii
            exit
        endif
    enddo

    if (pos .ne. -1) then
        obs(pos)%peak%mpeak = mh_obs
        obs(pos)%peak%mu = mu_obs
        usetoys = .True.
    else
        write (*, *) "WARNING in assign_toyvalues_to_peak: ID unknown."
    endif
end subroutine assign_toyvalues_to_peak

!> Modify the efficiencies of an observable.
!! You can first employ the subroutine io::get_peak_channels()
!! to obtain the relevant channel information of the observable.
!! @param ID observable ID
!! @param Nc number of channels of this observable
!! @param eff_ratios scaling factors for the efficiencies in each channel
subroutine assign_modelefficiencies_to_peak(ID, Nc, eff_ratios)
    use usefulbits_hs, only: obs
    implicit none

    integer, intent(in) :: ID, Nc
    double precision, dimension(Nc), intent(in) :: eff_ratios
    integer :: pos, ii

    pos = -1
    do ii = lbound(obs, dim=1), ubound(obs, dim=1)
        if (obs(ii)%id .eq. ID) then
            pos = ii
            exit
        endif
    enddo

    if (pos .ne. -1) then
        if (size(eff_ratios, dim=1) .ne. obs(pos)%table%Nc) then
            write (*, *) "WARNING in assign modelefficiencies_to_peak: Number of channels (", &
                size(eff_ratios, dim=1), "!=", obs(pos)%table%Nc, "does not match for observable ID = ", ID
        else
            obs(pos)%table%channel_eff_ratios = eff_ratios
        endif
    else
        write (*, *) "WARNING in assign_modelefficiencies_to_peak: ID unknown."
    endif
end subroutine assign_modelefficiencies_to_peak

!> Runs HiggsSignals for the LHC Run1 combination dataset.
!! You should also call run_higgssignals() and run_higgssignals_stxs() for
!! the complete HiggsSignals analysis.
!! @param Chisq_mu \f$\chi^2\f$ contribution of the rate observables in the LHC Run1 dataset
!! @param Chisq_mh \f$\chi^2\f$ contribution of the mass observables in the LHC Run1 dataset
!! @param Chisq combined \f$\chi^2\f$ contribution of the LHC Run1 dataset
!! @param nobs number of observables involved
!! @param Pvalue P-value calculated from the total \f$\chi^2\f$ and nobs - Np (see setup_nparam())
subroutine run_HiggsSignals_LHC_Run1_combination(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    use usefulbits, only: theo, just_after_run, ndat
    use theo_manip, only: HB5_complete_theo
    use usefulbits_HS, only: HSres, output_level, Nparam
    use evaluate, only: evaluate_LHC_Run1_combination

    implicit none
    !----------------------------------------output
    integer, intent(out) ::           nobs
    double precision, intent(out) ::  Pvalue, Chisq, Chisq_mu, Chisq_mh
    !-------------------------------------internal
    integer :: n, nobs_mu, nobs_mh
    !---------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine HiggsSignals_initialize must be called first'
    endif

    call HB5_complete_theo

    do n = 1, ndat

        call evaluate_LHC_Run1_combination(theo(n), n)

        Pvalue = HSres(n)%Pvalue_LHCRun1
        Chisq = HSres(n)%Chisq_LHCRun1
        Chisq_mu = HSres(n)%Chisq_LHCRun1_mu
        Chisq_mh = HSres(n)%Chisq_LHCRun1_mh
        nobs_mu = HSres(n)%nobs_LHCRun1_mu
        nobs_mh = HSres(n)%nobs_LHCRun1_mh
        nobs = nobs_mu + nobs_mh

        if (output_level .ne. 0) then
            write (*, *)
            write (*, *) '#*************************************************************************#'
            write (*, *) '#         HIGGSSIGNALS RESULTS (LHC ATLAS + CMS Run1 combination)         #'
            write (*, *) '#*************************************************************************#'
            write (*, '(A55,F21.8)') 'chi^2 from signal rate observables = ', Chisq_mu
            write (*, '(A55,F21.8)') 'chi^2 from Higgs mass observables = ', Chisq_mh
            write (*, '(A55,F21.8)') 'chi^2 (total) = ', Chisq
            write (*, '(A55,I21)') 'Number of rate observables = ', nobs_mu
            write (*, '(A55,I21)') 'Number of mass observables = ', nobs_mh
            write (*, '(A55,I21)') 'Number of observables (total) = ', nobs
            write (*, '(A48,I3,A4,F21.8)') 'Probability (ndf =', nobs - Nparam, ') = ', Pvalue
            write (*, *) '#*************************************************************************#'
            write (*, *)
        endif

    enddo

    just_after_run = .True.
end subroutine run_HiggsSignals_LHC_Run1_combination

!> Runs HiggsSignals for the STXS observables.
!! You should also call run_higgssignals() and run_higgssignals_lhc_run1_combination() for
!! the complete HiggsSignals analysis.
!! @param Chisq_STXS_rates \f$\chi^2\f$ contribution of the STXS rate observables
!! @param Chisq_STXS_mh \f$\chi^2\f$ contribution of the STXS mass observables
!! @param Chisq_STXS combined \f$\chi^2\f$ of the STXS observables
!! @param nobs_STXS number of observables involved
!! @param Pvalue_STXS P-value calculated from the total \f$\chi^2\f$ and nobs - Np (see setup_nparam())
subroutine run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)
    use STXS, only: evaluate_model_for_STXS, get_chisq_from_STXS, print_STXS, &
                    get_number_of_STXS_observables, STXSlist, print_STXS_to_file
    use usefulbits, only: theo, ndat, vsmall
    use usefulbits_hs, only: HSres, output_level, Nparam
    use theo_manip, only: HB5_complete_theo
    use numerics, only: gammp

    implicit none
    double precision, intent(out) :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, Pvalue_STXS
    integer, intent(out) :: nobs_STXS
    double precision :: Pvalue
    integer :: nobs_STXS_rates, nobs_STXS_mh, i, n

    call HB5_complete_theo

    Chisq_STXS_mh = 0.0D0

    do n = 1, ndat

        do i = lbound(STXSlist, dim=1), ubound(STXSlist, dim=1)
            call evaluate_model_for_STXS(STXSlist(i), theo(n))
        enddo
        call get_chisq_from_STXS(Chisq_STXS_rates, Pvalue_STXS)
        call get_number_of_STXS_observables(nobs_STXS_rates, nobs_STXS_mh)

        nobs_STXS = nobs_STXS_rates + nobs_STXS_mh
        ! Add routine for possible mh-observable in STXS here!

        Chisq_STXS = Chisq_STXS_rates + Chisq_STXS_mh

        HSres(n)%Chisq_STXS_rates = Chisq_STXS_rates
        HSres(n)%Chisq_STXS_mh = Chisq_STXS_mh
        HSres(n)%Chisq_STXS = Chisq_STXS
        HSres(n)%nobs_STXS_rates = nobs_STXS_rates
        HSres(n)%nobs_STXS_mh = nobs_STXS_mh
        HSres(n)%nobs_STXS = nobs_STXS

        Pvalue = 1.0D0
        if (Chisq_STXS .gt. vsmall .and. (nobs_STXS - Nparam) .gt. 0) then
            Pvalue = 1 - gammp(dble(nobs_STXS - Nparam)/2, Chisq_STXS/2)
        endif

        HSres(n)%Pvalue_STXS = Pvalue

    enddo

    if (output_level .eq. 1) call print_STXS
    if (output_level .eq. 3) then
        call print_STXS_to_file
    endif

    if (output_level .ne. 0) then
        write (*, *)
        write (*, *) '#*************************************************************************#'
        write (*, *) '#                HIGGSSIGNALS RESULTS (STXS observables)                  #'
        write (*, *) '#*************************************************************************#'
        write (*, '(A55,F21.8)') 'chi^2 (signal rate) from STXS observables = ', Chisq_STXS_rates
        write (*, '(A55,F21.8)') 'chi^2 (Higgs mass) from STXS observables = ', Chisq_STXS_mh
        write (*, '(A55,F21.8)') 'chi^2 (total) = ', Chisq_STXS
        write (*, '(A55,I21)') 'Number of STXS rate observables = ', nobs_STXS_rates
        write (*, '(A55,I21)') 'Number of STXS mass observables = ', nobs_STXS_mh
        write (*, '(A55,I21)') 'Number of STXS observables (total) = ', nobs_STXS
        write (*, '(A48,I3,A4,F21.8)') 'Probability (ndf =', nobs_STXS - Nparam, ') = ', Pvalue
        write (*, *) '#*************************************************************************#'
        write (*, *)
    endif

end subroutine run_HiggsSignals_STXS

!> Runs Higgssignals for the peak observables.
!! @param Chisq_mu \f$\chi^2\f$ contribution of the rate observables
!! @param Chisq_mh \f$\chi^2\f$ contribution of the mass observables
!! @param Chisq combined \f$\chi^2\f$ of the peak observables
!! @param nobs number of observables involved
!! @param Pvalue P-value calculated from the total \f$\chi^2\f$ and nobs - Np (see setup_nparam())
subroutine run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    use usefulbits, only: theo, just_after_run, ndat
    use usefulbits_HS, only: HSres, output_level, usescalefactor, Nparam, Exptdir
    use channels, only: check_channels
    use theo_manip, only: HB5_complete_theo
    use evaluate, only: evaluate_model

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------output
    integer, intent(out) ::           nobs
    double precision, intent(out) ::  Pvalue, Chisq, Chisq_mu, Chisq_mh
    !-------------------------------------internal
    integer :: n
    logical :: debug = .False.
    !---------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine HiggsSignals_initialize must be called first'
    endif
    if (debug) write (*, *) 'manipulating input...'; call flush (6)

    call HB5_complete_theo

    if (debug) write (*, *) 'compare each model to the experimental data...'; call flush (6)

    do n = 1, ndat

!  call recalculate_theo_for_datapoint(n)

        call evaluate_model(theo(n), n)

        Pvalue = HSres(n)%Pvalue_peak
        Chisq = HSres(n)%Chisq_peak
        Chisq_mu = HSres(n)%Chisq_peak_mu
        Chisq_mh = HSres(n)%Chisq_peak_mh
        nobs = HSres(n)%nobs_peak

        if (output_level .ne. 0) then
            write (*, *)
            write (*, *) '#*************************************************************************#'
            write (*, *) '#   HIGGSSIGNALS RESULTS (', trim(adjustl(Exptdir)), ') -- peak observables                 #'
            write (*, *) '#*************************************************************************#'
            write (*, '(A55,F21.8)') 'chi^2 (signal strength) from peak observables = ', &
                HSres(n)%Chisq_peak_mu
            write (*, '(A55,F21.8)') 'chi^2 (Higgs mass) from peak observables = ', HSres(n)%Chisq_peak_mh
            write (*, '(A55,F21.8)') 'chi^2 (total) from peak observables = ', HSres(n)%Chisq
            write (*, '(A55,I21)') 'Number of signal strength peak observables = ', &
                HSres(n)%nobs_peak_mu
            write (*, '(A55,I21)') 'Number of Higgs mass peak observables = ', HSres(n)%nobs_peak_mh
            write (*, '(A55,I21)') 'Number of peak observables (total) = ', HSres(n)%nobs_peak
            write (*, '(A48,I3,A4,F21.8)') 'Probability (ndf =', HSres(n)%nobs - Nparam, &
                ') using peak observables = ', HSres(n)%Pvalue_peak
            write (*, *) '#*************************************************************************#'
            write (*, *)
        endif

    enddo

    just_after_run = .True.
    usescalefactor = .False.

end subroutine run_HiggsSignals

!> Runs HiggsSignals with all available measurements.
!! Calls run_higgssignals(), run_higgssignals_stxs(),
!! run_higgssignals_lhc_run1_combination() and complete_hs_results().
!! Returns the combined results.
!! @param Chisq_mu \f$\chi^2\f$ contribution of all rate observables
!! @param Chisq_mh \f$\chi^2\f$ contribution of all mass observables
!! @param Chisq combined \f$\chi^2\f$ of all observables
!! @param nobs number of observables involved
!! @param Pvalue P-value calculated from the total \f$\chi^2\f$ and nobs - Np (see setup_nparam())
subroutine run_HiggsSignals_full(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    use usefulbits_HS, only: HSres
    use usefulbits, only: ndat

    implicit none
    double precision, intent(out) :: Chisq_mu, Chisq_mh, Chisq, Pvalue
    integer, intent(out) :: nobs

    call run_HiggsSignals(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    call run_HiggsSignals_LHC_Run1_combination(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    call run_HiggsSignals_STXS(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
    call complete_HS_results()

    Chisq_mu = HSres(ndat)%Chisq_mu
    Chisq_mh = HSres(ndat)%Chisq_mh
    Chisq = HSres(ndat)%Chisq
    nobs = HSres(ndat)%nobs
    Pvalue = HSres(ndat)%Pvalue
end subroutine run_HiggsSignals_full

#if defined(__GFORTRAN__) && __GNUC__>=7
!> Sets an upper limit on the scale factor s of the model prediction
!! as described in section 3.3.2 of the HiggsSignals manual.
!! Like run_higgssignals_full, this uses all available measurements.
!!
!! The implementation needs some Fortran 2007 features (local procedures as arguments)
!! and is therefore only available with newer compilers.
!!
!! @param sHat best fit value for the scalar factor s
!! @param sUp CL_s+b based upper limit on s
!! @param CLs the CL_s value for s = 1
subroutine set_HiggsSignals_limit(sHat, sUp, CLs)
    use usefulbits, only: theo, just_after_run, ndat
    use usefulbits_HS, only: HSres, wald_approximate_chisq
    use theo_manip, only: HB5_complete_theo
    use likelihoods, only: cumnor

    implicit none
    double precision, intent(out) :: sHat, sUp, CLs
    double precision:: sigmaS, CLsb, CLb, dummy
    double precision, parameter :: lowerSBound = -2D0, upperSBound = 2D0
    integer n
    if (.not. allocated(theo)) then
        stop 'subroutine HiggsSignals_initialize must be called first'
    endif

    call HB5_complete_theo

    do n = 1, ndat
        HSres(n)%limitSetting = .true.
        call wald_approximate_chisq(lowerSBound, upperSBound, chisq, HSres(n)%sHat, sigmaS)
        HSres(n)%sUp = HSres(n)%sHat + 1.64D0*sigmaS
        if (HSres(n)%sHat .ge. 1D0) then
            HSres(n)%CLs = 1D0
        else
            call cumnor((1 - HSres(n)%sHat)/sigmaS, dummy, CLsb)
            call cumnor((HSres(n)%sHat - 1)/sigmaS, dummy, CLb)
            HSres(n)%CLs = CLsb/CLb
        endif
    enddo

    sHat = HSres(ndat)%sHat
    sUp = HSres(ndat)%sUp
    CLs = HSres(ndat)%CLs

    just_after_run = .True.

contains
    function chisq(s)
        use STXS, only: evaluate_model_for_STXS, get_chisq_from_STXS, STXSlist
        use usefulbits, only: theo
        use evaluate, only: evaluate_model, evaluate_LHC_Run1_combination
        use usefulbits_HS, only: HSres
        implicit none
        double precision, intent(in) :: s
        double precision :: chisq
        double precision :: dummy
        integer :: i
        call evaluate_model(theo(n), n, s)
        call evaluate_LHC_Run1_combination(theo(n), n, s)
        do i = lbound(STXSlist, dim=1), ubound(STXSlist, dim=1)
            call evaluate_model_for_STXS(STXSlist(i), theo(n), s)
        enddo
        call get_chisq_from_STXS(chisq, dummy)
        chisq = chisq + HSres(n)%Chisq_LHCRun1_mu + HSres(n)%Chisq_peak_mu
    end function
end subroutine set_HiggsSignals_limit
#endif

!> Returns SM normalized signal rates of some relevant channels.
!! Does not include efficiencies.
!! @param ii which Higgs boson
!! @param collider which collider
!! value | collider
!! ------|---------
!!   1   | Tevatron
!!   2   | LHC7
!!   3   | LHC8
!!   4   | LHC13
!! @param R_H_WW signal rate in the \f$h\to WW\f$ channel
!! @param R_H_ZZ signal rate in the \f$h\to ZZ\f$ channel
!! @param R_H_gaga signal rate in the \f$h\to \gamma\gamma\f$ channel
!! @param R_H_tautau signal rate in the \f$h\to \tau\tau\f$ channel
!! @param R_H_bb signal rate in the \f$h\to bb\f$ channel
!! @param R_VH_bb signal rate in the \f$Vh\to bb\f$ channel
subroutine get_Rvalues(ii, collider, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
    integer, intent(in) :: ii, collider
    double precision, intent(out) :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb

    call get_rates(ii, collider, 5, (/"1.2", "2.2", "3.2", "4.2", "5.2"/), R_H_WW)
    call get_rates(ii, collider, 5, (/"1.3", "2.3", "3.3", "4.3", "5.3"/), R_H_ZZ)
    call get_rates(ii, collider, 5, (/"1.1", "2.1", "3.1", "4.1", "5.1"/), R_H_gaga)
    call get_rates(ii, collider, 5, (/"1.4", "2.4", "3.4", "4.4", "5.4"/), R_H_tautau)
    call get_rates(ii, collider, 5, (/"1.5", "2.5", "3.5", "4.5", "5.5"/), R_H_bb)
    call get_rates(ii, collider, 2, (/"3.5", "4.5"/), R_VH_bb)
end subroutine get_Rvalues

!> Get a SM normalized signal rate.
!! Does not include efficiencies
!! @param ii which Higgs boson
!! @param collider which collider
!! value | collider
!! ------|---------
!!   1   | Tevatron
!!   2   | LHC7
!!   3   | LHC8
!!   4   | LHC13
!! @param Nchannels number of subchanneld to include in the rate
!! @param IDchannels_str IDs of the [subchannels](doc/channels.md) to include in the form "p.d" (array of size Nchannels)
!! @param rate the resulting rate
subroutine get_rates(ii, collider, Nchannels, IDchannels_str, rate)
    use usefulbits, only: theo, np, Hneut
    use usefulbits_HS, only: mutable
    use evaluate, only: get_channelrates

    integer, intent(in) :: ii, collider, Nchannels
    character(LEN=*), dimension(Nchannels), intent(in) :: IDchannels_str
    double precision, intent(out) :: rate
!-Internal
    type(mutable) :: dummytable
    integer :: i, id, posperiod

!-Initialize a dummy mutable in order to run get_channelrates for the channels we want.
    if (collider .eq. 1) then
        dummytable%collider = 'TEV'
    else if (collider .eq. 2) then
        dummytable%collider = 'LHC'
        dummytable%energy = 7.0D0
    else if (collider .eq. 3) then
        dummytable%collider = 'LHC'
        dummytable%energy = 8.0D0
    else if (collider .eq. 4) then
        dummytable%collider = 'LHC'
        dummytable%energy = 13.0D0
    else
        write (*, *) 'WARNING: collider experiment for get_rates unknown.'
        continue
    endif

    dummytable%id = 999999
    dummytable%particle_x = 1
    dummytable%Nc = Nchannels
    allocate (dummytable%mass(1))
    dummytable%mass(1) = theo(1)%particle(dummytable%particle_x)%M(ii)
!  allocate(dummytable%channel_id(Nchannels))
    allocate (dummytable%channel_p_id(Nchannels))
    allocate (dummytable%channel_d_id(Nchannels))
    allocate (dummytable%channel_eff(Nchannels))
    allocate (dummytable%channel_eff_ratios(Nchannels))
!-Set all efficiencies equal:
    dummytable%channel_eff = 1.0D0
    dummytable%channel_eff_ratios = 1.0D0
    allocate (dummytable%channel_description(Nchannels, 2))
    allocate (dummytable%channel_w(Nchannels, np(Hneut)))
    allocate (dummytable%channel_w_corrected_eff(Nchannels, np(Hneut)))
    allocate (dummytable%channel_systSM(Nchannels, np(Hneut)))
    allocate (dummytable%channel_syst(Nchannels, np(Hneut)))
    allocate (dummytable%channel_mu(Nchannels, np(Hneut)))

! do i = 1,Nchannels
!  write(*,*) i, IDchannels_str(i)
! enddo

    do i = 1, Nchannels
        posperiod = index(IDchannels_str(i), '.')
!   write(*,*) IDchannels_str(i)
        if (posperiod .eq. 0) then
            if (len(trim(adjustl(IDchannels_str(i)))) .eq. 2) then
                read (IDchannels_str(i), *) id
                dummytable%channel_p_id(i) = int((id - modulo(id, 10))/dble(10))
                dummytable%channel_d_id(i) = modulo(id, 10)
            else
                stop " Error in get_rates: Cannot handle channel IDs!"
            endif
        else
!    write(*,*) dummytable%channel_p_id(i), dummytable%channel_d_id(i)
            read (IDchannels_str(i) (:posperiod - 1), *) dummytable%channel_p_id(i)
            read (IDchannels_str(i) (posperiod + 1:), *) dummytable%channel_d_id(i)
        endif
    enddo

!  write(*,*) dummytable%channel_p_id, dummytable%channel_d_id
    call get_channelrates(ii, theo(1), dummytable)
    rate = 0.0D0
    do i = lbound(dummytable%channel_mu, dim=1), ubound(dummytable%channel_mu, dim=1)
        rate = rate + dummytable%channel_mu(i, ii)*dummytable%channel_w(i, ii)
    enddo

    deallocate (dummytable%channel_p_id, dummytable%channel_d_id, dummytable%channel_eff, &
                dummytable%channel_w, dummytable%channel_systSM, dummytable%channel_syst, &
                dummytable%channel_mu, dummytable%channel_eff_ratios, dummytable%channel_description, &
                dummytable%channel_w_corrected_eff, dummytable%mass)

end subroutine get_rates

!> Calculates the P-value for the total \f$\chi^2\f$ value and the given number of free parameters.
!! *Always call complete_hs_results() before using this.*
!! @param nparam number of free parameters, overrides any value set by setup_nparam()
!! @param Pvalue the resulting P-value
subroutine get_Pvalue(nparam, Pvalue)
    use usefulbits, only: vsmall
    use usefulbits_hs, only: HSres
    use numerics
    implicit none
    integer, intent(in) :: nparam
    double precision, intent(out) :: Pvalue

    if (allocated(HSres)) then
        if (HSres(1)%Chisq .gt. vsmall .and. (HSres(1)%nobs - nparam) .gt. 0) then
            HSres(1)%Pvalue = 1 - gammp(dble(HSres(1)%nobs - nparam)/2, HSres(1)%Chisq/2)
        endif
    else
        write (*, *) "Warning: subroutine get_Pvalue should be called after run_HiggsSignals."
    endif

    Pvalue = HSres(1)%Pvalue

end subroutine get_Pvalue

!> Return the neutral Higgs masses and uncertainties.
!! May be useful for more detailed output.
!! @param Mh array of the Higgs masses
!! @param dMh array of the Higgs mass uncertainties
subroutine get_neutral_Higgs_masses(Mh, dMh)
    use usefulbits, only: theo, np, Hneut

    implicit none
    double precision, intent(out) :: Mh(np(Hneut)), dMh(np(Hneut))

    if (.not. allocated(theo)) then
        stop 'No model information given!'
    endif
    if (np(Hneut) .eq. 0) then
        write (*, *) 'Cannot access the neutral Higgs boson masses'
        write (*, *) 'because np(Hneut) == 0.'
        stop 'error in subroutine get_neutral_Higgs_masses'
    endif

    Mh = theo(1)%particle(Hneut)%M
    dMh = theo(1)%particle(Hneut)%dM

end subroutine get_neutral_Higgs_masses

!> Combines the results of the separate runs.
!! Call this after calling run_higgssignals(), run_higgssignals_stxs(),
!! and/or run_higgssignals_lhc_run1_combination() if you want to use the
!! builtin output routines or get combined results (e.g. through get_pvalue()).
subroutine complete_HS_results()
    use usefulbits, only: just_after_run, ndat, vsmall
    use usefulbits_HS, only: HSres, Nparam
    use numerics, only: gammp

    implicit none
    integer :: n

    if (just_after_run) then
        do n = 1, ndat

            HSres(n)%Chisq_mu = HSres(n)%Chisq_peak_mu + & !HSres(n)%Chisq_mpred + &
                                HSres(n)%Chisq_STXS_rates + HSres(n)%Chisq_LHCRun1_mu

            HSres(n)%Chisq_mh = HSres(n)%Chisq_peak_mh + HSres(n)%Chisq_LHCRun1_mh + &
                                HSres(n)%Chisq_STXS_mh

            HSres(n)%Chisq_STXS = HSres(n)%Chisq_STXS_rates + HSres(n)%Chisq_STXS_mh

            HSres(n)%Chisq_peak = HSres(n)%Chisq_peak_mu + HSres(n)%Chisq_peak_mh

            HSres(n)%Chisq_LHCRun1 = HSres(n)%Chisq_LHCRun1_mu + HSres(n)%Chisq_LHCRun1_mh

            HSres(n)%Chisq = HSres(n)%Chisq_mu + HSres(n)%Chisq_mh

            HSres(n)%nobs_mu = HSres(n)%nobs_peak_mu + &!HSres(n)%nobs_mpred + &
                               HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_STXS_rates

            HSres(n)%nobs_mh = HSres(n)%nobs_peak_mh + HSres(n)%nobs_LHCRun1_mh + &
                               HSres(n)%nobs_STXS_mh

            HSres(n)%nobs_peak = HSres(n)%nobs_peak_mu + HSres(n)%nobs_peak_mh

            HSres(n)%nobs_STXS = HSres(n)%nobs_STXS_rates + HSres(n)%nobs_STXS_mh

            HSres(n)%nobs_LHCRun1 = HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_LHCRun1_mh

            HSres(n)%nobs = HSres(n)%nobs_mu + HSres(n)%nobs_mh

            if (HSres(n)%Chisq .gt. vsmall .and. (HSres(n)%nobs - Nparam) .gt. 0) then
                HSres(n)%Pvalue = 1 - gammp(dble(HSres(n)%nobs - Nparam)/2.0D0, HSres(n)%Chisq/2.0D0)
            endif

            if (HSres(n)%Chisq_peak .gt. vsmall .and. (HSres(n)%nobs_peak - Nparam) .gt. 0) then
                HSres(n)%Pvalue_peak = 1 - gammp(dble(HSres(n)%nobs_peak - Nparam)/2.0D0, HSres(n)%Chisq_peak/2.0D0)
            endif

            if (HSres(n)%Chisq_LHCRun1 .gt. vsmall .and. (HSres(n)%nobs_LHCRun1 - Nparam) .gt. 0) then
                HSres(n)%Pvalue_LHCRun1 = 1 - gammp(dble(HSres(n)%nobs_LHCRun1 - Nparam)/2.0D0, HSres(n)%Chisq_LHCRun1/2.0D0)
            endif

            if (HSres(n)%Chisq_STXS .gt. vsmall .and. (HSres(n)%nobs_STXS - Nparam) .gt. 0) then
                HSres(n)%Pvalue_STXS = 1 - gammp(dble(HSres(n)%nobs_STXS - Nparam)/2.0D0, HSres(n)%Chisq_STXS/2.0D0)
            endif

        enddo
    else
        write (*, *) "Warning: complete_HS_results was called but just_after_run is", just_after_run
    endif
end subroutine complete_HS_results

!> Cleans up HiggsSignals.
!! Call this when you are done using HiggsSignals, deallocates array, closes files, etc.
subroutine finish_HiggsSignals
    use usefulbits, only: deallocate_usefulbits, debug, theo, debug, &
                          file_id_debug1, file_id_debug2
    use S95tables, only: deallocate_Exptranges
    use theory_BRfunctions, only: deallocate_BRSM
    use datatables, only: deallocate_observables
    use usefulbits_HS, only: deallocate_usefulbits_HS, analyses
    use mc_chisq, only: deallocate_mc_observables
    use install_data_HS

    if (debug) then
        close (file_id_debug2)
        close (file_id_debug1)
    endif

    if (debug) write (*, *) 'finishing off...'; call flush (6)
    if (.not. allocated(theo)) then
        if (debug) write (*, *) "HiggsBounds/HiggsSignals internal structure already deallocated!"
    else
        call deallocate_BRSM
        call deallocate_Exptranges
        call deallocate_usefulbits
    endif
    call deallocate_mc_observables
    call deallocate_observables
    if (allocated(analyses)) deallocate (analyses)
    call deallocate_usefulbits_HS
    call system('rm -f HS_analyses.txt')

    if (debug) write (*, *) 'finished'; call flush (6)

end subroutine finish_HiggsSignals
