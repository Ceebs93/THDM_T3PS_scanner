!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------

!> Contains HiggsSignals configuration parameters and core routines.
module usefulbits_HS
!--------------------------------------------------------------------
    implicit none

    integer, parameter :: f_dmh = 94
    character(LEN=100) :: Exptdir
!--------------------------------------------------------------------
!------------------    User Control parameters    -------------------
!--------------------------------------------------------------------
! Note: values can be changed with specific user subroutines.

    logical :: usescalefactor = .False.
    logical :: symmetricerrors = .False.
    logical :: useaveragemass = .True.
    logical :: EWcorr = .True. ! Whether to include EW-corrections in SM reference cross section (in STXS)
    logical :: correlations_mu = .True.
    logical :: correlations_mh = .True.
    logical :: normalize_rates_to_reference_position = .False.
    logical :: normalize_rates_to_reference_position_outside_dmtheo = .False.
! If false, the reference mass value will be mobs +- dmtheo (depending on predicted mass)
    logical :: LHC_combination_run1_SMXS_from_paper = .False.

!> If true, this sets all central values of the measurement to the SM prediction
    logical :: set_SM_expected = .False.

    !> The mass range (in standard deviations), in which the Higgs is forced
    !! to be assigned to a peak observable.
    !! Can be modified using setup_assignmentrange().
    double precision :: assignmentrange = 1.0D0

    !> The mass range (in standard deviations), in which the Higgs is forced
    !! to be assigned to to peak observables, which have a mass measurement.
    !! Can be modified using setup_assignmentrange_massobservables().
    double precision :: assignmentrange_massobs = 2.0D0

    !> The mass range (in standard deviations), in which the Higgs is forced
    !! to be assigned to an STXS observable.
    !! Can be modified using setup_assignmentrange_stxs().
    ! double precision :: assignmentrange_STXS = 1.0D0

    !> The minimal signal strength to be considered for an assignment to an observable.
    double precision :: mu_cutoff_for_assignment = 1.0D-03

    !> In case of multiple Higgs assignments in the box pdf case, this is
    !! the relative contribution of a Higgs boson to be counted as contributing in the
    !! mass separation chi^2.
    double precision :: significantcontribution = 0.01D0

    integer :: output_level = 0
    integer :: iterations = 0  ! default value: 0
    ! 1 -> try to assign as many Higgs bosons as possible to
    ! the observable, Higgs-to-peak assignment is based on
    ! Higg mass covariance matrices with maximal
    ! correlations.
    ! >1 -> use the covariance matrix of previous iteration.
    integer :: pdf = 2  ! default value: 2
    ! will automatically be set to 2 if not changed by the user
    ! via using subroutine set_pdf before.
    ! (1,2,3) = (box, gaussian, theory-box + exp-gaussian)
    integer :: Nparam = 0  ! Number of free model parameters (entering the Pvalue)
    ! Can be specified directly here or via the subroutine
    ! setup_nparam

!------------------- For internal debugging/testing -----------------
    logical :: withcorrexpsyst = .True. !(correlated experimental systematics)
    logical :: anticorrmu = .True.
    logical :: anticorrmh = .True.
    logical :: STXS_addtheoryuncertainty_corr = .False. ! Adds theory uncertainty naively using provided correlation matrix (only off-diagonal elements)
    logical :: STXS_addtheoryuncertainty_uncorr = .False. ! Adds theory uncertainty naively using provided correlation matrix (only diagonal elements)
    logical :: STXS_addtheoryuncertainty_full = .True. ! Using internal correlation information on inclusive XS to be folded into STXS chi^2 test.
    logical :: STXS_rescaleSMuncertainty = .True.  ! Rescale SM theory uncertainty by model-predicted signal strength modifier
    logical :: STXS_convert_abs_to_SMnorm = .False. ! Converts absolute measurements to signal strength modifiers before chi^2 calculation (default: .False.)
!-- sleeping features --
    logical :: useSMweights = .False.
    logical :: minimalchisq = .False.
    logical :: maximalchisq = .False.
    logical :: additional_output = .False. !(outdated, to be removed!)
    logical :: useSMtest = .False.
    double precision :: eps
!  logical :: use_SMrate_at_reference_position_for_STXS = .False.

!-------------------- Internal Control parameters --------------------

    logical :: usetoys = .False. !(outdated, to be removed!)
    logical :: absolute_errors = .False. ! Errors treated to original mu value, or toy-value
    logical :: SLHAdetailed = .False.
    logical :: newSLHAfile = .False.
    logical :: THU_included = .True.

!--------------------------------------------------------------------
    integer :: nanalys   !Total number of relevant analyses
    double precision, parameter :: vlarge = 1000000000000.0D0
!--------------------- Default rate uncertainties -------------------

    type rate_uncertainties
!- dCS_SM and dBR_SM for the SM
!- (from LHC HXSWG Yellow Report 3, arXiv:1307.1347)
!- dCS and dBR hold the model's rate uncertainties. Can be changed by user
!- with subroutine setup_rate_uncertainties. Default values are those of the SM.
        double precision :: dCS_SM(5) = (/0.147D0, 0.028D0, 0.037D0, 0.060D0, 0.12D0/)
        double precision :: dCS(5) = (/0.147D0, 0.028D0, 0.037D0, 0.060D0, 0.12D0/)
!- EDIT (TS 21/06/2013): Add new decay modes:
!- Channels: gammagamma, WW, ZZ, tautau, bb, Zgamma, cc, mumu, gg
        double precision :: dBR_SM(9) = (/0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0, &
                                          0.090D0, 0.122D0, 0.060D0, 0.100D0/)
        double precision :: dBR(9) = (/0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0, &
                                       0.090D0, 0.122D0, 0.060D0, 0.100D0/)
!--- IMPORTANT NOTE:
!-
!- The arrays dCS_SM, dCS, dBR_SM, dBR have been introduced in HiggsSignals-1.0.0
!- to hold the estimated theoretical uncertainties. These do not include correlations
!- via parametric uncertainties (e.g. scale, PDFs,...) or correlations in the BRs introduced
!- by the uncertainty of the total widths.
!-
!- Since HiggsSignals-1.1.0 the theoretical uncertainties for the cross sections and
!- branching ratios are evaluated with toy MC scripts including the correlations of
!- parametric error sources. The resulting covariance matrices are included per default
!- if the files "BRcov.in" and "XScov.in" are present in the main HiggsSignals directory.
!- If not, HiggsSignals will give a warning and use the old method.
!- The covariance matrices can also be re-evaluated by the user with the scripts
!- "smearErrorsBR.cpp" and "smearErrorsXS.cpp", which can be found in the directory
!- <HiggsSignals-main-directory>/supplements/
!-
!---
        logical :: BRcov_ok = .False.
        logical :: CScov_ok = .False.
        logical :: CScorr_ok = .False.
        logical :: usecov = .True.

        double precision, dimension(9, 9) :: BRcovSM = 0.0D0
        double precision, dimension(9, 9) :: BRcov = 0.0D0
        double precision, dimension(12, 12) :: CScovSM = 0.0D0
        double precision, dimension(12, 12) :: CScov = 0.0D0
        double precision, dimension(12, 12) :: CS13covSM = 0.0D0
        double precision, dimension(12, 12) :: CS13cov = 0.0D0
        double precision, dimension(12, 12) :: CScorrSM = 0.0D0
        double precision, dimension(12, 12) :: CScorr = 0.0D0
        double precision, dimension(12, 12) :: CS13corrSM = 0.0D0
        double precision, dimension(12, 12) :: CS13corr = 0.0D0

!--- ILC cross section uncertainties (under development)
!--- (none, none, WBF, ZH, ttH)
        double precision :: dCS_ILC_SM(5) = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
        double precision :: dCS_ILC(5) = (/0.0D0, 0.0D0, 0.01D0, 0.005D0, 0.01D0/)

    end type
    type(rate_uncertainties), save :: delta_rate

    type LHCrun1
        ! rate measurements
        integer :: channel_id
        double precision :: r, r_low, dr_low, r_up, dr_up
        double precision :: r_pred, dr, dr0
    end type

    type(LHCrun1), dimension(20) :: LHCrun1_rates

    double precision, dimension(20, 20) :: LHCrun1_correlationmatrix = 0.0D0

    type correlation_info
        integer :: obsID1, obsID2
        double precision :: corr
        double precision :: THUcorr = 0.0
    end type

    type(correlation_info), allocatable :: corrlist(:)

!-------------- Type definitions of internal structures --------------
    type neutHiggs
        double precision :: m, dm, mu
        integer :: mp_test  ! This variable is set to 1 (0) if the Higgs is (not) being tested in the m-pred chi^2 method.
        integer :: id
    end type
    !-Will contain info about all neutral Higgs for every considered observable, i.e.
!-neutHiggses has dimensions (number(observables),nH)
    type(neutHiggs), allocatable :: neutHiggses(:, :)

    type mutable
        integer :: id, nx, particle_x !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
        character(LEN=100) :: label
        character(LEN=100) :: desc
        character(LEN=3) :: expt
        character(LEN=10) :: collider
        character(LEN=10) :: collaboration
        double precision :: lumi, dlumi, energy  ! dlumi in %
!--TESTING correlated experimental systematics:
!  double precision, dimension(4) :: correxpsyst
!--END
        double precision :: xmax, xmin, sep, deltax
        double precision :: deltam
        character(LEN=100) :: assignmentgroup
        integer :: mhchisq
        double precision, allocatable :: mass(:)
        double precision, allocatable :: mu(:, :) ! in mu(a,b), a=row, b=1,2,3 for low,obs,up
        integer :: Nc    ! Number of channels
        character(LEN=5), allocatable :: channel_id_str(:)    ! Channels array as string, dim(Nc)
!   integer, allocatable :: channel_id(:)    ! Channels array, dim(Nc)
        integer, allocatable :: channel_p_id(:)    ! Production channels array, dim(Nc)
        integer, allocatable :: channel_d_id(:)    ! Decay channels array, dim(Nc)
        character(LEN=10), allocatable :: channel_description(:, :)
        double precision, allocatable :: channel_eff(:)     ! Channel efficiencies, dim(Nc)
        double precision, allocatable :: channel_eff_ratios(:)     ! Channel efficiency ratios (model vs. SM), dim(Nc)
        double precision, allocatable :: channel_w(:, :) ! Channel weights, dim(Nc, NHiggs)
        double precision, allocatable :: channel_w_corrected_eff(:, :) ! Channel weights, dim(Nc, NHiggs)
        double precision, allocatable :: channel_systSM(:, :) ! Channel systematics of SM, dim(Nc, NHiggs)
        double precision, allocatable :: channel_syst(:, :) ! Channel systematics, dim(Nc, NHiggs)
        double precision, allocatable :: channel_mu(:, :) ! SM normalized channel rates, dim(Nc, NHiggs)
        double precision :: eff_ref_mass    ! Reference Higgs mass for quoted efficiency
        integer :: npeaks
        double precision, allocatable :: Toys_mhobs(:)
        double precision, allocatable :: Toys_muobs(:)
        double precision :: scale_mu
    end type

    type mupeak
        integer :: id
        integer :: ilow, iup, ipeak
        double precision :: mpeak
        double precision :: dm
        double precision :: mu
        double precision :: mu_original
        double precision :: scale_mu!, scale_mh
        double precision :: dmuup, dmulow    ! Upper and lower cyan band
        double precision :: dmuup0sq, dmulow0sq    ! Cyan band squared subtracted by correlated uncertainties
        !-Peak object should contain everything needed for constructing the covariance matrices
        integer :: Nc    ! Number of channels
        integer, allocatable :: channel_id(:)    ! Channels array, dim(Nc)
        integer, allocatable :: channel_p_id(:)    ! Production channels array, dim(Nc)
        integer, allocatable :: channel_d_id(:)    ! Decay channels array, dim(Nc)
        double precision, allocatable :: channel_eff(:)     ! Channel efficiencies, dim(Nc)
        integer, allocatable :: Higgs_comb(:)    ! Assigned Higgs combination, dim(NHiggs)
        character(LEN=100) :: assignmentgroup
        type(neutHiggs), allocatable :: Higgses(:)
        integer :: domH    !  index of dominantly contributing Higgs
        integer :: NHiggs_comb    ! Number of combined Higgses
        integer :: Higgs_assignment_forced
        integer :: undo_assignment
        !--These arrays contain only the information about all Higgs bosons
        !--(need to have this for every peak separately because it can depend on the efficiencies
        !-- which are given for each peak separately)
        double precision, allocatable :: channel_w_allH(:, :)      ! Channel weights, dim(Nc, NHiggs)
        double precision, allocatable :: channel_w_corrected_eff_allH(:, :) ! Channel weights with corrected efficiencies, dim(Nc, NHiggs)
        double precision, allocatable :: channel_systSM_allH(:, :) ! Channel systematics of SM, dim(Nc, NHiggs)
        double precision, allocatable :: channel_syst_allH(:, :)   ! Channel systematics, dim(Nc, NHiggs)
        double precision, allocatable :: channel_mu_allH(:, :)     ! SM normalized channel rates, dim(Nc, NHiggs)
        !--These arrays contain only the information about the chosen Higgs combination:
        double precision, allocatable :: channel_w(:)             ! Channel weights, dim(Nc)
        double precision, allocatable :: channel_w_corrected_eff(:) ! Channel weights with corrected efficencies, dim(Nc)
        double precision, allocatable :: channel_systSM(:)        ! Channel systematics, dim(Nc)
        double precision, allocatable :: channel_syst(:)          ! Channel systematics, dim(Nc)
        double precision, allocatable :: channel_mu(:)            ! SM normalized channel rates, dim(Nc)
        double precision, allocatable :: channel_w_model(:)
        double precision :: total_mu
        double precision :: dlumi
        !-- Chisq values (mu and mh parts, total) after taking into account correlations with
        !-- other peaks:
        double precision :: chisq_mu
        double precision :: chisq_mh
        double precision :: chisq_tot
        double precision :: chisq_max
        integer :: internalnumber
    end type

    type mp_neutHiggs
!-This object is a Higgs or Higgscluster which are separately
!-tested with the predicted mass chi^2 method.
        type(neutHiggs), allocatable :: Higgses(:)
        double precision :: m, dm, mu
        integer :: mp_test
        double precision :: mu_obs, dmu_low_obs, dmu_up_obs, dmu_low0_obs, dmu_up0_obs, m_obs
        double precision, allocatable :: channel_w_model(:)
        double precision, allocatable :: channel_mu(:)
        double precision :: total_mu
        integer :: Higgscomb
        integer :: domH
        double precision :: chisq
!-n.b. these are the smeared observed signal strengths for this Higgs boson
    end type

    type mpred
        type(mp_neutHiggs), allocatable :: mp_Higgses(:)
        double precision :: mupred
    end type

    type observable
        integer :: id
        integer :: obstype
        type(mupeak) :: peak
        type(mutable) :: table
        type(neutHiggs), allocatable :: Higgses(:)
    end type
    type(observable), allocatable :: obs(:)

    type tablelist
        integer :: Npeaks
        integer :: id
        type(mutable) :: table
        type(mupeak), allocatable :: peaks(:)
        type(mpred) :: mpred
        type(neutHiggs), allocatable :: Higgses(:)
        ! This object holds primarily the Higgs boson predictions
        ! corresponding to tablelist%table. It corresponds to the full
        ! muplot if it is implemented (to enable the mpred-method).
    end type
    type(tablelist), allocatable :: analyses(:)

    type HSresults
        double precision :: Pvalue = -1.0D0
        double precision :: Pvalue_peak = -1.0D0
        double precision :: Pvalue_LHCRun1 = -1.0D0
        double precision :: Pvalue_STXS = -1.0D0
        double precision :: Chisq, Chisq_mu, Chisq_mh
        double precision :: Chisq_peak, Chisq_mpred, Chisq_peak_mu, Chisq_peak_mh
        double precision :: Chisq_LHCRun1, Chisq_LHCRun1_mu, Chisq_LHCRun1_mh, Chisq_LHCRun1_mhsep
        double precision :: Chisq_STXS, Chisq_STXS_rates, Chisq_STXS_mh
        double precision, allocatable :: mupred(:)
        integer, allocatable :: domH(:)
        integer, allocatable :: nH(:)
        integer, allocatable :: obsID(:)
        integer :: nobs, nobs_mu, nobs_mh
        integer :: nobs_peak, nobs_mpred, nobs_peak_mu, nobs_peak_mh, nanalysis
        integer :: nobs_LHCRun1, nobs_LHCRun1_mu, nobs_LHCRun1_mh
        integer :: nobs_STXS, nobs_STXS_rates, nobs_STXS_mh
        logical :: limitSetting = .false.
        double precision :: sHat, sUp, CLs
    end type
    type(HSresults), allocatable :: HSres(:)

!----------------------- Covariance matrices ----------------------
    double precision, allocatable :: cov(:, :)
    double precision, allocatable :: cov_mhneut(:, :, :)
    double precision, allocatable :: cov_mhneut_max(:, :, :)
    double precision, allocatable :: cov_mh_av(:, :)
    double precision, allocatable :: cov_mh_av_max(:, :)
    double precision, allocatable :: cov_mp(:, :)
    double precision, allocatable :: cov_mu_tot(:, :)
    double precision, allocatable :: mu_vector(:)
!--------------------------------------------------------------------
contains
!--------------------------------------------------------------------
    subroutine HiggsSignals_info
!--------------------------------------------------------------------
        use install_data_HS, only: version
        implicit none

        write (*, *)
        write (*, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write (*, *) "~                                                        ~"
        write (*, *) "~                  HiggsSignals "//adjustl(version)//"               ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~             Philip Bechtle, Sven Heinemeyer,           ~"
        write (*, *) "~              Tobias Klingl, Tim Stefaniak,             ~"
        write (*, *) "~             Georg Weiglein, Jonas Wittbrodt            ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~       arXiv:1305.1933 (Manual)                         ~"
        write (*, *) "~       arXiv:1403.1582 (application + more details)     ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~ It is based on the HiggsBounds-5 Fortran library.      ~"
        write (*, *) "~ Please consult and cite also the following references  ~"
        write (*, *) "~ for the HiggsBounds program                            ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~           arXiv:0811.4169, arXiv:1102.1898,            ~"
        write (*, *) "~           arXiv:1311.0055, arXiv:2006.06007            ~"
        write (*, *) "~      https://gitlab.com/higgsbounds/higgssignals       ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write (*, *)
        write (*, *) " HiggsSignals collects together results from "
        write (*, *)
        write (*, *) "    * the ATLAS and CMS Collaborations"
        write (*, *) "    * the CDF and D0 Collaborations"
        write (*, *) "    * the program HDECAY (arXiv:hep-ph/9704448)"
        write (*, *) "    * LHC Higgs Cross Section Working Group"
        write (*, *) "      (arXiv:1101.0593, arXiv:1201.3084, arXiv:1307.1347,"
        write (*, *) "       arXiv:1610.07922 and ref. therein)"
        write (*, *)

    end subroutine HiggsSignals_info
!--------------------------------------------------------------------
    subroutine print_dble_matrix(mat, title)
        !--------------------------------------------------------------------
        implicit none
        double precision, dimension(:, :), intent(in) :: mat(:, :)
        character(LEN=*), intent(in), optional :: title
        integer :: i

        if (present(title)) then
            write (*, *) "#*************************************************************************#"
            write (*, *) "# ", trim(title)
        endif
        write (*, *) "#*************************************************************************#"
        do i = lbound(mat, dim=1), ubound(mat, dim=1)
            write (*, *) mat(i, :)
        enddo
        write (*, *) "#*************************************************************************#"
    end subroutine print_dble_matrix
!--------------------------------------------------------------------
    subroutine deallocate_usefulbits_HS
!--------------------------------------------------------------------
        implicit none
        integer :: i
!  deallocate(neutHiggses)
        if (allocated(HSres)) then
            do i = lbound(HSres, dim=1), ubound(HSres, dim=1)
                if (allocated(HSres(i)%mupred)) deallocate (HSres(i)%mupred)
                if (allocated(HSres(i)%domH)) deallocate (HSres(i)%domH)
                if (allocated(HSres(i)%nH)) deallocate (HSres(i)%nH)
            enddo
            deallocate (HSres)
        endif

        if (allocated(corrlist)) deallocate (corrlist)

        call deallocate_covariance_matrices

    end subroutine deallocate_usefulbits_HS
!--------------------------------------------------------------------
    subroutine deallocate_covariance_matrices
!--------------------------------------------------------------------
        implicit none

        if (allocated(cov)) deallocate (cov)
        if (allocated(cov_mhneut)) deallocate (cov_mhneut)
        if (allocated(cov_mhneut_max)) deallocate (cov_mhneut_max)
        if (allocated(cov_mh_av)) deallocate (cov_mh_av)
        if (allocated(cov_mh_av_max)) deallocate (cov_mh_av_max)
        if (allocated(cov_mp)) deallocate (cov_mp)
        if (allocated(cov_mu_tot)) deallocate (cov_mu_tot)
        if (allocated(mu_vector)) deallocate (mu_vector)

    end subroutine deallocate_covariance_matrices
!--------------------------------------------------------------------
    function int_in_array(number, array)
        integer, intent(in) :: number
        integer, dimension(:), intent(in) :: array
        logical :: int_in_array

        integer :: i

        int_in_array = .False.

        do i = lbound(array, dim=1), ubound(array, dim=1)
            if (number .eq. array(i)) int_in_array = .True.
        enddo

    end function int_in_array
!------------------------------------------------------------------------------------

!--------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       StrCompress
!
! PURPOSE:
!       Subroutine to return a copy of an input string with all whitespace
!       (spaces and tabs) removed.
!
! CALLING SEQUENCE:
!       Result = StrCompress( String,  &  ! Input
!                             n = n    )  ! Optional Output
!
! INPUT ARGUMENTS:
!       String:         Character string to be compressed.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUT ARGUMENTS:
!       n:              Number of useful characters in output string
!                       after compression. From character n+1 -> LEN(Input_String)
!                       the output is padded with blanks.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT), OPTIONAL
!
! FUNCTION RESULT:
!       Result:         Input string with all whitespace removed before the
!                       first non-whitespace character, and from in-between
!                       non-whitespace characters.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(LEN(String))
!                       DIMENSION:  Scalar
!
! EXAMPLE:
!       Input_String = '  This is a string with spaces in it.'
!       Output_String = StrCompress( Input_String, n=n )
!       WRITE( *, '( a )' ) '>',Output_String( 1:n ),'<'
!   >Thisisastringwithspacesinit.<
!
!       or
!
!       WRITE( *, '( a )' ) '>',TRIM( Output_String ),'<'
!   >Thisisastringwithspacesinit.<
!
! PROCEDURE:
!       Definitions of a space and a tab character are made for the
!       ASCII collating sequence. Each single character of the input
!       string is checked against these definitions using the IACHAR()
!       intrinsic. If the input string character DOES NOT correspond
!       to a space or tab, it is not copied to the output string.
!
!       Note that for input that ONLY has spaces or tabs BEFORE the first
!       useful character, the output of this function is the same as the
!       ADJUSTL() instrinsic.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!
!:sdoc-:
!--------------------------------------------------------------------
    FUNCTION StrCompress(Input_String, n) RESULT(Output_String)
        ! Arguments
        CHARACTER(*), INTENT(IN)  :: Input_String
        INTEGER, OPTIONAL, INTENT(OUT) :: n
        ! Function result
        CHARACTER(LEN(Input_String)) :: Output_String
        ! Local parameters
        INTEGER, PARAMETER :: IACHAR_SPACE = 32
        INTEGER, PARAMETER :: IACHAR_TAB = 9
        ! Local variables
        INTEGER :: i, j
        INTEGER :: IACHAR_Character

        ! Setup
        ! -----
        ! Initialise output string
        Output_String = ' '
        ! Initialise output string "useful" length counter
        j = 0

        ! Loop over string contents character by character
        ! ------------------------------------------------
        DO i = 1, LEN(Input_String)

            ! Convert the current character to its position
            ! in the ASCII collating sequence
            IACHAR_Character = IACHAR(Input_String(i:i))

            ! If the character is NOT a space ' ' or a tab '->|'
            ! copy it to the output string.
            IF (IACHAR_Character /= IACHAR_SPACE .AND. &
                IACHAR_Character /= IACHAR_TAB) THEN
                j = j + 1
                Output_String(j:j) = Input_String(i:i)
            END IF

        END DO

        ! Save the non-whitespace count
        ! -----------------------------
        IF (PRESENT(n)) n = j

    END FUNCTION StrCompress

#if defined(__GFORTRAN__) && __GNUC__>=7
    !! Minimize the given chisq(s) between in `a<s<b` and return the
    !! parameters sHat and sigmaS that specify the corresponding
    !! parabola around the minimum.
    !!
    !! This only works if the compiler supports Fortran 2008 since
    !! we need to pass an internal function as argument.
    !!
    !! @param a lower bound on s
    !! @param b upper bound on s
    !! @param chisq a function of one parameter `s` which returns the chisq
    !! @param shat value of `s` that minimizes chisq
    !! @param sigmas width of the parabola around sHat in the Wald approximation and assuming s>0
    subroutine wald_approximate_chisq(a, b, chisq, sHat, sigmaS)

        interface
          function chisqFunc(s)
            double precision, intent(in) :: s
            double precision :: chisqFunc
          end function
        end interface
        double precision, intent(in) :: a, b
        procedure(chisqFunc) :: chisq
        double precision, intent(out) :: sHat, sigmaS

        double precision :: chisqMin, delta1

        call local_min(a, b, chisq, sHat, chisqMin)
        if (sHat < 0) then
            chisqMin = chisq(0D0)
        endif
        delta1 = zero(sHat, 2*b, deltaChisq1)
        sigmaS = sqrt(delta1 - sHat)
    contains
        function deltaChisq1(s)
            double precision, intent(in) :: s
            double precision deltaChisq1
            deltaChisq1 = chisq(s) - chisqMin - 1
            return
        end function

    end subroutine

#endif
!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    If the function F is defined on the interval (A,B), then local_min
!    finds an approximation X to the point at which F attatains its minimum
!    (or the appropriate limit point), and returns the value of F at X.
!
!    T and EPS define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is delta-unimodal for some delta less than TOL, the X approximates
!    the global minimum of F with an error less than 3*TOL.
!
!    If F is not delta-unimodal, then X may approximate a local, but
!    perhaps non-global, minimum.
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then, ignoring rounding errors, convergence is superlinear, and
!    usually of the order of about 1.3247.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    @param A lower bound of the interval
!    @param B upper bound of the interval
!    @param F the user-supplied function of the form `FUNCTION F ( X )` to minimise
!    @param xMin location of a local minimum in [A,B].
!    @param fxMin the value at the minimum F(X)
    subroutine local_min(a, b, f, xMin, fxMin)
        implicit none

        interface
            function func(x)
                double precision, intent(in) :: x
                double precision :: func
            end function
        end interface

        double precision, intent(in) :: a, b
        procedure(func) :: f
        double precision, intent(out) :: xMin, fxMin

        double precision, parameter :: vvsmall = 1.0D-20, vsmall = 1.0D-5, CGOLD = .3819660

        double precision d
        double precision e
        double precision fu
        double precision fv
        double precision fw
        double precision m
        double precision p
        double precision q
        double precision r
        double precision sa
        double precision sb
        double precision t2
        double precision tol
        double precision u
        double precision v
        double precision w

        sa = a
        sb = b
        xMin = sa + cgold*(b - a)
        w = xMin
        v = w
        e = 0.0D+00
        fxMin = f(xMin)
        fw = fxMin
        fv = fw

        do

            m = 0.5D+00*(sa + sb)
            tol = vvsmall*abs(xMin) + vsmall
            t2 = 2.0D+00*tol
            !
            !  Check the stopping criterion.
            !
            if (abs(xMin - m) <= t2 - 0.5D+00*(sb - sa)) then
                exit
            end if
            !
            !  Fit a parabola.
            !
            r = 0.0D+00
            q = r
            p = q

            if (tol < abs(e)) then

                r = (xMin - w)*(fxMin - fv)
                q = (xMin - v)*(fxMin - fw)
                p = (xMin - v)*q - (xMin - w)*r
                q = 2.0D+00*(q - r)

                if (0.0D+00 < q) then
                    p = -p
                end if

                q = abs(q)

                r = e
                e = d

            end if

            if (abs(p) < abs(0.5D+00*q*r) .and. &
                q*(sa - xMin) < p .and. &
                p < q*(sb - xMin)) then
                !
                !  Take the parabolic interpolation step.
                !
                d = p/q
                u = xMin + d
                !
                !  F must not be evaluated too close to A or B.
                !
                if ((u - sa) < t2 .or. (sb - u) < t2) then

                    if (xMin < m) then
                        d = tol
                    else
                        d = -tol
                    end if

                end if
                !
                !  A golden-section step.
                !
            else

                if (xMin < m) then
                    e = sb - xMin
                else
                    e = sa - xMin
                end if

                d = cgold*e

            end if
            !
            !  F must not be evaluated too close to X.
            !
            if (tol <= abs(d)) then
                u = xMin + d
            else if (0.0D+00 < d) then
                u = xMin + tol
            else
                u = xMin - tol
            end if

            fu = f(u)
            !
            !  Update A, B, V, W, and X.
            !
            if (fu <= fxMin) then

                if (u < xMin) then
                    sb = xMin
                else
                    sa = xMin
                end if

                v = w
                fv = fw
                w = xMin
                fw = fxMin
                xMin = u
                fxMin = fu

            else

                if (u < xMin) then
                    sa = u
                else
                    sb = u
                end if

                if (fu <= fw .or. abs(w - xMin) .lt. vvsmall) then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                else if (fu <= fv .or. abs(v - xMin) .lt. vvsmall .or. abs(w - v) .lt. vvsmall) then
                    v = u
                    fv = fu
                end if

            end if
        end do
    end subroutine

    function zero(a, b, f)
        !*****************************************************************************80
        !
        !! ZERO seeks the root of a function F(X) in an interval [A,B].
        !
        !  Discussion:
        !
        !    The interval [A,B] must be a change of sign interval for F.
        !    That is, F(A) and F(B) must be of opposite signs.  Then
        !    assuming that F is continuous implies the existence of at least
        !    one value C between A and B for which F(C) = 0.
        !
        !    The location of the zero is determined to within an accuracy
        !    of 6 * MACHEPS * abs ( C ) + 2 * T.
        !
        !    Thanks to Thomas Secretin for pointing out a transcription error in the
        !    setting of the value of P, 11 February 2013.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    11 February 2013
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Richard Brent.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Richard Brent,
        !    Algorithms for Minimization Without Derivatives,
        !    Dover, 2002,
        !    ISBN: 0-486-41998-3,
        !    LC: QA402.5.B74.
        !
        !  Parameters:
        !
        !    Input, double precision A, B, the endpoints of the change of
        !    sign interval.
        !
        !    Input, real ( kind = 8 ) T, a positive error tolerance.
        !
        !    Input, external real ( kind = 8 ) F, the name of a user-supplied
        !    function, of the form "FUNCTION F ( X )", which evaluates the
        !    function whose zero is being sought.
        !
        !    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
        !    the function F.
        !
        implicit none

        interface
            function func(x)
                double precision, intent(in) :: x
                double precision :: func
            end function
        end interface

        double precision, intent(in) :: a, b
        procedure(func) :: f
        double precision zero

        double precision c
        double precision d
        double precision e
        double precision fa
        double precision fb
        double precision fc
        double precision m
        double precision p
        double precision q
        double precision r
        double precision s
        double precision sa
        double precision sb
        double precision tol

        double precision, parameter :: machep = 1e-20, t = 1e-5

        !
        !  Make local copies of A and B.
        !
        sa = a
        sb = b
        fa = f(sa)
        fb = f(sb)

        c = sa
        fc = fa
        e = sb - sa
        d = e

        do

            if (abs(fc) < abs(fb)) then

                sa = sb
                sb = c
                c = sa
                fa = fb
                fb = fc
                fc = fa

            end if

            tol = 2.0D+00*machep*abs(sb) + t
            m = 0.5D+00*(c - sb)

            if (abs(m) <= tol .or. abs(fb) .lt. machep) then
                exit
            end if

            if (abs(e) < tol .or. abs(fa) <= abs(fb)) then

                e = m
                d = e

            else

                s = fb/fa

                if (abs(sa - c) .lt. machep) then

                    p = 2.0D+00*m*s
                    q = 1.0D+00 - s

                else

                    q = fa/fc
                    r = fb/fc
                    p = s*(2.0D+00*m*q*(q - r) - (sb - sa)*(r - 1.0D+00))
                    q = (q - 1.0D+00)*(r - 1.0D+00)*(s - 1.0D+00)

                end if

                if (0.0D+00 < p) then
                    q = -q
                else
                    p = -p
                end if

                s = e
                e = d

                if (2.0D+00*p < 3.0D+00*m*q - abs(tol*q) .and. &
                    p < abs(0.5D+00*s*q)) then
                    d = p/q
                else
                    e = m
                    d = e
                end if

            end if

            sa = sb
            fa = fb

            if (tol < abs(d)) then
                sb = sb + d
            else if (0.0D+00 < m) then
                sb = sb + tol
            else
                sb = sb - tol
            end if

            fb = f(sb)

            if ((0.0D+00 < fb .and. 0.0D+00 < fc) .or. &
                (fb <= 0.0D+00 .and. fc <= 0.0D+00)) then
                c = sa
                fc = fa
                e = sb - sa
                d = e
            end if

        end do

        zero = sb

        return
    end function
!--------------------------------------------------------------------
end module usefulbits_HS
!--------------------------------------------------------------------
