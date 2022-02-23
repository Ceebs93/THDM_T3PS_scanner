!******************************************************************
!> Contains global parameters, switches, objects.
module usefulbits
!******************************************************************
    implicit none
    logical :: debug = .False.

!> Switch for full mass uncertainty treatment.
!! If activated, HiggsBounds is re-run on dmhsteps number of mass variations
!! (for each Higgs boson), and takes the most conservative outcome.
    logical :: full_dmth_variation = .True.

!> Number of mass variations (within theoretical uncertainty),
!! if full_dmth_variation is set to True.
    integer :: dmhsteps = 3

!> Cut-off mass value. Lower mass uncertainties than this are not
!! considered fo variation.
    double precision :: small_mh = 0.1D0

!> Switches to the old HiggsBounds-3 run mode, where the global HiggsBounds result is taken only
!! from the most sensitive analysis among all Higgs bosons (and not for each Higgs boson individually).
    logical :: run_HB_classic = .False.

!> A file `Key.dat' that includes all applied analyses is generated if wantkey == True.
    logical :: wantkey = .True.

!> Allows to apply width-dependent limits even in case where the predicted width exceeds
!! the experimentally considered range.
    logical :: extrapolatewidth = .True.

!> If set to a non-zero value, experimental analyses are excluded from the
!! "standard" HiggsBounds run if there is a corresponding exclusion likelihood available.
!! This may be convenient in a parameter fit, where the likelihood is used instead.
    integer :: using_likelihood = 0

!> Number of HiggsBounds channel results (ranked by sensitivity) saved in output objects.
    integer, parameter :: numres = 3

!> Holds analyses-IDs those searches that should be considered in the standard HiggsBounds run.
!! Note: this is only used if whichanalyses=='list ', which can be set by
!! #initialize_higgsbounds.
    integer, allocatable :: analysislist(:)
!> Holds analyses-IDs for searches that should not be considered in the standard HiggsBounds run.
!! This can be specified by hand, or by using #channels::higgsbounds_deactivate_analyses().
    integer, allocatable :: analysis_exclude_list(:)

! For the LEP chisq extension:
    logical :: chisqcut_at_mumax = .False.
    logical :: BRdirectinput = .False.

! Global parameters (to be set during initialization)
    character(LEN=5) :: whichanalyses
    character(LEN=4) :: whichinput
    character(LEN=7) :: inputmethod = 'subrout'
    integer :: n_additional

! Internal parameters
    character(len=300) :: infile1, infile2
    integer, parameter :: file_id_common = 10
    integer, parameter :: file_id_common2 = 12
    integer, parameter :: file_id_common3 = 133
    integer, parameter :: file_id_common4 = 134
    integer, parameter :: file_id_debug1 = 444
    integer, parameter :: file_id_debug2 = 45

    !read from http://pdg.lbl.gov/ 22.10.2009
    double precision, parameter :: mt = 173.2D0
    double precision, parameter :: mc = 1.27D0
    !double precision, parameter :: ms = 0.105D0
    double precision, parameter :: mbmb = 4.20D0
    double precision, parameter :: mmu = 105.7D-3
    double precision, parameter :: mtau = 1.777D0
    !read in from https://pdg.lbl.gov/2019/tables/rpp2019-sum-quarks.pdf, 15.10.2020
    double precision, parameter :: mu = 2.16D-3
    double precision, parameter :: md = 4.67D-3
    double precision, parameter :: ms = 0.093D0
    double precision, parameter :: me = 0.511D-3

    ! double precision, parameter :: vev = 246.22D0
    double precision, parameter :: vev = 174.104D0
!-->
! CKM elements taken from PDG 2019 (fit):
    double precision, dimension(3, 3) :: CKM = reshape((/0.97446, 0.22438, 0.00896, &
                                                         0.22452, 0.97359, 0.04133, &
                                                         0.00365, 0.04214, 0.999105/), &
                                                       shape(CKM))
    ! Need to test correct index assignments!
    ! CKM(1,1) = 0.97446 ! V_ud
    ! CKM(1,2) = 0.22452 ! V_us
    ! CKM(1,3) = 0.00365 ! V_ub
    ! CKM(2,1) = 0.22438 ! V_cd
    ! CKM(2,2) = 0.97359 ! V_cs
    ! CKM(2,3) = 0.04214 ! V_cb
    ! CKM(3,1) = 0.00896 ! V_td
    ! CKM(3,2) = 0.04133 ! V_ts
    ! CKM(3,3) = 0.999105 ! V_tb
!<--

    double precision, parameter :: MZ = 91.1876D0 !PDG 2009
    double precision, parameter :: MW = 80.398D0 !PDG 2009
    double precision, parameter :: GF = 1.16637D-5
    double precision, parameter :: pi = 3.14159265358979323846264338328D0
    double precision, parameter :: alphas = 0.118D0

    double precision, parameter :: small = 1.0D-6
    double precision, parameter :: vsmall = 1.0D-16
    double precision, parameter :: vvsmall = 1.0D-100

    type particledescriptions
        character(LEN=10) :: short
        character(LEN=30) :: long
    end type

! particle codes: (n.b. these are NOT pdg)
    integer, parameter :: not_a_particle = 0
    integer, parameter :: Hneut = 1 !either Mhi, Mh2 or Mh3 (says nothing about CP properties)
    integer, parameter :: Hplus = 2 !single charged Higgs
    integer, parameter :: Chineut = 3 !either neutralino1, neutralino2, neutralino3 or neutralino4
    integer, parameter :: Chiplus = 4 !either chargino1 or chargino2
    integer, parameter :: anyH = 5 ! Considers both neutral and charged Higgs bosons
    integer :: np(0:5) = 1 !e.g np(Hneut) holds number of neutral Higgs considered
    type(particledescriptions), allocatable :: pdesc(:)

    ! HB-5.2: Needed for the channelrates_matrix
!  integer, parameter :: Nprod = 7
!  integer, parameter :: Ndecay = 9
!> Number of included Higgs boson production modes at hadron colliders
    integer, parameter :: Nprod = 11
!> Number of included Higgs boson decay modes
    integer, parameter :: Ndecay = 11

    !for subroutine version-------------------- (HB5: Removed!)
!  type inputsubroutineinfo
!   integer :: stat
!   character(LEN=40) :: desc
!   integer :: req
!  end type
!  type(inputsubroutineinfo),allocatable :: inputsub(:)

    logical :: just_after_run

    !associated with 'channels'----------------
    integer :: ntot

    type listprocesses
        integer :: tlist, ttype
        integer :: findi, findj
        integer :: corresponding_clsb_table_element
    end type

    type(listprocesses), allocatable :: pr(:)
    type(listprocesses), allocatable :: prsep(:, :)
    !-------------------------------------------

    !associated with 'input'--------------------

    type particlemasses
!> Central mass value (in GeV)
        double precision, allocatable  :: M(:)
!> Central mass value (used in mass variation procedure)
        double precision, allocatable  :: Mc(:)
!> Total decay width (in GeV)
        double precision, allocatable  :: GammaTot(:)
        double precision, allocatable  :: GammaTot_SM(:)
!> Mass uncertainties (chi-2 test) used in HiggsSignals
        double precision, allocatable  :: dM(:)
!> Mass uncertainties (variation) used in HiggsBounds
        double precision, allocatable  :: dMh(:)
    end type

    double precision, allocatable :: diffMhneut(:, :)
    double precision, allocatable :: diffMhch(:, :)
    double precision, allocatable :: dmn(:)
    double precision, allocatable :: dmch(:)

    integer ndmh

    integer ndat

    type lepdataset
        double precision, allocatable :: XS_hjZ_ratio(:)
        double precision, allocatable :: XS_bbhj_ratio(:)
        double precision, allocatable :: XS_tautauhj_ratio(:)
        double precision, allocatable :: XS_hjhi_ratio(:, :)
        double precision, allocatable :: XS_HpjHmj_ratio(:)
        double precision, allocatable :: XS_CpjCmj(:)
        double precision, allocatable :: XS_NjNi(:, :)
    end type

    type hadroncolliderdataset
        double precision, allocatable :: XS_hj_ratio(:)
        double precision, allocatable :: XS_gg_hj_ratio(:)  ! HB-5: for gluon fusion
        double precision, allocatable :: XS_bb_hj_ratio(:)  ! HB-5: for bb+Higgs production
        double precision, allocatable :: XS_hjZ_ratio(:)
        double precision, allocatable :: XS_gg_hjZ_ratio(:) ! HB-5 (TS 6.4.2018)
        double precision, allocatable :: XS_qq_hjZ_ratio(:) ! HB-5 (TS 6.4.2018)
        double precision, allocatable :: XS_hjW_ratio(:)
        double precision, allocatable :: XS_hjb_ratio(:)  ! still needed?
        double precision, allocatable :: XS_tthj_ratio(:)
        double precision, allocatable :: XS_vbf_ratio(:)
        double precision, allocatable :: XS_thj_tchan_ratio(:)  ! HB-5
        double precision, allocatable :: XS_thj_schan_ratio(:)  ! HB-5
        double precision, allocatable :: XS_tWhj_ratio(:)       ! HB-5
        double precision, allocatable :: XS_hjhi(:, :)           ! HB-5
! SM reference cross section holders:
        double precision, allocatable :: XS_HZ_SM(:)
        double precision, allocatable :: XS_gg_HZ_SM(:) ! HB-5 (TS 6.4.2018)
        double precision, allocatable :: XS_qq_HZ_SM(:) ! HB-5 (TS 6.4.2018)
        double precision, allocatable :: XS_HW_SM(:)
        double precision, allocatable :: XS_H_SM(:)
        double precision, allocatable :: XS_gg_H_SM(:)  ! HB-5
        double precision, allocatable :: XS_bb_H_SM(:)  ! HB-5
        !double precision, allocatable :: XS_H_SM_9713(:),XS_H_SM_9674(:)
        double precision, allocatable :: XS_ttH_SM(:)
        double precision, allocatable :: XS_tH_tchan_SM(:) ! HB-5
        double precision, allocatable :: XS_tH_schan_SM(:) ! HB-5
        double precision, allocatable :: XS_tWH_SM(:)  ! HB-5
        double precision, allocatable :: XS_vbf_SM(:)
        ! Higgs produced in association with b, where b is tagged, comes uncut and with various cuts
        ! see subroutines in theory_XS_SM_functions.f90 for details
        double precision, allocatable :: XS_Hb_SM(:)
        double precision, allocatable :: XS_Hb_c1_SM(:), XS_Hb_c2_SM(:), XS_Hb_c3_SM(:), XS_Hb_c4_SM(:)
        ! HB-5: Charged Higgs production cross sections (in pb)
        double precision, allocatable :: XS_vbf_Hpmj(:) ! for Hpm_j production in VBF
        double precision, allocatable :: XS_Hpmjtb(:)   ! for Hpm_j + t + b production
        double precision, allocatable :: XS_Hpmjcb(:)   ! for Hpm_j + c + b production
        double precision, allocatable :: XS_Hpmjbjet(:)   ! for Hpm_j + b + jet production
        double precision, allocatable :: XS_Hpmjcjet(:)   ! for Hpm_j + b + jet production
        double precision, allocatable :: XS_qq_Hpmj(:)   ! for Hpm_j + jet + jet production
        double precision, allocatable :: XS_HpmjW(:)   ! for Hpm_j + W production
        double precision, allocatable :: XS_HpmjZ(:)   ! for Hpm_j + Z production
        double precision, allocatable :: XS_HpjHmj(:) ! (j,i), for Hp_j Hm_j production
        double precision, allocatable :: XS_Hpmjhi(:, :) ! (j,i), for Hpm_j h_i production
! HB-5.2 beyond the narrow-width approximation matrix: holds the SM normalized channel rates
! with the dimensions (N_H, N_production-modes, N_decay-modes) = (N_H, 7, 9), where the
! ordering is the following
! 1: singleH, 2: VBF, 3: WH, 4: ZH, 5: ttH, 6: gg->phi, 7: bb->phi
! 1: gaga, 2: WW, 3: ZZ, 4: tautau, 5:bb, 6: Zga, 7: cc, 8: mumu, 9: gg
        double precision, allocatable :: channelrates(:, :, :)
! We need a temporary copy for the interface (will be copied in complete_theo)
        double precision, allocatable :: channelrates_tmp(:, :, :)
! This one holds the corresponding SM rates (in pb), assuming the NWA:
        double precision, allocatable :: channelrates_SM(:, :, :)
    end type

    type dataset
        logical :: gooddataset
        integer, allocatable :: CP_value(:)
        double precision, allocatable :: additional(:)
        type(particlemasses), allocatable :: particle(:)
        double precision, allocatable :: BR_hjss(:), BR_hjcc(:)
        double precision, allocatable :: BR_hjbb(:), BR_hjtt(:)
        double precision, allocatable :: BR_hjmumu(:), BR_hjtautau(:)
        double precision, allocatable :: BR_hjuu(:), BR_hjdd(:), BR_hjee(:) ! beyondHB extension
        double precision, allocatable :: BR_hjuc(:), BR_hjds(:) ! beyondHB extension
        double precision, allocatable :: BR_hjut(:), BR_hjdb(:) ! beyondHB extension
        double precision, allocatable :: BR_hjct(:), BR_hjsb(:) ! beyondHB extension
        double precision, allocatable :: BR_hjinvisible(:)
        logical :: full_BR_inv ! does BR_hjinvisible contain input or full BRinv values
        double precision, allocatable :: BR_hjhihi(:, :)     ! legacy HB-4
        double precision, allocatable :: BR_hkhjhi(:, :, :)   ! HB-5: for the decay h_k -> h_j h_i
        double precision, allocatable :: BR_hjhiZ(:, :)      ! HB-5: for the decay h_j -> h_i Z
        double precision, allocatable :: BR_hjemu(:), BR_hjetau(:), BR_hjmutau(:) ! HB-5
        double precision, allocatable :: BR_hjHpiW(:, :)   ! HB-5: for the decay h_j -> Hp_i W
        type(lepdataset) :: lep
        !-------------------------------------------
        double precision, allocatable :: BR_hjWW(:), BR_hjgaga(:)
        double precision, allocatable :: BR_hjZga(:)
        double precision, allocatable :: BR_hjZZ(:), BR_hjgg(:)

        double precision :: BR_tWpb
        double precision, allocatable :: BR_tHpjb(:)

        double precision, allocatable :: BR_Hpjcs(:)
        double precision, allocatable :: BR_Hpjcb(:)
        double precision, allocatable :: BR_Hpjtaunu(:)

! TS 2020-11-02: New BRs for (exotic) decays to light fermions:
        double precision, allocatable :: BR_Hpjud(:)
        double precision, allocatable :: BR_Hpjus(:)
        double precision, allocatable :: BR_Hpjcd(:)
        double precision, allocatable :: BR_Hpjub(:)
        double precision, allocatable :: BR_Hpjenu(:)
        double precision, allocatable :: BR_Hpjmunu(:)

        double precision, allocatable :: BR_Hpjtb(:)    ! HB-5: for the decay Hp_j -> t b
        double precision, allocatable :: BR_HpjWZ(:)   ! HB-5: for the decay Hp_j -> W Z
        double precision, allocatable :: BR_HpjhiW(:, :)   ! HB-5: for the decay Hp_j -> h_i W

        double precision, allocatable :: BR_CjqqNi(:, :)
        double precision, allocatable :: BR_CjlnuNi(:, :)
        double precision, allocatable :: BR_CjWNi(:, :)
        double precision, allocatable :: BR_NjqqNi(:, :)
        double precision, allocatable :: BR_NjZNi(:, :)

        type(hadroncolliderdataset) :: tev
        type(hadroncolliderdataset) :: lhc7
        type(hadroncolliderdataset) :: lhc8
        type(hadroncolliderdataset) :: lhc13 ! HB-5

! NEW(24/09/2014, TS):
!   double precision, allocatable :: gg_hj_ratio(:)
!   double precision, allocatable :: bb_hj_ratio(:)

        double precision, allocatable :: BR_Htt_SM(:), BR_Hbb_SM(:) !HB-5 new H->tt
        double precision, allocatable :: BR_Hcc_SM(:), BR_Hss_SM(:)
        double precision, allocatable :: BR_Hmumu_SM(:), BR_Htautau_SM(:)
        double precision, allocatable :: BR_HWW_SM(:), BR_HZZ_SM(:), BR_HZga_SM(:), BR_Hgaga_SM(:), BR_Hgg_SM(:)
        double precision, allocatable :: BR_Hjets_SM(:)
!         double precision, allocatable :: GammaTot_SM(:)
        !-------------------------------------------
    end type

    type(dataset), allocatable :: theo(:)

    type sqcouplratio
        double precision, allocatable :: hjss_s(:), hjss_p(:)
        double precision, allocatable :: hjcc_s(:), hjcc_p(:)
        double precision, allocatable :: hjbb_s(:), hjbb_p(:)
        double precision, allocatable :: hjtoptop_s(:), hjtoptop_p(:)  ! ToDo: Change name top -> t !
        double precision, allocatable :: hjmumu_s(:), hjmumu_p(:)
        double precision, allocatable :: hjtautau_s(:), hjtautau_p(:)

        double precision, allocatable :: hjWW(:), hjZZ(:)
        double precision, allocatable :: hjZga(:)
        double precision, allocatable :: hjgaga(:), hjgg(:), hjggZ(:)
        double precision, allocatable :: hjhiZ(:, :)
    end type

    type(sqcouplratio), allocatable :: g2(:)

    type couplratio
!-->
! TS 15-10-2020: Add first generation and FV coupling combinations
        double precision, allocatable :: hjee_s(:), hjee_p(:)
        double precision, allocatable :: hjuu_s(:), hjuu_p(:)
        double precision, allocatable :: hjdd_s(:), hjdd_p(:)
        double precision, allocatable :: hjuc_s(:), hjuc_p(:) ! n.b.: absolute (dimensionless) coupling
        double precision, allocatable :: hjds_s(:), hjds_p(:) ! n.b.: absolute (dimensionless) coupling
        double precision, allocatable :: hjut_s(:), hjut_p(:) ! n.b.: absolute (dimensionless) coupling
        double precision, allocatable :: hjdb_s(:), hjdb_p(:) ! n.b.: absolute (dimensionless) coupling
        double precision, allocatable :: hjct_s(:), hjct_p(:) ! n.b.: absolute (dimensionless) coupling
        double precision, allocatable :: hjsb_s(:), hjsb_p(:) ! n.b.: absolute (dimensionless) coupling
!<--
        double precision, allocatable :: hjcc_s(:), hjcc_p(:)
        double precision, allocatable :: hjss_s(:), hjss_p(:)
        double precision, allocatable :: hjtt_s(:), hjtt_p(:)
        double precision, allocatable :: hjbb_s(:), hjbb_p(:)
        double precision, allocatable :: hjmumu_s(:), hjmumu_p(:)
        double precision, allocatable :: hjtautau_s(:), hjtautau_p(:)

        double precision, allocatable :: hjWW(:), hjZZ(:)
        double precision, allocatable :: hjZga(:)
        double precision, allocatable :: hjgaga(:), hjgg(:) !,hjggZ(:)
        double precision, allocatable :: hjhiZ(:, :)
    end type

    type(couplratio), allocatable :: effC(:)
!-->
! TS 15-10-2020: Add charged Higgs effective couplings,
!                normalized to MFV 2HDM (see e.g. Eq.(8) in 1210.2465)
    type couplratio_Hc
        double precision, allocatable :: hcjud_L(:), hcjud_R(:)
        double precision, allocatable :: hcjcs_L(:), hcjcs_R(:)
        double precision, allocatable :: hcjtb_L(:), hcjtb_R(:)
        double precision, allocatable :: hcjus_L(:), hcjus_R(:)
        double precision, allocatable :: hcjub_L(:), hcjub_R(:)
        double precision, allocatable :: hcjcd_L(:), hcjcd_R(:)
        double precision, allocatable :: hcjcb_L(:), hcjcb_R(:)
        double precision, allocatable :: hcjtd_L(:), hcjtd_R(:)
        double precision, allocatable :: hcjts_L(:), hcjts_R(:)
    end type

    type(couplratio_Hc), allocatable :: effC_Hc(:)
!<--

    type hadroncolliderextras
        !nq_hjWp,nq_hjWm,nq_hj,nq_hjZ are set in allocate_hadroncolliderextras_parts below
        double precision, allocatable :: qq_hjWp(:, :)
        integer :: nq_hjWp!=2 i.e. (u dbar), (c sbar)  e.g. allocate(tR%qq_hjWp(tR%nq_hjWp,np(Hneut)))
        double precision, allocatable :: qq_hjWm(:, :)
        integer :: nq_hjWm!=2 i.e. (ubar d), (cbar s)

        double precision, allocatable :: gg_hj(:)
        double precision, allocatable :: qq_hj(:, :)
        integer :: nq_hj!=5 i.e.(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)

        double precision, allocatable :: gg_hjZ(:)
        double precision, allocatable :: qq_hjZ(:, :)
        integer :: nq_hjZ!=5 i.e.(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)

        double precision, allocatable :: bg_hjb(:)
    end type

    type(hadroncolliderextras), allocatable :: partR(:)
    !-------------------------------------------

    !associated with 'output'--------------------
    integer rep

    type results
        integer, allocatable :: chan(:)
        double precision, allocatable :: obsratio(:)
        double precision, allocatable :: predratio(:)
        double precision, allocatable :: sfactor(:)
        double precision, allocatable :: axis_i(:)
        double precision, allocatable :: axis_j(:)
        integer, allocatable :: allowed95(:)
        integer, allocatable :: ncombined(:)
        character(LEN=4), allocatable :: channelselection(:)
    end type

    type(results), allocatable :: res(:)

    !--new in HB-4:
    type fullresults
        integer :: chan = 0
        integer :: ncombined = 0
        integer :: allowed95 = 1
        double precision :: obsratio = 0.0D0
    end type

    type(fullresults), allocatable :: fullHBres(:)

    integer, allocatable :: allocate_if_stats_required(:)

! Needed to store relevant information on next-to-most sensitive channels:
    integer, allocatable ::  HBresult_all(:, :), chan_all(:, :), ncombined_all(:, :)
    double precision, allocatable :: obsratio_all(:, :), predratio_all(:, :)
    !-------------------------------------------

contains

    subroutine HiggsBounds_info
        use install_data, only: version
        implicit none

        write (*, *)
        write (*, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write (*, *) "~                                                        ~"
        write (*, *) "~                   HiggsBounds "//adjustl(version)//"               ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~             Philip Bechtle, Sven Heinemeyer,           ~"
        write (*, *) "~              Tobias Klingl, Tim Stefaniak,             ~"
        write (*, *) "~             Georg Weiglein, Jonas Wittbrodt            ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~            arXiv:0811.4169, arXiv:1102.1898,           ~"
        write (*, *) "~            arXiv:1301.2345, arXiv:1311.0055,           ~"
        write (*, *) "~            arXiv:1507.06706                            ~"
        write (*, *) "~       https://gitlab.com/higgsbounds/higgsbounds       ~"
        write (*, *) "~                                                        ~"
        write (*, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write (*, *)
        write (*, *) "HiggsBounds collects together results from "
        write (*, *)
        write (*, *) "    * the LEP collaborations and LEP Higgs Working Group"
        write (*, *) "    * the CDF and D0 Collaborations"
        write (*, *) "    * the ATLAS and CMS Collaborations"
        write (*, *) "    * the program HDECAY (arXiv:hep-ph/9704448)"
        write (*, *) "    * the program VH@NNLO"
        write (*, *) "      (arXiv:1210.5347,arXiv:1802.04817)"
        write (*, *) "    * TeV4LHC Higgs Working Group report"
        write (*, *) "      (see arXiv:hep-ph/0612172 and refs. therein)"
        write (*, *) "    * LHC Higgs Cross Section Working Group"
        write (*, *) "      (arXiv:1101.0593, arXiv:1201.3084, arXiv:1307.1347,"
        write (*, *) "       arXiv:1610.07922 and refs. therein, including the "
        write (*, *) "       gluon fusion N3LO prediction (arXiv:1602.00695).)"

    end subroutine HiggsBounds_info

    !**********************************************************
    function div(a, b, divlimit, div0res)
        !**********************************************************
        ! be careful about using this - not a mathematical limit
        double precision :: div
        !--------------------------------------input
        double precision :: a, b, divlimit, div0res
        !-----------------------------------internal
        double precision :: small1, small2
        !-------------------------------------------
        small1 = 1.0D-28
        small2 = 1.0D-20

        if (abs(b) .gt. small1) then
            div = a / b
        elseif (abs(a) .lt. small2) then
            div = divlimit
            if (div .lt. 0) stop 'error type divA (see function div in module usefulbits)'
        else
            div = div0res
            if (div .lt. 0) stop 'error type divB (see function div in module usefulbits)'
        endif

    end function

!--TESTING
    !**********************************************************
    subroutine iselementofarray(value, array, output)
        !**********************************************************
        implicit none
        !-------------------------------------input and output
        double precision, intent(in) :: value
        double precision, allocatable, dimension(:), intent(in) :: array
        integer, intent(out) :: output
        !---------------------------------------------internal
        integer :: i
        !-----------------------------------------------------
        output = -1

        if (allocated(array)) then
            do i = lbound(array, dim=1), ubound(array, dim=1)
                if (abs(value - array(i)) .le. vsmall) output = 1
            enddo
        else
            stop 'error: Passing an unallocated array to subroutine iselementofarray!'
        endif

    end subroutine iselementofarray
    !**********************************************************
    subroutine fill_pdesc
        !**********************************************************
        integer :: x

        if (ubound(np, dim=1) .ne. 5) stop 'error: have made a mistake in subroutine fill_pdesc (1)'

        x = 0
        allocate (pdesc(ubound(np, dim=1)))

        x = x + 1
        pdesc(x)%short = 'h'
        pdesc(x)%long = 'neutral Higgs boson'

        x = x + 1
        pdesc(x)%short = 'hplus'
        pdesc(x)%long = 'charged Higgs boson'

        x = x + 1
        pdesc(x)%short = 'N'
        pdesc(x)%long = 'neutralino'

        x = x + 1
        pdesc(x)%short = 'C'
        pdesc(x)%long = 'chargino'

        x = x + 1
        pdesc(x)%short = 'h/h+'
        pdesc(x)%long = 'neutral or charged scalar'

        if (x .ne. ubound(np, dim=1)) stop 'error: have made a mistake in subroutine fill_pdesc (2)'

    end subroutine fill_pdesc
    !**********************************************************
    subroutine allocate_dataset_parts(d, n_addit)
        !**********************************************************
        implicit none
        !-------------------------------------------
        type(dataset) :: d(:)
        !--------------------------------------input
        integer, intent(in) :: n_addit
        !-----------------------------------internal
        integer :: n_add, x, y
        integer, allocatable :: np_t(:)
        !-------------------------------------------
        allocate (np_t(lbound(np, dim=1):ubound(np, dim=1)))

        np_t = np
        do x = lbound(np_t, dim=1), ubound(np_t, dim=1)
            if (np(x) > 0) then
                np_t(x) = np(x)
            elseif (np(x) .eq. 0) then
                np_t(x) = 1
            else
                write (*, *) 'np=', np
                stop 'error in subroutine allocate_dataset_parts (1)'
            endif
        enddo

        if (n_addit > 0) then
            n_add = n_addit
        elseif (n_addit .eq. 0) then
            n_add = 1
        else
            stop 'error in subroutine allocate_dataset_parts (2)'
        endif

        do x = lbound(d, dim=1), ubound(d, dim=1)

            allocate (d(x)%additional(n_add))

            allocate (d(x)%particle(ubound(np_t, dim=1)))
            do y = 1, ubound(np_t, dim=1)
                allocate (d(x)%particle(y)%M(np_t(y)))
                allocate (d(x)%particle(y)%Mc(np_t(y)))
                allocate (d(x)%particle(y)%GammaTot(np_t(y)))
                allocate (d(x)%particle(y)%GammaTot_SM(np_t(y)))
                allocate (d(x)%particle(y)%dM(np_t(y)))
                allocate (d(x)%particle(y)%dMh(np_t(y)))
            enddo

            allocate (d(x)%lep%XS_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lep%XS_bbhj_ratio(np_t(Hneut)))
            allocate (d(x)%lep%XS_tautauhj_ratio(np_t(Hneut)))
            allocate (d(x)%lep%XS_hjhi_ratio(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%lep%XS_HpjHmj_ratio(np_t(Hplus)))
            allocate (d(x)%lep%XS_CpjCmj(np_t(Chiplus)))
            allocate (d(x)%lep%XS_NjNi(np_t(Chineut), np_t(Chineut)))

            allocate (d(x)%BR_hjuu(np_t(Hneut)))
            allocate (d(x)%BR_hjdd(np_t(Hneut)))
            allocate (d(x)%BR_hjss(np_t(Hneut)))
            allocate (d(x)%BR_hjcc(np_t(Hneut)))
            allocate (d(x)%BR_hjbb(np_t(Hneut)))
            allocate (d(x)%BR_hjtt(np_t(Hneut)))
            allocate (d(x)%BR_hjuc(np_t(Hneut)))
            allocate (d(x)%BR_hjds(np_t(Hneut)))
            allocate (d(x)%BR_hjut(np_t(Hneut)))
            allocate (d(x)%BR_hjdb(np_t(Hneut)))
            allocate (d(x)%BR_hjct(np_t(Hneut)))
            allocate (d(x)%BR_hjsb(np_t(Hneut)))
            allocate (d(x)%BR_hjee(np_t(Hneut)))
            allocate (d(x)%BR_hjmumu(np_t(Hneut)))
            allocate (d(x)%BR_hjtautau(np_t(Hneut)))

            allocate (d(x)%BR_hkhjhi(np_t(Hneut), np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%BR_hjhihi(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%BR_hjhiZ(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%BR_hjHpiW(np_t(Hneut), np_t(Hplus)))
            allocate (d(x)%BR_hjWW(np_t(Hneut)))
            allocate (d(x)%BR_hjZZ(np_t(Hneut)))
            allocate (d(x)%BR_hjZga(np_t(Hneut)))
            allocate (d(x)%BR_hjgaga(np_t(Hneut)))
            allocate (d(x)%BR_hjgg(np_t(Hneut)))
            d(x)%full_BR_inv = .false.
            allocate (d(x)%BR_hjinvisible(np_t(Hneut)))
            allocate (d(x)%BR_hjemu(np_t(Hneut)))
            allocate (d(x)%BR_hjetau(np_t(Hneut)))
            allocate (d(x)%BR_hjmutau(np_t(Hneut)))

            allocate (d(x)%BR_tHpjb(np_t(Hplus)))
            allocate (d(x)%BR_Hpjcs(np_t(Hplus)))
            allocate (d(x)%BR_Hpjcb(np_t(Hplus)))
            allocate (d(x)%BR_Hpjtaunu(np_t(Hplus)))
            allocate (d(x)%BR_Hpjud(np_t(Hplus)))
            allocate (d(x)%BR_Hpjus(np_t(Hplus)))
            allocate (d(x)%BR_Hpjcd(np_t(Hplus)))
            allocate (d(x)%BR_Hpjub(np_t(Hplus)))
            allocate (d(x)%BR_Hpjenu(np_t(Hplus)))
            allocate (d(x)%BR_Hpjmunu(np_t(Hplus)))
            allocate (d(x)%BR_Hpjtb(np_t(Hplus)))
            allocate (d(x)%BR_HpjWZ(np_t(Hplus)))
            allocate (d(x)%BR_HpjhiW(np_t(Hplus), np_t(Hneut)))

            allocate (d(x)%BR_CjqqNi(np_t(Chiplus), np_t(Chineut)))
            allocate (d(x)%BR_CjlnuNi(np_t(Chiplus), np_t(Chineut)))
            allocate (d(x)%BR_CjWNi(np_t(Chiplus), np_t(Chineut)))
            allocate (d(x)%BR_NjqqNi(np_t(Chineut), np_t(Chineut)))
            allocate (d(x)%BR_NjZNi(np_t(Chineut), np_t(Chineut)))

            allocate (d(x)%tev%XS_hjb_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_tthj_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_vbf_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_gg_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_qq_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_hjW_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_hj_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_gg_hj_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_bb_hj_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_thj_tchan_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_thj_schan_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_tWhj_ratio(np_t(Hneut)))
            allocate (d(x)%tev%XS_hjhi(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%tev%XS_vbf_Hpmj(np_t(Hplus)))
            allocate (d(x)%tev%XS_Hpmjtb(np_t(Hplus)))
            allocate (d(x)%tev%XS_Hpmjcb(np_t(Hplus)))
            allocate (d(x)%tev%XS_Hpmjbjet(np_t(Hplus)))
            allocate (d(x)%tev%XS_Hpmjcjet(np_t(Hplus)))
            allocate (d(x)%tev%XS_qq_Hpmj(np_t(Hplus)))
            allocate (d(x)%tev%XS_HpmjW(np_t(Hplus)))
            allocate (d(x)%tev%XS_HpmjZ(np_t(Hplus)))
            allocate (d(x)%tev%XS_HpjHmj(np_t(Hplus)))
            allocate (d(x)%tev%XS_Hpmjhi(np_t(Hplus), np_t(Hneut)))
            allocate (d(x)%tev%channelrates(np_t(Hneut), Nprod, Ndecay))
            allocate (d(x)%tev%channelrates_tmp(np_t(Hneut), Nprod, Ndecay))

            allocate (d(x)%lhc7%XS_hjb_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_tthj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_vbf_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_gg_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_qq_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_hjW_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_gg_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_bb_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_thj_tchan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_thj_schan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_tWhj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc7%XS_hjhi(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%lhc7%XS_vbf_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_Hpmjtb(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_Hpmjcb(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_Hpmjbjet(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_Hpmjcjet(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_qq_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_HpmjW(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_HpmjZ(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_HpjHmj(np_t(Hplus)))
            allocate (d(x)%lhc7%XS_Hpmjhi(np_t(Hplus), np_t(Hneut)))
            allocate (d(x)%lhc7%channelrates(np_t(Hneut), Nprod, Ndecay))
            allocate (d(x)%lhc7%channelrates_tmp(np_t(Hneut), Nprod, Ndecay))

            allocate (d(x)%lhc8%XS_hjb_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_tthj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_vbf_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_gg_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_qq_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_hjW_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_gg_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_bb_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_thj_tchan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_thj_schan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_tWhj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc8%XS_hjhi(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%lhc8%XS_vbf_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_Hpmjtb(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_Hpmjcb(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_Hpmjbjet(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_Hpmjcjet(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_qq_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_HpmjW(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_HpmjZ(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_HpjHmj(np_t(Hplus)))
            allocate (d(x)%lhc8%XS_Hpmjhi(np_t(Hplus), np_t(Hneut)))
            allocate (d(x)%lhc8%channelrates(np_t(Hneut), Nprod, Ndecay))
            allocate (d(x)%lhc8%channelrates_tmp(np_t(Hneut), Nprod, Ndecay))

            allocate (d(x)%lhc13%XS_hjb_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_tthj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_vbf_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_gg_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_qq_hjZ_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_hjW_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_gg_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_bb_hj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_thj_tchan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_thj_schan_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_tWhj_ratio(np_t(Hneut)))
            allocate (d(x)%lhc13%XS_hjhi(np_t(Hneut), np_t(Hneut)))
            allocate (d(x)%lhc13%XS_vbf_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_Hpmjtb(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_Hpmjcb(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_Hpmjbjet(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_Hpmjcjet(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_qq_Hpmj(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_HpmjW(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_HpmjZ(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_HpjHmj(np_t(Hplus)))
            allocate (d(x)%lhc13%XS_Hpmjhi(np_t(Hplus), np_t(Hneut)))
            allocate (d(x)%lhc13%channelrates(np_t(Hneut), Nprod, Ndecay))
            allocate (d(x)%lhc13%channelrates_tmp(np_t(Hneut), Nprod, Ndecay))

            allocate (d(x)%CP_value(np_t(Hneut)))

            do y = 1, ubound(np_t, dim=1)
                d(x)%particle(y)%M = -1.0D0
                d(x)%particle(y)%Mc = -1.0D0
                d(x)%particle(y)%GammaTot = 0.0D0
                d(x)%particle(y)%GammaTot_SM = 0.0D0
                d(x)%particle(y)%dM = 0.0D0
                d(x)%particle(y)%dMh = 0.0D0
            enddo

            d(x)%lep%XS_hjZ_ratio = 0.0D0
            d(x)%lep%XS_bbhj_ratio = 0.0D0
            d(x)%lep%XS_tautauhj_ratio = 0.0D0
            d(x)%lep%XS_hjhi_ratio = 0.0D0
            d(x)%lep%XS_HpjHmj_ratio = 0.0D0
            d(x)%lep%XS_CpjCmj = 0.0D0
            d(x)%lep%XS_NjNi = 0.0D0

            d(x)%BR_hjuu = 0.0D0
            d(x)%BR_hjdd = 0.0D0
            d(x)%BR_hjss = 0.0D0
            d(x)%BR_hjcc = 0.0D0
            d(x)%BR_hjbb = 0.0D0
            d(x)%BR_hjtt = 0.0D0
            d(x)%BR_hjuc = 0.0D0
            d(x)%BR_hjds = 0.0D0
            d(x)%BR_hjut = 0.0D0
            d(x)%BR_hjdb = 0.0D0
            d(x)%BR_hjct = 0.0D0
            d(x)%BR_hjsb = 0.0D0
            d(x)%BR_hjee = 0.0D0
            d(x)%BR_hjmumu = 0.0D0
            d(x)%BR_hjtautau = 0.0D0
            d(x)%BR_hjWW = 0.0D0
            d(x)%BR_hjZZ = 0.0D0
            d(x)%BR_hjZga = 0.0D0
            d(x)%BR_hjgaga = 0.0D0
            d(x)%BR_hjgg = 0.0D0
            d(x)%BR_hjinvisible = 0.0D0
            d(x)%BR_hjhihi = 0.0D0
            d(x)%BR_hjhiZ = 0.0D0
            d(x)%BR_hkhjhi = 0.0D0
            d(x)%BR_hjHpiW = 0.0D0
            d(x)%BR_hjemu = 0.0D0
            d(x)%BR_hjetau = 0.0D0
            d(x)%BR_hjmutau = 0.0D0

            d(x)%BR_tWpb = 0.0D0
            d(x)%BR_tHpjb = 0.0D0
            d(x)%BR_Hpjcs = 0.0D0
            d(x)%BR_Hpjcb = 0.0D0
            d(x)%BR_Hpjtaunu = 0.0D0
            d(x)%BR_Hpjud = 0.0D0
            d(x)%BR_Hpjus = 0.0D0
            d(x)%BR_Hpjcd = 0.0D0
            d(x)%BR_Hpjub = 0.0D0
            d(x)%BR_Hpjenu = 0.0D0
            d(x)%BR_Hpjmunu = 0.0D0

            d(x)%BR_Hpjtb = 0.0D0
            d(x)%BR_HpjWZ = 0.0D0
            d(x)%BR_HpjhiW = 0.0D0

            d(x)%BR_CjqqNi = 0.0D0
            d(x)%BR_CjlnuNi = 0.0D0
            d(x)%BR_CjWNi = 0.0D0
            d(x)%BR_NjqqNi = 0.0D0
            d(x)%BR_NjZNi = 0.0D0

            d(x)%tev%XS_hjb_ratio = 0.0D0
            d(x)%tev%XS_tthj_ratio = 0.0D0
            d(x)%tev%XS_vbf_ratio = 0.0D0
            d(x)%tev%XS_hj_ratio = 0.0D0
            d(x)%tev%XS_hjW_ratio = 0.0D0
            d(x)%tev%XS_hjZ_ratio = 0.0D0
            d(x)%tev%XS_gg_hj_ratio = 0.0D0
            d(x)%tev%XS_bb_hj_ratio = 0.0D0
            d(x)%tev%XS_thj_tchan_ratio = 0.0D0
            d(x)%tev%XS_thj_schan_ratio = 0.0D0
            d(x)%tev%XS_tWhj_ratio = 0.0D0
            d(x)%tev%XS_hjhi = 0.0D0
            d(x)%tev%XS_vbf_Hpmj = 0.0D0
            d(x)%tev%XS_Hpmjtb = 0.0D0
            d(x)%tev%XS_Hpmjcb = 0.0D0
            d(x)%tev%XS_Hpmjbjet = 0.0D0
            d(x)%tev%XS_Hpmjcjet = 0.0D0
            d(x)%tev%XS_qq_Hpmj = 0.0D0
            d(x)%tev%XS_HpmjW = 0.0D0
            d(x)%tev%XS_HpmjZ = 0.0D0
            d(x)%tev%XS_HpjHmj = 0.0D0
            d(x)%tev%XS_Hpmjhi = 0.0D0
            d(x)%tev%channelrates = 0.0D0
            d(x)%tev%channelrates_tmp = -1.0D0

            d(x)%lhc7%XS_hjb_ratio = 0.0D0
            d(x)%lhc7%XS_tthj_ratio = 0.0D0
            d(x)%lhc7%XS_vbf_ratio = 0.0D0
            d(x)%lhc7%XS_hj_ratio = 0.0D0
            d(x)%lhc7%XS_hjW_ratio = 0.0D0
            d(x)%lhc7%XS_hjZ_ratio = 0.0D0
            d(x)%lhc7%XS_gg_hj_ratio = 0.0D0
            d(x)%lhc7%XS_bb_hj_ratio = 0.0D0
            d(x)%lhc7%XS_thj_tchan_ratio = 0.0D0
            d(x)%lhc7%XS_thj_schan_ratio = 0.0D0
            d(x)%lhc7%XS_tWhj_ratio = 0.0D0
            d(x)%lhc7%XS_hjhi = 0.0D0
            d(x)%lhc7%XS_vbf_Hpmj = 0.0D0
            d(x)%lhc7%XS_Hpmjtb = 0.0D0
            d(x)%lhc7%XS_Hpmjcb = 0.0D0
            d(x)%lhc7%XS_Hpmjbjet = 0.0D0
            d(x)%lhc7%XS_Hpmjcjet = 0.0D0
            d(x)%lhc7%XS_qq_Hpmj = 0.0D0
            d(x)%lhc7%XS_HpmjW = 0.0D0
            d(x)%lhc7%XS_HpmjZ = 0.0D0
            d(x)%lhc7%XS_HpjHmj = 0.0D0
            d(x)%lhc7%XS_Hpmjhi = 0.0D0
            d(x)%lhc7%channelrates = 0.0D0
            d(x)%lhc7%channelrates_tmp = -1.0D0

            d(x)%lhc8%XS_hjb_ratio = 0.0D0
            d(x)%lhc8%XS_tthj_ratio = 0.0D0
            d(x)%lhc8%XS_vbf_ratio = 0.0D0
            d(x)%lhc8%XS_hj_ratio = 0.0D0
            d(x)%lhc8%XS_hjW_ratio = 0.0D0
            d(x)%lhc8%XS_hjZ_ratio = 0.0D0
            d(x)%lhc8%XS_gg_hj_ratio = 0.0D0
            d(x)%lhc8%XS_bb_hj_ratio = 0.0D0
            d(x)%lhc8%XS_thj_tchan_ratio = 0.0D0
            d(x)%lhc8%XS_thj_schan_ratio = 0.0D0
            d(x)%lhc8%XS_tWhj_ratio = 0.0D0
            d(x)%lhc8%XS_hjhi = 0.0D0
            d(x)%lhc8%XS_vbf_Hpmj = 0.0D0
            d(x)%lhc8%XS_Hpmjtb = 0.0D0
            d(x)%lhc8%XS_Hpmjcb = 0.0D0
            d(x)%lhc8%XS_Hpmjbjet = 0.0D0
            d(x)%lhc8%XS_Hpmjcjet = 0.0D0
            d(x)%lhc8%XS_qq_Hpmj = 0.0D0
            d(x)%lhc8%XS_HpmjW = 0.0D0
            d(x)%lhc8%XS_HpmjZ = 0.0D0
            d(x)%lhc8%XS_HpjHmj = 0.0D0
            d(x)%lhc8%XS_Hpmjhi = 0.0D0
            d(x)%lhc8%channelrates = 0.0D0
            d(x)%lhc8%channelrates_tmp = -1.0D0

            d(x)%lhc13%XS_hjb_ratio = 0.0D0
            d(x)%lhc13%XS_tthj_ratio = 0.0D0
            d(x)%lhc13%XS_vbf_ratio = 0.0D0
            d(x)%lhc13%XS_hj_ratio = 0.0D0
            d(x)%lhc13%XS_hjW_ratio = 0.0D0
            d(x)%lhc13%XS_hjZ_ratio = 0.0D0
            d(x)%lhc13%XS_gg_hj_ratio = 0.0D0
            d(x)%lhc13%XS_bb_hj_ratio = 0.0D0
            d(x)%lhc13%XS_thj_tchan_ratio = 0.0D0
            d(x)%lhc13%XS_thj_schan_ratio = 0.0D0
            d(x)%lhc13%XS_tWhj_ratio = 0.0D0
            d(x)%lhc13%XS_hjhi = 0.0D0
            d(x)%lhc13%XS_vbf_Hpmj = 0.0D0
            d(x)%lhc13%XS_Hpmjtb = 0.0D0
            d(x)%lhc13%XS_Hpmjcb = 0.0D0
            d(x)%lhc13%XS_Hpmjbjet = 0.0D0
            d(x)%lhc13%XS_Hpmjcjet = 0.0D0
            d(x)%lhc13%XS_qq_Hpmj = 0.0D0
            d(x)%lhc13%XS_HpmjW = 0.0D0
            d(x)%lhc13%XS_HpmjZ = 0.0D0
            d(x)%lhc13%XS_HpjHmj = 0.0D0
            d(x)%lhc13%XS_Hpmjhi = 0.0D0
            d(x)%lhc13%channelrates = 0.0D0
            d(x)%lhc13%channelrates_tmp = -1.0D0

            d(x)%additional = 0.0D0

            d(x)%CP_value = 0
        enddo

        select case (whichanalyses)
        case ('onlyH', 'LandH', 'onlyP', 'list ')
            do x = lbound(d, dim=1), ubound(d, dim=1)
                allocate (d(x)%tev%XS_HZ_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_gg_HZ_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_qq_HZ_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_HW_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_H_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_ttH_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_vbf_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_gg_H_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_bb_H_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_tH_tchan_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_tH_schan_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_tWH_SM(np_t(Hneut)))
                allocate (d(x)%tev%channelrates_SM(np_t(Hneut), Nprod, Ndecay))

                allocate (d(x)%tev%XS_Hb_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_Hb_c1_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_Hb_c2_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_Hb_c3_SM(np_t(Hneut)))
                allocate (d(x)%tev%XS_Hb_c4_SM(np_t(Hneut)))

                allocate (d(x)%lhc7%XS_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_gg_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_qq_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_HW_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_ttH_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_vbf_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_gg_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_bb_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_tH_tchan_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_tH_schan_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%XS_tWH_SM(np_t(Hneut)))
                allocate (d(x)%lhc7%channelrates_SM(np_t(Hneut), Nprod, Ndecay))

                allocate (d(x)%lhc7%XS_Hb_SM(np_t(Hneut)))
!     allocate(d(x)%lhc7%XS_Hb_c1_SM(      np_t(Hneut)    ))
!     allocate(d(x)%lhc7%XS_Hb_c2_SM(      np_t(Hneut)    ))
!     allocate(d(x)%lhc7%XS_Hb_c3_SM(      np_t(Hneut)    ))

                allocate (d(x)%lhc8%XS_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_gg_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_qq_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_HW_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_ttH_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_vbf_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_gg_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_bb_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_tH_tchan_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_tH_schan_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%XS_tWH_SM(np_t(Hneut)))
                allocate (d(x)%lhc8%channelrates_SM(np_t(Hneut), Nprod, Ndecay))

                allocate (d(x)%lhc8%XS_Hb_SM(np_t(Hneut)))
!     allocate(d(x)%lhc8%XS_Hb_c1_SM(      np_t(Hneut)    ))
!     allocate(d(x)%lhc8%XS_Hb_c2_SM(      np_t(Hneut)    ))
!     allocate(d(x)%lhc8%XS_Hb_c3_SM(      np_t(Hneut)    ))

                allocate (d(x)%lhc13%XS_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_gg_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_qq_HZ_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_HW_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_ttH_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_vbf_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_gg_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_bb_H_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_tH_tchan_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_tH_schan_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%XS_tWH_SM(np_t(Hneut)))
                allocate (d(x)%lhc13%channelrates_SM(np_t(Hneut), Nprod, Ndecay))
!     allocate(d(x)%lhc8%XS_Hb_SM(         np_t(Hneut)    ))

                allocate (d(x)%BR_Hbb_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hcc_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hss_SM(np_t(Hneut)))
                allocate (d(x)%BR_Htt_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hmumu_SM(np_t(Hneut)))
                allocate (d(x)%BR_Htautau_SM(np_t(Hneut)))
                allocate (d(x)%BR_HWW_SM(np_t(Hneut)))
                allocate (d(x)%BR_HZZ_SM(np_t(Hneut)))
                allocate (d(x)%BR_HZga_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hgaga_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hgg_SM(np_t(Hneut)))
                allocate (d(x)%BR_Hjets_SM(np_t(Hneut)))
!                 allocate (d(x)%GammaTot_SM(np_t(Hneut)))
            enddo
        case ('onlyL')
        case default
            stop 'error in allocate_dataset_parts (3)'
        end select
        deallocate (np_t)
    end subroutine allocate_dataset_parts

    !**********************************************************
    subroutine allocate_sqcouplratio_parts(gsq)
        ! to use this, gsq must be an array
        !**********************************************************
        implicit none
        !-------------------------------------------
        type(sqcouplratio) :: gsq(:)
        !-----------------------------------internal
        integer :: x
        integer :: nHiggsneut
        !-------------------------------------------

        if (np(Hneut) > 0) then
            nHiggsneut = np(Hneut)
        elseif (np(Hneut) .eq. 0) then
            nHiggsneut = 1
        else
            stop 'error in subroutine allocate_sqcouplratio_parts (1)'
        endif

        do x = lbound(gsq, dim=1), ubound(gsq, dim=1)
            allocate (gsq(x)%hjss_s(nHiggsneut), gsq(x)%hjss_p(nHiggsneut))
            allocate (gsq(x)%hjcc_s(nHiggsneut), gsq(x)%hjcc_p(nHiggsneut))
            allocate (gsq(x)%hjbb_s(nHiggsneut), gsq(x)%hjbb_p(nHiggsneut))
            allocate (gsq(x)%hjtoptop_s(nHiggsneut), gsq(x)%hjtoptop_p(nHiggsneut))
            allocate (gsq(x)%hjmumu_s(nHiggsneut), gsq(x)%hjmumu_p(nHiggsneut))
            allocate (gsq(x)%hjtautau_s(nHiggsneut), gsq(x)%hjtautau_p(nHiggsneut))

            allocate (gsq(x)%hjWW(nHiggsneut), gsq(x)%hjZZ(nHiggsneut))
            allocate (gsq(x)%hjZga(nHiggsneut))
            allocate (gsq(x)%hjgaga(nHiggsneut), gsq(x)%hjgg(nHiggsneut))
            allocate (gsq(x)%hjggZ(nHiggsneut))
            allocate (gsq(x)%hjhiZ(nHiggsneut, nHiggsneut))

            gsq(x)%hjss_s = 0.0D0
            gsq(x)%hjss_p = 0.0D0
            gsq(x)%hjcc_s = 0.0D0
            gsq(x)%hjcc_p = 0.0D0
            gsq(x)%hjbb_s = 0.0D0
            gsq(x)%hjbb_p = 0.0D0
            gsq(x)%hjtoptop_s = 0.0D0
            gsq(x)%hjtoptop_p = 0.0D0
            gsq(x)%hjmumu_s = 0.0D0
            gsq(x)%hjmumu_p = 0.0D0
            gsq(x)%hjtautau_s = 0.0D0
            gsq(x)%hjtautau_p = 0.0D0

            gsq(x)%hjWW = 0.0D0
            gsq(x)%hjZZ = 0.0D0
            gsq(x)%hjZga = 0.0D0
            gsq(x)%hjgaga = 0.0D0
            gsq(x)%hjgg = 0.0D0
            gsq(x)%hjggZ = 0.0D0
            gsq(x)%hjhiZ = 0.0D0
        enddo

    end subroutine allocate_sqcouplratio_parts
    !**********************************************************
    subroutine allocate_couplratio_parts(g)
        ! to use this, gsq must be an array
        !**********************************************************
        implicit none
        !-------------------------------------------
        type(couplratio) :: g(:)
        !-----------------------------------internal
        integer :: x
        integer :: nHiggsneut
        !-------------------------------------------

        if (np(Hneut) > 0) then
            nHiggsneut = np(Hneut)
        elseif (np(Hneut) .eq. 0) then
            nHiggsneut = 1
        else
            stop 'error in subroutine allocate_couplratio_parts (1)'
        endif

        do x = lbound(g, dim=1), ubound(g, dim=1)

            allocate (g(x)%hjee_s(nHiggsneut), g(x)%hjee_p(nHiggsneut))
            allocate (g(x)%hjuu_s(nHiggsneut), g(x)%hjuu_p(nHiggsneut))
            allocate (g(x)%hjdd_s(nHiggsneut), g(x)%hjdd_p(nHiggsneut))
            allocate (g(x)%hjuc_s(nHiggsneut), g(x)%hjuc_p(nHiggsneut))
            allocate (g(x)%hjds_s(nHiggsneut), g(x)%hjds_p(nHiggsneut))
            allocate (g(x)%hjut_s(nHiggsneut), g(x)%hjut_p(nHiggsneut))
            allocate (g(x)%hjdb_s(nHiggsneut), g(x)%hjdb_p(nHiggsneut))
            allocate (g(x)%hjct_s(nHiggsneut), g(x)%hjct_p(nHiggsneut))
            allocate (g(x)%hjsb_s(nHiggsneut), g(x)%hjsb_p(nHiggsneut))

            allocate (g(x)%hjss_s(nHiggsneut), g(x)%hjss_p(nHiggsneut))
            allocate (g(x)%hjcc_s(nHiggsneut), g(x)%hjcc_p(nHiggsneut))
            allocate (g(x)%hjbb_s(nHiggsneut), g(x)%hjbb_p(nHiggsneut))
            allocate (g(x)%hjtt_s(nHiggsneut), g(x)%hjtt_p(nHiggsneut))
            allocate (g(x)%hjmumu_s(nHiggsneut), g(x)%hjmumu_p(nHiggsneut))
            allocate (g(x)%hjtautau_s(nHiggsneut), g(x)%hjtautau_p(nHiggsneut))

            allocate (g(x)%hjWW(nHiggsneut), g(x)%hjZZ(nHiggsneut))
            allocate (g(x)%hjZga(nHiggsneut))
            allocate (g(x)%hjgaga(nHiggsneut), g(x)%hjgg(nHiggsneut))
!    allocate(g(x)%hjggZ(nHiggsneut)    )
            allocate (g(x)%hjhiZ(nHiggsneut, nHiggsneut))

            g(x)%hjee_s = 0.0D0
            g(x)%hjee_p = 0.0D0
            g(x)%hjuu_s = 0.0D0
            g(x)%hjuu_p = 0.0D0
            g(x)%hjdd_s = 0.0D0
            g(x)%hjdd_p = 0.0D0
            g(x)%hjuc_s = 0.0D0
            g(x)%hjuc_p = 0.0D0
            g(x)%hjut_s = 0.0D0
            g(x)%hjut_p = 0.0D0
            g(x)%hjds_s = 0.0D0
            g(x)%hjds_p = 0.0D0
            g(x)%hjdb_s = 0.0D0
            g(x)%hjdb_p = 0.0D0
            g(x)%hjct_s = 0.0D0
            g(x)%hjct_p = 0.0D0
            g(x)%hjsb_s = 0.0D0
            g(x)%hjsb_p = 0.0D0

            g(x)%hjss_s = 0.0D0
            g(x)%hjss_p = 0.0D0
            g(x)%hjcc_s = 0.0D0
            g(x)%hjcc_p = 0.0D0
            g(x)%hjbb_s = 0.0D0
            g(x)%hjbb_p = 0.0D0
            g(x)%hjtt_s = 0.0D0
            g(x)%hjtt_p = 0.0D0
            g(x)%hjmumu_s = 0.0D0
            g(x)%hjmumu_p = 0.0D0
            g(x)%hjtautau_s = 0.0D0
            g(x)%hjtautau_p = 0.0D0

            g(x)%hjWW = 0.0D0
            g(x)%hjZZ = 0.0D0
            g(x)%hjZga = 0.0D0
            g(x)%hjgaga = 0.0D0
            g(x)%hjgg = 0.0D0
!    g(x)%hjggZ          =0.0D0
            g(x)%hjhiZ = 0.0D0
        enddo

    end subroutine allocate_couplratio_parts

    !**********************************************************
    subroutine allocate_couplratio_Hc_parts(g)
        !**********************************************************
        implicit none
        !-------------------------------------------
        type(couplratio_Hc) :: g(:)
        !-----------------------------------internal
        integer :: x
        integer :: nHiggsplus
        !-------------------------------------------

        if (np(Hplus) > 0) then
            nHiggsplus = np(Hplus)
        elseif (np(Hplus) .eq. 0) then
            nHiggsplus = 1
        else
            stop 'error in subroutine allocate_couplratio_Hc_parts (1)'
        endif

        do x = lbound(g, dim=1), ubound(g, dim=1)

            allocate (g(x)%hcjud_L(nHiggsplus), g(x)%hcjud_R(nHiggsplus))
            allocate (g(x)%hcjcs_L(nHiggsplus), g(x)%hcjcs_R(nHiggsplus))
            allocate (g(x)%hcjtb_L(nHiggsplus), g(x)%hcjtb_R(nHiggsplus))
            allocate (g(x)%hcjus_L(nHiggsplus), g(x)%hcjus_R(nHiggsplus))
            allocate (g(x)%hcjub_L(nHiggsplus), g(x)%hcjub_R(nHiggsplus))
            allocate (g(x)%hcjcd_L(nHiggsplus), g(x)%hcjcd_R(nHiggsplus))
            allocate (g(x)%hcjcb_L(nHiggsplus), g(x)%hcjcb_R(nHiggsplus))
            allocate (g(x)%hcjtd_L(nHiggsplus), g(x)%hcjtd_R(nHiggsplus))
            allocate (g(x)%hcjts_L(nHiggsplus), g(x)%hcjts_R(nHiggsplus))

            g(x)%hcjud_L = 0.0D0
            g(x)%hcjud_R = 0.0D0
            g(x)%hcjcs_L = 0.0D0
            g(x)%hcjcs_R = 0.0D0
            g(x)%hcjtb_L = 0.0D0
            g(x)%hcjtb_R = 0.0D0
            g(x)%hcjus_L = 0.0D0
            g(x)%hcjus_R = 0.0D0
            g(x)%hcjub_L = 0.0D0
            g(x)%hcjub_R = 0.0D0
            g(x)%hcjcd_L = 0.0D0
            g(x)%hcjcd_R = 0.0D0
            g(x)%hcjcb_L = 0.0D0
            g(x)%hcjcb_R = 0.0D0
            g(x)%hcjtd_L = 0.0D0
            g(x)%hcjtd_R = 0.0D0
            g(x)%hcjts_L = 0.0D0
            g(x)%hcjts_R = 0.0D0
        enddo

    end subroutine allocate_couplratio_Hc_parts

    !**********************************************************
    subroutine allocate_hadroncolliderextras_parts(tR)
        !**********************************************************
        implicit none
        !-------------------------------------------
        type(hadroncolliderextras) :: tR(:)
        !-----------------------------------internal
        integer :: x
        integer :: nHiggsneut
        !-------------------------------------------

        if (np(Hneut) > 0) then
            nHiggsneut = np(Hneut)
        elseif (np(Hneut) .eq. 0) then
            nHiggsneut = 1
        else
            stop 'error in subroutine allocate_hadroncolliderextras_parts (1)'
        endif

        tR%nq_hjWp = 2 ! (u dbar), (c sbar)  e.g
        tR%nq_hjWm = 2 ! (ubar d), (cbar s)
        tR%nq_hj = 5   !(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)
        tR%nq_hjZ = 5  !(d dbar), (u ubar), (s sbar), (c cbar), (b bbar)

        do x = lbound(tR, dim=1), ubound(tR, dim=1)
            allocate (tR(x)%qq_hjWp(tR(x)%nq_hjWp, nHiggsneut))
            allocate (tR(x)%qq_hjWm(tR(x)%nq_hjWm, nHiggsneut))
            allocate (tR(x)%gg_hj(nHiggsneut))
            allocate (tR(x)%qq_hj(tR(x)%nq_hj, nHiggsneut))
            allocate (tR(x)%gg_hjZ(nHiggsneut))
            allocate (tR(x)%qq_hjZ(tR(x)%nq_hjZ, nHiggsneut))
            allocate (tR(x)%bg_hjb(nHiggsneut))

            tR(x)%qq_hjWp = 0.0D0
            tR(x)%qq_hjWm = 0.0D0
            tR(x)%gg_hj = 0.0D0
            tR(x)%qq_hj = 0.0D0
            tR(x)%gg_hjZ = 0.0D0
            tR(x)%qq_hjZ = 0.0D0
            tR(x)%bg_hjb = 0.0D0
        enddo

    end subroutine allocate_hadroncolliderextras_parts

    !**********************************************************
    subroutine deallocate_hadroncolliderextras_parts(tR)
        !**********************************************************
        implicit none
        !--------------------------------------input
        type(hadroncolliderextras) :: tR(:)
        !-----------------------------------internal
        integer :: x
        !-------------------------------------------

        do x = lbound(tR, dim=1), ubound(tR, dim=1)
            deallocate (tR(x)%qq_hjWp)
            deallocate (tR(x)%qq_hjWm)
            deallocate (tR(x)%gg_hj)
            deallocate (tR(x)%qq_hj)
            deallocate (tR(x)%gg_hjZ)
            deallocate (tR(x)%qq_hjZ)
            deallocate (tR(x)%bg_hjb)
        enddo

    end subroutine deallocate_hadroncolliderextras_parts

    !**********************************************************
    subroutine deallocate_usefulbits
        !**********************************************************
        ! deallocates theo,res (and everything inside)
        ! deallocates c,predratio,fact
        !************************************************************
        implicit none
        !-----------------------------------internal
        integer x, y
        !-------------------------------------------
        deallocate (pdesc)!allocated in fill_pdesc

        !these are allocated in subroutine do_input
        do x = lbound(theo, dim=1), ubound(theo, dim=1)
            deallocate (theo(x)%additional)

            do y = 1, ubound(np, dim=1)
                deallocate (theo(x)%particle(y)%M)
                deallocate (theo(x)%particle(y)%GammaTot)
                deallocate (theo(x)%particle(y)%GammaTot_SM)
                deallocate (theo(x)%particle(y)%dM)
                deallocate (theo(x)%particle(y)%dMh)
            enddo
            deallocate (theo(x)%particle)

            deallocate (theo(x)%lep%XS_hjZ_ratio)
            deallocate (theo(x)%lep%XS_bbhj_ratio)
            deallocate (theo(x)%lep%XS_tautauhj_ratio)
            deallocate (theo(x)%lep%XS_hjhi_ratio)
            deallocate (theo(x)%lep%XS_HpjHmj_ratio)
            deallocate (theo(x)%lep%XS_CpjCmj)
            deallocate (theo(x)%lep%XS_NjNi)

            deallocate (theo(x)%BR_hjuu)
            deallocate (theo(x)%BR_hjdd)
            deallocate (theo(x)%BR_hjss)
            deallocate (theo(x)%BR_hjcc)
            deallocate (theo(x)%BR_hjbb)
            deallocate (theo(x)%BR_hjtt)
            deallocate (theo(x)%BR_hjuc)
            deallocate (theo(x)%BR_hjds)
            deallocate (theo(x)%BR_hjut)
            deallocate (theo(x)%BR_hjdb)
            deallocate (theo(x)%BR_hjct)
            deallocate (theo(x)%BR_hjsb)
            deallocate (theo(x)%BR_hjee)
            deallocate (theo(x)%BR_hjmumu)
            deallocate (theo(x)%BR_hjtautau)
            deallocate (theo(x)%BR_hjhihi)
            deallocate (theo(x)%BR_hjhiZ)
            deallocate (theo(x)%BR_hkhjhi)
            deallocate (theo(x)%BR_hjHpiW)
            deallocate (theo(x)%BR_hjWW)
            deallocate (theo(x)%BR_hjZZ)
            deallocate (theo(x)%BR_hjZga)
            deallocate (theo(x)%BR_hjgaga)
            deallocate (theo(x)%BR_hjgg)
            deallocate (theo(x)%BR_hjinvisible)

            deallocate (theo(x)%BR_tHpjb)
            deallocate (theo(x)%BR_Hpjcs)
            deallocate (theo(x)%BR_Hpjcb)
            deallocate (theo(x)%BR_Hpjtaunu)
            deallocate (theo(x)%BR_Hpjud)
            deallocate (theo(x)%BR_Hpjus)
            deallocate (theo(x)%BR_Hpjcd)
            deallocate (theo(x)%BR_Hpjub)
            deallocate (theo(x)%BR_Hpjenu)
            deallocate (theo(x)%BR_Hpjmunu)
            deallocate (theo(x)%BR_Hpjtb)
            deallocate (theo(x)%BR_HpjWZ)
            deallocate (theo(x)%BR_HpjhiW)

            deallocate (theo(x)%BR_CjqqNi)
            deallocate (theo(x)%BR_CjlnuNi)
            deallocate (theo(x)%BR_CjWNi)
            deallocate (theo(x)%BR_NjqqNi)
            deallocate (theo(x)%BR_NjZNi)

            deallocate (theo(x)%tev%XS_hjb_ratio)
            deallocate (theo(x)%tev%XS_tthj_ratio)
            deallocate (theo(x)%tev%XS_vbf_ratio)
            deallocate (theo(x)%tev%XS_hjZ_ratio)
            deallocate (theo(x)%tev%XS_hjW_ratio)
            deallocate (theo(x)%tev%XS_hj_ratio)
            deallocate (theo(x)%tev%XS_gg_hj_ratio)
            deallocate (theo(x)%tev%XS_bb_hj_ratio)
            deallocate (theo(x)%tev%XS_thj_tchan_ratio)
            deallocate (theo(x)%tev%XS_thj_schan_ratio)
            deallocate (theo(x)%tev%XS_tWhj_ratio)
            deallocate (theo(x)%tev%XS_hjhi)
            deallocate (theo(x)%tev%XS_vbf_Hpmj)
            deallocate (theo(x)%tev%XS_Hpmjtb)
            deallocate (theo(x)%tev%XS_Hpmjcb)
            deallocate (theo(x)%tev%XS_Hpmjbjet)
            deallocate (theo(x)%tev%XS_Hpmjcjet)
            deallocate (theo(x)%tev%XS_qq_Hpmj)
            deallocate (theo(x)%tev%XS_HpmjW)
            deallocate (theo(x)%tev%XS_HpmjZ)
            deallocate (theo(x)%tev%XS_HpjHmj)
            deallocate (theo(x)%tev%XS_Hpmjhi)
            deallocate (theo(x)%tev%channelrates)
            deallocate (theo(x)%tev%channelrates_tmp)

            deallocate (theo(x)%lhc7%XS_hjb_ratio)
            deallocate (theo(x)%lhc7%XS_tthj_ratio)
            deallocate (theo(x)%lhc7%XS_vbf_ratio)
            deallocate (theo(x)%lhc7%XS_hjZ_ratio)
            deallocate (theo(x)%lhc7%XS_qq_hjZ_ratio)
            deallocate (theo(x)%lhc7%XS_gg_hjZ_ratio)
            deallocate (theo(x)%lhc7%XS_hjW_ratio)
            deallocate (theo(x)%lhc7%XS_hj_ratio)
            deallocate (theo(x)%lhc7%XS_gg_hj_ratio)
            deallocate (theo(x)%lhc7%XS_bb_hj_ratio)
            deallocate (theo(x)%lhc7%XS_thj_tchan_ratio)
            deallocate (theo(x)%lhc7%XS_thj_schan_ratio)
            deallocate (theo(x)%lhc7%XS_tWhj_ratio)
            deallocate (theo(x)%lhc7%XS_hjhi)
            deallocate (theo(x)%lhc7%XS_vbf_Hpmj)
            deallocate (theo(x)%lhc7%XS_Hpmjtb)
            deallocate (theo(x)%lhc7%XS_Hpmjcb)
            deallocate (theo(x)%lhc7%XS_Hpmjbjet)
            deallocate (theo(x)%lhc7%XS_Hpmjcjet)
            deallocate (theo(x)%lhc7%XS_qq_Hpmj)
            deallocate (theo(x)%lhc7%XS_HpmjW)
            deallocate (theo(x)%lhc7%XS_HpmjZ)
            deallocate (theo(x)%lhc7%XS_HpjHmj)
            deallocate (theo(x)%lhc7%XS_Hpmjhi)
            deallocate (theo(x)%lhc7%channelrates)
            deallocate (theo(x)%lhc7%channelrates_tmp)

            deallocate (theo(x)%lhc8%XS_hjb_ratio)
            deallocate (theo(x)%lhc8%XS_tthj_ratio)
            deallocate (theo(x)%lhc8%XS_vbf_ratio)
            deallocate (theo(x)%lhc8%XS_hjZ_ratio)
            deallocate (theo(x)%lhc8%XS_qq_hjZ_ratio)
            deallocate (theo(x)%lhc8%XS_gg_hjZ_ratio)
            deallocate (theo(x)%lhc8%XS_hjW_ratio)
            deallocate (theo(x)%lhc8%XS_hj_ratio)
            deallocate (theo(x)%lhc8%XS_gg_hj_ratio)
            deallocate (theo(x)%lhc8%XS_bb_hj_ratio)
            deallocate (theo(x)%lhc8%XS_thj_tchan_ratio)
            deallocate (theo(x)%lhc8%XS_thj_schan_ratio)
            deallocate (theo(x)%lhc8%XS_tWhj_ratio)
            deallocate (theo(x)%lhc8%XS_hjhi)
            deallocate (theo(x)%lhc8%XS_vbf_Hpmj)
            deallocate (theo(x)%lhc8%XS_Hpmjtb)
            deallocate (theo(x)%lhc8%XS_Hpmjcb)
            deallocate (theo(x)%lhc8%XS_Hpmjbjet)
            deallocate (theo(x)%lhc8%XS_Hpmjcjet)
            deallocate (theo(x)%lhc8%XS_qq_Hpmj)
            deallocate (theo(x)%lhc8%XS_HpmjW)
            deallocate (theo(x)%lhc8%XS_HpmjZ)
            deallocate (theo(x)%lhc8%XS_HpjHmj)
            deallocate (theo(x)%lhc8%XS_Hpmjhi)
            deallocate (theo(x)%lhc8%channelrates)
            deallocate (theo(x)%lhc8%channelrates_tmp)

            deallocate (theo(x)%lhc13%XS_hjb_ratio)
            deallocate (theo(x)%lhc13%XS_tthj_ratio)
            deallocate (theo(x)%lhc13%XS_vbf_ratio)
            deallocate (theo(x)%lhc13%XS_hjZ_ratio)
            deallocate (theo(x)%lhc13%XS_qq_hjZ_ratio)
            deallocate (theo(x)%lhc13%XS_gg_hjZ_ratio)
            deallocate (theo(x)%lhc13%XS_hjW_ratio)
            deallocate (theo(x)%lhc13%XS_hj_ratio)
            deallocate (theo(x)%lhc13%XS_gg_hj_ratio)
            deallocate (theo(x)%lhc13%XS_bb_hj_ratio)
            deallocate (theo(x)%lhc13%XS_thj_tchan_ratio)
            deallocate (theo(x)%lhc13%XS_thj_schan_ratio)
            deallocate (theo(x)%lhc13%XS_tWhj_ratio)
            deallocate (theo(x)%lhc13%XS_hjhi)
            deallocate (theo(x)%lhc13%XS_vbf_Hpmj)
            deallocate (theo(x)%lhc13%XS_Hpmjtb)
            deallocate (theo(x)%lhc13%XS_Hpmjcb)
            deallocate (theo(x)%lhc13%XS_Hpmjbjet)
            deallocate (theo(x)%lhc13%XS_Hpmjcjet)
            deallocate (theo(x)%lhc13%XS_qq_Hpmj)
            deallocate (theo(x)%lhc13%XS_HpmjW)
            deallocate (theo(x)%lhc13%XS_HpmjZ)
            deallocate (theo(x)%lhc13%XS_HpjHmj)
            deallocate (theo(x)%lhc13%XS_Hpmjhi)
            deallocate (theo(x)%lhc13%channelrates)
            deallocate (theo(x)%lhc13%channelrates_tmp)

            !deallocate(theo(x)%inLEPrange_Hpj)
            !deallocate(theo(x)%inTEVrange_Hpj)

            deallocate (theo(x)%CP_value)
        enddo

        select case (whichanalyses)
        case ('onlyH', 'LandH', 'onlyP', 'list ')
            do x = lbound(theo, dim=1), ubound(theo, dim=1)
                deallocate (theo(x)%BR_Hbb_SM)
                deallocate (theo(x)%BR_Hss_SM)
                deallocate (theo(x)%BR_Hcc_SM)
                deallocate (theo(x)%BR_Hmumu_SM)
                deallocate (theo(x)%BR_Htautau_SM)
                deallocate (theo(x)%BR_HWW_SM)
                deallocate (theo(x)%BR_HZZ_SM)
                deallocate (theo(x)%BR_HZga_SM)
                deallocate (theo(x)%BR_Hgaga_SM)
                deallocate (theo(x)%BR_Hgg_SM)
                deallocate (theo(x)%BR_Hjets_SM)
!                 deallocate (theo(x)%GammaTot_SM)

                deallocate (theo(x)%tev%XS_HZ_SM)
                deallocate (theo(x)%tev%XS_gg_HZ_SM)
                deallocate (theo(x)%tev%XS_qq_HZ_SM)
                deallocate (theo(x)%tev%XS_HW_SM)
                deallocate (theo(x)%tev%XS_H_SM)
                deallocate (theo(x)%tev%XS_gg_H_SM)
                deallocate (theo(x)%tev%XS_bb_H_SM)
                deallocate (theo(x)%tev%XS_ttH_SM)
                deallocate (theo(x)%tev%XS_vbf_SM)
                !deallocate(theo(x)%tev%XS_H_SM_9713)
                !deallocate(theo(x)%tev%XS_H_SM_9674)
                deallocate (theo(x)%tev%XS_tH_tchan_SM)
                deallocate (theo(x)%tev%XS_tH_schan_SM)
                deallocate (theo(x)%tev%XS_tWH_SM)
                deallocate (theo(x)%tev%channelrates_SM)

                deallocate (theo(x)%tev%XS_Hb_SM)
                deallocate (theo(x)%tev%XS_Hb_c1_SM)
                deallocate (theo(x)%tev%XS_Hb_c2_SM)
                deallocate (theo(x)%tev%XS_Hb_c3_SM)
                deallocate (theo(x)%tev%XS_Hb_c4_SM)

                deallocate (theo(x)%lhc7%XS_HZ_SM)
                deallocate (theo(x)%lhc7%XS_gg_HZ_SM)
                deallocate (theo(x)%lhc7%XS_qq_HZ_SM)
                deallocate (theo(x)%lhc7%XS_HW_SM)
                deallocate (theo(x)%lhc7%XS_H_SM)
                deallocate (theo(x)%lhc7%XS_gg_H_SM)
                deallocate (theo(x)%lhc7%XS_bb_H_SM)
                deallocate (theo(x)%lhc7%XS_ttH_SM)
                deallocate (theo(x)%lhc7%XS_vbf_SM)
                deallocate (theo(x)%lhc7%XS_tH_tchan_SM)
                deallocate (theo(x)%lhc7%XS_tH_schan_SM)
                deallocate (theo(x)%lhc7%XS_tWH_SM)
                deallocate (theo(x)%lhc7%XS_Hb_SM)
                deallocate (theo(x)%lhc7%channelrates_SM)
!     deallocate(theo(x)%lhc7%XS_Hb_c1_SM)
!     deallocate(theo(x)%lhc7%XS_Hb_c2_SM)
!     deallocate(theo(x)%lhc7%XS_Hb_c3_SM)

                deallocate (theo(x)%lhc8%XS_HZ_SM)
                deallocate (theo(x)%lhc8%XS_gg_HZ_SM)
                deallocate (theo(x)%lhc8%XS_qq_HZ_SM)
                deallocate (theo(x)%lhc8%XS_HW_SM)
                deallocate (theo(x)%lhc8%XS_H_SM)
                deallocate (theo(x)%lhc8%XS_gg_H_SM)
                deallocate (theo(x)%lhc8%XS_bb_H_SM)
                deallocate (theo(x)%lhc8%XS_ttH_SM)
                deallocate (theo(x)%lhc8%XS_vbf_SM)
                deallocate (theo(x)%lhc8%XS_tH_tchan_SM)
                deallocate (theo(x)%lhc8%XS_tH_schan_SM)
                deallocate (theo(x)%lhc8%XS_tWH_SM)
                deallocate (theo(x)%lhc8%XS_Hb_SM)
                deallocate (theo(x)%lhc8%channelrates_SM)
!     deallocate(theo(x)%lhc8%XS_Hb_c1_SM)
!     deallocate(theo(x)%lhc8%XS_Hb_c2_SM)
!     deallocate(theo(x)%lhc8%XS_Hb_c3_SM)

                deallocate (theo(x)%lhc13%XS_HZ_SM)
                deallocate (theo(x)%lhc13%XS_gg_HZ_SM)
                deallocate (theo(x)%lhc13%XS_qq_HZ_SM)
                deallocate (theo(x)%lhc13%XS_HW_SM)
                deallocate (theo(x)%lhc13%XS_H_SM)
                deallocate (theo(x)%lhc13%XS_gg_H_SM)
                deallocate (theo(x)%lhc13%XS_bb_H_SM)
                deallocate (theo(x)%lhc13%XS_ttH_SM)
                deallocate (theo(x)%lhc13%XS_vbf_SM)
                deallocate (theo(x)%lhc13%XS_tH_tchan_SM)
                deallocate (theo(x)%lhc13%XS_tH_schan_SM)
                deallocate (theo(x)%lhc13%XS_tWH_SM)
                deallocate (theo(x)%lhc13%channelrates_SM)

            enddo
        case ('onlyL')
        case default
            stop 'error in deallocate_usefulbits'
        end select

        deallocate (theo) !allocated in subroutine do_input

        !allocated in subroutine setup_output
        if (allocated(res)) then
            do x = lbound(res, dim=1), ubound(res, dim=1)
                deallocate (res(x)%chan)
                deallocate (res(x)%obsratio)
                deallocate (res(x)%predratio)
                deallocate (res(x)%axis_i)
                deallocate (res(x)%axis_j)
                deallocate (res(x)%sfactor)
                deallocate (res(x)%allowed95)
                deallocate (res(x)%ncombined)
            enddo
            deallocate (res) !allocated in subroutine setup_output
        endif

        if (allocated(fullHBres)) then
            deallocate (fullHBres)
        endif

!  call deallocate_sqcouplratio_parts(g2)
        do x = lbound(g2, dim=1), ubound(g2, dim=1)

            deallocate (g2(x)%hjss_s)
            deallocate (g2(x)%hjss_p)
            deallocate (g2(x)%hjcc_s)
            deallocate (g2(x)%hjcc_p)
            deallocate (g2(x)%hjbb_s)
            deallocate (g2(x)%hjbb_p)
            deallocate (g2(x)%hjtoptop_s)
            deallocate (g2(x)%hjtoptop_p)
            deallocate (g2(x)%hjmumu_s)
            deallocate (g2(x)%hjmumu_p)
            deallocate (g2(x)%hjtautau_s)
            deallocate (g2(x)%hjtautau_p)

            deallocate (g2(x)%hjWW)
            deallocate (g2(x)%hjZZ)
            deallocate (g2(x)%hjZga)
            deallocate (g2(x)%hjgaga)
            deallocate (g2(x)%hjgg)
            deallocate (g2(x)%hjggZ)
            deallocate (g2(x)%hjhiZ)
        enddo
        deallocate (g2)

        do x = lbound(effC, dim=1), ubound(effC, dim=1)

            deallocate (effC(x)%hjee_s)
            deallocate (effC(x)%hjee_p)
            deallocate (effC(x)%hjuu_s)
            deallocate (effC(x)%hjuu_p)
            deallocate (effC(x)%hjdd_s)
            deallocate (effC(x)%hjdd_p)
            deallocate (effC(x)%hjuc_s)
            deallocate (effC(x)%hjuc_p)
            deallocate (effC(x)%hjut_s)
            deallocate (effC(x)%hjut_p)
            deallocate (effC(x)%hjct_s)
            deallocate (effC(x)%hjct_p)
            deallocate (effC(x)%hjds_s)
            deallocate (effC(x)%hjds_p)
            deallocate (effC(x)%hjdb_s)
            deallocate (effC(x)%hjdb_p)
            deallocate (effC(x)%hjsb_s)
            deallocate (effC(x)%hjsb_p)

            deallocate (effC(x)%hjss_s)
            deallocate (effC(x)%hjss_p)
            deallocate (effC(x)%hjcc_s)
            deallocate (effC(x)%hjcc_p)
            deallocate (effC(x)%hjbb_s)
            deallocate (effC(x)%hjbb_p)
            deallocate (effC(x)%hjtt_s)
            deallocate (effC(x)%hjtt_p)
            deallocate (effC(x)%hjmumu_s)
            deallocate (effC(x)%hjmumu_p)
            deallocate (effC(x)%hjtautau_s)
            deallocate (effC(x)%hjtautau_p)

            deallocate (effC(x)%hjWW)
            deallocate (effC(x)%hjZZ)
            deallocate (effC(x)%hjZga)
            deallocate (effC(x)%hjgaga)
            deallocate (effC(x)%hjgg)
!    deallocate(effC(x)%hjggZ)
            deallocate (effC(x)%hjhiZ)
        enddo
        deallocate (effC)

        !these are allocated in subroutine do_input
        call deallocate_hadroncolliderextras_parts(partR)
        deallocate (partR) !allocated in subroutine do_input

        if (allocated(pr)) deallocate (pr)  !allocated in subroutine fill_pr or fill_pr_select
        if (allocated(prsep)) deallocate (prsep)  !allocated in subroutine fill_pr or fill_pr_select

        if (allocated(diffMhneut)) deallocate (diffMhneut)
        if (allocated(diffMhch)) deallocate (diffMhch)

        if (allocated(dmn)) deallocate (dmn)
        if (allocated(dmch)) deallocate (dmch)

        if (allocated(analysislist)) deallocate (analysislist)
        if (allocated(analysis_exclude_list)) deallocate (analysis_exclude_list)

        if (allocated(HBresult_all)) deallocate (HBresult_all)
        if (allocated(chan_all)) deallocate (chan_all)
        if (allocated(ncombined_all)) deallocate (ncombined_all)
        if (allocated(obsratio_all)) deallocate (obsratio_all)
        if (allocated(predratio_all)) deallocate (predratio_all)

    end subroutine deallocate_usefulbits
    !**********************************************************

    function BR_hjdijet(t, j, maxnbjet, gluons)
        implicit none

        type(dataset), intent(in) :: t
        integer, intent(in) :: j, maxnbjet
        double precision BR_hjdijet
        logical, intent(in) :: gluons

        BR_hjdijet = t%BR_hjuu(j) + t%BR_hjdd(j) + t%BR_hjcc(j) + t%BR_hjss(j) + &
                     t%BR_hjuc(j) + t%BR_hjds(j)

        if (gluons) BR_hjdijet = BR_hjdijet + t%BR_hjgg(j)
        if (maxnbjet >= 1) then
            BR_hjdijet = BR_hjdijet + t%BR_hjdb(j) + t%BR_hjsb(j)
        endif
        if (maxnbjet >= 2) then
            BR_hjdijet = BR_hjdijet + t%BR_hjbb(j)
        endif
        return
    end function

    function BR_Hpjdijet(t, j, maxnbjet)
        implicit none
        type(dataset), intent(in) :: t
        integer, intent(in) :: j, maxnbjet
        double precision BR_Hpjdijet

        BR_Hpjdijet = t%BR_Hpjud(j) + t%BR_Hpjus(j) + t%BR_Hpjcd(j) + t%BR_Hpjcs(j)
        if (maxnbjet >= 1) then
            BR_Hpjdijet = BR_Hpjdijet + t%BR_Hpjub(j) + t%BR_Hpjcb(j)
        endif
        return
    end function

    !> Returns an array of indices `d` such that `r(d)` is sorted in ascending order.
    !! Uses mergesort, slightly adapted from github.com/Astrokiwi/simple_fortran_argsort
    subroutine argsort(r, d)
        double precision, intent(in), dimension(:) :: r
        integer, intent(out), dimension(size(r)) :: d

        integer, dimension(size(r)) :: il

        integer :: stepsize
        integer :: i, j, n, left, k, ksize

        n = size(r)

        do i = 1, n
            d(i) = i
        end do

        if (n == 1) return

        stepsize = 1
        do while (stepsize < n)
            do left = 1, n - stepsize, stepsize * 2
                i = left
                j = left + stepsize
                ksize = min(stepsize * 2, n - left + 1)
                k = 1

                do while (i < left + stepsize .and. j < left + ksize)
                    if (r(d(i)) < r(d(j))) then
                        il(k) = d(i)
                        i = i + 1
                        k = k + 1
                    else
                        il(k) = d(j)
                        j = j + 1
                        k = k + 1
                    endif
                enddo

                if (i < left + stepsize) then
                    ! fill up remaining from left
                    il(k:ksize) = d(i:left + stepsize - 1)
                else
                    ! fill up remaining from right
                    il(k:ksize) = d(j:left + ksize - 1)
                endif
                d(left:left + ksize - 1) = il(1:ksize)
            end do
            stepsize = stepsize * 2
        end do

        return
    end subroutine

end module usefulbits
!******************************************************************
