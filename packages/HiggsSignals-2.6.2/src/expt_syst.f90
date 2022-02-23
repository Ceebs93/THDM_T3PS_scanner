module expt_syst

! use usefulbits_hs
    implicit none

    integer, parameter :: Nprod = 12
! ggH, VBF, WH, ZH, ttH
    integer, parameter :: Ndecay = 9
    integer, parameter :: Nsyst = 57
!  1: CMS H->gaga untagged 0-1 7 TeV event migration
!  2: CMS H->gaga untagged 1-2 7 TeV event migration
!  3: CMS H->gaga untagged 2-3 7 TeV event migration
!  4: CMS H->gaga untagged 0-1 8 TeV event migration
!  5: CMS H->gaga untagged 1-2 8 TeV event migration
!  6: CMS H->gaga untagged 2-3 8 TeV event migration
!  7: CMS H->gaga dijet 8 TeV event migration 0-1
!  8: Dijet tagging efficiency in dijet selection of CMS H->gaga analyses
!  9: ETmiss cut efficiency in ETmiss selection of CMS H->gaga analyses

! 10: ATLAS H->tautau ggH differential pT distribution and QCD scale
! 11: ATLAS H->tautau Top and Z->ll BG normalization (for hadlep and leplep channels)
! 12: ATLAS H->tautau hadronic tau identification and energy scale
! 13: ATLAS H->tautau JES eta calibration
! 14: ATLAS H->tautau Z->tautau normalization (for hadlep)
! 15: ATLAS H->tautau fake backgrounds (for leplep)
! 16: ATLAS H->tautau ditau(had) tagging efficiency

! 17: ATLAS H->gaga gg->H + (>2jets) cross section (affects VBF and VH(hadronic) channel)
! 18: ATLAS H->gaga gg->H + (3 jets) cross section (affect only VBF channels)
! 19: ATLAS H->gaga Underlying Event on gg->H yield
! 20: ATLAS H->gaga Underlying Event on ttH yield
! 21: ATLAS H->gaga pT spectrum modelling: migration between high- and low-pT categories of gg->H yield (central)
! 22: ATLAS H->gaga pT spectrum modelling: migration between high- and low-pT categories of gg->H yield (forward)
! 23: ATLAS H->gaga 2-jet Delta(Phi) angular distribution of gg->H in VBF categories
! 24: ATLAS H->gaga 2-jet Delta(Eta*) angular distribution of gg->H in VBF categories
! 25: ATLAS H->gaga gg->H contribution to ttH categories
! 26: ATLAS H->gaga heavy flavor fraction of gg->H, VBF, WH contribution to ttH categories
! 27: ATLAS H->gaga experimental: jet energy scale / resolution and vertex fraction
! 28: ATLAS H->gaga experimental: ETmiss energy scale and resolution
! 29: ATLAS H->gaga diphoton mass resolution (Tab. XII): Constant term
! 30: ATLAS H->gaga diphoton mass resolution (Tab. XII): Sampling term
! 31: ATLAS H->gaga diphoton mass resolution (Tab. XII): Material modeling
! 32: ATLAS H->gaga diphoton mass resolution (Tab. XII): Noise term
! 33: ATLAS H->gaga experimental: Photon Isolation (Tab. IX)
! 34: ATLAS H->gaga experimental: Photon ID (Tab. IX)
! 35: ATLAS H->gaga experimental: b-tagging efficiency
! 36: ATLAS H->gaga migration (not official): between VBF categories
! 37: ATLAS H->gaga migration (not official): between gg->H(forward) and VBF categories
! 38: ATLAS H->gaga unofficial VH error

! 39: CMS H->gaga untagged 3-4 8 TeV event migration
! 40: CMS H->gaga dijet 7 TeV event migration 0-1
! 41: CMS H->gaga dijet 8 TeV event migration 1-2
! 42: CMS H->gaga VH tight-loose migration 8 TeV
! 43: CMS H->gaga VH tight-loose common systematic 8 TeV
! 44: CMS H->gaga ttH(lepton)-VH tight migration
! 45: CMS H->gaga VH(loose)-VBF(dijet) migration
! 46: CMS H->gaga VBF(dijet2)-VH(ETmiss) migration
! 47: CMS H->gaga VH(ETmiss)-ttH multijet migration
! 48: CMS H->gaga ttH(multijet)-VH(dijet) migration
! 49: CMS H->gaga VH(dijet)-untagged migration
! 50: CMS H->gaga VBF uncertainty in untagged

! NEW 13 TeV results:
! CMS-16-020, H->gaga
! 51: UE and parton shower, jet energy scale/smearing: VBF-VBF Tag migration, 7% and 4-15%, respectively.
! 52: UE and parton shower, jet energy scale/smearing: VBF-untagged Tag migration, 9% and 4-15%, respectively.
! 53: Event migration untagged 0-1
! 54: Event migration untagged 1-2
! 55: Event migration untagged 2-3
! n.b: QCD scale uncertainty and energy scale/resolution, all categories, ~5-10% and ~6%, respectively. -> makes fit worse!
! 56: ggF contamination in VBF, ttH categories, ~39%
! 57: ggF contamination: VBF-VBF Tag migration, ~10%

    double precision :: rel_corr_err(2, 0:Nprod, 0:Ndecay, Nsyst)

! scaletype determines whether the systematic uncertainty is scaled with the
! observed (typical for BG uncertainties) [0] or with the predicted mu (typical
! for signal uncertainties) [1].
    integer :: scaletype(Nsyst)

contains

    !------------------------------------------------
    subroutine fill_scaletype
        !------------------------------------------------
        ! Set scaletypes (default is scaling with predicted)
        scaletype(:) = 1
        ! Set to 0 for mostly background-affecting systematics
        scaletype(11) = 0
        scaletype(12) = 0
        scaletype(13) = 0
        scaletype(14) = 0
        scaletype(15) = 0
        scaletype(16) = 0
        scaletype(27) = 0
        scaletype(28) = 0
        scaletype(29) = 0
        scaletype(30) = 0
        scaletype(31) = 0
        scaletype(32) = 0
        scaletype(33) = 0
        scaletype(34) = 0
        scaletype(35) = 0
!  scaletype(53)=0

    end subroutine fill_scaletype
!------------------------------------------------
    subroutine fill_rel_corr_err(ID, N)
        !------------------------------------------------
        implicit none
        integer, intent(in) :: ID, N

        rel_corr_err(N, :, :, :) = 0.0D0

        select case (ID)
! ----------------------------------------------------------------------
! This is for an outdated CMS CONF-Note for H->gamma gamma, which should not
! be used in conjunction with the updated CMS results for this analysis (0558...)
! ----------------------------------------------------------------------
        case (13001107) ! untagged 0 7TeV
            rel_corr_err(N, :, 1, 1) = +0.125D0
        case (13001108) ! untagged 1 7TeV
            rel_corr_err(N, :, 1, 1) = -0.125D0
            rel_corr_err(N, :, 1, 2) = +0.125D0
        case (13001109) ! untagged 2 7TeV
            rel_corr_err(N, :, 1, 2) = -0.125D0
            rel_corr_err(N, :, 1, 3) = +0.125D0
        case (13001110) ! untagged 3 7TeV
            rel_corr_err(N, :, 1, 3) = -0.125D0
        case (13001111) ! untagged 0 8TeV
            rel_corr_err(N, :, 1, 4) = +0.125D0
        case (13001112) ! untagged 1 8TeV
            rel_corr_err(N, :, 1, 4) = -0.125D0
            rel_corr_err(N, :, 1, 5) = +0.125D0
        case (13001113) ! untagged 2 8TeV
            rel_corr_err(N, :, 1, 5) = -0.125D0
            rel_corr_err(N, :, 1, 6) = +0.125D0
        case (13001114) ! untagged 3 8TeV
            rel_corr_err(N, :, 1, 6) = -0.125D0
        case (13001105) ! CMS H->gaga dijet loose tagged categories (8 TeV)
            rel_corr_err(N, :, 1, 7) = +0.15D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
        case (13001106) ! CMS H->gaga dijet tight tagged categories (8 TeV)
            rel_corr_err(N, :, 1, 7) = -0.15D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.15D0
        case (12015103) ! CMS H->gaga dijet tagged category (7 TeV)
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
        case (13001102) ! CMS H->gaga ETmiss tagged categories
            rel_corr_err(N, 1, 1, 9) = 0.15D0
            rel_corr_err(N, 2, 1, 9) = 0.15D0
            rel_corr_err(N, 3, 1, 9) = 0.04D0
            rel_corr_err(N, 4, 1, 9) = 0.04D0
            rel_corr_err(N, 5, 1, 9) = 0.04D0
! ----------------------------------------------------------------------
! updated (full LHC I run) CMS H->gamma gamma results
! partly using the same systematics than in older results
! ----------------------------------------------------------------------
        case (0558101) ! untagged 0 7TeV
            rel_corr_err(N, :, 1, 1) = +0.2D0
        case (0558102) ! untagged 1 7TeV
            rel_corr_err(N, :, 1, 1) = -0.2D0
            rel_corr_err(N, :, 1, 2) = +0.2D0
        case (0558103) ! untagged 2 7TeV
            rel_corr_err(N, :, 1, 2) = -0.2D0
            rel_corr_err(N, :, 1, 3) = +0.2D0
        case (0558104) ! untagged 3 7TeV
            rel_corr_err(N, :, 1, 3) = -0.2D0
        case (0558111) ! untagged 0 8TeV
            rel_corr_err(N, :, 1, 4) = +0.2D0
        case (0558112) ! untagged 1 8TeV
            rel_corr_err(N, :, 1, 4) = -0.2D0
            rel_corr_err(N, :, 1, 5) = +0.2D0
        case (0558113) ! untagged 2 8TeV
            rel_corr_err(N, :, 1, 5) = -0.2D0
            rel_corr_err(N, :, 1, 6) = +0.2D0
        case (0558114) ! untagged 3 8TeV
            rel_corr_err(N, :, 1, 6) = -0.2D0
            rel_corr_err(N, :, 1, 39) = +0.2D0
        case (0558115) !  untagged 4 8TeV
            rel_corr_err(N, :, 1, 39) = -0.2D0
        case (0558116) ! CMS H->gaga VBF dijet 0 categories (8 TeV)
            rel_corr_err(N, :, 1, 7) = +0.3D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
            rel_corr_err(N, :, 1, 45) = -0.2D0
        case (0558117) ! CMS H->gaga VBF dijet 1 categories (8 TeV)
            rel_corr_err(N, :, 1, 7) = -0.3D0
            rel_corr_err(N, :, 1, 41) = +0.3D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.15D0
        case (0558118) ! CMS H->gaga VBF dijet 2 categories (8 TeV)
            rel_corr_err(N, :, 1, 41) = -0.3D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.15D0
            rel_corr_err(N, :, 1, 46) = +0.2D0
        case (0558105) ! CMS H->gaga VBF dijet 0 categories (7 TeV)
            rel_corr_err(N, :, 1, 40) = +0.15D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.05D0
            rel_corr_err(N, :, 1, 44) = +0.2D0
        case (0558106) ! CMS H->gaga VBF dijet 1 categories (7 TeV)
            rel_corr_err(N, :, 1, 40) = -0.15D0
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.05D0
            rel_corr_err(N, :, 1, 44) = +0.2D0
        case (0558108) ! CMS H->gaga VH ETmiss 7 TeV
            rel_corr_err(N, 1, 1, 9) = 0.15D0
            rel_corr_err(N, 2, 1, 9) = 0.15D0
            rel_corr_err(N, 3, 1, 9) = 0.04D0
            rel_corr_err(N, 4, 1, 9) = 0.04D0
            rel_corr_err(N, 5, 1, 9) = 0.04D0
        case (0558121) ! CMS H->gaga VH ETmiss 8 TeV
            rel_corr_err(N, 1, 1, 9) = 0.15D0
            rel_corr_err(N, 2, 1, 9) = 0.15D0
            rel_corr_err(N, 3, 1, 9) = 0.04D0
            rel_corr_err(N, 4, 1, 9) = 0.04D0
            rel_corr_err(N, 5, 1, 9) = 0.04D0
            rel_corr_err(N, :, 1, 46) = -0.2D0
            rel_corr_err(N, :, 1, 47) = +0.2D0
        case (0558109) ! CMS H->gaga VH dijet 7 TeV
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
        case (0558110) ! CMS H->gaga ttH multijets 7 TeV
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
        case (0558122) ! CMS H->gaga VH dijet 8 TeV
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
            rel_corr_err(N, :, 1, 48) = +0.2D0
        case (0558124) ! CMS H->gaga ttH multijets 8 TeV
            rel_corr_err(N, 1, 1, 8) = -0.3D0
            rel_corr_err(N, 2, 1, 8) = +0.1D0
            rel_corr_err(N, :, 1, 47) = -0.2D0
            rel_corr_err(N, :, 1, 48) = +0.2D0
        case (0558119) ! CMS H->gaga VH tight
            rel_corr_err(N, :, 1, 42) = +0.2D0
            rel_corr_err(N, :, 1, 43) = +0.15D0
            rel_corr_err(N, :, 1, 44) = -0.2D0
        case (0558120) ! CMS H->gaga VH loose
            rel_corr_err(N, :, 1, 42) = -0.2D0
            rel_corr_err(N, :, 1, 43) = +0.15D0
            rel_corr_err(N, :, 1, 45) = +0.2D0
        case (0558123) ! CMS H->gaga ttH lepton
            rel_corr_err(N, :, 1, 44) = +0.2D0
! ----------------------------------------------------------------------
        case (2012160101, 2014061101) ! ATL H->tautau leplep boosted category
            rel_corr_err(N, 1, 4, 10) = 0.32D0
            rel_corr_err(N, :, 4, 11) = 0.15D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
            rel_corr_err(N, :, 4, 15) = 0.12D0
        case (2012160102, 2014061102) ! ATL H->tautau leplep VBF category
            rel_corr_err(N, 1, 4, 10) = 0.08D0
            rel_corr_err(N, :, 4, 11) = 0.15D0
            rel_corr_err(N, :, 4, 13) = 0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
            rel_corr_err(N, :, 4, 15) = 0.12D0
        case (2012160103, 2014061103) ! ATL H->tautau hadlep boosted category
            rel_corr_err(N, 1, 4, 10) = 0.32D0
            rel_corr_err(N, :, 4, 11) = 0.15D0
            rel_corr_err(N, :, 4, 12) = 0.04D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
            rel_corr_err(N, :, 4, 14) = 0.10D0
        case (2012160104, 2014061104) ! ATL H->tautau hadlep VBF category
            rel_corr_err(N, 1, 4, 10) = 0.08D0
            rel_corr_err(N, :, 4, 11) = 0.15D0
            rel_corr_err(N, :, 4, 12) = 0.04D0!
            rel_corr_err(N, :, 4, 13) = 0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
            rel_corr_err(N, :, 4, 14) = 0.10D0
        case (2012160105, 2014061105) ! ATL H->tautau hadhad boosted category
            rel_corr_err(N, 1, 4, 10) = 0.22D0
            rel_corr_err(N, :, 4, 12) = 0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
            rel_corr_err(N, :, 4, 16) = 0.07D0
        case (2012160106, 2014061106) ! ATL H->tautau hadhad VBF category
            rel_corr_err(N, 1, 4, 10) = 0.05D0
            rel_corr_err(N, :, 4, 12) = 0.12D0
            rel_corr_err(N, :, 4, 13) = 0.12D0
!    rel_corr_err(N,2,4,13)= -0.12D0
            rel_corr_err(N, :, 4, 16) = 0.07D0
!---NEW SUMMER 2014 results
        case (7084101) ! ATLAS H->gaga central low-pT
            rel_corr_err(N, 1, 1, 17) = -0.10D0 ! n.b.: not official, guess of migration to VBF/VH(had) categories (>2jet)
            rel_corr_err(N, 1, 1, 19) = -0.05D0 ! n.b.: not official, guess of migration due to UE
            rel_corr_err(N, 1, 1, 21) = 0.24D0
            rel_corr_err(N, 1, 1, 27) = 0.001D0
            rel_corr_err(N, 2, 1, 27) = 0.029D0
            rel_corr_err(N, 3, 1, 27) = 0.001D0
            rel_corr_err(N, 4, 1, 27) = 0.001D0
            rel_corr_err(N, 5, 1, 27) = 0.04D0
            rel_corr_err(N, 3, 1, 28) = 0.001D0
            rel_corr_err(N, 4, 1, 28) = 0.002D0
            rel_corr_err(N, 5, 1, 28) = 0.002D0
            rel_corr_err(N, :, 1, 29) = 0.075D0
            rel_corr_err(N, :, 1, 30) = 0.026D0
            rel_corr_err(N, :, 1, 31) = 0.049D0
            rel_corr_err(N, :, 1, 32) = 0.026D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 37) = -0.50D0 ! n.B.: not official, trying migration between VBF categories
!    rel_corr_err(N,1,1,39)=  0.15D0 ! n.B.: overall gg->H uncertainty in untagged categories
        case (7084102) ! ATLAS H->gaga central high-pT
            rel_corr_err(N, 1, 1, 17) = -0.10D0 ! n.b.: not official, guess of migration to VBF/VH(had) categories (>2jet)
            rel_corr_err(N, 1, 1, 19) = -0.05D0 ! n.b.: not official, guess of migration due to UE
            rel_corr_err(N, 1, 1, 21) = -0.24D0
            rel_corr_err(N, 1, 1, 27) = 0.011D0
            rel_corr_err(N, 2, 1, 27) = 0.045D0
            rel_corr_err(N, 3, 1, 27) = 0.014D0
            rel_corr_err(N, 4, 1, 27) = 0.014D0
            rel_corr_err(N, 5, 1, 27) = 0.035D0
            rel_corr_err(N, 3, 1, 28) = 0.001D0
            rel_corr_err(N, 4, 1, 28) = 0.002D0
            rel_corr_err(N, 5, 1, 28) = 0.002D0
            rel_corr_err(N, :, 1, 29) = 0.096D0
            rel_corr_err(N, :, 1, 30) = 0.056D0
            rel_corr_err(N, :, 1, 31) = 0.062D0
            rel_corr_err(N, :, 1, 32) = 0.017D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 37) = -0.500D0 ! n.B.: not official, trying migration between VBF categories
!    rel_corr_err(N,1,1,39)=  0.15D0 ! n.B.: overall gg->H uncertainty in untagged categories
        case (7084103) ! ATLAS H->gaga forward low-pT
            rel_corr_err(N, 1, 1, 17) = -0.15D0 ! n.b.: not official, guess of migration to VBF/VH(had) categories  (>2jet)
            rel_corr_err(N, 1, 1, 18) = -0.10D0   ! n.b.: not official, guess of migration to VBF categories (3jet)
            rel_corr_err(N, 1, 1, 19) = -0.05D0 ! n.b.: not official, guess of migration due to UE
            rel_corr_err(N, 1, 1, 22) = 0.24D0
            rel_corr_err(N, 1, 1, 27) = 0.001D0
            rel_corr_err(N, 2, 1, 27) = 0.029D0
            rel_corr_err(N, 3, 1, 27) = 0.001D0
            rel_corr_err(N, 4, 1, 27) = 0.001D0
            rel_corr_err(N, 5, 1, 27) = 0.04D0
            rel_corr_err(N, 3, 1, 28) = 0.001D0
            rel_corr_err(N, 4, 1, 28) = 0.002D0
            rel_corr_err(N, 5, 1, 28) = 0.002D0
            rel_corr_err(N, :, 1, 29) = 0.099D0
            rel_corr_err(N, :, 1, 30) = 0.013D0
            rel_corr_err(N, :, 1, 31) = 0.060D0
            rel_corr_err(N, :, 1, 32) = 0.021D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 37) = -0.300D0 ! n.B.: not official, trying migration between VBF categories
!    rel_corr_err(N,1,1,39)=  0.15D0 ! n.B.: overall gg->H uncertainty in untagged categories
        case (7084104) ! ATLAS H->gaga forward high-pT
            rel_corr_err(N, 1, 1, 17) = -0.15D0 ! n.b.: not official, guess of migration to VBF/VH(had) categories (>2jet)
            rel_corr_err(N, 1, 1, 18) = -0.10D0   ! n.b.: not official, guess of migration to VBF categories (3jet)
            rel_corr_err(N, 1, 1, 19) = -0.05D0 ! n.b.: not official, guess of migration due to UE
            rel_corr_err(N, 1, 1, 22) = -0.24D0
            rel_corr_err(N, 1, 1, 27) = 0.011D0
            rel_corr_err(N, 2, 1, 27) = 0.045D0
            rel_corr_err(N, 3, 1, 27) = 0.014D0
            rel_corr_err(N, 4, 1, 27) = 0.014D0
            rel_corr_err(N, 5, 1, 27) = 0.035D0
            rel_corr_err(N, 3, 1, 28) = 0.001D0
            rel_corr_err(N, 4, 1, 28) = 0.002D0
            rel_corr_err(N, 5, 1, 28) = 0.002D0
            rel_corr_err(N, :, 1, 29) = 0.120D0
            rel_corr_err(N, :, 1, 30) = 0.028D0
            rel_corr_err(N, :, 1, 31) = 0.078D0
            rel_corr_err(N, :, 1, 32) = 0.019D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 37) = -0.300D0 ! n.B.: not official, trying migration between VBF categories
!    rel_corr_err(N,1,1,39)=  0.15D0 ! n.B.: overall gg->H uncertainty in untagged categories
        case (7084105) ! ATLAS H->gaga VBF loose
            rel_corr_err(N, 1, 1, 17) = 0.20D0
            rel_corr_err(N, 1, 1, 18) = 0.25D0
            rel_corr_err(N, 1, 1, 19) = 0.06D0
            rel_corr_err(N, 1, 1, 23) = 0.089D0
            rel_corr_err(N, 1, 1, 24) = 0.048D0
            rel_corr_err(N, 1, 1, 27) = 0.120D0
            rel_corr_err(N, 2, 1, 27) = 0.044D0
            rel_corr_err(N, 3, 1, 27) = 0.130D0
            rel_corr_err(N, 4, 1, 27) = 0.130D0
            rel_corr_err(N, 5, 1, 27) = 0.076D0
            rel_corr_err(N, 3, 1, 28) = 0.001D0
            rel_corr_err(N, 4, 1, 28) = 0.002D0
            rel_corr_err(N, 5, 1, 28) = 0.010D0
            rel_corr_err(N, :, 1, 29) = 0.094D0
            rel_corr_err(N, :, 1, 30) = 0.026D0
            rel_corr_err(N, :, 1, 31) = 0.060D0
            rel_corr_err(N, :, 1, 32) = 0.021D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 36) = 0.500D0 ! n.B.: not official, add theory error for VBF
            rel_corr_err(N, 2, 1, 37) = 0.250D0 ! n.B.: not official, trying migration between VBF categories
        case (7084106) ! ATLAS H->gaga VBF tight
            rel_corr_err(N, 1, 1, 17) = 0.20D0
            rel_corr_err(N, 1, 1, 18) = 0.52D0
            rel_corr_err(N, 1, 1, 19) = 0.06D0
            rel_corr_err(N, 1, 1, 23) = 0.112D0
            rel_corr_err(N, 1, 1, 24) = 0.066D0
            rel_corr_err(N, 1, 1, 27) = 0.130D0
            rel_corr_err(N, 2, 1, 27) = 0.091D0
            rel_corr_err(N, 3, 1, 27) = 0.170D0
            rel_corr_err(N, 4, 1, 27) = 0.170D0
            rel_corr_err(N, 5, 1, 27) = 0.063D0
            rel_corr_err(N, 3, 1, 28) = 0.011D0
            rel_corr_err(N, 5, 1, 28) = 0.027D0
            rel_corr_err(N, :, 1, 29) = 0.100D0
            rel_corr_err(N, :, 1, 30) = 0.038D0
            rel_corr_err(N, :, 1, 31) = 0.065D0
            rel_corr_err(N, :, 1, 32) = 0.021D0
            rel_corr_err(N, :, 1, 33) = 0.023D0
            rel_corr_err(N, :, 1, 34) = 0.010D0
            rel_corr_err(N, 2, 1, 36) = 0.500D0 ! n.B.: not official, add theory error for VBF
            rel_corr_err(N, 2, 1, 37) = 0.250D0 ! n.B.: not official, trying migration between VBF categories
        case (7084107) ! ATLAS H->gaga VH hadronic
            rel_corr_err(N, 1, 1, 17) = 0.20D0
            rel_corr_err(N, 1, 1, 27) = 0.028D0
            rel_corr_err(N, 2, 1, 27) = 0.041D0
            rel_corr_err(N, 3, 1, 27) = 0.025D0
            rel_corr_err(N, 4, 1, 27) = 0.025D0
            rel_corr_err(N, 5, 1, 27) = 0.095D0
            rel_corr_err(N, 4, 1, 28) = 0.001D0
            rel_corr_err(N, 5, 1, 28) = 0.007D0
            rel_corr_err(N, :, 1, 29) = 0.110D0
            rel_corr_err(N, :, 1, 30) = 0.040D0
            rel_corr_err(N, :, 1, 31) = 0.072D0
            rel_corr_err(N, :, 1, 32) = 0.016D0
            rel_corr_err(N, :, 1, 33) = 0.038D0
            rel_corr_err(N, :, 1, 34) = 0.041D0
            rel_corr_err(N, 1, 1, 35) = -0.15D0 ! n.b.: not official, guess of migration due to b-tagging uncertainty
            rel_corr_err(N, 3, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
            rel_corr_err(N, 4, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
        case (7084108) ! ATLAS H->gaga VH ETmiss
            rel_corr_err(N, 1, 1, 27) = 0.026D0
            rel_corr_err(N, 2, 1, 27) = 0.090D0
            rel_corr_err(N, 3, 1, 27) = 0.002D0
            rel_corr_err(N, 4, 1, 27) = 0.002D0
            rel_corr_err(N, 5, 1, 27) = 0.012D0
            rel_corr_err(N, 1, 1, 28) = -0.35D0
            rel_corr_err(N, 2, 1, 28) = -0.35D0
            rel_corr_err(N, 3, 1, 28) = -0.013D0
            rel_corr_err(N, 4, 1, 28) = -0.009D0
            rel_corr_err(N, 5, 1, 28) = -0.011D0
            rel_corr_err(N, :, 1, 29) = 0.110D0
            rel_corr_err(N, :, 1, 30) = 0.036D0
            rel_corr_err(N, :, 1, 31) = 0.074D0
            rel_corr_err(N, :, 1, 32) = 0.017D0
            rel_corr_err(N, :, 1, 33) = 0.038D0
            rel_corr_err(N, :, 1, 34) = 0.041D0
            rel_corr_err(N, 3, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
            rel_corr_err(N, 4, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
        case (7084109) ! ATLAS H->gaga VH one-lepton
            rel_corr_err(N, 1, 1, 27) = 0.049D0
            rel_corr_err(N, 2, 1, 27) = 0.062D0
            rel_corr_err(N, 3, 1, 27) = 0.005D0
            rel_corr_err(N, 4, 1, 27) = 0.005D0
            rel_corr_err(N, 5, 1, 27) = 0.028D0
            rel_corr_err(N, 1, 1, 28) = 0.045D0
            rel_corr_err(N, 2, 1, 28) = 0.045D0
            rel_corr_err(N, 3, 1, 28) = 0.004D0
            rel_corr_err(N, 4, 1, 28) = 0.040D0
            rel_corr_err(N, 5, 1, 28) = 0.006D0
            rel_corr_err(N, :, 1, 29) = 0.098D0
            rel_corr_err(N, :, 1, 30) = 0.028D0
            rel_corr_err(N, :, 1, 31) = 0.063D0
            rel_corr_err(N, :, 1, 32) = 0.021D0
            rel_corr_err(N, :, 1, 33) = 0.038D0
            rel_corr_err(N, :, 1, 34) = 0.041D0
            rel_corr_err(N, 3, 1, 35) = -0.05D0 ! n.b.: not official, guess of migration due to b-tagging uncertainty
            rel_corr_err(N, 3, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
            rel_corr_err(N, 4, 1, 38) = 0.50D0  !n.b.: not official, added theory error on VH
        case (7084110) ! ATLAS H->gaga ttH hadronic
            rel_corr_err(N, 1, 1, 19) = 0.60D0
            rel_corr_err(N, 5, 1, 20) = 0.11D0
            rel_corr_err(N, 1, 1, 25) = 0.50D0
            rel_corr_err(N, 1, 1, 26) = 1.00D0
            rel_corr_err(N, 1, 2, 26) = 1.00D0
            rel_corr_err(N, 1, 3, 26) = 1.00D0
            rel_corr_err(N, 1, 1, 27) = 0.11D0
            rel_corr_err(N, 2, 1, 27) = 0.21D0
            rel_corr_err(N, 3, 1, 27) = 0.22D0
            rel_corr_err(N, 4, 1, 27) = 0.22D0
            rel_corr_err(N, 5, 1, 27) = 0.073D0
            rel_corr_err(N, :, 1, 29) = 0.096D0
            rel_corr_err(N, :, 1, 30) = 0.036D0
            rel_corr_err(N, :, 1, 31) = 0.063D0
            rel_corr_err(N, :, 1, 32) = 0.019D0
            rel_corr_err(N, :, 1, 33) = 0.038D0
            rel_corr_err(N, :, 1, 34) = 0.041D0
            rel_corr_err(N, 1, 1, 35) = 0.30D0
        case (7084111) ! ATLAS H->gaga ttH leptonic
            rel_corr_err(N, 5, 1, 20) = 0.03D0
            rel_corr_err(N, 1, 1, 25) = 0.50D0
            rel_corr_err(N, 1, 1, 26) = 1.00D0
            rel_corr_err(N, 1, 2, 26) = 1.00D0
            rel_corr_err(N, 1, 3, 26) = 1.00D0
            rel_corr_err(N, 1, 1, 27) = 0.37D0
            rel_corr_err(N, 2, 1, 27) = 0.077D0
            rel_corr_err(N, 3, 1, 27) = 0.074D0
            rel_corr_err(N, 4, 1, 27) = 0.074D0
            rel_corr_err(N, 5, 1, 27) = 0.005D0
            rel_corr_err(N, 1, 1, 28) = 0.019D0
            rel_corr_err(N, 2, 1, 28) = 0.019D0
            rel_corr_err(N, 3, 1, 28) = 0.010D0
            rel_corr_err(N, 4, 1, 28) = 0.030D0
            rel_corr_err(N, 5, 1, 28) = 0.001D0
            rel_corr_err(N, :, 1, 29) = 0.095D0
            rel_corr_err(N, :, 1, 30) = 0.034D0
            rel_corr_err(N, :, 1, 31) = 0.062D0
            rel_corr_err(N, :, 1, 32) = 0.021D0
            rel_corr_err(N, :, 1, 33) = 0.038D0
            rel_corr_err(N, :, 1, 34) = 0.041D0
            rel_corr_err(N, 3, 1, 35) = 0.07D0
! LHC 13 TeV results
! CMS H-gaga, CMS-16-020:
! 51: UE and parton shower, jet energy scale/smearing: VBF-VBF Tag migration, 7% and 4-15%, respectively.
! 52: UE and parton shower, jet energy scale/smearing: VBF-untagged Tag migration, 9% and 4-15%, respectively.
! 53: Event migration untagged 0-1
! 54: Event migration untagged 1-2
! 55: Event migration untagged 2-3
! n.b: QCD scale uncertainty and energy scale/resolution, all categories, ~5-10% and ~6%, respectively. -> makes fit worse!
! 56: ggF contamination in VBF, ttH categories, ~39%
! 57: ggF contamination: VBF-VBF Tag migration, ~10%
        case (1602011) ! untagged 0
!     rel_corr_err(N,1,1,52)=  -0.24D0
            rel_corr_err(N, 2, 1, 52) = +0.10D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, :, 1, 53) = +0.10D0
        case (1602012) ! untagged 1
!     rel_corr_err(N,1,1,52)=  -0.24D0
            rel_corr_err(N, 2, 1, 52) = +0.05D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, :, 1, 53) = -0.10D0
            rel_corr_err(N, :, 1, 54) = +0.10D0
        case (1602013) ! untagged 2
!     rel_corr_err(N,1,1,52)=  -0.24D0
            rel_corr_err(N, 2, 1, 52) = +0.05D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, :, 1, 54) = -0.10D0
            rel_corr_err(N, :, 1, 55) = +0.10D0
        case (1602014) ! untagged 3
!     rel_corr_err(N,1,1,52)=  -0.24D0
            rel_corr_err(N, 2, 1, 52) = +0.05D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, :, 1, 55) = -0.10D0
        case (1602015) ! VBF tag 0
            rel_corr_err(N, :, 1, 51) = +0.10D0
!     rel_corr_err(N,2,1,51)=  -0.22D0
            rel_corr_err(N, 2, 1, 52) = -0.10D0
!     rel_corr_err(N,2,1,52)=  -0.24D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, 1, 1, 56) = +0.39D0
            rel_corr_err(N, 1, 1, 57) = +0.10D0
        case (1602016) ! VBF tag 1
            rel_corr_err(N, :, 1, 51) = -0.10D0
!     rel_corr_err(N,2,1,51)=  +0.22D0
!    rel_corr_err(N,1,1,52)=  +0.24D0
            rel_corr_err(N, :, 1, 52) = -0.15D0
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, 1, 1, 56) = +0.39D0
            rel_corr_err(N, 1, 1, 57) = -0.10D0
        case (1602017) ! TTH tag hadr
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, 1, 1, 56) = +0.39D0
        case (1602018) ! TTH tag lept
!     rel_corr_err(N,:,1,53)=  +0.10D0
            rel_corr_err(N, 1, 1, 56) = +0.39D0
        case default

        end select

    end subroutine fill_rel_corr_err
!------------------------------------------------
    subroutine get_expt_syst_corr_for_peaks(value, peak1, mu1, peak2, mu2, model)
!------------------------------------------------
        use usefulbits_hs, only: mupeak, print_dble_matrix
        implicit none

        type(mupeak), intent(in) :: peak1, peak2
        integer, intent(in) :: model
        double precision, intent(in) :: mu1, mu2 ! observed mu
        double precision, intent(out) :: value

        integer :: p1, d1, p2, d2
        integer :: i, j, k

        call fill_rel_corr_err(peak1%id, 1)
        call fill_rel_corr_err(peak2%id, 2)

!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) then
!  write(*,*)'#-------------- ',peak1%id,' --------------#'
!  do k=1,Nprod
!   write(*,*) rel_corr_err(1,k,1,:)
!  enddo
!  write(*,*)'#-------------- ',peak2%id,' --------------#'
!  do k=1,Nprod
!   write(*,*) rel_corr_err(2,k,1,:)
!  enddo
!  endif

        value = 0.0D0

        do i = lbound(peak1%channel_p_id, dim=1), ubound(peak1%channel_p_id, dim=1)
            do j = lbound(peak2%channel_p_id, dim=1), ubound(peak2%channel_p_id, dim=1)

                p1 = peak1%channel_p_id(i)
                d1 = peak1%channel_d_id(i)
                p2 = peak2%channel_p_id(j)
                d2 = peak2%channel_d_id(j)

!     id1 = peak1%channel_id(i)
!     p1 = int((id1-modulo(id1,10))/dble(10))
!     d1 = modulo(id1,10)
!     id2 = peak2%channel_id(j)
!     p2 = int((id2-modulo(id2,10))/dble(10))
!     d2 = modulo(id2,10)

!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) then
!   write(*,*) id1, p1, d1, id2, p2, d2
!   write(*,*) value
!  endif

                do k = 1, Nsyst
                    if (model .eq. 1) then
                        if (scaletype(k) .eq. 1) then
                            value = value + &
                                    peak1%channel_w_model(i)*rel_corr_err(1, p1, d1, k)*peak1%total_mu* &
                                    peak2%channel_w_model(j)*rel_corr_err(2, p2, d2, k)*peak2%total_mu
                        elseif (scaletype(k) .eq. 0) then
                            value = value + &
                                    peak1%channel_w(i)*rel_corr_err(1, p1, d1, k)*mu1* &
                                    peak2%channel_w(j)*rel_corr_err(2, p2, d2, k)*mu2
                        else
                            write (*, *) "WARNING in get_expt_syst_corr_for peaks: Unknown scaletype of ", k
                        endif
                    else
                        value = value + &
                                peak1%channel_w(i)*rel_corr_err(1, p1, d1, k)*mu1* &
                                peak2%channel_w(j)*rel_corr_err(2, p2, d2, k)*mu2
                    endif
                enddo
            enddo
        enddo

!  if(abs(value).ge.0.0000001D0) then
!   write(*,*) "Non-zero correlated systematics:", peak1%id, peak2%id, mu1, mu2, value
!   write(*,*) "1st weights: ",peak1%channel_w_model
!   do k=1,Nprod
!    write(*,*) rel_corr_err(1,k,1,:)
!   enddo
!   write(*,*) "2nd weights: ",peak2%channel_w_model
!   do k=1,Nprod
!    write(*,*) rel_corr_err(2,k,1,:)
!   enddo
!
!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001105) write(22,*) value
!  if(peak1%id.eq.13001106.and.peak2%id.eq.13001106) write(23,*) value, peak1%channel_w_model(1)*mu1
!  if(peak1%id.eq.12015103.and.peak2%id.eq.12015103) write(24,*) value, peak1%channel_w_model(1)*mu1
!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) write(25,*) value
!  if(peak1%id.eq.12015103.and.peak2%id.eq.13001106) write(26,*) value

    end subroutine get_expt_syst_corr_for_peaks
!------------------------------------------------
end module
