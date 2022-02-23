!> @file
!> External access to rescaled SM Higgs quantities in the effective coupling approximation (aka kappa framework).
!*********************************
!> Cross section (in pb) of inclusive \f$HZ\f$ production, as a function of mass and effective Higgs couplings.
!! Both scalar and pseudoscalar coupling to fermions are taken into account (see manual for more information).
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!!
!! @param Mh mass of the Higgs boson (allowed mass range \f$[10,3000]~\mathrm{GeV}\f$)
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param ghZZ (SM-normalized) effective Higgs coupling to Z bosons
!! @param ghtt_s Scalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghtt_p Pseudoscalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_p Pseudoscalar (SM-normalized) effective Higgs coupling to bottom quarks
function SMCS_effC_HZ(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    integer, intent(in) :: collider
    double precision, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p

    double precision :: SMCS_effC_HZ

    call testBRSM

    SMCS_effC_HZ = ZH_cpmix_nnlo_ggqqbb(Mh, collider_name(collider), ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function SMCS_effC_HZ
!*********************************
!> Cross section (in pb) of \f$gg\to HZ\f$ production, as a function of mass and effective Higgs couplings (aka scale factors).
!! Both scalar and pseudoscalar coupling to fermions are taken into account (see manual for more information).
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!!
!! @param Mh mass of the Higgs boson (allowed mass range \f$[10,3000]~\mathrm{GeV}\f$)
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param ghZZ (SM-normalized) effective Higgs coupling to Z bosons
!! @param ghtt_s Scalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghtt_p Pseudoscalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_p Pseudoscalar (SM-normalized) effective Higgs coupling to bottom quarks
function SMCS_effC_gg_HZ(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    integer, intent(in) :: collider
    double precision, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p
    double precision :: SMCS_effC_gg_HZ

    call testBRSM

    SMCS_effC_gg_HZ = ZH_cpmix_nnlo_gg(Mh, collider_name(collider), ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function SMCS_effC_gg_HZ
!*********************************
!> Cross section (in pb) of \f$q\bar{q}\to HZ\f$ production (including \f$b\bar{b} \to HZ\f$), as a function of mass and effective Higgs couplings (aka scale factors).
!! Both scalar and pseudoscalar coupling to fermions are taken into account (see manual for more information).
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!!
!! @param Mh mass of the Higgs boson (allowed mass range \f$[10,3000]~\mathrm{GeV}\f$)
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param ghZZ (SM-normalized) effective Higgs coupling to Z bosons
!! @param ghtt_s Scalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghtt_p Pseudoscalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_p Pseudoscalar (SM-normalized) effective Higgs coupling to bottom quarks
function SMCS_effC_qq_HZ(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    integer, intent(in) :: collider
    double precision, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p
    double precision :: SMCS_effC_qq_HZ

    call testBRSM

    SMCS_effC_qq_HZ = ZH_cpmix_nnlo_qqbb(Mh, collider_name(collider), ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function SMCS_effC_qq_HZ
!*********************************
!> Cross section (in pb) of \f$b\bar{b}\to HZ\f$ production, as a function of mass and effective Higgs couplings (aka scale factors).
!! Both scalar and pseudoscalar coupling to fermions are taken into account (see manual for more information).
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!!
!! @param Mh mass of the Higgs boson (allowed mass range \f$[10,3000]~\mathrm{GeV}\f$)
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param ghbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghbb_p Pseudoscalar (SM-normalized) effective Higgs coupling to bottom quarks
function SMCS_effC_bb_HZ(Mh, collider, ghbb_s, ghbb_p)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    integer, intent(in) :: collider
    double precision, intent(in) :: ghbb_s, ghbb_p
    double precision :: SMCS_effC_bb_HZ

    call testBRSM

    SMCS_effC_bb_HZ = ZH_cpmix_nnlo_bb(Mh, collider_name(collider), ghbb_s, ghbb_p, .True.)
end function SMCS_effC_bb_HZ
!> Cross section (in pb) of inclusive \f$HZ\f$ production, as a function of mass and effective Higgs couplings.
!! Both scalar and pseudoscalar coupling to fermions are taken into account (see manual for more information).
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!!
!! @param Mh mass of the Higgs boson (allowed mass range \f$[10,3000]~\mathrm{GeV}\f$)
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13    | LHC at 13 TeV|
!! @param ghWW (SM-normalized) effective Higgs coupling to Z bosons
!! @param ghtt_s Scalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
function SMCS_effC_HW(Mh, collider, ghWW, ghtt_s, ghbb_s)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    integer, intent(in) :: collider
    double precision, intent(in) :: ghWW, ghtt_s,ghbb_s

    double precision :: SMCS_effC_HW

    call testBRSM

    SMCS_effC_HW = WH_nnlo(Mh, collider_name(collider), ghWW, ghtt_s, ghbb_s, .True., .False.)
end function SMCS_effC_HW
!> Cross section (in pb) of \f$ p p \to t H\pm\f$ production at the 13TeV LHC.
!!
!! Evaluated on a coefficient grid extracted from the 2HDM results of
!! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMCharged.
!! See references on the website. If \f$ \mathrm{BR}(t\to H_j^+ b)>0 \f$ the relation
!! \f$ \mathrm{BR}(t\to W^+ b) + \mathrm{BR}(t\to H^+_j b)\approx 1\f$ is assumed.
!!
!! @param MHc mass of \f$ H^\pm_j \f$ (allowed mass range \f$[145,2000]~\mathrm{GeV}\f$)
!! @param gHcjt scaling factor \f$ \kappa^{j\pm}_t \f$ of the charged-Higgs top coupling (\f$=1/\tan\beta\f$ in the 2HDM type-II)
!! @param gHcjb scaling factor \f$ \kappa^{j\pm}_b \f$ of the charged-Higgs bottom coupling (\f$=\tan\beta\f$ in the 2HDM type-II)
!! @param BR_tHpjb \f$ \mathrm{BR}(t\to H^+_j b) \f$
function HCCS_tHc(MHc, gHcjt, gHcjb, BR_tHpjb)
    use theory_XS_SM_functions, only: tHc_cxn
    implicit none
    double precision, intent(in) :: MHc, gHcjt, gHcjb, BR_tHpjb
    double precision :: HCCS_tHc

    HCCS_tHc = tHc_cxn(Mhc, gHcjt, gHcjb, BR_tHpjb)
end function HCCS_tHc
