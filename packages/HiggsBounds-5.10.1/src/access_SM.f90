!> @file
!> External access to SM Higgs predictions.
!! Provides SM cross sections, branching ratios, total width as a function of mass.
!! Values are mostly based on the values provided in the CERN Yellow Report 4.
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to WW^{(*)}\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_HWW(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_HWW

    call testBRSM

    SMBR_HWW = BRSM_HWW(Mh, .True.)

end function SMBR_HWW
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to ZZ^{(*)}\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_HZZ(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_HZZ

    call testBRSM

    SMBR_HZZ = BRSM_HZZ(Mh, .True.)

end function SMBR_HZZ
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to bb\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hbb(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hbb

    call testBRSM

    SMBR_Hbb = BRSM_Hbb(Mh, .True.)

end function SMBR_Hbb
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to \tau\tau\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Htautau(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Htautau

    call testBRSM

    SMBR_Htautau = BRSM_Htautau(Mh, .True.)

end function SMBR_Htautau
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to \gamma\gamma\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hgamgam(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hgamgam

    call testBRSM

    SMBR_Hgamgam = BRSM_Hgaga(Mh, .True.)

end function SMBR_Hgamgam
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to gg\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hgg(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hgg

    call testBRSM

    SMBR_Hgg = BRSM_Hgg(Mh, .True.)
end function SMBR_Hgg
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to tt\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Htoptop(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Htoptop

    call testBRSM

    SMBR_Htoptop = BRSM_Htoptop(Mh, .True.)
end function SMBR_Htoptop
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to cc\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hcc(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hcc

    call testBRSM

    SMBR_Hcc = BRSM_Hcc(Mh, .True.)
end function SMBR_Hcc
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to ss\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hss(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hss

    call testBRSM

    SMBR_Hss = BRSM_Hss(Mh, .True.)
end function SMBR_Hss
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to \mu\mu\f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_Hmumu(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_Hmumu

    call testBRSM

    SMBR_Hmumu = BRSM_Hmumu(Mh, .True.)
end function SMBR_Hmumu
!*********************************
!> SM Higgs branching ratio for the decay \f$ H\to Z\gamma \f$ as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMBR_HZgam(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMBR_HZgam

    call testBRSM

    SMBR_HZgam = BRSM_HZga(Mh, .True.)
end function SMBR_HZgam
!*********************************
!> SM Higgs total decay width (in GeV) as a function of mass.
!! Prediction is taken from the CERN Yellow Report 4 by the LHC Higgs Cross Section
!! Working Group.
!! @param Mh mass of the SM Higgs boson (allowed range \f$[0,1000]~\mathrm{GeV}\f$)
function SMGamma_h(Mh)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMGamma_h

    call testBRSM

    SMGamma_h = BRSM_GammaTot(Mh, .True.)

end function SMGamma_h
!*********************************
!> Decay width (in GeV) of the top quark into a W-boson and a b-quark, as a function of mass.
!! Numbers and equation read from http://pdg.lbl.gov/2018/reviews/rpp2018-rev-top-quark.pdf, Eq. (67.1).
!! @param mt mass of the top quark
function SMGamma_tWpb(mt)
    use theory_BRfunctions
    implicit none
    double precision, intent(in) :: mt
    double precision :: SMGamma_tWpb

    call testBRSM

    SMGamma_tWpb = BRSM_Gamma_tWpb(mt, .True.)

end function SMGamma_tWpb
!*********************************
!> SM cross section (in pb) of \f$HW\f$ production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_HW(Mh)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_HW

    call testBRSM

    SMCS_tev_HW = 1.0D-3*XS_tev_HW_SM(Mh, .True.)
end function SMCS_tev_HW
!*********************************
!> SM cross section (in pb) of \f$HZ\f$ production at the Tevatron, as a function of mass.
!! This is the reference function used in HiggsBounds-4. It is now superceded by
!! \ref smcs_tev_hz.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_HZ_HB4(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_HZ_HB4

    call testBRSM

    SMCS_tev_HZ_HB4 = 1.0D-3*XS_tev_HZ_SM(Mh, .True.)
end function SMCS_tev_HZ_HB4
!*********************************
!> SM cross section (in pb) of inclusive \f$HZ\f$ production at the Tevatron, as a function of mass.
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_tev_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_HZ

    call testBRSM

    SMCS_tev_HZ = ZH_cpmix_nnlo_ggqqbb(Mh, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_tev_HZ
!*********************************
!> SM cross section (in pb) of \f$ q\bar{q} \to HZ \f$ production (including \f$b\bar{b} \to HZ\f$) at the Tevatron, as a function of mass.
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!! @param Mh mass of the SM Higgs boson (mass range \f$ [10,3000]~\mathrm{GeV} \f$)
function SMCS_tev_qq_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_qq_HZ

    call testBRSM

    SMCS_tev_qq_HZ = ZH_cpmix_nnlo_qqbb(Mh, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_tev_qq_HZ
!*********************************
!> SM cross section (in pb) of \f$gg \to HZ\f$ production at the Tevatron, as a function of mass.
!! Prediction is evaluated from fits to `VH@NNLO-2` data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_tev_gg_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_gg_HZ

    call testBRSM

    SMCS_tev_gg_HZ = ZH_cpmix_nnlo_gg(Mh, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_tev_gg_HZ
!*********************************
!> SM cross section (in pb) of \f$b\bar{b} \to HZ\f$ production at the Tevatron, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_tev_bb_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_bb_HZ

    call testBRSM

    SMCS_tev_bb_HZ = ZH_cpmix_nnlo_bb(Mh, 'TEV  ', 1.0D0, 0.0D0, .True.)
end function SMCS_tev_bb_HZ
!*********************************
!*********************************
!> SM cross section (in pb) of \f$gg\to H\f$ production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_gg_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_gg_H

    call testBRSM

!  SMCS_tev_gg_H=1.0D-3*XS_tev_gg_H_SM(Mh) !use this for outdated CS
    SMCS_tev_gg_H = 1.0D-3*XS_tev_gg_H_SM(Mh, .True.) !use this for updated (31/03/2011) CS
end function SMCS_tev_gg_H
!*********************************
!> SM cross section (in pb) of \f$b\bar{b}\to H\f$ production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_bb_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_bb_H

    call testBRSM

    SMCS_tev_bb_H = 1.0D-3*XS_tev_bb_H_SM(Mh, .True.)
end function SMCS_tev_bb_H
!*********************************
!> SM cross section (in pb) of \f$Hqq\f$ (VBF) production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_vbf_H(Mh)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_vbf_H

    call testBRSM

    SMCS_tev_vbf_H = 1.0D-3*XS_tev_vbf_SM(Mh, .True.)
end function SMCS_tev_vbf_H
!*********************************
!> SM cross section (in pb) of \f$bg \to Hb\f$ production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_bg_Hb(Mh)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_bg_Hb

    call testBRSM

    SMCS_tev_bg_Hb = 1.0D-3*XS_tev_bg_Hb_SM(Mh, .True.)
end function SMCS_tev_bg_Hb
!*********************************
!> SM cross section (in pb) of \f$Ht\bar{t}\f$ production at the Tevatron, as a function of mass.
!! @param Mh mass of the SM Higgs boson (fitted mass range \f$[60,360]~\mathrm{GeV}\f$)
function SMCS_tev_ttH(Mh)
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_ttH

    call testBRSM

    SMCS_tev_ttH = 1.0D-3*XS_tev_ttH_SM(Mh, .True.)
end function SMCS_tev_ttH
!*********************************
!> SM cross section (in pb) of \f$HW\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_HW(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_HW

    call testBRSM

    SMCS_lhc7_HW = XS_lhc7_HW_SM(Mh, .True.)
end function SMCS_lhc7_HW
!*********************************
! > SM cross section (in pb) of \f$HZ\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
! ! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
! function SMCS_lhc7_HZ(Mh)
! note: in pb
!  use theory_XS_SM_functions
!   implicit none
!   double precision,intent(in) :: Mh
!   double precision :: SMCS_lhc7_HZ
!
!   call testBRSM
!
!   SMCS_lhc7_HZ=XS_lhc7_HZ_SM(Mh,.True.)
! end function SMCS_lhc7_HZ
!*********************************
!> SM cross section (in pb) of inclusive \f$HZ\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_HZ

    call testBRSM

    SMCS_lhc7_HZ = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc7_HZ
!*********************************
!> SM cross section (in pb) of \f$ q\bar{q} \to HZ \f$ production (including \f$b\bar{b} \to HZ\f$) at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$ [10,3000]~\mathrm{GeV} \f$)
function SMCS_lhc7_qq_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_qq_HZ

    call testBRSM

    SMCS_lhc7_qq_HZ = ZH_cpmix_nnlo_qqbb(Mh, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc7_qq_HZ
!*********************************
!> SM cross section (in pb) of \f$gg \to HZ\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_gg_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_gg_HZ

    call testBRSM

    SMCS_lhc7_gg_HZ = ZH_cpmix_nnlo_gg(Mh, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc7_gg_HZ
!*********************************
!> SM cross section (in pb) of \f$b\bar{b} \to HZ\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_bb_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_bb_HZ

    call testBRSM

    SMCS_lhc7_bb_HZ = ZH_cpmix_nnlo_bb(Mh, 'LHC7 ', 1.0D0, 0.0D0, .True.)
end function SMCS_lhc7_bb_HZ
!*********************************
!> SM cross section (in pb) of \f$gg\to H\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_gg_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_gg_H

    call testBRSM

    SMCS_lhc7_gg_H = XS_lhc7_gg_H_SM(Mh, .True.)
end function SMCS_lhc7_gg_H
!*********************************
!> SM cross section (in pb) of \f$b\bar{b}\to H\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_bb_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_bb_H

    call testBRSM

    SMCS_lhc7_bb_H = XS_lhc7_bb_H_SM(Mh, .True.)
end function SMCS_lhc7_bb_H
!*********************************
!> SM cross section (in pb) of \f$Hqq\f$ (VBF) production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_vbf_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_vbf_H

    call testBRSM

    SMCS_lhc7_vbf_H = XS_lhc7_vbf_SM(Mh, .True.)
end function SMCS_lhc7_vbf_H
!*********************************
!> SM cross section (in pb) of \f$Ht\bar{t}\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_ttH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_ttH

    call testBRSM

    SMCS_lhc7_ttH = XS_lhc7_ttH_SM(Mh, .True.)
end function SMCS_lhc7_ttH
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (s-channel) at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_tH_schan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_tH_schan

    call testBRSM

    SMCS_lhc7_tH_schan = XS_lhc7_tH_schan_SM(Mh, .True., .True.)

end function SMCS_lhc7_tH_schan
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (t-channel) at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc7_tH_tchan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_tH_tchan

    call testBRSM

    SMCS_lhc7_tH_tchan = XS_lhc7_tH_tchan_SM(Mh, .True., .True.)

end function SMCS_lhc7_tH_tchan
!*********************************
!> SM cross section (in pb) of \f$tWH\f$ production at the LHC with \f$\sqrt{s}=7~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (only available for 125~\mathrm{GeV}\f$)
function SMCS_lhc7_tWH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc7_tWH
    call testBRSM

    SMCS_lhc7_tWH = XS_lhc7_gb_tWH_SM(Mh)
end function SMCS_lhc7_tWH
!*********************************
!*********************************
!> SM cross section (in pb) of \f$HW\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_HW(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_HW

    call testBRSM

    SMCS_lhc8_HW = XS_lhc8_HW_SM(Mh, .True.)
end function SMCS_lhc8_HW
! *********************************
! > SM cross section (in pb) of \f$HZ\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
! ! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
! function SMCS_lhc8_HZ(Mh)
! note: in pb
!  use theory_XS_SM_functions
!   implicit none
!   double precision,intent(in) :: Mh
!   double precision :: SMCS_lhc8_HZ
!
!   call testBRSM
!
!   SMCS_lhc8_HZ=XS_lhc8_HZ_SM(Mh,.True.)
! end function SMCS_lhc8_HZ
!*********************************
!> SM cross section (in pb) of inclusive \f$HZ\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_HZ

    call testBRSM

    SMCS_lhc8_HZ = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc8_HZ
!*********************************
!> SM cross section (in pb) of \f$ q\bar{q} \to HZ \f$ production (including \f$b\bar{b} \to HZ\f$) at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$ [10,3000]~\mathrm{GeV} \f$)
function SMCS_lhc8_qq_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_qq_HZ

    call testBRSM

    SMCS_lhc8_qq_HZ = ZH_cpmix_nnlo_qqbb(Mh, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc8_qq_HZ
!*********************************
!> SM cross section (in pb) of \f$gg \to HZ\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_gg_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_gg_HZ

    call testBRSM

    SMCS_lhc8_gg_HZ = ZH_cpmix_nnlo_gg(Mh, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc8_gg_HZ
!*********************************
!> SM cross section (in pb) of \f$b\bar{b} \to HZ\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_bb_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_bb_HZ

    call testBRSM

    SMCS_lhc8_bb_HZ = ZH_cpmix_nnlo_bb(Mh, 'LHC8 ', 1.0D0, 0.0D0, .True.)
end function SMCS_lhc8_bb_HZ
!*********************************
!> SM cross section (in pb) of \f$gg\to H\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_gg_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_gg_H

    call testBRSM

    SMCS_lhc8_gg_H = XS_lhc8_gg_H_SM(Mh, .True.)
end function SMCS_lhc8_gg_H
!*********************************
!> SM cross section (in pb) of \f$b\bar{b}\to H\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_bb_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_bb_H

    call testBRSM

    SMCS_lhc8_bb_H = XS_lhc8_bb_H_SM(Mh, .True.)
end function SMCS_lhc8_bb_H
!*********************************
!> SM cross section (in pb) of \f$Hqq\f$ (VBF) production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_vbf_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_vbf_H

    call testBRSM

    SMCS_lhc8_vbf_H = XS_lhc8_vbf_SM(Mh, .True.)
end function SMCS_lhc8_vbf_H
!*********************************
!> SM cross section (in pb) of \f$Ht\bar{t}\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_ttH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_ttH

    call testBRSM

    SMCS_lhc8_ttH = XS_lhc8_ttH_SM(Mh, .True.)
end function SMCS_lhc8_ttH
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (s-channel) at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_tH_schan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_tH_schan

    call testBRSM

    SMCS_lhc8_tH_schan = XS_lhc8_tH_schan_SM(Mh, .True., .True.)

end function SMCS_lhc8_tH_schan
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (t-channel) at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc8_tH_tchan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_tH_tchan

    call testBRSM

    SMCS_lhc8_tH_tchan = XS_lhc8_tH_tchan_SM(Mh, .True., .True.)

end function SMCS_lhc8_tH_tchan
!*********************************
!> SM cross section (in pb) of \f$tWH\f$ production at the LHC with \f$\sqrt{s}=8~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (only available for 125~\mathrm{GeV}\f$)
function SMCS_lhc8_tWH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc8_tWH
    call testBRSM

    SMCS_lhc8_tWH = XS_lhc8_gb_tWH_SM(Mh)
end function SMCS_lhc8_tWH
!*********************************
!> SM cross section (in pb) of \f$HW\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_HW(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_HW

    call testBRSM

    SMCS_lhc13_HW = XS_lhc13_HW_SM(Mh, .True.)
end function SMCS_lhc13_HW
!*********************************
! > SM cross section (in pb) of \f$HZ\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
! ! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
! function SMCS_lhc13_HZ(Mh)
! note: in pb
!  use theory_XS_SM_functions
!   implicit none
!   double precision,intent(in) :: Mh
!   double precision :: SMCS_lhc13_HZ
!
!   call testBRSM
!
!   SMCS_lhc13_HZ=XS_lhc13_HZ_SM(Mh,.True.)
! end function SMCS_lhc13_HZ
!*********************************
!> SM cross section (in pb) of inclusive \f$HZ\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_HZ

    call testBRSM

    SMCS_lhc13_HZ = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc13_HZ
!*********************************
!> SM cross section (in pb) of \f$q\bar{q} \to HZ\f$ production (including \f$b\bar{b} \to HZ\f$) at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_qq_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_qq_HZ

    call testBRSM

    SMCS_lhc13_qq_HZ = ZH_cpmix_nnlo_qqbb(Mh, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc13_qq_HZ
!*********************************
!> SM cross section (in pb) of \f$gg \to HZ\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_gg_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_gg_HZ

    call testBRSM

    SMCS_lhc13_gg_HZ = ZH_cpmix_nnlo_gg(Mh, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function SMCS_lhc13_gg_HZ
!*********************************
!> SM cross section (in pb) of \f$b\bar{b} \to HZ\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! Prediction is evaluated from fits to`VH@NNLO-2`data.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_bb_HZ(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_bb_HZ

    call testBRSM

    SMCS_lhc13_bb_HZ = ZH_cpmix_nnlo_bb(Mh, 'LHC13', 1.0D0, 0.0D0, .True.)
end function SMCS_lhc13_bb_HZ
!*********************************
!> SM cross section (in pb) of \f$gg\to H\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_gg_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_gg_H

    call testBRSM

    SMCS_lhc13_gg_H = XS_lhc13_gg_H_SM(Mh, .True.)
end function SMCS_lhc13_gg_H
!*********************************
!> SM cross section (in pb) of \f$b\bar{b}\to H\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_bb_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_bb_H

    call testBRSM

    SMCS_lhc13_bb_H = XS_lhc13_bb_H_SM(Mh, .True.)
end function SMCS_lhc13_bb_H
!*********************************
!> SM cross section (in pb) of \f$Hqq\f$ (VBF) production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_vbf_H(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_vbf_H

    call testBRSM

    SMCS_lhc13_vbf_H = XS_lhc13_vbf_SM(Mh, .True.)
end function SMCS_lhc13_vbf_H
!*********************************
!> SM cross section (in pb) of \f$Ht\bar{t}\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_ttH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_ttH

    call testBRSM

    SMCS_lhc13_ttH = XS_lhc13_ttH_SM(Mh, .True.)
end function SMCS_lhc13_ttH
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (s-channel) at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_tH_schan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_tH_schan

    call testBRSM

    SMCS_lhc13_tH_schan = XS_lhc13_tH_schan_SM(Mh, .True., .True.)

end function SMCS_lhc13_tH_schan
!*********************************
!> SM cross section (in pb) of \f$Ht\f$ production (t-channel) at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (mass range \f$[10,3000]~\mathrm{GeV}\f$)
function SMCS_lhc13_tH_tchan(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_tH_tchan

    call testBRSM

    SMCS_lhc13_tH_tchan = XS_lhc13_tH_tchan_SM(Mh, .True., .True.)

end function SMCS_lhc13_tH_tchan
!*********************************
!> SM cross section (in pb) of \f$tWH\f$ production at the LHC with \f$\sqrt{s}=13~\mathrm{TeV}\f$, as a function of mass.
!! @param Mh mass of the SM Higgs boson (only available for 125~\mathrm{GeV}\f$)
function SMCS_lhc13_tWH(Mh)
!note: in pb
    use theory_XS_SM_functions
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_lhc13_tWH
    call testBRSM

    SMCS_lhc13_tWH = XS_lhc13_gb_tWH_SM(Mh)
end function SMCS_lhc13_tWH
!*********************************



subroutine testBRSM()
    use theory_BRfunctions
    implicit none

    if (.not. allocated(BRSM)) then
        write (*, *) 'You can only use this function between'
        write (*, *) 'calling the subroutines'
        write (*, *) 'initialize_HiggsBounds and'
        write (*, *) 'finish_HiggsBounds'
        stop 'error (see standard output for more info)'
    endif
end subroutine testBRSM
!*********************************
! OLD FUNCTIONS
!*********************************
! function SMCS_tev_qq_HZ(Mh)
!   implicit none
!   double precision,intent(in) :: Mh
!   double precision :: SMCS_tev_qq_HZ
!   double precision :: SMCS_tev_HZ
!
!   write(*,*)'Note: function SMCS_tev_qq_HZ has been'
!   write(*,*)'superceded by function SMCS_tev_HZ'
!   SMCS_tev_qq_HZ=SMCS_tev_HZ(Mh)
!
! end function SMCS_tev_qq_HZ
!*********************************
function SMCS_tev_qq_HW(Mh)
    implicit none
    double precision, intent(in) :: Mh
    double precision :: SMCS_tev_qq_HW
    double precision :: SMCS_tev_HW

    write (*, *) 'Note: function SMCS_tev_qq_HW has been'
    write (*, *) 'superceded by function SMCS_tev_HW'
    SMCS_tev_qq_HW = SMCS_tev_HW(Mh)

end function SMCS_tev_qq_HW
!*********************************
