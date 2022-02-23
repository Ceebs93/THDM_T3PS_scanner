! ---- Branching ratios ----
function SMBR_HWW_C(Mh) &
    bind(C, name='SMBR_HWW')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_HWW_C

    call testBRSM
    SMBR_HWW_C = BRSM_HWW(Mh, .True.)
end function

function SMBR_HZZ_C(Mh) &
    bind(C, name='SMBR_HZZ')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_HZZ_C

    call testBRSM
    SMBR_HZZ_C = BRSM_HZZ(Mh, .True.)
end function

function SMBR_Hbb_C(Mh) &
    bind(C, name='SMBR_Hbb')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hbb_C

    call testBRSM
    SMBR_Hbb_C = BRSM_Hbb(Mh, .True.)
end function

function SMBR_Htautau_C(Mh) &
    bind(C, name='SMBR_Htautau')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Htautau_C

    call testBRSM
    SMBR_Htautau_C = BRSM_Htautau(Mh, .True.)
end function

function SMBR_Hgamgam_C(Mh) &
    bind(C, name='SMBR_Hgamgam')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hgamgam_C

    call testBRSM
    SMBR_Hgamgam_C = BRSM_Hgaga(Mh, .True.)
end function

function SMBR_Hgg_C(Mh) &
    bind(C, name='SMBR_Hgg')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hgg_C

    call testBRSM
    SMBR_Hgg_C = BRSM_Hgg(Mh, .True.)
end function

function SMBR_Htoptop_C(Mh) &
    bind(C, name='SMBR_Htoptop')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Htoptop_C

    call testBRSM
    SMBR_Htoptop_C = BRSM_Htoptop(Mh, .True.)
end function

function SMBR_Hcc_C(Mh) &
    bind(C, name='SMBR_Hcc')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hcc_C

    call testBRSM
    SMBR_Hcc_C = BRSM_Hcc(Mh, .True.)
end function

function SMBR_Hss_C(Mh) &
    bind(C, name='SMBR_Hss')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hss_C

    call testBRSM
    SMBR_Hss_C = BRSM_Hss(Mh, .True.)
end function

function SMBR_Hmumu_C(Mh) &
    bind(C, name='SMBR_Hmumu')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_Hmumu_C

    call testBRSM
    SMBR_Hmumu_C = BRSM_Hmumu(Mh, .True.)
end function

function SMBR_HZgam_C(Mh) &
    bind(C, name='SMBR_HZgam')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMBR_HZgam_C

    call testBRSM
    SMBR_HZgam_C = BRSM_HZga(Mh, .True.)
end function

! ---- Widths ----
function SMGamma_h_C(Mh) &
    bind(C, name='SMGamma_H')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMGamma_h_C

    call testBRSM
    SMGamma_h_C = BRSM_GammaTot(Mh, .True.)
end function

function SMGamma_tWpb_C(mt) &
    bind(C, name='SMGamma_tWpb')
    use, intrinsic :: iso_c_binding
    use theory_BRfunctions

    implicit none
    real(kind=c_double), value, intent(in) :: mt
    real(kind=c_double) :: SMGamma_tWpb_C

    call testBRSM
    SMGamma_tWpb_C = BRSM_Gamma_tWpb(mt, .True.)
end function

! ---- Cross sections - LEP ----
function SMCS_tev_HW_C(Mh) &
    bind(C, name='SMCS_tev_HW')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_HW_C

    call testBRSM
    SMCS_tev_HW_C = 1.0D-3*XS_tev_HW_SM(Mh, .True.)
end function

function SMCS_tev_HZ_C(Mh) &
    bind(C, name='SMCS_tev_HZ')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_HZ_C

    call testBRSM
    SMCS_tev_HZ_C = ZH_cpmix_nnlo_ggqqbb(Mh, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function

function SMCS_tev_gg_H_C(Mh) &
    bind(C, name='SMCS_tev_gg_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_gg_H_C

    call testBRSM
    SMCS_tev_gg_H_C = 1.0D-3*XS_tev_gg_H_SM(Mh, .True.)
end function

function SMCS_tev_bb_H_C(Mh) &
    bind(C, name='SMCS_tev_bb_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_bb_H_C

    call testBRSM
    SMCS_tev_bb_H_C = 1.0D-3*XS_tev_bb_H_SM(Mh, .True.)
end function

function SMCS_tev_vbf_H_C(Mh) &
    bind(C, name='SMCS_tev_vbf_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_vbf_H_C

    call testBRSM
    SMCS_tev_vbf_H_C = 1.0D-3*XS_tev_vbf_SM(Mh, .True.)
end function

function SMCS_tev_bg_Hb_C(Mh) &
    bind(C, name='SMCS_tev_bg_Hb')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_bg_Hb_C

    call testBRSM
    SMCS_tev_bg_Hb_C = 1.0D-3*XS_tev_bg_Hb_SM(Mh, .True.)
end function

function SMCS_tev_ttH_C(Mh) &
    bind(C, name='SMCS_tev_ttH')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_tev_ttH_C

    call testBRSM
    SMCS_tev_ttH_C = 1.0D-3*XS_tev_ttH_SM(Mh, .True.)
end function

! ---- Cross sections - LHC 7TeV ----
function SMCS_lhc7_HW_C(Mh) &
    bind(C, name='SMCS_lhc7_HW')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_HW_C

    call testBRSM
    SMCS_lhc7_HW_C = XS_lhc7_HW_SM(Mh, .True.)
end function

function SMCS_lhc7_HZ_C(Mh) &
    bind(C, name='SMCS_lhc7_HZ')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_HZ_C

    call testBRSM
    SMCS_lhc7_HZ_C = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function

function SMCS_lhc7_gg_H_C(Mh) &
    bind(C, name='SMCS_lhc7_gg_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_gg_H_C

    call testBRSM
    SMCS_lhc7_gg_H_C = XS_lhc7_gg_H_SM(Mh, .True.)
end function

function SMCS_lhc7_bb_H_C(Mh) &
    bind(C, name='SMCS_lhc7_bb_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_bb_H_C

    call testBRSM
    SMCS_lhc7_bb_H_C = XS_lhc7_bb_H_SM(Mh, .True.)
end function

function SMCS_lhc7_vbf_H_C(Mh) &
    bind(C, name='SMCS_lhc7_vbf_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_vbf_H_C

    call testBRSM
    SMCS_lhc7_vbf_H_C = XS_lhc7_vbf_SM(Mh, .True.)
end function

function SMCS_lhc7_ttH_C(Mh) &
    bind(C, name='SMCS_lhc7_ttH')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc7_ttH_C

    call testBRSM
    SMCS_lhc7_ttH_C = XS_lhc7_ttH_SM(Mh, .True.)
end function

! ---- Cross sections - LHC 8TeV ----
function SMCS_lhc8_HW_C(Mh) &
    bind(C, name='SMCS_lhc8_HW')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_HW_C

    call testBRSM
    SMCS_lhc8_HW_C = XS_lhc8_HW_SM(Mh, .True.)
end function

function SMCS_lhc8_HZ_C(Mh) &
    bind(C, name='SMCS_lhc8_HZ')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_HZ_C

    call testBRSM
    SMCS_lhc8_HZ_C = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function

function SMCS_lhc8_gg_H_C(Mh) &
    bind(C, name='SMCS_lhc8_gg_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_gg_H_C

    call testBRSM
    SMCS_lhc8_gg_H_C = XS_lhc8_gg_H_SM(Mh, .True.)
end function

function SMCS_lhc8_bb_H_C(Mh) &
    bind(C, name='SMCS_lhc8_bb_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_bb_H_C

    call testBRSM
    SMCS_lhc8_bb_H_C = XS_lhc8_bb_H_SM(Mh, .True.)
end function

function SMCS_lhc8_vbf_H_C(Mh) &
    bind(C, name='SMCS_lhc8_vbf_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_vbf_H_C

    call testBRSM
    SMCS_lhc8_vbf_H_C = XS_lhc8_vbf_SM(Mh, .True.)
end function

function SMCS_lhc8_ttH_C(Mh) &
    bind(C, name='SMCS_lhc8_ttH')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc8_ttH_C

    call testBRSM
    SMCS_lhc8_ttH_C = XS_lhc8_ttH_SM(Mh, .True.)
end function

! ---- Cross sections - LHC 13TeV ----
function SMCS_lhc13_HW_C(Mh) &
    bind(C, name='SMCS_lhc13_HW')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_HW_C

    call testBRSM
    SMCS_lhc13_HW_C = XS_lhc13_HW_SM(Mh, .True.)
end function

function SMCS_lhc13_HZ_C(Mh) &
    bind(C, name='SMCS_lhc13_HZ')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_HZ_C

    call testBRSM
    SMCS_lhc13_HZ_C = ZH_cpmix_nnlo_ggqqbb(Mh, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
end function

function SMCS_lhc13_gg_H_C(Mh) &
    bind(C, name='SMCS_lhc13_gg_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_gg_H_C

    call testBRSM
    SMCS_lhc13_gg_H_C = XS_lhc13_gg_H_SM(Mh, .True.)
end function

function SMCS_lhc13_bb_H_C(Mh) &
    bind(C, name='SMCS_lhc13_bb_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_bb_H_C

    call testBRSM
    SMCS_lhc13_bb_H_C = XS_lhc13_bb_H_SM(Mh, .True.)
end function

function SMCS_lhc13_vbf_H_C(Mh) &
    bind(C, name='SMCS_lhc13_vbf_H')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_vbf_H_C

    call testBRSM
    SMCS_lhc13_vbf_H_C = XS_lhc13_vbf_SM(Mh, .True.)
end function

function SMCS_lhc13_ttH_C(Mh) &
    bind(C, name='SMCS_lhc13_ttH')
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    real(kind=c_double) :: SMCS_lhc13_ttH_C

    call testBRSM
    SMCS_lhc13_ttH_C = XS_lhc13_ttH_SM(Mh, .True.)
end function
