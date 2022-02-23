subroutine initialize_HiggsBounds_C(nHneut, nHplus, analyses_flag) &
    bind(c, name='initialize_HiggsBounds')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_size_t), intent(in), value :: nHneut, nHplus
    integer(kind=c_int), intent(in), value :: analyses_flag

    call initialize_HiggsBounds_int(int(nHneut), int(nHplus), analyses_flag)
end subroutine

! subroutine HiggsBounds_input_SLHA(infile)

subroutine HiggsBounds_neutral_input_properties_C(Mh, GammaTotal_hj, CP_value) &
    bind(c, name='HiggsBounds_neutral_input_properties')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: Mh(np(Hneut)), GammaTotal_hj(np(Hneut))
    integer(kind=c_int), intent(in) :: CP_value(np(Hneut))

    call HiggsBounds_neutral_input_properties(Mh, GammaTotal_hj, CP_value)
end subroutine

subroutine HiggsBounds_neutral_input_effC_C(ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                                            ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                                            ghjmumu_s, ghjmumu_p, &
                                            ghjtautau_s, ghjtautau_p, &
                                            ghjWW, ghjZZ, ghjZga, &
                                            ghjgaga, ghjgg, ghjhiZ) &
    bind(C, name='HiggsBounds_neutral_input_effC')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: ghjss_s(np(Hneut)), ghjss_p(np(Hneut)), &
                                       ghjcc_s(np(Hneut)), ghjcc_p(np(Hneut)), &
                                       ghjbb_s(np(Hneut)), ghjbb_p(np(Hneut)), &
                                       ghjtt_s(np(Hneut)), ghjtt_p(np(Hneut)), &
                                       ghjmumu_s(np(Hneut)), ghjmumu_p(np(Hneut)), &
                                       ghjtautau_s(np(Hneut)), ghjtautau_p(np(Hneut)), &
                                       ghjWW(np(Hneut)), ghjZZ(np(Hneut)), ghjZga(np(Hneut)), &
                                       ghjgaga(np(Hneut)), ghjgg(np(Hneut)), &
                                       ghjhiZ(np(Hneut), np(Hneut))

    call HiggsBounds_neutral_input_effC(ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
                                        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
                                        ghjmumu_s, ghjmumu_p, &
                                        ghjtautau_s, ghjtautau_p, &
                                        ghjWW, ghjZZ, ghjZga, &
                                        ghjgaga, ghjgg, ghjhiZ)
end subroutine

subroutine HiggsBounds_neutral_input_SMBR_C(BR_hjss, BR_hjcc, BR_hjbb, &
                                            BR_hjtt, BR_hjmumu, &
                                            BR_hjtautau, BR_hjWW, &
                                            BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                            BR_hjgg) &
    bind(C, name='HiggsBounds_neutral_input_SMBR')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: BR_hjss(np(Hneut)), BR_hjcc(np(Hneut)), &
                                       BR_hjbb(np(Hneut)), BR_hjtt(np(Hneut)), &
                                       BR_hjmumu(np(Hneut)), BR_hjtautau(np(Hneut)), &
                                       BR_hjWW(np(Hneut)), BR_hjZZ(np(Hneut)), &
                                       BR_hjZga(np(Hneut)), BR_hjgaga(np(Hneut)), &
                                       BR_hjgg(np(Hneut))

    call HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, &
                                        BR_hjtt, BR_hjmumu, &
                                        BR_hjtautau, BR_hjWW, &
                                        BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                        BR_hjgg)
end subroutine

subroutine HiggsBounds_neutral_input_nonSMBR_C(BR_hjinvisible, BR_hkhjhi_C, BR_hjhiZ_C, &
                                               BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW_C) &
    bind(C, name='HiggsBounds_neutral_input_nonSMBR')

    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    real(kind=c_double), intent(in) :: BR_hjinvisible(np(Hneut)), &
                                       BR_hkhjhi_C(np(Hneut), np(Hneut), np(Hneut)), &
                                       BR_hjhiZ_C(np(Hneut), np(Hneut)), &
                                       BR_hjemu(np(Hneut)), &
                                       BR_hjetau(np(Hneut)), &
                                       BR_hjmutau(np(Hneut))
    real(kind=c_double), intent(in) :: BR_hjHpiW_C(np(Hplus), np(Hneut))
    double precision :: BR_Hkhjhi(np(Hneut), np(Hneut), np(Hneut)), &
        BR_hjhiZ(np(Hneut), np(Hneut)), &
        BR_hjHpiW(np(Hneut), np(Hplus))

    BR_Hkhjhi = RESHAPE(BR_Hkhjhi_C, (/np(Hneut), np(Hneut), np(Hneut)/), ORDER=(/3, 2, 1/))
    BR_hjhiZ = TRANSPOSE(BR_hjhiZ_C)
    BR_hjHpiW = TRANSPOSE(BR_hjHpiW_C)
    call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi, BR_hjhiZ, &
                                           BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW)
end subroutine

subroutine HiggsBounds_neutral_input_effC_firstgen_C( &
    ghjuu_s, ghjuu_p, ghjdd_s, ghjdd_p, ghjee_s, ghjee_p) &
    bind(C, name='HiggsBounds_neutral_input_effC_firstgen')

    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: ghjuu_s(np(Hneut)), ghjuu_p(np(Hneut)), &
                                       ghjdd_s(np(Hneut)), ghjdd_p(np(Hneut)), &
                                       ghjee_s(np(Hneut)), ghjee_p(np(Hneut))
    call HiggsBounds_neutral_input_effC_firstgen(ghjuu_s, ghjuu_p, ghjdd_s, ghjdd_p, ghjee_s, ghjee_p)
end subroutine

subroutine HiggsBounds_neutral_input_effC_FV_C( &
    ghjuc_s, ghjuc_p, ghjut_s, ghjut_p, ghjct_s, ghjct_p, &
    ghjds_s, ghjds_p, ghjdb_s, ghjdb_p, ghjsb_s, ghjsb_p) &
    bind(C, name='HiggsBounds_neutral_input_effC_FV')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: ghjuc_s(np(Hneut)), ghjuc_p(np(Hneut)), &
                                       ghjut_s(np(Hneut)), ghjut_p(np(Hneut)), &
                                       ghjct_s(np(Hneut)), ghjct_p(np(Hneut)), &
                                       ghjds_s(np(Hneut)), ghjds_p(np(Hneut)), &
                                       ghjdb_s(np(Hneut)), ghjdb_p(np(Hneut)), &
                                       ghjsb_s(np(Hneut)), ghjsb_p(np(Hneut))

    call HiggsBounds_neutral_input_effC_FV( &
        ghjuc_s, ghjuc_p, ghjut_s, ghjut_p, ghjct_s, ghjct_p, &
        ghjds_s, ghjds_p, ghjdb_s, ghjdb_p, ghjsb_s, ghjsb_p)

end subroutine

subroutine HiggsBounds_neutral_input_firstgenBR_C(BR_hjuu, BR_hjdd, BR_hjee) &
    bind(C, name='HiggsBounds_neutral_input_firstgenBR')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: BR_hjuu(np(Hneut)), BR_hjdd(np(Hneut)), BR_hjee(np(Hneut))
    call HiggsBounds_neutral_input_firstgenBR(BR_hjuu, BR_hjdd, BR_hjee)
end subroutine

subroutine HiggsBounds_neutral_input_FVBR_C(BR_hjuc, BR_hjds, BR_hjut, BR_hjdb, BR_hjct, BR_hjsb) &
    bind(C, name='HiggsBounds_neutral_input_FVBR')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: BR_hjuc(np(Hneut)), BR_hjds(np(Hneut)), &
                                       BR_hjut(np(Hneut)), BR_hjdb(np(Hneut)), &
                                       BR_hjct(np(Hneut)), BR_hjsb(np(Hneut))
    call HiggsBounds_neutral_input_FVBR(BR_hjuc, BR_hjds, BR_hjut, BR_hjdb, BR_hjct, BR_hjsb)
end subroutine

subroutine HiggsBounds_neutral_input_LEP_C(XS_ee_hjZ_ratio, XS_ee_bbhj_ratio, &
                                           XS_ee_tautauhj_ratio, XS_ee_hjhi_ratio) &
    bind(C, name='HiggsBounds_neutral_input_LEP')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    real(kind=c_double), intent(in) :: XS_ee_hjZ_ratio(np(Hneut)), &
                                       XS_ee_bbhj_ratio(np(Hneut)), XS_ee_tautauhj_ratio(np(Hneut)), &
                                       XS_ee_hjhi_ratio(np(Hneut), np(Hneut))

    call HiggsBounds_neutral_input_LEP(XS_ee_hjZ_ratio, XS_ee_bbhj_ratio, &
                                       XS_ee_tautauhj_ratio, XS_ee_hjhi_ratio)
end subroutine HiggsBounds_neutral_input_LEP_C

subroutine HiggsBounds_neutral_input_hadr_C(collider, CS_hj_ratio, &
                                            CS_gg_hj_ratio, CS_bb_hj_ratio, &
                                            CS_hjW_ratio, CS_hjZ_ratio, &
                                            CS_vbf_ratio, CS_tthj_ratio, &
                                            CS_thj_tchan_ratio, CS_thj_schan_ratio, &
                                            CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, &
                                            CS_tWhj_ratio, &
                                            CS_hjhi) &
    bind(C, name='HiggsBounds_neutral_input_hadr')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut

    implicit none
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), intent(in) :: CS_hj_ratio(np(Hneut)), &
                                       CS_gg_hj_ratio(np(Hneut)), CS_bb_hj_ratio(np(Hneut)), &
                                       CS_hjW_ratio(np(Hneut)), CS_hjZ_ratio(np(Hneut)), &
                                       CS_vbf_ratio(np(Hneut)), CS_tthj_ratio(np(Hneut)), &
                                       CS_thj_tchan_ratio(np(Hneut)), CS_thj_schan_ratio(np(Hneut)), &
                                       CS_qq_hjZ_ratio(np(Hneut)), CS_gg_hjZ_ratio(np(Hneut)), &
                                       CS_tWhj_ratio(np(Hneut)), &
                                       CS_hjhi(np(Hneut), np(Hneut))

    call HiggsBounds_neutral_input_hadr(collider, CS_hj_ratio, &
                                        CS_gg_hj_ratio, CS_bb_hj_ratio, &
                                        CS_hjW_ratio, CS_hjZ_ratio, &
                                        CS_vbf_ratio, CS_tthj_ratio, &
                                        CS_thj_tchan_ratio, CS_thj_schan_ratio, &
                                        CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, &
                                        CS_tWhj_ratio, &
                                        CS_hjhi)
end subroutine

! subroutine HiggsBounds_neutral_input_hadr_channelrates(collider,channelrates)

subroutine HiggsBounds_charged_input_C(Mhplus, GammaTotal_Hpj, &
                                       CS_ee_HpjHmj_ratio, &
                                       BR_tWpb, BR_tHpjb, &
                                       BR_Hpjcs, BR_Hpjcb, BR_Hpjtaunu, BR_Hpjtb, &
                                       BR_HpjWZ, BR_HpjhiW_C) &
    bind(C, name='HiggsBounds_charged_input')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    real(kind=c_double), intent(in) :: Mhplus(np(Hplus)), GammaTotal_Hpj(np(Hplus)), &
                                       CS_ee_HpjHmj_ratio(np(Hplus)), &
                                       BR_tHpjb(np(Hplus)), &
                                       BR_Hpjcs(np(Hplus)), BR_Hpjcb(np(Hplus)), BR_Hpjtaunu(np(Hplus)), &
                                       BR_Hpjtb(np(Hplus)), BR_HpjWZ(np(Hplus))
    real(kind=c_double), intent(in), value :: BR_tWpb
    real(kind=c_double), intent(in) :: BR_HpjhiW_C(np(Hneut), np(Hplus))
    double precision :: BR_HpjhiW(np(Hplus), np(Hneut))

    BR_HpjhiW = TRANSPOSE(BR_HpjhiW_C)
    call HiggsBounds_charged_input(Mhplus, GammaTotal_Hpj, &
                                   CS_ee_HpjHmj_ratio, &
                                   BR_tWpb, BR_tHpjb, &
                                   BR_Hpjcs, BR_Hpjcb, BR_Hpjtaunu, BR_Hpjtb, &
                                   BR_HpjWZ, BR_HpjhiW)
end subroutine

subroutine HiggsBounds_charged_input_exoticBR_C(BR_Hpjud, BR_Hpjus, BR_Hpjcd, BR_Hpjub, BR_Hpjenu, BR_Hpjmunu) &
    bind(C, name='HiggsBounds_charged_input_exoticBR')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hplus

    implicit none
    real(kind=c_double), intent(in) ::BR_Hpjud(np(Hplus)), BR_Hpjus(np(Hplus)), &
                                       BR_Hpjcd(np(Hplus)), BR_Hpjub(np(Hplus)), &
                                       BR_Hpjenu(np(Hplus)), BR_Hpjmunu(np(Hplus))

    call HiggsBounds_charged_input_exoticBR(BR_Hpjud, BR_Hpjus, BR_Hpjcd, BR_Hpjub, BR_Hpjenu, BR_Hpjmunu)
end subroutine

subroutine HiggsBounds_charged_input_hadr_C(collider, CS_Hpmjtb, CS_Hpmjcb, &
                                            CS_Hpmjbjet, CS_Hpmjcjet, CS_Hpmjjetjet, CS_HpmjW, &
                                            CS_HpmjZ, CS_vbf_Hpmj, CS_HpjHmj, CS_Hpmjhi_C) &
    bind(C, name='HiggsBounds_charged_input_hadr')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    real(kind=c_double), intent(in) :: CS_Hpmjtb(np(Hplus)), CS_Hpmjcb(np(Hplus)), &
                                       CS_Hpmjbjet(np(Hplus)), CS_Hpmjcjet(np(Hplus)), &
                                       CS_Hpmjjetjet(np(Hplus)), &
                                       CS_HpmjW(np(Hplus)), CS_HpmjZ(np(Hplus)), &
                                       CS_vbf_Hpmj(np(Hplus)), CS_HpjHmj(np(Hplus))
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), intent(in) :: CS_Hpmjhi_C(np(Hneut), np(Hplus))
    double precision :: CS_HpmjHi(np(Hplus), np(Hneut))

    CS_HpmjHi = TRANSPOSE(CS_HpmjHi_C)
    call HiggsBounds_charged_input_hadr(collider, CS_Hpmjtb, CS_Hpmjcb, &
                                        CS_Hpmjbjet, CS_Hpmjcjet, CS_Hpmjjetjet, CS_HpmjW, &
                                        CS_HpmjZ, CS_vbf_Hpmj, CS_HpjHmj, CS_Hpmjhi)
end subroutine

subroutine HiggsBounds_charged_input_effC_fermions_C( &
    hcjud_L, hcjud_R, hcjcs_L, hcjcs_R, hcjtb_L, hcjtb_R, &
    hcjus_L, hcjus_R, hcjub_L, hcjub_R, hcjcd_L, hcjcd_R, &
    hcjcb_L, hcjcb_R, hcjtd_L, hcjtd_R, hcjts_L, hcjts_R) &
    bind(C, name='HiggsBounds_charged_input_effC_fermions')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hplus

    implicit none
    real(kind=c_double), intent(in) :: hcjud_L(np(Hplus)), hcjud_R(np(Hplus)), &
                                       hcjcs_L(np(Hplus)), hcjcs_R(np(Hplus)), &
                                       hcjtb_L(np(Hplus)), hcjtb_R(np(Hplus)), &
                                       hcjus_L(np(Hplus)), hcjus_R(np(Hplus)), &
                                       hcjub_L(np(Hplus)), hcjub_R(np(Hplus)), &
                                       hcjcd_L(np(Hplus)), hcjcd_R(np(Hplus)), &
                                       hcjcb_L(np(Hplus)), hcjcb_R(np(Hplus)), &
                                       hcjtd_L(np(Hplus)), hcjtd_R(np(Hplus)), &
                                       hcjts_L(np(Hplus)), hcjts_R(np(Hplus))
    call HiggsBounds_charged_input_effC_fermions( &
        hcjud_L, hcjud_R, hcjcs_L, hcjcs_R, hcjtb_L, hcjtb_R, &
        hcjus_L, hcjus_R, hcjub_L, hcjub_R, hcjcd_L, hcjcd_R, &
        hcjcb_L, hcjcb_R, hcjtd_L, hcjtd_R, hcjts_L, hcjts_R)
end subroutine

subroutine HiggsBounds_get_neutral_hadr_CS_C(i, collider, &
                                             singleH, ggH, bbH, VBF, WH, ZH, ttH, &
                                             tH_tchan, tH_schan, qqZH, ggZH) &
    bind(C, name='HiggsBounds_get_neutral_hadr_CS')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: i, collider
    real(kind=c_double), intent(out) :: singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH, ggZH

    call HiggsBounds_get_neutral_hadr_CS(i, collider, &
                                         singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH, ggZH)
end subroutine

subroutine HiggsBounds_get_neutral_BR_C(i, BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, &
                                        BR_hjmumu, BR_hjtautau, &
                                        BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg) &
    bind(C, name='HiggsBounds_get_neutral_BR')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: i
    real(kind=c_double), intent(out) :: BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, &
                                        BR_hjmumu, BR_hjtautau, &
                                        BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg

    call HiggsBounds_get_neutral_BR(i, BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt, &
                                    BR_hjmumu, BR_hjtautau, &
                                    BR_hjWW, BR_hjZZ, BR_hjZga, BR_hjgaga, BR_hjgg)
end subroutine

subroutine HiggsBounds_set_mass_uncertainties_C(dMhneut, dMhch) &
    bind(C, name='HiggsBounds_set_mass_uncertainties')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    real(kind=c_double), intent(in) :: dMhneut(np(Hneut)), dMhch(np(Hplus))

    call HiggsBounds_set_mass_uncertainties(dMhneut, dMhch)
end subroutine

subroutine run_HiggsBounds_C(HBresult, chan, obsratio, ncombined) &
    bind(C, name='run_HiggsBounds')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: HBresult, chan, ncombined
    real(kind=c_double), intent(out) :: obsratio

    call run_HiggsBounds(HBresult, chan, obsratio, ncombined)
end subroutine

subroutine run_HiggsBounds_single_C(h, HBresult, chan, obsratio, ncombined) &
    bind(C, name='run_HiggsBounds_single')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: h
    integer(kind=c_int), intent(out) :: HBresult, chan, ncombined
    real(kind=c_double), intent(out) :: obsratio

    call run_HiggsBounds_single(h, HBresult, chan, obsratio, ncombined)
end subroutine

subroutine run_higgsBounds_full_C(HBresult, chan, obsratio, ncombined) &
    bind(C, name='run_HiggsBounds_full')
    use, intrinsic :: iso_c_binding
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    integer(kind=c_int), intent(out)::  HBresult(0:np(Hneut) + np(Hplus))
    integer(kind=c_int), intent(out)::  chan(0:np(Hneut) + np(Hplus))
    integer(kind=c_int), intent(out)::  ncombined(0:np(Hneut) + np(Hplus))
    real(kind=c_double), intent(out) :: obsratio(0:np(Hneut) + np(Hplus))

    call run_higgsBounds_full(HBresult, chan, obsratio, ncombined)
end subroutine

subroutine run_HiggsBounds_classic_C(HBresult, chan, obsratio, ncombined) &
    bind(C, name='run_HiggsBounds_classic')
    use, intrinsic :: iso_c_binding

    implicit none
    integer(kind=c_int), intent(out) :: HBresult, chan, ncombined
    real(kind=c_double), intent(out) :: obsratio

    call run_HiggsBounds_classic(HBresult, chan, obsratio, ncombined)
end subroutine

subroutine HiggsBounds_get_most_sensitive_channels_per_Higgs_C(nH, pos, &
                                                               HBresult, &
                                                               chan, &
                                                               obsratio, &
                                                               predratio, &
                                                               ncombined) &
    bind(C, name='HiggsBounds_get_most_sensitive_channels_per_Higgs')
    use, intrinsic:: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: nH, pos
    integer(kind=c_int), intent(out) :: HBresult, chan, ncombined
    real(kind=c_double), intent(out) :: obsratio, predratio

    call HiggsBounds_get_most_sensitive_channels_per_Higgs(nH, pos, HBresult, &
                                                           chan, obsratio, &
                                                           predratio, ncombined)
end subroutine

subroutine HiggsBounds_get_most_sensitive_channels_C(pos, HBresult, chan, &
                                                     obsratio, predratio, &
                                                     ncombined) &
    bind(C, name='HiggsBounds_get_most_sensitive_channels')
    use, intrinsic:: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: pos
    integer(kind=c_int), intent(out) :: HBresult, chan, ncombined
    real(kind=c_double), intent(out) :: obsratio, predratio

    call HiggsBounds_get_most_sensitive_channels(pos, HBresult, chan, &
                                                 obsratio, predratio, &
                                                 ncombined)
end subroutine

! subroutine run_HiggsBounds_classic( HBresult,chan,obsratio,ncombined)

subroutine HiggsBounds_get_likelihood_C(analysisID, Hindex, nc, cbin, M, &
                                        llh, obspred_flag) &
    bind(C, name='HiggsBounds_get_likelihood')
    use, intrinsic:: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: analysisID, obspred_flag
    integer(kind=c_int), intent(out) :: Hindex, nc, cbin
    real(kind=c_double), intent(out) :: llh, M

    select case (obspred_flag)
    case (0)
        call HiggsBounds_get_likelihood(analysisID, Hindex, nc, cbin, M, llh, &
                                        'pred')
    case (1)
        call HiggsBounds_get_likelihood(analysisID, Hindex, nc, cbin, M, llh, &
                                        'obs')
    case default
        print *, 'The function HiggsBounds_get_likelihood received an invalid obspred_flag'
        stop 1
    end select
end subroutine

subroutine HiggsBounds_get_likelihood_for_Higgs_C(analysisID, cbin_in, Hindex, &
                                                  nc, cbin, M, llh, obspred_flag) &
    bind(C, name='HiggsBounds_get_likelihood_for_Higgs')
    use, intrinsic:: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: analysisID, Hindex, cbin_in, obspred_flag
    integer(kind=c_int), intent(out) ::  nc, cbin
    real(kind=c_double), intent(out) :: llh, M

    select case (obspred_flag)
    case (0)
        call HiggsBounds_get_likelihood_for_Higgs(analysisID, cbin_in, Hindex, &
                                                  nc, cbin, M, llh, 'pred')
    case (1)
        call HiggsBounds_get_likelihood_for_Higgs(analysisID, cbin_in, Hindex, &
                                                  nc, cbin, M, llh, 'obs')
    case default
        print *, 'The function HiggsBounds_get_likelihood_for_Higgs received an invalid obspred_flag'
        stop 1
    end select
end subroutine

subroutine HiggsBounds_get_likelihood_for_comb_C(analysisID, cbin_in, Hindex, &
                                                 nc, cbin, M, llh, obspred_flag) &
    bind(C, name='HiggsBounds_get_likelihood_for_comb')
    use, intrinsic:: iso_c_binding

    implicit none
    integer(kind=c_int), value, intent(in) :: analysisID, cbin_in, obspred_flag
    integer(kind=c_int), intent(out) ::  Hindex, nc, cbin
    real(kind=c_double), intent(out) :: llh, M

    select case (obspred_flag)
    case (0)
        call HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, 'pred')
    case (1)
        call HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, 'obs')
    case default
        print *, 'The function HiggsBounds_get_likelihood_for_comb received an invalid obspred_flag'
        stop 1
    end select
end subroutine

! subroutine HiggsBounds_SLHA_output

#ifdef enableCHISQ
subroutine initialize_HiggsBounds_chisqtables_C() &
    bind(C, name='initialize_HiggsBounds_chisqtables')
    use, intrinsic :: iso_c_binding

    call initialize_HiggsBounds_chisqtables
end subroutine

subroutine HiggsBounds_get_LEPCHisq_C(theory_uncertainty_1s, chisq_withouttheory, chisq_withtheory, channel) &
    bind(C, name='HiggsBounds_get_LEPChisq')
    use, intrinsic :: iso_c_binding

    implicit none
    real(kind=c_double), intent(in), value :: theory_uncertainty_1s
    integer(kind=c_int), intent(out) :: channel
    real(kind=c_double), intent(out) :: chisq_withouttheory, chisq_withtheory

    call HiggsBounds_get_LEPCHisq(theory_uncertainty_1s, chisq_withouttheory, chisq_withtheory, channel)
end subroutine

subroutine finish_HiggsBounds_chisqtables_C() &
    bind(C, name='finish_HiggsBounds_chisqtables')
    use, intrinsic :: iso_c_binding

    call finish_HiggsBounds_chisqtables
end subroutine

! subroutine HB_calc_stats(theory_uncertainty_1s,chisq_withouttheory,chisq_withtheory,chan2)
#endif

subroutine finish_HiggsBounds_C() &
    bind(C, name='finish_HiggsBounds')
    use, intrinsic :: iso_c_binding

    call finish_HiggsBounds
end subroutine

! subroutine HiggsBounds_neutral_input_effC_single(quantity,val)

! subroutine HiggsBounds_neutral_input_effC_double(quantity,val)

! subroutine HiggsBounds_neutral_input_LEP_single(quantity,val)

! subroutine HiggsBounds_neutral_input_LEP_double(quantity,val)

! subroutine HiggsBounds_neutral_input_hadr_single(collider,quantity,val)

! subroutine HiggsBounds_neutral_input_hadr_double(collider,quantity,val)

subroutine HiggsBounds_neutral_input_hadr_channelrates_single_C(collider, &
                                                                nHiggs, p, d, &
                                                                val) &
    bind(C, name='HiggsBounds_neutral_input_hadr_channelrates_single')
    use, intrinsic :: iso_c_binding

    implicit none
    real(kind=c_double), value, intent(in) :: val
    integer(kind=c_int), value, intent(in) :: collider, p, d, nHiggs

    call HiggsBounds_neutral_input_hadr_channelrates_single(collider, nHiggs, &
                                                            p, d, val)
end subroutine

subroutine HiggsBounds_neutral_input_hadr_channelrates_clean_C() &
    bind(C, name='HiggsBounds_neutral_input_hadr_channelrates_clean')
    use, intrinsic :: iso_c_binding

    call HiggsBounds_neutral_input_hadr_channelrates_clean
end subroutine
