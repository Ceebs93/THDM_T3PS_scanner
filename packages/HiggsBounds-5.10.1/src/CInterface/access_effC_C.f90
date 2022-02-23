function SMCS_effC_HZ_C(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p) &
    bind(C, name="SMCS_effC_HZ")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: collider_name, ZH_cpmix_nnlo_ggqqbb

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), value, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p
    real(kind=c_double) :: SMCS_effC_HZ_C

    call testBRSM

    SMCS_effC_HZ_C = ZH_cpmix_nnlo_ggqqbb(Mh, collider_name(collider), &
                                          ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function

function SMCS_effC_gg_HZ_C(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p) &
    bind(C, name="SMCS_effC_gg_HZ")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: collider_name, ZH_cpmix_nnlo_gg

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), value, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p
    real(kind=c_double) :: SMCS_effC_gg_HZ_C

    call testBRSM

    SMCS_effC_gg_HZ_C = ZH_cpmix_nnlo_gg(Mh, collider_name(collider), &
                                         ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function

function SMCS_effC_qq_HZ_C(Mh, collider, ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p) &
    bind(C, name="SMCS_effC_qq_HZ")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: collider_name, ZH_cpmix_nnlo_qqbb

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), value, intent(in) :: ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p
    real(kind=c_double) :: SMCS_effC_qq_HZ_C

    call testBRSM

    SMCS_effC_qq_HZ_C = ZH_cpmix_nnlo_qqbb(Mh, collider_name(collider), &
                                           ghZZ, ghtt_s, ghbb_s, ghtt_p, ghbb_p, .True.)
end function

function SMCS_effC_bb_HZ_C(Mh, collider, ghbb_s, ghbb_p) &
    bind(C, name="SMCS_effC_bb_HZ")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: collider_name, ZH_cpmix_nnlo_bb

    implicit none
    real(kind=c_double), value, intent(in) :: Mh
    integer(kind=c_int), intent(in) :: collider
    real(kind=c_double), value, intent(in) :: ghbb_s, ghbb_p
    real(kind=c_double) :: SMCS_effC_bb_HZ_C

    call testBRSM

    SMCS_effC_bb_HZ_C = ZH_cpmix_nnlo_bb(Mh, collider_name(collider), &
                                         ghbb_s, ghbb_p, .True.)
end function

function SMCS_effC_HW_C(Mh, collider, ghWW, ghtt_s, ghbb_s) &
    bind(C, name="SMCS_effC_HW")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: collider_name, WH_nnlo
    implicit none
    real(kind = c_double),value, intent(in) :: Mh
    integer(kind=c_int), value, intent(in) :: collider
    real(kind=c_double), value, intent(in) :: ghWW, ghtt_s,ghbb_s
    real(kind=c_double) :: SMCS_effC_HW_C

    call testBRSM

    SMCS_effC_HW_C = WH_nnlo(Mh, collider_name(collider), ghWW, ghtt_s, ghbb_s, .True., .False.)
end function SMCS_effC_HW_C

function HCCS_tHc_C(MHc, gHcjt, gHcjb, BR_tHpjb) &
    bind(C, name="HCCS_tHc")
    use, intrinsic :: iso_c_binding
    use theory_XS_SM_functions, only: tHc_cxn

    implicit none
    real(kind=c_double), value, intent(in) :: mHc, gHcjt, gHcjb, BR_tHpjb
    real(kind=c_double) :: HCCS_tHc_C
    HCCS_tHc_C = tHc_cxn(MHc, gHcjt, gHcjb, BR_tHpjb)
end function
