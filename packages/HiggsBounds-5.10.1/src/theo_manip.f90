!> @file
!! Contains the module #theo_manip

!> Performs internal manipulations of the theory input.
module theo_manip
    use usefulbits, only: ndat, np, Hneut, Hplus, anyH, theo, partR, hadroncolliderextras, pdesc

    implicit none

contains
    !> Completes the theory input
    subroutine HB5_complete_theo
        use usefulbits, only: whichanalyses, whichinput, BRdirectinput
        implicit none

        if (np(Hneut) > 0) then

            select case (whichinput)
            case ('effC')
                call HB5_csratios_from_effC
                call HB5_cp_from_effC
                if (.not. BRdirectinput) then
                    call HB5_br_from_effC
                end if
            case ('SLHA')
                call HB5_csratios_from_effC
                call HB5_cp_from_effC
            case ('hadr')
            case default
                stop 'error in subroutine complete_theo (2): unknown whichinput!'
            end select
        end if

        call complete_BRs

        call check_dataset

        if (np(Hneut) > 0) then
            select case (whichanalyses)
            case ('onlyH', 'LandH', 'onlyP', 'list ')
!  n.b. there's no LEP SM cross sections at the moment
                call fill_theo_SM
! HB-5.2: calculate hadronic channelrates (need SM reference values)
                call complete_channelrates
            case ('onlyL')
            case default
                stop 'error in subroutine complete_theo (2): unknown whichinput!'
            end select
        end if

    end subroutine HB5_complete_theo

    !> Does the same as complete_theo but just for a single datapoint.
    !! @param n index of the datapoint
    subroutine HB5_recalculate_theo_for_datapoint(n)

        use usefulbits, only: whichanalyses, whichinput, BRdirectinput
        implicit none
        integer, intent(in) :: n

        if (np(Hneut) > 0) then
            select case (whichinput)
            case ('effC')
                call HB5_csratios_from_effC_for_datapoint(n)
                if (.not. BRdirectinput) then
                    call HB5_br_from_effC_for_datapoint(n)
                end if
            case ('SLHA')
                call HB5_csratios_from_effC_for_datapoint(n)
            case ('hadr', 'part')
            case default
                stop 'error in subroutine recalculate_theo_for_datapoint (1)'
            end select
        end if

        if (np(Hplus) > 0 .and. whichinput .eq. 'effC') then
            call HB5_charged_csratios_from_effC_for_datapoint(n)
        end if

        call check_dataset ! Checks consistency in BRs and total width

        if (np(Hneut) > 0) then
            select case (whichanalyses)
            case ('onlyH', 'LandH', 'onlyP', 'list ')
                ! n.b. there's no LEP SM cross sections at the moment
                call fill_theo_SM_for_datapoint(n)
            case ('onlyL')
            case default
                stop 'error in subroutine recalculate_theo_for_datapoint (2)'
            end select
        end if

        ! write (*, *) '# --------- complete_theo debugging --------- #'
        ! write (*, *) 'XS(ggH)_norm at TeV: ', theo(1)%tev%XS_gg_hj_ratio
        ! write (*, *) 'XS(bbH)_norm at TeV: ', theo(1)%tev%XS_bb_hj_ratio
        ! write (*, *) 'XS(VBF)_norm at Tev: ', theo(1)%tev%XS_vbf_ratio
        ! write (*, *) 'XS(HZ)_norm at Tev: ', theo(1)%tev%XS_hjZ_ratio
        ! write (*, *) 'XS(HW)_norm at TeV: ', theo(1)%tev%XS_hjW_ratio
        ! write (*, *) 'XS(ttH)_norm at TeV: ', theo(1)%tev%XS_tthj_ratio
        ! write (*, *) 'XS(tH)_norm at TeV: ', theo(1)%tev%XS_thj_tchan_ratio
        !
        ! write (*, *) 'XS(ggH)_norm at 7 TeV: ', theo(1)%lhc7%XS_gg_hj_ratio
        ! write (*, *) 'XS(bbH)_norm at 7 TeV: ', theo(1)%lhc7%XS_bb_hj_ratio
        ! write (*, *) 'XS(VBF)_norm at 7 Tev: ', theo(1)%lhc7%XS_vbf_ratio
        ! write (*, *) 'XS(HZ)_norm at 7 TeV: ', theo(1)%lhc7%XS_hjZ_ratio
        ! write (*, *) 'XS(HW)_norm at 7 TeV: ', theo(1)%lhc7%XS_hjW_ratio
        ! write (*, *) 'XS(ttH)_norm at 7 TeV: ', theo(1)%lhc7%XS_tthj_ratio
        ! write (*, *) 'XS(tH)_norm at 7 TeV: ', theo(1)%lhc7%XS_thj_tchan_ratio
        !
        ! write (*, *) 'XS(ggH)_norm at 8 TeV: ', theo(1)%lhc8%XS_gg_hj_ratio
        ! write (*, *) 'XS(bbH)_norm at 8 TeV: ', theo(1)%lhc8%XS_bb_hj_ratio
        ! write (*, *) 'XS(VBF)_norm at 8 Tev: ', theo(1)%lhc8%XS_vbf_ratio
        ! write (*, *) 'XS(HZ)_norm at 8 TeV: ', theo(1)%lhc8%XS_hjZ_ratio
        ! write (*, *) 'XS(HW)_norm at 8 TeV: ', theo(1)%lhc8%XS_hjW_ratio
        ! write (*, *) 'XS(ttH)_norm at 8 TeV: ', theo(1)%lhc8%XS_tthj_ratio
        ! write (*, *) 'XS(tH)_norm at 8 TeV: ', theo(1)%lhc8%XS_thj_tchan_ratio
        !
        ! write (*, *) 'XS(ggH)_norm at 13 TeV: ', theo(1)%lhc13%XS_gg_hj_ratio
        ! write (*, *) 'XS(bbH)_norm at 13 TeV: ', theo(1)%lhc13%XS_bb_hj_ratio
        ! write (*, *) 'XS(VBF)_norm at 13 Tev: ', theo(1)%lhc13%XS_vbf_ratio
        ! write (*, *) 'XS(HZ)_norm at 13 TeV: ', theo(1)%lhc13%XS_hjZ_ratio
        ! write (*, *) 'XS(HW)_norm at 13 TeV: ', theo(1)%lhc13%XS_hjW_ratio
        ! write (*, *) 'XS(ttH)_norm at 13 TeV: ', theo(1)%lhc13%XS_tthj_ratio
        ! write (*, *) 'XS(tH)_norm at 13 TeV: ', theo(1)%lhc13%XS_thj_tchan_ratio
        ! write (*, *) '# ---------      end debugging      --------- #'

    end subroutine HB5_recalculate_theo_for_datapoint

    !> Calls hb5_csratios_from_effc_for_datapoint for each datapoint
    subroutine HB5_csratios_from_effC
        use usefulbits, only: ndat
        implicit none
        integer :: jj

        if (np(Hneut) < 1) stop 'error in csratios_from_g2 (np(Hneut))'

        do jj = 1, ndat
            call HB5_csratios_from_effC_for_datapoint(jj)
            call HB5_charged_csratios_from_effC_for_datapoint(jj)
        end do

    end subroutine HB5_csratios_from_effC

    !> Uses the effective couplings to calculate the hadronic cross section ratios
    !! @param jj dataset index
    subroutine HB5_csratios_from_effC_for_datapoint(jj)
        use usefulbits, only: effC
        use theory_colliderSfunctions
        use theory_XS_SM_functions
        use S95tables, only: inrange
        implicit none
        integer, intent(in) :: jj
        ! internal
        integer :: i, j
        double precision :: TEVSM_ZZ_contrib_to_VBF, TEVSM_WW_contrib_to_VBF
        double precision :: Mhi

! relative contributuion of WW- and ZZ-fusion to VBF (in LO) for
! p p-bar collisions at SqrtS=1.96 TeV (calcuated by T. Figy with VBFNLO):s
        TEVSM_ZZ_contrib_to_VBF = 0.23D0
        TEVSM_WW_contrib_to_VBF = 0.77D0

        do i = 1, np(Hneut)

            Mhi = theo(jj)%particle(Hneut)%M(i)

            !---------------------------------------!
            !                   LEP                 !
            !---------------------------------------!
            theo(jj)%lep%XS_hjZ_ratio(i) = effC(jj)%hjZZ(i)**2
            theo(jj)%lep%XS_bbhj_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2
            !n.b.: LEP tables with bbhj at the moment can not be applied to mixed CP Higgs
            theo(jj)%lep%XS_tautauhj_ratio(i) = effC(jj)%hjtautau_s(i)**2 + effC(jj)%hjtautau_p(i)**2
            !n.b.: LEP tables with tautauhj at the moment can not be applied to mixed CP Higgs

            !---------------------------------------!
            !                TEVATRON               !
            !---------------------------------------!
            theo(jj)%tev%XS_gg_hj_ratio(i) = effC(jj)%hjgg(i)**2
            theo(jj)%tev%XS_bb_hj_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2
            theo(jj)%tev%XS_hjb_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2  ! still needed?
            ! calculate inclusive single Higgs production:
            ! n.b.: neglect cc,ss->hj for inclusive single Higgs production (in effC approximation)
            if (inrange(Mhi, 'TEV  ')) then
                theo(jj)%tev%XS_hj_ratio(i) = (theo(jj)%tev%XS_gg_hj_ratio(i)*XS_tev_gg_H_SM(Mhi) + &
                                               theo(jj)%tev%XS_bb_hj_ratio(i)*XS_tev_bb_H_SM(Mhi))/ &
                                              (XS_tev_gg_H_SM(Mhi) + XS_tev_bb_H_SM(Mhi))
            else
                theo(jj)%tev%XS_hj_ratio(i) = 0.0D0
            end if
            theo(jj)%tev%XS_vbf_ratio(i) = effC(jj)%hjWW(i)**2*TEVSM_WW_contrib_to_VBF &
                                           + effC(jj)%hjZZ(i)**2*TEVSM_ZZ_contrib_to_VBF

! TS 13/09/2019: Approximate this by fitted polynomial to MadGraph-5 predictions (for now, strictly valid only for LHC13 and Mh=125 GeV!)
            theo(jj)%tev%XS_tthj_ratio(i) = ttH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i))
            theo(jj)%tev%XS_thj_tchan_ratio(i) = tH_tchan_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
!            theo(jj)%tev%XS_tthj_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
!            theo(jj)%tev%XS_thj_tchan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
            ! THIS IS NOT ACCURATE, as the CP-odd component contributes differently (see ttH, tH tchan).
            theo(jj)%tev%XS_thj_schan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
            ! n.b.: Tev tables for tthj at the moment can only use CP even Higgs
            if (inrange(Mhi, 'TEV  ')) then

                theo(jj)%tev%XS_hjW_ratio(i) = WH_nnlo(Mhi, 'TEV  ', effC(jj)%hjWW(i), &
                                                       effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), .True., .True.)/ &
                                               WH_nnlo(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

                theo(jj)%tev%XS_hjZ_ratio(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'TEV  ', effC(jj)%hjZZ(i), &
                                                                    effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                    effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                               ZH_cpmix_nnlo_ggqqbb(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%tev%XS_gg_hjZ_ratio(i) = ZH_cpmix_nnlo_gg(Mhi, 'TEV  ', effC(jj)%hjZZ(i), &
                                                                   effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                   effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                  ZH_cpmix_nnlo_gg(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%tev%XS_qq_hjZ_ratio(i) = ZH_cpmix_nnlo_qqbb(Mhi, 'TEV  ', effC(jj)%hjZZ(i), &
                                                                     effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                     effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                  ZH_cpmix_nnlo_qqbb(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

            else
                theo(jj)%tev%XS_hjW_ratio(i) = 0.0D0
                theo(jj)%tev%XS_hjZ_ratio(i) = 0.0D0
                theo(jj)%tev%XS_gg_hjZ_ratio(i) = 0.0D0
                theo(jj)%tev%XS_qq_hjZ_ratio(i) = 0.0D0
            end if
            !---------------------------------------!
            !                  LHC 7                !
            !---------------------------------------!
            theo(jj)%lhc7%XS_gg_hj_ratio(i) = effC(jj)%hjgg(i)**2
            theo(jj)%lhc7%XS_bb_hj_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2
            theo(jj)%lhc7%XS_hjb_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2  ! still needed?
            ! calculate inclusive single Higgs production:
            ! n.b.: neglect cc,ss->hj for inclusive single Higgs production (in effC approximation)
            if (inrange(Mhi, 'LHC7 ')) then
                theo(jj)%lhc7%XS_hj_ratio(i) = (theo(jj)%lhc7%XS_gg_hj_ratio(i)*XS_lhc7_gg_H_SM(Mhi) + &
                                                theo(jj)%lhc7%XS_bb_hj_ratio(i)*XS_lhc7_bb_H_SM(Mhi))/ &
                                               (XS_lhc7_gg_H_SM(Mhi) + XS_lhc7_bb_H_SM(Mhi))
            else
                theo(jj)%lhc7%XS_hj_ratio(i) = 0.0D0
            end if
            if (inrange(Mhi, 'LHC7 ')) then
                theo(jj)%lhc7%XS_vbf_ratio(i) = effC(jj)%hjWW(i)**2*lhc7_rHVBF_WW(Mhi) + &
                                                effC(jj)%hjZZ(i)**2*lhc7_rHVBF_ZZ(Mhi)
            else
                theo(jj)%lhc7%XS_vbf_ratio(i) = 0.0D0
            end if
! TS 13/09/2019: Approximate this by fitted polynomial to MadGraph-5 predictions (for now, strictly valid only for LHC13 and Mh=125 GeV!)
            theo(jj)%lhc7%XS_tthj_ratio(i) = ttH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i))
            theo(jj)%lhc7%XS_thj_tchan_ratio(i) = tH_tchan_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
            theo(jj)%lhc7%XS_tWhj_ratio(i) = tWH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
!            theo(jj)%lhc7%XS_tthj_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
!            theo(jj)%lhc7%XS_thj_tchan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
! THIS IS NOT ACCURATE, as the CP-odd component contributes differently (see ttH, tH tchan).
            theo(jj)%lhc7%XS_thj_schan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
            ! n.b.: lhc7 tables for tthj at the moment can only use CP even Higgs
            if (inrange(Mhi, 'LHC7 ')) then

                theo(jj)%lhc7%XS_hjW_ratio(i) = WH_nnlo(Mhi, 'LHC7 ', effC(jj)%hjWW(i), &
                                                        effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), .True., .True.)/ &
                                                WH_nnlo(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

                theo(jj)%lhc7%XS_hjZ_ratio(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC7 ', effC(jj)%hjZZ(i), &
                                                                     effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                     effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc7%XS_gg_hjZ_ratio(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC7 ', effC(jj)%hjZZ(i), &
                                                                    effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                    effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                   ZH_cpmix_nnlo_gg(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc7%XS_qq_hjZ_ratio(i) = ZH_cpmix_nnlo_qqbb(Mhi, 'LHC7 ', effC(jj)%hjZZ(i), &
                                                                      effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                      effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                   ZH_cpmix_nnlo_qqbb(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

            else
                theo(jj)%lhc7%XS_hjW_ratio(i) = 0.0D0
                theo(jj)%lhc7%XS_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc7%XS_gg_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc7%XS_qq_hjZ_ratio(i) = 0.0D0
            end if

            !---------------------------------------!
            !                  LHC 8                !
            !---------------------------------------!
            theo(jj)%lhc8%XS_gg_hj_ratio(i) = effC(jj)%hjgg(i)**2
            theo(jj)%lhc8%XS_bb_hj_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2
            theo(jj)%lhc8%XS_hjb_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2  ! still needed?
            ! calculate inclusive single Higgs production:
            ! n.b.: neglect cc,ss->hj for inclusive single Higgs production (in effC approximation)
            if (inrange(Mhi, 'LHC8 ')) then
                theo(jj)%lhc8%XS_hj_ratio(i) = (theo(jj)%lhc8%XS_gg_hj_ratio(i)*XS_lhc8_gg_H_SM(Mhi) + &
                                                theo(jj)%lhc8%XS_bb_hj_ratio(i)*XS_lhc8_bb_H_SM(Mhi))/ &
                                               (XS_lhc8_gg_H_SM(Mhi) + XS_lhc8_bb_H_SM(Mhi))
            else
                theo(jj)%lhc8%XS_hj_ratio(i) = 0.0D0
            end if
            if (inrange(Mhi, 'LHC8 ')) then
                theo(jj)%lhc8%XS_vbf_ratio(i) = effC(jj)%hjWW(i)**2*lhc8_rHVBF_WW(Mhi) + &
                                                effC(jj)%hjZZ(i)**2*lhc8_rHVBF_ZZ(Mhi)
            else
                theo(jj)%lhc8%XS_vbf_ratio(i) = 0.0D0
            end if
! TS 13/09/2019: Approximate this by fitted polynomial to MadGraph-5 predictions (for now, strictly valid only for LHC13 and Mh=125 GeV!)
            theo(jj)%lhc8%XS_tthj_ratio(i) = ttH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i))
            theo(jj)%lhc8%XS_thj_tchan_ratio(i) = tH_tchan_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
            theo(jj)%lhc8%XS_tWhj_ratio(i) = tWH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
!            theo(jj)%lhc8%XS_tthj_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
!            theo(jj)%lhc8%XS_thj_tchan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
! THIS IS NOT ACCURATE, as the CP-odd component contributes differently (see ttH, tH tchan):
            theo(jj)%lhc8%XS_thj_schan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
            ! n.b.: lhc8 tables for tthj at the moment can only use CP even Higgs
            if (inrange(Mhi, 'LHC8 ')) then

                theo(jj)%lhc8%XS_hjW_ratio(i) = WH_nnlo(Mhi, 'LHC8 ', effC(jj)%hjWW(i), &
                                                        effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), .True., .True.)/ &
                                                WH_nnlo(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

                theo(jj)%lhc8%XS_hjZ_ratio(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC8 ', effC(jj)%hjZZ(i), &
                                                                     effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                     effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc8%XS_gg_hjZ_ratio(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC8 ', effC(jj)%hjZZ(i), &
                                                                    effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                    effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                   ZH_cpmix_nnlo_gg(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc8%XS_qq_hjZ_ratio(i) = ZH_cpmix_nnlo_qqbb(Mhi, 'LHC8 ', effC(jj)%hjZZ(i), &
                                                                      effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                      effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                   ZH_cpmix_nnlo_qqbb(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

            else
                theo(jj)%lhc8%XS_hjW_ratio(i) = 0.0D0
                theo(jj)%lhc8%XS_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc8%XS_gg_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc8%XS_qq_hjZ_ratio(i) = 0.0D0

            end if
            !---------------------------------------!
            !                  LHC 13                !
            !---------------------------------------!
            theo(jj)%lhc13%XS_gg_hj_ratio(i) = effC(jj)%hjgg(i)**2
            theo(jj)%lhc13%XS_bb_hj_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2
            theo(jj)%lhc13%XS_hjb_ratio(i) = effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2  ! still needed?
            ! calculate inclusive single Higgs production:
            ! n.b.: neglect cc,ss->hj for inclusive single Higgs production (in effC approximation)
            if (inrange(Mhi, 'LHC13')) then
                theo(jj)%lhc13%XS_hj_ratio(i) = (theo(jj)%lhc13%XS_gg_hj_ratio(i)*XS_lhc13_gg_H_SM(Mhi) + &
                                                 theo(jj)%lhc13%XS_bb_hj_ratio(i)*XS_lhc13_bb_H_SM(Mhi))/ &
                                                (XS_lhc13_gg_H_SM(Mhi) + XS_lhc13_bb_H_SM(Mhi))
                if ((Mhi .ge. 250D0) .and. (Mhi .le. 1050D0)) then
                    theo(jj)%lhc13%XS_hj_ratio(i) = theo(jj)%lhc13%XS_hj_ratio(i) + &
                                                    singleH_from_quarks_at_13TeV(jj, i, Mhi)/ &
                                                    (XS_lhc13_gg_H_SM(Mhi) + XS_lhc13_bb_H_SM(Mhi))
                end if
            else
                theo(jj)%lhc13%XS_hj_ratio(i) = 0.0D0
            end if
            if (inrange(Mhi, 'LHC13')) then
                theo(jj)%lhc13%XS_vbf_ratio(i) = effC(jj)%hjWW(i)**2*lhc13_rHVBF_WW(Mhi) + &
                                                 effC(jj)%hjZZ(i)**2*lhc13_rHVBF_ZZ(Mhi)
            else
                theo(jj)%lhc13%XS_vbf_ratio(i) = 0.0D0
            end if
! TS 13/09/2019: Approximate this by fitted polynomial to MadGraph-5 predictions (for now, strictly valid only for LHC13 and Mh=125 GeV!)
            theo(jj)%lhc13%XS_tthj_ratio(i) = ttH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i))
            theo(jj)%lhc13%XS_thj_tchan_ratio(i) = tH_tchan_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
            theo(jj)%lhc13%XS_tWhj_ratio(i) = tWH_modifier(effC(jj)%hjtt_s(i), effC(jj)%hjtt_p(i), effC(jj)%hjWW(i))
!            theo(jj)%lhc13%XS_tthj_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
!            theo(jj)%lhc13%XS_thj_tchan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
! THIS IS NOT ACCURATE, as the CP-odd component contributes differently (see ttH, tH tchan):
            theo(jj)%lhc13%XS_thj_schan_ratio(i) = effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2
            ! n.b.: lhc13 tables for tthj at the moment can only use CP even Higgs
            ! write (*, *) 'inrange(Mhi,LHC13) = ', inrange(Mhi, 'LHC13')
            if (inrange(Mhi, 'LHC13')) then

                theo(jj)%lhc13%XS_hjW_ratio(i) = WH_nnlo(Mhi, 'LHC13', effC(jj)%hjWW(i), &
                                                         effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), .True., .True.)/ &
                                                 WH_nnlo(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

                ! write (*, *) "theo(jj)%lhc13%XS_hjW_ratio(i) = ", theo(jj)%lhc13%XS_hjW_ratio(i)
                ! write (*, *) WH_nnlo(Mhi, 'LHC13', effC(jj)%hjWW(i), &
                !                      effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), .True., .True.)
                ! write (*, *) WH_nnlo(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                ! write (*, *) WH_nnlo_SM(Mhi, 'LHC13', .True., .True.)

                theo(jj)%lhc13%XS_hjZ_ratio(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC13', effC(jj)%hjZZ(i), &
                                                                      effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                      effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                 ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc13%XS_gg_hjZ_ratio(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC13', effC(jj)%hjZZ(i), &
                                                                     effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                     effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                    ZH_cpmix_nnlo_gg(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

                theo(jj)%lhc13%XS_qq_hjZ_ratio(i) = ZH_cpmix_nnlo_qqbb(Mhi, 'LHC13', effC(jj)%hjZZ(i), &
                                                                       effC(jj)%hjtt_s(i), effC(jj)%hjbb_s(i), &
                                                                       effC(jj)%hjtt_p(i), effC(jj)%hjbb_p(i), .True.)/ &
                                                    ZH_cpmix_nnlo_qqbb(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)

            else
                theo(jj)%lhc13%XS_hjW_ratio(i) = 0.0D0
                theo(jj)%lhc13%XS_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc13%XS_gg_hjZ_ratio(i) = 0.0D0
                theo(jj)%lhc13%XS_qq_hjZ_ratio(i) = 0.0D0
            end if

        end do

        theo(jj)%lep%XS_hjhi_ratio = effC(jj)%hjhiZ**2! note only half of XS_hjhi_ratio is filled here

        do j = 2, np(Hneut)
            do i = 1, j - 1
                theo(jj)%lep%XS_hjhi_ratio(i, j) = theo(jj)%lep%XS_hjhi_ratio(j, i)
            end do
        end do

    contains
        function singleH_from_quarks_at_13TeV(jj, i, Mhi)
            use theory_XS_SM_functions, only: neutral_scalar_from_quarks_lhc13
            integer, intent(in) :: jj, i
            double precision, intent(in) :: Mhi
            double precision :: singleH_from_quarks_at_13TeV, coupling
            integer :: iii

            singleH_from_quarks_at_13TeV = 0.0D0
            do iii = 1, 9 ! Loop over relevant couplings
                select case (iii) ! Select relevant coupling strength
                    ! n.b. make use effC directly.
                case (1)
                    coupling = sqrt(effC(jj)%hjuu_s(i)**2.0D0 + effC(jj)%hjuu_p(i)**2.0D0)
                case (2)
                    coupling = sqrt(effC(jj)%hjdd_s(i)**2.0D0 + effC(jj)%hjdd_p(i)**2.0D0)
                case (3)
                    coupling = sqrt(effC(jj)%hjcc_s(i)**2.0D0 + effC(jj)%hjcc_p(i)**2.0D0)
                case (4)
                    coupling = sqrt(effC(jj)%hjss_s(i)**2.0D0 + effC(jj)%hjss_p(i)**2.0D0)
                case (5)
                    !coupling = sqrt(effC(jj)%hjbb_s(i)**2.0D0 + effC(jj)%hjbb_p(i)**2.0D0)
                    coupling = 0.0D0 ! skip the bb->phi process, as it is already included.
                case (6)
                    coupling = sqrt(effC(jj)%hjuc_s(i)**2.0D0 + effC(jj)%hjuc_p(i)**2.0D0)
                case (7)
                    coupling = sqrt(effC(jj)%hjds_s(i)**2.0D0 + effC(jj)%hjds_p(i)**2.0D0)
                case (8)
                    coupling = sqrt(effC(jj)%hjdb_s(i)**2.0D0 + effC(jj)%hjdb_p(i)**2.0D0)
                case (9)
                    coupling = sqrt(effC(jj)%hjsb_s(i)**2.0D0 + effC(jj)%hjsb_p(i)**2.0D0)
                end select

                singleH_from_quarks_at_13TeV = singleH_from_quarks_at_13TeV + &
                                               neutral_scalar_from_quarks_lhc13(Mhi, coupling, iii, .True.)
            end do
        end function singleH_from_quarks_at_13TeV

    end subroutine HB5_csratios_from_effC_for_datapoint

    subroutine HB5_charged_csratios_from_effC_for_datapoint(jj)
        use usefulbits, only: np, Hplus
        implicit none
        integer, intent(in) :: jj

        integer :: i

        do i = 1, np(Hplus)
            theo(jj)%lhc13%XS_qq_Hpmj(i) = qq_Hplus_from_quarks_at_13TeV(theo(jj)%particle(Hplus)%M(i))
        end do

    contains
        double precision function qq_Hplus_from_quarks_at_13TeV(Mhi)
            use theory_XS_SM_functions, only: positive_charged_scalar_from_quarks_lhc13, negative_charged_scalar_from_quarks_lhc13
            use usefulbits, only: effC_Hc
            double precision, intent(in) :: Mhi
            double precision :: cL, cR
            integer :: iii

            qq_Hplus_from_quarks_at_13TeV = 0.0D0
            if (allocated(effC_Hc)) then
            do iii = 1, 6 ! Loop over relevant couplings
                select case (iii) ! Select relevant coupling strength
                    ! n.b. make use effC directly.
                case (1)
                    cL = effC_Hc(jj)%hcjud_L(i)
                    cR = effC_Hc(jj)%hcjud_R(i)
                case (2)
                    cL = effC_Hc(jj)%hcjcs_L(i)
                    cR = effC_Hc(jj)%hcjcs_R(i)
                case (3)
                    cL = effC_Hc(jj)%hcjus_L(i)
                    cR = effC_Hc(jj)%hcjus_R(i)
                case (4)
                    cL = effC_Hc(jj)%hcjcd_L(i)
                    cR = effC_Hc(jj)%hcjcd_R(i)
                case (5)
                    cL = effC_Hc(jj)%hcjub_L(i)
                    cR = effC_Hc(jj)%hcjub_R(i)
                case (6)
                    cL = effC_Hc(jj)%hcjcb_L(i)
                    cR = effC_Hc(jj)%hcjcb_R(i)
                end select

                qq_Hplus_from_quarks_at_13TeV = qq_Hplus_from_quarks_at_13TeV &
                                                + positive_charged_scalar_from_quarks_lhc13(Mhi, cL, cR, iii, .True.) &
                                                + negative_charged_scalar_from_quarks_lhc13(Mhi, cL, cR, iii, .True.)
            end do
            end if
        end function qq_Hplus_from_quarks_at_13TeV

    end subroutine HB5_charged_csratios_from_effC_for_datapoint

    !> Uses the effective couplings to calculate the cp property of neutral higgs
    subroutine HB5_cp_from_effC

        use usefulbits, only: effC, ndat, vsmall
        implicit none
        ! internal
        integer :: i, jj
        double precision :: max_hjff_s, max_hjff_p

        if (np(Hneut) < 1) stop 'error in cp_from_effC (np(Hneut))'

        do jj = 1, ndat

            do i = 1, np(Hneut)
                max_hjff_s = max(effC(jj)%hjss_s(i), effC(jj)%hjcc_s(i), effC(jj)%hjbb_s(i), &
                                 effC(jj)%hjtt_s(i), effC(jj)%hjmumu_s(i), effC(jj)%hjtautau_s(i))

                max_hjff_p = max(effC(jj)%hjss_p(i), effC(jj)%hjcc_p(i), effC(jj)%hjbb_p(i), &
                                 effC(jj)%hjtt_p(i), effC(jj)%hjmumu_p(i), effC(jj)%hjtautau_p(i))

                if (max_hjff_p .lt. vsmall) then !CP even
                    theo(jj)%CP_value(i) = 1
                elseif (max_hjff_s .lt. vsmall) then !CP odd
                    theo(jj)%CP_value(i) = -1
                else                              !mixed CP
                    theo(jj)%CP_value(i) = 0
                end if
            end do

        end do

    end subroutine HB5_cp_from_effC

    !> Calls the subroutine br_from_effc_for_datapoint for each datapoint
    subroutine HB5_br_from_effC
        use usefulbits, only: np, Hneut, ndat
        implicit none
        ! internal
        integer :: jj

        if (np(Hneut) < 1) stop 'error in br_from_effC (np(Hneut))'

        do jj = 1, ndat
            call HB5_br_from_effC_for_datapoint(jj)
        end do

    end subroutine HB5_br_from_effC

    !> Uses the effective couplings to calculate the branching ratios
    !! @param jj index of the dataset
    subroutine HB5_br_from_effC_for_datapoint(jj)
        use theory_BRfunctions
        use S95tables, only: inrange
        use usefulbits, only: effC, ms, mc, mt, mbmb, mmu, mtau, small, debug
        implicit none
        integer, intent(in) :: jj
        ! internal
        integer :: i, k, kk
        double precision :: Mhi, GammaRat, GammaTotSM, BRNP

        do i = 1, np(Hneut)
            Mhi = theo(jj)%particle(Hneut)%M(i)
            if (theo(jj)%particle(Hneut)%Mc(i) .ge. small) Mhi = theo(jj)%particle(Hneut)%Mc(i)

            theo(jj)%BR_hjss(i) = 0.0D0
            theo(jj)%BR_hjcc(i) = 0.0D0
            theo(jj)%BR_hjbb(i) = 0.0D0
            theo(jj)%BR_hjmumu(i) = 0.0D0
            theo(jj)%BR_hjtautau(i) = 0.0D0
            theo(jj)%BR_hjWW(i) = 0.0D0
            theo(jj)%BR_hjZZ(i) = 0.0D0
            theo(jj)%BR_hjZga(i) = 0.0D0
            theo(jj)%BR_hjgaga(i) = 0.0D0
            theo(jj)%BR_hjgg(i) = 0.0D0
            theo(jj)%BR_hjtt(i) = 0.0D0

            if (inrange(Mhi, 'SMBR')) then

                if (theo(jj)%particle(Hneut)%GammaTot(i) .lt. 0.0D0) then
                    if (debug) write (*, *) 'Derive total width from effective coupling input and BR(Higgs-to-NP)...'

                    GammaTotSM = (effC(jj)%hjss_s(i)**2 + effC(jj)%hjss_p(i)**2*invbsq(ms, Mhi)) &
                                 *BRSM_Hss(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + (effC(jj)%hjcc_s(i)**2 + effC(jj)%hjcc_p(i)**2*invbsq(mc, Mhi)) &
                                 *BRSM_Hcc(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + (effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2*invbsq(mbmb, Mhi)) &
                                 *BRSM_Hbb(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + (effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2*invbsq(mt, Mhi)) &
                                 *BRSM_Htoptop(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + (effC(jj)%hjmumu_s(i)**2 + effC(jj)%hjmumu_p(i)**2*invbsq(mmu, Mhi)) &
                                 *BRSM_Hmumu(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + (effC(jj)%hjtautau_s(i)**2 + effC(jj)%hjtautau_p(i)**2*invbsq(mtau, Mhi)) &
                                 *BRSM_Htautau(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + effC(jj)%hjWW(i)**2*BRSM_HWW(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + effC(jj)%hjZZ(i)**2*BRSM_HZZ(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + effC(jj)%hjZga(i)**2*BRSM_HZga(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + effC(jj)%hjgaga(i)**2*BRSM_Hgaga(Mhi)*BRSM_GammaTot(Mhi)
                    GammaTotSM = GammaTotSM + effC(jj)%hjgg(i)**2*BRSM_Hgg(Mhi)*BRSM_GammaTot(Mhi)

!                     call complete_BRs

                    BRNP = theo(jj)%BR_hjinvisible(i)
                    do k = lbound(theo(jj)%BR_hkhjhi, dim=2), ubound(theo(jj)%BR_hkhjhi, dim=2)
                        do kk = 1, k
                            BRNP = BRNP + theo(jj)%BR_hkhjhi(i, k, kk)
                        end do
                        BRNP = BRNP + theo(jj)%BR_hjhiZ(i, k)
                    end do
                    BRNP = BRNP + theo(jj)%BR_hjemu(i)
                    BRNP = BRNP + theo(jj)%BR_hjetau(i)
                    BRNP = BRNP + theo(jj)%BR_hjmutau(i)

                    BRNP = BRNP + theo(jj)%BR_hjuu(i) + theo(jj)%BR_hjdd(i) &
                           + theo(jj)%BR_hjee(i) + theo(jj)%BR_hjuc(i) &
                           + theo(jj)%BR_hjds(i) + theo(jj)%BR_hjut(i) &
                           + theo(jj)%BR_hjdb(i) + theo(jj)%BR_hjct(i) &
                           + theo(jj)%BR_hjsb(i) ! beyondHB extension

                    do k = lbound(theo(jj)%BR_hjHpiW, dim=2), ubound(theo(jj)%BR_hjHpiW, dim=2)
                        BRNP = BRNP + theo(jj)%BR_hjHpiW(i, k)
                    end do

                    theo(jj)%particle(Hneut)%GammaTot(i) = GammaTotSM/(1.0D0 - BRNP)
                    theo(jj)%particle(anyH)%GammaTot(i) = theo(jj)%particle(Hneut)%GammaTot(i)
                end if

                GammaRat = theo(jj)%particle(Hneut)%GammaTot(i)/BRSM_GammaTot(Mhi)

                if (theo(jj)%particle(Hneut)%GammaTot(i) .gt. 0.0D0) then
                    theo(jj)%BR_hjss(i) = (effC(jj)%hjss_s(i)**2 + effC(jj)%hjss_p(i)**2*invbsq(ms, Mhi)) &
                                          *BRSM_Hss(Mhi)/GammaRat
                    theo(jj)%BR_hjcc(i) = (effC(jj)%hjcc_s(i)**2 + effC(jj)%hjcc_p(i)**2*invbsq(mc, Mhi)) &
                                          *BRSM_Hcc(Mhi)/GammaRat
                    theo(jj)%BR_hjbb(i) = (effC(jj)%hjbb_s(i)**2 + effC(jj)%hjbb_p(i)**2*invbsq(mbmb, Mhi)) &
                                          *BRSM_Hbb(Mhi)/GammaRat
                    theo(jj)%BR_hjtt(i) = (effC(jj)%hjtt_s(i)**2 + effC(jj)%hjtt_p(i)**2*invbsq(mt, Mhi)) &
                                          *BRSM_Htoptop(Mhi)/GammaRat
                    theo(jj)%BR_hjmumu(i) = (effC(jj)%hjmumu_s(i)**2 + effC(jj)%hjmumu_p(i)**2*invbsq(mmu, Mhi)) &
                                            *BRSM_Hmumu(Mhi)/GammaRat
                    theo(jj)%BR_hjtautau(i) = (effC(jj)%hjtautau_s(i)**2 + effC(jj)%hjtautau_p(i)**2*invbsq(mtau, Mhi)) &
                                              *BRSM_Htautau(Mhi)/GammaRat

                    theo(jj)%BR_hjWW(i) = effC(jj)%hjWW(i)**2*BRSM_HWW(Mhi)/GammaRat
                    theo(jj)%BR_hjZZ(i) = effC(jj)%hjZZ(i)**2*BRSM_HZZ(Mhi)/GammaRat
                    theo(jj)%BR_hjZga(i) = effC(jj)%hjZga(i)**2*BRSM_HZga(Mhi)/GammaRat
                    theo(jj)%BR_hjgaga(i) = effC(jj)%hjgaga(i)**2*BRSM_Hgaga(Mhi)/GammaRat
                    theo(jj)%BR_hjgg(i) = effC(jj)%hjgg(i)**2*BRSM_Hgg(Mhi)/GammaRat

                    ! write (*, *) 'BR h->ss = ', theo(jj)%BR_hjss(i), 'SM =', BRSM_Hss(Mhi)
                    ! write (*, *) 'BR h->cc = ', theo(jj)%BR_hjcc(i), 'SM =', BRSM_Hcc(Mhi)
                    ! write (*, *) 'BR h->bb = ', theo(jj)%BR_hjbb(i), 'SM =', BRSM_Hbb(Mhi)
                    ! write (*, *) 'BR h->tt = ', theo(jj)%BR_hjtt(i), 'SM =', BRSM_Htoptop(Mhi)
                    ! write (*, *) 'BR h->mumu = ', theo(jj)%BR_hjmumu(i), 'SM =', BRSM_Hmumu(Mhi)
                    ! write (*, *) 'BR h->tautau = ', theo(jj)%BR_hjtautau(i), 'SM =', BRSM_Htautau(Mhi)
                    ! write (*, *) 'BR h->WW = ', theo(jj)%BR_hjWW(i), 'SM =', BRSM_HWW(Mhi)
                    ! write (*, *) 'BR h->ZZ = ', theo(jj)%BR_hjZZ(i), 'SM =', BRSM_HZZ(Mhi)
                    ! write (*, *) 'BR h->gaga = ', theo(jj)%BR_hjgaga(i), 'SM =', BRSM_Hgaga(Mhi)
                    ! write (*, *) 'BR h->gg = ', theo(jj)%BR_hjgg(i), 'SM =', BRSM_Hgg(Mhi)
                    ! write (*, *) 'BR h->Zga = ', theo(jj)%BR_hjZga(i), 'SM =', BRSM_HZga(Mhi)

                elseif (theo(jj)%particle(Hneut)%GammaTot(i) .lt. 0.0D0) then

                    write (*, *) 'at jj=', jj, 'i=', i
                    write (*, *) 'total decay width is less than zero:', theo(jj)%particle(Hneut)%GammaTot(i)
                end if
            end if
        end do

    end subroutine HB5_br_from_effC_for_datapoint

    !> Calls the subroutine complete_channelrates_for_datapoint for each datapoint.
    subroutine complete_channelrates
        use usefulbits, only: np, Hneut, ndat
        implicit none
        ! internal
        integer :: jj

        if (np(Hneut) < 1) stop 'error in complete_channelrates (np(Hneut))'

        do jj = 1, ndat
            call complete_channelrates_for_datapoint(jj)
        end do

    end subroutine complete_channelrates

    !> Obtains the channelrates either from the XS and BR input.
    !! Assumes the narrow width approximation), or, if provided directly, uses the user's input
    !! The important case of 0.0D0 will be taken over! (To enable treating the interference
    !! of several Higgs bosons in one slot and de-activate (i.e. set to zero) the other slot)
    !! @param jj index of the dataset
    subroutine complete_channelrates_for_datapoint(jj)
        use usefulbits, only: Nprod, Ndecay
        implicit none
        ! internal
        integer :: i, jj, p, d
        double precision :: sigma, BR

        do i = 1, np(Hneut)

            do p = 1, Nprod
                do d = 1, Ndecay
                    if (theo(jj)%tev%channelrates_tmp(i, p, d) .ge. 0.0D0) then
                        theo(jj)%tev%channelrates(i, p, d) = theo(jj)%tev%channelrates_tmp(i, p, d)
                    else
                        select case (p)
                        case (1)
                            sigma = theo(jj)%tev%XS_hj_ratio(i)
                        case (2)
                            sigma = theo(jj)%tev%XS_vbf_ratio(i)
                        case (3)
                            sigma = theo(jj)%tev%XS_hjW_ratio(i)
                        case (4)
                            sigma = theo(jj)%tev%XS_hjZ_ratio(i)
                        case (5)
                            sigma = theo(jj)%tev%XS_tthj_ratio(i)
                        case (6)
                            sigma = theo(jj)%tev%XS_gg_hj_ratio(i)
                        case (7)
                            sigma = theo(jj)%tev%XS_bb_hj_ratio(i)
                        case (8)
                            sigma = theo(jj)%tev%XS_thj_tchan_ratio(i)
                        case (9)
                            sigma = theo(jj)%tev%XS_thj_schan_ratio(i)
                        case (10)
                            sigma = theo(jj)%tev%XS_qq_hjZ_ratio(i)
                        case (11)
                            sigma = theo(jj)%tev%XS_gg_hjZ_ratio(i)
                        end select
                        select case (d)
                        case (1)
                            BR = theo(jj)%BR_hjgaga(i)
                        case (2)
                            BR = theo(jj)%BR_hjWW(i)
                        case (3)
                            BR = theo(jj)%BR_hjZZ(i)
                        case (4)
                            BR = theo(jj)%BR_hjtautau(i)
                        case (5)
                            BR = theo(jj)%BR_hjbb(i)
                        case (6)
                            BR = theo(jj)%BR_hjZga(i)
                        case (7)
                            BR = theo(jj)%BR_hjcc(i)
                        case (8)
                            BR = theo(jj)%BR_hjmumu(i)
                        case (9)
                            BR = theo(jj)%BR_hjgg(i)
                        case (10)
                            BR = theo(jj)%BR_hjss(i)
                        case (11)
                            BR = theo(jj)%BR_hjtt(i)
                        end select
                        theo(jj)%tev%channelrates(i, p, d) = sigma*BR
                    end if

                    if (theo(jj)%lhc7%channelrates_tmp(i, p, d) .ge. 0.0D0) then
                        theo(jj)%lhc7%channelrates(i, p, d) = theo(jj)%lhc7%channelrates_tmp(i, p, d)
                    else
                        select case (p)
                        case (1)
                            sigma = theo(jj)%lhc7%XS_hj_ratio(i)
                        case (2)
                            sigma = theo(jj)%lhc7%XS_vbf_ratio(i)
                        case (3)
                            sigma = theo(jj)%lhc7%XS_hjW_ratio(i)
                        case (4)
                            sigma = theo(jj)%lhc7%XS_hjZ_ratio(i)
                        case (5)
                            sigma = theo(jj)%lhc7%XS_tthj_ratio(i)
                        case (6)
                            sigma = theo(jj)%lhc7%XS_gg_hj_ratio(i)
                        case (7)
                            sigma = theo(jj)%lhc7%XS_bb_hj_ratio(i)
                        case (8)
                            sigma = theo(jj)%lhc7%XS_thj_tchan_ratio(i)
                        case (9)
                            sigma = theo(jj)%lhc7%XS_thj_schan_ratio(i)
                        case (10)
                            sigma = theo(jj)%lhc7%XS_qq_hjZ_ratio(i)
                        case (11)
                            sigma = theo(jj)%lhc7%XS_gg_hjZ_ratio(i)
                        end select
                        select case (d)
                        case (1)
                            BR = theo(jj)%BR_hjgaga(i)
                        case (2)
                            BR = theo(jj)%BR_hjWW(i)
                        case (3)
                            BR = theo(jj)%BR_hjZZ(i)
                        case (4)
                            BR = theo(jj)%BR_hjtautau(i)
                        case (5)
                            BR = theo(jj)%BR_hjbb(i)
                        case (6)
                            BR = theo(jj)%BR_hjZga(i)
                        case (7)
                            BR = theo(jj)%BR_hjcc(i)
                        case (8)
                            BR = theo(jj)%BR_hjmumu(i)
                        case (9)
                            BR = theo(jj)%BR_hjgg(i)
                        case (10)
                            BR = theo(jj)%BR_hjss(i)
                        case (11)
                            BR = theo(jj)%BR_hjtt(i)
                        end select
                        theo(jj)%lhc7%channelrates(i, p, d) = sigma*BR
                    end if

                    if (theo(jj)%lhc8%channelrates_tmp(i, p, d) .ge. 0.0D0) then
                        theo(jj)%lhc8%channelrates(i, p, d) = theo(jj)%lhc8%channelrates_tmp(i, p, d)
                    else
                        select case (p)
                        case (1)
                            sigma = theo(jj)%lhc8%XS_hj_ratio(i)
                        case (2)
                            sigma = theo(jj)%lhc8%XS_vbf_ratio(i)
                        case (3)
                            sigma = theo(jj)%lhc8%XS_hjW_ratio(i)
                        case (4)
                            sigma = theo(jj)%lhc8%XS_hjZ_ratio(i)
                        case (5)
                            sigma = theo(jj)%lhc8%XS_tthj_ratio(i)
                        case (6)
                            sigma = theo(jj)%lhc8%XS_gg_hj_ratio(i)
                        case (7)
                            sigma = theo(jj)%lhc8%XS_bb_hj_ratio(i)
                        case (8)
                            sigma = theo(jj)%lhc8%XS_thj_tchan_ratio(i)
                        case (9)
                            sigma = theo(jj)%lhc8%XS_thj_schan_ratio(i)
                        case (10)
                            sigma = theo(jj)%lhc8%XS_qq_hjZ_ratio(i)
                        case (11)
                            sigma = theo(jj)%lhc8%XS_gg_hjZ_ratio(i)
                        end select
                        select case (d)
                        case (1)
                            BR = theo(jj)%BR_hjgaga(i)
                        case (2)
                            BR = theo(jj)%BR_hjWW(i)
                        case (3)
                            BR = theo(jj)%BR_hjZZ(i)
                        case (4)
                            BR = theo(jj)%BR_hjtautau(i)
                        case (5)
                            BR = theo(jj)%BR_hjbb(i)
                        case (6)
                            BR = theo(jj)%BR_hjZga(i)
                        case (7)
                            BR = theo(jj)%BR_hjcc(i)
                        case (8)
                            BR = theo(jj)%BR_hjmumu(i)
                        case (9)
                            BR = theo(jj)%BR_hjgg(i)
                        case (10)
                            BR = theo(jj)%BR_hjss(i)
                        case (11)
                            BR = theo(jj)%BR_hjtt(i)
                        end select
                        theo(jj)%lhc8%channelrates(i, p, d) = sigma*BR
                    end if

                    if (theo(jj)%lhc13%channelrates_tmp(i, p, d) .ge. 0.0D0) then
                        theo(jj)%lhc13%channelrates(i, p, d) = theo(jj)%lhc13%channelrates_tmp(i, p, d)
                    else
                        select case (p)
                        case (1)
                            sigma = theo(jj)%lhc13%XS_hj_ratio(i)
                        case (2)
                            sigma = theo(jj)%lhc13%XS_vbf_ratio(i)
                        case (3)
                            sigma = theo(jj)%lhc13%XS_hjW_ratio(i)
                        case (4)
                            sigma = theo(jj)%lhc13%XS_hjZ_ratio(i)
                        case (5)
                            sigma = theo(jj)%lhc13%XS_tthj_ratio(i)
                        case (6)
                            sigma = theo(jj)%lhc13%XS_gg_hj_ratio(i)
                        case (7)
                            sigma = theo(jj)%lhc13%XS_bb_hj_ratio(i)
                        case (8)
                            sigma = theo(jj)%lhc13%XS_thj_tchan_ratio(i)
                        case (9)
                            sigma = theo(jj)%lhc13%XS_thj_schan_ratio(i)
                        case (10)
                            sigma = theo(jj)%lhc13%XS_qq_hjZ_ratio(i)
                        case (11)
                            sigma = theo(jj)%lhc13%XS_gg_hjZ_ratio(i)
                        end select
                        select case (d)
                        case (1)
                            BR = theo(jj)%BR_hjgaga(i)
                        case (2)
                            BR = theo(jj)%BR_hjWW(i)
                        case (3)
                            BR = theo(jj)%BR_hjZZ(i)
                        case (4)
                            BR = theo(jj)%BR_hjtautau(i)
                        case (5)
                            BR = theo(jj)%BR_hjbb(i)
                        case (6)
                            BR = theo(jj)%BR_hjZga(i)
                        case (7)
                            BR = theo(jj)%BR_hjcc(i)
                        case (8)
                            BR = theo(jj)%BR_hjmumu(i)
                        case (9)
                            BR = theo(jj)%BR_hjgg(i)
                        case (10)
                            BR = theo(jj)%BR_hjss(i)
                        case (11)
                            BR = theo(jj)%BR_hjtt(i)
                        end select
                        theo(jj)%lhc13%channelrates(i, p, d) = sigma*BR
                        ! write (*, *) "debug: getting 13 TeV from XS x BR:", i, p, d, sigma, BR
                    end if
                end do
            end do
        end do

    end subroutine complete_channelrates_for_datapoint

    !> Calls the subroutine clean_channelrates_for_datapoint for each datapoint.
    subroutine clean_channelrates
        use usefulbits, only: np, Hneut, ndat
        implicit none
        ! internal
        integer :: jj

        if (np(Hneut) < 1) stop 'error in complete_channelrates (np(Hneut))'
        do jj = 1, ndat
            call clean_channelrates_for_datapoint(jj)
        end do
    end subroutine clean_channelrates

    !> Fills all channelrates matrices with -1.0D0
    !! @param jj index of the datapoint
    subroutine clean_channelrates_for_datapoint(jj)
        implicit none
        ! internal
        integer :: jj

        theo(jj)%tev%channelrates = -1.0D0
        theo(jj)%tev%channelrates_tmp = -1.0D0
        theo(jj)%lhc7%channelrates = -1.0D0
        theo(jj)%lhc7%channelrates_tmp = -1.0D0
        theo(jj)%lhc8%channelrates = -1.0D0
        theo(jj)%lhc8%channelrates_tmp = -1.0D0
        theo(jj)%lhc13%channelrates = -1.0D0
        theo(jj)%lhc13%channelrates_tmp = -1.0D0
    end subroutine clean_channelrates_for_datapoint

    !> Completing calculation of Higgs branching ratios.
    !! This routine does some internal calculations for Higgs-to-Higgs decays and
    !! Higgs decays to invisible final states. For the latter, the invisible
    !! decay rate contribution from cascade decays through other Higgs-Higgs and Higgs-Z bosons
    !! is added to the value for the invisible decay rate provided by the user.
    subroutine complete_BRs
        use usefulbits, only: argsort
        implicit none
        integer :: jj, i, j, k, massOrder(np(Hneut))

        do jj = 1, ndat
            ! copying over the (k,i,i) elements of BR_hkhjhi into BR_hjhihi.
            do i = 1, np(Hneut)
                theo(jj)%BR_hjhihi(:, i) = theo(jj)%BR_hkhjhi(:, i, i)
            end do

            if (.not. theo(jj)%full_BR_inv) then
                ! get the full BR(hj->inv)
                ! traverse from lightest to heaviest to account for cascades
                call argsort(theo(jj)%particle(Hneut)%M(:), massOrder)
                do i = 2, np(Hneut)
                    do j = 1, i - 1
                        ! contribution from hi -> (hj -> inv) (Z -> vv)
                        theo(jj)%BR_hjinvisible(massOrder(i)) = theo(jj)%BR_hjinvisible(massOrder(i)) &
                                                                + theo(jj)%BR_hjhiZ(massOrder(i), massOrder(j)) &
                                                                *theo(jj)%BR_hjinvisible(massOrder(j))*0.2D0
                        do k = j, i - 1
                            ! contribution from hi -> hj hk -> inv
                            theo(jj)%BR_hjinvisible(massOrder(i)) = theo(jj)%BR_hjinvisible(massOrder(i)) &
                                                                    + theo(jj)%BR_hkhjhi(massOrder(i), massOrder(j), massOrder(k)) &
                                                                    *theo(jj)%BR_hjinvisible(massOrder(j)) &
                                                                    *theo(jj)%BR_hjinvisible(massOrder(k))
                        end do
                    end do
                end do
                theo(jj)%full_BR_inv = .true.
            end if
        end do
    end subroutine complete_BRs

    function invbsq(mf, mh)
        implicit none
        double precision, intent(in) :: mf, mh
        double precision :: invbsq
        if (mh > 2.0D0*mf) then
            invbsq = 1.0D0/(1.0D0 - 4.0D0*(mf/mh)**2.0D0)
        else
            invbsq = 0.0D0
        end if
    end function invbsq

    !> Checks for each datapoint if the Higgs masses and branching ratios make sense
    !! Sets theo(jj)%gooddataset accordingly
    subroutine check_dataset
        use usefulbits, only: theo, ndat, debug, np, vsmall
        implicit none
        ! internal
        integer :: jj, kk, ll, mm, x
        double precision :: testsumBR, testsumBR_t
        double precision, allocatable :: testBR(:)
        double precision, parameter :: fuzziness = 0.01D0
        double precision, allocatable :: sumhjhi(:), sumhjHpi(:)

        testsumBR = 0.0D0
        testsumBR_t = 0.0D0

        if (np(Hneut) > 0) then
            allocate (testBR(np(Hneut)))
            allocate (sumhjhi(np(Hneut)))
            allocate (sumhjHpi(np(Hneut)))

            ! testing to see if the dataset is ok
            do jj = 1, ndat

                do kk = 1, np(Hneut)
                    do ll = 1, np(Hneut)
                        do mm = 1, np(Hneut)
                            !   write(*,'(a,I2,a,I2,a,I2,a,1F10.8)') "BR_hkhjhi(",kk,",",ll,",",mm,")=",theo(jj)%BR_hkhjhi(kk,ll,mm)

                            if (abs(theo(jj)%BR_hkhjhi(kk, ll, mm) - theo(jj)%BR_hkhjhi(kk, mm, ll)) .gt. vsmall) then
                                if (theo(jj)%BR_hkhjhi(kk, ll, mm) .lt. vsmall) then
                                    theo(jj)%BR_hkhjhi(kk, ll, mm) = theo(jj)%BR_hkhjhi(kk, mm, ll)
!           write(*,'(a,I2,a,I2,a,I2,a)') "WARNING: BR_hkhjhi is not symmetric. Correcting BR_hkhjhi(",&
! &          kk,",",ll,",",mm,") element..."
                                else if (theo(jj)%BR_hkhjhi(kk, mm, ll) .lt. vsmall) then
                                    theo(jj)%BR_hkhjhi(kk, mm, ll) = theo(jj)%BR_hkhjhi(kk, ll, mm)
!           write(*,'(a,I2,a,I2,a,I2,a)') "WARNING: BR_hkhjhi is not symmetric. Correcting BR_hkhjhi(",&
! &          kk,",",mm,",",ll,") element..."
                                else
                                    write (*, *) "WARNING: BR_hkhjhi is not symmetric."
                                end if
                            end if
                        end do
                    end do
                end do

                sumhjhi = 0.0D0
                do kk = lbound(theo(jj)%BR_hkhjhi, dim=1), ubound(theo(jj)%BR_hkhjhi, dim=1)
                    do ll = lbound(theo(jj)%BR_hkhjhi, dim=2), ubound(theo(jj)%BR_hkhjhi, dim=2)
                        do mm = lbound(theo(jj)%BR_hkhjhi, dim=3), ll
                            sumhjhi(kk) = sumhjhi(kk) + theo(jj)%BR_hkhjhi(kk, ll, mm)
!        write(*,*) "kk,ll,mm, sumhjhi = ", kk, ll, mm, sumhjhi
                        end do
                    end do
                end do

                sumhjHpi = 0.0D0
                if (np(Hplus) .gt. 0) then
                    sumhjHpi = sum(theo(jj)%BR_hjHpiW, dim=2)
                end if

                testBR = theo(jj)%BR_hjss &
                         + theo(jj)%BR_hjcc &
                         + theo(jj)%BR_hjbb &
                         + theo(jj)%BR_hjtt &
                         + theo(jj)%BR_hjmumu &
                         + theo(jj)%BR_hjtautau &
                         + theo(jj)%BR_hjemu &
                         + theo(jj)%BR_hjetau &
                         + theo(jj)%BR_hjmutau &
                         + theo(jj)%BR_hjWW &
                         + theo(jj)%BR_hjZZ &
                         + theo(jj)%BR_hjZga &
                         + theo(jj)%BR_hjgg &
                         + theo(jj)%BR_hjgaga &
                         + sum(theo(jj)%BR_hjhiZ, dim=2) &
                         + sumhjhi + sumhjHpi
!         write(*,*) 'sumhjhi = ',sumhjhi
!     write(*,*) 'sum(theo(jj)%BR_hjhiZ,dim=2) = ', sum(theo(jj)%BR_hjhiZ,dim=2)
!     write(*,*) 'debug: testBR = ', testBR
                testsumBR = maxval(testBR)
                if (testsumBR .gt. 1.0D0 + fuzziness) then
                    write (*, *) 'warning: sum of BR for '//trim(adjustl(pdesc(Hneut)%long))// &
                        ' ', maxloc(testBR), ' at line number=', jj, 'is', testsumBR
                    write (*, *) 'BR(h', maxloc(testBR), '->WW)=', theo(jj)%BR_hjWW(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->ZZ)=', theo(jj)%BR_hjZZ(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->gg)=', theo(jj)%BR_hjgg(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->gaga)=', theo(jj)%BR_hjgaga(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->bb)=', theo(jj)%BR_hjbb(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->tautau)=', theo(jj)%BR_hjtautau(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->tt)=', theo(jj)%BR_hjtt(maxloc(testBR))
                    write (*, *) 'BR(h', maxloc(testBR), '->hiZ)=', theo(jj)%BR_hjhiZ(maxloc(testBR), :)
                    write (*, *) 'sum(BR(h', maxloc(testBR), '->hjhi))=', sumhjhi(maxloc(testBR))
                    write (*, *) 'sum(BR(h', maxloc(testBR), '->HpjW))=', sumhjHpi(maxloc(testBR))

                end if

            end do
            deallocate (testBR)
        end if

        if (np(Hplus) > 0) then
            allocate (testBR(np(Hplus)))

            do jj = 1, ndat
                testBR = theo(jj)%BR_Hpjcs &
                         + theo(jj)%BR_Hpjcb &
                         + theo(jj)%BR_Hpjtaunu
                testsumBR = maxval(testBR)

                testsumBR_t = theo(jj)%BR_tWpb &
                              + sum(theo(jj)%BR_tHpjb, dim=1)

                if (testsumBR .gt. 1.0D0 + fuzziness) then
                    if (debug) write (*, *) 'warning: sum of BR for '//trim(adjustl(pdesc(Hplus)%long)) &
                        //' at line number=', jj, 'is', testsumBR
                elseif (testsumBR_t .gt. 1.0D0 + fuzziness) then
                    if (debug) write (*, *) 'warning: sum of BR for the top quark at jj=', jj, 'is', testsumBR_t
                end if

            end do
            deallocate (testBR)
        end if

        do jj = 1, ndat
            theo(jj)%gooddataset = .True.
        end do

        do x = 1, ubound(np, dim=1)
            if (np(x) > 0) then
                do jj = 1, ndat
                    if (minval(theo(jj)%particle(x)%M) .lt. 0.0D0) then
                        theo(jj)%gooddataset = .False.
                        write (*, *) 'warning: negative mass for '//trim(adjustl(pdesc(x)%long)) &
                            //' at line number=', jj, theo(jj)%particle(x)%M
                    elseif (.not. (sum(theo(jj)%particle(x)%M) .ge. 0.0D0)) then
                        theo(jj)%gooddataset = .False.
                        write (*, *) 'warning: mass is NaN for '//trim(adjustl(pdesc(x)%long)) &
                            //' at line number=', jj, theo(jj)%particle(x)%M

                    elseif (minval(theo(jj)%particle(x)%GammaTot) .lt. 0.0D0) then
                        theo(jj)%gooddataset = .False.
                        write (*, *) 'warning: negative total decay width for '//trim(adjustl(pdesc(x)%long))// &
                            ' at line number=', jj, theo(jj)%particle(x)%GammaTot
                    elseif (.not. (sum(theo(jj)%particle(x)%GammaTot) .ge. 0.0D0)) then
                        theo(jj)%gooddataset = .False.
                        write (*, *) 'warning: total decay width is NaN for '//trim(adjustl(pdesc(x)%long))// &
                            ' at line number=', jj, theo(jj)%particle(x)%GammaTot
                    end if
                end do
            end if
        end do

    end subroutine check_dataset

    !> Fills the Standard Model part of theo.
    !! We do this here to save computational time - these  quantities will be
    !! needed a few times in subroutine calcfact_t1, so don't want to calculate them each time
    subroutine fill_theo_SM
        use theory_BRfunctions
        use theory_XS_SM_functions
        use usefulbits, only: ndat
        use S95tables, only: inrange
        implicit none
        ! internal
        integer :: n

        if (np(Hneut) < 1) stop 'error in subroutine fill_theo_SM (np(Hneut))'

        do n = 1, ndat
            call fill_theo_SM_for_datapoint(n)
        end do

    end subroutine fill_theo_SM

    !> Fills the Standard Model part of theo for a single datapoint
    !! We do this here to save computational time - these  quantities will be
    !! needed a few times in subroutine calcfact_t1, so don't want to calculate them each time
    !! @param n datapoint index
    subroutine fill_theo_SM_for_datapoint(n)
        use theory_BRfunctions
        use theory_XS_SM_functions
        use usefulbits, only: theo, small
        use S95tables, only: inrange
        implicit none
        integer, intent(in) :: n
        ! internal
        integer :: i
        double precision :: Mhi

        if (theo(n)%gooddataset) then
            do i = 1, np(Hneut)

                Mhi = theo(n)%particle(Hneut)%M(i)
                if (theo(n)%particle(Hneut)%Mc(i) .ge. small) Mhi = theo(n)%particle(Hneut)%Mc(i)

!   write(*,*) 'DEBUG HB - running fill_theo_SM_for_datapoint for theo_Mh = ', theo(n)%particle(Hneut)%M, &
!   i, theo(n)%particle(Hneut)%M(i)
!   write(*,*) 'DEBUG HB - running fill_theo_SM_for_datapoint for Mh = ', Mhi
!   write(*,*) 'DEBUG HB - running fill_theo_SM_for_datapoint BRs = ', BRSM_HWW(Mhi)

                if (inrange(Mhi, 'SMBR')) then
                    theo(n)%BR_HWW_SM(i) = BRSM_HWW(Mhi)
                    theo(n)%BR_HZZ_SM(i) = BRSM_HZZ(Mhi)
                    theo(n)%BR_Hbb_SM(i) = BRSM_Hbb(Mhi)
                    theo(n)%BR_Htt_SM(i) = BRSM_Htoptop(Mhi)     !HB-5 new
                    theo(n)%BR_Hcc_SM(i) = BRSM_Hcc(Mhi)
                    theo(n)%BR_Hss_SM(i) = BRSM_Hss(Mhi)
                    theo(n)%BR_Hmumu_SM(i) = BRSM_Hmumu(Mhi)
                    theo(n)%BR_Htautau_SM(i) = BRSM_Htautau(Mhi)
                    theo(n)%BR_HZga_SM(i) = BRSM_HZga(Mhi)
                    theo(n)%BR_Hgaga_SM(i) = BRSM_Hgaga(Mhi)
                    theo(n)%BR_Hgg_SM(i) = BRSM_Hgg(Mhi)
                    theo(n)%BR_Hjets_SM(i) = BRSM_Hss(Mhi) + BRSM_Hcc(Mhi) + BRSM_Hbb(Mhi) + BRSM_Hgg(Mhi)
                    theo(n)%particle(Hneut)%GammaTot_SM(i) = BRSM_GammaTot(Mhi)
                else
                    theo(n)%BR_HWW_SM(i) = 0.0D0
                    theo(n)%BR_HZZ_SM(i) = 0.0D0
                    theo(n)%BR_Hbb_SM(i) = 0.0D0
                    theo(n)%BR_Hcc_SM(i) = 0.0D0
                    theo(n)%BR_Hss_SM(i) = 0.0D0
                    theo(n)%BR_Hmumu_SM(i) = 0.0D0
                    theo(n)%BR_Htautau_SM(i) = 0.0D0
                    theo(n)%BR_HZga_SM(i) = 0.0D0
                    theo(n)%BR_Hgaga_SM(i) = 0.0D0
                    theo(n)%BR_Hgg_SM(i) = 0.0D0
                    theo(n)%BR_Hjets_SM(i) = 0.0D0
                    theo(n)%particle(Hneut)%GammaTot_SM(i) = 0.0D0
                end if

                if (inrange(Mhi, 'TEV  ')) then ! n.b.: in fb
!     theo(n)%tev%XS_HZ_SM(i) = XS_tev_HZ_SM(Mhi)
!     theo(n)%tev%XS_HW_SM(i) = XS_tev_HW_SM(Mhi)
                    theo(n)%tev%XS_HZ_SM(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%tev%XS_gg_HZ_SM(i) = ZH_cpmix_nnlo_gg(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%tev%XS_qq_HZ_SM(i) = theo(n)%tev%XS_HZ_SM(i) - theo(n)%tev%XS_gg_HZ_SM(i)
                    theo(n)%tev%XS_HW_SM(i) = WH_nnlo(Mhi, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                    theo(n)%tev%XS_H_SM(i) = XS_tev_gg_H_SM(Mhi) + XS_tev_bb_H_SM(Mhi)
                    theo(n)%tev%XS_gg_H_SM(i) = XS_tev_gg_H_SM(Mhi)     !HB-5 new
                    theo(n)%tev%XS_bb_H_SM(i) = XS_tev_bb_H_SM(Mhi)     !HB-5 new
                    theo(n)%tev%XS_vbf_SM(i) = XS_tev_vbf_SM(Mhi)
                    theo(n)%tev%XS_ttH_SM(i) = XS_tev_ttH_SM(Mhi)
                    theo(n)%tev%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%tev%XS_tH_schan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%tev%XS_tWH_SM(i) = 0.0D0     !HB-5 new

                    theo(n)%tev%XS_Hb_SM(i) = XS_tev_bg_Hb_SM(Mhi)
                    theo(n)%tev%XS_Hb_c1_SM(i) = XS_tev_bg_Hb_c1_SM(Mhi)
                    theo(n)%tev%XS_Hb_c2_SM(i) = XS_tev_bg_Hb_c2_SM(Mhi)
                    theo(n)%tev%XS_Hb_c3_SM(i) = XS_tev_bg_Hb_c3_SM(Mhi)
                    theo(n)%tev%XS_Hb_c4_SM(i) = XS_tev_bg_Hb_c4_SM(Mhi)
                else
                    theo(n)%tev%XS_HW_SM(i) = 0.0D0
                    theo(n)%tev%XS_H_SM(i) = 0.0D0
                    theo(n)%tev%XS_gg_H_SM(i) = 0.0D0
                    theo(n)%tev%XS_bb_H_SM(i) = 0.0D0
                    theo(n)%tev%XS_HZ_SM(i) = 0.0D0
                    theo(n)%tev%XS_gg_HZ_SM(i) = 0.0D0
                    theo(n)%tev%XS_qq_HZ_SM(i) = 0.0D0
                    theo(n)%tev%XS_vbf_SM(i) = 0.0D0
                    theo(n)%tev%XS_ttH_SM(i) = 0.0D0
                    theo(n)%tev%XS_tH_tchan_SM(i) = 0.0D0
                    theo(n)%tev%XS_tH_tchan_SM(i) = 0.0D0
                    theo(n)%tev%XS_tWH_SM(i) = 0.0D0

                    theo(n)%tev%XS_Hb_SM(i) = 0.0D0
                    theo(n)%tev%XS_Hb_c1_SM(i) = 0.0D0
                    theo(n)%tev%XS_Hb_c2_SM(i) = 0.0D0
                    theo(n)%tev%XS_Hb_c3_SM(i) = 0.0D0
                    theo(n)%tev%XS_Hb_c4_SM(i) = 0.0D0
                end if

                if (inrange(Mhi, 'LHC7 ')) then

                    theo(n)%lhc7%XS_H_SM(i) = XS_lhc7_gg_H_SM(Mhi) + XS_lhc7_bb_H_SM(Mhi)
                    theo(n)%lhc7%XS_gg_H_SM(i) = XS_lhc7_gg_H_SM(Mhi)     !HB-5 new
                    theo(n)%lhc7%XS_bb_H_SM(i) = XS_lhc7_bb_H_SM(Mhi)     !HB-5 new

                    theo(n)%lhc7%XS_HZ_SM(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc7%XS_gg_HZ_SM(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc7%XS_qq_HZ_SM(i) = theo(n)%lhc7%XS_HZ_SM(i) - theo(n)%lhc7%XS_gg_HZ_SM(i)
                    theo(n)%lhc7%XS_HW_SM(i) = WH_nnlo(Mhi, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
!     theo(n)%lhc7%XS_HW_SM(i) = XS_lhc7_HW_SM(Mhi)
!     theo(n)%lhc7%XS_HZ_SM(i) = XS_lhc7_HZ_SM(Mhi)
                    theo(n)%lhc7%XS_vbf_SM(i) = XS_lhc7_vbf_SM(Mhi)
                    theo(n)%lhc7%XS_ttH_SM(i) = XS_lhc7_ttH_SM(Mhi)
                    theo(n)%lhc7%XS_tH_tchan_SM(i) = XS_lhc7_tH_tchan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc7%XS_tH_schan_SM(i) = XS_lhc7_tH_schan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc7%XS_tWH_SM(i) = XS_lhc7_gb_tWH_SM(Mhi)
                else
                    theo(n)%lhc7%XS_HW_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_H_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_gg_H_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_bb_H_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_HZ_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_gg_HZ_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_qq_HZ_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_vbf_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_ttH_SM(i) = 0.0D0
                    theo(n)%lhc7%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc7%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc7%XS_tWH_SM(i) = 0.0D0
                end if

                if (inrange(Mhi, 'LHC8 ')) then
!     theo(n)%lhc8%XS_HW_SM(i) = XS_lhc8_HW_SM(Mhi)
                    theo(n)%lhc8%XS_H_SM(i) = XS_lhc8_gg_H_SM(Mhi) + XS_lhc8_bb_H_SM(Mhi)
                    theo(n)%lhc8%XS_gg_H_SM(i) = XS_lhc8_gg_H_SM(Mhi)     !HB-5 new
                    theo(n)%lhc8%XS_bb_H_SM(i) = XS_lhc8_bb_H_SM(Mhi)     !HB-5 new

                    theo(n)%lhc8%XS_HZ_SM(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc8%XS_gg_HZ_SM(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc8%XS_qq_HZ_SM(i) = theo(n)%lhc8%XS_HZ_SM(i) - theo(n)%lhc8%XS_gg_HZ_SM(i)
                    theo(n)%lhc8%XS_HW_SM(i) = WH_nnlo(Mhi, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

!     theo(n)%lhc8%XS_HZ_SM(i) = XS_lhc8_HZ_SM(Mhi)
                    theo(n)%lhc8%XS_vbf_SM(i) = XS_lhc8_vbf_SM(Mhi)
                    theo(n)%lhc8%XS_ttH_SM(i) = XS_lhc8_ttH_SM(Mhi)
                    theo(n)%lhc8%XS_tH_tchan_SM(i) = XS_lhc8_tH_tchan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc8%XS_tH_schan_SM(i) = XS_lhc8_tH_schan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc8%XS_tWH_SM(i) = XS_lhc8_gb_tWH_SM(Mhi)

                else
                    theo(n)%lhc8%XS_HW_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_H_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_HZ_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_gg_HZ_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_qq_HZ_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_vbf_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_ttH_SM(i) = 0.0D0
                    theo(n)%lhc8%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc8%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc8%XS_tWH_SM(i) = 0.0D0
                end if

                if (inrange(Mhi, 'LHC13')) then
!     theo(n)%lhc13%XS_HW_SM(i) = XS_lhc13_HW_SM(Mhi)
                    theo(n)%lhc13%XS_H_SM(i) = XS_lhc13_gg_H_SM(Mhi) + XS_lhc13_bb_H_SM(Mhi)
                    theo(n)%lhc13%XS_gg_H_SM(i) = XS_lhc13_gg_H_SM(Mhi)     !HB-5 new
                    theo(n)%lhc13%XS_bb_H_SM(i) = XS_lhc13_bb_H_SM(Mhi)     !HB-5 new
                    theo(n)%lhc13%XS_HZ_SM(i) = ZH_cpmix_nnlo_ggqqbb(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc13%XS_gg_HZ_SM(i) = ZH_cpmix_nnlo_gg(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    theo(n)%lhc13%XS_qq_HZ_SM(i) = theo(n)%lhc13%XS_HZ_SM(i) - theo(n)%lhc13%XS_gg_HZ_SM(i)
                    theo(n)%lhc13%XS_HW_SM(i) = WH_nnlo(Mhi, 'LHC13', 1.0D0, 1.0D0, 1.0D0, .True., .True.)

!     theo(n)%lhc13%XS_HZ_SM(i) = XS_lhc13_HZ_SM(Mhi)
                    theo(n)%lhc13%XS_vbf_SM(i) = XS_lhc13_vbf_SM(Mhi)
                    theo(n)%lhc13%XS_ttH_SM(i) = XS_lhc13_ttH_SM(Mhi)
                    theo(n)%lhc13%XS_tH_tchan_SM(i) = XS_lhc13_tH_tchan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc13%XS_tH_schan_SM(i) = XS_lhc13_tH_schan_SM(Mhi)     !HB-5 new
                    theo(n)%lhc13%XS_tWH_SM(i) = XS_lhc13_gb_tWH_SM(Mhi)

                else
                    theo(n)%lhc13%XS_HW_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_H_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_HZ_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_gg_HZ_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_qq_HZ_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_vbf_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_ttH_SM(i) = 0.0D0
                    theo(n)%lhc13%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc13%XS_tH_tchan_SM(i) = 0.0D0     !HB-5 new
                    theo(n)%lhc13%XS_tWH_SM(i) = 0.0D0
                end if

            end do
        end if

    end subroutine fill_theo_SM_for_datapoint

end module theo_manip
