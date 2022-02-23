module evaluate
    implicit none

contains

    subroutine evaluate_model(t, n, signalScaling)
        ! internal routine
        !------------------------------------------------------------
        ! This subroutine evaluates the signal strength modifier for every Higgs boson and
        ! considered analysis. It fills a matrix neutHiggs(:,:) of type neutHiggs with dimensions
        ! (number(considered analyses),nH).
        !------------------------------------------------------------
        use usefulbits, only: np, Hneut, dataset, results, vsmall
        use usefulbits_hs, only: neutHiggses, nanalys, HSresults, obs, analyses, &
                                 iterations, deallocate_covariance_matrices, &
                                 output_level, Nparam, nanalys, HSres
        use datatables, only: setup_tablelist, check_available_Higgses
        use pc_chisq
        use mc_chisq
        use all_chisq
        use numerics
        implicit none
        !--------------------------------------input
        type(dataset), intent(in) :: t
        integer, intent(in) :: n
        double precision, intent(in), optional :: signalScaling
        !-------------------------------------output
        !  type(HSresults), intent(out) :: r

        integer :: ii, jj, iii, jjj

        double precision :: totchisq, muchisq, mhchisq
        integer :: nobs, Nmu, Nmh
        character(LEN=100), allocatable :: assignmentgroups(:)
        integer, allocatable :: assignmentgroups_domH(:)
        integer, allocatable :: assignmentgroups_Higgs_comb(:, :)

        allocate (assignmentgroups(nanalys), assignmentgroups_domH(nanalys))
        allocate (assignmentgroups_Higgs_comb(nanalys, np(Hneut)))

        assignmentgroups = ''

        !---Initialize assignmentgroups arrays with default values
        do ii = lbound(assignmentgroups_domH, dim=1), ubound(assignmentgroups_domH, dim=1)
            assignmentgroups_domH(ii) = 0
            assignmentgroups_Higgs_comb(ii, :) = 0
        end do

        !---First, evaluate the model predictions
        allocate (neutHiggses(nanalys, np(Hneut)))
        !-Loop over considered analyses
        do ii = lbound(neutHiggses, dim=1), ubound(neutHiggses, dim=1)
            !-Loop over the neutral Higgs bosons of the model
            do jj = lbound(neutHiggses, dim=2), ubound(neutHiggses, dim=2)
                !!   write(*,*) "hello evaluate model:", ii, jj
                call calc_mupred(jj, t, obs(ii)%table, neutHiggses(ii, jj))
            end do
            if (.not. allocated(obs(ii)%Higgses)) allocate (obs(ii)%Higgses(np(Hneut)))
            obs(ii)%Higgses(:) = neutHiggses(ii, :)
        end do

        !-Pass the observables and their predicted Higgs properties (obs%Higgses)
        !-to the tablelist "analyses"
        call setup_tablelist

        jjj = 0
        do ii = lbound(analyses, dim=1), ubound(analyses, dim=1)
            call deallocate_covariance_matrices
            call assign_Higgs_to_peaks(analyses(ii)%table, analyses(ii)%peaks, 0)
            do iii = lbound(analyses(ii)%peaks, dim=1), ubound(analyses(ii)%peaks, dim=1)
                if (analyses(ii)%table%mhchisq .eq. 1 .and. &
                    len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))) .ne. 0) then
                    jjj = jjj + 1
                    assignmentgroups(jjj) = analyses(ii)%peaks(iii)%assignmentgroup
                    assignmentgroups_Higgs_comb(jjj, :) = analyses(ii)%peaks(iii)%Higgs_comb
                    assignmentgroups_domH(jjj) = analyses(ii)%peaks(iii)%domH
                    !     write(*,*) "Found leader of group ",assignmentgroups(jjj)
                    !     write(*,*) "ID ",analyses(ii)%peaks(iii)%id
                    !     write(*,*) "with Higgs combination ",assignmentgroups_Higgs_comb(jjj,:)
                    !     write(*,*) "and dominant Higgs boson ",assignmentgroups_domH(jjj)
                end if
            end do
        end do
        do ii = lbound(analyses, dim=1), ubound(analyses, dim=1)
            do iii = lbound(analyses(ii)%peaks, dim=1), ubound(analyses(ii)%peaks, dim=1)
                if (analyses(ii)%table%mhchisq .eq. 0 .and. &
                    len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))) .ne. 0) then
                    !SELECT ASSIGNMENT GROUP FOLLOWERS
                    do jjj = lbound(assignmentgroups, dim=1), ubound(assignmentgroups, dim=1)
                        if (analyses(ii)%peaks(iii)%assignmentgroup .eq. assignmentgroups(jjj)) then
                            !TAKE OVER THE HIGGS ASSIGNMENT OF THE LEADING PEAK
                            analyses(ii)%peaks(iii)%Higgs_comb = assignmentgroups_Higgs_comb(jjj, :)
                            analyses(ii)%peaks(iii)%domH = assignmentgroups_domH(jjj)
                            if (assignmentgroups_domH(jjj) .ne. 0) then
                                analyses(ii)%peaks(iii)%Higgs_assignment_forced = 1
                            end if
                            call evaluate_peak(analyses(ii)%peaks(iii), analyses(ii)%table)
                        end if
                    end do
                end if
            end do
        end do

        !   write(*,*) "Starting assignment procedure..."

        ! Do the iterative Higgs-to-peak-assignment here:
        call assign_Higgs_to_peaks_with_correlations(iterations)

        if(present(signalScaling)) then
            call calculate_total_pc_chisq(totchisq, muchisq, mhchisq, nobs, Nmu, Nmh,signalScaling)
        else
            call calculate_total_pc_chisq(totchisq, muchisq, mhchisq, nobs, Nmu, Nmh)
        endif
        !   write(*,*) "...done."

        if (output_level .eq. 1) call print_peakinformation
        if (output_level .eq. 2) call print_peakinformation_essentials
        if (output_level .eq. 3) then
            call print_peaks_to_file
            call print_peaks_signal_rates_to_file
        end if

        call add_peaks_to_HSresults(HSres(n))

        HSres(n)%Chisq_peak = totchisq
        HSres(n)%Chisq_peak_mu = muchisq
        HSres(n)%Chisq_mpred = 0.0D0
        HSres(n)%Chisq_peak_mu = muchisq
        HSres(n)%Chisq_peak_mh = mhchisq
        HSres(n)%nobs_mpred = 0
        HSres(n)%nobs_peak_mu = Nmu
        HSres(n)%nobs_peak_mh = Nmh
        HSres(n)%nanalysis = size(analyses)
        HSres(n)%nobs_peak = nobs
        !
        !   if(HSres(n)%Chisq.gt.vsmall.and.(HSres(n)%nobs-Nparam).gt.0) then
        !    HSres(n)%Pvalue_peak = 1 - gammp(dble(HSres(n)%nobs-Nparam)/2,HSres(n)%Chisq/2)
        !   endif

        if (HSres(n)%Chisq_peak .gt. vsmall .and. (HSres(n)%nobs_peak - Nparam) .gt. 0) then
            HSres(n)%Pvalue_peak = 1 - gammp(dble(HSres(n)%nobs_peak - Nparam)/2, HSres(n)%Chisq_peak/2)
        end if

        deallocate (neutHiggses)
        deallocate (assignmentgroups, assignmentgroups_domH, assignmentgroups_Higgs_comb)
    end subroutine evaluate_model

    subroutine calc_mupred(j, t, mutab, Higgs)
        ! internal routine
        ! Calculates the model-predicted signal strength modifier
        !------------------------------------------------------------
        use usefulbits, only: dataset, div
        use usefulbits_HS, only: neutHiggs, mutable, useSMtest, eps
        implicit none

        integer, intent(in) :: j                                        ! Higgs index
        type(dataset), intent(in) :: t
        type(mutable), intent(inout) :: mutab
        type(neutHiggs), intent(inout) :: Higgs

        integer :: i
        double precision :: c, dcbyc
        integer :: testSMratios
        logical :: correct_properties

        Higgs%m = t%particle(mutab%particle_x)%M(j)
        Higgs%dm = t%particle(mutab%particle_x)%dM(j)
        Higgs%id = j

        call get_channelrates(j, t, mutab)

        correct_properties = .True.

        !--Evaluate the predicted signal strength modifier c of the model
        c = 0.
        do i = 1, mutab%Nc
            !----use a weighted average of the channel rate ratios
            c = c + mutab%channel_w(i, j)*mutab%channel_mu(i, j)
        end do

        !--Evaluate the deviation of each channel rate ratio to the signal
        !--strength modifier c and test SM likeness criterium, if this is
        !--activated.
        testSMratios = 1  !passes the SM-like ratios test
        do i = 1, mutab%Nc
            dcbyc = div((mutab%channel_mu(i, j) - c), c, 0.0D0, 1.0D9)
            if (dcbyc*mutab%channel_w(i, j) .gt. eps .and. useSMtest) then
                testSMratios = -1  !fails the SM-like ratios test
            end if
        end do

        if (testSMratios .lt. 0) correct_properties = .False.

        if (correct_properties) then
            Higgs%mu = c
        else
            Higgs%mu = 0.0D0
        end if

    end subroutine calc_mupred

    subroutine get_channelrates(j, t, mutab)
        ! internal routine
        !
        ! This subroutine assignes the rates, weights and systematic rate uncertainty of
        ! the Higgs boson (j) for the channels considered by the analysis (mutab).
        !
        ! WARNING: if normalize_rates_to_reference_position is true
        ! The rates are normalized w.r.t. a reference rate at the (peak) mass position.
        ! This does not work with the mass-centered chi^2 method.
        ! Also, theoretical mass uncertainties are problematic!
        !------------------------------------------------------------
        use usefulbits, only: dataset, div, small
        use usefulbits_HS, only: neutHiggs, mutable, delta_rate, normalize_rates_to_reference_position, &
                                 normalize_rates_to_reference_position_outside_dmtheo
        use theory_XS_SM_functions
        use theory_BRfunctions

        integer, intent(in) :: j
        type(dataset), intent(in) :: t
        type(mutable), intent(inout) :: mutab

        integer :: i, p, d ! id
        integer :: ii, p1, p2, d1, d2 !id1, id2
        double precision :: rate, SMrate, drsq_SM, drsq, dBR, dBRSM, drcov, drcovSM

        double precision :: rate_SMref, refmass, BR_SMref!,BR_SMref_mpeak
        double precision :: dynamicalmass, rate_SMdyn

        ! TS (17/10/2018: dynamicalmass is the default reference mass position for the SM normalization)

        if (size(mutab%mass, dim=1) .eq. 1) then
            refmass = mutab%mass(1)

            ! TS (17/10/2018): Take dynamical reference mass for SM-normalization at mobs+dmtheo box boundary.
            if (t%particle(mutab%particle_x)%M(j) .ge. (mutab%mass(1) + t%particle(mutab%particle_x)%dM(j))) then
                dynamicalmass = mutab%mass(1) + t%particle(mutab%particle_x)%dM(j)
            else if (t%particle(mutab%particle_x)%M(j) .le. (mutab%mass(1) - t%particle(mutab%particle_x)%dM(j))) then
                dynamicalmass = mutab%mass(1) - t%particle(mutab%particle_x)%dM(j)
            else
                dynamicalmass = t%particle(mutab%particle_x)%M(j)
            end if

            !   write(*,*) "HS debug, dynamicalmass, refmass = ",dynamicalmass, refmass
            !---

        else
            ! Only relevant for the mass-centered chi^2 method
            refmass = t%particle(mutab%particle_x)%M(j)
        end if

        !write(*,*) 'DEBUG HS: id = ', mutab%id
        !write(*,*) 'DEBUG HS, m = ', t%particle(mutab%particle_x)%M(j)

        do i = 1, mutab%Nc
            !   id = mutab%channel_id(i)
            !   p = int((id-modulo(id,10))/dble(10))
            !   d = modulo(id,10)
            p = mutab%channel_p_id(i)
            d = mutab%channel_d_id(i)
            !--Do the production rate for the relevant experiment and cms-energy
            if (mutab%collider .eq. 'LHC') then
                if (abs(mutab%energy - 7.0D0) .le. small) then
                    if (p .eq. 1) then
                        rate = t%lhc7%XS_hj_ratio(j)
                        SMrate = t%lhc7%XS_H_SM(j)
                        rate_SMdyn = XS_lhc7_gg_H_SM(dynamicalmass) + XS_lhc7_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_gg_H_SM(refmass) + XS_lhc7_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'singleH'
                    else if (p .eq. 2) then
                        rate = t%lhc7%XS_vbf_ratio(j)
                        SMrate = t%lhc7%XS_vbf_SM(j)
                        rate_SMdyn = XS_lhc7_vbf_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_vbf_SM(refmass)
                        mutab%channel_description(i, 1) = 'VBF'
                    else if (p .eq. 3) then
                        rate = t%lhc7%XS_hjW_ratio(j)
                        SMrate = t%lhc7%XS_HW_SM(j)
                        rate_SMdyn = WH_nnlo(dynamicalmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        rate_SMref = WH_nnlo(refmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        mutab%channel_description(i, 1) = 'HW'
                    else if (p .eq. 4) then
                        rate = t%lhc7%XS_hjZ_ratio(j)
                        SMrate = t%lhc7%XS_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_ggqqbb(dynamicalmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_ggqqbb(refmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'HZ'
                    else if (p .eq. 5) then
                        rate = t%lhc7%XS_tthj_ratio(j)
                        SMrate = t%lhc7%XS_ttH_SM(j)
                        rate_SMdyn = XS_lhc7_ttH_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_ttH_SM(refmass)
                        mutab%channel_description(i, 1) = 'ttH'
                    else if (p .eq. 6) then
                        rate = t%lhc7%XS_gg_hj_ratio(j)
                        SMrate = t%lhc7%XS_gg_H_SM(j)
                        rate_SMdyn = XS_lhc7_gg_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_gg_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'ggH'
                    else if (p .eq. 7) then
                        rate = t%lhc7%XS_bb_hj_ratio(j)
                        SMrate = t%lhc7%XS_bb_H_SM(j)
                        rate_SMdyn = XS_lhc7_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'bbH'
                    else if (p .eq. 8) then
                        rate = t%lhc7%XS_thj_tchan_ratio(j)
                        SMrate = t%lhc7%XS_tH_tchan_SM(j)
                        rate_SMdyn = XS_lhc7_tH_tchan_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_tH_tchan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH t-chan'
                    else if (p .eq. 9) then
                        rate = t%lhc7%XS_thj_schan_ratio(j)
                        SMrate = t%lhc7%XS_tH_schan_SM(j)
                        rate_SMdyn = XS_lhc7_tH_schan_SM(dynamicalmass)
                        rate_SMref = XS_lhc7_tH_schan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH s-chan'
                    else if (p .eq. 10) then
                        rate = t%lhc7%XS_qq_hjZ_ratio(j)
                        SMrate = t%lhc7%XS_qq_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_qqbb(dynamicalmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_qqbb(refmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'qq-HZ'
                    else if (p .eq. 11) then
                        rate = t%lhc7%XS_gg_hjZ_ratio(j)
                        SMrate = t%lhc7%XS_gg_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_gg(dynamicalmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_gg(refmass, 'LHC7 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'gg-HZ'
                    else if (p .eq. 12) then
                        rate = t%lhc7%XS_tWhj_ratio(j)
                        SMrate = t%lhc7%XS_tWH_SM(j)
                        rate_SMdyn = t%lhc7%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        rate_SMref = t%lhc7%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        mutab%channel_description(i, 1) = 'tWH'
                    else if (p .eq. 0) then
                        rate = 1.0D0
                        SMrate = 1.0D0
                        rate_SMdyn = 1.0D0
                        rate_SMref = 1.0D0
                        mutab%channel_description(i, 1) = 'none'
                    end if
                else if (abs(mutab%energy - 8.0D0) .le. small) then
                    if (p .eq. 1) then
                        rate = t%lhc8%XS_hj_ratio(j)
                        SMrate = t%lhc8%XS_H_SM(j)
                        rate_SMdyn = XS_lhc8_gg_H_SM(dynamicalmass) + XS_lhc8_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_gg_H_SM(refmass) + XS_lhc8_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'singleH'
                    else if (p .eq. 2) then
                        rate = t%lhc8%XS_vbf_ratio(j)
                        SMrate = t%lhc8%XS_vbf_SM(j)
                        rate_SMdyn = XS_lhc8_vbf_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_vbf_SM(refmass)
                        mutab%channel_description(i, 1) = 'VBF'
                    else if (p .eq. 3) then
                        rate = t%lhc8%XS_hjW_ratio(j)
                        SMrate = t%lhc8%XS_HW_SM(j)
                        rate_SMdyn = WH_nnlo(dynamicalmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        rate_SMref = WH_nnlo(refmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        mutab%channel_description(i, 1) = 'HW'
                    else if (p .eq. 4) then
                        rate = t%lhc8%XS_hjZ_ratio(j)
                        SMrate = t%lhc8%XS_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_ggqqbb(dynamicalmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_ggqqbb(refmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'HZ'
                    else if (p .eq. 5) then
                        rate = t%lhc8%XS_tthj_ratio(j)
                        SMrate = t%lhc8%XS_ttH_SM(j)
                        rate_SMdyn = XS_lhc8_ttH_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_ttH_SM(refmass)
                        mutab%channel_description(i, 1) = 'ttH'
                    else if (p .eq. 6) then
                        rate = t%lhc8%XS_gg_hj_ratio(j)
                        SMrate = t%lhc8%XS_gg_H_SM(j)
                        rate_SMdyn = XS_lhc8_gg_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_gg_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'ggH'
                    else if (p .eq. 7) then
                        rate = t%lhc8%XS_bb_hj_ratio(j)
                        SMrate = t%lhc8%XS_bb_H_SM(j)
                        rate_SMdyn = XS_lhc8_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'bbH'
                    else if (p .eq. 8) then
                        rate = t%lhc8%XS_thj_tchan_ratio(j)
                        SMrate = t%lhc8%XS_tH_tchan_SM(j)
                        rate_SMdyn = XS_lhc8_tH_tchan_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_tH_tchan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH t-chan'
                    else if (p .eq. 9) then
                        rate = t%lhc8%XS_thj_schan_ratio(j)
                        SMrate = t%lhc8%XS_tH_schan_SM(j)
                        rate_SMdyn = XS_lhc8_tH_schan_SM(dynamicalmass)
                        rate_SMref = XS_lhc8_tH_schan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH s-chan'
                    else if (p .eq. 10) then
                        rate = t%lhc8%XS_qq_hjZ_ratio(j)
                        SMrate = t%lhc8%XS_qq_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_qqbb(dynamicalmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_qqbb(refmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'qq-HZ'
                    else if (p .eq. 11) then
                        rate = t%lhc8%XS_gg_hjZ_ratio(j)
                        SMrate = t%lhc8%XS_gg_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_gg(dynamicalmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_gg(refmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'gg-HZ'
                    else if (p .eq. 12) then
                        rate = t%lhc8%XS_tWhj_ratio(j)
                        SMrate = t%lhc8%XS_tWH_SM(j)
                        rate_SMdyn = t%lhc8%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        rate_SMref = t%lhc8%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        mutab%channel_description(i, 1) = 'tWH'
                    else if (p .eq. 0) then
                        rate = 1.0D0
                        SMrate = 1.0D0
                        rate_SMdyn = 1.0D0
                        rate_SMref = 1.0D0
                        mutab%channel_description(i, 1) = 'none'
                    end if
                else if (abs(mutab%energy - 13.0D0) .le. small) then
                    if (p .eq. 1) then
                        rate = t%lhc13%XS_hj_ratio(j)
                        SMrate = t%lhc13%XS_H_SM(j)
                        rate_SMdyn = XS_lhc13_gg_H_SM(dynamicalmass) + XS_lhc13_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_gg_H_SM(refmass) + XS_lhc13_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'singleH'
                    else if (p .eq. 2) then
                        rate = t%lhc13%XS_vbf_ratio(j)
                        SMrate = t%lhc13%XS_vbf_SM(j)
                        rate_SMdyn = XS_lhc13_vbf_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_vbf_SM(refmass)
                        mutab%channel_description(i, 1) = 'VBF'
                    else if (p .eq. 3) then
                        rate = t%lhc13%XS_hjW_ratio(j)
                        SMrate = t%lhc13%XS_HW_SM(j)
                        rate_SMdyn = WH_nnlo(dynamicalmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        rate_SMref = WH_nnlo(refmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                        mutab%channel_description(i, 1) = 'HW'
                    else if (p .eq. 4) then
                        rate = t%lhc13%XS_hjZ_ratio(j)
                        SMrate = t%lhc13%XS_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_ggqqbb(dynamicalmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_ggqqbb(refmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'HZ'
                    else if (p .eq. 5) then
                        rate = t%lhc13%XS_tthj_ratio(j)
                        SMrate = t%lhc13%XS_ttH_SM(j)
                        rate_SMdyn = XS_lhc13_ttH_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_ttH_SM(refmass)
                        mutab%channel_description(i, 1) = 'ttH'
                    else if (p .eq. 6) then
                        rate = t%lhc13%XS_gg_hj_ratio(j)
                        SMrate = t%lhc13%XS_gg_H_SM(j)
                        rate_SMdyn = XS_lhc13_gg_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_gg_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'ggH'
                    else if (p .eq. 7) then
                        rate = t%lhc13%XS_bb_hj_ratio(j)
                        SMrate = t%lhc13%XS_bb_H_SM(j)
                        rate_SMdyn = XS_lhc13_bb_H_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_bb_H_SM(refmass)
                        mutab%channel_description(i, 1) = 'bbH'
                    else if (p .eq. 8) then
                        rate = t%lhc13%XS_thj_tchan_ratio(j)
                        SMrate = t%lhc13%XS_tH_tchan_SM(j)
                        rate_SMdyn = XS_lhc13_tH_tchan_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_tH_tchan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH t-chan'
                    else if (p .eq. 9) then
                        rate = t%lhc13%XS_thj_schan_ratio(j)
                        SMrate = t%lhc13%XS_tH_schan_SM(j)
                        rate_SMdyn = XS_lhc13_tH_schan_SM(dynamicalmass)
                        rate_SMref = XS_lhc13_tH_schan_SM(refmass)
                        mutab%channel_description(i, 1) = 'tH s-chan'
                    else if (p .eq. 10) then
                        rate = t%lhc13%XS_qq_hjZ_ratio(j)
                        SMrate = t%lhc13%XS_qq_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_qqbb(dynamicalmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_qqbb(refmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'qq-HZ'
                    else if (p .eq. 11) then
                        rate = t%lhc13%XS_gg_hjZ_ratio(j)
                        SMrate = t%lhc13%XS_gg_HZ_SM(j)
                        rate_SMdyn = ZH_cpmix_nnlo_gg(dynamicalmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        rate_SMref = ZH_cpmix_nnlo_gg(refmass, 'LHC13', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                        mutab%channel_description(i, 1) = 'gg-HZ'
                    else if (p .eq. 12) then
                        rate = t%lhc13%XS_tWhj_ratio(j)
                        SMrate = t%lhc13%XS_tWH_SM(j)
                        rate_SMdyn = t%lhc13%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        rate_SMref = t%lhc13%XS_tWH_SM(j) ! Implemented function is not mass-dependent
                        ! write(*,*) "get_channelrates: tWH rate, SMrate: ", rate, SMrate
                        mutab%channel_description(i, 1) = 'tWH'
                    else if (p .eq. 0) then
                        rate = 1.0D0
                        SMrate = 1.0D0
                        rate_SMdyn = 1.0D0
                        rate_SMref = 1.0D0
                        mutab%channel_description(i, 1) = 'none'
                    end if
                end if
            else if (mutab%collider .eq. 'TEV') then
                if (p .eq. 1) then
                    rate = t%tev%XS_hj_ratio(j)
                    SMrate = t%tev%XS_H_SM(j)
                    rate_SMdyn = XS_tev_gg_H_SM(dynamicalmass) + XS_tev_bb_H_SM(dynamicalmass)
                    rate_SMref = XS_tev_gg_H_SM(refmass) + XS_tev_bb_H_SM(refmass)
                    mutab%channel_description(i, 1) = 'singleH'
                else if (p .eq. 2) then
                    rate = t%tev%XS_vbf_ratio(j)
                    SMrate = t%tev%XS_vbf_SM(j)
                    rate_SMdyn = XS_tev_vbf_SM(dynamicalmass)
                    rate_SMref = XS_tev_vbf_SM(refmass)
                    mutab%channel_description(i, 1) = 'VBF'
                else if (p .eq. 3) then
                    rate = t%tev%XS_hjW_ratio(j)
                    SMrate = t%tev%XS_HW_SM(j)
                    rate_SMdyn = WH_nnlo(dynamicalmass, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                    rate_SMref = WH_nnlo(refmass, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                    mutab%channel_description(i, 1) = 'HW'
                else if (p .eq. 4) then
                    rate = t%tev%XS_hjZ_ratio(j)
                    SMrate = t%tev%XS_HZ_SM(j)
                    rate_SMdyn = ZH_cpmix_nnlo_ggqqbb(dynamicalmass, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    rate_SMref = ZH_cpmix_nnlo_ggqqbb(refmass, 'TEV  ', 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, .True.)
                    mutab%channel_description(i, 1) = 'HZ'
                else if (p .eq. 5) then
                    rate = t%tev%XS_tthj_ratio(j)
                    SMrate = t%tev%XS_ttH_SM(j)
                    rate_SMdyn = XS_tev_ttH_SM(dynamicalmass)
                    rate_SMref = XS_tev_ttH_SM(refmass)
                    mutab%channel_description(i, 1) = 'ttH'
                else if (p .eq. 0) then
                    rate = 1.0D0
                    SMrate = 1.0D0
                    rate_SMdyn = 1.0D0
                    rate_SMref = 1.0D0
                    mutab%channel_description(i, 1) = 'none'
                end if
            else if (mutab%collider .eq. 'ILC') then
                !--n.B.: As a first attempt, we use the LHC8 normalized cross sections for ZH, VBF, ttH.
                !        In order to do this properly, a separate input for the ILC cross sections
                !        has to be provided! It works only for single production mode observables (no
                !        correct weighting of channels included!)Then, at least in the effective coupling
                !        approximation, there is no difference to a full implementation.
                !        The theoretical uncertainty of the ILC production modes will are defined in
                !        usefulbits_HS.f90.
                if (p .eq. 1 .or. p .eq. 2) then
                    write (*, *) 'Warning: Unknown ILC production mode (', p, ') in table ', mutab%id
                    rate = 0.0D0
                    SMrate = 1.0D0
                    rate_SMref = 1.0D0
                    mutab%channel_description(i, 1) = 'unknown'
                else if (p .eq. 3) then
                    rate = t%lhc8%XS_hjW_ratio(j)
                    SMrate = t%lhc8%XS_HW_SM(j)
                    rate_SMdyn = WH_nnlo(dynamicalmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                    rate_SMref = WH_nnlo(refmass, 'LHC8 ', 1.0D0, 1.0D0, 1.0D0, .True., .True.)
                    mutab%channel_description(i, 1) = 'WBF'
                else if (p .eq. 4) then
                    rate = t%lhc8%XS_hjZ_ratio(j)
                    SMrate = t%lhc8%XS_HZ_SM(j)
                    rate_SMdyn = XS_lhc8_HZ_SM(dynamicalmass)
                    rate_SMref = XS_lhc8_HZ_SM(refmass)
                    mutab%channel_description(i, 1) = 'HZ'
                else if (p .eq. 5) then
                    rate = t%lhc8%XS_tthj_ratio(j)
                    SMrate = t%lhc8%XS_ttH_SM(j)
                    rate_SMdyn = XS_lhc8_ttH_SM(dynamicalmass)
                    rate_SMref = XS_lhc8_ttH_SM(refmass)
                    mutab%channel_description(i, 1) = 'ttH'
                else if (p .eq. 0) then
                    rate = 1.0D0
                    SMrate = 1.0D0
                    rate_SMdyn = 1.0D0
                    rate_SMref = 1.0D0
                    mutab%channel_description(i, 1) = 'none'
                end if
            end if

            !    write(*,*) "DEBUG, after production mode, rate_SMdyn = ",rate_SMdyn
            !--Multiply now by the decay rate
            if (d .eq. 1) then
                rate = rate*div(t%BR_hjgaga(j), t%BR_Hgaga_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hgaga_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hgaga(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hgaga(refmass)
                mutab%channel_description(i, 2) = 'gammagamma'
            else if (d .eq. 2) then
                rate = rate*div(t%BR_hjWW(j), t%BR_HWW_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_HWW_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_HWW(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_HWW(refmass)
                mutab%channel_description(i, 2) = 'WW'
            else if (d .eq. 3) then
                rate = rate*div(t%BR_hjZZ(j), t%BR_HZZ_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_HZZ_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_HZZ(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_HZZ(refmass)
                mutab%channel_description(i, 2) = 'ZZ'
            else if (d .eq. 4) then
                rate = rate*div(t%BR_hjtautau(j), t%BR_Htautau_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Htautau_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Htautau(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Htautau(refmass)
                mutab%channel_description(i, 2) = 'tautau'
            else if (d .eq. 5) then
                rate = rate*div(t%BR_hjbb(j), t%BR_Hbb_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hbb_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hbb(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hbb(refmass)
                mutab%channel_description(i, 2) = 'bb'
            else if (d .eq. 6) then
                rate = rate*div(t%BR_hjZga(j), t%BR_HZga_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_HZga_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_HZga(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_HZga(refmass)
                mutab%channel_description(i, 2) = 'Zgamma'
            else if (d .eq. 7) then
                rate = rate*div(t%BR_hjcc(j), t%BR_Hcc_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hcc_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hcc(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hcc(refmass)
                mutab%channel_description(i, 2) = 'cc'
            else if (d .eq. 8) then
                rate = rate*div(t%BR_hjmumu(j), t%BR_Hmumu_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hmumu_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hmumu(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hmumu(refmass)
                mutab%channel_description(i, 2) = 'mumu'
            else if (d .eq. 9) then
                rate = rate*div(t%BR_hjgg(j), t%BR_Hgg_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hgg_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hgg(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hgg(refmass)
                mutab%channel_description(i, 2) = 'gg'
            else if (d .eq. 10) then
                rate = rate*div(t%BR_hjss(j), t%BR_Hss_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Hss_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Hss(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Hss(refmass)
                mutab%channel_description(i, 2) = 'ss'
            else if (d .eq. 11) then
                rate = rate*div(t%BR_hjtt(j), t%BR_Htt_SM(j), 0.0D0, 1.0D0)
                SMrate = SMrate*t%BR_Htt_SM(j)
                rate_SMdyn = rate_SMdyn*BRSM_Htoptop(dynamicalmass)
                rate_SMref = rate_SMref*BRSM_Htoptop(refmass)
                mutab%channel_description(i, 2) = 'tt'
            else if (d .eq. 0) then
                rate = rate*1.0D0
                SMrate = SMrate*1.0D0
                rate_SMdyn = rate_SMdyn*1.0D0
                rate_SMref = rate_SMref*1.0D0
                mutab%channel_description(i, 2) = 'none'
            end if

            !-------------------------
            ! NEW FEATURE (since HB-5.2): Enable to set channelrates directly.
            if (p .ne. 0 .and. d .ne. 0 .and. p .ne. 12) then ! Need to extend channelrates to tWH production (p==12)
                select case (d)
                case (1)
                    BR_SMref = t%BR_Hgaga_SM(j)
                    !      BR_SMref_mpeak = BRSM_Hgaga(refmass)
                case (2)
                    BR_SMref = t%BR_HWW_SM(j)
                    !      BR_SMref_mpeak = BRSM_HWW(refmass)
                case (3)
                    BR_SMref = t%BR_HZZ_SM(j)
                    !      BR_SMref_mpeak = BRSM_HZZ(refmass)
                case (4)
                    BR_SMref = t%BR_Htautau_SM(j)
                    !      BR_SMref_mpeak = BRSM_Htautau(refmass)
                case (5)
                    BR_SMref = t%BR_Hbb_SM(j)
                    !      BR_SMref_mpeak = BRSM_Hbb(refmass)
                case (6)
                    BR_SMref = t%BR_HZga_SM(j)
                    !      BR_SMref_mpeak = BRSM_HZga(refmass)
                case (7)
                    BR_SMref = t%BR_Hcc_SM(j)
                    !      BR_SMref_mpeak = BRSM_Hcc(refmass)
                case (8)
                    BR_SMref = t%BR_Hmumu_SM(j)
                    !      BR_SMref_mpeak = BRSM_Hmumu(refmass)
                case (9)
                    BR_SMref = t%BR_Hgg_SM(j)
                    !      BR_SMref_mpeak = BRSM_Hgg(refmass)
                case (10)
                    BR_SMref = t%BR_Hss_SM(j)
                case (11)
                    BR_SMref = t%BR_Htt_SM(j)
                end select
                if (mutab%collider .eq. 'LHC') then
                    if (abs(mutab%energy - 7.0D0) .le. small) then
                        if (t%lhc7%channelrates(j, p, d) .ge. 0.0d0) then
                            rate = div(t%lhc7%channelrates(j, p, d), BR_SMref, 0.0D0, 1.0D0)
                        end if
                    else if (abs(mutab%energy - 8.0D0) .le. small) then
                        if (t%lhc8%channelrates(j, p, d) .ge. 0.0d0) then
                            rate = div(t%lhc8%channelrates(j, p, d), BR_SMref, 0.0D0, 1.0D0)
                        end if
                    else if (abs(mutab%energy - 13.0D0) .le. small) then
                        if (t%lhc13%channelrates(j, p, d) .ge. 0.0d0) then
                            rate = div(t%lhc13%channelrates(j, p, d), BR_SMref, 0.0D0, 1.0D0)
                        end if
                    end if
                else if (mutab%collider .eq. 'TEV') then
                    if (t%tev%channelrates(j, p, d) .ge. 0.0d0) then
                        rate = div(t%tev%channelrates(j, p, d), BR_SMref, 0.0D0, 1.0D0)
                    end if
                end if
            end if
            !-------------------------

            ! write(*,*) 'DEBUG HS: SM BRs = ', t%BR_HWW_SM(j), t%BR_HZZ_SM(j), t%BR_Hgaga_SM(j)
            ! write(*,*) 'DEBUG HS: rate, SMrate(i), rate_SMdyn = ', rate, SMrate, rate_SMdyn
            ! write(*,*) 'DEBUG HS: eff(i) = ', mutab%channel_eff(i)

            if (normalize_rates_to_reference_position) then
                !! THIS IS STILL IN TESTING PHASE !!
                mutab%channel_mu(i, j) = rate*div(SMrate, rate_SMref, 0.0D0, 0.0D0)
                !   write(*,*) "DEBUG HS: normalize_rates_to_reference_position, mu = ",mutab%channel_mu(i,j)
            else
                !   mutab%channel_mu(i,j)=rate  !! OLD WAY
                mutab%channel_mu(i, j) = rate*div(SMrate, rate_SMdyn, 0.0D0, 0.0D0)  !! NEW WAY  (TS 17/10/2018)
                !   write(*,*) "DEBUG HS: not normalize_rates_to_reference_position, mu = ",mutab%channel_mu(i,j)
            end if

            if (normalize_rates_to_reference_position_outside_dmtheo) then
                if (abs(refmass - t%particle(mutab%particle_x)%M(j)) .ge. t%particle(mutab%particle_x)%dM(j)) then
                    mutab%channel_mu(i, j) = rate*div(SMrate, rate_SMref, 0.0D0, 0.0D0)
                endif
            endif

            mutab%channel_w(i, j) = mutab%channel_eff(i)*SMrate
            !  mutab%channel_w_corrected_eff(i,j)=mutab%channel_eff_ratios(i)*mutab%channel_eff(i)*SMrate
        end do

        if(trim(mutab%label).eq."CMS-PAS-HIG-18-019") then
        !     print*,"mutab%channel_eff(i)",mutab%channel_eff
        ! print*,"mutab%channel_w",mutab%channel_w
        endif

        ! write(*,*) 'DEBUG HS: BRs = ', t%BR_hjWW, t%BR_hjZZ, t%BR_hjgaga
        ! write(*,*) 'DEBUG HS: LHC13 = ', t%lhc13%XS_hj_ratio, t%lhc13%XS_vbf_ratio, t%lhc13%XS_hjW_ratio,&
        !                                 t%lhc13%XS_hjZ_ratio, t%lhc13%XS_tthj_ratio

        SMrate = sum(mutab%channel_w(:, j))
        ! write(*,*) 'DEBUG HS: SMrate = ', SMrate
        ! modelrate=sum(mutab%channel_w_corrected_eff(:,j))

        do i = 1, mutab%Nc
            mutab%channel_w(i, j) = div(mutab%channel_w(i, j), SMrate, 0.0D0, 1.0D9)
            !  mutab%channel_w_corrected_eff(i,j)=div(mutab%channel_w_corrected_eff(i,j),modelrate,0.0D0,1.0D9)
        end do

        ! (TS 30/10/2013):
        ! write(*,*) "get_channelrates (mu, w, weff), id = ", mutab%id
        ! write(*,*) mutab%channel_mu
        ! write(*,*)  mutab%channel_w
        ! write(*,*) mutab%channel_eff_ratios
        do i = 1, mutab%Nc
            mutab%channel_w_corrected_eff(i, j) = mutab%channel_eff_ratios(i)*mutab%channel_w(i, j)
            ! n.b.: model weights are not normalized to 1!
        end do

        ! write(*,*) j,mutab%id, "SM         = ", mutab%channel_w(:,j)
        ! write(*,*) j,mutab%id, "SM effcorr = ",mutab%channel_w_corrected_eff(:,j)

        do i = 1, mutab%Nc
            drsq_SM = 0.0D0
            drsq = 0.0D0

            !   id1 = mutab%channel_id(i)
            !   p1 = int((id1-modulo(id1,10))/dble(10))
            !   d1 = modulo(id1,10)
            p1 = mutab%channel_p_id(i)
            d1 = mutab%channel_d_id(i)
            if (mutab%collider .ne. 'ILC') then
                do ii = 1, mutab%Nc
                    p2 = mutab%channel_p_id(ii)
                    d2 = mutab%channel_d_id(ii)
                    !     id2 = mutab%channel_id(ii)
                    !     p2 = int((id2-modulo(id2,10))/dble(10))
                    !     d2 = modulo(id2,10)
                    if (p1 .eq. p2 .and. p1 .ne. 0) then
                        if (delta_rate%CScov_ok .and. delta_rate%usecov) then
                            !-- TS 29/03/2017: Add 13 TeV XS covariance matrix here
                            if (abs(mutab%energy - 13.0D0) .le. small) then
                                drcov = delta_rate%CS13cov(p1, p1)
                                drcovSM = delta_rate%CS13covSM(p1, p1)
                            else
                                drcov = delta_rate%CScov(p1, p1)
                                drcovSM = delta_rate%CScovSM(p1, p1)
                            end if
                            drsq = drsq + drcov*mutab%channel_w_corrected_eff(i, j)*mutab%channel_w_corrected_eff(ii, j)
                            drsq_SM = drsq_SM + drcovSM*mutab%channel_w(i, j)*mutab%channel_w(ii, j)
                        else
                            drsq = drsq + delta_rate%dCS(p1)**2 &
                                   *mutab%channel_w_corrected_eff(i, j)*mutab%channel_w_corrected_eff(ii, j)
                            drsq_SM = drsq_SM + delta_rate%dCS_SM(p1)**2*mutab%channel_w(i, j)*mutab%channel_w(ii, j)
                        end if
                    end if
                    if (d1 .eq. d2 .and. d1 .ne. 0) then
                        if (delta_rate%BRcov_ok .and. delta_rate%usecov) then
                            dBRSM = delta_rate%BRcovSM(d1, d1)
                            dBR = delta_rate%BRcov(d1, d1)
                        else
                            dBRSM = delta_rate%dBR_SM(d1)**2
                            dBR = delta_rate%dBR(d1)**2
                        end if
                        drsq = drsq + dBR*mutab%channel_w_corrected_eff(i, j)*mutab%channel_w_corrected_eff(ii, j)
                        drsq_SM = drsq_SM + dBRSM*mutab%channel_w(i, j)*mutab%channel_w(ii, j)
                    end if
                end do
            end if
            mutab%channel_syst(i, j) = sqrt(drsq)
            mutab%channel_systSM(i, j) = sqrt(drsq_SM)
        end do

        ! write(*,*) 'DEBUG HS: mu = ', mutab%channel_mu
        ! write(*,*) 'DEBUG HS: w = ', mutab%channel_w
        ! write(*,*) 'DEBUG HS: eff = ', mutab%channel_eff

    end subroutine get_channelrates

    subroutine evaluate_LHC_Run1_combination(t, n, signalScaling)
        use usefulbits, only: np, Hneut, dataset, results, vsmall
        use usefulbits_hs, only: HSresults, Nparam, pdf, &
                                 LHCrun1_rates, LHCrun1_correlationmatrix, &
                                 assignmentrange, assignmentrange_massobs, HSres, normalize_rates_to_reference_position, &
                                 normalize_rates_to_reference_position_outside_dmtheo
        use pc_chisq, only: csq_mh
        use numerics, only: invmatrix, matmult, gammp

        implicit none
        !--------------------------------------input
        type(dataset), intent(in) :: t
        integer, intent(in) :: n
        double precision, intent(in), optional :: signalScaling
        !--------------------------------------internal
        integer :: p, d, id, i, j, k, ncomb
        double precision, allocatable :: covmat(:, :), invcovmat(:, :)
        double precision, allocatable :: covmatzero(:, :), invcovmatzero(:, :)
        integer, allocatable :: index_assigned_Higgs(:)
        double precision, allocatable :: relative_signal_strength(:, :)
        double precision, dimension(20) :: v, v2, csq_mu, vzero, vzero2, csq_mu_max
        double precision, dimension(20, 1) :: vmat, vzeromat
        double precision :: mobs = 125.09D0
        double precision :: dmobs = 0.24D0
        double precision :: dmbbtautau = 20.0D0
        double precision :: dmWW = 5.0D0
        double precision :: expmassrange, allowed_massrange
        double precision :: Higgs_signal_k
        double precision :: num1, num2, dnum1, dnum2, denom1, denom2, mav, dmav, csq_sep
        double precision :: globalScaling, SM_rate

        allocate (covmat(20, 20), invcovmat(20, 20))
        allocate (covmatzero(20, 20), invcovmatzero(20, 20))
        allocate (index_assigned_Higgs(np(Hneut)))
        allocate (relative_signal_strength(2, np(Hneut)))

        if (present(signalScaling)) then
            globalScaling = signalScaling
        else
            globalScaling = 1D0
        endif
        mav = 0.0D0
        dmav = 0.0D0

        denom1 = 0.0D0
        denom2 = 0.0D0
        num1 = 0.0D0
        num2 = 0.0D0
        dnum1 = 0.0D0
        dnum2 = 0.0D0

        ! initialize all entries with zero; set k-th entry to 1 if Higgs h_k is assigned
        index_assigned_Higgs = 0
        relative_signal_strength = 0.0D0 !(first index --> gaga or ZZ, second index -> Higgs)

        do i = lbound(LHCrun1_rates, dim=1), ubound(LHCrun1_rates, dim=1)
            id = LHCrun1_rates(i)%channel_id
            p = int((id - modulo(id, 10))/dble(10))
            d = modulo(id, 10)

            if (d .eq. 4 .or. d .eq. 5) then
                expmassrange = dmbbtautau
            elseif (d .eq. 2) then
                expmassrange = dmWW
            else
                ! if (pdf .eq. 1) then
                expmassrange = dmobs
                ! else
                ! expmassrange = assignmentrange_LHCrun1*dmobs
                ! endif
            end if

            LHCrun1_rates(i)%r_pred = 0.0D0

            ncomb = 0
            do k = 1, np(Hneut)
                if (d .eq. 1 .or. d .eq. 3) then ! Apply assignment range for Gaussian pdf for gaga and ZZ final states.
                    if ((pdf .eq. 1) .or. (pdf .eq. 3)) then ! box and box+Gaussian cases
                        allowed_massrange = assignmentrange_massobs*expmassrange + t%particle(Hneut)%dM(k)
                    else
                        allowed_massrange = assignmentrange_massobs*sqrt(expmassrange**2.0D0 + t%particle(Hneut)%dM(k)**2.0D0)
                    endif
                else
                    if ((pdf .eq. 1) .or. (pdf .eq. 3)) then ! box and box+Gaussian cases
                        allowed_massrange = assignmentrange*expmassrange + t%particle(Hneut)%dM(k)
                    else
                        allowed_massrange = assignmentrange*sqrt(expmassrange**2.0D0 + t%particle(Hneut)%dM(k)**2.0D0)
                    end if
                    ! allowed_massrange = sqrt(expmassrange**2.0D0 + t%particle(Hneut)%dM(k)**2.0D0)
                end if

                if (abs(t%particle(Hneut)%M(k) - mobs) .le. allowed_massrange) then
                    index_assigned_Higgs(k) = 1
                    Higgs_signal_k = signalrate(k, p, d, mobs, t%particle(Hneut)%M(k), t%particle(Hneut)%dM(k))
                    LHCrun1_rates(i)%r_pred = LHCrun1_rates(i)%r_pred + Higgs_signal_k
                    if (id .eq. 11) then ! gg -> h_k -> gaga weighted mass average
                        num1 = num1 + Higgs_signal_k*t%particle(Hneut)%M(k)
                        dnum1 = dnum1 + Higgs_signal_k*t%particle(Hneut)%dM(k)
                        relative_signal_strength(1, k) = Higgs_signal_k
                    else if (id .eq. 13) then ! gg -> h_k -> ZZ -> 4l weighted mass average
                        num2 = num2 + Higgs_signal_k*t%particle(Hneut)%M(k)
                        dnum2 = dnum2 + Higgs_signal_k*t%particle(Hneut)%dM(k)
                        relative_signal_strength(2, k) = Higgs_signal_k
                    end if
                    ncomb = ncomb + 1
                end if
            end do

            SM_rate = SM_production_rate(mobs,p)*SM_decay_rate(mobs,d)
            LHCrun1_rates(i)%r_pred = SM_rate + globalScaling*(LHCrun1_rates(i)%r_pred - SM_rate)

            if (id .eq. 11) then
                denom1 = LHCrun1_rates(i)%r_pred
                relative_signal_strength(1, :) = relative_signal_strength(1, :)/denom1
            else if (id .eq. 13) then
                denom2 = LHCrun1_rates(i)%r_pred
                relative_signal_strength(2, :) = relative_signal_strength(2, :)/denom2
            end if

            if (LHCrun1_rates(i)%r_pred .gt. LHCrun1_rates(i)%r) then
                LHCrun1_rates(i)%dr = LHCrun1_rates(i)%dr_up
            else
                LHCrun1_rates(i)%dr = LHCrun1_rates(i)%dr_low
            end if

            if (LHCrun1_rates(i)%r .lt. 0.0D0) then
                LHCrun1_rates(i)%dr0 = LHCrun1_rates(i)%dr_up
            else
                LHCrun1_rates(i)%dr0 = LHCrun1_rates(i)%dr_low
            end if

            v(i) = LHCrun1_rates(i)%r_pred - LHCrun1_rates(i)%r
            vmat(i, 1) = v(i)

            vzero(i) = LHCrun1_rates(i)%r
            vzeromat(i, 1) = vzero(i)
        end do

        if (denom1 .gt. vsmall .and. denom2 .gt. vsmall) then
            mav = 0.5D0*(num1/denom1 + num2/denom2)
            dmav = 0.5D0*(dnum1/denom1 + dnum2/denom2)
            !   write(*,*) "Averaged mass is ",mav, " +- ",dmav
            !  else
            !    write(*,*) "denom1 and denom2 are ",denom1, denom2
        end if

        do i = lbound(LHCrun1_rates, dim=1), ubound(LHCrun1_rates, dim=1)
            do j = lbound(LHCrun1_rates, dim=1), ubound(LHCrun1_rates, dim=1)
                covmat(i, j) = LHCrun1_correlationmatrix(i, j)* &
                               LHCrun1_rates(i)%dr*LHCrun1_rates(j)%dr
                covmatzero(i, j) = LHCrun1_correlationmatrix(i, j)* &
                                   LHCrun1_rates(i)%dr0*LHCrun1_rates(j)%dr0
            end do
        end do

        call invmatrix(covmat, invcovmat)
        call matmult(invcovmat, vmat, v2, 20, 1)

        call invmatrix(covmatzero, invcovmatzero)
        call matmult(invcovmatzero, vzeromat, vzero2, 20, 1)

        do i = 1, 20
            csq_mu(i) = v(i)*v2(i)
        end do

        do i = 1, 20
            csq_mu_max(i) = vzero(i)*vzero2(i)
        end do

        if (mav .lt. vsmall) then
            HSres(n)%Chisq_LHCRun1_mh = 0.0D0
            HSres(n)%Chisq_LHCRun1_mhsep = 0.0D0
        else
            csq_sep = 0.0D0
            do k = 1, np(Hneut)
                csq_sep = csq_sep + csq_sep_mh_per_Higgs(t%particle(Hneut)%M(k), mav, &
                                                         t%particle(Hneut)%dM(k), dmobs, relative_signal_strength(:, k))
            end do
            HSres(n)%Chisq_LHCRun1_mhsep = csq_sep
            HSres(n)%Chisq_LHCRun1_mh = csq_mh(mav, mobs, dmav, dmobs) + csq_sep
        end if

        if ((HSres(n)%Chisq_LHCRun1_mh + sum(csq_mu)) .gt. sum(csq_mu_max).and..not.present(signalScaling)) then
            HSres(n)%Chisq_LHCRun1_mu = sum(csq_mu_max)
            HSres(n)%Chisq_LHCRun1_mh = 0.0D0
        else
            HSres(n)%Chisq_LHCRun1_mu = sum(csq_mu)
        end if

        HSres(n)%Chisq_LHCRun1 = HSres(n)%Chisq_LHCRun1_mu + HSres(n)%Chisq_LHCRun1_mh
        HSres(n)%nobs_LHCRun1_mu = 20
        HSres(n)%nobs_LHCRun1_mh = 1
        if (HSres(n)%Chisq_LHCRun1 .gt. vsmall .and. (HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_LHCRun1_mh - Nparam) .gt. 0) then
            HSres(n)%Pvalue_LHCRun1 = 1 - gammp(dble(HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_LHCRun1_mh - Nparam)/2, &
                                                HSres(n)%Chisq_LHCRun1/2)
        end if

        deallocate (covmat, invcovmat)
        deallocate (covmatzero, invcovmatzero)

    contains

        function SM_production_rate(mass,p)
            use usefulbits_hs, only: LHC_combination_run1_SMXS_from_paper
            !--------------------------------------external functions
            double precision, external :: SMCS_lhc8_gg_H, SMCS_lhc8_bb_H, SMCS_lhc8_vbf_H, &
                SMCS_lhc8_HW, SMCS_lhc8_HZ, SMCS_lhc8_ttH, SMCS_lhc8_tWH, SMCS_lhc8_tH_schan, SMCS_lhc8_tH_tchan
            integer, intent(in) :: p
            double precision, intent(in) :: mass
            double precision :: SM_production_rate

            select case (p)
            case (1)
                if (LHC_combination_run1_SMXS_from_paper) then
                    SM_production_rate = 19.2D0 + 0.203D0
                else
                    SM_production_rate = SMCS_lhc8_gg_H(mass) + SMCS_lhc8_bb_H(mass)
                endif
            case (2)
                if (LHC_combination_run1_SMXS_from_paper) then
                    SM_production_rate = 1.58D0
                else
                    SM_production_rate = SMCS_lhc8_vbf_H(mass)
                endif
            case (3)
                if (LHC_combination_run1_SMXS_from_paper) then
                    SM_production_rate = 0.703D0
                else
                    SM_production_rate = SMCS_lhc8_HW(mass)
                endif
            case (4)
                if (LHC_combination_run1_SMXS_from_paper) then
                    SM_production_rate = 0.446D0
                else
                    SM_production_rate = SMCS_lhc8_HZ(mass)
                endif
            case (5)
                if (LHC_combination_run1_SMXS_from_paper) then
                    SM_production_rate = 0.129D0
                else
                    SM_production_rate = SMCS_lhc8_ttH(mass) + &
                                         SMCS_lhc8_tH_schan(mass) + &
                                         SMCS_lhc8_tH_tchan(mass) + &
                                         SMCS_lhc8_tWH(mass)

                endif
            case DEFAULT
                write (*, *) "Invalid production channel ", p, " in evaluate_LHC_run1_combination."
                stop 1
            end select
        end function SM_production_rate

        function SM_decay_rate(mass, d)
            double precision :: SMBR_Hgamgam, SMBR_HWW, SMBR_HZZ, SMBR_Htautau, SMBR_Hbb
            double precision, intent(in) :: mass
            integer, intent(in) :: d
            double precision :: SM_decay_rate

            select case (d)
            case (1)
                SM_decay_rate = SMBR_Hgamgam(mass)
            case (2)
                SM_decay_rate = SMBR_HWW(mass)
            case (3)
                SM_decay_rate = SMBR_HZZ(mass)
            case (4)
                SM_decay_rate = SMBR_Htautau(mass)
            case (5)
                SM_decay_rate = SMBR_Hbb(mass)
            case DEFAULT
                write (*, *) "Invalid decay channel ", d, " in evaluate_LHC_run1_combination."
                stop 1
            end select
        end function SM_decay_rate
        !------------------------------------------------------------
        function signalrate(k, p, d, mobs, m, dmtheo)
            !------------------------------------------------------------
            use usefulbits_hs, only: LHC_combination_run1_SMXS_from_paper
            !--------------------------------------external functions
            double precision :: SMCS_lhc8_gg_H, SMCS_lhc8_bb_H, SMCS_lhc8_vbf_H, &
                SMCS_lhc8_HW, SMCS_lhc8_HZ, SMCS_lhc8_ttH, SMBR_Hgamgam, SMBR_HWW, &
                SMBR_HZZ, SMBR_Htautau, SMBR_Hbb, SMCS_lhc8_tH_schan, SMCS_lhc8_tH_tchan, &
                SMCS_lhc8_tWH
            double precision, intent(in) :: mobs, m, dmtheo
            integer, intent(in) :: k, p, d
            double precision :: signalrate, production_rate, decay_rate, mass, refmass
            double precision :: production_rate_scalefactor, decay_rate_scalefactor
            mass = t%particle(Hneut)%M(k)

            ! TS (17/10/2018): Take reference mass for SM-normalization at mobs+dmtheo box boundary.
            if (mass .ge. (mobs + dmtheo)) then
                refmass = mobs + dmtheo
            else if (mass .le. (mobs - dmtheo)) then
                refmass = mobs - dmtheo
            else
                refmass = mass
            end if
            !---

            select case (p)
            case (1)
                if (LHC_combination_run1_SMXS_from_paper) then
                    production_rate = t%lhc8%XS_gg_hj_ratio(k)*19.2D0 &
                                      + t%lhc8%XS_bb_hj_ratio(k)*0.203D0
                else
                    production_rate = t%lhc8%XS_gg_hj_ratio(k)*SMCS_lhc8_gg_H(mass) &
                                      + t%lhc8%XS_bb_hj_ratio(k)*SMCS_lhc8_bb_H(mass)
                end if
                ! NOTE: Here we make a small error in the scalefactor. Correct would be to rescale
                !       the gg and bb contributions separately.
                production_rate_scalefactor = (SMCS_lhc8_gg_H(mobs) + SMCS_lhc8_bb_H(mobs))/ &
                                              (SMCS_lhc8_gg_H(refmass) + SMCS_lhc8_bb_H(refmass))
            case (2)
                production_rate = t%lhc8%XS_vbf_ratio(k)*SM_production_rate(mass, 2)
                production_rate_scalefactor = SMCS_lhc8_vbf_H(mobs)/SMCS_lhc8_vbf_H(refmass)
            case (3)
                production_rate = t%lhc8%XS_hjW_ratio(k)*SM_production_rate(mass,3)
                production_rate_scalefactor = SMCS_lhc8_HW(mobs)/SMCS_lhc8_HW(refmass)
            case (4)
                production_rate = t%lhc8%XS_hjZ_ratio(k)*SM_production_rate(mass,4)
                production_rate_scalefactor = SMCS_lhc8_HZ(mobs)/SMCS_lhc8_HZ(refmass)
            case (5)
                if (LHC_combination_run1_SMXS_from_paper) then
                    production_rate = t%lhc8%XS_tthj_ratio(k)*0.129D0 ! TODO: Still need to include the other channels!
                else
                    production_rate = t%lhc8%XS_tthj_ratio(k)*SMCS_lhc8_ttH(mass) + &
                                      t%lhc8%XS_thj_schan_ratio(k)*SMCS_lhc8_tH_schan(mass) + &
                                      t%lhc8%XS_thj_tchan_ratio(k)*SMCS_lhc8_tH_tchan(mass) + &
                                      t%lhc8%XS_tWhj_ratio(k)*SMCS_lhc8_tWH(mass)
                end if
                production_rate_scalefactor = SMCS_lhc8_ttH(mobs)/SMCS_lhc8_ttH(refmass)
            case DEFAULT
                write (*, *) "Invalid production channel ", p, " in evaluate_LHC_run1_combination."
                stop 1
            end select
            select case (d)
            case (1)
                decay_rate = t%BR_hjgaga(k)
                decay_rate_scalefactor = SMBR_Hgamgam(mobs)/SMBR_Hgamgam(refmass)
            case (2)
                decay_rate = t%BR_hjWW(k)
                decay_rate_scalefactor = SMBR_HWW(mobs)/SMBR_HWW(refmass)
            case (3)
                decay_rate = t%BR_hjZZ(k)
                decay_rate_scalefactor = SMBR_HZZ(mobs)/SMBR_HZZ(refmass)
            case (4)
                decay_rate = t%BR_hjtautau(k)
                decay_rate_scalefactor = SMBR_Htautau(mobs)/SMBR_Htautau(refmass)
            case (5)
                decay_rate = t%BR_hjbb(k)
                decay_rate_scalefactor = SMBR_Hbb(mobs)/SMBR_Hbb(refmass)
            case DEFAULT
                write (*, *) "Invalid decay channel ", d, " in evaluate_LHC_run1_combination."
                stop 1
            end select

            if (normalize_rates_to_reference_position) then
                signalrate = production_rate*decay_rate
            else
                ! This is the default:
                signalrate = production_rate*production_rate_scalefactor* &
                             decay_rate*decay_rate_scalefactor
            end if

            if (normalize_rates_to_reference_position_outside_dmtheo) then
                if (abs(mobs - m) .ge. dmtheo) then
                    signalrate = production_rate*decay_rate
                end if
            end if

        end function signalrate

        function csq_sep_mh_per_Higgs(mh, mav, dmh, dm_exp, relative_mu)
            use usefulbits_hs, only: significantcontribution

            double precision, intent(in) :: mh, mav, dmh, dm_exp
            double precision, intent(in) :: relative_mu(2)
            double precision :: csq_sep_mh_per_Higgs, relative_mu_average

            csq_sep_mh_per_Higgs = 0.0D0

            relative_mu_average = 0.5D0*(relative_mu(1) + relative_mu(2))

            select case (pdf)
            case (1) ! box pdf
                if (relative_mu_average .ge. significantcontribution) then
                    if (abs(mh - mav) .gt. (dmh + dm_exp)) then
                        csq_sep_mh_per_Higgs = 100000.0D0  !Large number
                    else
                        continue
                    end if
                end if
            case (2) ! Gaussian pdf
                csq_sep_mh_per_Higgs = relative_mu_average*(mh - mav)**2.0D0/ &
                                       (dmh**2.0D0 + dm_exp**2.0D0)
            case (3)
                if ((mh + dmh) .le. mav) then
                    csq_sep_mh_per_Higgs = (mh + dmh - mav)**2.0D0/dm_exp**2.0D0
                else if ((mh - dmh) .gt. mav) then
                    csq_sep_mh_per_Higgs = (mh - dmh - mav)**2.0D0/dm_exp**2.0D0
                else
                    csq_sep_mh_per_Higgs = 0.0D0
                end if
                csq_sep_mh_per_Higgs = relative_mu_average*csq_sep_mh_per_Higgs
            end select

        end function csq_sep_mh_per_Higgs

    end subroutine evaluate_LHC_Run1_combination

end module evaluate
