!--------------------------------------------------------------------------------------
! This example program is part of HiggsBounds-5 (TS 22/02/2017).
!--------------------------------------------------------------------------------------
program HBeffC
! Example program to demonstrate the effective coupling input for HiggsBounds.
! The test point is a complete toy model with 3 neutral Higgs bosons and 1 charged
! Higgs boson, it does not correspond to any sensible new physics model.
!--------------------------------------------------------------------------------------
    use usefulbits, only: analysislist
    implicit none

    integer :: nHzero, nHplus
    integer :: HBresult, chan, ncombined, chan2
    integer :: HBresult_tmp, chan_tmp, ncombined_tmp
    integer, parameter :: fileid = 78
    double precision :: obsratio, obsratio_tmp, predratio_tmp
    double precision :: SMGamma_h, GammaTotal(3)
    double precision :: Mh(3), ghjss_s(3), ghjss_p(3), ghjcc_s(3), ghjcc_p(3), &
        ghjbb_s(3), ghjbb_p(3), ghjtt_s(3), ghjtt_p(3), &
        ghjmumu_s(3), ghjmumu_p(3), ghjtautau_s(3), ghjtautau_p(3), &
        ghjWW(3), ghjZZ(3), ghjZga(3), ghjgaga(3), ghjgg(3), &
        ghjhiZ(3), BR_hkhjhi(3, 3, 3), BR_hjhiZ(3, 3), BR_hjinvisible(3), &
        BR_hjemu(3), BR_hjetau(3), BR_hjmutau(3)
    double precision :: Mhplus, GammaTotal_Hpj, &
        CS_lep_HpjHmj_ratio, &
        BR_tWpb, BR_tHpjb, &
        BR_Hpjcs, BR_Hpjcb, BR_Hpjtaunu, &
        BR_Hpjtb, BR_HpjWZ, BR_HpjhiW(1, 3), &
        BR_hjHpiW(3, 1)
    double precision :: dMhneut(3), dMhch(1)
    double precision :: CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, CS_Hpjcjet, &
        CS_Hpjjetjet, CS_HpjW, CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi(1, 3)
    double precision :: nonSMBR(3), sumiBR_hkhihj(3), sumiBR_hihjZ(3)
    character(len=100)::filename
    integer :: i, j, k
    double precision :: singleH_tmp, ggH_tmp, bbH_tmp, VBF_tmp, WH_tmp, &
        ZH_tmp, ttH_tmp, tH_tchan_tmp, tH_schan_tmp, &
        BR_hjss_tmp, BR_hjcc_tmp, BR_hjbb_tmp, &
        BR_hjtt_tmp, BR_hjmumu_tmp, BR_hjtautau_tmp, BR_hjWW_tmp, &
        BR_hjZZ_tmp, BR_hjZga_tmp, BR_hjgaga_tmp, BR_hjgg_tmp, &
        qqZh_tmp, ggZH_tmp

!! We could use a subset of experimental analysis by employing the analysislist:
!
!! Allocate list of analyses ID numbers to include
!! and initialize with bogus (-1) numbers:
!
!   allocate(analysislist(15))
!   do i=1,ubound(analysislist,1)
!    analysislist(i) = -1
!   enddo
!
!   analysislist(5) = 14011
!   analysislist(2) = 160312
!   analysislist(3) = 20160851
!   analysislist(4) = 20160852
!   analysislist(5) = 20160151
!   analysislist(6) = 20160152
!   analysislist(7) = 16002
!   analysislist(8) = 16029
!   analysislist(9) = 14013
!   analysislist(10) = 2016071
!   analysislist(11) = 20160741
!   analysislist(12) = 20160742
!   analysislist(13) = 2016062
!   analysislist(14) = 2016049
!   analysislist(3) = 044781
!   analysislist(4) = 044782
!   analysislist(1) = 011811
!   analysislist(2) = 011812
!   analysislist(15) = 6065
!   analysislist(3) = 2016059
!   analysislist(4) = 14022

    nHzero = 3
    nHplus = 1

!!  call initialize_HiggsBounds(nHzero, nHplus, "list ")
    call initialize_HiggsBounds(nHzero, nHplus, "LandH")

!! set the neutral Higgs masses and uncertainties
    Mh = (/25.0D0, 125.0D0, 340.0D0/)
    dMhneut = (/0.0D0, 0.0D0, 0.0D0/)

!! set some example non-SM neutral Higgs decays
    BR_hkhjhi = 0.0D0
    BR_hkhjhi(3, 1, 2) = 0.01D0
    BR_hkhjhi(3, 2, 1) = 0.01D0
    BR_hkhjhi(3, 1, 1) = 0.1D0
    BR_hkhjhi(2, 1, 1) = 0.4D0
    BR_hjhiZ = 0.0D0
    BR_hjhiZ(3, 2) = 0.1D0
    BR_hjhiZ(3, 1) = 0.1D0
    BR_hjhiZ(2, 1) = 0.1D0
    BR_hjHpiW(1, 1) = 0.0D0
    BR_hjHpiW(2, 1) = 0.1D0
    BR_hjHpiW(3, 1) = 0.0D0
    BR_hjemu = 0.0D0
    BR_hjetau = 0.0D0
    BR_hjmutau = 0.0D0
    BR_hjinvisible = 0d0

    do i = 1, 3
        sumiBR_hkhihj(i) = 0.0D0
        sumiBR_hihjZ(i) = 0.0D0
        do j = 1, 3
            do k = 1, j
                sumiBR_hkhihj(i) = sumiBR_hkhihj(i) + BR_hkhjhi(i, j, k)
            enddo
            sumiBR_hihjZ(i) = sumiBR_hihjZ(i) + BR_hjhiZ(i, j)
        enddo

        nonSMBR(i) = BR_hjemu(i) + BR_hjetau(i) + BR_hjmutau(i) + &
                     BR_hjinvisible(i) + sumiBR_hihjZ(i) + sumiBR_hkhihj(i) + &
                     BR_hjHpiW(i, 1)
    enddo

!!  Set some example values for the neutral Higgs effective couplings. Here,
!!  as a toy example, set all Higgs couplings to the SM Higgs couplings.
    do i = 1, 3

        ghjss_s(i) = 1.0d0
        ghjss_p(i) = 0.0d0
        ghjcc_s(i) = 1.0d0
        ghjcc_p(i) = 0.0d0
        ghjbb_s(i) = 1.0d0
        ghjbb_p(i) = 0.0d0
        ghjtt_s(i) = 1.0d0
        ghjtt_p(i) = 0.0d0
        ghjmumu_s(i) = 1.0d0
        ghjmumu_p(i) = 0.0d0
        ghjtautau_s(i) = 1.0d0
        ghjtautau_p(i) = 0.0d0
        ghjWW(i) = 1.0d0
        ghjZZ(i) = 1.0d0
        ghjZga(i) = 1d0
        ghjgaga(i) = 1.0d0
        ghjgg(i) = 1.0d0
        ghjhiZ(i) = 0d0

!         GammaTotal(i) = SMGamma_h(Mh(i))/(1 - nonSMBR(i))
        GammaTotal(i) = -1.0D0 ! Derive total decay width internally from provided input!
    enddo

    ghjWW(1) = 0.0d0
    ghjZZ(1) = 0.0d0

!! Charged Higgs input (example values)
    Mhplus = 55.0D0
    dMhch = (/0.0D0/)
    GammaTotal_Hpj = 1.0D0 ! n.b.: this is not needed at the moment
    CS_lep_HpjHmj_ratio = 1.0D0
    BR_tWpb = 1.0D0
    BR_tHpjb = 0.0D0
    BR_Hpjcs = 0.0D0
    BR_Hpjcb = 0.0D0
    BR_Hpjtaunu = 0.2D0
    BR_Hpjtb = 0.00D0
    BR_HpjWZ = 0.00D0
    BR_HpjhiW(1, 1) = 0.8D0
    BR_HpjhiW(1, 2) = 0.0D0
    BR_HpjhiW(1, 3) = 0.0D0
!! Cross sections (CS) are in pb:
    CS_Hpjtb = 1.0D-02
!! n.b.: for the following quantities there are currently no experimental searches:
    CS_Hpjcb = 0.5D-02
    CS_Hpjbjet = 1.0D-02
    CS_Hpjcjet = 1.0D-02
    CS_Hpjjetjet = 1.0D-03
    CS_HpjW = 0.01D0
    CS_HpjZ = 0.01D0
    CS_vbf_Hpj = 0.01D0
    CS_HpjHmj = 0.01D0
    CS_Hpjhi(1, 1) = 0.01D0
    CS_Hpjhi(1, 2) = 0.005D0
    CS_Hpjhi(1, 3) = 0.001D0

!! Calling the HiggsBounds input routines:

    call HiggsBounds_neutral_input_properties(Mh, GammaTotal, (/0, 0, 0/))

    call HiggsBounds_set_mass_uncertainties(dMhneut, dMhch)

    call HiggsBounds_neutral_input_effC( &
        ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
        ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
        ghjmumu_s, ghjmumu_p, &
        ghjtautau_s, ghjtautau_p, &
        ghjWW, ghjZZ, ghjZga, &
        ghjgaga, ghjgg, ghjhiZ)

    call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi, &
                                           BR_hjhiZ, BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW)

    call HiggsBounds_charged_input(Mhplus, GammaTotal_Hpj, &
                                   CS_lep_HpjHmj_ratio, &
                                   BR_tWpb, BR_tHpjb, &
                                   BR_Hpjcs, BR_Hpjcb, BR_Hpjtaunu, BR_Hpjtb, &
                                   BR_HpjWZ, BR_HpjhiW)

    call HiggsBounds_charged_input_hadr(13, CS_Hpjtb, CS_Hpjcb, &
                                        CS_Hpjbjet, CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
                                        CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)

!! Run HiggsBounds and get some more details (most sensitive channels) from the output:

    open (fileid, file="HBeffC-output.dat")

    call run_HiggsBounds(HBresult, chan, obsratio, ncombined)
    write (*, *) '**********************************************************************'
    write (*, *) 'HiggsBounds main results'
    write (*, *) '**********************************************************************'
    write (*, '(A,1I5,A)') ' HBresult   = ', HBresult, '          (0: excluded, 1: allowed)'
    write (*, '(A,1I5)') ' channel ID = ', chan
    write (*, '(A,1F8.3)') ' obsratio   = ', obsratio
    write (*, '(A,1I5)') ' ncombined  = ', ncombined
    write (*, *) '**********************************************************************'
    write (*, *) 'Control output for some (SM normalized) cross sections and BRs '
    write (*, *) '**********************************************************************'

    write (fileid, '(A,1I5,A)') ' HBresult   = ', HBresult, '          (0: excluded, 1: allowed)'
    write (fileid, '(A,1I5)') ' channel ID = ', chan
    write (fileid, '(A,1F8.3)') ' obsratio   = ', obsratio
    write (fileid, '(A,1I5)') ' ncombined  = ', ncombined

    do i = 1, 3
        call HiggsBounds_get_neutral_hadr_CS(i, 13, singleH_tmp, ggH_tmp, &
                                             bbH_tmp, VBF_tmp, WH_tmp, ZH_tmp, ttH_tmp, tH_tchan_tmp, tH_schan_tmp, &
                                             qqZh_tmp, ggZH_tmp)

        call HiggsBounds_get_neutral_BR(i, BR_hjss_tmp, BR_hjcc_tmp, BR_hjbb_tmp, &
                                        BR_hjtt_tmp, BR_hjmumu_tmp, BR_hjtautau_tmp, BR_hjWW_tmp, &
                                        BR_hjZZ_tmp, BR_hjZga_tmp, BR_hjgaga_tmp, BR_hjgg_tmp)

        write (*, "(A,I1,A,1F5.1,A)") ' Higgs boson h', i, ' (m = ', Mh(i), ' GeV) at 13 TeV LHC:'
        write (*, *) ' '
        write (*, "(A,I1,A,1F8.3)") ' pp -> h', i, ' (inclusive) :', singleH_tmp
        write (*, "(A,I1,A,1F8.3)") ' gg -> h', i, ':             ', ggH_tmp
        write (*, "(A,I1,A,1F8.3)") ' bb -> h', i, ':             ', bbH_tmp
        write (*, "(A,I1,A,1F8.3)") ' pp -> qqh', i, ':           ', VBF_tmp
        write (*, "(A,I1,A,1F8.3)") ' pp -> Wh', i, ':            ', WH_tmp
        write (*, "(A,I1,A,1F8.3)") ' pp -> Zh', i, ':            ', ZH_tmp
        write (*, "(A,I1,A,1F8.3)") ' pp -> tth', i, ':           ', ttH_tmp
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->ss):           ', 100.0D0*BR_hjss_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->cc):           ', 100.0D0*BR_hjcc_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->bb):           ', 100.0D0*BR_hjbb_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->tt):           ', 100.0D0*BR_hjtt_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->mumu):         ', 100.0D0*BR_hjmumu_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->tautau):       ', 100.0D0*BR_hjtautau_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->WW):           ', 100.0D0*BR_hjWW_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->ZZ):           ', 100.0D0*BR_hjZZ_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->Zga):          ', 100.0D0*BR_hjZga_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->gaga):         ', 100.0D0*BR_hjgaga_tmp, ' %'
        write (*, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->gg):           ', 100.0D0*BR_hjgg_tmp, ' %'
        write (*, *) ' '

        write (fileid, "(A,I1,A,1F5.1,A)") ' Higgs boson h', i, ' (m = ', Mh(i), ' GeV) at 13 TeV LHC:'
        write (fileid, *) ' '
        write (fileid, "(A,I1,A,1F8.3)") ' pp -> h', i, ' (inclusive) :', singleH_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' gg -> h', i, ':             ', ggH_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' bb -> h', i, ':             ', bbH_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' pp -> qqh', i, ':           ', VBF_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' pp -> Wh', i, ':            ', WH_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' pp -> Zh', i, ':            ', ZH_tmp
        write (fileid, "(A,I1,A,1F8.3)") ' pp -> tth', i, ':           ', ttH_tmp
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->ss):           ', 100.0D0*BR_hjss_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->cc):           ', 100.0D0*BR_hjcc_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->bb):           ', 100.0D0*BR_hjbb_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->tt):           ', 100.0D0*BR_hjtt_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->mumu):         ', 100.0D0*BR_hjmumu_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->tautau):       ', 100.0D0*BR_hjtautau_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->WW):           ', 100.0D0*BR_hjWW_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->ZZ):           ', 100.0D0*BR_hjZZ_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->Zga):          ', 100.0D0*BR_hjZga_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->gaga):         ', 100.0D0*BR_hjgaga_tmp, ' %'
        write (fileid, "(A,I1,A,1F8.3,A)") ' BR(h', i, '->gg):           ', 100.0D0*BR_hjgg_tmp, ' %'
        write (fileid, *) ' '
    enddo

    write (*, *) '**********************************************************************'
    write (*, *) 'Ranking of most sensitive channels for each Higgs boson'
    write (*, *) '**********************************************************************'
    write (*, *) ' Higgs-no      rank  HBresult   channel  obsratio  predratio ncombined'
    write (fileid, *) ' Higgs-no      rank  HBresult   channel  obsratio  predratio ncombined'
    do i = 1, 4
        do j = 1, 3
            call HiggsBounds_get_most_sensitive_channels_per_Higgs(i, j, HBresult_tmp, chan_tmp, &
                                                                   obsratio_tmp, predratio_tmp, ncombined_tmp)
            write (*, '(4I10,2F10.3,1I10)') i, j, HBresult_tmp, chan_tmp, obsratio_tmp, &
                predratio_tmp, ncombined_tmp
            write (fileid, '(4I10,2F10.3,1I10)') i, j, HBresult_tmp, chan_tmp, obsratio_tmp, &
                predratio_tmp, ncombined_tmp
        enddo
    enddo
    write (*, *) '**********************************************************************'
    write (*, *) 'Ranking of most sensitive channels (overall)'
    write (*, *) '**********************************************************************'
    write (*, *) '     rank  HBresult   channel  obsratio  predratio ncombined'
    write (fileid, *) '     rank  HBresult   channel  obsratio  predratio ncombined'
    do j = 1, 3
        call HiggsBounds_get_most_sensitive_channels(j, HBresult_tmp, chan_tmp, obsratio_tmp, &
                                                     predratio_tmp, ncombined_tmp)
        write (*, '(3I10,2F10.3,1I10)') j, HBresult_tmp, chan_tmp, obsratio_tmp, &
            predratio_tmp, ncombined_tmp
        write (fileid, '(3I10,2F10.3,1I10)') j, HBresult_tmp, chan_tmp, obsratio_tmp, &
            predratio_tmp, ncombined_tmp
    enddo
    write (*, *) '**********************************************************************'

    call finish_HiggsBounds

end program HBeffC
