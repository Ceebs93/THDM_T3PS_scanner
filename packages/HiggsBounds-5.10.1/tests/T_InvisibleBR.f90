!> Unit test for the cascade decay to invisible
program Test_InvisibleBR
    use usefulbits, only: theo, whichinput
    use theo_manip, only: HB5_complete_theo

    implicit none
    double precision :: SMGamma_h, GammaTotal(3)
    double precision :: Mh(3), BR_hkhjhi(3, 3, 3), BR_hjhiZ(3, 3), BR_hjinvisible(3), &
        BR_hjemu(3), BR_hjetau(3), BR_hjmutau(3), BR_HjHpiW(3, 0)
    double precision :: sumiBR_hkhihj(3), sumiBR_hihjZ(3), nonSMBr(3)
    integer :: i, j, k
    double precision :: testval

    call initialize_HiggsBounds(3, 0, 'LandH')
    whichinput = 'effC'
    Mh = (/125.0D0, 25.0D0, 340.0D0/)

!! set some example non-SM neutral Higgs decays, not mass ordered
    BR_hkhjhi = 0.0D0
    BR_hkhjhi(3, 2, 1) = 0.1D0
    BR_hkhjhi(3, 1, 2) = 0.1D0
    BR_hkhjhi(3, 2, 2) = 0.15D0
    BR_hkhjhi(3, 1, 1) = 0.2D0
    BR_hkhjhi(1, 2, 2) = 0.25D0
    BR_hjhiZ = 0.0D0
    BR_hjhiZ(3, 1) = 0.05D0
    BR_hjhiZ(3, 2) = 0.1D0
    BR_hjhiZ(1, 2) = 0.15D0
    BR_hjemu = 0.0D0
    BR_hjetau = 0.0D0
    BR_hjmutau = 0.0D0
    BR_hjinvisible = (/0.4, 0.5, 0.3/)

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
                     BR_hjinvisible(i) + sumiBR_hihjZ(i) + sumiBR_hkhihj(i)
        GammaTotal(i) = SMGamma_h(Mh(i))/(1 - nonSMBR(i))
    enddo

    call HiggsBounds_neutral_input_properties(Mh, GammaTotal)
    call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi, BR_hjhiZ, &
                                           BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW)
    call HB5_complete_theo

    if (abs(BR_hjinvisible(2) - theo(1)%BR_hjinvisible(2)) > 1e-6) then
        print *, "Test 1: ", BR_hjinvisible(2), "!=", theo(1)%BR_hjinvisible(2)
        stop 1
    endif

    testval = BR_hjinvisible(1) + & ! direct
              BR_Hkhjhi(1, 2, 2)*theo(1)%BR_hjinvisible(2)**2 + & ! 2->11
              BR_hjhiZ(1, 2)*theo(1)%BR_hjinvisible(2)*0.2 ! 2->1Z
    if (abs(testval - theo(1)%BR_hjinvisible(1)) > 1e-6) then
        print *, "Test 2: ", testval, "!=", theo(1)%BR_hjinvisible(1)
        stop 1
    endif

    testval = BR_hjinvisible(3) + & ! direct
              BR_Hkhjhi(3, 1, 1)*theo(1)%BR_hjinvisible(1)**2 + & ! 3->11
              BR_Hkhjhi(3, 2, 2)*theo(1)%BR_hjinvisible(2)**2 + & ! 3->22
              BR_hkhjhi(3, 1, 2)*theo(1)%BR_hjinvisible(1)*theo(1)%BR_hjinvisible(2) + & ! 3->12
              BR_hjhiZ(3, 1)*theo(1)%BR_hjinvisible(1)*0.2 + & ! 3->1Z
              BR_hjhiZ(3, 2)*theo(1)%BR_hjinvisible(2)*0.2 ! 3->2Z
    if (abs(testval - theo(1)%BR_hjinvisible(3)) > 1e-6) then
        print *, "Test 3: ", testval, "!=", theo(1)%BR_hjinvisible(3)
        stop 1
    endif

!! now check the trivial case
    BR_hkhjhi = 0.0D0
    BR_hkhjhi(3, 1, 2) = 0.1D0
    BR_hkhjhi(3, 2, 1) = 0.1D0
    BR_hkhjhi(3, 1, 1) = 0.15D0
    BR_hkhjhi(3, 2, 2) = 0.2D0
    BR_hkhjhi(2, 1, 1) = 0.25D0
    BR_hjhiZ = 0.0D0
    BR_hjhiZ(3, 2) = 0.05D0
    BR_hjhiZ(3, 1) = 0.1D0
    BR_hjhiZ(2, 1) = 0.15D0
    BR_hjemu = 0.0D0
    BR_hjetau = 0.0D0
    BR_hjmutau = 0.0D0
    BR_hjinvisible = (/0.0, 0.0, 0.0/)

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
                     BR_hjinvisible(i) + sumiBR_hihjZ(i) + sumiBR_hkhihj(i)
        GammaTotal(i) = SMGamma_h(Mh(i))/(1 - nonSMBR(i))
    enddo

    call HiggsBounds_neutral_input_properties(Mh, GammaTotal)
    call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi, BR_hjhiZ, &
                                           BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW)
    call HB5_complete_theo

    if (abs(BR_hjinvisible(1)) > 1e-6) then
        print *, "Test 4: ", BR_hjinvisible(1), "!= 0"
        stop 1
    endif

    testval = BR_hjinvisible(2) + & ! direct
              BR_Hkhjhi(2, 1, 1)*BR_hjinvisible(1)**2 + & ! 2->11
              BR_hjhiZ(2, 1)*BR_hjinvisible(1)*0.2 ! 2->1Z
    if (abs(testval) > 1e-6) then
        print *, "Test 5: ", testval, "!= 0"
        stop 1
    endif

    testval = BR_hjinvisible(3) + & ! direct
              BR_Hkhjhi(3, 1, 1)*BR_hjinvisible(1)**2 + & ! 3->11
              BR_Hkhjhi(3, 2, 2)*BR_hjinvisible(2)**2 + & ! 3->22
              BR_hkhjhi(3, 1, 2)*BR_hjinvisible(1)*BR_hjinvisible(2) + & ! 3->12
              BR_hjhiZ(3, 1)*BR_hjinvisible(1)*0.2 + & ! 3->1Z
              BR_hjhiZ(3, 2)*BR_hjinvisible(2)*0.2 ! 3->2Z
    if (abs(testval) > 1e-6) then
        print *, "Test 6: ", testval, "!= 0"
        stop 1
    endif

    call finish_HiggsBounds
end program
