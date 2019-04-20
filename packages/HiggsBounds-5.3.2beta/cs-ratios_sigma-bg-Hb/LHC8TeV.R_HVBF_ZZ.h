!******************************************************
        real*8 function lhc8_rHVBF_ZZ(x)
!******************************************************
!* Used HAWK-2.0.1.
!* (PDF: NNPDF23_nlo_as_0118_qed)
!* Calculated Born CS
!* (a) with only SM HWW coupling,
!* (b) with only SM HZZ coupling,
!* (c) with both SM HWW and HZZ coupling.
!* Calculated ratios (a)/((a)+(b)) and (b)/((a)+(b)), and
!* performed simultaneous fit, such that the sum is 1.
!* Cross checked that (a)+(b) ~= (c), therefore interference
!* is negligible.
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [10:1000].
!* 01/06/2016, TS
!******************************
        implicit none
        real*8 am1,a0,a1,a2,x

        am1 = 9.99844541408024e-06
        a0 = 0.250107980787387
        a1 = 6.52202639729497e-05
        a2 = -3.20376357582592e-08

        lhc8_rHVBF_ZZ=0d0

        if(x .ge. 0d0) then
        lhc8_rHVBF_ZZ= am1/x + a0 + a1*x + a2*x**2.0d0
        endif
        
        if(x .gt. 1050d0) then
!        write(*,*)'function lhc8_rHVBF_ZZ might not be a good fit (m_H > 1050 GeV)'
        endif

        end function
        
