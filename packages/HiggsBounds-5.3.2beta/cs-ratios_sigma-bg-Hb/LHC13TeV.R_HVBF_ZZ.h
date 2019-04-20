!******************************************************
        real*8 function lhc13_rHVBF_ZZ(x)
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
!* fit: strictly valid in range [10:3000].
!* 01/06/2016, TS
!******************************
        implicit none
        real*8 a0,a1,a2,x

        a0 = 0.263427271363775
        a1 = 2.72728035702793e-05
        a2 = -3.89549044847535e-09
        
        lhc13_rHVBF_ZZ=0d0

        if(x .ge. 0d0) then
        lhc13_rHVBF_ZZ= a0 + a1*x + a2*x**2.0d0
        endif
        
        if(x .gt. 3050d0) then
!        write(*,*)'function lhc13_rHVBF_ZZ might not be a good fit (m_H > 3050 GeV)'
        endif

        end function
        
