!******************************************************
        real*8 function lhc8_rHZ_bbb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.2%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,x,log_lhc8_rHZ_bbb

        a1 = -3.82303073904333D0
        b1 = -0.00155195945974061D0
        c1 = -3.43240362660468D-05
        d1 = 1.5238122765212D-07
        e1 = -1.96173770125554D-10


        lhc8_rHZ_bbb=0d0
        
        if(x .lt. 80d0) then
!         write(*,*)'function lhc8_rHZ_bbb might not be a good fit (m_H < 80 GeV)'
!--Update (16/09/2014, TS): set value to value at lowest allowed mass
        x=80d0
        log_lhc8_rHZ_bbb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc8_rHZ_bbb=exp(log_lhc8_rHZ_bbb)
!--End update        
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc8_rHZ_bbb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc8_rHZ_bbb=exp(log_lhc8_rHZ_bbb)
        endif
        if(x .gt. 310d0) then
!         write(*,*)'function lhc8_rHZ_bbb might not be a good fit (m_H > 310 GeV)'
!--Update (16/09/2014, TS): set value to value at highest allowed mass
        x=310d0
        log_lhc8_rHZ_bbb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc8_rHZ_bbb=exp(log_lhc8_rHZ_bbb)
!--End update                
        endif
                                                                                                                                        
        end function
        
