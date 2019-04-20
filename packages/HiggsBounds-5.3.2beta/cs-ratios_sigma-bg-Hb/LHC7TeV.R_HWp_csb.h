!******************************************************
        real*8 function lhc7_rHWp_csb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.4%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,x,log_lhc7_rHWp_csb

        a1 = -2.6362541789099
        b1 = -0.00532565578631144
        c1 = 4.40029490664571e-06
        d1 = -2.32475135442569e-09

        log_lhc7_rHWp_csb=0d0
        
        if(x .lt. 80d0) then
!         write(*,*)'function lhc7_rHWp_csb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc7_rHWp_csb=a1+b1*x+c1*x**2+d1*x**3
        lhc7_rHWp_csb=exp(log_lhc7_rHWp_csb)
        endif
        if(x .gt. 310d0) then
!         write(*,*)'function lhc7_rHWp_csb might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                        
        end function
        
