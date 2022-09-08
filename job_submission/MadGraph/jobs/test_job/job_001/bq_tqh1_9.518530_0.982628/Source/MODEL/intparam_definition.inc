ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        MDL_TB__EXP__2 = MDL_TB**2

        MDL_CB = SQRT(DCMPLX(1.000000D+00/(1.000000D+00+MDL_TB__EXP__2)
     $   ))

        MDL_SINBMA__EXP__2 = MDL_SINBMA**2

        MDL_COSBMA = SQRT(DCMPLX(1.000000D+00-MDL_SINBMA__EXP__2))

        MDL_L7 = MDL_LR7

        MDL_LI5 = 0.000000D+00

        MDL_LI6 = 0.000000D+00

        MDL_MZ__EXP__2 = MDL_MZ**2

        MDL_MZ__EXP__4 = MDL_MZ**4

        MDL_SQRT__2 = SQRT(DCMPLX(2.000000D+00))

        MDL_CB__EXP__2 = MDL_CB**2

        MDL_SB = SQRT(DCMPLX(1.000000D+00-MDL_CB__EXP__2))

        MDL_CA = MDL_CB*MDL_COSBMA+MDL_SB*MDL_SINBMA

        MDL_CA__EXP__2 = MDL_CA**2

        MDL_SA = SQRT(DCMPLX(1.000000D+00-MDL_CA__EXP__2))

        MDL_MHC__EXP__2 = MDL_MHC**2

        MDL_MH1__EXP__2 = MDL_MH1**2

        MDL_MH2__EXP__2 = MDL_MH2**2

        MDL_MH3__EXP__2 = MDL_MH3**2

        MDL_ASIN__SINBMA = ASIN(MDL_SINBMA)

        MDL_COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

        MDL_COSBMA__EXP__2 = MDL_COSBMA**2

        MDL_COSBMA__EXP__4 = MDL_COSBMA**4

        MDL_COSBMA__EXP__3 = MDL_COSBMA**3

        MDL_SINBMA__EXP__3 = MDL_SINBMA**3

        MDL_SINBMA__EXP__4 = MDL_SINBMA**4

        MDL_AEW = 1.000000D+00/AEWM1

        MDL_MW = SQRT(DCMPLX(MDL_MZ__EXP__2/2.000000D+00
     $   +SQRT(DCMPLX(MDL_MZ__EXP__4/4.000000D+00-(MDL_AEW*PI
     $   *MDL_MZ__EXP__2)/(MDL_GF*MDL_SQRT__2)))))

        MDL_SQRT__AEW = SQRT(DCMPLX(MDL_AEW))

        MDL_EE = 2.000000D+00*MDL_SQRT__AEW*SQRT(DCMPLX(PI))

        MDL_MW__EXP__2 = MDL_MW**2

        MDL_SW2 = 1.000000D+00-MDL_MW__EXP__2/MDL_MZ__EXP__2

        MDL_CW = SQRT(DCMPLX(1.000000D+00-MDL_SW2))

        MDL_SQRT__SW2 = SQRT(DCMPLX(MDL_SW2))

        MDL_SW = MDL_SQRT__SW2

        MDL_SW__EXP__2 = MDL_SW**2

        MDL_EE__EXP__2 = MDL_EE**2

        MDL_MU2 = MDL_MHC__EXP__2-(2.000000D+00*MDL_L3*MDL_MW__EXP__2
     $   *MDL_SW__EXP__2)/MDL_EE__EXP__2

        MDL_L1 = -(MDL_EE__EXP__2*(-(MDL_MH1__EXP__2
     $   *MDL_SINBMA__EXP__2)-MDL_MH2__EXP__2*(1.000000D+00
     $   -MDL_SINBMA__EXP__2)))/(8.000000D+00*MDL_MW__EXP__2
     $   *MDL_SW__EXP__2)

        MDL_L4 = (MDL_EE__EXP__2*(2.000000D+00*MDL_MH1__EXP__2
     $   +2.000000D+00*MDL_MH2__EXP__2+4.000000D+00*MDL_MH3__EXP__2
     $   -8.000000D+00*MDL_MHC__EXP__2+2.000000D+00*(-MDL_MH1__EXP__2
     $   +MDL_MH2__EXP__2)*COS(2.000000D+00*(PI/2.000000D+00
     $   -MDL_ASIN__SINBMA))))/(1.600000D+01*MDL_MW__EXP__2
     $   *MDL_SW__EXP__2)

        MDL_LR5 = (MDL_EE__EXP__2*(2.000000D+00*(MDL_MH1__EXP__2
     $   +MDL_MH2__EXP__2-2.000000D+00*MDL_MH3__EXP__2)-2.000000D+00
     $   *(MDL_MH1-MDL_MH2)*(MDL_MH1+MDL_MH2)*COS(2.000000D+00*(PI
     $   /2.000000D+00-MDL_ASIN__SINBMA))))/(3.200000D+01
     $   *MDL_MW__EXP__2*MDL_SW__EXP__2)

        MDL_LR6 = (MDL_EE__EXP__2*(-MDL_MH1__EXP__2+MDL_MH2__EXP__2)
     $   *MDL_SINBMA*SQRT(DCMPLX(1.000000D+00-MDL_SINBMA__EXP__2)))
     $   /(4.000000D+00*MDL_MW__EXP__2*MDL_SW__EXP__2)

        MDL_YB = (MDL_EE*MDL_YMB)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YC = (MDL_EE*MDL_YMC)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YDO = (MDL_EE*MDL_YMDO)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YE = (MDL_EE*MDL_YME)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YM = (MDL_EE*MDL_YMM)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YS = (MDL_EE*MDL_YMS)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YT = (MDL_EE*MDL_YMT)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YTAU = (MDL_EE*MDL_YMTAU)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_YUP = (MDL_EE*MDL_YMUP)/(MDL_MW*MDL_SW*MDL_SQRT__2)

        MDL_MU1 = (-4.000000D+00*MDL_L1*MDL_MW__EXP__2*MDL_SW__EXP__2)
     $   /MDL_EE__EXP__2

        MDL_L5 = MDL_COMPLEXI*MDL_LI5+MDL_LR5

        MDL_L6 = MDL_COMPLEXI*MDL_LI6+MDL_LR6

        MDL_MU3 = (-2.000000D+00*MDL_L6*MDL_MW__EXP__2*MDL_SW__EXP__2)
     $   /MDL_EE__EXP__2

        MDL_CW__EXP__2 = MDL_CW**2

        MDL_SW__EXP__3 = MDL_SW**3

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      AS = G**2/4/PI

      MDL_SQRT__AS = SQRT(DCMPLX(AS))

      MDL_G__EXP__2 = G**2

C     
C     Parameters that should be updated for the loops.
C     
C     
C     Definition of the EW coupling used in the write out of aqed
C     
      GAL(1) = 3.5449077018110318D0 / DSQRT(ABS(AEWM1))
      GAL(2) = 1D0

