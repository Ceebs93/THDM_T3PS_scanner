! This file is part of SusHi.
! 
! It includes the real corrections to gluon fusion,
! which are called within sigma.f .
! 
C-{{{ coefficients

C-{{{ C1SM

      subroutine C1SMh(s,t,u,m2,c1)
      implicit none
      double complex c1
      double precision s,t,u,m2
      double complex Integral3,Integral4

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      stu4 = Integral4(s,t,u,m2)
      sut4 = Integral4(s,u,t,m2)
      tsu4 = Integral4(t,s,u,m2)
      s3   = Integral3(s,s + t + u,m2)
      t3   = Integral3(t,s + t + u,m2)
      u3   = Integral3(u,s + t + u,m2)

      c1 = 8.d0*(s + t + u) - (-4.d0*m2 + s + t + u) * (
     &     2.d0*((t + u)*s3 + (s + u)*t3 + (s + t)*u3)
     &     + s*u*stu4 + s*t*sut4 + t*u*tsu4 )
      end

      subroutine C1SMA(s,t,u,m2,c1)
      implicit none
      double complex c1
      double precision s,t,u,m2
      double complex Integral3,Integral4

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      stu4 = Integral4(s,t,u,m2)
      sut4 = Integral4(s,u,t,m2)
      tsu4 = Integral4(t,s,u,m2)
      s3   = Integral3(s,s + t + u,m2)
      t3   = Integral3(t,s + t + u,m2)
      u3   = Integral3(u,s + t + u,m2)

      c1 = - (s + t + u) * (
     &       2*( (t + u)*s3 + (s + u)*t3 + (s + t)*u3 )
     &     + s*u*stu4 + s*t*sut4 + t*u*tsu4 )
      end

C-}}}
C-{{{ C1SUSY

      subroutine C1SUSYh(s,t,u,m2,c1)
      implicit none
      double complex c1
      double precision s,t,u,m2
      double complex Integral3,Integral4

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      stu4 = Integral4(s,t,u,m2)
      sut4 = Integral4(s,u,t,m2)
      tsu4 = Integral4(t,s,u,m2)
      s3   = Integral3(s,s + t + u,m2)
      t3   = Integral3(t,s + t + u,m2)
      u3   = Integral3(u,s + t + u,m2)

      c1 = -2.d0*(s + t + u) - m2 * (2.d0 * (
     &     (t + u)*s3 + (s + u)*t3 + (s + t)*u3 )
     &     + s*u*stu4 + s*t*sut4 + t*u*tsu4)
      end

C-}}}
C-{{{ C2SM

      subroutine C2SMh(s,t,u,m2,c2,c3,c4)
      implicit none
      double complex c2,c3,c4,Integral3,C2SMhcalc
      double precision s,t,u,m2

      double complex lstu4,lsut4,ltsu4,ls3,lt3,lu3,ls03,lt03,lu03
      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      lstu4 = stu4
      lsut4 = sut4
      ltsu4 = tsu4
      ls3   = s3
      lt3   = t3
      lu3   = u3
      ls03  = Integral3(0.d0,s,m2)
      lt03  = Integral3(0.d0,t,m2)
      lu03  = Integral3(0.d0,u,m2)

      t03   = lt03
      u03   = lu03
      c2    = C2SMhcalc(s,t,u,m2)

      stu4  = ltsu4
      tsu4  = lstu4
      s3    = lt3
      t3    = ls3
      t03   = ls03
      c3    = C2SMhcalc(t,s,u,m2)

      sut4  = lstu4
      tsu4  = lsut4
      s3    = lu3
      u3    = lt3
      u03   = lt03
      c4    = C2SMhcalc(u,s,t,m2)
      end

      function C2SMhcalc(s,t,u,m2)
      implicit none
      double complex C2SMhcalc,Integral2!,Integral3
      double precision s,t,u,m2

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03,t2,u2,stu2
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      t2 = Integral2(t,m2)
      u2 = Integral2(u,m2)
      stu2 = Integral2(s + t + u,m2)

c$$$      t03 = Integral3(0.d0, t, m2)
c$$$      u03 = Integral3(0.d0, u, m2)
c$$$
c$$$      s3 = Integral3(s,s + t + u, m2)
c$$$      t3 = Integral3(t,s + t + u, m2)
c$$$      u3 = Integral3(u,s + t + u, m2)

      C2SMhcalc=(8.d0*(-((s*t)/(s + t)) + s**2/(s + u) + 
     &     t*u*(((2.d0*s + u)*t2)/(s + u)**2 + 
     &     ((2.d0*s + t)*u2)/(s + t)**2 - 
     &     ((2.d0*s + t)/(s + t)**2 + 
     &        (2.d0*s + u)/(s + u)**2)*
     &      stu2 + ((t*t03 + u*u03))/(2.d0*s))) + 
     &  2.d0*((4.d0*m2 - s)*(t + u)*s3 + 
     &  (-s**2 - 2.d0*(2.d0*m2 + t)*u - (2.d0*t*u**2)/s + 
     &     s*(4.d0*m2 + u - (8.d0*m2*u)/(s + u)))*
     &   t3 + ((-2.d0*t**2*u)/s - (-4.d0*m2 + s)*s + 
     &     t*(-4.d0*m2 - 2.d0*u + s - (8.d0*m2*s)/(s + t)))*
     &   u3) + (4.d0*m2 - s)*s*u*stu4 + (4.d0*m2 - s)*s*t*sut4 + 
     & t*u*(-12.d0*m2 + s - (4.d0*t*u)/s)*tsu4)
      end


      subroutine C2SMA(s,t,u,m2,c2,c3,c4)
      implicit none
      double complex c2,c3,c4,C2SMAcalc
      double precision s,t,u,m2

      double complex lstu4,lsut4,ltsu4,ls3,lt3,lu3
      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      lstu4 = stu4
      lsut4 = sut4
      ltsu4 = tsu4
      ls3   = s3
      lt3   = t3
      lu3   = u3

      c2    = C2SMAcalc(s,t,u,m2)

      stu4  = ltsu4
      tsu4  = lstu4
      s3    = lt3
      t3    = ls3
      c3    = C2SMAcalc(t,s,u,m2)

      sut4  = lstu4
      tsu4  = lsut4
      s3    = lu3
      u3    = lt3
      c4    = C2SMAcalc(u,s,t,m2)
      end

      function C2SMAcalc(s,t,u,m2)
      implicit none
      double complex C2SMAcalc!,Integral3
      double precision s,t,u,m2

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

c$$$      s3 = Integral3(s,s + t + u, m2)
c$$$      t3 = Integral3(t,s + t + u, m2)
c$$$      u3 = Integral3(u,s + t + u, m2)

      C2SMAcalc = -s * (
     &     2 * ( (t + u)*s3 + (s - u)*t3 + (s - t)*u3 )
     &     + s*u*stu4 + s*t*sut4 - t*u*tsu4 )
      end

C-}}}
C-{{{ C2SUSY

      subroutine C2SUSYh(s,t,u,m2,c2,c3,c4)
      implicit none
      double complex c2,c3,c4,Integral3,C2SUSYhcalc
      double precision s,t,u,m2

      double complex lstu4,lsut4,ltsu4,ls3,lt3,lu3,ls03,lt03,lu03
      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      lstu4 = stu4
      lsut4 = sut4
      ltsu4 = tsu4
      ls3   = s3
      lt3   = t3
      lu3   = u3
      ls03  = Integral3(0.d0,s,m2)
      lt03  = Integral3(0.d0,t,m2)
      lu03  = Integral3(0.d0,u,m2)

      t03   = lt03
      u03   = lu03
      c2    = C2SUSYhcalc(s,t,u,m2)

      stu4  = ltsu4
      tsu4  = lstu4
      s3    = lt3
      t3    = ls3
      t03   = ls03
      c3    = C2SUSYhcalc(t,s,u,m2)

      sut4  = lstu4
      tsu4  = lsut4
      s3    = lu3
      u3    = lt3
      u03   = lt03
      c4    = C2SUSYhcalc(u,s,t,m2)
      end

      function C2SUSYhcalc(s,t,u,m2)
      implicit none
      double complex C2SUSYhcalc
      double precision s,t,u,m2
      double complex Integral2!,Integral3

      double complex stu4,sut4,tsu4,s3,t3,u3,t03,u03,u2,t2,stu2
      common /stu4/ stu4,sut4,tsu4,s3,t3,u3,t03,u03

      u2 = Integral2(u,m2)
      t2 = Integral2(t,m2)
      stu2 = Integral2(s + t + u,m2)

c$$$      t03 = Integral3(0.d0, t, m2)
c$$$      u03 = Integral3(0.d0, u, m2)
c$$$
c$$$      s3 = Integral3(s,s + t + u, m2)
c$$$      t3 = Integral3(t,s + t + u, m2)
c$$$      u3 = Integral3(u,s + t + u, m2)

      C2SUSYhcalc=-(2.d0*(-((s*t)/(s + t)) + s**2/(s + u) + 
     &     t*u*(((2.d0*s + u)*t2)/(s + u)**2 + 
     &     ((2.d0*s + t)*u2)/(s + t)**2 - 
     &     ((2.d0*s + t)/(s + t)**2 + 
     &        (2.d0*s + u)/(s + u)**2)*
     &      stu2 + ((t*t03 + u*u03))/(2.d0*s))) + 
     &  (4.d0*m2*(t + u)*s3 + 
     &  (- 2.d0*t*u - (2.d0*t*u**2)/s + 
     &     s*4.d0*m2 - s*(8.d0*m2*u)/(s + u)-4.d0*m2*u)*
     &   t3 + ((-2.d0*t**2*u)/s + 4.d0*m2*s + 
     &     t*(-4.d0*m2 - 2.d0*u - (8.d0*m2*s)/(s + t)))*
     &   u3)/2.d0 + (4.d0*m2*s*u*stu4 + 4.d0*m2*s*t*sut4 + 
     & t*u*(-12.d0*m2 - (4.d0*t*u)/s)*tsu4)/4.d0)
      end

C-}}}
C-{{{ C1eff

      subroutine C1effh(s,t,u,c1)
      implicit none
      double precision s,t,u
      double complex c1
      c1 = 4/3.d0*(s + t + u)**2
      end

      subroutine C1effA(s,t,u,c1)
      implicit none
      double precision s,t,u
      double complex c1
      c1 = 2*(s + t + u)**2
      end

C-}}}
C-{{{ C2eff

      subroutine C2effh(s,t,u,c2,c3,c4)
      implicit none
      double precision s,t,u
      double complex c2,c3,c4
      c2 = 4/3.d0*s**2
      c3 = 4/3.d0*t**2
      c4 = 4/3.d0*u**2
      end

      subroutine C2effA(s,t,u,c2,c3,c4)
      implicit none
      double precision s,t,u
      double complex c2,c3,c4
      c2 = 2*s**2
      c3 = 2*t**2
      c4 = 2*u**2
      end

C-}}}
C-{{{ ASM

      function ASMh(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASMh
      double complex Integral2,Integral3
      ASMh=-3.d0*(2.d0*s*Integral2(s,m2)-2.d0*s*Integral2(s+tpu,m2)+
     &   (tpu)*(-2.d0 +
     &   (tpu - 4.d0*m2)*Integral3(s,s + tpu,m2)))/tpu**2
      end

      function ASMA(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASMA
      double complex Integral3
      ASMA=3.d0*Integral3(s,s + tpu,m2)
      end

C-}}}
C-{{{ ASUSY

      function ASUSYh(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASUSYh
      double complex Integral2,Integral3 
      ASUSYh=3/4.d0*(2*s*Integral2(s,m2)-2*s*Integral2(s+tpu,m2)-
     &   (tpu)*(2 +
     &   4*m2*Integral3(s,s + tpu,m2)))/tpu**2
      end

C-}}}
C-{{{ ASMeff

      function ASMeffh()
      implicit none
      double precision ASMeffh
      ASMeffh=1.d0
      end

      function ASMeffA()
      implicit none
      double precision ASMeffA
      ASMeffA=-3/2.d0
      end

C-}}}
C-{{{ ASUSYeff

      function ASUSYeffh()
      implicit none
      double precision ASUSYeffh
      ASUSYeffh=1/8.d0
      end

C-}}}
C-{{{ ALOSM

      function ALOSMh(mh2,m2)
      implicit none
      double precision mh2,m2
      double complex ALOSMh
      double complex Integral3
      ALOSMh=3*(2+(4*m2-
     &mh2)*Integral3(0.d0,mh2,m2))/mh2
      end

      function ALOSMA(mh2,m2)
      implicit none
      double precision mh2,m2
      double complex ALOSMA
      double complex Integral3
      ALOSMA=-3*Integral3(0.d0,mh2,m2)
      end

C-}}}
C-{{{ ALOSMeff

      function ALOSMeffh()
      implicit none
      double precision ALOSMeffh
      ALOSMeffh=1.d0
      end

      function ALOSMeffA()
      implicit none
      double precision ALOSMeffA
      ALOSMeffA=3/2.d0
      end

C-}}}
C-{{{ ALOSUSY

      function ALOSUSYh(mh2,mb2)
      implicit none
      double complex ALOSUSYh
      double precision mh2,mb2
      double complex Integral3
      ALOSUSYh=-3/2.d0*(1+2*mb2*Integral3(0.d0,mh2,mb2))/mh2
      end

C-}}}
C-{{{ ALOSUSYeff

      function ALOSUSYeffh()
      implicit none
      double precision ALOSUSYeffh
      ALOSUSYeffh=1/8.d0
      end

C-}}}

C-}}}
C-{{{ squared amplidute for gg -> gh

      function AMPgg(s,t,u)
c     squared amplidute for gg -> gh,
c     normalized to the LO amplitude
      implicit none
      double precision AMPgg,AMPggpure,AMPLO,s,t,u
      double complex AMPLOeff

      include '../commons/common-quark.f'

      if(htl) then
         AMPgg = AMPggpure(s,t,u)/cdabs(AMPLOeff(s+t+u))**2
      else
         AMPgg = AMPggpure(s,t,u)/AMPLO(s+t+u)
      endif
      end

      function AMPggpure(s,t,u)
c     squared amplidute for gg -> gh
c     without normalization
      implicit none
      double precision AMPggpure,s,t,u
      double complex c1,c2,c3,c4

      double complex charm1,charm2,charm3,charm4,top1,top2,top3,top4,
     &bottom1,bottom2,bottom3,bottom4,stop11,stop12,stop13,stop14,
     &stop21,stop22,stop23,stop24,sbot11,sbot12,sbot13,sbot14,
     &sbot21,sbot22,sbot23,sbot24,botp1,botp2,botp3,botp4,
     &topp1,topp2,topp3,topp4,stopp11,stopp12,stopp13,stopp14,
     &stopp21,stopp22,stopp23,stopp24,sbotp11,sbotp12,sbotp13,sbotp14,
     &sbotp21,sbotp22,sbotp23,sbotp24,dim51,dim52,dim53,dim54

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'

      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            call C1SMA(s,t,u,mc2,charm1)
            call C2SMA(s,t,u,mc2,charm2,charm3,charm4)
            c1 = c1 + mc2*gc*charm1
            c2 = c2 + mc2*gc*charm2
            c3 = c3 + mc2*gc*charm3
            c4 = c4 + mc2*gc*charm4
         endif

         if(gb.ne.0.d0) then
            call C1SMA(s,t,u,mb2,bottom1)
            call C2SMA(s,t,u,mb2,bottom2,bottom3,bottom4)
            c1 = c1 + mb2*gb*bottom1
            c2 = c2 + mb2*gb*bottom2
            c3 = c3 + mb2*gb*bottom3
            c4 = c4 + mb2*gb*bottom4
         endif

         if(gt.ne.0.d0) then
            if(htl) then
               call C1effA(s,t,u,top1)
               call C2effA(s,t,u,top2,top3,top4)
               c1 = c1 + gt*top1
               c2 = c2 + gt*top2
               c3 = c3 + gt*top3
               c4 = c4 + gt*top4
            else
               call C1SMA(s,t,u,mt2,top1)
               call C2SMA(s,t,u,mt2,top2,top3,top4)
               c1 = c1 + mt2*gt*top1
               c2 = c2 + mt2*gt*top2
               c3 = c3 + mt2*gt*top3
               c4 = c4 + mt2*gt*top4
            endif
         endif

         if(gbp.ne.0.d0) then
            call C1SMA(s,t,u,mbp2,botp1)
            call C2SMA(s,t,u,mbp2,botp2,botp3,botp4)
            c1 = c1 + mbp2*gbp*botp1
            c2 = c2 + mbp2*gbp*botp2
            c3 = c3 + mbp2*gbp*botp3
            c4 = c4 + mbp2*gbp*botp4
         endif

         if(gtp.ne.0.d0) then
            call C1SMA(s,t,u,mtp2,topp1)
            call C2SMA(s,t,u,mtp2,topp2,topp3,topp4)
            c1 = c1 + mtp2*gtp*topp1
            c2 = c2 + mtp2*gtp*topp2
            c3 = c3 + mtp2*gtp*topp3
            c4 = c4 + mtp2*gtp*topp4
         endif

         if(c5sp(1,0).ne.0.d0) then
            call C1effA(s,t,u,dim51)
            call C2effA(s,t,u,dim52,dim53,dim54)
            c1 = c1 + c5sp(1,0)*dim51
            c2 = c2 + c5sp(1,0)*dim52
            c3 = c3 + c5sp(1,0)*dim53
            c4 = c4 + c5sp(1,0)*dim54
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            call C1SMh(s,t,u,mc2,charm1)
            call C2SMh(s,t,u,mc2,charm2,charm3,charm4)
            c1 = c1 + mc2*gc*charm1
            c2 = c2 + mc2*gc*charm2
            c3 = c3 + mc2*gc*charm3
            c4 = c4 + mc2*gc*charm4
         endif

         if(gb.ne.0.d0) then
            call C1SMh(s,t,u,mb2,bottom1)
            call C2SMh(s,t,u,mb2,bottom2,bottom3,bottom4)
            c1 = c1 + mb2*gb*bottom1
            c2 = c2 + mb2*gb*bottom2
            c3 = c3 + mb2*gb*bottom3
            c4 = c4 + mb2*gb*bottom4
         endif

         if(gt.ne.0.d0) then
            if(htl) then
               call C1effh(s,t,u,top1)
               call C2effh(s,t,u,top2,top3,top4)
               c1 = c1 + gt*top1
               c2 = c2 + gt*top2
               c3 = c3 + gt*top3
               c4 = c4 + gt*top4
            else
               call C1SMh(s,t,u,mt2,top1)
               call C2SMh(s,t,u,mt2,top2,top3,top4)
               c1 = c1 + mt2*gt*top1
               c2 = c2 + mt2*gt*top2
               c3 = c3 + mt2*gt*top3
               c4 = c4 + mt2*gt*top4
            endif
         endif

         if(gbp.ne.0.d0) then
            call C1SMh(s,t,u,mbp2,botp1)
            call C2SMh(s,t,u,mbp2,botp2,botp3,botp4)
            c1 = c1 + mbp2*gbp*botp1
            c2 = c2 + mbp2*gbp*botp2
            c3 = c3 + mbp2*gbp*botp3
            c4 = c4 + mbp2*gbp*botp4
         endif

         if(gtp.ne.0.d0) then
            call C1SMh(s,t,u,mtp2,topp1)
            call C2SMh(s,t,u,mtp2,topp2,topp3,topp4)
            c1 = c1 + mtp2*gtp*topp1
            c2 = c2 + mtp2*gtp*topp2
            c3 = c3 + mtp2*gtp*topp3
            c4 = c4 + mtp2*gtp*topp4
         endif

         if(gb1.ne.0.d0) then
            call C1SUSYh(s,t,u,msb12,sbot11)
            call C2SUSYh(s,t,u,msb12,sbot12,sbot13,sbot14)
            c1 = c1 + mbsb2*gb1*sbot11
            c2 = c2 + mbsb2*gb1*sbot12
            c3 = c3 + mbsb2*gb1*sbot13
            c4 = c4 + mbsb2*gb1*sbot14
            call C1SUSYh(s,t,u,msb22,sbot21)
            call C2SUSYh(s,t,u,msb22,sbot22,sbot23,sbot24)
            c1 = c1 + mbsb2*gb2*sbot21
            c2 = c2 + mbsb2*gb2*sbot22
            c3 = c3 + mbsb2*gb2*sbot23
            c4 = c4 + mbsb2*gb2*sbot24
         endif

         if(gt1.ne.0.d0) then
            call C1SUSYh(s,t,u,mst12,stop11)
            call C2SUSYh(s,t,u,mst12,stop12,stop13,stop14)
            c1 = c1 + mt2*gt1*stop11
            c2 = c2 + mt2*gt1*stop12
            c3 = c3 + mt2*gt1*stop13
            c4 = c4 + mt2*gt1*stop14
            call C1SUSYh(s,t,u,mst22,stop21)
            call C2SUSYh(s,t,u,mst22,stop22,stop23,stop24)
            c1 = c1 + mt2*gt2*stop21
            c2 = c2 + mt2*gt2*stop22
            c3 = c3 + mt2*gt2*stop23
            c4 = c4 + mt2*gt2*stop24
         endif

         if(gbp1.ne.0.d0) then
            call C1SUSYh(s,t,u,msbp12,sbotp11)
            call C2SUSYh(s,t,u,msbp12,sbotp12,sbotp13,sbotp14)
            c1 = c1 + mbp2*gt1*sbotp11
            c2 = c2 + mbp2*gt1*sbotp12
            c3 = c3 + mbp2*gt1*sbotp13
            c4 = c4 + mbp2*gt1*sbotp14
            call C1SUSYh(s,t,u,msbp22,sbotp21)
            call C2SUSYh(s,t,u,msbp22,sbotp22,sbotp23,sbotp24)
            c1 = c1 + mbp2*gt2*sbotp21
            c2 = c2 + mbp2*gt2*sbotp22
            c3 = c3 + mbp2*gt2*sbotp23
            c4 = c4 + mbp2*gt2*sbotp24
         endif

         if(gtp1.ne.0.d0) then
            call C1SUSYh(s,t,u,mstp12,stopp11)
            call C2SUSYh(s,t,u,mstp12,stopp12,stopp13,stopp14)
            c1 = c1 + mtp2*gtp1*stopp11
            c2 = c2 + mtp2*gtp1*stopp12
            c3 = c3 + mtp2*gtp1*stopp13
            c4 = c4 + mtp2*gtp1*stopp14
            call C1SUSYh(s,t,u,mstp22,stopp21)
            call C2SUSYh(s,t,u,mstp22,stopp22,stopp23,stopp24)
            c1 = c1 + mtp2*gtp2*stopp21
            c2 = c2 + mtp2*gtp2*stopp22
            c3 = c3 + mtp2*gtp2*stopp23
            c4 = c4 + mtp2*gtp2*stopp24
         endif

         if(c5sp(0,0).ne.0.d0) then
            call C1effh(s,t,u,dim51)
            call C2effh(s,t,u,dim52,dim53,dim54)
            c1 = c1 + c5sp(0,0)*dim51
            c2 = c2 + c5sp(0,0)*dim52
            c3 = c3 + c5sp(0,0)*dim53
            c4 = c4 + c5sp(0,0)*dim54
         endif

C-}}}

      endif
      AMPggpure = (cdabs(c1)**2+cdabs(c2)**2+cdabs(c3)**2+cdabs(c4)**2)
     &     /s/t/u*9/32.d0/(s+t+u)
      end

C-}}}
C-{{{ squared amplidute for qg -> qh

      function AMPqg(s,t,u)
c     squared amplidute for qg -> qh,
c     normalized to the LO amplitude
      implicit none
      double precision AMPqg,AMPqgpure,AMPLO,s,t,u
      double complex AMPLOeff

      include '../commons/common-quark.f'

      if(htl) then
         AMPqg=AMPqgpure(s,t,u)/cdabs(AMPLOeff(s+t+u))**2
      else
         AMPqg=AMPqgpure(s,t,u)/AMPLO(s+t+u)
      endif
      end

      function AMPqgpure(s,t,u)
c     squared amplidute for qg -> qh
c     without normalization
      implicit none
      double precision AMPqgpure,s,t,u
      double complex A
      double precision ASMeffh,ASMeffA
      double complex ASMh,ASMA,ASUSYh

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'      

      A = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMA(t,s+u,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMA(t,s+u,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffA()
            else
               A = A + mt2*gt*ASMA(t,s+u,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMA(t,s+u,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMA(t,s+u,mtp2)
         endif
         if(c5sp(1,0).ne.0.d0) then
               A = A + c5sp(1,0)*ASMeffA()
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMh(t,s+u,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMh(t,s+u,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffh()
            else
               A = A + mt2*gt*ASMh(t,s+u,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMh(t,s+u,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMh(t,s+u,mtp2)
         endif
         if(gb1.ne.0.d0) then
            A = A + mbsb2*(gb1*ASUSYh(t,s+u,msb12)
     &                 + gb2*ASUSYh(t,s+u,msb22))
         endif
         if(gt1.ne.0.d0) then
            A = A + mt2*(gt1*ASUSYh(t,s+u,mst12)
     &                 + gt2*ASUSYh(t,s+u,mst22))
         endif
         if(gbp1.ne.0.d0) then
            A = A + mbp2*(gbp1*ASUSYh(t,s+u,msbp12)
     &                  + gbp2*ASUSYh(t,s+u,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            A = A + mtp2*(gt1*ASUSYh(t,s+u,mstp12)
     &                  + gt2*ASUSYh(t,s+u,mstp22))
         endif
         if(c5sp(0,0).ne.0.d0) then
               A = A + c5sp(0,0)*ASMeffh()
         endif

C-}}}

      endif
      AMPqgpure = -(s**2+u**2)/t*cdabs(A)**2/(s+t+u)
      end

C-}}}
C-{{{ squared amplidute for qq -> gh

      function AMPqq(s,tpu)
c     squared amplidute for qq -> gh,
c     normalized to the LO amplitude
      implicit none
      double precision AMPqq,AMPqqpure,AMPLO,s,tpu
      double complex AMPLOeff

      include '../commons/common-quark.f'

      if(htl) then
         AMPqq = AMPqqpure(s,tpu)/cdabs(AMPLOeff(s+tpu))**2
      else
         AMPqq = AMPqqpure(s,tpu)/AMPLO(s+tpu)
      endif
      end

      function AMPqqpure(s,tpu)
c     squared amplidute for qq -> gh
c     without normalization
      implicit none
      double precision AMPqqpure,s,tpu
      double complex A
      double precision ASMeffh,ASMeffA
      double complex ASMh,ASMA,ASUSYh

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'

      A = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMA(s,tpu,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMA(s,tpu,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffA()
            else
               A = A + mt2*gt*ASMA(s,tpu,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMA(s,tpu,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMA(s,tpu,mtp2)
         endif
         if(c5sp(1,0).ne.0.d0) then
            A = A + c5sp(1,0)*ASMeffA()
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMh(s,tpu,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMh(s,tpu,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffh()
            else
               A = A + mt2*gt*ASMh(s,tpu,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMh(s,tpu,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMh(s,tpu,mtp2)
         endif
         if(gb1.ne.0.d0) then
            A = A + mbsb2*(gb1*ASUSYh(s,tpu,msb12)
     &                   + gb2*ASUSYh(s,tpu,msb22))
         endif
         if(gt1.ne.0.d0) then
            A = A + mt2*(gt1*ASUSYh(s,tpu,mst12)
     &                 + gt2*ASUSYh(s,tpu,mst22))
         endif
         if(gbp1.ne.0.d0) then
            A = A + mbp2*(gbp1*ASUSYh(s,tpu,msbp12)
     &                  + gbp2*ASUSYh(s,tpu,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            A = A + mtp2*(gtp1*ASUSYh(s,tpu,mstp12)
     &                  + gtp2*ASUSYh(s,tpu,mstp22))
         endif
         if(c5sp(0,0).ne.0.d0) then
            A = A + c5sp(0,0)*ASMeffh()
         endif

C-}}}

      endif
      AMPqqpure = cdabs(A)**2
      end

C-}}}
C-{{{ squared amplidute for gg -> h

      function AMPLO(mhiggs2)
      implicit none
      double precision AMPLO,mhiggs2
      double complex AMP
      double precision ALOSMeffh,ALOSMeffA
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMA(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMA(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
c            if(htl) then
c               AMP = AMP + gt*ALOSMeffA()
c            else
               AMP = AMP + mt2*gt*ALOSMA(mhiggs2,mt2)
c            endif
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMA(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMA(mhiggs2,mtp2)
         endif
         if(c5sp(1,0).ne.0.d0) then
            AMP = AMP + c5sp(1,0)*ALOSMeffA()
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMh(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMh(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
c            if(htl) then
c               AMP = AMP + gt*ALOSMeffh()
c            else
               AMP = AMP + mt2*gt*ALOSMh(mhiggs2,mt2)
c            endif
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMh(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMh(mhiggs2,mtp2)
         endif
         if(gb1.ne.0.d0) then
            AMP = AMP + mbsb2*(gb1*ALOSUSYh(mhiggs2,msb12)
     &                       + gb2*ALOSUSYh(mhiggs2,msb22))
         endif
         if(gt1.ne.0.d0) then
            AMP = AMP + mt2*(gt1*ALOSUSYh(mhiggs2,mst12)
     &                     + gt2*ALOSUSYh(mhiggs2,mst22))
         endif
         if(gbp1.ne.0.d0) then
            AMP = AMP + mbp2*(gbp1*ALOSUSYh(mhiggs2,msbp12)
     &                      + gbp2*ALOSUSYh(mhiggs2,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            AMP = AMP + mtp2*(gtp1*ALOSUSYh(mhiggs2,mstp12)
     &                      + gtp2*ALOSUSYh(mhiggs2,mstp22))
         endif
         if(c5sp(0,0).ne.0.d0) then
            AMP = AMP + c5sp(0,0)*ALOSMeffh()
         endif

C-}}}

      endif

      AMPLO=cdabs(AMP)**2

      end

C-}}}
C-{{{ squared amplidute for gg -> h

      function AMPLOpure(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMP,AMPLOpure
      double precision ALOSMeffh,ALOSMeffA
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMA(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMA(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               AMP = AMP + gt*ALOSMeffA()
            else
               AMP = AMP + mt2*gt*ALOSMA(mhiggs2,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMA(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMA(mhiggs2,mtp2)
         endif
         if(c5sp(1,0).ne.0.d0) then
            AMP = AMP + c5sp(1,0)*ALOSMeffA()
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMh(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMh(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               AMP = AMP + gt*ALOSMeffh()
            else
               AMP = AMP + mt2*gt*ALOSMh(mhiggs2,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMh(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMh(mhiggs2,mtp2)
         endif
         if(gb1.ne.0.d0) then
            AMP = AMP + mbsb2*(gb1*ALOSUSYh(mhiggs2,msb12)
     &                       + gb2*ALOSUSYh(mhiggs2,msb22))
         endif
         if(gt1.ne.0.d0) then
            AMP = AMP + mt2*(gt1*ALOSUSYh(mhiggs2,mst12)
     &                     + gt2*ALOSUSYh(mhiggs2,mst22))
         endif
         if(gbp1.ne.0.d0) then
            AMP = AMP + mbp2*(gbp1*ALOSUSYh(mhiggs2,msbp12)
     &                      + gbp2*ALOSUSYh(mhiggs2,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            AMP = AMP + mtp2*(gtp1*ALOSUSYh(mhiggs2,mstp12)
     &                      + gtp2*ALOSUSYh(mhiggs2,mstp22))
         endif
         if(c5sp(0,0).ne.0.d0) then
            AMP = AMP + c5sp(0,0)*ALOSMeffh()
         endif

C-}}}

      endif

      AMPLOpure=AMP

      end

C-}}}
C-{{{ amplidute for gg -> h eff

      function AMPLOeff(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMPLOeff
      double precision ALOSUSYeffh,ALOSMeffh,ALOSMeffA
      double complex ALOSMh,ALOSMA

      include '../commons/common-quark.f'
      include '../commons/common-vars.f'

      if(pseudo.eq.1) then

         if(htl) then
            AMPLOeff = mc2 * gc * ALOSMA(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMA(mhiggs2,mb2)
     &            +          gt * ALOSMeffA()
     &            +    mbp2* gbp* ALOSMA(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMA(mhiggs2,mtp2)
     &            +          c5sp(1,0)* ALOSMeffA()
         else
            AMPLOeff = mc2 * gc * ALOSMA(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMA(mhiggs2,mb2)
     &            +    mt2 * gt * ALOSMA(mhiggs2,mt2)
     &            +    mbp2* gbp* ALOSMA(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMA(mhiggs2,mtp2)
     &            +          c5sp(1,0)* ALOSMeffA()
         endif

      else

         if(htl) then
            AMPLOeff = mc2 * gc * ALOSMh(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMh(mhiggs2,mb2)
     &            +          gt * ALOSMeffh()
     &            +    mbp2* gbp* ALOSMh(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMh(mhiggs2,mtp2)
     &            + gtpe + gbpe
     &            +          c5sp(0,0)* ALOSMeffh()
            if(gt1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mt2 * (gt1/mst12+gt2/mst22)*ALOSUSYeffh()
            endif
            if(gb1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mbsb2 * (gb1/msb12+gb2/msb22)*ALOSUSYeffh()
            endif
         else
            AMPLOeff = mc2 * gc * ALOSMh(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMh(mhiggs2,mb2)
     &            +    mt2 * gt * ALOSMh(mhiggs2,mt2)
     &            +    mbp2* gbp* ALOSMh(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMh(mhiggs2,mtp2)
     &            + gtpe + gbpe
     &            +          c5sp(0,0)* ALOSMeffh()
            if(gt1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mt2 * (gt1/mst12+gt2/mst22)*ALOSUSYeffh()
            endif
            if(gb1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mbsb2 *(gb1/msb12+gb2/msb22)*ALOSUSYeffh()
            endif
         endif
      endif
      end

C-}}}
