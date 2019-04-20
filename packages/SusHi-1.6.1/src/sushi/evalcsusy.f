! This file is part of SusHi.
! 
! It includes the squark-quark-gluino contributions
! to gluon fusion in the limit of a vanishing Higgs mass
! in form of the published code evalcsusy.f.
!
C-{{{ subroutine quarkhiggscouplings

      subroutine quarkhiggscoup(beta,Hmx,Amx,gth,gbh
     &,model,thdmv,lpseudo,Sind)
      !This subroutine calculates the Higgs-quark couplings
      !in the MSSM or in the 2-Higgs-Doublet models as defined
      !in 1207.4835, Table 3
      implicit none
      double precision beta, gth, gbh, Hmx(3,3), Amx(3,3)
      integer model, thdmv, lpseudo, Sind
      double precision tanb, salpha,calpha,sbeta,cbeta

      tanb = dtan(beta)
      sbeta = dsin(beta)
      cbeta = dcos(beta)

      !SM
      if (model.eq.0) then
         gth = 1.d0
         gbh = 1.d0
      endif

      !MSSM or NMSSM
      if ((model.eq.1).or.(model.eq.3)) then
      if (lpseudo.eq.0) then         !light and heavy Higgs
         gth = Hmx(Sind,2)/sbeta
         gbh = Hmx(Sind,1)/cbeta
      else if (lpseudo.eq.1) then    !pseudoscalar Higgs
         gth = Amx(Sind,2)*1.d0/tanb
         gbh = Amx(Sind,2)*tanb
      endif
      endif

      !2HDM
      if (model.eq.2) then
       if (lpseudo.eq.0) then
         gth = Hmx(Sind,2)/sbeta
       else if (lpseudo.eq.1) then
         gth = 1.d0/tanb
       endif

       if ((thdmv.eq.1).or.(thdmv.eq.4)) then
        if (lpseudo.eq.0) then
           gbh = Hmx(Sind,2)/sbeta
        else if (lpseudo.eq.1) then
           gbh = -1.d0/tanb
        endif
       else if ((thdmv.eq.2).or.(thdmv.eq.3)) then
        if (lpseudo.eq.0) then
           gbh = Hmx(Sind,1)/cbeta
        else if (lpseudo.eq.1) then
           gbh = tanb
        endif
       else
         write(*,*) "2HDM version unknown!"
         stop
       endif
      endif

      end

C-}}}
C-{{{ subroutine higgssquarkcouplingsNMSSM

      subroutine squarkhiggscoupNMSSM(Mbsb,Msbot1,Msbot2,
     &Mt,Mstop1,Mstop2,beta,Mw,Mz,GF,Hmx,Amx,lam,muSUSY,
     &cthetat,sthetat,cthetab,sthetab,gbh11,gbh12,gbh21,gbh22,
     &gth11,gth12,gth21,gth22,lpseudo,Sind)
      !This subroutine calculates the Higgs-squark couplings
      !in the NMSSM
      !NB: Amx is prerotated!
      !2014.09.29: Stefan Liebler
      !input
      double precision Mbsb, Msbot1, Msbot2, Mt, Mstop1, Mstop2
      double precision beta, Mw, Mz, GF, lam, muSUSY
      double precision Hmx(3,3), Amx(3,3)
      double precision cthetat, sthetat, cthetab, sthetab
      integer lpseudo, Sind
      !output
      double precision gbh11, gbh12, gbh21, gbh22
      double precision gth11, gth12, gth21, gth22
      !internal
      double precision stht, ctht, s2tht, c2tht, stht2, ctht2
      double precision stht3, ctht3, sthb, cthb
      double precision s2thb, c2thb, sthb2, cthb2
      double precision SiB, CoB, TaB, Si2B, SiB2, CoB2
      double precision msb12, msb22, mst12, mst22, ctW, stW, vevS
      double precision mbsb2, mt2, s2W, sW2, cW2, c2W, mZ2
      double precision Ab11Amx(3), Ab12Amx(3), Ab21Amx(3), Ab22Amx(3)
      double precision At11amx(3), At12Amx(3), At21Amx(3), At22Amx(3)
      double precision Hb11Hmx(3), Hb12Hmx(3), Hb21Hmx(3), Hb22Hmx(3)
      double precision Ht11Hmx(3), Ht12Hmx(3), Ht21Hmx(3), Ht22Hmx(3)

      double precision Sqrt2
      double precision vev 

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      mbsb2 = mbsb*mbsb
      mt2 = mt*mt
      mZ2 = mZ*mZ

      stht = sthetat
      ctht = cthetat
      s2tht = 2.d0*stht*ctht
      c2tht = ctht**2 - stht**2
      stht2 = stht*stht
      ctht2 = ctht*ctht
      stht3 = stht2*stht
      ctht3 = ctht2*ctht
      sthb = sthetab
      cthb = cthetab
      s2thb = 2.d0*sthb*cthb
      c2thb = cthb**2 - sthb**2
      sthb2 = sthb*sthb
      cthb2 = cthb*cthb

      SiB = dsin(beta)
      CoB = dcos(beta)
      TaB = dtan(beta)
      Si2B = dsin(2.d0*beta)
      SiB2 = SiB*SiB
      CoB2 = CoB*CoB

      ctW = Mw/Mz
      stW = dsqrt(1.d0-ctW**2)
      sW2 = stW*stW
      cW2 = ctW*ctW
      c2W = dcos(2.d0*dacos(ctW))
      s2W = dsin(2.d0*dacos(ctW))

      mst12 = Mstop1*Mstop1
      mst22 = Mstop2*Mstop2
      msb12 = Msbot1*Msbot1
      msb22 = Msbot2*Msbot2

!      write(*,*) 'NMSSM Higgs Squark couplings'
!      write(*,*) 'vev,vevS',vev,vevS
!      write(*,*) 'stht,ctht,sthb,cthb',stht,ctht,sthb,cthb
!      write(*,*) 'SiB,ctW',SiB,ctW
!      write(*,*) 'Stop masses',Mstop1,Mstop2
!      write(*,*) 'Sbottom masses',Msbot1,Msbot2
!      write(*,*) 'Hmx',Hmx
!      write(*,*) 'Amx',Amx

      if (lpseudo.eq.0) then

      Hb11Hmx(1) = (-1/3.d0/CoB*(cthb2*(-6.*mbsb2 
     & + mZ2*CoB2*(2. + c2W)) + 2.*sthb2*(-3.*mbsb2 
     & + mZ2*CoB2*sW2) - 3.*cthb*sthb*((msb12 - msb22)*s2thb 
     & + lam*mbsb*Sqrt2*vevS*TaB)))/mbsb2
      Hb11Hmx(2) = (1/3.d0*(mZ2*cthb2*(2. + c2W)*SiB 
     & - 3.*lam*mbsb*Sqrt2*vevS*cthb/CoB*sthb
     & + 2.*mZ2*SiB*sthb2*sW2))/mbsb2
      Hb11Hmx(3) = (-1.*lam*Sqrt2*vev*cthb*sthb*TaB)/mbsb

      Hb12Hmx(1) = (1/6.d0/CoB*(mZ2*CoB2*(1. + 2.*c2W)*s2thb 
     & + 3.*c2thb*((msb12 - msb22)*s2thb 
     & + lam*mbsb*Sqrt2*vevS*TaB)))/mbsb2
      Hb12Hmx(2) = (1/6.d0/CoB*(-3.*lam*mbsb*Sqrt2*vevS*cthb2 
     & - mZ2*cthb*(1. + 2.*c2W)*Si2B*sthb 
     & + 3.*lam*mbsb*Sqrt2*vevS*sthb2))/mbsb2
      Hb12Hmx(3) = (-1/Sqrt2*lam*vev*c2thb*TaB)/mbsb

      Hb21Hmx(1) = (1/6.d0/CoB*(mZ2*CoB2*(1. + 2.*c2W)*s2thb 
     & + 3.*c2thb*((msb12 - msb22)*s2thb 
     & + lam*mbsb*Sqrt2*vevS*TaB)))/mbsb2
      Hb21Hmx(2) = (1/6.d0/CoB*(-3.*lam*mbsb*Sqrt2*vevS*cthb2 
     & - mZ2*cthb*(1. + 2.*c2W)*Si2B*sthb
     & + 3.*lam*mbsb*Sqrt2*vevS*sthb2))/mbsb2
      Hb21Hmx(3) = (-1/Sqrt2*lam*vev*c2thb*TaB)/mbsb

      Hb22Hmx(1) = (1/3.d0*(CoB*sthb2*(-3.*mZ2*cW2 
     & + 6.*mbsb2/CoB2 - mZ2*sW2) 
     & + cthb2*(6.*mbsb2/CoB - 2.*mZ2*CoB*sW2) 
     & - 3.*cthb/CoB*sthb*((msb12 - msb22)*s2thb
     & + lam*mbsb*Sqrt2*vevS*TaB)))/mbsb2
      Hb22Hmx(2) = (1/3.d0*(3.*lam*mbsb*Sqrt2*vevS*cthb/CoB*sthb
     & + mZ2*(2. + c2W)*SiB*sthb2 + 2.*mZ2*cthb2*SiB*sW2))/mbsb2
      Hb22Hmx(3) = (lam*Sqrt2*vev*cthb*sthb*TaB)/mbsb

      Ht11Hmx(1) = (1/6.d0*(-mZ2*CoB*(-3. + c2tht - 4.*c2tht*c2W) 
     & - 6.*lam*mt*Sqrt2*vevS*ctht/SiB*stht))/mt2
      Ht11Hmx(2) = (1/3.d0*(3.*lam*mt*Sqrt2*vevS*ctht/TaB/SiB*stht 
     & + ctht2*(-mZ2*(1. + 2.*c2W)*SiB + 6./SiB*(mt2
     & + (mst12 - mst22)*stht2)) 
     & + 2.*SiB*stht2*(3.*mt2/SiB2 - 2.*mZ2*sW2)))/mt2
      Ht11Hmx(3) = (-1.*lam*Sqrt2*vev*ctht/TaB*stht)/mt

      Ht12Hmx(1) = (1/6.d0/SiB*(-3.*lam*mt*Sqrt2*vevS*ctht2
     & - mZ2*ctht*(-1. + 4.*c2W)*Si2B*stht
     & + 3.*lam*mt*Sqrt2*vevS*stht2))/mt2
      Ht12Hmx(2) = (1/6.d0/SiB*(3.*lam*mt*Sqrt2*vevS*ctht2/TaB
     & + 6.*(mst12 - mst22)*ctht3*stht 
     & - 3.*lam*mt*Sqrt2*vevS/TaB*stht2
     & + 6.*ctht*(mZ2*cW2*SiB2*stht
     & + (mst22 - mst12)*stht3)
     & - 5.*mZ2*SiB2*s2tht*sW2))/mt2
      Ht12Hmx(3) = (-1/Sqrt2*lam*vev*c2tht/TaB)/mt

      Ht21Hmx(1) = (1/6.d0/SiB*(-3.*lam*mt*Sqrt2*vevS*ctht2
     & - mZ2*ctht*(-1. + 4.*c2W)*Si2B*stht
     & + 3.*lam*mt*Sqrt2*vevS*stht2))/mt2
      Ht21Hmx(2) = (1/6.d0/SiB*(3.*lam*mt*Sqrt2*vevS*ctht2/TaB
     & + 6.*(mst12 - mst22)*ctht3*stht
     & - 3.*lam*mt*Sqrt2*vevS/TaB*stht2
     & + 6.*ctht*(mZ2*cW2*SiB2*stht
     & + (mst22 - mst12)*stht3)
     & - 5.*mZ2*SiB2*s2tht*sW2))/mt2
      Ht21Hmx(3) = (-1/Sqrt2*lam*vev*c2tht/TaB)/mt

      Ht22Hmx(1) = (1/3.d0*(3.*lam*mt*Sqrt2*vevS*ctht/SiB*stht
     & + mZ2*CoB*(1. + 2.*c2W)*stht2 + 4.*mZ2*CoB*ctht2*sW2))/mt2
      Ht22Hmx(2) = (1/3.d0*(-3.*ctht/SiB*stht*(lam*mt*Sqrt2*vevS/TaB
     & + (mst12 - mst22)*s2tht) + SiB*stht2*(-3.*mZ2*cW2
     & + 6.*mt2/SiB2 + mZ2*sW2)
     & + ctht2*(6.*mt2/SiB - 4.*mZ2*SiB*sW2)))/mt2
      Ht22Hmx(3) = (lam*Sqrt2*vev*ctht/TaB*stht)/mt

      gth11 = Hmx(Sind,1)*Ht11Hmx(1)
     & + Hmx(Sind,2)*Ht11Hmx(2) + Hmx(Sind,3)*Ht11Hmx(3)
      gth12 = Hmx(Sind,1)*Ht12Hmx(1)
     & + Hmx(Sind,2)*Ht12Hmx(2) + Hmx(Sind,3)*Ht12Hmx(3)
      gth21 = Hmx(Sind,1)*Ht21Hmx(1)
     & + Hmx(Sind,2)*Ht21Hmx(2) + Hmx(Sind,3)*Ht21Hmx(3)
      gth22 = Hmx(Sind,1)*Ht22Hmx(1)
     & + Hmx(Sind,2)*Ht22Hmx(2) + Hmx(Sind,3)*Ht22Hmx(3)

      gbh11 = Hmx(Sind,1)*Hb11Hmx(1)
     & + Hmx(Sind,2)*Hb11Hmx(2) + Hmx(Sind,3)*Hb11Hmx(3)
      gbh12 = Hmx(Sind,1)*Hb12Hmx(1)
     & + Hmx(Sind,2)*Hb12Hmx(2) + Hmx(Sind,3)*Hb12Hmx(3)
      gbh21 = Hmx(Sind,1)*Hb21Hmx(1)
     & + Hmx(Sind,2)*Hb21Hmx(2) + Hmx(Sind,3)*Hb21Hmx(3)
      gbh22 = Hmx(Sind,1)*Hb22Hmx(1)
     & + Hmx(Sind,2)*Hb22Hmx(2) + Hmx(Sind,3)*Hb22Hmx(3)

      else

      Ab11Amx(1) = 0.d0
      Ab11Amx(2) = 0.d0
      Ab11Amx(3) = 0.d0

      Ab12Amx(1) = (1/2.d0/CoB*((msb12 - msb22)*s2thb 
     &  + lam*mbsb*Sqrt2*vevS*TaB))/mbsb2
      Ab12Amx(2) = (1/Sqrt2*lam*vevS*1/CoB)/mbsb
      Ab12Amx(3) = (1/Sqrt2*lam*vev*TaB)/mbsb

      Ab21Amx(1) = (-1/2.d0/CoB*((msb12 - msb22)*s2thb
     & + lam*mbsb*Sqrt2*vevS*TaB))/mbsb2
      Ab21Amx(2) = (-1/Sqrt2*lam*vevS*1/CoB)/mbsb
      Ab21Amx(3) = (-1/Sqrt2*lam*vev*TaB)/mbsb

      Ab22Amx(1) = 0.d0
      Ab22Amx(2) = 0.d0
      Ab22Amx(3) = 0.d0

      At11Amx(1) = 0.d0
      At11Amx(2) = 0.d0
      At11Amx(3) = 0.d0

      At12Amx(1) = (1/Sqrt2*lam*vevS/SiB)/mt
      At12Amx(2) = (1/2.d0/SiB*(lam*mt*Sqrt2*vevS/TaB 
     &  + (mst12 - mst22)*s2tht))/mt2
      At12Amx(3) = (1/Sqrt2*lam*vev/TaB)/mt

      At21Amx(1) = (-1/Sqrt2*lam*vevS/SiB)/mt
      At21Amx(2) = (-1/2.d0/SiB*(lam*mt*Sqrt2*vevS/TaB
     &  + (mst12 - mst22)*s2tht))/mt2
      At21Amx(3) = (-1/Sqrt2*lam*vev/TaB)/mt

      At22Amx(1) = 0.d0
      At22Amx(2) = 0.d0
      At22Amx(3) = 0.d0

!      gth11 = Amx(Sind,1)*At11Amx(1)
!     & + Amx(Sind,2)*At11Amx(2) + Amx(Sind,3)*At11Amx(3)
!      gth12 = Amx(Sind,1)*At12Amx(1)
!     & + Amx(Sind,2)*At12Amx(2) + Amx(Sind,3)*At12Amx(3)
!      gth21 = Amx(Sind,1)*At21Amx(1)
!     & + Amx(Sind,2)*At21Amx(2) + Amx(Sind,3)*At21Amx(3)
!      gth22 = Amx(Sind,1)*At22Amx(1)
!     & + Amx(Sind,2)*At22Amx(2) + Amx(Sind,3)*At22Amx(3)

!      gbh11 = Amx(Sind,1)*Ab11Amx(1)
!     & + Amx(Sind,2)*Ab11Amx(2) + Amx(Sind,3)*Ab11Amx(3)
!      gbh12 = Amx(Sind,1)*Ab12Amx(1)
!     & + Amx(Sind,2)*Ab12Amx(2) + Amx(Sind,3)*Ab12Amx(3)
!      gbh21 = Amx(Sind,1)*Ab21Amx(1)
!     & + Amx(Sind,2)*Ab21Amx(2) + Amx(Sind,3)*Ab21Amx(3)
!      gbh22 = Amx(Sind,1)*Ab22Amx(1)
!     & + Amx(Sind,2)*Ab22Amx(2) + Amx(Sind,3)*Ab22Amx(3)

      !Amx starts from the basis (G,A,S), thus prerotation
      !from (H_d,H_u) to (G,A) with Beta
      gth11 = Amx(Sind,1)*(CoB*At11Amx(1) - SiB*At11Amx(2))
     & + Amx(Sind,2)*(SiB*At11Amx(1) + CoB*At11Amx(2))
     & + Amx(Sind,3)*At11Amx(3)
      gth12 = Amx(Sind,1)*(CoB*At12Amx(1) - SiB*At12Amx(2))
     & + Amx(Sind,2)*(SiB*At12Amx(1) + CoB*At12Amx(2))
     & + Amx(Sind,3)*At12Amx(3)
      gth21 = Amx(Sind,1)*(CoB*At21Amx(1) - SiB*At21Amx(2))
     & + Amx(Sind,2)*(SiB*At21Amx(1) + CoB*At21Amx(2))
     & + Amx(Sind,3)*At21Amx(3)
      gth22 = Amx(Sind,1)*(CoB*At22Amx(1) - SiB*At22Amx(2))
     & + Amx(Sind,2)*(SiB*At22Amx(1) + CoB*At22Amx(2))
     & + Amx(Sind,3)*At22Amx(3)
     
      gbh11 = Amx(Sind,1)*(CoB*Ab11Amx(1) - SiB*Ab11Amx(2))
     & + Amx(Sind,2)*(SiB*Ab11Amx(1) + CoB*Ab11Amx(2))
     & + Amx(Sind,3)*Ab11Amx(3)
      gbh12 = Amx(Sind,1)*(CoB*Ab12Amx(1) - SiB*Ab12Amx(2))
     & + Amx(Sind,2)*(SiB*Ab12Amx(1) + CoB*Ab12Amx(2))
     & + Amx(Sind,3)*Ab12Amx(3)
      gbh21 = Amx(Sind,1)*(CoB*Ab21Amx(1) - SiB*Ab21Amx(2))
     & + Amx(Sind,2)*(SiB*Ab21Amx(1) + CoB*Ab21Amx(2))
     & + Amx(Sind,3)*Ab21Amx(3)
      gbh22 = Amx(Sind,1)*(CoB*Ab22Amx(1) - SiB*Ab22Amx(2))
     & + Amx(Sind,2)*(SiB*Ab22Amx(1) + CoB*Ab22Amx(2))
     & + Amx(Sind,3)*Ab22Amx(3)

      end if

      end

C-}}}
C-{{{ header:

c.. evalcsusy.f  --  Version 2.00
c.. ------------------------------------------------------------
c.. Program: evalcsusy
c.. Authors: R. Harlander (U Karlsruhe) and M. Steinhauser (U Hamburg)
c.. Date   : August 2004
c.. v2.00  : June 2005
c.. ------------------------------------------------------------
c.. 
c.. for documentation, see
c.. 
c.. ***********
c.. "Supersymmetric Higgs production in gluon fusion at
c..  next-to-leading order"
c..  by Robert V. Harlander and Matthias Steinhauser
c..  JHEP 09 (2004) 066 [hep-ph/0409010].
c.. 
c..  and
c.. 
c.. "Pseudo-scalar Higgs production at next-to-leading order SUSY-QCD"
c..  by Robert V. Harlander and Franziska Hofmann
c..  TTP05-07, SFB/CPP-05-17, hep-ph/0507041
c.. ***********
c.. 
c.. compilation:
c.. > make
c.. 
c.. or:
c.. > f77 -c functions.f readslha.f slhablocks.f \
c..   odeint8.f rkck8.f rkqs8.f.f obj/*.f
c.. > mv gghmix*.o obj/
c.. > f77 -o evalcsusy evalcsusy.f functions.o readslha.o slhablocks.o \
c..   odeint8.o rkck8.o rkqs8.o \
c..   obj/*.o -L`cernlib -v pro kernlib,mathlib`
c.. 
c.. execution:
c.. 
c.. > ./x.evalcsusy example.in test.out
c.. 
c.. subsequent check:
c.. > diff test.out example.out
c.. 
c.. 
c.. **************************
c.. 
c.. Remarks about compilers:
c.. * The current version has been tested on several platforms
c..   using the g77 compiler (version 2.95 or higher).
c..   Problems due to long continuation lines may arise with other compilers.
c.. * Compile time ranges from a few minutes to a few hours, depending
c..   on the platform.
c.. 
c.. **************************
c.. 
c.. CHANGES:
c.. --------
c.. * v2 includes pseudo-scalar Higgs production
c.. 
c.. **************************

C-}}}
C-{{{ subroutine evalcsusy 

      subroutine evalcsusy(Mbot,Msbot1,Msbot2,Mtop,Mstop1,Mstop2,
     &Mgluino,beta,Mw,Mz,alpha,muR,muSUSY,cthetat,sthetat,
     &cthetab,sthetab,c1sm1,c1sm2,c1sm3,c1susy1t,c1susy2t,
     &c1susy1b,c1susy2b,gbh11,gbh12,gbh21,gbh22,
     &gth11,gth12,gth21,gth22,lpseudo,Sind)
      implicit real*8 (a-h,k-z)
      implicit integer(i)
      integer lpseudo, Sind
      data pi/3.141592653589793238462643d0/

      !2013.02.01: quark-Higgs couplings are input

      tanb = dtan(beta)

      salpha = dsin(alpha)
      calpha = dcos(alpha)
      sbeta = dsin(beta)
      cbeta = dcos(beta)

      salphambeta = dsin(alpha-beta)
      calphambeta = dcos(alpha-beta)

      s2thetat = 2*sthetat*cthetat
      c2thetat = cthetat**2 - sthetat**2

      s2thetab = 2*sthetab*cthetab
      c2thetab = cthetab**2 - sthetab**2

      ctW = Mw/Mz
      stW = dsqrt(1.d0-ctW**2)

      c1ewt = 0.d0
      c2ewt = 0.d0
      c1ewb = 0.d0
      c2ewb = 0.d0

      if (lpseudo.eq.0) then

       if (Sind.eq.1) then

         c1EWt = -2.d0*Mz**2/Mtop**2 * (1.d0/2.d0-2.d0/3.d0*stW**2) *
     .     dsin( alpha + beta )
         c2EWt = -2.d0*Mz**2/Mtop**2 * 2.d0/3.d0*stW**2 *
     .     dsin( alpha + beta )
         c1EWmc2EWt = Mtop**2*(c1EWt - c2EWt)
         c1EWpc2EWt = Mtop**2*(c1EWt + c2EWt)

         c1EWb =  2.d0*Mz**2/Mbot**2 * (1.d0/2.d0-1.d0/3.d0*stW**2) *
     .     dsin( alpha + beta )
         c2EWb =  2.d0*Mz**2/Mbot**2 * 1.d0/3.d0*stW**2 *
     .     dsin( alpha + beta )
         c1EWmc2EWb = Mbot**2*(c1EWb - c2EWb)
         c1EWpc2EWb = Mbot**2*(c1EWb + c2EWb)

       else if (Sind.eq.2) then

         c1EWt = 2.d0*Mz**2/Mtop**2 * (1.d0/2.d0-2.d0/3.d0*stW**2) *
     .     dcos( alpha + beta )
         c2EWt = 2.d0*Mz**2/Mtop**2 * 2.d0/3.d0*stW**2 *
     .     dcos( alpha + beta )
         c1EWmc2EWt = Mtop**2*(c1EWt - c2EWt)
         c1EWpc2EWt = Mtop**2*(c1EWt + c2EWt)

         c1EWb = - 2.d0*Mz**2/Mbot**2 * (1.d0/2.d0-1.d0/3.d0*stW**2) *
     .     dcos( alpha + beta )
         c2EWb = - 2.d0*Mz**2/Mbot**2 * 1.d0/3.d0*stW**2 *
     .     dcos( alpha + beta )
         c1EWmc2EWb = Mbot**2*(c1EWb - c2EWb)
         c1EWpc2EWb = Mbot**2*(c1EWb + c2EWb)

       end if

      else if (lpseudo.eq.1) then

         c1EWt = 0.d0
         c2EWt = 0.d0
         c1EWmc2EWt = 0.d0
         c1EWpc2EWt = 0.d0

         c1EWb =  0.d0
         c2EWb =  0.d0
         c1EWmc2EWb = 0.d0
         c1EWpc2EWb = 0.d0

      end if

      mfact = (Mstop1**2 - Mstop2**2)/2.d0/Mtop**2
      mfacb = (Msbot1**2 - Msbot2**2)/2.d0/Mbot**2

c..   top, stop and bottom, sbottom couplings:
      if (lpseudo.eq.0) then

       if (Sind.eq.1) then
         !light Higgs
         !gth = calpha/sbeta
         !gbh = - salpha/cbeta

         gth11ew = c1EWt*cthetat**2 + c2EWt*sthetat**2
         gth22ew = c1EWt*sthetat**2 + c2EWt*cthetat**2
         gth12ew = 1/2.d0*(c2EWt - c1EWt)*s2thetat
         gth21ew = gth12ew
         gth11mu = muSUSY/Mtop*calphambeta/sbeta**2*s2thetat
         gth22mu = -gth11mu
         gth12mu = muSUSY/Mtop*calphambeta/sbeta**2*c2thetat
         gth21mu = gth12mu
         gth11al = calpha/sbeta*( 2.d0 + mfact*s2thetat**2 )
         gth22al = calpha/sbeta*( 2.d0 - mfact*s2thetat**2 )
         gth12al = mfact*calpha/sbeta*s2thetat*c2thetat
         gth21al = gth12al

         gth11 = gth11ew + gth11mu + gth11al
         gth22 = gth22ew + gth22mu + gth22al
         gth12 = gth12ew + gth12mu + gth12al
         gth21 = gth21ew + gth21mu + gth21al

         gbh11ew = c1EWb*cthetab**2 + c2EWb*sthetab**2
         gbh22ew = c1EWb*sthetab**2 + c2EWb*cthetab**2
         gbh12ew = 1/2.d0*(c2EWb - c1EWb)*s2thetab
         gbh21ew = gbh12ew
         gbh11mu = -muSUSY/Mbot*calphambeta/cbeta**2*s2thetab
         gbh22mu = -gbh11mu
         gbh12mu = -muSUSY/Mbot*calphambeta/cbeta**2*c2thetab
         gbh21mu = gbh12mu
         gbh11al = -salpha/cbeta*( 2.d0 + mfacb*s2thetab**2 )
         gbh22al = -salpha/cbeta*( 2.d0 - mfacb*s2thetab**2 )
         gbh12al = -mfacb*salpha/cbeta*s2thetab*c2thetab
         gbh21al = gbh12al

         gbh11 = gbh11ew + gbh11mu + gbh11al
         gbh22 = gbh22ew + gbh22mu + gbh22al
         gbh12 = gbh12ew + gbh12mu + gbh12al
         gbh21 = gbh21ew + gbh21mu + gbh21al

        else if (Sind.eq.2) then
         !heavy Higgs
         !changes compared to light Higgs: cos -> sin, sin -> -cos
         !gth = salpha/sbeta
         !gbh = calpha/cbeta

         gth11ew = c1EWt*cthetat**2 + c2EWt*sthetat**2
         gth22ew = c1EWt*sthetat**2 + c2EWt*cthetat**2
         gth12ew = 1/2.d0*(c2EWt - c1EWt)*s2thetat
         gth21ew = gth12ew
         gth11mu = muSUSY/Mtop*salphambeta/sbeta**2*s2thetat
         gth22mu = -gth11mu
         gth12mu = muSUSY/Mtop*salphambeta/sbeta**2*c2thetat
         gth21mu = gth12mu
         gth11al = salpha/sbeta*( 2.d0 + mfact*s2thetat**2 )
         gth22al = salpha/sbeta*( 2.d0 - mfact*s2thetat**2 )
         gth12al = mfact*salpha/sbeta*s2thetat*c2thetat
         gth21al = gth12al

         gth11 = gth11ew + gth11mu + gth11al
         gth22 = gth22ew + gth22mu + gth22al
         gth12 = gth12ew + gth12mu + gth12al
         gth21 = gth21ew + gth21mu + gth21al

         gbh11ew = c1EWb*cthetab**2 + c2EWb*sthetab**2
         gbh22ew = c1EWb*sthetab**2 + c2EWb*cthetab**2
         gbh12ew = 1/2.d0*(c2EWb - c1EWb)*s2thetab
         gbh21ew = gbh12ew
         gbh11mu = -muSUSY/Mbot*salphambeta/cbeta**2*s2thetab
         gbh22mu = -gbh11mu
         gbh12mu = -muSUSY/Mbot*salphambeta/cbeta**2*c2thetab
         gbh21mu = gbh12mu
         gbh11al = calpha/cbeta*( 2.d0 + mfacb*s2thetab**2 )
         gbh22al = calpha/cbeta*( 2.d0 - mfacb*s2thetab**2 )
         gbh12al = mfacb*calpha/cbeta*s2thetab*c2thetab
         gbh21al = gbh12al

         gbh11 = gbh11ew + gbh11mu + gbh11al
         gbh22 = gbh22ew + gbh22mu + gbh22al
         gbh12 = gbh12ew + gbh12mu + gbh12al
         gbh21 = gbh21ew + gbh21mu + gbh21al

       end if

      else if (lpseudo.eq.1) then
         !pseudoscalar Higgs
         !gth = 1.d0/tanb
         !gbh = tanb

         gth12m0 = mfact*s2thetat/tanb
         gth12m1 = muSUSY/Mtop*(1.d0+1.d0/tanb**2)

         gth11 = 0.d0
         gth22 = 0.d0
         gth12 = gth12m0 + gth12m1
         gth21 = - gth12

         gbh12m0 = mfacb*s2thetab*tanb
         gbh12m1 = muSUSY/Mbot*(1.d0+tanb**2)

         gbh11 = 0.d0
         gbh22 = 0.d0
         gbh12 = gbh12m0 + gbh12m1
         gbh21 = - gbh12

      endif

C-{{{ csusy:

      rnf = 5.d0

      if (lpseudo.eq.0) then
       if (Sind.eq.1) then
c..   scalar coefficient functions:
         c1susy1t = 3.d0 * gghmixlo(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &        cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EWt,c1EWpc2EWt)

         c1susy1b = 3.d0 * gghmixlo(Mbot,Msbot1,Msbot2,Mgluino,muR,
     &        cthetab,sthetab,alpha+Pi/2.d0,beta+Pi/2.d0,
     &       -muSUSY,c1EWmc2EWb,c1EWpc2EWb)

         c1sm1 = 1.d0

         c1susy2t = 3.d0 * gghmixnlosusy(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &        cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EWt,c1EWpc2EWt)

         c1susy2b = 3.d0 * gghmixnlosusy(Mbot,Msbot1,Msbot2,Mgluino,
     &        muR,cthetab,sthetab,alpha+Pi/2.d0,beta+Pi/2.d0,
     &        muSUSY,c1EWmc2EWb,c1EWpc2EWb)

         c1sm2 = 11.d0/4.d0
         c1sm3 = gghqcdnnlo(Mtop,muR,rnf)

       else if (Sind.eq.2) then
c..   heavy scalar coefficient functions:
         c1susy1t = 3.d0 * gghmixlo(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &        cthetat,sthetat,alpha-Pi/2.d0,beta,
     &        muSUSY,c1EWmc2EWt,c1EWpc2EWt)

         c1susy1b = 3.d0 * gghmixlo(Mbot,Msbot1,Msbot2,Mgluino,muR,
     &        cthetab,sthetab,alpha,beta+Pi/2.d0,
     &       -muSUSY,c1EWmc2EWb,c1EWpc2EWb)

         c1sm1 = 1.d0

         c1susy2t = 3.d0 * gghmixnlosusy(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &        cthetat,sthetat,alpha-Pi/2.d0,beta,
     &        muSUSY,c1EWmc2EWt,c1EWpc2EWt)

         c1susy2b = 3.d0 * gghmixnlosusy(Mbot,Msbot1,Msbot2,Mgluino,
     &        muR,cthetab,sthetab,alpha,beta+Pi/2.d0,
     &        muSUSY,c1EWmc2EWb,c1EWpc2EWb)

         c1sm2 = 11.d0/4.d0
         c1sm3 = gghqcdnnlo(Mtop,muR,rnf)

       end if

      else if (lpseudo.eq.1) then
c..   pseudo-scalar coefficient functions:
         c1sm1 = 1.d0
         c1sm2 = 0.d0
         c1sm3 = 0.d0
         c2sm1 = -1/2.d0 + 2*dlog(muR/Mtop)

         c1susy1t = 1.d0/tanb
         c1susy2t = ggamixnlosusy(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &        cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EWt,c1EWpc2EWt)

         c1susy1b = tanb
!         c1susy2b = ggamixnlosusy(Mbot,Msbot1,Msbot2,Mgluino,muR,
!     &        sthetab,alpha+Pi/2.d0,beta+Pi/2.d0,
!     &        muSUSY,c1EWmc2EWb,c1EWpc2EWb)
         c1susy2b = 0.d0

      endif

C-}}}

      end

C-}}}
C-{{{ function b0fin:

      real*8 function b0fin(q2,m1,m2,mu)
c..
c..   This is the finite part of the B0 function.
c..   Only the real part is implemented; the imaginary part
c..   is set to zero.
c..   Note: not all cases are covered yet; corresponding 
c..   'if-endif' still missing.
c..   
c..   Is it necessary to implement also special cases (mi=0, q2=0,  ...)?
c..   
c..   q2 == q^2
c..   m1 == m_1
c..   m2 == m_2
c..

      implicit real*8 (a-z)

      mu2 = mu*mu

      tmp = m1
      m1 = dabs(m1)

      if ( (q2.le.(m1-m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        + dsqrt((m1+m2)**2-q2)*dsqrt((m1-m2)**2-q2)/q2 *
     .        dlog(  
     .        (dsqrt((m1+m2)**2-q2)+dsqrt((m1-m2)**2-q2)) /
     .        (dsqrt((m1+m2)**2-q2)-dsqrt((m1-m2)**2-q2)) 
     .        ) )
      elseif ( (q2.ge.(m1+m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - dsqrt(q2-(m1+m2)**2)*dsqrt(q2-(m1-m2)**2)/q2 *
     .        dlog(-  
     .        (dsqrt(q2-(m1+m2)**2)+dsqrt(q2-(m1-m2)**2)) /
     .        (dsqrt(q2-(m1+m2)**2)-dsqrt(q2-(m1-m2)**2)) 
     .        ) )
      elseif ( (q2.lt.(m1+m2)**2).and.(q2.gt.(m1-m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - 2.d0/q2 *
     .        dsqrt( dabs(q2-(m1+m2)**2) )*
     .        dsqrt( dabs(q2-(m1-m2)**2) )*
     .        datan( dsqrt(dabs(q2-(m1-m2)**2)) / 
     .        dsqrt(dabs(q2-(m1+m2)**2)) )
     .        )
      else
         write(6,*) '<function b0fin>: ',
     .        'arguments (',q2,',',m1,',',m2,') not implemented.'
         stop
      endif
      m1 = tmp
      end

C-}}}
C-{{{ function b0finim:

      function b0finim(q2,m1,m2,mu)
c..
c..   This is the finite part of the B0 function.
c..   Note: not all cases are covered yet; corresponding 
c..   'if-endif' still missing.
c..   
c..   Is it necessary to implement also special cases (mi=0, q2=0,  ...)?
c..   
c..   q2 == q^2
c..   m1 == m_1
c..   m2 == m_2
c..

      implicit real*8 (c-z)
      double complex b0finim
      parameter (pi=3.1415926535897932385d0)

      mu2 = mu*mu

      if ( (q2.le.(m1-m2)**2) ) then
         b0finim = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        + dsqrt((m1+m2)**2-q2)*dsqrt((m1-m2)**2-q2)/q2 *
     .        dlog(  
     .        (dsqrt((m1+m2)**2-q2)+dsqrt((m1-m2)**2-q2)) /
     .        (dsqrt((m1+m2)**2-q2)-dsqrt((m1-m2)**2-q2)) 
     .        ) )
      elseif ( (q2.ge.(m1+m2)**2) ) then
         b0finim = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - dsqrt(q2-(m1+m2)**2)*dsqrt(q2-(m1-m2)**2)/q2 *
     .        dlog(-  
     .        (dsqrt(q2-(m1+m2)**2)+dsqrt(q2-(m1-m2)**2)) /
     .        (dsqrt(q2-(m1+m2)**2)-dsqrt(q2-(m1-m2)**2)) 
     .        ) ) 
     .    + (0.d0,1.d0)*Pi*dsqrt(q2-(m1+m2)**2)*dsqrt(q2-(m1-m2)**2)/q2
      elseif ( (q2.lt.(m1+m2)**2).and.(q2.gt.(m1-m2)**2) ) then
         b0finim = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - 2.d0/q2 *
     .        dsqrt( dabs(q2-(m1+m2)**2) )*
     .        dsqrt( dabs(q2-(m1-m2)**2) )*
     .        datan( dsqrt(dabs(q2-(m1-m2)**2)) / 
     .        dsqrt(dabs(q2-(m1+m2)**2)) )
     .        )
      else
         write(6,*) '<function b0fin>: ',
     .        'arguments (',q2,',',m1,',',m2,') not implemented.'
         stop
      endif
      end

C-}}}
C-{{{ function gghmixlo:

      real*8 function gghmixlo(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
      
      implicit real*8 (a-z)
c..   
c..   LO ggH coefficient function in the MSSM
c..   In the notation of the paper gghmixlo corresponds to
c..   1/3*( c^{(0)}_{1,t} + c^{(0)}_{1,\tilde t} )
c..   Note: C_1 = -\alpha_s/\pi * gghmixlo + h.o.
c..
      sbeta = dsin(beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      calphambeta = dcos(alpha-beta)
      !cthetat = dsqrt(1.d0-sthetat*sthetat)

      dum1 = calpha*(1/(3.*sbeta) + Mtop**2/(12.*Mstop1**2*sbeta) + Mtop
     &     **2/(12. *Mstop2**2*sbeta) + (cthetat**2*sthetat**2)/(6.
     &     *sbeta) - (cthetat**2*Mstop1**2*sthetat**2)/(12.*Mstop2**2
     &     *sbeta) - (cthetat**2*Mstop2**2*sthetat**2)/(12.*Mstop1**2
     &     *sbeta))

      dum2 = muSUSY*( (calphambeta*cthetat*Mtop*sthetat)/(12.*Mstop1**2
     &     *sbeta**2)-(calphambeta*cthetat*Mtop*sthetat)/(12.*Mstop2**2
     &     *sbeta**2) )
      
      dum3 = c1EWpc2EW*( 1/(48.*Mstop1**2) + 1/(48.*Mstop2**2) )

      dum4 = c1EWmc2EW*(
     &     cthetat**2/(48.*Mstop1**2) - cthetat**2/(48.*Mstop2**2) - 
     &     sthetat**2/(48.*Mstop1**2) + sthetat**2/(48.*Mstop2**2) )

      gghmixlo = dum1 + dum2 + dum3 + dum4

      end

C-}}}
C-{{{ function gghmixnlosusy:

      real*8 function gghmixnlosusy(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)

c..
c..

      implicit real*8 (a-z)

c..   allowed difference between mass-differences:
      del = .01d0

      Mtoptmp1 = Mtop*(1.d0+15.d0*del)
      Mtoptmp2 = Mtop*(1.d0-15.d0*del)
      Mgluinotmp1 = Mgluino*(1.d0+15.d0*del)
      Mgluinotmp2 = Mgluino*(1.d0-15.d0*del)
      Mstop1tmp1 = Mstop1*(1.d0+del)
      Mstop1tmp2 = Mstop1*(1.d0-del)

      Mstop1def = Mstop1
      if ( dabs(Mstop1-Mstop2)/Mstop1.le.del ) then 
         Mstop1=Mstop2*(1.+del/10.d0)
      endif

      dumbare = 0.d0
      dumCT = 0.d0

      if ( ( (dabs(Mstop1-(Mgluino-Mtop))/Mstop1).le.del ).or.
     &     ( (dabs(Mstop2-(Mgluino-Mtop))/Mstop2).le.del ).or.
     &     ( (dabs(Mstop1-(Mgluino+Mtop))/Mstop1).le.del ).or.
     &     ( (dabs(Mstop2-(Mgluino+Mtop))/Mstop2).le.del ) ) then
         dumbare1= (gghnlobare(Mtoptmp1,Mstop1,Mstop2,Mgluino,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtoptmp1,Mstop1,Mstop2,Mgluino,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
         dumbare2= (gghnlobare(Mtoptmp2,Mstop1,Mstop2,Mgluino,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtoptmp2,Mstop1,Mstop2,Mgluino,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
 
         dumbare = (dumbare1+dumbare2)/2.d0

       !high SUSY masses case
       if ( ( (dabs(Mstop1-Mgluino)/Mstop1).le.del ).or.
     &     ( (dabs(Mstop2-Mgluino)/Mstop2).le.del ).or.
     &     ( (dabs(Mstop1-Mgluino)/Mstop1).le.del ).or.
     &     ( (dabs(Mstop2-Mgluino)/Mstop2).le.del ) ) then
         dumbare1= (gghnlobare(Mtop,Mstop1,Mstop2,Mgluinotmp1,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtop,Mstop1,Mstop2,Mgluinotmp1,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
         dumbare2= (gghnlobare(Mtop,Mstop1,Mstop2,Mgluinotmp2,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtop,Mstop1,Mstop2,Mgluinotmp2,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
          dumbare = (dumbare1+dumbare2)/2.d0
       end if

      elseif ( dabs(Mstop1-Mstop2)/Mstop1.le.del ) then 
         dumbare1= (gghnlobare(Mtop,Mstop1tmp1,Mstop2,Mgluino,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtop,Mstop1tmp1,Mstop2,Mgluino,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
         dumbare2= (gghnlobare(Mtop,Mstop1tmp2,Mstop2,Mgluino,muR,
     &        Mstop1,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtop,Mstop1tmp2,Mstop2,Mgluino,muR,
     &        Mstop2,cthetat,sthetat,alpha,beta,
     &        muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0

         dumbare = (dumbare1+dumbare2)/2.d0
      else
         dumbare = (gghnlobare(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &              gghnlobare(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0
      endif

      dumCT =  (gghmix_CT1(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT1(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_CT2(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT2(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_CT3(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT3(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_CT4(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT4(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_CT5(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT5(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_CT6(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop1,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     gghmix_CT6(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     Mstop2,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0

      gghmixnlosusy = dumbare + dumCT

      Mstop1 = Mstop1def

      end

C-}}}
C-{{{ function ggamixnlosusy:

      real*8 function ggamixnlosusy(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
c..
c..
      implicit real*8 (a-z)
      real*8 masses(4),mmax(4)
      integer imax(4),i,j,k,npoles,tmp

c..   allowed difference between masses:
      del = 1d-4
      delfac = 0.d0

      if ((Mgluino.lt.1.d0).or.(Mtop.lt.1.d0).or.(Mstop1.lt.1.d0).or
     &     .(Mstop2.lt.1.d0))  then
         write(6,*) 'evalcsusy.f (fatal): ',
     &        'One of the masses is too small.',' Stopped.'
c         call printdie('One of the masses is too small.')
         stop
      endif

      masses(1) = Mstop1
      masses(2) = Mstop2
      masses(3) = Mtop
      masses(4) = Mgluino

      do i=1,4
         mmax(i) = masses(i)
         imax(i) = i
      enddo

      do i=1,4
         do j=i+1,4
            if (mmax(j).gt.mmax(i)) then
               tmp = mmax(j)
               mmax(j) = mmax(i)
               mmax(i) = tmp
               tmp = imax(i)
               imax(i) = imax(j)
               imax(j) = tmp
            endif
         enddo
      enddo
      
 101  npoles = 0
      do k=1,4
         do i=k+1,4
            do j=i+1,4
C               print*,imax(k),imax(i),imax(j)
               if ( dabs((masses(imax(i))+masses(imax(j)))/masses(imax(k
     &              )) - 1.d0).lt.del )then
                  masses(imax(k)) = (masses(imax(i))+masses(imax(j)))
     &                 * (1 + delfac*del)
C                  print*,masses(imax(i)),masses(imax(j)),masses(imax(k))
C     &                 ,mmax(k)
                  npoles = npoles+1
               endif
            enddo
         enddo
      enddo
      if (npoles.gt.0) then
         delfac = delfac + 1.d0
c         print*,'delfac = ',delfac,npoles
         goto 101
      endif

      if (dabs(masses(1)/masses(2) - 1.d0).lt.del) then
         masses(1) = masses(2) * ( 1 + 2*del )
      endif


      Mstop1l = masses(1)
      Mstop2l = masses(2)
      Mtopl = masses(3)
      Mgluinol = masses(4)

      ggamixnlosusy = (ggamix_bare(Mtopl,Mstop1l,Mstop2l,Mgluinol,muR,
     &     Mstop1l,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW) +
     &     ggamix_bare(Mtopl,Mstop1l,Mstop2l,Mgluinol,muR,
     &     Mstop2l,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW))/2.d0

      end

C-}}}
C-{{{ function gghnlobare:

      real*8 function gghnlobare(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)

      implicit real*8 (a-z)
c..   
c..
c..

         gghnlobare = gghmix_bare1(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_bare2(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_bare3(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_bare4(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_bare5(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)
     &     + gghmix_bare6(Mtop,Mstop1,Mstop2,Mgluino,muR,
     &     q0,cthetat,sthetat,alpha,beta,
     &     muSUSY,c1EWmc2EW,c1EWpc2EW)

      end

C-}}}
C-{{{ function gghqcdnnlo:

      real*8 function gghqcdnnlo(Mtop,muR,nl)

      implicit real*8 (a-z)
c..   
c..   NNLO (pure QCD) coefficient of C_1, i.e. the coefficient of
c..   api5^2 of the expression (we have nl=5):
c..
c..   1 + (11*api5)/4 + api5^2*(2777/288 + (19*lmm)/16 + (-67/96 + lmm/3)*nl)  
c..

      lmm = dlog(muR**2/Mtop**2)

      gghqcdnnlo = 2777.d0/288.d0 + (19.d0*lmm)/16.d0 
     .     + (-67.d0/96.d0 + lmm/3.d0)*nl

      end

C-}}}
C-{{{ function davytausk:

      real*8 function davytausk(xin,yin,m3in)
c..
c..   Version 2.0
c..
c..   This is the function
c..   
c..   \lambda^2(x,y) * \Phi(x,y)
c..   
c..   as defined in
c..   A.I. Davydychev and J.B. Tausk,  Nucl. Phys. B 397 (1993) 123
c..   Eq. (4.10). 
c..
c..   It uses CERNLIB for the Clausen function, and the dilogarithm from
c..   other sources.
c..   Thus, compile e.g. with
c..
c..   setenv LIBS `cernlib -v pro kernlib,mathlib`
c..   f77 -o davytausk davytausk.f ../ggah/v2.0/functions.f -L$LIBS
c..
c..   Note that the printed version of Eq.(4.10) has a typo:
c..   the arguments
c..   1/2 *  1 + x - y - \lambda
c..   and
c..   1/2 *  1 - x + y - \lambda
c..   should read
c..   1/2 * (1 + x - y - \lambda)
c..   and
c..   1/2 * (1 + x - y - \lambda)
c..
c..   
      implicit real*8 (a-z)
      real*8 m32,m3in
      data pi/3.14159265358979323846264338328d0/
      complex*16 dtphixyc,cli2

      sqx = dsqrt(xin)
      sqy = dsqrt(yin)

      m32 = 0.d0
      xx = 0.d0
      yy = 0.d0

      dtlam2 = (1.d0 - xin - yin)**2 - 4*xin*yin

c..   First the case lambda^2 >= 0, which is given by Eq.(4.10).
c..   But before, we need to do some mapping to cover all
c..   kinematical cases (see discussion below Eq.(4.10)).
      if (dtlam2.gt.0) then
c..   m1 + m2 <= m3:
         if     ( (sqx+sqy).le.1.d0 ) then
            xx = xin
            yy = yin
            m32 = m3in**2
c..   m2 + m3 <= m1  (m3 <-> m1):
         elseif ( (sqx-sqy).gt.1.d0 ) then
            xx = 1/xin
            yy = yin/xin
            m32 = xin*m3in**2
c..   m1 + m3 <= m2  (m2 <-> m3):
         elseif ( (sqy-sqx).gt.1.d0 ) then
            xx = xin/yin
            yy = 1/yin
            m32 = yin*m3in**2
         endif

c..   the sign of lambda^2 is invariant under the mappings above, but
c..   lambda itself is not:
         dtlam2xy = (1.d0 - xx - yy)**2 - 4*xx*yy
         dtlambda = dsqrt(dtlam2xy)

c..   Eq.(4.10) times lambda^2:
         dtphixyc = 
     &        dtlambda * m32 * ( 2*dlog( dtarg(xx,yy,dtlambda) )
     &        *dlog( dtarg(yy,xx,dtlambda) )
     &        - dlog(xx)*dlog(yy)
     &        - 2*cli2( (1.d0,0.d0) * dtarg(xx,yy,dtlambda) )
     &        - 2*cli2( (1.d0,0.d0) * dtarg(yy,xx,dtlambda) ) 
     &        + pi*pi/3.d0 )

c..   Now the case lambda^2 < 0:
      elseif (dtlam2.lt.0) then
c..   and again the mappings, but this time   
c..   m1 + m2 >= m3  is the default:
         if     ( (sqx+sqy).ge.1.d0 ) then
            xx = xin
            yy = yin
            m32 = m3in**2
         elseif ( (sqx-sqy).lt.1.d0 ) then
            xx = 1/xin
            yy = yin/xin
            m32 = xin*m3in**2
         elseif ( (sqy-sqx).lt.1.d0 ) then
            xx = xin/yin
            yy = 1/yin
            m32 = yin*m3in**2
         endif

         dtlam2xy = (1.d0 - xx - yy)**2 - 4*xx*yy
         dtlambda = dsqrt(-dtlam2xy)

c..   Eq.(4.15) times lambda^2:
         dtphixyc = - 2.d0 * m32 * dtlambda *
     &        ( dtclaus(1.d0,xx,yy)
     &        + dtclaus(yy,xx,1.d0)
     &        + dtclaus(xx,yy,1.d0) )

      else
         dtphixyc = 0.d0
      endif
         
      if (dabs(dimag(dtphixyc)).gt.1d-10) stop 1003
      davytausk = dreal(dtphixyc)

      end

C-}}}
C-{{{ function dtclaus:

      real*8 function dtclaus(arg1,arg2,arg3)
c..
c..   This is Clausen's function with the particular argument as 
c..   defined in Eq.(4.15) of
c..   Davydychev, Tausk,  NPB 397 (1993) 123.
c..
c..   dclaus(x) is taken from CERNLIB which has to be compiled in.
c..   
      implicit real*8 (a-z)
      external dclaus

      argu = 2*acos( (- arg1 + arg2 + arg3)/2.d0/dsqrt(arg2*arg3) )
      dtclaus = dclaus(argu)

      end

C-}}}
C-{{{ function dclaus

      function dclaus(x)
      implicit none
      double precision dclaus,x
      double complex arg,cli2
      arg = cdexp((0.d0,1.d0)*x)
      dclaus = Aimag(cli2(arg))
      end

C-}}}
C-{{{ function dtarg:

      real*8 function dtarg(xin,yin,xylam)

      implicit real*8 (a-z)

      dtarg = (1 + xin - yin - xylam)/2.d0

      end

C-}}}
