! This file is part of SusHi.
!
! It includes the call of FeynHiggs.
!

#include "versions.h"

      subroutine runfeynhiggs(alfasmz,mt,mbmb,mw,mz,mh,pseudo
     &,msb12,msb22,mst12,mst22)

      implicit none
      integer nmfv, pseudo
      double precision alfasmz,mt,mbmb,mw,mz,mh
      include '../commons/common-inputoutput.f'
      include '../commons/common-ren.f'
      include '../commons/common-consts.f' !Pi
      include '../commons/common-keys.f'
      include '../commons/common-citations.f'
      !FH specific variables
      double precision MCha(2),MNeu(4),MHtree(4)
      double precision MHiggs(4),SAtree,Mst(2),MTrun,MSb(2),SUMbSL2
      double complex UCha(2,2),VCha(2,2),ZNeu(4,4)
      double complex SAeff,UHiggs(3,3),ZHiggs(3,3),USt(2,2),USb(2,2)
      double complex DeltaMB,Deltab
      double precision msb12,msb22,mst12,mst22
      character fhvers*20,hbvers*20,hsvers*20,thdmcvers*20
      common /vers/ fhvers,hbvers,hsvers,thdmcvers
      integer error
#ifndef FEYNHIGGS2100
      double precision MASf(6,4)
      double complex UASf(6,6,4)
      double precision MSf(2,4,3)
      double complex USf(2,2,4,3)
!extension of mass matrices and vectors in FH 2.10.0
#else
      double precision MASf(6,5)
      double complex UASf(6,6,5)
      double precision MSf(2,5,3)
      double complex USf(2,2,5,3)
#endif

#ifdef HIGGSBOUNDS
#include "HBandHSwithFH.h"
      citations(31) = 1
      citations(32) = 1
      citations(33) = 1
#endif

      fhvers = FHVERSION
      hbvers = HBVERSION
      hsvers = HSVERSION

      error = 0
#ifndef FEYNHIGGS2120
      call FHSetFlags(error,
     &     mssmpart, fieldren, tanbren, higgsmix, p2approx,
     &     looplevel, runningMT, botResum, tlCplxApprox)
#else
      call FHSetFlags(error,
     &     mssmpart, fieldren, tanbren, higgsmix, p2approx,
     &     looplevel, loglevel, runningMT, botResum, tlCplxApprox)
#endif

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHSetFlags ', error
         write(*,105)
         stop
      endif

      !to avoid a wrong setting of alpha_EM in the SusHi to FeynHiggs
      !link FH is asked to use its default value - changed 1.6.0, 27.03.2015
      !Therefore appears a -1 in FHSetSMPara.
      !However note that invAlfaMZ has no effect on the masses calculations in FH!
      call FHSetSMPara(error,
     &     -1, AlfasMZ, GF,
     &     ME, MU, MD, MM, MC, MS, ML, mbmb,
     &     MW, MZ,
     &     CKMlambda, CKMA, CKMrhobar, CKMetabar)

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHSetSMPara ', error
         write(*,105)
         stop
      endif

      call FHSetPara(error,scalefactor,MT,tanb,MA0,MHp,
     &     M3SL,M3SE,M3SQ,M3SU,M3SD,
     &     M2SL,M2SE,M2SQ,M2SU,M2SD,
     &     M1SL,M1SE,M1SQ,M1SU,M1SD,
     &     MUE,cAtau,cAt,cAb,
     &     cAmu,cAc,cAs,cAe,cAu,cAd,M_1,M_2,M_3,Qtau,Qt,Qb)

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHSetPara ', error
         write(*,105)
         stop
      endif

#ifndef FEYNHIGGS295
      call FHGetPara(error, nmfv, MASf, UASf,
     &     MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &     MHtree, SAtree)
#else
      call FHGetPara(error, nmfv, MSf, USf, MASf, UASf,
     &     MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &     MHtree, SAtree)
#endif

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHGetPara ', error
         write(*,105)
         stop
      endif

      call FHHiggsCorr(error,MHiggs,SAeff,UHiggs,ZHiggs)

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHHiggsCorr ', error
         write(*,105)
         stop
      endif

#ifndef FEYNHIGGS2113
      call FHGetTLPara(error,MSb,USb,SUMbSL2,Deltab)
#else
      call FHGetTLPara(error,MSt,USt,MTrun,MSb,USb,SUMbSL2,Deltab)
#endif

      if(error.ne.0) then
         write(*,105)
         write(*,*) 'Error in FHGetTLPara ', error
         write(*,105)
         stop
      endif

#ifdef HIGGSBOUNDS
!Call HiggsBounds and HiggsSignals
#include "HBandHSwithFH.F"
#endif

      dMLb = SUMbSL2-M3SQ**2

!Set delmb equal to DeltaMB from FH - own calculation in renscheme.f
      delmb = Real(DeltaMB)

      msb12 = Msb(1)**2
      msb22 = Msb(2)**2

      cthetab = Real(Usb(1,1))
      sthetab = Real(Usb(1,2))

      !only real MSSM implemented, SAeff expected to be real
      alpha = dasin(Real(SAeff))
      mst12 = MASf(3,3)**2
      mst22 = MASf(6,3)**2
      cthetat = Real(UASf(3,3,3))
      sthetat = Real(UASf(3,6,3))
      muSUSY = Real(MUE)

      Hmx = 0.d0
      Hmx(1,1) = -Real(SAeff)
      Hmx(1,2) = dsqrt(1.d0-Real(SAeff)*Real(SAeff))
      Hmx(2,1) = dsqrt(1.d0-Real(SAeff)*Real(SAeff))
      Hmx(2,2) = Real(SAeff)

      Amx = 0.d0
      Amx(2,2) = 1.d0

!thetab from OS-sbottom masses in renormalize()

      if (pseudo.eq.0) then
       if (Sind.eq.1) then
         Mh = MHiggs(1)         !h
       else if (Sind.eq.2) then
         Mh = MHiggs(2)         !H
       end if
      else if (pseudo.eq.1) then
         Mh = MHiggs(3)         !A
      endif

 105  format('#--------------------------------------------------#')

      end



#ifdef HIGGSBOUNDS
!************************************************************** 
	subroutine get_singleH_uncertainty(dCS,dggh,dbbh,g2hgg,g2hbb,mh)
	use theory_colliderSfunctions

	implicit none	
	double precision dCS
 	double precision dggh, dbbh, g2hgg, g2hbb, mh
 	double precision vsmall

 	vsmall=1.0D-10

	if(g2hgg.le.vsmall.and.g2hbb.le.vsmall) then
		dCS = 0.0D0
	else 
		dCS=(g2hgg*LHC8_rH_gg(mh)*dggh+
     & 		g2hbb*LHC8_rH_bb(mh)*dbbh)/
     &      (g2hgg*LHC8_rH_gg(mh)+g2hbb*LHC8_rH_bb(mh))
	endif
	end
!************************************************************** 
#endif
