! This file is part of SusHi.
!
! This file contains the 2 loop amplitudes from squark-quark-gluino
! diagrams taken from G. Degrassi, P. Slavich and S. Di Vita.
! Implementation by S.Liebler and H.Mantler.
! References are: arXiv:1007.3465, 1107.0914, 1204.1016
! October 2014: Added NMSSM contributions - S. Liebler
! February 2016: Making the routines more modular, i.e.
! no common blocks are used, instead all input parameters
! have to be specified. - S. Liebler
!
#define debug 0
#define DRbartoOSt 1
#define DRbartoOSb 1

C-{{{ B1SUSYDS

      function B1SUSYDS(mt2,mb2,mbsb2,mbmb,mbyuk,mh2,mgl
     & ,beta,muSUSY,GF,lam,delta_mb,mw,mz
     & ,Ab,msb12,msb22,s2thetab,c2thetab
     & ,At,mst12,mst22,s2thetat,c2thetat
     & ,Hmx,Amx
     & ,tanbresum,mbrunyuk,yukfac,muRggh
     & ,dmsb1,dmsb2,dthetab,dmbsb,dAbds
     & ,stopcont,sbotcont,pseudo,Sindin)
      implicit none
      !relevant input:
      double precision stopcont,sbotcont
      integer pseudo,Sindin
      double precision mb2,mh2,msb12,msb22,mgl,muSUSY,Ab,beta
      double precision delta_mb,GF,mt2,mst12,mst22,lam,mbyuk
      double precision At,s2thetat,c2thetat,s2thetab,c2thetab
      double precision yukfac(9),mbsb2,mbmb
      double precision Hmx(3,3), Amx(3,3),mw,mz,muRggh
      integer tanbresum,mbrunyuk
      !input counterterms:
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)
      !internal:
      double complex B1SUSYDS,H1DS,H2DS,H3DS,Kbbg,Kttg
      double complex KbbgS, KttgS
      double precision Kfac,Kfacexp
      double precision Sqrt2, vev, vevS,mb

!      write(*,*) "1 ",stopcont,sbotcont
!      write(*,*) "2 ",pseudo,Sindin
!      write(*,*) "3 ",mb2,mh2,msb12,msb22,mgl,muSUSY,Ab,beta
!      write(*,*) "4 ",delta_mb,GF,mt2,mst12,mst22,lam,mbyuk
!      write(*,*) "5 ",At,s2thetat,c2thetat,s2thetab,c2thetab
!      write(*,*) "6 ",yukfac,mbsb2,mbmb
!      write(*,*) "7 ",Hmx, Amx,mw,mz,muRggh
!      write(*,*) "8 ",tanbresum,mbrunyuk
!      write(*,*) "9 ",dmsb1,dmsb2,dthetab,dmbsb,dAbds

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      mb = dsqrt(mb2)

      !Multiplying all Higgs-bottom couplings within the sbottom
      !corrections with an additional factor K
      !1. In case of the usage of mb^DRbar_MSSM with:
      Kfac = mbyuk/mb
      Kfacexp = Kfac
      !2. In case of the usage of mb^OS with tan(beta)-resummation with:
      if (tanbresum.eq.1) then
         Kfac = 1.d0/(1.d0 + delta_mb)
         Kfacexp = Kfac
      else if (tanbresum.eq.2) then
      !Delta_b approximation:
      if (pseudo.eq.0) then
         Kfac = 1.d0/(1.d0 + delta_mb)
     &    * (1 + Hmx(Sindin,2)*delta_mb/(Hmx(Sindin,1)*dtan(beta))
     &   + Hmx(Sindin,3)*vev*dcos(beta)/(Hmx(Sindin,1)*vevS)*delta_mb)
         Kfacexp = Kfac - Hmx(Sindin,2)/(Hmx(Sindin,1)*dtan(beta))
     &   - Hmx(Sindin,3)*vev*dcos(beta)/(Hmx(Sindin,1)*vevS)
      else if (pseudo.eq.1) then
         Kfac = 1.d0/(1.d0 + delta_mb)
     &    * (1 - delta_mb/((dtan(beta))**2) 
     &    - delta_mb*Amx(Sindin,3)*vev/(Amx(Sindin,2)*vevS*dtan(beta)))
         Kfacexp = Kfac + 1.d0/((dtan(beta))**2)
     &    + Amx(Sindin,3)*vev/(Amx(Sindin,2)*vevS*dtan(beta))
      endif      
      endif

      if(pseudo.eq.0) then
         B1SUSYDS = -3/2.d0 * (
     &     Hmx(Sindin,1)
     &      * H1DS(muRggh,mt2,mb2,mbmb,mbsb2,mh2,mgl,msb12,msb22
     &       ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &       ,lam,GF,muSUSY,beta,mw,mz,At,Ab
     &       ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac
     &       ,mbrunyuk,tanbresum,stopcont,sbotcont,Kfac,Kfacexp)
     &    + Hmx(Sindin,2)
     &      * H2DS(muRggh,mt2,mb2,mbsb2,mh2,mgl,msb12,msb22
     &       ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &       ,lam,GF,muSUSY,beta,mw,mz,At
     &       ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac
     &       ,stopcont,sbotcont,Kfac,Kfacexp)
     &    + Hmx(Sindin,3)
     &      * H3DS(muRggh,mt2,mb2,mbsb2,mh2,mgl,msb12,msb22
     &       ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &       ,lam,GF,muSUSY,beta
     &       ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac,stopcont,sbotcont))
      else if(pseudo.eq.1) then
         B1SUSYDS = -3/2.d0 * Amx(Sindin,2) * (
     &       yukfac(7)*sbotcont * dtan(beta)
     &       * Kbbg(mb2,mbsb2,mh2,msb12,msb22,s2thetab,c2thetab
     &              ,mgl,muSUSY,Ab,beta,Kfac,Kfacexp,mbrunyuk,tanbresum)
     &     + yukfac(6)*stopcont / dtan(beta)
     &       * Kttg(mt2,mh2,mst12,mst22,s2thetat,c2thetat
     &              ,mgl,At,muSUSY,beta))
     &     - 3/2.d0 * Amx(Sindin,3) * (
     &       yukfac(7)*sbotcont * dtan(beta)
     &       * KbbgS(mb2,mbsb2,mh2,msb12,msb22,s2thetab,c2thetab
     &              ,mgl,muSUSY,lam,GF,beta)
     &     + yukfac(6)*stopcont / dtan(beta) 
     &       * KttgS(mt2,mh2,mst12,mst22,s2thetat,c2thetat
     &              ,mgl,At,muSUSY,beta,lam,GF))
      endif
      end

C-}}}

C-{{{ scalar higgs

C-{{{ H1DS

      function H1DS(muRggh,mt2,mb2,mbmb,mbsb2,mh2,mgl,msb12,msb22
     &  ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &  ,lam,GF,muSUSY,beta,mw,mz,At,Ab
     &  ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac
     &  ,mbrunyuk,tanbresum,stopcont,sbotcont,Kfac,Kfacexp)
      implicit none
      double complex H1DS,F2lb,G2lb,D2lb,F2lt,D2lt
      double precision stopcont,sbotcont,Kfac,Kfacexp
      double precision yukfac(9)
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mt2,mst12,mst22,s2thetat,c2thetat
      double precision lam,GF,muSUSY,beta
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)
      double precision At,mw,mz,mbmb,Ab
      integer mbrunyuk,tanbresum

      H1DS = 
c     top:
     &(-dsqrt(mt2)*muSUSY 
     &   * F2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
     & + mz**2*dsin(2*beta)
     & * D2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,mw,mz))
     & / dsin(beta)
     & * yukfac(6) * stopcont
c     bottom:
     & +(dsqrt(mbsb2)*Ab*s2thetab
     &   * F2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab
     &     ,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
     & + 2*mb2 * G2lb(muRggh,mb2,mbsb2,mbmb,mh2,msb12,msb22,mgl
     &    ,s2thetab,c2thetab,muSUSY,dtan(beta),mw,mz,dmsb1,dmsb2,dthetab
     &    ,dmbsb,dAbds,tanbresum,mbrunyuk,Kfac,Kfacexp)
     & + 2*mz**2*dcos(beta)**2 
     &   * D2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &    ,s2thetab,c2thetab,mw,mz,dmsb1,dmsb2,dthetab,dmbsb,dAbds)) 
     & / dcos(beta) * yukfac(7) * sbotcont

#if(DRbartoOSb>0)
      H1DS = H1DS - dsqrt(mbsb2)/6.d0/dcos(beta)*s2thetab * (
     &       dAbds(3)*(1/msb12-1/msb22)
     &     + dAbds(2)*(1/msb12-1/msb22
     &     +2/15.d0*mh2/msb12**2-2/15.d0*mh2/msb22**2)
     &     + dAbds(1)*(1/msb12-1/msb22
     &     +2/15.d0*mh2/msb12**2-2/15.d0*mh2/msb22**2
     &     +3/140.d0*mh2**2/msb12**3-3/140.d0*mh2**2/msb22**3) )
     & * yukfac(7) * Kfac * sbotcont
#endif
#if( debug > 0)
      write(*,*) 'H1DS: ',H1DS
#endif
      end

C-}}}
C-{{{ H2DS

      function H2DS(muRggh,mt2,mb2,mbsb2,mh2,mgl,msb12,msb22
     &  ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &  ,lam,GF,muSUSY,beta,mw,mz,At
     &  ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac
     &  ,stopcont,sbotcont,Kfac,Kfacexp)
      implicit none
      double complex H2DS,F2lb,D2lb,F2lt,G2lt,D2lt
      double precision dAtosfins2,Kfac,Kfacexp
      double precision stopcont,sbotcont
      double precision yukfac(9)
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mt2,mst12,mst22,s2thetat,c2thetat
      double precision lam,GF,muSUSY,beta
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)
      double precision At,mw,mz

      H2DS = 
c     bottom:
     &   (-dsqrt(mbsb2)*muSUSY*s2thetab
     &    * F2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab
     &     ,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
     & - mz**2*dsin(2*beta)
     &   *D2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &    ,s2thetab,c2thetab,mw,mz,dmsb1,dmsb2,dthetab,dmbsb,dAbds))
     & / dcos(beta) * yukfac(7) * sbotcont
c     top:
     & +(((mst12-mst22)/2.d0*s2thetat+muSUSY/dtan(beta)*dsqrt(mt2))
     & * F2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
     & + 2*mt2 * G2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
     & - 2*mz**2*dsin(beta)**2
     & * D2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,mw,mz))
     & / dsin(beta)
     & * yukfac(6) * stopcont

c     DRbar -> OS
#if(DRbartoOSt>0)
      H2DS = H2DS - dsqrt(mt2)/6.d0/dsin(beta)*(
     &   dAtosfins2(1,At,muSUSY,beta,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,muRggh) * (1/mst12-1/mst22)
     &    + dAtosfins2(-1,At,muSUSY,beta,mt2,mst12,mst22
     &     ,mgl,s2thetat,c2thetat,muRggh)*2*mh2/15.d0*(1/mst12**2-1
     &     /mst22**2))* yukfac(6) * stopcont
#endif
#if(DRbartoOSt>1)
      write(*,*) 'H2DS: ',-dsqrt(mt2)/6.d0/dsin(beta)*(
     &   dAtosfins2(1,At,muSUSY,beta,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,muRggh) * (1/mst12-1/mst22)
     &    + dAtosfins2(-1,At,muSUSY,beta,mt2,mst12,mst22
     &     ,mgl,s2thetat,c2thetat,muRggh) *2*mh2/15.d0*(1/mst12**2-1
     &     /mst22**2))* yukfac(6) * stopcont
#endif
#if(debug > 0)
      write(*,*) 'H2DS: ',H2DS
#endif
      end

C-}}}
C-{{{ H3DS

      function H3DS(muRggh,mt2,mb2,mbsb2,mh2,mgl,msb12,msb22
     &  ,s2thetab,c2thetab,mst12,mst22,s2thetat,c2thetat
     &  ,lam,GF,muSUSY,beta
     &  ,dmsb1,dmsb2,dthetab,dmbsb,dAbds,yukfac,stopcont,sbotcont)
      implicit none
      double complex H3DS,F2lb,F2lt
      double precision stopcont,sbotcont,Sqrt2,vev,vevS
      double precision yukfac(9)
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mt2,mst12,mst22,s2thetat,c2thetat
      double precision lam,GF,muSUSY,beta
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam

      H3DS = 
c     bottom:
     &   (-dsqrt(mbsb2)*lam*vev*dsin(beta)*s2thetab/Sqrt2
     &    * F2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab
     &     ,dmsb1,dmsb2,dthetab,dmbsb,dAbds))
     &    * yukfac(7) * sbotcont / dcos(beta)
c     top: factor s2thetat included in F2lt!
     & + (-dsqrt(mt2)*lam*vev*dcos(beta)/Sqrt2
     &    * F2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat))
     &    * yukfac(6) * stopcont / dsin(beta)

#if(debug > 0)
      write(*,*) 'H3DS: ',H3DS
#endif
      end

C-}}}

C-{{{ bottom

C-{{{ F2lb

      function F2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
      implicit none
      double complex F2lb,YBS1,YBS2,YC2B
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)

      F2lb =YBS1(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)
     & - YBS2(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)
     & - 4.d0*c2thetab**2/(msb12-msb22)
     & *YC2B(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)

#if(DRbartoOSb>0)
      F2lb = F2lb + (
     &     2*dmsb1(1)/msb12
     &   - 2*dmsb2(1)/msb22
     & - (dmbsb(2)+2*c2thetab/s2thetab*dthetab(2))
     & * (1/msb12-1/msb22)
     & - (dmbsb(1)+2*c2thetab/s2thetab*dthetab(1))
     & * (1/msb12+2/15.d0*mh2/msb12**2-1/msb22-2/15.d0*mh2/msb22**2) )
     & / 6.d0
#endif
#if( debug > 0)
      write(*,*) 'F2lb: ',F2lb
#endif
      end

C-}}}
C-{{{ G2lb

      function G2lb(muRggh,mb2,mbsb2,mbmb,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,muSUSY,tanb,mw,mz,dmsb1,dmsb2,dthetab
     &,dmbsb,dAbds,tanbresum,mbrunyuk,Kfac,Kfacexp)
      implicit none
      double complex G2lb,G1l,I
      double precision Q,x1,x2,Kfac,Kfacexp,dli2,pi
      parameter (I=(0.d0,1.d0))
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mw,mz,mbmb,muSUSY,tanb
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)
      integer tanbresum,mbrunyuk

      pi = 3.14159265358979323846264338328d0

      Q = muRggh

      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2

      G2lb = 0.d0

      if (mbrunyuk.ge.3) then
         G2lb = G2lb 
     & -1/2.d0/mb2*G1l(mb2,mh2,.true.,0)*4/3.d0*( - 1.d0 *(
     & + (x1-3)/4/(1-x1)+x1*(x1-2)/2/(1-x1)**2*dlog(x1)
     & + (x2-3)/4/(1-x2)+x2*(x2-2)/2/(1-x2)**2*dlog(x2)
     &  )/4.d0) * Kfac

       if (mbrunyuk.eq.3) then
         G2lb = G2lb - 1/2.d0/mb2*G1l(mb2,mh2,.true.,0)*4/3.d0*(
     & - (2*dlog(Mgl/Q))/4.d0) * Kfac
       else if (mbrunyuk.eq.4) then
         G2lb = G2lb - 1/2.d0/mb2*G1l(mb2,mh2,.true.,0)*4/3.d0*(
     & - (2*dlog(Mgl/mbmb))/4.d0) * Kfac
       end if

       !Order in G1l can be changed
         G2lb = G2lb 
     & -1/2.d0/mb2*G1l(mb2,mh2,.true.,0)*4/3.d0*(
     & - (Mgl/dsqrt(mbsb2)*s2thetab*
     & (x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2)))/4.d0) * Kfac

      endif

      if ((mbrunyuk.eq.1).or.(mbrunyuk.eq.2)) then
      !same as tanbresum.ne.0 some lines below
         G2lb = G2lb 
     & -1/2.d0/mb2*G1l(mb2,mh2,.true.,0)*4/3.d0*(
     & (Mgl/dsqrt(mbsb2)*2*dsqrt(mbsb2)/(msb12-msb22)*muSUSY*tanb*
     & (x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2)))/4.d0) * Kfac
      endif

      G2lb = G2lb - 1/6.d0*G1l(mb2,mh2,.true.,0)*Mgl/mb2/dsqrt(mbsb2)* !/dsqrt(mb2)*
     &        s2thetab
     &        * (x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2)) * Kfac

      if (tanbresum.ne.0) then
       G2lb = G2lb - 1/6.d0*G1l(mb2,mh2,.true.,0)*Mgl/mb2/dsqrt(mbsb2)* !/dsqrt(mb2)*
     &        2*dsqrt(mbsb2)/(msb12-msb22)*muSUSY*tanb
     &        * (x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     &        * Kfacexp
      endif

      G2lb = G2lb - 2/3.d0/Mgl/dsqrt(mb2)*s2thetab*(dsqrt(mbsb2/mb2))*(
     &     ( 2*(1-x1)**3*dlog(Mgl/Q)+2*(x1**3+2*x1**2)*dlog(x1)
     &     - 3*(x1**3-x1-2*x1**2*dlog(x1))*(dlog(mh2/Mgl**2)-I*Pi)
     &     + 5*x1**3-5*x1**2+x1-1
     &     - 12*x1**2*dli2(1-1/x1)-6*x1**2*dlog(x1)**2
     &     ) /6/x1/(1-x1)**3-
     &     ( 2*(1-x2)**3*dlog(Mgl/Q)+2*(x2**3+2*x2**2)*dlog(x2)
     &     - 3*(x2**3-x2-2*x2**2*dlog(x2))*(dlog(mh2/Mgl**2)-I*Pi)
     &     + 5*x2**3-5*x2**2+x2-1
     &     - 12*x2**2*dli2(1-1/x2)-6*x2**2*dlog(x2)**2
     &     ) /6/x2/(1-x2)**3
     &     ) * Kfac

      G2lb = G2lb+3.d0/2.d0/Mgl/dsqrt(mb2)*s2thetab*(dsqrt(mbsb2/mb2))*(
     &     ( 2*x1*(1+x1)*dlog(x1)+2*x1-2-6*x1*dli2(1-1/x1)
     &     - 3*x1*dlog(x1)**2
     &     + 3*(1-x1+x1*dlog(x1))*(dlog(mh2/Mgl**2)-I*Pi)
     &     ) /6/(1-x1)**2 -
     &     ( 2*x2*(1+x2)*dlog(x2)+2*x2-2-6*x2*dli2(1-1/x2)
     &     - 3*x2*dlog(x2)**2
     &     + 3*(1-x2+x2*dlog(x2))*(dlog(mh2/Mgl**2)-I*Pi)
     &     ) /6/(1-x2)**2
     &     ) * Kfac

#if(DRbartoOSb>0)
      G2lb = G2lb-dmbsb(1)*(1/msb12+1/msb22)/3.d0*mbsb2/mb2
#endif
#if( debug > 0)
      write(*,*) "G2lb: ",G2lb
#endif
      end

C-}}}
C-{{{ D2lb

      function D2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mw,mz,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
      implicit none
      double complex D2lb,GS2lb,FS2lb
      double precision sthetaw2
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mw,mz
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)

      sthetaw2 = 1 - (Mw/Mz)**2

      D2lb = - 1/4.d0 * GS2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &         ,s2thetab,c2thetab,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
     &       + c2thetab*(sthetaw2/3.d0-1/4.d0) 
     &  * FS2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab
     &         ,dmsb1,dmsb2,dthetab,dmbsb,dAbds)

c      write(*,*) 'thetaw: ',dasin(dsqrt(sthetaw2))
#if( debug > 0)
      write(*,*) 'D2lb: ', D2lb
#endif
      end

C-}}}
C-{{{ GS2lb

      function GS2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
      implicit none
      double complex GS2lb,YBS1,YBS2
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)

      GS2lb=YBS1(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)
     & + YBS2(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)

#if(DRbartoOSb>0)
      GS2lb = GS2lb+(2*dmsb1(1)/msb12
     &             + 2*dmsb2(1)/msb22)/6.d0
#endif
      end

C-}}}
C-{{{ FS2lb

      function FS2lb(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,dmsb1,dmsb2,dthetab,dmbsb,dAbds)
      implicit none
      double complex FS2lb,YBS1,YBS2,YC2B
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab
      double precision dmsb1(2),dmsb2(2),dthetab(3),dmbsb(3),dAbds(3)

      FS2lb=YBS1(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)
     & - YBS2(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)
     & + 4.d0*s2thetab**2/(msb12-msb22)
     & *YC2B(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl,s2thetab,c2thetab)

#if(DRbartoOSb>0)
      FS2lb = FS2lb + ( 2*dmsb1(1)/msb12
     &                - 2*dmsb2(1)/msb22
     & + 2*s2thetab/c2thetab*dthetab(2)
     &     *(1/msb12-1/msb22)
     & + 2*s2thetab/c2thetab*dthetab(1)
     &     *(1/msb12+2/15.d0*mh2/msb12**2-1/msb22-2/15.d0*mh2/msb22**2)
     &     ) / 6.d0
#endif
      end

C-}}}
C-{{{ YBS1

      function YBS1(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab)
      implicit none
      double complex YBS1,G1l
      double precision x1,Q
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab

      Q = muRggh

      x1 = msb12 / Mgl**2

      YBS1 = -3/4.d0/msb12

      YBS1 = YBS1 - (
     &     (c2thetab**2*msb12+s2thetab**2*msb22)/msb12**2
     &    + s2thetab**2/msb12**2/msb22
     &    * (msb12**2*dlog(msb12/Q**2) - msb22**2*dlog(msb22/Q**2))
     &     ) / 18.d0

      YBS1 = YBS1 + s2thetab/3.d0/dsqrt(mbsb2)/Mgl*G1l(mb2,mh2,.true.,0)
     &  * (1/(1-x1)+1/(1-x1)**2*dlog(x1))
     &  - (1/(1-x1)+1/(1-x1)**2*dlog(x1)-1/x1**2*(1-2*dlog(Mgl/Q)))
     &    *2/9.d0/Mgl**2

      YBS1 = YBS1 - (1/(1-x1)+1/(1-x1)**2*dlog(x1))/4.d0/Mgl**2

#if( debug > 0)
      write(*,*) 'YBS1: ', YBS1
#endif
      end

C-}}}
C-{{{ YBS2

      function YBS2(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab)
      implicit none
      double complex YBS2,G1l
      double precision x2,Q
      double precision muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab

      Q = muRggh

      x2 = msb22 / Mgl**2

      YBS2 = -3/4.d0/msb22

      YBS2 = YBS2 - (
     &     (c2thetab**2*msb22+s2thetab**2*msb12)/msb22**2
     &    + s2thetab**2/msb22**2/msb12
     &    * (msb22**2*dlog(msb22/Q**2) - msb12**2*dlog(msb12/Q**2))
     &     ) / 18.d0

      YBS2 = YBS2 - s2thetab/3.d0/dsqrt(mbsb2)/Mgl*G1l(mb2,mh2,.true.,0)
     &  * (1/(1-x2)+1/(1-x2)**2*dlog(x2))
     &  - (1/(1-x2)+1/(1-x2)**2*dlog(x2)-1/x2**2*(1-2*dlog(Mgl/Q)))
     &    *2/9.d0/Mgl**2

      YBS2 = YBS2 - (1/(1-x2)+1/(1-x2)**2*dlog(x2))/4.d0/Mgl**2

#if( debug > 0)
      write(*,*) 'YBS2: ', YBS2
#endif
      end

C-}}}
C-{{{ YC2B

      function YC2B(muRggh,mb2,mbsb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab)
      implicit none
      double complex YC2B,G1l
      double precision x1,x2,Q
      double precision muRggh,mb2,mh2,msb12,msb22,mgl
     &,s2thetab,c2thetab,mbsb2

      Q = muRggh

      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2

      YC2B = - ( (msb12-msb22)**2/msb12/msb22
     &        - (msb12-msb22)/msb22*dlog(msb12/Q**2)
     &        - (msb22-msb12)/msb12*dlog(msb22/Q**2)
     &        ) / 18.d0

      YC2B = YC2B - Mgl/6.d0/dsqrt(mbsb2)/s2thetab*G1l(mb2,mh2,.true.,0)
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))

#if( debug > 0)
      write(*,*) 'YC2B: ', YC2B
#endif
      end

C-}}}

C-}}}
C-{{{ top

C-{{{ D2lt

      function D2lt(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat,mw,mz)
      implicit none
      double complex D2lt,GS2lt,FS2lt
      double precision sthetaw2
      double precision muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat,mw,mz

      sthetaw2 = 1 - (Mw/Mz)**2

      D2lt = 1/4.d0 * GS2lt(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat)
     &       + c2thetat*(-2*sthetaw2/3.d0+1/4.d0)
     & * FS2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)

#if( debug > 0)
      write(*,*) 'D2lt: ', D2lt
#endif
      end

C-}}}
C-{{{ F2lt

      function F2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
      implicit none
      double complex F2lt,YTS1,YTS2,YC2Ts2
      double precision dmst12osfin,dmst22osfin,dmtosfin,dthetatosfin
      double precision muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat

      !be careful due to additional factor s2thetat, also for CT
      F2lt = s2thetat*(
     &   YTS1(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,1)
     & - YTS2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,1))
     & - 4.d0*c2thetat**2/(mst12-mst22)
     &   *YC2Ts2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,1)

c     DRbar -> OS
#if(DRbartoOSt>0)
      F2lt = F2lt + 1.d0/6.d0 * s2thetat * (
     &  dmst12osfin(muRggh,mt2,mst12,mst22
     &       ,mgl,s2thetat,c2thetat,1)/mst12**2
     &  - dmst22osfin(muRggh,mt2,mst12,mst22
     &       ,mgl,s2thetat,c2thetat,1)/mst22**2
     &  - (dmtosfin(.true.,2,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2)
     &    + 2.d0*c2thetat/s2thetat*dthetatosfin(muRggh,mt2,mst12,mst22
     &                               ,mgl,s2thetat,c2thetat,1))
     &  * (1/mst12-1/mst22) 
     &  - 2.d0*mh2/15.d0 * (1.d0/mst12**2-1.d0/mst22**2)
     &    * dmtosfin(.true.,0,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2))
#endif
#if(DRbartoOSt>1)
      write(*,*) 'F2: ',1.d0/6.d0 * (
     &  dmst12osfin(muRggh,mt2,mst12,mst22
     &         ,mgl,s2thetat,c2thetat,1)/mst12**2
     &  - dmst22osfin(muRggh,mt2,mst12,mst22
     &         ,mgl,s2thetat,c2thetat,1)/mst22**2
     &  - (dmtosfin(.true.,2,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2)
     &    + 2.d0*c2thetat/s2thetat*dthetatosfin(muRggh,mt2,mst12,mst22
     &                                ,mgl,s2thetat,c2thetat,1))
     &  * (1/mst12-1/mst22) 
     &  - 2.d0*mh2/15.d0 * (1.d0/mst12**2-1.d0/mst22**2)
     &    * dmtosfin(.true.,0,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2))
#endif
#if( debug > 0)
      write(*,*) 'F2lt: ',F2lt
#endif
      end

C-}}}
C-{{{ G2lt

      function G2lt(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
      implicit none
      double complex G2lt,YTS1,YTS2,F1l,YTh
      double precision dmst12osfin,dmst22osfin,dmtosfin
      double precision muRggh,mt2,mst12,mst22,mgl,s2thetat,c2thetat,mh2

      G2lt = YTS1(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)
     &   + YTS2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)
     &   + YTh(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)

c     DRbar -> OS
#if(DRbartoOSt>0)
      G2lt = G2lt + (dmst12osfin(muRggh,mt2,mst12,mst22
     &          ,mgl,s2thetat,c2thetat,0)/mst12**2
     &   + dmst22osfin(muRggh,mt2,mst12,mst22
     &          ,mgl,s2thetat,c2thetat,0)/mst22**2
     & - 2*dmtosfin(.true.,1,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2)*(1/mst12+1/mst22) ) /6.d0
     & - 2/3.d0*F1l(mt2,mh2)*dmtosfin(.false.,3,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/mt2/dsqrt(mt2)
#endif
#if(DRbartoOSt>1)
      write(*,*) 'G2: ', (dmst12osfin(muRggh,mt2,mst12,mst22
     &          ,mgl,s2thetat,c2thetat,0)/mst12**2
     & + dmst22osfin(muRggh,mt2,mst12,mst22
     &          ,mgl,s2thetat,c2thetat,0)/mst22**2
     & - 2*dmtosfin(.true.,1,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/dsqrt(mt2)*(1/mst12+1/mst22) ) /6.d0
     & - 2/3.d0*F1l(mt2,mh2)*dmtosfin(.false.,3,mt2,mgl
     &,s2thetat,mst12,mst22,muRggh)/mt2/dsqrt(mt2)
#endif
#if( debug > 0)
      write(*,*) "G2lt: ",G2lt
#endif
      end

C-}}}
C-{{{ FS2lt

      function FS2lt(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat)
      implicit none
      double complex FS2lt,YTS1,YTS2,YC2Ts2
      double precision dmst12osfin,dmst22osfin,dthetatosfin
      double precision muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat

      FS2lt = YTS1(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)
     & - YTS2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)
     & + 4.d0*s2thetat/(mst12-mst22)
     &   *YC2Ts2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)

c     DRbar -> OS
#if(DRbartoOSt>0)
      FS2lt = FS2lt + (dmst12osfin(muRggh,mt2,mst12,mst22
     &       ,mgl,s2thetat,c2thetat,0)/mst12**2
     & - dmst22osfin(muRggh,mt2,mst12,mst22
     &       ,mgl,s2thetat,c2thetat,0)/mst22**2
     & + 2*s2thetat*dthetatosfin(muRggh,mt2,mst12,mst22
     &    ,mgl,s2thetat,c2thetat,0)/c2thetat*(1/mst12-1/mst22)
     &   ) / 6.d0
#endif
#if(DRbartoOSt>1)
      write(*,*) 'FS2: ',(dmst12osfin(muRggh,mt2,mst12,mst22
     &        ,mgl,s2thetat,c2thetat,0)/mst12**2
     & - dmst22osfin(muRggh,mt2,mst12,mst22
     &        ,mgl,s2thetat,c2thetat,0)/mst22**2
     & + 2*s2thetat*dthetatosfin(muRggh,mt2,mst12,mst22
     &    ,mgl,s2thetat,c2thetat,0)/c2thetat*(1/mst12-1/mst22)
     &   ) / 6.d0
#endif
#if( debug > 0)
      write(*,*) 'FS2lt: ',FS2lt
#endif
      end

C-}}}
C-{{{ GS2lt

      function GS2lt(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat)
      implicit none
      double complex GS2lt,YTS1,YTS2
      double precision dmst12osfin,dmst22osfin
      double precision muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat

      GS2lt = YTS1(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)
     & + YTS2(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat,0)

c     DRbar -> OS
#if(DRbartoOSt>0)
      GS2lt = GS2lt + (dmst12osfin(muRggh,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,0)/mst12**2
     &  + dmst22osfin(muRggh,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,0)/mst22**2)/6.d0
#endif
#if(DRbartoOSt>1)
      write(*,*) 'GS2: ', 
     &(dmst12osfin(muRggh,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,0)/mst12**2
     &+dmst22osfin(muRggh,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,0)/mst22**2)/6.d0
#endif
#if( debug > 0)
      write(*,*) 'GS2lt: ',GS2lt
#endif
      end

C-}}}

C-{{{ YTS1

      function YTS1(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat,flag)
      implicit none
      double complex YTS1,G1l
      double precision x1,Q,dli2
      integer flag
      double precision muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat

      Q = muRggh

      x1 = mst12 / Mgl**2

      YTS1 = -3/4.d0/mst12

#if( debug > 1)
      write(*,*) 'YTS1a: ', YTS1
#endif

      YTS1 = YTS1 - (
     &     (c2thetat**2*mst12+s2thetat**2*mst22)/mst12**2
     &    + s2thetat**2/mst12**2/mst22
     &    * (mst12**2*dlog(mst12/Q**2) - mst22**2*dlog(mst22/Q**2))
     &     ) / 18.d0

#if( debug > 1)
      write(*,*) 'YTS1b: ', - (
     &     (c2thetat**2*mst12+s2thetat**2*mst22)/mst12**2
     &    + s2thetat**2/mst12**2/mst22
     &    * (mst12**2*dlog(mst12/Q**2) - mst22**2*dlog(mst22/Q**2))
     &     ) / 18.d0
#endif

      YTS1 = YTS1 + (s2thetat/3.d0/dsqrt(mt2)/Mgl*G1l(mt2,mh2,.false.,0)
     & - 17/36.d0/Mgl**2)*(1/(1-x1)+1/(1-x1)**2*dlog(x1))
     & + 1/18.d0/(Mgl*x1)**2/(1-x1)**3*( 4*(1-x1)**3*(1-2*dlog(Mgl/Q))
     &   - 3*x1**2*G1l(mt2,mh2,.false.,0)*((1-x1)*(3-x1)+2*dlog(x1)) )

      if(flag.eq.1) then
         YTS1 = YTS1 + s2thetat*dsqrt(mt2)*2/9.d0/Mgl**3/x1**2/(1-x1)**4
     & * ( 3*x1**2 * ((1-x1)*(x1+5)+2*(2*x1+1)*dlog(x1))
     & * ( G1l(mt2,mh2,.false.,0)/4.d0-dlog(mt2/Mgl**2) ) 
     & + (1-x1)**4*2*dlog(Mgl/Q)
     & - 12*x1**2*(2*x1+1)*dli2(1-x1)-(1-x1)*(14*x1**2-3*x1+1)-2*x1**2
     & * (x1**2+18*x1+5)*dlog(x1) )
     & + s2thetat*dsqrt(mt2)/Mgl**3/2.d0/(1-x1)**3*(
     &   3*(2-2*x1+(x1+1)*dlog(x1))*(1+dlog(mt2/Mgl**2))
     & + 6*(1+x1)*dli2(1-x1)+2*x1*(1-x1)+2*(6*x1+1)*dlog(x1) )
     & + s2thetat*mh2*G1l(mt2,mh2,.false.,0)/36.d0/dsqrt(mt2)
     &   /Mgl**3/x1/(1-x1)**4*(
     &   (1-x1)*(x1**2-5*x1-2)-6*x1*dlog(x1) )
      endif

#if( debug > 1)
      write(*,*) 'YTS1c: ',(s2thetat/3.d0/dsqrt(mt2)
     & /Mgl*G1l(mt2,mh2,.false.,0)
     & - 17/36.d0/Mgl**2)*(1/(1-x1)+1/(1-x1)**2*dlog(x1))
     & + 1/18.d0/(Mgl*x1)**2/(1-x1)**3*( 4*(1-x1)**3*(1-2*dlog(Mgl/Q))
     &   - 3*x1**2*G1l(mt2,mh2,.false.,0)*((1-x1)*(3-x1)+2*dlog(x1)))
     & + s2thetat*dsqrt(mt2)*2/9.d0/Mgl**3/x1**2/(1-x1)**4*( 3*x1**2
     & * ( (1-x1)*(x1+5)+2*(2*x1+1)*dlog(x1) )
     & * ( G1l(mt2,mh2,.false.,0)/4.d0-dlog(mt2/Mgl**2) )
     & + (1-x1)**4*2*dlog(Mgl/Q)
     & - 12*x1**2*(2*x1+1)*dli2(1-x1)-(1-x1)*(14*x1**2-3*x1+1)-2*x1**2
     & * (x1**2+18*x1+5)*dlog(x1) )
     & + s2thetat*dsqrt(mt2)/Mgl**3/2.d0/(1-x1)**3*(
     &   3*(2-2*x1+(x1+1)*dlog(x1))*(1+dlog(mt2/Mgl**2))
     & + 6*(1+x1)*dli2(1-x1)+2*x1*(1-x1)+2*(6*x1+1)*dlog(x1) )
     & + s2thetat*mh2*G1l(mt2,mh2,.false.,0)/36.d0/dsqrt(mt2)
     & /Mgl**3/x1/(1-x1)**4*(
     &   (1-x1)*(x1**2-5*x1-2)-6*x1*dlog(x1) )
#endif
#if( debug > 0)
      write(*,*) 'YTS1: ', YTS1
#endif
      end

C-}}}
C-{{{ YTS2

      function YTS2(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat,flag)
      implicit none
      double complex YTS2,G1l
      double precision x2,Q,dli2
      integer flag
      double precision muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat

      Q = muRggh

      x2 = mst22 / Mgl**2

      YTS2 = -3/4.d0/mst22

#if( debug > 1)
      write(*,*) 'YTS2a: ', YTS2
#endif

      YTS2 = YTS2 - (
     &     (c2thetat**2*mst22+s2thetat**2*mst12)/mst22**2
     &    + s2thetat**2/mst22**2/mst12
     &    * (mst22**2*dlog(mst22/Q**2) - mst12**2*dlog(mst12/Q**2))
     &     ) / 18.d0

#if( debug > 1)
      write(*,*) 'YTS2b: ', - (
     &     (c2thetat**2*mst22+s2thetat**2*mst12)/mst22**2
     &    + s2thetat**2/mst22**2/mst12
     &    * (mst22**2*dlog(mst22/Q**2) - mst12**2*dlog(mst12/Q**2))
     &     ) / 18.d0
#endif

      YTS2 = YTS2 +(-s2thetat/3.d0/dsqrt(mt2)/Mgl*G1l(mt2,mh2,.false.,0)
     & - 17/36.d0/Mgl**2)*(1/(1-x2)+1/(1-x2)**2*dlog(x2))
     & + 1/18.d0/(Mgl*x2)**2/(1-x2)**3*( 4*(1-x2)**3*(1-2*dlog(Mgl/Q))
     &   - 3*x2**2*G1l(mt2,mh2,.false.,0)*((1-x2)*(3-x2)+2*dlog(x2)))

      if(flag.eq.1) then
         YTS2 = YTS2 - s2thetat*dsqrt(mt2)*2/9.d0/Mgl**3/x2**2/(1-x2)**4
     & * ( 3*x2**2 * ((1-x2)*(x2+5)+2*(2*x2+1)*dlog(x2))
     & * ( G1l(mt2,mh2,.false.,0)/4.d0-dlog(mt2/Mgl**2) ) 
     & + (1-x2)**4*2*dlog(Mgl/Q)
     & - 12*x2**2*(2*x2+1)*dli2(1-x2)-(1-x2)*(14*x2**2-3*x2+1)-2*x2**2
     & * (x2**2+18*x2+5)*dlog(x2) )
     & - s2thetat*dsqrt(mt2)/Mgl**3/2.d0/(1-x2)**3*(
     &   3*(2-2*x2+(x2+1)*dlog(x2))*(1+dlog(mt2/Mgl**2))
     & + 6*(1+x2)*dli2(1-x2)+2*x2*(1-x2)+2*(6*x2+1)*dlog(x2) )
     & - s2thetat*mh2*G1l(mt2,mh2,.false.,0)/36.d0/dsqrt(mt2)
     & /Mgl**3/x2/(1-x2)**4*(
     &   (1-x2)*(x2**2-5*x2-2)-6*x2*dlog(x2) )
         endif

#if( debug > 1)
      write(*,*) 'YTS2c: ',(-s2thetat/3.d0/dsqrt(mt2)
     & /Mgl*G1l(mt2,mh2,.false.,0)
     & - 17/36.d0/Mgl**2)*(1/(1-x2)+1/(1-x2)**2*dlog(x2))
     & + 1/18.d0/(Mgl*x2)**2/(1-x2)**3*( 4*(1-x2)**3*(1-2*dlog(Mgl/Q))
     &   - 3*x2**2*G1l(mt2,mh2,.false.,0)*((1-x2)*(3-x2)+2*dlog(x2)))
     & - s2thetat*dsqrt(mt2)*2/9.d0/Mgl**3/x2**2/(1-x2)**4*( 3*x2**2
     & * ( (1-x2)*(x2+5)+2*(2*x2+1)*dlog(x2) )
     & * ( G1l(mt2,mh2,.false.,0)/4.d0-dlog(mt2/Mgl**2) ) 
     & + (1-x2)**4*2*dlog(Mgl/Q)
     & - 12*x2**2*(2*x2+1)*dli2(1-x2)-(1-x2)*(14*x2**2-3*x2+1)-2*x2**2
     & * (x2**2+18*x2+5)*dlog(x2) )
     & - s2thetat*dsqrt(mt2)/Mgl**3/2.d0/(1-x2)**3*(
     &   3*(2-2*x2+(x2+1)*dlog(x2))*(1+dlog(mt2/Mgl**2))
     & + 6*(1+x2)*dli2(1-x2)+2*x2*(1-x2)+2*(6*x2+1)*dlog(x2) )
     & - s2thetat*mh2*G1l(mt2,mh2,.false.,0)/36.d0/dsqrt(mt2)
     & /Mgl**3/x2/(1-x2)**4*(
     &   (1-x2)*(x2**2-5*x2-2)-6*x2*dlog(x2) )
#endif
#if( debug > 0)
      write(*,*) 'YTS2: ', YTS2
#endif
      end

C-}}}
C-{{{ YC2Ts2

      function YC2Ts2(muRggh,mt2,mh2,mst12,mst22,mgl
     &,s2thetat,c2thetat,flag)
      implicit none
      double complex YC2Ts2,G1l
      double precision x1,x2,Q,dli2
      integer flag
      double precision muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat

      Q = muRggh

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      YC2Ts2 = - ( (mst12-mst22)**2/mst12/mst22
     &        - (mst12-mst22)/mst22*dlog(mst12/Q**2)
     &        - (mst22-mst12)/mst12*dlog(mst22/Q**2)
     &        ) / 18.d0 * s2thetat

      YC2Ts2 = YC2Ts2 - Mgl/6.d0/dsqrt(mt2)*G1l(mt2,mh2,.false.,0)
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))

      if(flag.eq.1) then
         YC2Ts2 = YC2Ts2 + (dsqrt(mt2)/9.d0/Mgl/x1/(1-x1)**3 * (
     &   3*x1*(x1**2-2*x1*dlog(x1)-1)*(G1l(mt2,mh2,.false.,0)/4.d0
     & -dlog(mt2/Mgl**2))
     & + (1-x1)**3*2*dlog(Mgl/Q)+12*x1**2*dli2(1-x1)-(1-2*x1)*(1-x1)**2
     & + 2*x1**2*(x1+5)*dlog(x1) )
     & + dsqrt(mt2)*x1/4.d0/Mgl/(1-x1)**2* ( (x1-1-dlog(x1))
     & * (3*dlog(mt2/Mgl**2)+1)-6*dli2(1-x1)-2*(x1+2)*dlog(x1) )
     & + mh2*G1l(mt2,mh2,.false.,0)/24.d0/dsqrt(mt2)/Mgl/(1-x1)**2 * (
     &   (1-x1)*(x1+x2-2*x1*x2)/(1-x2)/(x1-x2)
     & + 2*x1*(x1**2+x1*x2-2*x2)*dlog(x1)/(x1-x2)**2 )
     & - dsqrt(mt2)/9.d0/Mgl/x2/(1-x2)**3 * (
     &   3*x2*(x2**2-2*x2*dlog(x2)-1)*(G1l(mt2,mh2,.false.,0)/4.d0
     & -dlog(mt2/Mgl**2))
     & + (1-x2)**3*2*dlog(Mgl/Q)+12*x2**2*dli2(1-x2)-(1-2*x2)*(1-x2)**2
     & + 2*x2**2*(x2+5)*dlog(x2) )
     & - dsqrt(mt2)*x2/4.d0/Mgl/(1-x2)**2* ( (x2-1-dlog(x2))
     & * (3*dlog(mt2/Mgl**2)+1)-6*dli2(1-x2)-2*(x2+2)*dlog(x2) )
     & - mh2*G1l(mt2,mh2,.false.,0)/24.d0/dsqrt(mt2)/Mgl/(1-x2)**2 * (
     &   (1-x2)*(x2+x1-2*x2*x1)/(1-x1)/(x2-x1)
     & + 2*x2*(x2**2+x2*x1-2*x1)*dlog(x2)/(x2-x1)**2 )
     &)
      endif

#if( debug > 0)
      write(*,*) 'YC2Ts2: ', YC2Ts2,YC2Ts2/s2thetat
#endif
      end

C-}}}
C-{{{ YTh

      function YTh(muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat)
      implicit none
      double complex YTh,R1T,R2T,F1l,G1l
      double precision x1,x2,Q
      double precision muRggh,mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat

      Q = muRggh

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      YTh = -2/9.d0/mt2*F1l(mt2,mh2)*(2*dlog(Mgl/Q)
     & + Mgl/dsqrt(mt2)*s2thetat*x1/(1-x1)*dlog(x1)
     & + (x1-3)/4.d0/(1-x1)+x1*(x1-2)/2.d0/(1-x1)**2*dlog(x1)
     & + dsqrt(mt2)/Mgl*s2thetat/2.d0/(1-x1)**3*(1-x1**2+2*x1*dlog(x1))
     & + mt2/Mgl**2/6.d0/(1-x1)**3*(x1**2-5*x1-2-6*x1/(1-x1)*dlog(x1))
     & - Mgl/dsqrt(mt2)*s2thetat*x2/(1-x2)*dlog(x2)
     & + (x2-3)/4.d0/(1-x2)+x2*(x2-2)/2.d0/(1-x2)**2*dlog(x2)
     & - dsqrt(mt2)/Mgl*s2thetat/2.d0/(1-x2)**3*(1-x2**2+2*x2*dlog(x2))
     & + mt2/Mgl**2/6.d0/(1-x2)**3*(x2**2-5*x2-2-6*x2/(1-x2)*dlog(x2)))

      YTh = YTh - 1/6.d0*G1l(mt2,mh2,.false.,0)
     & * Mgl/mt2/dsqrt(mt2)*s2thetat
     & * (x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     & + R1T(muRggh,mt2,mh2,mst12,mst22,mgl)
     &   * s2thetat/Mgl/dsqrt(mt2)/2.d0
     & + R2T(muRggh,mt2,mh2,mst12,mst22,mgl) /Mgl**2/2.d0

#if( debug > 0)
      write(*,*) 'YT: ', Yth
#endif
      end
C-}}}
C-{{{ R1T

c$$$      eq. (24) from 1204.1016
      function R1T(muRggh,mt2,mh2,mst12,mst22,mgl)
      implicit none
      double complex R1T,BDS,G1l,K1l
      double precision x1,x2,Q,dli2
      double precision muRggh,mt2,mh2,mst12,mst22,mgl

      Q = muRggh

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      R1T = 1/2.d0/(1-x1)**2*(3*(1-x1+x1*dlog(x1))
     & * (dlog(mt2/Mgl**2)-BDS(mt2,mh2)-K1l(mt2,mh2)/2.d0+2)
     & + 6*x1*dli2(1-x1)+2*x1+2*x1*(1+x1)*dlog(x1)-2)
     & - 2/9.d0/x1/(1-x1)**3*(3*(x1-x1**3+2*x1**2*dlog(x1))
     & * (dlog(mt2/Mgl**2)-BDS(mt2,mh2)-G1l(mt2,mh2,.false.,0)/4.d0
     & -K1l(mt2,mh2)/2.d0+2)
     & + (1-x1)**3*2*dlog(Mgl/Q)+12*x1**2*dli2(1-x1)+5*x1**3-5*x1**2+x1
     & - 1+2*(x1**3+2*x1**2)*dlog(x1))
     & - 1/2.d0/(1-x2)**2*(3*(1-x2+x2*dlog(x2))
     & * (dlog(mt2/Mgl**2)-BDS(mt2,mh2)-K1l(mt2,mh2)/2.d0+2)
     & + 6*x2*dli2(1-x2)+2*x2+2*x2*(1+x2)*dlog(x2)-2)
     & + 2/9.d0/x2/(1-x2)**3*(3*(x2-x2**3+2*x2**2*dlog(x2))
     & * (dlog(mt2/Mgl**2)-BDS(mt2,mh2)-G1l(mt2,mh2,.false.,0)/4.d0
     & -K1l(mt2,mh2)/2.d0+2)
     & + (1-x2)**3*2*dlog(Mgl/Q)+12*x2**2*dli2(1-x2)+5*x2**3-5*x2**2+x2
     & - 1+2*(x2**3+2*x2**2)*dlog(x2))

      end

C-}}}
C-{{{ R2T

c$$$      eq. (25) from 1204.1016
      function R2T(muRggh,mt2,mh2,mst12,mst22,mgl)
      implicit none
      double complex R2T,BDS,G1l,K1l
      double precision x1,x2,Q,dli2
      double precision muRggh,mt2,mh2,mst12,mst22,mgl

      Q = muRggh

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      R2T = -1/4.d0/(1-x1)**3*(3*(1-x1**2+2*x1*dlog(x1))
     & * (2*dlog(mt2/Mgl**2)-BDS(mt2,mh2)-K1l(mt2,mh2)/2.d0+2)
     & + 24*x1*dli2(1-x1)+1-x1**2+2*x1*(3*x1+10)*dlog(x1))
     & + 2/27.d0/x1/(1-x1)**4 * (
     &   3*x1*((1-x1)*(5*x1-x1**2+2)+6*x1*dlog(x1))
     & * (2*dlog(mt2/Mgl**2)-BDS(mt2,mh2)-G1l(mt2,mh2,.false.,0)/2.d0
     & -K1l(mt2,mh2)/2.d0+2)
     & + 12*(1-x1)**4*dlog(Mgl/Q)+72*x1**2*dli2(1-x1)
     & - x1*(1-x1)**2*(11*x1-26)-6*(1-x1)+6*x1**2*(2*x1+9)*dlog(x1))
     & -1/4.d0/(1-x2)**3*(3*(1-x2**2+2*x2*dlog(x2))
     & * (2*dlog(mt2/Mgl**2)-BDS(mt2,mh2)-K1l(mt2,mh2)/2.d0+2)
     & + 24*x2*dli2(1-x2)+1-x2**2+2*x2*(3*x2+10)*dlog(x2))
     & + 2/27.d0/x2/(1-x2)**4 * (
     &   3*x2*((1-x2)*(5*x2-x2**2+2)+6*x2*dlog(x2))
     & * (2*dlog(mt2/Mgl**2)-BDS(mt2,mh2)-G1l(mt2,mh2,.false.,0)/2.d0
     & -K1l(mt2,mh2)/2.d0+2)
     & + 12*(1-x2)**4*dlog(Mgl/Q)+72*x2**2*dli2(1-x2)
     & - x2*(1-x2)**2*(11*x2-26)-6*(1-x2)+6*x2**2*(2*x2+9)*dlog(x2))

      end

C-}}}

C-}}}

C-{{{ F1l

      function F1l(mq2,mh2)
      implicit none
      double complex F1l,Integral3,ftauq
      double precision tau,mq2,mh2

      tau = mh2/mq2/4.d0
      ftauq = - Integral3(0.d0,mh2,mq2) * mh2 / 2.d0
      if(Real(cdsqrt(ftauq)*(1-tau)).gt.0.d0) then
         F1l = 3/tau**2*(tau+(tau-2)*ftauq
     &        +cdsqrt((1.d0,0.d0)*(1-tau)*tau)*cdsqrt(ftauq))
      else
         F1l = 3/tau**2*(tau+(tau-2)*ftauq
     &        -cdsqrt((1.d0,0.d0)*(1-tau)*tau)*cdsqrt(ftauq))
      endif

#if( debug > 0)
      write(*,*) "F1l: ", F1l
#endif
      end

C-}}}
C-{{{ G1l

      function G1l(mq2,mh2,expflag,exporder)
      implicit none
      double complex G1l,ALOSMh, I
      double precision mq2, taub, mh2, pi
      logical expflag
      integer exporder

      pi = 3.14159265358979323846264338328d0

      if (expflag.eqv.(.true.)) then
         I = (0.d0,1.d0)
         taub = 4.d0 * mq2 / mh2
         G1l = -2.d0 * taub + taub / 2.d0 * (dLog(4.d0/taub)-I*Pi)**2
      if (exporder.eq.1) then
         G1l = G1l
     &-taub**2/2.d0*(dLog(4.d0/taub)-I*Pi+(dLog(4.d0/taub)-I*Pi)**2)
      endif
      else
         G1l = -4/3.d0*mq2*ALOSMh(mh2,mq2)
      endif

#if( debug > 0)
      write(*,*) "G1l: ", G1l
#endif
      end

C-}}}
C-{{{ B

      function BDS(mq2,mh2)
      implicit none
      double complex BDS,Integral2
      double precision mq2,mh2

      BDS = 2.d0+Integral2(mh2,mq2)

      end

C-}}}
C-{{{ K1l

      function K1l(mq2,mh2)
      implicit none
      double complex K1l,Integral3
      double precision mq2,mh2

      K1l = 4.d0*mq2*Integral3(0.d0,mh2,mq2)

      end

C-}}}

C-}}}
C-{{{ pseudoscalar higgs

C-{{{ bottom

C-{{{ Kbbg

      function Kbbg(mb2,mbsb2,mh2,msb12,msb22,s2thetab,c2thetab
     &,mgl,muSUSY,Ab,beta,Kfac,Kfacexp,mbrunyuk,tanbresum)
      implicit none
      double complex Kbbg,K1b,R1b
      double precision x1,x2,Yb
      double precision Kfac, Kfacexp
      double precision s2thetab,c2thetab,mbsb2,tanb
      double precision mb2,mh2,msb12,msb22,mgl,muSUSY,Ab,beta
      integer mbrunyuk,tanbresum

      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2
      Yb = Ab + muSUSY/dtan(beta)
      tanb = dtan(beta)

!      Kbbg = 2/3.d0*K1b(mb2,mh2,.true.)*Mgl*muSUSY/(msb12-msb22)*(tanb+1/tanb)
      Kbbg = -2/3.d0*K1b(mb2,mh2,.true.)*Mgl/dsqrt(mbsb2)
     &     *(s2thetab/2.d0 * Kfac - dsqrt(mbsb2)*Yb/(msb12-msb22))
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     &     -dsqrt(mbsb2)/Mgl*s2thetab
     &      *R1b(mh2,msb12,msb22,mgl,Ab,muSUSY,beta,Kfac)

      if ((mbrunyuk.eq.1).or.(mbrunyuk.eq.2)) then
      Kbbg = Kbbg - 2/3.d0*K1b(mb2,mh2,.true.)
     &     *Mgl*muSUSY/(msb12-msb22)*tanb
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2)) * Kfac
      endif

      if (tanbresum.ne.0) then
      Kbbg = Kbbg - 2/3.d0*K1b(mb2,mh2,.true.)
     &     *Mgl*muSUSY/(msb12-msb22)*tanb
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2)) * Kfacexp
      end if

#if( debug > 0)
      write(*,*) 'Kbbg: ',Kbbg
#endif
      end

C-}}}
C-{{{ KbbgS

      function KbbgS(mb2,mbsb2,mh2,msb12,msb22,s2thetab,c2thetab
     &,mgl,muSUSY,lam,GF,beta)
      implicit none
      double complex KbbgS,K1b,R1bS
      double precision x1,x2,Ysave
      double precision vev,vevS,Sqrt2,mbsb2,s2thetab,c2thetab
      double precision mb2,mh2,msb12,msb22,mgl,muSUSY,lam,GF,beta

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam

      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2
      Ysave = lam*vev/Sqrt2

      KbbgS = 2/3.d0*K1b(mb2,mh2,.true.)*Mgl/dsqrt(mbsb2)
     &     *(dsqrt(mbsb2)*Ysave/(msb12-msb22))
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     &     -dsqrt(mbsb2)/Mgl*s2thetab
     &       *R1bS(msb12,msb22,mgl,muSUSY,lam,GF)

#if( debug > 0)
      write(*,*) 'KbbgS: ',KbbgS
#endif
      end

C-}}}
C-{{{ K1b

      function K1b(mb2,mh2,expflag)
      implicit none
      logical expflag
      double complex K1b,ALOSMA
      double precision taubA,pi,mb2,mh2

      pi = 3.14159265358979323846264338328d0

      if (expflag.eqv.(.true.)) then
      taubA = 4.d0*mb2/mh2
      K1b = taubA/2d0*(dcmplx(log(4/taubA),-pi))**2
      else
      K1b = -4/3.d0*mb2*ALOSMA(mh2,mb2)
      end if

      end

C-}}}
C-{{{ R1b

      function R1b(mh2,msb12,msb22,mgl,Ab,muSUSY,beta,Kfac)
      implicit none
      double complex R1b,I
      double precision x1,x2,Yb,Kfac,dli2
      double precision msb12,msb22,mgl,Ab,muSUSY,beta,mh2,pi

      pi = 3.14159265358979323846264338328d0

      I = (0.d0,1.d0)
      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2
      Yb = Ab + muSUSY/dtan(beta)

      R1b = Kfac*1/3.d0/(1-x1)**3 * (
     &       (1-x1**2+2*x1*dlog(x1))*(1-2*(dlog(mh2/Mgl**2)-I*Pi))
     &      - 8*x1*dli2(1-x1)-2*x1*(3+x1)*dlog(x1) )
     &   + Kfac*3/2.d0/(1-x1)**2 * (
     &       (1-x1+x1*dlog(x1))*((dlog(mh2/Mgl**2)-I*Pi)-1)
     &      + 2*x1*dli2(1-x1)+x1*(1+x1)*dlog(x1) )
     &   + 4/3.d0/(x1-x2)**2*Yb/Mgl * (
     &      x1**2*(1-2*x2)/2.d0/(1-x1)/(1-x2)+x1/2.d0/(1-x1)**2
     &     *(x1**2-2*x2+x1*x2)*dlog(x1) )
     &   - Kfac*1/3.d0/(1-x2)**3 * (
     &       (1-x2**2+2*x2*dlog(x2))*(1-2*(dlog(mh2/Mgl**2)-I*Pi))
     &      - 8*x2*dli2(1-x2)-2*x2*(3+x2)*dlog(x2) )
     &   - Kfac*3/2.d0/(1-x2)**2 * (
     &       (1-x2+x2*dlog(x2))*((dlog(mh2/Mgl**2)-I*Pi)-1)
     &      + 2*x2*dli2(1-x2)+x2*(1+x2)*dlog(x2) )
     &   - 4/3.d0/(x2-x1)**2*Yb/Mgl * (
     &      x2**2*(1-2*x1)/2.d0/(1-x2)/(1-x1)+x2/2.d0/(1-x2)**2
     &     *(x2**2-2*x1+x2*x1)*dlog(x2) )
      end

C-}}}
C-{{{ R1bS

      function R1bS(msb12,msb22,mgl,muSUSY,lam,GF)
      implicit none
      double complex R1bS
      double precision x1,x2,Ysave,dli2
      double precision vev,vevS,Sqrt2
      double precision muSUSY,lam,GF,msb12,msb22,mgl

      x1 = msb12 / Mgl**2
      x2 = msb22 / Mgl**2

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      Ysave = lam*vev/Sqrt2

      R1bS = 4/3.d0/(x1-x2)**2*Ysave/Mgl * (
     &      x1**2*(1-2*x2)/2.d0/(1-x1)/(1-x2)+x1/2.d0/(1-x1)**2
     &     *(x1**2-2*x2+x1*x2)*dlog(x1) )
     &   - 4/3.d0/(x2-x1)**2*Ysave/Mgl * (
     &      x2**2*(1-2*x1)/2.d0/(1-x2)/(1-x1)+x2/2.d0/(1-x2)**2
     &     *(x2**2-2*x1+x2*x1)*dlog(x2) )

      end

C-}}}

C-}}}
C-{{{ top

C-{{{ Kttg

      function Kttg(mt2,mh2,mst12,mst22,s2thetat,c2thetat
     &,mgl,At,muSUSY,beta)
      implicit none
      double complex Kttg,K1t,R1tA,R2tA,R3tA,R4tA
      double precision x1,x2,Yt
      double precision mt2,mh2,mst12,mst22,s2thetat,c2thetat
     &,mgl,At,muSUSY,beta

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2
      Yt = At + muSUSY*dtan(beta)

c      Kttg = 2/3.d0*K1t(mt2,mh2)*Mgl*muSUSY/(mst12-mst22)*(tanb+1/tanb)
      Kttg = -2/3.d0*K1t(mt2,mh2)*Mgl/dsqrt(mt2)
     &     *(s2thetat/2.d0-dsqrt(mt2)*Yt/(mst12-mst22))
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     &     -dsqrt(mt2)/Mgl*s2thetat
     &    *R1tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
     &  + 2.d0*mt2*Yt/Mgl/(mst12-mst22)
     &    *R2tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
     &  + mt2/(Mgl**2)*R3tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
     &  - K1t(mt2,mh2)/2.d0*mh2/(mst12-mst22)
     &    *R4tA(mst12,mst22,mgl,At,muSUSY,beta)

#if( debug > 0)
      write(*,*) 'Kttg: ',Kttg
#endif
      end

C-}}}
C-{{{ KttgS

      function KttgS(mt2,mh2,mst12,mst22,s2thetat,c2thetat
     &,mgl,At,muSUSY,beta,lam,GF)
      implicit none
      double complex KttgS,K1t,R2tA,R1tAS,R4tAS
      double precision x1,x2,Ysave
      double precision vev,vevS,Sqrt2
      double precision mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta,lam,GF
      double precision s2thetat,c2thetat

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2
      Ysave = lam*vev/Sqrt2

      KttgS = 2/3.d0*K1t(mt2,mh2)*Mgl/dsqrt(mt2)
     &     *(dsqrt(mt2)*Ysave/(mst12-mst22))
     &     *(x1/(1-x1)*dlog(x1)-x2/(1-x2)*dlog(x2))
     &  -dsqrt(mt2)/Mgl*s2thetat
     &     *R1tAS(mt2,mh2,mst12,mst22,mgl,muSUSY,lam,GF)
     &  + 2.d0*mt2*Ysave/Mgl/(mst12-mst22)
     &    *R2tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
     &  - K1t(mt2,mh2)/2.d0*mh2/(mst12-mst22)
     &     *R4tAS(mst12,mst22,mgl,muSUSY,lam,GF)

#if( debug > 0)
      write(*,*) 'Kttg: ',KttgS
#endif
      end

C-}}}
C-{{{ R1tA

      function R1tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
      implicit none
      double complex R1tA,K1t,BA
      double precision x1,x2,Yt,dli2
      double precision mst12,mst22,mgl,At,muSUSY,beta,mt2,mh2

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2
      Yt = At + muSUSY*dtan(beta)

      R1tA = 1/3.d0/(1-x1)**3 * (
     &     (1-x1**2+2*x1*dlog(x1))*
     &     (2*dlog(Mgl**2/mt2)-3-3/2.d0*K1t(mt2,mh2)+2*BA(mt2,mh2))
     &     - 8*x1*dli2(1-x1)-2*x1*(3+x1)*dlog(x1) )
     &   + 3/2.d0/(1-x1)**2 * (
     &     (1-x1+x1*dlog(x1))*(dlog(mt2/Mgl**2)
     &        +1+1/2.d0*K1t(mt2,mh2)-BA(mt2,mh2))
     &     + 2*x1*dli2(1-x1)+x1*(1+x1)*dlog(x1) )
     &   + 4/3.d0/(x1-x2)**2*Yt/Mgl * (1+1/2.d0*K1t(mt2,mh2)) * (
     &     x1**2*(1-2*x2)/2/(1-x1)/(1-x2)+x1/2/(1-x1)**2
     &     *(x1**2-2*x2+x1*x2)*dlog(x1) )
     &   - 1/3.d0/(1-x2)**3 * (
     &     (1-x2**2+2*x2*dlog(x2))*
     &     (2*dlog(Mgl**2/mt2)-3-3/2.d0*K1t(mt2,mh2)+2*BA(mt2,mh2))
     &     - 8*x2*dli2(1-x2)-2*x2*(3+x2)*dlog(x2) )
     &   - 3/2.d0/(1-x2)**2 * (
     &     (1-x2+x2*dlog(x2))*(dlog(mt2/Mgl**2)
     &         +1+1/2.d0*K1t(mt2,mh2)-BA(mt2,mh2))
     &     + 2*x2*dli2(1-x2)+x2*(1+x2)*dlog(x2) )
     &   - 4/3.d0/(x2-x1)**2*Yt/Mgl * (1+1/2.d0*K1t(mt2,mh2)) * (
     &     x2**2*(1-2*x1)/2/(1-x2)/(1-x1)+x2/2/(1-x2)**2
     &     *(x2**2-2*x1+x2*x1)*dlog(x2) )

      end

C-}}}
C-{{{ R1tAS

      function R1tAS(mt2,mh2,mst12,mst22,mgl,muSUSY,lam,GF)
      implicit none
      double complex R1tAS,K1t
      double precision x1,x2,Ysave,BA,dli2
      double precision vev,vevS,Sqrt2,mh2,mt2
      double precision mst12,mst22,mgl,muSUSY,lam,GF

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      Ysave = lam*vev/Sqrt2

      R1tAS = 4/3.d0/(x1-x2)**2*Ysave/Mgl * (1+1/2.d0*K1t(mt2,mh2)) * (
     &     x1**2*(1-2*x2)/2/(1-x1)/(1-x2)+x1/2/(1-x1)**2
     &     *(x1**2-2*x2+x1*x2)*dlog(x1) )
     &   - 4/3.d0/(x2-x1)**2*Ysave/Mgl * (1+1/2.d0*K1t(mt2,mh2)) * (
     &     x2**2*(1-2*x1)/2/(1-x2)/(1-x1)+x2/2/(1-x2)**2
     &     *(x2**2-2*x1+x2*x1)*dlog(x2) )

      end

C-}}}
C-{{{ R2tA

      function R2tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
      implicit none
      double complex R2tA,K1t
      double precision x1,x2,dli2
      double precision mst12,mst22,mgl,At,muSUSY,beta,mt2,mh2

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      R2tA = 1/3.d0/(1-x1)**3 * (
     &    2*(1-x1**2+2*x1*dlog(x1))*dlog(Mgl**2/mt2)-8*x1*dli2(1-x1)
     &    +(1-x1**2)*(1+1/2.d0*K1t(mt2,mh2))
     &  - 2*x1*(2+x1-1/2.d0*K1t(mt2,mh2))*dlog(x1))
     &  + 3/2.d0/(1-x1)**2 * ( 
     &    (1-x1+x1*dlog(x1))*dlog(mt2/Mgl**2)
     &    +2*x1*dli2(1-x1)+x1*(1+x1)*dlog(x1) )
     &  - 1/3.d0/(1-x2)**3 * (
     &    2*(1-x2**2+2*x2*dlog(x2))*dlog(Mgl**2/mt2)-8*x2*dli2(1-x2)
     &    +(1-x2**2)*(1+1/2.d0*K1t(mt2,mh2))
     &  - 2*x2*(2+x2-1/2.d0*K1t(mt2,mh2))*dlog(x2))
     &  - 3/2.d0/(1-x2)**2 * ( 
     &    (1-x2+x2*dlog(x2))*dlog(mt2/Mgl**2)
     &    +2*x2*dli2(1-x2)+x2*(1+x2)*dlog(x2) )
      end

C-}}}
C-{{{ R3tA

      function R3tA(mt2,mh2,mst12,mst22,mgl,At,muSUSY,beta)
      implicit none
      double complex R3tA,K1t,BA
      double precision x1,x2
      double precision mst12,mst22,mgl,At,muSUSY,beta,mt2,mh2

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      R3tA = 2/9.d0/(1-x1)**4 * (
     &     -2-3*x1+6*x1**2-x1**3-6*x1*dlog(x1) )
     &     * (2+K1t(mt2,mh2)-BA(mt2,mh2))
     &  + 3/8.d0/(1-x1)**3 * (
     &     1-x1**2+2*x1*dlog(x1) ) * (2+K1t(mt2,mh2)-2*BA(mt2,mh2))
     &  + 2/9.d0/(1-x2)**4 * (
     &     -2-3*x2+6*x2**2-x2**3-6*x2*dlog(x2) )
     &     * (2+K1t(mt2,mh2)-BA(mt2,mh2))
     &  + 3/8.d0/(1-x2)**3 * (
     &     1-x2**2+2*x2*dlog(x2) ) * (2+K1t(mt2,mh2)-2*BA(mt2,mh2))
      
      end

C-}}}
C-{{{ R4tA

      function R4tA(mst12,mst22,mgl,At,muSUSY,beta)
      implicit none
      double complex R4tA
      double precision x1,x2,Yt
      double precision mst12,mst22,mgl
      double precision At,muSUSY,beta

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2
      Yt = At + muSUSY*dtan(beta)

      R4tA = 4/3.d0/(x1-x2)**2*Yt/Mgl * (
     &     x1**2*(1-2*x2)/2.d0/(1-x1)/(1-x2)
     &   + x1/2.d0/(1-x1)**2*(x1**2-2*x2+x1*x2)*dlog(x1)
     &   - x2**2*(1-2*x1)/2.d0/(1-x2)/(1-x1)
     &   - x2/2.d0/(1-x2)**2*(x2**2-2*x1+x2*x1)*dlog(x2) )

      end

C-}}}
C-{{{ R4tAS

      function R4tAS(mst12,mst22,mgl,muSUSY,lam,GF)
      implicit none
      double complex R4tAS
      double precision x1,x2,Ysave
      double precision vev,vevS,Sqrt2
      double precision mst12,mst22,mgl,muSUSY,lam,GF

      x1 = mst12 / Mgl**2
      x2 = mst22 / Mgl**2

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      Ysave = lam*vev/Sqrt2

      R4tAS = 4/3.d0/(x1-x2)**2*Ysave/Mgl * (
     &     x1**2*(1-2*x2)/2.d0/(1-x1)/(1-x2)
     &   + x1/2.d0/(1-x1)**2*(x1**2-2*x2+x1*x2)*dlog(x1)
     &   - x2**2*(1-2*x1)/2.d0/(1-x2)/(1-x1)
     &   - x2/2.d0/(1-x2)**2*(x2**2-2*x1+x2*x1)*dlog(x2) )

      end

C-}}}
C-{{{ K1t

      function K1t(mt2,mh2)
      implicit none
      double precision mt2,mh2
      double complex K1t,ALOSMA

      K1t = -4/3.d0*mt2*ALOSMA(mh2,mt2)

      end

C-}}}
C-{{{ BA

      function BA(mt2,mh2)
      implicit none
      double precision mt,mt2,mh2
      double complex b0finim
      double complex BA

      mt = dsqrt(mt2)
      BA = b0finim(mh2,mt,mt,mt)

      end

C-}}}
C-{{{ Kttgeval

      function Kttgeval(mt2,mh2,mst12,mst22,mgl,
     &At,muSUSY,beta,s2thetat,c2thetat)
      implicit none
      double complex Kttgeval,FDSA
      double precision Yt,At,muSUSY,beta
      double precision mt2,mh2,mst12,mst22,mgl,s2thetat,c2thetat

      Yt = At + muSUSY*dtan(beta)

      Kttgeval = (s2thetat/2.d0-dsqrt(mt2)*Yt/(mst12-mst22))
     &     * (FDSA(Mgl**2,mt2,mst12) - FDSA(Mgl**2,mt2,mst22))

      end

C-}}}
C-{{{ FDSA

      function FDSA(m1,m2,m3)
      !squared masses to be entered
      implicit none
      double precision m1,m2,m3
      double complex FDSA,DeltaA,PhiA

      FDSA = 4/3.D0*dsqrt(m1)/dsqrt(m2)/DeltaA(m1,m2,m3)
     & *(m2*(m1-m2+m3)*dlog(m2/m1)+m3*(m1+m2-m3)*dlog(m3/m1)
     &  + 2*m1*m2*PhiA(m1,m2,m3))
     & + 3.D0*dsqrt(m2)/dsqrt(m1)/DeltaA(m1,m2,m3) 
     & *(m3*(m3-m2-m1)*dlog(m2/m1)+m3*(m2-m3-m1)*dlog(m3/m1)
     &  + m1*(m2+m3-m1)*PhiA(m1,m2,m3))

      end

C-}}}
C-{{{ DeltaA

      function DeltaA(m1,m2,m3)
      !squared masses to be entered - Kaehlen function
      implicit none
      double precision m1,m2,m3
      double complex DeltaA

      DeltaA = m2**2+m1**2+m3**2-2*(m1*m2+m1*m3+m2*m3)

      end

C-}}}
C-{{{ PhiDSA

      function PhiA(m1,m2,m3)
      !not the general definition - only for m1<m3 and m2<m3
      implicit none
      double precision m1,m2,m3,phiu,phiv,phixpl,phixmi,pi,dli2
      double complex PhiA,phila

      pi = 3.14159265358979323846264338328d0

      phiu = m1/m3
      phiv = m2/m3
      phila = dsqrt((1-phiu-phiv)**2-4*phiu*phiv)
      phixpl = 0.5D0*(1 + (phiu-phiv) - phila)
      phixmi = 0.5D0*(1 - (phiu-phiv) - phila)

      PhiA = 1/phila*(2*dlog(phixpl)*dlog(phixmi)-dlog(phiu)*dlog(phiv)
     & -2*(dli2(phixpl) + dli2(phixmi)) + pi**2/3)

      end

C-}}}

C-}}}

C-}}}
