! This file is part of SusHi.
! 
! This file contains various routines related to the renormalization
! of the squark sectors and the determination of the bottom Yukawa
! coupling including resummation.
!
C-{{{ subroutine renormalizeSM

      subroutine renormalizeSM
      !calculates the counterterms for mb and Yb in SM and 2HDM
      implicit none

      include '../commons/common-inputoutput.f'
      include '../commons/common-quark.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-consts.f'

      gb1 = 0.d0
      gb2 = 0.d0
      gt1 = 0.d0
      gt2 = 0.d0
      gte = 0.d0
      gte2 = 0.d0
      gbp1 = 0.d0
      gbp2 = 0.d0
      gbpe = 0.d0
      gbpe2 = 0.d0
      gtp1 = 0.d0
      gtp2 = 0.d0
      gtpe = 0.d0
      gtpe2 = 0.d0

      !if mbrunloop should be activated, just set dmb(3)!
      if ((mbrunloop.eq.0).and.(mbrunyuk.eq.0)) then
          mb = mbos
          mbyuk = mbos
          dmb = 0.d0
          dgb = 0.d0
          muranalytic = .true.
      else if ((mbrunloop.eq.0).and.(mbrunyuk.eq.1)) then
          mb = mbos
          mbyuk = mbmb
          dmb = 0.d0
          gb = mbyuk/mb * yukfac(3)
          !dgb = (mbmb - mbos)/ mbyuk * Pi/alphaggh
          dgb = - cf * alphaggh / Pi * Pi / alphaggh !1-loop
      else if ((mbrunloop.eq.0).and.(mbrunyuk.eq.2)) then
          mb = mbos
          mbyuk = mbMSbarmuR
          dmb = 0.d0
          gb = mbyuk/mb * yukfac(3)
          !dgb = (mbMSbarmuR - mbos)/ mbyuk * Pi/alphaggh
          dgb = (- cf + 2.d0 * dlog(mbmb/muRggh)) 
     &    * alphaggh / Pi * Pi / alphaggh !1-loop
      else
          write(*,105)
          write(*,*) "Renormalization scheme for bottom"
          write(*,*) "sector not known in the SM/2HDM."
          write(*,105)
          stop
      end if

      mb2 = mb**2

 105  format('#--------------------------------------------------#')

      end

C-}}}
C-{{{ subroutine renormalize

      subroutine renormalize()
      !calculates all counterterms relevant for gluon fusion in the MSSM
      implicit none
      double precision b0fin,dAb(3),dmbfin,alphasren,
     &deltaML,SUSHI_alphas,sthetaw,
     &dAbosfinmbHM,dAbosfintbHM,dAbosfinmbDS,dAbosfintbDS!,tmp
      double precision dmbstd,dthetabstd,dAbstd
      double precision dmbnontanb, dmbtanb, mbMSDRtrans
      double precision dmbtb, dmbntanbwosm, dthetabos(2)
      double precision delAb, ifunc
      double precision vev, vevS, Sqrt2
      logical expflag
      integer dep

      include '../commons/common-inputoutput.f'
      include '../commons/common-quark.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-consts.f'
      include '../commons/common-flags.f'

      !determine whether to use expanded CTs (true) or exact ones (false)
c      expflag = .true.
      expflag = .false.

      sthetaw = dsqrt(1.d0-(MW/MZ)**2)
      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam

!Various mb calculations

!Reconstruction of mb from FH - equals direct extraction according
!to s2thetab / (Ab - mu * tan(beta)) * (msb12 - msb22) / 2
      if (fhflag.eqv.(.true.)) then
         mbsb = dsqrt(cthetab**2 * msb12 + sthetab**2 * msb22 - dMLb
     &- M3SQ**2 - MZ**2 * dcos(2.d0*beta) * (-1/2.d0+1/3.d0*sthetaw**2))
      end if
!Reconstruction of mbsb
      if (onshellflag) then
         s2thetab = 2*sthetab*cthetab
         mbsb = s2thetab /(Ab - muSUSY*tan(beta))*(msb12 - msb22)/2.d0
         if (mbsb.le.0.d0) then
         call printerrorsushi(1,6
     &        ,'Your sbottom sector input is inconsistent!')
         end if
      thetab = dasin(sthetab)
      if (cthetab.lt.0.d0) then
      thetab = Pi/2.d0 - thetab
      end if
      c2thetab = cthetab**2 - sthetab**2
      s2thetab = 2*sthetab*cthetab
      t2thetab = s2thetab/c2thetab
      c2thetat = cthetat**2 - sthetat**2
      s2thetat = 2*sthetat*cthetat
      end if

!Calculation of mb_DRbar_SM
      !at scale muR
      alphasren = SUSHI_alphas(muRggh)
      mbDRSMmuR = mbMSDRtrans(mbMSbarmuR,1,alphasren) !transition to DRbar
      !at scale mT
      alphasren = SUSHI_alphas(mt)
      mbDRSMmt = mbMSDRtrans(mbMSbarmt,1,alphasren)
      !at scale mbmb
      alphasren = SUSHI_alphas(mbmb)
      mbDRSMmb = mbMSDRtrans(mbmb,1,alphasren)
      !at scale muB
      alphasren = SUSHI_alphas(muB)
      mbDRSMmuB = mbMSDRtrans(mbMSbarmuB,1,alphasren)

!Calculation of the tree-level-sbottom-masses:
      if (.not.onshellflag) then
      if (.not.fhflag) then
      call calcsquarkmass(M3SU,M3SQ,M3SD,mw,mz,mt,mbos
     &,muSUSY,beta,At,Ab,0.d0
     &,msb12,msb22,mst12,mst22
     &,cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)
      end if
      end if

!Calculation of delta_b with tree-level-sbottom-masses:
!Setting of muD - decoupling scale supersymmetry
      if (model.gt.0) then
         muD = (mgl + dsqrt(msb12) + dsqrt(msb22))/3.d0

         if (delmbfh.eq.2) then
         muD = 0.5 * (mgl + dsqrt(msb12) + dsqrt(msb22))/3.d0
         else if (delmbfh.eq.3) then
         muD = 2 * (mgl + dsqrt(msb12) + dsqrt(msb22))/3.d0
         end if

      else
         muD = muRggh
      end if
      !Calculation of mass m_b^MSbar(muD)
      call runmass(mbmb,SUSHI_alphas(mbmb)/Pi,SUSHI_alphas(muD)/Pi,nf,4,
     &mbMSbarmuD)
      mbDRSMmuD = mbMSDRtrans(mbMSbarmuD,1,SUSHI_alphas(muD))

      if (delmbfh.eq.1) then
      !Take delta_b from FH - includes ew corrections
         delta_mb = delmb
      !write(*,*) "delmb",delmb
      else if (delmbfh.eq.4) then
         delta_mb = delmb - 0.1*delmb
         delmb = delta_mb
      else if (delmbfh.eq.5) then
         delta_mb = delmb + 0.1*delmb
         delmb = delta_mb
      else
      !Calculate delta_b on our own
         alphasren = SUSHI_alphas(muD)
         delmb = ifunc(msb12,msb22,mgl**2)/(2.d0*Pi)
         delAb = - cf*alphasren*mgl*Ab*delmb
         delmb = cf*alphasren*mgl*muSUSY*tanb*delmb
         delta_mb = delmb
      !write(*,*) "delmb",delmb
      end if

!Calculation of running bottom mass at scale muB
      mbDRMSmuB = mbMSbarmuB / (1.d0 + delmb)

!NB: alphas/Pi * dmbtb/mb gives the same as delmb above
!Calculation of mb_DRbar_MSSM in accordance to PS
      !at scale mT
      alphasren = SUSHI_alphas(mt)
      dmbntanbwosm = dmbnontanb(.false.,.true.,mbDRSMmt
     &,mgl,msb12,msb22,Ab,mt)
      dmbtb = dmbtanb(mbDRSMmt,mgl,msb12,msb22,beta,muSUSY,mt)
      mbDRMSmt = mbDRSMmt
     & * (1.d0 + alphasren/Pi * dmbntanbwosm/mbDRSMmt)
     & / (1.d0 - alphasren/Pi * dmbtb/mbDRSMmt)
      !at scale muRggh
      alphasren = SUSHI_alphas(muRggh)
      dmbntanbwosm = dmbnontanb(.false.,.true.,mbDRSMmuR
     &,mgl,msb12,msb22,Ab,muRggh)
      dmbtb = dmbtanb(mbDRSMmuR,mgl,msb12,msb22,beta,muSUSY
     &     ,muRggh)
      mbDRMSmuR = mbDRSMmuR
     & * (1.d0 + alphasren/Pi * dmbntanbwosm/mbDRSMmuR)
     & / (1.d0 - alphasren/Pi * dmbtb/mbDRSMmuR)
      !at scale mbmb
      alphasren = SUSHI_alphas(mbmb)
      dmbntanbwosm = dmbnontanb(.false.,.true.,mbDRSMmb
     &,mgl,msb12,msb22,Ab,mbmb)
      dmbtb = dmbtanb(mbDRSMmb,mgl,msb12,msb22,beta,muSUSY,mbmb)
      mbDRMSmb = mbDRSMmb
     & * (1.d0 + alphasren/Pi * dmbntanbwosm/mbDRSMmb)
     & / (1.d0 - alphasren/Pi * dmbtb/mbDRSMmb)

!mbOSsbot calculation similar to PS
!without any Yukawa corrections, just shifting mbDRMSmt to on-shell value
!      alphasren = SUSHI_alphas(mt)
!      mbOSsbot = mbDRMSmt -dmbfin(.true.,.true.,mbDRMSmt
!     &,mgl,msb12,msb22,s2thetab,mt)*alphasren/Pi

!Calculation of the OS-sbottom-masses:
      if ((fhflag.eqv.(.false.)).and.(onshellflag.eqv.(.false.))) then
         alphasren = SUSHI_alphas(muD)
         dMLb = alphasren/Pi*deltaML(mt,dsqrt(mst12),dsqrt(mst22),
     &        cthetat,sthetat,
     &        mbos,dsqrt(msb12),dsqrt(msb22),cthetab,sthetab,mgl,muD) 
         mb = mbos
         mb2 = mb**2
         mt2 = mt**2
         call sbottomrenormHM(mbsb,mb,mgl,msb12,msb22,sthetab,cthetab
     &        ,dmbsbos,dthetabos,dmsb1,dmsb2,muD)

         if(dabs(s2thetab).gt.1.d-8) then
            mbsb = mbDRSMmuD - ( 2.d0*mb/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        - (c2thetab*(msb12-msb22)*dcos(beta)*dsin(beta))/muSUSY
     &        * dthetabos(1)
     &        - 2.d0*mb2/(msb12-msb22)/s2thetab
     &        * dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &    sthetab,cthetab,Ab,muSUSY,dmsb1(1),dmsb2(1),0.d0,muD))
     &        * alphasren/Pi
         else
c     limit for sin(2*thetab) -> 0
            mbsb = mbDRSMmuD - ( 2.d0*mb/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        + (mb*((-6*msb12*dmsb1(1) +
     -        6*msb22*dmsb2(1) -
     -        8*msb12 - 4*lb1*msb12 + 8*msb22 + 
     -        4*lb2*msb22)*muSUSY + 
     -        2*b0fin(msb22,mb,mgl,muD)*
     -        ((mb2 + mgl**2 - msb22)*muSUSY + 
     -        mgl*(msb12 - msb22)*dcos(beta)*dsin(beta)) + 
     -        b0fin(msb12,mb,mgl,muD)*
     -        (-2*(mb2 + mgl**2 - msb12)*muSUSY + mgl*(msb12 - msb22)
     -        *dsin(2*beta))))/(3.d0*(msb12 - msb22)*muSUSY)
     -        ) * alphasren/Pi
            if(thetab.le.1.d0) then
c     thetab -> 0
               mbsb =  mbsb + ((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetabos(1) * alphasren/Pi
            else
c     thetab -> Pi/2
               mbsb =  mbsb - ((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetabos(1) * alphasren/Pi
            endif
         endif
      endif

      if (.not.onshellflag) then
      call calcsquarkmass(M3SU,M3SQ,M3SD,mw,mz,mt,mbsb
     &,muSUSY,beta,At,Ab,dMLb
     &,msb12,msb22,mst12,mst22
     &,cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)
      end if

!Various mb MSbar and DRbar values:
!      write(*,*)
!      write(*,*) 'Various internal (s)bottom masses/parameters:'
!      write(*,*) 'm_b^MSbar(mb):         ',mbmb
!      write(*,*) 'm_b^MSbar(muR):        ',mbMSbarmuR
!      write(*,*) 'm_b^OS:                ',mbOS
!      write(*,*) 'm_b in sbottom sector: ',mbsb
!      write(*,*) 'Delta_b:               ',delmb
!      write(*,*)
!End of mb/sbottom calculations

!Calculation of counterterms in mbsb, Ab and thetab.
!In case of on-shell input we dont allow for the various options,
!but mbsb is always considered to be a dependent parameter.
      mb = mbsb
      mb2 = mb**2

      if(expflag.eqv..true.) then
         call sbottomrenormDS(mbsb,thetab,mgl,msb12,msb22,dmbsbos
     &        ,dthetabos,dmsb1,dmsb2,muRggh)
      else
         call sbottomrenormHM(mbsb,mb,mgl,msb12,msb22,sthetab,cthetab
     &        ,dmbsbos,dthetabos,dmsb1,dmsb2,muRggh)
      endif

      dep = 3
      dAb = 0.d0
      if(mbsbch.eq.0) then      !mb OS
         dmbsb = dmbsbos
      else if(mbsbch.eq.1) then !mb DRbar
         dmbsb = 0.d0
      else
         dep = 0
      endif
      if(thetabch.eq.0) then    !thetab OS
         dthetab(1) = 0.d0
         dthetab(2) = dthetabos(1)
         dthetab(3) = dthetabos(2)
      else if(thetabch.eq.1) then !thetab DRbar
         dthetab = 0.d0
      else
         if(dep.eq.0) then
            write(*,*) 'error: two dependent quantities.'
            stop
         endif
         dep = 1
      endif
      if(Abch.eq.0) then        !Ab OS
         if(dep.eq.0) then
            if(expflag.eqv..true.) then
               dAb(3) = dAbosfintbDS(mgl,msb12,msb22
     &,thetab,Ab,muSUSY,beta,dthetabos(1),dthetab(2),muRggh)
            else
               dAb(3) = dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &sthetab,cthetab,Ab,muSUSY,dmsb1(1),dmsb2(1),dthetab(2),muRggh)
            endif
         else if(dep.eq.1) then
            if(expflag.eqv..true.) then
               dAb(3) = dAbosfinmbDS(beta,mb,mgl,msb12,msb22,Ab,muSUSY
     &,dmbsbos(1),dmbsbos(2),dmbsb(2),muRggh)
            else
               dAb(3) = dAbosfinmbHM(beta,mb,mgl,msb12
     &,msb22,Ab,muSUSY,dmbsb(2),muRggh)
            endif
         else
            write(*,*) 'error: no dependent quantity.'
            stop
         endif
      else if(Abch.eq.1) then   !Ab DRbar
         dAb(3) = 0.d0
      else
         if(dep.lt.2) then
            write(*,*) 'error: two dependent quantities.'
            stop
         endif
         dep = 2
      endif

      if(dep.eq.0) then
      !mb dependent
         dmbsb(1) = 0.d0
         dmbsb(2) = 2.d0*mb/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
         if(expflag.eqv..true.) then
            dmbsb(2) = dmbsb(2) + 2.d0/t2thetab*mb * dthetab(2)
         elseif (Abch.ne.0) then
            dmbsb(2) = dmbsb(2) + 2.d0/t2thetab*mb * dthetab(2)
     &           - 2.d0*mb2/(msb12-msb22)/s2thetab * dAb(3)
         elseif(dabs(s2thetab).gt.1.d-8) then
            dmbsb(2) = dmbsb(2) 
     &        - (c2thetab*(msb12-msb22)*dcos(beta)*dsin(beta))/muSUSY
     &        * dthetab(2)
     &        - 2.d0*mb2/(msb12-msb22)/s2thetab
     &        * dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &sthetab,cthetab,Ab,muSUSY,dmsb1(1),dmsb2(1),0.d0,muRggh)
         else
c     limit for sin(2*thetab) -> 0
            dmbsb(2) = dmbsb(2)
     -        + (mb*((-6*msb12*dmsb1(1) +
     -        6*msb22*dmsb2(1) -
     -        8*msb12 - 4*lb1*msb12 + 8*msb22 + 
     -        4*lb2*msb22)*muSUSY + 
     -        2*b0fin(msb22,mb,mgl,muRggh)*
     -        ((mb2 + mgl**2 - msb22)*muSUSY + 
     -        mgl*(msb12 - msb22)*dcos(beta)*dsin(beta)) + 
     -        b0fin(msb12,mb,mgl,muRggh)*
     -        (-2*(mb2 + mgl**2 - msb12)*muSUSY + mgl*(msb12 - msb22)
     -        *dsin(2*beta))))/(3.d0*(msb12 - msb22)*muSUSY)
            if(thetab.le.1.d0) then
c     thetab -> 0
               dmbsb(2) = dmbsb(2)-((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetab(2)
            else
c     thetab -> Pi/2
               dmbsb(2) = dmbsb(2)+((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetab(2)
            endif
         endif
         dmbsb(3) = 0.d0
      else if(dep.eq.1) then
      !thetab dependent
      !not offered as option
         write(*,*) 'error: theta_b cannot be chosen',
     &' as a dependent quantity. '
         stop
         dthetab(1) = t2thetab/mb/2.d0 * dmbsb(1)
         dthetab(2) = -t2thetab/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        + t2thetab/mb/2.d0 * dmbsb(2)
         if(expflag.eqv..false.) then
            dthetab(2) = dthetab(2) + mb/(msb12-msb22)/c2thetab * dAb(3)
         else
            dthetab(3) = 0.d0
         endif
      else if(dep.eq.2) then
      !Ab dependent

         dAb(1) = - s2thetab*(msb12-msb22)/2.d0/mb2 * dmbsb(1)

         dAb(2) = s2thetab/mb
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        - s2thetab*(msb12-msb22)/2.d0/mb2 * dmbsb(2)
     &        + c2thetab*(msb12-msb22)/mb * dthetab(2)

         dAb(3) = s2thetab/mb
     &        * (dmsb1(2)*msb12 - dmsb2(2)*msb22)
     &        - s2thetab*(msb12-msb22)/2.d0/mb2 * dmbsb(3)
     &        + c2thetab*(msb12-msb22)/mb * dthetab(3)

      else
         write(*,*) 'error: no dependent quantity.'
         stop
      endif

      dthetabstd = dthetabos(1)
      dmbstd = 2.d0*mb/(msb12-msb22)
     &     * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
      if(expflag.eqv..true.) then
         dmbstd = dmbstd + 2.d0/t2thetab*mb * dthetabstd
         dAbstd = dAbosfintbDS(mgl,msb12,msb22,
     &thetab,Ab,muSUSY,beta,dthetabos(1),dthetabstd,muRggh)
      else
         dAbstd = dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &sthetab,cthetab,Ab,muSUSY,dmsb1(1),dmsb2(1),dthetabstd,muRggh)
         if(dabs(s2thetab).gt.1.d-8) then
            dmbstd = dmbstd 
     &        - (c2thetab*(msb12-msb22)*dcos(beta)*dsin(beta))/muSUSY
     &        * dthetabstd
     &        - 2.d0*mb2/(msb12-msb22)/s2thetab
     &        * dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &sthetab,cthetab,Ab,muSUSY,dmsb1(1),dmsb2(1),0.d0,muRggh)
         else
c     limit for sin(2*thetab) -> 0
            dmbstd = dmbstd
     -        + (mb*((-6*msb12*dmsb1(1) +
     -        6*msb22*dmsb2(1) -
     -        8*msb12 - 4*lb1*msb12 + 8*msb22 + 
     -        4*lb2*msb22)*muSUSY + 
     -        2*b0fin(msb22,mb,mgl,muRggh)*
     -        ((mb2 + mgl**2 - msb22)*muSUSY + 
     -        mgl*(msb12 - msb22)*dcos(beta)*dsin(beta)) + 
     -        b0fin(msb12,mb,mgl,muRggh)*
     -        (-2*(mb2 + mgl**2 - msb12)*muSUSY + mgl*(msb12 - msb22)
     -        *dsin(2*beta))))/(3.d0*(msb12 - msb22)*muSUSY)
            if(thetab.le.1.d0) then
c     thetab -> 0
               dmbstd = dmbstd - ((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetabstd
            else
c     thetab -> Pi/2
               dmbstd = dmbstd + ((msb12-msb22)*dcos(beta)*dsin(beta))
     &   / muSUSY * dthetabstd
            endif
         endif
      endif

c$$$      write(*,*) mgl,Ab,dsqrt(msb12),dsqrt(msb22),mh,mb,thetab,
c$$$     &     muSUSY,beta,dasin(sthetaw),mz,alpha,muRggh
c$$$      write(*,*) 'mb_OS: ',dmbsbos(1),dmbsbos(1)+dmbsbos(2),
c$$$     &     dmbsbos(1)+dmbsbos(2)+dmbsbos(3)
c$$$      write(*,*) 'tb_OS: ',dthetabos(1),dthetabos(1)+dthetabos(2)
c$$$      write(*,*) 'msb_OS: ',dmsb1(1),dmsb1(1)+dmsb1(2),
c$$$     &     dmsb2(1),dmsb2(1)+dmsb2(2)
c$$$      write(*,*) 'mb_RS: ',dmbsb(1),dmbsb(1)+dmbsb(2),
c$$$     &     dmbsb(1)+dmbsb(2)+dmbsb(3)
c$$$      write(*,*) 'tb_RS: ',dthetab(1),dthetab(1)+dthetab(2),
c$$$     &     dthetab(1)+dthetab(2)+dthetab(3)
c$$$      write(*,*) 'Ab_RS: ',dAb(1),dAb(1)+dAb(2),dAb(1)+dAb(2)+dAb(3)
c$$$      write(*,*) 'delstd:',dmbstd, dthetabstd, dAbstd

!CT in Ab and thetab
      Ab = Ab + alphaggh/Pi * (dAbstd - (dAb(1)+dAb(2)+dAb(3)))
      thetab = thetab+alphaggh/Pi*(dthetabstd-(dthetab(1)+dthetab(2)))
      mbsb = mbsb + alphaggh/Pi * (dmbstd - (dmbsb(1)+dmbsb(2)))

      mbsb2 = mbsb**2
      dmbsb = dmbsb / mbsb

!Settings for internal bottom mass/Yukawa coupling
      call sbottomrenormDS(mbos,thetab,mgl,msb12,msb22,dmbsbos
     &        ,dthetabos,dmsb1,dmsb2,muRggh)

      if(mbrunloop.eq.0) then
      !mb_OS:
         mb = mbos
         dmb(1) = 0.d0
         dmb(2) = 0.d0
         dmb(3) = 0.d0
      elseif(mbrunloop.eq.1) then
      !mb_DRbar_MSSM(muRggh):
         mb = mbDRMSmuR
         dmb(1) = dmbsbos(1)
         dmb(3) = -mbos*(5/3.d0+dlog(muRggh**2/mbos**2))
         dmb(2) = dmbsbos(2)-dmb(3)
      elseif(mbrunloop.eq.2) then
      !mb_DRbar_MSSM(mb):
         mb = mbDRMSmb
         dmb(1) = dmbsbos(1)
         dmb(3) = -mbos*(5/3.d0+dlog(mbmb**2/mbos**2))
         dmb(2) = dmbsbos(2) + mbos*2/3.d0*dlog(muRggh**2/mbmb**2)
     &        -dmb(3)
      elseif(mbrunloop.eq.3) then
      !mb_DRbar_MSSM(mub):
         mb = mbDRMSmub
         dmb(1) = dmbsbos(1)
         dmb(3) = -mbos*(5/3.d0+dlog(mub**2/mbos**2))
         dmb(2) = dmbsbos(2) + mbos*2/3.d0*dlog(muRggh**2/mub**2)-dmb(3)
      else
         stop
      endif
      mb2 = mb**2
      dmb = dmb / mbos

      mbyuk = mbos
      dgb = 0.d0
      !NOTE: dgb only contains the SM part, the SUSY part is in virtds.
      !The sign of dgb is different from the sign in virtds.
      if ((mbrunyuk.eq.1).or.(mbrunyuk.eq.2)) then
      !difference between 1 and 2: 1 muB=mbmb, 2 muB = muRggh
         mbyuk = mbDRMSmuB
         dgb = 2.d0 * dlog(muD/muB)
     &   + 2.d0 * dlog(mb/muD) - 5/3.d0 + 1/3.d0
!          dgb = 2.d0 * dlog(muD/muB) * alphasPDF(muB)/alphasPDF(muRggh)
!     &   + 2.d0 * dlog(mb/muD) - 5/3.d0
!     &   + 1/3.d0 * alphasPDF(muB)/alphasPDF(muRggh)
         gb = mbDRMSmuB / mbos * yukfac(3) 
         tanbresum = 0
      else if(mbrunyuk.eq.3) then
         mbyuk = mbDRMSmuR
         dgb = (dmbfin(.true.,.true.,mbos,mgl,msb12,msb22,s2thetab
     &        ,muRggh)- dmbfin(.false.,.true.,mbos,mgl,msb12,msb22
     &        ,s2thetab,muRggh))/ mbDRMSmuR * mbDRMSmuR / mbos !last factor: Kfac
         gb = mbDRMSmuR / mbos * yukfac(3) 
         tanbresum = 0
      else if(mbrunyuk.eq.4) then
         mbyuk = mbDRMSmb
         dgb = (dmbfin(.true.,.true.,mbos,mgl,msb12,msb22,s2thetab,mbmb) 
     &      - dmbfin(.false.,.true.,mbos,mgl,msb12,msb22,s2thetab,mbmb))
     & / mbDRMSmb * mbDRMSmb / mbos !last factor: Kfac
         gb = mbDRMSmb / mbos * yukfac(3) 
         tanbresum = 0
      else
      !Delta_b resummation:
         if(tanbresum.eq.1) then
            gb = yukfac(3) / ( 1.d0 + delmb)
         else if(tanbresum.eq.2) then   
            if (pseudo.eq.1) then
               gb = yukfac(3) * (1.d0 - (1.d0/tanb**2)*delmb
     &           - Amx(Sind,3)*vev/(Amx(Sind,2)*vevS*tanb)*delmb)
     &              * (1.d0/(1.d0+delmb))
            else if (pseudo.eq.0) then
               gb = yukfac(3) * (1.d0/(1.d0+delmb))
     &         * (1.d0 + delmb * Hmx(Sind,2)/(tanb*Hmx(Sind,1))
     &       + delmb*Hmx(Sind,3)*vev*Cos(Atan(tanb))/(Hmx(Sind,1)*vevS))
            end if 
         else
            gb = yukfac(3)
         end if
      end if

      !Final settings of renormalize
      cthetab = dcos(thetab)
      sthetab = dsin(thetab)
      t2thetab = dtan(2*thetab)
      c2thetab = cthetab**2 - sthetab**2
      s2thetab = 2*sthetab*cthetab

!Screen output of relevant data in the (s)bottom sector:
!      write(*,103) 'Bottom masses used for:              '
!      write(*,101) 'internal masses:      ',mb
!      if (mbrunyuk.gt.0) then
!      write(*,101) 'bottom Yuk:           ',mbyuk
!      else
!      write(*,101) 'bottom Yuk:           ',mbyuk*gb/yukfac(3)
!      end if
!      if ((model.eq.1).or.(model.eq.3)) then
!      write(*,101) 'sbottom sector:       ',mbsb
!      write(*,105)
!      write(*,103) 'Renormalized sbottom sector:         '
!      write(*,101) 'm_b:                   ',mbsb
!      write(*,102) 'sin(theta_b):          ',sthetab
!      write(*,102) 'cos(theta_b):          ',cthetab
!      write(*,101) 'A_b:                   ',Ab
!      write(*,101) 'OS-sbottom mass 1:    ',dsqrt(msb12)
!      write(*,101) 'OS-sbottom mass 2:    ',dsqrt(msb22)
!      end if
!      write(*,105)

!Recalculation of counterterms using the new parameters
      call sbottomrenormDS(mbsb,thetab,mgl,msb12,msb22,dmbsbos,dthetabos
     &        ,dmsb1,dmsb2,muRggh)

      dep = 3
      dAb = 0.d0
      if(mbsbch.eq.0) then      !mb OS
         dmbsb = dmbsbos
      else if(mbsbch.eq.1) then !mb DRbar
         dmbsb = 0.d0
      else
         dep = 0
      endif
      if(thetabch.eq.0) then    !thetab OS
         dthetab(1) = 0.d0
         dthetab(2) = dthetabos(1)
         dthetab(3) = dthetabos(2)
      else if(thetabch.eq.1) then !thetab DRbar
         dthetab = 0.d0
      else
         dep = 1
      endif
      if(Abch.eq.0) then        !Ab OS
         if(dep.eq.0) then
            dAb(3) = dAbosfintbDS(mgl,msb12,msb22,
     &thetab,Ab,muSUSY,beta,dthetabos(1),dthetab(2),muRggh)
         else if(dep.eq.1) then
            dAb(3) = dAbosfinmbDS(beta,mb,mgl,msb12,msb22,Ab,muSUSY
     &,dmbsbos(1),dmbsbos(2),dmbsb(2),muRggh)
         endif
      else if(Abch.eq.1) then   !Ab DRbar
         dAb(3) = 0.d0
      else
         dep = 2
      endif

      if(dep.eq.0) then
!     mb dependent
         dmbsb(1) = 0.d0
         dmbsb(2) = 2.d0*mbsb/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        + 2.d0/t2thetab*mbsb * dthetab(2)
         dmbsb(3) = 0.d0
      else if(dep.eq.1) then
!     thetab dependent
         dthetab(1) = t2thetab/mbsb/2.d0 * dmbsb(1)
         dthetab(2) = -t2thetab/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        + t2thetab/mbsb/2.d0 * dmbsb(2)
         dthetab(3) = 0.d0
      else if(dep.eq.2) then
!     Ab dependent

         dAb(1) = - s2thetab*(msb12-msb22)/2.d0/mbsb2 * dmbsb(1)

         dAb(2) = s2thetab/mbsb
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        - s2thetab*(msb12-msb22)/2.d0/mbsb2 * dmbsb(2)
     &        + c2thetab*(msb12-msb22)/mbsb * dthetab(2)

         dAb(3) = s2thetab/mbsb
     &        * (dmsb1(2)*msb12 - dmsb2(2)*msb22)
     &        - s2thetab*(msb12-msb22)/2.d0/mbsb2 * dmbsb(3)
     &        + c2thetab*(msb12-msb22)/mbsb * dthetab(3)
      endif

      dAbds = dAb
      dmbsb = dmbsb / mbsb

      if ((mbrunloop.eq.0).and.(mbrunyuk.eq.0)) then
       muranalytic = .true.
      end if
      if ((mbsbch+thetabch+Abch).gt.2) then
       muranalytic = .false.
      end if

! 100  format(a45)
! 103  format(a30)
! 101  format(a22,'    ',f10.5, '   GeV')
! 102  format(a22,'    ',f10.5)
! 105  format('#--------------------------------------------------#')
      end

C-}}}
C-{{{ subroutine renormalizesimpleFH

      subroutine renormalizesimpleFH(M3SU,M3SQ,M3SD
     &,msb12,msb22,mst12,mst22,dMLb,beta,mw,mz,mgl,yukfac,muRggh
     &,mbos,mbsb,mbsb2,delmb,delta_mb,muSUSY,Ab,At,tanbresum
     &,gb,dthetab,dmbsb,dAb,dmb,dmsb1,dmsb2,dgb
     &,mb,mb2,mt,mbyuk
     &,cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)
      !calculates all counterterms relevant for gluon fusion in the MSSM 
      !minimal version of the above routine - Feb 2016 S. Liebler
      implicit none
      !input parameters:
      double precision msb12,msb22,mst12,mst22,dMLb
      double precision beta,mw,mz,mgl,yukfac(9),muRggh
      double precision mbos,delmb,muSUSY,Ab,At,mt
      double precision M3SU,M3SQ,M3SD
      integer tanbresum
      !output parameters:
      double precision gb,dthetab(3),dmbsb(3),dAb(3),dmb(3),mbsb
      double precision mbsb2,dgb,delta_mb
      double precision dmsb1(2),dmsb2(2),mb,mb2,mbyuk
      double precision cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab
      !internal parameters:
      double precision sthetaw,dthetabos(2),dmbsbos(3)
      double precision pi,dAbosfintbDS
      integer dep

      pi = 3.14159265358979323846264338328d0

      sthetaw = dsqrt(1.d0-(mw/mz)**2)

      !Reconstruction of mb from FH - equals direct extraction according
      !to s2thetab / (Ab - mu * tan(beta)) * (msb12 - msb22) / 2
      mbsb = dsqrt(cthetab**2 * msb12 + sthetab**2 * msb22 - dMLb
     &- M3SQ**2 - MZ**2 * dcos(2.d0*beta) * (-1/2.d0+1/3.d0*sthetaw**2))
      thetab = dasin(sthetab)
      if (cthetab.lt.0.d0) then
      thetab = Pi/2.d0 - thetab
      end if
      c2thetab = cthetab**2 - sthetab**2
      s2thetab = 2*sthetab*cthetab
      t2thetab = dtan(2*thetab)

      call calcsquarkmass(M3SU,M3SQ,M3SD,mw,mz,mt,mbsb
     &,muSUSY,beta,At,Ab,dMLb
     &,msb12,msb22,mst12,mst22
     &,cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)

      call sbottomrenormDS(mbsb,thetab,mgl,msb12,msb22,dmbsbos
     &        ,dthetabos,dmsb1,dmsb2,muRggh)
  
      !a few standard settings:
      dmb = 0.d0 !OS counterterm for bottom mass
      mb = mbos
      mb2 = mb**2
      mbyuk = mbos
      dgb = 0.d0
      delta_mb = delmb

      !Delta_b resummation:
      if(tanbresum.eq.1) then
         gb = yukfac(3) / ( 1.d0 + delmb)
      else
         gb = yukfac(3)
      end if

      dAb = 0.d0
      dep = 0 !only mbsb is accepted to be dependent in this routine
      dthetab(1) = 0.d0
      dthetab(2) = dthetabos(1)
      dthetab(3) = dthetabos(2)
      dAb = 0.d0
      dAb(3) = dAbosfintbDS(mgl,msb12,msb22,
     &thetab,Ab,muSUSY,beta,dthetabos(1),dthetab(2),muRggh)

      if(dep.eq.0) then
         dmbsb(1) = 0.d0
         dmbsb(2) = 2.d0*mbsb/(msb12-msb22)
     &        * (dmsb1(1)*msb12 - dmsb2(1)*msb22)
     &        + 2.d0/t2thetab*mbsb * dthetab(2)
         dmbsb(3) = 0.d0
      end if

      dmbsb = dmbsb / mbsb
      mbsb2 = mbsb*mbsb

      end

C-}}}
C-{{{ subroutine calcsquarkmass

      subroutine calcsquarkmass(M3SU,M3SQ,M3SD,mw,mz,mt,mb
     &,muSUSY,beta,At,Ab,dMLb
     &,msb12,msb22,mst12,mst22
     &,cthetat,sthetat,c2thetat,s2thetat
     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)
      implicit none
      !input parameters 
      double precision M3SU,M3SQ,M3SD,mw,mz,mt,mb,muSUSY,beta,At,Ab
      double precision dMLb
      !output parameters
      double precision msb12,msb22,mst12,mst22
      double precision cthetat,sthetat,c2thetat,s2thetat
      double precision cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab
      !internal parameters
      double precision tanb,pi
      double precision sthetaW
      double precision MsbotEig(2,2), MstopEig(2,2)
      double precision M3SU2, M3SQ2, M3SD2, thetat
      double precision msbot1,msbot2,mstop1,mstop2,mb2,mt2
c      double precision StopMix(2,2), SbotMix(2,2)
c      double precision TestMatrix(2,2)

      pi = 3.14159265358979323846264338328d0

      sthetaw = dsqrt(1.d0-(MW/MZ)**2)
      mb2 = mb**2
      mt2 = mt**2
      tanb = dtan(beta)

      M3SU2 = M3SU**2
      M3SQ2 = M3SQ**2
      M3SD2 = M3SD**2
      if (M3SU.lt.0.d0) M3SU2 = - M3SU2
      if (M3SQ.lt.0.d0) M3SQ2 = - M3SQ2
      if (M3SD.lt.0.d0) M3SD2 = - M3SD2

      MstopEig(1,1) = M3SQ2 + mt2 + MZ**2 
     &* dcos(2.d0*beta) * (1/2.d0-2/3.d0*sthetaw**2)
      MstopEig(2,2) = M3SU2 + mt2 + MZ**2 
     &* dcos(2.d0*beta) * (2/3.d0*sthetaw**2)
      MstopEig(1,2) = mt * (At - muSUSY/tanb)
      MstopEig(2,1) = MstopEig(1,2)

      !calculate mstop1 and mstop2 with mstop1 < mstop2
      mstop1=dsqrt(M3SQ2/2.d0 + M3SU2/2.d0 + mt2
     &       +MZ**2/4.d0*dcos(2.d0*beta)
     &     -dsqrt((M3SQ2 - M3SU2 + MZ**2*dcos(2.d0*beta)
     &       *(1/2.d0-4/3.d0*sthetaw**2))**2
     &     + 4.d0*mt2*(At-muSUSY/tanb)**2 )/2.d0)
      mstop2=dsqrt(M3SQ2/2.d0 + M3SU2/2.d0 + mt2
     &       +MZ**2/4.d0*dcos(2.d0*beta)
     &     +dsqrt((M3SQ2 - M3SU2 + MZ**2 *dcos(2.d0*beta)
     &       *(1/2.d0-4/3.d0*sthetaw**2))**2
     &     + 4.d0*mt2*(At-muSUSY/tanb)**2 )/2.d0)
      mst12 = mstop1**2
      mst22 = mstop2**2

      thetat = dasin(2.d0*mt/(mst12-mst22)*(At-muSUSY/tanb))/2.d0
      if (((2.d0*mt/(mst12-mst22)*(At-muSUSY/tanb)).gt.(1.d0)).and.
     & ((2.d0*mt/(mst12-mst22)*(At-muSUSY/tanb)).lt.(1.0001d0))) then
         thetat = dasin(1.d0)/2.d0
      end if
      if (((2.d0*mt/(mst12-mst22)*(At-muSUSY/tanb)).lt.(-1.d0)).and.
     & ((2.d0*mt/(mst12-mst22)*(At-muSUSY/tanb)).gt.(-1.0001d0))) then
         thetat = dasin(-1.d0)/2.d0
      end if
      sthetat = dsin(thetat)
      cthetat = dcos(thetat)

      !change sin <-> cos
      if(dabs(cthetat**2 *mst12 +sthetat**2 *mst22 
     & - MstopEig(1,1)).gt.(1.d-3)) then
         thetat = Pi/2.d0 - thetat
         sthetat = dsin(thetat)
         cthetat = dcos(thetat)
      end if

      c2thetat = cthetat**2 - sthetat**2
      s2thetat = 2*sthetat*cthetat

!      write(*,*) "---------------------------Test diagonalization"
!      write(*,*) "thetat",thetat,sthetat,cthetat
!      write(*,*) "11",cthetat**2 * mst12 + sthetat**2 * mst22
!     & ,MstopEig(1,1)
!      write(*,*) "12",cthetat*sthetat* (mst12 - mst22)
!     & ,MstopEig(1,2)
!      write(*,*) "22",sthetat**2 * mst12 + cthetat**2 * mst22
!     & ,MstopEig(2,2)
!      StopMix(1,1) = cthetat
!      StopMix(1,2) = sthetat
!      StopMix(2,1) = -sthetat
!      StopMix(2,2) = cthetat
!      TestMatrix = Matmul(Matmul(StopMix,MstopEig),Transpose(StopMix))
!      write(*,*) "Matrix-1",Sqrt(TestMatrix(1,1)),TestMatrix(1,2)
!      write(*,*) "Matrix-2",TestMatrix(2,1),Sqrt(TestMatrix(2,2))
!      write(*,*) "Masses",mstop1,mstop2

      MsbotEig(1,1) = M3SQ2 + dMLb + mb2 + MZ**2 
     &* dcos(2.d0*beta) * (-1/2.d0+1/3.d0*sthetaw**2)
!      MsbotEig(2,2) = M3SD2 + mb2 + MZ**2 
!     &* dcos(2.d0*beta) * (-1/3.d0*sthetaw**2)
!      MsbotEig(1,2) = mb * (Ab - muSUSY*tanb)
!      MsbotEig(2,1) = MsbotEig(1,2)

      !calculate msbot1 and msbot2 with msbot1 < msbot2
      msbot1=dsqrt(M3SQ2/2.d0 + M3SD2/2.d0 + dMLb/2.d0 + mb2
     &     -MZ**2/4.d0*dcos(2*beta)-dsqrt(
     &     (M3SQ2 - M3SD2 + dMLb
     &     + MZ**2*dcos(2*beta)*(-1/2.d0+2/3.d0*sthetaw**2))**2
     &     + 4.d0*mb2*(Ab-muSUSY*tanb)**2 )/2.d0)
      msbot2=dsqrt(M3SQ2/2.d0 + M3SD2/2.d0 + dMLb/2.d0 + mb2
     &     -MZ**2/4.d0*dcos(2*beta)+dsqrt(
     &     (M3SQ2 - M3SD2 + dMLb
     &     + MZ**2*dcos(2*beta)*(-1/2.d0+2/3.d0*sthetaw**2))**2
     &     + 4.d0*mb2*(Ab-muSUSY*tanb)**2 )/2.d0)

      msb12 = msbot1**2
      msb22 = msbot2**2

      thetab = dasin(2.d0*mb/(msb12-msb22)*(Ab-muSUSY*tanb))/2.d0
      if (((2.d0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)).gt.(1.d0)).and.
     & ((2.d0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)).lt.(1.0001d0))) then
         thetab = dasin(1.d0)/2.d0
      end if
      if (((2.d0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)).lt.(-1.d0)).and.
     & ((2.d0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)).gt.(-1.0001d0))) then
         thetab = dasin(-1.d0)/2.d0
      end if

      sthetab = dsin(thetab)
      cthetab = dcos(thetab)

      if(dabs(cthetab**2 * msb12 + sthetab**2 * msb22 
     & - MsbotEig(1,1)).gt.(1.d-3)) then
         thetab = Pi/2.d0 - thetab
         sthetab = dsin(thetab)
         cthetab = dcos(thetab)
      end if

      t2thetab = dtan(2*thetab)
      c2thetab = cthetab**2 - sthetab**2
      s2thetab = 2*sthetab*cthetab

!      write(*,*) "---------------------------Test diagonalization"
!      write(*,*) "thetab",thetab,sthetab,cthetab
!      write(*,*) "11",cthetab**2 * msb12 + sthetab**2 * msb22
!     & ,MsbotEig(1,1)
!      write(*,*) "12",cthetab*sthetab* (msb12 - msb22)
!     & ,MsbotEig(1,2)
!      write(*,*) "22",sthetab**2 * msb12 + cthetab**2 * msb22
!     & ,MsbotEig(2,2)
!      SbotMix(1,1) = cthetab
!      SbotMix(1,2) = sthetab
!      SbotMix(2,1) = -sthetab
!      SbotMix(2,2) = cthetab
!      TestMatrix = Matmul(Matmul(SbotMix,MsbotEig),Transpose(SbotMix))
!      write(*,*) "Matrix-1",Sqrt(TestMatrix(1,1)),TestMatrix(1,2)
!      write(*,*) "Matrix-2",TestMatrix(2,1),Sqrt(TestMatrix(2,2))
!      write(*,*) "Masses",msbot1,msbot2

      end

C-}}}
C-{{{ function deltaML

      function deltaML(mtop,mstop1,mstop2,cthetat,sthetat,mb,msb1,
     &msb2,cthetab,sthetab,mgluino,mu)
      !This routine encodes the shift in the (1,1) element of the sbottom mass matrix
      implicit none
      double precision mtop,mstop1,mstop2,cthetat,mb,msb1,
     &msb2,cthetab,mgluino,lb,lb1,lb2,lt,lt1,lt2,lg,b0fin,
     &mu,sthetab,sthetat,c2thetab,s2thetab,c2thetat,s2thetat
     &,deltaML
c      double complex b0finim

      c2thetat = cthetat**2 - sthetat**2
      s2thetat = 2*sthetat*cthetat
      c2thetab = cthetab**2 - sthetab**2
      s2thetab = 2*sthetab*cthetab

      lb = 2.d0*dlog(mu/mb)
      lb1 = 2.d0*dlog(mu/msb1)
      lb2 = 2.d0*dlog(mu/msb2)
      lt = 2.d0*dlog(mu/mtop)
      lt1 = 2.d0*dlog(mu/mstop1)
      lt2 = 2.d0*dlog(mu/mstop2)
      lg = 2.d0*dlog(mu/mgluino)

      deltaML =(-8*mb**2 - 4*lb*mb**2 + 4*msb1**2 + 3*c2thetab*msb1**2 + 
     -     2*lb1*msb1**2 + c2thetab*lb1*msb1**2 + 4*msb2**2 - 
     -     3*c2thetab*msb2**2 + 2*lb2*msb2**2 - c2thetab*lb2*msb2**2 - 
     -     4*mstop1**2 - 3*c2thetat*mstop1**2 - 2*lt1*mstop1**2 - 
     -     c2thetat*lt1*mstop1**2 - 4*mstop2**2 + 3*c2thetat*mstop2**2 - 
     -     2*lt2*mstop2**2 + c2thetat*lt2*mstop2**2 + 8*mtop**2 + 
     -     4*lt*mtop**2 + (mb**2 + mgluino**2 - msb1**2 - 
     -     2*mb*mgluino*s2thetab)*b0fin(mb**2,mgluino,msb1,mu) + 
     -     (mb**2 + mgluino**2 - msb2**2 + 2*mb*mgluino*s2thetab)*
     -     b0fin(mb**2,mgluino,msb2,mu) + 
     -     mb**2*b0fin(msb1**2,mb,mgluino,mu) + 
     -     c2thetab*mb**2*b0fin(msb1**2,mb,mgluino,mu) + 
     -     mgluino**2*b0fin(msb1**2,mb,mgluino,mu) + 
     -     c2thetab*mgluino**2*b0fin(msb1**2,mb,mgluino,mu) - 
     -     msb1**2*b0fin(msb1**2,mb,mgluino,mu) - 
     -     c2thetab*msb1**2*b0fin(msb1**2,mb,mgluino,mu) - 
     -     2*mb*mgluino*s2thetab*b0fin(msb1**2,mb,mgluino,mu) + 
     -     mb**2*b0fin(msb2**2,mb,mgluino,mu) - 
     -     c2thetab*mb**2*b0fin(msb2**2,mb,mgluino,mu) + 
     -     mgluino**2*b0fin(msb2**2,mb,mgluino,mu) - 
     -     c2thetab*mgluino**2*b0fin(msb2**2,mb,mgluino,mu) - 
     -     msb2**2*b0fin(msb2**2,mb,mgluino,mu) + 
     -     c2thetab*msb2**2*b0fin(msb2**2,mb,mgluino,mu) + 
     -     2*mb*mgluino*s2thetab*b0fin(msb2**2,mb,mgluino,mu) - 
     -     mgluino**2*b0fin(mstop1**2,mtop,mgluino,mu) - 
     -     c2thetat*mgluino**2*b0fin(mstop1**2,mtop,mgluino,mu) + 
     -     mstop1**2*b0fin(mstop1**2,mtop,mgluino,mu) + 
     -     c2thetat*mstop1**2*b0fin(mstop1**2,mtop,mgluino,mu) - 
     -     mtop**2*b0fin(mstop1**2,mtop,mgluino,mu) - 
     -     c2thetat*mtop**2*b0fin(mstop1**2,mtop,mgluino,mu) + 
     -     2*mgluino*mtop*s2thetat*b0fin(mstop1**2,mtop,mgluino,mu) - 
     -     mgluino**2*b0fin(mstop2**2,mtop,mgluino,mu) + 
     -     c2thetat*mgluino**2*b0fin(mstop2**2,mtop,mgluino,mu) + 
     -     mstop2**2*b0fin(mstop2**2,mtop,mgluino,mu) - 
     -     c2thetat*mstop2**2*b0fin(mstop2**2,mtop,mgluino,mu) - 
     -     mtop**2*b0fin(mstop2**2,mtop,mgluino,mu) + 
     -     c2thetat*mtop**2*b0fin(mstop2**2,mtop,mgluino,mu) - 
     -     2*mgluino*mtop*s2thetat*b0fin(mstop2**2,mtop,mgluino,mu) - 
     -     mgluino**2*b0fin(mtop**2,mgluino,mstop1,mu) + 
     -     mstop1**2*b0fin(mtop**2,mgluino,mstop1,mu) - 
     -     mtop**2*b0fin(mtop**2,mgluino,mstop1,mu) + 
     -     2*mgluino*mtop*s2thetat*
     -     b0fin(mtop**2,mgluino,mstop1,mu) - 
     -     mgluino**2*b0fin(mtop**2,mgluino,mstop2,mu) + 
     -     mstop2**2*b0fin(mtop**2,mgluino,mstop2,mu) - 
     -     mtop**2*b0fin(mtop**2,mgluino,mstop2,mu) - 
     -     2*mgluino*mtop*s2thetat*b0fin(mtop**2,mgluino,mstop2,mu)
     -     )/3.d0
      end

C-}}}

C-{{{ sbottomrenormHM

      subroutine sbottomrenormHM(mbsb,mb,mgl,msb12,msb22,sthetab,cthetab
     &,dhbhb_out,dthetab_out,dmsb1_out,dmsb2_out,mu_in)
      !Calculation of DR->OS counterterms in the sbottom sector
      implicit none
      !output
      double precision dhbhb_out(3), dthetab_out(2)
      double precision dmsb1_out(2), dmsb2_out(2)
      !input
      double precision mbsb,mb,mgl,msb12,msb22,sthetab,cthetab,mu_in
      !other routines
      double precision dmsb1osfin,dmsb2osfin,dthetabosfin,
     &dmbfin
      !local
      double precision s2thetab,c2thetab,msbot1,msbot2,mb2

      s2thetab = 2.d0*sthetab*cthetab
      c2thetab = cthetab**2 - sthetab**2
      msbot1 = dsqrt(msb12)
      msbot2 = dsqrt(msb22)
      mb2 = mb**2

      dmsb1_out(1) 
     & = dmsb1osfin(mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in)/msbot1
      dmsb1_out(2) = 0.d0
      dmsb2_out(1) 
     & = dmsb2osfin(mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in)/msbot2
      dmsb2_out(2) = 0.d0
      dthetab_out(1) = dthetabosfin(mb2,mgl,msb12,msb22
     &,s2thetab,c2thetab,mu_in)
      dthetab_out(2) = 0.d0
      dhbhb_out(1) = 0.d0
      dhbhb_out(2) = dmbfin(.true.,.true.,mbsb
     &,mgl,msb12,msb22,s2thetab,mu_in)
      dhbhb_out(3) = 0.d0

      end

C-}}}
C-{{{ sbottomrenormDS

      subroutine sbottomrenormDS(mb,thetab,mgl,msb12,msb22,dhbhb_out,
     &dthetab_out,dmsb1_out,dmsb2_out,mu)
      !Calculation of DR->OS counterterms in the sbottom sector
      !formulas from 1007.3465 (B1)-(B4): Degrassi, Slavich
      implicit none
      !output
      double precision dhbhb_out(3), dthetab_out(2)
      double precision thetab, dmsb1_out(2), dmsb2_out(2)
      !input and local
      double precision mb, mb2, mu
      double precision cthetab2, sthetab2,x1,x2
      double precision msb12,msb22,mgl

      cthetab2 = dcos(2.d0*thetab)
      sthetab2 = dsin(2.d0*thetab)
      mb2 = mb**2

      dmsb1_out(1) = msb12 * 1/3.d0 *
     & (3.d0 * dlog(msb12/mu**2) 
     & - 3.d0 - cthetab2**2 * (dlog(msb12/mu**2) - 1.d0)
     & - sthetab2**2 * msb22/msb12 * (dlog(msb22/mu**2) - 1.d0)
     & - 6.d0 * mgl**2/msb12
     & - 4.d0 * (1.d0 - 2.d0 * mgl**2/msb12)*dlog(mgl/mu)
     & - 2.d0 * (1.d0 - mgl**2/msb12)**2
     &   * dlog(dabs(1.d0 - msb12/mgl**2))  )

      dmsb1_out(2) = (-4*sthetab2*mgl
     -     *(-mgl**2*dlog(dabs(1.d0 - msb12/mgl**2)) + 
     -       msb12*(-2 + dlog(dabs(mgl**2 - msb12)/mu**2)))*
     -     dsqrt(mb2))/(3.d0*msb12)

      dmsb2_out(1) = msb22 * 1/3.d0 *
     & (3.d0 * dlog(msb22/mu**2) 
     & - 3.d0 - cthetab2**2 * (dlog(msb22/mu**2) - 1.d0)
     & - sthetab2**2 * msb12/msb22 * (dlog(msb12/mu**2) - 1.d0)
     & - 6.d0 * mgl**2/msb22
     & - 4.d0 * (1.d0 - 2.d0 * mgl**2/msb22)*dlog(mgl/mu)
     & - 2.d0 * (1.d0 - mgl**2/msb22)**2
     &   * dlog(dabs(1.d0 - msb22/mgl**2))  )

      dmsb2_out(2) = (4*sthetab2*mgl*
     -     (-mgl**2*dlog(dabs(1.d0 - msb22/mgl**2)) + 
     -       msb22*(-2 + dlog(dabs(mgl**2 - msb22)/mu**2)))*
     -     dsqrt(mb2))/(3.d0*msb22)

      dthetab_out(1) = sthetab2 * 1/3.d0 *
     & (-2.d0 * cthetab2**2 + 2.d0 * cthetab2**2 /(msb12 - msb22)
     &  * (msb12 * dlog(msb12/mu**2) - msb22 * dlog(msb22/mu**2)))

      dthetab_out(2) = (4*mgl*cthetab2**2
     -     *(mgl**2*msb22*dlog(dabs(1.d0 - msb12/mgl**2)) + 
     -      msb12*(mgl**2*dlog(dabs(1.d0 - msb22/mgl**2)) + 
     -      msb22*(4 - dlog(dabs((mgl**2 - msb12)*(mgl**2 - msb22))/
     -     mu**4))))*dsqrt(mb2))/
     -     (3.d0*msb12*(msb12 - msb22)*msb22)

      x1 = msb12/mgl**2
      x2 = msb22/mgl**2

!     dhbhb_pole from (32) and (29)
      dhbhb_out(1) = (sthetab2*mgl*((x1*dlog(x1))/(-1 + x1) - 
     -     (x2*dlog(x2))/(-1 + x2)))/(3.*dsqrt(mb2))
      dhbhb_out(2) = dlog(mb2/mu**2) + 
     -   (-5 + (-3 + x1)/(4.*(-1 + x1)) + (-3 + x2)/(4.*(-1 + x2)) - 
     -      dlog(mgl**2/mu**2) - 
     -      ((-2 + x1)*x1*dlog(x1))/(2.*(-1 + x1)**2) - 
     -      ((-2 + x2)*x2*dlog(x2))/(2.*(-1 + x2)**2))/3.d0
      dhbhb_out(3) = sthetab2 * (
     -     ((-3 + x1)*x1)/(2.*(-1 + x1)**2) - 
     -     ((-3 + x2)*x2)/(2.*(-1 + x2)**2) + 
     -     (x1*dlog(x1))/(-1 + x1)**3 - 
     -     (x2*dlog(x2))/(-1 + x2)**3 ) * dsqrt(mb2)/(3.d0*mgl)

      dhbhb_out = dhbhb_out*mb
      dmsb1_out = dmsb1_out/2.d0/msb12
      dmsb2_out = dmsb2_out/2.d0/msb22
      dthetab_out = dthetab_out/2.d0/dcos(2*thetab)

      end

C-}}}
C-{{{ function dAbosfinmb

      function dAbosfinmbHM(beta,mb,mgl,msb12,msb22
     &,Ab,muSUSY,deltamb,mu_in)
      implicit none
      !output
      double precision dAbosfinmbHM
      !input
      double precision beta,deltamb,mu_in,mb,mgl,msb12,msb22,Ab,muSUSY
      !other routines
      double precision b0fin
      !local
      double precision tanb,lb1,lb2,mb2

      tanb = dtan(beta)
      lb1 = dlog(mu_in**2/msb12)
      lb2 = dlog(mu_in**2/msb22)
      mb2 = mb**2

      dAbosfinmbHM = - ((4*lb1*mb*msb12 + (3*deltamb+ 8*mb)
     - *(msb12 - msb22) - 4*lb2*mb*msb22)*(Ab + muSUSY/tanb) + 
     - 2*mb*(Ab*(mb2 + mgl**2 - msb12) + mgl*(-msb12 + msb22) + 
     - ((mb2 + mgl**2 - msb12)*muSUSY)/tanb)*b0fin(msb12,mb,mgl,mu_in)
     - - 2*mb*(Ab*(mb2 + mgl**2 - msb22) + mgl*(msb12 - msb22) + 
     - ((mb2 + mgl**2 - msb22)*muSUSY)/tanb)*b0fin(msb22,mb,mgl,mu_in))/
     - (3.d0*mb*(msb12 - msb22))

      end

      function dAbosfinmbDS(beta,mb,mgl,msb12,msb22,Ab,muSUSY
     &,dmbsbos1,dmbsbos2,deltamb,mu_in)
      implicit none
      !output
      double precision dAbosfinmbDS
      !input
      double precision beta,mb,mgl,msb12,msb22,Ab,muSUSY
      double precision dmbsbos1,dmbsbos2,deltamb,mu_in
      !local
      double precision tanb

      tanb = dtan(beta)

      dAbosfinmbDS = (2*mgl*(4
     -     + ((mgl**2 - msb12)*dlog(dabs(1 - msb12/mgl**2)))/msb12
     -     + ((mgl**2 - msb22)*dlog(dabs(1 - msb22/mgl**2)))/msb22
     -     - (3*(deltamb-(dmbsbos1+dmbsbos2))
     -     /mb*(Ab+muSUSY/dtan(beta)))/(2.d0*mgl)
     -     + 4*dlog(mu_in/mgl) ))/3.d0

      end

C-}}}
C-{{{ function dAbosfintb

      function dAbosfintbHM(beta,mb,mgl,msb12,msb22,
     &sthetab,cthetab,Ab,muSUSY,dmsb1_1,dmsb2_1,deltatb,mu_in)
      implicit none
      !output
      double precision dAbosfintbHM
      !input
      double precision tanb,deltatb,mu_in,dmsb1_1,dmsb2_1,mb
      double precision mgl,msb12,msb22,Ab,muSUSY,beta,sthetab,cthetab
      !other routines
      double precision b0fin
      !local
      double precision dmsb12,dmsb22,s2thetab,c2thetab,lb1,lb2,mb2

      tanb = dtan(beta)
      dmsb12 = 2.d0*dmsb1_1*msb12
      dmsb22 = 2.d0*dmsb2_1*msb22
      lb1 = dlog(mu_in**2/msb12)
      lb2 = dlog(mu_in**2/msb22)
      mb2 = mb**2

      s2thetab = 2*sthetab*cthetab
      c2thetab = cthetab**2 - sthetab**2

      dAbosfintbHM = ((6*c2thetab*deltatb*(msb12 - msb22) + 
     -   (3*dmsb12 - 3*dmsb22 + 8*msb12 + 4*lb1*msb12 - 8*msb22 - 
     -   4*lb2*msb22)*s2thetab)*(Ab + muSUSY/tanb) - 
     -   (2*s2thetab*b0fin(msb22,mb,mgl,mu_in)*
     -   ((mb2 + mgl**2 - msb22)*muSUSY*dcos(beta) + 
     -   (Ab*(mb2 + mgl**2 - msb22) + mgl*(msb12 - msb22))*dsin(beta)))/
     -   dsin(beta) + (2*s2thetab*b0fin(msb12,mb,mgl,mu_in)*
     -   ((mb2 + mgl**2 - msb12)*muSUSY*dcos(beta) + 
     -   (Ab*(mb2 + mgl**2 - msb12) + mgl*(-msb12 + msb22))*dsin(beta)))
     -   /dsin(beta))/
     -   (6*Ab*mb - 3*(msb12 - msb22)*s2thetab + (6*mb*muSUSY)/tanb) 

      end

      function dAbosfintbDS(mgl,msb12,msb22,
     &thetab,Ab,muSUSY,beta,dthetabos1,deltatb,mu_in)
      implicit none
      !output
      double precision dAbosfintbDS
      !input
      double precision mgl,msb12,msb22,Ab,muSUSY,beta,thetab
      double precision dthetabos1,deltatb,mu_in

      dAbosfintbDS = (2*mgl*(4
     -  + ((mgl**2 - msb12)*dlog(dabs(1 - msb12/mgl**2)))/msb12
     -  + ((mgl**2 - msb22)*dlog(dabs(1 - msb22/mgl**2)))/msb22
     -  - (3*(deltatb - dthetabos1)/dtan(2*thetab)*
     -    (Ab+muSUSY/dtan(beta)))/mgl
     -  + 4*dlog(mu_in/mgl) ))/3.d0

      end

C-}}}
C-{{{ function mbMSDRtrans

      function mbMSDRtrans(mb,orderin,alphas)
      implicit none
      integer orderin
      double precision mbMSDRtrans, mb, alphas

      include '../commons/common-consts.f'

      mbMSDRtrans = mb * (1.d0 - 1.d0/3.d0*alphas/Pi)

      if (orderin.eq.1) then
      mbMSDRtrans = mbMSDRtrans - mb*(29.d0/72.d0 *(alphas/Pi)**2)
      end if

      end

C-}}}
C-{{{ function dmbfin

      function dmbfin(smflag,drflag,mb,mgl,msb12,msb22,s2thetab,mu)
      implicit none
      !output
      double precision dmbfin
      !input
      double precision mb,mgl,msb12,msb22,s2thetab,mu
      logical smflag,drflag
      !local
      double precision lb, lb1, lb2, lg, mb2
      !other routines
      double precision b0fin

      mb2 = mb**2
      lb = dlog(mu**2/mb2)
      lb1 = dlog(mu**2/msb12)
      lb2 = dlog(mu**2/msb22)
      lg = 2.d0*dlog(mu/mgl)

      dmbfin = (- 2.d0*(1.d0 +lg)*mgl**2+
     -     msb12 + lb1*msb12 + msb22 + lb2*msb22 +
     -     (mb2 + mgl**2 - msb12 - 2.d0*mb*mgl*s2thetab)*
     -     b0fin(mb2,mgl,dsqrt(msb12),mu) 
     -    + (mb2 + mgl**2 - msb22 +
     -     2*mb*mgl*s2thetab)*b0fin(mb2,mgl,dsqrt(msb22),mu))
     -    /(6.d0*mb)

      if (smflag.eqv.(.true.)) then
      if (drflag.eqv.(.true.)) then
         dmbfin = dmbfin - 2.d0*(5.d0+ 3.d0*lb)*mb2/(6.d0*mb)
      else
         dmbfin = dmbfin - 2.d0*(4.d0+ 3.d0*lb)*mb2/(6.d0*mb)
      end if
      end if
      end

C-}}}
C-{{{ function dmbnontanb - non-tanbeta enhanced terms only

      function dmbnontanb(smflag,drflag,mb,mgl,msb12,msb22,Ab,mu)
      implicit none
      !output
      double precision dmbnontanb
      !input
      double precision mb,mgl,msb12,msb22,Ab,mu
      logical smflag,drflag
      !local
      double precision lb, lb1, lb2, lg, mb2
      !other routines
      double precision b0fin

      mb2 = mb**2
      lb = dlog(mu**2/mb2)
      lb1 = dlog(mu**2/msb12)
      lb2 = dlog(mu**2/msb22)
      lg = 2.d0*dlog(mu/mgl)

      dmbnontanb = (- 2.d0*(1.d0+lg)*mgl**2+
     -     msb12 + lb1*msb12 + msb22 + lb2*msb22 +
     -     (mb2 + mgl**2 - msb12 - 4.d0*mb2*mgl*Ab/(msb12-msb22))*
     -     b0fin(mb2,mgl,dsqrt(msb12),mu) 
     -      + (mb2 + mgl**2 - msb22 + 4*mb2*mgl*Ab/(msb12-msb22))
     -     *b0fin(mb2,mgl,dsqrt(msb22),mu))
     -    /(6.d0*mb)

      if (smflag.eqv.(.true.)) then
      if (drflag.eqv.(.true.)) then
         dmbnontanb = dmbnontanb - 2.d0*(5.d0+ 3.d0*lb)*mb2/(6.d0*mb)
      else
         dmbnontanb = dmbnontanb - 2.d0*(4.d0+ 3.d0*lb)*mb2/(6.d0*mb)
      end if
      end if

      end

C-}}}
C-{{{ function dmbtanb - tanbeta enhanced terms only

      function dmbtanb(mb,mgl,msb12,msb22,beta,muSUSY,mu)
      implicit none
      !output
      double precision dmbtanb
      !input
      double precision mb,mgl,msb12,msb22,beta,muSUSY,mu
      !local
      double precision mb2
      !other routines
      double precision b0fin

      mb2 = mb**2

      dmbtanb = (4*mb2*mgl*muSUSY*tan(beta)/(msb12-msb22)
     -     *b0fin(mb2,mgl,dsqrt(msb12),mu)  -
     -     4*mb2*mgl*muSUSY*tan(beta)/(msb12-msb22)
     -     *b0fin(mb2,mgl,dsqrt(msb22),mu))
     -    /(6.d0*mb)

      end

C-}}}
C-{{{ function dthetabosfin

      function dthetabosfin(mb2,mgl,msb12,msb22
     &,s2thetab,c2thetab,mu_in)
      implicit none
      !output
      double precision dthetabosfin
      !input
      double precision mb2,mgl,msb12,msb22,s2thetab,c2thetab,mu_in
      !local
      double precision lb1,lb2
      !other routines
      double precision b0fin

      lb1 = dlog(mu_in**2/msb12)
      lb2 = dlog(mu_in**2/msb22)

      dthetabosfin = (c2thetab*((-((1 + lb1)*msb12) + 
     - (1 + lb2)*msb22)*s2thetab + 
     - 2.d0*dsqrt(mb2)*mgl*(b0fin(msb12,dsqrt(mb2),mgl,mu_in)
     -    + b0fin(msb22,dsqrt(mb2),mgl,mu_in))))/
     - (3.d0*(msb12 - msb22))

      end

C-}}}
C-{{{ function dmsb1osfin

      function dmsb1osfin(mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in)
      implicit none
      !output
      double precision dmsb1osfin
      !input
      double precision mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in
      !local
      double precision lb,lb1,lb2,lg,s2thetab
      !other routines
      double precision b0fin

      lb = dlog(mu_in**2/mb2)
      lb1 = dlog(mu_in**2/msb12)
      lb2 = dlog(mu_in**2/msb22)
      lg = 2.d0*dlog(mu_in/mgl)

      s2thetab = 2*sthetab*cthetab

      dmsb1osfin = -((1 + lb)*mb2 + (1 + lg)*mgl**2 + (3 + lb1)*msb12 + 
     -     2*cthetab**2*((1 + lb1)*msb12 - (1 + lb2)*msb22)*sthetab**2 + 
     -(mb2 + mgl**2 - msb12 
     -   - 2*dsqrt(mb2)*mgl*s2thetab)*b0fin(msb12,dsqrt(mb2),mgl,mu_in)
     -      )/(3.d0*dsqrt(msb12))

      end

C-}}}
C-{{{ function dmsb2osfin

      function dmsb2osfin(mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in)
      implicit none
      !output
      double precision dmsb2osfin
      !input
      double precision mb2,mgl,msb12,msb22,sthetab,cthetab,mu_in
      !local
      double precision lb,lb1,lb2,lg,s2thetab
      !other routines
      double precision b0fin

      lb = dlog(mu_in**2/mb2)
      lb1 = dlog(mu_in**2/msb12)
      lb2 = dlog(mu_in**2/msb22)
      lg = 2.d0*dlog(mu_in/mgl)

      s2thetab = 2*sthetab*cthetab

      dmsb2osfin = -((1 + lb)*mb2 + (1 + lg)*mgl**2 + (3 + lb2)*msb22 + 
     -  2*cthetab**2*(-((1 + lb1)*msb12) + (1 + lb2)*msb22)*sthetab**2 + 
     -(mb2 + mgl**2 - msb22 
     -   + 2*dsqrt(mb2)*mgl*s2thetab)*b0fin(msb22,dsqrt(mb2),mgl,mu_in)
     -     )/(3.d0*dsqrt(msb22))

      end

C-}}}
C-{{{ dAtosfins2

      function dAtosfins2(mtord,At,muSUSY,beta,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,muRggh)
      !finite part of delta A_t^OS * sin(2 thetat)
      implicit none
      !output
      double precision dAtosfins2
      !input
      double precision mt2,mst12,mst22,mgl,s2thetat,c2thetat,muRggh
      double precision At,muSUSY,beta
      integer mtord
      !other routines
      double precision dmst12osfin,dmst22osfin
      double precision dthetatosfin,dmtosfin

      dAtosfins2 = (s2thetat*(dmst12osfin(muRggh,mt2,mst12,mst22
     &        ,mgl,s2thetat,c2thetat,mtord)
     &   - dmst22osfin(muRggh,mt2,mst12,mst22
     &        ,mgl,s2thetat,c2thetat,mtord))/(mst12-mst22)
     & + 2*c2thetat*dthetatosfin(muRggh,mt2,mst12,mst22
     &        ,mgl,s2thetat,c2thetat,mtord)
     &   - dmtosfin(.true.,mtord+1,mt2,mgl
     &       ,s2thetat,mst12,mst22,muRggh)*s2thetat/dsqrt(mt2)
     & ) * (At - muSUSY/dtan(beta))

      end

C-}}}
C-{{{ dmtosfin

      function dmtosfin(smflag,mtord,mt2,mgl
     &,s2thetat,mst12,mst22,mu_in)
      implicit none
      !output
      double precision dmtosfin
      !input
      double precision mu_in,mt2,mgl,mst12,mst22,s2thetat
      integer mtord
      logical smflag
      !other routines
      double precision dmtossusy
      !local
      double precision lt

      dmtosfin = 0.d0

      lt = dlog(mu_in**2/mt2)
      
      dmtosfin = dmtossusy(mtord,mt2,mgl,s2thetat,mst12/mgl**2,mu_in)
     &  + dmtossusy(mtord,mt2,mgl,-s2thetat,mst22/mgl**2,mu_in)

      if ((smflag.eqv.(.true.)).and.(mtord.gt.0)) then
         dmtosfin = dmtosfin - (5.d0 + 3.d0*lt)*dsqrt(mt2)/3.d0
      endif

      end

C-}}}
C-{{{ dmtossusy

      function dmtossusy(mtord,mt2,mgl,s2thetat,x_in,mu_in)
      implicit none
      !output
      double precision dmtossusy
      !input
      double precision mu_in,x_in,s2thetat,mt2,mgl
      integer mtord

      dmtossusy = 0.d0  

      if (mtord.gt.(-1)) then
      dmtossusy = -1/3.d0 *
     & s2thetat * mgl * x_in /(1.d0 - x_in) *dlog(x_in)
      end if

      if (mtord.gt.0) then
      dmtossusy = dmtossusy - 1/3.d0 * dsqrt(mt2)
     & * (dlog(mgl/mu_in) + (x_in-3.d0)/(4.d0*(1.d0-x_in))
     & + x_in * (x_in - 2.d0) / (2.d0 * (1.d0 - x_in)**2) * dlog(x_in))
      end if

      if (mtord.gt.1) then
      dmtossusy = dmtossusy - 1/3.d0 * mt2 * s2thetat
     & / (2.d0 * mgl * (1.d0-x_in)**3)
     & * (1.d0 - x_in**2 + 2.d0 * x_in * dlog(x_in))
      end if

      if (mtord.gt.2) then
      dmtossusy = dmtossusy - 1/3.d0 * mt2 * dsqrt(mt2)
     & / (6.d0 * mgl**2 * (1.d0-x_in)**3)
     & * (x_in**2 - 5.d0 * x_in - 2.d0
     &    - 6.d0 * x_in/(1.d0-x_in) * dlog(x_in))
      end if

      end

C-}}}
C-{{{ dthetatosfin

      function dthetatosfin(mu_in,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,mtord)
      implicit none
      !output
      double precision dthetatosfin
      !input
      double precision mu_in,mst12,mst22,mgl,s2thetat,c2thetat,mt2
      integer mtord

      dthetatosfin = 0.d0

      if (mtord.gt.(-1)) then
      dthetatosfin = s2thetat * c2thetat * 1/3.d0 / (mst12 - mst22)
     & * (mst12 * (dlog(mst12/mu_in**2) - 1.d0))
     &   + s2thetat * c2thetat * 1/3.d0 / (mst22 - mst12)
     & * (mst22 * (dlog(mst22/mu_in**2) - 1.d0))
      end if

      if (mtord.gt.0) then
      dthetatosfin = dthetatosfin 
     &  - c2thetat * dsqrt(mt2) * mgl * 2/3.d0 / (mst12 - mst22)
     &     *(2.d0 * dlog(mgl/mu_in) + (1.d0-mgl**2/mst12)
     &        *dlog(dabs(1.d0 - mst12/mgl**2)) - 2.d0)
     &  + c2thetat * dsqrt(mt2) * mgl * 2/3.d0 / (mst22 - mst12)
     &     *(2.d0 * dlog(mgl/mu_in) + (1.d0-mgl**2/mst22)
     &        *dlog(dabs(1.d0 - mst22/mgl**2)) - 2.d0)
       end if

      end

C-}}}
C-{{{ dmst12osfin

      function dmst12osfin(mu_in,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,mtord)
      implicit none
      !output
      double precision dmst12osfin
      !input
      double precision mu_in,mst12,mst22,mgl,s2thetat,c2thetat,mt2
      integer mtord

      dmst12osfin = 0.d0

      if (mtord.gt.(-1)) then
      dmst12osfin = mst12 * 1/3.d0 *
     & (3.d0 * dlog(mst12/mu_in**2) 
     & - 3.d0 - c2thetat**2 * (dlog(mst12/mu_in**2) - 1.d0)
     & - s2thetat**2 * mst22/mst12 * (dlog(mst22/mu_in**2) - 1.d0)
     & - 6.d0 * mgl**2/mst12
     & - 4.d0 * (1.d0 - 2.d0 * mgl**2/mst12)*dlog(mgl/mu_in)
     & - 2.d0 * (1.d0 - mgl**2/mst12)**2
     &   * dlog(dabs(1.d0 - mst12/mgl**2))  )
      end if

      if (mtord.gt.0) then
      dmst12osfin = dmst12osfin - 4/3.d0*s2thetat*mgl*dsqrt(mt2)
     &     *(2.d0 * dlog(mgl/mu_in) + (1.d0-mgl**2/mst12)
     &        *dlog(dabs(1.d0 - mst12/mgl**2)) - 2.d0)
      end if

      end

C-}}}
C-{{{ dmst22osfin

      function dmst22osfin(mu_in,mt2,mst12,mst22
     &,mgl,s2thetat,c2thetat,mtord)
      implicit none
      !output
      double precision dmst22osfin
      !input
      double precision mu_in,mst12,mst22,mgl,s2thetat,c2thetat,mt2
      integer mtord

      dmst22osfin = 0.d0

      if (mtord.gt.(-1)) then
      dmst22osfin = mst22 * 1/3.d0 *
     & (3.d0 * dlog(mst22/mu_in**2) 
     & - 3.d0 - c2thetat**2 * (dlog(mst22/mu_in**2) - 1.d0)
     & - s2thetat**2 * mst12/mst22 * (dlog(mst12/mu_in**2) - 1.d0)
     & - 6.d0 * mgl**2/mst22
     & - 4.d0 * (1.d0 - 2.d0 * mgl**2/mst22)*dlog(mgl/mu_in)
     & - 2.d0 * (1.d0 - mgl**2/mst22)**2
     &   * dlog(dabs(1.d0 - mst22/mgl**2))  )
      end if

      if (mtord.gt.0) then
      dmst22osfin = dmst22osfin + 4/3.d0*s2thetat*mgl*dsqrt(mt2)
     &     *(2.d0 * dlog(mgl/mu_in) + (1.d0-mgl**2/mst22)
     &        *dlog(dabs(1.d0 - mst22/mgl**2)) - 2.d0)
      end if

      end

C-}}}

C-{{{ alphas

      real*8 function SUSHI_alphas(rmurl)
      implicit none
      double precision rmurl,api0,apiout

      include '../commons/common-consts.f'      
      include '../commons/common-vars.f' 

      api0 = AlfasMZ/Pi

      call runalpha(api0,mz,rmurl,nf,4,0,apiout)

      SUSHI_alphas = apiout*Pi

      end

C-}}}
C-{{{ subroutine polemass:

      subroutine polemasshm(mqmq,apimq,nloop,mqpole)
c..
c..   Computes the pole mass mqpole from the MS-bar mass mqmq = mq(mq).
c..   apimq = alpha_s(mqmq)/pi.
c..
c      implicit real*8 (a-z)
      implicit none
      integer nloop
      double precision mqmq,apimq,mqpole,lmum

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'

      lmum = 0.d0

      if (nloop.eq.0) then
         mqpole = mqmq
      elseif (nloop.eq.1) then
         mqpole = mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.d0) )
      elseif (nloop.eq.2) then
         mqpole = mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.d0)
     &        + apimq**2*(
     &        + cf**2*( -71/128.d0 - (9*lmum)/32.d0 + (9*lmum**2)/32.d0
     &        +(15*z2)/8.d0  - 3*ln2*z2 + (3*z3)/4.d0 )+ cf*ca*( 1111
     &        /384.d0 + (185*lmum)/96.d0 + (11*lmum**2)/32.d0  - z2/2.d0
     &        + (3*ln2*z2)/2.d0 - (3*z3)/8.d0 )+ cf*tr*( -3/4.d0 - (71
     &        *nf)/96.d0 - (13*lmum*nf)/24.d0 -(lmum**2*nf)/8.d0 + (3*z2
     &        )/2.d0 - (nf*z2)/2.d0) ) )
      elseif (nloop.eq.3) then
         mqpole = mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.d0)
     &        + apimq**2*(
     &        + cf**2*( -71/128.d0 - (9*lmum)/32.d0 + (9*lmum**2)/32.d0
     &        +(15*z2)/8.d0  - 3*ln2*z2 + (3*z3)/4.d0 )+ cf*ca*( 1111
     &        /384.d0 + (185*lmum)/96.d0 + (11*lmum**2)/32.d0  - z2/2.d0
     &        + (3*ln2*z2)/2.d0 - (3*z3)/8.d0 )+ cf*tr*( -3/4.d0 - (71
     &        *nf)/96.d0 - (13*lmum*nf)/24.d0 -(lmum**2*nf)/8.d0 + (3*z2
     &        )/2.d0 - (nf*z2)/2.d0) )
     &        + apimq**3*(.6527d0*(nf-1.d0)**2 - 26.655d0*(nf-1.d0)
     &        +190.595d0) ) 
      else
         write(6,*) '<function polemass>: nloop = ',nloop
     &        ,' not implemented'
      endif

      end

C-}}}
C-{{{ subroutine getmass

      subroutine getmass()
c      from ggh@nnlo
      implicit none
      double precision apim, apimt, apimuB!, apimuD
      double precision SUSHI_alphas

      include '../commons/common-consts.f'
      include '../commons/common-quark.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'
      include '../commons/common-ren.f'

c      external runmass
c      apimuR = alphaggh/Pi
      apimuR = SUSHI_alphas(muRggh)/Pi

      apim  = SUSHI_alphas(mcmc)/Pi
      call polemasshm(mcmc,apim,3,mcos) !3-loop conversion

!mass mbos
      apim  = SUSHI_alphas(mbmb)/Pi
!      call polemasshm(mbmb,apim,1,mbos) !1-loop conversion
!      call polemasshm(mbmb,apim,2,mbos) !2-loop conversion
      call polemasshm(mbmb,apim,3,mbos) !3-loop conversion


!mass mbmuR
!      call runmass(mbmb,apim,apimuR,nf,norderggh+1,mbMSbarmuR)
      call runmass(mbmb,apim,apimuR,nf,4,mbMSbarmuR) !4-loop running


!mass mbmuB
      apimuB  = SUSHI_alphas(muB)/Pi
      call runmass(mbmb,apim,apimuB,nf,4,mbMSbarmuB) !4-loop running

!mass mbmt
      apimt  = SUSHI_alphas(mt)/Pi
      call runmass(mbmb,apim,apimt,nf,4,mbMSbarmt) !4-loop running

      mc2 = mcos**2

      end

C-}}}
