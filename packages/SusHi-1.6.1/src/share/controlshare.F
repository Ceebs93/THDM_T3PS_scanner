C-{{{ function amplosq:

      real*8 function amplosq(nexp)

      implicit none
      integer nexp
      real*8 atauexph,amplore,amploim
      complex*16 atau,tildeatau
      complex*16 amplo,amplot,amplost1,amplost2
      complex*16 ALOSUSYh
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'

      include '../commons/common-ren.f'
      include '../commons/common-quark.f'

      errnorm = 0

      if (nexp.eq.-2) then
         amplot = gt * atau( 4*mt**2/mh**2 )
      else
         amplot = gt * dsqrt(atauexph( 4*mt**2/mh**2, nexp ))*(1.d0
     &        ,0.d0)
      endif
      amplost1 = (0.d0,0.d0)
      amplost2 = (0.d0,0.d0)

      if ((.not.lpseudo).and.nnlostop) then
         amplost1 = mt**2*gt1*ALOSUSYh(mh**2,mst12)
         amplost2 = mt**2*gt2*ALOSUSYh(mh**2,mst22)
      end if

      amplo = amplot + amplost1 + amplost2

      if (lpseudo) then
         amplo = amplo + c5sp(1,0)
      else
         amplo = amplo + c5sp(0,0)
      end if
      
      amplore = dreal( amplo )
      amploim = dimag( amplo )
      amplosq = amplore**2 + amploim**2
      
      end

C-}}}

C-{{{ function atauexph(tau,n):

      real*8 function atauexph(tau,nn)
c..
c..   Expansion of LO cross section (square of atau(tau)) in mh/mt:
c..
      implicit real*8 (a-h,o-z)
      integer nn
      include '../commons/common-keys.f'

      taum1 = 1/tau

      atauexph = 1.d0
      if (.not.lpseudo) then
         if (nn.ge.2) then
            atauexph = atauexph + 7/15.d0*taum1
         endif
         if (nn.ge.4) then
            atauexph = atauexph + 1543/6300.d0*taum1**2
         endif
         if (nn.ge.6) then
            atauexph = atauexph + 226/1575.d0*taum1**3
         endif
         if (nn.ge.8) then
            atauexph = atauexph + 55354/606375.d0*taum1**4
         endif
         if (nn.ge.10) then
            atauexph = atauexph + 1461224/23648625.d0*taum1**5 
         endif
      else
         if (nn.ge.2) then
            atauexph = atauexph + taum1/3.d0
         endif
         if (nn.ge.4) then
            atauexph = atauexph + taum1**2*8/45.d0
         endif
         if (nn.ge.6) then
            atauexph = atauexph + taum1**3*4/35.d0
         endif
         if (nn.ge.8) then
            atauexph = atauexph + taum1**4*128/1575.d0
         endif
         if (nn.ge.10) then
            atauexph = atauexph + taum1**5*128/2079.d0
         endif
      endif

      if (nn.ge.12) stop 'stopped in atauexph()'
      end

C-}}}
C-{{{ function atau(tau):

      complex*16 function atau(tau)
c..
c..   Coefficient of LO cross section including full mass dependence.
c..
      implicit real*8 (a-h,o-z)
      complex*16 ctau,ftau
      include '../commons/common-keys.f'

c..   ctau = complex(tau,0.d0)  ! for g77 compiler
      ctau = (1.d0,0.d0)*tau

      if (lpseudo) then
c..   pseudo-scalar:
         atau = tau*ftau(ctau)
      else
c..   scalar:
         atau = 3/2.d0 * tau*( (1.d0,0) + (1.d0-tau)*ftau(ctau) )
      endif

      end

C-}}}
C-{{{ function tildeatau(tau):

      complex*16 function tildeatau(tau)
c..
c..   Coefficient of LO cross section including full mass dependence.
c..
      implicit real*8 (a-h,o-z)
      complex*16 ctau,ftau
      include '../commons/common-keys.f'

c..   ctau = complex(tau,0.d0)  ! for g77 compiler
      ctau = (1.d0,0.d0)*tau

      tildeatau = -3/4.d0*tau*( (1.d0,0) - tau*ftau(ctau) )

      end

C-}}}
C-{{{ function ftau(tau):

      complex*16 function ftau(ctau)
c..
c..   used by atau(tau).
c..
      implicit real*8 (a-h,o-z)
      complex*16 ctau
      include '../commons/common-consts.f'

      tau = dreal(ctau)

      if (tau.ge.1.d0) then
c         ftau = complex(asin(dsqrt(1/tau))**2,0)  ! for g77
         ftau = (1.d0,0)*asin(dsqrt(1/tau))**2 
      else
         snum = 1 + dsqrt(1-tau)
         sden = 1 - dsqrt(1-tau)
         dl = dlog(snum/sden)
         ftau  = -1/4.d0 * ( (1.d0,0.d0)*dl - (0,1)*pi )**2
      endif

      end

C-}}}
C-{{{ function vartrans:

      real*8 function vartrans(rlumi,partsig,xx)
c..
c..   Doing a transformation of variables.
c..   Taken from the zwprod.f code by W. van Neerven.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external rlumi,partsig

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      vartrans = tauh * pmeas * rlumi(ww/zt,tauh/ww)*partsig(zt)/zt**2

      end

C-}}}
C-{{{ function vartrans1:

      real*8 function vartrans1(rlumi,partsig,xx)
c..
c..   Similar to vartrans, but does an additional
c..   tranformation  z -> 1/z^2, so that 1/z terms are
c..   taken into consideration more properly.
c..   Used only for checks.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external rlumi,partsig

      vt = xx(1)
      xt = xx(2)

      zt=tauh/vt

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      vartrans1 = pmeas * rlumi(ww/zt,tauh/ww)*partsig(zt)

      end

C-}}}
C-{{{ subroutine convolute:

      subroutine convolute(ikey,fun,rint,del,chi2)
c..
c..   Integrate fun.
c..
      implicit real*8 (a-h,o-z)
      integer ikey
      logical lveg1
      include '../commons/common-vegpar.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'
      common/bveg1/xl(10),xu(10),acc1,ndim,ncall,itmx,nprn
      external fun

      ndim=2
      acc1=acc

      do iv=3,10
         xl(iv)=0.d0
         xu(iv)=0.d0
      enddo

      xl(1) = tauh
      xu(1) = 1.d0
      xl(2) = tauh
      xu(2) = 1.d0

      if (ikey.eq.1) then
         nprn=nprnvggh
         lveg1=lveg1ggh
         itmx=itmx1ggh
         ncall=ncall1ggh
         nprn=nprnvggh
      elseif (ikey.eq.2) then
         nprn=nprnvbbh
         lveg1=lveg1bbh
         itmx=itmx1bbh
         ncall=ncall1bbh
         nprn=nprnvbbh
      else
         call printdie('unknown ikey in convolute()')
      endif

      call vegas(fun,rintl,dell,chi2l)
      if (lveg1) then
         if (ikey.eq.1) then
            itmx=itmx2ggh
            ncall=ncall2ggh
         elseif (ikey.eq.2) then
            itmx=itmx2bbh
            ncall=ncall2bbh
         else
            call printdie('unknown ikey in convolute()')
         endif
         call vegas1(fun,rintl,dell,chi2l)
      endif

      chi2 = chi2l
      rint = rintl
      del = dell

      return
      end

C-}}}
C-{{{ function stepfun(xx):

      real*8 function stepfun(xin)

      implicit real*8 (a-h,o-z)

      if (xin.gt.0.d0) then
         stepfun = 1.d0
      else
         stepfun = 0.d0
      endif
      
      end

C-}}}

