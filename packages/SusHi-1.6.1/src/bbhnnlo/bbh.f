
C-{{{ subroutine intdel:

      subroutine intdelbbh(del,errdel)
c..
c..   Integrating the delta(1-x) part over PDFs.
c..   
c..   del:    result
c..   errdel: uncertainty
c..   
      implicit none
      real*8 acc1,xl,xu,del,errdel,chi2a,bbqfun
      integer ndim,ncall,itmx,nprn,iv
      include '../commons/common-vars.f'
      include '../commons/common-vegpar.f'
      common/bveg1/xl(10),xu(10),acc1,ndim,ncall,itmx,nprn
      external bbqfun

      ndim=1
      nprn=nprnvbbh
      acc1=acc

      do iv=2,10
         xl(iv)=0.d0
         xu(iv)=0.d0
      enddo

      itmx=itmx1bbh
      ncall=ncall1bbh
      xl(1) = tauh
      xu(1) = 1.d0

      call vegas(bbqfun,del,errdel,chi2a)

      if (lveg1bbh) then
         itmx=itmx2bbh
         ncall=ncall2bbh
         call vegas1(bbqfun,del,errdel,chi2a)
      endif

      end

C-}}}
C-{{{ function bbqfun(yy):

      real*8 function bbqfun(xt)
c..
c..   integrand for intdel
c..
      implicit none
      real*8 bbqpdf,xt
      include '../commons/common-vars.f'

      bbqfun = tauh * bbqpdf(xt,tauh/xt)/xt

      return
      end

C-}}}
C-{{{ subroutine soft1bbh(ddsoft1):

      subroutine soft1bbh(ddsoft1,errsoft1)
c..
c..   Integrating the D-terms at NLO.
c..
      implicit none
      real*8 ppdt1bbh,ddsoft1,errsoft1,chi2a
      external ppdt1bbh

      call convolute(2,ppdt1bbh,ddsoft1,errsoft1,chi2a)

      end

C-}}}
C-{{{ subroutine soft2bbh(ddsoft1):

      subroutine soft2bbh(ddsoft2,errsoft2)
c..
c..   Integrating the D-terms at NNLO.
c..
      implicit none
      real*8 ppdt2bbh,ddsoft2,errsoft2,chi2a
      external ppdt2bbh

      call convolute(2,ppdt2bbh,ddsoft2,errsoft2,chi2a)

      end

C-}}}
C-{{{ function ppdt1bbh(yy)

      real*8 function ppdt1bbh(xx,wgt)
c..
c..   Integrand for the D-terms at NLO.
c..   
      implicit none
      real*8 xx(10),zt,xt,ww,wgt,pmeas,bbqpdf,dterms1bbh
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      ppdt1bbh = tauh * ( pmeas/zt**2 * bbqpdf(ww/zt,tauh/ww)
     &     - bbqpdf(xt,tauh/xt)/xt )*dterms1bbh(zt)

      end

C-}}}
C-{{{ function ppdt2bbh(yy)

      real*8 function ppdt2bbh(xx,wgt)
C
C     Integrand for D-terms at NNLO.
C
      implicit none
      real*8 xx(10),zt,xt,ww,wgt,pmeas,bbqpdf,dterms2bbh
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      external bbqpdf,dterms2bbh

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      ppdt2bbh = tauh * ( pmeas/zt**2 * bbqpdf(ww/zt,tauh/ww)
     &     - bbqpdf(xt,tauh/xt)/xt )*dterms2bbh(zt)

      end

C-}}}
C-{{{ function delta1bbh():

      real*8 function delta1bbh()
c..
c..   Coefficient of the delta-function at NLO.
c..   
      implicit none
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      delta1bbh = -1.3333333333333333d0 - 2*lfr + (8*z2)/3.d0

      end

C-}}}
C-{{{ function delta2bbh():

      real*8 function delta2bbh()
c..
c..   Coefficient of the delta-function at NNLO.
c..   
      implicit none
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      delta2bbh = 6.388888888888889d0 + 3*lfh - (25*lfr)/12.d0 +
     &     (19*lfr**2)/4.d0+(2*nf)/27.d0 + (lfr*nf)/18.d0 -
     &     (lfr**2*nf)/6.d0+ (58*z2)/9.d0 +(8*lfh*z2)/3.d0 - (32
     &     *lfh**2*z2)/9.d0 - (38*lfr*z2)/3.d0 -(10*nf*z2)/27.d0 +
     &     (4*lfr*nf*z2)/9.d0 - (26*z3)/3.d0 - (122*lfh*z3)/9.d0
     &     +(2*nf*z3)/3.d0 - (19*z4)/18.d0

      end

C-}}}
C-{{{ function dterms1bbh(...):

      real*8 function dterms1bbh(xt)
c..
c..   The plus-distributions at NLO.
c..   
      implicit none
      real*8 dd0,dd1,xt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'

      dd0 = 1.d0/(1.d0 - xt)
      dd1 = dd0*dlog(1.d0 - xt)

      dterms1bbh =  (16*dd1)/3.d0 - (8*dd0*lfh)/3.d0

      end

C-}}}
C-{{{ function dterms2bbh(...):

      real*8 function dterms2bbh(xt)

      implicit none
      real*8 dd0,dd1,dd2,dd3,xt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      dd0 = 1.d0/(1.d0 - xt)
      dd1 = dd0*dlog(1.d0 - xt)
      dd2 = dd1*dlog(1.d0 - xt)
      dd3 = dd2*dlog(1.d0 - xt)

      dterms2bbh = (128*dd3)/9.d0 + dd2*(-14.666666666666666d0 
     &     - (64*lfh)/3.d0 +(8*nf)/9.d0) +dd1*(22.666666666666668d0 
     &     + (44*lfh)/3.d0
     &     + (64*lfh**2)/9.d0 -(76*lfr)/3.d0 + (
     &     -1.4814814814814814d0 -(8*lfh)/9.d0 + (8*lfr)/9.d0)*nf
     &     - (200*z2)/9.d0) +dd0*(-14.962962962962964d0 - (11*lfh**2)
     &     /3.d0 +nf*(0.691358024691358d0 + (2*lfh**2)/9.d0 +lfh
     &     *(0.7407407407407407d0 - (4*lfr)/9.d0) - (8*z2)/9.d0) +(44
     &     *z2)/3.d0 +lfh*(-11.333333333333334d0 + (38*lfr)/3.d0
     &     +(100*z2)/9.d0) + (382*z3)/9.d0)

      end

C-}}}
C-{{{ function dtsub1bbh(...):

      real*8 function dtsub1bbh()
c..
c..   Contributions arising from the fact that the integrals over
c..   plus-distributions do not run from 0 to 1, but from z to 1.
c..   
      implicit none
      real*8 ddz0,ddz1
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      ddz0 = dlog(1.d0 - tauh)
      ddz1 = ddz0**2/2.d0

      dtsub1bbh =  (16*ddz1)/3.d0 - (8*ddz0*lfh)/3.d0

      end

C-}}}
C-{{{ function dtsub2bbh(...):

      real*8 function dtsub2bbh()
c..
c..   Contributions arising from the fact that the integrals over
c..   plus-distributions do not run from 0 to 1, but from z to 1.
c..   
      implicit none
      real*8 ddz0,ddz1,ddz2,ddz3
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      ddz0 = dlog(1.d0 - tauh)
      ddz1 = ddz0**2/2.d0
      ddz2 = ddz0**3/3.d0
      ddz3 = ddz0**4/4.d0

      dtsub2bbh = (128*ddz3)/9.d0 + ddz2*(-14.666666666666666d0 
     &     - (64*lfh)/3.d0 +(8*nf)/9.d0) +ddz1*(22.666666666666668d0 
     &     + (44*lfh)/3.d0 + (64*lfh**2)/9.d0 -(76*lfr)/3.d0 
     &     + (-1.4814814814814814d0 - (8*lfh)/9.d0 
     &     + (8*lfr)/9.d0)*nf - (200
     &     *z2)/9.d0) +ddz0*(-14.962962962962964d0 - (11*lfh**2)/3.d0
     &     +nf*(0.691358024691358d0 + (2*lfh**2)/9.d0 +lfh*(0
     &     .7407407407407407d0 - (4*lfr)/9.d0) - (8*z2)/9.d0) +(44
     &     *z2)/3.d0 +lfh*(-11.333333333333334d0 + (38*lfr)/3.d0
     &     +(100*z2)/9.d0) + (382*z3)/9.d0)

      end

C-}}}
C-{{{ subroutine evalhard1bbh(...):

      subroutine evalhard1bbh(hard1,err1)
c..
c..   Hard contributions at NLO.
c..   
      implicit none
      real*8 hard1,err1,chi2a,ppall1bbh,ppbbq1bbh,ppbg1bbh
      include '../commons/common-keys.f'
      external ppall1bbh,ppbbq1bbh,ppbg1bbh

      if (nsubprocbbh.eq.0) then
         call convolute(2,ppall1bbh,hard1,err1,chi2a)
      elseif (nsubprocbbh.eq.1) then
         call convolute(2,ppbbq1bbh,hard1,err1,chi2a)
      elseif (nsubprocbbh.eq.2) then
         call convolute(2,ppbg1bbh,hard1,err1,chi2a)
      else
         hard1 = 0.d0
         err1 = 0.d0
      endif

      end
      
C-}}}
C-{{{ subroutine evalhard2bbh(...):

      subroutine evalhard2bbh(hard2,err2)
c..
c..   Hard contributions at NNLO.
c..   hardall2 contains the result obtained by summing the individual
c..   contributions BEFORE integration.
c..   It should be the same as the sum of hardgg2, hardqg2, etc.
c..   
      implicit none
      real*8 hard2,err2,chi2a,ppall2bbh,ppbbq2bbh,ppbg2bbh,ppgg2bbh
     &     ,ppbb2bbh,ppbq2bbh,ppqqb2bbh
      include '../commons/common-keys.f'
      external ppall2bbh,ppbbq2bbh,ppbg2bbh,ppgg2bbh
     &     ,ppbb2bbh,ppbq2bbh,ppqqb2bbh
      
      if (nsubprocbbh.eq.0) then
         call convolute(2,ppall2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.1) then
         call convolute(2,ppbbq2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.2) then
         call convolute(2,ppbg2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.3) then
         call convolute(2,ppgg2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.4) then
         call convolute(2,ppbb2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.5) then
         call convolute(2,ppbq2bbh,hard2,err2,chi2a)
      elseif (nsubprocbbh.eq.6) then
         call convolute(2,ppqqb2bbh,hard2,err2,chi2a)
      else
         hard2 = 0.d0
         err2 = 0.d0
      endif

      end
      
C-}}}
C-{{{ function sbbq1bbh(xt)

      real*8 function sbbq1bbh(xt)
C..
C..   sbbq1bbh(xt) is the one-loop result minus the purely soft terms
C..
      implicit none
      real*8 dlxm1,xm1,dlnx,xt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'

      dlxm1 = dlog(1-xt)
      xm1 = 1.d0-xt
      dlnx = dlog(xt)

      sbbq1bbh = (16*dlnx)/3.d0 - (32*dlxm1)/3.d0 + (16*lfh)/3.d0 
     & - (8*dlnx)/(3.d0*xm1) +(1.3333333333333333d0 - 4*dlnx + 8*dlxm1 
     &     - 4*lfh)*xm1 +(-1.3333333333333333d0 + (4*dlnx)/3.d0 - (8
     &     *dlxm1)/3.d0 + (4*lfh)/3.d0)*xm1**2

      end

C-}}}
C-{{{ function sbg1bbh(yy)

      real*8 function sbg1bbh(xt)
C..
C..   bg contribution at NLO, exact.
C..   
      implicit none
      real*8 xm1,dlxm1,dlnx,xt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)
      dlnx = dlog(xt)

      sbg1bbh = -dlnx/4.d0 + dlxm1/2.d0 - lfh/4.d0 + (0.5d0 + (3
     &     *dlnx)/4.d0 - (3*dlxm1)/2.d0 + (3*lfh)/4.d0)*xm1 + (
     &     -1.375d0 -dlnx + 2*dlxm1 - lfh)*xm1**2 + (0.875d0 + dlnx
     &     /2.d0 -dlxm1 + lfh/2.d0)*xm1**3

      end

C-}}}
C-{{{ function ppbbq1bbh(yy)

      real*8 function ppbbq1bbh(xx,wgt)
c..
c..   Integrand for hard bb-bar-contribution at NLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bbqpdf,sbbq1bbh
      external bbqpdf,sbbq1bbh

      ppbbq1bbh = vartrans(bbqpdf,sbbq1bbh,xx)

      end

C-}}}
C-{{{ function ppbg1bbh(yy)

      real*8 function ppbg1bbh(xx,wgt)
c..
c..   Integrand for hard bg-contribution at NLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bgpdf,sbg1bbh
      external bgpdf,sbg1bbh

      ppbg1bbh = vartrans(bgpdf,sbg1bbh,xx)

      end

C-}}}
C-{{{ function ppall1bbh(yy)

      real*8 function ppall1bbh(xx,wgt)
c..
c..   Sum of all exact integrands of the sub-processes at NLO.
c..
      implicit none
      real*8 xx(10),wgt,ppbbq1bbh,ppbg1bbh

      ppall1bbh = ppbbq1bbh(xx,wgt) + ppbg1bbh(xx,wgt)

      end

C-}}}
C-{{{ function ppbbq2bbh(yy)

      real*8 function ppbbq2bbh(xx,wgt)
c..
c..   Integrand for hard bb-bar-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bbqpdf,sbbq2bbh
      external bbqpdf,sbbq2bbh

      ppbbq2bbh = vartrans(bbqpdf,sbbq2bbh,xx)

      end

C-}}}
C-{{{ function ppbg2bbh(yy)

      real*8 function ppbg2bbh(xx,wgt)
c..
c..   Integrand for hard bg-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bgpdf,sbg2bbh
      external bgpdf,sbg2bbh

      ppbg2bbh = vartrans(bgpdf,sbg2bbh,xx)

      end

C-}}}
C-{{{ function ppgg2bbh(yy)

      real*8 function ppgg2bbh(xx,wgt)
c..
c..   Integrand for hard gg-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,ggpdf,sgg2bbh
      external ggpdf,sgg2bbh

      ppgg2bbh = vartrans(ggpdf,sgg2bbh,xx)

      end

C-}}}
C-{{{ function ppbb2bbh(yy)

      real*8 function ppbb2bbh(xx,wgt)
c..
c..   Integrand for hard bb-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bbpdf,sbb2bbh
      external bbpdf,sbb2bbh

      ppbb2bbh = vartrans(bbpdf,sbb2bbh,xx)

      end

C-}}}
C-{{{ function ppbq2bbh(yy)

      real*8 function ppbq2bbh(xx,wgt)
c..
c..   Integrand for hard bq-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,bqpdf,sbq2bbh
      external bqpdf,sbq2bbh

      ppbq2bbh = vartrans(bqpdf,sbq2bbh,xx)

      end

C-}}}
C-{{{ function ppqqb2bbh(yy)

      real*8 function ppqqb2bbh(xx,wgt)
c..
c..   Integrand for hard q-qbar contribution at NNLO.
c..
      implicit none
      real*8 xx(10),vartrans,wgt,qqbpdf,sqqb2bbh
      external qqbpdf,sqqb2bbh

      ppqqb2bbh = vartrans(qqbpdf,sqqb2bbh,xx)

      end

C-}}}
C-{{{ function ppall2bbh(yy)

      real*8 function ppall2bbh(xx,wgt)
c..
c..   Sum of all integrands of the sub-processes at NNLO.
c..
      implicit none
      real*8 xx(10),wgt,vartrans,ppbbq2bbh,ppbg2bbh,ppgg2bbh,ppbb2bbh
     &     ,ppbq2bbh,ppqqb2bbh

      ppall2bbh = ppbbq2bbh(xx,wgt) + ppbg2bbh(xx,wgt) 
     &     + ppgg2bbh(xx,wgt) + ppbb2bbh(xx,wgt) 
     &     + ppbq2bbh(xx,wgt) + ppqqb2bbh(xx,wgt)

      end

C-}}}
