C
C-{{{ subroutine intdel:

      subroutine intdel(del,errdel)
c..
c..   Integrating the delta(1-x) part over PDFs.
c..   
c..   del:    result
c..   errdel: uncertainty
c..   
      implicit real*8 (a-h,o-z)
      integer ndim,ncall,itmx,nprn
      include '../commons/common-vars.f'
      include '../commons/common-vegpar.f'
      common/bveg1/xl(10),xu(10),acc1,ndim,ncall,itmx,nprn
      external ggfun

      ndim=1
      nprn=nprnvggh
      acc1=acc

      do iv=2,10
         xl(iv)=0.d0
         xu(iv)=0.d0
      enddo

      itmx=itmx1ggh
      ncall=ncall1ggh
      xl(1) = tauh
      xu(1) = 1.d0
      call vegas(ggfun,del,errdel,chi2a)

      if (lveg1ggh) then
         itmx=itmx2ggh
         ncall=ncall2ggh
         call vegas1(ggfun,del,errdel,chi2a)
      endif

      end

C-}}}
C-{{{ function ggfun(yy):

      real*8 function ggfun(xt)
c..
c..   integrand for intdel
c..
      implicit real*8 (a-h,o-z)
      include '../commons/common-vars.f'
      external ggpdf

      ggfun = tauh * ggpdf(xt,tauh/xt)/xt
     
      return
      end

C-}}}
C-{{{ subroutine soft1(ddsoft1):

      subroutine soft1(ddsoft1,errsoft1)
c..
c..   Integrating the D-terms at NLO.
c..
      implicit real*8 (a-h,o-z)
      external ppdt1,ppdt1mt
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'

      call convolute(1,ppdt1,ddsoft1,errsoft1,chi2a)

      return
      end

C-}}}
C-{{{ subroutine soft2(ddsoft2):

      subroutine soft2(ddsoft2,errsoft2)
c..
c..   Integrating the D-terms at NNLO.
c..
      implicit real*8 (a-h,o-z)
      external ppdt2,ppdt2mt
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'
      include '../commons/common-citations.f'

      call convolute(1,ppdt2,ddsoft2,errsoft2,chi2a)

      return
      end

C-}}}
C-{{{ subroutine soft3(ddsoft3):

      subroutine soft3(ddsoft3,errsoft3)
c..
c..   Integrating the D-terms at N3LO.
c..
      implicit real*8 (a-h,o-z)
      external ppdt3,ppdt3mt
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'

      call convolute(1,ppdt3,ddsoft3,errsoft3,chi2a)

      return
      end

C-}}}
C-{{{ function ppdt1(yy)

      real*8 function ppdt1(xx)
c..
c..   Integrand for the D-terms at NLO.
c..   
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      external ggpdf,dterms1

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      ppdt1 = tauh * ( pmeas/zt**2 * ggpdf(ww/zt,tauh/ww) *zt**ncat1
     &     - ggpdf(xt,tauh/xt)/xt )*dterms1(zt)

      return
      end

C-}}}
C-{{{ function ppdt2(yy)

      real*8 function ppdt2(xx)
C
C     Integrand for D-terms at NNLO.
C
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      external ggpdf,dterms2

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      ppdt2 = tauh * ( pmeas/zt**2 * ggpdf(ww/zt,tauh/ww) *zt**ncat2
     &     - ggpdf(xt,tauh/xt)/xt )*dterms2(zt)

      return
      end

C-}}}
C-{{{ function ppdt3(yy)

      real*8 function ppdt3(xx)
C
C     Integrand for D-terms at N3LO.
C
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      external ggpdf,dterms3

      zt = xx(1)
      xt = xx(2)

      ww = ( (zt - tauh)*xt + tauh*(1.d0 - zt) )/(1.d0 - tauh)
      pmeas = (zt - tauh)/( (zt - tauh)*xt + tauh*(1.d0 - zt) )

      ppdt3 = tauh * ( pmeas/zt**2 * ggpdf(ww/zt,tauh/ww) *zt**ncat3
     &     - ggpdf(xt,tauh/xt)/xt )*dterms3(zt)

      return
      end

C-}}}
C-{{{ function delta1():

      real*8 function delta1(nmt)
c..
c..   Coefficient of the delta-function at NLO.
c..   
      implicit none
      integer nmt
      real*8 delta1mt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   pseudo-scalar:
      delta1 = 6.d0 + 6*z2

c..   lfr-terms:
      delta1 = delta1 - (11*lfr)/2.d0 + (lfr*nf)/3.d0

c..   switch to scalar:
      if (.not.lpseudo) then
         delta1 = delta1 - 0.5d0
      endif

c--   checked against checks.m [17/06/09,rh]

      if (nmt.gt.0) then
c..   add 1/mt terms:
         delta1 = delta1 + delta1mt(nmt)
      endif

c..   factor out SM Wilson coefficient:
      delta1 = delta1 - 2*c1sm1/c1sm0

      return
      end

C-}}}
C-{{{ function delta1mt():

      real*8 function delta1mt(nmt)
c..
c..   Coefficient of the delta-function at NLO, 1/mt-terms.
c..   
      implicit none
      integer nmt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   only scalar so far, and muF=muR
      if ((nmt.gt.0).and.lpseudo) stop 'delta1mt'
      if (nmt.gt.10) stop 'delta1mt'

      delta1mt = 0.d0

      if (nmt.ge.2) then
         delta1mt = delta1mt
     &        - (mh/mt)**2*(
     &        -193/216.d0 - (7*z2)/10.d0
     &        )
         delta1mt = delta1mt
     &        - lfr*(mh/mt)**2*(
     &        77/120.d0 - (7*nf)/180.d0)
      endif
      if (nmt.ge.4) then
         delta1mt = delta1mt
     &        + (mh/mt)**4*(
     &        29213/201600.d0 + (1543*z2)/16800.d0
     &        )
         delta1mt = delta1mt
     &        + lfr*(mh/mt)**4*(
     &        -16973/201600.d0 + (1543*nf)/302400.d0)
      endif
      if (nmt.ge.6) then
         delta1mt = delta1mt
     &        - (mh/mt)**6*(
     &        -4697767/190512000.d0 - (113*z2)/8400.d0
     &        )
         delta1mt = delta1mt
     &        + lfr*(mh/mt)**6*(
     &        -1243/100800.d0 + (113*nf)/151200.d0)
      endif
      if (nmt.ge.8) then
         delta1mt = delta1mt
     &        + (mh/mt)**8*(
     &        24571559/5588352000.d0 + (27677*z2)/12936000.d0
     &        )
         delta1mt = delta1mt
     &        + lfr*(mh/mt)**8*(
     &        -27677/14112000.d0 + (27677*nf)/232848000.d0)
      endif
      if (nmt.ge.10) then
         delta1mt = delta1mt
     &        - (mh/mt)**10*(
     &        -1086176593/1331890560000.d0 - (182653*z2)/504504000.d0
     &        )
         delta1mt = delta1mt
     &        + lfr*(mh/mt)**10*(
     &        -182653/550368000.d0 + (182653*nf)/9081072000.d0)
      endif

c--   checked against checks1l.m [19/06/09,rh]

      end

C-}}}
C-{{{ function delta2(nmt):

      real*8 function delta2(nmt)
c..
c..   Coefficient of the delta-function at NNLO.
c..   
      implicit none
      real*8 delta2mt,delta1
      integer nmt
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   pseudo-scalar:
      delta2 =  741/8.d0 - (689*nf)/72.d0 + lft*nf + (139*z2)/2.d0 - (5
     &     *nf*z2)/3.d0 -(165*z3)/4.d0 + (5*nf*z3)/6.d0 - (9*z4)/8.d0
      
      delta2 = delta2 -18*lfh**2*z2 + lfh*( 69/4.d0 - (9*nf)/4.d0 + (33
     &     *z2) /2.d0 - nf*z2 -(171*z3)/2.d0)

      delta2 = delta2 
     &     + lfr**2 * ( 363/16.d0 - 11/4.d0*nf + nf**2/12.d0 )
     &     + lfr * (-249/4.d0 + (55*nf)/12.d0 - (99*z2)/2.d0 + 3*nf*z2)

c..   switch to scalar:
      if (.not.lpseudo) then 
         delta2 = delta2 - (1939/144.d0 + 3*z2 + (15*lfh)/4.d0
     &        - (19*lft) /8.d0 
     &        + nf*( -21/16.d0 - 5/12.d0*lfh + 1/3.d0*lft) )
         delta2 = delta2 + lfr*(33/8.d0 - nf/4.d0)
      endif

c-- checked against checks.m [17 Jun 2009,rh]

      if (nmt.gt.0) then
         delta2 = delta2 + delta2mt(nmt)
      endif

c..   factor out SM Wilson coefficient:
      delta2 = delta2 - 2*c1sm1/c1sm0*delta1(nmt) - ( (c1sm1/c1sm0)
     &     **2 + 2*c1sm2/c1sm0 ) - c2sm1/c1sm0*nf
c..   last term is for pseudo-scalar

      return
      end

C-}}}
C-{{{ function delta2mt():

      real*8 function delta2mt(nmt)
c..
c..   Coefficient of the delta-function at NNLO, 1/mt-terms.
c..   
      implicit none
      integer i,nmt
      real*8 help(10),nl
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      if ((nmt.gt.0).and.lpseudo) stop 'delta2mt'
      if (nmt.gt.6) stop 'delta2mt'

      do i=1,10
         help(i) = 0.d0
      enddo

      nl=nf

c..   results in terms of on-shell top mass:

C-{{{ 1/mt^2:

      help(1) = -35947007/1244160.d0 - (4729*lth)/4320.d0 - (258031*nl)
     &     /388800.d0 +(19*lth*nl)/960.d0 +       (701*z2)/72.d0 + (7
     &     *ln2*z2)/45.d0- (49*nl*z2)/180.d0 - (21*z4)/160.d0 + (1643069
     &     *z3)/55296.d0 + (7*nl*z3)/72.d0
      help(1) = help(1) +(10993*lfh)/4320.d0 - (12757*lfr)/1440.d0 +
     &     (847*lfr**2)/320.d0 - (577*lfh*nl)/3240.d0 +  (341*lfr*nl)
     &     /540.d0- (77*lfr**2*nl)/240.d0 + (7*lfr**2*nl**2)/720.d0 +
     &     (77*lfh*z2) /40.d0 -  (21*lfh**2*z2)/10.d0 - (231*lfr*z2)
     &     /40.d0 - (7*lfh *nl*z2)/60.d0 + (7*lfr*nl*z2)/20.d0 -(399*lfh
     &     *z3)/40.d0

C-}}}
C-{{{ 1/mt^4:

      help(2) = -18366412471.d0/1560674304.d0 - (3757*lth)/17920.d0 -
     &     (115004203*nl)/1524096000.d0 +  (219577*lth*nl)/21772800.d0 +
     &     (303347*z2)/201600.d0 + (1543*ln2*z2)/37800.d0 -  (1543*nl
     &     *z2)/33600.d0 + (8135101621.d0*z3)/743178240.d0 + (1543*nl
     &     *z3) /120960.d0 -  (1543*z4)/89600.d0
      help(2) = help(2) + (330601*lfh)/806400.d0 - (74761*lfr)/53760.d0
     &     +(186703*lfr**2)/537600.d0 -  (16921*lfh*nl)/604800.d0 +
     &     (4177*lfr*nl)/43200.d0 - (16973*lfr**2*nl)/403200.d0 +  (1543
     &     *lfr**2*nl**2)/1209600.d0 + (16973*lfh*z2)/67200.d0 - (1543
     &     *lfh**2*z2)/5600.d0 -(16973*lfr*z2)/22400.d0 - (1543*lfh*nl
     &     *z2)/100800.d0 + (1543*lfr*nl*z2)/33600.d0 -  (29317*lfh*z3)
     &     /22400.d0

C-}}}
C-{{{ 1/mt^6:

      help(3) = -44.95627198396834d0 - 0.03848809510686991d0*lth -
     &     0.010305245049139095d0*nl + 0.002259751441728955d0*lth*nl +
     &     0.2499657344419249d0*z2 + 0.008968253968253969d0*ln2*z2 -
     &     0.00822089947089947d0*nl*z2 + 37.62443662184186d0*z3 +
     &     0.0018683862433862433d0*nl*z3 - 0.0025223214285714285d0*z4
      help(3) = help(3) + (0.06949281016418914d0*lfh -
     &     0.23202009715923405d0*lfr + 0.05086681547619048d0*lfr**2 -
     &     0.004670289360600207d0*lfh*nl + 0.015879254325186866d0*lfr*nl
     &     - 0.0061656746031746035d0*lfr**2*nl +
     &     0.00018683862433862433d0*lfr**2*nl**2 + 0.03699404761904762d0
     &     *lfh*z2 - 0.040357142857142855d0*lfh**2*z2 -
     &     0.11098214285714286d0*lfr*z2 - 0.0022420634920634923d0*lfh*nl
     &     *z2 + 0.006726190476190476d0*lfr*nl*z2 -
     &     0.19169642857142857d0*lfh*z3)

C-}}}

      delta2mt = 0.d0

      do i=1,nmt/2
         delta2mt = delta2mt + (mh/mt)**(2*i)* help(i)
      enddo

      return
      end

C-}}}
C-{{{ function delta3():

      real*8 function delta3()
c..
c..   Coefficient of the delta-function at NNNLO.
c..   
      implicit none
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

c.. all formulas obtained with ~/math/ggh/rgfull-n3lo.m:         

c..   mu-independent terms:
      delta3 = 1120.4739583333333d0 - 217.74054783950618d0*nf +
     &     6.671360596707819d0*nf**2 + 336.4791666666667d0*z2 -
     &     15.666666666666666d0*nf*z2 - 1.6574074074074074d0*nf**2*z2 -
     &     1382.0625d0*z3 + 82.77314814814815d0*nf*z3 -
     &     0.3611111111111111d0 *nf**2*z3 + 1101.375d0*z2*z3 - 81.75d0
     &     *nf*z2*z3 + 1858.5d0*z3**2 - 476.78125d0*z4 +
     &     55.38194444444444d0*nf*z4 - 1.6203703703703705d0 *nf**2*z4 -
     &     1421.0625d0*z5 + 109.73611111111111d0*nf*z5 - 1126.6875d0*z6
      
c..   mu-dependent terms
      delta3 = delta3 + (393.0625d0*lfh + 37.125d0*lfh**2 - 714.109375d0
     &     *lfr -148.5d0*lfh*lfr + 227.90625d0*lfr**2 - 83.1875d0*lfr**3
     &     -74.79513888888889d0*lfh*nf - 7.291666666666667d0*lfh**2*nf
     &     +135.22743055555554d0*lfr*nf + 29.166666666666668d0*lfh*lfr
     &     *nf -42.114583333333336d0*lfr**2*nf + 15.125d0*lfr**3*nf
     &     +2.6493055555555554d0*lfh*nf**2 + 0.3055555555555556d0*lfh**2
     &     *nf**2 - 5.138310185185185d0*lfr*nf**2 - 1.2222222222222223d0
     &     *lfh*lfr*nf**2 + 1.7152777777777777d0*lfr**2*nf**2
     &     -0.9166666666666666d0*lfr**3*nf**2 + 0.018518518518518517d0
     &     *lfr**3*nf**3 + 99.25d0*lfh*z2 - 155.625d0*lfh**2*z2 - 49.5d0
     &     *lfh**3*z2 -559.75d0*lfr*z2 - 181.5d0*lfh*lfr*z2 + 198.d0*lfh
     &     **2*lfr*z2 +272.25d0*lfr**2*z2 - 26.d0*lfh*nf*z2 + 4.5d0*lfh
     &     **2*nf*z2 + 3.d0*lfh**3*nf*z2 + 64.41666666666667d0*lfr*nf*z2
     &     + 22.d0*lfh*lfr*nf*z2 -12.d0*lfh**2*lfr*nf*z2 - 33.d0*lfr**2
     &     *nf*z2 + 0.5555555555555556d0*lfh*nf**2*z2 +
     &     0.16666666666666666d0*lfh**2*nf**2*z2 -1.1111111111111112d0
     &     *lfr*nf**2*z2 - 0.6666666666666666d0*lfh*lfr*nf**2*z2 + lfr
     &     **2*nf**2*z2 - 1181.625d0*lfh*z3 - 334.125d0*lfh**2*z3 -
     &     72.d0*lfh**3*z3 + 453.75d0*lfr*z3 + 940.5d0*lfh*lfr*z3
     &     +65.83333333333333d0*lfh*nf*z3 + 20.25d0*lfh**2*nf*z3
     &     -36.666666666666664d0*lfr*nf*z3 - 57.d0*lfh*lfr*nf*z3
     &     -0.2777777777777778d0*lfh*nf**2*z3 + 0.5555555555555556d0*lfr
     &     *nf**2*z3 + 1633.5d0*lfh*z2*z3 + 77.34375d0*lfh*z4 - 243.d0
     &     *lfh**2*z4+ 12.375d0*lfr*z4 - 4.6875d0*lfh*nf*z4 - 0.75d0*lfr
     &     *nf*z4 - 2524.5d0*lfh*z5)

c..   only scalar so far:
      if (lpseudo) 
     &     call printdieggh('pseudo scalar Higgs not available at N3LO')

      end

C-}}}
C-{{{ function delta3mt():

      real*8 function delta3mt()
c..
c..   Coefficient of the delta-function at NNLO, 1/mt-terms.
c..   
      implicit real*8 (a-h,o-z)
      real*8 help(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      
      write(6,*) 'not implemented yet'
      stop
      delta3mt = 0.d0

      return
      end

C-}}}
C-{{{ function dterms1(...):

      real*8 function dterms1(xt)
c..
c..   The plus-distributions at NLO.
c..   
      implicit none
      include '../commons/common-expand.f'
      real*8 dd0,dd1,xt,dcoefs1

      dd0 = 1.d0/(1.d0 - xt)
      dd1 = dd0*dlog(1.d0 - xt)

      dterms1 = dcoefs1(nmtlim(1,1),dd0,dd1)

      end

C-}}}
C-{{{ function dtsub1(...):

      real*8 function dtsub1(nmt)
c..
c..   Contributions arising from the fact that the integrals over
c..   plus-distributions do not run from 0 to 1, but from z to 1.
c..   
      implicit none
      real*8 ddz0,ddz1,dcoefs1
      integer nmt
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'

      ddz0 = dlog(1.d0 - tauh)
      ddz1 = ddz0**2/2.d0

      dtsub1 = dcoefs1(nmt,ddz0,ddz1)

      return
      end

C-}}}
C-{{{ function dcoefs1(...):

      real*8 function dcoefs1(nmt,dd0,dd1)
c..
c..   The plus-distributions at NLO.
c..   
      implicit none
      integer nmt
      real*8 dd0,dd1,dcoefs1mt
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

c..   pseudo-scalar and scalar are the same:
      dcoefs1 = - 6.d0*dd0*lfh + 12.d0*dd1

      if (nmt.gt.0) dcoefs1 = dcoefs1 + dcoefs1mt(nmt,dd0,dd1)

c--   checked against checks.m [17/06/09,rh]

      end

C-}}}
C-{{{ function dcoefs1mt(...):

      real*8 function dcoefs1mt(nmt,dd0,dd1)
c..
c..   The plus-distributions at NLO.
c..   
      implicit none
      integer nmt
      real*8 dd0,dd1
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'

c..   only scalar so far
      if ((nmt.gt.0).and.lpseudo) stop 'dcoefs1mt'
      if (nmt.gt.10) stop 'dcoefs1mt'

      dcoefs1mt = 0.d0
      if (nmt.ge.2) then
         dcoefs1mt = dcoefs1mt 
     &        - (mh/mt)**2*(
     &        (-7.d0*dd1)/5.d0 + (7.d0*dd0*lfh)/10.d0
     &        )
      endif
      if (nmt.ge.4) then
         dcoefs1mt = dcoefs1mt 
     &        + (mh/mt)**4*(
     &         (1543*dd1)/8400.d0 - (1543*dd0*lfh)/16800.d0
     &        )
      endif
      if (nmt.ge.6) then
         dcoefs1mt = dcoefs1mt 
     &        - (mh/mt)**6*(
     &        (-113*dd1)/4200.d0+ (113*dd0*lfh)/8400.d0
     &        )
      endif
      if (nmt.ge.8) then
         dcoefs1mt = dcoefs1mt 
     &        + (mh/mt)**8*(
     &        (27677*dd1)/6468000.d0- (27677*dd0*lfh)/12936000.d0
     &        )
      endif
      if (nmt.ge.10) then
         dcoefs1mt = dcoefs1mt 
     &        - (mh/mt)**10*(
     &        (-182653*dd1)/252252000.d0+ (182653*dd0*lfh)/504504000.d0
     &        )
      endif

c--   checked against checks.m [18/06/09,rh]

      end

C-}}}
C-{{{ function dterms2(...):

      real*8 function dterms2(xt)

      implicit none
      include '../commons/common-expand.f'
      real*8 dd0,dd1,dd2,dd3,xt,dcoefs2

      dd0 = 1.d0/(1.d0 - xt)
      dd1 = dd0*dlog(1.d0 - xt)
      dd2 = dd1*dlog(1.d0 - xt)
      dd3 = dd2*dlog(1.d0 - xt)

      dterms2 = dcoefs2(nmtlim(2,1),dd0,dd1,dd2,dd3)

      return
      end

C-}}}
C-{{{ function dtsub2(...):

      real*8 function dtsub2(nmt)

      implicit none
      integer nmt
      include '../commons/common-vars.f'
      real*8 ddz0,ddz1,ddz2,ddz3,dcoefs2

      ddz0 = dlog(1.d0 - tauh)
      ddz1 = ddz0**2/2.d0
      ddz2 = ddz0**3/3.d0
      ddz3 = ddz0**4/4.d0

      dtsub2 = dcoefs2(nmt,ddz0,ddz1,ddz2,ddz3)
      
      return
      end

C-}}}
C-{{{ function dcoefs2(...):

      real*8 function dcoefs2(nmt,dd0,dd1,dd2,dd3)

      implicit none
      integer nmt
      real*8 dd0,dd1,dd2,dd3,dcoefs2mt,dcoefs1
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   pseudo-scalar:
      dcoefs2 = (-101/3.d0 + (14*nf)/9.d0 + 33*z2 - 2*nf*z2
     &     +(351*z3)/2.d0)*dd0 + (139 - (10*nf)/3.d0 - 90*z2)*dd1 +
     &     (-33+ 2*nf)*dd2 +72*dd3
      
      dcoefs2 = dcoefs2 + (lfh**2*(-33/4.d0 + nf/2.d0) + lfh*(-139/2.d0
     &     +(5*nf)/3.d0 + 45*z2))*dd0 + (36*lfh**2 + lfh*(33 - 2*nf))
     &     *dd1- 108*lfh*dd2

      dcoefs2 = dcoefs2 + lfh*lfr*(99/2.d0 - 3*nf)*dd0 + lfr*(-99 + 6
     &     *nf)*dd1

c..   switch to scalar:
      if (.not.lpseudo) then
         dcoefs2 = dcoefs2 - ( -3*lfh*dd0 + 6*dd1 )
      endif

      if (nmt.gt.0) dcoefs2 = dcoefs2 + dcoefs2mt(nmt,dd0,dd1,dd2
     &     ,dd3)

c..   factor out SM Wilson coefficient:
      dcoefs2 = dcoefs2 - 2*c1sm1/c1sm0*dcoefs1(nmt,dd0,dd1)

c--   checked against checks.m [17/06/2009,rh]

      end

C-}}}
C-{{{ function dcoefs2mt(...):

      real*8 function dcoefs2mt(nmt,dd0,dd1,dd2,dd3)

      implicit none
      integer i,nmt
      real*8 dd0,dd1,dd2,dd3,nl,help(10)      
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'

c..   only scalar so far
      if (lpseudo) stop 'delta2mt'
      if (nmt.gt.6) stop 'delta2mt'

      do i=1,10
         help(i) = 0.d0
      enddo

      nl=nf

C-{{{ 1/mt^2:

      help(1) = (42*dd3)/5.d0 + dd2*(-77/20.d0 + (7*nl)/30.d0) +  dd1
     &     *(3337/180.d0 -(7*nl)/18.d0 - (21*z2)/2.d0) +  dd0*(-707
     &     /180.d0 + (49*nl)/270.d0 + (77*z2)/20.d0 - (7*nl*z2)/30.d0 +
     &     (819*z3)/40.d0)
      help(1) = help(1) +(-63*dd2*lfh)/5.d0 + dd1*((77*lfh)/20.d0 + (21
     &     *lfh**2)/5.d0 - (231*lfr)/20.d0 -    (7*lfh*nl)/30.d0 + (7
     &     *lfr *nl)/10.d0) +  dd0*((-3337*lfh)/360.d0 - (77*lfh**2)
     &     /80.d0 + (231*lfh*lfr)/40.d0 + (7*lfh*nl)/36.d0 +    (7*lfh
     &     **2*nl) /120.d0 - (7*lfh*lfr*nl)/20.d0 + (21*lfh*z2)/4.d0)

C-}}}
C-{{{ 1/mt^4:

      help(2) = (1543*dd3)/1400.d0 + dd2*(-16973/33600.d0 + (1543*nl)
     &     /50400.d0) +dd1*(278659/100800.d0 - (1543*nl)/30240.d0 -
     &     (1543*z2)/1120.d0) +  dd0*(-155843/302400.d0 + (1543*nl)
     &     /64800.d0 + (16973*z2)/33600.d0 - (1543*nl*z2)/50400.d0 +
     &     (60177*z3)/22400.d0)
      help(2) = help(2) +(-4629*dd2*lfh)/2800.d0 + dd1*((16973*lfh)
     &     /33600.d0 +(1543*lfh**2)/2800.d0 -    (16973*lfr)/11200.d0 -
     &     (1543*lfh*nl)/50400.d0 + (1543*lfr*nl)/16800.d0) +  dd0*((
     &     -278659*lfh)/201600.d0 -(16973*lfh**2)/134400.d0 + (16973*lfh
     &     *lfr)/22400.d0 +    (1543*lfh*nl)/60480.d0 + (1543*lfh**2*nl)
     &     /201600.d0 - (1543*lfh*lfr*nl)/33600.d0+    (1543*lfh*z2)
     &     /2240.d0)

C-}}}
C-{{{ 1/mt^6:

      help(3) = -0.07548280423280423d0*dd0 + 0.446121945074326d0*dd1 -
     &     0.07398809523809524d0*dd2 + 0.16142857142857142d0*dd3 +
     &     0.003487654320987654d0*dd0*nl - 0.007473544973544973d0*dd1*nl
     &     + 0.0044841269841269845d0*dd2*nl + 0.07398809523809524d0*dd0
     &     *z2 - 0.2017857142857143d0*dd1*z2 - 0.0044841269841269845d0
     &     *dd0*nl*z2 + 0.39348214285714284d0*dd0*z3
      help(3) = help(3) + (-0.223060972537163d0*dd0*lfh +
     &     0.07398809523809524d0*dd1*lfh - 0.24214285714285713d0*dd2*lfh
     &     - 0.01849702380952381d0*dd0*lfh**2 + 0.08071428571428571d0
     &     *dd1*lfh**2 - 0.22196428571428573d0*dd1*lfr +
     &     0.11098214285714286d0*dd0*lfh*lfr + 0.0037367724867724867d0
     &     *dd0*lfh*nl - 0.0044841269841269845d0*dd1*lfh*nl +
     &     0.0011210317460317461d0*dd0*lfh**2*nl +
     &     0.013452380952380952d0*dd1*lfr*nl - 0.006726190476190476d0
     &     *dd0*lfh*lfr*nl + 0.10089285714285715d0*dd0*lfh*z2)

C-}}}

      dcoefs2mt = 0.d0
      do i=1,nmt/2
         dcoefs2mt = dcoefs2mt + (mh/mt)**(2*i) * help(i)
      enddo

      return
      end

C-}}}
C-{{{ function dterms3(...):

      real*8 function dterms3(xt)

      implicit real*8 (a-h,o-z)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      dd0 = 1.d0/(1.d0 - xt)
      dd1 = dd0*dlog(1.d0 - xt)
      dd2 = dd1*dlog(1.d0 - xt)
      dd3 = dd2*dlog(1.d0 - xt)
      dd4 = dd3*dlog(1.d0 - xt)
      dd5 = dd4*dlog(1.d0 - xt)  

      dterms3 = dcoefs3(dd0,dd1,dd2,dd3,dd4,dd5)

      end

C-}}}
C-{{{ function dtsub3(...):

      real*8 function dtsub3()

      implicit real*8 (a-h,o-z)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      ddz0 = dlog(1.d0 - tauh)
      ddz1 = ddz0**2/2.d0
      ddz2 = ddz0**3/3.d0
      ddz3 = ddz0**4/4.d0
      ddz4 = ddz0**5/5.d0
      ddz5 = ddz0**6/6.d0

      dtsub3 = dcoefs3(ddz0,ddz1,ddz2,ddz3,ddz4,ddz5)

      end

C-}}}
C-{{{ function dcoefs3(...):

c..
c..   The plus-distributions at NNNLO.
c..   
      real*8 function dcoefs3(dd0,dd1,dd2,dd3,dd4,dd5)

      implicit none
      real*8 dd0,dd1,dd2,dd3,dd4,dd5
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c.. all formulas obtained with ~/math/ggh/rgfull-n3lo.m:         
c..   mu-independent terms:
      dcoefs3 = 0.d0

      dcoefs3 = dcoefs3 + dd0*(-343.78356481481484d0 +
     &     32.08256172839506d0*nf - 0.23868312757201646d0*nf**2 +
     &     713.5833333333334d0*z2 -62.361111111111114d0*nf*z2 +
     &     1.1111111111111112d0*nf**2*z2 +2235.25d0*z3 -
     &     122.97222222222223d0*nf*z3 + 0.5555555555555556d0*nf**2*z3 -
     &     3262.5d0*z2*z3 + 284.625d0*z4 - 14.75d0*nf*z4 + 5022.d0*z5)

      dcoefs3 = dcoefs3 + dd1*(1273.7083333333333d0 -
     &     148.47222222222223d0*nf + 0.9259259259259259d0*nf**2 -
     &     1368.d0*z2 + 94.d0*nf*z2 -1.3333333333333333d0*nf**2*z2 -
     &     3168.d0*z3 + 162.d0*nf*z3 - 2079.d0*z4)

      dcoefs3 = dcoefs3 + dd2*(-1051.d0 + 78.16666666666667d0*nf -
     &     1.1111111111111112d0*nf**2 + 1683.d0*z2 - 102.d0*nf*z2 +
     &     4887.d0*z3)

      dcoefs3 = dcoefs3 + dd3*(925.d0 - 54.666666666666664d0*nf +
     &     0.4444444444444444d0*nf**2 - 1512.d0*z2)

      dcoefs3 = dcoefs3 + dd4*(-330.d0 + 20.d0*nf)

      dcoefs3 = dcoefs3 + dd5*(216.d0)

c.. mu-dependent terms:
      dcoefs3 = dcoefs3 + dd0*(lfh*lfr**2*(-272.25d0 + 33.d0*nf - 1.d0
     &     *nf**2)+ lfh**3*(-15.125d0 + 1.8333333333333333d0*nf
     &     -0.05555555555555555d0*nf**2 + 108.d0*z2) + lfh**2*(-192.25d0
     &     +23.541666666666668d0*nf - 0.2777777777777778d0*nf**2 +
     &     173.25d0*z2- 10.5d0*nf*z2 + 945.d0*z3) + lfr
     &     *(370.3333333333333d0 -39.55555555555556d0*nf +
     &     1.037037037037037d0*nf**2 + lfh**2*(90.75d0 - 11.d0*nf +
     &     0.3333333333333333d0*nf**2) - 363.d0*z2 + 44.d0*nf*z2 -
     &     1.3333333333333333d0*nf**2*z2 + lfh*(559.75d0
     &     -64.41666666666667d0*nf + 1.1111111111111112d0*nf**2 - 495.d0
     &     *z2 +30.d0*nf*z2) - 1930.5d0*z3 + 117.d0*nf*z3) + lfh*(
     &     -636.8541666666666d0 + 74.23611111111111d0*nf
     &     -0.46296296296296297d0*nf**2 + 684.d0*z2 - 47.d0*nf*z2
     &     +0.6666666666666666d0*nf**2*z2 + 1584.d0*z3 - 81.d0*nf*z3 +
     &     1039.5d0*z4))

      dcoefs3 = dcoefs3 + dd1*(lfh**3*(99.d0 - 6.d0*nf) + lfr**2
     &     *(544.5d0 -66.d0*nf + 2.d0*nf**2) + lfh**2*(492.75d0 - 31.d0
     &     *nf +0.3333333333333333d0*nf**2 - 972.d0*z2) + lfr*(-1119.5d0
     &     +128.83333333333334d0*nf - 2.2222222222222223d0*nf**2 + lfh
     &     **2*(-396.d0 + 24.d0*nf) + lfh*(-363.d0 + 44.d0*nf -
     &     1.3333333333333333d0*nf**2) + 990.d0*z2 - 60.d0*nf*z2) + lfh
     &     *(1011.d0 - 90.83333333333333d0*nf + 1.1111111111111112d0*nf
     &     **2 - 1287.d0*z2 + 78.d0*nf*z2 - 4860.d0*z3))

      dcoefs3 = dcoefs3 + dd2*(-108.d0*lfh**3 + lfh**2*(-445.5d0 + 27.d0
     &     *nf) +lfr*(363.d0 + lfh*(1188.d0 - 72.d0*nf) - 44.d0*nf
     &     +1.3333333333333333d0*nf**2) + lfh*(-1387.5d0 + 82.d0*nf
     &     -0.6666666666666666d0*nf**2 + 2268.d0*z2))

      dcoefs3 = dcoefs3 + dd3*(432.d0*lfh**2 + lfh*(660.d0 - 40.d0*nf) +
     &     lfr*(-792.d0 + 48.d0*nf))

      dcoefs3 = dcoefs3 + dd4*(-540.d0*lfh)
      
      dcoefs3 = dcoefs3 + dd5*(0.d0)

c.. for checking:
c$$$      print*,'ddrep = {Dd[0] -> Rationalize[',dd0,',0],'
c$$$      print*,'Dd[1] -> Rationalize[',dd1,',0],'
c$$$      print*,'Dd[2] -> Rationalize[',dd2,',0],'
c$$$      print*,'Dd[3] -> Rationalize[',dd3,',0],'
c$$$      print*,'Dd[4] -> Rationalize[',dd4,',0],'
c$$$      print*,'Dd[5] -> Rationalize[',dd5,',0]}'
c$$$      print*,dcoefs3
c$$$      stop

c..   only scalar so far:
      if (lpseudo) then
         call printdieggh('pseudo scalar Higgs not available at N3LO')
      endif

      end

C-}}}
C-{{{ function dcoefs3mt(...):

      real*8 function dcoefs3mt(dd0,dd1,dd2,dd3,dd4,dd5)

      implicit real*8 (a-h,o-z)
      real*8 help(10)
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'

      write(6,*) 'not implemented yet'
      stop
      dcoefs3mt = 0.d0

      return
      end

C-}}}
C-{{{ subroutine evalhard1(...):

      subroutine evalhard1(hard1,err1)
c..
c..   Hard contributions at NLO.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard1(-1:5),err1(-1:5)
      external ppall1,ppgg1,ppqg1,ppqqb1
      external ppall1x0,ppgg1x0,ppqg1x0,ppqqb1x0
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'
      include '../commons/common-citations.f'
      
      hard1=0.d0
      err1=0.d0

      if (lball(1)) then
         citations(38) = 1
         call dohard1(ppall1x0,ppgg1x0,ppqg1x0,ppqqb1x0,hard1
     &        ,err1)
      else
         call dohard1(ppall1,ppgg1,ppqg1,ppqqb1,hard1,err1)
      endif

      return
      end
      
C-}}}
C-{{{ subroutine dohard1(...):

      subroutine dohard1(fppall1,fppgg1,fppqg1,fppqqb1,hard1,err1)
c..
c..   Hard contributions at NLO.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard1(-1:5),err1(-1:5)
      external fppall1,fppgg1,fppqg1,fppqqb1
      include '../commons/common-keys.f'

      do i=-1,5
         hard1(i) = 0.d0
         err1(i) = 0.d0
      enddo

      if (nsubprocggh.eq.0) then
         call convolute(1,fppall1,hard1(0),err1(0),chi2a)
      else
         if (nsubprocggh.eq.1) then
            call convolute(1,fppgg1,hard1(1),err1(1),chi2a)
         elseif (nsubprocggh.eq.2) then
            call convolute(1,fppqg1,hard1(2),err1(2),chi2a)
         elseif (nsubprocggh.eq.3) then
            call convolute(1,fppqqb1,hard1(3),err1(3),chi2a)
         elseif ((nsubprocggh.eq.4).or.(nsubprocggh.eq.5)) then
         elseif (nsubprocggh.eq.10) then
            call convolute(1,fppgg1,hard1(1),err1(1),chi2a)
            call convolute(1,fppqg1,hard1(2),err1(2),chi2a)
            call convolute(1,fppqqb1,hard1(3),err1(3),chi2a)
         else
            call printdieggh('no such subprocess')
         endif
         hard1(0) = hard1(1) + hard1(2) + hard1(3)
         err1(0) = dsqrt(err1(1)**2 + err1(2)**2 + err1(3)**2)
      endif

      return
      end
      
C-}}}
C-{{{ subroutine evalhard2(...):

      subroutine evalhard2(hard2,err2)
c..
c..   Hard contributions at NNLO.
c..   hardall2 contains the result obtained by summing the individual
c..   contributions BEFORE integration.
c..   It should be the same as the sum of hardgg2, hardqg2, etc.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard2(-1:5),err2(-1:5)
      external ppall2,ppgg2,ppqg2,ppqqb2,ppqq2,ppqu2
      external ppall2x0,ppgg2x0,ppqg2x0,ppqqb2x0,ppqq2x0,ppqu2x0
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      include '../commons/common-citations.f'

      hard2 = 0.d0
      err2 = 0.d0

      if (lball(2)) then
         citations(38) = 1
         call dohard2(ppall2x0,ppgg2x0,ppqg2x0,ppqqb2x0,ppqq2x0,ppqu2x0
     &        ,hard2,err2)
      else
         call dohard2(ppall2,ppgg2,ppqg2,ppqqb2,ppqq2,ppqu2,hard2
     &        ,err2)
      endif

      return
      end
      
C-}}}
C-{{{ subroutine dohard2(...):

      subroutine dohard2(fppall2,fppgg2,fppqg2,fppqqb2,fppqq2,fppqu2
     &     ,hard2,err2)
c..
c..   Hard contributions at NNLO.
c..   hardall2 contains the result obtained by summing the individual
c..   contributions BEFORE integration.
c..   It should be the same as the sum of hardgg2, hardqg2, etc.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard2(-1:5),err2(-1:5)
      external fppall2,fppgg2,fppqg2,fppqqb2,fppqq2,fppqu2
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      do i=-1,5
         hard2(i) = 0.d0
         err2(i) = 0.d0
      enddo

      if (nsubprocggh.eq.0) then
         call convolute(1,fppall2,hard2(0),err2(0),chi2a)
      else
         if (nsubprocggh.eq.1) then
            call convolute(1,fppgg2,hard2(1),err2(1),chi2a)
         elseif (nsubprocggh.eq.2) then
            call convolute(1,fppqg2,hard2(2),err2(2),chi2a)
         elseif (nsubprocggh.eq.3) then
            call convolute(1,fppqqb2,hard2(3),err2(3),chi2a)
         elseif (nsubprocggh.eq.4) then
            call convolute(1,fppqq2,hard2(4),err2(4),chi2a)
         elseif (nsubprocggh.eq.5) then
            call convolute(1,fppqu2,hard2(5),err2(5),chi2a)
         elseif (nsubprocggh.eq.10) then
            call convolute(1,fppgg2,hard2(1),err2(1),chi2a)
            call convolute(1,fppqg2,hard2(2),err2(2),chi2a)
            call convolute(1,fppqqb2,hard2(3),err2(3),chi2a)
            call convolute(1,fppqq2,hard2(4),err2(4),chi2a)
            call convolute(1,fppqu2,hard2(5),err2(5),chi2a)
         else
            call printdieggh('no such subprocess')
         endif
         hard2(0) = hard2(1) + hard2(2) + hard2(3) + hard2(4) +
     &        hard2(5)
         err2(0) = dsqrt( err2(1)**2 + err2(2)**2 + err2(3)**2 +
     &        err2(4)**2 + err2(5)**2 )
      endif

      return
      end
      
C-}}}
C-{{{ subroutine evalhard3(...):

      subroutine evalhard3(hard3,err3)
c..
c..   Hard contributions at N3LO.
c..   hardall3 contains the result obtained by summing the individual
c..   contributions BEFORE integration.
c..   It should be the same as the sum of hardgg3, hardqg3, etc.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard3(-1:5),err3(-1:5)
      external ppall3,ppgg3,ppqg3,ppqqb3,ppqq3,ppqu3
      external ppall3x0,ppgg3x0,ppqg3x0,ppqqb3x0,ppqq3x0,ppqu3x0
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      include '../commons/common-citations.f'

      hard3 = 0.d0
      err3 = 0.d0

      if (lball(3)) then
c$$$         call dohard3(ppall3x0,ppgg3x0,ppqg3x0,ppqqb3x0,ppqq3x0,ppqu3x0
c$$$     &        ,hard3,err3)

         call dohard3(ppall3,ppgg3x0,ppqg3,ppqqb3,ppqq3,ppqu3
     &        ,hard3,err3)
      else
         call dohard3(ppall3,ppgg3,ppqg3,ppqqb3,ppqq3,ppqu3,hard3
     &        ,err3)
      endif

      return
      end
      
C-}}}
C-{{{ subroutine dohard3(...):

      subroutine dohard3(fppall3,fppgg3,fppqg3,fppqqb3,fppqq3,fppqu3
     &     ,hard3,err3)
c..
c..   Hard contributions at NNLO.
c..   hardall3 contains the result obtained by summing the individual
c..   contributions BEFORE integration.
c..   It should be the same as the sum of hardgg3, hardqg3, etc.
c..   
      implicit real*8 (a-h,o-z)
      real*8 hard3(-1:5),err3(-1:5)
      external fppall3,fppgg3,fppqg3,fppqqb3,fppqq3,fppqu3
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      do i=-1,5
         hard3(i) = 0.d0
         err3(i) = 0.d0
      enddo

      if (nsubprocggh.eq.0) then
         call convolute(1,fppall3,hard3(0),err3(0),chi2a)
      else
         if (nsubprocggh.eq.1) then
            call convolute(1,fppgg3,hard3(1),err3(1),chi2a)
         elseif (nsubprocggh.eq.2) then
            call convolute(1,fppqg3,hard3(2),err3(2),chi2a)
         elseif (nsubprocggh.eq.3) then
            call convolute(1,fppqqb3,hard3(3),err3(3),chi2a)
         elseif (nsubprocggh.eq.4) then
            call convolute(1,fppqq3,hard3(4),err3(4),chi2a)
         elseif (nsubprocggh.eq.5) then
            call convolute(1,fppqu3,hard3(5),err3(5),chi2a)
         elseif (nsubprocggh.eq.10) then
            call convolute(1,fppgg3,hard3(1),err3(1),chi2a)
            call convolute(1,fppqg3,hard3(2),err3(2),chi2a)
            call convolute(1,fppqqb3,hard3(3),err3(3),chi2a)
            call convolute(1,fppqq3,hard3(4),err3(4),chi2a)
            call convolute(1,fppqu3,hard3(5),err3(5),chi2a)
         else
            call printdieggh('no such subprocess')
         endif
         hard3(0) = hard3(1) + hard3(2) + hard3(3) + hard3(4) +
     &        hard3(5)
         err3(0) = dsqrt( err3(1)**2 + err3(2)**2 + err3(3)**2 +
     &        err3(4)**2 + err3(5)**2 )
      endif

      return
      end
      
C-}}}
C-{{{ function sgg1(xt)

      real*8 function sgg1(xt)
C..
C..   sgg1(xt) is the one-loop result minus the purely soft terms, i.e.
C..
C..   full = sig0*( delta(1-z)*( 11/2 + 6*z2 )  +  12*D1  +  sgg1 )
C..
      implicit none
      include '../commons/common-expand.f'
      real*8 xt,sgg1exact,sgg1exp

      if (nexpand(1).eq.0) then
         sgg1 = sgg1exact(ncat1,nmtlim(1,1),xt)
      else
         sgg1 = sgg1exp(nmtlim(1,1),nsoft(1,1),ncat1,xt)
      endif
      
      return
      end

C-}}}
C-{{{ function sgg1exact(xt)

      real*8 function sgg1exact(ncat,nmt,xt)
C..
C..   sgg1exact(xt) is the one-loop result minus the purely soft terms, i.e.
C..
C..   full = sig0*( delta(1-z)*( 11/2 + 6*z2 )  +  12*D1  +  sgg1exact )
C..
      implicit none
      real*8 dlxm1,xm1,xt,dlnx,sgg1mtexact,dcoefs1,factor
      integer k,nmt,ncat
c..      implicit real*8 (a-h,o-z)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      dlxm1 = dlog(1-xt)
      xm1 = 1.d0-xt
      dlnx = dlog(xt)

c..   pseudo-scalar and scalar are the same:
      sgg1exact = 12*dlnx - 24*dlxm1 - (6*dlnx)/xm1 - 18*dlnx*xm1 + 36
     &     *dlxm1*xm1 +12*dlnx*xm1**2 - 24*dlxm1*xm1**2 - (11*xm1**3)
     &     /2.d0 - 6*dlnx*xm1**3 +12*dlxm1*xm1**3
      
      sgg1exact = sgg1exact + lfh*(12 - 18*xm1 + 12*xm1**2 - 6*xm1**3)
c..   there are no lfr-terms

      if (nmt.gt.0) sgg1exact = sgg1exact + sgg1mtexact(nmt,xt)

c..   factor = (1+xt+xt^2+...+xt^(n-1)) = (1-xt^n)/(1-xt)
      factor = 0.d0
      do k=0,ncat-1
         factor = factor + xt**k
      enddo
      sgg1exact = sgg1exact + factor*dcoefs1(nmt,1.d0,dlxm1)

c--   checked against checks.m [17/06/09,rh]

      return
      end

C-}}}
C-{{{ function sgg1exp(xt)

      real*8 function sgg1exp(nmt,nsft,ncat,xt)
C..
C..   sgg1exp(xt) is sgg1(xt) expanded in xm1
C..
      implicit none
      integer i,k,nmt,nsft,ncat
      real*8 xt,xm1,dlxm1,sum,dcoefs1,bininv,sgg1mtexp
      real*8 help(0:20),helpx(0:20),helpd(0:20)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)
      
      if (nsft.gt.16) stop 'sgg1exp(): nsft too large'
c..
c..   input produced with ~/math/gghmt/rgfull.m,
c..   function printitnlo[a,b]
c..
      help(0) =  6 - 24*dlxm1 + 12*lfh
      help(1) = -9 + 36*dlxm1 - 18*lfh
      help(2) = 14 - 24*dlxm1 + 12*lfh
      help(3) = -11 + 12*dlxm1 - 6*lfh
      help(4) = 21/5.d0
      help(5) = 21/10.d0
      help(6) = 51/35.d0
      help(7) = 159/140.d0
      help(8) = 197/210.d0
      help(9) = 0.8023809523809524d0
      help(10) = 0.7025974025974026d0
      help(11) = 0.6257575757575757d0
      help(12) = 0.5645687645687646d0
      help(13) = 0.5145854145854146d0
      help(14) = 0.4729270729270729d0
      help(15) = 0.43763736263736264d0
      help(16) = 0.4073367808661926d0
      
c..   we want  x^n * expand( 1/x^n *( (1-x^n) * dcoefs(1,dlxm1)/(1-x) + sgg1) )
c..   where n=ncat

c..   these are the coefficients of (1-x)^i in expand( (1-x^n)/(1-x)*dcoefs ):
      helpd=0.d0
      do i=0,ncat-1
         helpd(i) = -bininv(-ncat,i+1)*dcoefs1(nmt,1.d0,dlxm1)
      enddo
      
      modified = 0
c..   this is  G(x) = sigma/xt^ncat :
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * (help(k-i) + helpd(k-i))
         enddo
         helpx(k) = sum
      enddo

      sgg1exp = 0.d0
      do i=0,nsft
         sgg1exp = sgg1exp + xm1**i * ( helpx(i) )
      enddo

      sgg1exp = xt**ncat * sgg1exp
      
      if (nmt.gt.0) sgg1exp = sgg1exp + sgg1mtexp(nmt,nsft
     &     ,ncat,xt)

      return
      end

C-}}}
C-{{{ function sgg1mtexact(xt)

      real*8 function sgg1mtexact(nmt,xt)
C..
C..   sgg1mtexact(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      integer nmt
      real*8 xt,xm1,dlxm1,dlnx
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      dlxm1 = dlog(1-xt)
      xm1 = 1.d0-xt
      dlnx = dlog(xt)

c..   only scalar so far
      if ((nmt.gt.0).and.lpseudo) stop 'sgg1mtexact'
      if (nmt.gt.10) stop 'sgg1mtexact'

c..   1/mt expansion:
      sgg1mtexact = 0.d0
      if (nmt.ge.2) then
         sgg1mtexact = sgg1mtexact - (mh/mt)**2* 1/64.d0 * (+ 308/45.d0
     &        *cf*ca **2*tr - 284/15.d0*cf*xt*ca**2*tr + 448/15.d0*cf*xt
     &        *ca **2*dlxm1*tr - 224/15.d0*cf*xt*dlnx*ca**2*tr + 284
     &        /15.d0*cf*xt**2*ca**2*tr - 224/15.d0*cf*xt**2*ca**2 *dlxm1
     &        *tr + 112/15.d0*cf*xt**2*dlnx*ca**2*tr - 308 /45.d0*cf*xt
     &        **3*ca**2*tr + 224/15.d0*cf*xt**3*ca**2 *dlxm1*tr - 112
     &        /15.d0*cf*xt**3*dlnx*ca**2*tr + 112 /15.d0*cf*dlnx/xm1*ca
     &        **2*tr)

         sgg1mtexact = sgg1mtexact - (mh/mt)**2 * lfh/64.d0 * ( - 224/15
     &        *cf*xt *ca**2*tr +112/15*cf*xt**2*ca**2*tr- 112/15*cf*xt
     &        **3*ca **2*tr)
      endif
      if (nmt.ge.4) then
         sgg1mtexact = sgg1mtexact + (mh/mt)**4*(-19/2800.d0 + (1543
     &        *dlnx)/8400.d0- (1543*dlxm1)/4200.d0 + (1543 *lfh)/8400.d0
     &        + 37/11200.d0/xt**2 + 39/11200.d0/xt - (1543 *dlnx)
     &        /(16800.d0*xm1) -(211*xm1)/11200.d0 - (1543*dlnx*xm1)
     &        /5600.d0 +(1543*dlxm1*xm1) /2800.d0 - (1543*lfh*xm1)
     &        /5600.d0 + (2*xm1**2)/175.d0+ (1543*dlnx*xm1**2) /8400.d0
     &        - (1543*dlxm1*xm1**2)/4200.d0+ (1543*lfh*xm1 **2)/8400.d0
     &        - (16973 *xm1**3)/201600.d0 -(1543*dlnx *xm1**3)/16800.d0
     &        + (1543*dlxm1*xm1**3)/8400.d0 -(1543 *lfh*xm1**3)
     &        /16800.d0)
      endif
      if (nmt.ge.6) then
         sgg1mtexact = sgg1mtexact - (mh/mt)**6*(5729/2016000.d0 - (113
     &        *dlnx)/4200.d0 +(113*dlxm1)/2100.d0 - (113 *lfh)/4200.d0 -
     &        23/12000.d0/xt**3 + 473/134400.d0/xt**2 - 1/225.d0/xt
     &        +(113*dlnx)/(8400.d0*xm1) + (1921*xm1)/1008000.d0 + (113
     &        *dlnx*xm1)/2800.d0 - (113*dlxm1*xm1)/1400.d0 + (113 *lfh
     &        *xm1)/2800.d0 +(569*xm1**2) /2016000.d0 - (113*dlnx*xm1
     &        **2)/4200.d0 + (113*dlxm1*xm1**2)/2100.d0 - (113*lfh *xm1
     &        **2)/4200.d0 + (1243*xm1**3)/100800.d0 + (113*dlnx*xm1
     &        **3)/8400.d0 - (113*dlxm1 *xm1**3)/4200.d0 + (113*lfh*xm1
     &        **3)/8400.d0)
      endif 
      if (nmt.ge.8) then
         sgg1mtexact = sgg1mtexact + (mh/mt)**8*(-2602757/2328480000.d0
     &        + (27677*dlnx)/6468000.d0 - (27677*dlxm1) /3234000.d0 +
     &        (27677*lfh) /6468000.d0 + 3347/(77616000.d0*xt**4) +
     &        104143/(141120000.d0 *xt**3) - 1219/(689920.d0*xt**2) +
     &        391849/(186278400.d0*xt) - (27677*dlnx)/(12936000.d0 *xm1)
     &        -(14159 *xm1)/47520000.d0 - (27677*dlnx*xm1) /4312000.d0
     &        +(27677 *dlxm1*xm1)/2156000.d0 - (27677*lfh *xm1)
     &        /4312000.d0- (1756999 *xm1**2)/4656960000.d0 + (27677*dlnx
     &        *xm1**2)/6468000.d0 - (27677 *dlxm1*xm1**2) /3234000.d0
     &        +(27677*lfh *xm1**2)/6468000.d0 - (27677 *xm1**3)
     &        /14112000.d0 -(27677 *dlnx*xm1**3)/12936000.d0 + (27677
     &        *dlxm1*xm1**3)/6468000.d0 - (27677*lfh*xm1**3)
     &        /12936000.d0)
      endif
      if (nmt.ge.10) then
         sgg1mtexact = sgg1mtexact - (mh/mt)**10*(284553923
     &        /847566720000.d0 -(182653*dlnx)/252252000.d0 + (182653
     &        *dlxm1) /126126000.d0- (182653*lfh)/252252000.d0 - 2597663
     &        /(94174080000.d0*xt**5)+2696483/(40360320000.d0*xt**4) -
     &        11868379/(40360320000.d0*xt**3) + 508611 /(896896000.d0*xt
     &        **2)-2241643/(3459456000.d0 *xt) + (182653*dlnx)
     &        /(504504000.d0*xm1)+(12341353*xm1) /169513344000.d0 +
     &        (182653*dlnx*xm1)/168168000.d0- (182653 *dlxm1*xm1)
     &        /84084000.d0 + (182653*lfh*xm1) /168168000.d0+ (100690649
     &        *xm1**2)/847566720000.d0 -(182653*dlnx*xm1**2)
     &        /252252000.d0+(182653*dlxm1*xm1 **2)/126126000.d0 -(182653
     &        *lfh*xm1**2)/252252000.d0 + (182653*xm1**3)/550368000.d0
     &        +(182653*dlnx*xm1**3) /504504000.d0 -(182653*dlxm1*xm1**3)
     &        /252252000.d0+(182653*lfh*xm1**3)/504504000.d0)
      endif
      
c$$$      factor = 0.d0
c$$$      do k=0,ncat1-1
c$$$         factor = factor + xt**k
c$$$      enddo
c..   factor = (1+xt+xt^2+...+xt^(n-1)) = (1-xt^n)/(1-xt)
c$$$      sgg1mtexact = sgg1mtexact + factor*dcoefs1mt(1.d0,dlxm1)

c--   checked against checks1l.m [18/06/09,rh]

c..   no lfr-terms

      end

C-}}}
C-{{{ function sgg1mtexp(xt)

      real*8 function sgg1mtexp(nmt,nsft,ncat,xt)
C..
C..   sgg1mt(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      real*8 xt,xm1,dlxm1,sum
      integer i,j,k,nmt,nsft,ncat
      real*8 helpmt(0:20),helpx(0:20),helpd(0:20),help(10,0:20)
      real*8 bininv
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)

      do i=1,10
         do j=0,20
            help(i,j)  = 0.d0
         enddo
      enddo
      
c..   only scalar so far, and muF=muR
      if ((nmt.gt.0).and.lpseudo) stop 'sgg1mtexp'
      if (nmt.gt.10) stop 'sgg1mtexp'
      if (nsft.gt.16) stop 'sgg1mtexp'
c..
c..   input produced with ~/math/gghmt/rgfull.m,
c..   function printitnlo[a,b]
c..
      help(1,0) = 7/10.d0 - (14*dlxm1)/5.d0 + (7*lfh)/5.d0
      help(1,1) = -6/5.d0 + (21*dlxm1)/5.d0 - (21*lfh)/10.d0
      help(1,2) = 107/60.d0 - (14*dlxm1)/5.d0 + (7*lfh)/5.d0
      help(1,3) = -77/60.d0 + (7*dlxm1)/5.d0 - (7*lfh)/10.d0
      help(1,4) = 49/100.d0
      help(1,5) = 49/200.d0
      help(1,6) = 17/100.d0
      help(1,7) = 53/400.d0
      help(1,8) = 197/1800.d0
      help(1,9) = 0.09361111111111112d0
      help(1,10) = 0.08196969696969697d0
      help(1,11) = 0.0730050505050505d0
      help(1,12) = 0.06586635586635586d0
      help(1,13) = 0.060034965034965034d0
      help(1,14) = 0.05517482517482517d0
      help(1,15) = 0.05105769230769231d0
      help(1,16) = 0.04752262443438914d0
      
      help(2,0) = 1543/16800.d0 - (1543*dlxm1)/4200.d0 + (1543*lfh)
     &     /8400.d0
      help(2,1) = -1641/11200.d0 + (1543*dlxm1)/2800.d0 - (1543*lfh)
     &     /5600.d0
      help(2,2) = 3013/12600.d0 - (1543*dlxm1)/4200.d0 + (1543*lfh)
     &     /8400.d0
      help(2,3) = -1529/10080.d0 + (1543*dlxm1)/8400.d0 - (1543*lfh)
     &     /16800.d0
      help(2,4) = 2023/24000.d0
      help(2,5) = 18631/336000.d0
      help(2,6) = 57521/1176000.d0
      help(2,7) = 222479/4704000.d0
      help(2,8) = 1007051/21168000.d0
      help(2,9) = 0.04880033541194256d0
      help(2,10) = 0.05057646619253762d0
      help(2,11) = 0.05270380892255892d0
      help(2,12) = 0.055070730195730194d0
      help(2,13) = 0.057609179510965225d0
      help(2,14) = 0.06027506422149279d0
      help(2,15) = 0.06303843700941915d0
      help(2,16) = 0.06587818108012436d0

      help(3,0) = 113/8400.d0 - (113*dlxm1)/2100.d0 + (113*lfh)/4200.d0
      help(3,1) = -53/2800.d0 + (113*dlxm1)/1400.d0 - (113*lfh)/2800.d0
      help(3,2) = 1051/28800.d0 - (113*dlxm1)/2100.d0 + (113*lfh)
     &     /4200.d0
      help(3,3) = -61/4032.d0 + (113*dlxm1)/4200.d0 - (113*lfh)/8400.d0
      help(3,4) = 50429/2016000.d0
      help(3,5) = 28513/1008000.d0
      help(3,6) = 518513/14112000.d0
      help(3,7) = 337531/7056000.d0
      help(3,8) = 2587729/42336000.d0
      help(3,9) = 0.07646664777021919d0
      help(3,10) = 0.0938069148113791d0
      help(3,11) = 0.1131152898027898d0
      help(3,12) = 0.13437542203167202d0
      help(3,13) = 0.15757734428270143d0
      help(3,14) = 0.1827145983976341d0
      help(3,15) = 0.20978279805511948d0
      help(3,16) = 0.23877885032218016d0

      help(4,0) = 27677/12936000.d0 - (27677*dlxm1)/3234000.d0 + (27677
     &     *lfh)/6468000.d0
      help(4,1) = -1/392.d0 + (27677*dlxm1)/2156000.d0 - (27677*lfh)
     &     /4312000.d0
      help(4,2) = 974389/155232000.d0 - (27677*dlxm1)/3234000.d0 +
     &     (27677*lfh)/6468000.d0
      help(4,3) = -600001/931392000.d0 + (27677*dlxm1)/6468000.d0 -
     &     (27677*lfh)/12936000.d0
      help(4,4) = 4276133/582120000.d0
      help(4,5) = 3380789/332640000.d0
      help(4,6) = 474006733/32598720000.d0
      help(4,7) = 163933873/8149680000.d0
      help(4,8) = 72967679/2716560000.d0
      help(4,9) = 0.03479654535515505d0
      help(4,10) = 0.04395790307382956d0
      help(4,11) = 0.05438288766608572d0
      help(4,12) = 0.06611202118588483d0
      help(4,13) = 0.07918684095014104d0
      help(4,14) = 0.09364944237353212d0
      help(4,15) = 0.10954225038590472d0
      help(4,16) = 0.1269078956239783d0

      help(5,0) = 182653/504504000.d0 - (182653*dlxm1)/126126000.d0 +
     &     (182653*lfh)/252252000.d0
      help(5,1) = -3559/10192000.d0 + (182653*dlxm1)/84084000.d0 -
     &     (182653*lfh)/168168000.d0
      help(5,2) = 14320547/12108096000.d0 -(182653*dlxm1)/126126000.d0 +
     &     (182653*lfh)/252252000.d0
      help(5,3) = 1153249/4036032000.d0 + (182653*dlxm1)/252252000.d0 -
     &     (182653*lfh)/504504000.d0
      help(5,4) = 125283397/60540480000.d0
      help(5,5) = 18061163/5503680000.d0
      help(5,6) = 274428569/52972920000.d0
      help(5,7) = 6654603503.d0/847566720000.d0
      help(5,8) = 29149008179.d0/2542700160000.d0
      help(5,9) = 0.01622297335561579d0
      help(5,10) = 0.02236398142148084d0
      help(5,11) = 0.03015023468544843d0
      help(5,12) = 0.03987307004455689d0
      help(5,13) = 0.05185157982550857d0
      help(5,14) = 0.06643253442469092d0
      help(5,15) = 0.08399034362818704d0
      help(5,16) = 0.10492703566134588d0
      
      do k=0,nsft
         sum = 0.d0
         do j=2,nmt,2
            sum = sum + (mh/mt)**j * help(j/2,k)
         enddo
         helpmt(k) = sum
      enddo
      
c..   these are the coefficients of (1-x)^i in 
c..   expand( (1-x^n)/(1-x)*dcoefs ):
c.. (commented out because this is now taken care of in sgg1exp)
c$$$      helpd=0.d0
c$$$      do i=0,ncat-1
c$$$         helpd(i) = -bininv(-ncat,i+1)*dcoefs1mt(1.d0,dlxm1)
c$$$      enddo

      modified = 0
c..   this is  G(x) = sigma/xt^ncat :
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * helpmt(k-i)
         enddo
         helpx(k) = sum
      enddo
      
      sgg1mtexp = 0.d0
      do k=0,nsft
         sgg1mtexp = sgg1mtexp + xm1**k * helpx(k)
      enddo

      sgg1mtexp = xt**ncat * sgg1mtexp
      
      return
      end

C-}}}
C-{{{ function sqg1(yy)

      real*8 function sqg1(xt)
C..
C..   qg contribution at NLO, exact.
C..   
      implicit none
      real*8 xt,sqg1exact,sqg1exp
      include '../commons/common-expand.f'

      if (nexpand(1).eq.0) then
         sqg1 = sqg1exact(nmtlim(1,2),xt)
      else
         sqg1 = sqg1exp(nmtlim(1,2),nsoft(1,2),ncat1,xt)
      endif

      return
      end

C-}}}
C-{{{ function sqg1exact(yy)

      real*8 function sqg1exact(nmt,xt)
C..
C..   qg contribution at NLO, exact.
C..   
      implicit none
      integer nmt
      real*8 sqg1mtexact
      real*8 xt,xm1,dlxm1,dlnx
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)
      dlnx = dlog(xt)

c..   pseudo-scalar and scalar are the same:
      sqg1exact = 2/3.d0 - (2*dlnx)/3.d0 + (4*dlxm1)/3.d0 - (4
     &     *xm1)/3.d0 -xm1**2/3.d0 - (2*dlnx*xm1**2)/3.d0 + (4*dlxm1
     &     *xm1**2)/3.d0
      
      sqg1exact = sqg1exact +lfh*( -2/3.d0 - (2*xm1**2)/3.d0 )

c..   there are no lfr-terms
         
c--   checked against checks.m [17/06/09,rh]

      if (nmt.gt.0) sqg1exact = sqg1exact + sqg1mtexact(nmt,xt)

      return
      end

C-}}}
C-{{{ function sqg1exp(nsft,ncat,xt)

      real*8 function sqg1exp(nmt,nsft,ncat,xt)
C..
C..   sqg1exp(xt) is sqg1(xt) expanded in xm1
C..
      implicit none
      integer i,k,nmt,nsft,ncat
      real*8 help(0:20),helpx(0:20)
      real*8 xt,xm1,dlxm1,sum,bininv,sqg1mtexp
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)

      if (nsft.gt.16) stop 'sqg1exp(): nsft too large'

c..
c..   input produced with ~/math/gghmt/rgfull.m,
c..   function printitnlo[a,b]
c..
      help(0) = 2/3.d0 + (4*dlxm1)/3.d0 - (2*lfh)/3.d0
      help(1) = -2/3.d0
      help(2) = (4*dlxm1)/3.d0 - (2*lfh)/3.d0
      help(3) = 8/9.d0
      help(4) = 1/2.d0
      help(5) = 16/45.d0
      help(6) = 5/18.d0
      help(7) = 8/35.d0
      help(8) = 7/36.d0
      help(9) = 0.1693121693121693d0
      help(10) = 0.15d0
      help(11) = 0.13468013468013468d0
      help(12) = 0.12222222222222222d0
      help(13) = 0.11188811188811189d0
      help(14) = 0.10317460317460317d0
      help(15) = 0.09572649572649573d0
      help(16) = 0.08928571428571429d0
      
      modified = 0
c..   this is  G(x) = sigma/xt^ncat :
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * help(k-i)
         enddo
         helpx(k) = sum
      enddo

      sqg1exp = 0.d0
      do i=0,nsft
         sqg1exp = sqg1exp + xm1**i * ( helpx(i) )
      enddo

      sqg1exp = xt**ncat * sqg1exp

      if (nmt.gt.0) sqg1exp = sqg1exp + sqg1mtexp(nmt,nsft
     &     ,ncat,xt)

      return
      end

C-}}}
C-{{{ function sqg1mtexact(xt)

      real*8 function sqg1mtexact(nmt,xt)
C..
C..   sqg1mtexact(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      integer nmt
      real*8 xt,xm1,dlxm1,dlnx
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      dlxm1 = dlog(1-xt)
      xm1 = 1.d0-xt
      dlnx = dlog(xt)

c..   only scalar so far
      if ((nmt.gt.0).and.lpseudo) stop 'sgg1mt'
      if (nmt.gt.10) stop 'sqg1mtexact'

      sqg1mtexact = 0.d0 
      if (nmt.ge.2) then
         sqg1mtexact = sqg1mtexact - (mh/mt)**2* 1/64.d0 * (+ 176/45.d0
     &        *cf**2 /xt*ca*tr - 46/15.d0*cf**2*ca*tr - 112/15.d0*cf**2
     &        *ca *dlxm1*tr + 56/15.d0*cf**2*dlnx*ca*tr - 8/3.d0*cf**2
     &        *xt *ca*tr+112/15.d0*cf**2*xt*ca*dlxm1*tr - 56/15.d0*cf**2
     &        *xt*dlnx*ca*tr - 2/45.d0*cf**2*xt**2*ca*tr - 56/15.d0 *cf
     &        **2*xt**2*ca*dlxm1*tr + 28/15.d0*cf**2*xt**2*dlnx*ca *tr)

         sqg1mtexact = sqg1mtexact - (mh/mt)**2 * lfh/64.d0 * (+ 56/15
     &        *cf**2*ca*tr - 56/15*cf**2*xt*ca*tr + 28/15*cf**2*xt**2*ca
     &        *tr) 
      endif
      if (nmt.ge.4) then
         sqg1mtexact = sqg1mtexact + (mh/mt)**4* 1/64.d0 * (+ 3487
     &        /10800.d0*cf **2*1/xt**2*cA*tr - 1297/945.d0*cf**2*1/xt*cA
     &        *tr + 1593/1400.d0*cf**2*cA*tr + 1543/1575.d0*cf**2*cA
     &        *dlxm1*tr-1543/3150.d0*cf**2*dlnx*cA*tr + 703/9450.d0*cf
     &        **2*xt *cA*tr- 1543/1575.d0*cf**2*xt*cA*dlxm1*tr + 1543
     &        /3150.d0*cf**2*xt*dlnx*cA*tr + 6221/75600.d0*cf**2*xt **2
     &        *cA*tr +1543/3150.d0*cf**2*xt**2*cA*dlxm1*tr - 1543
     &        /6300.d0*cf**2*xt**2*dlnx*cA*tr - 1543/3150.d0*lfh*cf **2
     &        *cA*tr + 1543/3150.d0*lfh*cf**2*xt*cA*tr - 1543 /6300.d0
     &        *lfh*cf**2*xt**2*cA*tr)
      endif
      if (nmt.ge.6) then
         sqg1mtexact = sqg1mtexact - (mh/mt)**6* 1/64.d0 * (+ 539
     &        /13500.d0*cf **2*1/xt**3*cA*tr - 41/216.d0*cf**2*1/xt**2
     &        *cA*tr + 1838/4725.d0*cf**2*1/xt*cA*tr - 5233/18900.d0*cf
     &        **2*cA *tr-226/1575.d0*cf**2*cA*dlxm1*tr + 113/1575.d0*cf
     &        **2 *dlnx*cA*tr + 409/18900.d0*cf**2*xt*cA*tr + 226
     &        /1575.d0 *cf**2*xt*cA*dlxm1*tr - 113/1575.d0*cf**2*xt*dlnx
     &        *cA*tr - 533/27000.d0*cf**2*xt**2*cA*tr - 113/1575.d0*cf
     &        **2*xt **2*cA*dlxm1*tr + 113/3150.d0*cf**2*xt**2*dlnx*cA
     &        *tr + 113/1575.d0*lfh*cf**2*cA*tr - 113/1575.d0*lfh*cf**2
     &        *xt *cA*tr+ 113/3150.d0*lfh*cf**2*xt**2*cA*tr)
      endif
      if (nmt.ge.8) then
         sqg1mtexact = sqg1mtexact + (mh/mt)**8* 1/64.d0 * (+ 107369
     &        /18191250.d0*cf**2*1/xt**4*cA*tr - 20833/630000.d0*cf **2
     &        *1/xt**3*cA*tr+ 92251/1188000.d0*cf**2*1/xt**2*cA*tr -
     &        375883/3638250.d0*cf**2*1/xt*cA*tr + 302471 /4851000.d0*cf
     &        **2*cA*tr + 27677/1212750.d0*cf**2*cA *dlxm1*tr - 27677
     &        /2425500.d0*cf**2*dlnx*cA*tr - 1144387 /145530000.d0*cf**2
     &        *xt*cA*tr - 27677/1212750.d0*cf**2 *xt*cA*dlxm1*tr + 27677
     &        /2425500.d0*cf**2*xt*dlnx*cA*tr + 392407/97020000.d0*cf**2
     &        *xt**2*cA*tr +27677 /2425500.d0*cf**2*xt**2*cA*dlxm1*tr
     &        -27677/4851000.d0 *cf**2*xt**2*dlnx*cA*tr - 27677
     &        /2425500.d0*lfh*cf**2*cA *tr + 27677/2425500.d0*lfh*cf**2
     &        *xt*cA*tr - 27677 /4851000.d0*lfh*cf**2*xt**2*cA*tr)
      endif
      if (nmt.ge.10) then
         sqg1mtexact = sqg1mtexact - (mh/mt)**10* 1/64.d0 * (+ 872827
     &        /902947500.d0*cf**2*1/xt**5*cA*tr - 896803/141891750.d0
     &        *cf**2*1/xt**4*cA*tr + 9084427/515970000.d0*cf**2*1/xt **3
     &        *cA*tr -8857573/324324000.d0*cf**2*1/xt**2*cA*tr + 3790307
     &        /141891750.d0*cf**2*1/xt*cA*tr - 19593521 /1418917500.d0
     &        *cf**2*cA*tr - 182653/47297250.d0*cf**2 *cA*dlxm1*tr +
     &        182653/94594500.d0*cf**2*dlnx*cA*tr + 452219/227026800.d0
     &        *cf**2*xt*cA*tr + 182653/47297250.d0 *cf**2*xt*cA*dlxm1*tr
     &        -182653/94594500.d0*cf**2*xt*dlnx *cA*tr - 63821123
     &        /79459380000.d0*cf**2*xt**2*cA*tr - 182653/94594500.d0*cf
     &        **2*xt**2*cA*dlxm1*tr + 182653 /189189000.d0*cf**2*xt**2
     &        *dlnx*cA*tr + 182653 /94594500.d0*lfh*cf**2*cA*tr - 182653
     &        /94594500.d0*lfh *cf**2*xt*cA*tr + 182653/189189000.d0*lfh
     &        *cf**2*xt**2 *cA*tr)
      endif

c--   checked against checks1l.m [18/06/09,rh]

c..   no lfr-terms

      return
      end

C-}}}
C-{{{ function sqg1mtexp(xt)

      real*8 function sqg1mtexp(nmt,nsft,ncat,xt)
C..
C..   sqg1mt(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      integer i,j,k,nmt,nsft,ncat
      real*8 xt,xm1,dlxm1,sum,bininv
      real*8 help(10,0:20),helpmt(0:20),helpx(0:20)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'
      include '../commons/common-expand.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)

      do i=1,10
         do j=0,20
            help(i,j) = 0.d0
         enddo
      enddo

c..   only scalar so far, and muF=muR
      if ((nmt.gt.0).and.lpseudo) stop 'sqg1mtexp'
      if (nmt.gt.10) stop 'sqg1mtexp'
      if (nsft.gt.16) stop 'sqg1mtexp'

c..
c..   input produced with ~/math/gghmt/rgfull.m,
c..   function printitnlo[a,b]
c..
      help(1,0) = 7/90.d0 + (7*dlxm1)/45.d0 - (7*lfh)/90.d0
      help(1,1) = -1/5.d0
      help(1,2) = -11/90.d0 + (7*dlxm1)/45.d0 - (7*lfh)/90.d0
      help(1,3) = -8/135.d0
      help(1,4) = -113/1080.d0
      help(1,5) = -82/675.d0
      help(1,6) = -47/360.d0
      help(1,7) = -92/675.d0
      help(1,8) = -101/720.d0
      help(1,9) = -0.14320987654320988d0
      help(1,10) = -0.14546296296296296d0
      help(1,11) = -0.14725028058361392d0
      help(1,12) = -0.1487037037037037d0
      help(1,13) = -0.1499093499093499d0
      help(1,14) = -0.15092592592592594d0
      help(1,15) = -0.1517948717948718d0
      help(1,16) = -0.1525462962962963d0

      help(2,0) = 1543/151200.d0 + (1543*dlxm1)/75600.d0 - (1543*lfh)
     &     /151200.d0
      help(2,1) = -4541/151200.d0
      help(2,2) = -2509/302400.d0 + (1543*dlxm1)/75600.d0 - (1543*lfh)
     &     /151200.d0
      help(2,3) = 221/21600.d0
      help(2,4) = 383/21600.d0
      help(2,5) = 131423/4536000.d0
      help(2,6) = 37409/907200.d0
      help(2,7) = 15857/294000.d0
      help(2,8) = 80881/1209600.d0
      help(2,9) = 0.07993412488452172d0
      help(2,10) = 0.09309143518518519d0
      help(2,11) = 0.10630985783763562d0
      help(2,12) = 0.11957208994708994d0
      help(2,13) = 0.13286683224183224d0
      help(2,14) = 0.14618638196019149d0
      help(2,15) = 0.1595253018586352d0
      help(2,16) = 0.1728796414399093d0

      help(3,0) = 113/75600.d0 + (113*dlxm1)/37800.d0 - (113*lfh)
     &     /75600.d0
      help(3,1) = -1/216.d0
      help(3,2) = -1/1120.d0 + (113*dlxm1)/37800.d0 - (113*lfh)/75600.d0
      help(3,3) = 89/113400.d0
      help(3,4) = -1/2016.d0
      help(3,5) = -41/14175.d0
      help(3,6) = -5143/756000.d0
      help(3,7) = -181/14700.d0
      help(3,8) = -2353/120960.d0
      help(3,9) = -0.02823591164861006d0
      help(3,10) = -0.038669642857142854d0
      help(3,11) = -0.05075800331355887d0
      help(3,12) = -0.06450352733686067d0
      help(3,13) = -0.07990786990786991d0
      help(3,14) = -0.09697215923406399d0
      help(3,15) = -0.1156971916971917d0
      help(3,16) = -0.13608354591836735d0

      help(4,0) = 27677/116424000.d0 + (27677*dlxm1)/58212000.d0 -
     &     (27677*lfh)/116424000.d0
      help(4,1) = -4187/5544000.d0
      help(4,2) = -1973/16632000.d0 + (27677*dlxm1)/58212000.d0 - (27677
     &     *lfh)/116424000.d0
      help(4,3) = 107/1134000.d0
      help(4,4) = -17/1862784.d0
      help(4,5) = 42071/582120000.d0
      help(4,6) = 1819999/3492720000.d0
      help(4,7) = 27256/17364375.d0
      help(4,8) = 3221767/931392000.d0
      help(4,9) = 0.006433153067223362d0
      help(4,10) = 0.010736639710598043d0
      help(4,11) = 0.0166148226579684d0
      help(4,12) = 0.024313225222749032d0
      help(4,13) = 0.03407751057264044d0
      help(4,14) = 0.046153425655976674d0
      help(4,15) = 0.060786770196294006d0
      help(4,16) = 0.07822337855085516d0

      help(5,0) = 182653/4540536000.d0 + (182653*dlxm1)/2270268000.d0 -
     &     (182653*lfh)/4540536000.d0
      help(5,1) = -294313/2270268000.d0
      help(5,2) = -16249/908107200.d0 + (182653*dlxm1)/2270268000.d0 -
     &     (182653*lfh)/4540536000.d0
      help(5,3) = 893/65488500.d0
      help(5,4) = 629/1009008000.d0
      help(5,5) = 25777/8513505000.d0
      help(5,6) = -1153501/136216080000.d0
      help(5,7) = -12101/127338750.d0
      help(5,8) = -2325041/6519744000.d0
      help(5,9) = -0.0009330118385296012d0
      help(5,10) = -0.0020041080006853815d0
      help(5,11) = -0.0037899879964334027d0
      help(5,12) = -0.006550970619621413d0
      help(5,13) = -0.010587627749465912d0
      help(5,14) = -0.016240793839149847d0
      help(5,15) = -0.023891571162459005d0
      help(5,16) = -0.03396133288362242d0

c..   if lzinv=.true., reconstruct the 1/xt-terms:
C$$$      lzinv = .false.
C$$$      if (lzinv) then
C$$$         do i=0,nsft1qgmt
C$$$            help(1,i) = help(1,i) + 22/135.d0
C$$$         enddo
C$$$         sqg1mtexp = -22/135.d0/xt*(mh/mt)**2
C$$$      else
C$$$         sqg1mtexp = 0.d0
C$$$      endif

      do k=0,nsft
         sum = 0.d0
         do j=2,nmt,2
            sum = sum + (mh/mt)**j * help(j/2,k)
         enddo
         helpmt(k) = sum
      enddo

      modified = 0
c..   this is  G(x) = sigma/xt^ncat:
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * helpmt(k-i)
         enddo
         helpx(k) = sum
      enddo

      sqg1mtexp = 0.d0
      do k=0,nsft
         sqg1mtexp = sqg1mtexp + xm1**k * helpx(k)
      enddo

      sqg1mtexp = xt**ncat * sqg1mtexp

c--   checked against checks1l.m [19/06/09,rh]

      return
      end

C-}}}
C-{{{ function sqqb1(yy)

      real*8 function sqqb1(xt)
C..
C..   q-qbar contribution at NLO, exact.
C..   
      implicit none
      real*8 xt,sqqb1exact,sqqb1exp
      include '../commons/common-expand.f'
      
      if (nexpand(1).eq.0) then
         sqqb1 = sqqb1exact(nmtlim(1,3),xt)
      else
         sqqb1 = sqqb1exp(nmtlim(1,3),nsoft(1,3),ncat1,xt)
      endif

      return
      end

C-}}}
C-{{{ function sqqb1exact(yy)

      real*8 function sqqb1exact(nmt,xt)
C..
C..   q-qbar contribution at NLO, exact.
C..   
      implicit none
      integer nmt
      real*8 xt,xm1,sqqb1mtexact
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      xm1 = 1.d0 - xt

c..   pseudo-scalar and scalar are the same:
      sqqb1exact = 32/27.d0*xm1**3

c..   there are no lfh- and no lfr-terms

c--   checked against checks.m [17/06/09,rh]

      if (nmt.gt.0) sqqb1exact = sqqb1exact + sqqb1mtexact(nmt,xt)

      return
      end

C-}}}
C-{{{ function sqqb1exp(yy)

      real*8 function sqqb1exp(nmt,nsft,ncat,xt)
C..
C..   q-qbar contribution at NLO, exact.
C..   
      implicit none
      integer i,k,nmt,nsft,ncat
      real*8 xt,xm1,dlxm1,sqqb1mtexp,sum,bininv
      real*8 help(0:20),helpx(0:20)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0 - xt
      dlxm1 = dlog(xm1)

      do i=0,20
         help(i) = 0.d0
      enddo

      help(3) = 32/27.d0

      modified = 0
c..   this is  G(x) = sigma/xt^ncat :
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * help(k-i)
         enddo
         helpx(k) = sum
      enddo

      sqqb1exp = 0.d0
      do i=0,nsft
         sqqb1exp = sqqb1exp + xm1**i * ( helpx(i) )
      enddo

      sqqb1exp = xt**ncat * sqqb1exp

      if (nmt.gt.1) sqqb1exp = sqqb1exp + sqqb1mtexp(nmt
     &     ,nsft,ncat,xt)

      return
      end

C-}}}
C-{{{ function sqqb1mtexact(xt)

      real*8 function sqqb1mtexact(nmt,xt)
C..
C..   sqqb1mtexact(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      integer nmt
      real*8 xt,xm1,dlxm1,dlnx
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      dlxm1 = dlog(1-xt)
      xm1 = 1.d0-xt
      dlnx = dlog(xt)

c..   only scalar so far
      if ((nmt.gt.0).and.lpseudo) stop 'sqqb1mtexact'
      if (nmt.gt.10) stop 'sqqb1mtexact'

      sqqb1mtexact = 0.d0 
      if (nmt.ge.2) then
         sqqb1mtexact = sqqb1mtexact - (mh/mt)**2* 1/64.d0 * (- 88/45.d0
     &        *cf**3*1/xt*cA + 208/45.d0*cf**3*cA - 32/15.d0*cf**3*xt*cA
     &        -16/9.d0*cf**3*xt**2*cA + 56/45.d0*cf**3*xt**3*cA)
      endif
      if (nmt.ge.4) then
         sqqb1mtexact = sqqb1mtexact + (mh/mt)**4* 1/64.d0 * (+ 3487
     &        /9450.d0*cf**3*1/xt**2*cA - 7463/9450.d0*cf**3*1/xt*cA +
     &        43/135.d0*cf**3*cA + 439/4725.d0*cf**3*xt*cA + 233/1350.d0
     &        *cf**3*xt**2*cA- 1543/9450.d0*cf**3*xt**3*cA)
      endif
      if (nmt.ge.6) then
         sqqb1mtexact = sqqb1mtexact - (mh/mt)**6* 1/64.d0 * (- 49
     &        /675.d0*cf**3 *1/xt**3*cA + 46/315.d0*cf**3*1/xt**2*cA -
     &        83/1575.d0 *cf**3*1/xt*cA - 76/4725.d0*cf**3*cA - 11
     &        /1575.d0*cf**3 *xt*cA-34/1575.d0*cf**3*xt**2*cA + 113
     &        /4725.d0*cf**3 *xt**3*cA)
      endif
      if (nmt.ge.8) then
         sqqb1mtexact = sqqb1mtexact + (mh/mt)**8* 1/64.d0 * (+ 107369
     &        /7276500.d0*cf**3*1/xt**4*cA - 22969/808500.d0*cf**3*1 /xt
     &        **3*cA+6257/661500.d0*cf**3*1/xt**2*cA + 53 /18900.d0*cf
     &        **3*1/xt*cA + 37/26460.d0*cf**3*cA + 4841 /7276500.d0*cf
     &        **3*xt*cA+2071/661500.d0*cf**3*xt**2*cA - 27677/7276500.d0
     &        *cf**3*xt**3*cA)
      endif
      if (nmt.ge.10) then
         sqqb1mtexact = sqqb1mtexact - (mh/mt)**10* 1/64.d0 * (- 872827
     &        /283783500.d0*cf**3*1/xt**5*cA + 812887/141891750.d0*cf
     &        **3*1/xt**4*cA - 9802/5457375.d0*cf**3*1/xt**3*cA - 19
     &        /36750.d0*cf**3*1/xt**2*cA - 17/66150.d0*cf**3*1/xt*cA -17
     &        /110250.d0*cf**3*cA - 466/6449625.d0*cf**3*xt*cA - 5461
     &        /10914750.d0*cf**3*xt**2*cA + 182653/283783500.d0 *cf**3
     &        *xt**3*cA)
      endif

c--   checked against checks1l.m [18/06/09,rh]

c..   no lfr-terms

      return
      end

C-}}}
C-{{{ function sqqb1mtexp(xt)

      real*8 function sqqb1mtexp(nmt,nsft,ncat,xt)
C..
C..   sqqb1mt(xt) are the 1/mt-terms of the hard one-loop result
C..
      implicit none
      integer i,j,k,nmt,nsft,ncat
      real*8 xt,xm1,dlxm1,sum,bininv
      real*8 help(10,0:20),helpmt(0:20),helpx(0:20)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'
      include '../commons/common-keys.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      dlxm1 = dlog(xm1)

      do i=1,10
         do j=0,20
            help(i,j) = 0.d0
         enddo
      enddo

c..   only scalar so far, and muF=muR
      if ((nmt.gt.0).and.lpseudo) stop 'sqqb1mtexp'
      if (nmt.gt.10) stop 'sqqb1mtexp'
      if (nsft.gt.16) stop 'sqqb1mtexp'

c..
c..   input produced with ~/math/gghmt/rgfull.m,
c..   function printitnlo[a,b]
c..
      help(1,0) = 0.d0
      help(1,1) = 0.d0
      help(1,2) = 0.d0
      help(1,3) = 16/45.d0
      help(1,4) = 88/405.d0
      help(1,5) = 88/405.d0
      help(1,6) = 88/405.d0
      help(1,7) = 88/405.d0
      help(1,8) = 88/405.d0
      help(1,9) = 88/405.d0
      help(1,10) = 88/405.d0
      help(1,11) = 88/405.d0
      help(1,12) = 88/405.d0
      help(1,13) = 88/405.d0
      help(1,14) = 88/405.d0
      help(1,15) = 88/405.d0
      help(1,16) = 88/405.d0

      help(2,0) = 0.d0
      help(2,1) = 0.d0
      help(2,2) = 0.d0
      help(2,3) = 446/4725.d0
      help(2,4) = 554/4725.d0
      help(2,5) = 13459/85050.d0
      help(2,6) = 8473/42525.d0
      help(2,7) = 973/4050.d0
      help(2,8) = 2392/8505.d0
      help(2,9) = 0.32224573780129334d0
      help(2,10) = 0.36324514991181656d0
      help(2,11) = 0.4042445620223398d0
      help(2,12) = 0.44524397413286304d0
      help(2,13) = 0.48624338624338626d0
      help(2,14) = 0.5272427983539094d0
      help(2,15) = 0.5682422104644327d0
      help(2,16) = 0.609241622574956d0

      help(3,0) = 0.d0
      help(3,1) = 0.d0
      help(3,2) = 0.d0
      help(3,3) = 344/14175.d0
      help(3,4) = 8/175.d0
      help(3,5) = 368/4725.d0
      help(3,6) = 5023/42525.d0
      help(3,7) = 337/2025.d0
      help(3,8) = 3158/14175.d0
      help(3,9) = 0.2872192827748383d0
      help(3,10) = 0.35971781305114636d0
      help(3,11) = 0.4402821869488536d0
      help(3,12) = 0.5289124044679601d0
      help(3,13) = 0.6256084656084656d0
      help(3,14) = 0.7303703703703703d0
      help(3,15) = 0.8431981187536743d0
      help(3,16) = 0.9640917107583774d0

      help(4,0) = 0.d0
      help(4,1) = 0.d0
      help(4,2) = 0.d0
      help(4,3) = 2242/363825.d0
      help(4,4) = 51082/3274425.d0
      help(4,5) = 70163/2182950.d0
      help(4,6) = 29627/519750.d0
      help(4,7) = 1202669/13097700.d0
      help(4,8) = 2263322/16372125.d0
      help(4,9) = 0.19790039472579155d0
      help(4,10) = 0.2724363819601915d0
      help(4,11) = 0.3634898799025783d0
      help(4,12) = 0.47270039777976286d0
      help(4,13) = 0.601707444818556d0
      help(4,14) = 0.7521505302457684d0
      help(4,15) = 0.9256691632882109d0
      help(4,16) = 1.1239028531726944d0

      help(5,0) = 0.d0
      help(5,1) = 0.d0
      help(5,2) = 0.d0
      help(5,3) = 22072/14189175.d0
      help(5,4) = 478/96525.d0
      help(5,5) = 1529494/127702575.d0
      help(5,6) = 31050161/1277025750.d0
      help(5,7) = 9378004/212837625.d0
      help(5,8) = 62708909/851350500.d0
      help(5,9) = 0.11588445064635541d0
      help(5,10) = 0.17386327370454355d0
      help(5,11) = 0.2510589849108368d0
      help(5,12) = 0.3512776848861505d0
      help(5,13) = 0.4786672163815021d0
      help(5,14) = 0.6377171642780108d0
      help(5,15) = 0.8332588555868979d0
      help(5,16) = 1.0704653594494864d0

      do k=0,nsft
         sum = 0.d0
         do j=2,nmt,2
            sum = sum + (mh/mt)**j * help(j/2,k)
         enddo
         helpmt(k) = sum
      enddo

      modified = 0
c..   this is  G(x) = sigma/xt^ncat :
      do k=0,nsft
         sum = 0.d0
         do i=0,k
            sum = sum + bininv(ncat,i) * helpmt(k-i)
         enddo
         helpx(k) = sum
      enddo

      sqqb1mtexp = 0.d0
      do k=0,nsft
         sqqb1mtexp = sqqb1mtexp + xm1**k * helpx(k)
      enddo

      sqqb1mtexp = xt**ncat * sqqb1mtexp

c--   checked against checks1l.m [19/06/09,rh]

      return
      end

C-}}}
C-{{{ function ppgg1(yy)

      real*8 function ppgg1(xx)
c..
c..   Integrand for hard gg-contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external ggpdf,sgg1

      ppgg1 = vartrans(ggpdf,sgg1,xx)

      end

C-}}}
C-{{{ function ppgg1x0(yy)

      real*8 function ppgg1x0(xx)
c..
c..   Integrand for hard gg-contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external ggpdf,sgg1x0

      ppgg1x0 = vartrans(ggpdf,sgg1x0,xx)

      end

C-}}}
C-{{{ function ppqg1(yy)

      real*8 function ppqg1(xx)
c..
c..   Integrand for hard qg-contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external qgpdf,sqg1

      ppqg1 = vartrans(qgpdf,sqg1,xx)

      end

C-}}}
C-{{{ function ppqg1x0(yy)

      real*8 function ppqg1x0(xx)
c..
c..   Integrand for hard qg-contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external qgpdf,sqg1x0

      ppqg1x0 = vartrans(qgpdf,sqg1x0,xx)

      return
      end

C-}}}
C-{{{ function ppqqb1(yy)

      real*8 function ppqqb1(xx)
c..
c..   Integrand for hard q-qbar contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external qqbpdf,sqqb1

      ppqqb1 = vartrans(qqbpdf,sqqb1,xx)

      return
      end

C-}}}
C-{{{ function ppqqb1x0(yy)

      real*8 function ppqqb1x0(xx)
c..
c..   Integrand for hard q-qbar contribution at NLO.
c..
      implicit none
      real*8 xx(10),vartrans
      external qqbpdf,sqqb1x0

      ppqqb1x0 = vartrans(qqbpdf,sqqb1x0,xx)

      return
      end

C-}}}
C-{{{ function ppall1(yy)

      real*8 function ppall1(xx)
c..
c..   Sum of all exact integrands of the sub-processes at NLO.
c..
      implicit none
      real*8 xx(10),ppgg1,ppqg1,ppqqb1

      ppall1 = ppgg1(xx) + ppqg1(xx) + ppqqb1(xx)

      return
      end

C-}}}
C-{{{ function ppall1x0(yy)

      real*8 function ppall1x0(xx)
c..
c..   Sum of all exact integrands of the sub-processes at NLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external ppgg1x0,ppqg1x0,ppqqb1x0

      ppall1x0 = ppgg1x0(xx) + ppqg1x0(xx) + ppqqb1x0(xx)

      return
      end

C-}}}

C-{{{ function ppgg2(yy)

      real*8 function ppgg2(xx)
c..
c..   Integrand for hard gg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external ggpdf,sgg2

      ppgg2 = vartrans(ggpdf,sgg2,xx)
      
      return
      end

C-}}}
C-{{{ function ppgg2x0(yy)

      real*8 function ppgg2x0(xx)
c..
c..   Integrand for hard gg-contribution at NNLO.
c..
      implicit none
      real*8 xx(10),sgg2x0,ggpdf,vartrans
      external ggpdf,sgg2x0

      ppgg2x0 = vartrans(ggpdf,sgg2x0,xx)

      return
      end

C-}}}
C-{{{ function ppqg2(yy)

      real*8 function ppqg2(xx)
c..
c..   Integrand for hard qg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qgpdf,sqg2

      ppqg2 = vartrans(qgpdf,sqg2,xx)

      end

C-}}}
C-{{{ function ppqg2x0(yy)

      real*8 function ppqg2x0(xx)
c..
c..   Integrand for hard qg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qgpdf,sqg2x0

      ppqg2x0 = vartrans(qgpdf,sqg2x0,xx)

      return
      end

C-}}}
C-{{{ function ppqqb2(yy)

      real*8 function ppqqb2(xx)
c..
c..   Integrand for hard qqb-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqbpdf,sqqb2

      ppqqb2 = vartrans(qqbpdf,sqqb2,xx)

      end

C-}}}
C-{{{ function ppqqb2x0(yy)

      real*8 function ppqqb2x0(xx)
c..
c..   Integrand for hard qqb-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqbpdf,sqqb2x0

      ppqqb2x0 = vartrans(qqbpdf,sqqb2x0,xx)

      return
      end

C-}}}
C-{{{ function ppqq2(yy)

      real*8 function ppqq2(xx)
c..
c..   Integrand for hard qq-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqpdf,sqq2

      ppqq2 = vartrans(qqpdf,sqq2,xx)

      end

C-}}}
C-{{{ function ppqq2x0(yy)

      real*8 function ppqq2x0(xx)
c..
c..   Integrand for hard qq-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqpdf,sqq2x0

      ppqq2x0 = vartrans(qqpdf,sqq2x0,xx)

      return
      end

C-}}}
C-{{{ function ppqu2(yy)

      real*8 function ppqu2(xx)
c..
c..   Integrand for hard qu-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qupdf,qubpdf,squ2

      ppqu2 = vartrans(qupdf,squ2,xx) + vartrans(qubpdf,squ2,xx)

      end

C-}}}
C-{{{ function ppqu2x0(yy)

      real*8 function ppqu2x0(xx)
c..
c..   Integrand for hard qu-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qupdf,qubpdf,squ2x0

      ppqu2x0 = vartrans(qupdf,squ2x0,xx) + vartrans(qubpdf,squ2x0,xx)

      return
      end

C-}}}
C-{{{ function ppall2(yy)

      real*8 function ppall2(xx)
c..
c..   Sum of all integrands of the sub-processes at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external ppgg2,ppqg2,ppqqb2,ppqq2,ppqu2

      ppall2 = ppgg2(xx) + ppqg2(xx) + ppqqb2(xx)+ ppqq2(xx) +
     &     ppqu2(xx)

      end

C-}}}
C-{{{ function ppall2x0(yy)

      real*8 function ppall2x0(xx)
c..
c..   Sum of all integrands of the sub-processes at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external ppgg2x0,ppqg2x0,ppqqb2x0,ppqq2x0,ppqu2x0

      ppall2x0 = ppgg2x0(xx) + ppqg2x0(xx) + ppqqb2x0(xx)+ ppqq2x0(xx) +
     &     ppqu2x0(xx)

      return
      end

C-}}}

C-{{{ function ppgg3(yy)

      real*8 function ppgg3(xx)
c..
c..   Integrand for hard gg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external ggpdf,sgg3

      ppgg3 = vartrans(ggpdf,sgg3,xx)

      return
      end

C-}}}
C-{{{ function ppgg3x0(yy)

      real*8 function ppgg3x0(xx)
c..
c..   Integrand for hard gg-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external ggpdf,sgg3x0

      ppgg3x0 = vartrans(ggpdf,sgg3x0,xx)

      return
      end

C-}}}
C-{{{ function ppqg3(yy)

      real*8 function ppqg3(xx)
c..
c..   Integrand for hard qg-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qgpdf,sqg3

      ppqg3 = vartrans(qgpdf,sqg3,xx)

      end

C-}}}
C-{{{ function ppqg3x0(yy)

      real*8 function ppqg3x0(xx)
c..
c..   Integrand for hard qg-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qgpdf,sqg3x0

      ppqg3x0 = vartrans(qgpdf,sqg3x0,xx)

      return
      end

C-}}}
C-{{{ function ppqqb3(yy)

      real*8 function ppqqb3(xx)
c..
c..   Integrand for hard qqb-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqbpdf,sqqb3

      ppqqb3 = vartrans(qqbpdf,sqqb3,xx)

      end

C-}}}
C-{{{ function ppqqb3x0(yy)

      real*8 function ppqqb3x0(xx)
c..
c..   Integrand for hard qqb-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqbpdf,sqqb3x0

      ppqqb3x0 = vartrans(qqbpdf,sqqb3x0,xx)

      return
      end

C-}}}
C-{{{ function ppqq3(yy)

      real*8 function ppqq3(xx)
c..
c..   Integrand for hard qq-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqpdf,sqq3

      ppqq3 = vartrans(qqpdf,sqq3,xx)

      end

C-}}}
C-{{{ function ppqq3x0(yy)

      real*8 function ppqq3x0(xx)
c..
c..   Integrand for hard qq-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqpdf,sqq3x0

      ppqq3x0 = vartrans(qqpdf,sqq3x0,xx)

      return
      end

C-}}}
C-{{{ function ppqu3(yy)

      real*8 function ppqu3(xx)
c..
c..   Integrand for hard qu-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qupdf,qubpdf,squ3

      ppqu3 = vartrans(qupdf,squ3,xx) + vartrans(qubpdf,squ3,xx)

      end

C-}}}
C-{{{ function ppqu3x0(yy)

      real*8 function ppqu3x0(xx)
c..
c..   Integrand for hard qu-contribution at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qupdf,qubpdf,squ3x0

      ppqu3x0 = vartrans(qupdf,squ3x0,xx) + vartrans(qubpdf,squ3x0,xx)

      return
      end

C-}}}
C-{{{ function ppall3(yy)

      real*8 function ppall3(xx)
c..
c..   Sum of all integrands of the sub-processes at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external ppgg3,ppqg3,ppqqb3,ppqq3,ppqu3

      ppall3 = ppgg3(xx) + ppqg3(xx) + ppqqb3(xx)+ ppqq3(xx) +
     &     ppqu3(xx)

      end

C-}}}
C-{{{ function ppall3x0(yy)

      real*8 function ppall3x0(xx)
c..
c..   Sum of all integrands of the sub-processes at N3LO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external ppgg3x0,ppqg3x0,ppqqb3x0,ppqq3x0,ppqu3x0

      ppall3x0 = ppgg3x0(xx) + ppqg3x0(xx) + ppqqb3x0(xx)+ ppqq3x0(xx) +
     &     ppqu3x0(xx)

      return
      end

C-}}}

C-{{{ function fmatch

      real*8 function fmatch(xt,xt0,ww)

      implicit none
      real*8 xt,xt0,ww

      fmatch = 1/3.d0*( 1 + dtanh( ( xt0 - xt )/ww ) )

      return
      end

C-}}}
C-{{{ function rfact:

      real*8 function rfact(nn)
c..
c..   rfact(n) = n!
c..
c..   note:  rfact(0) = 1
c..
      implicit none
      integer ii,nn

      if (nn.lt.0) then
         call printdieggh('Negativ argument in rfact().')
      endif

      rfact = 1.d0
      do ii=1,nn
         rfact = rfact*ii
      enddo

      return
      end

C-}}}
C-{{{ function bininv:

      real*8 function bininv(nn,kk)
c..
c..   bininv(nn,kk) = Gamma(nn+kk)/Gamma(kk+1)/Gamma(nn)
c..   
c..   1/(1-x)^n = sum( bininv(n,k) * x^k, k=0,...,infinity )
c..   
      implicit none
      real*8 denom,rfact
      integer nn,mnn,kk

      if (nn.gt.0) then
         denom = 1.d0*rfact(kk)*rfact(nn-1)
         bininv = rfact(nn+kk-1)/denom
      elseif (nn.lt.0) then
         mnn = -nn
         if (kk.gt.mnn) then
            bininv = 0.d0
         else
            denom = 1.d0*rfact(kk)*rfact(mnn-kk)
            bininv = (-1)**kk * rfact(mnn)/denom
         endif
      elseif (nn.eq.0) then
         if (kk.eq.0) then
            bininv = 1.d0
         else
            bininv = 0.d0
         endif
      else
         bininv = 0.d0
         call printdieggh('parameters out of range')
      endif

      return
      end

C-}}}
