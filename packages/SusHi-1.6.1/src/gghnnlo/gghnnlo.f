c..
c..   exact NNLO results for gg -> A  and  gg -> H
c..   
C-{{{ RCS:

c..   
c..   $Id: gghnnloexact.f,v 3.3 2005/02/13 15:21:16 rharland Exp $
c..   $Log: gghnnloexact.f,v $
c..   Revision 3.3  2005/02/13 15:21:16  rharland
c..   bug fix.
c..
c..   Revision 3.1  2004/11/26 15:15:58  rharland
c..   G- and sigma-expansion re-implemented
c..
c..   Revision 3.0  2004/11/18 13:41:48  rharland
c..   soft expansion included
c..
c..   Revision 2.1  2004/10/29 14:30:14  rharland
c..   pseudo-scalar included
c..   error reporting improved
c..
c..   Revision 2.0  2004/09/22 19:40:27  rharland
c..   version used for hep-ph/0409010
c..
c..   Revision 1.1  2004/03/23 16:44:31  rharland
c..   Initial revision
c..
c..   Revision 1.2  2004/02/20 03:06:09  rharland
c..   removed unused commons.
c..
c..   Revision 1.1  2004/02/19 22:59:52  rharland
c..   Initial revision
c..

C-}}}
C-{{{ function sgg2(xt)

      real*8 function sgg2(xt)
C..
C..   sgg2(xt) is the two-loop result minus the purely soft terms
C..
C..
      implicit none
      include '../commons/common-expand.f'
      real*8 xt,sgg2exact,sgg2exp

      if (nexpand(2).eq.0) then
         sgg2 = sgg2exact(xt)
      else
         sgg2 = sgg2exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sgg2exact(yy)

      real*8 function sgg2exact(xt)
C
C
      implicit real*8 (a-h,o-z)
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      sgg2exact = sggaexact(xt) + nf*sggfexact(xt)

c..   factor = (1+xt+xt^2+...+xt^(n-1)) = (1-xt^n)/(1-xt)
      factor = 0.d0
      do k=0,ncat2-1
         factor = factor + xt**k
      enddo
      sgg2exact = sgg2exact + factor*dcoefs2(nmtlim(2,1),1.d0,dlxm1
     &     ,dlxm1**2,dlxm1**3)

      if (nmtlim(2,1).gt.0) then
         stop 'Error in sgg2exact.'
      endif

c..   factor out SM Wilson coefficient:      
      sgg2exact = sgg2exact - 2*c1sm1/c1sm0*sgg1exact(ncat2,nmtlim(2,1)
     &     ,xt)

      return
      end

C-}}}
C-{{{ function sqg2(yy)

      real*8 function sqg2(xt)
C
      implicit none
      real*8 xt,sqgaexact,sqgfexact,sqg2exp,sqg1exact
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
C
      if (nexpand(2).eq.0) then
         sqg2 = sqgaexact(xt) + nf*sqgfexact(xt)
         if (nmtlim(2,2).gt.0) stop 'Error in sqg2.'
c..   factor out SM Wilson coefficient:
         sqg2 = sqg2 - 2*c1sm1/c1sm0*sqg1exact(nmtlim(2,2),xt)
      else
         sqg2 = sqg2exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sqqb2(yy)

      real*8 function sqqb2(xt)
C
      implicit none
      real*8 xt,sqqbaexact,sqqbfexact,sqqb2exp,sqqb1exact
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      if (nexpand(2).eq.0) then
         sqqb2 = sqqbaexact(xt) + nf*sqqbfexact(xt)
         if (nmtlim(2,3).gt.0) stop 'Error in sqqb2.'
c..   factor out SM Wilson coefficient:
         sqqb2 = sqqb2 - 2*c1sm1/c1sm0*sqqb1exact(nmtlim(2,3),xt)
      else
         sqqb2 = sqqb2exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sqq2(yy)

      real*8 function sqq2(xt)
C
      implicit none
      real*8 xt,sqqaexact,sqqfexact,sqq2exp
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
C
      if (nexpand(2).eq.0) then
         sqq2 = sqqaexact(xt) + nf*sqqfexact(xt)
         if (nmtlim(2,4).gt.0) stop 'Error in sqq2.'
      else
         sqq2 = sqq2exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function squ2(yy)

      real*8 function squ2(xt)
C
C
      implicit none
      real*8 xt,squaexact,squfexact,squ2exp
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      if (nexpand(2).eq.0) then
         squ2 = squaexact(xt) + nf*squfexact(xt)
         if (nmtlim(2,5).gt.0) stop 'Error in squ2.'
      else
         squ2 = squ2exp(xt)
      endif

      return
      end

C-}}}

c-{{{ function sggaexact:

      real*8 function sggaexact(xx)
      implicit real*8 (a-h,o-z)
      external sudilog,trilog
      complex*16 sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      sggaexact = 67.33333333333333d0 + (1983*dlnx)/4.d0 + 57*dlnx**3 -
     &     (1109*dlxm1)/4.d0 -198*dlnx*dlxm1 + 36*dlnx**2*dlxm1 + 66
     &     *dlxm1**2 - 9*dlnx*dlxm1**2 -144*dlxm1**3 + (21*dlnx**3)/(2
     &     .d0*(2 - xm1)) -(9*dlnx**2*dlxm1)/(2 - xm1) - (139*dlnx)/(2
     &     .d0*xm1) -(33*dlnx**2)/(8.d0*xm1) + (3*dlnx**3)/xm1 + (33
     &     *dlnx*dlxm1)/xm1 +(72*dlnx**2*dlxm1)/xm1 - (279*dlnx*dlxm1**2
     &     )/(2.d0*xm1) -(3435*xm1)/4.d0 - (5823*dlnx*xm1)/4.d0 - (657
     &     *dlnx**2*xm1)/4.d0 -(327*dlnx**3*xm1)/4.d0 + 1530*dlxm1*xm1 +
     &     1017*dlnx*dlxm1*xm1 -45*dlnx**2*dlxm1*xm1 - 675*dlxm1**2*xm1
     &     +(27*dlnx*dlxm1**2*xm1)/2.d0 + 216*dlxm1**3*xm1 + (11107*xm1
     &     **2)/12.d0 +(10365*dlnx*xm1**2)/8.d0 + (2007*dlnx**2*xm1**2)
     &     /8.d0 +24*dlnx**3*xm1**2 - (5567*dlxm1*xm1**2)/4.d0 -1149
     &     *dlnx*dlxm1*xm1**2 - 72*dlnx**2*dlxm1*xm1**2 

      sggaexact = sggaexact + 642*dlxm1**2*xm1**2 + 135*dlnx*dlxm1**2
     &     *xm1**2 - 144*dlxm1**3*xm1**2 - (7583*xm1**3)/16.d0 - (2171
     &     *dlnx*xm1**3)/4.d0 - (561*dlnx**2*xm1**3)/4.d0 + 3*dlnx**3
     &     *xm1**3+ (2583*dlxm1*xm1**3)/4.d0 + 561*dlnx*dlxm1*xm1**3 +
     &     81*dlnx**2*dlxm1*xm1**3 - 330*dlxm1**2*xm1**3 - (279*dlnx
     &     *dlxm1**2*xm1**3)/2.d0 + 72*dlxm1**3*xm1**3 - 66*z2 + 126
     &     *dlnx*z2 +180*dlxm1*z2 + (81*dlnx*z2)/xm1 + 675*xm1*z2 - 189
     &     *dlnx*xm1*z2 - 270*dlxm1*xm1*z2 - 642*xm1**2*z2 - 18*dlnx*xm1
     &     **2*z2 +180*dlxm1*xm1**2*z2 + (1089*xm1**3*z2)/4.d0 + 81*dlnx
     &     *xm1**3*z2 - 90*dlxm1*xm1**3*z2 - 351*z3 + (1053*xm1*z3)/2.d0
     &     - 351*xm1**2*z3 + (351*xm1**3*z3)/2.d0 + (-117 + 333*dlnx -
     &     432*dlxm1 - (27*dlnx)/(2 - xm1) + (36*dlxm1)/(2 - xm1) - 33
     &     /(4.d0*xm1) + (99*dlnx)/(2.d0*xm1) - (261*xm1)/4.d0 - (999
     &     *dlnx*xm1)/2.d0 + 612*dlxm1*xm1 - 138*xm1**2 + 144*dlnx*xm1
     &     **2 -144*dlxm1*xm1**2 + (363*xm1**3)/4.d0 + (99*dlnx*xm1**3)
     &     /2.d0- 36*dlxm1*xm1**3) *dli2a + (-4.5d0 + (45*dlnx)/2.d0 -
     &     72*dlxm1 + (18*dlnx) /(2 - xm1) - (18*dlxm1)/(2 - xm1) + (189
     &     *xm1)/4.d0 - (99*dlnx *xm1)/4.d0 + 126*dlxm1*xm1 - (81*xm1
     &     **2)/2.d0 - 9*dlnx*xm1**2 - 72*dlxm1*xm1**2 + (33*xm1**3)
     &     /4.d0 +(9*dlnx*xm1**3)/2.d0 + 18*dlxm1*xm1**3)*dli2b 

      sggaexact = sggaexact + (189 + 27/(2 - xm1) - 171 /(2.d0*xm1) -
     &     234*xm1 + 36*xm1**2 - (81*xm1**3)/2.d0)*dli3a + (99 + 45 /(2
     &     -xm1) - (333*xm1)/2.d0 + 72*xm1**2 - 18*xm1 **3)* dli3b + (
     &     -99- 45/(2 - xm1) + (333*xm1) /2.d0 - 72*xm1**2 + 18*xm1 **3)
     &     * dli3c + (-441 - 27/(2 - xm1) - 81/xm1 + (1449*xm1)/2.d0
     &     -270*xm1**2 - 27*xm1 **3)*dli3d + (29.25d0 - 27/(4.d0*(2 -
     &     xm1)) - (531 *xm1)/8.d0 + 63*xm1**2 - 18*xm1**3)* dli3e +
     &     (6.75d0 - 9/(4.d0*(2 - xm1)) - (189*xm1)/8.d0 + 27*xm1**2 - 9
     &     *xm1**3) *dli3f


c..   lfh-terms:
      sggaexact = sggaexact + lfh**2*(16.5d0 - 36*dlnx - 72*dlxm1 - (18
     &     *dlnx)/xm1 - (675*xm1)/4.d0 +54*dlnx*xm1 + 108*dlxm1*xm1 +
     &     (321*xm1**2)/2.d0 - 72*dlxm1*xm1**2 -(297*xm1**3)/4.d0 - 18
     &     *dlnx*xm1**3 + 36*dlxm1*xm1**3) +lfh*(139 + 102*dlnx - 45
     &     *dlnx**2 - 66*dlxm1 + 36*dlnx*dlxm1 +216*dlxm1**2 + (9*dlnx
     &     **2)/(2.d0*(2 - xm1)) - (33*dlnx)/(2.d0*xm1) -(45*dlnx**2)/(2
     &     .d0*xm1) + (126*dlnx*dlxm1)/xm1 - (3051*xm1)/4.d0 -513*dlnx
     &     *xm1 + 63*dlnx**2*xm1 + 675*dlxm1*xm1 -54*dlnx*dlxm1*xm1 -
     &     324*dlxm1**2*xm1 + (2773*xm1**2)/4.d0 +(1185*dlnx*xm1**2)/2
     &     .d0 + 9*dlnx**2*xm1**2 - 642*dlxm1*xm1**2 -108*dlnx*dlxm1*xm1
     &     **2 + 216*dlxm1**2*xm1**2 - (2449*xm1**3)/8.d0 -(561*dlnx*xm1
     &     **3)/2.d0 - 27*dlnx**2*xm1**3 + 330*dlxm1*xm1**3 +126*dlnx
     &     *dlxm1*xm1**3 - 108*dlxm1**2*xm1**3 - 90*z2 +135*xm1*z2 - 90
     &     *xm1**2*z2 + 45*xm1**3*z2 +(216 - 18/(2 - xm1) - 306*xm1 + 72
     &     *xm1**2 + 18*xm1**3)*dli2a +(36 + 9/(2 - xm1) - 63*xm1 +
     &     36*xm1**2 - 9*xm1**3)*dli2b)

c..   lfr-terms:
      sggaexact = sggaexact + lfr*(-99*dlnx + 198*dlxm1 + (99*dlnx)
     &     /(2.d0*xm1) + (297*dlnx*xm1)/2.d0 -297*dlxm1*xm1 - 99*dlnx
     &     *xm1**2 +198*dlxm1*xm1**2 +(363*xm1**3)/8.d0 + (99*dlnx*xm1
     &     **3)/2.d0 -99 *dlxm1*xm1**3 +lfh*(-99 + (297*xm1)/2.d0 - 99
     &     *xm1**2 + (99*xm1 **3)/2.d0))

c..   turn to scalar:
      if (.not.lpseudo) then
         sggaexact = sggaexact - ((69*dlnx)/2.d0 - 9*dlnx**2 - 12*dlxm1
     &        +6*lfh - (3*dlnx)/xm1 +27*xm1 - (39*dlnx*xm1)/2.d0 + 9
     &        *dlnx**2*xm1 + 18*dlxm1*xm1 -9*lfh*xm1 + (57*xm1**2)/4.d0
     &        + 6*dlnx*xm1**2 - 12*dlxm1*xm1 **2 +6*lfh*xm1**2 - (11*xm1
     &        **3)/4.d0 - 3*dlnx*xm1**3 + 6 *dlxm1*xm1**3 -3*lfh*xm1**3
     &        )

c..   no additional lfr-terms
         
      endif

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function sggfexact:

      real*8 function sggfexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      sggfexact = -3.111111111111111d0 - (265*dlnx)/216.d0 + (287*dlnx
     &     **2)/72.d0 +(17*dlnx**3)/72.d0 + (77*dlxm1)/12.d0 - (68*dlnx
     &     *dlxm1)/9.d0 -(16*dlnx**2*dlxm1)/3.d0 - 4*dlxm1**2 + (16*dlnx
     &     *dlxm1**2)/3.d0 +(5*dlnx)/(3.d0*xm1) + dlnx**2/(4.d0*xm1) -(2
     &     *dlnx*dlxm1)/xm1 +(8441*xm1)/216.d0 + (883*dlnx*xm1)/36.d0-
     &     (8*dlnx**2*xm1)/3.d0 -(dlnx**3*xm1)/3.d0 - (751*dlxm1*xm1)
     &     /18.d0 + (8*dlnx*dlxm1*xm1)/3.d0 +8*dlnx**2*dlxm1*xm1 + (38
     &     *dlxm1**2*xm1)/3.d0 - 8*dlnx*dlxm1**2*xm1 -(17227*xm1**2)
     &     /432.d0 - (2165*dlnx*xm1**2)/72.d0 -(59*dlnx**2*xm1**2)/24.d0
     &     +(dlnx**3*xm1**2)/8.d0 +(1487*dlxm1*xm1**2)/36.d0 + (26*dlnx
     &     *dlxm1*xm1**2)/3.d0 -(8*dlnx**2*dlxm1*xm1**2)/3.d0 - (32
     &     *dlxm1**2*xm1**2)/3.d0 +(8*dlnx*dlxm1**2*xm1**2)/3.d0 + (8887
     &     *xm1**3)/648.d0 

      sggfexact = sggfexact + (298*dlnx*xm1**3)/27.d0 + (11*dlnx**2*xm1
     &     **3)/6.d0 - (785*dlxm1*xm1**3)/54.d0 - (50*dlnx*dlxm1*xm1**3)
     &     /9.d0 + (34*dlxm1**2*xm1**3)/9.d0 + 4*z2 - (16*dlnx*z2)/3.d0
     &     -(38*xm1*z2)/3.d0 + 8*dlnx*xm1*z2 + (32*xm1**2*z2)/3.d0 - (8
     &     *dlnx*xm1**2*z2)/3.d0 - (34*xm1**3*z2)/9.d0 + (
     &     -12.13888888888889d0 - (95*dlnx)/12.d0 + (32*dlxm1)/3.d0 + 1
     &     /(2.d0*xm1) + (121*xm1)/6.d0 + 12*dlnx*xm1 - 16*dlxm1*xm1 -
     &     (27*xm1**2)/4.d0 - (47*dlnx*xm1**2)/12.d0 + (16*dlxm1*xm1**2)
     &     /3.d0 + xm1**3/9.d0)*dli2a + (-5.5d0 + 8*xm1 - (17*xm1**2)
     &     /6.d0)*dli3a + (5.25d0 - 8*xm1 + (31*xm1**2)/12.d0)*dli3d


c..   lfh-terms:
      sggfexact = sggfexact + lfh**2*(-1 + (4*dlnx)/3.d0 + (19*xm1)/6.d0
     &     -2*dlnx*xm1 -(8*xm1**2)/3.d0 + (2*dlnx*xm1**2)/3.d0 + (17*xm1
     &     **3)/18.d0) +lfh*(-3.3333333333333335d0 + (34*dlnx)/9.d0 + (8
     &     *dlnx**2)/3.d0 + 4*dlxm1 -(16*dlnx*dlxm1)/3.d0 + dlnx/xm1
     &     +(190*xm1)/9.d0 -(4*dlnx*xm1)/3.d0 -4*dlnx**2*xm1 - (38*dlxm1
     &     *xm1)/3.d0 + 8*dlnx*dlxm1*xm1 -(187*xm1**2)/9.d0 - (13*dlnx
     &     *xm1**2)/3.d0 + (4*dlnx**2*xm1**2)/3.d0 +(32*dlxm1*xm1**2)
     &     /3.d0 - (8*dlnx*dlxm1*xm1**2)/3.d0 +(785*xm1**3)/108.d0 + (25
     &     *dlnx*xm1**3)/9.d0 - (34*dlxm1*xm1**3)/9.d0 +(
     &     -5.333333333333333d0 + 8*xm1 - (8*xm1**2)/3.d0)*dli2a)

c..   lfr-terms:
      sggfexact = sggfexact + lfr*(6*dlnx - 12*dlxm1 - (3*dlnx)/xm1 - 9
     &     *dlnx*xm1 + 18*dlxm1*xm1+6*dlnx*xm1**2 - 12*dlxm1*xm1**2 -
     &     (11*xm1**3)/4.d0 - 3*dlnx*xm1**3 +6*dlxm1*xm1**3 + lfh*(6 - 9
     &     *xm1 + 6*xm1**2 - 3*xm1**3))

c..   turn to scalar:
      if (.not.lpseudo) then
         sggfexact = sggfexact - (dlnx + (2*dlnx**2)/3.d0 + (3*xm1)/2.d0
     &        -dlnx*xm1 - (2*dlnx**2 *xm1)/3.d0 -(5*xm1**2)/3.d0)
c..   no additional lfh- and lfr-terms
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function sqgaexact:

      real*8 function sqgaexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      sqgaexact = 0.7407407407407407d0 + (6103*dlnx)/108.d0 + (451*dlnx
     &     **2)/54.d0 +(50*dlnx**3)/3.d0 + (353*dlxm1)/18.d0 - (601*dlnx
     &     *dlxm1)/27.d0 +(392*dlnx**2*dlxm1)/9.d0 + (85*dlxm1**2)/36.d0
     &     -(1385*dlnx*dlxm1**2)/18.d0 + (367*dlxm1**3)/54.d0 - (11119
     &     *xm1)/108.d0 -(8521*dlnx*xm1)/54.d0 - (251*dlnx**2*xm1)/9.d0
     &     - 16*dlnx**3*xm1 +(4136*dlxm1*xm1)/27.d0 + (964*dlnx*dlxm1
     &     *xm1)/9.d0 -44*dlnx**2*dlxm1*xm1 - (841*dlxm1**2*xm1)/9.d0
     &     +72*dlnx*dlxm1**2*xm1 + (1153*xm1**2)/72.d0 + (1471*dlnx*xm1
     &     **2)/36.d0 +(74*dlnx**2*xm1**2)/9.d0 + (115*dlnx**3*xm1**2)
     &     /27.d0 -(815*dlxm1*xm1**2)/27.d0 - 38*dlnx*dlxm1*xm1**2 +(52
     &     *dlnx**2*dlxm1*xm1**2)/3.d0 + (325*dlxm1**2*xm1**2)/12.d0
     &     -(553*dlnx*dlxm1**2*xm1**2)/18.d0 + (367*dlxm1**3*xm1**2)/54
     &     .d0 -(1537*xm1**3)/486.d0 - (616*dlnx*xm1**3)/81.d0 -(113
     &     *dlnx**2*xm1**3)/27.d0 + (392*dlxm1*xm1**3)/81.d0 +(400*dlnx
     &     *dlxm1*xm1**3)/27.d0 - 8*dlxm1**2*xm1**3 + (29*z2)/6.d0 

      sqgaexact = sqgaexact + (629*dlnx*z2)/9.d0 - (50*dlxm1*z2)/9.d0 +
     &     (701*xm1*z2)/9.d0 - 72*dlnx*xm1*z2 - (281*xm1**2*z2)/9.d0 +
     &     (71*dlnx*xm1**2*z2)/3.d0 - (50*dlxm1*xm1**2*z2)/9.d0 + 8*xm1
     &     **3*z2 + (311*z3)/18.d0 + (311*xm1**2*z3)/18.d0 + (-37
     &     .333333333333336d0 + (761*dlnx)/9.d0 - (322*dlxm1)/3.d0 -
     &     (209*xm1)/9.d0 - 96*dlnx*xm1 + 128*dlxm1*xm1 + (59*xm1**2)/18
     &     .d0 + (245*dlnx*xm1**2)/9.d0 - (278*dlxm1*xm1**2)/9.d0 + (26
     &     *xm1**3)/9.d0)*dli2a + (7.87037037037037d0 + 10*dlnx - 10
     &     *dlxm1 - (50*xm1)/9.d0 - 8*dlnx*xm1 + 8*dlxm1*xm1 + (5*xm1**2
     &     )/6.d0 + 2*dlnx*xm1**2 - 2*dlxm1*xm1**2 - (2*xm1**3)/27.d0)
     &     *dli2b + (22.444444444444443d0 - 32*xm1 + (2*xm1**2) /9.d0)
     &     *dli3a + (20 - 16*xm1 + 4*xm1**2)*dli3b + (-20 + 16*xm1 - 4
     &     *xm1**2)*dli3c + (-123 .88888888888889d0 + 128*xm1 - (113*xm1
     &     **2)/3.d0)*dli3d + (-2.5d0 + 2*xm1 - xm1**2/2.d0)*dli3e + (-2
     &     .5d0 + 2*xm1 - xm1**2/2.d0)*dli3f

c..   lfh-terms:
      sqgaexact = sqgaexact + lfh**2*(-1.5d0 - (160*dlnx)/9.d0 + (31
     &     *dlxm1)/9.d0- (191*xm1)/9.d0 +18*dlnx*xm1 + (50*xm1**2)/9.d0
     &     - (56*dlnx*xm1**2)/9.d0 +(31*dlxm1*xm1**2)/9.d0 - 2*xm1**3)
     &     +lfh*(-10.444444444444445d0 + (263*dlnx)/27.d0 - (185*dlnx
     &     **2)/9.d0-dlxm1/3.d0 + 76*dlnx*dlxm1 - (31*dlxm1**2)/3.d0 -
     &     (2092*xm1)/27.d0 -(152*dlnx*xm1)/3.d0 + 22*dlnx**2*xm1 + (836
     &     *dlxm1*xm1)/9.d0 -72*dlnx*dlxm1*xm1 + (743*xm1**2)/54.d0 +
     &     (52*dlnx*xm1**2)/3.d0 -(67*dlnx**2*xm1**2)/9.d0 - (218*dlxm1
     &     *xm1**2)/9.d0+(268*dlnx*dlxm1*xm1**2)/9.d0 - (31*dlxm1**2*xm1
     &     **2)/3.d0-(196*xm1**3)/81.d0 - (200*dlnx*xm1**3)/27.d0 + 8
     &     *dlxm1*xm1**3 +(25*z2)/9.d0 + (25*xm1**2*z2)/9.d0
     &     +(54.22222222222222d0- 64*xm1 + 16*xm1**2)*dli2a+(5 - 4*xm1 +
     &     xm1**2)*dli2b)

c..   lfr-terms:
      sqgaexact = sqgaexact + lfr*(-5.5d0 + (11*dlnx)/2.d0 - 11*dlxm1 +
     &     11*xm1 + (11*xm1**2)/4.d0 +(11*dlnx*xm1**2)/2.d0 - 11*dlxm1
     &     *xm1**2 + lfh*(5.5d0 + (11*xm1**2)/2.d0))

c..   turn to scalar:
      if (.not.lpseudo) then
         sqgaexact = sqgaexact - (
     &        0.3333333333333333d0 + 17*dlnx - (28*dlnx**2)/9.d0 + (2
     &        *dlxm1)/3.d0 -lfh/3.d0 + (140*xm1)/9.d0 - (28*dlnx*xm1)/3
     &        .d0 + (28*dlnx**2*xm1)/9.d0 +(17*xm1**2)/6.d0 - (dlnx*xm1
     &        **2)/3.d0 + (2*dlxm1*xm1**2)/3.d0 -(lfh*xm1**2)/3.d0 )
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function sqgfexact:

      real*8 function sqgfexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

c..   pseudo-scalar:
      sqgfexact = 0.16049382716049382d0 + (10*dlnx)/27.d0 + dlnx**2/9.d0
     &     - (2*dlxm1)/3.d0 -(2*dlnx*dlxm1)/9.d0 + dlxm1**2/18.d0 + (10
     &     *xm1)/27.d0 + (2*dlxm1*xm1)/9.d0 +(179*xm1**2)/162.d0 + (19
     &     *dlnx*xm1**2)/27.d0 + (dlnx**2*xm1**2)/9.d0 -dlxm1*xm1**2 -
     &     (2*dlnx*dlxm1*xm1**2)/9.d0 + (dlxm1**2*xm1**2)/18.d0

c..   lfh-terms:
      sqgfexact = sqgfexact + lfh**2*(0.1111111111111111d0 + xm1**2
     &     /9.d0)+lfh*(0.37037037037037035d0 + (2*dlnx)/9.d0 - (2*dlxm1)
     &     /9.d0+ (19*xm1**2)/27.d0 + (2*dlnx*xm1**2)/9.d0 - (2*dlxm1
     &     *xm1**2)/9.d0)

c..   lfr-terms:
      sqgfexact = sqgfexact + lfr*(0.3333333333333333d0 - dlnx/3.d0 + (2
     &     *dlxm1)/3.d0 - (2*xm1)/3.d0 -xm1**2/6.d0 - (dlnx*xm1**2)/3.d0
     &     + (2*dlxm1*xm1**2)/3.d0 +lfh*(-0.3333333333333333d0 - xm1**2
     &     /3.d0))

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function sqqaexact:

      real*8 function sqqaexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      sqqaexact = (304*dlnx)/27.d0 + (4*dlnx**2)/27.d0 + (328*dlnx**3)
     &     /81.d0 -(8*dlnx*dlxm1)/9.d0 + 8*dlnx**2*dlxm1 - 16*dlnx*dlxm1
     &     **2-(476*xm1)/27.d0 - (700*dlnx*xm1)/27.d0 - (112*dlnx**2
     &     *xm1)/27.d0 -(8*dlnx**3*xm1)/3.d0 + (88*dlxm1*xm1)/3.d0 + 16
     &     *dlnx*dlxm1*xm1 -(16*dlnx**2*dlxm1*xm1)/3.d0 - (128*dlxm1**2
     &     *xm1)/9.d0 +(32*dlnx*dlxm1**2*xm1)/3.d0 + (44*xm1**2)/9.d0 +
     &     (46*dlnx*xm1**2)/9.d0 +(4*dlnx**2*xm1**2)/3.d0 + (40*dlnx**3
     &     *xm1**2)/81.d0 -(20*dlxm1*xm1**2)/3.d0 - (40*dlnx*dlxm1*xm1
     &     **2)/9.d0 +(8*dlnx**2*dlxm1*xm1**2)/9.d0 + (32*dlxm1**2*xm1
     &     **2)/9.d0 -(16*dlnx*dlxm1**2*xm1**2)/9.d0 + 16*dlnx*z2 + (128
     &     *xm1*z2)/9.d0 -(32*dlnx*xm1*z2)/3.d0 - (32*xm1**2*z2)/9.d0 +
     &     (16*dlnx*xm1**2*z2)/9.d0 +(-0.8888888888888888d0 + (656*dlnx)
     &     /27.d0 - 32*dlxm1 - (16*xm1)/3.d0 -16*dlnx*xm1 + (64*dlxm1
     &     *xm1)/3.d0 + (8*xm1**2)/9.d0 +(80*dlnx*xm1**2)/27.d0 - (32
     &     *dlxm1*xm1**2)/9.d0)*dli2a +(-0.2962962962962963d0 - (8*xm1
     &     **2)/27.d0)*dli3a +(-32.2962962962963d0 + (64*xm1)/3.d0 -
     &     (104*xm1**2)/27.d0)*dli3d

c..   lfh-terms:
      sqqaexact = sqqaexact + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8*dlnx
     &     *xm1)/3.d0 + (8*xm1**2) /9.d0 -(4*dlnx*xm1**2)/9.d0) +lfh*((4
     &     *dlnx)/9.d0 - 4*dlnx**2 + 16 *dlnx*dlxm1 - (44*xm1)/3.d0 -8
     &     *dlnx*xm1 + (8*dlnx**2*xm1)/3.d0 + (128*dlxm1*xm1)/9.d0 -(32
     &     *dlnx*dlxm1*xm1)/3.d0 + (10*xm1**2)/3.d0 + (20*dlnx*xm1**2)/9
     &     .d0 -(4*dlnx**2*xm1**2)/9.d0 - (32*dlxm1*xm1**2)/9.d0 +(16
     &     *dlnx*dlxm1*xm1**2)/9.d0 +(16 - (32*xm1)/3.d0 + (16*xm1 **2)
     &     /9.d0)*dli2a)

c..   there are no lfr-terms

c..   turn to scalar:
      if (.not.lpseudo) then
         sqqaexact = sqqaexact - ((272*dlnx)/27.d0 - (64*dlnx**2)/27.d0
     &        +(272*xm1)/27.d0 - (176*dlnx*xm1)/27.d0 + (64*dlnx**2*xm1)
     &        /27.d0 + (8*xm1**2)/9.d0 )
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function sqqfexact:

      real*8 function sqqfexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   pseudo-scalar
      sqqfexact = 0.d0 * xx

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function squaexact:

      real*8 function squaexact(xx)
c..   
c..   qq' contribution (different quarks)
c..   
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      squaexact = 12*dlnx + 4*dlnx**3 - (8*dlnx*dlxm1)/9.d0 + 8*dlnx**2
     &     *dlxm1 -16*dlnx*dlxm1**2 - (152*xm1)/9.d0 - 28*dlnx*xm1 -(32
     &     *dlnx**2*xm1)/9.d0 - (8*dlnx**3*xm1)/3.d0 + (88*dlxm1*xm1)/3
     &     .d0 +16 *dlnx*dlxm1*xm1 - (16*dlnx**2*dlxm1*xm1)/3.d0 -(128
     &     *dlxm1**2 *xm1)/9.d0 + (32*dlnx*dlxm1**2*xm1)/3.d0 + (10*xm1
     &     **2)/3.d0 +(58 *dlnx*xm1**2)/9.d0 + (8*dlnx**2*xm1**2)/9.d0
     &     +(4*dlnx**3*xm1**2) /9.d0 - (20*dlxm1*xm1**2)/3.d0 -(40*dlnx
     &     *dlxm1*xm1**2)/9.d0 + (8 *dlnx**2*dlxm1*xm1**2)/9.d0 +(32
     &     *dlxm1**2*xm1**2)/9.d0 - (16*dlnx *dlxm1**2*xm1**2)/9.d0 +16
     &     *dlnx*z2 + (128*xm1*z2)/9.d0 - (32*dlnx *xm1*z2)/3.d0 -(32
     &     *xm1**2*z2)/9.d0 + (16*dlnx*xm1**2*z2)/9.d0 +(-0
     &     .8888888888888888d0 + 24*dlnx - 32*dlxm1 - (16*xm1)/3.d0 -16
     &     *dlnx *xm1 + (64*dlxm1*xm1)/3.d0 + (8*xm1**2)/9.d0 +(8*dlnx
     &     *xm1**2)/3.d0 - (32*dlxm1*xm1**2)/9.d0)*dli2a +(-32 + (64
     &     *xm1)/3.d0 - (32 *xm1**2)/9.d0)*dli3d

c..   lfh-terms:
      squaexact = squaexact + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8
     &     *dlnx*xm1)/3.d0 + (8*xm1**2) /9.d0 -(4*dlnx*xm1**2)/9.d0)
     &     +lfh*((4*dlnx)/9.d0 - 4*dlnx**2 + 16 *dlnx*dlxm1 - (44*xm1)
     &     /3.d0 -8*dlnx*xm1 + (8*dlnx**2*xm1)/3.d0 + (128*dlxm1*xm1)
     &     /9.d0 -(32*dlnx*dlxm1*xm1)/3.d0 + (10*xm1**2)/3.d0 + (20*dlnx
     &     *xm1**2)/9.d0 -(4*dlnx**2*xm1**2)/9.d0 - (32*dlxm1*xm1**2)
     &     /9.d0 +(16*dlnx*dlxm1*xm1**2)/9.d0 +(16 - (32*xm1)/3.d0 + (16
     &     *xm1 **2)/9.d0)*dli2a)
      
c..   there are no lfr-terms

c..   turn to scalar:
      if (.not.lpseudo) then
         squaexact = squaexact - ((80*dlnx)/9.d0 - (16*dlnx**2)/9.d0 +
     &        (80*xm1)/9.d0 - (16*dlnx *xm1)/3.d0+(16*dlnx**2*xm1)/9.d0
     &        + (8*xm1**2)/9.d0 )
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function squfexact:

      real*8 function squfexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

c..   pseudo-scalar:
      squfexact = 0.d0 * xx

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function sqqbaexact:

      real*8 function sqqbaexact(xx)
c..
c..   q q-bar contribution
c..
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = sudilog(xm1)
      dli2b = sudilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar
      sqqbaexact = (64*dlnx)/9.d0 - (136*dlnx**2)/81.d0 + (88*dlnx**3)
     &     /27.d0 +(184*dlnx*dlxm1)/81.d0 + 8*dlnx**2*dlxm1 - 16*dlnx
     &     *dlxm1**2 -(2020*xm1)/81.d0 - (1708*dlnx*xm1)/81.d0 - (8*dlnx
     &     **2*xm1)/27.d0 -(56*dlnx**3*xm1)/27.d0 + (2632*dlxm1*xm1)
     &     /81.d0+(304*dlnx*dlxm1*xm1)/27.d0 - (16*dlnx**2*dlxm1*xm1)
     &     /3.d0-(128*dlxm1**2*xm1)/9.d0 + (32*dlnx*dlxm1**2*xm1)/3.d0
     &     +(860*xm1**2)/81.d0 - (386*dlnx*xm1**2)/81.d0 -(136*dlnx**2
     &     *xm1**2)/27.d0 + (8*dlnx**3*xm1**2)/27.d0 -(796*dlxm1*xm1**2)
     &     /81.d0+ (8*dlnx*dlxm1*xm1**2)/27.d0 +(8*dlnx**2*dlxm1*xm1**2)
     &     /9.d0+ (32*dlxm1**2*xm1**2)/9.d0 -(16*dlnx*dlxm1**2*xm1**2)
     &     /9.d0 +(1060*xm1**3)/27.d0 +(512*dlnx*xm1**3)/27.d0 + (352
     &     *dlnx**2*xm1**3)/81.d0 -(512*dlxm1*xm1**3)/27.d0 - (512*dlnx
     &     *dlxm1*xm1**3)/81.d0 +(416*dlxm1**2*xm1**3)/81.d0 + 16*dlnx
     &     *z2 +(128*xm1*z2)/9.d0 -(32*dlnx*xm1*z2)/3.d0 - (32*xm1**2
     &     *z2)/9.d0 + (16*dlnx*xm1**2*z2)/9.d0 -(688*xm1**3*z2)/81.d0 

      sqqbaexact = sqqbaexact + (7 .604938271604938d0 + (256*dlnx)/9.d0
     &     -32*dlxm1 - (208*xm1)/9 .d0 - (176*dlnx*xm1)/9.d0 + (64*dlxm1
     &     *xm1)/3.d0 + (152*xm1**2)/9.d0 + (32*dlnx*xm1**2)/9.d0 - (32
     &     *dlxm1*xm1**2)/9.d0 - (208*xm1**3)/81.d0)*dli2a + (-2
     &     .6666666666666665d0 - (20*dlnx)/9.d0 + (176*xm1)/27.d0 + (16
     &     *dlnx*xm1)/9.d0 - (152*xm1**2)/27.d0 - (4*dlnx*xm1**2)/9.d0 +
     &     (16*xm1**3)/9.d0)*dli2b + (-7.407407407407407d0 + (160*xm1)
     &     /27.d0 - (40*xm1**2)/27.d0)*dli3a + (-1 .4814814814814814d0 +
     &     (32*xm1)/27.d0 - (8*xm1**2)/27.d0) * dli3b + (1
     &     .4814814814814814d0 - (32*xm1)/27 .d0 + (8*xm1**2)/27.d0)
     &     * dli3c + (-36 .44444444444444d0 + (224*xm1)/9.d0 - (40*xm1
     &     **2)/9.d0)*dli3d + (1.8518518518518519d0 - (40*xm1)/27.d0 +
     &     (10*xm1 **2)/27.d0)* dli3e + (1.1111111111111112d0 - (8*xm1)
     &     /9.d0 + (2*xm1**2)/9.d0)* dli3f

c..   lfh-terms:
      sqqbaexact = sqqbaexact + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8
     &     *dlnx*xm1)/3.d0 + (8*xm1**2)/9.d0 -(4*dlnx*xm1**2)/9.d0) +lfh
     &     *((-220 *dlnx)/81.d0 - 4*dlnx**2 + 16*dlnx*dlxm1 - (1444*xm1)
     &     /81.d0 -(88 *dlnx*xm1)/27.d0 + (8*dlnx**2*xm1)/3.d0 + (128
     &     *dlxm1*xm1)/9.d0 -(32 *dlnx*dlxm1*xm1)/3.d0 + (526*xm1**2)
     &     /81.d0-(68*dlnx*xm1**2)/27.d0 - (4*dlnx**2*xm1**2)/9.d0 -(32
     &     *dlxm1*xm1**2)/9.d0 + (16*dlnx *dlxm1*xm1**2)/9.d0 +(88*xm1
     &     **3)/9.d0 + (256*dlnx*xm1**3)/81.d0 - (256*dlxm1*xm1**3)
     &     /81.d0 +(16- (32*xm1)/3.d0 + (16*xm1**2)/9.d0) *dli2a)

c..   lfr-terms:
      sqqbaexact = sqqbaexact + (-88*lfr*xm1**3)/9.d0

c..   turn to scalar:
      if (.not.lpseudo) then
         sqqbaexact = sqqbaexact - ((352*dlnx)/27.d0 + (32*dlnx**2)
     &        /27.d0+(352*xm1)/27.d0 - (256*dlnx*xm1)/27.d0 - (32*dlnx
     &        **2*xm1)/27.d0- (64*xm1**2)/9.d0 + (16*xm1**3)/27.d0 )
c..   no additional lfr- and lfh-terms
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function sqqbfexact:

      real*8 function sqqbfexact(xx)
      implicit real*8 (a-h,o-z)
      complex*16 sudilog,trilog
      external sudilog,trilog
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

c..   pseudo-scalar:
      sqqbfexact = (-32*dlnx)/81. - (32*xm1)/81. + (16*dlnx*xm1)/27. +
     &     (32*xm1**2)/81. - (328*xm1**3)/243. - (64*dlnx*xm1**3)/81. +
     &     (32*dlxm1*xm1**3)/81.

c..   lfh-terms:
      sqqbfexact = sqqbfexact + (-16*lfh*xm1**3)/27.d0

c..   lfr-terms:
      sqqbfexact = sqqbfexact + (16*lfr*xm1**3)/27.d0

c..   turn to scalar:
      if (.not.lpseudo) then
         sqqbfexact = sqqbfexact - ((-32*dlnx)/27.d0 - (32*xm1)/27.d0 +
     &        (32*dlnx*xm1)/27.d0 + (16 *xm1**2)/27.d0 )
c..   no additional lfh- and lfr-terms
      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}





