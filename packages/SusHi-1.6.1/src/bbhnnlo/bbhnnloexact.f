C-{{{ DESC:

c..   bbhnnloexact.f  is part of bbh@nnlo
c..   See README for a description of this package.

C-}}}
C-{{{ RCS:

c..   
c..   $Id: bbhnnloexact.f,v 1.3 2005/02/01 10:17:27 rharland Exp $
c..   $Log: bbhnnloexact.f,v $
c..   Revision 1.3  2005/02/01 10:17:27  rharland
c..   before taking out LHAPDF again
c..
c..   Revision 1.1  2004/11/16 11:54:37  rharland
c..   Initial revision
c..
c..

C-}}}
C-{{{ function sbbq2bbh(yy)

      real*8 function sbbq2bbh(xt)
C
C
      implicit none
      real*8 xt,deltabbqabbh,deltabbqfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sbbq2bbh = deltabbqabbh(xt) + nf*deltabbqfbbh(xt)

      end

C-}}}
C-{{{ function sbg2bbh(yy)

      real*8 function sbg2bbh(xt)
C
C
      implicit none
      real*8 xt,deltabgabbh,deltabgfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sbg2bbh = deltabgabbh(xt) + nf*deltabgfbbh(xt)

      end

C-}}}
C-{{{ function sgg2bbh(yy)

      real*8 function sgg2bbh(xt)
C
C
      implicit none
      real*8 xt,deltaggabbh,deltaggfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sgg2bbh = deltaggabbh(xt) + nf*deltaggfbbh(xt)

      end

C-}}}
C-{{{ function sbb2bbh(yy)

      real*8 function sbb2bbh(xt)
C
C
      implicit none
      real*8 xt,deltabbabbh,deltabbfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sbb2bbh = deltabbabbh(xt) + nf*deltabbfbbh(xt)

      end

C-}}}
C-{{{ function sbq2bbh(yy)

      real*8 function sbq2bbh(xt)
C
C
      implicit none
      real*8 xt,deltabqabbh,deltabqfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sbq2bbh = deltabqabbh(xt) + nf*deltabqfbbh(xt)

      end

C-}}}
C-{{{ function sqqb2bbh(yy)

      real*8 function sqqb2bbh(xt)
C
C
      implicit none
      real*8 xt,deltaqqbabbh,deltaqqbfbbh
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      sqqb2bbh = deltaqqbabbh(xt) + nf*deltaqqbfbbh(xt)

      end

C-}}}

c-{{{ function deltabbqabbh:

      real*8 function deltabbqabbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabbqabbh = 28.814814814814813d0 + (74*dli3a)/3.d0 
     &     + (2*dli3b)/3.d0
     &     - (2*dli3c)/3.d0 +(82*dli3d)/3.d0 - dli3f/18.d0 + dli2b
     &     *(0.3333333333333333d0 + dlnx/3.d0) +(1607*dlnx)/54.d0 + (157
     &     *dlnx**2)/9.d0 + (10*dlnx**3)/3.d0 -(403*dlxm1)/9.d0 - (154
     &     *dlnx*dlxm1)/3.d0 -(368*dlnx**2*dlxm1)/9.d0 +(88*dlxm1**2)/3
     &     .d0 + (584*dlnx*dlxm1**2)/9.d0 - (256*dlxm1**3)/9.d0 +dli2a*(
     &     -4.888888888888889d0 - (236*dlnx)/9.d0 + (136*dlxm1)/9.d0
     &     -(88*lfh)/9.d0) + (7.333333333333333d0+ (86*dlnx)/9.d0
     &     -(128*dlxm1)/9.d0)*lfh**2 + ((-76*dlnx)/3.d0 +(152*dlxm1)
     &     /3.d0)*lfr +lfh*(22.666666666666668d0 + (77*dlnx)/3.d0
     &     +(136*dlnx**2)/9.d0 -(88*dlxm1)/3.d0 - (536*dlnx*dlxm1)/9.d0
     &     + (128*dlxm1**2)/3.d0 -(76*lfr)/3.d0 - (200*z2)/9.d0) -
     &     (88*z2)/3.d0 - (416*dlnx*z2)/9.d0 +(400*dlxm1*z2)/9.d0 + ((
     &     -142*dli3a)/9.d0 - (64*dli3d)/9.d0 -(146*dlnx)/9.d0 - (23
     &     *dlnx**2)/3.d0- (44*dlnx**3)/27.d0 +24*dlnx*dlxm1 + (148*dlnx
     &     **2*dlxm1)/9.d0 - (248*dlnx*dlxm1**2)/9.d0 +dli2a*(1 + (58
     &     *dlnx)/9.d0 +(20*dlxm1)/9.d0) +(-12*dlnx- (50*dlnx**2)/9.d0 +
     &     (224*dlnx*dlxm1)/9.d0)*lfh -(32*dlnx*lfh**2)/9.d0 + (38
     &     *dlnx*lfr)/3.d0+ (164*dlnx*z2)/9.d0)/xm1 -(764*z3)/9.d0 

      deltabbqabbh = deltabbqabbh + xm1**2*(15.805555555555555d0 
     &     + (55*dli3a)/9.d0
     &     + 2*dli3b -2*dli3c + (76*dli3d)/9.d0 - dli3f/6.d0 + (31
     &     *dlnx)/3.d0+(79*dlnx**2)/18.d0 + (11*dlnx**3)/9.d0 +dli2b*(0
     &     .8888888888888888d0 + dlnx) - (185*dlxm1)/9.d0 -(88*dlnx
     &     *dlxm1)/9.d0 - (110*dlnx**2*dlxm1)/9.d0 + (68*dlxm1**2)/9.d0
     &     +(56*dlnx*dlxm1**2)/3.d0 - (64*dlxm1**3)/9.d0 +dli2a*(2
     &     .2222222222222223d0 - (94*dlnx)/9.d0 + (26*dlxm1)/3.d0 -(44
     &     *lfh)/9.d0) + (2.7777777777777777d0+ 3*dlnx - (32*dlxm1)
     &     /9.d0)*lfh**2 + (6.333333333333333d0 - (19*dlnx)/3.d0 +
     &     (38*dlxm1)/3.d0)*lfr + lfh*(10.555555555555555d0 +(16
     &     *dlnx)/3.d0 + (43*dlnx**2)/9.d0 -(68*dlxm1)/9.d0 - (52*dlnx
     &     *dlxm1)/3.d0 + (32*dlxm1**2)/3.d0 -(19*lfr)/3.d0 - (50*z2)
     &     /9.d0) -(41*z2)/3.d0- 14*dlnx*z2 +(100*dlxm1*z2)/9.d0 - (191
     &     *z3)/9.d0 - (7*dli3e)/6.d0)

      deltabbqabbh = deltabbqabbh + xm1**3*(-2.2253086419753085d0 
     &     - (2*dli3a)/3.d0
     &     - (2*dli3b)/3.d0 + (2*dli3c)/3.d0 + (2*dli3d)/3.d0 +
     &     dli3f/18.d0 +dli2b*(-0.2222222222222222d0 - dlnx/3.d0) +
     &     dli2a*(-1.1111111111111112d0 + (2*dlnx)/9.d0) - (19*dlnx)/27
     &     .d0 + (2*dlnx**2)/3.d0 - (4*dlnx**3)/27.d0 + (44*dlxm1)/27.d0
     &     - (8*dlnx*dlxm1)/3.d0 +(8*dlxm1**2)/9.d0 + (-0
     &     .8148148148148148d0 + (4*dlnx)/3.d0 - (8*dlxm1)/9.d0)*lfh
     &     +(2*lfh**2)/9.d0 - (8*z2)/9.d0 + (7*dli3e)/18.d0)- (7
     &     *dli3e)/18.d0 + xm1*(-28.074074074074073d0 - (43*dli3a)/3.d0
     &     - 2*dli3b + 2*dli3c - (88*dli3d)/3.d0 + dli3f/6.d0 + dli2b*(
     &     -1 -dlnx) - (139*dlnx)/6.d0 -(89*dlnx**2)/6.d0 - (25*dlnx**3)
     &     /9.d0 + (124*dlxm1)/3.d0 + (358*dlnx*dlxm1)/9.d0 + (110*dlnx
     &     **2*dlxm1)/3.d0 - (200*dlxm1**2)/9.d0- 56*dlnx*dlxm1**2 + (64
     &     *dlxm1**3)/3.d0 + (-6.444444444444445d0 -9*dlnx + (32*dlxm1)
     &     /3.d0)*lfh**2 + dli2a*(4.555555555555555d0 +30*dlnx - 26
     &     *dlxm1 + (44*lfh)/3.d0) + (-6.333333333333333d0 + 19*dlnx
     &     -38*dlxm1)*lfr + (85*z2)/3.d0 + 42*dlnx*z2 - (100*dlxm1
     &     *z2)/3.d0 + lfh*(-21.22222222222222d0 - (61*dlnx)/3.d0 -
     &     (43*dlnx**2)/3.d0 + (200*dlxm1)/9.d0 + 52*dlnx*dlxm1 - 32
     &     *dlxm1**2 +19*lfr + (50*z2)/3.d0) + (191*z3)/3.d0 + (7
     &     *dli3e)/6.d0)

      end

c-}}}
c-{{{ function deltabbqfbbh:

      real*8 function deltabbqfbbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabbqfbbh = -1.382716049382716d0 
     &     + (2*dli2a)/9.d0 - (20*dlnx)/9.d0
     &     - (13*dlnx **2)/9.d0 +(80*dlxm1)/27.d0 + (32*dlnx*dlxm1)/9.d0
     &     - (16*dlxm1**2) /9.d0 -(4*lfh**2)/9.d0 + lfh*(-1
     &     .4814814814814814d0 -(16*dlnx)/9.d0 +(16*dlxm1)/9.d0 + (8
     &     *lfr)/9.d0) + ((8*dlnx)/9.d0 -(16*dlxm1)/9.d0) *lfr +((
     &     -2*dli2a)/9.d0 + (10*dlnx)/9.d0 + (2*dlnx**2)/3.d0 -(16 *dlnx
     &     *dlxm1)/9.d0 + (8*dlnx*lfh)/9.d0 - (4*dlnx*lfr)/9.d0)
     &     /xm1 +xm1 *(1.3703703703703705d0 + (17*dlnx)/9.d0 + (7*dlnx
     &     **2)/6.d0 -(8*dlxm1)/3.d0 - (8*dlnx*dlxm1)/3.d0 + (4*dlxm1
     &     **2)/3.d0+ lfh**2/3.d0 +lfh*(1.3333333333333333d0 + (4
     &     *dlnx)/3.d0 - (4*dlxm1)/3.d0 -(2*lfr)/3.d0) +
     &     (0.2222222222222222d0 - (2*dlnx)/3.d0 +(4*dlxm1)/3.d0)*lfr
     &     - (4*z2)/3.d0) +xm1**2*(-0.6790123456790124d0 -(7*dlnx)/9.d0
     &     - (7 *dlnx**2)/18.d0 +(32*dlxm1)/27.d0 + (8*dlnx*dlxm1)/9.d0
     &     - (4*dlxm1 **2)/9.d0 -lfh**2/9.d0 + lfh*(
     &     -0.5925925925925926d0 - (4*dlnx)/9.d0 +(4*dlxm1)/9.d0 + (2
     &     *lfr)/9.d0)+(-0.2222222222222222d0 + (2*dlnx) /9.d0 - (4
     &     *dlxm1)/9.d0)*lfr + (4*z2)/9.d0) + (16*z2)/9.d0

      end

c-}}}
c-{{{ function deltabgabbh:

      real*8 function deltabgabbh(xx)

      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabgabbh = -0.7916666666666666d0 
     &     - (31*dli3a)/3.d0 - (15*dli3b)/2.d0
     &     +(15*dli3c)/2.d0 + (113*dli3d)/12.d0 + (13*dli3e)/48.d0 +
     &     (13*dli3f)/48.d0 - (425*dlnx)/48.d0 - (157*dlnx**2)/32.d0 -
     &     (193*dlnx**3)/144.d0 + (13*dlxm1)/16.d0 + (331*dlnx*dlxm1)/24
     &     .d0 - (79*dlnx**2*dlxm1)/24.d0 - dlxm1**2/2.d0 + (21*dlnx
     &     *dlxm1**2)/8.d0 + (257*dlxm1**3)/144.d0 + dli2a*(16.375d0 -
     &     (9*dlnx)/2.d0 + (125*dlxm1)/24.d0 - (5*lfh)/2.d0) + dli2b
     &     *(-1.6666666666666667d0 - (37*dlnx)/12.d0 + (15*dlxm1)/4.d0 -
     &     (15*lfh)/8.d0) + (-0.125d0 + (9*dlnx)/8.d0 + (7*dlxm1)
     &     /8.d0)*lfh**2 + ((19*dlnx)/16.d0 - (19*dlxm1)/8.d0)*lfr
     &     + z2/2.d0 -(89*dlnx*z2)/24.d0 - (35*dlxm1*z2)/12.d0 + lfh*
     &     (-0.4583333333333333d0 - (163*dlnx)/24.d0 + 2*dlnx**2 + dlxm1
     &     /2.d0 -(37*dlnx*dlxm1)/12.d0 - (21*dlxm1**2)/8.d0 + (19
     &     *lfr)/16.d0 + (35*z2)/24.d0) 

      deltabgabbh = deltabgabbh + xm1*(7.75d0 
     &     + 19*dli3a + (33*dli3b)/2.d0 -
     &     (33*dli3c)/2.d0 - (65*dli3d)/4.d0 - dli3e/16.d0 - dli3f/16.d0
     &     + (625*dlnx)/16.d0 + (641*dlnx**2)/32.d0 + (121*dlnx**3)/48
     &     .d0 - (31*dlxm1)/2.d0 - (457*dlnx*dlxm1)/8.d0 + (19*dlnx**2
     &     *dlxm1)/8.d0 +(113*dlxm1**2)/8.d0 + (9*dlnx*dlxm1**2)/8.d0 -
     &     (257*dlxm1**3)/48.d0+ (3.25d0 - (9*dlnx)/8.d0 - (21*dlxm1)/8
     &     .d0)*lfh**2 + dli2a*(-40.375d0 + (9*dlnx)/2.d0 - (29
     &     *dlxm1)/8.d0 + (3*lfh)/2.d0) + dli2b*(3.375d0 + (25*dlnx)
     &     /4.d0 - (33*dlxm1)/4.d0 + (33*lfh)/8.d0) + (-2.375d0- (57
     &     *dlnx)/16.d0 +(57*dlxm1)/8.d0)*lfr + lfh
     &     *(8.166666666666666d0+ (227*dlnx)/8.d0 - (9*dlnx**2)/4.d0 -
     &     (43*dlxm1)/3.d0 + (dlnx*dlxm1)/4.d0+ (63*dlxm1**2)/8.d0 - (57
     &     *lfr)/16.d0 - (35*z2)/8.d0) -(93*z2)/8.d0 + (17*dlnx*z2)
     &     /8.d0 + (35*dlxm1*z2)/4.d0 - (161*z3)/16.d0)

      deltabgabbh = deltabgabbh + xm1**3*(6.184027777777778d0 
     &     + 3*dli3a + 3*dli3b - 3*dli3c 
     &     + (3*dli3d)/2.d0 + (7*dli3e)/24.d0 + (7
     &     *dli3f) /24.d0 + (331*dlnx)/16.d0 + (257*dlnx**2)/24.d0 +
     &     dlnx**3/12.d0 - (121*dlxm1)/8.d0 - (129*dlnx*dlxm1)/4.d0 -
     &     (49*dlnx**2*dlxm1)/12.d0 + (237*dlxm1**2)/16.d0 + (31*dlnx
     &     *dlxm1**2)/4.d0 - (257*dlxm1**3) /72.d0 + dli2a*(-9
     &     .833333333333334d0 - (17*dlnx)/6.d0 + (43*dlxm1) /12.d0 - 2
     &     *lfh) + dli2b*(0.4375d0 + (5*dlnx)/6.d0 - (3*dlxm1)/2.d0 +
     &     (3*lfh)/4.d0) + (2.9375d0 + dlnx - (7*dlxm1)/4.d0)*lfh
     &     **2+ (-4.15625d0 - (19*dlnx)/8.d0 + (19*dlxm1)/4.d0)*lfr
     &     +lfh*(7.270833333333333d0 + (46*dlnx)/3.d0 + (4*dlnx**2)
     &     /3.d0 - (353*dlxm1)/24.d0 - (41*dlnx*dlxm1)/6.d0 + (21*dlxm1
     &     **2)/4.d0 -(19*lfr)/8.d0 - (35*z2)/12.d0) - (169*z2)/16.d0
     &     - (67*dlnx*z2)/12.d0 + (35*dlxm1*z2)/6.d0 - (161*z3)/24.d0)
     &     +(161*z3)/48.d0 

      deltabgabbh = deltabgabbh + xm1**2*(-13.864583333333334d0 
     &     - (35*dli3a)/3.d0- 12*dli3b + 12*dli3c 
     &     + (16*dli3d)/3.d0 - dli3e/2.d0 -
     &     dli3f/2.d0 -(2443*dlnx)/48.d0 - (155*dlnx**2)/6.d0 - (91*dlnx
     &     **3)/72.d0 + (1447*dlxm1)/48.d0 + (907*dlnx*dlxm1)/12.d0 + 5
     &     *dlnx**2*dlxm1 - (439*dlxm1**2)/16.d0 - (23*dlnx*dlxm1**2)/2
     &     .d0 + (257*dlxm1**3)/36.d0 +dli2b*(-2.1458333333333335d0 - 4
     &     *dlnx + 6*dlxm1 - 3*lfh) + (-5.8125d0 - dlnx + (7*dlxm1)
     &     /2.d0)*lfh**2 + dli2a*(35.833333333333336d0 + (17*dlnx)
     &     /6.d0 - (31*dlxm1)/6.d0 + 3*lfh) + (6.53125d0 + (19*dlnx)
     &     /4.d0 - (19*dlxm1)/2.d0)*lfr + (331*z2)/16.d0 +(43*dlnx
     &     *z2)/6.d0 - (35*dlxm1*z2)/3.d0 + lfh*(
     &     -15.145833333333334d0- (443*dlnx)/12.d0 - (13*dlnx**2)/12.d0
     &     + (661*dlxm1)/24.d0 + (29*dlnx*dlxm1)/3.d0 - (21*dlxm1**2)
     &     /2.d0 + (19*lfr)/4.d0 + (35*z2)/6.d0) +(161*z3)/12.d0)

      end

c-}}}
c-{{{ function deltabgfbbh:

      real*8 function deltabgfbbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabgfbbh = (-dlnx/24.d0 + dlxm1/12.d0)*lfr 
     &     - (lfh*lfr)/24.d0 + ((0.08333333333333333d0 
     &     + dlnx/8.d0 - dlxm1/4.d0)*lfr + (lfh*lfr
     &     )/8.d0) * xm1 + ((-0.22916666666666666d0 - dlnx/6.d0 + dlxm1
     &     /3.d0)*lfr - (lfh*lfr)/6.d0)*xm1**2 +
     &     ((0.14583333333333334d0+ dlnx/12.d0 - dlxm1/6.d0)*lfr +
     &     (lfh*lfr)/12.d0)* xm1**3

      end

c-}}}
c-{{{ function deltaggabbh:

      real*8 function deltaggabbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltaggabbh = (11*dli3a)/4.d0 + (5*dli3b)/32.d0 
     &     - (5*dli3c)/32.d0 -
     &     (105*dli3d)/16.d0 -(89*dli3e)/128.d0 - (79*dli3f)/128.d0
     &     +dli2b*(-0.15625d0 + (47*dlnx)/64.d0) - (143*dlnx)/64.d0 -(21
     &     *dlnx**2)/16.d0 +(11*dlnx**3)/8.d0 + (17*dlnx*dlxm1)/4.d0 +(9
     &     *dlnx**2*dlxm1)/4.d0 -(9*dlnx*dlxm1**2)/2.d0 +((-17*dlnx)/8
     &     .d0 - (9*dlnx**2)/8.d0 + (9*dlnx*dlxm1)/2.d0)*lfh -(9*dlnx
     &     *lfh**2)/8.d0 + dli2a*(4.5625d0 +(85*dlnx)/16.d0 - 9*dlxm1
     &     +(9*lfh)/2.d0) + (9*dlnx*z2)/2.d0 +xm1*(-14.484375d0 - (61
     &     *dli3a)/8.d0 - (11*dli3b)/32.d0 + (11*dli3c)/32.d0+(225
     &     *dli3d)/16.d0 + (247*dli3e)/128.d0 + (225*dli3f)/128.d0
     &     +dli2b*(0.5d0 - (129*dlnx)/64.d0) + (9*dlnx)/64.d0 + (109
     &     *dlnx**2)/64.d0-(53*dlnx**3)/16.d0 + (49*dlxm1)/4.d0 - (25
     &     *dlnx*dlxm1)/4.d0 -(21*dlnx**2*dlxm1)/4.d0 - 4*dlxm1**2 + (21
     &     *dlnx*dlxm1**2)/2.d0+dli2a*(-13.25d0 - (189*dlnx)/16.d0 + 21
     &     *dlxm1 - (21*lfh)/2.d0) +(-6.125d0 + (25*dlnx)/8.d0 + (21
     &     *dlnx**2)/8.d0 +4*dlxm1 -(21*dlnx*dlxm1)/2.d0)*lfh + (-1 +
     &     (21*dlnx)/8.d0)*lfh**2 + 4*z2 -(21*dlnx*z2)/2.d0) 

      deltaggabbh = deltaggabbh + xm1**3* (-13.2109375d0 
     &     - (9*dli3a)/4.d0 -
     &     dli3b/16.d0 + dli3c/16.d0 + (15*dli3d)/8.d0 + (37*dli3e)/64
     &     .d0 + (35*dli3f)/64.d0 + dli2b*(0.1875d0 - (19*dlnx)/32.d0) -
     &     (223*dlnx)/64.d0 - (43*dlnx**2)/64.d0 - (17*dlnx**3)/24.d0 +
     &     (75*dlxm1)/8.d0 + dlnx*dlxm1- dlnx**2*dlxm1 - 3*dlxm1**2 + 2
     &     *dlnx*dlxm1**2 + dli2a*(-3.875d0 - (15*dlnx)/8.d0 + 4*dlxm1 -
     &     2*lfh) + (-4.6875d0 - dlnx/2.d0 +dlnx**2/2.d0 + 3*dlxm1 -
     &     2*dlnx*dlxm1)*lfh + (-0.75d0 + dlnx/2.d0)*lfh**2 + 3*z2
     &     - 2*dlnx*z2)

      deltaggabbh = deltaggabbh + xm1**2*(27.6953125d0 
     &     + (57*dli3a)/8.d0 +
     &     dli3b/4.d0 - dli3c/4.d0 -(75*dli3d)/8.d0 - (29*dli3e)/16.d0 -
     &     (27*dli3f)/16.d0+(357*dlnx)/64.d0 + (9*dlnx**2)/32.d0 + (127
     &     *dlnx**3)/48.d0 +dli2b*(-0.53125d0 + (15*dlnx)/8.d0) - (173
     &     *dlxm1)/8.d0 + dlnx*dlxm1 +4*dlnx**2*dlxm1 + 7*dlxm1**2 - 8
     &     *dlnx*dlxm1**2 +(10.8125d0 -dlnx/2.d0 - 2*dlnx**2 - 7*dlxm1 +
     &     8*dlnx*dlxm1)*lfh +(1.75d0 - 2*dlnx)*lfh**2 +dli2a
     &     *(12.5625d0+ (67*dlnx)/8.d0 - 16*dlxm1 + 8*lfh) - 7*z2 +8
     &     *dlnx*z2)

      end

c-}}}
c-{{{ function deltaggfbbh:

      real*8 function deltaggfbbh(xx)

      implicit none
      real*8 xx

      deltaggfbbh = 0.d0

      end

c-}}}
c-{{{ function deltabbabbh:

      real*8 function deltabbabbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabbabbh = (31*dli2a)/9.d0 + dli2b/9.d0 - (10*dli3a)/9.d0 + (20
     &     *dli3b)/9.d0 -(20*dli3c)/9.d0 + (5*dli3d)/3.d0 - dli3e/3.d0 -
     &     dli3f/9.d0 - (71*dlnx)/27.d0 -(49*dli2a*dlnx)/9.d0 + (8*dli2b
     &     *dlnx)/9.d0 - (23*dlnx**2)/18.d0 +(11*dlnx**3)/18.d0 + (64
     &     *dli2a*dlxm1)/9.d0 - (8*dli2b*dlxm1)/9.d0 +(34*dlnx*dlxm1)/9
     &     .d0 - (28*dlnx**2*dlxm1)/9.d0 + (8*dlnx*dlxm1**2)/3.d0 -(32
     &     *dli2a*lfh)/9.d0 + (4*dli2b*lfh)/9.d0 -(17*dlnx*lfh)
     &     /9.d0+(14*dlnx**2*lfh)/9.d0 - (8*dlnx*dlxm1*lfh)/3.d0+
     &     (2*dlnx*lfh**2)/3.d0 -(4*dli3a)/(3.d0*(2 - xm1)) - (20
     &     *dli3b)/(9.d0*(2 -xm1)) +(20*dli3c)/(9.d0*(2 - xm1)) + (4
     &     *dli3d)/(3.d0*(2- xm1)) +dli3e/(3.d0*(2 - xm1)) + dli3f/(9.d0
     &     *(2 - xm1)) +(4*dli2a*dlnx)/(3.d0*(2 - xm1)) - (8*dli2b*dlnx)
     &     /(9.d0*(2 - xm1))-(14*dlnx**3)/(27.d0*(2 - xm1)) - (16*dli2a
     &     *dlxm1)/(9.d0*(2 - xm1)) +(8*dli2b*dlxm1)/(9.d0*(2 - xm1)) +
     &     (4*dlnx**2*dlxm1)/(9.d0*(2- xm1)) +(8*dli2a*lfh)/(9.d0*(2
     &     - xm1)) - (4*dli2b*lfh)/(9.d0*(2- xm1)) -(2*dlnx**2
     &     *lfh)/(9.d0*(2 - xm1)) + (16*xm1)/27.d0 - (53*dli2a*xm1)
     &     /9.d0 -(dli2b*xm1)/9.d0 + (17*dli3a*xm1)/9.d0 - (13*dli3b
     &     *xm1)/9.d0 +(13*dli3c*xm1)/9.d0 - (40*dli3d*xm1)/9.d0 + (13
     &     *dli3e*xm1)/36.d0 +(dli3f*xm1)/12.d0 

      deltabbabbh = deltabbabbh + (127*dlnx*xm1)/18.d0 
     &     + (22*dli2a*dlnx*xm1)/3
     &     .d0-(11*dli2b*dlnx*xm1)/18.d0 + (125*dlnx**2*xm1)/36.d0 -(11
     &     *dlnx**3*xm1)/27.d0 - (26*dlxm1*xm1)/9.d0 - (80*dli2a*dlxm1
     &     *xm1)/9.d0 +(4*dli2b*dlxm1*xm1)/9.d0 - (34*dlnx*dlxm1*xm1)/3
     &     .d0 +(38*dlnx**2*dlxm1*xm1)/9.d0 + (10*dlxm1**2*xm1)/3.d0 -4
     &     *dlnx*dlxm1**2*xm1 +(13*lfh*xm1)/9.d0 + (40*dli2a*lfh
     &     *xm1)/9.d0 -(2*dli2b*lfh*xm1)/9.d0+ (17*dlnx*lfh*xm1)
     &     /3.d0 -(19*dlnx**2*lfh*xm1)/9.d0 - (10*dlxm1*lfh*xm1)
     &     /3.d0 +4*dlnx*dlxm1*lfh*xm1 + (5*lfh**2*xm1)/6.d0 -
     &     dlnx*lfh**2*xm1 +(35*xm1**2)/27.d0 + (49*dli2a*xm1**2)
     &     /9.d0 + (5*dli3a*xm1**2)/9.d0 +(13*dli3b*xm1**2)/9.d0 - (13
     &     *dli3c*xm1**2)/9.d0 + (5*dli3d*xm1**2)/3.d0 -(13*dli3e*xm1
     &     **2)/36.d0 - (dli3f*xm1**2)/12.d0 -(73*dlnx*xm1**2)/18.d0
     &     -(29*dli2a*dlnx*xm1**2)/9.d0 +(11*dli2b*dlnx*xm1**2)/18.d0
     &     -(113*dlnx**2*xm1**2)/36.d0 + (17*dlnx**3*xm1**2)/54.d0 +(14
     &     *dlxm1*xm1**2)/9.d0 + (32*dli2a*dlxm1*xm1**2)/9.d0 -(4*dli2b
     &     *dlxm1*xm1**2)/9.d0 + (92*dlnx*dlxm1*xm1**2)/9.d0 -(14*dlnx
     &     **2*dlxm1*xm1**2)/9.d0 - (10*dlxm1**2*xm1**2)/3.d0 +(4*dlnx
     &     *dlxm1**2*xm1**2)/3.d0 - (7*lfh*xm1**2)/9.d0 -(16*dli2a
     &     *lfh*xm1**2)/9.d0 

      deltabbabbh = deltabbabbh + (2*dli2b*lfh*xm1**2)/9.d0 
     &     - (46*dlnx*lfh*xm1
     &     **2)/9.d0 + (7*dlnx**2*lfh*xm1**2)/9.d0 + (10*dlxm1*lfh
     &     *xm1**2)/3 .d0 - (4*dlnx*dlxm1*lfh*xm1**2)/3.d0 - (5
     &     *lfh**2*xm1**2)/6.d0 + (dlnx*lfh**2*xm1**2)/3.d0 - (205
     &     *xm1**3)/81.d0 -(11*dli2a*xm1**3)/9.d0 - (2*dli3d*xm1**3)
     &     /9.d0 - (10*dlnx*xm1**3)/27.d0 + (17*dlnx **2*xm1**3)/18.d0 +
     &     (44*dlxm1*xm1**3)/27.d0 - (8*dlnx*dlxm1*xm1 **3)/3.d0 + (8
     &     *dlxm1**2*xm1**3)/9.d0 - (22*lfh*xm1**3)/27.d0 + (4 *dlnx
     &     *lfh*xm1**3)/3.d0 - (8*dlxm1*lfh*xm1**3)/9.d0 + (2
     &     *lfh**2 *xm1**3)/9.d0 - (8*dlnx*z2)/3.d0 - (10*xm1*z2)
     &     /3.d0 + 4*dlnx*xm1*z2 + (10*xm1**2*z2)/3.d0 - (4*dlnx*xm1**2
     &     *z2)/3.d0 - (8*xm1**3*z2)/9 .d0


      end

c-}}}
c-{{{ function deltabbfbbh:

      real*8 function deltabbfbbh(xx)

      implicit none
      real*8 xx

      deltabbfbbh = 0.d0

      end

c-}}}
c-{{{ function deltabqabbh:

      real*8 function deltabqabbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltabqabbh = (-4*dli3a)/3.d0 
     &     + (4*dli3d)/3.d0 - (151*dlnx)/108.d0 -
     &     (29*dlnx**2)/36.d0 +dlnx**3/18.d0 + (7*dlnx*dlxm1)/3.d0 - (4
     &     *dlnx**2*dlxm1)/3.d0 +(4*dlnx*dlxm1**2)/3.d0 +dli2a*(2
     &     .3333333333333335d0 -2*dlnx + (8*dlxm1)/3.d0 - (4*lfh)
     &     /3.d0)+((-7*dlnx)/6.d0 + (2*dlnx**2)/3.d0 - (4*dlnx*dlxm1)
     &     /3.d0)*lfh +(dlnx*lfh**2)/3.d0 + xm1**3*(
     &     -1.0848765432098766d0 - (4*dli2a)/9.d0 - (5*dlnx)/27.d0 +(5
     &     *dlnx**2)/9.d0 + (22*dlxm1)/27.d0 - (4*dlnx*dlxm1)/3.d0 +(4
     &     *dlxm1**2)/9.d0 +(-0.4074074074074074d0 + (2*dlnx)/3.d0 -(4
     &     *dlxm1)/9.d0)*lfh +lfh**2/9.d0 - (4*z2)/9.d0) - (4*dlnx
     &     *z2)/3.d0 +xm1**2*(1.0787037037037037d0 - (2*dli3a)/3.d0 + (2
     &     *dli3d)/3.d0 -(17*dlnx)/9.d0- (15*dlnx**2)/8.d0 + dlnx**3
     &     /36.d0 + dlxm1/3.d0+(16*dlnx*dlxm1)/3.d0 - (2*dlnx**2*dlxm1)
     &     /3.d0 - (5*dlxm1**2)/3.d0 +(2*dlnx*dlxm1**2)/3.d0 +dli2a
     &     *(2.5d0 - dlnx + (4*dlxm1)/3.d0 - (2*lfh)/3.d0) +(
     &     -0.16666666666666666d0 - (8*dlnx)/3.d0 + dlnx**2/3.d0 + (5
     &     *dlxm1)/3.d0-(2*dlnx*dlxm1)/3.d0)*lfh+(
     &     -0.4166666666666667d0 + dlnx/6.d0)*lfh**2 + (5*z2)/3.d0
     &     -(2*dlnx*z2)/3.d0) 

      deltabqabbh = deltabqabbh + xm1* (-0.3148148148148148d0 
     &     + 2*dli3a - 2*dli3d 
     &     + (125*dlnx)/36.d0 + (17*dlnx**2)/8.d0 - dlnx**3/12.d0
     &     - dlxm1 - (19*dlnx*dlxm1)/3.d0 + 2*dlnx**2*dlxm1 + (5*dlxm1
     &     **2)/3 .d0 - 2*dlnx*dlxm1**2 + (0.5d0 + (19*dlnx)/6.d0 - dlnx
     &     **2 - (5 *dlxm1)/3.d0 + 2*dlnx*dlxm1)* lfh + (0
     &     .4166666666666667d0 - dlnx /2.d0)*lfh**2 + dli2a*(-3.5d0 +
     &     3*dlnx - 4*dlxm1 + 2*lfh) - (5*z2)/3.d0 + 2*dlnx*z2)

      end

c-}}}
c-{{{ function deltabqfbbh:

      real*8 function deltabqfbbh(xx)

      implicit none
      real*8 xx

      deltabqfbbh = 0.d0

      end

c-}}}
c-{{{ function deltaqqbabbh:

      real*8 function deltaqqbabbh(xx)

      implicit none
      real*8 xm1,xx,dlnx,dlxm1
      complex*16 sudilog,trilog
      complex*16 dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,dli3e,dli3f
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

      deltaqqbabbh = (-4*dli2a)/9.d0 + (2*dli2b)/9.d0 + (4*dlnx)/9.d0 +
     &     dlnx**2/9.d0 +(0.4444444444444444d0 + (4*dli2a)/3.d0 - (2
     &     *dli2b)/3.d0 - (4*dlnx)/3.d0 -dlnx**2/3.d0)*xm1 + (-1
     &     .1111111111111112d0 - (4*dli2a)/3.d0 +(2*dli2b)/3.d0 + (11
     &     *dlnx)/9.d0 + dlnx**2/3.d0)*xm1**2 +(0.6666666666666666d0 +
     &     (4*dli2a)/9.d0 - (2*dli2b)/9.d0 - dlnx/3.d0-dlnx**2/9.d0)*xm1
     &     **3

      end

c-}}}
c-{{{ function deltaqqbfbbh:

      real*8 function deltaqqbfbbh(xx)

      implicit none
      real*8 xx

      deltaqqbfbbh = 0.d0

      end

c-}}}

