C-{{{ function sgg1x0(xt)

      real*8 function sgg1x0(xt)
C..
C..  sgg1x0(xt) = sgg1(xt) - sgg1(xt->0)
C..   (see R. Harlander et al, arxiv:0912.2104)
C..
      implicit none
      real*8 tau,xm1,xt,xeps,b1gg,b1ggc,dcoefs1,sgg1,c1ggx0
      real*8 c1
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      tau = 4*(mt/mh)**2
      xeps = 0.d0

      c1 = c1ggx0( (mt/mh)**2 )
      c1 = c1 - 2*lfh

      sgg1x0 = sgg1(xt)
     &     - xm1**(nsoftmt(1,1)+1)*(sgg1(xeps)
     &     + dcoefs1(nmtlim(1,1),1.d0,0.d0)
     &     - ca*c1*ataut2)

      return
      end

C-}}}
C-{{{ function sqg1x0(yy)

      real*8 function sqg1x0(xt)
C..
C..   qg contribution at NLO, exact.
C..   
      implicit none
      real*8 tau,xt,xm1,xeps,sqg1,b1qg,b1qgc,c1qgx0
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'

      xm1 = 1.d0-xt
      tau = 4*(mt/mh)**2

      xeps = 0.d0

c$$$      sqg1x0 = sqg1(xt) - xm1**(nsoft1qgmt+1)*(sqg1(xeps) - ( 3*c1qgx0(
c$$$     &     (mt/mh)**2 ) )*ataut2)

      b1qgc = 3*c1qgx0( (mt/mh)**2 )

      sqg1x0 = sqg1(xt) - xm1**(nsoftmt(1,2)+1)*(sqg1(xeps) - ( b1qgc
     &     - cf*lfh )*ataut2)

      return
      end

C-}}}
C-{{{ function sqqb1x0(yy)

      real*8 function sqqb1x0(xt)
C..
C..   q-qbar contribution at NLO, exact.
C..   
      implicit none
      real*8 xm1,xt,xeps,sqqb1
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'

      xm1 = 1.d0-xt
      xeps = 0.d0
c$$$      xeps = xt

      sqqb1x0 = sqqb1(xt) - xm1**(nsoftmt(1,3)+1)*(sqqb1(xeps))

      return
      end

C-}}}
C-{{{ function sgg2x0balbe(yy)

      real*8 function sgg2x0balbe(xin)
C
C
      implicit real*8 (a-h,o-z)
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'
      common/tmp/balbe

      xt = xin
      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      sgg2x0balbe = sgg2exp(xt) + sgg2mtexp(xt)
     &     - 9*c2ggx0( (mt/mh)**2 )*( sub + dlog(xt) )
     &     + (balbe-1.d0)*xm1**(nsoft(2,0)+1)*sgg2exp(0.d0)
      modified = 1

      return
      end

C-}}}
C-{{{ function sgg2x0(yy)

      real*8 function sgg2x0(xtin)
C
C
      implicit none
      integer i
      real*8 tau,xt,xtin,xm1,sub,den,xeps,sgg2,c1ggx0,c2ggx0,a2gg,b1gg
     &     ,b2gg,b1ggc,a2ggc
      real*8 dcoefs2,beta0,c1,c2
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

      xt=xtin
      xm1 = 1.d0 - xt
      tau = 4*(mt/mh)**2
      xeps = 0.d0
      beta0 = (33 - 2*nf)/12.d0      
      
c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

c$$$      sgg2x0 = sgg2(xt)
c$$$     &     - (   9*c2ggx0( (mt/mh)**2 ) )*ataut2*( sub + dlog(xt) )

      c1 = c1ggx0( (mt/mh)**2 )
      c2 = c2ggx0( (mt/mh)**2 ) - 2*c1*lfh + 2*lfh**2
      
c.. this is equation (24) of arxiv:0912.2104:      
      sgg2x0 = sgg2(xt)
     &     - ca**2*c2*ataut2*( sub + dlog(xt) )
     &     - xm1**(nsoft(2,0)+1)*(sgg2(xeps)
     &     + dcoefs2(nmtlim(2,1),1.d0,0.d0,0.d0,0.d0)
     &     - ( b2gg(nmtlim(2,1),tau) - 3*beta0*lfr*(b1ggc - 2*ca*lfh) )
     &     *ataut2)
      
      return
      end

C-}}}
C-{{{ function sqg2x0(yy)

      real*8 function sqg2x0(xt)
C
      implicit none
      integer i
      real*8 tau,xt,xm1,den,sub,c2qgx0,sqg2,a2qg,b1gg,b1qg,b2qg,xeps
      real*8 b1ggc,a2qgc,b1qgc,c1ggx0,c1qgx0
      real*8 beta0
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0 - xt
      tau = 4*(mt/mh)**2
      xeps = 0.d0
      beta0 = (33 - 2*nf)/12.d0      

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      b1ggc = 3*c1ggx0( (mt/mh)**2 )
      b1qgc = 3*c1qgx0( (mt/mh)**2 )
      a2qgc = 9*c2qgx0( (mt/mh)**2 )

c$$$      sqg2x0 = sqg2(xt)
c$$$     &     - (   9*c2qgx0( (mt/mh)**2 ) )*ataut2*( sub + dlog(xt) )
      sqg2x0 = sqg2(xt)
     &     - (   a2qgc - 3/2.d0*cf*b1ggc
     &     *lfh + 3/2.d0*ca*cf*lfh**2 )*ataut2*( sub + dlog(xt) )
     &     - xm1**(nsoft(2,0)+1)*(sqg2(xeps)
     &     - ( b2qg(tau) - 3*beta0*lfr*(b1qgc - cf*lfh) )*ataut2)

      return
      end

C-}}}
C-{{{ function sqqb2x0(yy)

      real*8 function sqqb2x0(xt)
C
      implicit none
      integer i
      real*8 xt,xm1,tau,den,sub,c2qqx0,sqqb2,a2qq,a2qqc,b1gg,b1ggc
     &     ,c1ggx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      xm1 = 1.d0 - xt
      tau = 4*(mt/mh)**2

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

c$$$      sqqb2x0 = sqqb2(xt)
c$$$     &     - (9*c2qqx0( (mt/mh)**2 ) )*ataut2*( sub + dlog(xt) )

      b1ggc = 3*c1ggx0( (mt/mh)**2 )
      a2qqc = 9*c2qqx0( (mt/mh)**2 )

      sqqb2x0 = sqqb2(xt)
     &     - ( a2qqc - cf**2/ca*b1ggc*lfh + cf**2*lfh**2 )
     &     *ataut2*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function sqq2x0(yy)

      real*8 function sqq2x0(xt)
C
      implicit none
      integer i
      real*8 xt,xm1,tau,den,sub,c2qqx0,sqq2,a2qq,a2qqc,b1gg,b1ggc,c1ggx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      xm1 = 1.d0 - xt
      tau = 4*(mt/mh)**2

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

c$$$      sqq2x0 = sqq2(xt)
c$$$     &     - 9*ataut2*c2qqx0( (mt/mh)**2 )*( sub + dlog(xt) )

      b1ggc = 3*c1ggx0( (mt/mh)**2 )
      a2qqc = 9*c2qqx0( (mt/mh)**2 )

      sqq2x0 = sqq2(xt)
     &     - ( a2qqc - cf**2/ca*b1ggc*lfh + cf**2*lfh**2 )
     &     *ataut2*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function squ2x0(yy)

      real*8 function squ2x0(xt)
C
      implicit none
      integer i
      real*8 tau,xt,xm1,den,sub,c2qqx0,squ2,a2qq,a2qqc,b1gg,b1ggc,c1ggx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'

      xm1 = 1.d0 - xt
      tau = 4*(mt/mh)**2

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(2,0)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

c$$$      squ2x0 = squ2(xt)
c$$$     &     - 9*ataut2*c2qqx0( (mt/mh)**2 )*( sub + dlog(xt) )

      b1ggc = 3*c1ggx0( (mt/mh)**2 )
      a2qqc = 9*c2qqx0( (mt/mh)**2 )

      squ2x0 = squ2(xt)
     &     - ( a2qqc - cf**2/ca*b1ggc*lfh + cf**2*lfh**2 )
     &     *ataut2*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function sgg3x0(yy)

      real*8 function sgg3x0(xin)
C
C
      implicit none
      integer i
      real*8 xin,xt,sub,den,xm1,sgg3,c1ggx0,c2ggx0,c3ggx0
      real*8 c1,c2,c3
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'
      include '../commons/common-errors.f'

      xt = xin
      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(3,1)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      c1 = c1ggx0( (mt/mh)**2 )
      c2 = c2ggx0( (mt/mh)**2 )
      c3 = c3ggx0( (mt/mh)**2 ) - 2*c2*lfh + 2*c1*lfh**2 - 4/3.d0*lfh**3
      
      sgg3x0 = sgg3(xt) + ca**3*c3*ataut2*( sub + dlog(xt) )**2

      return
      end

C-}}}
C-{{{ function sqg3x0(yy)

      real*8 function sqg3x0(xt)
C
      implicit none
      integer i
      real*8 xin,xt,sub,den,xm1,sqg3exp,sqg3mtexp,c3qgx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(3,2)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      sqg3x0 = sqg3exp(xt) + sqg3mtexp(xt)
     &     - (   9*c3qgx0( (mt/mh)**2 ) )*ataut2*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function sqqb3x0(yy)

      real*8 function sqqb3x0(xt)
C
      implicit none
      integer i
      real*8 xin,xt,sub,den,xm1,sqqb3exp,sqqb3mtexp,c3qqx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(3,3)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      sqqb3x0 = sqqb3exp(xt) + sqqb3mtexp(xt)
     &     - (9*c3qqx0( (mt/mh)**2 ) )*ataut2*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function sqq3x0(yy)

      real*8 function sqq3x0(xt)
C
      implicit none
      integer i
      real*8 xin,xt,sub,den,xm1,sqq3exp,sqq3mtexp,c3qqx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(3,4)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      sqq3x0 = sqq3exp(xt) + sqq3mtexp(xt)
     &     - 9*ataut2*c3qqx0( (mt/mh)**2 )*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function squ3x0(yy)

      real*8 function squ3x0(xt)
C
      implicit none
      integer i
      real*8 xin,xt,sub,den,xm1,squ3exp,squ3mtexp,c3qqx0
      include '../commons/common-expand.f'
      include '../commons/common-vars.f'

      xm1 = 1.d0 - xt

c..   sub = -ln(xt) = -ln(1-xm1) = xm1 + xm1^2/2 + xm1^3/3 + ...
      sub = 0.d0
      do i=1,nsoft(3,5)
         den = 1.d0*i
         sub = sub + xm1**i/den
      enddo

      squ3x0 = squ3exp(xt) + squ3mtexp(xt)
     &     - 9*ataut2*c3qqx0( (mt/mh)**2 )*( sub + dlog(xt) )

      return
      end

C-}}}
C-{{{ function c1ggx0:

      real*8 function c1ggx0(yt)
c..
c..   fit function that produces the numbers for C^(1) of Table 1
c..   in Marzani et al, arxiv:0801.2544
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

c     result based on Table 1 in Marzani et al, arxiv:0801.2544
c$$$      c1ggx0 = -194.43579253947573d0 - 8.125759764115632d0/yt +
c$$$     &     59.03542270053724d0/dsqrt(yt) +  355.5075911094033d0
c$$$     &     *dsqrt(yt) - 374.98404105792184d0*yt +  231.08927255626415d0
c$$$     &     *yt*dsqrt(yt) - 68.99988288612144d0*yt**2 +
c$$$     &     4.370202819588193d0*yt**3 - 0.24775707823408039d0*yt**4

c     result based on Table 5.1 in Marzani's PhD thesis
      c1ggx0 = -259.714805697563d0 -13.6306626897191d0/yt +
     &     91.4647502300768d0/dsqrt(yt) +  384.855847213461d0
     &     *dsqrt(yt) - 297.774741403335d0*yt + 108.902043512131d0
     &     *yt*sqrt(yt) - 8.55355565126968d0*yt**2 -
     &     2.67596027512186d0*yt**3 + 0.337145313477596d0*yt**4

      return
      end

C-}}}
C-{{{ function c1qgx0:

      real*8 function c1qgx0(yt)
c..
c..   fit function that produces the limit x->0 for the qg channel;
c..   constructed from the data sent my S. Marzani,
c..   see ~/math/gghmt/ball-etal.m
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

c..   NOTE: that fit is not too great...
C$$$      c1qgx0 = -1.4271402145785763d0 - 0.2331135871187728d0/yt +
C$$$     &     0.64968570960198d0/dsqrt(yt) + 2.3826032537358794d0*dsqrt(yt)
C$$$     &     - 0.6927746610892669d0*yt + 0.03384028847048702d0*yt**2

      c1qgx0 = -0.1270955178104202d0/yt + 0.8354121865107461d0*dsqrt(yt)
     &     +0.05848253416462733d0*yt - 0.06286376780909438d0*yt**2
     &     +0.009083141325159723d0*yt**3 - 0.030418362990246223d0
     &     *dlog(yt)

      return
      end

C-}}}
C-{{{ function c2ggx0:

      real*8 function c2ggx0(yt)
c..
c..   fit function that produces the numbers for C^(2) of Table 1
c..   in Marzani et al, arxiv:0801.2544
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

c     result based on Table 1 in Marzani et al, arxiv:0801.2544
c$$$      c2ggx0 = 250.41843107662382d0 + 11.427544031360895d0/yt -
c$$$     &     79.67241527262105d0/dsqrt(yt) -  450.40406507303146d0
c$$$     &     *dsqrt(yt) + 500.54694144679763d0*yt -  311.2604723909205d0
c$$$     &     *yt*dsqrt(yt) + 93.18547848041179d0*yt**2 -
c$$$     &     5.908094620578248d0*yt**3 + 0.334798688860235d0*yt**4

c     result based on Table 5.1 in Marzani's PhD thesis
      c2ggx0 = 1151.20911427526d0 + 42.3120547985503d0/yt -
     &     336.299203232444d0/dsqrt(yt) - 2179.27887828242d0
     &     *sqrt(yt) + 2448.77933591726d0*yt - 1572.9396792799d0
     &     *yt*dsqrt(yt) + 485.263143689644d0*yt**2 -
     &     32.2679690234154d0*yt**3 + 1.89100722479237d0*yt**4

      return
      end

C-}}}
C-{{{ function c2gqx0:

      real*8 function c2qgx0(yt)
c..
c..   fit function that produces the numbers for C^(2) of Table 1
c..   in Harlander, ..., Marzani et al, arxiv:???????
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

      c2qgx0 = 0.5247689919369154d0/yt + 1.5156042013259359d0
     &     *dsqrt(yt) +0.6023061837965867d0*yt - 0.03436330601053507d0
     &     *yt**2-0.001949391979200509d0*yt**3 + 1.5159876785925033d0
     &     *dlog(yt)
      
      return
      end

C-}}}
C-{{{ function c2qqx0:

      real*8 function c2qqx0(yt)
c..
c..   fit function that produces the numbers for C^(2) of Table 1
c..   in Harlander, ..., Marzani et al, arxiv:???????
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

      c2qqx0 = 0.17089821452225454d0/yt + 0.2169100856420175d0
     &     *dsqrt(yt) +0.22740881402993324d0*yt - 0.009837319004985617d0
     &     *yt**2 -0.0009404078304873317d0*yt**3 + 0.5292193557994358d0
     &     *dlog(yt)

      return
      end

C-}}}
C-{{{ function c3ggx0:

      real*8 function c3ggx0(yt)
c..
c..   fit function that produces the numbers for C^(3) of Table 5.1
c..   in Marzani's PhD thesis
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

      c3ggx0 = -1711.30435539634d0 - 73.874077873869d0/yt +
     &     545.034273959081d0/dsqrt(yt) + 2906.91452861361d0*dsqrt(yt) -
     &     2860.15928059394d0*yt + 1626.32810730606d0*yt*dsqrt(yt) -
     &     437.112650546066d0*yt**2 + 20.6085309678596d0*yt**3 -
     &     0.734780850785562d0*yt**4

      return
      end

C-}}}
C-{{{ function c3qgux0:

      real*8 function c3qgx0(yt)
c..
c..   fit function that produces the numbers for C^(2) of Table 1
c..   in Harlander, ..., Marzani et al, arxiv:???????
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

      call printdieggh('function c3qgx0 not functional')

      c3qgx0 = -1000000000.d0
      
      return
      end

C-}}}
C-{{{ function c3qqx0:

      real*8 function c3qqx0(yt)
c..
c..   fit function that produces the numbers for C^(2) of Table 1
c..   in Harlander, ..., Marzani et al, arxiv:???????
c..
      implicit none
      real*8 yt
      include '../commons/common-vars.f'

      call printdieggh('function c3qqx0 not functional')
      c3qqx0 = -1000000000.d0

      return
      end

C-}}}

C-{{{ function b1gg:

      real*8 function b1gg(tauin)
c..
c..   b1gg( 4*mt^2/mh^2 ) = 3*c1ggx0( mt^2/mh^2 ),
c..   but rather than an interpolating function, it uses
c..   linear pointwise interpolation.
c..
      implicit real*8 (a-z)
      integer ii
      real*8 tau(23),b1ggtab(23)

      data tau/1.0,1.5,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0,
     &     5.5d0,6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0,
     &     9.5d0,10.0d0,10.5d0,11.0d0,11.5d0,12.0d0/

      data b1ggtab/-0.8821d0,2.9212d0,5.0234d0,6.5538d0,7.7650d0,
     &     8.7693d0, 9.6279d0,10.3781d0,11.0444d0,11.6437d0,12.1883d0
     &     ,12.6875d0,13.1482d0,13.5760d0,13.9752d0,14.3495d0,14.7018d0
     &     ,15.0345d0,15.3497d0,15.6491d0,15.9343d0,16.2065d0,16.4670d0/

      if ((tauin.lt.tau(1)).or.(tauin.gt.tau(23))) then
         write(6,*) 'b1gg: input out of range'
         stop
      endif

      ii = 1
      do while (tauin.gt.tau(ii))
         ii = ii+1
      enddo

      tval = (tauin - tau(ii-1))/(tau(ii)-tau(ii-1))
      dum = tau(ii-1)+tval*tau(ii)
      b1gg = (1.d0-tval)*b1ggtab(ii-1) + tval*b1ggtab(ii)
      
      return
      end

C-}}}
C-{{{ function b1qg:

      real*8 function b1qg(tauin)
c..
c..   b1qg( 4*mt^2/mh^2 ) = 3*c1qgx0( mt^2/mh^2 ),
c..   but rather than an interpolating function, it uses
c..   linear pointwise interpolation.
c..
      implicit real*8 (a-z)
      integer ii
      real*8 tau(23),b1qgtab(23)

      data tau/1.0,1.5,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0,
     &     5.5d0,6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0,
     &     9.5d0,10.0d0,10.5d0,11.0d0,11.5d0,12.0d0/

      data b1qgtab/-0.1960d0,0.6492d0,1.1163d0,1.4564d0,1.7255d0,1.9487,
     &     2.1395d0,2.3062d0,2.4543d0,2.5875d0,2.7085d0,2.8194d0,2.9218,
     &     3.0169d0,3.1056d0,3.1888d0,3.2671d0,3.3410d0,3.4110d0,3.4776,
     &     3.5410d0,3.6015d0,3.6593d0/

      if ((tauin.lt.tau(1)).or.(tauin.gt.tau(23))) then
         write(6,*) 'b1qg: input out of range'
         stop
      endif

      ii = 1
      do while (tauin.gt.tau(ii))
         ii = ii+1
      enddo

      tval = (tauin - tau(ii-1))/(tau(ii)-tau(ii-1))
      dum = tau(ii-1)+tval*tau(ii)
      b1qg = (1.d0-tval)*b1qgtab(ii-1) + tval*b1qgtab(ii)
      
      return
      end

C-}}}
C-{{{ function b2gg:

      real*8 function b2gg(nmtlim,tau)
c..
c..   Not known yet. We set it to the value chosen (implicitely) in
c..   R. Harlander et al., arXiv:0912.2104.
c..
      implicit none
      real*8 tau,xeps,sgg2,dcoefs2
      integer nmtlim
      include '../commons/common-vars.f'

      xeps = 0.d0

      b2gg = (sgg2(xeps) + dcoefs2(nmtlim,1.d0,0.d0,0.d0,0.d0))/ataut2
      
      return
      end

C-}}}
C-{{{ function b2qg:

      real*8 function b2qg(tau)
c..
c..   Not known yet. We set it to the value chosen (implicitely) in
c..   R. Harlander et al., arXiv:0912.2104.
c..
      implicit none
      real*8 tau,xeps,sqg2
      include '../commons/common-vars.f'

      xeps = 0.d0

      b2qg = sqg2(xeps)/ataut2
      
      return
      end

C-}}}
C-{{{ function a2gg:

      real*8 function a2gg(tauin)
c..
c..   a2gg( 4*mt^2/mh^2 ) = 9*c2ggx0( mt^2/mh^2 ),
c..   but rather than an interpolating function, it uses
c..   linear pointwise interpolation.
c..
      implicit real*8 (a-z)
      integer ii
      real*8 tau(23),a2ggtab(23)

      data tau/1.0,1.5,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0, 5.5d0
     &     ,6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0, 9.5d0,10.0d0
     &     ,10.5d0,11.0d0,11.5d0,12.0d0/
      
      data a2ggtab/33.0465d0,  35.9907d0,  44.2884d0,  53.1336d0,
     &     61.8029d0,70.1088d0,78.0127d0,85.5245d0,  92.6698d0,
     &     99.4782d0,105.9788d0,112.1985d0,118.1616d0,123.8897d0,
     &     129.4021d0,134.7158d0,139.8461d0,144.8064d0,149.6090d0
     &     ,154.2646d0,158.7829d0,163.1728d0,167.4422d0/

      if ((tauin.lt.tau(1)).or.(tauin.gt.tau(23))) then
         write(6,*) 'a2gg: input out of range'
         stop
      endif

      ii = 1
      do while (tauin.gt.tau(ii))
         ii = ii+1
      enddo

      tval = (tauin - tau(ii-1))/(tau(ii)-tau(ii-1))
      dum = tau(ii-1)+tval*tau(ii)
      a2gg = (1.d0-tval)*a2ggtab(ii-1) + tval*a2ggtab(ii)
      
      return
      end

C-}}}
C-{{{ function a2qg:

      real*8 function a2qg(tauin)
c..
c..   a2qg( 4*mt^2/mh^2 ) = 9*c2qgx0( mt^2/mh^2 ),
c..   but rather than an interpolating function, it uses
c..   linear pointwise interpolation.
c..
      implicit real*8 (a-z)
      integer ii
      real*8 tau(23),a2qgtab(23)

      data tau/1.0,1.5,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0, 5.5d0
     &     ,6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0, 9.5d0,10.0d0
     &     ,10.5d0,11.0d0,11.5d0,12.0d0/
      
      data a2qgtab/8.7703d0,9.5484d0,12.2677d0,15.1924d0,18.0679d0
     &     ,20.8272d0,23.4553d0,25.9547d0,28.3331d0,30.6002d0,32.7654d0
     &     ,34.8374d0,36.8244d0,38.7333d0,40.5706d0,42.3420d0,44.0524d0
     &     ,45.7063d0,47.3077d0,48.8603d0,50.3672d0,51.8315d0,53.2556d0/

      if ((tauin.lt.tau(1)).or.(tauin.gt.tau(23))) then
         write(6,*) 'a2qg: input out of range'
         stop
      endif

      ii = 1
      do while (tauin.gt.tau(ii))
         ii = ii+1
      enddo

      tval = (tauin - tau(ii-1))/(tau(ii)-tau(ii-1))
      dum = tau(ii-1)+tval*tau(ii)
      a2qg = (1.d0-tval)*a2qgtab(ii-1) + tval*a2qgtab(ii)
      
      return
      end

C-}}}
C-{{{ function a2qq:

      real*8 function a2qq(tauin)
c..
c..   a2qq( 4*mt^2/mh^2 ) = 9*c2qqx0( mt^2/mh^2 ),
c..   but rather than an interpolating function, it uses
c..   linear pointwise interpolation.
c..
      implicit real*8 (a-z)
      integer ii
      real*8 tau(23),a2qqtab(23)

      data tau/1.0,1.5,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0, 5.5d0
     &     ,6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0, 9.5d0,10.0d0
     &     ,10.5d0,11.0d0,11.5d0,12.0d0/

      data a2qqtab/1.2681d0,1.3782d0,2.1563d0,3.0088d0,3.8524d0,4.6644d0
     &     ,5.4393d0,6.1771d0,6.8798d0,7.5501d0,8.1907d0,8.8039d0
     &     ,9.3922d0,9.9576d0,10.5019d0,11.0268d0,11.5337d0,12.0240d0
     &     ,12.4989d0,12.9594d0,13.4064d0,13.8408d0,14.2633d0/


      if ((tauin.lt.tau(1)).or.(tauin.gt.tau(23))) then
         write(6,*) 'a2qq: input out of range'
         stop
      endif

      ii = 1
      do while (tauin.gt.tau(ii))
         ii = ii+1
      enddo

      tval = (tauin - tau(ii-1))/(tau(ii)-tau(ii-1))
      dum = tau(ii-1)+tval*tau(ii)
      a2qq = (1.d0-tval)*a2qqtab(ii-1) + tval*a2qqtab(ii)
      
      return
      end

C-}}}
