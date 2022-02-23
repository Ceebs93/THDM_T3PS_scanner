c..   NLO and NNLO coefficients of the correlator <O1O1>.
c..   LO term is normalized to 1.
c..   These results are building blocks for the construction of the 
c..   decay width of the Higgs boson.

      program deltaPS

      implicit real*8(a-z)

      mh   = 100.d0
      mur  = 240.d0
      nl   = 5.d0

      dum1 = deltaPS1(mh,mur,nl)
      dum2 = deltaPS2(mh,mur,nl)

      print*,dum1,dum2

      end


C-{{{ function deltaPS1(mh,mur,nl):

      real*8 function deltaPS1(mh,mur,nl)
c..
c..   NLO coefficient of the correlator <O1O1>
c..
      implicit real*8 (a-z)

      LQQ = dlog(mh**2/mur**2)
      
      deltaPS1 = ( 73.d0/4.d0 - 11.d0/2.d0*LQQ ) 
     &     + nl*(  - 7.d0/6.d0 + 1.d0/3.d0*LQQ )

      end

C-}}}
C-{{{ function deltaPS2(mh,mur,nl):

      real*8 function deltaPS2(mh,mur,nl)
c..
c..   NNLO coefficient of the correlator <O1O1>
c..
      implicit real*8 (a-z)

      LQQ = dlog(mh**2/mur**2)
      z2 = 1.64493406684822643647241516665d0
      z3 = 1.20205690315959428539973816151d0
      
      deltaPS2 = 37631.d0/96.d0 - 495.d0/8.d0*z3 - 363.d0/8.d0*z2 
     &     - 2817.d0/16.d0*LQQ + 363.d0/16.d0*LQQ**2
     &     + nl*(  - 7189.d0/144.d0 + 5.d0/4.d0*z3 + 11.d0/2.d0*z2 
     &     + 263.d0/12.d0*LQQ - 11.d0/4.d0*LQQ**2 )
     &     + nl**2*( 127.d0/108.d0 - 1.d0/6.d0*z2 - 7.d0/12.d0*LQQ 
     &     + 1.d0/12.d0*LQQ**2 )

      end

C-}}}


c.. result for 'hggnum.mat':
c   resG2G2 = (
c       + nl*api * (  - 7/6 + 1/3*LQQ )
c       + nl*api^2 * (  - 7189/144 + 5/4*z3 + 11/2*z2 + 263/12*LQQ - 11/4*LQQ^2
c          )
c       + nl^2*api^2 * ( 127/108 - 1/6*z2 - 7/12*LQQ + 1/12*LQQ^2 )
c       + api * ( 73/4 - 11/2*LQQ )
c       + api^2 * ( 37631/96 - 495/8*z3 - 363/8*z2 - 2817/16*LQQ + 363/16*LQQ^2
c          )
c       + 1
c);


