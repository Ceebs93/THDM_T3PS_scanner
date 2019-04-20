! This file is part of SusHi.
! 
! It includes the calculation of loop integrals
! as well as logarithms and polylogarithms.
! 
#define integrals 2

#if (integrals==1)

c-{{{ with    QCD-Loop

c-{{{ 2-point-function

      function Integral2(m1,mb2)
      implicit none
      double complex Integral2
      double precision m1, mb2
      integer i
      i=0
      call ffxb0(Integral2, dlog(mb2), 1.d0, m1, mb2, mb2, i)
      end

c-}}}
c-{{{ 3-point-function

      function Integral3(m1,m2,mb2)
      implicit none
      double complex Integral3
      double precision m1,m2,mb2,m(6)
      integer i
      i=0
      m(1)=mb2
      m(2)=mb2
      m(3)=mb2
      m(4)=m1
      m(5)=0.d0
      m(6)=m2
      call ffxc0(Integral3, m, i)
      end

c-}}}
c-{{{ 4-point-function

      function Integral4(m1,m2,m3,mb2)
      implicit none
      double complex Integral4
      double precision m1,m2,m3,mb2,m(13)
      integer i
      i=0
      m(1)=mb2
      m(2)=mb2
      m(3)=mb2
      m(4)=mb2
      m(5)=0.d0
      m(6)=0.d0
      m(7)=0.d0
      m(8)=m1+m2+m3
      m(9)=m1
      m(10)=m3
      m(11)=0.d0
      m(12)=0.d0
      m(13)=0.d0
      call ffxd0(Integral4, m, i)
      end

c-}}}

c-}}}

#elif (integrals>=2)

c-{{{ without QCD-Loop

c-{{{ 2-point-function

      function Integral2(ma,mq2)
      implicit none
      double complex Integral2
      double precision ma,m1,mq2,re,im,Pi
      parameter (Pi=3.141592653589793238462643d0)
      m1=ma/mq2
      if(m1.eq.0.d0) then
         re = -2.d0
         im = 0.d0
      elseif(m1.lt.0.d0) then
         re = -2.d0*dsqrt(1.d0-4.d0/m1)*
     &        dlog((dsqrt(-m1)+dsqrt(-m1+4.d0))/2.d0)
         im = 0.d0
      elseif(m1.lt.4.d0) then
         re = -2.d0*dsqrt(4.d0/m1-1.d0)*dasin(dsqrt(m1)/2.d0)
         im = 0.d0
      elseif(m1.gt.4.d0) then
         re = -2.d0*dsqrt(1-4.d0/m1)*(dlog((dsqrt(m1)+dsqrt(m1-4.d0))/
     &(2.d0)))
         im = 2.d0*dsqrt(1.d0-4.d0/m1)*Pi/2.d0
      else 
         write(*,*) 'SusHi error in Integral2'
         stop
      endif
      Integral2 = (1.d0,0.d0)*re + (0.d0,1.d0)*im
      end

c-}}}
c-{{{ 3-point-function

      function Integral3(ma,mb,mq2)
      implicit none
      double complex Integral3
      double precision ma,mb,mq2,m1,m2,re,im,Pi
      parameter (Pi=3.141592653589793238462643D0)
      m1=ma/mq2
      m2=mb/mq2

      if(m1.eq.0.d0) then
         re=0.d0
         im=0.d0
      elseif(m1.lt.0.d0) then
         re=(dlog(2.d0/(2.d0-m1+dsqrt(m1*(-4+m1))))**2)/4.d0
         im=0.d0
      elseif(m1.lt.4.d0) then
         re=-dasin(dsqrt(m1/4.d0))**2
         im=0.d0
      elseif(m1.gt.4.d0) then
         re=(dlog((-2.d0+m1+dsqrt(m1*(-4.d0+m1)))/(2.d0))**2-Pi**2)/4.d0
         im=-Pi*dlog((-2.d0+m1+dsqrt(m1*(-4.d0+m1)))/2.d0)/(m1-m2)
      endif

      if(m2.eq.0.d0) then
      elseif(m2.lt.0.d0) then
         re=re-(dlog(2.d0/(2.d0-m2+dsqrt(m2*(-4+m2))))**2)/4.d0
      elseif(m2.le.4.d0) then
         re=re+dasin(dsqrt(m2/4.d0))**2
      elseif(m2.gt.4.d0) then
        re=re-(dlog((-2.d0+m2+dsqrt(m2*(-4.d0+m2)))/2.d0)**2-Pi**2)/4.d0
        im=im+Pi*dlog((-2.d0+m2+dsqrt(m2*(-4.d0+m2)))/(2.d0))/(m1-m2)
      else 
         write(*,*) 'SusHi error in Integral3'
         stop
      endif
      Integral3=((1.d0,0.d0)*2.d0/(m1-m2)*re+(0.d0,1.d0)*im)/mq2
      end

c-}}}
c-{{{ 4-point-function

      function Integral4(ma,mb,mc,mq2)
      implicit none
      double complex Integral4,wrz,p(1:3)
      double precision ma,mb,mc,mq2,m1,m2,m3,x(1:2),pn(1:3),tmp,
     &dli2,arg,zeta2,zeta3,pi
      integer i,k
      double complex cli2
      COMMON/ZETACONST/ZETA2,ZETA3
      parameter (Pi=3.141592653589793238462643D0)

      m1=ma/mq2
      m2=mb/mq2
      m3=mc/mq2
      pn(1)=m1
      pn(2)=m1+m2+m3
      pn(3)=m3

      tmp=1+m1*m3/m2/2.d0
      if(tmp**2.eq.1.d0) then
         Integral4=1/6.d0/mq2**2
         return
      endif

      x(2)=-tmp-dsqrt(tmp**2-1)
      x(1)=1/x(2)

      Integral4=(2.d0,0.d0)*(2*dli2(1+x(2))+dlog(-x(2))**2/2)
      do i=1,3
         wrz=cdsqrt((1.d0,0.d0)*(1-pn(i)/2.d0)**2-1)
         if(pn(i).eq.0.d0) then
            p(i) = 1.d0
         elseif(pn(i).gt.0.d0) then
            p(i)=1/(1-pn(i)/2.d0-wrz)
         else
            p(i)=1-pn(i)/2.d0+wrz
         endif

         do k=1,2
            if(pn(i).le.0.d0) then
               arg = real(p(i)*x(k))
               Integral4=Integral4+(-1)**(i+k)*
     &              (-2*dli2(1+1/arg)-dlog(-arg)**2/2)
            elseif(pn(i).ge.4.d0) then
               arg = real(p(i)*x(k))
               Integral4=Integral4+(-1)**(i+k) * (
     &              dli2(-1/arg)-zeta2-dlog(1+1/arg)*dlog(arg)
     &              -(dlog(arg)-(0.d0,1.d0)*Pi)**2/4.d0
     &              -(0.d0,1.d0)*Pi*dlog(1+arg))*2.d0
            else
               Integral4=Integral4+(-1)**(i+k)*(cli2(1+p(i)*x(k))
     &              +cli2(1+x(k)/p(i)))
            endif
         enddo
      enddo
      Integral4=Integral4/dsqrt(tmp**2-1)/m2/mq2**2/2.d0
      end

c-}}}

c-}}}

#endif

c-{{{ complex dilogarithm (from HIGLU)

      COMPLEX*16 FUNCTION CLI2(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,SUSHI_CLI2
      COMMON/ZETACONST/ZETA2,ZETA3
      CLI2 = 0.D0
      ZERO=1.D-40
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      IF(R2.LE.ZERO)THEN
        CLI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          CLI2=DCMPLX(ZETA2)
        ELSE
          CLI2=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
       CLI2=SUSHI_CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)+0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        CLI2=-SUSHI_CLI2(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        CLI2=-SUSHI_CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        CLI2=SUSHI_CLI2(Y)
        RETURN
      ENDIF
      END

      COMPLEX*16 FUNCTION SUSHI_CLI2(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      SUSHI_CLI2=B2(NBER)
      DO 111 I=N,1,-1
        SUSHI_CLI2=Z*SUSHI_CLI2+B2(I)
111   CONTINUE
      SUSHI_CLI2=Z**2*SUSHI_CLI2+Z
      RETURN
      END

      SUBROUTINE SUSHI_BERNINI(N)
C--INITIALIZATION OF COEFFICIENTS FOR POLYLOGARITHMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(18),PB(19)
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/ZETACONST/ZETA2,ZETA3
      COMMON/POLY/NBER
 
      NBER=N
      PI=4.D0*DATAN(1.D0)
 
      B(1)=-1.D0/2.D0
      B(2)=1.D0/6.D0
      B(3)=0.D0
      B(4)=-1.D0/30.D0
      B(5)=0.D0
      B(6)=1.D0/42.D0
      B(7)=0.D0
      B(8)=-1.D0/30.D0
      B(9)=0.D0
      B(10)=5.D0/66.D0
      B(11)=0.D0
      B(12)=-691.D0/2730.D0
      B(13)=0.D0
      B(14)=7.D0/6.D0
      B(15)=0.D0
      B(16)=-3617.D0/510.D0
      B(17)=0.D0
      B(18)=43867.D0/798.D0
      ZETA2=PI**2/6.D0
      ZETA3=1.202056903159594D0
 
      DO 995 I=1,18
        B2(I)=B(I)/SUSHI_FACULT(I+1)
        B12(I)=DFLOAT(I+1)/SUSHI_FACULT(I+2)*B(I)/2.D0
        PB(I+1)=B(I)
        B3(I)=0.D0
995   CONTINUE
      PB(1)=1.D0
      DO 996 I=1,18
      DO 996 J=0,I
       B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/SUSHI_FACULT(I-J)/SUSHI_FACULT(J+1)
     .                                            /DFLOAT(I+1)
996   CONTINUE
 
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SUSHI_FACULT(N)
C--DOUBLE PRECISION VERSION OF FACULTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SUSHI_FACULT=1.D0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        SUSHI_FACULT=SUSHI_FACULT*DFLOAT(I)
999   CONTINUE
      RETURN
      END

c-}}}
C-{{{ real dilogarithm (from HqT)

        FUNCTION dli2(x)
        implicit none
*      !! Dilogarithm for arguments x < = 1.0
        real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
        real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD
        real*8 dli2
        integer  i

        DATA ZERO /0.0d0/, ONE /1.0d0/
        DATA HALF /0.5d0/, MALF /-0.5d0/
        DATA MONE /-1.0d0/, MTWO /-2.0d0/
        DATA PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/

        DATA C( 0) / 0.4299669356081370d0/                              
        DATA C( 1) / 0.4097598753307711d0/                              
        DATA C( 2) /-0.0185884366501460d0/                              
        DATA C( 3) / 0.0014575108406227d0/                              
        DATA C( 4) /-0.0001430418444234d0/                              
        DATA C( 5) / 0.0000158841554188d0/                              
        DATA C( 6) /-0.0000019078495939d0/                              
        DATA C( 7) / 0.0000002419518085d0/                              
        DATA C( 8) /-0.0000000319334127d0/                              
        DATA C( 9) / 0.0000000043454506d0/                              
        DATA C(10) /-0.0000000006057848d0/                              
        DATA C(11) / 0.0000000000861210d0/                              
        DATA C(12) /-0.0000000000124433d0/                              
        DATA C(13) / 0.0000000000018226d0/                              
        DATA C(14) /-0.0000000000002701d0/                              
        DATA C(15) / 0.0000000000000404d0/                              
        DATA C(16) /-0.0000000000000061d0/                              
        DATA C(17) / 0.0000000000000009d0/                              
        DATA C(18) /-0.0000000000000001d0/                              
                                                                           
        if(x .gt. 1.00000000001d0) then                                    
          write(6,*)'problems in LI2'
          write(6,*)'x=',x 
          stop                                               
        elseif(x .gt. 1.0d0) then                                          
          x = 1.d0                                                      
        endif                                                              
        IF(X .EQ. ONE) THEN                                                
         LI2OLD=PI6
         dli2=LI2OLD                                                       
         RETURN                                                            
        ELSE IF(X .EQ. MONE) THEN                                          
         LI2OLD=MALF*PI6
         dli2=LI2OLD                                                  
         RETURN                                                            
        END IF                                                             
        T=-X                                                               
        IF(T .LE. MTWO) THEN                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
        ELSE IF(T .LT. MONE) THEN                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
        ELSE IF(T .LE. MALF) THEN                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
        ELSE IF(T .LT. ZERO) THEN                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
        ELSE IF(T .LE. ONE) THEN                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
        ELSE                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
        END IF                                                             
                                                                           
        H=Y+Y-ONE                                                          
        ALFA=H+H                                                           
        B1=ZERO                                                            
        B2=ZERO                                                            
        DO  I = 18,0,-1                                                    
          B0=C(I)+ALFA*B1-B2                                               
          B2=B1                                                            
          B1=B0                                                            
        ENDDO                                                              
        LI2OLD=-(S*(B0-H*B2)+A) 
         dli2=LI2OLD
        end     

C-}}}
