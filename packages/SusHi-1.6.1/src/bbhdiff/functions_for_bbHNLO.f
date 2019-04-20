com-- 
com-- function to create Momenta for 2->2 Process - begin
com-- 

C-{{{ subroutine: xtoP

      subroutine xtoP(x,PiA,Pf)
      implicit none
      double precision x(10),PiA(4,3),Pf(4,3),E3,E4,costheta,
     f sintheta,E,k5,sinchi,coschi,x1,x2
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      x1 = x(1)
      x2 = x(2)
      E = dsqrt(x1*x2*Shad)/2d0
** initial state momenta
      PiA(4,1) = E
      PiA(1,1) = 0d0
      PiA(2,1) = 0d0
      PiA(3,1) = E
      PiA(4,2) = E
      PiA(1,2) = 0d0
      PiA(2,2) = 0d0
      PiA(3,2) = -E

      costheta=x(3)
      sintheta=dsin(dacos(costheta))
      E3=(4*E**2-Hmass**2)/(4*E)
      E4=(4*E**2+Hmass**2)/(4*E)

      Pf(4,1) = E3
      Pf(1,1) = 0d0
      Pf(2,1) = E3*sintheta
      Pf(3,1) = E3*costheta

      Pf(4,2)=E4
      Pf(1,2)=-Pf(1,1)
      Pf(2,2)=-Pf(2,1)
      Pf(3,2)=-Pf(3,1)
      end

C-}}}      

com-- 
com-- end - function to create Momenta for 2->2 Process
com-- 

com-- 
com--  matrix elements for  bb->gH, bg->bH and gb->bH
com-- 

C-{{{ FUNCTION: Mtree

      DOUBLE PRECISION FUNCTION Mtree(p12,p34,cA,cF)
      implicit none
      DOUBLE PRECISION p12(4,3),p34(4,3)
      DOUBLE PRECISION s,t,u,Mhi,dot,cA,cF
      s=2.d0*dot(p12,1,p12,2)
      t=-2.d0*dot(p12,1,p34,1)
      u=-2.d0*dot(p12,2,p34,1)
      Mhi=dot(p34,2,p34,2)
      Mtree=cA*cF*4*((t+u)**2+2*Mhi*s)/(t*u)
      RETURN
      END

C-}}}
C-{{{ FUNCTION: Mtreegb

      DOUBLE PRECISION FUNCTION Mtreegb(P12,P34)
      implicit none
      DOUBLE PRECISION p12(4,3),p34(4,3)
      DOUBLE PRECISION s,t,u,Mhi,dot
      include '../commons/common-consts.f'
      s=-2.d0*dot(p34,1,p12,2)
      t=-2.d0*dot(p34,1,p12,1)
      u=2.d0*dot(p12,2,p12,1)
      Mhi=dot(p34,2,p34,2)
      Mtreegb=-cA*cF*4*((t+u)**2+2*Mhi*s)/(t*u)
      RETURN
      END

C-}}}      
C-{{{ FUNCTION: Mtreebg

      DOUBLE PRECISION FUNCTION Mtreebg(P12,P34)
      implicit none
      DOUBLE PRECISION p12(4,3),p34(4,3)
      DOUBLE PRECISION s,t,u,Mhi,dot
      include '../commons/common-consts.f'
      s=-2.d0*dot(p34,1,p12,1)
      t=-2.d0*dot(p34,1,p12,2)
      u=2.d0*dot(p12,2,p12,1)
      Mhi=dot(p34,2,p34,2)
      Mtreebg=-cA*cF*4*((t+u)**2+2*Mhi*s)/(t*u)
      RETURN
      END

C-}}}

com-- 
com-- end - matrix elements for  bb->gH, bg->bH and gb->bH
com-- 

com-- 
com-- runs Couplings and fills common-block /Couplings/ 
com-- 

C-{{{ SUBROUTINE: runCouplings

      SUBROUTINE runCouplings(muR_in,muY,orderin)
      implicit none
      DOUBLE PRECISION muY,alphamb,alphamZ,alphamuY,mbmuY,muR_in
      INTEGER orderin
      include '../commons/common-vars.f'
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-cplngs.f'
      include '../commons/common-inputoutput.f'

      gev2pb = .38937966d+9 !  conversion  1/GeV^2 -> pb

      call runalpha(apimz,mZ,muR_in,nf,orderin+1,0,alphabbh)
      alphabbh = alphabbh * Pi

      call runalpha(apimz,mZ,muY,nf,orderin+1,0,alphamuY)
      call runalpha(apimz,mZ,qqhscale,nf,orderin+1,0,alphamb)
      call runmass(qqhyuk,alphamb,alphamuY,nf,orderin+1,mbmuY)

      HbbCoupl = mbmuY**2 * Sqrt(2.d0) * GF ! /246.220569d0**2
      !HbbCoupl = 1.d0 ! /246.220569d0**2
      mbmuY = 0.d0

      !call polemass(mbmb,alphamb,nf,orderin,mbpole)

      END

C-}}}

com-- 
com-- functions to get couplings in calculation
com-- 

C-{{{ FUNCTION: alpha_s

      DOUBLE PRECISION FUNCTION alpha_s(order)
      implicit none
      include '../commons/common-bbhdiff-cplngs.f'
      INTEGER order
      alpha_s=0d0
      IF(order==0) alpha_s=alphabbh
      IF(order==1) alpha_s=alphabbh
      IF(order==2) alpha_s=alphabbh
      IF(alpha_s==0d0) THEN
      print*, 'ACHTUNG keine Ordnung gesetzt -- Programm gestoppt !!!'
      stop
      ENDIF
      END

C-}}}
C-{{{ FUNCTION: Hbb

      DOUBLE PRECISION FUNCTION Hbb(order)
      implicit none
      include '../commons/common-bbhdiff-cplngs.f'
      INTEGER order
      Hbb=0d0
      IF(order==0) Hbb=HbbCoupl
      IF(order==1) Hbb=HbbCoupl
      IF(order==2) Hbb=HbbCoupl
      IF(Hbb==0d0) THEN
      print*, 'ACHTUNG keine Ordnung oder Kopplung gesetzt',
     f'-- Programm gestoppt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop
      ENDIF
      END

C-}}}

com-- 
com-- end - functions to calculate mb(muR) and alphaS(muR)
com-- 


com-- 
com-- Fjet functions - begin
com-- 

C-{{{ FUNCTION: Fjet2

      INTEGER FUNCTION Fjet2(PfnT,pTquest)
      implicit none
      DOUBlE PRECISION PfnT(4,3),PinT(4,3),Pf(4,3),Pi(4,3),y,pTquest
     f ,pTH,yH,etaH,pTjet,yjet,CALCpT,CALCy,CALCeta
      include '../commons/common-bbhdiff-cuts.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-btag.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      Fjet2=1

      CALL LorentzTrafo3(Pi,Pf,PinT,PfnT)

!      write(*,*) "Pi",Pi
!      write(*,*) "Pf",Pf

      pTH = CALCpT(Pf,2)
      yH  = CALCy(Pf,2)
      etaH  = CALCeta(Pf,2)
       
!      write(*,*) "kinematical variables",pTH,yH,pT_Hmin
!     &,pT_Hmax,y_Hmax,y_Hmin

      IF(pTH < pT_Hmin .or. pTH >= pT_Hmax) THEN
      Fjet2=0
      return
      ENDIF

      if (pseudorapflag) then
      IF (Abs(etaH) >= y_Hmax .or. Abs(etaH) < y_Hmin) THEN
      Fjet2=0
      return
      ENDIF
      else
      IF (Abs(yH) >= y_Hmax .or. Abs(yH) < y_Hmin) THEN
      Fjet2=0
      return
      ENDIF
      end if

      IF(btag == -1) THEN
      IF(jets == -1) return

      pTjet=CALCpT(Pf,1)
      yjet=CALCy(Pf,1)
      IF(pTjet < pT_JETmin .or. pTjet >= pT_JETmax) Fjet2=0
      IF(yjet >= y_JETmax .or. yjet < y_JETmin) Fjet2=0

      ELSEIF(btag == 1 .or. btag == 3) THEN

      IF(bnr == 0) THEN
      Fjet2=0
      return
      ENDIF

      pTjet=CALCpT(Pf,1)
      yjet=CALCy(Pf,1)

      IF(pTjet < pT_bmin .or. pTjet >= pT_bmax) Fjet2=0
      IF(yjet >= y_bmax .or. yjet < y_bmin) Fjet2=0
      ELSE
      Fjet2=0
      ENDIF
      END

C-}}}      
C-{{{ FUNCTION: Fjet1

      INTEGER FUNCTION Fjet1(PfnT,pTquest)
      implicit none
      DOUBlE PRECISION PfnT(4,3),PinT(4,3),Pf(4,3),Pi(4,3),y,pTquest
     f ,CALCy,CALCeta,eta
      include '../commons/common-bbhdiff-cuts.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      Fjet1=1

      IF(jets == -1) THEN
       
      CALL LorentzTrafo3(Pi,Pf,PinT,PfnT)

      IF(dsqrt(Pf(1,2)**2+Pf(2,2)**2) < pT_Hmin) Fjet1=0
      IF(dsqrt(Pf(1,2)**2+Pf(2,2)**2) >= pT_Hmax) Fjet1=0

    !  write(*,*) "TEST",dsqrt(Pf(1,2)**2+Pf(2,2)**2),
    ! &pT_Hmin,pT_Hmax,y_Hmin,y_Hmax,Fjet1

!      IF(y_Hmax==9d20 .AND. y_Hmin==-9d20) RETURN
!      y = CALCy(Pf,2)
!      IF(y >= y_Hmax) Fjet1=0
!      IF(y < y_Hmin) Fjet1=0

      if (pseudorapflag) then
      IF((y_Hmax==(1.d100)).AND.(y_Hmin==0d0)) RETURN
      else
      IF((y_Hmax==(dlog(sqrt(shad)/hmass))).AND.(y_Hmin==0d0)) RETURN
      end if

      y = CALCy(Pf,2)
      eta = CALCeta(Pf,2)   

      if (pseudorapflag) then
      IF(Abs(eta) >= y_Hmax) Fjet1=0
      IF(Abs(eta) < y_Hmin) Fjet1=0      
      else
      IF(Abs(y) >= y_Hmax) Fjet1=0
      IF(Abs(y) < y_Hmin) Fjet1=0
      end if

      ELSE
      Fjet1=0
      ENDIF
      END

C-}}}      
C-{{{ FUNCTION: arctan2

      DOUBLE PRECISION FUNCTION arctan2(y,x)
      implicit none
      DOUBLE PRECISION x,y
      include '../commons/common-consts.f'
      
      IF(x == 0d0) THEN
      IF(y == 0d0) THEN
       print *, '######### !!! Error !!! ########   ----   arctan(0/0)'
      stop
      ENDIF
      ENDIF
      
      IF(y >= 0d0) arctan2=dacos(x/(dsqrt(x**2+y**2)))
      IF(y < 0d0) arctan2=2d0*Pi-dacos(x/(dsqrt(x**2+y**2)))
     
      END

C-}}}
C-{{{ FUNCTION: CALCpT

      DOUBLE PRECISION FUNCTION CALCpT(P,nr)
      implicit none
      DOUBLE PRECISION P(4,3)
      INTEGER nr
      CALCpT = dsqrt(P(1,nr)**2+P(2,nr)**2)
      END

C-}}}
C-{{{ FUNCTION: CALCy

      DOUBLE PRECISION FUNCTION CALCy(P,nr)
      implicit none
      DOUBLE PRECISION P(4,3)
      INTEGER nr
      CALCy = 1d0/2d0*dlog((P(4,nr)+P(3,nr))/(P(4,nr)-P(3,nr)))
      END

C-}}}
C-{{{ FUNCTION: CALCeta

      DOUBLE PRECISION FUNCTION CALCeta(P,nr)
      implicit none
      DOUBLE PRECISION P(4,3)
      DOUBLE PRECISION Pabs
      INTEGER nr
      Pabs = Sqrt(P(1,nr)**2 + P(2,nr)**2 + P(3,nr)**2)
      CALCeta = 1d0/2d0*dlog((Pabs+P(3,nr))/(Pabs-P(3,nr)))
      END

C-}}}
com-- 
com-- end - Fjet functions
com-- 

com-- 
com-- PDF functions - begin
com-- 

C-{{{ FUNCTION: PDFsbbbar

      DOUBLE PRECISION FUNCTION PDFsbbbar(eta1,eta2,muF,order)
      implicit none
      DOUBLE PRECISION eta1,eta2,muF
      real*8 ff1(-6:6),ff2(-6:6)
      INTEGER order
      CHARACTER prefix*50
      include '../commons/common-lumis.f'

      PDFsbbbar=0d0
      IF(eta1>1d0) RETURN
      IF(eta2>1d0) RETURN

      call pdfs(1,order,eta1,muf,ff1)
      call pdfs(ppbarsign,order,eta2,muf,ff2)

      PDFsbbbar = ff1(ptn1)*ff2(ptn2) + ff1(ptn2)*ff2(ptn1)

      END

C-}}}      
C-{{{ FUNCTION: PDFsgb

      DOUBLE PRECISION FUNCTION PDFsgb(eta1,eta2,muF,order)
      implicit none
      DOUBLE PRECISION eta1,eta2,muF
      INTEGER order
      real*8 ff1(-6:6),ff2(-6:6)
      include '../commons/common-lumis.f'

      PDFsgb=0d0
      IF(eta1>1d0) RETURN
      IF(eta2>1d0) RETURN

      call pdfs(1,order,eta1,muf,ff1)
      call pdfs(ppbarsign,order,eta2,muf,ff2)

      PDFsgb = ff1(0) * ( ff2(ptn1) + ff2(ptn2) )

      END

C-}}}
C-{{{ FUNCTION: PDFsbg

      DOUBLE PRECISION FUNCTION PDFsbg(eta1,eta2,muF,order)
      implicit none
      DOUBLE PRECISION eta1,eta2,muF
      INTEGER order
      real*8 ff1(-6:6),ff2(-6:6)
      include '../commons/common-lumis.f'

      PDFsbg=0d0
      IF(eta1>1d0) RETURN
      IF(eta2>1d0) RETURN

      call pdfs(1,order,eta1,muf,ff1)
      call pdfs(ppbarsign,order,eta2,muf,ff2)

      PDFsbg = ( ff1(ptn1) + ff1(ptn2) ) * ff2(0) 

      END

C-}}}
com-- 
com-- end - PDF functions
com-- 
C-{{{ FUNCTION: dot

      DOUBLE PRECISION FUNCTION dot(T,x,K,y)
      DOUBLE PRECISION T(4,3)
      DOUBLE PRECISION K(4,3)
      INTEGER x, y
      dot=T(4,x)*K(4,y)-T(1,x)*K(1,y)-T(2,x)*K(2,y)-T(3,x)*K(3,y)
      RETURN
      END

C-}}}
C-{{{ FUNCTION: theta

      DOUBLE PRECISION FUNCTION theta(x,eta)
      DOUBLE PRECISION x,eta
      theta=0d0
      IF(x>=eta) theta=1d0
      END

C-}}}
C-{{{ SUBROUTINE: LorentzTrafo2

      SUBROUTINE LorentzTrafo2(PiTrafo,PfTrafo,PiA,Pf)
      DOUBLE PRECISION PiTrafo(4,3),PfTrafo(4,3),PiA(4,3),Pf(4,3)
     f ,beta,Eprime,gamma
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-fraction.f'

      beta=(eta2-eta1)/(eta1+eta2)
      Eprime=E/dsqrt(eta1*eta2) !!! E is undefined!?!?
      gamma=1/dsqrt(1-beta**2)

      PiTrafo(4,1)=eta1*Eprime
      PiTrafo(1,1)=0d0
      PiTrafo(2,1)=0d0
      PiTrafo(3,1)=eta1*Eprime

      PiTrafo(4,2)=eta2*Eprime
      PiTrafo(1,2)=0d0
      PiTrafo(2,2)=0d0
      PiTrafo(3,2)=-eta2*Eprime

      PfTrafo(4,1)= gamma *(Pf(4,1)-beta*Pf(3,1))
      PfTrafo(1,1)= Pf(1,1)
      PfTrafo(2,1)= Pf(2,1)
      PfTrafo(3,1)= gamma *(Pf(3,1)-beta*Pf(4,1))
      
      PfTrafo(4,2)= gamma *(Pf(4,2)-beta*Pf(3,2))
      PfTrafo(1,2)= Pf(1,2)
      PfTrafo(2,2)= Pf(2,2)
      PfTrafo(3,2)= gamma *(Pf(3,2)-beta*Pf(4,2))
      END

C-}}}
C-{{{ SUBROUTINE: LorentzTrafo3

      SUBROUTINE LorentzTrafo3(PiTrafo,PfTrafo,PiA,Pf)
      DOUBLE PRECISION PiTrafo(4,3),PfTrafo(4,3),PiA(4,3),Pf(4,3)
     f ,beta,Eprime,gamma
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-fraction.f'

      beta=(eta2-eta1)/(eta1+eta2)
      Eprime=dsqrt(Shad)/2d0
      gamma=1d0/dsqrt(1-beta**2)

      PiTrafo(4,1)=eta1*Eprime
      PiTrafo(1,1)=0d0
      PiTrafo(2,1)=0d0
      PiTrafo(3,1)=eta1*Eprime
      
      PiTrafo(4,2)=eta2*Eprime
      PiTrafo(1,2)=0d0
      PiTrafo(2,2)=0d0
      PiTrafo(3,2)=-eta2*Eprime

      PiTrafo(1,3)=0d0
      PiTrafo(2,3)=0d0
      PiTrafo(3,3)=0d0
      PiTrafo(4,3)=0d0

      PfTrafo(4,1)= gamma *(Pf(4,1)-beta*Pf(3,1))
      PfTrafo(1,1)= Pf(1,1)
      PfTrafo(2,1)= Pf(2,1)
      PfTrafo(3,1)= gamma *(Pf(3,1)-beta*Pf(4,1))

      PfTrafo(4,2)= gamma *(Pf(4,2)-beta*Pf(3,2))
      PfTrafo(1,2)= Pf(1,2)
      PfTrafo(2,2)= Pf(2,2)
      PfTrafo(3,2)= gamma *(Pf(3,2)-beta*Pf(4,2))

      PfTrafo(4,3)= gamma *(Pf(4,3)-beta*Pf(3,3))
      PfTrafo(1,3)= Pf(1,3)
      PfTrafo(2,3)= Pf(2,3)
      PfTrafo(3,3)= gamma *(Pf(3,3)-beta*Pf(4,3))

      END

C-}}}
C-{{{ FUNCTION: P_qq

      DOUBLE PRECISION FUNCTION P_qq(z)
      implicit none
      DOUBLE PRECISION z
      include '../commons/common-consts.f'
      P_qq=cF*(1+z**2)/(1-z)  ! +3d0/2d0*delta(1-z)
      END

C-}}}
C-{{{ FUNCTION: Psplitgq

      DOUBLE PRECISION FUNCTION Psplitgq(x)
      implicit none
      DOUBLE PRECISION x
      include '../commons/common-consts.f'
      Psplitgq=0d0
      Psplitgq= tR*(x**2 + (1d0 - x)**2)
      END

C-}}}
C-{{{ FUNCTION: Kbargq

      DOUBLE PRECISION FUNCTION Kbargq(x)
      implicit none
      DOUBLE PRECISION x,Psplitgq
      include '../commons/common-consts.f'
      Kbargq=0d0
      Kbargq=Psplitgq(x) * dlog((1d0 - x)/x) + tR*2d0*x*(1d0 - x)
      END

C-}}}
C-{{{ FUNCTION: Ktildegq

      DOUBLE PRECISION FUNCTION Ktildegq(x)
      implicit none
      DOUBLE PRECISION x,Psplitgq
      Ktildegq=0d0
      Ktildegq=Psplitgq(x) * dlog(1 - x)
      END

C-}}}

C-{{{ SUBROUTINE: xp

      SUBROUTINE xp(p12t,y,z,p12,x)
      DOUBLE PRECISION p12t(4,3),p12(4,3)
      DOUBLE PRECISION x
      INTEGER y,z
      p12t(1,y)=x*p12(1,y)
      p12t(2,y)=x*p12(2,y)
      p12t(3,y)=x*p12(3,y)
      p12t(4,y)=x*p12(4,y)
      p12t(1,z)=p12(1,z)
      p12t(2,z)=p12(2,z)
      p12t(3,z)=p12(3,z)
      p12t(4,z)=p12(4,z)
      RETURN
      END

C-}}}
C-{{{ SUBROUTINE: kkk

      SUBROUTINE kkk(p12t,p34t,y,z,p12,a,b,p345)
      DOUBLE PRECISION p34t(4,3),p12t(4,3),p12(4,3),p345(4,3),
     %K(4,3)
      DOUBLE PRECISION dot
      INTEGER y,z,a,b

      K(1,1)=p12(1,1)+p12(1,2)-p345(1,a)    ! K^mu
      K(2,1)=p12(2,1)+p12(2,2)-p345(2,a)
      K(3,1)=p12(3,1)+p12(3,2)-p345(3,a)
      K(4,1)=p12(4,1)+p12(4,2)-p345(4,a)
      K(1,2)=p12t(1,y)+p12(1,z)             ! Ktilde^mu
      K(2,2)=p12t(2,y)+p12(2,z)
      K(3,2)=p12t(3,y)+p12(3,z)
      K(4,2)=p12t(4,y)+p12(4,z)
      K(1,3)=K(1,1)+K(1,2)                  ! K^mu+Ktilde^mu
      K(2,3)=K(2,1)+K(2,2)
      K(3,3)=K(3,1)+K(3,2)
      K(4,3)=K(4,1)+K(4,2)

      p34t(1,1)=p345(1,b)-2.d0*dot(p345,b,K,3)/dot(K,3,K,3)*K(1,3)+2.d0*    ! parton momenutum
     %dot(p345,b,K,1)/dot(K,1,K,1)*K(1,2)
      p34t(2,1)=p345(2,b)-2.d0*dot(p345,b,K,3)/dot(K,3,K,3)*K(2,3)+2.d0*
     %dot(p345,b,K,1)/dot(K,1,K,1)*K(2,2)
      p34t(3,1)=p345(3,b)-2.d0*dot(p345,b,K,3)/dot(K,3,K,3)*K(3,3)+2.d0*
     %dot(p345,b,K,1)/dot(K,1,K,1)*K(3,2)
      p34t(4,1)=p345(4,b)-2.d0*dot(p345,b,K,3)/dot(K,3,K,3)*K(4,3)+2.d0*
     %dot(p345,b,K,1)/dot(K,1,K,1)*K(4,2)
      p34t(1,2)=p345(1,2)-2.d0*dot(p345,2,K,3)/dot(K,3,K,3)*K(1,3)+2.d0*    ! Higgs momentum
     %dot(p345,2,K,1)/dot(K,1,K,1)*K(1,2)
      p34t(2,2)=p345(2,2)-2.d0*dot(p345,2,K,3)/dot(K,3,K,3)*K(2,3)+2.d0*
     %dot(p345,2,K,1)/dot(K,1,K,1)*K(2,2)
      p34t(3,2)=p345(3,2)-2.d0*dot(p345,2,K,3)/dot(K,3,K,3)*K(3,3)+2.d0*
     %dot(p345,2,K,1)/dot(K,1,K,1)*K(3,2)
      p34t(4,2)=p345(4,2)-2.d0*dot(p345,2,K,3)/dot(K,3,K,3)*K(4,3)+2.d0*
     %dot(p345,2,K,1)/dot(K,1,K,1)*K(4,2)
      RETURN
      END

C-}}}
C-{{{ subroutine: xtoP_1

      subroutine xtoP_1(x,PiA,Pf)
      implicit none
      double precision x(10),PiA(4,3),Pf(4,3),E3,E4,costheta,
     f sintheta,E,k5,sinchi,coschi,x1,x2
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
*      print *,x(1),x(2),x(3),x(4),x(5),x(6)
      x1 = x(1)
      x2 = x(2)
C      x1=0.1d0
C      x2=0.2d0
      E = dsqrt(x1*x2*Shad)/2d0
** initial state momenta
      PiA(4,1) = E
      PiA(1,1) = 0d0
      PiA(2,1) = 0d0
      PiA(3,1) = E
      PiA(4,2) = E
      PiA(1,2) = 0d0
      PiA(2,2) = 0d0
      PiA(3,2) = -E

*      costheta= 2d0*x(3)-1d0
      costheta=x(3)
C      costheta=0.3d0
      sintheta=dsin(dacos(costheta))
*      print*, costheta
      E3=(4*E**2-Hmass**2)/(4*E)
      E4=(4*E**2+Hmass**2)/(4*E)

      Pf(4,1) = E3
      Pf(1,1) = 0d0
      Pf(2,1) = E3*sintheta
      Pf(3,1) = E3*costheta

      Pf(4,2)=E4
      Pf(1,2)=-Pf(1,1)
      Pf(2,2)=-Pf(2,1)
      Pf(3,2)=-Pf(3,1)
C      print*, x(1),x(2),x(3),dsqrt(Pf(1,2)**2+Pf(2,2)**2),sintheta,E
      end

C-}}}
C-{{{ subroutine: xtoP_12

      subroutine xtoP_12(x1,x2,x3,PiA,Pf)
      implicit none
      double precision x1,x2,x3,PiA(4,3),Pf(4,3),E3,E4,costheta,
     f sintheta,E,k5,sinchi,coschi
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
*      print *,x(1),x(2),x(3),x(4),x(5),x(6)
C      x1=0.1d0
C      x2=0.2d0
      E = dsqrt(x1*x2*Shad)/2d0
** initial state momenta
      PiA(4,1) = E
      PiA(1,1) = 0d0
      PiA(2,1) = 0d0
      PiA(3,1) = E
      PiA(4,2) = E
      PiA(1,2) = 0d0
      PiA(2,2) = 0d0
      PiA(3,2) = -E

*      costheta= 2d0*x3-1d0
      costheta=x3
C      costheta=0.3d0
      sintheta=dsin(dacos(costheta))
*      print*, costheta
      E3=(4*E**2-Hmass**2)/(4*E)
      E4=(4*E**2+Hmass**2)/(4*E)

      Pf(4,1) = E3
      Pf(1,1) = 0d0
      Pf(2,1) = E3*sintheta
      Pf(3,1) = E3*costheta

      Pf(4,2)=E4
      Pf(1,2)=-Pf(1,1)
      Pf(2,2)=-Pf(2,1)
      Pf(3,2)=-Pf(3,1)
C      print*, x(1),x(2),x(3),dsqrt(Pf(1,2)**2+Pf(2,2)**2),sintheta,E
      end

C-}}}
