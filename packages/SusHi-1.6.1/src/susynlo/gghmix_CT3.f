      real*8 function gghmix_CT3(Mtop,Mstop1,Mstop2,Mgluino,mu,
     & q0,cthetat,sthetat,alpha,beta,
     & muSUSY,c1EWmc2EW,c1EWpc2EW)
c..
c..
      implicit real*8 (a-z)
      data pi/3.14159265358979323846264338328d0/

      mu2 = mu*mu
      lMtop    = dlog(mu2/Mtop/Mtop)
      lMtg = dlog(Mtop*Mtop/Mgluino/Mgluino)
      lMt1 = dlog(Mtop*Mtop/Mstop1/Mstop1)
      lMt2 = dlog(Mtop*Mtop/Mstop2/Mstop2)
c      cthetat = dsqrt(1.d0 - sthetat**2)
c      beta = datan(tanb)
      sbeta = dsin(beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      calphambeta = dcos(alpha - beta)
      q02 = q0*q0

      gghmix_CT3 =
     & (calphambeta*muSUSY*(-(lMtop*Mgluino*Mtop**2)/
     &         (9.*Mstop1**2*Mstop2**2) + 
     &        (lMt1*Mgluino*Mtop**2)/
     &         (9.*Mstop1**2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)) + 
     &        (lMt2*Mgluino*Mtop**2)/
     &         (9.*Mstop2**2*(-Mstop1 + Mstop2)*(Mstop1 + Mstop2)) + 
     &        cthetat*((2*lMtop*(Mstop1 - Mstop2)**3*
     &              (Mstop1 + Mstop2)**3*Mtop)/(9.*Mstop1**4*Mstop2**4)
     &             + (2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (Mstop1**4 - Mstop1**2*Mstop2**2 + Mstop2**4)*Mtop)
     &             /(9.*Mstop1**4*Mstop2**4) + 
     &           (lMt2*(Mstop1**6 - 4*Mstop1**4*Mstop2**2 + 
     &                2*Mstop1**2*Mstop2**4 - Mstop2**6)*Mtop)/
     &            (9.*Mstop1**4*Mstop2**4) + 
     &           (lMt1*(Mstop1**6 - 2*Mstop1**4*Mstop2**2 + 
     &                4*Mstop1**2*Mstop2**4 - Mstop2**6)*Mtop)/
     &            (9.*Mstop1**4*Mstop2**4))*sthetat**3 + 
     &        cthetat*((-2*lMtop*(Mstop1 - Mstop2)**3*
     &              (Mstop1 + Mstop2)**3*Mtop)/(9.*Mstop1**4*Mstop2**4)
     &             - (2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (Mstop1**4 - Mstop1**2*Mstop2**2 + Mstop2**4)*Mtop)
     &             /(9.*Mstop1**4*Mstop2**4) - 
     &           (lMt2*(Mstop1**6 - 4*Mstop1**4*Mstop2**2 + 
     &                2*Mstop1**2*Mstop2**4 - Mstop2**6)*Mtop)/
     &            (9.*Mstop1**4*Mstop2**4) - 
     &           (lMt1*(Mstop1**6 - 2*Mstop1**4*Mstop2**2 + 
     &                4*Mstop1**2*Mstop2**4 - Mstop2**6)*Mtop)/
     &            (9.*Mstop1**4*Mstop2**4))*sthetat**5 + 
     &        cthetat*sthetat*((lMtg*Mgluino**2*(Mstop1 - Mstop2)*
     &              (Mstop1 + Mstop2)*
     &              (Mstop1**2*Mstop2**2 - 2*Mstop1**2*Mtop**2 - 
     &                2*Mstop2**2*Mtop**2))/
     &            (36.*Mstop1**4*Mstop2**4*Mtop) - 
     &           (lMt2*(Mstop1**2*Mstop2**4 - Mstop2**6 + 
     &                8*Mgluino**2*Mstop1**2*Mtop**2 - 
     &                4*Mstop1**2*Mstop2**2*Mtop**2 + 
     &                4*Mstop2**4*Mtop**2 + 8*Mstop1**2*Mtop**4))/
     &            (72.*Mstop1**2*Mstop2**4*Mtop) + 
     &           ((Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (2*Mgluino**2*Mstop1**2*Mstop2**2 - 
     &                Mstop1**4*Mstop2**2 - Mstop1**2*Mstop2**4 - 
     &                12*Mgluino**2*Mstop1**2*Mtop**2 - 
     &                12*Mgluino**2*Mstop2**2*Mtop**2 + 
     &                2*Mstop1**2*Mstop2**2*Mtop**2 - 
     &                12*Mstop1**2*Mtop**4 - 12*Mstop2**2*Mtop**4))/
     &            (72.*Mstop1**4*Mstop2**4*Mtop) + 
     &           (lMtop*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (2*Mgluino**2*Mstop1**2*Mstop2**2 - 
     &                Mstop1**4*Mstop2**2 - Mstop1**2*Mstop2**4 - 
     &                12*Mgluino**2*Mstop1**2*Mtop**2 - 
     &                12*Mgluino**2*Mstop2**2*Mtop**2 + 
     &                14*Mstop1**2*Mstop2**2*Mtop**2 - 
     &                12*Mstop1**2*Mtop**4 - 12*Mstop2**2*Mtop**4))/
     &            (72.*Mstop1**4*Mstop2**4*Mtop) - 
     &           (lMt1*(Mstop1**6 - Mstop1**4*Mstop2**2 - 
     &                4*Mstop1**4*Mtop**2 - 
     &                8*Mgluino**2*Mstop2**2*Mtop**2 + 
     &                4*Mstop1**2*Mstop2**2*Mtop**2 - 
     &                8*Mstop2**2*Mtop**4))/
     &            (72.*Mstop1**4*Mstop2**2*Mtop) - 
     &           (Mtop*(-Mgluino**2 + Mstop1**2 - Mtop**2)*
     &              B0fin(Mstop1**2,Mtop,Mgluino,mu))/(18.*Mstop1**4)
     &            + (Mtop*(-Mgluino**2 + Mstop2**2 - Mtop**2)*
     &              B0fin(Mstop2**2,Mtop,Mgluino,mu))/(18.*Mstop2**4)
     &            - ((Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (Mgluino**2 - Mstop1**2 + Mtop**2)*
     &              B0fin(Mtop**2,Mgluino,Mstop1,mu))/
     &            (72.*Mstop1**2*Mstop2**2*Mtop) - 
     &           ((Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              (Mgluino**2 - Mstop2**2 + Mtop**2)*
     &              B0fin(Mtop**2,Mgluino,Mstop2,mu))/
     &            (72.*Mstop1**2*Mstop2**2*Mtop)) - 
     &        (Mgluino*Mtop**2*B0fin(q0**2,Mtop,Mgluino,mu))/
     &         (9.*Mstop1**2*Mstop2**2) + 
     &        sthetat**4*((2*lMtop*Mgluino*(Mstop1 - Mstop2)**2*
     &              (Mstop1 + Mstop2)**2*Mtop**2)/
     &            (9.*Mstop1**4*Mstop2**4) + 
     &           (2*lMt2*Mgluino*(Mstop1**2 - 3*Mstop2**2)*Mtop**2)/
     &            (9.*(Mstop1 - Mstop2)*Mstop2**4*(Mstop1 + Mstop2)) + 
     &           (2*lMt1*Mgluino*(3*Mstop1**2 - Mstop2**2)*Mtop**2)/
     &            (9.*Mstop1**4*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)) + 
     &           (2*Mgluino*(Mstop1**4 + Mstop2**4)*Mtop**2)/
     &            (9.*Mstop1**4*Mstop2**4) + 
     &           (2*Mgluino*Mtop**2*B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop1**4) + 
     &           (2*Mgluino*Mtop**2*B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop2**4) - 
     &           (Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              B0fin(Mtop**2,Mgluino,Mstop1,mu))/
     &            (18.*Mstop1**2*Mstop2**2) + 
     &           (Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              B0fin(Mtop**2,Mgluino,Mstop2,mu))/
     &            (18.*Mstop1**2*Mstop2**2) - 
     &           (4*Mgluino*Mtop**2*B0fin(q0**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop1**2*Mstop2**2)) + 
     &        sthetat**2*((-2*lMtop*Mgluino*(Mstop1 - Mstop2)**2*
     &              (Mstop1 + Mstop2)**2*Mtop**2)/
     &            (9.*Mstop1**4*Mstop2**4) - 
     &           (2*lMt2*Mgluino*(Mstop1**2 - 3*Mstop2**2)*Mtop**2)/
     &            (9.*(Mstop1 - Mstop2)*Mstop2**4*(Mstop1 + Mstop2)) - 
     &           (2*lMt1*Mgluino*(3*Mstop1**2 - Mstop2**2)*Mtop**2)/
     &            (9.*Mstop1**4*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)) - 
     &           (2*Mgluino*(Mstop1**4 + Mstop2**4)*Mtop**2)/
     &            (9.*Mstop1**4*Mstop2**4) - 
     &           (2*Mgluino*Mtop**2*B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop1**4) - 
     &           (2*Mgluino*Mtop**2*B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop2**4) + 
     &           (Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              B0fin(Mtop**2,Mgluino,Mstop1,mu))/
     &            (18.*Mstop1**2*Mstop2**2) - 
     &           (Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &              B0fin(Mtop**2,Mgluino,Mstop2,mu))/
     &            (18.*Mstop1**2*Mstop2**2) + 
     &           (4*Mgluino*Mtop**2*B0fin(q0**2,Mtop,Mgluino,mu))/
     &            (9.*Mstop1**2*Mstop2**2))))/sbeta**2
      end
