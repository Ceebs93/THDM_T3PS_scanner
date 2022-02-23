      real*8 function gghmix_CT6(Mtop,Mstop1,Mstop2,Mgluino,mu,
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

      gghmix_CT6 =
     & (-(calpha*lMtg*Mgluino**2*(Mstop1**4*Mstop2**2 + 
     &            Mstop1**2*Mstop2**4 - Mstop1**4*Mtop**2 - 
     &            Mstop2**4*Mtop**2))/(18.*Mstop1**4*Mstop2**4) + 
     &      (calpha*lMt2*(Mstop1**2*Mstop2**4 + Mstop2**6 + 
     &           4*Mgluino**2*Mstop1**2*Mtop**2 - 
     &           2*Mstop1**2*Mstop2**2*Mtop**2 + 4*Mstop1**2*Mtop**4))/
     &       (36.*Mstop1**2*Mstop2**4) + 
     &      (calpha*lMt1*(Mstop1**6 + Mstop1**4*Mstop2**2 + 
     &           4*Mgluino**2*Mstop2**2*Mtop**2 - 
     &           2*Mstop1**2*Mstop2**2*Mtop**2 + 4*Mstop2**2*Mtop**4))/
     &       (36.*Mstop1**4*Mstop2**2) + 
     &      (calpha*lMtop*(-2*Mgluino**2*Mstop1**4*Mstop2**2 + 
     &           Mstop1**6*Mstop2**2 - 
     &           2*Mgluino**2*Mstop1**2*Mstop2**4 + 
     &           2*Mstop1**4*Mstop2**4 + Mstop1**2*Mstop2**6 + 
     &           6*Mgluino**2*Mstop1**4*Mtop**2 - 
     &           8*Mstop1**4*Mstop2**2*Mtop**2 + 
     &           6*Mgluino**2*Mstop2**4*Mtop**2 - 
     &           8*Mstop1**2*Mstop2**4*Mtop**2 + 6*Mstop1**4*Mtop**4 + 
     &           6*Mstop2**4*Mtop**4))/(36.*Mstop1**4*Mstop2**4) + 
     &      (calpha*(-2*Mgluino**2*Mstop1**4*Mstop2**2 + 
     &           Mstop1**6*Mstop2**2 - 
     &           2*Mgluino**2*Mstop1**2*Mstop2**4 + 
     &           18*Mstop1**4*Mstop2**4 + Mstop1**2*Mstop2**6 + 
     &           6*Mgluino**2*Mstop1**4*Mtop**2 - 
     &           4*Mstop1**4*Mstop2**2*Mtop**2 + 
     &           6*Mgluino**2*Mstop2**4*Mtop**2 - 
     &           4*Mstop1**2*Mstop2**4*Mtop**2 + 6*Mstop1**4*Mtop**4 + 
     &           6*Mstop2**4*Mtop**4))/(36.*Mstop1**4*Mstop2**4) + 
     &      ((-4*calpha*lMtop*(Mstop1 - Mstop2)**4*
     &            (Mstop1 + Mstop2)**4)/(9.*Mstop1**4*Mstop2**4) - 
     &         (4*calpha*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (Mstop1**4 - Mstop1**2*Mstop2**2 + Mstop2**4))/
     &          (9.*Mstop1**4*Mstop2**4) - 
     &         (2*calpha*lMt2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**6 - 4*Mstop1**4*Mstop2**2 + 
     &              2*Mstop1**2*Mstop2**4 - Mstop2**6))/
     &          (9.*Mstop1**4*Mstop2**4) - 
     &         (2*calpha*lMt1*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**6 - 2*Mstop1**4*Mstop2**2 + 
     &              4*Mstop1**2*Mstop2**4 - Mstop2**6))/
     &          (9.*Mstop1**4*Mstop2**4))*sthetat**6 + 
     &      ((2*calpha*lMtop*(Mstop1 - Mstop2)**4*
     &            (Mstop1 + Mstop2)**4)/(9.*Mstop1**4*Mstop2**4) + 
     &         (2*calpha*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (Mstop1**4 - Mstop1**2*Mstop2**2 + Mstop2**4))/
     &          (9.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMt2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**6 - 4*Mstop1**4*Mstop2**2 + 
     &              2*Mstop1**2*Mstop2**4 - Mstop2**6))/
     &          (9.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMt1*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**6 - 2*Mstop1**4*Mstop2**2 + 
     &              4*Mstop1**2*Mstop2**4 - Mstop2**6))/
     &          (9.*Mstop1**4*Mstop2**4))*sthetat**8 - 
     &      (calpha*Mtop**2*(-Mgluino**2 + Mstop1**2 - Mtop**2)*
     &         B0fin(Mstop1**2,Mtop,Mgluino,mu))/(18.*Mstop1**4) - 
     &      (calpha*Mtop**2*(-Mgluino**2 + Mstop2**2 - Mtop**2)*
     &         B0fin(Mstop2**2,Mtop,Mgluino,mu))/(18.*Mstop2**4) + 
     &      sthetat**2*(-(calpha*lMtg*Mgluino**2*(Mstop1 - Mstop2)**2*
     &             (Mstop1 + Mstop2)**2*(Mstop1**2 + Mstop2**2))/
     &          (18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (-3*Mgluino**2*Mstop1**2 - 3*Mgluino**2*Mstop2**2 + 
     &              2*Mstop1**2*Mstop2**2 - 7*Mstop1**2*Mtop**2 - 
     &              7*Mstop2**2*Mtop**2))/(18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMtop*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (-3*Mgluino**2*Mstop1**2 - 3*Mgluino**2*Mstop2**2 + 
     &              4*Mstop1**2*Mstop2**2 - 7*Mstop1**2*Mtop**2 - 
     &              7*Mstop2**2*Mtop**2))/(18.*Mstop1**4*Mstop2**4) - 
     &         (calpha*lMt2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (2*Mgluino**2*Mstop1**4 - Mstop1**4*Mstop2**2 + 
     &              3*Mstop1**2*Mstop2**4 + 4*Mstop1**4*Mtop**2 - 
     &              2*Mstop1**2*Mstop2**2*Mtop**2 - 2*Mstop2**4*Mtop**2
     &              ))/(18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMt1*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (3*Mstop1**4*Mstop2**2 + 2*Mgluino**2*Mstop2**4 - 
     &              Mstop1**2*Mstop2**4 - 2*Mstop1**4*Mtop**2 - 
     &              2*Mstop1**2*Mstop2**2*Mtop**2 + 4*Mstop2**4*Mtop**2
     &              ))/(18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*
     &            (Mgluino**2 - Mstop1**2 + Mtop**2)*
     &            B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &          (18.*Mstop1**4*Mstop2**2) - 
     &         (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*
     &            (Mgluino**2 - Mstop2**2 + Mtop**2)*
     &            B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &          (18.*Mstop1**2*Mstop2**4)) + 
     &      sthetat**4*((calpha*lMtg*Mgluino**2*(Mstop1 - Mstop2)**2*
     &            (Mstop1 + Mstop2)**2*(Mstop1**2 + Mstop2**2))/
     &          (18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMtop*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (3*Mgluino**2*Mstop1**2 + 4*Mstop1**4 + 
     &              3*Mgluino**2*Mstop2**2 - 12*Mstop1**2*Mstop2**2 + 
     &              4*Mstop2**4 + 7*Mstop1**2*Mtop**2 + 
     &              7*Mstop2**2*Mtop**2))/(18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*(Mstop1 - Mstop2)**2*(Mstop1 + Mstop2)**2*
     &            (3*Mgluino**2*Mstop1**2 + 4*Mstop1**4 + 
     &              3*Mgluino**2*Mstop2**2 - 6*Mstop1**2*Mstop2**2 + 
     &              4*Mstop2**4 + 7*Mstop1**2*Mtop**2 + 
     &              7*Mstop2**2*Mtop**2))/(18.*Mstop1**4*Mstop2**4) + 
     &         (calpha*lMt2*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (2*Mgluino**2*Mstop1**4 + 2*Mstop1**6 - 
     &              9*Mstop1**4*Mstop2**2 + 7*Mstop1**2*Mstop2**4 - 
     &              2*Mstop2**6 + 4*Mstop1**4*Mtop**2 - 
     &              2*Mstop1**2*Mstop2**2*Mtop**2 - 2*Mstop2**4*Mtop**2
     &              ))/(18.*Mstop1**4*Mstop2**4) - 
     &         (calpha*lMt1*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (-2*Mstop1**6 + 7*Mstop1**4*Mstop2**2 + 
     &              2*Mgluino**2*Mstop2**4 - 9*Mstop1**2*Mstop2**4 + 
     &              2*Mstop2**6 - 2*Mstop1**4*Mtop**2 - 
     &              2*Mstop1**2*Mstop2**2*Mtop**2 + 4*Mstop2**4*Mtop**2
     &              ))/(18.*Mstop1**4*Mstop2**4) - 
     &         (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*
     &            (Mgluino**2 - Mstop1**2 + Mtop**2)*
     &            B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &          (18.*Mstop1**4*Mstop2**2) + 
     &         (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*
     &            (Mgluino**2 - Mstop2**2 + Mtop**2)*
     &            B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &          (18.*Mstop1**2*Mstop2**4)) + 
     &      (calpha*(Mstop1**2 + Mstop2**2)*
     &         (Mgluino**2 - Mstop1**2 + Mtop**2)*
     &         B0fin(Mtop**2,Mgluino,Mstop1,mu))/
     &       (36.*Mstop1**2*Mstop2**2) - 
     &      (calpha*(Mstop1**2 + Mstop2**2)*
     &         (-Mgluino**2 + Mstop2**2 - Mtop**2)*
     &         B0fin(Mtop**2,Mgluino,Mstop2,mu))/
     &       (36.*Mstop1**2*Mstop2**2) + 
     &      cthetat*sthetat**5*((2*calpha*lMtop*Mgluino*
     &            (Mstop1 - Mstop2)**3*(Mstop1 + Mstop2)**3*Mtop)/
     &          (9.*Mstop1**4*Mstop2**4) + 
     &         (2*calpha*lMt2*Mgluino*(Mstop1**2 - 3*Mstop2**2)*Mtop)/
     &          (9.*Mstop2**4) + 
     &         (2*calpha*lMt1*Mgluino*(3*Mstop1**2 - Mstop2**2)*Mtop)/
     &          (9.*Mstop1**4) + 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**4 + Mstop2**4)*Mtop)/
     &          (9.*Mstop1**4*Mstop2**4) + 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**4*Mstop2**2) + 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**2*Mstop2**4) - 
     &         (8*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            Mtop*B0fin(q0**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**2*Mstop2**2)) + 
     &      cthetat*sthetat*((2*calpha*Mgluino*(Mstop1 - Mstop2)*
     &            (Mstop1 + Mstop2)*(Mstop1**2 + Mstop2**2)*Mtop**3)/
     &          (9.*Mstop1**4*Mstop2**4) + 
     &         (2*calpha*lMt1*Mgluino*(Mstop1 - Mtop)*Mtop*
     &            (Mstop1 + Mtop))/(9.*Mstop1**4) - 
     &         (2*calpha*lMt2*Mgluino*(Mstop2 - Mtop)*Mtop*
     &            (Mstop2 + Mtop))/(9.*Mstop2**4) - 
     &         (2*calpha*lMtop*Mgluino*(Mstop1 - Mstop2)*
     &            (Mstop1 + Mstop2)*Mtop*
     &            (Mstop1**2*Mstop2**2 - Mstop1**2*Mtop**2 - 
     &              Mstop2**2*Mtop**2))/(9.*Mstop1**4*Mstop2**4) - 
     &         (2*calpha*Mgluino*Mtop**3*
     &            B0fin(Mstop1**2,Mtop,Mgluino,mu))/(9.*Mstop1**4) + 
     &         (2*calpha*Mgluino*Mtop**3*
     &            B0fin(Mstop2**2,Mtop,Mgluino,mu))/(9.*Mstop2**4) - 
     &         (calpha*Mgluino*(Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mtop**2,Mgluino,Mstop1,mu))/
     &          (9.*Mstop1**2*Mstop2**2) + 
     &         (calpha*Mgluino*(Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mtop**2,Mgluino,Mstop2,mu))/
     &          (9.*Mstop1**2*Mstop2**2) - 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            Mtop*B0fin(q0**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**2*Mstop2**2)) + 
     &      cthetat*sthetat**3*((-2*calpha*lMtop*Mgluino*
     &            (Mstop1 - Mstop2)**3*(Mstop1 + Mstop2)**3*Mtop)/
     &          (9.*Mstop1**4*Mstop2**4) - 
     &         (2*calpha*lMt2*Mgluino*(Mstop1**2 - 3*Mstop2**2)*Mtop)/
     &          (9.*Mstop2**4) - 
     &         (2*calpha*lMt1*Mgluino*(3*Mstop1**2 - Mstop2**2)*Mtop)/
     &          (9.*Mstop1**4) - 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**4 + Mstop2**4)*Mtop)/
     &          (9.*Mstop1**4*Mstop2**4) - 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mstop1**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**4*Mstop2**2) - 
     &         (2*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (Mstop1**2 + Mstop2**2)*Mtop*
     &            B0fin(Mstop2**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**2*Mstop2**4) + 
     &         (8*calpha*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            Mtop*B0fin(q0**2,Mtop,Mgluino,mu))/
     &          (9.*Mstop1**2*Mstop2**2)))/sbeta
      end
