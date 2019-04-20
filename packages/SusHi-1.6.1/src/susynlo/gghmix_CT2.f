      real*8 function gghmix_CT2(Mtop,Mstop1,Mstop2,Mgluino,mu,
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

      gghmix_CT2 =
     & 0
      end
