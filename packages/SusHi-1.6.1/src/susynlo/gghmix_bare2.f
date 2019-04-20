      real*8 function gghmix_bare2(Mtop,Mstop1,Mstop2,Mgluino,mu,
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

      gghmix_bare2 =
     & c1EWpc2EW*(-(Mtop**2*(-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &             13*Mgluino**4*Mstop2**4 - 21*Mgluino**2*Mstop2**6 + 
     &             18*Mstop2**8 + 10*Mgluino**6*Mtop**2 - 
     &             19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &             32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &             9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &             71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &             63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &             27*Mstop2**2*Mtop**6 + 9*Mtop**8)*
     &           davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &            Mgluino))/
     &        (48.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        + (cthetat*Mtop*(Mgluino**10 - 19*Mgluino**8*Mstop2**2 + 
     &            42*Mgluino**6*Mstop2**4 - 22*Mgluino**4*Mstop2**6 - 
     &            11*Mgluino**2*Mstop2**8 + 9*Mstop2**10 - 
     &            11*Mgluino**8*Mtop**2 + 
     &            52*Mgluino**6*Mstop2**2*Mtop**2 - 
     &            46*Mgluino**4*Mstop2**4*Mtop**2 - 
     &            4*Mgluino**2*Mstop2**6*Mtop**2 + 
     &            9*Mstop2**8*Mtop**2 + 20*Mgluino**6*Mtop**4 - 
     &            92*Mgluino**4*Mstop2**2*Mtop**4 + 
     &            24*Mgluino**2*Mstop2**4*Mtop**4 - 
     &            72*Mstop2**6*Mtop**4 - 20*Mgluino**4*Mtop**6 - 
     &            28*Mgluino**2*Mstop2**2*Mtop**6 + 
     &            72*Mstop2**4*Mtop**6 + 19*Mgluino**2*Mtop**8 - 
     &            9*Mstop2**2*Mtop**8 - 9*Mtop**10)*sthetat*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (24.*Mgluino**3*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)))
     &      + c1EWmc2EW*((Mtop**2*
     &          (-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &            13*Mgluino**4*Mstop2**4 - 21*Mgluino**2*Mstop2**6 + 
     &            18*Mstop2**8 + 10*Mgluino**6*Mtop**2 - 
     &            19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &            32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &            9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &            71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &            63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &            27*Mstop2**2*Mtop**6 + 9*Mtop**8)*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (48.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        - (cthetat*Mtop*(Mgluino**10*Mstop1**2 + 
     &            Mgluino**10*Mstop2**2 - 
     &            19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &            7*Mgluino**8*Mstop2**4 + 
     &            42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &            42*Mgluino**6*Mstop2**6 - 
     &            22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &            94*Mgluino**4*Mstop2**8 - 
     &            11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &            85*Mgluino**2*Mstop2**10 + 9*Mstop1**2*Mstop2**10 - 
     &            27*Mstop2**12 - 11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &            11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &            52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &            56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &            46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &            50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &            4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &            40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &            9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &            45*Mstop2**10*Mtop**2 + 
     &            20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &            20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &            92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &            56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &            24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &            120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &            72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &            36*Mstop2**8*Mtop**4 - 
     &            20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &            20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &            28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &            56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &            72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &            108*Mstop2**6*Mtop**6 + 
     &            19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &            19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &            9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &            63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &            9*Mstop2**2*Mtop**10)*sthetat*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (24.*Mgluino**3*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &          (-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        - (Mtop**2*(-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &            13*Mgluino**4*Mstop2**4 - 21*Mgluino**2*Mstop2**6 + 
     &            18*Mstop2**8 + 10*Mgluino**6*Mtop**2 - 
     &            19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &            32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &            9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &            71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &            63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &            27*Mstop2**2*Mtop**6 + 9*Mtop**8)*sthetat**2*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (24.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        + (cthetat*Mtop*(Mgluino**10*Mstop1**2 + 
     &            Mgluino**10*Mstop2**2 - 
     &            19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &            7*Mgluino**8*Mstop2**4 + 
     &            42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &            42*Mgluino**6*Mstop2**6 - 
     &            22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &            94*Mgluino**4*Mstop2**8 - 
     &            11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &            85*Mgluino**2*Mstop2**10 + 9*Mstop1**2*Mstop2**10 - 
     &            27*Mstop2**12 - 11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &            11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &            52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &            56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &            46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &            50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &            4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &            40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &            9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &            45*Mstop2**10*Mtop**2 + 
     &            20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &            20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &            92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &            56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &            24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &            120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &            72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &            36*Mstop2**8*Mtop**4 - 
     &            20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &            20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &            28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &            56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &            72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &            108*Mstop2**6*Mtop**6 + 
     &            19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &            19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &            9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &            63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &            9*Mstop2**2*Mtop**10)*sthetat**3*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (12.*Mgluino**3*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &          (-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)))
     &      + (calphambeta*muSUSY*
     &       (-(Mstop2**2*Mtop**2*
     &             (-Mgluino**2 + 9*Mstop2**2 + 9*Mtop**2)*
     &             (Mgluino**4 - 2*Mgluino**2*Mstop2**2 + Mstop2**4 - 
     &               2*Mstop2**2*Mtop**2 + Mtop**4)*
     &             davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &              Mgluino))/
     &          (12.*Mgluino**3*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (-Mgluino + Mstop2 - Mtop)**2*
     &            (Mgluino + Mstop2 - Mtop)**2*
     &            (-Mgluino + Mstop2 + Mtop)**2*
     &            (Mgluino + Mstop2 + Mtop)**2*
     &            ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &              (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)
     &            ) + (cthetat*Mtop**3*
     &            (-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &              13*Mgluino**4*Mstop2**4 - 
     &              21*Mgluino**2*Mstop2**6 + 18*Mstop2**8 + 
     &              10*Mgluino**6*Mtop**2 - 
     &              19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &              32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &              9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &              71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &              63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &              27*Mstop2**2*Mtop**6 + 9*Mtop**8)*sthetat*
     &            davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &             Mgluino))/
     &          (12.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &            (Mgluino + Mstop2 - Mtop)**3*
     &            (-Mgluino + Mstop2 + Mtop)**3*
     &            (Mgluino + Mstop2 + Mtop)**3*
     &            ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &              (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)
     &            ) - (Mtop**2*(Mgluino**10*Mstop1**2 + 
     &              Mgluino**10*Mstop2**2 - 
     &              19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &              7*Mgluino**8*Mstop2**4 + 
     &              42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &              42*Mgluino**6*Mstop2**6 - 
     &              22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &              94*Mgluino**4*Mstop2**8 - 
     &              11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &              85*Mgluino**2*Mstop2**10 + 
     &              9*Mstop1**2*Mstop2**10 - 27*Mstop2**12 - 
     &              11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &              11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &              52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &              56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &              46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &              50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &              4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &              40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &              9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &              45*Mstop2**10*Mtop**2 + 
     &              20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &              20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &              92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &              56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &              24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &              120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &              72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &              36*Mstop2**8*Mtop**4 - 
     &              20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &              20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &              28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &              56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &              72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &              108*Mstop2**6*Mtop**6 + 
     &              19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &              19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &              9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &              63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &              9*Mstop2**2*Mtop**10)*sthetat**2*
     &            davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &             Mgluino))/
     &          (6.*Mgluino**3*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (-Mgluino + Mstop2 - Mtop)**3*
     &            (Mgluino + Mstop2 - Mtop)**3*
     &            (-Mgluino + Mstop2 + Mtop)**3*
     &            (Mgluino + Mstop2 + Mtop)**3*
     &            ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &              (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)
     &            ) + (Mtop**2*(Mgluino**10*Mstop1**2 + 
     &              Mgluino**10*Mstop2**2 - 
     &              19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &              7*Mgluino**8*Mstop2**4 + 
     &              42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &              42*Mgluino**6*Mstop2**6 - 
     &              22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &              94*Mgluino**4*Mstop2**8 - 
     &              11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &              85*Mgluino**2*Mstop2**10 + 
     &              9*Mstop1**2*Mstop2**10 - 27*Mstop2**12 - 
     &              11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &              11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &              52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &              56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &              46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &              50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &              4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &              40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &              9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &              45*Mstop2**10*Mtop**2 + 
     &              20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &              20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &              92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &              56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &              24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &              120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &              72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &              36*Mstop2**8*Mtop**4 - 
     &              20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &              20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &              28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &              56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &              72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &              108*Mstop2**6*Mtop**6 + 
     &              19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &              19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &              9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &              63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &              9*Mstop2**2*Mtop**10)*sthetat**4*
     &            davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &             Mgluino))/
     &          (6.*Mgluino**3*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*
     &            (-Mgluino + Mstop2 - Mtop)**3*
     &            (Mgluino + Mstop2 - Mtop)**3*
     &            (-Mgluino + Mstop2 + Mtop)**3*
     &            (Mgluino + Mstop2 + Mtop)**3*
     &            ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &              (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)
     &            )))/sbeta**2 + 
     &    ((calpha*Mtop**2*(-Mgluino**2 + Mstop2**2 - Mtop**2)*
     &          (-(Mgluino**6*Mstop2**2) + 11*Mgluino**4*Mstop2**4 - 
     &            19*Mgluino**2*Mstop2**6 + 9*Mstop2**8 - 
     &            Mgluino**6*Mtop**2 + 
     &            30*Mgluino**4*Mstop2**2*Mtop**2 + 
     &            Mgluino**2*Mstop2**4*Mtop**2 + 
     &            18*Mstop2**6*Mtop**2 + 11*Mgluino**4*Mtop**4 + 
     &            Mgluino**2*Mstop2**2*Mtop**4 - 
     &            54*Mstop2**4*Mtop**4 - 19*Mgluino**2*Mtop**6 + 
     &            18*Mstop2**2*Mtop**6 + 9*Mtop**8)*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (12.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        - (calpha*cthetat*Mtop*
     &          (-(Mgluino**10*Mstop2**2) + 13*Mgluino**8*Mstop2**4 - 
     &            42*Mgluino**6*Mstop2**6 + 58*Mgluino**4*Mstop2**8 - 
     &            37*Mgluino**2*Mstop2**10 + 9*Mstop2**12 - 
     &            Mgluino**10*Mtop**2 + 
     &            34*Mgluino**8*Mstop2**2*Mtop**2 - 
     &            80*Mgluino**6*Mstop2**4*Mtop**2 + 
     &            62*Mgluino**4*Mstop2**6*Mtop**2 - 
     &            15*Mgluino**2*Mstop2**8*Mtop**2 + 
     &            11*Mgluino**8*Mtop**4 - 
     &            42*Mgluino**6*Mstop2**2*Mtop**4 - 
     &            26*Mgluino**4*Mstop2**4*Mtop**4 + 
     &            106*Mgluino**2*Mstop2**6*Mtop**4 - 
     &            81*Mstop2**8*Mtop**4 - 20*Mgluino**6*Mtop**6 + 
     &            102*Mgluino**4*Mstop2**2*Mtop**6 - 
     &            38*Mgluino**2*Mstop2**4*Mtop**6 + 
     &            144*Mstop2**6*Mtop**6 + 20*Mgluino**4*Mtop**8 + 
     &            3*Mgluino**2*Mstop2**2*Mtop**8 - 
     &            81*Mstop2**4*Mtop**8 - 19*Mgluino**2*Mtop**10 + 
     &            9*Mtop**12)*sthetat*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (6.*Mgluino**3*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        + (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*Mtop**2*
     &          (-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &            13*Mgluino**4*Mstop2**4 - 21*Mgluino**2*Mstop2**6 + 
     &            18*Mstop2**8 + 10*Mgluino**6*Mtop**2 - 
     &            19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &            32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &            9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &            71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &            63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &            27*Mstop2**2*Mtop**6 + 9*Mtop**8)*sthetat**2*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (12.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        - (calpha*cthetat*Mtop*
     &          (Mgluino**10*Mstop1**2 + Mgluino**10*Mstop2**2 - 
     &            19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &            7*Mgluino**8*Mstop2**4 + 
     &            42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &            42*Mgluino**6*Mstop2**6 - 
     &            22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &            94*Mgluino**4*Mstop2**8 - 
     &            11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &            85*Mgluino**2*Mstop2**10 + 9*Mstop1**2*Mstop2**10 - 
     &            27*Mstop2**12 - 11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &            11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &            52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &            56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &            46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &            50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &            4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &            40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &            9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &            45*Mstop2**10*Mtop**2 + 
     &            20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &            20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &            92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &            56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &            24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &            120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &            72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &            36*Mstop2**8*Mtop**4 - 
     &            20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &            20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &            28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &            56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &            72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &            108*Mstop2**6*Mtop**6 + 
     &            19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &            19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &            9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &            63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &            9*Mstop2**2*Mtop**10)*sthetat**3*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (6.*Mgluino**3*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        - (calpha*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)*Mtop**2*
     &          (-Mgluino**8 + 17*Mgluino**6*Mstop2**2 - 
     &            13*Mgluino**4*Mstop2**4 - 21*Mgluino**2*Mstop2**6 + 
     &            18*Mstop2**8 + 10*Mgluino**6*Mtop**2 - 
     &            19*Mgluino**4*Mstop2**2*Mtop**2 + 
     &            32*Mgluino**2*Mstop2**4*Mtop**2 + 
     &            9*Mstop2**6*Mtop**2 - 8*Mgluino**4*Mtop**4 + 
     &            71*Mgluino**2*Mstop2**2*Mtop**4 - 
     &            63*Mstop2**4*Mtop**4 - 10*Mgluino**2*Mtop**6 + 
     &            27*Mstop2**2*Mtop**6 + 9*Mtop**8)*sthetat**4*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (12.*Mgluino**2*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2))
     &        + (calpha*cthetat*Mtop*
     &          (Mgluino**10*Mstop1**2 + Mgluino**10*Mstop2**2 - 
     &            19*Mgluino**8*Mstop1**2*Mstop2**2 - 
     &            7*Mgluino**8*Mstop2**4 + 
     &            42*Mgluino**6*Mstop1**2*Mstop2**4 + 
     &            42*Mgluino**6*Mstop2**6 - 
     &            22*Mgluino**4*Mstop1**2*Mstop2**6 - 
     &            94*Mgluino**4*Mstop2**8 - 
     &            11*Mgluino**2*Mstop1**2*Mstop2**8 + 
     &            85*Mgluino**2*Mstop2**10 + 9*Mstop1**2*Mstop2**10 - 
     &            27*Mstop2**12 - 11*Mgluino**8*Mstop1**2*Mtop**2 - 
     &            11*Mgluino**8*Mstop2**2*Mtop**2 + 
     &            52*Mgluino**6*Mstop1**2*Mstop2**2*Mtop**2 + 
     &            56*Mgluino**6*Mstop2**4*Mtop**2 - 
     &            46*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**2 - 
     &            50*Mgluino**4*Mstop2**6*Mtop**2 - 
     &            4*Mgluino**2*Mstop1**2*Mstop2**6*Mtop**2 - 
     &            40*Mgluino**2*Mstop2**8*Mtop**2 + 
     &            9*Mstop1**2*Mstop2**8*Mtop**2 + 
     &            45*Mstop2**10*Mtop**2 + 
     &            20*Mgluino**6*Mstop1**2*Mtop**4 + 
     &            20*Mgluino**6*Mstop2**2*Mtop**4 - 
     &            92*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**4 + 
     &            56*Mgluino**4*Mstop2**4*Mtop**4 + 
     &            24*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**4 - 
     &            120*Mgluino**2*Mstop2**6*Mtop**4 - 
     &            72*Mstop1**2*Mstop2**6*Mtop**4 + 
     &            36*Mstop2**8*Mtop**4 - 
     &            20*Mgluino**4*Mstop1**2*Mtop**6 - 
     &            20*Mgluino**4*Mstop2**2*Mtop**6 - 
     &            28*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**6 + 
     &            56*Mgluino**2*Mstop2**4*Mtop**6 + 
     &            72*Mstop1**2*Mstop2**4*Mtop**6 - 
     &            108*Mstop2**6*Mtop**6 + 
     &            19*Mgluino**2*Mstop1**2*Mtop**8 + 
     &            19*Mgluino**2*Mstop2**2*Mtop**8 - 
     &            9*Mstop1**2*Mstop2**2*Mtop**8 + 
     &            63*Mstop2**4*Mtop**8 - 9*Mstop1**2*Mtop**10 - 
     &            9*Mstop2**2*Mtop**10)*sthetat**5*
     &          davytausk(Mtop**2/Mgluino**2,Mstop2**2/Mgluino**2,
     &           Mgluino))/
     &        (6.*Mgluino**3*(-Mgluino + Mstop2 - Mtop)**3*
     &          (Mgluino + Mstop2 - Mtop)**3*
     &          (-Mgluino + Mstop2 + Mtop)**3*
     &          (Mgluino + Mstop2 + Mtop)**3*
     &          ((-4*Mstop2**2*Mtop**2)/Mgluino**4 + 
     &            (1 - Mstop2**2/Mgluino**2 - Mtop**2/Mgluino**2)**2)))
     &      /sbeta
      end
