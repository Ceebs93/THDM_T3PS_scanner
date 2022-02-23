
      real*8 function ggamix_bare(Mtop,Mstop1,Mstop2,Mgluino,mu,
     & q0,cthetat,sthetat,alpha,beta,
     & muSUSY,c1EWmc2EW,c1EWpc2EW)
c..
c..
      implicit real*8 (a-z)
      data pi/3.14159265358979323846264338328d0/


      L12 = 2.d0*dlog(Mstop1/Mstop2)
      L2gl = 2.d0*dlog(Mstop2/Mgluino)
      Ltgl = 2.d0*dlog(Mtop/Mgluino)
      tanb = dtan(beta)
      ctbeta = 1.d0/tanb

      dtstop1 = 1.d0/Mgluino**2 * davytausk(Mtop**2/Mgluino**2,Mstop1**2
     &     /Mgluino**2,Mgluino)

      dtstop2 = 1.d0/Mgluino**2 * davytausk(Mtop**2/Mgluino**2,Mstop2**2
     &     /Mgluino**2,Mgluino)

      amp2ldeno = (6.*Mgluino*(Mstop1 - Mstop2)*(Mstop1 + Mstop2)
     &     * (Mgluino - Mstop1 - Mtop)**2*(Mgluino + Mstop1 - Mtop)**2
     &     * (Mgluino - Mstop2 - Mtop)**2*(Mgluino + Mstop2 - Mtop)**2
     &     * (Mgluino - Mstop1 + Mtop)**2*(Mgluino + Mstop1 + Mtop)**2
     &     * (Mgluino - Mstop2 + Mtop)**2*(Mgluino + Mstop2 + Mtop)**2)


      amp2lfac = (1 + ctbeta**2)*muSUSY
      
      amp2lnum = 4*L12*Mgluino**16*Mstop1**2 + 4*L2gl*Mgluino**16*Mstop1
     &     **2 - 12*L12*Mgluino**14*Mstop1**4 - 12*L2gl*Mgluino**14
     &     *Mstop1**4 + 12*L12*Mgluino**12*Mstop1**6 + 12*L2gl*Mgluino
     &     **12*Mstop1**6 - 4*L12*Mgluino**10*Mstop1**8 - 4*L2gl*Mgluino
     &     **10*Mstop1**8 - 4*L2gl*Mgluino**16*Mstop2**2 - 16*L12
     &     *Mgluino**14*Mstop1**2*Mstop2**2 + 48*L12*Mgluino**12*Mstop1
     &     **4*Mstop2**2 + 24*L2gl*Mgluino**12*Mstop1**4*Mstop2**2 - 48
     &     *L12*Mgluino**10*Mstop1**6*Mstop2**2 - 32*L2gl*Mgluino**10
     &     *Mstop1**6*Mstop2**2 + 16*L12*Mgluino**8*Mstop1**8*Mstop2**2
     &     + 12*L2gl*Mgluino**8*Mstop1**8*Mstop2**2 + 12*L2gl*Mgluino
     &     **14*Mstop2**4 + 24*L12*Mgluino**12*Mstop1**2*Mstop2**4 - 24
     &     *L2gl*Mgluino**12*Mstop1**2*Mstop2**4 - 72*L12*Mgluino**10
     &     *Mstop1**4*Mstop2**4 

      amp2lnum = amp2lnum + 72*L12*Mgluino**8*Mstop1**6*Mstop2**4 + 24
     &     *L2gl*Mgluino**8*Mstop1**6*Mstop2**4 - 24*L12*Mgluino**6
     &     *Mstop1**8*Mstop2**4 - 12*L2gl*Mgluino**6*Mstop1**8*Mstop2**4
     &     - 12*L2gl*Mgluino**12*Mstop2**6 - 16*L12*Mgluino**10*Mstop1
     &     **2*Mstop2**6 + 32*L2gl*Mgluino**10*Mstop1**2*Mstop2**6 + 48
     &     *L12*Mgluino**8*Mstop1**4*Mstop2**6 - 24*L2gl*Mgluino**8
     &     *Mstop1**4*Mstop2**6 - 48*L12*Mgluino**6*Mstop1**6*Mstop2**6
     &     + 16*L12*Mgluino**4*Mstop1**8*Mstop2**6 + 4*L2gl*Mgluino**4
     &     *Mstop1**8*Mstop2**6 + 4*L2gl*Mgluino**10*Mstop2**8 + 4*L12
     &     *Mgluino**8*Mstop1**2*Mstop2**8 - 12*L2gl*Mgluino**8*Mstop1
     &     **2*Mstop2**8 - 12*L12*Mgluino**6*Mstop1**4*Mstop2**8 + 12
     &     *L2gl*Mgluino**6*Mstop1**4*Mstop2**8 + 12*L12*Mgluino**4
     &     *Mstop1**6*Mstop2**8 - 4*L2gl*Mgluino**4*Mstop1**6*Mstop2**8
     &     - 4*L12*Mgluino**2*Mstop1**8*Mstop2**8 - dtstop1*Mgluino**14
     &     *Mstop1**2*Mtop**2 - 29*L12*Mgluino**14*Mstop1**2*Mtop**2 -
     &     29*L2gl*Mgluino**14*Mstop1**2*Mtop**2

      amp2lnum = amp2lnum + 3*Ltgl*Mgluino**14*Mstop1**2*Mtop**2 + 9
     &     *dtstop1*Mgluino**12*Mstop1**4*Mtop**2 + 49*L12*Mgluino**12
     &     *Mstop1**4*Mtop**2 + 49*L2gl*Mgluino**12*Mstop1**4*Mtop**2 -
     &     Ltgl*Mgluino**12*Mstop1**4*Mtop**2 - 27*L12*Mgluino**10
     &     *Mstop1**6*Mtop**2 - 27*L2gl*Mgluino**10*Mstop1**6*Mtop**2 -
     &     7*Ltgl*Mgluino**10*Mstop1**6*Mtop**2 + 7*L12*Mgluino**8
     &     *Mstop1**8*Mtop**2 + 7*L2gl*Mgluino**8*Mstop1**8*Mtop**2 + 5
     &     *Ltgl*Mgluino**8*Mstop1**8*Mtop**2 + dtstop2*Mgluino**14
     &     *Mstop2**2*Mtop**2 + 29*L2gl*Mgluino**14*Mstop2**2*Mtop**2 -
     &     3*Ltgl*Mgluino**14*Mstop2**2*Mtop**2 + 4*dtstop1*Mgluino**12
     &     *Mstop1**2*Mstop2**2*Mtop**2 - 4*dtstop2*Mgluino**12*Mstop1
     &     **2*Mstop2**2*Mtop**2 + 68*L12*Mgluino**12*Mstop1**2*Mstop2
     &     **2*Mtop**2 - 36*dtstop1*Mgluino**10*Mstop1**4*Mstop2**2*Mtop
     &     **2

      amp2lnum = amp2lnum + 6*dtstop2*Mgluino**10*Mstop1**4*Mstop2**2
     &     *Mtop**2 - 52*L12*Mgluino**10*Mstop1**4*Mstop2**2*Mtop**2 +
     &     10*L2gl*Mgluino**10*Mstop1**4*Mstop2**2*Mtop**2 - 14*Ltgl
     &     *Mgluino**10*Mstop1**4*Mstop2**2*Mtop**2 - 4*dtstop2*Mgluino
     &     **8*Mstop1**6*Mstop2**2*Mtop**2 - 36*L12*Mgluino**8*Mstop1**6
     &     *Mstop2**2*Mtop**2 - 72*L2gl*Mgluino**8*Mstop1**6*Mstop2**2
     &     *Mtop**2 + 40*Ltgl*Mgluino**8*Mstop1**6*Mstop2**2*Mtop**2 +
     &     dtstop2*Mgluino**6*Mstop1**8*Mstop2**2*Mtop**2 + 20*L12
     &     *Mgluino**6*Mstop1**8*Mstop2**2*Mtop**2 + 33*L2gl*Mgluino**6
     &     *Mstop1**8*Mstop2**2*Mtop**2 - 23*Ltgl*Mgluino**6*Mstop1**8
     &     *Mstop2**2*Mtop**2 - 9*dtstop2*Mgluino**12*Mstop2**4*Mtop**2
     &     - 49*L2gl*Mgluino**12*Mstop2**4*Mtop**2 + Ltgl*Mgluino**12
     &     *Mstop2**4*Mtop**2 - 6*dtstop1*Mgluino**10*Mstop1**2*Mstop2
     &     **4*Mtop**2 + 36*dtstop2*Mgluino**10*Mstop1**2*Mstop2**4*Mtop
     &     **2 - 62*L12*Mgluino**10*Mstop1**2*Mstop2**4*Mtop**2 - 10
     &     *L2gl*Mgluino**10*Mstop1**2*Mstop2**4*Mtop**2

      amp2lnum = amp2lnum + 14*Ltgl*Mgluino**10*Mstop1**2*Mstop2**4*Mtop
     &     **2 + 54*dtstop1*Mgluino**8*Mstop1**4*Mstop2**4*Mtop**2 - 54
     &     *dtstop2*Mgluino**8*Mstop1**4*Mstop2**4*Mtop**2 - 42*L12
     &     *Mgluino**8*Mstop1**4*Mstop2**4*Mtop**2 + 36*dtstop2*Mgluino
     &     **6*Mstop1**6*Mstop2**4*Mtop**2 + 174*L12*Mgluino**6*Mstop1
     &     **6*Mstop2**4*Mtop**2 + 130*L2gl*Mgluino**6*Mstop1**6*Mstop2
     &     **4*Mtop**2 - 46*Ltgl*Mgluino**6*Mstop1**6*Mstop2**4*Mtop**2
     &     - 9*dtstop2*Mgluino**4*Mstop1**8*Mstop2**4*Mtop**2 - 70*L12
     &     *Mgluino**4*Mstop1**8*Mstop2**4*Mtop**2 - 71*L2gl*Mgluino**4
     &     *Mstop1**8*Mstop2**4*Mtop**2 + 31*Ltgl*Mgluino**4*Mstop1**8
     &     *Mstop2**4*Mtop**2 + 27*L2gl*Mgluino**10*Mstop2**6*Mtop**2 +
     &     7*Ltgl*Mgluino**10*Mstop2**6*Mtop**2 + 4*dtstop1*Mgluino**8
     &     *Mstop1**2*Mstop2**6*Mtop**2 + 36*L12*Mgluino**8*Mstop1**2
     &     *Mstop2**6*Mtop**2

      amp2lnum = amp2lnum + 72*L2gl*Mgluino**8*Mstop1**2*Mstop2**6*Mtop
     &     **2 - 40*Ltgl*Mgluino**8*Mstop1**2*Mstop2**6*Mtop**2 - 36
     &     *dtstop1*Mgluino**6*Mstop1**4*Mstop2**6*Mtop**2 + 44*L12
     &     *Mgluino**6*Mstop1**4*Mstop2**6*Mtop**2 - 130*L2gl*Mgluino**6
     &     *Mstop1**4*Mstop2**6*Mtop**2 + 46*Ltgl*Mgluino**6*Mstop1**4
     &     *Mstop2**6*Mtop**2 - 132*L12*Mgluino**4*Mstop1**6*Mstop2**6
     &     *Mtop**2 + 52*L12*Mgluino**2*Mstop1**8*Mstop2**6*Mtop**2 + 31
     &     *L2gl*Mgluino**2*Mstop1**8*Mstop2**6*Mtop**2 - 13*Ltgl
     &     *Mgluino**2*Mstop1**8*Mstop2**6*Mtop**2 - 7*L2gl*Mgluino**8
     &     *Mstop2**8*Mtop**2 - 5*Ltgl*Mgluino**8*Mstop2**8*Mtop**2 -
     &     dtstop1*Mgluino**6*Mstop1**2*Mstop2**8*Mtop**2 - 13*L12
     &     *Mgluino**6*Mstop1**2*Mstop2**8*Mtop**2 - 33*L2gl*Mgluino**6
     &     *Mstop1**2*Mstop2**8*Mtop**2 + 23*Ltgl*Mgluino**6*Mstop1**2
     &     *Mstop2**8*Mtop**2 + 9*dtstop1*Mgluino**4*Mstop1**4*Mstop2**8
     &     *Mtop**2 + L12*Mgluino**4*Mstop1**4*Mstop2**8*Mtop**2 + 71
     &     *L2gl*Mgluino**4*Mstop1**4*Mstop2**8*Mtop**2 - 31*Ltgl
     &     *Mgluino**4*Mstop1**4*Mstop2**8*Mtop**2 + 21*L12*Mgluino**2
     &     *Mstop1**6*Mstop2**8*Mtop**2 - 31*L2gl*Mgluino**2*Mstop1**6
     &     *Mstop2**8*Mtop**2

      amp2lnum = amp2lnum + 13*Ltgl*Mgluino**2*Mstop1**6*Mstop2**8*Mtop
     &     **2 - 9*L12*Mstop1**8*Mstop2**8*Mtop**2 + 13*dtstop1*Mgluino
     &     **12*Mstop1**2*Mtop**4 + 99*L12*Mgluino**12*Mstop1**2*Mtop**4
     &     + 99*L2gl*Mgluino**12*Mstop1**2*Mtop**4 - 11*Ltgl*Mgluino**12
     &     *Mstop1**2*Mtop**4 - 36*dtstop1*Mgluino**10*Mstop1**4*Mtop**4
     &     - 70*L12*Mgluino**10*Mstop1**4*Mtop**4 - 70*L2gl*Mgluino**10
     &     *Mstop1**4*Mtop**4 - 30*Ltgl*Mgluino**10*Mstop1**4*Mtop**4 +
     &     15*L12*Mgluino**8*Mstop1**6*Mtop**4 + 15*L2gl*Mgluino**8
     &     *Mstop1**6*Mtop**4 + 33*Ltgl*Mgluino**8*Mstop1**6*Mtop**4 +
     &     12*L12*Mgluino**6*Mstop1**8*Mtop**4 + 12*L2gl*Mgluino**6
     &     *Mstop1**8*Mtop**4 - 24*Ltgl*Mgluino**6*Mstop1**8*Mtop**4 -
     &     13*dtstop2*Mgluino**12*Mstop2**2*Mtop**4 - 99*L2gl*Mgluino
     &     **12*Mstop2**2*Mtop**4

      amp2lnum = amp2lnum + 11*Ltgl*Mgluino**12*Mstop2**2*Mtop**4 - 40
     &     *dtstop1*Mgluino**10*Mstop1**2*Mstop2**2*Mtop**4 + 40*dtstop2
     &     *Mgluino**10*Mstop1**2*Mstop2**2*Mtop**4 - 128*L12*Mgluino
     &     **10*Mstop1**2*Mstop2**2*Mtop**4 + 36*dtstop1*Mgluino**8
     &     *Mstop1**4*Mstop2**2*Mtop**4 - 50*dtstop2*Mgluino**8*Mstop1
     &     **4*Mstop2**2*Mtop**4 - 68*L12*Mgluino**8*Mstop1**4*Mstop2**2
     &     *Mtop**4 - 178*L2gl*Mgluino**8*Mstop1**4*Mstop2**2*Mtop**4 +
     &     114*Ltgl*Mgluino**8*Mstop1**4*Mstop2**2*Mtop**4 + 32*dtstop2
     &     *Mgluino**6*Mstop1**6*Mstop2**2*Mtop**4 + 24*L12*Mgluino**6
     &     *Mstop1**6*Mstop2**2*Mtop**4 + 64*L2gl*Mgluino**6*Mstop1**6
     &     *Mstop2**2*Mtop**4 - 32*Ltgl*Mgluino**6*Mstop1**6*Mstop2**2
     &     *Mtop**4 - 9*dtstop2*Mgluino**4*Mstop1**8*Mstop2**2*Mtop**4 -
     &     52*L12*Mgluino**4*Mstop1**8*Mstop2**2*Mtop**4 - 75*L2gl
     &     *Mgluino**4*Mstop1**8*Mstop2**2*Mtop**4 + 35*Ltgl*Mgluino**4
     &     *Mstop1**8*Mstop2**2*Mtop**4 + 36*dtstop2*Mgluino**10*Mstop2
     &     **4*Mtop**4 + 70*L2gl*Mgluino**10*Mstop2**4*Mtop**4 + 30*Ltgl
     &     *Mgluino**10*Mstop2**4*Mtop**4

      amp2lnum = amp2lnum + 50*dtstop1*Mgluino**8*Mstop1**2*Mstop2**4
     &     *Mtop**4 - 36*dtstop2*Mgluino**8*Mstop1**2*Mstop2**4*Mtop**4
     &     + 110*L12*Mgluino**8*Mstop1**2*Mstop2**4*Mtop**4 + 178*L2gl
     &     *Mgluino**8*Mstop1**2*Mstop2**4*Mtop**4 - 114*Ltgl*Mgluino**8
     &     *Mstop1**2*Mstop2**4*Mtop**4 + 36*dtstop1*Mgluino**6*Mstop1
     &     **4*Mstop2**4*Mtop**4 - 36*dtstop2*Mgluino**6*Mstop1**4
     &     *Mstop2**4*Mtop**4 - 32*L12*Mgluino**6*Mstop1**4*Mstop2**4
     &     *Mtop**4 + 36*dtstop2*Mgluino**4*Mstop1**6*Mstop2**4*Mtop**4
     &     + 318*L12*Mgluino**4*Mstop1**6*Mstop2**4*Mtop**4 + 346*L2gl
     &     *Mgluino**4*Mstop1**6*Mstop2**4*Mtop**4 - 138*Ltgl*Mgluino**4
     &     *Mstop1**6*Mstop2**4*Mtop**4 - 60*L12*Mgluino**2*Mstop1**8
     &     *Mstop2**4*Mtop**4 - 66*L2gl*Mgluino**2*Mstop1**8*Mstop2**4
     &     *Mtop**4 + 30*Ltgl*Mgluino**2*Mstop1**8*Mstop2**4*Mtop**4 -
     &     15*L2gl*Mgluino**8*Mstop2**6*Mtop**4 - 33*Ltgl*Mgluino**8
     &     *Mstop2**6*Mtop**4 - 32*dtstop1*Mgluino**6*Mstop1**2*Mstop2
     &     **6*Mtop**4 - 40*L12*Mgluino**6*Mstop1**2*Mstop2**6*Mtop**4 -
     &     64*L2gl*Mgluino**6*Mstop1**2*Mstop2**6*Mtop**4

      amp2lnum = amp2lnum + 32*Ltgl*Mgluino**6*Mstop1**2*Mstop2**6*Mtop
     &     **4 - 36*dtstop1*Mgluino**4*Mstop1**4*Mstop2**6*Mtop**4 - 28
     &     *L12*Mgluino**4*Mstop1**4*Mstop2**6*Mtop**4 - 346*L2gl
     &     *Mgluino**4*Mstop1**4*Mstop2**6*Mtop**4 + 138*Ltgl*Mgluino**4
     &     *Mstop1**4*Mstop2**6*Mtop**4 - 192*L12*Mgluino**2*Mstop1**6
     &     *Mstop2**6*Mtop**4 + 36*L12*Mstop1**8*Mstop2**6*Mtop**4 + 9
     &     *L2gl*Mstop1**8*Mstop2**6*Mtop**4 - 9*Ltgl*Mstop1**8*Mstop2
     &     **6*Mtop**4 - 12*L2gl*Mgluino**6*Mstop2**8*Mtop**4 + 24*Ltgl
     &     *Mgluino**6*Mstop2**8*Mtop**4 + 9*dtstop1*Mgluino**4*Mstop1
     &     **2*Mstop2**8*Mtop**4 + 23*L12*Mgluino**4*Mstop1**2*Mstop2**8
     &     *Mtop**4 + 75*L2gl*Mgluino**4*Mstop1**2*Mstop2**8*Mtop**4 -
     &     35*Ltgl*Mgluino**4*Mstop1**2*Mstop2**8*Mtop**4 + 6*L12
     &     *Mgluino**2*Mstop1**4*Mstop2**8*Mtop**4 + 66*L2gl*Mgluino**2
     &     *Mstop1**4*Mstop2**8*Mtop**4 - 30*Ltgl*Mgluino**2*Mstop1**4
     &     *Mstop2**8*Mtop**4 + 27*L12*Mstop1**6*Mstop2**8*Mtop**4 - 9
     &     *L2gl*Mstop1**6*Mstop2**8*Mtop**4 + 9*Ltgl*Mstop1**6*Mstop2
     &     **8*Mtop**4 - 42*dtstop1*Mgluino**10*Mstop1**2*Mtop**6 - 209
     &     *L12*Mgluino**10*Mstop1**2*Mtop**6 - 209*L2gl*Mgluino**10
     &     *Mstop1**2*Mtop**6

      amp2lnum = amp2lnum + 19*Ltgl*Mgluino**10*Mstop1**2*Mtop**6 + 54
     &     *dtstop1*Mgluino**8*Mstop1**4*Mtop**6 + 3*L12*Mgluino**8
     &     *Mstop1**4*Mtop**6 + 3*L2gl*Mgluino**8*Mstop1**4*Mtop**6 + 93
     &     *Ltgl*Mgluino**8*Mstop1**4*Mtop**6 - 30*L12*Mgluino**6*Mstop1
     &     **6*Mtop**6 - 30*L2gl*Mgluino**6*Mstop1**6*Mtop**6 - 30*Ltgl
     &     *Mgluino**6*Mstop1**6*Mtop**6 - 38*L12*Mgluino**4*Mstop1**8
     &     *Mtop**6 - 38*L2gl*Mgluino**4*Mstop1**8*Mtop**6 + 42*Ltgl
     &     *Mgluino**4*Mstop1**8*Mtop**6 + 42*dtstop2*Mgluino**10*Mstop2
     &     **2*Mtop**6 + 209*L2gl*Mgluino**10*Mstop2**2*Mtop**6 - 19
     &     *Ltgl*Mgluino**10*Mstop2**2*Mtop**6 + 32*dtstop1*Mgluino**8
     &     *Mstop1**2*Mstop2**2*Mtop**6 - 32*dtstop2*Mgluino**8*Mstop1
     &     **2*Mstop2**2*Mtop**6 + 116*L12*Mgluino**8*Mstop1**2*Mstop2
     &     **2*Mtop**6 + 36*dtstop1*Mgluino**6*Mstop1**4*Mstop2**2*Mtop
     &     **6 - 30*dtstop2*Mgluino**6*Mstop1**4*Mstop2**2*Mtop**6 + 184
     &     *L12*Mgluino**6*Mstop1**4*Mstop2**2*Mtop**6

      amp2lnum = amp2lnum + 308*L2gl*Mgluino**6*Mstop1**4*Mstop2**2*Mtop
     &     **6 - 44*Ltgl*Mgluino**6*Mstop1**4*Mstop2**2*Mtop**6 + 36
     &     *dtstop2*Mgluino**4*Mstop1**6*Mstop2**2*Mtop**6 + 144*L12
     &     *Mgluino**4*Mstop1**6*Mstop2**2*Mtop**6 + 144*L2gl*Mgluino**4
     &     *Mstop1**6*Mstop2**2*Mtop**6 - 112*Ltgl*Mgluino**4*Mstop1**6
     &     *Mstop2**2*Mtop**6 - 20*L12*Mgluino**2*Mstop1**8*Mstop2**2
     &     *Mtop**6 + 3*L2gl*Mgluino**2*Mstop1**8*Mstop2**2*Mtop**6 + 15
     &     *Ltgl*Mgluino**2*Mstop1**8*Mstop2**2*Mtop**6 - 54*dtstop2
     &     *Mgluino**8*Mstop2**4*Mtop**6 - 3*L2gl*Mgluino**8*Mstop2**4
     &     *Mtop**6 - 93*Ltgl*Mgluino**8*Mstop2**4*Mtop**6 + 30*dtstop1
     &     *Mgluino**6*Mstop1**2*Mstop2**4*Mtop**6 - 36*dtstop2*Mgluino
     &     **6*Mstop1**2*Mstop2**4*Mtop**6 - 124*L12*Mgluino**6*Mstop1
     &     **2*Mstop2**4*Mtop**6 - 308*L2gl*Mgluino**6*Mstop1**2*Mstop2
     &     **4*Mtop**6 + 44*Ltgl*Mgluino**6*Mstop1**2*Mstop2**4*Mtop**6
     &     + 54*dtstop1*Mgluino**4*Mstop1**4*Mstop2**4*Mtop**6 - 54
     &     *dtstop2*Mgluino**4*Mstop1**4*Mstop2**4*Mtop**6 - 132*L12
     &     *Mgluino**4*Mstop1**4*Mstop2**4*Mtop**6 + 234*L12*Mgluino**2
     &     *Mstop1**6*Mstop2**4*Mtop**6 + 150*L2gl*Mgluino**2*Mstop1**6
     &     *Mstop2**4*Mtop**6 - 114*Ltgl*Mgluino**2*Mstop1**6*Mstop2**4
     &     *Mtop**6 - 54*L12*Mstop1**8*Mstop2**4*Mtop**6 - 27*L2gl
     &     *Mstop1**8*Mstop2**4*Mtop**6 + 27*Ltgl*Mstop1**8*Mstop2**4
     &     *Mtop**6 + 30*L2gl*Mgluino**6*Mstop2**6*Mtop**6

      amp2lnum = amp2lnum + 30*Ltgl*Mgluino**6*Mstop2**6*Mtop**6 - 36
     &     *dtstop1*Mgluino**4*Mstop1**2*Mstop2**6*Mtop**6 - 144*L2gl
     &     *Mgluino**4*Mstop1**2*Mstop2**6*Mtop**6 + 112*Ltgl*Mgluino**4
     &     *Mstop1**2*Mstop2**6*Mtop**6 + 84*L12*Mgluino**2*Mstop1**4
     &     *Mstop2**6*Mtop**6 - 150*L2gl*Mgluino**2*Mstop1**4*Mstop2**6
     &     *Mtop**6 + 114*Ltgl*Mgluino**2*Mstop1**4*Mstop2**6*Mtop**6 -
     &     108*L12*Mstop1**6*Mstop2**6*Mtop**6 + 38*L2gl*Mgluino**4
     &     *Mstop2**8*Mtop**6 - 42*Ltgl*Mgluino**4*Mstop2**8*Mtop**6 -
     &     23*L12*Mgluino**2*Mstop1**2*Mstop2**8*Mtop**6 - 3*L2gl
     &     *Mgluino**2*Mstop1**2*Mstop2**8*Mtop**6 - 15*Ltgl*Mgluino**2
     &     *Mstop1**2*Mstop2**8*Mtop**6 - 27*L12*Mstop1**4*Mstop2**8
     &     *Mtop**6 + 27*L2gl*Mstop1**4*Mstop2**8*Mtop**6 - 27*Ltgl
     &     *Mstop1**4*Mstop2**8*Mtop**6 + 58*dtstop1*Mgluino**8*Mstop1
     &     **2*Mtop**8 + 295*L12*Mgluino**8*Mstop1**2*Mtop**8 + 295*L2gl
     &     *Mgluino**8*Mstop1**2*Mtop**8 - 35*Ltgl*Mgluino**8*Mstop1**2
     &     *Mtop**8 - 36*dtstop1*Mgluino**6*Mstop1**4*Mtop**8 + 128*L12
     &     *Mgluino**6*Mstop1**4*Mtop**8 + 128*L2gl*Mgluino**6*Mstop1**4
     &     *Mtop**8 - 116*Ltgl*Mgluino**6*Mstop1**4*Mtop**8 + 90*L12
     &     *Mgluino**4*Mstop1**6*Mtop**8

      amp2lnum = amp2lnum + 90*L2gl*Mgluino**4*Mstop1**6*Mtop**8 - 38
     &     *Ltgl*Mgluino**4*Mstop1**6*Mtop**8 + 32*L12*Mgluino**2*Mstop1
     &     **8*Mtop**8 + 32*L2gl*Mgluino**2*Mstop1**8*Mtop**8 - 32*Ltgl
     &     *Mgluino**2*Mstop1**8*Mtop**8 - 58*dtstop2*Mgluino**8*Mstop2
     &     **2*Mtop**8 - 295*L2gl*Mgluino**8*Mstop2**2*Mtop**8 + 35*Ltgl
     &     *Mgluino**8*Mstop2**2*Mtop**8 + 40*dtstop1*Mgluino**6*Mstop1
     &     **2*Mstop2**2*Mtop**8 - 40*dtstop2*Mgluino**6*Mstop1**2
     &     *Mstop2**2*Mtop**8 + 16*L12*Mgluino**6*Mstop1**2*Mstop2**2
     &     *Mtop**8 - 36*dtstop1*Mgluino**4*Mstop1**4*Mstop2**2*Mtop**8
     &     - 54*dtstop2*Mgluino**4*Mstop1**4*Mstop2**2*Mtop**8 - 88*L12
     &     *Mgluino**4*Mstop1**4*Mstop2**2*Mtop**8 - 188*L2gl*Mgluino**4
     &     *Mstop1**4*Mstop2**2*Mtop**8 + 4*Ltgl*Mgluino**4*Mstop1**4
     &     *Mstop2**2*Mtop**8 + 24*L12*Mgluino**2*Mstop1**6*Mstop2**2
     &     *Mtop**8 - 32*L2gl*Mgluino**2*Mstop1**6*Mstop2**2*Mtop**8 +
     &     32*Ltgl*Mgluino**2*Mstop1**6*Mstop2**2*Mtop**8 + 36*L12
     &     *Mstop1**8*Mstop2**2*Mtop**8 + 27*L2gl*Mstop1**8*Mstop2**2
     &     *Mtop**8 - 27*Ltgl*Mstop1**8*Mstop2**2*Mtop**8 + 36*dtstop2
     &     *Mgluino**6*Mstop2**4*Mtop**8 - 128*L2gl*Mgluino**6*Mstop2**4
     &     *Mtop**8 + 116*Ltgl*Mgluino**6*Mstop2**4*Mtop**8 + 54*dtstop1
     &     *Mgluino**4*Mstop1**2*Mstop2**4*Mtop**8 + 36*dtstop2*Mgluino
     &     **4*Mstop1**2*Mstop2**4*Mtop**8

      amp2lnum = amp2lnum + 100*L12*Mgluino**4*Mstop1**2*Mstop2**4*Mtop
     &     **8 + 188*L2gl*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**8 - 4
     &     *Ltgl*Mgluino**4*Mstop1**2*Mstop2**4*Mtop**8 - 72*L12*Mgluino
     &     **2*Mstop1**4*Mstop2**4*Mtop**8 + 162*L12*Mstop1**6*Mstop2**4
     &     *Mtop**8 + 54*L2gl*Mstop1**6*Mstop2**4*Mtop**8 - 54*Ltgl
     &     *Mstop1**6*Mstop2**4*Mtop**8 - 90*L2gl*Mgluino**4*Mstop2**6
     &     *Mtop**8 + 38*Ltgl*Mgluino**4*Mstop2**6*Mtop**8 + 56*L12
     &     *Mgluino**2*Mstop1**2*Mstop2**6*Mtop**8 + 32*L2gl*Mgluino**2
     &     *Mstop1**2*Mstop2**6*Mtop**8 - 32*Ltgl*Mgluino**2*Mstop1**2
     &     *Mstop2**6*Mtop**8 + 108*L12*Mstop1**4*Mstop2**6*Mtop**8 - 54
     &     *L2gl*Mstop1**4*Mstop2**6*Mtop**8 + 54*Ltgl*Mstop1**4*Mstop2
     &     **6*Mtop**8 - 32*L2gl*Mgluino**2*Mstop2**8*Mtop**8 + 32*Ltgl
     &     *Mgluino**2*Mstop2**8*Mtop**8 + 9*L12*Mstop1**2*Mstop2**8
     &     *Mtop**8 - 27*L2gl*Mstop1**2*Mstop2**8*Mtop**8 + 27*Ltgl
     &     *Mstop1**2*Mstop2**8*Mtop**8 - 37*dtstop1*Mgluino**6*Mstop1
     &     **2*Mtop**10 - 279*L12*Mgluino**6*Mstop1**2*Mtop**10 - 279
     &     *L2gl*Mgluino**6*Mstop1**2*Mtop**10 + 65*Ltgl*Mgluino**6
     &     *Mstop1**2*Mtop**10 + 9*dtstop1*Mgluino**4*Mstop1**4*Mtop**10
     &     - 185*L12*Mgluino**4*Mstop1**4*Mtop**10 - 185*L2gl*Mgluino**4
     &     *Mstop1**4*Mtop**10 + 105*Ltgl*Mgluino**4*Mstop1**4*Mtop**10
     &     - 87*L12*Mgluino**2*Mstop1**6*Mtop**10 - 87*L2gl*Mgluino**2
     &     *Mstop1**6*Mtop**10 + 69*Ltgl*Mgluino**2*Mstop1**6*Mtop**10 -
     &     9*L12*Mstop1**8*Mtop**10 - 9*L2gl*Mstop1**8*Mtop**10

      amp2lnum = amp2lnum + 9*Ltgl*Mstop1**8*Mtop**10 + 37*dtstop2
     &     *Mgluino**6*Mstop2**2*Mtop**10 + 279*L2gl*Mgluino**6*Mstop2
     &     **2*Mtop**10 - 65*Ltgl*Mgluino**6*Mstop2**2*Mtop**10 - 36
     &     *dtstop1*Mgluino**4*Mstop1**2*Mstop2**2*Mtop**10 + 36*dtstop2
     &     *Mgluino**4*Mstop1**2*Mstop2**2*Mtop**10 - 148*L12*Mgluino**4
     &     *Mstop1**2*Mstop2**2*Mtop**10 - 132*L12*Mgluino**2*Mstop1**4
     &     *Mstop2**2*Mtop**10 - 30*L2gl*Mgluino**2*Mstop1**4*Mstop2**2
     &     *Mtop**10 - 6*Ltgl*Mgluino**2*Mstop1**4*Mstop2**2*Mtop**10 -
     &     108*L12*Mstop1**6*Mstop2**2*Mtop**10 - 72*L2gl*Mstop1**6
     &     *Mstop2**2*Mtop**10 + 72*Ltgl*Mstop1**6*Mstop2**2*Mtop**10 -
     &     9*dtstop2*Mgluino**4*Mstop2**4*Mtop**10 + 185*L2gl*Mgluino**4
     &     *Mstop2**4*Mtop**10 - 105*Ltgl*Mgluino**4*Mstop2**4*Mtop**10
     &     - 102*L12*Mgluino**2*Mstop1**2*Mstop2**4*Mtop**10 + 30*L2gl
     &     *Mgluino**2*Mstop1**2*Mstop2**4*Mtop**10 + 6*Ltgl*Mgluino**2
     &     *Mstop1**2*Mstop2**4*Mtop**10 - 162*L12*Mstop1**4*Mstop2**4
     &     *Mtop**10 + 87*L2gl*Mgluino**2*Mstop2**6*Mtop**10 - 69*Ltgl
     &     *Mgluino**2*Mstop2**6*Mtop**10 - 36*L12*Mstop1**2*Mstop2**6
     &     *Mtop**10 + 72*L2gl*Mstop1**2*Mstop2**6*Mtop**10 - 72*Ltgl
     &     *Mstop1**2*Mstop2**6*Mtop**10 + 9*L2gl*Mstop2**8*Mtop**10 - 9
     &     *Ltgl*Mstop2**8*Mtop**10 + 9*dtstop1*Mgluino**4*Mstop1**2
     &     *Mtop**12 + 169*L12*Mgluino**4*Mstop1**2*Mtop**12 + 169*L2gl
     &     *Mgluino**4*Mstop1**2*Mtop**12 - 73*Ltgl*Mgluino**4*Mstop1**2
     &     *Mtop**12 + 114*L12*Mgluino**2*Mstop1**4*Mtop**12 + 114*L2gl
     &     *Mgluino**2*Mstop1**4*Mtop**12 - 78*Ltgl*Mgluino**2*Mstop1**4
     &     *Mtop**12 + 27*L12*Mstop1**6*Mtop**12 + 27*L2gl*Mstop1**6
     &     *Mtop**12 - 27*Ltgl*Mstop1**6*Mtop**12 - 9*dtstop2*Mgluino**4
     &     *Mstop2**2*Mtop**12 - 169*L2gl*Mgluino**4*Mstop2**2*Mtop**12

      amp2lnum = amp2lnum + 73*Ltgl*Mgluino**4*Mstop2**2*Mtop**12 + 128
     &     *L12*Mgluino**2*Mstop1**2*Mstop2**2*Mtop**12 + 108*L12*Mstop1
     &     **4*Mstop2**2*Mtop**12 + 54*L2gl*Mstop1**4*Mstop2**2*Mtop**12
     &     - 54*Ltgl*Mstop1**4*Mstop2**2*Mtop**12 - 114*L2gl*Mgluino**2
     &     *Mstop2**4*Mtop**12 + 78*Ltgl*Mgluino**2*Mstop2**4*Mtop**12 +
     &     54*L12*Mstop1**2*Mstop2**4*Mtop**12 - 54*L2gl*Mstop1**2
     &     *Mstop2**4*Mtop**12 + 54*Ltgl*Mstop1**2*Mstop2**4*Mtop**12 -
     &     27*L2gl*Mstop2**6*Mtop**12 + 27*Ltgl*Mstop2**6*Mtop**12 - 59
     &     *L12*Mgluino**2*Mstop1**2*Mtop**14 - 59*L2gl*Mgluino**2
     &     *Mstop1**2*Mtop**14 + 41*Ltgl*Mgluino**2*Mstop1**2*Mtop**14 -
     &     27*L12*Mstop1**4*Mtop**14 - 27*L2gl*Mstop1**4*Mtop**14 + 27
     &     *Ltgl*Mstop1**4*Mtop**14 + 59*L2gl*Mgluino**2*Mstop2**2*Mtop
     &     **14 - 41*Ltgl*Mgluino**2*Mstop2**2*Mtop**14 - 36*L12*Mstop1
     &     **2*Mstop2**2*Mtop**14 + 27*L2gl*Mstop2**4*Mtop**14 - 27*Ltgl
     &     *Mstop2**4*Mtop**14 + 9*L12*Mstop1**2*Mtop**16 + 9*L2gl
     &     *Mstop1**2*Mtop**16 - 9*Ltgl*Mstop1**2*Mtop**16 - 9*L2gl
     &     *Mstop2**2*Mtop**16 + 9*Ltgl*Mstop2**2*Mtop**16


c      print*,amp2lfac,amp2lnum,amp2ldeno

      ggamix_bare = amp2lfac * amp2lnum / amp2ldeno


      return
      end
