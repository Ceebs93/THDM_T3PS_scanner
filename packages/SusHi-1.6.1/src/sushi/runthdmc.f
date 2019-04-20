! This file is part of SusHi.
! 
! Interface to 2HDMC
! 

#include "versions.h"

      subroutine runthdmc()

      implicit none
      integer slha,i,mapintdef
      real*8 tanbtest
c$$$      real*8 thdmlam(7)!,thdmcres(15)
      real*8 dumpar(10)
      real*8 m122
      real*8 mapreal8
      real*8 smparam(6)
      character fhvers*20,hbvers*20,hsvers*20,thdmcvers*20
      common /vers/ fhvers,hbvers,hsvers,thdmcvers
      integer thdmc_set_param
      external thdmc_set_param
      external mapreal8
      include '../commons/common-inputoutput.f'
      include '../commons/common-slha.f'
      include '../commons/common-ren.f'
      include '../commons/common-vars.f'
      include '../commons/common-quark.f'

      thdmcvers = THDMCVERSION

      do i=1,10
         dumpar(i) = 0.d0
      enddo

      smparam(1) = SUinvAlfa
      smparam(2) = gf
      smparam(3) = AlfasMZ
      smparam(4) = mz
      smparam(5) = mbmb
      smparam(6) = mt

      thdmckey = mapreal8(thdmcpar(1),lthdmcpar(1),
     &     '2HDMC key undefined')
      twohdmver = mapreal8(thdmcpar(2),lthdmcpar(2),
     &     '2HDM version undefined')
      thdmcbreak = mapintdef(int(thdmcpar(10)),lthdmcpar(10),0
     &     ,'2HDM theory validity checks will be ignored')

      do i=1,15
         tparam(i) = 0.d0
      enddo

      if (thdmckey.eq.1) then
         tanb = mapreal8(thdmcpar(3),lthdmcpar(3)
     &        ,'tanb for 2HDMC undefined.')
         tparam(1) = mapreal8(thdmcpar(11),lthdmcpar(11),
     &        'lambda1 for 2HDMC undefined.')
         tparam(2) = mapreal8(thdmcpar(12),lthdmcpar(12),
     &        'lambda2 for 2HDMC undefined.')
         tparam(3) = mapreal8(thdmcpar(13),lthdmcpar(13),
     &        'lambda3 for 2HDMC undefined.')
         tparam(4) = mapreal8(thdmcpar(14),lthdmcpar(14),
     &        'lambda4 for 2HDMC undefined.')
         tparam(5) = mapreal8(thdmcpar(15),lthdmcpar(15),
     &        'lambda5 for 2HDMC undefined.')
         tparam(6) = mapreal8(thdmcpar(16),lthdmcpar(16),
     &        'lambda6 for 2HDMC undefined.')
         tparam(7) = mapreal8(thdmcpar(17),lthdmcpar(17),
     &        'lambda7 for 2HDMC undefined.')
         tparam(8) =  mapreal8(thdmcpar(4),lthdmcpar(4)
     &        ,'m12_2 for 2HDMC undefined.')
      tparam(9) = tanb
      tparam(10) = twohdmver
c Warning/WARNING/Comment/comment
c Why on earth do you made m12 as an input instead of m12_2???
c What a piece of shit.
c        tanb = mapreal8(thdmcpar(3),lthdmcpar(3)
c    &        ,'tanb for 2HDMC undefined.')
c        thdmcm12 = mapreal8(thdmcpar(4),lthdmcpar(4)
c    &        ,'m12 for 2HDMC undefined.')
c        tparam(1) = mapreal8(thdmcpar(11),lthdmcpar(11),
c    &        'lambda1 for 2HDMC undefined.')
c        tparam(2) = mapreal8(thdmcpar(12),lthdmcpar(12),
c    &        'lambda2 for 2HDMC undefined.')
c        tparam(3) = mapreal8(thdmcpar(13),lthdmcpar(13),
c    &        'lambda3 for 2HDMC undefined.')
c        tparam(4) = mapreal8(thdmcpar(14),lthdmcpar(14),
c    &        'lambda4 for 2HDMC undefined.')
c        tparam(5) = mapreal8(thdmcpar(15),lthdmcpar(15),
c    &        'lambda5 for 2HDMC undefined.')
c        tparam(6) = mapreal8(thdmcpar(16),lthdmcpar(16),
c    &        'lambda6 for 2HDMC undefined.')
c        tparam(7) = mapreal8(thdmcpar(17),lthdmcpar(17),
c    &        'lambda7 for 2HDMC undefined.')
c     tparam(8) = thdmcm12*thdmcm12
c     tparam(9) = tanb
c     tparam(10) = twohdmver
      else if (thdmckey.eq.2) then
         tanb = mapreal8(thdmcpar(3),lthdmcpar(3)
     &        ,'tanb for 2HDMC undefined.')
         thdmcm12 = mapreal8(thdmcpar(4),lthdmcpar(4)
     &        ,'m12 for 2HDMC undefined.')
         tparam(1) = mapreal8(thdmcpar(21),lthdmcpar(21),
     &        'mh for 2HDMC undefined.')
         tparam(2) = mapreal8(thdmcpar(22),lthdmcpar(22),
     &        'mH for 2HDMC undefined.')
         tparam(3) = mapreal8(thdmcpar(23),lthdmcpar(23),
     &        'mA for 2HDMC undefined.')
         tparam(4) = mapreal8(thdmcpar(24),lthdmcpar(24),
     &        'mC for 2HDMC undefined.')
         tparam(5) = mapreal8(thdmcpar(25),lthdmcpar(25),
     &        'sin(beta-alpha) for 2HDMC undefined')
         tparam(6) = mapreal8(thdmcpar(26),lthdmcpar(26),
     &        'lambda6 for 2HDMC undefined.')
         tparam(7) = mapreal8(thdmcpar(27),lthdmcpar(27),
     &        'lambda7 for 2HDMC undefined.')
      tparam(8) = thdmcm12*thdmcm12
      tparam(9) = tanb
      tparam(10) = twohdmver
      else if (thdmckey.eq.3) then
         tanb = mapreal8(thdmcpar(3),lthdmcpar(3)
     &        ,'tanb for 2HDMC undefined.')
         tparam(1) = mapreal8(thdmcpar(31),lthdmcpar(31),
     &        'mh for 2HDMC undefined.')
         tparam(2) = mapreal8(thdmcpar(32),lthdmcpar(32),
     &        'mH for 2HDMC undefined.')
         tparam(3) = mapreal8(thdmcpar(33),lthdmcpar(33),
     &        'sin(beta-alpha) for 2HDMC undefined.')
         tparam(4) = tanb
         tparam(5) = mapreal8(thdmcpar(34),lthdmcpar(34),
     &        'Z4 for 2HDMC undefined.')
         tparam(6) = mapreal8(thdmcpar(35),lthdmcpar(35),
     &        'Z5 for 2HDMC undefined.')
         tparam(7) = mapreal8(thdmcpar(36),lthdmcpar(36),
     &        'Z7 for 2HDMC undefined.')
      tparam(10) = twohdmver
      else if (thdmckey.eq.4) then
         tanb = mapreal8(thdmcpar(3),lthdmcpar(3)
     &        ,'tanb for 2HDMC undefined.')
         tparam(1) = mapreal8(thdmcpar(31),lthdmcpar(31),
     &        'mh for 2HDMC undefined.')
         tparam(2) = mapreal8(thdmcpar(32),lthdmcpar(32),
     &        'mH for 2HDMC undefined.')
         tparam(3) = mapreal8(thdmcpar(33),lthdmcpar(33),
     &        'cos(beta-alpha) for 2HDMC undefined.')
         tparam(4) = tanb
         tparam(5) = mapreal8(thdmcpar(34),lthdmcpar(34),
     &        'Z4 for 2HDMC undefined.')
         tparam(6) = mapreal8(thdmcpar(35),lthdmcpar(35),
     &        'Z5 for 2HDMC undefined.')
         tparam(7) = mapreal8(thdmcpar(36),lthdmcpar(36),
     &        'Z7 for 2HDMC undefined.')
      tparam(10) = twohdmver
      else
         call printdie('thdmckey not implemented')
      endif

      slha=1

      if (thdmc_set_param(thdmckey,smparam,tparam
     &     ,thdmcres,slha).ne.0) then
         call printdie('2HDMC run time error.')
      endif

      write(6,*) 'Stability:      ',int(thdmcres(7))
      write(6,*) 'Perturbativity: ',int(thdmcres(8))
      write(6,*) 'Unitarity:      ',int(thdmcres(9))
      if ((thdmcbreak.eq.1).and.thdmcres(7)*thdmcres(8)
     &     *thdmcres(9).eq.0.d0) then
         call printinfosushi(6,'2HDM theoretically unviable.')
         call printinfosushi(6
     &        ,'To ignore this, set Block 2HDMC, entry 10 to 0.')
         call printdie('')
      endif

      if (pseudo.eq.0) then
       if (Sind.eq.1) then
         mh = thdmcres(1)
       else if (Sind.eq.2) then
         mh = thdmcres(2)
       end if
      elseif (pseudo.eq.1) then
         mh = thdmcres(3)
      endif
      
      alpha = thdmcres(4)
      beta = datan(thdmcres(5))

      Hmx = 0.d0
      Hmx(1,1) = -dsin(alpha)
      Hmx(1,2) = dcos(alpha)
      Hmx(2,1) = dcos(alpha)
      Hmx(2,2) = dsin(alpha)

      Amx = 0.d0
      Amx(2,2) = 1.d0

      end
