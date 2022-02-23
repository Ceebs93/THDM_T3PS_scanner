! This file is part of SusHi.
! 
! It includes relevant references for SusHi, which are automatically listed.
! 
C-{{{ subroutine citeinit:

      subroutine citeinit()
c..   
c..   initialize citations list
c..      
      integer i
      include '../commons/common-citations.f'

      do i=1,100
         citations(i) = 0
      enddo

      end

C-}}}
C-{{{ subroutine printcite:

      subroutine printcite(unit,key,instring)
c..
c..   print texkeys of all papers used
c..   
      integer i,unit,key
      include '../commons/common-citations.f'
      character(*) instring
      character*20 refs,string
      character reflist*100

      do i=1,100
         if (citations(i).eq.key) then
            write(unit,101) instring,refs(i)
         endif
      enddo
      
 101  format(2a)

      end

C-}}}
C-{{{ function refs:

      character*20 function refs(entry)
c..
c..   define citations list
c..
      integer entry
      character*20 references(100)

      do i=1,100
         references(i) = ''
      enddo

      !Inspire identifiers of relevant papers:
      references(1) = 'Harlander:2012pb' !SusHi - always to be cited
      references(2) = 'Harlander:2016hcx' ! SusHi 1.6.0 - always to be cited
      references(3) = 'Liebler:2015bka' ! NMSSM version of SusHi
      references(4) = 'Harlander:2015xur' ! Heavy quark annihilation
      !free numbers - leave free for future publications
      references(7) = 'Harlander:2002wh' !ggh@nnlo - before (2)
      references(8) = 'Harlander:2003ai' !bbh@nnlo
      references(9) = 'Harlander:2010cz' !bbh diff.
      references(10) = 'Actis:2008ug' !SM EW factor
      references(11) = 'Aglietti:2004nj' !EW factor lq
      references(12) = 'Bonciani:2010ms' !EW factor lq
      references(13) = 'Degrassi:2010eu' !SUSY cont.
      references(14) = 'Degrassi:2011vq' !SUSY cont.
      references(15) = 'Degrassi:2012vt' !SUSY cont.
      references(16) = 'Heinemeyer:1998yj' !FeynHiggs
      references(17) = 'Heinemeyer:1998np' !FeynHiggs
      references(18) = 'Degrassi:2002fi' !FeynHiggs
      references(19) = 'Frank:2006yh' !FeynHiggs
      references(20) = 'Djouadi:1991tka' !NLO gluon fusion
      references(21) = 'Spira:1995rr' !NLO gluon fusion
      references(22) = 'Dawson:1990zj' !NLO gluon fusion
      references(23) = 'Anastasiou:2002yz' !NNLO gluon fusion
      references(24) = 'Ravindran:2003um' !NNLO gluon fusion
      references(25) = 'Eriksson:2009ws' !2HDMC
      references(26) = 'Harlander:2005rq' !Harlander/Kant
      references(27) = 'Georgi:1977gs' !gg->H LO
      references(28) = 'Harlander:2009bw' ! NNLO 1/mt virtual
      references(29) = 'Harlander:2009mq' ! NNLO 1/mt no matching
      references(30) = 'Harlander:2009my' ! NNLO 1/mt matching
      references(31) = 'Bechtle:2008jh' ! HiggsBounds
      references(32) = 'Bechtle:2011sb' ! HiggsBounds
      references(33) = 'Bechtle:2013wla' ! HiggsBounds
      references(34) = 'Bechtle:2013xfa' ! HiggsSignals
      !free number
      references(36) = 'Harlander:2002vv' !ggh@nnlo pseudoscalar
      references(37) = 'Anastasiou:2002wq' !ggh@nnlo pseudoscalar
      references(38) = 'Marzani:2008az' ! lim x->0
      references(39) = 'Harlander:2009mq' ! 1/mt @ NNLO
      references(40) = 'Harlander:2009my' ! 1/mt @ NNLO with Marzani
      references(41) = 'Harlander:2009bw' ! 1/mt @ NNLO virtual
      references(42) = 'Pak:2009dg' ! 1/mt by Steinhauser et al.
      references(43) = 'Pak:2009bx' ! 1/mt by Steinhauser et al.
      references(44) = 'Pak:2011hs' ! 1/mt by Steinhauser et al.
      references(45) = 'Anastasiou:2014lda' ! N3LO soft expansion
      references(46) = 'Anastasiou:2015ema' ! N3LO soft expansion - does not contain formulas
      references(47) = 'Anastasiou:2015yha' ! N3LO soft expansion
      references(48) = 'Anastasiou:2016cez' ! N3LO soft expansion
      !free numbers
      references(50) = 'Harlander:2003bb' !evalcsusy - before (3)
      references(51) = 'Harlander:2004tp' !evalcsusy - before (4)
      references(52) = 'Harlander:2003kf' !evalcsusy - before (5)
      references(53) = 'Harlander:2005if' !evalcsusy pseudo - before (6)
      references(54) = 'Harlander:2010wr' !bottom/sbottom - before (7)
      references(55) = 'Chetyrkin:2000yt' !RunDec
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 
c$$$      references() = 

      refs = references(entry)

      end

C-}}}
