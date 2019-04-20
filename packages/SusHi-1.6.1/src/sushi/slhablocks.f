! This file is part of SusHi.
! 
! The routines of this file are used by the input routine
! to be found in inputoutput.f .
! 
C-{{{ subroutine slhablocks:

      subroutine slhablocks(blocktype,cline,ierr)
c..
c..   Read the data in string CLINE according to the format
c..   of block BLOCKTYPE.
c..   Note that BLOCKTYPE must ALWAYS be in capital letters, and
c..   of length *15 (trailing blanks).
c..
c..   Here the COMMON blocks are actually filled with data.
c..
c..   Example:
c..   BLOCKTYPE = 'MASS           '
c..   CLINE = ' 25  1.14272938e02  '
c..
c..   Then smassbo(25) would be given the value 1.14272938e02.
c..   In other words, M_h = 114 GeV.
c..   For the particle codes, see SUBROUTINE SLHAMAP.
c..
c..   Note that no checking of the correct format of CLINE is done.
c..   
      implicit real*8(a-h,o-z)
      real*8 rval,rval1
      integer ival
      character blocktype*15,cline*200
      include '../commons/common-slha.f'

      ierr=0

      if (blocktype.eq.'MASS          ') then
         ifound = 1
         read(cline,*,ERR=510) i,rval
         if (i.le.100) then
            smassbo(i) = rval
            lsmassbo(i) = .true.
         elseif (i.lt.2000000) then
            smassfe1(i-1000000) = rval
            lsmassfe1(i-1000000) = .true.
         elseif (i.lt.3000000) then
            smassfe2(i-2000000) = rval
            lsmassfe2(i-2000000) = .true.
         elseif (i.lt.4000000) then
            smassfe3(i-3000000) = rval
            lsmassfe3(i-3000000) = .true.
         elseif (i.lt.5000000) then
            smassfe4(i-4000000) = rval
            lsmassfe4(i-4000000) = .true.
         elseif (i.lt.6000000) then
            smassfe5(i-5000000) = rval
            lsmassfe5(i-5000000) = .true.
         else 
            smassfe6(i-6000000) = rval
            lsmassfe6(i-6000000) = .true.
         endif
      endif
      
      if (blocktype.eq.'TMP           ') then
         read(cline,*,ERR=510) i
         if (i.le.10) then
            read(cline,*,ERR=510) i,ival
            tmpint(i) = ival
            ltmpint(i) = .true.
         else
            read(cline,*,ERR=510) i,rval
            tmpreal(i) = rval
            ltmpreal(i) = .true.
         endif
      endif

      if (blocktype.eq.'SMINPUTS      ') then
         read(cline,*,ERR=510) i,rval
         sminputs(i) = rval
         lsminputs(i) = .true.
      endif

      if (blocktype.eq.'QQH           ') then
         read(cline,*,ERR=510) i
         if (i.le.10) then
            read(cline,*,ERR=510) i,ival
            qqhint(i) = ival
            lqqhint(i) = .true.
         else
            read(cline,*,ERR=510) i,rval
            qqhreal(i-10) = rval
            lqqhreal(i-10) = .true.
         endif
      endif

      if (blocktype.eq.'STOPMIX       ') then
         read(cline,*,ERR=510) i,j,rval
         stopmix(i,j) = rval
         lstopmix(i,j) = .true.
      endif

      if (blocktype.eq.'SBOTMIX       ') then
         read(cline,*,ERR=510) i,j,rval
         sbotmix(i,j) = rval         
         lsbotmix(i,j) = .true.
      endif

      if (blocktype.eq.'NMHMIX        ') then
         read(cline,*,ERR=510) i,j,rval
         nmhmix(i,j) = rval         
         lnmhmix(i,j) = .true.
      endif

      if (blocktype.eq.'NMAMIX        ') then
         read(cline,*,ERR=510) i,j,rval
         nmamix(i,j) = rval         
         lnmamix(i,j) = .true.
      endif

      if (blocktype.eq.'NMAMIXR       ') then
         read(cline,*,ERR=510) i,j,rval
         nmamixr(i,j) = rval         
         lnmamixr(i,j) = .true.
      endif

      if (blocktype.eq.'MINPAR        ') then
         read(cline,*,ERR=510) i,rval
         minpar(i) = rval
         lminpar(i) = .true.
      endif

      if (blocktype.eq.'ALPHA         ') then
         read(cline,*,ERR=510) rval
         slhaalpha = rval
         lslhaalpha = .true.
      endif

      if (blocktype.eq.'2HDM          ') then
         read(cline,*,ERR=510) rval
         twohdmpar = rval
         ltwohdmpar = .true.
      endif
      
      if (blocktype.eq.'GAUGE         ') then
         read(cline,*,ERR=510) i,rval
         gauge(i) = rval
         lgauge(i) = .true.
      endif

      if (blocktype.eq.'SUSHI         ') then
         read(cline,*,ERR=510) i
         if (i.eq.4) then
            read(cline,*,ERR=510) i,rval
            sushipar(i) = rval
            lsushipar(i) = .true.
         elseif (i.eq.8) then
            read(cline,*,ERR=510) i,pdfname(2)
            lpdfname(2)=.true.
         elseif (i.eq.9) then
            read(cline,*,ERR=510) i,pdfname(3) 
            lpdfname(3)=.true.
         else
            read(cline,*,ERR=510) i,j
            sushipar(i) = j
            lsushipar(i) = .true.
         endif
      endif

      if (blocktype.eq.'PDFSPEC       ') then
         read(cline,*,ERR=510) i
         if (i.eq.1) then
            read(cline,*,ERR=510) i,pdfname(1)
            lpdfname(1)=.true.
         elseif (i.eq.2) then
            read(cline,*,ERR=510) i,pdfname(2)
            lpdfname(2)=.true.
         elseif (i.eq.3) then
            read(cline,*,ERR=510) i,pdfname(3)
            lpdfname(3)=.true.
         elseif (i.eq.4) then
            read(cline,*,ERR=510) i,pdfname(4)
            lpdfname(4)=.true.
         else
            read(cline,*,ERR=510) i,j
            pdfspecpar(i) = j
            lpdfspecpar(i) = .true.
         endif
      endif

      if (blocktype.eq.'VEGAS         ') then
         read(cline,*,ERR=510) i,j
         vegaspar(i) = j
         lvegaspar(i) = .true.
      endif

      if (blocktype.eq.'FACTORS       ') then
         read(cline,*,ERR=510) i,rval
         yukawafac(i) = rval
         lyukawafac(i) = .true.
      endif
      
      if (blocktype.eq.'4GEN          ') then
         read(cline,*,ERR=510) i,rval
         fourgen(i) = rval
         lfourgen(i) = .true.
      endif

      if (blocktype.eq.'SCALES        ') then
         read(cline,*,ERR=510) i
         if (i.eq.101) then
            read(cline,*,ERR=510) i,rval1,rval2 !,ival
            frgghmin = rval1
            frgghmax = rval2
            lmurscanpar(1) = .true.
         else if (i.eq.102) then
            read(cline,*,ERR=510) i,rval1,rval2,ival
            frgghmintab = rval1
            frgghmaxtab = rval2
            frgghstep = ival
            lmurscanpar(2) = .true.
         else
            read(cline,*,ERR=510) i,rval
            scalespar(i) = rval
            lscalespar(i) = .true.
         endif
      endif

      if (blocktype.eq.'RENORMBOT     ') then
         read(cline,*,ERR=510) i,rval
         renormbotpar(i) = rval
         lrenormbotpar(i) = .true.
      endif

      if (blocktype.eq.'RENORMSBOT    ') then
         read(cline,*,ERR=510) i,rval
         renormsbotpar(i) = rval
         lrenormsbotpar(i) = .true.
      endif

      if (blocktype.eq.'DISTRIB       ') then
         read(cline,*,ERR=510) i
         if (i.le.5) then
         read(cline,*,ERR=510) i,j
         distribpar(i) = j
         else
         read(cline,*,ERR=510) i,rval
         distribpar(i) = rval
         endif
         ldistribpar(i) = .true.
      endif

      if (blocktype.eq.'FEYNHIGGS     ') then
         read(cline,*,ERR=510) i,rval
         feynhiggspar(i) = rval
         lfeynhiggspar(i) = .true.
      endif

      if (blocktype.eq.'FEYNHIGGSFLAGS') then
         read(cline,*,ERR=510) i,j
         fhflagspar(i) = j
         lfhflags(i) = .true.
      endif

      if (blocktype.eq.'2HDMC         ') then
         read(cline,*,ERR=510) i,rval
         thdmcpar(i) = rval
         lthdmcpar(i) = .true.
      endif

      if (blocktype.eq.'EXTPAR        ') then
         read(cline,*,ERR=510) i,rval
         extpar(i) = rval
         lextpar(i) = .true.
      endif

      if (blocktype.eq.'GGHSOFT       ') then
         read(cline,*,ERR=510) i,ival1,ival2,ival3
         gghsoftpar(i,1) = ival1
         gghsoftpar(i,2) = ival2
         gghsoftpar(i,3) = ival3
         lgghsoftpar(i) = .true.
      endif

      if (blocktype.eq.'GGHMT         ') then
         read(cline,*,ERR=510) i,rval
         if (i.eq.9) then
            gghmtreal(i) = rval
            lgghmtreal(i) = .true.
         else
            gghmtpar(i) = rval
            lgghmtpar(i) = .true.
         endif
      endif

      if (blocktype.eq.'DIM5          ') then
         read(cline,*,ERR=510) i
         if (i.eq.0) then
         read(cline,*,ERR=510) i,j
         else
         read(cline,*,ERR=510) i,rval
         end if
!         if (i.le.10) then
!            call printerrorsushi
!     &           (1,6,'Block DIM5: specify Higgs by 11,12,...')
!         endif
         if (i.eq.0) then
            dim5runpar = j
            ldim5runpar = .true.
         elseif (i.lt.100) then
            dim5par(i,0) = rval
            ldim5par(i,0) = .true.
         elseif (i.lt.200) then
            dim5par(i-100,1) = rval
            ldim5par(i-100,1) = .true.
         elseif (i.lt.300) then
            dim5par(i-200,2) = rval
            ldim5par(i-200,2) = .true.
         elseif (i.lt.400) then
            dim5par(i-300,3) = rval
            ldim5par(i-300,3) = .true.
         else
            write(6,*) 'Block DIM5, entry',i
            call printerrorsushi(1,6,'Block DIM5: wrong format')
         endif
      endif

      return

 510  ierr=1
      
      end

C-}}}
