! This file is part of SusHi.
! 
! It includes the input and output routines of SusHi.
! 

C-{{{ subroutine readinsushi

      subroutine readinsushi()

      implicit none
c..
c..   extflag: .true., if called from external program
c..
      logical ifilein,ppcoll,int2bool
      double precision mapreal8,mapreal8def,Amxpre(10,10)
      integer mapint,mapintdef,mapintdef2,i,ndum
      integer count1, count2, count3, blockflag, rensbotc
      integer ntmp
      character blocks(20)*15, bltype*15
      character tostring*10
      double precision mstop1,mstop2,msbot1,msbot2

      include '../commons/common-inputoutput.f'
      !internal definitions
      include '../commons/common-keys.f'
      include '../commons/common-slha.f'
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-quark.f'
      include '../commons/common-int.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-flags.f'
      include '../commons/common-expand.f'
      include '../commons/common-lumis.f'

      common /coll/ ppcoll

      inquire(file=jfilein,exist=ifilein)
      if(.not.ifilein)then
         call printerrorsushi(1,6,'ERROR: FILE '//jfilein/
     &        /' DOES NOT EXIST!')
      endif
      open(unit=12,file=jfilein)

      fhflag  = .false.
      thdmcflag = .false.
      onshellflag = .false.
      lqqhkey = .false.
      nnlostop=.false.

!     setting read-in-flags to false
      lqqhint = .false.
      lqqhreal = .false.
      ltmpint = .false.
      ltmpreal = .false.
      lpdfname = .false.
      lsmassbo   = .false.
      lsmassfe1  = .false.
      lsmassfe2  = .false.
      lsminputs  = .false.
      lminpar    = .false.
      lsushipar  = .false.
      ldistribpar    = .false.
      lextpar        = .false.
      lstopmix       = .false.
      lvegaspar      = .false.
      lyukawafac     = .false.
      lfourgen       = .false.
      lfeynhiggspar  = .false.
      lscalespar     = .false.
      lmurscanpar    = .false.
      lrenormbotpar  = .false.
      lrenormsbotpar = .false.
      lfhflags       = .false.
      lslhaalpha     = .false.
      lexpandc1o1    = .true.

      !read in of the main parameters for SUSHI
      blocks = ' '
      blocks(1) = 'SUSHI'
      call readblocks(12,blocks)

      verbose = mapreal8def(sushipar(19),lsushipar(19),1.d0
     & ,'SUSHI(19) defaults to 1: Normal output.')

!     SUSHI parameters
      model = mapint(int(sushipar(1)),lsushipar(1),'SUSHI(1) '//
     &     'undefined.')
      htype= mapint(int(sushipar(2)),lsushipar(2),'SUSHI(2) undefined.')

c..   Currently, we do not allow for factorizing out C1, because
c..   it is not compatible with 'ewreweight' [rh, 22Apr2016]:
      ndum = mapintdef2(0,int(sushipar(11)),lsushipar(11),1,1
     &     ,1,'SUSHI(11) defaults to 1: C1 will be perturbatively'/
     &     /'expanded.','SUSHI(11) out of range.')
      lexpandc1o1 = int2bool(ndum)
      
c..   translate from old notation for Higgs type to NMSSM-style:
      select case (htype)
      case (0)
         htype=11
      case (1)
         htype=21
      case (2)
         htype=12
      end select

      if (((model.eq.0).and.((htype.ne.11).and.(htype.ne.21))).or.
     &     (((model.eq.1).or.(model.eq.2)).and.((htype.eq.13).or.
     &     (htype.eq.22)))) then
         write(6,*) 'Block SUSHI, entry 1 =',int(sushipar(1))
         write(6,*) 'Block SUSHI, entry 2 =',int(sushipar(2))
         call printerrorsushi(1,6,'Incompatible input.')
      endif

      select case (htype)
      case (11)                 ! SM Higgs, light Higgs, or H1
         pseudo=0
         sind=1
      case (12)                 ! heavy Higgs or H2
         pseudo=0
         sind=2
      case (13)                 ! H3
         pseudo=0
         sind=3
      case (21)                 ! pseudo-scalar Higgs or A1
         pseudo=1
         sind=2             
      case (22)                 ! A2
         pseudo=1
         sind=3
      case default
         write(6,*) 'Block SUSHI, entry 1 =',int(sushipar(1))
         write(6,*) 'Block SUSHI, entry 2 =',int(sushipar(2))
         call printerrorsushi(1,6,'Incompatible input.')
      end select

      if(.not.extflag) then
         if (.not.lsushipar(7)) call printdie('SUSHI(7) undefined. ')
         call printcopysushi(1,6,verbose)
         nsubprocggh = mapreal8def(sushipar(20),lsushipar(20),0.d0
     &        ,'SUSHI(20) defaults to 0: '/
     &        /'All ggh@nnlo subprocesses included.')
         nsubprocbbh = mapreal8def(sushipar(21),lsushipar(21),0.d0
     &        ,'SUSHI(21) defaults to 0: '/
     &        /'All bbh@nnlo subprocesses included.')
      else
         call printcopysushi(0,6,verbose)
      endif

      if(sushipar(7).eq.1) then
         ew = 1
      else if (sushipar(7).eq.2) then
         ew = 2
      else
         ew = 0
      endif

      ldim5 = .false.
      c5sp = 0.d0
      bltype = 'DIM5'
      call readslha(12,bltype,blockflag)
      if (blockflag.eq.1) then
         ldim5 = .true.
            dim5run = mapintdef2(0,dim5runpar,ldim5runpar,0,2
     &        ,1,'','DIM5(0) out of range.')
         if (pseudo.eq.0) then
          c5spin(0,0)=mapreal8def(dim5par(htype,0),ldim5par(htype,0)
     &           ,0.d0,'DIM5(0,0) defaults to 0.d0: '/
     &           /'no extra dim-5 operator included')
          c5spin(0,1)=mapreal8def(dim5par(htype,1),ldim5par(htype,1)
     &           ,0.d0,'DIM5(0,1) defaults to 0.d0: '/
     &           /'no NLO dim-5 operator included')
          c5spin(0,2)=mapreal8def(dim5par(htype,2),ldim5par(htype,2)
     &           ,0.d0,'DIM5(0,2) defaults to 0.d0: '/
     &           /'no NNLO dim-5 operator included')
          c5spin(0,3)=mapreal8def(dim5par(htype,3),ldim5par(htype,3)
     &           ,0.d0,'DIM5(0,3) defaults to 0.d0: '/
     &           /'no N3LO dim-5 operator included')
         else
          c5spin(1,0)=mapreal8def(dim5par(htype,0),ldim5par(htype,0)
     &           ,0.d0,'DIM5(1,0) defaults to 0.d0: '/
     &           /'no extra dim-5 operator included')
          c5spin(1,1)=mapreal8def(dim5par(htype,1),ldim5par(htype,1)
     &           ,0.d0,'DIM5(1,1) defaults to 0.d0: '/
     &           /'no NLO dim-5 operator included')
          c5spin(1,2)=mapreal8def(dim5par(htype,2),ldim5par(htype,2)
     &           ,0.d0,'DIM5(1,2) defaults to 0.d0: '/
     &           /'no NNLO dim-5 operator included')
         endif
         c5sp = c5spin
      endif


      !---NMSSM
      if (model.eq.3) then
      !generic NMSSM - mixing matrices provided by hand
         blocks = ' '
         !on-shell input
         bltype = 'STOPMIX'
         call readslha(12,bltype,blockflag)
         if(blockflag.eq.1) then
          call printinfosushi(6,
     &    "SusHi was called with on-shell squark sector input:")
          call printinfosushi(6,
     &    "Please guarantee the consistency of the renormalization")
          call printinfosushi(6,
     &    "scheme using on-shell parameter input.")
          blocks(11) = 'SBOTMIX'
          blocks(12) = 'STOPMIX'
          onshellflag = .true.
         end if

         blocks(1) = 'SUSHI'
         blocks(2) = 'SMINPUTS'
         blocks(3) = 'MINPAR'
         blocks(4) = 'EXTPAR'
         blocks(5) = 'RENORMBOT'
         blocks(6) = 'RENORMSBOT'
         blocks(7) = 'FACTORS'
         blocks(8) = 'MASS'
         blocks(9) = 'NMHMIX'

         bltype = 'NMAMIX'
         call readslha(12,bltype,blockflag)
         if(blockflag.eq.1) then
            call printinfosushi(6,
     &    "CP-odd Higgs mixing matrix read in without prerotation.")
         blocks(10) = 'NMAMIX'
         else
            call printinfosushi(6,
     &    "CP-odd Higgs mixing matrix read in with prerotation.")
         blocks(10) = 'NMAMIXR'
         end if
         call readblocks(12,blocks)

      !---MSSM
      else if (model.eq.1) then
      !looking for the block FEYNHIGGS
         bltype = 'FEYNHIGGS'
         call readslha(12,bltype,blockflag)

         if(blockflag.eq.1) then
            call printinfosushi(6,
     &           "SusHi was called with the SLHA-block 'FEYNHIGGS':")
            call printinfosushi(6,
     &           "Thus FeynHiggs is used for the calculation of SUSY")
            call printinfosushi(6,
     &           "Higgs masses, mu and stop, sbottom masses/angles.")
            fhflag  = .true.
         end if

      !read-in not dependent on FH-use
         blocks = ' '
         blocks(1) = 'SUSHI'
         blocks(2) = 'SMINPUTS'
         blocks(3) = 'MINPAR'
         blocks(4) = 'EXTPAR'
         blocks(5) = 'RENORMBOT'
         blocks(6) = 'RENORMSBOT'
         blocks(7) = 'FACTORS'
         call readblocks(12,blocks)

         !read-in dependent on FH-use
         blocks = ' '
         if (.not.fhflag) then
            bltype = 'STOPMIX'
            call readslha(12,bltype,blockflag)
            if(blockflag.eq.1) then
            call printinfosushi(6,
     &      "SusHi was called with on-shell squark sector input:")
            call printinfosushi(6,
     &      "Please guarantee the consistency of the renormalization")
            call printinfosushi(6,
     &      "scheme using on-shell parameter input.")
            blocks(3) = 'SBOTMIX'
            blocks(4) = 'STOPMIX'
            onshellflag = .true.
            end if
            blocks(1) = 'ALPHA'
            blocks(2) = 'MASS'
            call readblocks(12,blocks)
         else
            blocks(1) = 'FEYNHIGGS'
            blocks(2) = 'FEYNHIGGSFLAGS'
            call readblocks(12,blocks)
            !standard flags for FH 2.9 and higher
            mssmpart = 4
            fieldren = 0
            tanbren = 0
            higgsmix = 2
            p2approx = 4
            looplevel = 2
            loglevel = 0
            runningMT = 1
            botResum = 1
            tlCplxApprox = 0
            if (lfhflags(1)) mssmpart = fhflagspar(1)
            if (lfhflags(2)) fieldren = fhflagspar(2)
            if (lfhflags(3)) tanbren = fhflagspar(3)
            if (lfhflags(4)) higgsmix = fhflagspar(4)
            if (lfhflags(5)) p2approx = fhflagspar(5)
            if (lfhflags(6)) looplevel = fhflagspar(6)
            if (lfhflags(7)) runningMT = fhflagspar(7)
            if (lfhflags(8)) botResum = fhflagspar(8)
            if (lfhflags(9)) tlCplxApprox = fhflagspar(9)
            if (lfhflags(10)) loglevel = fhflagspar(10)
         end if

      !---2HDM
      else if (model.eq.2) then
!     looking for the block 2HDMC
         bltype = '2HDMC'
         call readslha(12,bltype,blockflag)
         if(blockflag.eq.1) then
            call printinfosushi(6,
     &           "The input file contains the block '2HDMC':")
            call printinfosushi(6,
     &           "Thus 2HDMC is used for the calculation of the")
            call printinfosushi(6,
     &           "Higgs masses and the quark Yukawa couplings.")
            thdmcflag  = .true.
         else
            thdmcflag  = .false.
         end if

      !read-in not dependent on 2HDMC-use:
         blocks = ' '
         blocks(1) = 'SUSHI'
         blocks(2) = 'SMINPUTS'
         blocks(3) = 'RENORMBOT'
         blocks(4) = 'FACTORS'
         call readblocks(12,blocks)

         !read-in dependent on 2HDMC-use
         blocks = ' '
         if (.not.thdmcflag) then
            blocks(1) = 'ALPHA'
            blocks(2) = 'MASS'
            blocks(3) = '2HDM'
            blocks(4) = 'MINPAR'
            call readblocks(12,blocks)
         else
            blocks(1) = '2HDMC'
            call readblocks(12,blocks)
         end if

      !---SM
      else if (model.eq.0) then

         bltype = 'FEYNHIGGS'
         call readslha(12,bltype,blockflag)

         if(blockflag.eq.1) then
            call printinfosushi(6,
     &           "SusHi was called with the SLHA-block"//
     &           " FEYNHIGGS in the SM:")
            call printinfosushi(6,"Thus FeynHiggs is used for"//
     &           " the calculation of a Higgs mass,")
            call printinfosushi(6,
     &           "however the Higgs is the SM Higgs!")
            fhflag  = .true.
         end if

         blocks = ' '
         if (.not.fhflag) then
            blocks(1) = 'SUSHI'
            blocks(2) = 'SMINPUTS'
            blocks(3) = 'MASS'
            blocks(4) = 'RENORMBOT'
            blocks(5) = 'FACTORS'
            call readblocks(12,blocks)
         else
            blocks(1) = 'SUSHI'
            blocks(2) = 'SMINPUTS'
            blocks(3) = 'MINPAR'
            blocks(4) = 'EXTPAR'
            blocks(5) = 'RENORMBOT'
            blocks(6) = 'RENORMSBOT'
            blocks(7) = 'FACTORS'
            blocks(8) = 'FEYNHIGGS'
            blocks(9) = 'FEYNHIGGSFLAGS'
            call readblocks(12,blocks)
            !standard flags for FH 2.9 and higher
            mssmpart = 4
            fieldren = 0
            tanbren = 0
            higgsmix = 2
            p2approx = 4
            looplevel = 2
            loglevel = 0
            runningMT = 1
            botResum = 1
            tlCplxApprox = 0
            if (lfhflags(1)) mssmpart = fhflagspar(1)
            if (lfhflags(2)) fieldren = fhflagspar(2)
            if (lfhflags(3)) tanbren = fhflagspar(3)
            if (lfhflags(4)) higgsmix = fhflagspar(4)
            if (lfhflags(5)) p2approx = fhflagspar(5)
            if (lfhflags(6)) looplevel = fhflagspar(6)
            if (lfhflags(7)) runningMT = fhflagspar(7)
            if (lfhflags(8)) botResum = fhflagspar(8)
            if (lfhflags(9)) tlCplxApprox = fhflagspar(9)
            if (lfhflags(10)) loglevel = fhflagspar(10)
!to be used for a fixed top mass within FH
!          runningMT = 0
         end if

      else
         write(*,105)
         write(*,*) 
     &"Error: We are sorry, but the specified model"
         write(*,*) 
     &"is no reasonable SusHi dish."
         write(*,105)
         stop
      end if

      !---Setting parameters

      !Standard model input parameters
      SUinvAlfa=mapreal8(sminputs(1),lsminputs(1)
     &     ,'SMINPUTS(1): alpha_em undefined.')
      GF = mapreal8(sminputs(2),lsminputs(2)
     &     ,'SMINPUTS(2): G_Fermi undefined. ')
      AlfasMZ=mapreal8(sminputs(3),lsminputs(3),
     &     'SMINPUTS(3): alphas(MZ) undefined.')
      MZ = mapreal8(sminputs(4),lsminputs(4)
     &     ,'SMINPUTS(4): Mz undefined. ')
      mbmb=mapreal8(sminputs(5),lsminputs(5)
     &     ,'SMINPUTS(5): bottom mass undefined. ')
      MT = mapreal8(sminputs(6),lsminputs(6)
     &     ,'SMINPUTS(6): top mass undefined. ')
      mt2 = mt**2
      MS = -1
      MC = -1
      ME = -1
      MU = -1
      MD = -1
      MM = -1
      MS = -1
      ML = -1

      Hmx = 0.d0
      Amx = 0.d0
      Hmx(1,1) = 1.d0
      Amx(2,2) = 1.d0

      !calculating mW from mZ, GF and alpha_EM
      mW = dsqrt(0.25D0- (1.D0/SUinvAlfa) * Pi / (dsqrt(2.D0)*GF*mZ*mZ))
      mW = mZ * mZ * (0.5D0 + mW) / 0.985D0
      mW = dsqrt(mW)

      lam = 1.d0
      kap = 1.d0
      Alam = 1.d0
      Akap = 1.d0

c..   shall we calculate qq'->H instead of bb->H ?
      ptn1=5
      ptn2=-5
      qqhyuk = mbmb
      qqhscale = mbmb
      qqhint = 0
      qqhreal = 0.d0
      bltype = 'QQH'
      qqhstring='bbh'
      call readslha(12,bltype,blockflag)
      if (blockflag.eq.1) then
         lqqhkey = .true.
         qqhstring='qqh'
         ptn1 = mapintdef2(0,qqhint(1),lqqhint(1),-5,5,5,''
     &        ,'QQH(1) out of range')
         ptn2 = mapintdef2(0,qqhint(2),lqqhint(2),-5,5,-ptn1,''
     &        ,'QQH(2) out of range')
         qqhyuk = mapreal8def(qqhreal(1),lqqhreal(1),mbmb
     &        ,'QQH(11) defaults to mb(mb)')
         qqhscale = mapreal8def(qqhreal(2),lqqhreal(2),mbmb
     &        ,'QQH(12) defaults to mb(mb)')
      endif

      bltype = 'TMP'
      ltmpblock = .false.
      call readslha(12,bltype,blockflag)
      if (blockflag.eq.1) then
         ltmpblock = .true.
         hxswg = mapintdef2(0,tmpint(1),ltmpint(1),0,1,0,''
     &        ,'TMP(1) out of range')
      endif

      !---NMSSM
      if (model.eq.3) then
         !Minpar parameters
         tanb = mapreal8(minpar(3),lminpar(3)
     &        ,'MINPAR(3): tanb undefined. ')
         !SUSY parameters from extpar or external masses
         Mgl = mapreal8(dabs(extpar(3)),lextpar(3)
     &     ,'EXTPAR(3): gluino mass undefined. ')
         !muSUSY = mapreal8(extpar(23),lextpar(23),'mu undefined. ')
            if (onshellflag) then
             !either read in on-shell masses
             Msbot1 = mapreal8(smassfe1(5),lsmassfe1(5)
     &       ,'Sbottom mass 1 undefined. ')
             Msbot2 = mapreal8(smassfe2(5),lsmassfe2(5)
     &       ,'Sbottom mass 2 undefined. ')
             Mstop1 = mapreal8(smassfe1(6),lsmassfe1(6)
     &       ,'Stop mass 1 undefined. ')
             Mstop2 = mapreal8(smassfe2(6),lsmassfe2(6)
     &       ,'Stop mass 2 undefined. ')
             cthetat = mapreal8(stopmix(1,1),lstopmix(1,1)
     &       ,'Stop mixing undefined. ')
             sthetat = mapreal8(stopmix(1,2),lstopmix(1,2)
     &       ,'Stop mixing undefined. ')
             cthetab = mapreal8(sbotmix(1,1),lsbotmix(1,1)
     &       ,'Sbottom mixing undefined. ')
             sthetab = mapreal8(sbotmix(1,2),lsbotmix(1,2)
     &       ,'Sbottom mixing undefined. ')
            msb12 = msbot1**2
            msb22 = msbot2**2
            mst12 = mstop1**2
            mst22 = mstop2**2
            else
             !or calculate them internally later in renscheme.f
             !original option of SusHi
             M3SQ = mapreal8(extpar(43),lextpar(43),'M_SQ undefined. ')
             M3SU = mapreal8(extpar(46),lextpar(46),'M_SU undefined. ')
             M3SD = mapreal8(extpar(49),lextpar(49),'M_SD undefined. ')
            end if
         cAt = mapreal8(extpar(11),lextpar(11),'At undefined. ')
         cAb = mapreal8(extpar(12),lextpar(12),'Ab undefined. ')
         lam = mapreal8(extpar(61),lextpar(61),'lambda undefined. ')
         !kap = mapreal8(extpar(62),lextpar(62),'kappa undefined. ')
         !Alam = mapreal8(extpar(63),lextpar(63),'Alambda undefined. ')
         !Akap = mapreal8(extpar(64),lextpar(64),'Akappa undefined. ')
         if (lextpar(23)) then
           write(*,105)
           write(*,*) "Entry 23 of EXTPAR is ignored in the NMSSM."
           write(*,*) "Please use entry 65 for an eff. mu instead."
           write(*,105)
           stop
         end if
         muSUSY = mapreal8(extpar(65),lextpar(65),'effmu undefined. ')
         !read in mass spectrum
         if (pseudo.eq.0) then
            Mh = mapreal8(smassbo(25+10*(Sind-1))
     &       ,lsmassbo(25+10*(Sind-1)),'Higgs mass undefined. ')
         else if (pseudo.eq.1) then
            Mh = mapreal8(smassbo(26+10*(Sind-1))
     &       ,lsmassbo(26+10*(Sind-1)),'Higgs mass undefined. ')
         end if
         !read in mixing matrices
         Hmx = 0.d0
         Amx = 0.d0
         if (lnmamix(2,2)) then
           Amxpre = 0.d0
           Amxpre(1,1) = dcos(atan(tanb))
           Amxpre(1,2) = -dsin(atan(tanb))
           Amxpre(2,1) = dsin(atan(tanb))
           Amxpre(2,2) = dcos(atan(tanb))
           Amxpre(3,3) = 1.d0
           do count1 = 2,3
           do count2 = 1,3
            !misuse of Hmx, to be correctly set later
            !changed to count1-1 on 18/7/2015
            Hmx(count1,count2) = mapreal8(nmamix(count1-1,count2),
     &       lnmamix(count1-1,count2),'NMAMIX partially undefined. ')
            !write(*,*) 'Hmxcheck1',count1,count2,Hmx(count1,count2)
           end do
           end do
           do count1 = 2,3
           do count2 = 1,3
            Amx(count1,count2) = 0.d0
            do count3 = 1,3
            Amx(count1,count2) = Amx(count1,count2)
     &        + Hmx(count1,count3)*Amxpre(count2,count3)
            end do
           end do
           end do
!           !testing by rotating back:
!           Hmx = 0.d0
!           do count1=2,3
!           do count2=1,3
!           do count3=1,3
!             Hmx(count1,count2) = Hmx(count1,count2)
!     &             + Amx(count1,count3)*Amxpre(count3,count2)
!           end do
!            write(*,*) 'Hmxcheck2',count1,count2,Hmx(count1,count2)
!           end do
!           end do
         else
           do count1 = 2,3
           do count2 = 2,3
           !changed to count1-1,count2-1 on 18/7/2015
           Amx(count1,count2) = mapreal8(nmamixr(count1-1,count2-1),
     &      lnmamixr(count1-1,count2-1),'NMAMIXR partially undefined. ')
           end do
           end do
         end if

         do count1 = 1,3
         do count2 = 1,3
         Hmx(count1,count2) = mapreal8(nmhmix(count1,count2),
     &     lnmhmix(count1,count2),'NMHMIX partially undefined. ')
         end do
         end do

      !---MSSM
      else if ((model.eq.1).or.((model.eq.0).and.fhflag)) then
         !Minpar parameters
         tanb = mapreal8(minpar(3),lminpar(3),
     &        'MINPAR(3): tanb undefined. ')

         !SUSY parameters from extpar or external masses
         if (.not.fhflag) then
            Mgl = mapreal8(dabs(extpar(3)),lextpar(3)
     &     ,'gluino mass undefined. ')
            alpha = mapreal8(slhaalpha,lslhaalpha,'alpha undefined. ')
            muSUSY = mapreal8(extpar(23),lextpar(23),'mu undefined. ')
            if (onshellflag) then
             !either read in on-shell masses
             Msbot1 = mapreal8(smassfe1(5),lsmassfe1(5)
     &       ,'Sbottom mass 1 undefined. ')
             Msbot2 = mapreal8(smassfe2(5),lsmassfe2(5)
     &       ,'Sbottom mass 2 undefined. ')
             Mstop1 = mapreal8(smassfe1(6),lsmassfe1(6)
     &       ,'Stop mass 1 undefined. ')
             Mstop2 = mapreal8(smassfe2(6),lsmassfe2(6)
     &       ,'Stop mass 2 undefined. ')
             cthetat = mapreal8(stopmix(1,1),lstopmix(1,1)
     &       ,'Stop mixing undefined. ')
             sthetat = mapreal8(stopmix(1,2),lstopmix(1,2)
     &       ,'Stop mixing undefined. ')
             cthetab = mapreal8(sbotmix(1,1),lsbotmix(1,1)
     &       ,'Sbottom mixing undefined. ')
             sthetab = mapreal8(sbotmix(1,2),lsbotmix(1,2)
     &       ,'Sbottom mixing undefined. ')
            msb12 = msbot1**2
            msb22 = msbot2**2
            mst12 = mstop1**2
            mst22 = mstop2**2
            else
             !or calculate them internally later in renscheme.f
             !original option of SusHi
             M3SQ = mapreal8(extpar(43),lextpar(43),'M_SQ undefined. ')
             M3SU = mapreal8(extpar(46),lextpar(46),'M_SU undefined. ')
             M3SD = mapreal8(extpar(49),lextpar(49),'M_SD undefined. ')
            end if
            cAt = mapreal8(extpar(11),lextpar(11),'At undefined. ')
            cAb = mapreal8(extpar(12),lextpar(12),'Ab undefined. ')
            if (pseudo.eq.0) then
             if (Sind.eq.1) then
               Mh = mapreal8(smassbo(25),lsmassbo(25)
     &              ,'light Higgs mass undefined. ')
             else if (Sind.eq.2) then
               Mh = mapreal8(smassbo(35),lsmassbo(35)
     &              ,'heavy Higgs mass undefined. ')
             end if
            else if (pseudo.eq.1) then
               Mh = mapreal8(smassbo(36),lsmassbo(36)
     &              ,'pseudoscalar mass undefined. ')
            end if
            Hmx = 0.d0
            Hmx(1,1) = -dsin(alpha)
            Hmx(1,2) = dcos(alpha)
            Hmx(2,1) = dcos(alpha)
            Hmx(2,2) = dsin(alpha)
            Amx = 0.d0
            Amx(2,2) = 1.d0
         else
            MUE = 
     &dcmplx(mapreal8(extpar(23),lextpar(23),'mu undefined. '))
            M3SQ = mapreal8(extpar(43),lextpar(43),'M_SQ undefined. ')
            M3SU = mapreal8(extpar(46),lextpar(46),'M_SU undefined. ')
            M3SD = mapreal8(extpar(49),lextpar(49),'M_SD undefined. ')
            cAt = mapreal8(extpar(11),lextpar(11),'At undefined. ')
            cAb = mapreal8(extpar(12),lextpar(12),'Ab undefined. ')
            M_3=dcmplx(mapreal8(extpar(3),lextpar(3),'M_3 undefined. '))
            MA0 = mapreal8(extpar(26),lextpar(26),'MA0 undefined. ')

            M_1=dcmplx(mapreal8(feynhiggspar(1),lfeynhiggspar(1),
     &'M_1 undefined. '))
            M_2=dcmplx(mapreal8(feynhiggspar(2),lfeynhiggspar(2),
     &'M_2 undefined. '))

            !Setting of trilinear couplings in FeynHiggs
            if (lfeynhiggspar(3)) then
            cAu = mapreal8(feynhiggspar(3),lfeynhiggspar(3),
     &'A for FH undefined. ')
            cAmu = cAu
            cAtau = cAu
            cAc = cAu
            cAs = cAu
            cAe = cAu
            cAd = cAu
            else
            cAtau = mapreal8(feynhiggspar(13),lfeynhiggspar(13),
     &'A_TAU for FH undefined. ')
            cAc = mapreal8(feynhiggspar(14),lfeynhiggspar(14),
     &'A_C for FH undefined. ')
            cAs = mapreal8(feynhiggspar(15),lfeynhiggspar(15),
     &'A_S for FH undefined. ')
            cAmu = mapreal8(feynhiggspar(16),lfeynhiggspar(16),
     &'A_MU for FH undefined. ')
            cAu = mapreal8(feynhiggspar(17),lfeynhiggspar(17),
     &'A_U for FH undefined. ')
            cAd = mapreal8(feynhiggspar(18),lfeynhiggspar(18),
     &'A_D for FH undefined. ')
            cAe = mapreal8(feynhiggspar(19),lfeynhiggspar(19),
     &'A_E for FH undefined. ')
            end if

            !Setting of soft-breaking masses in FeynHiggs
            if (lfeynhiggspar(4)) then
            MSUSY = mapreal8(feynhiggspar(4),lfeynhiggspar(4),
     &'MSUSY for FH undefined. ')
            M3SL = MSUSY
            M3SE = MSUSY
            M2SL = MSUSY
            M2SE = MSUSY
            M2SQ = MSUSY
            M2SU = MSUSY
            M2SD = MSUSY
            M1SL = MSUSY
            M1SE = MSUSY
            M1SQ = MSUSY
            M1SU = MSUSY
            M1SD = MSUSY
            else
            M3SL = mapreal8(feynhiggspar(33),lfeynhiggspar(33),
     &'M3SL for FH undefined. ')
            M3SE = mapreal8(feynhiggspar(36),lfeynhiggspar(36),
     &'M3SE for FH undefined. ')
            M2SL = mapreal8(feynhiggspar(32),lfeynhiggspar(32),
     &'M2SL for FH undefined. ')
            M2SE = mapreal8(feynhiggspar(35),lfeynhiggspar(35),
     &'M2SE for FH undefined. ')
            M2SQ = mapreal8(feynhiggspar(42),lfeynhiggspar(42),
     &'M2SQ for FH undefined. ')
            M2SU = mapreal8(feynhiggspar(45),lfeynhiggspar(45),
     &'M2SU for FH undefined. ')
            M2SD = mapreal8(feynhiggspar(48),lfeynhiggspar(48),
     &'M2SD for FH undefined. ')
            M1SL = mapreal8(feynhiggspar(31),lfeynhiggspar(31),
     &'M1SL for FH undefined. ')
            M1SE = mapreal8(feynhiggspar(34),lfeynhiggspar(34),
     &'M1SE for FH undefined. ')
            M1SQ = mapreal8(feynhiggspar(41),lfeynhiggspar(41),
     &'M1SQ for FH undefined. ')
            M1SU = mapreal8(feynhiggspar(44),lfeynhiggspar(44),
     &'M1SU for FH undefined. ')
            M1SD = mapreal8(feynhiggspar(47),lfeynhiggspar(47),
     &'M1SD for FH undefined. ')
            end if

         end if

      !---SM
      else if ((model.eq.0).and.(.not.fhflag)) then
         
         Mh = mapreal8(smassbo(25),lsmassbo(25)
     &        ,'scalar Higgs mass undefined. ')

      end if

!---  2HDM
      if (model.eq.2) then
         if (.not.thdmcflag) then
            twohdmver = mapreal8(twohdmpar,ltwohdmpar,
     &           '2HDM version undefined. ')
            tanb = mapreal8(minpar(3),lminpar(3),'tanb undefined. ')
            alpha = mapreal8(slhaalpha,lslhaalpha,'alpha undefined. ')
            if (pseudo.eq.0) then
             if (Sind.eq.1) then
               Mh = mapreal8(smassbo(25),lsmassbo(25)
     &              ,'light Higgs mass undefined. ')
             else if (Sind.eq.2) then
               Mh = mapreal8(smassbo(35),lsmassbo(35)
     &              ,'heavy Higgs mass undefined. ')
             end if
            else if (pseudo.eq.1) then
               Mh = mapreal8(smassbo(36),lsmassbo(36)
     &              ,'pseudoscalar mass undefined. ')
            end if

            Hmx = 0.d0
            Hmx(1,1) = -dsin(alpha)
            Hmx(1,2) = dcos(alpha)
            Hmx(2,1) = dcos(alpha)
            Hmx(2,2) = dsin(alpha)
            Amx = 0.d0
            Amx(2,2) = 1.d0

         endif
      end if
         
!     some settings for FeynHiggs if used
      if ((model.eq.1).or.((model.eq.0).and.fhflag)) then
         CKMlambda = -1
         CKMA      = -1
         CKMrhobar = -1
         CKMetabar = -1
         
         scalefactor = 1.d0

         Qtau = 0.d0
         Qt = 0.d0
         Qb = 0.d0

         MHp = 0.d0
      end if

      !Yukawa factors
      yukfac(1)=mapreal8def(yukawafac(1),lyukawafac(1),1.d0,
     &     'charm yukawa-coupling undefined. ')
      yukfac(2)=mapreal8def(yukawafac(2),lyukawafac(2),1.d0,
     &     'top yukawa-coupling undefined. ')
      yukfac(3)=mapreal8def(yukawafac(3),lyukawafac(3),1.d0,
     &     'bottom yukawa-coupling undefined. ')

      if ((model.eq.1).or.(model.eq.3)) then
         yukfac(6)=mapreal8def(yukawafac(4),lyukawafac(4),1.d0,
     &'stop yukawa-coupling undefined. ')
         yukfac(7)=mapreal8def(yukawafac(5),lyukawafac(5),1.d0,
     &'sbottom yukawa-coupling undefined. ')
      else
         yukfac(6) = 0.d0
         yukfac(7) = 0.d0
      end if

      Mcharm = 1.28d0
      Mbp = 1000.d0
      Mtp = 1000.d0

      if(yukfac(1).ne.0.d0)
     & Mcharm = mapreal8(sminputs(8),lsminputs(8),
     & 'charm mass undefined. ')

      !4.Generation parameters
!Form of the corresponding block in the input files:
!Block 4GEN
!  1   3.00000000e+02   # Mb'
!  2   3.00000000e+02   # Mt'
!  3   0.d0             # factor for yukawa-couplings: t'
!  4   0.d0             # b'
!  5   0.d0             # st'
!  6   0.d0             # sb'
      bltype = '4GEN'
      call readslha(12,bltype,blockflag)
      if (blockflag.eq.1) then
         blocks = ' '
         blocks(1) = '4GEN'
         call readblocks(12,blocks)
         yukfac(4) = mapreal8(fourgen(3),lfourgen(3),
     &        'top yukawa-coupling undefined. ')
         yukfac(5) = mapreal8(fourgen(4),lfourgen(4),
     &        'bottom yukawa-coupling undefined. ')
         if ((model.eq.1).or.(model.eq.3)) then
            yukfac(8) = mapreal8(fourgen(5),lfourgen(5),
     &           'stop yukawa-coupling undefined. ')
            yukfac(9) = mapreal8(fourgen(6),lfourgen(6),
     &           'sbottom yukawa-coupling undefined. ')
         else
            yukfac(8) = 0.d0
            yukfac(9) = 0.d0
         end if
         if(yukfac(5).ne.0.d0)
     & Mbp = mapreal8(fourgen(1),lfourgen(1),'b mass undefined. ')
         if(yukfac(4).ne.0.d0)
     & Mtp = mapreal8(fourgen(2),lfourgen(2),'t mass undefined. ')
      else
         yukfac(4) = 0.d0
         yukfac(5) = 0.d0
         yukfac(8) = 0.d0
         yukfac(9) = 0.d0
      end if

      !bottom sector
      if (lrenormbotpar(1)) then
         mbrunyuk = renormbotpar(1)
      else
         mbrunyuk = 0
      endif

      if (lrenormbotpar(2)) then
         tanbresum = renormbotpar(2)
      if (model.eq.0) then
         write(*,105)
         write(*,*)
     &"Tan(b)-Resummation for the SM is not reasonable!"
         write(*,105)
      else
      if (mbrunyuk.gt.0) then
         write(*,105)
         write(*,*)
     &"Running bottom mass in the MSSM is always resummed!"
         write(*,105)
      endif
      endif
      else
         tanbresum = 0
      endif

      if (lrenormbotpar(3)) then
         delmbfh = renormbotpar(3)
      else
         delmbfh = 0
      endif

      if (lrenormbotpar(4)) then
         mbossave = renormbotpar(4)
      else
         mbossave = 0.D0
      endif

      if (lrenormbotpar(5)) then
         mbrunloop = renormbotpar(5)
      else
         mbrunloop = 0
      endif

      !sbottom sector
      if (((model.eq.1).or.(model.eq.3)).and.(.not.onshellflag)) then
      !first check whether all flags are set
         if (.not.lrenormsbotpar(1)) then
            write(*,105)
            write(*,*)
     &"Please choose a scheme for mb in the sbottom sector!"
            write(*,105)
            stop
         end if
         if (.not.lrenormsbotpar(2)) then
            write(*,105)
            write(*,*)
     &"Please choose a scheme for Ab in the sbottom sector!"
            write(*,105)
            stop
         end if
         if (.not.lrenormsbotpar(3)) then
            write(*,105)
            write(*,*)
     &"Please choose a scheme for thetab in the sbottom sector!"
            write(*,105)
            stop
         end if

         rensbotc = 0
         do count1=1,3
            if (renormsbotpar(count1).eq.2) rensbotc = rensbotc + 1
         enddo
         if (rensbotc.ne.1) then
            write(*,105)
            write(*,*)
     &"Exactly one sbottom parameter has to be choosen -dependent-!"
            write(*,105)
            stop
         endif
         mbsbch = renormsbotpar(1)
         Abch = renormsbotpar(2)
         thetabch = renormsbotpar(3)
      endif

      !standard settings for 2HDM and onshell input:
      if (model.eq.2) then
         tanbresum = 0
         mbsbch = 2
         Abch = 0
         thetabch = 0
      end if
      if (onshellflag) then
         mbsbch = 2
         Abch = 0
         thetabch = 0
      end if

      if(extflag) then
         norderggh = 1
         norderbbh = -1
         close(12)
      else
         norderggh = mapint(int(sushipar(5)),lsushipar(5)
     &        ,'SUSHI(5) undefined.')
         if (norderggh.ge.12) then
            norderggh=norderggh-10
            nnlostop=.true.
         endif
         norderbbh = mapint(int(sushipar(6)),lsushipar(6)
     &        ,'SUSHI(6) undefined.')
         if (((norderbbh.ge.3).and.(norderbbh.le.9))
     &        .or.(norderbbh.ge.13).or.(norderbbh.le.-2)) then
            call printerrorsushi(1,6
     &           ,'Block SUSHI, entry 6 out of range.')
         endif

         blocks(1) = 'DISTRIB'
         blocks(2) = 'SCALES'
         blocks(3) = 'VEGAS'
         blocks(4) = 'PDFSPEC'
         call readblocks(12,blocks)

         do i=1,max(norderggh+1,norderbbh+1)
            if (lpdfname(i)) then
                SUpdfs(i) = pdfname(i)
            else
               write(6,*) 'SusHi: PDFSPEC(',i,') undefined.'//
     &              ' Stopped.'
               stop
            endif
         enddo

         if (lpdfspecpar(10)) then
            if (lpdfspecpar(11).or.lpdfspecpar(12).or.lpdfspecpar(13)
     &           .or.lpdfspecpar(14)) then
               call printwarnsushi(verbose,6
     &              ,'Block PDFSPEC, entry 10 defined:')
               call printwarnsushi(verbose,6,
     &              ' => overriding possible entries 11-14!')
            endif
            SUiset(1) = pdfspecpar(10)
            SUiset(2) = pdfspecpar(10)
            SUiset(3) = pdfspecpar(10)
            SUiset(4) = pdfspecpar(10)
         else
            SUiset(1) = mapintdef(pdfspecpar(11),lpdfspecpar(11),0
     &           ,'PDF number set to 0.')
            SUiset(2) = mapintdef(pdfspecpar(12),lpdfspecpar(12),0
     &           ,'PDF number set to 0.')
            SUiset(3) = mapintdef(pdfspecpar(13),lpdfspecpar(13),0
     &           ,'PDF number set to 0.')
            SUiset(4) = mapintdef(pdfspecpar(14),lpdfspecpar(14),0
     &           ,'PDF number set to 0.')
         endif

         if (.not.lsushipar(3)) call printdie('SUSHI(3) undefined ')
         if(sushipar(3).eq.0) then
            ppcoll = .true.
            ppbarsign = +1
         elseif (sushipar(3).eq.1) then
            ppcoll = .false.
            ppbarsign = -1
         else
            call printdie('SUSHI(3) out of range. ')
         endif
         cme = 
     &        mapreal8(sushipar(4),lsushipar(4),'cm-energy undefined. ')

         !scale choices
         murfacggh = 
     &        mapreal8def(scalespar(1),lscalespar(1),1.d0,
     &        'Using default renormalization scale muR=mHiggs.')
         muffacggh = 
     &        mapreal8def(scalespar(2),lscalespar(2),murfacggh,
     &        'Using default factorization scale muF=muR.')
         murfacbbh = 
     &        mapreal8def(scalespar(11),lscalespar(11),murfacggh,
     &        'Using common renormalization scale for ggh and '
     &        //qqhstring//'.')
         muffacbbh = 
     &        mapreal8def(scalespar(12),lscalespar(12),muffacggh,
     &        'Using default muF=muR in '//qqhstring//'.')
         facmurgghmin = .1d0
         facmurgghmax = 10.d0
         facmurgghstep = 100
         facmurint1 = 0.5d0
         facmurint2 = 2.d0

         if (lmurscanpar(1)) then
            facmurint1 = frgghmin
            facmurint2 = frgghmax
         endif
         if (lmurscanpar(2)) then
            facmurgghmin = frgghmintab
            facmurgghmax = frgghmaxtab
            facmurgghstep = frgghstep
         endif
         if((facmurint1.le.0.d0).or.(facmurint2.le.0.d0))
     &        call printerrorsushi(1,6
     &        ,'Block SCALES, entry 101 out of range.')
         if((facmurgghmin.le.0.d0).or.(facmurgghmax.le.0.d0)
     &        .or.(facmurgghstep.le.0.d0))
     &        call printerrorsushi(1,6
     &        ,'Block SCALES, entry 102 out of range.')
         if(facmurint1.gt.facmurint2) then
            facmurint1 = frgghmax
            facmurint2 = frgghmin
         endif
c$$$         if(facmurgghmin.gt.facmurgghmax) then
c$$$            facmurgghmin = frgghmaxtab
c$$$            facmurgghmax = frgghmintab
c$$$         endif

         !scalespt determines whether muR = factor * Sqrt(mh^2+p_T^2)
         !should be used in case of p_T distributions.
         if (lscalespar(3)) then
            if (scalespar(3).eq.1) scalespt = .true.
         else
            scalespt = .false.
         endif

         if (lscalespar(4)) then
            muBfactor = scalespar(4)
         else
            muBfactor = 1.d0
         endif

         !distribution parameters
         dist = 0
         ptc = 0
         rapc = 0
         if (ldistribpar(1)) dist = distribpar(1)
         if (ldistribpar(2)) ptc = distribpar(2)
         if (ldistribpar(3)) rapc = distribpar(3)

c..   ltotinc = .true. : calculate total inclusive cross section
         ltotinc = .false.
         if ((dist.eq.0).and.(ptc.eq.0).and.(rapc.eq.0)) 
     &        ltotinc = .true.
         
         minptc = 0.d0
         minrapc = 0.d0
         if (ldistribpar(21)) minptc = distribpar(21)
         if (ldistribpar(22)) maxptc = distribpar(22)
         if (ldistribpar(31)) minrapc = distribpar(31)
         if (ldistribpar(32)) maxrapc = distribpar(32)

         if (ptc.ge.2) then
            maxptc = mapreal8(distribpar(22),ldistribpar(22),
     &           'max. value of pt undefined. ')
         endif
         if ((rapc.eq.1).or.(rapc.eq.3)) then
            maxrapc = mapreal8(distribpar(32),ldistribpar(32),
     &           'max. value of rapidity undefined. ')
         end if

         if ((dist.eq.1).or.(dist.eq.3)) then
            maxptc = mapreal8(distribpar(22),ldistribpar(22),
     &           'value of pt undefined. ')
         end if
         if ((dist.eq.2).or.(dist.eq.3)) then
            maxrapc = mapreal8(distribpar(32),ldistribpar(32),
     &           'value of rapidity undefined. ')
         end if
         
         pseudorap = .true.
         if((distribpar(4).eq.0).and.((dist.ge.2).or.(rapc.ne.0))) then
            pseudorap = .false.
         endif

         juser = .false.
c$$$         if(distribpar(5).eq.1) then
c$$$            juser = .true.
c$$$         endif

         acc=1.d-8
         lveg1sushi=.true.
         ncall1sushi=mapintdef(vegaspar(1),lvegaspar(1),10000,
     &        'Using default VEGAS ncall=10000.')
         itmx1sushi=mapintdef(vegaspar(2),lvegaspar(2),5,
     &        'Using default VEGAS itmx=5.')
         nprnvsushi=mapintdef(vegaspar(3),lvegaspar(3),10,
     &        'Using default VEGAS nprn=10.')
         ncall2sushi=mapintdef(vegaspar(11),lvegaspar(11),0,
     &        'VEGAS will be called only once.')
         if (ncall2sushi.gt.0) then
            itmx2sushi=mapintdef(vegaspar(12),lvegaspar(12),5,
     &           'Using default VEGAS itmx=5.')
         else
            lveg1sushi=.false.
            itmx2sushi=0
         endif

         lveg1ggh=.true.
         ncall1ggh=mapintdef(vegaspar(4),lvegaspar(4),2000,
     &        'Using default VEGAS ncall1=2000 for ggh.')
         itmx1ggh=mapintdef(vegaspar(5),lvegaspar(5),5,
     &        'Using default VEGAS itmx1=5 for ggh.')
         nprnvggh=mapintdef(vegaspar(6),lvegaspar(6),0,
     &        'Using default VEGAS nprnggh=0.')
         ncall2ggh=mapintdef(vegaspar(14),lvegaspar(14),5000,
     &        'Using default VEGAS ncall2=5000 for ggh.')
         if (ncall2ggh.eq.0) then
            lveg1ggh=.false.
            itmx2ggh=0
         else
            itmx2ggh=mapintdef(vegaspar(15),lvegaspar(15),2,
     &           'Using default VEGAS itmx2=2 for ggh.')
         endif
         
         lveg1bbh=.true.
         ncall1bbh=mapintdef(vegaspar(7),lvegaspar(7),2000,
     &        'Using default VEGAS ncall1=2000 for '//qqhstring//'.')
         itmx1bbh=mapintdef(vegaspar(8),lvegaspar(8),5,
     &        'Using default VEGAS itmx1=5 for '//qqhstring//'.')
         nprnvbbh=mapintdef(vegaspar(9),lvegaspar(9),0,
     &        'Using default VEGAS nprnbbh=0.')
         ncall2bbh=mapintdef(vegaspar(17),lvegaspar(17),5000,
     &        'Using default VEGAS ncall2=5000 for '//qqhstring//'.')
         if (ncall2bbh.eq.0) then
            lveg1bbh=.false.
            itmx2bbh=0
         else
            itmx2bbh=mapintdef(vegaspar(18),lvegaspar(18),2,
     &           'Using default VEGAS itmx2=2 for '//qqhstring//'.')
         endif

         if (verbose.eq.0) then
            nprnvsushi = 0
            nprnvggh = 0
            nprnvbbh = 0
         end if

c..   Variables for soft and 1/mt expansion:
c..
c..   nsoft(k,l): expand through xm1^nsoft in NkLO result for gg,qg,...
c..   nsoftmt(k,l): expand through xm1^nsoft in NkLO 1/mt terms for gg,qg,...
c..   nmtlim(k,l):  expand through 1/mt^nmtlim in NkLO result for gg,qg,...
c..   ncat1: do soft expansion of sigma1(x)/x^ncat1 in leading 1/mt terms
c..   ncat1mt: do soft expansion of sigma1(x)/x^ncat1 in 1/mt terms
c..   lball(1)=.true.: do matching to x->0 limit
c..
c..   nexpand(1)=0: no soft expansion
c..   nexpand(1)=1: needs all parameters defined above
c..   nexpand(1)=2: needs nsoft(1,0) and sets nsoft(1,1),... equal to nsoft(1,0)
c..   nexpand(1)=3: like nexpand(1)=2, 
c..               but sets 1/mt terms for pure quark channels =0
c..               (they do not converge well)
c.. 
c..   analogue parameters for NNLO and N3LO
      htlflag = 0
      gghmtpar=0
      gghmtreal=0.d0
      lgghmtpar=.false.
      lgghmtreal=.false.
      nfaclo=3
      nmtlim=0
      bltype = 'GGHMT'
      call readslha(12,bltype,blockflag)

      if (norderggh.ge.0) then
         nfaclo = mapintdef2(0,gghmtpar(-1),lgghmtpar(-1),-1,3,3
     &        ,'factor out exact LO result through'//
     &        ' LO(=0)/NLO(=1)/etc. [def=3]'
     &        ,'GGHMT(-1) not within allowed range. ')
         nmtlim(0,0) = mapintdef2(0,gghmtpar(0),lgghmtpar(0),-1,80000,-1
     &        ,'expand LO result through 1/mt^n [def=-1]'
     &        ,'GGHMT(0) not within allowed range. ')
      endif
      if (norderggh.ge.1) then
         nmtlim(1,0) = mapintdef2(0,gghmtpar(1),lgghmtpar(1),0,10,0
     &        ,'expand NLO result through 1/mt^n [def=-1]'
     &        ,'GGHMT(1) not within allowed range. ')
         nmtlim(1,1) = mapintdef2(0,gghmtpar(11),lgghmtpar(11),0,10
     &        ,nmtlim(1,0),'expand NLO gg result through 1/mt^n'
     &        ,'GGHMT(11) not within allowed range. ')
         nmtlim(1,2) = mapintdef2(0,gghmtpar(12),lgghmtpar(12),0,10
     &        ,nmtlim(1,0),'expand NLO qg result through 1/mt^n'
     &        ,'GGHMT(12) not within allowed range. ')
         nmtlim(1,3) = mapintdef2(0,gghmtpar(13),lgghmtpar(13),0,10
     &        ,nmtlim(1,0),'expand NLO qqbar result through 1/mt^n'
     &        ,'GGHMT(13) not within allowed range. ')
         nmtlim(1,0) = max(nmtlim(1,0),nmtlim(1,1),nmtlim(1,2),nmtlim(1
     &        ,3))
         lball(1)=.false.
         ndum = mapintdef2(0,gghmtpar(10),lgghmtpar(10),0,1,0
     &        ,'do not match NNLO 1/mt expansion to high energy limit'
     &        ,'GGHMT(10) not within allowed range. ')
         if (ndum.eq.1) lball(1)=.true.
      endif
      if (norderggh.ge.2) then
         nmtlim(2,0) = mapintdef2(0,gghmtpar(2),lgghmtpar(2),0,6,0
     &        ,'expand NNLO result through 1/mt^n [def=0]'
     &        ,'GGHMT(2) not within allowed range. ')
         nmtlim(2,1) = mapintdef2(0,gghmtpar(21),lgghmtpar(21),0,6
     &        ,nmtlim(2,0),'expand NNLO gg result through 1/mt^n'
     &        ,'GGHMT(21) not within allowed range. ')
         nmtlim(2,2) = mapintdef2(0,gghmtpar(22),lgghmtpar(22),0,6
     &        ,nmtlim(2,0),'expand NNLO qg result through 1/mt^n'
     &        ,'GGHMT(22) not within allowed range. ')
         nmtlim(2,3) = mapintdef2(0,gghmtpar(23),lgghmtpar(23),0,6
     &        ,nmtlim(2,0),'expand NNLO qqbar result through 1/mt^n'
     &        ,'GGHMT(23) not within allowed range. ')
         nmtlim(2,4) = mapintdef2(0,gghmtpar(24),lgghmtpar(24),0,6
     &        ,nmtlim(2,0),'expand NNLO qq result through 1/mt^n'
     &        ,'GGHMT(24) not within allowed range. ')
         nmtlim(2,5) = mapintdef2(0,gghmtpar(25),lgghmtpar(25),0,6
     &        ,nmtlim(2,0),'expand NNLO qqprime result through 1/mt^n'
     &        ,'GGHMT(25) not within allowed range. ')
         nmtlim(2,0) = max(nmtlim(2,0),nmtlim(2,1),nmtlim(2,2),nmtlim(2
     &        ,3),nmtlim(2,4),nmtlim(2,5))
         lball(2)=.false.
         ndum = mapintdef2(0,gghmtpar(20),lgghmtpar(20),0,1,0
     &        ,'do not match NNLO 1/mt expansion to high energy limit'
     &        ,'GGHMT(20) not within allowed range. ')
         if (ndum.eq.1) lball(2)=.true.
      endif
      if (norderggh.ge.3) then
         lball(3)=.false.
         ndum = mapintdef2(0,gghmtpar(30),lgghmtpar(30),0,1,0
     &        ,'do not match N3LO result to high energy limit'
     &        ,'GGHMT(30) not within allowed range. ')
         if (ndum.eq.1) lball(3)=.true.
      endif
      htlflag = 0
      delmatch = mapreal8def(gghmtreal(9),lgghmtreal(9),5d-2,'')
         
      nsoft=0
      nsoftmt=0

      gghsoftpar=0
      lgghsoftpar=.false.
      nexpand(1)=0
      ncat1=0
      nexpandmt(1)=0
      ncat1mt=0
      nexpand(2)=0
      ncat2=0
      nexpandmt(2)=0
      ncat2mt=0
      nexpand(3)=0
      ncat3=0
      nexpandmt(3)=0
      ncat3mt=0

      bltype = 'GGHSOFT'
      call readslha(12,bltype,blockflag)
      if (norderggh.ge.1) then
         nexpand(1)=mapintdef2(0,gghsoftpar(1,1),lgghsoftpar(1),
     &        0,1,0,'[0/1]: soft expansion at NLO [no/yes]'
     &        ,'GGHSOFT(1,1) not within allowed range. ')
c$$$         if (nmtlim1.gt.0) then
c$$$            nexpandmt(1)=mapintdef2(0,gghsoftpar(11,1),lgghsoftpar(11),
c$$$     &           0,1,nexpand(1),'[0/1]: soft exp of NLO 1/mt [no/yes]'
c$$$     &           ,'GGHSOFT(11,1) not within allowed range. ')
c$$$         endif
      endif
      if (nexpand(1).gt.0) then
         nsoft(1,0)=mapintdef2(0,gghsoftpar(1,2),lgghsoftpar(1),-1,16,16
     &        ,'soft exp of NLO through (1-x)^16'
     &        ,'GGHSOFT(1,2) not within allowed range. ')
         ncat1=mapintdef2(0,gghsoftpar(1,3),lgghsoftpar(1),0,1000,0
     &        ,'expand sigma(x)/x^0 at NLO'
     &        ,'GGHSOFT(1,3) not within allowed range. ')
         nsoftmt(1,0)=nsoft(1,0)
         ncat1mt=ncat1
      endif
c$$$      if (nexpandmt(1).gt.0) then
c$$$         nsoft1mt=mapintdef2(0,gghsoftpar(11,2),lgghsoftpar(11),-1,16
c$$$     &        ,nsoft1,'Soft exp of NLO 1/mt through (1-x)^'
c$$$     &        //tostring(nsoft1)
c$$$     &        ,'GGHSOFT(11,2) not within allowed range. ')
c$$$         ncat1mt=mapintdef2(0,gghsoftpar(11,3),lgghsoftpar(11),0,1000
c$$$     &        ,ncat1,'expand sigma(x)/x^'//tostring(ncat1)/
c$$$     &        /'at NLO 1/mt','GGHSOFT(11,3)not within allowed range. ')
c$$$      endif
      if (norderggh.ge.2) then
         nexpand(2)=mapintdef2(0,gghsoftpar(2,1),lgghsoftpar(2),
     &        0,1,0,'[0/1]: soft exp of NNLO [no/yes]'
     &        ,'GGHSOFT(2,1) not within allowed range. ')
c$$$         if (nmtlim2.gt.0) then
c$$$            nexpandmt(2)=mapintdef2(0,gghsoftpar(12,1),lgghsoftpar(12),
c$$$     &           1,1,1,'[0/1]: soft exp of NNLO 1/mt [no/yes]'
c$$$     &           ,'GGHSOFT(12,1) not within allowed range. ')
c$$$         endif
      endif
      if (nexpand(2).gt.0) then
         nsoft(2,0)=mapintdef2(0,gghsoftpar(2,2),lgghsoftpar(2),-1,16,16
     &        ,'Soft exp of NNLO through (1-x)^16'
     &        ,'GGHSOFT(2,2) not within allowed range. ')
         ncat2=mapintdef2(0,gghsoftpar(2,3),lgghsoftpar(2),0,1000,0
     &        ,'expand sigma(x)/x^0 at NNLO'
     &        ,'GGHSOFT(2,3) not within allowed range. ')
         nsoftmt(2,0)=nsoft(2,0)
         ncat2mt=ncat2
      endif
c$$$      if (nexpandmt(2).gt.0) then
c$$$         nsoft2mt=mapintdef2(0,gghsoftpar(12,2),lgghsoftpar(12),-1,16
c$$$     &        ,nsoft2,'Soft exp of NNLO 1/mt through (1-x)^'
c$$$     &        //tostring(nsoft2)
c$$$     &        ,'GGHSOFT(12,2) not within allowed range. ')
c$$$         ncat2mt=mapintdef2(0,gghsoftpar(12,3),lgghsoftpar(12),0,1000
c$$$     &        ,ncat2,'expand sigma(x)/x^'//tostring(ncat2)/
c$$$     &        /' at NNLO 1/mt'
c$$$     &        ,'GGHSOFT(12,3) not within allowed range. ')
c$$$      endif
      if (norderggh.ge.3) then
         nexpand(3)=mapintdef2(0,gghsoftpar(3,1),lgghsoftpar(3),
     &        1,1,1,'[0/1]: soft exp of NNLO [no/yes]'
     &        ,'GGHSOFT(3,1) not within allowed range. ')
c$$$         if (nmtlim3.gt.0) then
c$$$            nexpandmt(3)=mapintdef2(0,gghsoftpar(13,1),lgghsoftpar(13),
c$$$     &           0,1,nexpand(3),'[0/1]: soft exp of NNLO 1/mt [no/yes]'
c$$$     &           ,'GGHSOFT(13,1) not within allowed range. ')
c$$$         endif
      endif
      if (nexpand(3).gt.0) then
         nsoft(3,0)=mapintdef2(0,gghsoftpar(3,2),lgghsoftpar(3),-1,16,16
     &        ,'Soft exp of NNLO through (1-x)^16'
     &        ,'GGHSOFT(3,2) not within allowed range. ')
         ncat3=mapintdef2(0,gghsoftpar(3,3),lgghsoftpar(3),0,1000,0
     &        ,'expand sigma(x)/x^0 at NNLO'
     &        ,'GGHSOFT(3,3) not within allowed range. ')
         nsoftmt(3,0)=nsoft(3,0)
         ncat3mt=ncat3
      endif
c$$$      if (nexpandmt(3).gt.0) then
c$$$c..   nsoft3mt has no effect yet:
c$$$         nsoft3mt=mapintdef2(0,gghsoftpar(13,2),lgghsoftpar(13),-1,16
c$$$     &        ,nsoft3,'Soft exp of NNLO 1/mt through (1-x)^'
c$$$     &        //tostring(nsoft3)
c$$$     &        ,'GGHSOFT(13,2) not within allowed range. ')
c$$$c..   ncat3mt has no effect yet:
c$$$         ncat3mt=mapintdef2(0,gghsoftpar(13,3),lgghsoftpar(13),0,1000
c$$$     &        ,ncat3,'expand sigma(x)/x^'//tostring(ncat3)/
c$$$     &        /' at NNLO 1/mt'
c$$$     &        ,'GGHSOFT(13,3) not within allowed range. ')
c$$$      endif

      endif
      close(12)

 105  format('#--------------------------------------------------#')
      end

C-}}}
C-{{{ subroutine outputsushi

      subroutine outputsushi()

      implicit none

      integer i,count1,count2
      logical ppcoll
      double precision effres, effreserr, effres2, effres2err
      character effresstring*30
      character effresstring2*30
      integer date_time(8)
      character realcl(3)*12,tmpstring*4,orderstring(4)*4,
     &     fhvers*20,hbvers*20,hsvers*20,thdmcvers*20
      common /vers/ fhvers,hbvers,hsvers,thdmcvers

      include '../commons/common-inputoutput.f'
      !internal definitions

!defined at other places:
!AlfasMZ,muRfacggh,muFfacggh,cme,dist,prefix,prefix2,pty,iset
!norderggh,mw,mc,ew,mz,mbmb,renscheme,mt,pseudo,tanbresum,pseudorap
      include '../commons/common-slha.f'
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-quark.f'
      include '../commons/common-int.f'
      include '../commons/common-citations.f'
      include '../commons/common-keys.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-bbhnnlo.f'
      include '../commons/common-flags.f'
      include '../commons/common-ptnames.f'
      include '../commons/common-lumis.f'

      common /coll/ ppcoll

      open(unit=13,file=jfnameact,status='unknown')

      orderstring(1) = 'LO  '
      orderstring(2) = 'NLO '
      orderstring(3) = 'NNLO'
      orderstring(4) = 'N3LO'

      call printcopysushi(0,13,verbose)

      call date_and_time(realcl(1),realcl(2),realcl(3),date_time)

      write (13,106) date_time(3), date_time(2), date_time(1),
     &date_time(5), date_time(6), date_time(7)

      write(13,103)
     &"# For the reconstructed input file, check below.   #"
      write(13,103)
     &"#--------------------------------------------------#"
      write(13,"(a)")
     &"# Please cite the following papers (for this run):"
      call printcite(13,1,'# ')
      write(13,103)
     &"#--------------------------------------------------#"

      !most important results
      if ((dist.eq.0).or.(dist.eq.2)) then
         effresstring = "XS in pb"
      else
         effresstring = "XS in pb/GeV"
      endif

      effres = xeff
      effreserr = xefferr
      if (norderggh.eq.0) then
         effres = xsec(1)
         effreserr = xsecerr(1)
      endif

!xeffd,xeffu,murgghd,murgghu

      write(13,103) 'Block SUSHIinfo'
      write(13,1023) 1,SUversion,'SusHi version'
      if(fhflag.and.(.not.extflag)) then
         write(13,1023) 2,fhvers,'FeynHiggs version'
         if(citations(31).gt.0) then
            write(13,1023) 3,hbvers,'HiggsBounds version'
            if(citations(34).gt.0) then
               write(13,1023) 4,hsvers,'HiggsSignals version'
            endif
         endif
      endif
      if(thdmcflag) then
         write(13,1023) 5,thdmcvers,'2HDMC version'
      endif

      if (norderggh.ge.0) then
         write(13,103) 
     &        "Block SUSHIggh # Bon appetit "
          write(13,102) 1,effres,'ggh '//effresstring
          write(13,102) 101,effreserr,
     &    '+/- integ. error: ggh '//effresstring
          if (muranalytic) then
          write(13,102) 102,gghXSvar(1),
     &    '- from muR variation: ggh '//effresstring
          write(13,102) 103,gghXSvar(2),
     &    '+ from muR variation: ggh '//effresstring
          write(13,1030) '# within [',dexp(xeffvar(1))/murggh,
     &  ',',dexp(xeffvar(2))/murggh,']*SCALES(1)*MASSOUT(1).'
          end if

          if ( ((dist.eq.1).or.(dist.eq.3).or.ptcut) .and.
     &           (ptref.lt.20.d0)) then
               write(13,103) 
     &        "# Note: Fixed order p_T distribution (no resummation)!"
               write(13,103) "# p_T < 20 GeV not reasonable!"
          endif
      endif
      
      if (norderbbh.eq.0) then
         effres2 = SUbbh(1)
         effres2err = SUbbherr(1)
      else if (norderbbh.ge.0) then
         effres2 = SUbbh(2)
         effres2err = SUbbherr(2)
      end if
      if (norderbbh.eq.2) then
         if ((.not.((ptc.gt.0).or.(rapc.gt.0))).and.(dist.eq.0)) then
            effres2 = SUbbh(3)
            effres2err = SUbbherr(3)
         end if
      end if
      effresstring2 = qqhstring//' '//effresstring

      if ((norderbbh.ge.0).and.(dist.le.1)) then
         write(13,103) 
     &        "Block SUSHIbbh # Bon appetit "
         if (lqqhkey) write(13,103) '# calculated '//ptname(ptn1)
     &        //ptname(ptn2)//'->H xsec - see Block QQH below'
         write(13,102) 1,effres2,effresstring2
         write(13,102) 101,effres2err,
     &   '+/- integr. error: '//effresstring2
         if (((dist.eq.1).or.(dist.eq.3).or.ptcut) .and.
     &        (ptref.lt.20.d0)) then
            write(13,103) 
     &        "# Note: Fixed order p_T distribution (no resummation)!"
            write(13,103) "# p_T < 20 GeV not reasonable!"
         endif
      endif
     
      !detailled output
      !SM ggh
      if (norderggh.ge.0) then
         if ((dist.eq.0).or.(dist.eq.2)) then
            write(13,103) 'Block XSGGH # ggh xsec in pb (w/o EW) '//
     &           '(+/- integ.err.)'
         else
            write(13,103)
     &           'Block XSGGH # ggh xsec in pb/GeV  (w/o EW) '//
     &           '(+/- integ.err.)'
         endif
         if (norderggh.eq.0) then
            write(13,102) 1,xsec(1),"LO"
         elseif (norderggh.ge.1) then
            write(13,102) 1,xsec(1),"LO w/ NLO PDFs"
            write(13,102) 2,xsec(2),"NLO"
            write(13,102) 21,ggchan,"NLO gg"
            write(13,102) 22,qgchan,"NLO qg"
            write(13,102) 23,qqchan,"NLO qq"
            write(13,102) 101,xsecerr(1),"+/-: LO w/ NLO PDFs"
            write(13,102) 102,xsecerr(2),"+/-: NLO"
            write(13,102) 121,ggchanerr,"+/-: NLO gg"
            write(13,102) 122,qgchanerr,"+/-: NLO qg"
            write(13,102) 123,qqchanerr,"+/-: NLO qq"
         endif
         write(13,103) 'Block XSGGHEFF # ggh xsec in heavy top limit '
     &        //'(from ggh@nnlo) (+/- integ.err.)'
         if (ltotinc) then
            if (norderggh.ge.1) then
               write(13,102) 1,SUggh(2,2,0),"NLO"
               write(13,102) 101,SUggherr(2,2,0),"+/-: NLO"
               write(13,102) 19,SUggh(norderggh+1,2,0),"NLO'"
               write(13,102) 191,SUggherr(norderggh+1,2,0),
     &              "+/-: NLO'"
               if (nsubprocggh.eq.10) then
!     ggh@nnlo computed all subchannels separately:
                  write(13,102) 10,SUggh(2,2,-1),"NLO: soft+virt"
                  write(13,102) 11,SUggh(2,2,1),"NLO: gg-channel"
                  write(13,102) 12,SUggh(2,2,2),"NLO: qg-channel"
                  write(13,102) 13,SUggh(2,2,3),"NLO: qqb-channel"
                  write(13,102) 110,SUggherr(2,2,-1)
     &                 ,"+/- : NLO soft+virt"
                  write(13,102) 111,SUggherr(2,2,1),"+/- : NLO gg"
                  write(13,102) 112,SUggherr(2,2,2),"+/- : NLO qg"
                  write(13,102) 113,SUggherr(2,2,3),"+/- : NLO qqb"
               endif
            endif
            do i=3,4
               if (norderggh.ge.(i-1)) then
                tmpstring = orderstring(i)
                write(13,102) i-1,SUggh(i,i,0),tmpstring
                write(13,102) 100+i-1,SUggherr(i,i,0),"+/-: "//tmpstring
                if (nsubprocggh.eq.10) then
!     ggh@nnlo computed all subchannels separately:
                     write(13,102) 10*(i-1),SUggh(i,i,-1),tmpstring/
     &                    /": soft+virt"
                     write(13,102) 10*(i-1)+1,SUggh(i,i,1),tmpstring/
     &                    /": gg-channel"
                     write(13,102) 10*(i-1)+2,SUggh(i,i,2),tmpstring/
     &                    /": qg-channel"
                     write(13,102) 10*(i-1)+3,SUggh(i,i,3),tmpstring/
     &                    /": qqb-channel"
                     write(13,102) 10*(i-1)+4,SUggh(i,i,4),tmpstring/
     &                    /": qq-channel"
                     write(13,102) 10*(i-1)+5,SUggh(i,i,5),tmpstring/
     &                    /": qu-channel"
                     write(13,102) 100+10*(i-1),SUggherr(i,i,-1)
     &                    ,"+/- : soft+virt"
                     write(13,102) 100+10*(i-1)+1,SUggherr(i,i,1),
     &                    "+/- : "//tmpstring//" gg"
                     write(13,102) 100+10*(i-1)+2,SUggherr(i,i,2),
     &                    "+/- : "//tmpstring//" qg"
                     write(13,102) 100+10*(i-1)+3,SUggherr(i,i,3),
     &                    "+/- : "//tmpstring//" qqb"
                     write(13,102) 100+10*(i-1)+4,SUggherr(i,i,4),
     &                    "+/- : "//tmpstring//" qq"
                     write(13,102) 100+10*(i-1)+5,SUggherr(i,i,5),
     &                    "+/- : "//tmpstring//" qu"
                  endif
               endif
            enddo
         endif
         write(13,102) 4,elw,"electroweak factor"
      endif
      
      !bbh
      if (norderbbh.ge.0) then
         if ((dist.eq.0).or.(dist.eq.1)) then
            if (dist.eq.0) then
            write(13,103) 'Block XSBBH # '//qqhstring/
     &           /' xsec in pb (+/- integ.err.)'
            else
            write(13,103) 'Block XSBBH # '//qqhstring/
     &           /' xsec in pb/GeV (+/- integ.err.)'
            end if
            if (lqqhkey) write(13,103) '# calculated '//ptname(ptn1)
     &           //ptname(ptn2)//'->H xsec - see Block QQH below'
            write(13,102) 1,SUbbh(1),"LO"
            write(13,102) 101,SUbbherr(1),"+/-: LO"
            if (norderbbh.gt.0) then
               write(13,102) 2,SUbbh(2),"NLO"
               write(13,102) 102,SUbbherr(2),"+/-: NLO"
            end if
            if (norderbbh.gt.1) then
               if ((ptc.eq.0).and.(rapc.eq.0).and.(dist.eq.0)) then
                  write(13,102) 3,SUbbh(3),"NNLO"
                  write(13,102) 103,SUbbherr(3),"+/-: NNLO"
               end if
            end if
         endif
      endif

      call writecouplings()

      write(13,103) 'Block ALPHASOUT # values from LHAPDF set'
      if (scalespt) then
         if (norderggh.ge.1) then
            write(13,102) 2,awithpt,
     &           'alpha_s(muR)@NLO'
            write(13,102) 20,awopt,
     &           'alpha_s(muR(w/o)pt)@NLO'
         endif
      else
         if (norderggh.eq.0) then
            write(13,102) 1,pi*apiggh(1),'alpha_s(muRggh) @ LO'
         endif
         if (norderggh.ge.1) then
            write(13,102) 2,pi*apiggh(2),'alpha_s(muRggh) @ NLO'
         endif
         if (norderbbh.ge.0) then
            write(13,102) 11,apimurbbh*pi,'alpha_s(muRbbh)'
         endif
         
         if (norderggh.ge.2) then
            write(13,102) 3,pi*apiggh(3),'alpha_s(muRggh) @ NNLO'
         end if
         if (norderggh.ge.3) then
            write(13,102) 4,pi*apiggh(4),'alpha_s(muRggh) @ N3LO'
         end if
      end if
      do i=0,max(norderggh,norderbbh)
         write(13,102) 101+i,apimzlist(i+1)*pi,'alpha_s(Mz) @ '
     &        //orderstring(i+1)
      enddo

      if(citations(31).gt.0) then
         write(13,103) 'Block HiggsBoundsResults # HiggsBounds output'
         write(13,1021) 1,HB_chan,'most sensitive channel '//
     &        '(see Key.dat)'
         write(13,1021) 2,HB_result,'result flag (1: allowed, 0: '//
     &        'excluded, -1: unphysical)'
         write(13,102) 3,HB_obsratio,'ratio [sig x BR] model/[sig x '//
     &        'BR] limit (< 1: allowed, > 1: excluded)'
         write(13,1021) 4,HB_ncombined,'number of Higgs bosons '//
     &        'combined for most sensitive channel'
      endif

      if(citations(34).gt.0) then
         write(13,103) 'Block HiggsSignalsResults # HiggsSignals output'
         write(13,1021) 2,1,'Chi-squared method (1:peak-c, '//
     &        '2:mass-c, 3:both)'
         write(13,1021) 3,2,'Higgs mass pdf (1:box, 2:Gaussian, '//
     &        '3:box+Gaussian)'
         write(13,1021) 7,HS_nobs,'Total number of observables'
         write(13,102) 8,HS_Chisq_mu,'Chi^2 contr. from signal '//
     &        'strength observables'
         write(13,102) 9,HS_Chisq_mh,'Chi^2 contr. from mass '//
     &        'observables'
         write(13,102) 12,HS_Chisq,'Total Chi^2'
         write(13,102) 13,HS_Pvalue,'Probability (total Chi^2, '//
     &        'total number observables)'
      endif

      write(13,103)
     &     "#--------------------------------------------------#"
      write(13,103)
     &     "# Reconstructed input file of this run             #"
      write(13,103)
     &     "#--------------------------------------------------#"
      call writeinpfile(13,cme,scalespt,0.d0,muFfacggh,SUpdfs,SUiset)

      !call filedump(12,13)
      write(13,103)
     &"#--------------------------------------------------#"
      write(13,103) 
     &"# End of SusHi output                              #"
      write(13,103)
     &"#--------------------------------------------------#"

      close(13)

 102  format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
c 1020 format(1x,I9,3x,1P,i2.1,0P,3x,'#',1x,A)
 1021 format(1x,I9,3x,1P,i6.1,0P,3x,'#',1x,A)
c 1022 format(1x,I9,3x,1P,A45,0P,3x,'#',1x,A)
 1023 format(1x,I9,3x,1P,A10,0P,3x,'#',1x,A)
c 1024 format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,I1,A)
 103  format(a)
 1030  format(A,E16.8,A,E16.8,A)
c 1030 format(a,I1,a)
 106  format ('# Date: ',i2.2,'.',i2.2,'.',i4.4,' at ',
     &i2.2,':',i2.2,':',i2.2,'                     #')
c 107  format(9x,1P,E16.8,0P,3x,'#',1x,A)
c 1070 format(9x,1P,i2.1,0P,3x,'#',1x,A)
c 108  format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)

      end

C-}}}
C-{{{ subroutine writecouplings

      subroutine writecouplings()
      implicit none
      double precision vev,vevS,Sqrt2
      character psstr*7
      integer Sind2

      include '../commons/common-ren.f'
      include '../commons/common-inputoutput.f'
      include '../commons/common-flags.f'
      include '../commons/common-vars.f'
      include '../commons/common-quark.f'
      include '../commons/common-keys.f'

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam

      if (pseudo.eq.0) then
         psstr = 'CP-even'
         Sind2 = Sind
      else
         psstr = ' CP-odd'
         Sind2 = Sind - 1 
      end if

!     effective couplings in (N)MSSM
      if (model.ge.1) then
         write(13,1030) 
     &        "Block HGGSUSY # couplings of ",Sind2,". "//psstr/
     &        /" H to 3. gen."
         write(13,102) 101,gt3,'g_t^phi'
         if ((model.eq.1).or.(model.eq.3)) then
            write(13,102) 111,gt13,'g_st11^phi'
            write(13,102) 122,gt23,'g_st22^phi'
            write(13,102) 112,gth123,'g_st12^phi'
            write(13,102) 121,gth213,'g_st21^phi'
         end if
         write(13,102) 201,gb3,'g_b^phi'
         if ((model.eq.1).or.(model.eq.3)) then
            write(13,102) 211,gb13,'g_sb11^phi'
            write(13,102) 222,gb23,'g_sb22^phi'
            write(13,102) 212,gbh123,'g_sb12^phi'
            write(13,102) 221,gbh213,'g_sb21^phi'
         end if
         if((yukfac(4).ne.0.d0).or.(yukfac(5).ne.0.d0)) then
            write(13,1030) 
     &           "Block HGGSUSY # couplings of ",Sind2,". "//psstr/
     &           /" H to 4. gen."
            write(13,102) 101,gt4,'g_t^phi'
            if ((model.eq.1).or.(model.eq.3)) then
               write(13,102) 111,gt14,'g_st11^phi'
               write(13,102) 122,gt24,'g_st22^phi'
               write(13,102) 112,gth124,'g_st12^phi'
               write(13,102) 121,gth214,'g_st21^phi'
            end if
            write(13,102) 201,gb4,'g_b^phi'
            if ((model.eq.1).or.(model.eq.3)) then
               write(13,102) 211,gb14,'g_sb11^phi'
               write(13,102) 222,gb24,'g_sb22^phi'
               write(13,102) 212,gbh124,'g_sb12^phi'
               write(13,102) 221,gbh214,'g_sb21^phi'
            end if
         endif
         if ((model.eq.1).or.(model.eq.3)) then
            write(13,103)
     &"Block BBHREWEIGHT # top and bottom reweighting factors"
         write(13,102) 1,gt3,'g_t'
         if (pseudo.eq.1) then
            write(13,102) 2,gb3* (1.d0 - (1.d0/tanb
     &           **2 - Amx(Sind,3)*vev/(Amx(Sind,2)*vevS*tanb))*delmb)
     &           ,'g_b incl. 1/(1+Delta_b)*(1+C*Delta_b)'
         else if (pseudo.eq.0) then
            write(13,102) 2,gb3* (1.d0 + (Hmx(Sind,2)/(tanb
     &           *Hmx(Sind,1)) + Hmx(Sind,3)*vev*dcos(datan(tanb))
     &             /(Hmx(Sind,1)*vevS))*delmb)* (1.d0/(1.d0+delmb))
     &           ,'g_b incl. 1/(1+Delta_b)*(1+C*Delta_b)'
         end if
         end if
         if (model.eq.2) then
         write(13,103) 
     &"Block BBHREWEIGHT # top and bottom reweighting factors"
         write(13,102) 1,gt3,'g_t'
         write(13,102) 2,gb3,'g_b'
         end if
      endif
!additional output
      if (ldim5) then
         write(13,103) 'Block DIM5OUT'
            write(13,102) htype,c5sp(pseudo,0)
     &           ,'LO coeff of dim-5 operator'
            write(13,102) 100+htype,c5sp(pseudo,1)
     &           ,'NLO coeff of dim-5 operator'
            write(13,102) 200+htype,c5sp(pseudo,2)
     &           ,'NNLO coeff of dim-5 operator'
         if (pseudo.eq.0) then
            write(13,102) 300+htype,c5sp(pseudo,3)
     &           ,'N3LO coeff of dim-5 operator'
         end if
      endif
      write(13,103) 'Block MASSOUT'
      write(13,102) 1,Mh,'Mphi in GeV'
      if(yukfac(1).ne.0.d0) write(13,102) 4,mcmc,'m_c(m_c), MSbar'
      write(13,102) 5,mbmb,'m_b(m_b), MSbar'
      write(13,102) 6,mt,'m_t(pole)'
      write(13,102) 23,mZ,'m_Z'
      write(13,102) 24,mW,'m_W'
      write(13,1024) 25+10*(Sind-1)+pseudo,Mh,
     &       Sind2,'. '//psstr//' Higgs mass in GeV'
      if ((model.eq.1).or.(model.eq.3)) then
         write(13,102) 1000005,dsqrt(msb12),'sbottom1 mass in GeV'
         write(13,102) 2000005,dsqrt(msb22),'sbottom2 mass in GeV'
         write(13,102) 1000006,dsqrt(mst12),'stop1 mass in GeV'
         write(13,102) 2000006,dsqrt(mst22),'stop2 mass in GeV'
         write(13,102) 1000021,Mgl,'gluino mass in GeV'
      end if

      if ((model.eq.2).and.(thdmcflag)) then
         write(13,103) 'Block 2HDMCOUT # 2HDMC results for'
         write(13,1020) 1,int(thdmcres(7))
     &        ,'stability, 0=failed, 1=ok'
         write(13,1020) 2,int(thdmcres(8)),
     &        'perturbativity, 0=failed, 1=ok'
         write(13,1020) 3,int(thdmcres(9))
     &        ,'unitarity, 0=failed, 1=ok'
         write(13,102) 7,thdmcres(11),'h width in GeV'
         write(13,102) 8,thdmcres(12),'H width in GeV'
         write(13,102) 9,thdmcres(13),'A width in GeV'      
      end if
      if ((model.eq.1).and.(fhflag)) then
         write(13,103) 'Block ALPHA # Higgs mixing parameter from FH'
         write(13,107) alpha,'alpha'
      end if

      if ((model.eq.1).or.(model.eq.3)) then
         write(13,103) 'Block STOPMIXOUT # stop mixing matrix'
         write(13,108) 1,1,cthetat,'V_11'
         write(13,108) 1,2,sthetat,'V_12'
         write(13,108) 2,1,-sthetat,'V_21'
         write(13,108) 2,2,cthetat,'V_22'
         write(13,103) 'Block SBOTMIXOUT # sbottom mixing matrix'
         write(13,108) 1,1,cthetab,'V_11'
         write(13,108) 1,2,sthetab,'V_12'
         write(13,108) 2,1,-sthetab,'V_21'
         write(13,108) 2,2,cthetab,'V_22'
         write(13,103) 'Block AD'
         write(13,108) 3,3,Ab
     &        ,'used Ab in GeV - def. accord. to scheme'
         write(13,103) 'Block AU'
         write(13,108) 3,3,At,'used At in GeV'
      endif
      
      write(13,103) 'Block INTERNALMASSES # Masses in GeV'
      if(yukfac(1).ne.0.d0) write(13,102) 40,mcos,'m_c(pole)'
      write(13,102) 50,mbmb,'m_b(m_b), MSbar'
      write(13,102) 51,mbMSbarmuR ,'m_b(mu_R) MSbar'
      write(13,102) 52,mbOS       ,'m_b(pole)'
      write(13,102) 53,mb
     &     ,'m_b used for internal masses'
      if (tanbresum.eq.1) then
         write(13,102) 54,mbyuk/(1+delmb)
     &        ,'m_b for bottom Yukawa 1'
      else if (tanbresum.eq.2) then
         if (pseudo.eq.1) then
            write(13,102) 54,mbyuk* (1.d0 - (1.d0/tanb
     &           **2 - Amx(Sind,3)*vev/(Amx(Sind,2)*vevS*tanb))*delmb)
     &             * (1.d0/(1.d0+delmb))
     &           ,'m_b for bottom Yukawa 2'
         else if (pseudo.eq.0) then
            write(13,102) 54,mbyuk* (1.d0 + (Hmx(Sind,2)/(tanb
     &           *Hmx(Sind,1)) + Hmx(Sind,3)*vev*dcos(datan(tanb))
     &             /(Hmx(Sind,1)*vevS))*delmb)* (1.d0/(1.d0+delmb))
     &           ,'m_b for bottom Yukawa 2'
         end if 
      else
         write(13,102) 54,mbyuk
     &           ,'m_b for bottom Yukawa'
      endif
      if (model.eq.1) then
         write(13,102) 55,mbsb     ,'m_b for sbottom sector'
      endif

 102  format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 1020 format(1x,I9,3x,1P,i2.1,0P,3x,'#',1x,A)
c 1021 format(1x,I9,3x,1P,i6.1,0P,3x,'#',1x,A)
c 1022 format(1x,I9,3x,1P,A45,0P,3x,'#',1x,A)
 1024 format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,I1,A)
 103  format(a)
 1030 format(a,I1,a)
 107  format(9x,1P,E16.8,0P,3x,'#',1x,A)
c 1070 format(9x,1P,i2.1,0P,3x,'#',1x,A)
 108  format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
      end

C-}}}
C-{{{ subroutine writeinpfile

      subroutine writeinpfile(file,
     &     shad,runscale,qres,mufout,pdfnames,pdfset)
      implicit none
      integer file,i,bool2int
      logical runscale
      double precision shad,qres,mufout
      character pdfnames(1:4)*50
      character ordstr(4)*4
      integer pdfset(1:4)
      logical ppcoll
      integer count1,count2,Sind2
      character htypestring(0:3)*30
      character psstr*7,resstring*30

      include '../commons/common-inputoutput.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-int.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-quark.f'
      include '../commons/common-flags.f'
      include '../commons/common-expand.f'
      include '../commons/common-lumis.f'
      include '../commons/common-ptnames.f'

      common /coll/ ppcoll

      ordstr(1) = 'LO  '
      ordstr(2) = 'NLO '
      ordstr(3) = 'NNLO'
      ordstr(4) = 'N3LO'

      if (pseudo.eq.0) then
         psstr = 'CP-even'
         Sind2 = Sind
      else
         psstr = ' CP-odd'
         Sind2 = Sind - 1 
      end if

!     Reconstruction of the input file
      write(file,103) 'Block SUSHI'
      write(file,1020) 1,model,
     &     'Chosen model: 0=SM, 1=MSSM, 2=2HDM, 3=NMSSM'
      htypestring(0) = '11=scalar, 21=pseudo-scalar'
      htypestring(1) = '11=h, 12=H, 21=A'
      htypestring(2) = '11=h, 12=H, 21=A'
      htypestring(3) = '11=H1,12=H2,13=H3,21=A1,22=A2'
      write(file,1020) 2,htype,htypestring(model)

      if (ppcoll) then
         write(file,1020) 3,0,'Particle collider: 0=pp, 1=ppbar'
      else
         write(file,1020) 3,1,'Particle collider: 0=pp, 1=ppbar'
      end if
      write(file,102)  4,shad,'center-of-mass energy in GeV'
      if(.not.extflag) then
         write(file,1020) 5,norderggh,'Order for ggh'
         write(file,1020) 6,norderbbh,'Order for '//qqhstring
      else
         write(file,1020) 5,1,'Order for ggh'
         write(file,1020) 6,-1,'Order for '//qqhstring
      endif
      write(file,1020) 7,ew,'Electroweak contributions to ggh'
      write(file,1020) 11,bool2int(lexpandc1o1),
     &     '[1/0] do/not expand C1 perturbatively'
      write(file,1020) 19,verbose
     &     ,'0 = silent mode of SusHi, 1 = normal output'
      if(.not.extflag) then
         write(file,1020) 20,nsubprocggh
     &     ,'ggh@nnlo subprocesses: 0=all, 10=ind. contributions'
         write(file,1020) 21,nsubprocbbh,'bbh@nnlo subprocesses: 0=all'
      endif

      if ((norderggh.ge.0).and.(.not.extflag)) then
         call writegghmtblock(13)
      endif
         
      if ((model.eq.2).and.(thdmcflag.eqv.(.false.))) then
         write(file,103) 'Block 2HDM'
         write(file,1070) twohdmver,'2HDM type'
      end if
      if (thdmcflag) then
         write(file,103) 'Block 2HDMC'
         write(file,1020) 1,thdmckey,'parametrization of 2HDM'
         write(file,1020) 2,twohdmver,'2HDM type'
         write(file,1020) 10,thdmcbreak
     &        ,'[0/1] ignore validity checks [yes/no]'
         if (thdmckey.eq.1) then
            write(file,102) 3,tanb,'tan(beta)'
            write(file,102) 11,tparam(1),'lambda1'
            write(file,102) 12,tparam(2),'lambda2'
            write(file,102) 13,tparam(3),'lambda3'
            write(file,102) 14,tparam(4),'lambda4'
            write(file,102) 15,tparam(5),'lambda5'
            write(file,102) 16,tparam(6),'lambda6'
            write(file,102) 17,tparam(7),'lambda7'
c Warning/WARNING/ changed here too
            write(file,102) 4,tparam(8),'m12_2'
         else if (thdmckey.eq.2) then
            write(file,102) 3,tanb,'tan(beta)'
            write(file,102) 4,thdmcm12,'m12'
            write(file,102) 21,tparam(1),'mh'
            write(file,102) 22,tparam(2),'mH'
            write(file,102) 23,tparam(3),'mA'
            write(file,102) 24,tparam(4),'mC'
            write(file,102) 25,tparam(5),'sin(beta-alpha)'
            write(file,102) 26,tparam(6),'lambda6'
            write(file,102) 27,tparam(7),'lambda7'
         else if (thdmckey.eq.3) then
            write(file,102) 3,tanb,'tan(beta)'
            write(file,102) 31,tparam(1),'mh'
            write(file,102) 32,tparam(2),'mH'
            write(file,102) 33,tparam(3),'sin(beta-alpha)'
            write(file,102) 34,tparam(5),'Z4'
            write(file,102) 35,tparam(6),'Z5'
            write(file,102) 36,tparam(7),'Z7'
         else if (thdmckey.eq.4) then
            write(file,102) 3,tanb,'tan(beta)'
            write(file,102) 31,tparam(1),'mh'
            write(file,102) 32,tparam(2),'mH'
            write(file,102) 33,tparam(3),'cos(beta-alpha)'
            write(file,102) 34,tparam(5),'Z4'
            write(file,102) 35,tparam(6),'Z5'
            write(file,102) 36,tparam(7),'Z7'
         endif
      endif

      if (lqqhkey) then
         write(file,103) 'Block QQH # parameters for qq->H process'
         write(file,1020) 1,ptn1,'parton 1 = '//ptname(ptn1)
         write(file,1020) 2,ptn2,'parton 2 = '//ptname(ptn2)
         write(file,102) 11,qqhyuk,'Yukawa coupling'
         write(file,102) 12,qqhscale
     &        ,'renorm.-scale for input Yuk.-coupl.'
      endif

      if (ltmpblock) then
         write(file,103) 'Block TMP # temporary parameters'
         write(file,1020) 1,hxswg,'1: qqh running as required by LHXSWG'
      endif

      if (ldim5) then
         write(file,103) 'Block DIM5'
            write(file,1020) 0,dim5run
     &         ,'Running of DIM5 operators 0=off/1=pert./2=res.'
            write(file,102) htype,c5spin(pseudo,0)
     &           ,'LO coeff of dim-5 operator'
            write(file,102) 100+htype,c5spin(pseudo,1)
     &           ,'NLO coeff of dim-5 operator'
            write(file,102) 200+htype,c5spin(pseudo,2)
     &           ,'NNLO coeff of dim-5 operator'
         if (pseudo.eq.0) then
            write(file,102) 300+htype,c5spin(pseudo,3)
     &           ,'N3LO coeff of dim-5 operator'
         end if
      endif

      write(file,103) 'Block SMINPUTS'
      write(file,102) 1,SUinvAlfa,'1/alpha_em(MZ) SM MSbar'
      write(file,102) 2,GF,'G_F'
      write(file,102) 3,AlfasMZ,'alpha_s(m_Z)'
      write(file,102) 4,mZ,'m_Z(pole)'
      write(file,102) 5,mbmb,'m_b(m_b)'
      write(file,102) 6,mt,'m_t(pole)'
      if (yukfac(1).ne.(0.d0)) write(file,102) 8,mcmc,'m_c(m_c)'
      if (((model.eq.2).and.(thdmcflag.eqv.(.false.))).or.
     &     (model.eq.1).or.(model.eq.3)) then
         write(file,103) 'Block MINPAR'
         write(file,102) 3,tanb,'tan(beta)'
      end if
!Specific SUSY parameters
      if (model.eq.1) then
         write(file,103) 'Block EXTPAR'
         write(file,102) 3,mgl,'Gluino mass'
         write(file,102) 11,Real(cAt),'A_t'
         write(file,102) 12,Real(cAb),'A_b'
         write(file,102) 23,muSUSY,'mu'
         if (fhflag) write(file,102) 26,mA0,'mA'
         if (.not.onshellflag) then
            write(file,102) 43,M3SQ,'M3SQ'
            write(file,102) 46,M3SU,'M3SU'
            write(file,102) 49,M3SD,'M3SD'
         end if
         if (fhflag) then
            write(file,103) 'Block FEYNHIGGSFLAGS'
            write(file,1020) 1,mssmpart,'mssmpart'
            write(file,1020) 2,fieldren,'fieldren'
            write(file,1020) 3,tanbren,'tanbren'
            write(file,1020) 4,higgsmix,'higgsmix'
            write(file,1020) 5,p2approx,'p2approx'
            write(file,1020) 6,looplevel,'looplevel'
            write(file,1020) 7,runningMT,'runningMT'
            write(file,1020) 8,botResum,'botResum'
            write(file,1020) 9,tlCplxApprox,'tlCplxApprox'
            write(file,1020) 10,loglevel,'loglevel'
            write(file,103) 'Block FEYNHIGGS'
            write(file,102) 1,Real(M_1),'M1'
            write(file,102) 2,Real(M_2),'M2'
            write(file,102) 13,Real(cAtau),'A_tau'
            write(file,102) 14,Real(cAc),'A_c'
            write(file,102) 15,Real(cAs),'A_s'
            write(file,102) 16,Real(cAmu),'A_mu'
            write(file,102) 17,Real(cAu),'A_u'
            write(file,102) 18,Real(cAd),'A_d'
            write(file,102) 19,Real(cAe),'A_e'
            write(file,102) 31,M1SL,'M1SL'
            write(file,102) 32,M2SL,'M2SL'
            write(file,102) 33,M3SL,'M3SL'
            write(file,102) 34,M1SE,'M1SE'
            write(file,102) 35,M2SE,'M2SE'
            write(file,102) 36,M3SE,'M3SE'
            write(file,102) 41,M1SQ,'M1SQ'
            write(file,102) 42,M2SQ,'M2SQ'
            write(file,102) 44,M1SU,'M1SU'
            write(file,102) 45,M2SU,'M2SU'
            write(file,102) 47,M1SD,'M1SD'
            write(file,102) 48,M2SD,'M2SD'
         end if
      end if
!NMSSM specific
      if (model.eq.3) then
         write(file,103) 'Block EXTPAR'
         write(file,102) 3,mgl,'Gluino mass'
         write(file,102) 11,Real(cAt),'A_t'
         write(file,102) 12,Real(cAb),'A_b'
         if (.not.onshellflag) then
            write(file,102) 43,M3SQ,'M3SQ'
            write(file,102) 46,M3SU,'M3SU'
            write(file,102) 49,M3SD,'M3SD'
         end if
         write(file,102) 61,lam,'lambda'
!write(file,102) 62,kap,'kappa'
!write(file,102) 63,Alam,'Alambda'
!write(file,102) 64,Akap,'Akappa'
         write(file,102) 65,muSUSY,'eff. mu'
      end if
      if (((model.eq.2).and.(thdmcflag.eqv.(.false.))).or.
     &     ((model.eq.1).and.(fhflag.eqv.(.false.)))) then
         write(file,103) 'Block ALPHA # Higgs mixing parameter'
         write(file,107) alpha,'alpha'
      end if
      if (model.eq.0) then
         write(file,103) 'Block MASS'
         write(file,102) 25,Mh,'Higgs mass'
      end if
      if (((model.eq.2).and.(thdmcflag.eqv.(.false.))).or.
     &   ((model.eq.1).and.(fhflag.eqv.(.false.))).or.(model.eq.3)) then
         write(file,103) 'Block MASS'
         write(file,1024) 25+10*(Sind-1)+pseudo,Mh,
     &        Sind2,'. '//psstr//' Higgs mass in GeV'
!for on-shell flag case
         if (((model.eq.1).or.(model.eq.3)).and.(onshellflag)) then
            write(file,102) 1000005,dsqrt(msb12),'Sbottom1 mass in GeV'
            write(file,102) 2000005,dsqrt(msb22),'Sbottom2 mass in GeV'
            write(file,102) 1000006,dsqrt(mst12),'Stop1 mass in GeV'
            write(file,102) 2000006,dsqrt(mst22),'Stop2 mass in GeV'
         end if
         if (((model.eq.1).or.(model.eq.3)).and.(onshellflag)) then
            write(file,103) 'Block STOPMIX # stop mixing matrix'
            write(file,108) 1,1,cthetat,'V_11'
            write(file,108) 1,2,sthetat,'V_12'
            write(file,103) 'Block SBOTMIX # sbottom mixing matrix'
            write(file,108) 1,1,cthetab,'V_11'
            write(file,108) 1,2,sthetab,'V_12'
         endif

      end if
      if (model.eq.3) then
         write(file,103) 'Block NMHMIX # CP-even Higgs mixing matrix'
         do count1=1,3
            do count2=1,3
               write(file,108) count1,count2,Hmx(count1,count2),
     &              'Mixing matrix'
            end do
         end do
         write(file,103) 'Block NMAMIXR # Prerot.'
     &        //' CP-odd Higgs mixing matrix'
         do count1=2,3
            do count2=2,3
               !changed to count1-1, 18/7/2015
               write(file,108) count1-1,count2-1,Amx(count1,count2),
     &              'Mixing matrix'
            end do
         end do
      end if
!only write relevant DISTRIBUTION input
      if(.not.extflag) then
         write(file,103) 'Block DISTRIB'
         if (dist.eq.0) then
            resstring = 'total XS'
         else if (dist.eq.1) then
            resstring = 'diff. XS dsigma/dpt'
         else if (dist.eq.2) then
            resstring = 'diff. XS dsigma/dy'
            if (pseudorap) resstring = 'diff. XS dsigma/deta'
         else if (dist.eq.3) then
            resstring = 'diff. XS dsigma/dy/dpt'
            if (pseudorap) resstring = 'diff. XS dsigma/deta/dpt'
         end if
         write(file,1020) 1,dist,resstring
         if (ptc.eq.0) then
            if ((dist.eq.1).or.(dist.eq.3))
     &           write(file,102) 22,pt,
     &           'pt in GeV (for diff. XS evaluation)'
         else if (ptc.ge.1) then
            write(file,1020) 2,ptc,'pt cut'
            if (ptc.eq.1) then
               write(file,102) 21,minptc,'Minimal pt value in GeV'
            else if (ptc.eq.2) then
               write(file,102) 22,maxptc,'Maximal pt value in GeV'
            else if (ptc.eq.3) then
               write(file,102) 21,minptc,'Minimal pt value in GeV'
               write(file,102) 22,maxptc,'Maximal pt value in GeV'
            end if
         end if
         if (pseudorap) then
            if (rapc.eq.0) then
               if ((dist.eq.2).or.(dist.eq.3))
     &              write(file,102) 32,y,'eta (for diff. XS evaluation)'
            else if (rapc.ge.1) then
               write(file,1020) 3,rapc,'eta cut'
               if (rapc.eq.2) then
                  write(file,102) 31,minrap,'Minimal eta value'
               else if (rapc.eq.1) then
                  write(file,102) 32,maxrap,'Maximal eta value'
               else if (rapc.eq.3) then
                  write(file,102) 31,minrap,'Minimal eta value'
                  write(file,102) 32,maxrap,'Maximal eta value'
               end if
            end if
         else
            if (rapc.eq.0) then
               if ((dist.eq.2).or.(dist.eq.3))
     &              write(file,102) 32,y,'y (for diff. XS evaluation)'
            else if (rapc.ge.1) then
               write(file,1020) 3,rapc,'y cut'
               if (rapc.eq.2) then
                  write(file,102) 31,minrap,'Minimal y value'
               else if (rapc.eq.1) then
                  write(file,102) 32,maxrap,'Maximal y value'
               else if (rapc.eq.3) then
                  write(file,102) 31,minrap,'Minimal y value'
                  write(file,102) 32,maxrap,'Maximal y value'
               end if
            end if
         end if

         if((dist.ge.2).or.(rapc.ne.0)) then
            if (pseudorap) then
               write(file,1020) 4,1,'Pseudorapidity'
            else
               write(file,1020) 4,0,'Rapidity'
            endif
         endif
      endif

      write(file,103) 'Block SCALES'
      if (runscale) then
         write(file,102) 1,muRfacggh,
     &        'Renormalization scale muR/Sqrt(mh^2+pt^2)'
         write(file,102) 2,muFout,
     &        'Factorization scale muF/Sqrt(mh^2+pt^2)'
         write(file,1020) 3,1,'Use p_t dependent scale choice'
         if(extflag) write(file,102) 4,Qres,
     &        'Resummation scale Qres in GeV'
      else
         write(file,102) 1,muRggh/mh,'Renormalization scale muR/mh'
         write(file,102) 2,muFout,'Factorization scale muF/mh'
         if(.not.extflag) then
            write(file,109) 101,facmurint1,facmurint2
     &        ,'min and max for muR scale uncertainty around SCALES(1)'
            write(file,110) 102,facmurgghmin,facmurgghmax,facmurgghstep
     &        ,'min/max/steps for table output of muR variation'
            if(norderbbh.ge.0) then
               write(file,102) 11,murfacbbh
     &              ,'Renormalization scale muR/mh for '//qqhstring
               write(file,102) 12,muffacbbh,'Factorization scale'
     &              //' muF/mh for '//qqhstring
            endif
         else
            write(file,1020) 3,0,'Use p_t dependent scale choice'
            write(file,102) 4,Qres,'Resummation scale Qres in GeV'
         endif
      end if
      write(file,103) 'Block RENORMBOT'
      write(file,1020) 1,mbrunyuk,'mb used for bottom Yukawa'
      if ((model.eq.1).or.(model.eq.3)) then
         write(file,1020) 2,tanbresum,'Resummation of sbottom effects'
         if (fhflag) write(file,1020) 3,delmbfh,'Delta_b from FeynHiggs'
      end if
      if (mbossave.gt.0.D0)
     &     write(file,102) 4,mbossave,'Fixed m_b(pole) in GeV'
      if ((model.eq.1).or.(model.eq.3)) then
         write(file,103) 'Block RENORMSBOT'
         write(file,1020) 1,mbsbch,'Renormalization of m_b'
         write(file,1020) 2,Abch,'Renormalization of A_b'
         write(file,1020) 3,thetabch,'Renormalization of theta_b'
      end if
      write(file,103) 'Block PDFSPEC'
      do i=1,1+max(norderggh,norderbbh)
         write(file,1022) i,pdfnames(i),ordstr(i)//' PDF set'
         write(file,1020) i+10,pdfset(i),ordstr(i)//' PDF set number'
      enddo

      if(.not.extflag) then
         write(file,103) 'Block VEGAS'
         write(file,103) '# parameters for NLO SusHi'
         write(file,1021) 1,ncall1sushi,'Number of points'
         write(file,1021) 2,itmx1sushi,'Number of iterations'
         if (lveg1sushi) then
            write(file,1021) 11,ncall2sushi
     &           ,'Number of points in second run'
            write(file,1021) 12,itmx2sushi
     &           ,'Number of iterations in second run'
         endif
         write(file,1021) 3,nprnvsushi,'Output format of VEGAS'
     &        //' integration'
         write(file,103) '# parameters for ggh@nnlo:'
         write(file,1021) 4,ncall1ggh,'Number of points'
         write(file,1021) 5,itmx1ggh,'Number of iterations'
         if (lveg1ggh) then
            write(file,1021) 14,ncall2ggh,'Number of points in second'
     &           //' run'
            write(file,1021) 15,itmx2ggh
     &           ,'Number of iterations in second run'
         endif
         write(file,1021) 6,nprnvggh,'Output format of VEGAS'
     &        //' integration'
         write(file,103) '# parameters for bbh@nnlo:'
         write(file,1021) 7,ncall1bbh,'Number of points'
         write(file,1021) 8,itmx1bbh,'Number of iterations'
         if (lveg1bbh) then
            write(file,1021) 17,ncall2bbh,'Number of points in'
     &           //' second run'
            write(file,1021) 18,itmx2bbh
     &           ,'Number of iterations in second run'
         endif
         write(file,1021) 9,nprnvbbh,'Output format of VEGAS'
     &        //' integration'
      endif

      write(file,103) 'Block FACTORS'
      write(file,102) 1,yukfac(1),
     &     'Factor multiplied with Yukawa of c quark'
      write(file,102) 2,yukfac(2),'t quark'
      write(file,102) 3,yukfac(3),'b quark'
      if ((model.eq.1).or.(model.eq.3)) then
         write(file,102) 4,yukfac(6),'top squark'
         write(file,102) 5,yukfac(7),'bottom squark'
      end if

 102  format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 109  format(1x,I9,3x,1P,2E10.2,3x,'#',1x,A)
 110  format(1x,I9,3x,1P,2E10.2,1x,I5,3x,'#',1x,A)
 1020 format(1x,I9,3x,1P,i4.1,0P,3x,'#',1x,A)
 1021 format(1x,I9,3x,1P,i6.1,0P,3x,'#',1x,A)
 1022 format(1x,I9,3x,1P,A45,0P,3x,'#',1x,A)
 1024 format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,I1,A)
 1025 format(1x,I9,3x,i2,1x,i2,1x,i2,3x,'#',1x,A)
 103  format(a)
 107  format(9x,1P,E16.8,0P,3x,'#',1x,A)
 1070 format(9x,1P,i2.1,0P,3x,'#',1x,A)
 108  format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
      end

C-}}}
C-{{{ subroutine writegghmtblock:

      subroutine writegghmtblock(file)

      implicit none
      integer ndum,file

      include '../commons/common-inputoutput.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-int.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-quark.f'
      include '../commons/common-flags.f'
      include '../commons/common-expand.f'
      include '../commons/common-lumis.f'
      include '../commons/common-ptnames.f'

      write(file,103) 'Block GGHMT'
      write(file,1020) -1,nfaclo,'factor out exact LO result'//
     &     ' at LO(=0)/NLO(=1)/etc.'
      write(file,1020) 0,nmtlim(0,0)
     &     ,'expand through 1/mt^n at LO [-1=exact]'

      if (norderggh.ge.1) then
         write(file,1020) 1,nmtlim(1,0),'expand through 1/mt^n at NLO'
         write(file,1020) 11,nmtlim(1,1)
     &        ,'expand gg through 1/mt^n at NLO'
         write(file,1020) 12,nmtlim(1,2)
     &        ,'expand qg through 1/mt^n at NLO'
         write(file,1020) 13,nmtlim(1,3)
     &        ,'expand qqbar through 1/mt^n at NLO'
      endif
      if (norderggh.ge.2) then
c$$$  write(file,102), 9,delmatch,'parameter for matching to x->0'
         write(file,1020) 2,nmtlim(2,0),'expand through 1/mt^n at NNLO'
         write(file,1020) 21,nmtlim(2,1)
     &        ,'expand gg through 1/mt^n at NNLO'
         write(file,1020) 22,nmtlim(2,2)
     &        ,'expand qg through 1/mt^n at NNLO'
         write(file,1020) 23,nmtlim(2,3)
     &        ,'expand qqbar through 1/mt^n at NNLO'
         write(file,1020) 24,nmtlim(2,4)
     &        ,'expand qq through 1/mt^n at NNLO'
         write(file,1020) 25,nmtlim(2,5)
     &        ,'expand qqprime through 1/mt^n at NNLO'
      endif
      if (norderggh.ge.3) then
c$$$  write(file,102), 9,delmatch,'parameter for matching to x->0'
         write(file,1020) 3,nmtlim(3,0),'expand through 1/mt^n at N3LO'
      endif

      ndum=0
      if ((norderggh.ge.1).and.lball(1)) ndum=1
      write(file,1020) 10,ndum
     &     ,'[0/1]: do not/match'/
     &     /' to high energy limit at NLO'
      ndum=0
      if ((norderggh.ge.2).and.lball(2)) ndum=1
      write(file,1020) 20,ndum,'[0/1]: do not/match'//
     &     ' to high energy limit at NNLO'
      ndum=0
c.. this is not available yet:
      if ((norderggh.ge.3).and.lball(3)) ndum=1
      write(file,1020) 30,ndum,'[0/1]: do not/match'//
     &     ' to high energy limit at N3LO'

      if (norderggh.ge.1) then
      write(file,103) 'Block GGHSOFT # parameters for soft expansion'
         write(file,1025) 1,nexpand(1),nsoft(1,0),ncat1
     &        ,'NLO: [0/1=n/y] [order] sig(x)/x^[n]'
c$$$         write(file,1025) 11,nexpandmt(1),nsoft1mt,ncat1mt
c$$$     &        ,'NLO  1/mt:   [0/1=n/y] [order] sig(x)/x^[n]'
      endif
      if (norderggh.ge.2) then
         write(file,1025) 2,nexpand(2),nsoft(2,0),ncat2
     &        ,'NNLO: [0/1=n/y] [order] sig(x)/x^[n]'
c$$$         write(file,1025) 12,nexpandmt(2),nsoft2mt,ncat2mt
c$$$     &        ,'NNLO 1/mt:   [0/1=n/y] [order] sig(x)/x^[n]'
      endif
c$$$      if (norderggh.ge.3) then
      if (norderggh.ge.3) then
         write(file,1025) 3,nexpand(3),nsoft(3,0),ncat3
     &        ,'N3LO: [0/1=n/y] [order] sig(x)/x^[n]'
c$$$         write(file,1025) 13,nexpandmt(3),nsoft3mt,ncat3mt
c$$$     &        ,'N3LO 1/mt:   [0/1=n/y] [order] sig(x)/x^[n]'
      endif

 103  format(a)
 1020 format(1x,I9,3x,1P,i4.1,0P,3x,'#',1x,A)
 1025 format(1x,I9,3x,i2,1x,i2,1x,i2,3x,'#',1x,A)

      end

C-}}}

C-{{{ subroutine screenoutputsushi

      subroutine screenoutputsushi()

      implicit none
      character hepunit*10

      include '../commons/common-inputoutput.f'
      !internal definitions

!defined at other places:
!AlfasMZ,muRfacggh,muFfacggh,cme,dist,prefix,prefix2,pty,iset
!norder,mw,mc,ew,mz,mbmb,renscheme,mt,pseudo,tanbresum,pseudorap
      include '../commons/common-keys.f'
      include '../commons/common-slha.f'
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-quark.f'
      include '../commons/common-int.f'

      if ((norderggh.eq.0).and.(ew.ge.1)) then
            write(*,105)
            write(*,*) 
     &"Note: At LO electroweak contributions for ggh are not given."
            write(*,105)
      endif

      if (((dist.ge.1).or.(ptcut.or.rapcut)).and.(ew.eq.2)) then
            write(*,105)
            write(*,*) 
     &"Note: The SM electroweak factor is not applicable"
            write(*,*) "for cuts and distributions."
            write(*,105)
      endif

      if (((dist.eq.1).or.(dist.eq.3).or.ptcut) .and.
     &     (ptref.lt.20.d0)) then
            write(*,105)
            write(*,*) 
     &"Note: Fixed order p_T distribution (no resummation)!"
            write(*,*) "p_T < 20 GeV not reasonable."
            write(*,105)
      endif

      if (norderggh.ge.2) then
         if (dist.ne.0) then
            write(*,105)
            write(*,*) 
     &"Note: Differential cross sections for ggh are not available ",
     &"at NNLO."
            write(*,*) 
     &"Thus only NLO results will be given."
            write(*,105)
         elseif (ptcut.or.rapcut) then
            write(*,105)
            write(*,*) 
     &"Note: Cuts can only be applied at NLO for ggh."
            write(*,*) 
     &"Thus NNLO weighted results will not be given."
            write(*,105)
         elseif (model.eq.1) then
            if(norderggh.eq.2) then
               write(*,105)
               write(*,*) 
     &"Note: NNLO SUSY contribution not taken into account."
               write(*,*) 
     &"Use order ggh = 3 for approximate NNLO stop cont. for h."
               write(*,105)
            else
               write(*,105)
               write(*,*) 
     &"Note: Using approximate NNLO top+stop contributions."
               write(*,*) 
            endif
         endif
      endif
      if (norderbbh.ge.0) then
         if (dist.gt.1) then
            write(*,105)
            write(*,*) 
     &           "Note: Differential cross sections for "//qqhstring/
     &           /" are only available for pT."
            write(*,*) 'In addition, cuts for the bbh total cross'//
     &           ' sections can be set:'
            write(*,*) 'use DISTRIB(1)=0, DISTRIB(2)>0, DISTRIB(3)>0'
            write(*,*) 'or DISTRIB(1)=1, DISTRIB(3)>0'
            write(*,105)
         endif
         if ((dist.eq.0).and.(norderbbh.eq.2)) then
         if (ptcut.or.rapcut) then
            write(*,105)
            write(*,*) 
     &"Note: Cuts can only be applied at NLO for "//qqhstring//"."
            write(*,*) 
     &"Thus NNLO results will not be given."
            write(*,105)
         endif
         endif
      end if

      if ((dist.eq.0).or.(dist.eq.2)) then
         hepunit = "   pb"
      else
         hepunit = "   pb/GeV"
      end if

      if (norderbbh.ge.0) then
         if ((dist.eq.0).or.(dist.eq.1)) then
            write(*,105)
            write(*,107) " The "//qqhstring//" results are given by:"
            write(*,106) "LO: ",SUbbh(1),SUbbherr(1),hepunit
            if (norderbbh.gt.0) then
               write(*,106) "NLO: ",SUbbh(2),SUbbherr(2),hepunit
            end if
            if (norderbbh.gt.1) then
            if ((.not.((ptc.gt.0).or.(rapc.gt.0))).and.(dist.eq.0)) then
               write(*,106) "NNLO: ",SUbbh(3),SUbbherr(3),hepunit
            end if
            end if
         end if
      end if

      if (norderggh.ge.0) then
         write(*,105)
         write(*,107) " The ggh results are given by:"
         write(*,106) "LO: ",xsec(1),xsecerr(1),hepunit
         if (norderggh.gt.0) then
            write(*,106) "gg: ",ggchan,ggchanerr,hepunit
            write(*,106) "qg: ",qgchan, qgchanerr,hepunit
            write(*,106) "qq: ",qqchan, qqchanerr,hepunit
            write(*,106) "NLO: ",xsec(2),xsecerr(2),hepunit
!            if ((dist.eq.0).and.((ew.ge.1).or.(norderggh.ge.2))) then
            if (((ew.ge.1).or.(norderggh.ge.2))
     &.and.(((dist.eq.0).and.(rapc.eq.0).and.(ptc.eq.0)).or.(ew.eq.1)))
     & then
         write(*,106) "xsec: ",xeff,xefferr,hepunit
!               write(*,107)
!     &" Details with respect to the weighted result in output-file."
            end if
         end if
         if (muranalytic) then
         write(*,109)
     &   "  muR unc.: ",dabs(gghXSvar(1)),gghXSvar(2),hepunit
         write(*,107)
     &   " The scale uncertainty is obtained by a variation"
         write(*,1030) " of muR within [",dexp(xeffvar(1))/murggh,
     &  ',',dexp(xeffvar(2))/murggh,']*',murggh,' GeV.'
         end if
      end if

      write(*,105)

 105  format('#--------------------------------------------------#')
 106  format(a11,'    ',f12.5,'   +-',f12.5,a10)
 109  format(a11,'  - ',f12.5,'   + ',f12.5,a10)
 107  format(a)
 108  format(a11,'    ',f16.5, '   pb')
 1030 format(A,E12.5,A,E12.5,A,E12.5,A)

      end
C-}}}
C-{{{ function mapreal8:

      real*8 function mapreal8(val,lval,message)

      implicit real*8(a-z)
      character(*) message
      logical lval

      if (lval) then
         mapreal8 = val
      else
         mapreal8 = -12345678901234567890.d0
         call printdie(message)
      endif

      return
      end

C-}}}
C-{{{ function mapreal8def:

      real*8 function mapreal8def(val,lval,default,message)

      real*8 default,val
      character(*) message
      logical lval

      if (lval) then
         mapreal8def = val
      else
         mapreal8def = default
         call printinfosushi(6,message)
      endif

      return
      end

C-}}}
C-{{{ function mapint:

      integer function mapint(val,lval,message)

      implicit real*8(a-z)
      integer val
      character(*) message
      logical lval

      if (lval) then
         mapint = val
      else
         call printdie(message)
      endif

      return
      end

C-}}}
C-{{{ function mapintdef:

      integer function mapintdef(val,lval,default,message)

      integer default,val
      character(*) message
      logical lval

      if (lval) then
         mapintdef = val
      else
         mapintdef = default
         call printinfosushi(6,message)
      endif

      return
      end

C-}}}
C-{{{ function mapintdef2:

      integer function mapintdef2(verb,val,lval,min,max,default,defmess
     &     ,errmess)

      integer default,val,verb,min,max
      character(*) defmess,errmess
      logical lval

      if (lval) then
         if ((val.ge.min).and.(val.le.max)) then
            mapintdef2 = val
         else
            call printdie(errmess)
         endif
      else
         mapintdef2 = default
         call printinfosushi2(verb,6,defmess)
      endif

      return
      end

C-}}}
C-{{{ subroutine filedump:

      subroutine filedump(iunitin,iunitout)
c..
c..   Write file unit=iunitin into file unit=iunitout.
c..   Both files must be open before calling filedump,
c..   and both files will be at EOF after calling filedump.
c..
      integer iunitin,iunitout
      character*200 achar

      rewind(iunitin)
 201  read(iunitin,1112,END=202) achar
      length = lnblnk(achar)
      write(iunitout,1112) achar(1:length)
      goto 201
 202  continue

 1112 format(a)

      end

C-}}}
C-{{{ subroutine printdie:

      subroutine printdie(strng)

      character*(*) strng

      write(6,*) char(27)//'[31mSusHi (fatal):'//char(27)//'[0m ',strng
     &     ,char(27),'[31mStopped.',char(27),'[0m'
      stop

      end

C-}}}
C-{{{ subroutine printcopysushi:

      subroutine printcopysushi(fullflag,iu,verbose)

      integer fullflag, iu, verbose

      include '../commons/common-inputoutput.f'

      if (verbose.ne.0) then
      write(iu,1) "#--------------------------------------------------#"
      write(iu,1) "# SusHi: (Supersymmetric) Higgs production through #"
      write(iu,1) "#  __      __        gluon fusion and bottom-quark #"
      write(iu,1) "# [_  | | [_  |_| |  annihilation                  #"
      write(iu,1) "# __] |_| __] | | |                                #"
      write(iu,1) '#                    Version '//SUversion//SUdate//
     & '  #'
      if (fullflag.eq.1) then
      write(iu,1) "#        R. Harlander, S. Liebler, and H. Mantler  #"
      write(iu,1) "#               (harlander@physik.rwth-aachen.de)  #"
      write(iu,1) "#                        (stefan.liebler@desy.de)  #"
      write(iu,1) "#                       (hendrik.mantler@kit.edu)  #"
      write(iu,1) "#            BU Wuppertal, RWTH Aachen, DESY, KIT  #"
      write(iu,1) "#--------------------------------------------------#"
      write(iu,1) "# SusHi is based on a number of calculations       #"
      write(iu,1) "# due to various groups. Please acknowledge these  #"
      write(iu,1) "# efforts by citing the list of references which   #"
      write(iu,1) "# are included in the output file of every run.    #"
      end if
      write(iu,1) "#--------------------------------------------------#"
      endif

 1    format(a)
      
      end

C-}}}
C-{{{ subroutine printwarnsushi:

      subroutine printwarnsushi(verb,unit,strng)

      integer unit,verb
      character*(*) strng

      if (verb.ne.0) write(unit,*) char(27)//'[33mSusHi (warning): '
     &     //char(27)//'[0m',strng

      end

C-}}}
C-{{{ subroutine printinfosushi:

      subroutine printinfosushi(unit,strng)

      integer unit,verb
      character*(*) strng
      include '../commons/common-keys.f'
      
      if ((verbose.ne.0).and.(strng.ne.'')) write(unit,*)
     &     'SusHi (info): ',strng

      end

C-}}}
C-{{{ subroutine printinfosushi2:

      subroutine printinfosushi2(verb,unit,strng)

      integer unit,verb
      character*(*) strng
      include '../commons/common-keys.f'
      
      if (verb.ne.0) write(unit,*) 'SusHi (info): ',strng

      end

C-}}}
C-{{{ subroutine printerrorsushi:

      subroutine printerrorsushi(verb,unit,strng)

      integer unit,verb
      character*(*) strng

      if (verb.eq.1) write(unit,*) 'SusHi (error): ',strng
      stop

      end

C-}}}
C-{{{ function bool2int:

      integer function bool2int(bool)

      implicit none
      logical bool

      if (bool) then
         bool2int=1
      else
         bool2int=0
      endif

      return
      end
      
C-}}}
C-{{{ function int2bool:

      logical function int2bool(int)

      implicit none
      integer int

      if (int.eq.0) then
         int2bool=.false.
      else
         int2bool=.true.
      endif

      return
      end
      
C-}}}
