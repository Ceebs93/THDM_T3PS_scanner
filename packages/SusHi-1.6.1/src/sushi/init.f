! This file is part of SusHi.
!
! It includes initialization routines.
!

C-{{{ initsushi

      subroutine initsushi(htl_in)
      implicit none
      double precision ms42
      integer htl_in,i

      double precision thdmcparams(15)
      integer slha
      integer thdmc_set_param
      character newline*17

      include '../commons/common-expand.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-citations.f'
      include '../commons/common-inputoutput.f'
      include '../commons/common-quark.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'
      include '../commons/common-ren.f'
      include '../commons/common-flags.f'
      include '../commons/common-lumis.f'

      SUversion='1.6.1       '
      SUdate=   ' Jul 2016 '

c..   run alphas at 3-loop or 4-loop order for N3LO:
      runasn3lo = 4
      
      call citeinit()
      citations(1) = 1
      citations(2) = 1

      !read-in of parameters using readinsushi from inputoutput.f
      call readinsushi()

      if ((ew.eq.1).and.(norderggh.ge.1)) then
         citations(11) = 1
         citations(12) = 1
      endif
      if ((ew.eq.2).and.(norderggh.ge.1).and.
     &     ((dist.eq.0).and.(rapc.eq.0).and.(ptc.eq.0))) then
         citations(10) = 1
      endif
      if(norderggh.ge.1) then
         citations(26) = 1
         citations(20) = 2
         citations(21) = 2
         citations(22) = 2
      endif
      if (model.eq.3) citations(3) = 1
      if (lqqhkey) citations(4) = 1

      !previously u and d excluded through (abs(ptn1/2).lt.3)
      if ((ptn1*ptn2.ge.0)) then      ! only qq-bar allowed
         write(6,*) 'SusHi: parton 1 = ',ptn1
         write(6,*) 'SusHi: parton 2 = ',ptn2
         call printdie('Settings in Block QQH not allowed. ')
      endif
      if (lqqhkey) norderggh = -1

      if (nmtlim(2,0).gt.0) then
         citations(39) = 1
         citations(40) = 2
         citations(41) = 2
         citations(42) = 2
         citations(43) = 2
      endif

      apimzlist=0.d0

c$$$         if (lball(1).or.lball(2).or.lball(3)) then
c$$$            write(6,*) 'lball(1) = ',lball(1)
c$$$            write(6,*) 'lball(2) = ',lball(2)
c$$$            write(6,*) 'lball(3) = ',lball(3)
c$$$            call printdie('Matching to s->infty limit not possible'
c$$$     &           //' for Mh > 2*Mt')
c$$$         endif

!     run FeynHiggs, if required:
      if (fhflag) then
#ifndef FEYNHIGGS
         write(*,*) "SusHi: usage error"
         write(*,*) "SusHi is currently not linked to FeynHiggs,"
         write(*,*) "but input-file contains FEYNHIGGS block:"
         write(*,*) "please modify input file, or re-compile with"
         write(*,*) "make predef=FH"
         stop
#else
         citations(16) = 1
         citations(17) = 1
         citations(18) = 1
         citations(19) = 1
         call runfeynhiggs(alfasmz,mt,mbmb,mw,mz,mh,pseudo
     &,msb12,msb22,mst12,mst22)
#endif
      end if

!     run 2HDMC, if required:
      if (thdmcflag) then
#ifndef THDMC
         write(*,*) "Input file contains block 2DHMC,"
         write(*,*) "but SusHi was compiled without 2HDMC support."
         write(*,*) "Either re-compile with 'make predef=2HDMC',"
         write(*,*) "or remove 2HDMC Block from input file."
         stop
#else
         citations(25) = 1
         call runthdmc()
#endif
      endif

      if (Mh.Le.20.d0) then
         call printwarnsushi(1,6,'You called SusHi'//
     &           ' with a small Higgs mass.')
         call printwarnsushi(1,6,'Results might not be reasonable.')
         call printwarnsushi(1,6,'Electroweak contributions'//
     &           ' not implemented.')
      endif
      if (Mh.Ge.1500.d0) then
         call printwarnsushi(1,6,'You called SusHi'//
     &           ' with a large Higgs mass.')
         call printwarnsushi(1,6,'Results might not be reasonable.')
         call printwarnsushi(1,6,'Electroweak contributions'//
     &           ' not implemented.')
      end if

      q0 = q0c + q01*dsqrt(mst12) + q02*dsqrt(mst22)

      !determine whether to use the heavy top limit 
      !m_t -> inf in the top-cont.
      if(htl_in.eq.0) then
         htl = .false.
      else
         htl = .true.
      endif

      !call ffini

      !settings in case of a fourth generation
      if((yukfac(4).ne.0.d0).or.(yukfac(5).ne.0.d0)) then

         if(mbp.ne.mtp) then
            write(*,105)
            write(*,*) 'Warning: masses of the 4th generation quarks '
            write(*,*) 'are not equal. Setting m_b` = m_t`. '
            write(*,105)
            mbp = mtp
         endif

         ms42 = MSUSY**2+mtp**2
         mstp12 = ms42 - MSUSY**2
         mstp22 = ms42 + MSUSY**2
         msbp12 = mstp12
         msbp22 = mstp22
         mstp1 = dsqrt(mstp12)
         mstp2 = dsqrt(mstp22)
         msbp1 = dsqrt(msbp12)
         msbp2 = dsqrt(msbp22)

         Mh = dsqrt(
     &          Mh**2
     &        + (3/2.d0 * dlog(ms42/Mbp**2) * Mbp**4
     &        +  3/2.d0 * dlog(ms42/Mtp**2) * Mtp**4 ) 
     &        / (Pi/dsqrt(dsqrt(2.d0)*GF))**2 )

         if(Mh.gt.350.d0) then
            write(*,105)
            write(*,*) "Error: Higgs mass too large.", Mh
            write(*,105)
            stop
         endif
      endif

!      mst12 = Mstop1**2
!      mst22 = Mstop2**2     
      mcmc = Mcharm
      mbp2 = Mbp**2
      mtp2 = Mtp**2
!      msb12 = Msbot1**2
!      msb22 = Msbot2**2     

      !definition of renormalization scale and variables
      muRggh = Mh * muRfacggh
      muFggh = Mh * muFfacggh
      muRbbh = Mh * muRfacbbh
      muFbbh = Mh * muFfacbbh

      if (mbrunyuk.eq.1) then
         muB = mbmb
      else if (mbrunyuk.eq.2) then
         muB = muRggh
      else
         muB = Mh
      end if

      mh2 = Mh**2
      z = mh2/cme**2
      sigmanull = GF/288.d0/dsqrt(2.d0)/Pi
     &     * 0.38937966d9 ! last factor: conversion gev -> pb
      nf = 5.d0 ! note: N3LO-results only available for nf=5!
      betanull = 11/12.d0*ca-nf/6.d0
      ln2 = dlog(2.d0)
      lfh = dlog(muFggh**2/mh2)
      lrh = dlog(muRggh**2/mh2)
      lfr = dlog(muFggh**2/muRggh**2)
      lrt = dlog(muRggh**2/mt**2)

      !some routines which should capture ill-defined SusHi calls
      !needs to be at the end of init to have a well-defined Higgs mass
      do i=1,2
         if ( (nmtlim(i,0).ge.2).and.(nexpand(i).eq.0)) then
         call printdie
     &           ('1/mt expansion only available with soft expansion. ')
         endif
      enddo

      newline = NEW_LINE('a')//'                '
      do i=1,3
         if ( (lball(i).and.(norderggh.ge.i))
     &        .and. (nexpand(i).eq.0) ) then
            call printdie
     &           ('Matching to high energy limit only available with'
     &           //newline//'soft expansion. '//
     &           'Change GGHMT(10/20/30) or Block GGHSOFT.'//newline)
         endif
      enddo

      do i=2,3
         if((norderggh.ge.i).and.(nexpand(i-1).gt.0)) then
c     check orders of soft expansion 
            if((nexpand(i).eq.0).or.(nsoft(i,0).gt.nsoft(i-1,0))) then
               call printdie
     &              ('Order of the soft expansion must not increase'
     &              // ' with ascending perturbativ order. ')
            endif
c     check orders of expansion in 1/mt
            if(nmtlim(i,0).gt.nmtlim(i-1,0)) then
               call printdie
     &              ('Order of the expansion in 1/mt must not increase'
     &              // ' with ascending perturbativ order. ')
            endif
         endif
      enddo

      if (ldim5) then
         if(.not.all(nmtlim.le.1)) then         ! first 1/mt term appears at O(1/mt^2)
!     do not allow for 1/mt terms in case dim5 operators are included
            call printdie
     &           ('1/mt terms cannot be used together with'
     &           // ' dim5 operators. ')            
         endif
!         if(lball(1).or.lball(2).or.lball(3)) then
!!     do not allow for matching in case dim5 operators are included - doable, if Higgs is light enough
!            call printdie
!     &           ('Matching to high energy limit must be'
!     &           // ' deactivated when including dim5 operators. ')            
!         endif
         if(nfaclo.ne.-1) then
!     do not factor out LO top mass dependence in case dim5 operators are included
            call printdie
     &           ('Factoring out LO top mass dependence must be'
     &           // ' deactivated when including dim5 operators. ')            
         endif
      end if

      do i=1,5
         nsoft(1,i)  = nsoft(1,0)
         nsoftmt(1,i)  = nsoftmt(1,0)
         nsoft(2,i)  = nsoft(2,0)
         nsoftmt(2,i)  = nsoftmt(2,0)
         nsoft(3,i)  = nsoft(3,0)
         nsoftmt(3,i)  = nsoftmt(3,0)
      enddo
       
      if (mh.gt.(2*mt)) then
         if (norderggh.gt.2) then
            write(6,*) 'mh = ',mh
            write(6,*) 'mt = ',mt
            call printwarnsushi(1,6,'mh > 2*mt: results beyond'//
     &           ' NNLO are outside validity range.')
            call printwarnsushi(1,6,'The result needs to be'//
     &           ' taken with extreme care!')
            if (ldim5) then
            call printwarnsushi(1,6,'Results induced through'//
     &           ' DIM5 operators only are ok!')
            end if
         endif
         if ((nmtlim(1,0).gt.0).or.(nmtlim(2,0).gt.0)) then
            call printdie('1/mt terms (set in GGHMT) are '//
     &           'divergent for Higgs masses beyond 2*mt.'//newline/
     &           /'Switch them off! ')
         endif
      endif

      if ((mh.le.(100.d0)).or.(mh.gt.(300.d0))) then
      if (lball(1).or.lball(2).or.lball(3)) then
        call printdie(' The matching to the high energy limit, x->0, '//
     &        'is only implemented for'//newline/
     &        /'Higgs masses between 100GeV and 300GeV. ')
      end if
      end if       
      
      if (pseudo.eq.1) then
         if (norderggh.gt.2) then
            call printdie(' Results beyond NNLO are not available'//
     &           ' for pseudoscalar.'//newline/
     &           /'Change SUSHI(5). ')
         endif
         if ((nmtlim(1,0).gt.0).or.(nmtlim(2,0).gt.0)) then
            call printdie(' 1/mt terms are not available'//
     &           ' for pseudoscalar.'//newline//'Delete the Block
     &           GGHMT. ')
         endif
         if ((nexpand(1).ne.0).or.(nexpand(2).ne.0)) then
            call printdie(' soft expansion not available'//
     &           ' for pseudoscalar.'//newline/
     &           /'Change Block GGHSOFT. ')
         endif
      endif
         
!     calculate _murdep file, only set to true in case of on-shell
!     parameters
      !this is done in the renormalize routines in renscheme.f
      muranalytic = .false.

      call SUSHI_BERNINI(18)

 105  format('#--------------------------------------------------#')
      end

C-}}}
C-{{{ initcouplings

      subroutine initcouplings()
      implicit none
      double precision cthetab4,sthetab4,cthetat4,sthetat4,z3
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include '../commons/common-inputoutput.f'
      include '../commons/common-quark.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-ren.f'
      include '../commons/common-readdata.f'
      include '../commons/common-flags.f'
      
      !Naming of couplings:
      !gt, gb, gc = Relative Yukawa couplings to top, bottom and charm quark
      !gt3, gb3: Within the routine and the inputoutput we also use those
      !  for the top- and bottom relative Yukawa. They just differ to
      !  the above Yukawa couplings by the numbers of the block FACTORS.
      !gt1,gt2,gb1,gb2: The couplings of the Higgs-bosons to stops and
      !  sbottoms, where the numbers are understood as gt1 = H-stop1-stop1.
      !gt13,gt23,gb13,gb23: Again they equal gt1,gt2,gb1,gb2 except from
      !  the numbers given in block FACTORS and occur here and in
      !  inputoutput.f
      !4-th generation: gbp,gtp,gbp1,gbp2,gtp1,gtp2 are (s)top and (s)bottoms
      !  couplings for a 4th generation of particles.
      !gte, gte2: Wilson coefficients of stop contributions for the virtual amplitude

      z3 = 1.2020569031595942853997d0
      gt3 = 1.d0
      gb3 = 1.d0

      !settings for a fourth generation
      if((yukfac(4).ne.0.d0).or.(yukfac(5).ne.0.d0)) then

         cthetat4 = 0.d0
         sthetat4 = 1.d0
         cthetab4 = 0.d0
         sthetab4 = 1.d0

         call quarkhiggscoup(beta,Hmx,Amx,gt4,gb4,model,
     &twohdmver,pseudo,Sind)

         call evalcsusy(Mbp,msbp1,msbp2,Mtp,Mstp1,Mstp2,Mgl,beta,
     &        Mw,Mz,alpha,muRggh,muSUSY,cthetat4,sthetat4,cthetab4
     &        ,sthetab4,c1sm14,c1sm24,c1sm34,c1susy1t4,c1susy2t4
     &        ,c1susy1b4,c1susy2b4,gb14,gbh124,gbh214,gb24,gt14,gth124
     &        ,gth214,gt24,pseudo,Sind)

         gbpe  = (c1susy1b4-gb4*c1sm14) * yukfac(9)
         gbpe2 = (c1susy2b4-gb4*c1sm24) * yukfac(9)

         gtpe  = (c1susy1t4-gt4*c1sm14) * yukfac(8)
         gtpe2 = (c1susy2t4-gt4*c1sm24) * yukfac(8)

         gbp = gb4 * yukfac(5)
         gtp = gt4 * yukfac(4)
         gbp1 = gb14 * yukfac(9)
         gbp2 = gb24 * yukfac(9)
         gtp1 = gt14 * yukfac(8)
         gtp2 = gt24 * yukfac(8)

      else

         gbp = 0.d0
         gtp = 0.d0
         gbp1 = 0.d0
         gbp2 = 0.d0
         gtp1 = 0.d0
         gtp2 = 0.d0
         gbpe = 0.d0
         gbpe2 = 0.d0
         gtpe = 0.d0
         gtpe2 = 0.d0

      endif

      call getmass()

      if (mbossave.gt.0.D0) then
         mbos = mbossave
      end if

      comc1eff0rd = (1.d0,0.d0)
      c1eff1rd = 2.75d0
      c1eff2rd = (2777/288.d0 + 19/16.d0*lrt + nf*( -67/96.d0 +lrt
     &     /3.d0))
      c1eff3rd = (897943/9216.d0*z3 + 209/64.d0*lrt**2 
     & + 1733/288.d0*lrt - 2892659/41472.d0 + nf*(-110779/13824.d0*z3
     &  + 23/32.d0*lrt**2 + 55/54.d0*lrt + 40291/20736.d0) + nf**2*( 
     &     -1/18.d0*lrt**2 + 77/1728.d0*lrt - 6865/31104.d0))
      c2eff1rd = 0.d0
      if (pseudo.eq.1) then
         c1eff1rd = 0.d0
         c1eff2rd = 0.d0
         c2eff1rd = (-1/2.d0 + lrt)
      endif

      !different running options of dimension 5 operators
      if (pseudo.eq.1) dim5run = 0 !pseudoscalar coef. does not run
      if (dim5run.eq.1) then !default
      call runc5pert(pseudo,c5sp,lrh)
      else if (dim5run.eq.2) then
      call runc5(pseudo,3,AlfasMZ,mZ,c5sp,mh,murggh) !always at order=3
      c5sp(0,1:3) = 0.d0 !set higher orders to 0
      end if
      
      !SM - avoid wrong gt multiplication in subsequent models
      gt = gt3*yukfac(2)
      if (model.eq.0) then
        comc1eff0rd = gt*comc1eff0rd + c5sp(pseudo,0)
        c1eff1rd = gt*c1eff1rd + c5sp(pseudo,1)
        c1eff2rd = gt*c1eff2rd + c5sp(pseudo,2)
        c1eff3rd = gt*c1eff3rd + c5sp(pseudo,3)
        c2eff1rd = gt*c2eff1rd
      end if

      !NMSSM
      if (model.eq.3) then
         Ab = Real(cAb)
         At = Real(cAt)
         beta = datan(tanb)
         call renormalize()
         q0 = q0c + q01*dsqrt(mst12) + q02*dsqrt(mst22)
         call quarkhiggscoup(beta,Hmx,Amx,gt3,gb3,
     &        model,twohdmver,pseudo,Sind)

	 call squarkhiggscoupNMSSM(mbsb,dsqrt(msb12),dsqrt(msb22),
     &mt,dsqrt(mst12),dsqrt(mst22),beta,Mw,Mz,GF,Hmx,Amx,lam,muSUSY,
     &cthetat,sthetat,cthetab,sthetab,gb13,gbh123,gbh213,gb23,
     &gt13,gth123,gth213,gt23,pseudo,Sind)

         gb = gb3 * gb
         gc = gt3 * yukfac(1)
         gt = gt3 * yukfac(2)
         gb1 = gb13 * yukfac(7)
         gb2 = gb23 * yukfac(7)
         gt1 = gt13 * yukfac(6)
         gt2 = gt23 * yukfac(6)
         if (pseudo.eq.1) then
            comc1eff0rd = mt2*gt*ALOSMA(mh**2,mt2) + c5sp(pseudo,0)
         else
            comc1eff0rd = mt2*gt*ALOSMh(mh**2,mt2) + c5sp(pseudo,0)
         end if
         c1eff1rd = gt*c1eff1rd + c5sp(pseudo,1)
         c1eff2rd = gt*c1eff2rd + c5sp(pseudo,2)
         c1eff3rd = gt*c1eff3rd + c5sp(pseudo,3)
         c2eff1rd = gt*c2eff1rd
      endif
      
      !MSSM
      if (model.eq.1) then
!Ab either taken from input file or FH
         Ab = Real(cAb)
         At = Real(cAt)
         beta = datan(tanb)
         call renormalize()

!         call renormalizesimpleFH(M3SU,M3SQ,M3SD
!     &,msb12,msb22,mst12,mst22,dMLb,beta,mw,mz,mgl,yukfac,muRggh
!     &,mbos,mbsb,mbsb2,delmb,delta_mb,muSUSY,Ab,At,tanbresum
!     &,gb,dthetab,dmbsb,dAbds,dmb,dmsb1,dmsb2,dgb
!     &,mb,mb2,mt,mbyuk
!     &,cthetat,sthetat,c2thetat,s2thetat
!     &,cthetab,sthetab,c2thetab,s2thetab,t2thetab,thetab)

         q0 = q0c + q01*dsqrt(mst12) + q02*dsqrt(mst22)
         call quarkhiggscoup(beta,Hmx,Amx,gt3,gb3,
     &        model,twohdmver,pseudo,Sind)
         call evalcsusy(mbsb,dsqrt(msb12),dsqrt(msb22),mt,
     &        dsqrt(mst12),dsqrt(mst22),Mgl,beta,
     &        Mw,Mz,alpha,muRggh,muSUSY,cthetat,sthetat,cthetab,sthetab,
     &        c1sm13,c1sm23,c1sm33,c1susy1t3,c1susy2t3,c1susy1b3,
     &        c1susy2b3,gb13,gbh123,gbh213,gb23,gt13,gth123,gth213,gt23,
     &        pseudo,Sind)

         gb = gb3 * gb
         gc = gt3 * yukfac(1)
         gt = gt3 * yukfac(2)
         gb1 = gb13 * yukfac(7)
         gb2 = gb23 * yukfac(7)
         gt1 = gt13 * yukfac(6)
         gt2 = gt23 * yukfac(6)
         gte  = (c1susy1t3-gt3*c1sm13) * yukfac(6)
         gte2 = (c1susy2t3-gt3*c1sm23) * yukfac(6)
c.. this is input for ggh@nnlo:         
         !c1eff0rd = gt*c1sm13
         !replaced by exakt amplitude:
         if (pseudo.eq.1) then
            comc1eff0rd = mt2*gt*ALOSMA(mh**2,mt2) + c5sp(pseudo,0)
         else
            comc1eff0rd = mt2*gt*ALOSMh(mh**2,mt2) + c5sp(pseudo,0)
         end if
         !write(*,*) "c1eff0",gt*c1sm13,comc1eff0rd
         !write(*,*) "check",dsqrt(mt2),dsqrt(mst12)
         c1eff1rd = gt*c1sm23   + c5sp(pseudo,1)
         c1eff2rd = gt*c1sm33   + c5sp(pseudo,2)
         c1eff3rd = gt*c1eff3rd + c5sp(pseudo,3)
         c2eff1rd = gt*c2eff1rd
         if (nnlostop) then
            !c1eff0rd = c1eff0rd + gte
            if (pseudo.ne.1) then
            comc1eff0rd = comc1eff0rd
     &                + mt2*(gt1*ALOSUSYh(mh**2,mst12)
     &                     + gt2*ALOSUSYh(mh**2,mst22))
            end if
            !write(*,*) "c1eff0-3",gt*c1sm13+gte,comc1eff0rd
            c1eff1rd = c1eff1rd + gte2
         endif
      endif
      !2HDM
      if (model.eq.2) then
        beta = datan(tanb)
        call quarkhiggscoup(beta,Hmx,Amx,gt3,gb3,
     &     model,twohdmver,pseudo,Sind)
        gb = gb3 * yukfac(3)
        gc = gt3 * yukfac(1)
        gt = gt3 * yukfac(2)
        comc1eff0rd = gt*comc1eff0rd + c5sp(pseudo,0)
        c1eff1rd = gt*c1eff1rd + c5sp(pseudo,1)
        c1eff2rd = gt*c1eff2rd + c5sp(pseudo,2)
        c1eff3rd = gt*c1eff3rd + c5sp(pseudo,3)
        c2eff1rd = gt*c2eff1rd
      endif     
      end

C-}}}
