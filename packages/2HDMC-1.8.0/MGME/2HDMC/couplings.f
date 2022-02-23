c----------------------------------------------------------------------
C  couplings.f 
c----------------------------------------------------------------------
c  This files takes the inputs of the standard model from a Les Houches 
c  file (param_card.dat) and calculates all the couplings that end up
c  in the Feynman rules, i.e., in the HELAS routines.
c   
c  With the convention of the New Web Generation MadEvent in this
c  part no non-trivial computation is made. The SM model secondary
c  parameters have been all calculated by the SM calculator, SMCalc
c  which produced param_card.dat.
c
c  The only exception to the above rule is for alpha_S. In the case
c  where a pdf is used in the initial state, then the values as(MZ)
c  set in the les houches card is superseeded.
c  Scales and pdf's are set in the run_card.dat.
c
c This file contains the following routines:
c 
c- madgraph original call to set the parameters
c- lh_readin is called from here.
c  It talks to the lh_readin through the common block values.
c      subroutine setpara
c
c-routine to read the LH parameters in
c      subroutine lh_readin
c
c-to set
c      subroutine set_it(npara,ivalue,value,name,id,
c      subroutine case_trap(string,length)
c      subroutine no_spaces(buff,nchars)
c---------------------------------------------------------------------- 


      subroutine setpara(param_name,readlha)
c***********************************************************************
c This subroutine sets up the HELAS couplings of the STANDARD MODEL.
c***********************************************************************
      implicit none
c
c local
c
      character*(*) param_name
      logical readlha
      integer i
      real*8 dum
c
c     common file with the couplings
c
      include 'coupl.inc'
      include 'input.inc'
c
c     local
c
      double precision  v
      double precision  ee, ee2, ez, ey, sw, cw, sc2, sin2w, wm
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  alphas, alfa, alfaw, mfrun
      external          alphas, alfa, alfaw, mfrun
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
c
c     alfas and run
c************
c Uncomment the following lines in order to use alphas from the PDF
c      include '../alfas.inc'
c      include '../run.inc'
c***********
c------------------------------------------
c Start calculating the couplings for HELAS
c------------------------------------------
c
      if(readlha) then 
         call lh_readin(param_name)
         G = DSQRT(4d0*PI*ALFAS) ! use setting of the param_card.dat @ NLO
      endif
c     

      GG(1) = -G
      GG(2) = -G     
c
c auxiliary local values
c
      wm = sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      sin2w  = One-(wm/zmass)**2
      cw  = sqrt( One - sin2w )
      ee2 = alpha * Fourpi
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee   ! the wmass is used to calculate v
      lambda = hmass**2 / (Two * v**2)
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )

c
c     CKM matrix: 
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd Vcb=Vub=0
c
c     >>>>>>>>>>>>>>>***** NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter only when interaction_CKM.dat
c     is used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal
c     CKM matrix.

      do i=1,2
         gwfud(i) = gwf(i)*Vud
         gwfus(i) = gwf(i)*Vus
         gwfub(i) = gwf(i)*Vub
         gwfcd(i) = gwf(i)*Vcd
         gwfcs(i) = gwf(i)*Vcs
         gwfcb(i) = gwf(i)*Vcb
         gwftd(i) = gwf(i)*Vtd
         gwfts(i) = gwf(i)*Vts
         gwftb(i) = gwf(i)*Vtb
      enddo

c---------------------------------------------------------
c Set Photon Width to Zero, used by symmetry optimization
c---------------------------------------------------------

      awidth = 0d0
c************************************            
c UserMode couplings
c************************************

      GH1EE(1)=dcmplx(IMSGH1EE-IMPGH1EE,REPGH1EE-RESGH1EE)
      GH1EE(2)=dcmplx(IMSGH1EE+IMPGH1EE,-RESGH1EE-REPGH1EE)
      GH2EE(1)=dcmplx(IMSGH2EE-IMPGH2EE,REPGH2EE-RESGH2EE)
      GH2EE(2)=dcmplx(IMSGH2EE+IMPGH2EE,-RESGH2EE-REPGH2EE)
      GH3EE(1)=dcmplx(IMSGH3EE-IMPGH3EE,REPGH3EE-RESGH3EE)
      GH3EE(2)=dcmplx(IMSGH3EE+IMPGH3EE,-RESGH3EE-REPGH3EE)

      GH1MUMU(1)=dcmplx(IMSGH1MUMU-IMPGH1MUMU,REPGH1MUMU-RESGH1MUMU)
      GH1MUMU(2)=dcmplx(IMSGH1MUMU+IMPGH1MUMU,-RESGH1MUMU-REPGH1MUMU)
      GH2MUMU(1)=dcmplx(IMSGH2MUMU-IMPGH2MUMU,REPGH2MUMU-RESGH2MUMU)
      GH2MUMU(2)=dcmplx(IMSGH2MUMU+IMPGH2MUMU,-RESGH2MUMU-REPGH2MUMU)
      GH3MUMU(1)=dcmplx(IMSGH3MUMU-IMPGH3MUMU,REPGH3MUMU-RESGH3MUMU)
      GH3MUMU(2)=dcmplx(IMSGH3MUMU+IMPGH3MUMU,-RESGH3MUMU-REPGH3MUMU)

      GH1TATA(1)=dcmplx(IMSGH1TATA-IMPGH1TATA,REPGH1TATA-RESGH1TATA)
      GH1TATA(2)=dcmplx(IMSGH1TATA+IMPGH1TATA,-RESGH1TATA-REPGH1TATA)
      GH2TATA(1)=dcmplx(IMSGH2TATA-IMPGH2TATA,REPGH2TATA-RESGH2TATA)
      GH2TATA(2)=dcmplx(IMSGH2TATA+IMPGH2TATA,-RESGH2TATA-REPGH2TATA)
      GH3TATA(1)=dcmplx(IMSGH3TATA-IMPGH3TATA,REPGH3TATA-RESGH3TATA)
      GH3TATA(2)=dcmplx(IMSGH3TATA+IMPGH3TATA,-RESGH3TATA-REPGH3TATA)

      GH1UU(1)=dcmplx(IMSGH1UU-IMPGH1UU,REPGH1UU-RESGH1UU)
      GH1UU(2)=dcmplx(IMSGH1UU+IMPGH1UU,-RESGH1UU-REPGH1UU)
      GH2UU(1)=dcmplx(IMSGH2UU-IMPGH2UU,REPGH2UU-RESGH2UU)
      GH2UU(2)=dcmplx(IMSGH2UU+IMPGH2UU,-RESGH2UU-REPGH2UU)
      GH3UU(1)=dcmplx(IMSGH3UU-IMPGH3UU,REPGH3UU-RESGH3UU)
      GH3UU(2)=dcmplx(IMSGH3UU+IMPGH3UU,-RESGH3UU-REPGH3UU)

      GH1CC(1)=dcmplx(IMSGH1CC-IMPGH1CC,REPGH1CC-RESGH1CC)
      GH1CC(2)=dcmplx(IMSGH1CC+IMPGH1CC,-RESGH1CC-REPGH1CC)
      GH2CC(1)=dcmplx(IMSGH2CC-IMPGH2CC,REPGH2CC-RESGH2CC)
      GH2CC(2)=dcmplx(IMSGH2CC+IMPGH2CC,-RESGH2CC-REPGH2CC)
      GH3CC(1)=dcmplx(IMSGH3CC-IMPGH3CC,REPGH3CC-RESGH3CC)
      GH3CC(2)=dcmplx(IMSGH3CC+IMPGH3CC,-RESGH3CC-REPGH3CC)

      GH1TT(1)=dcmplx(IMSGH1TT-IMPGH1TT,REPGH1TT-RESGH1TT)
      GH1TT(2)=dcmplx(IMSGH1TT+IMPGH1TT,-RESGH1TT-REPGH1TT)
      GH2TT(1)=dcmplx(IMSGH2TT-IMPGH2TT,REPGH2TT-RESGH2TT)
      GH2TT(2)=dcmplx(IMSGH2TT+IMPGH2TT,-RESGH2TT-REPGH2TT)
      GH3TT(1)=dcmplx(IMSGH3TT-IMPGH3TT,REPGH3TT-RESGH3TT)
      GH3TT(2)=dcmplx(IMSGH3TT+IMPGH3TT,-RESGH3TT-REPGH3TT)

      GH1DD(1)=dcmplx(IMSGH1DD-IMPGH1DD,REPGH1DD-RESGH1DD)
      GH1DD(2)=dcmplx(IMSGH1DD+IMPGH1DD,-RESGH1DD-REPGH1DD)
      GH2DD(1)=dcmplx(IMSGH2DD-IMPGH2DD,REPGH2DD-RESGH2DD)
      GH2DD(2)=dcmplx(IMSGH2DD+IMPGH2DD,-RESGH2DD-REPGH2DD)
      GH3DD(1)=dcmplx(IMSGH3DD-IMPGH3DD,REPGH3DD-RESGH3DD)
      GH3DD(2)=dcmplx(IMSGH3DD+IMPGH3DD,-RESGH3DD-REPGH3DD)

      GH1SS(1)=dcmplx(IMSGH1SS-IMPGH1SS,REPGH1SS-RESGH1SS)
      GH1SS(2)=dcmplx(IMSGH1SS+IMPGH1SS,-RESGH1SS-REPGH1SS)
      GH2SS(1)=dcmplx(IMSGH2SS-IMPGH2SS,REPGH2SS-RESGH2SS)
      GH2SS(2)=dcmplx(IMSGH2SS+IMPGH2SS,-RESGH2SS-REPGH2SS)
      GH3SS(1)=dcmplx(IMSGH3SS-IMPGH3SS,REPGH3SS-RESGH3SS)
      GH3SS(2)=dcmplx(IMSGH3SS+IMPGH3SS,-RESGH3SS-REPGH3SS)

      GH1BB(1)=dcmplx(IMSGH1BB-IMPGH1BB,REPGH1BB-RESGH1BB)
      GH1BB(2)=dcmplx(IMSGH1BB+IMPGH1BB,-RESGH1BB-REPGH1BB)
      GH2BB(1)=dcmplx(IMSGH2BB-IMPGH2BB,REPGH2BB-RESGH2BB)
      GH2BB(2)=dcmplx(IMSGH2BB+IMPGH2BB,-RESGH2BB-REPGH2BB)
      GH3BB(1)=dcmplx(IMSGH3BB-IMPGH3BB,REPGH3BB-RESGH3BB)
      GH3BB(2)=dcmplx(IMSGH3BB+IMPGH3BB,-RESGH3BB-REPGH3BB)

      GHCUD(1)=dcmplx(IMSGHCUD+IMPGHCUD,RESGHCUD+REPGHCUD)
      GHCUD(2)=dcmplx(IMSGHCUD-IMPGHCUD,RESGHCUD-REPGHCUD)
      GHCUS(1)=dcmplx(IMSGHCUS+IMPGHCUS,RESGHCUS+REPGHCUS)
      GHCUS(2)=dcmplx(IMSGHCUS-IMPGHCUS,RESGHCUS-REPGHCUS)
      GHCUB(1)=dcmplx(IMSGHCUB+IMPGHCUB,RESGHCUB+REPGHCUB)
      GHCUB(2)=dcmplx(IMSGHCUB-IMPGHCUB,RESGHCUB-REPGHCUB)

      GHCCD(1)=dcmplx(IMSGHCCD+IMPGHCCD,RESGHCCD+REPGHCCD)
      GHCCD(2)=dcmplx(IMSGHCCD-IMPGHCCD,RESGHCCD-REPGHCCD)
      GHCCS(1)=dcmplx(IMSGHCCS+IMPGHCCS,RESGHCCS+REPGHCCS)
      GHCCS(2)=dcmplx(IMSGHCCS-IMPGHCCS,RESGHCCS-REPGHCCS)
      GHCCB(1)=dcmplx(IMSGHCCB+IMPGHCCB,RESGHCCB+REPGHCCB)
      GHCCB(2)=dcmplx(IMSGHCCB-IMPGHCCB,RESGHCCB-REPGHCCB)

      GHCTD(1)=dcmplx(IMSGHCTD+IMPGHCTD,RESGHCTD+REPGHCTD)
      GHCTD(2)=dcmplx(IMSGHCTD-IMPGHCTD,RESGHCTD-REPGHCTD)
      GHCTS(1)=dcmplx(IMSGHCTS+IMPGHCTS,RESGHCTS+REPGHCTS)
      GHCTS(2)=dcmplx(IMSGHCTS-IMPGHCTS,RESGHCTS-REPGHCTS)
      GHCTB(1)=dcmplx(IMSGHCTB+IMPGHCTB,RESGHCTB+REPGHCTB)
      GHCTB(2)=dcmplx(IMSGHCTB-IMPGHCTB,RESGHCTB-REPGHCTB)

      GHCDU(1)=dcmplx(IMSGHCUD-IMPGHCUD,REPGHCUD-RESGHCUD)
      GHCDU(2)=dcmplx(IMSGHCUD+IMPGHCUD,-RESGHCUD-REPGHCUD)
      GHCDC(1)=dcmplx(IMSGHCCD-IMPGHCCD,REPGHCCD-RESGHCCD)
      GHCDC(2)=dcmplx(IMSGHCCD+IMPGHCCD,-RESGHCCD-REPGHCCD)
      GHCDT(1)=dcmplx(IMSGHCTD-IMPGHCTD,REPGHCTD-RESGHCTD)
      GHCDT(2)=dcmplx(IMSGHCTD+IMPGHCTD,-RESGHCTD-REPGHCTD)

      GHCSU(1)=dcmplx(IMSGHCUS-IMPGHCUS,REPGHCUS-RESGHCUS)
      GHCSU(2)=dcmplx(IMSGHCUS+IMPGHCUS,-RESGHCUS-REPGHCUS)
      GHCSC(1)=dcmplx(IMSGHCCS-IMPGHCCS,REPGHCCS-RESGHCCS)
      GHCSC(2)=dcmplx(IMSGHCCS+IMPGHCCS,-RESGHCCS-REPGHCCS)
      GHCST(1)=dcmplx(IMSGHCTS-IMPGHCTS,REPGHCTS-RESGHCTS)
      GHCST(2)=dcmplx(IMSGHCTS+IMPGHCTS,-RESGHCTS-REPGHCTS)

      GHCBU(1)=dcmplx(IMSGHCUB-IMPGHCUB,REPGHCUB-RESGHCUB)
      GHCBU(2)=dcmplx(IMSGHCUB+IMPGHCUB,-RESGHCUB-REPGHCUB)
      GHCBC(1)=dcmplx(IMSGHCCB-IMPGHCCB,REPGHCCB-RESGHCCB)
      GHCBC(2)=dcmplx(IMSGHCCB+IMPGHCCB,-RESGHCCB-REPGHCCB)
      GHCBT(1)=dcmplx(IMSGHCTB-IMPGHCTB,REPGHCTB-RESGHCTB)
      GHCBT(2)=dcmplx(IMSGHCTB+IMPGHCTB,-RESGHCTB-REPGHCTB)

      GHCVEE(1)=dcmplx(IMSGHCVEE+IMPGHCVEE,RESGHCVEE+REPGHCVEE)
      GHCVEE(2)=dcmplx(IMSGHCVEE-IMPGHCVEE,RESGHCVEE-REPGHCVEE)
      GHCVMMU(1)=dcmplx(IMSGHCVMMU+IMPGHCVMMU,RESGHCVMMU+REPGHCVMMU)
      GHCVMMU(2)=dcmplx(IMSGHCVMMU-IMPGHCVMMU,RESGHCVMMU-REPGHCVMMU)
      GHCVTTA(1)=dcmplx(IMSGHCVTTA+IMPGHCVTTA,RESGHCVTTA+REPGHCVTTA)
      GHCVTTA(2)=dcmplx(IMSGHCVTTA-IMPGHCVTTA,RESGHCVTTA-REPGHCVTTA)

      GHCEVE(1)=dcmplx(IMSGHCVEE-IMPGHCVEE,REPGHCVEE-RESGHCVEE)
      GHCEVE(2)=dcmplx(IMSGHCVEE+IMPGHCVEE,-RESGHCVEE-REPGHCVEE)
      GHCMUVM(1)=dcmplx(IMSGHCVMMU-IMPGHCVMMU,REPGHCVMMU-RESGHCVMMU)
      GHCMUVM(2)=dcmplx(IMSGHCVMMU+IMPGHCVMMU,-RESGHCVMMU-REPGHCVMMU)
      GHCTAVT(1)=dcmplx(IMSGHCVTTA-IMPGHCVTTA,REPGHCVTTA-RESGHCVTTA)
      GHCTAVT(2)=dcmplx(IMSGHCVTTA+IMPGHCVTTA,-RESGHCVTTA-REPGHCVTTA)

      GWWH1=dcmplx(IMGWWH1,-REGWWH1)
      GWWH2=dcmplx(IMGWWH2,-REGWWH2)
      GZZH1=dcmplx(IMGZZH1,-REGZZH1)
      GZZH2=dcmplx(IMGZZH2,-REGZZH2)

      GAHCHC=dcmplx(-IMGAHCHC,REGAHCHC)

      GZH1H2=dcmplx(IMGZH1H2,-REGZH1H2)
      GZH1H3=dcmplx(IMGZH1H3,-REGZH1H3)
      GZH2H3=dcmplx(IMGZH2H3,-REGZH2H3)
      GZHCHC=dcmplx(-IMGZHCHC,REGZHCHC)

      GWPHCH1=dcmplx(-IMGWPHCH1,-REGWPHCH1)
      GWMH1HC=dcmplx(-IMGWPHCH1,REGWPHCH1)
      GWPHCH2=dcmplx(-IMGWPHCH2,-REGWPHCH2)
      GWMH2HC=dcmplx(-IMGWPHCH2,REGWPHCH2)
      GWPHCH3=dcmplx(-IMGWPHCH3,-REGWPHCH3)
      GWMH3HC=dcmplx(-IMGWPHCH3,REGWPHCH3)

      GH1H1H1=dcmplx(IMGH1H1H1,-REGH1H1H1)
      GH1H1H2=dcmplx(IMGH1H1H2,-REGH1H1H2)
      GH1H2H2=dcmplx(IMGH1H2H2,-REGH1H2H2)
      GH1H3H3=dcmplx(IMGH1H3H3,-REGH1H3H3)
      GH2H2H2=dcmplx(IMGH2H2H2,-REGH2H2H2)
      GH2H3H3=dcmplx(IMGH2H3H3,-REGH2H3H3)
      GH1HCHC=dcmplx(IMGH1HCHC,-REGH1HCHC)
      GH2HCHC=dcmplx(IMGH2HCHC,-REGH2HCHC)

      GAAHCHC=dcmplx(IMGAAHCHC,-REGAAHCHC)
      GAZHCHC=dcmplx(IMGAZHCHC,-REGAZHCHC)
      GAWPHCH1=dcmplx(IMGAWPHCH1,REGAWPHCH1)
      GAWPHCH2=dcmplx(IMGAWPHCH2,REGAWPHCH2)
      GAWPHCH3=dcmplx(IMGAWPHCH3,REGAWPHCH3)
      GAWMH1HC=dcmplx(IMGAWPHCH1,-REGAWPHCH1)
      GAWMH2HC=dcmplx(IMGAWPHCH2,-REGAWPHCH2)
      GAWMH3HC=dcmplx(IMGAWPHCH3,-REGAWPHCH3)

      GZZH1H1=dcmplx(IMGZZH1H1,-REGZZH1H1)
      GZZH2H2=dcmplx(IMGZZH2H2,-REGZZH2H2)
      GZZH3H3=dcmplx(IMGZZH3H3,-REGZZH3H3)
      GZZHCHC=dcmplx(IMGZZHCHC,-REGZZHCHC)
      GZWPHCH1=dcmplx(IMGZWPHCH1,REGZWPHCH1)
      GZWPHCH2=dcmplx(IMGZWPHCH2,REGZWPHCH2)
      GZWPHCH3=dcmplx(IMGZWPHCH3,REGZWPHCH3)
      GZWMH1HC=dcmplx(IMGZWPHCH1,-REGZWPHCH1)
      GZWMH2HC=dcmplx(IMGZWPHCH2,-REGZWPHCH2)
      GZWMH3HC=dcmplx(IMGZWPHCH3,-REGZWPHCH3)

      GWWH1H1=dcmplx(IMGWWH1H1,-REGWWH1H1)
      GWWH2H2=dcmplx(IMGWWH2H2,-REGWWH2H2)
      GWWH3H3=dcmplx(IMGWWH3H3,-REGWWH3H3)
      GWWHCHC=dcmplx(IMGWWHCHC,-REGWWHCHC)

c----------------------------
c end subroutine coupsm
c----------------------------


      return
      end

      
