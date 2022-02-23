! This file is part of SusHi.
! 
! Interface for external codes.
! 

C-{{{ subroutine SusHicall

      subroutine SusHicall(jfilein_in,muRfac_in,htl_in,model_out,GF_out,
     &mt_out,mb_out,mbos_out,mc_out,mh_out,hpdg_out,ew_out)

C$$$----------------------------------------------------------------------
c..
c..   Input:
c..
c..   jfilein_in = input file
c..   muRfac_in = renormalization scale muR/mh
c..   htl_in = heavy top limit [0=do not use htl, 1=use htl]
c..
c..   Output:
c..
c..   model_out = model [0=SM, 1=MSSM, 2=2HDM, 3=NMSSM]
c..   GF_out = Fermi coupling constant
c..   mt_out = top quark pole mass
c..   mb_out = bottom quark mass used in loop
c..   mbos_out = bottom quark pole mass
c..   mc_out = charm quark mass
c..   mh_out = Higgs mass
c..   hpdg_out = Higgs type [PDG code]
c..   ew_out = electroweak contributions [0=no, 1=light quarks, 2=SM EW factor]
c..
C$$$----------------------------------------------------------------------

      implicit none
      character jfilein_in*60
      integer htl_in,model_out,hpdg_out,ew_out
      double precision muRfac_in,GF_out,mt_out,mb_out,mbos_out,mc_out,
     &mh_out
      logical first

      !common variables used by various routines
      include '../commons/common-inputoutput.f'
      include '../commons/common-vars.f'
      include '../commons/common-quark.f'
      include '../commons/common-ren.f'
      include '../commons/common-flags.f'

      save first
      data first /.true./

      jfilein = jfilein_in

      if(first) then

         first = .false.
         extflag = .true.

         murfacggh = murfac_in
         call initsushi(htl_in)
         call initcouplings()

         if (model.eq.2) call renormalizeSM

         if (model.eq.0) then
            !set couplings to SM values:
            gc = yukfac(1)
            gt = yukfac(2)
            gb = yukfac(3)
            gtp = yukfac(4)
            gbp = yukfac(5)
            call renormalizeSM
         end if

c$$$         if (gt.ne.0.d0) then
c$$$            write(*,*) "Top Quark enabled"
c$$$         endif
c$$$         if (gb.ne.0.d0) then
c$$$            write(*,*) "Bottom Quark enabled"
c$$$         endif
c$$$         if (gc.ne.0.d0) then
c$$$            write(*,*) "Charm Quark enabled"
c$$$         endif

      endif

      if (pseudo.eq.0) then
         hpdg_out = 25 + 10*(Sind - 1)
      else if (pseudo.eq.1) then
         hpdg_out = 26 + 10*(Sind - 1)
      endif

      ew_out = ew
      model_out = model
      GF_out = GF
      if(gt.eq.0.d0) then
         mt_out = -1000000d0
      else
         mt_out = dsqrt(mt2)
c         mt_out = mt
      endif
      if(gb.eq.0.d0) then
         mb_out = -1000000d0
      else
         mb_out = dsqrt(mb2)
c         mb_out = mb
      endif
      if(gc.eq.0.d0) then
         mc_out = -1000000d0
      else
         mc_out = dsqrt(mc2)
      endif
      mh_out = mh
      mbos_out = mbos
      end

C-}}}

!end of file
