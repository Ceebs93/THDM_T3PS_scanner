C-{{{ subroutine gghcall:

      subroutine gghcall(model,sqrtsin,mhiggsin,rmurmhin,
     &     rmufmhin,norderin,nscalpseudin,ncolliderin,gfermiin,
     &     mzin,mbmbin,mtopin,sigout,errout)
C$$$----------------------------------------------------------------------
c..
c..   See example-gghcall.f for an example how to use this subroutine.
c..
c..   Input:
c..
c..   sqrtsin = c.m.s. energy
c..   mhiggsin = Higgs mass
c..   rmurmhin = renormalization scale in units of Higgs mass
c..   rmufmhin = factorization scale in units of Higgs mass
c..   norderin = order of calculation [0=LO, 1=NLO, 2=NNLO]
c..   nscalpseudin = scalar [0] / pseudo-scalar [1] Higgs production
c..   ncolliderin = pp [0] / pp-bar [1] collider
c..   gfermiin = Fermi coupling constant
c..   mzin = Z mass
c..   mbmbin = mb(mb) bottom quark mass
c..   mtopin = top quark pole mass
c..   
c..   Output:
c..   siglo =   LO total cross section w/o EW corrections
c..   signlo =  NLO total cross section w/o EW corrections
c..   signnlo = NNLO total cross section w/o EW corrections
c..
C$$$----------------------------------------------------------------------



      implicit real*8(a-h,o-z)
      real*8 mhiggsin,mzin,mbmbin,mtopin
      real*8 sigout(0:3,-1:5),errout(0:3,-1:5)
      integer model,i,j
      include '../commons/common-sigma.f'
c..
c..   call function in order to use ggh@nnlo as a library
c..
      include '../commons/common-readdata.f'
c      include '../commons/common-slha.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-expand.f'
      include '../commons/common-errors.f'
      include '../commons/common-citations.f'

      citations(7) = 1
      
c--   which subprocess to evaluate?  [0 = sum of all subprocesses]
c$$$      nsubprocggh = 0

c--   order of calculation:
      norder = norderin

      if (norder.ge.2) then
         if (lpseudo) then
            citations(36) = 1
            citations(37) = 2
         else
            citations(23) = 2
            citations(24) = 2
         endif
      endif
      if (norder.ge.3) then
         if (lpseudo) then
            stop 'gghcall.f 1'
         else
            citations(45) = 1
            citations(46) = 2
            citations(47) = 1
            citations(48) = 1
         endif
      endif

c--   scalar/pseudo-scalar?
      nscalpseud = nscalpseudin
      if (nscalpseud.eq.0) then
         lpseudo = .false.
      elseif (nscalpseud.eq.1) then
         lpseudo = .true.
      else
         call printdieggh('HNNLO1 (2) must be 0 or 1')
      endif

      ncollider = ncolliderin
c--   proton-proton or proton-antiproton?
      if (ncollider.eq.0) then
         ppbar = .false.
      else
         ppbar = .true.
      endif

      nmodel = model

      gfermird = gfermiin
      mzrd = mzin
      mbmbrd = mbmbin
      mtoprd = mtopin

      sqscmsrd = sqrtsin
      rmurmh  = rmurmhin
      rmufmh  = rmufmhin
      
      mhiggsrd = mhiggsin

c..   these are just to switch on/off certain contributions:
      gth = 1.d0
      gbh = 0.d0
      gth11 = 0.d0
      gth12 = 0.d0
      gth22 = 0.d0
      gbh11 = 0.d0
      gbh12 = 0.d0
      gbh22 = 0.d0

      call initggh()

c$$$      print*,'>>>>>>>>>> WARNING!!!!!!!!!!!! gghcall.f'
c$$$      open(unit=1122,file='heino.out')
c$$$      xtl=-3.d0
c$$$      do i=1,290
c$$$         xtl=xtl+.01d0
c$$$         write(1122,*) 10**xtl,sgg3x0(10**xtl),sgg3exp(10**xtl)
c$$$      enddo
c$$$      close(1122)
c$$$      print*,'<<<<<<<<<< WARNING!!!!!!!!!!!! gghcall.f'

      call evalsigmaggh()

      do i=0,3,1
         do j=-1,5,1
            sigout(i,j)=sigma(i,j)
            errout(i,j)=sigerr(i,j)
         enddo
      enddo
      
c$$$      print*,'>>>>',sigout

c$$$  siglo = sigma(0,0)
c$$$      signlo = sigma(1,0)
c$$$      signnlo = sigma(2,0)
c$$$      errlo = sigerr(0,0)
c$$$      errnlo = sigerr(1,0)
c$$$      errnnlo = sigerr(2,0)
         
      end

C-}}}
