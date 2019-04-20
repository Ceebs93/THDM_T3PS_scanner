C-{{{ subroutine evalsigma:

      subroutine evalsigmabbh()
c..
c..   Basically nothing is done here, only subroutine calls.
c..   
      implicit none
      integer i
      real*8 apimb,apiref,muref,mbref
      real*8 sigtotbbh(0:norder)
      logical lhxswg
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'
      
      if (lmyapimz.and.(apimz.ne.myapimz)) then
         call printwarnbbh('using alphas value incompatible with pdfs')
         write(6,*) '   correct alphas off by',(apimz-myapimz)/apimz
     &        *100,'%'
         apimz = myapimz
      endif
      
c..   running alpha_s with my own routines instead of LHAPDF:
      if (hxswg.eq.1) then
c..   first run everything at highest order to muref:
         muref=mh
         call runalpha(apimz,mz,muref,nfrun,4,0,apiref)
         call runalpha(apiref,muref,murfacbbh*mh,nfrun,norder+1,0
     &        ,apimurbbh)
         call runalpha(apimz,mz,qqhscale,nfrun,4,0,apimb)
         call runmass(qqhyuk,apimb,apiref,nfrun,4,mbref)
         call runmass(mbref,apiref,apimurbbh,nfrun,norder+1
     &        ,mbMSbarmuRbbh)
      else
         call runalpha(apimz,mz,murfacbbh*mh,nfrun,norder+1,0,apimurbbh)
         call runalpha(apimz,mz,qqhscale,nfrun,norder+1,0,apimb)
         call runmass(qqhyuk,apimb,apimurbbh,nfrun,norder+1
     &        ,mbMSbarmuRbbh)
      endif

      call normalizationbbh()
      call dointegralsbbh()
      
      do i=0,norder
         sigmabbh(i) = sallbbh(i)
      enddo

      if (errnorm.eq.1) then
         error = 1
      endif
      if (warnnorm.eq.1) then
         warning = 1
      endif

      end

C-}}}
C-{{{ subroutine dointegrals:

      subroutine dointegralsbbh()
c..
c..   NNLO cross section (sall).
c..   Partly fills common/sigma/.
c..
      implicit none
      real*8 ddsoft1,errsoft1
      real*8 ddsoft2,errsoft2,del2,errdel2,hard1,hard2,err1,err2
      real*8 sall0,sall1,sall2,sall3,sall4
      real*8 serr0,serr1,serr2,serr3,serr4
      real*8 dtsub1bbh,dtsub2bbh,delta1bbh,delta2bbh
      include '../commons/common-vars.f'
      include '../commons/common-sigma.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-vegpar.f'

c--   LO:   ---

      call intdelbbh(del2,errdel2)

c..   LO is non-zero only for the bb-bar subprocess:
      if (nsubprocbbh.le.1) then
         sall0 = del2
         serr0 = errdel2
      else
         sall0 = 0.d0
         serr0 = 0.d0
      endif

c--   NLO:  ---
      if (norder.gt.0) then

c--   D-terms:
         call soft1bbh(ddsoft1,errsoft1)

c--   hard contribution:
         call evalhard1bbh(hard1,err1)
         
         if (nsubprocbbh.le.1) then
            if (nsubprocbbh.ge.-2) then
               sall1 = ( delta1bbh() )*del2
               serr1 = ( delta1bbh()*errdel2 )**2
            endif
            if (nsubprocbbh.ge.-1) then
               sall1 = sall1 + ( dtsub1bbh() )*del2 + ddsoft1
               serr1 = serr1 + ( dtsub1bbh()*errdel2 )**2 + errsoft1**2
            endif
            if (nsubprocbbh.ge.0) then
               sall1 = sall1 + hard1
               serr1 = serr1 + err1**2
            endif
            err1 = dsqrt(err1)
         else
            sall1 = hard1
            serr1 = err1
         endif
      else
         sall1 = 0.d0
         serr1 = 0.d0
      endif

      if (norder.gt.1) then
c--   NNLO: ---

c--   D-terms:
         call soft2bbh(ddsoft2,errsoft2)

c--   hard contribution:
         call evalhard2bbh(hard2,err2)

         if (nsubprocbbh.le.1) then
            if (nsubprocbbh.ge.-2) then
               sall2 = ( delta2bbh() )*del2
               serr2 = ( delta2bbh() *errdel2 )**2
            endif
            if (nsubprocbbh.ge.-1) then
               sall2 = sall2 + ( dtsub2bbh() )*del2 + ddsoft2
               serr2 = serr2 + ( dtsub2bbh()*errdel2 )**2 + errsoft2**2
            endif
            if (nsubprocbbh.ge.0) then
               sall2 = sall2 + hard2
               serr2 = serr2 + err2**2
            endif
            serr2 = dsqrt(serr2)
         else
            sall2 = hard2
            serr2 = err2
         endif
      else
         sall2 = 0.d0
         serr2 = 0.d0
      endif

      sall3=0.d0
      serr3=0.d0
      sall4=0.d0
      serr4=0.d0

!C       if (norder.gt.2) then
!C c--   N^3LO: ---
!C
!C c--   D-terms:
!C          call soft3(ddsoft3,errsoft3)
!C
!C c--   hard contribution is not available yet:
!C          hard3 = 0.d0
!C
!C          if (nsubprocbbh.le.1) then
!C          sall3 = ( delta3() + dtsub3() )*del2
!C      &        + ddsoft3 + hard3
!C          else
!C             sall3 = hard3
!C          endif
!C       else
!C          sall3 = 0d0
!C       endif
!C
!C       if (norder.gt.3) then
!C c--   N^4LO: ---
!C
!C c--   D-terms:
!C          call soft4(ddsoft4,errsoft4)
!C
!C c--   hard contribution is not available yet:
!C          hard4 = 0.d0
!C
!C          if (nsubprocbbh.le.1) then
!C          sall4 = ( delta4() + dtsub4() )*del2
!C      &        + ddsoft4 + hard4
!C          else
!C             sall4 = hard4
!C          endif
!C       else
!C          sall4 = 0d0
!C       endif

      sallbbh(0) = 0.d0
      sallbbh(1) = 0.d0
      sallbbh(2) = 0.d0
      sallbbh(3) = 0.d0
      sallbbh(4) = 0.d0
      serrbbh(0) = 0.d0
      serrbbh(1) = 0.d0
      serrbbh(2) = 0.d0
      serrbbh(3) = 0.d0
      serrbbh(4) = 0.d0

      sallbbh(0) = prefac*( sall0 )
      serrbbh(0) = prefac*( serr0 )
      if (norder.gt.0) then
         sallbbh(1) = prefac*( sall0 + apimurbbh*sall1 )
         serrbbh(1) = prefac*dsqrt( serr0**2 + (apimurbbh*serr1)**2 )
      endif
      if (norder.gt.1) then
         sallbbh(2) = prefac*( sall0 + apimurbbh*sall1 + apimurbbh**2
     &        *sall2 )
         serrbbh(2) = prefac*dsqrt( serr0**2 + (apimurbbh*serr1)**2 
     &        + (apimurbbh**2*serr2)**2 )
      endif
      if (norder.gt.2) then
         sallbbh(3) = prefac*( sall0 
     &        + apimurbbh*sall1 
     &        + apimurbbh**2*sall2 
     &        + apimurbbh**3*sall3 )
         serrbbh(3) = prefac*dsqrt( serr0**2
     &        + (apimurbbh*serr1)**2
     &        + (apimurbbh**2*serr2)**2
     &        + (apimurbbh**3*serr3)**2 )
      endif
      if (norder.gt.3) then
         sallbbh(4) = prefac*( sall0 
     &        + apimurbbh*sall1 
     &        + apimurbbh**2*sall2 
     &        + apimurbbh**3*sall3 
     &        + apimurbbh**4*sall4 
     &        )
         serrbbh(4) = prefac*dsqrt( serr0**2
     &        + (apimurbbh*serr1)**2
     &        + (apimurbbh**2*serr2)**2
     &        + (apimurbbh**3*serr3)**2
     &        + (apimurbbh**4*serr4)**2
     &        )
      endif

      end

C-}}}
C-{{{ subroutine initbbh:

      subroutine initbbh()
c..
c..   Initialize some parameters.
c..   
      implicit none
c      include '../commons/common-slha.f'
      include '../commons/common-readdata.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'

      pi = 3.14159265358979323846264338328d0
      z2 = 1.6449340668482264364724d0
      z3 = 1.2020569031595942853997d0
      z4 = 1.0823232337111381915160d0
      z5 = 1.0369277551433699263314d0
      ca = 3.d0
      cf = 1.33333333333333333333d0
      tr = 0.5d0
     
      lsloppy = lsloppyrd

c$$$      pdfname = pdfnamerd
c$$$      npdfstart = npdfstartrd
c$$$      npdfend = npdfstartrd
c$$$      npdfend = npdfendrd
c$$$      if (npdfstart.ne.npdfend) then
c$$$         call printdie('currently only one PDF set per run allowed -')
c$$$      endif

c--   allow only for a single pdf-set for the moment:
c      if (npdfstart.ne.npdfend) then
c         call printdie('Averaging over PDF-sets temporarily disabled.')
c      endif

c--   initialize errors and warnings:
      error = 0
      errnorm = 0
      warning = 0
      warnnorm = 0

c--   order of calculation:
      norder = norderrd
      
c--   proton-proton or proton-antiproton?
      if (ncolliderrd.eq.0) then
         ppbar = .false.
      else
         ppbar = .true.
      endif

c--   scalar/pseudo-scalar?
      if (nscalpseudrd.eq.0) then
         lpseudo = .false.
      else
         lpseudo = .true.
      endif

c--   cms energy:
      sqrts = sqscmsrd

c--   Higgs mass:
      mh = mhiggsrd

c--   bottom mass:
      mbmb = mbottomrd

c--   bottom Yukawa couplings:
c      gb = ybottomrd

c--   number of massless flavors:
      nf=5.d0
      nfrun=5.d0

c--   Vegas parameters:
c$$$      acc = 1.d-8       !  accuracy for vegas
c$$$      itmx1=5         !  interations for vegas run 1
c$$$      ncall1=2000      !  calls for vegas run 1
c$$$      itmx2=2         !  interations for vegas run 2
c$$$      ncall2=5000     !  calls for vegas run 2
c$$$      nprnv=0          !  =1 -- verbose mode
c$$$      lveg1 = .true.  !  if .false., run vegas only once!

c--   constants and normalization:
      gfermi = gfermird
      mz = mzrd

      tauh = mh**2/sqrts**2

c--   logarithms:
      lfh = 2*dlog(mufmhrd)
      lfr = 2*dlog(mufmhrd/murmhrd)
      lrt = lft + lfr
      rmur = murmhrd*mh
      rmuf = mufmhrd*mh


c--   SM or non-SM?  (0: SM  --  1: other model)
      if (nc1rd.eq.1) then
         lstdmodel = .false.
      else
         lstdmodel = .true.
      endif

c--   top and bottom Yukawa couplings:
      gth = gthrd
      gbh = gbhrd
      gth11 = 0.d0
      gth22 = 0.d0
      gth12 = 0.d0
      gth21 = 0.d0

c--   which subprocess to evaluate?  [0 = sum of all subprocesses]
c$$$      nsubprocbbh = nsubprocrd

c      call rluxgo(3,12348271,0d0,0d0)    !  runlux initialization (optional)

      end

C-}}}
C-{{{ subroutine polemass:

      subroutine polemassbbh(mqmq,apimq,nfh,nloop,mqpole)
c..
c..   Computes the pole mass mqpole from the MS-bar mass mqmq = mq(mq).
c..   apimq = alpha_s(mqmq)/pi.
c..
      implicit none
      integer nloop,myloop
      real*8 nf,nfh,mqmq,apimq,mqpole
      include '../commons/common-consts.f'

      nf = nfh-1
      
      myloop=nloop

      if (myloop.gt.2) then
         write(6,*) '<function polemass>: nloop = ',myloop,
     &        ' not implemented '
         write(6,*) '      using 2-loop expression for mqpole'
     &        
         myloop=2
      endif

      if (myloop.eq.0) then
         mqpole = mqmq
      elseif (myloop.eq.1) then
         mqpole = mqmq*( 1 + apimq*4/3.d0 )
      elseif (myloop.eq.2) then
         mqpole = mqmq*( 1 + apimq*4/3.d0 + apimq**2 * (
     &        307/32.d0 - (71*nf)/144.d0 + 2*z2 + (2*dlog(2.d0)*z2)/3.d0
     &        - (nf*z2)/3.d0 -z3/6.d0 ) )
      else
         write(6,*) '<function polemass>: ERROR'
         stop
      endif

      end

C-}}}
C-{{{ subroutine printnotice:

      subroutine printnoticebbh(unit)

      implicit none
      integer unit

      write(unit,2011) '# -----------------------'
      write(unit,2011) '# IMPORTANT NOTE'
      write(unit,2011) '# -----------------------'
      write(unit,2011) 
     &     '# If you use bbh@nnlo (or parts of it, or a modified '//
     &     'version of it)'
      write(unit,2011)
     &     '# for a publication, you have to refer to the paper'
      write(unit,2011) '#' 
      write(unit,2011) 
     &     '#   Robert V. Harlander and William B. Kilgore,'
      write(unit,2011) 
     &     '#     "Higgs boson production in bottom quark fusion'
      write(unit,2011) 
     &     '#     at next-to-next-to-leading order",'
      write(unit,2011) 
     &     '#     Phys. Rev. D68, 013001 (2003), arXiv:hep-ph/0304035'
      write(unit,2011) '#'
      write(unit,2011) 
     &     '# If you re-distribute bbh@nnlo (or parts of it, '//
     &     'or a modified'
      write(unit,2011) 
     &     '# version of it), or distribute code that links bbh@nnlo '//
     &     'or is based'
      write(unit,2011) 
     &     '# on results of bbh@nnlo (e.g. interpolations), '//
     &     'you are required to:'
      write(unit,2011)
     &     '# * inform the author of bbh@nnlo '//
     &     '(robert.harlander@uni-wuppertal.de)'
      write(unit,2011)
     &     '# * clearly refer to the original source in your code '//
     &     'and its output'
      write(unit,2011) 
     &     '# * point out to any user that she or he must '//
     &     'refer to the above paper'
      write(unit,2011) '#   in publications produced with the help '//
     &     'of this code'
      write(unit,2011) '#'
      write(unit,2011) 
     &     '# If you disagree with any of the above, do not use '//
     &     'bbh@nnlo.'
      write(unit,2011) 
     &     '# '//
     &     ''

 2011 format(A)

      end

C-}}}
C-{{{ subroutine printwarn:

      subroutine printwarnbbh(strng)

      implicit none
      character*(*) strng

      write(6,*) 'bbh@nnlo (WARNING): ',strng

      end

C-}}}
C-{{{ subroutine normalization:

      subroutine normalizationbbh()

      implicit none
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'
!      data gev2pb/.38937966d+9/   

      gev2pb = .38937966d+9 !  conversion  1/GeV^2 -> pb

      errnorm = 0

      rho0 = pi/6.d0*dsqrt(2.d0)*gfermi/mh**2*gev2pb
      prefac = rho0*gbh**2*mbMSbarmuRbbh**2

      end

C-}}}

