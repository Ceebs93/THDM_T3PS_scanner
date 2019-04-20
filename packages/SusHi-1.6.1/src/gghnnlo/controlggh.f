C-{{{ subroutine evalsigmaggh:

      subroutine evalsigmaggh()
c..
c..   Basically nothing is done here, only subroutine calls.
c..   
      implicit none
      integer i,j
      real*8 apimb
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'
      real*8 sigtot(0:3,-1:5)

      lmyapimz=.false.
      if (lmyapimz.and.(apimz.ne.myapimz)) then
         apimz = myapimz
         call printwarnggh('using alphas value incompatible with pdfs')
      endif

c..   running alpha_s with my own routines instead of LHAPDF:
      call runalpha(apimz,mz,mbmb,nf,norder+1,0,apimb)
      call runalpha(apimz,mz,rmur,nf,norder+1,0,apimuR)
      call polemass(mbmb,apimb,nf,norder,mbpole)
c     call runmass(mbmb,apimb,apimuR,nf,norder+1,mbMSbarmuR)

      call normalizationggh()
      call dodointggh(.false.)

      sigma = 0.d0
      
      do i=0,norder
         do j=-1,5
            sigma(i,j) = sall(i,j)
         enddo
      enddo

      if (errnorm.eq.1) then
         error = 1
      endif

      return
      end

C-}}}
C-{{{ subroutine dodointggh()

      subroutine dodointggh(onlytopho)

      implicit none
      logical lbll(3),onlytopho
      integer nmt(0:3,0:5)
      integer nfclo
      real*8 c1re(0:3),c1im(0:3),c2re(0:3),c2im(0:3),prfc
      real*8 sallx(0:3,-1:5),serrx(0:3,-1:5)
      include '../commons/common-sigma.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-expand.f'
      
      c1re=0.d0
      c2re=0.d0
      c1im=0.d0
      c2im=0.d0
      sigerr=0.d0
      sall=0.d0
      nmt=nmtlim
      lbll=lball
      nfclo=nfaclo
      
c$$$  if (any(nmtlim.ne.0)) then
      if (onlytopho) then
         write(6,*) 'calculate SM result including 1/mt terms'
         c1re(0) = c1sm0
         c1re(1) = c1sm1
         c1re(2) = c1sm2
         c1re(3) = c1sm3
         c2re(0) = 0.d0
         c2re(1) = c2sm1
         
         lbll=lball
         
         nmt=nmtlim
         nfclo = nfaclo
         prfc = rho0*(apimuR/3.d0)**2

         call dointggh(prfc,c1re,c1im,c2re,c2im,nmt,nfclo,lbll,sallx
     &        ,serrx)
         
         sall=sallx
         sigerr=serrx
         
         write(6,*) 'now calculate SM result without 1/mt terms'
         nmt=0
         nfclo = nfaclo
      
         call dointggh(prfc,c1re,c1im,c2re,c2im,nmt,nfclo,lbll,sallx
     &        ,serrx)
         
         sall=sall-sallx
         sigerr=dsqrt(sigerr**2 + serrx**2)

         write(6,*) 'calculate full result without 1/mt terms'
         nmt=0.d0
         c1re(0) = dreal(comc1eff0)
         c1im(0) = dimag(comc1eff0)
         c1re(1) = c1eff1
         c1re(2) = c1eff2
         c1re(3) = c1eff3
         c2re(1) = c2eff1
         
         call dointggh(prefac/ataut2,c1re,c1im,c2re,c2im,nmt,nfclo,lbll
     &        ,sallx,serrx)
         
         sall = sall+sallx
         sigerr=dsqrt(sigerr**2 + serrx**2)
      else
         nmt=nmtlim
         c1re(0) = dreal(comc1eff0)
         c1im(0) = dimag(comc1eff0)
         c1re(1) = c1eff1
         c1re(2) = c1eff2
         c1re(3) = c1eff3
         c2re(1) = c2eff1
         
         call dointggh(prefac,c1re,c1im,c2re,c2im,nmt,nfclo,lbll
     &        ,sallx,serrx)
        
         sall = sallx
         sigerr=serrx
      endif
         
      return
      end

C-}}}
C-{{{ subroutine dointggh:

      subroutine dointggh(prfc,c1re,c1im,c2re,c2im,nmt,nfclo,lbll,sallx
     &     ,serrx)
c..
c..   NNLO cross section (sall).
c..   Partly fills common/sigma/.
c..
c..   del0, del1, del2 are the results for Wilson coefficient=1.
c..   sigma0, sigma1, sigma2  are the true cross sections
c..   (devided by the leading order), with proper coefficient functions
c..   included.
c..
      implicit none
      integer i,nfclo,nmt(0:3,0:5),nmtsave(0:3,0:5)
      logical lballsave(3),lbll(3)
      complex*16 atau
      real*8 amplosq
      real*8 ddsoft1mt0,prfc,taut,c1abs2
      real*8 c1re(0:3),c1im(0:3),c2re(0:3),c2im(0:3),ddsoft1,ddsoft2
     &     ,ddsoft3
      real*8 c1fac(1:3),c2fac(1:3),c1factor(0:3),dumre(0:3),dumim(0:3)
      real*8 errsoft1,errsoft2,errsoft3,delt2,errdelt2
      real*8 dtsub1,dtsub2,dtsub3,atauexph,psigma1,delta1,delta2,delta3
      real*8 hard1(-1:5),err1(-1:5),hard2(-1:5),err2(-1:5),hard3(-1:5)
     &     ,err3(-1:5)
      real*8 errdel0(-1:5),errdel1(-1:5),errdel2(-1:5),errdel3(-1:5)
      real*8 errsig0(-1:5),errsig1(-1:5),errsig2(-1:5),errsig3(-1:5)
      real*8 sigma0(-1:5),sigma1(-1:5),sigma2(-1:5),sigma3(-1:5)
      real*8 sallx(0:3,-1:5),serrx(0:3,-1:5)
      real*8 del0(-1:5),del1(-1:5),del2(-1:5), del3(-1:5)
      include '../commons/common-vars.f'
      include '../commons/common-sigma.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-expand.f'
      include '../commons/common-vegpar.f'

      nmtsave = nmtlim
      nmtlim=nmt
      lballsave=lball

c..   expand c1^2*<o1o1> perturbatively or not?
      c1abs2 = c1re(0)**2 + c1im(0)**2
      c1fac=0.d0
      c2fac=0.d0
      c2fac(1) = (c2re(1)*c1re(0) + c2im(1)*c1im(0))/c1abs2
      if (lexpandc1o1) then
         c1fac(1) = 2*(c1re(1)*c1re(0)+c1im(1)*c1im(0))/c1abs2
         c1fac(2) = (c1re(1)**2 + 2*(c1re(2)*c1re(0) + c1im(2)*c1im(0)))
     &        /c1abs2
         c1fac(3) = 2*(c1re(1)*c1re(2) + c1im(1)*c1im(2) + c1re(3)
     &        *c1re(0) + c1im(3)*c1im(0))/c1abs2
      endif
      
c--   LO:   ---

      del0 = 0.d0
      del1 = 0.d0
      del2 = 0.d0
      del3 = 0.d0
      sigma0 = 0.d0
      sigma1 = 0.d0
      sigma2 = 0.d0
      sigma3 = 0.d0
      hard1 = 0.d0
      hard2 = 0.d0
      hard3 = 0.d0
      errsig0 = 0.d0
      errsig1 = 0.d0
      errsig2 = 0.d0
      errsig3 = 0.d0
      errdel0 = 0.d0
      errdel1 = 0.d0
      errdel2 = 0.d0
      errdel3 = 0.d0

      call intdel(delt2,errdelt2)

c.. 0=full, 1=gg (real), 2=qg, 3=qqbar
      del0(-1) = delt2
      errdel0(-1) = errdelt2
      if ((nsubprocggh.eq.0).or.(nsubprocggh.eq.10)) then
         del0(0) = del0(-1)
         errdel0(0) = errdel0(-1)
      endif
      do i=-1,5
         sigma0(i) = del0(i)
         errsig0(i) = errdel0(i)
      enddo

c--   NLO:  ---
      if (norder.ge.1) then

c--   D-terms:
         call soft1(ddsoft1,errsoft1)
         
c--   hard contribution:
         call evalhard1(hard1,err1)

         del1(-1) = ( delta1(nmt(1,1)) + dtsub1(nmt(1,1)) )*delt2
     &        +ddsoft1
         errdel1(-1) = dsqrt((( delta1(nmt(1,1)) + dtsub1(nmt(1
     &        ,1)) )*errdelt2)**2+ errsoft1**2)
         if ((nsubprocggh.eq.0).or.(nsubprocggh.eq.10)) then
            del1(0) = del1(-1)
            errdel1(0) = errdel1(-1)
         endif
         
         do i=-1,5
            del1(i) = del1(i) + hard1(i)
            errdel1(i) = dsqrt(errdel1(i)**2 + err1(i)**2)
         enddo
         do i=-1,5
            sigma1(i) = del1(i) + c1fac(1)*del0(i)
            errsig1(i) = dsqrt(errdel1(i)**2 
     &           + (c1fac(1)*errdel0(i))**2)
         enddo
         
      endif

      if (norder.ge.2) then
c--   NNLO: ---

c--   D-terms:
         call soft2(ddsoft2,errsoft2)

c--   hard contribution:
         call evalhard2(hard2,err2)
         
         del2(-1) = ( delta2(nmt(2,1)) + dtsub2(nmt(2,1)) )*delt2
     &        +ddsoft2
         errdel2(-1) = dsqrt( (( delta2(nmt(2,1)) + dtsub2(nmt(2
     &        ,1)) )*errdelt2)**2+ errsoft2**2 )
         if ((nsubprocggh.eq.0).or.(nsubprocggh.eq.10)) then
            del2(0) = del2(-1)
            errdel2(0) = errdel2(-1)
         endif

         do i=-1,5
            del2(i) = del2(i) + hard2(i)
         enddo
          do i=-1,5
             sigma2(i) = del2(i) + c1fac(1)*del1(i) + c1fac(2)*del0(i)
             errsig2(i) = dsqrt( errdel2(i)**2 
     &            + (c1fac(1)*errdel1(i))**2
     &            + (c1fac(2)*errdel0(i))**2)
         enddo

         psigma1 = nf*sigma0(0)
         if (lpseudo) then
            sigma2(-1) = sigma2(-1) + c2fac(1)*psigma1
            sigma2(0) = sigma2(0) + c2fac(1)*psigma1
         endif
      endif

      if (norder.ge.3) then
c--   N3LO: ---

c--   D-terms:
         call soft3(ddsoft3,errsoft3)

c--   hard contribution:
         call evalhard3(hard3,err3)

         del3(-1) = ( delta3() + dtsub3() )*delt2 + ddsoft3
         errdel3(-1) = dsqrt( (( delta3() + dtsub3() )*errdelt2)**2
     &        + errsoft3**2 )
         if ((nsubprocggh.eq.0).or.(nsubprocggh.eq.10)) then
            del3(0) = del3(-1)
            errdel3(0) = errdel3(-1)
         endif

         do i=-1,5
            del3(i) = del3(i) + hard3(i)
         enddo

         do i=-1,5
            sigma3(i) = del3(i) + c1fac(1)*del2(i) + c1fac(2)*del1(i)
     &           +c1fac(3)*del0(i)
            errsig3(i) = dsqrt(errdel3(i)**2 + (c1fac(1)*errdel2(i))**2
     &           +(c1fac(2)*errdel1(i))**2 +  (c1fac(3)*errdel0(i))
     &           **2)
         enddo
      endif
      
      dumre(0) = c1re(0)
      dumim(0) = c1im(0)
      c1factor=1.d0
      if (.not.lexpandc1o1) then
         do i=1,norder
            dumre(i) = dumre(i-1) + apimur**i*c1re(i)
            dumim(i) = dumim(i-1) + apimur**i*c1im(i)
         enddo
         do i=0,norder
            c1factor(i) = (dumre(i)**2 + dumim(i)**2)/c1abs2
         enddo
      endif
      
      do i=-1,5
         sallx(0,i) = 0.d0
         sallx(1,i) = 0.d0
         sallx(2,i) = 0.d0
         sallx(3,i) = 0.d0

         taut = 4*mt**2/mh**2
         if (nfclo.eq.-1) then
c..   plain expansion in 1/mt:
            sallx(0,i) = sigma0(i)*amplosq( nmt(0,0) )
         else
c..   full mt-dependence:
            sallx(0,i) = sigma0(i)*amplosq( -2 )
         endif
         sallx(0,i) = prfc*c1factor(0)*sallx(0,i)
         serrx(0,i) = prfc*c1factor(0)*errsig0(i)*amplosq( -2 )

         if (norder.gt.0) then
            if (nfclo.eq.-1) then
c..   plain expansion in 1/mt:
               sallx(1,i) =
     &              sigma0(i)*amplosq( nmt(0,0) )
     &              + c1abs2*apimur*sigma1(i)
            elseif (nfclo.eq.0) then
c..   keep LO mt-dep, expand NLO in 1/mt:
               sallx(1,i) = sigma0(i)*amplosq( -2 )
     &              + c1abs2*apimur*sigma1(i)
            else
c..   factor out LO mt-dep:
               sallx(1,i) = amplosq( -2 )*(
     &              sigma0(i)
     &              + apimur*sigma1(i)/atauexph( taut, nmt(1,0)) )
            endif
            sallx(1,i) = prfc*c1factor(1)*sallx(1,i)
            serrx(1,i) = prfc*c1factor(1)*dsqrt(errsig0(i)
     &           **2 +(apimur*errsig1(i))**2)*amplosq( -2 )
         endif

         if (norder.gt.1) then
            if (nfclo.eq.-1) then
c..   plain expansion in 1/mt:
               sallx(2,i) =
     &              sigma0(i) * amplosq( nmt(0,0) )
     &              + c1abs2*(apimur*sigma1(i)
     &              + apimur**2*sigma2(i))
            elseif (nfclo.eq.0) then
c..   keep LO mt-dep, expand NLO and NNLO in 1/mt:
               sallx(2,i) = sigma0(i)*amplosq( -2 )
     &              + c1abs2*(apimur*sigma1(i)
     &              + apimur**2*sigma2(i))
            elseif (nfclo.eq.1) then
c..   factor out LO mt-dep through NLO, expand NNLO in 1/mt:
               sallx(2,i) = amplosq( -2 )*(sigma0(i)
     &              + apimur*sigma1(i)/atauexph( taut, nmt(1,0)))
     &              + c1abs2*apimur**2*sigma2(i)
            else
c..   factor out LO mt-dep through NNLO:
               sallx(2,i) = amplosq( -2 )*(sigma0(i)
     &              + apimur*sigma1(i) / atauexph( taut,nmt(1,0) )
     &              + apimur**2*sigma2(i) / atauexph( taut,nmt(2,0) ) )
            endif
            sallx(2,i) = prfc*c1factor(2)*sallx(2,i)
            serrx(2,i) = prfc*c1factor(2)*dsqrt(errsig0(i)**2 +
     &           (apimur*errsig1(i))**2 + (apimur**2*errsig2(i)))
     &           *amplosq( -2 )
         endif

         if (norder.ge.3) then
            if (nfclo.eq.-1) then
c..   plain expansion in 1/mt:
               sallx(3,i) = sigma0(i) * amplosq( nmt(0,0) )
     &              + c1abs2*(apimur*sigma1(i)
     &              + apimur**2*sigma2(i)
     &              + apimur**3*sigma3(i))
            elseif (nfclo.eq.0) then
c..   keep LO mt-dep, expand NLO,NNLO,N3LO in 1/mt:
               sallx(3,i) = sigma0(i)*amplosq( -2 )
     &              + c1abs2*(apimur*sigma1(i)
     &              + apimur**2*sigma2(i) 
     &              + apimur**3*sigma3(i))
            elseif (nfclo.eq.1) then
c..   factor out LO mt-dep through NLO, expand NNLO,N3LO in 1/mt:
               sallx(3,i) =  amplosq( -2 )*(sigma0(i)
     &              + apimur*sigma1(i) / atauexph( taut, nmt(1,0)))
     &              + c1abs2*(apimur**2*sigma2(i)
     &              + apimur**3*sigma3(i))
            elseif (nfclo.eq.2) then
c..   factor out LO mt-dep through NNLO, expand N3LO in 1/mt:
               sallx(3,i) = amplosq( -2 )*( sigma0(i)
     &              + apimur*sigma1(i) / atauexph( taut,nmt(1,0) )
     &              + apimur**2*sigma2(i) / atauexph( taut,nmt(2,0) ) )
     &              + c1abs2*apimur**3*sigma3(i)
            else
c..   factor out LO mt-dep through N3LO:
               sallx(3,i) = amplosq( -2 )*( sigma0(i)
     &              + apimur*sigma1(i) / atauexph( taut,nmt(1,0) )
     &              + apimur**2*sigma2(i) / atauexph( taut,nmt(2,0) )
     &              + apimur**3*sigma3(i) / atauexph( taut,nmt(3,0) ) )
            endif
            sallx(3,i) = prfc*c1factor(3)*sallx(3,i)
            serrx(3,i) = prfc*c1factor(3)*dsqrt(errsig0(i)**2 +
     &           (apimur*errsig1(i))**2 + (apimur**2*errsig2(i)) +
     &           (apimur**3*errsig3(i)))*amplosq( -2 )
         endif
      enddo
      
      nmtlim=nmtsave
      lball=lballsave
      
      return
      end

c-}}}
C-{{{ subroutine initggh:

      subroutine initggh()
c..
c..   Initialize some parameters.
c..   
      implicit real*8(a-h,o-z)
      complex*16 atau
      include '../commons/common-expand.f'
      include '../commons/common-readdata.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-consts.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'

      include '../commons/common-ren.f'

      pi = 3.14159265358979323846264338328d0
      z2 = 1.6449340668482264364724d0
      z3 = 1.2020569031595942853997d0
      z4 = 1.0823232337111381915160d0
      z5 = 1.0369277551433699263314d0
      z6 = 1.0173430619844491397145d0
      ca = 3.d0
      cf = 1.33333333333333333333d0
      tr = 0.5d0
      gev2pb = .38937966d+9   !  conversion  1/GeV^2 -> pb
      ln2 = dlog(2.d0)

c--   initialize errors and warnings:
      error = 0
      errnorm = 0
      do i=1,100
         warnings(i) = 0
      enddo

c--   if you modify the program, set this to 1 at the relevant places
      modified = 0

c--   cms energy:
      sqrts = sqscmsrd

c--   Higgs mass:
      mh = mhiggsrd

c--   renormalization/factorization scale:
      rmur = rmurmh * mh
      rmuf = rmufmh * mh

c--   mt in MSbar scheme?
      lmtms = .false.

c--   top mass:
      mt = mtoprd

c--   bottom mass:
      mbmb = mbmbrd

c--   number of massless flavors:
      nf = 5.d0

c--   Vegas parameters:
c$$$      acc = 1.d-8       !  accuracy for vegas
c--   sufficient for permille accuracy:
c$$$      itmx1=5         !  interations for vegas run 1
c$$$      ncall1=2000      !  calls for vegas run 1
c$$$      itmx2=2         !  interations for vegas run 2
c$$$      ncall2=5000     !  calls for vegas run 2
c--   higher accuracy:
C$$$      itmx1=10         !  interations for vegas run 1
C$$$      ncall1=5000      !  calls for vegas run 1
C$$$      itmx2=5         !  interations for vegas run 2
C$$$      ncall2=20000     !  calls for vegas run 2
c$$$      nprnv=0          !  =1 -- verbose mode
c$$$      lveg1 = .true.  !  if .false., run vegas only once!

c--   constants and normalization:
      gfermi = gfermird
      mz = mzrd

      tauh = mh**2/sqrts**2

c--   logarithms:
      lfh = dlog(rmuf**2/mh**2)
      lft = dlog(rmuf**2/mt**2)
      lfr = dlog(rmuf**2/rmur**2)
      lrt = lft - lfr
      lth = dlog(mt**2/mh**2)
      lrh = dlog(rmur**2/mh**2)

      ataut2 = abs(atau( 4*mt**2/mh**2 ))**2

c--   electro-weak correction factor:
      if (lpseudo) then
         elwfac = 1.d0
      else
         elwfac = 1.d0 + glgl_elw(mt,mh)
      endif

c--   SM or non-SM?  (0: SM  --  1: other model)
      !default choice lstdmodel = false, also for SM
      if (nmodel.eq.0) then
         lstdmodel = .true.
      else
         lstdmodel = .false.
      endif

c--   SM coefficient functions for effective vertex:
      if (.not.lpseudo) then
         c1sm0 = 1.d0
         c1sm1 = 11/4.d0
         c1sm2 = 2777/288.d0 + 19/16.d0*lrt + nf*( -67/96.d0 +lrt/3.d0)
         c1sm3 = (897943/9216.d0*z3 + 209/64.d0*lrt**2 + 1733/288.d0*lrt
     &        - 2892659/41472.d0 + nf*(-110779/13824.d0*z3 + 23/32.d0
     &        *lrt**2 + 55/54.d0*lrt + 40291/20736.d0) + nf**2*( -1
     &        /18.d0*lrt**2 + 77/1728.d0*lrt - 6865/31104.d0))
         c2sm1 = 0.d0
      else
         c1sm0 = 1.d0
         c1sm1 = 0.d0
         c1sm2 = 0.d0
         c2sm1 = -1/2.d0 + lrt
      endif

c--   initially, set coefficient functions to SM values:
      comc1eff0 = c1sm0*(1.d0,0.d0)
      c1eff1 = c1sm1
      c1eff2 = c1sm2
      c1eff3 = c1sm3
      c2eff1 = c2sm1

c      call rluxgo(3,12348271,0d0,0d0)    !  runlux initialization (optional)
c-- From this line: Obtain effective vertices defined in init.f
      if (.not.lstdmodel) then
            
         comc1eff0 = comc1eff0rd
         if (norder.ge.1) then
            c1eff1 = c1eff1rd
         endif
         if (norder.ge.2) then
            c1eff2 = c1eff2rd
            if (lpseudo) c2eff1 = c2eff1rd
         endif
         if (norder.ge.3) then
            c1eff3 = c1eff3rd
         end if
         
      endif

      return
      end

C-}}}
C-{{{ subroutine normalization:

      subroutine normalizationggh()

      implicit real*8 (a-h,o-z)
      include '../commons/common-consts.f'
      include '../commons/common-keys.f'
      include '../commons/common-vars.f'
      include '../commons/common-sigma.f'
      include '../commons/common-errors.f'

      errnorm = 0

c..   note: these factors are correct! (Also the apimuR/3.d0!)
      if (.not.lpseudo) then
         rho0 = 1/64.d0 * pi*sqrt(2.d0)*gfermi*gev2pb
         prefac = rho0*(apimuR/3.d0)**2
      else
         rho0 = 1/16.d0 * pi*sqrt(2.d0)*gfermi*gev2pb
         prefac = rho0*(apimuR/4.d0)**2
      endif

      end

C-}}}
C-{{{ subroutine printnotice:

      subroutine printnoticeggh(unit)

      integer unit

      write(unit,2011) '# -----------------------'
      write(unit,2011) '# IMPORTANT NOTE'
      write(unit,2011) '# -----------------------'
      write(unit,2011) 
     &     '# If you use ggh@nnlo (or parts of it, or a modified '//
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
     &     '# If you re-distribute ggh@nnlo (or parts of it, '//
     &     'or a modified'
      write(unit,2011) 
     &     '# version of it), or distribute code that links ggh@nnlo '//
     &     'or is based'
      write(unit,2011) 
     &     '# on results of ggh@nnlo (e.g. interpolations), '//
     &     'you are required to:'
      write(unit,2011)
     &     '# * inform the author of ggh@nnlo '//
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
     &     'ggh@nnlo.'
      write(unit,2011) 
     &     '# '//
     &     ''

 2011 format(A)

      end

C-}}}
C-{{{ subroutine printdie:

      subroutine printdieggh(strng)

      character*(*) strng

      write(6,*) 'ggh@nnlo (fatal): ',strng,' Stopped.'
      stop

      end

C-}}}
C-{{{ subroutine printwarn:

      subroutine printwarnggh(strng)

      character*(*) strng

      write(6,*) 'ggh@nnlo (warning): ',strng

      end

C-}}}
C-{{{ subroutine printinfo:

      subroutine printinfoggh(strng)

      character*(*) strng

      write(6,*) 'ggh@nnlo (info): ',strng

      end

C-}}}
