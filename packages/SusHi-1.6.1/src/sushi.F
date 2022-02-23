! This file is main file of SusHi.
! 
! written by Robert Harlander, Stefan Liebler, Hendrik Mantler
! Bergische Universitaet Wuppertal
! Version 1.6.1
! 
      program SUSHI

      implicit none
      integer count1,count2,ninp,noutp,i,j
      double precision gghdum(4,-1:5),ggherrdum(4,-1:5)
      double precision SigLO,errLO,Sig_bb,err_bb,Sig_bg,err_bg
      double precision alphasPDFm
      double precision wHex(2)
      logical ppcoll
      common /coll/ ppcoll

      !common variables used by various routines
      include 'commons/common-inputoutput.f'
      include 'commons/common-consts.f'
      include 'commons/common-vars.f'
      include 'commons/common-keys.f'
      include 'commons/common-quark.f'
      include 'commons/common-int.f'
      include 'commons/common-ren.f'
      include 'commons/common-citations.f'
      include 'commons/common-flags.f'

      common /vgasio/ ninp,noutp

      call get_command_argument(1,jfilein)
      call get_command_argument(2,jfnameact)

      if ((len_trim(jfilein).eq.0).or.(len_trim(jfnameact).eq.0)) then
       write(*,*) "SusHi must be called with two arguments:"
       write(*,*) "./SusHi in.in out.out"
       write(*,*) "The input-file in.in has to be provided"
       write(*,*) "in the SLHA-standard. Example files: SM.in, ..."
       stop
      end if

      extflag = .false.

      call initsushi(htlflag)
      
      if (verbose.ne.0) then
         write(*,105)
         write(*,106) 'Higgs mass:',Mh
         write(*,105)
      endif

      !VEGAS output
      noutp = 6

      
C-{{{ cuts

      !definition of cuts
      subtr = .true.
      ptcut = .false.
      rapcut = .false.

      minpt = 0.d0
      maxpt = cme*(1-z)/2.d0
      minrap = 0.d0
      maxrap = -dlog(z)/2.d0

      pt = maxptc
      y = maxrapc
c..   ptref is only introduced to check whether we are in the
c..   low-pt region or not (if yes, a warning is printed)
      ptref=pt

c     juser=.true. => apply cuts that the user can define in the function 'jetuser'
c     in src/sushi/integration.f
      if(juser) then
         call printinfosushi(6,
     &        'Warning: Applying cuts defined in "juser" in'
     &        //' src/sushi/integration.f.')
         call printinfosushi(6,
     &        'Hence, ignoring all cuts defined in the Block DISTRIB.')
         ptc = 0
         rapc = 0
      endif

      if(dist.ge.2) then
c..   dsigma/dy or d2sigma/(dpt dy)
         if((maxrapc.lt.0.d0).or.(
     &     (maxrapc.gt.maxrap).and.(.not.pseudorap) )) then
            write(*,105)
            write(*,*) 'Error: y out of range. '
            write(*,*) 'Abs[y] < ', maxrap
            write(*,105)
            stop
         endif         
      endif

      if((dist.eq.1).or.(dist.eq.3)) then
c..   dsigma/dpt or d2sigma/(dpt dy)
         subtr = .false.
         if ((norderggh.eq.0).or.(norderbbh.eq.0)) then
            write(6,*) 'SUSHI(5)=',norderggh
            write(6,*) 'SUSHI(6)=',norderbbh
            write(6,*) 'DISTRIB(1)=',dist
            call printinfosushi(6,'At LO, the pt of the Higgs '//
     &           'is always zero.')
            call printerrorsushi(1,6
     &           ,'Therefore, no distribution can be calulated. ')
         endif
         if((maxptc.le.0.d0).or.(maxptc.gt.maxpt)) then
            call printerrorsushi(1,6,'Error: pt out of range. ')
         endif
      else
         scalespt = .false.
      endif
      if(dist.eq.3) then
c..   d2sigma/(dpt dy)
         if(maxptc.gt.dsqrt(mh2
     &        *(dexp(2.d0*maxrapc)*(1.d0/z+z)-1.d0-dexp(4.d0*maxrapc)))/
     &        (1.d0+dexp(2.d0*maxrapc))) then
            write(*,105)
            write(*,*) 'Error: pt out of range. '
            write(*,*) 'for y = ', maxrapc, ', pt must satisfy'
            write(*,*) '0 < pt < ', dsqrt(mh2*(dexp(2.d0*maxrapc)
     &*(1.d0/z+z)-1.d0-dexp(4.d0*maxrapc)))/(1.d0+dexp(2.d0*maxrapc))
            write(*,105)
            stop
         endif
      endif

      if(ptc.eq.0) then
      !no pt-cut
         minptc = minpt
         maxptc = maxpt
      elseif(ptc.eq.1) then
      !pt > ptmin
         ptref=minptc
         if((minptc.lt.0.d0).or.(minptc.ge.maxpt)) then
            write(*,105)
            write(*,*) 'Error: ptmin out of range. '
            write(*,105)
            stop
         endif
         if(minptc.gt.0.d0) then
            subtr = .false.
            ptcut = .true.
            minpt = minptc
         endif
         maxptc = maxpt
      elseif(ptc.eq.2) then
         ptref=maxptc
      !pt < ptmax
         if(maxptc.lt.0.d0) then
            write(*,105)
            write(*,*) 'Error: ptmax out of range. '
            write(*,105)
            stop
         endif
         ptcut = .true.
         if(maxptc.gt.maxpt) then
            write(*,105)
            write(*,*) 'Warning: ptmax too large. '
            write(*,*) 'maximal value for ptmax is ', maxpt, ' GeV. '
            write(*,105)
            maxptc = maxpt
         endif
         minptc = minpt
      elseif(ptc.eq.3) then
      !ptmin < pt < ptmax
         ptref=minptc
         if((minptc.lt.0.d0).or.(minptc.ge.maxpt)) then
            write(*,105)
            write(*,*) 'Error: ptmin out of range. '
            write(*,105)
            stop
         endif
         if(maxptc.lt.0.d0) then
            write(*,105)
            write(*,*) 'Error: ptmax out of range. '
            write(*,105)
            stop
         endif
         if(maxptc.le.minptc) then
            write(*,105) 
            write(*,*) 'Error: ptmin larger then ptmax. '
            write(*,105) 
            stop
         endif
         ptcut = .true.
         if(minptc.gt.0.d0) then
            subtr = .false.
            minpt = minptc
         endif
         if(maxptc.gt.maxpt) then
            write(*,105)
            write(*,*) 'Warning: ptmax too large. '
            write(*,*) 'maximal value for ptmax is ', maxpt, ' GeV. '
            write(*,105)
            maxptc = maxpt
         endif
      else
         write(*,105)
         write(*,*) 'Error: bad value in Block DISTRIB. '
         write(*,105)
         stop
      endif

      if(rapc.eq.0) then
      !no ycut
      elseif(rapc.eq.1) then
      !Abs[y] < ymax
         if(maxrapc.lt.0.d0) then
            write(*,105)
            write(*,*) 'Error: ymax out of range. '
            write(*,105)
            stop
         endif            
         rapcut = .true.
         if((maxrapc.gt.maxrap).and.(.not.pseudorap)) then
            write(*,105)
            write(*,*) 'Warning: ymax too large. '
            write(*,*) 'maximal value for ymax is ', maxrap, '. '
            write(*,105)
         else
            maxrap = maxrapc
         endif
      elseif(rapc.eq.2) then
      !Abs[y] > ymin
         if(pseudorap) then
            if(minrapc.lt.0.d0) then
               write(*,105)
               write(*,*) 'Error: ymin out of range. '
               write(*,*) '0 < ymin'
               write(*,105)
               stop
            endif
            maxrap = 1.d100
         else
            if((minrapc.lt.0.d0).or.(minrapc.ge.maxrap)) then
               write(*,105)
               write(*,*) 'Error: ymin out of range. '
               write(*,*) '0 < ymin < ', maxrap
               write(*,105)
               stop
            endif
         endif
         rapcut = .true.
         minrap = minrapc
      elseif(rapc.eq.3) then
      !ymin < Abs[y] < ymax
         if(pseudorap) then
            if(minrapc.lt.0.d0) then
               write(*,105)
               write(*,*) 'Error: ymin out of range. '
               write(*,*) '0 < ymin'
               write(*,105)
               stop
            endif
            maxrap = maxrapc
         else
            if((minrapc.lt.0.d0).or.(minrapc.ge.maxrap)) then
               write(*,105)
               write(*,*) 'Error: ymin out of range. '
               write(*,*) '0 < ymin < ', maxrap
               write(*,105)
               stop
            endif
            if(maxrapc.gt.maxrap) then
               write(*,105)
               write(*,*) 'Warning: ymax too large. '
               write(*,*) 'maximal value for ymax is ', maxrap, '. '
               write(*,105)
            else
               maxrap = maxrapc
            endif
         endif
         if(maxrapc.lt.0.d0) then
            write(*,105)
            write(*,*) 'Error: ymax out of range. '
            write(*,105)
            stop
         endif            
         rapcut = .true.
         minrap = minrapc
      else
         write(*,105)
         write(*,*) 'Error: bad value in Block DISTRIB. '
         write(*,105)
         stop         
      endif

      if(rapcut.and.pseudorap.and.subtr) then
         write(*,105)
         write(*,*) 'Error: Pseudorapidity cuts cannot be ',
     &'applied when pt = 0 is allowed. '
         write(*,105)
         stop
      endif

      if((dist.eq.2).and.subtr) then
         write(*,105)
         write(*,*) 'Error: (Pseudo)rapidity distributions cannot be ',
     &'calculated when pt = 0 is allowed. '
         write(*,105)
         stop
      endif

      if(((dist.eq.1).or.(dist.eq.3)).and.ptcut) then
         write(*,105)
         write(*,*) 'Warning: Cuts in pt will be ignored ',
     &'when calculating pt-distribution. '
         write(*,105)
         ptcut = .false.
         minpt = 0.d0
c         maxpt = cme*(1-z)/2.d0
      endif

      if((dist.ge.2).and.rapcut) then
         write(*,105)
         write(*,*) 'Warning: (Pseudo)rapidity cuts will be ignored ',
     &'when calculating (pseudo)rapidity-distribution. '
         write(*,105)
         rapcut = .false.
         minrap = 0.d0
         maxrap = -dlog(z)/2.d0
      endif

c     regarding phase space integration and convolution with PDFs:
c     integration over pT, y and x1 (pty=.true.) or t_hat, x1 and x2 (pty=.false.)
      pty = .true.

C-}}}
C-{{{ calculation of cross sections - calling ggh@nnlo and bbh@nnlo

      subbh=0.d0
      subbherr=0.d0
      suggh=0.d0
      suggherr=0.d0

c$$$      if((norderggh.eq.0).or.((dist.eq.0).and.(norderbbh.ge.0))) then
c$$$         call SUinitpdf(1,SUpdfs(1),SUiset(1))
c$$$         !for bbh alpha_s(muRfacbbh*mh) is calculated internally in bbh@nnlo
c$$$      endif
c$$$      if((norderggh.ge.1).or.((dist.eq.0).and.(norderbbh.ge.1))) then
c$$$         call SUinitpdf(2,SUpdfs(2),SUiset(2))
c$$$      endif
c$$$      if(ltotinc.and.((norderggh.ge.2).or.(norderbbh.ge.2))) then
c$$$         call SUinitpdf(3,SUpdfs(3),SUiset(3))
c$$$      endif

c$$$            print*,dist,ptc,rapc,norderbbh,lqqhkey
c$$$            stop

      if (ptcut) then
         if(pseudorap) then
            write(*,*) "Applied p_t cuts:         ",minptc,maxptc
         else
            write(*,*) "Applied p_t cuts:    ",minptc,maxptc
         endif
      end if
      if (rapcut) then
         if(pseudorap) then
            write(*,*) "Applied pseudorapidity cuts:",minrap,maxrap
         else
            write(*,*) "Applied rapidity cuts:",minrap,maxrap
         endif
      end if

      !Differential distributions for bbh
      if ((dist.eq.0).and.((ptc.gt.0).or.(rapc.gt.0))
     &.or.(dist.eq.1)) then
         if (norderbbh.ge.0) then
!            if (lqqhkey) then
!               write(6,*) 'SusHi: Differential distributions'/
!     &              /' incompatible with Block QQH.'
!               write(6,*) '       Remove Block QQH from input'/
!     &              /' file, or change DISTRIB(1) to 0.'
!               stop 'SusHi stopped.'
!            endif

            if (verbose.ne.0) then
               write(*,105)
               write(*,*) "Performing (N)LO bbh calculation"
               write(*,105)
            endif

            citations(9) = 1

            call SUinitpdf(1,SUpdfs(1),SUiset(1))
            call runCouplings(murbbh,murbbh,0)
            call bbhdiffinit(0,cme,mh,0.d0,mt
     &           ,-1,-1,0,pseudorap,minptc,maxptc,minrap,maxrap
     &           ,dist,pt,murbbh,mufbbh)
            call Sigma_bbHNLO(SigLO,errLO,Sig_bb,err_bb,
     &           Sig_bg,err_bg,SUbbh(1),SUbbherr(1))

            if (norderbbh.ge.1) then
               call SUinitpdf(2,SUpdfs(2),SUiset(2))
               call runCouplings(murbbh,murbbh,1)
               call bbhdiffinit(1,cme,mh,0.d0,mt
     &           ,-1,-1,0,pseudorap,minptc,maxptc,minrap,maxrap
     &           ,dist,pt,murbbh,mufbbh)
               call Sigma_bbHNLO(SigLO,errLO,Sig_bb,err_bb,
     &              Sig_bg,err_bg,SUbbh(2),SUbbherr(2))
            end if

            if (norderbbh.eq.2) then
               write(*,105)
               write(*,*)
     &"Warning: Differential cross sections for bb->H are"
               write(*,*)
     &"not available at NNLO. NLO results are given."
               write(*,105)
            end if
         end if
      end if

C-}}}
C-{{{ settings for NLO calculation

      call initcouplings()

      if (ltotinc) then

         if (verbose.ne.0) then
            write(*,105)
            write(*,*) "Calling ggh@nnlo and bbh@nnlo"
            write(*,105)
         endif

         do i=0,max(norderggh,norderbbh)
            call XSgghbbh(i,gghdum,ggherrdum,SUbbh(i+1),SUbbherr(i+1))
            do j=1,i+1
               do count1=-1,5
                  SUggh(i+1,j,count1) = gghdum(j,count1)
                  SUggherr(i+1,j,count1) = ggherrdum(j,count1)
               enddo
            enddo
         enddo
      end if

      if (verbose.ne.0) then
      if (norderggh.eq.0) then
         write(*,105)
         write(*,*) "Performing LO ggh calculation"
         write(*,105)
      else if (norderggh.gt.0) then
         write(*,105)
         write(*,*) "Performing NLO ggh calculation"
         write(*,105)
      endif
      endif

C-}}}
C-{{{ weighting of bbh results
      wHex = 1.d0
      if ((dist.eq.0).or.(dist.eq.1)) then
      !weighting bbh result
         if (model.ge.1) then
            call bbhweights(model,twohdmver,Hmx,Amx,pseudo,Sind,tanb,
     &           muSUSY,lam,GF,delmb,wHex)
         end if
         !take into account BLOCK FACTORS for bottom annihilation
         wHex = wHex * yukfac(3) * yukfac(3)
         if (norderbbh.eq.0) then
            do count1=1,3
              SUbbh(count1) = wHex(1) * SUbbh(count1) 
              SUbbherr(count1) = wHex(1) * SUbbherr(count1) 
            end do
         else
            do count1=1,3
              SUbbh(count1) = wHex(2) * SUbbh(count1) 
              SUbbherr(count1) = wHex(2) * SUbbherr(count1) 
            end do
         end if

      end if
C-}}}
C-{{{ calculation of cross sections - NLO full contributions

      lfh = 2.d0*dlog(muffacggh)
      lfr = 2.d0*dlog(muffacggh/murfacggh)

      if (norderggh.ge.0) then
         
         call SUinitpdf(norderggh+1,SUpdfs(norderggh+1),
     &        SUiset(norderggh+1))
         call runalpha(apimz,mz,murggh,nf,norderggh+1,0,apimur)
         alphaggh=apimur*pi

         if (model.eq.2) then
            call renormalizeSM
         endif
         
         if (model.eq.0) then
            !set couplings to SM values:
            gc = yukfac(1)
            gt = yukfac(2)
            gb = yukfac(3)
            gtp = yukfac(4)
            gbp = yukfac(5)
            call renormalizeSM
         end if

         if ((ldim5).and.(dim5run.ne.1)) then
          muranalytic = .false. !analytic mur determination only useful for perturbative dim5 running
         end if

         i = len(jfnameact)
         do while (jfnameact(i:i).eq.' ')
            i = i - 1
         enddo
         if (muranalytic)
     &    open(unit=11,file=jfnameact(1:i)//'_murdep',status='unknown')
         call XSFull(SUggh,SUggherr,ggchan,qgchan,qqchan,ggchanerr,
     &        qgchanerr,qqchanerr,xsec,xsecerr,elw,xeff,xefferr,
     &        xeffvar,gghXSvar)
         if (muranalytic) close(11)

      end if


C-}}}
C-{{{ showing screen output and writing the output-file
      call screenoutputsushi()
      if (scalespt.and.(norderggh.ge.1)) then
         call SUinitpdf(2,SUpdfs(2),SUiset(2))
         awithpt = alphasPDFm(2,muRfacggh*dsqrt(mh**2+pt**2))
         awopt = alphasPDFm(2,murggh)
      endif
      call outputsushi()

      write(*,104) " Results written in output file: ",jfnameact
      write(*,*) achar(27)//'[32mMake sure to cite the references listed
     &'//' in output file '//trim(jfnameact)//'!'//achar(27)//'[0m'
      stop
      call printcopysushi(0,6,verbose)

 104  format(a,a)
 105  format('#--------------------------------------------------#')
 106  format(a12,'    ',f10.5, '   GeV')
C-}}}

      end

!end of file
