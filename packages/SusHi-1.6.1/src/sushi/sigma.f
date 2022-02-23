! This file is part of SusHi.
! 
! It includes the cross section calculations for gluon fusion
! at NLO and the routine calling ggh@nnlo and bbh@nnlo.
!
C-{{{ calculate cross section

      subroutine XSFull(sigggh,errggh,gg,qg,qq,errorgg,errorqg,
     &errorqq,sigma,errorges,elw,sigweightelw,errorweightelw,
     &xlrvar,XSvar)
      !This routine calculates the full NLO corrections
      !for gluon fusion in the SM and MSSM.
      implicit none
      !input from ggh@nnlo
      double precision sigggh(4,4,-1:5),errggh(4,4,-1:5)
      !output
      double precision gg,qg,qq,elw,elwp1,errorgg,errorqg,errorqq
      double precision sigweightelw,errorweightelw
      double precision sigma(2),errorges(2) !1=LO, 2=NLO
      !internal
      double precision F02,error,chi2,
     &delta,reell,virt,errord,errorr,errorv,prefac,sqrtz,Mtop,
     &sigmanlo,sigmannlo,errornlo,errornnlo,mpggh
      double precision AMPLO,glgl_elw,intgg,intqg,intqq,PDFgg
      double precision xlrvar(2),XSvar(2)
      integer hoint

      external deltagg

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'
      include '../commons/common-quark.f'

      Mtop = dsqrt(mt2)

      F02 = AMPLO(mh2)

      gg = 0.d0
      qg = 0.d0
      qq = 0.d0
      errorgg = 0.d0
      errorqg = 0.d0
      errorqq = 0.d0

      elw = 0.d0

      sigma=0.d0
      errorges=0.d0

      call SUinitpdf(1,SUpdfs(1),SUiset(1))
      call runalpha(apimz,mz,murggh,nf,1,0,apimur)
      alphaggh=apimur*pi
      apiggh(1) = apimur
      if (dist.eq.0) then
         if(subtr) then
            call integrate(1,deltagg,delta,error,chi2)
            sigma(1) = delta*sigmanull*alphaggh**2*F02
            errorges(1) = error*sigmanull*alphaggh**2*F02
         else
            sigma(1) = 0.d0
            errorges(1) = 0.d0
         endif
      elseif(dist.eq.2) then
         if(pseudorap) then
            write(*,105)
            write(*,*)
     &           'Warning: no pseudorapidity-distribution at LO, '
            write(*,*)
     &           'calculating rapidity-distribution instead. '
            write(*,105)
         endif
         sqrtz = dsqrt(z)

         if(subtr) then
            delta = PDFgg(sqrtz*dexp(y),sqrtz*dexp(-y))
         else
            delta = 0.d0
         endif

         sigma(1) = delta*sigmanull*alphaggh**2*F02
      endif
      mpggh = sigma(1)

      if (norderggh.ge.1) then
         call SUinitpdf(2,SUpdfs(2),SUiset(2))
         call runalpha(apimz,mz,murggh,nf,2,0,apimur)
         alphaggh=apimur*pi
         apiggh(2) = apimur

         gg = intgg(delta,reell,virt,errord,errorr,errorv)
         qg = intqg(errorqg)
         qq = intqq(errorqq)

         prefac = alphaggh**3 / Pi * sigmanull * F02

         gg = ( delta * Pi / alphaggh + gg ) * prefac
         qg = qg * prefac
         qq = qq * prefac

         errorgg = dsqrt(errorr**2+errorv**2+(Pi/alphaggh*errord)**2) 
     &        * prefac
         errorqg = errorqg * prefac
         errorqq = errorqq * prefac

         sigma(1) = delta*alphaggh**2*sigmanull*F02
         errorges(1) = errord*alphaggh**2*sigmanull*F02

         sigma(2) = gg + qg + qq
         errorges(2) = dsqrt(errorgg**2 + errorqg**2 + errorqq**2)

!     In the following the results are weighted with ggh@nnlo and/or
!     electroweak corrections, giving sigweightelw and errorweightelw!
!     sigggh and errggh are input from ggh@nnlo!
!     See SusHi manual for reweighting formula.
         if (ew.eq.1) then
            call elw_lq()
         elseif (ew.eq.2) then
            if (pseudo.eq.1) then
               elw = 0.d0
            else
               elw = glgl_elw(Mtop,Mh)
            endif
!glgl_elw provides unreasonable results for mh>1TeV, thus:
            if ((mh.lt.100.d0).or.(mh.gt.1000.d0)) then
               elw = 0.d0
               call printwarnsushi(1,6,'EW corrections only'//
     &              'available for 100 < mh/GeV < 1000. Using 0.')
            endif
         endif
         elwp1 = 1.d0+elw

         call ewreweight(norderggh,ew,sigma(2),errorges(2),sigggh(2,2
     &        ,0),errggh(2,2,0),sigggh(norderggh+1,norderggh+1,0),
     &        errggh(norderggh+1,norderggh+1,0),elwp1,sigweightelw
     &        ,errorweightelw)
         mpggh = sigweightelw

      endif

      !setting the highest output order
      hoint = norderggh
      if ((dist.gt.0).or.(ptcut).or.(rapcut)) hoint = 1

      if (muranalytic) then
      call murvariation(hoint,sigggh,errggh,sigma,errorges,mpggh
     &,elwp1,xlrvar,XSvar)
      end if

 105  format('#--------------------------------------------------#')

      end

C-}}}
C-{{{ subroutine murvariation:

      subroutine murvariation(hoint,sigggh,errggh,sigma,errorges,mpggh
     & ,elwp1,xlrvar,XSvar)
      !analytic determination of the muR dependence of the ggh XS
      implicit none
      !input from ggh@nnlo
      double precision sigggh(4,4,-1:5),errggh(4,4,-1:5)
      integer i,j,k,nstep,hoint !highest order of output
      !The introduction of hoint is necessary, since the highest
      !order output does not always match norderggh, e.g. norderggh
      !can be 3, but in case of a pT cut only the NLO prediction is given.
      double precision mpggh !most precise ggh prediction for output
      double precision sigma(2),errorges(2)
      double precision sigr(0:3),tmpin(0:3),tmpout0(0:3),tmpout1(0:3)
     &    ,tmpout2(0:3),tmp1,tmp2,tmpnlo,rbeg,rend,xlr
      double precision alphaspdfm,apitmp,sigmalo,elwp1
      double precision xlrvar(2),XSvar(2)

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'
      include '../commons/common-quark.f'

      rbeg = facmurgghmin
      rend = facmurgghmax
      nstep = facmurgghstep
      tmp1=0.d0
      tmp2=0.d0
      tmpout0=0.d0
      tmpout1=0.d0
      tmpout2=0.d0
      sigr=0.d0

      XSvar(1) = mpggh !initialization of min and max values
      XSvar(2) = mpggh
      xlrvar(1) = dlog(murggh)
      xlrvar(2) = dlog(murggh)
      write(11,103) ' mh = ',mh,'GeV'
      write(11,104)
     &     ' muR/GeV - sigma(LO)/pb - sigma(NLO)/pb'//
     &     '- sigma(NNLO)/pb - sigma(N3LO)/pb'

c..   
      do i=0,nstep+2
         xlr = dlog(murggh) + dlog(rbeg)*(1-i/(1.d0*nstep)) + i
     &        /(1.d0*nstep)*dlog(rend)
         if (i.eq.(nstep+1)) then
         xlr = dlog(facmurint1*murggh) !include exactly muR/2
         else if (i.eq.(nstep+2)) then
         xlr = dlog(facmurint2*murggh) !include exactly 2*muR
         endif

c..   At LO, tmpin is already defined:
         call rensca_as(apiggh(1),murggh,dexp(xlr),nf,2,0,sigma(1)
     &        ,tmpout0)
         sigr(0)=tmpout0(0)

c..   Now NLO:
         if (hoint.ge.1) then
c..   first mur-dep of exact NLO: tmpout1(1)
            tmpin=0.d0
            tmpin(0) = sigma(1)
            tmpin(1) = sigma(2)
            call rensca_as(apiggh(2),murggh,dexp(xlr),nf,2,1,tmpin
     &           ,tmpout1)

c..   Now NLO,NNLO,N3LO in the heavy-top limit: tmpout2(1,2,3)
            tmpnlo=0.d0
            do k=1,hoint
               do j=0,k
                  tmpin(j)=sigggh(k+1,j+1,0)
               enddo
               call rensca_as(apiggh(k+1),murggh,dexp(xlr),nf,2,k,tmpin
     &              ,tmpout2)
               if (k.eq.1) tmpnlo=tmpout2(1)
               call ewreweight(hoint,ew,tmpout1(1),errorges(2)
     &              ,tmpnlo,errggh(2,2,0),tmpout2(k),errggh(k+1,k+1,0)
     &              ,elwp1,tmp1,tmp2)
               sigr(k)=tmp1
            enddo

         endif

         if (i.le.nstep) then
         write(11,*) dexp(xlr),sigr(0),sigr(1),sigr(2),sigr(3)
         end if

          !obtain boundaries to write in the main output file
          !look at maximum in the interval 0.5<muR/mh<2
          !if the user sets a more narrow band, the obtained
          !result might actually be screwed up.
          if ((facmurint1*murggh-dexp(xlr)).le.(0.0001d0)
     &     .and.((facmurint2*murggh-dexp(xlr)).ge.(-0.0001d0))) then
           if (dexp(xlr).le.dexp(xlrvar(1))) xlrvar(1) = xlr
           if (dexp(xlr).ge.dexp(xlrvar(2))) xlrvar(2) = xlr
           if (sigr(hoint).le.XSvar(1)) XSvar(1) = sigr(hoint)
           if (sigr(hoint).ge.XSvar(2)) XSvar(2) = sigr(hoint)
           !write(*,*) 'mean,min,max',mpggh,XSvar(1),XSvar(2)
          end if

      enddo

      XSvar(1) = XSvar(1)-mpggh
      XSvar(2) = XSvar(2)-mpggh

      !choose the larger deviation as symmetric error
c$$$      if (dabs(XSvar(2)).ge.dabs(XSvar(1))) then
c$$$       XSvar(1) = - XSvar(2)
c$$$      else
c$$$       XSvar(2) = - XSvar(1)
c$$$      end if

 103  format('#',a,0p,e16.8,1x,a)
 104  format('#',a)

      end

C-}}}
C-{{{ subroutine ewreweight:

      subroutine ewreweight(nord,ew,sigma2,errorges2
     &     ,sigmanlo,errornlo,sigmannlo,errornnlo,elwp1,sigweightelw
     &     ,errorweightelw)

      implicit none
      integer nord,ew
      real*8 sigma2,errorges2,sigmanlo,errornlo,sigmannlo,errornnlo
     &     ,elwp1,sigweightelw,errorweightelw

      if (nord.eq.1) then
         if (ew.eq.1) then
            sigweightelw = sigma2*elwp1
            errorweightelw = errorges2*elwp1
         elseif (ew.eq.2) then
            sigweightelw = sigma2 + sigmanlo*(elwp1 - 1)
            errorweightelw = dsqrt(errorges2**2
     &           +(errornlo*(elwp1 - 1))**2)
         else
            sigweightelw = sigma2
            errorweightelw = errorges2
         endif
      else
         if (ew.eq.1) then
            sigweightelw = sigma2*elwp1
     &           + sigmannlo - sigmanlo
            errorweightelw = dsqrt((errorges2*elwp1)**2
     &           + errornnlo**2 + errornlo**2)
         elseif (ew.eq.2) then
            sigweightelw = sigma2 + sigmannlo*elwp1 - sigmanlo
            errorweightelw = dsqrt(errorges2**2
     &           +(errornnlo*elwp1)**2+errornlo**2)
         else
            sigweightelw = sigma2 + sigmannlo - sigmanlo
            errorweightelw = dsqrt(errorges2**2
     &           +errornnlo**2+errornlo**2)
         endif
      endif
      return
      end

C-}}}
C-{{{ calculate cross section ggh + bbh at various orders

      subroutine XSgghbbh(ordercall,sigggh,errggh,sigbbh,errbbh)
      !This routine calls ggh@nnlo and bbh@nnlo at order ordercall.
      implicit none
      integer i,j,ppcollin
      !input
      integer ordercall, pseudosc
      !output
      double precision sigggh(4,-1:5),errggh(4,-1:5),sigbbh,errbbh
!     internal
      double precision sigdum(0:3,-1:5),errdum(0:3,-1:5),Mtop
      logical ppcoll
      common /coll/ ppcoll

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'
      include '../commons/common-quark.f'
      include '../commons/common-inputoutput.f'
      include '../commons/common-keys.f'

      if (ppcoll) then
       ppcollin = 0
      else
       ppcollin = 1
      end if

      sigggh=0.d0
      errggh=0.d0
      sigdum=0.d0
      errdum=0.d0
      sigbbh = 0.D0

      if (pseudo.eq.1) then
       pseudosc = 1
      else
       pseudosc = 0
      end if

      Mtop = dsqrt(mt2)

      if (ordercall.eq.0) then
         call SUinitpdf(1,SUpdfs(1),SUiset(1))
         if (norderbbh.ge.0) then
!obtain a LO prediction for bb->H from bbh@nnlo using LO PDF set
            call bbhcall(0,ppcollin,cme,mh,murfacbbh,muffacbbh,
     &           mbmb,gf,Mz,sigbbh,errbbh)
         end if
      end if

      if (ordercall.eq.1) then
         call SUinitpdf(2,SUpdfs(2),SUiset(2))
         if ((gt.ne.0.d0).or.((gt1.ne.0.d0).and.nnlostop).or.ldim5) then
            !obtain a NLO prediction for gg->H from ggh@nnlo using NLO PDF set
            !The result is only needed if
            !1. the effective Higgs-top-quark gt coupling is not vanishing.
            !2. stop contributions are added to the Wilson coefficients, thus
            !   we check if the stop1-H-coupling gt1 is nonzero and
            !   stop contributions are activated within ggh@nnlo
            !3. dim5 operators are included
            call gghcall(1,cme,Mh,muRfacggh,muFfacggh,
     &           1,pseudosc,ppcollin,GF,Mz,mbmb,mtop,sigdum,errdum)
         endif
         if (norderbbh.ge.1) then
!     obtain a NLO prediction for bb->H from bbh@nnlo using NLO PDF set
            call bbhcall(1,ppcollin,cme,Mh,murfacbbh,muffacbbh,mbmb
     &           ,GF,Mz,sigbbh,errbbh)
         end if
      end if

      if (ordercall.eq.2) then
         call SUinitpdf(3,SUpdfs(3),SUiset(3))
         call runalpha(apimz,mz,murggh,nf,3,0,apimur)
         apiggh(3) = apimur
         if ((gt.ne.0.d0).or.((gt1.ne.0.d0).and.nnlostop).or.ldim5) then
            !obtain a NNLO prediction for gg->H from ggh@nnlo using NNLO PDF set
            !Check the NLO case for the if-statements
            call gghcall(1,cme,Mh,muRfacggh,muFfacggh,
     &           2,pseudosc,ppcollin,GF,Mz,mbmb,mtop,sigdum,errdum)
         endif
         if(norderbbh.ge.2) then
            !obtain a NNLO prediction for bb->H from bbh@nnlo using NNLO PDF set
            call bbhcall(norderbbh,ppcollin,cme,Mh,murfacbbh
     &           ,muffacbbh,mbmb,GF,Mz,sigbbh,errbbh)
         end if
      end if


!     N3LO
      if (ordercall.eq.3) then
c..   still using NNLO pdfs, as long as N3LO is not available:
         call SUinitpdf(4,SUpdfs(4),SUiset(4))
         call runalpha(apimz,mz,murggh,nf,runasn3lo,0,apimur)
         apiggh(4) = apimur
         if ((gt.ne.0.d0).or.((gt1.ne.0.d0).and.nnlostop).or.ldim5) then
            call gghcall(1,cme,Mh,muRfacggh,muFfacggh,
     &           3,pseudosc,0,GF,Mz,mbmb,mtop,sigdum,errdum)
         endif
      end if

      do j=1,ordercall+1
         do i=-1,5
            sigggh(j,i) = sigdum(j-1,i)
            errggh(j,i) = errdum(j-1,i)
         enddo
      enddo

      end
C-}}}
C-{{{ gg-channel

      function intgg(delta,reell,virt,errord,errorr,errorv)
      implicit none
      double precision intgg,error(3),chi2(3),delta,virtual,
     &hard,sub,errord,errorr,errorv,reell,virt
      double precision c_virt,dtermsz,PDFgg

      external deltagg,intsubgg,realgg,realggpt,realggy,realggpty

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'

      delta = 0.d0
      sub = 0.d0
      virtual = 0.d0
      error(1) = 0.d0
      error(3) = 0.d0

      if(dist.eq.0) then

         if(subtr) then

c     virtual contribution

            virtual = c_virt()
            call integrate(1,deltagg,delta,error(1),chi2(1))

c     subtraction term + collinear counterterm

            call integrate(2,intsubgg,sub,error(3),chi2(3))
         endif

c     real contribution

         call integrate(3,realgg,hard,error(2),chi2(2))

      elseif(dist.eq.1) then

c     pt - distribution

         if(maxpt.lt.30.d0) then
            write(*,105)
            write(*,*) '# warning: pt < 30 GeV. '
            write(*,105)
         endif
         call integrate(2,realggpt,hard,error(2),chi2(2))

         if(scalespt) then
            hard = hard / alphaggh**3
            error(2) = error(2) / alphaggh**3
         endif

      elseif(dist.eq.2) then

c     y - distribution

         call integrate(2,realggy,hard,error(2),chi2(2))

      elseif(dist.eq.3) then

c     pt - y - distribution

         if(maxpt.lt.30.d0) then
            write(*,105)
            write(*,*) '# warning: pt < 30 GeV. '
            write(*,105)
         endif
         call integrate(1,realggpty,hard,error(2),chi2(2))
         if(scalespt) then
            hard = hard / alphaggh**3
            error(2) = error(2) / alphaggh**3
         endif

      endif

      errord = error(1)

      reell = ca*(hard+sub)
      errorr = ca*dsqrt(error(2)**2+error(3)**2)

      if(subtr) then
         virt = delta*(virtual+Pi**2-2*betanull*lfr+ca*dtermsz())
         errorv = error(1)*(virtual+Pi**2-2*betanull*lfr+ca*dtermsz())
      else
         virt = 0.d0
         errorv = 0.d0
      endif

      intgg = reell + virt

 105  format('#--------------------------------------------------#')

      end

C-}}}
C-{{{ qg-channel

      function intqg(error1)
      implicit none
      double precision intqg,hard,sub,error(2),chi2(2),error1

      external realqg,intsubqg,realqgpt,realqgy,realqgpty

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'

      sub = 0.d0
      error(2) = 0.d0

      if(dist.eq.0) then

c     real contribution

         call integrate(3,realqg,hard,error(1),chi2(1))

c     subtraction term + collinear counterterm

         if(subtr) then
            call integrate(2,intsubqg,sub,error(2),chi2(2))
         endif

      elseif(dist.eq.1) then

c     pt - distribution

         call integrate(2,realqgpt,hard,error(1),chi2(1))
         if(scalespt) then
            hard = hard / alphaggh**3
            error(1) = error(1) / alphaggh**3
         endif
      elseif(dist.eq.2) then

c     y - distribution

         call integrate(2,realqgy,hard,error(1),chi2(1))

      elseif(dist.eq.3) then

c     pt - y - distribution

         call integrate(1,realqgpty,hard,error(1),chi2(1))
         if(scalespt) then
            hard = hard / alphaggh**3
            error(1) = error(1) / alphaggh**3
         endif
      endif

      intqg=(hard+sub)*cf/2.d0
      error1 = dsqrt(error(1)**2+error(2)**2)*cf/2.d0
      end

C-}}}
C-{{{ qq-channel

      function intqq(error)
      implicit none
      double precision intqq,hard,error,chi2

      external realqq,realqqpt,realqqy,realqqpty

      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      include '../commons/common-int.f'

      if(dist.eq.0) then

         call integrate(3,realqq,hard,error,chi2)

      elseif(dist.eq.1) then

c     pt - distribution

         call integrate(2,realqqpt,hard,error,chi2)
         if(scalespt) then
            hard = hard / alphaggh**3
            error = error / alphaggh**3
         endif

      elseif(dist.eq.2) then

c     y - distribution

         call integrate(2,realqqy,hard,error,chi2)

      elseif(dist.eq.3) then

c     pt - y - distribution

         call integrate(1,realqqpty,hard,error,chi2)
         if(scalespt) then
            hard = hard / alphaggh**3
            error = error / alphaggh**3
         endif
      endif

      intqq = hard*32/27.d0
      error = error*32/27.d0
      end

C-}}}
