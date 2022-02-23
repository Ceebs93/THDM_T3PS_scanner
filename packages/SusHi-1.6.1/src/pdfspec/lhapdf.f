C-{{{ pdfinit
      subroutine SUinitpdf(orderin,pdfnamein,isetin)
      implicit none
      double precision alphasPDFm
      character pdfnamein*50
      integer order,orderin,isetin,i

      include '../commons/common-consts.f'
      include '../commons/common-int.f'
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      order=orderin
      if (order.gt.3) then
         call printwarnsushi(verbose,6
     &        ,'N3LO PDF evolution not available.'//' Using NNLO.')
         order=3
      endif

c..   check if requested PDF set is currently active:

      if ((order.eq.orderc).and.(isetc.eq.isetin).and.
     &     (pdfnamec.eq.pdfnamein)) return

      if (verbose.eq.0) call SetLHAPARM('SILENT')
      call InitPDFsetByNamem(order,pdfnamein)
      call InitPDFm(order,isetin)

      orderc=order
      pdfnamec=pdfnamein
      isetc=isetin

      apimz = alphasPDFm(order,mz)/Pi
c..   We introduce the do-loop here, because otherwise currently
c..   apimzlist(4) would never be set. Now it is set to apimzlist(3):
      do i=order,4
         apimzlist(i)=apimz
      enddo
         
      pdforder = order

      return
      end
C-}}}
C-{{{ GetPDFs

      subroutine GetPDFs(x,up,dn,usea,dsea,
     &str,sbar,chm,cbar,bot,bbar,glu)
      implicit none
      double precision x,up,dn,usea,dsea,
     &str,sbar,chm,cbar,bot,bbar,glu,f(-6:6)

      include '../commons/common-int.f'
      include '../commons/common-vars.f'

      call evolvePDFm(pdforder,x,muFggh,f)
      up = f(2)
      dn = f(1)
      usea = f(-2)
      dsea = f(-1)
      str = f(3)
      sbar = f(-3)
      chm = f(4)
      cbar = f(-4)
      bot = f(5)
      bbar = f(-5)
      glu = f(0)
      end

C-}}}
C-{{{ function pdfsold(xx,rmuf):

      subroutine pdfsold(orderin,xx,rmuf,up,dn,usea,dsea,str,sbar,
     &     chm,cbar,bot,bbar,glu)
      !needed by ggh@nnlo and bbh@nnlo routines
      implicit real*8 (a-z)
      real*8 ff(-6:6)
      integer order,orderin,neigen
      parameter(neigen=20)      ! number of eigenvectors
      double precision up,dn,upv,dnv,usea,dsea,str,chm,bot

      xxl = xx
      rmufl = rmuf

c..   !!! WARNING: this should be changed once N3LO pdfs are available
      order=min(orderin,2)

      call evolvepdfm(order+1,xxl,rmufl,ff)

      up   = ff(1)/xx
      usea = ff(-1)/xx
      dn   = ff(2)/xx
      dsea = ff(-2)/xx
      str  = ff(3)/xx
      sbar = ff(-3)/xx
      chm  = ff(4)/xx
      cbar = ff(-4)/xx
      bot  = ff(5)/xx
      bbar = ff(-5)/xx
      glu  = ff(0)/xx

c$$$      upv = up-usea
c$$$      dnv = dn-dsea

      end

C-}}}
C-{{{ function pdfs(xx,rmuf):

      subroutine pdfs(sign,orderin,xx,rmuf,ffout)
      !needed by ggh@nnlo and bbh@nnlo routines
      implicit none
      integer i
      real*8 ff(-6:6),ffout(-6:6)
      real*8 xx,rmuf
      integer order,orderin,neigen,sign
      parameter(neigen=20)      ! number of eigenvectors

c..   !!! WARNING: this should be changed once N3LO pdfs are available
      order=min(orderin,2)

      call evolvepdfm(order+1,xx,rmuf,ff)

      do i=-6,6
         ffout(i)   = ff(sign*i)/xx
      enddo

      end

C-}}}


