!The following routines are from bbh (ggh below):

C-{{{ function bbqpdf(x1,x2):

      real*8 function bbqpdf(x1,x2)
C
C     b-bbar pdfs
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-lumis.f'

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(ppbarsign,norder,x2,rmuf,ff2)

      bbqpdf = ff1(ptn1)*ff2(ptn2) + ff1(ptn2)*ff2(ptn1)

      return
      end

C-}}}
C-{{{ function bgpdf(x1,x2):

      real*8 function bgpdf(x1,x2)
C
C     bottom-gluon pdfs, or more accurately:
C
C     b(x1)*g(x2) + g(x1)*b(x2) + bbar(x1)*g(x2) + g(x1)*bbar(x2)
C
C     this can be done because  
C     \sigma_{b g} = \sigma_{g b} = \sigma_{bbar g} = \sigma_{g bbar}
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-lumis.f'

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(ppbarsign,norder,x2,rmuf,ff2)

c..   we can use bot+bbar, because sigma(b g)=sigma(bbar g):
      bgpdf = ( ff1(ptn1) + ff1(ptn2) ) * ff2(0) 
     &     + ff1(0) * ( ff2(ptn1) + ff2(ptn2) )

      return
      end

C-}}}
C-{{{ function bbpdf(xx,yy):

      real*8 function bbpdf(x1,x2)
C
C     b-b pdfs, or more accurately:
C
C     b(x1)*b(x2) + bbar(x1)*bbar(x2).
C
C     This can be done because
C
C     \sigma_{q q} = \sigma_{qbar qbar}
C
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-lumis.f'

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(ppbarsign,norder,x2,rmuf,ff2)
C
C     note that sigma(b b) = sigma(bbar bbar):
C
      bbpdf = ff1(ptn1)*ff2(ptn1) + ff1(ptn2)*ff2(ptn2)

      return
      end

C-}}}
C-{{{ function bqpdf(xx,yy):

      real*8 function bqpdf(x1,x2)
C
C     b-q pdfs, or more accurately:
C
C     ( b(x1) + bbar(x1) )*( q(x2) + qbar(x2) )
C        + ( b(x2) + bbar(x2) )*( q(x1) + qbar(x1) )
C
C     This can be done because
C
C     \sigma_{b q} = \sigma_{b qbar} = \sigma_{bbar q} = \sigma_{bbar qbar}
C
C     Note that q and q' denote different flavors!
C     The case where q and q' have the same flavor is treated in
C     <function qqpdf>!
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      real*8 qrk1,qrk2,brk1,brk2
      integer i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'
      include '../commons/common-lumis.f'

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(ppbarsign,norder,x2,rmuf,ff2)

      qrk1=0.d0
      qrk2=0.d0
      brk1=0.d0
      brk2=0.d0

      do i=1,5
         qrk1 = qrk1 + ff1(i)
         qrk2 = qrk2 + ff2(i)
         brk1 = brk1 + ff1(-i)
         brk2 = brk2 + ff2(-i)
      enddo

      qrk1 = qrk1 - ff1(ptn1)
      qrk2 = qrk2 - ff2(ptn1)
      brk1 = brk1 - ff1(ptn2)
      brk2 = brk2 - ff2(ptn2)

      bqpdf = (ff1(ptn1) + ff1(ptn2)) * (qrk2 + brk2) 
     &     + (ff2(ptn1) + ff2(ptn2)) * (qrk1 + brk1)

      return
      end

C-}}}

!The following routines are from ggh:

C-{{{ function ggpdf(x1,x2):

      real*8 function ggpdf(x1,x2)
C
C     gluon-gluon pdfs
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      ggpdf = ff1(0)*ff2(0)

      return
      end

C-}}}
C-{{{ function qgpdf(x1,x2):

      real*8 function qgpdf(x1,x2)
C
C     quark-gluon pdfs, or more accurately:
C
C     q(x1)*g(x2) + g(x1)*q(x2) + qbar(x1)*g(x2) + g(x1)*qbar(x2)
C
C     this can be done because  
C     \sigma_{q g} = \sigma_{g q} = \sigma_{qbar g} = \sigma_{g qbar}
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      real*8 qrk1,qrk2,brk1,brk2
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      qrk1=0.d0
      qrk2=0.d0
      brk1=0.d0
      brk2=0.d0

      do i=1,5
         qrk1 = qrk1 + ff1(i)
         qrk2 = qrk2 + ff2(i)
         brk1 = brk1 + ff1(-i)
         brk2 = brk2 + ff2(-i)
      enddo

      qgpdf = ( qrk1 + brk1 )*ff2(0) + ff1(0)*( qrk2 + brk2 )

      return
      end

C-}}}
C-{{{ function qqbpdf(x1,x1):

      real*8 function qqbpdf(x1,x2)
C
C     q-qbar pdfs for equal quarks, or more accurately:
C
C     q(x1)*qbar(x2) + qbar(x1)*q(x2)
C
C     This can be done because
C
C     \sigma_{q qbar} = \sigma_{qbar q}
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      qqbpdf = 0.d0

      do i=1,5
         qqbpdf = qqbpdf + ff1(i)*ff2(-i) + ff1(-i)*ff2(i)
      enddo

      return
      end

C-}}}
C-{{{ function qqpdf(x1,x2):

      real*8 function qqpdf(x1,x2)
C
C     q-q pdfs, or more accurately:
C
C     q(x1)*q(x2) + qbar(x1)*qbar(x2).
C
C     This can be done because
C
C     \sigma_{q q} = \sigma_{qbar qbar}
C
C     Note that both quarks must be of the same flavor!
C     The case where the quarks have different flavors is treated in
C     <function qupdf>!
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      qqpdf = 0.d0

      do i=1,5
         qqpdf = qqpdf + ff1(i)*ff2(i) + ff1(-i)*ff2(-i)
      enddo

      return
      end

C-}}}
C-{{{ function qupdf(x1,x2):

      real*8 function qupdf(x1,x2)
C
C     q-q' pdfs, or more accurately:
C
C     q(x1)*q'(x2) + qbar(x1)*q'bar(x2).
C
C     This can be done because
C
C     \sigma_{q q'} = \sigma_{qbar q'bar}
C
C     Note that q and q' denote different flavors!
C     The case where q and q' have the same flavor is treated in
C     <function qqpdf>!
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      real*8 qrk1,qrk2,brk1,brk2
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      qrk1=0.d0
      qrk2=0.d0
      brk1=0.d0
      brk2=0.d0

      do i=1,5
         qrk1 = qrk1 + ff1(i)
         qrk2 = qrk2 + ff2(i)
         brk1 = brk1 + ff1(-i)
         brk2 = brk2 + ff2(-i)
      enddo

      qupdf = qrk1*qrk2 + brk1*brk2

      do i=1,5
         qupdf = qupdf - ff1(i)*ff2(i) - ff1(-i)*ff2(-i)
      enddo

      return
      end

C-}}}
C-{{{ function qubpdf(x1,x2):

      real*8 function qubpdf(x1,x2)
C
C     q-qbar pdfs for non-equal quarks, or more accurately:
C
C     q(x1)*q'bar(x2) + q'bar(x1)*q(x2)
C
C     This can be done because
C
C     \sigma_{q q'bar} = \sigma_{q'bar q}
C
C     Note that q and q' denote different flavors!
C     The case where q and q' have the same flavor is treated in
C     <function qbpdf>!
C
      implicit none
      real*8 x1,x2
      real*8 ff1(-6:6),ff2(-6:6)
      real*8 qrk1,qrk2,brk1,brk2
      integer sign,i
      include '../commons/common-vars.f'
      include '../commons/common-keys.f'

      sign=1
      if (ppbar) sign=-1

      call pdfs(1,norder,x1,rmuf,ff1)
      call pdfs(sign,norder,x2,rmuf,ff2)

      qrk1=0.d0
      qrk2=0.d0
      brk1=0.d0
      brk2=0.d0

      do i=1,5
         qrk1 = qrk1 + ff1(i)
         qrk2 = qrk2 + ff2(i)
         brk1 = brk1 + ff1(-i)
         brk2 = brk2 + ff2(-i)
      enddo

      qubpdf = qrk1*brk2 + brk1*qrk2

      do i=1,5
         qubpdf = qubpdf - ff1(i)*ff2(-i) - ff1(-i)*ff2(i)
      enddo

      return
      end

C-}}}
