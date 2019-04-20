C-{{{ RCS:

c..   $Id: ggh-exp.f,v 3.3 2005/02/13 15:21:16 rharland Exp $
c..   $Log: ggh-exp.f,v $
c..   Revision 3.3  2005/02/13 15:21:16  rharland
c..   bug fix.
c..
c..   Revision 3.1  2004/11/26 15:15:58  rharland
c..   G- and sigma-expansion re-implemented
c..
c..   Revision 3.0  2004/11/18 13:41:48  rharland
c..   soft expansion included
c..

C-}}}
C-{{{ function ppgg2exp(yy)

      real*8 function ppgg2exp(xx)
c..
c..   Integrand for hard gg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external ggpdf,sgg2exp

      ppgg2exp = vartrans(ggpdf,sgg2exp,xx)
      
      end

C-}}}
C-{{{ function ppqg2exp(yy)

      real*8 function ppqg2exp(xx)
c..
c..   Integrand for hard qg-contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qgpdf,sqg2exp

      ppqg2exp = vartrans(qgpdf,sqg2exp,xx)

      end

C-}}}
C-{{{ function ppqqb2exp(yy)

      real*8 function ppqqb2exp(xx)
c..
c..   Integrand for hard q-qbar contribution at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqbpdf,sqqb2exp

      ppqqb2exp = vartrans(qqbpdf,sqqb2exp,xx)

      end

C-}}}
C-{{{ function ppqq2exp(yy)

      real*8 function ppqq2exp(xx)
c..
c..   Integrand for hard qq contribution (equal quarks) at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qqpdf,sqq2exp

      ppqq2exp = vartrans(qqpdf,sqq2exp,xx)

      end

C-}}}
C-{{{ function ppqu2exp(yy)

      real*8 function ppqu2exp(xx)
c..
c..   Integrand for hard qq contribution (non-equal quarks) at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      external qupdf,qubpdf,squ2exp

      ppqu2exp = vartrans(qupdf,squ2exp,xx) + vartrans(qubpdf,squ2exp,xx
     &     )

      end

C-}}}
C-{{{ function ppall2exp(yy)

      real*8 function ppall2exp(xx)
c..
c..   Sum of all integrands of the sub-processes at NNLO.
c..
      implicit real*8 (a-h,o-z)
      real*8 xx(10)
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      external ppgg2exp,ppqg2exp,ppqqb2exp,ppqq2exp,ppqu2exp

      ppall2exp = ppgg2exp(xx) + ppqg2exp(xx) + ppqqb2exp(xx)+
     &     ppqq2exp(xx) + ppqu2exp(xx)

      end

C-}}}
