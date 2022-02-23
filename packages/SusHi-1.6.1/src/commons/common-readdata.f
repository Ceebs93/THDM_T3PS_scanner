
      double precision sqscmsrd,mhiggsrd,murmhrd,mufmhrd
      double precision c1eff0rd,c1eff1rd,c1eff2rd,c1eff3rd,c2eff1rd
      complex*16 comc1eff0rd
      double precision mtoprd,mbpolerd,mbmbrd
      double precision mzrd,gfermird

!von bbh
      character*80 pdfnamerd
      logical lsloppyrd
      integer norderrd,nscalpseudrd,nmodelrd,ncolliderrd,npdfrd
     &     ,npdfstartrd,npdfendrd,nsubprocrd,nc1rd
      double precision Mstop1rd,Mstop2rd
      double precision gthrd,gth11rd,gth22rd,gth12rd,gth21rd,gbhrd
      double precision mbottomrd,ytoprd,ybottomrd
      character*20 pdfstringrd
!Ende von bbh
      
      common/rddata/sqscmsrd,mhiggsrd,murmhrd,mufmhrd,mtoprd
     &     ,mbottomrd,ytoprd,ybottomrd,mzrd,gfermird,
     &     comc1eff0rd,c1eff0rd,c1eff1rd,c1eff2rd,c1eff3rd,c2eff1rd,
     &     gthrd,gth11rd,gth22rd,gth12rd,gth21rd,
     &     gbhrd,mbpolerd,mbmbrd,
     &     Mstop1rd,Mstop2rd,
     &     norderrd,nscalpseudrd,nmodelrd,ncolliderrd,npdfrd,npdfstartrd
     &     ,npdfendrd,nsubprocrd,nc1rd,lsloppyrd
      common/rdchars/ pdfnamerd


