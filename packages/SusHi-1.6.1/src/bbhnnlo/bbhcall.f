
      subroutine bbhcall(norderin,ncolliderin,sqrtsin,mhin,murfacin
     &     ,muffacin,mbottomin,gfermiin,mzin,sigtot,erroutbbh)
c..
      implicit none
      integer norderin,ncolliderin
      real*8 sqrtsin,mhin,murfacin,muffacin,mbottomin,gfermiin,mzin
     &     ,sigtot,erroutbbh
      include '../commons/common-readdata.f'
      include '../commons/common-sigma.f'
      include '../commons/common-keys.f'
      include '../commons/common-citations.f'

      citations(8) = 1
      norderrd = norderin
      ncolliderrd = ncolliderin
      sqscmsrd = sqrtsin
      mhiggsrd = mhin
      murmhrd = murfacin
      mufmhrd = muffacin
      mbottomrd = mbottomin
      gfermird = gfermiin
      mzrd = mzin

      !nscalpseudrd = pseudoscin

      lsloppyrd = .false.
      ybottomrd = 1.d0
      nc1rd = 0
      gthrd = 1.d0
      gbhrd = 1.d0
c$$$      nsubprocrd = 0

      call initbbh()
      call evalsigmabbh()

      sigtot = sigmabbh(norderin)
      erroutbbh = serrbbh(norderin)

      end

