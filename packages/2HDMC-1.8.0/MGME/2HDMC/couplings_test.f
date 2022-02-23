      subroutine testcoupl
      implicit none
c
c     local
c
      integer i,iformat
      character*2 ab(2)
      real*8 ene
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
c
c     include
c
      include 'coupl.inc'
      include 'input.inc'
c
c     Common to lh_readin and printout
c
      double precision  alpha, sin2w, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
      open(unit=1,file="couplings_check.txt")
c
c output all info
c
 10   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV' )
 11   format( 1x,a10,2x,f11.5,2x,f11.5,a3,f11.5,2x,f11.5 )
 12   format( 1x,a10,2x,f6.2,a )
 13   format( 1x,a10,2x,f6.4,a )
 14   format( 1x,2(a10,2x,f10.7,2x,f10.7) )
 15   format( 1x,a10,2x,f9.5,a )
 16   format( 1x,a10,2x,f7.5 )
 17   format( 1x,a10,2x,f8.4 )
 18   format( 1x,a10,2x,f8.4,' GeV' )
 19   format( 1x,a10,2x,f6.4,a13,2x,f6.4 )
 20   format( 1x,a10,2x,f11.5,1x,f11.5 )
 21   format( 1x,a10,2x,f8.4,' GeV',1x,a13,2x,f8.4,' GeV' )
 22   format( 1x,a10,2x,f10.8,a13,2x,f6.4 )
 23   format( 1x,a10,2x,f8.4)
 24   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV  (calc @ LO)')
 25   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV')


      write(1,11) 'GH1EE ', GH1EE(1),' ', GH1EE(2)
      write(1,11) 'GH2EE ', GH2EE(1),' ', GH2EE(2)
      write(1,11) 'GH3EE ', GH3EE(1),' ', GH3EE(2)
      write(1,11) 'GH1MUMU ', GH1MUMU(1),' ', GH1MUMU(2)
      write(1,11) 'GH2MUMU ', GH2MUMU(1),' ', GH2MUMU(2)
      write(1,11) 'GH3MUMU ', GH3MUMU(1),' ', GH3MUMU(2)
      write(1,11) 'GH1TATA ', GH1TATA(1),' ', GH1TATA(2)
      write(1,11) 'GH2TATA ', GH2TATA(1),' ', GH2TATA(2)
      write(1,11) 'GH3TATA ', GH3TATA(1),' ', GH3TATA(2)
      write(1,11) 'GH1UU ', GH1UU(1),' ', GH1UU(2)
      write(1,11) 'GH2UU ', GH2UU(1),' ', GH2UU(2)
      write(1,11) 'GH3UU ', GH3UU(1),' ', GH3UU(2)
      write(1,11) 'GH1CC ', GH1CC(1),' ', GH1CC(2)
      write(1,11) 'GH2CC ', GH2CC(1),' ', GH2CC(2)
      write(1,11) 'GH3CC ', GH3CC(1),' ', GH3CC(2)
      write(1,11) 'GH1TT ', GH1TT(1),' ', GH1TT(2)
      write(1,11) 'GH2TT ', GH2TT(1),' ', GH2TT(2)
      write(1,11) 'GH3TT ', GH3TT(1),' ', GH3TT(2)
      write(1,11) 'GH1DD ', GH1DD(1),' ', GH1DD(2)
      write(1,11) 'GH2DD ', GH2DD(1),' ', GH2DD(2)
      write(1,11) 'GH3DD ', GH3DD(1),' ', GH3DD(2)
      write(1,11) 'GH1SS ', GH1SS(1),' ', GH1SS(2)
      write(1,11) 'GH2SS ', GH2SS(1),' ', GH2SS(2)
      write(1,11) 'GH3SS ', GH3SS(1),' ', GH3SS(2)
      write(1,11) 'GH1BB ', GH1BB(1),' ', GH1BB(2)
      write(1,11) 'GH2BB ', GH2BB(1),' ', GH2BB(2)
      write(1,11) 'GH3BB ', GH3BB(1),' ', GH3BB(2)
      write(1,11) 'GHCUD ', GHCUD(1),' ', GHCUD(2)
      write(1,11) 'GHCUS ', GHCUS(1),' ', GHCUS(2)
      write(1,11) 'GHCUB ', GHCUB(1),' ', GHCUB(2)
      write(1,11) 'GHCCD ', GHCCD(1),' ', GHCCD(2)
      write(1,11) 'GHCCS ', GHCCS(1),' ', GHCCS(2)
      write(1,11) 'GHCCB ', GHCCB(1),' ', GHCCB(2)
      write(1,11) 'GHCTD ', GHCTD(1),' ', GHCTD(2)
      write(1,11) 'GHCTS ', GHCTS(1),' ', GHCTS(2)
      write(1,11) 'GHCTB ', GHCTB(1),' ', GHCTB(2)
      write(1,11) 'GHCDU ', GHCDU(1),' ', GHCDU(2)
      write(1,11) 'GHCDC ', GHCDC(1),' ', GHCDC(2)
      write(1,11) 'GHCDT ', GHCDT(1),' ', GHCDT(2)
      write(1,11) 'GHCSU ', GHCSU(1),' ', GHCSU(2)
      write(1,11) 'GHCSC ', GHCSC(1),' ', GHCSC(2)
      write(1,11) 'GHCST ', GHCST(1),' ', GHCST(2)
      write(1,11) 'GHCBU ', GHCBU(1),' ', GHCBU(2)
      write(1,11) 'GHCBC ', GHCBC(1),' ', GHCBC(2)
      write(1,11) 'GHCBT ', GHCBT(1),' ', GHCBT(2)
      write(1,11) 'GHCVEE ', GHCVEE(1),' ', GHCVEE(2)
      write(1,11) 'GHCVMMU ', GHCVMMU(1),' ', GHCVMMU(2)
      write(1,11) 'GHCVTTA ', GHCVTTA(1),' ', GHCVTTA(2)
      write(1,11) 'GHCEVE ', GHCEVE(1),' ', GHCEVE(2)
      write(1,11) 'GHCMUVM ', GHCMUVM(1),' ', GHCMUVM(2)
      write(1,11) 'GHCTAVT ', GHCTAVT(1),' ', GHCTAVT(2)
      write(1,20) 'GWWH1 ', GWWH1
      write(1,20) 'GWWH2 ', GWWH2
      write(1,20) 'GZZH1 ', GZZH1
      write(1,20) 'GZZH2 ', GZZH2
      write(1,20) 'GAHCHC ', GAHCHC
      write(1,20) 'GZH1H2 ', GZH1H2
      write(1,20) 'GZH1H3 ', GZH1H3
      write(1,20) 'GZH2H3 ', GZH2H3
      write(1,20) 'GZHCHC ', GZHCHC
      write(1,20) 'GWPHCH1 ', GWPHCH1
      write(1,20) 'GWMH1HC ', GWMH1HC
      write(1,20) 'GWPHCH2 ', GWPHCH2
      write(1,20) 'GWMH2HC ', GWMH2HC
      write(1,20) 'GWPHCH3 ', GWPHCH3
      write(1,20) 'GWMH3HC ', GWMH3HC
      write(1,20) 'GH1H1H1 ', GH1H1H1
      write(1,20) 'GH1H1H2 ', GH1H1H2
      write(1,20) 'GH1H2H2 ', GH1H2H2
      write(1,20) 'GH1H3H3 ', GH1H3H3
      write(1,20) 'GH2H2H2 ', GH2H2H2
      write(1,20) 'GH2H3H3 ', GH2H3H3
      write(1,20) 'GH1HCHC ', GH1HCHC
      write(1,20) 'GH2HCHC ', GH2HCHC
      write(1,20) 'GAAHCHC ', GAAHCHC
      write(1,20) 'GAZHCHC ', GAZHCHC
      write(1,20) 'GAWPHCH1 ', GAWPHCH1
      write(1,20) 'GAWPHCH2 ', GAWPHCH2
      write(1,20) 'GAWPHCH3 ', GAWPHCH3
      write(1,20) 'GAWMH1HC ', GAWMH1HC
      write(1,20) 'GAWMH2HC ', GAWMH2HC
      write(1,20) 'GAWMH3HC ', GAWMH3HC
      write(1,20) 'GZZH1H1 ', GZZH1H1
      write(1,20) 'GZZH2H2 ', GZZH2H2
      write(1,20) 'GZZH3H3 ', GZZH3H3
      write(1,20) 'GZZHCHC ', GZZHCHC
      write(1,20) 'GZWPHCH1 ', GZWPHCH1
      write(1,20) 'GZWPHCH2 ', GZWPHCH2
      write(1,20) 'GZWPHCH3 ', GZWPHCH3
      write(1,20) 'GZWMH1HC ', GZWMH1HC
      write(1,20) 'GZWMH2HC ', GZWMH2HC
      write(1,20) 'GZWMH3HC ', GZWMH3HC
      write(1,20) 'GWWH1H1 ', GWWH1H1
      write(1,20) 'GWWH2H2 ', GWWH2H2
      write(1,20) 'GWWH3H3 ', GWWH3H3
      write(1,20) 'GWWHCHC ', GWWHCHC
      write(1,20) 'gwwa ', gwwa
      write(1,20) 'gwwz ', gwwz
      write(1,20) 'gwwh ', gwwh
      write(1,20) 'gzzh ', gzzh
      write(1,20) 'ghhh ', ghhh
      write(1,20) 'gwwhh ', gwwhh
      write(1,20) 'gzzhh ', gzzhh
      write(1,20) 'ghhhh ', ghhhh
      write(1,11) 'gal ',gal(1),' ',gal(2)
      write(1,11) 'gau ',gau(1),' ',gau(2)
      write(1,11) 'gad ',gad(1),' ',gad(2)
      write(1,11) 'gwf ',gwf(1),' ',gwf(2)
      write(1,11) 'gwfud ',gwfud(1),' ',gwfud(2)
      write(1,11) 'gwfus ',gwfus(1),' ',gwfus(2)
      write(1,11) 'gwfub ',gwfub(1),' ',gwfub(2)
      write(1,11) 'gwfcd ',gwfcd(1),' ',gwfcd(2)
      write(1,11) 'gwfcs ',gwfcs(1),' ',gwfcs(2)
      write(1,11) 'gwfcb ',gwfcb(1),' ',gwfcb(2)
      write(1,11) 'gwftd ',gwftd(1),' ',gwftd(2)
      write(1,11) 'gwfts ',gwfts(1),' ',gwfts(2)
      write(1,11) 'gwftb ',gwftb(1),' ',gwftb(2)
      write(1,11) 'gzn ',gzn(1),' ',gzn(2)
      write(1,11) 'gzl ',gzl(1),' ',gzl(2)
      write(1,11) 'gzu ',gzu(1),' ',gzu(2)
      write(1,11) 'gzd ',gzd(1),' ',gzd(2)
      write(1,14) 'gHtop ',ghtop(1),' ',ghtop(2)
      write(1,14) 'gHbot ',ghbot(1),' ',ghbot(2)
      write(1,14) 'gHcha ',ghcha(1),' ',ghcha(2)
      write(1,14) 'gHtau ',ghtau(1),' ',ghtau(2)
      write(1,14) 'gg ',gg(1),' ',gg(2)

      return
      end
