      subroutine printout
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
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
c output all info
c
 10   format( 1x,a16,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV' )
 11   format( 1x,a13,2x,f11.5,2x,f11.5,2x,a13,2x,f11.5,2x,f11.5 )
 12   format( 1x,a13,2x,f6.2,a )
 13   format( 1x,a13,2x,f6.4,a )
 14   format( 1x,2(a13,2x,f10.7,2x,f10.7) )
 15   format( 1x,a13,2x,f9.5,a )
 16   format( 1x,a13,2x,f7.5 )
 17   format( 1x,a13,2x,f8.4 )
 18   format( 1x,a13,2x,f8.4,' GeV' )
 19   format( 1x,a13,2x,f6.4,a13,2x,f6.4 )
 20   format( 1x,a13,2x,f13.5,1x,f13.5 )
 21   format( 1x,a13,2x,f8.4,' GeV',1x,a13,2x,f8.4,' GeV' )
 22   format( 1x,a13,2x,f10.8,a13,2x,f6.4 )
 23   format( 1x,a13,2x,f8.4)
 24   format( 1x,a16,2x,f9.3,' GeV        ',a16,2x,f7.4,' GeV  (calc @ LO)')
 25   format( 1x,a16,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV')


      write(*,*) '*****************************************************'
      write(*,*) '*                    MadEvent                       *'
      write(*,*) '*        --------------------------------           *'
      write(*,*) '*          http://madgraph.hep.uiuc.edu             *'
      write(*,*) '*          http://madgraph.phys.ucl.ac.be           *'
      write(*,*) '*          http://madgraph.roma2.infn.it            *'
      write(*,*) '*        --------------------------------           *'
      write(*,*) '*                                                   *'
      write(*,*) '*         INTEGRATION CHANNEL LOG FILE              *'
      write(*,*) '*                                                   *'
      write(*,*) '*****************************************************'
      write(6,*)
      write(*,*) '*****************************************************'
      write(*,*) '*          SUMMARY OF THE SM PARAMETERS             *'
      write(*,*) '*****************************************************'
      write(6,*)
      write(6,*)  ' EW Parameters:'
      write(6,*)  '---------------'
      write(6,*)
      write(6,23) ' GF (10^-5*GeV^-2) = ',gfermi*1d5
      write(6,23) ' 1/alpha           = ',1d0/alpha
      write(6,23) ' M_Z   (GeV)       = ',zmass
      write(6,*)
      write(6,*)
      write(6,*)  'Boson masses and widths:'
      write(6,*)  '------------------------'
      write(6,*)
      write(6,24) 'Z mass  =  ',zmass, 'Z width  = ',zwidth
      write(6,24) 'W mass  =  ',wmass, 'W width  = ',wwidth
      write(6,24) 'H mass  =  ',hmass, 'H width  = ',hwidth
      write(6,*)
      write(6,*)  'Fermion masses and widths:'
      write(6,*)  '--------------------------'
      write(6,*)
      write(6,24) 'top    mass  =  ', tmass, 'top    width  = ', twidth
      write(6,10) 'bottom mass  =  ', bmass, 'bottom width  = ', Zero
      write(6,10) 'charm  mass  =  ', cmass, 'charm  width  = ', Zero
      write(6,10) 'tau    mass  =  ', lmass, 'tau    width  = ', Zero
      write(6,*)  'all other quark and lepton masses set to zero'
      write(6,*)
      write(6,*)  'User Model couplings:'
      write(6,*)  '---------------------'
      write(6,*)
      write(6,11) 'GH1EE(L)  =  ', GH1EE(1),'GH1EE(R)  =  ', GH1EE(2)
      write(6,11) 'GH2EE(L)  =  ', GH2EE(1),'GH2EE(R)  =  ', GH2EE(2)
      write(6,11) 'GH3EE(L)  =  ', GH3EE(1),'GH3EE(R)  =  ', GH3EE(2)
      write(6,11) 'GH1MUMU(L)  =  ', GH1MUMU(1),'GH1MUMU(R)  =  ', GH1MUMU(2)
      write(6,11) 'GH2MUMU(L)  =  ', GH2MUMU(1),'GH2MUMU(R)  =  ', GH2MUMU(2)
      write(6,11) 'GH3MUMU(L)  =  ', GH3MUMU(1),'GH3MUMU(R)  =  ', GH3MUMU(2)
      write(6,11) 'GH1TATA(L)  =  ', GH1TATA(1),'GH1TATA(R)  =  ', GH1TATA(2)
      write(6,11) 'GH2TATA(L)  =  ', GH2TATA(1),'GH2TATA(R)  =  ', GH2TATA(2)
      write(6,11) 'GH3TATA(L)  =  ', GH3TATA(1),'GH3TATA(R)  =  ', GH3TATA(2)
      write(6,11) 'GH1UU(L)  =  ', GH1UU(1),'GH1UU(R)  =  ', GH1UU(2)
      write(6,11) 'GH2UU(L)  =  ', GH2UU(1),'GH2UU(R)  =  ', GH2UU(2)
      write(6,11) 'GH3UU(L)  =  ', GH3UU(1),'GH3UU(R)  =  ', GH3UU(2)
      write(6,11) 'GH1CC(L)  =  ', GH1CC(1),'GH1CC(R)  =  ', GH1CC(2)
      write(6,11) 'GH2CC(L)  =  ', GH2CC(1),'GH2CC(R)  =  ', GH2CC(2)
      write(6,11) 'GH3CC(L)  =  ', GH3CC(1),'GH3CC(R)  =  ', GH3CC(2)
      write(6,11) 'GH1TT(L)  =  ', GH1TT(1),'GH1TT(R)  =  ', GH1TT(2)
      write(6,11) 'GH2TT(L)  =  ', GH2TT(1),'GH2TT(R)  =  ', GH2TT(2)
      write(6,11) 'GH3TT(L)  =  ', GH3TT(1),'GH3TT(R)  =  ', GH3TT(2)
      write(6,11) 'GH1DD(L)  =  ', GH1DD(1),'GH1DD(R)  =  ', GH1DD(2)
      write(6,11) 'GH2DD(L)  =  ', GH2DD(1),'GH2DD(R)  =  ', GH2DD(2)
      write(6,11) 'GH3DD(L)  =  ', GH3DD(1),'GH3DD(R)  =  ', GH3DD(2)
      write(6,11) 'GH1SS(L)  =  ', GH1SS(1),'GH1SS(R)  =  ', GH1SS(2)
      write(6,11) 'GH2SS(L)  =  ', GH2SS(1),'GH2SS(R)  =  ', GH2SS(2)
      write(6,11) 'GH3SS(L)  =  ', GH3SS(1),'GH3SS(R)  =  ', GH3SS(2)
      write(6,11) 'GH1BB(L)  =  ', GH1BB(1),'GH1BB(R)  =  ', GH1BB(2)
      write(6,11) 'GH2BB(L)  =  ', GH2BB(1),'GH2BB(R)  =  ', GH2BB(2)
      write(6,11) 'GH3BB(L)  =  ', GH3BB(1),'GH3BB(R)  =  ', GH3BB(2)
      write(6,11) 'GHCUD(L)  =  ', GHCUD(1),'GHCUD(R)  =  ', GHCUD(2)
      write(6,11) 'GHCUS(L)  =  ', GHCUS(1),'GHCUS(R)  =  ', GHCUS(2)
      write(6,11) 'GHCUB(L)  =  ', GHCUB(1),'GHCUB(R)  =  ', GHCUB(2)
      write(6,11) 'GHCCD(L)  =  ', GHCCD(1),'GHCCD(R)  =  ', GHCCD(2)
      write(6,11) 'GHCCS(L)  =  ', GHCCS(1),'GHCCS(R)  =  ', GHCCS(2)
      write(6,11) 'GHCCB(L)  =  ', GHCCB(1),'GHCCB(R)  =  ', GHCCB(2)
      write(6,11) 'GHCTD(L)  =  ', GHCTD(1),'GHCTD(R)  =  ', GHCTD(2)
      write(6,11) 'GHCTS(L)  =  ', GHCTS(1),'GHCTS(R)  =  ', GHCTS(2)
      write(6,11) 'GHCTB(L)  =  ', GHCTB(1),'GHCTB(R)  =  ', GHCTB(2)
      write(6,11) 'GHCDU(L)  =  ', GHCDU(1),'GHCDU(R)  =  ', GHCDU(2)
      write(6,11) 'GHCDC(L)  =  ', GHCDC(1),'GHCDC(R)  =  ', GHCDC(2)
      write(6,11) 'GHCDT(L)  =  ', GHCDT(1),'GHCDT(R)  =  ', GHCDT(2)
      write(6,11) 'GHCSU(L)  =  ', GHCSU(1),'GHCSU(R)  =  ', GHCSU(2)
      write(6,11) 'GHCSC(L)  =  ', GHCSC(1),'GHCSC(R)  =  ', GHCSC(2)
      write(6,11) 'GHCST(L)  =  ', GHCST(1),'GHCST(R)  =  ', GHCST(2)
      write(6,11) 'GHCBU(L)  =  ', GHCBU(1),'GHCBU(R)  =  ', GHCBU(2)
      write(6,11) 'GHCBC(L)  =  ', GHCBC(1),'GHCBC(R)  =  ', GHCBC(2)
      write(6,11) 'GHCBT(L)  =  ', GHCBT(1),'GHCBT(R)  =  ', GHCBT(2)
      write(6,11) 'GHCVEE(L)  =  ', GHCVEE(1),'GHCVEE(R)  =  ', GHCVEE(2)
      write(6,11) 'GHCVMMU(L)  =  ', GHCVMMU(1),'GHCVMMU(R)  =  ', GHCVMMU(2)
      write(6,11) 'GHCVTTA(L)  =  ', GHCVTTA(1),'GHCVTTA(R)  =  ', GHCVTTA(2)
      write(6,11) 'GHCEVE(L)  =  ', GHCEVE(1),'GHCEVE(R)  =  ', GHCEVE(2)
      write(6,11) 'GHCMUVM(L)  =  ', GHCMUVM(1),'GHCMUVM(R)  =  ', GHCMUVM(2)
      write(6,11) 'GHCTAVT(L)  =  ', GHCTAVT(1),'GHCTAVT(R)  =  ', GHCTAVT(2)
      write(6,20) 'GWWH1  =  ', GWWH1
      write(6,20) 'GWWH2  =  ', GWWH2
      write(6,20) 'GZZH1  =  ', GZZH1
      write(6,20) 'GZZH2  =  ', GZZH2
      write(6,20) 'GAHCHC  =  ', GAHCHC
      write(6,20) 'GZH1H2  =  ', GZH1H2
      write(6,20) 'GZH1H3  =  ', GZH1H3
      write(6,20) 'GZH2H3  =  ', GZH2H3
      write(6,20) 'GZHCHC  =  ', GZHCHC
      write(6,20) 'GWPHCH1  =  ', GWPHCH1
      write(6,20) 'GWMH1HC  =  ', GWMH1HC
      write(6,20) 'GWPHCH2  =  ', GWPHCH2
      write(6,20) 'GWMH2HC  =  ', GWMH2HC
      write(6,20) 'GWPHCH3  =  ', GWPHCH3
      write(6,20) 'GWMH3HC  =  ', GWMH3HC
      write(6,20) 'GH1H1H1  =  ', GH1H1H1
      write(6,20) 'GH1H1H2  =  ', GH1H1H2
      write(6,20) 'GH1H2H2  =  ', GH1H2H2
      write(6,20) 'GH1H3H3  =  ', GH1H3H3
      write(6,20) 'GH2H2H2  =  ', GH2H2H2
      write(6,20) 'GH2H3H3  =  ', GH2H3H3
      write(6,20) 'GH1HCHC  =  ', GH1HCHC
      write(6,20) 'GH2HCHC  =  ', GH2HCHC
      write(6,20) 'GAAHCHC  =  ', GAAHCHC
      write(6,20) 'GAZHCHC  =  ', GAZHCHC
      write(6,20) 'GAWPHCH1  =  ', GAWPHCH1
      write(6,20) 'GAWPHCH2  =  ', GAWPHCH2
      write(6,20) 'GAWPHCH3  =  ', GAWPHCH3
      write(6,20) 'GAWMH1HC  =  ', GAWMH1HC
      write(6,20) 'GAWMH2HC  =  ', GAWMH2HC
      write(6,20) 'GAWMH3HC  =  ', GAWMH3HC
      write(6,20) 'GZZH1H1  =  ', GZZH1H1
      write(6,20) 'GZZH2H2  =  ', GZZH2H2
      write(6,20) 'GZZH3H3  =  ', GZZH3H3
      write(6,20) 'GZZHCHC  =  ', GZZHCHC
      write(6,20) 'GZWPHCH1  =  ', GZWPHCH1
      write(6,20) 'GZWPHCH2  =  ', GZWPHCH2
      write(6,20) 'GZWPHCH3  =  ', GZWPHCH3
      write(6,20) 'GZWMH1HC  =  ', GZWMH1HC
      write(6,20) 'GZWMH2HC  =  ', GZWMH2HC
      write(6,20) 'GZWMH3HC  =  ', GZWMH3HC
      write(6,20) 'GWWH1H1  =  ', GWWH1H1
      write(6,20) 'GWWH2H2  =  ', GWWH2H2
      write(6,20) 'GWWH3H3  =  ', GWWH3H3
      write(6,20) 'GWWHCHC  =  ', GWWHCHC
      write(6,*)
      write(6,*) 'Boson couplings:'
      write(6,*) '----------------'
      write(6,*)
      write(6,20) 'gwwa  = ', gwwa
      write(6,20) 'gwwz  = ', gwwz
      write(6,20) 'gwwh  = ', gwwh
      write(6,20) 'gzzh  = ', gzzh
      write(6,20) 'ghhh  = ', ghhh
      write(6,*)
      write(6,20) 'gwwhh = ', gwwhh
      write(6,20) 'gzzhh = ', gzzhh
      write(6,20) 'ghhhh = ', ghhhh
      write(6,*)
      write(6,*) 'FFV couplings:'
      write(6,*) '--------------'
      write(6,*)
      write(6,11) 'gal(L)   =  ',gal(1), 'gal(R)   =  ',gal(2)
      write(6,11) 'gau(L)   =  ',gau(1), 'gau(R)   =  ',gau(2)
      write(6,11) 'gad(L)   =  ',gad(1), 'gad(R)   =  ',gad(2)
      write(6,*)
      write(6,11) 'gwf(L)   =  ',gwf(1), 'gwf(R)   =  ',gwf(2)
      write(6,*)
      write(6,11) 'gwfud(L) =  ',gwfud(1), 'gwfud(R) =  ',gwfud(2)
      write(6,11) 'gwfus(L) =  ',gwfus(1), 'gwfud(R) =  ',gwfus(2)
      write(6,11) 'gwfub(L) =  ',gwfub(1), 'gwfud(R) =  ',gwfub(2)
      write(6,11) 'gwfcd(L) =  ',gwfcd(1), 'gwfcd(R) =  ',gwfcd(2)
      write(6,11) 'gwfcs(L) =  ',gwfcs(1), 'gwfcd(R) =  ',gwfcs(2)
      write(6,11) 'gwfcb(L) =  ',gwfcb(1), 'gwfcd(R) =  ',gwfcb(2)
      write(6,11) 'gwftd(L) =  ',gwftd(1), 'gwftd(R) =  ',gwftd(2)
      write(6,11) 'gwfts(L) =  ',gwfts(1), 'gwftd(R) =  ',gwfts(2)
      write(6,11) 'gwftb(L) =  ',gwftb(1), 'gwftd(R) =  ',gwftb(2)
      write(6,*)
      write(6,11) 'gzn(L)   =  ',gzn(1), 'gzn(R)   =  ',gzn(2)
      write(6,11) 'gzl(L)   =  ',gzl(1), 'gzl(R)   =  ',gzl(2)
      write(6,11) 'gzu(L)   =  ',gzu(1), 'gzu(R)   =  ',gzu(2)
      write(6,11) 'gzd(L)   =  ',gzd(1), 'gzd(R)   =  ',gzd(2)
      write(6,*)
      write(6,*) 'FFH couplings:'
      write(6,*) '--------------'
      write(6,*)
      write(6,14) 'gHtop(L) =  ',ghtop(1), 'gHtop(R) =  ',ghtop(2)
      write(6,14) 'gHbot(L) =  ',ghbot(1), 'gHbot(R) =  ',ghbot(2)
      write(6,14) 'gHcha(L) =  ',ghcha(1), 'gHcha(R) =  ',ghcha(2)
      write(6,14) 'gHtau(L) =  ',ghtau(1), 'gHtau(R) =  ',ghtau(2)
      write(6,*)
      write(6,*) 'Strong couplings:'
      write(6,*) '-----------------'
      write(6,*)
      write(6,14) 'gg(1)    =  ',gg(1)   , 'gg(2)    =  ',gg(2)
      write(6,*)
      write(6,*) 'User Model masses and widths:'
      write(6,*) '-----------------------------'
      write(6,*)
      write(6,24) 'H1MASS  =  ', H1MASS, 'H1WIDTH    = ', H1WIDTH
      write(6,24) 'H2MASS  =  ', H2MASS, 'H2WIDTH    = ', H2WIDTH
      write(6,24) 'H3MASS  =  ', H3MASS, 'H3WIDTH    = ', H3WIDTH
      write(6,24) 'HCMASS  =  ', HCMASS, 'HCWIDTH    = ', HCWIDTH
      write(6,*)
      write(6,*) 'User Model Parameters:'
      write(6,*) '----------------------'
      write(6,24) 'RESGH1EE     =  ', RESGH1EE   
      write(6,24) 'IMSGH1EE     =  ', IMSGH1EE   
      write(6,24) 'REPGH1EE     =  ', REPGH1EE   
      write(6,24) 'IMPGH1EE     =  ', IMPGH1EE   
      write(6,24) 'RESGH2EE     =  ', RESGH2EE   
      write(6,24) 'IMSGH2EE     =  ', IMSGH2EE   
      write(6,24) 'REPGH2EE     =  ', REPGH2EE   
      write(6,24) 'IMPGH2EE     =  ', IMPGH2EE   
      write(6,24) 'RESGH3EE     =  ', RESGH3EE   
      write(6,24) 'IMSGH3EE     =  ', IMSGH3EE   
      write(6,24) 'REPGH3EE     =  ', REPGH3EE   
      write(6,24) 'IMPGH3EE     =  ', IMPGH3EE   
      write(6,24) 'RESGH1MUMU   =  ', RESGH1MUMU 
      write(6,24) 'IMSGH1MUMU   =  ', IMSGH1MUMU 
      write(6,24) 'REPGH1MUMU   =  ', REPGH1MUMU 
      write(6,24) 'IMPGH1MUMU   =  ', IMPGH1MUMU 
      write(6,24) 'RESGH2MUMU   =  ', RESGH2MUMU 
      write(6,24) 'IMSGH2MUMU   =  ', IMSGH2MUMU 
      write(6,24) 'REPGH2MUMU   =  ', REPGH2MUMU 
      write(6,24) 'IMPGH2MUMU   =  ', IMPGH2MUMU 
      write(6,24) 'RESGH3MUMU   =  ', RESGH3MUMU 
      write(6,24) 'IMSGH3MUMU   =  ', IMSGH3MUMU 
      write(6,24) 'REPGH3MUMU   =  ', REPGH3MUMU 
      write(6,24) 'IMPGH3MUMU   =  ', IMPGH3MUMU 
      write(6,24) 'RESGH1TATA   =  ', RESGH1TATA 
      write(6,24) 'IMSGH1TATA   =  ', IMSGH1TATA 
      write(6,24) 'REPGH1TATA   =  ', REPGH1TATA 
      write(6,24) 'IMPGH1TATA   =  ', IMPGH1TATA 
      write(6,24) 'RESGH2TATA   =  ', RESGH2TATA 
      write(6,24) 'IMSGH2TATA   =  ', IMSGH2TATA 
      write(6,24) 'REPGH2TATA   =  ', REPGH2TATA 
      write(6,24) 'IMPGH2TATA   =  ', IMPGH2TATA 
      write(6,24) 'RESGH3TATA   =  ', RESGH3TATA 
      write(6,24) 'IMSGH3TATA   =  ', IMSGH3TATA 
      write(6,24) 'REPGH3TATA   =  ', REPGH3TATA 
      write(6,24) 'IMPGH3TATA   =  ', IMPGH3TATA 
      write(6,24) 'RESGH1UU     =  ', RESGH1UU   
      write(6,24) 'IMSGH1UU     =  ', IMSGH1UU   
      write(6,24) 'REPGH1UU     =  ', REPGH1UU   
      write(6,24) 'IMPGH1UU     =  ', IMPGH1UU   
      write(6,24) 'RESGH2UU     =  ', RESGH2UU   
      write(6,24) 'IMSGH2UU     =  ', IMSGH2UU   
      write(6,24) 'REPGH2UU     =  ', REPGH2UU   
      write(6,24) 'IMPGH2UU     =  ', IMPGH2UU   
      write(6,24) 'RESGH3UU     =  ', RESGH3UU   
      write(6,24) 'IMSGH3UU     =  ', IMSGH3UU   
      write(6,24) 'REPGH3UU     =  ', REPGH3UU   
      write(6,24) 'IMPGH3UU     =  ', IMPGH3UU   
      write(6,24) 'RESGH1CC     =  ', RESGH1CC   
      write(6,24) 'IMSGH1CC     =  ', IMSGH1CC   
      write(6,24) 'REPGH1CC     =  ', REPGH1CC   
      write(6,24) 'IMPGH1CC     =  ', IMPGH1CC   
      write(6,24) 'RESGH2CC     =  ', RESGH2CC   
      write(6,24) 'IMSGH2CC     =  ', IMSGH2CC   
      write(6,24) 'REPGH2CC     =  ', REPGH2CC   
      write(6,24) 'IMPGH2CC     =  ', IMPGH2CC   
      write(6,24) 'RESGH3CC     =  ', RESGH3CC   
      write(6,24) 'IMSGH3CC     =  ', IMSGH3CC   
      write(6,24) 'REPGH3CC     =  ', REPGH3CC   
      write(6,24) 'IMPGH3CC     =  ', IMPGH3CC   
      write(6,24) 'RESGH1TT     =  ', RESGH1TT   
      write(6,24) 'IMSGH1TT     =  ', IMSGH1TT   
      write(6,24) 'REPGH1TT     =  ', REPGH1TT   
      write(6,24) 'IMPGH1TT     =  ', IMPGH1TT   
      write(6,24) 'RESGH2TT     =  ', RESGH2TT   
      write(6,24) 'IMSGH2TT     =  ', IMSGH2TT   
      write(6,24) 'REPGH2TT     =  ', REPGH2TT   
      write(6,24) 'IMPGH2TT     =  ', IMPGH2TT   
      write(6,24) 'RESGH3TT     =  ', RESGH3TT   
      write(6,24) 'IMSGH3TT     =  ', IMSGH3TT   
      write(6,24) 'REPGH3TT     =  ', REPGH3TT   
      write(6,24) 'IMPGH3TT     =  ', IMPGH3TT   
      write(6,24) 'RESGH1DD     =  ', RESGH1DD   
      write(6,24) 'IMSGH1DD     =  ', IMSGH1DD   
      write(6,24) 'REPGH1DD     =  ', REPGH1DD   
      write(6,24) 'IMPGH1DD     =  ', IMPGH1DD   
      write(6,24) 'RESGH2DD     =  ', RESGH2DD   
      write(6,24) 'IMSGH2DD     =  ', IMSGH2DD   
      write(6,24) 'REPGH2DD     =  ', REPGH2DD   
      write(6,24) 'IMPGH2DD     =  ', IMPGH2DD   
      write(6,24) 'RESGH3DD     =  ', RESGH3DD   
      write(6,24) 'IMSGH3DD     =  ', IMSGH3DD   
      write(6,24) 'REPGH3DD     =  ', REPGH3DD   
      write(6,24) 'IMPGH3DD     =  ', IMPGH3DD   
      write(6,24) 'RESGH1SS     =  ', RESGH1SS   
      write(6,24) 'IMSGH1SS     =  ', IMSGH1SS   
      write(6,24) 'REPGH1SS     =  ', REPGH1SS   
      write(6,24) 'IMPGH1SS     =  ', IMPGH1SS   
      write(6,24) 'RESGH2SS     =  ', RESGH2SS   
      write(6,24) 'IMSGH2SS     =  ', IMSGH2SS   
      write(6,24) 'REPGH2SS     =  ', REPGH2SS   
      write(6,24) 'IMPGH2SS     =  ', IMPGH2SS   
      write(6,24) 'RESGH3SS     =  ', RESGH3SS   
      write(6,24) 'IMSGH3SS     =  ', IMSGH3SS   
      write(6,24) 'REPGH3SS     =  ', REPGH3SS   
      write(6,24) 'IMPGH3SS     =  ', IMPGH3SS   
      write(6,24) 'RESGH1BB     =  ', RESGH1BB   
      write(6,24) 'IMSGH1BB     =  ', IMSGH1BB   
      write(6,24) 'REPGH1BB     =  ', REPGH1BB   
      write(6,24) 'IMPGH1BB     =  ', IMPGH1BB   
      write(6,24) 'RESGH2BB     =  ', RESGH2BB   
      write(6,24) 'IMSGH2BB     =  ', IMSGH2BB   
      write(6,24) 'REPGH2BB     =  ', REPGH2BB   
      write(6,24) 'IMPGH2BB     =  ', IMPGH2BB   
      write(6,24) 'RESGH3BB     =  ', RESGH3BB   
      write(6,24) 'IMSGH3BB     =  ', IMSGH3BB   
      write(6,24) 'REPGH3BB     =  ', REPGH3BB   
      write(6,24) 'IMPGH3BB     =  ', IMPGH3BB   
      write(6,24) 'RESGHCUD     =  ', RESGHCUD   
      write(6,24) 'IMSGHCUD     =  ', IMSGHCUD   
      write(6,24) 'REPGHCUD     =  ', REPGHCUD   
      write(6,24) 'IMPGHCUD     =  ', IMPGHCUD   
      write(6,24) 'RESGHCUS     =  ', RESGHCUS   
      write(6,24) 'IMSGHCUS     =  ', IMSGHCUS   
      write(6,24) 'REPGHCUS     =  ', REPGHCUS   
      write(6,24) 'IMPGHCUS     =  ', IMPGHCUS   
      write(6,24) 'RESGHCUB     =  ', RESGHCUB   
      write(6,24) 'IMSGHCUB     =  ', IMSGHCUB   
      write(6,24) 'REPGHCUB     =  ', REPGHCUB   
      write(6,24) 'IMPGHCUB     =  ', IMPGHCUB   
      write(6,24) 'RESGHCCD     =  ', RESGHCCD   
      write(6,24) 'IMSGHCCD     =  ', IMSGHCCD   
      write(6,24) 'REPGHCCD     =  ', REPGHCCD   
      write(6,24) 'IMPGHCCD     =  ', IMPGHCCD   
      write(6,24) 'RESGHCCS     =  ', RESGHCCS   
      write(6,24) 'IMSGHCCS     =  ', IMSGHCCS   
      write(6,24) 'REPGHCCS     =  ', REPGHCCS   
      write(6,24) 'IMPGHCCS     =  ', IMPGHCCS   
      write(6,24) 'RESGHCCB     =  ', RESGHCCB   
      write(6,24) 'IMSGHCCB     =  ', IMSGHCCB   
      write(6,24) 'REPGHCCB     =  ', REPGHCCB   
      write(6,24) 'IMPGHCCB     =  ', IMPGHCCB   
      write(6,24) 'RESGHCTD     =  ', RESGHCTD   
      write(6,24) 'IMSGHCTD     =  ', IMSGHCTD   
      write(6,24) 'REPGHCTD     =  ', REPGHCTD   
      write(6,24) 'IMPGHCTD     =  ', IMPGHCTD   
      write(6,24) 'RESGHCTS     =  ', RESGHCTS   
      write(6,24) 'IMSGHCTS     =  ', IMSGHCTS   
      write(6,24) 'REPGHCTS     =  ', REPGHCTS   
      write(6,24) 'IMPGHCTS     =  ', IMPGHCTS   
      write(6,24) 'RESGHCTB     =  ', RESGHCTB   
      write(6,24) 'IMSGHCTB     =  ', IMSGHCTB   
      write(6,24) 'REPGHCTB     =  ', REPGHCTB   
      write(6,24) 'IMPGHCTB     =  ', IMPGHCTB   
      write(6,24) 'RESGHCVEE    =  ', RESGHCVEE  
      write(6,24) 'IMSGHCVEE    =  ', IMSGHCVEE  
      write(6,24) 'REPGHCVEE    =  ', REPGHCVEE  
      write(6,24) 'IMPGHCVEE    =  ', IMPGHCVEE  
      write(6,24) 'RESGHCVMMU   =  ', RESGHCVMMU 
      write(6,24) 'IMSGHCVMMU   =  ', IMSGHCVMMU 
      write(6,24) 'REPGHCVMMU   =  ', REPGHCVMMU 
      write(6,24) 'IMPGHCVMMU   =  ', IMPGHCVMMU 
      write(6,24) 'RESGHCVTTA   =  ', RESGHCVTTA 
      write(6,24) 'IMSGHCVTTA   =  ', IMSGHCVTTA 
      write(6,24) 'REPGHCVTTA   =  ', REPGHCVTTA 
      write(6,24) 'IMPGHCVTTA   =  ', IMPGHCVTTA 
      write(6,24) 'REGWWH1       =  ', REGWWH1     
      write(6,24) 'IMGWWH1       =  ', IMGWWH1     
      write(6,24) 'REGWWH2       =  ', REGWWH2     
      write(6,24) 'IMGWWH2       =  ', IMGWWH2     
      write(6,24) 'REGZZH1       =  ', REGZZH1     
      write(6,24) 'IMGZZH1       =  ', IMGZZH1     
      write(6,24) 'REGZZH2       =  ', REGZZH2     
      write(6,24) 'IMGZZH2       =  ', IMGZZH2     
      write(6,24) 'REGAHCHC      =  ', REGAHCHC    
      write(6,24) 'IMGAHCHC      =  ', IMGAHCHC    
      write(6,24) 'REGZH1H2      =  ', REGZH1H2    
      write(6,24) 'IMGZH1H2      =  ', IMGZH1H2    
      write(6,24) 'REGZH1H3      =  ', REGZH1H3    
      write(6,24) 'IMGZH1H3      =  ', IMGZH1H3    
      write(6,24) 'REGZH2H3      =  ', REGZH2H3    
      write(6,24) 'IMGZH2H3      =  ', IMGZH2H3    
      write(6,24) 'REGZHCHC      =  ', REGZHCHC    
      write(6,24) 'IMGZHCHC      =  ', IMGZHCHC    
      write(6,24) 'REGWPHCH1     =  ', REGWPHCH1   
      write(6,24) 'IMGWPHCH1     =  ', IMGWPHCH1   
      write(6,24) 'REGWPHCH2     =  ', REGWPHCH2   
      write(6,24) 'IMGWPHCH2     =  ', IMGWPHCH2   
      write(6,24) 'REGWPHCH3     =  ', REGWPHCH3   
      write(6,24) 'IMGWPHCH3     =  ', IMGWPHCH3   
      write(6,24) 'REGH1H1H1     =  ', REGH1H1H1   
      write(6,24) 'IMGH1H1H1     =  ', IMGH1H1H1   
      write(6,24) 'REGH1H1H2     =  ', REGH1H1H2   
      write(6,24) 'IMGH1H1H2     =  ', IMGH1H1H2   
      write(6,24) 'REGH1H2H2     =  ', REGH1H2H2   
      write(6,24) 'IMGH1H2H2     =  ', IMGH1H2H2   
      write(6,24) 'REGH1H3H3     =  ', REGH1H3H3   
      write(6,24) 'IMGH1H3H3     =  ', IMGH1H3H3   
      write(6,24) 'REGH2H2H2     =  ', REGH2H2H2   
      write(6,24) 'IMGH2H2H2     =  ', IMGH2H2H2   
      write(6,24) 'REGH2H3H3     =  ', REGH2H3H3   
      write(6,24) 'IMGH2H3H3     =  ', IMGH2H3H3   
      write(6,24) 'REGH1HCHC     =  ', REGH1HCHC   
      write(6,24) 'IMGH1HCHC     =  ', IMGH1HCHC   
      write(6,24) 'REGH2HCHC     =  ', REGH2HCHC   
      write(6,24) 'IMGH2HCHC     =  ', IMGH2HCHC   
      write(6,24) 'REGAAHCHC     =  ', REGAAHCHC   
      write(6,24) 'IMGAAHCHC     =  ', IMGAAHCHC   
      write(6,24) 'REGAZHCHC     =  ', REGAZHCHC   
      write(6,24) 'IMGAZHCHC     =  ', IMGAZHCHC   
      write(6,24) 'REGAWPHCH1    =  ', REGAWPHCH1  
      write(6,24) 'IMGAWPHCH1    =  ', IMGAWPHCH1  
      write(6,24) 'REGAWPHCH2    =  ', REGAWPHCH2  
      write(6,24) 'IMGAWPHCH2    =  ', IMGAWPHCH2  
      write(6,24) 'REGAWPHCH3    =  ', REGAWPHCH3  
      write(6,24) 'IMGAWPHCH3    =  ', IMGAWPHCH3  
      write(6,24) 'REGZZH1H1     =  ', REGZZH1H1   
      write(6,24) 'IMGZZH1H1     =  ', IMGZZH1H1   
      write(6,24) 'REGZZH2H2     =  ', REGZZH2H2   
      write(6,24) 'IMGZZH2H2     =  ', IMGZZH2H2   
      write(6,24) 'REGZZH3H3     =  ', REGZZH3H3   
      write(6,24) 'IMGZZH3H3     =  ', IMGZZH3H3   
      write(6,24) 'REGZZHCHC     =  ', REGZZHCHC   
      write(6,24) 'IMGZZHCHC     =  ', IMGZZHCHC   
      write(6,24) 'REGZWPHCH1    =  ', REGZWPHCH1  
      write(6,24) 'IMGZWPHCH1    =  ', IMGZWPHCH1  
      write(6,24) 'REGZWPHCH2    =  ', REGZWPHCH2  
      write(6,24) 'IMGZWPHCH2    =  ', IMGZWPHCH2  
      write(6,24) 'REGZWPHCH3    =  ', REGZWPHCH3  
      write(6,24) 'IMGZWPHCH3    =  ', IMGZWPHCH3  
      write(6,24) 'REGWWH1H1     =  ', REGWWH1H1   
      write(6,24) 'IMGWWH1H1     =  ', IMGWWH1H1   
      write(6,24) 'REGWWH2H2     =  ', REGWWH2H2   
      write(6,24) 'IMGWWH2H2     =  ', IMGWWH2H2   
      write(6,24) 'REGWWH3H3     =  ', REGWWH3H3   
      write(6,24) 'IMGWWH3H3     =  ', IMGWWH3H3   
      write(6,24) 'REGWWHCHC     =  ', REGWWHCHC   
      write(6,24) 'IMGWWHCHC     =  ', IMGWWHCHC   
      write(6,*)

      return
      end

