*********************************************************************
* HBwithFH
*
* Updated demo program for HiggsBounds 4 using MSSM input from
* FeynHiggs (FH version > 2.14.x recommended)
*
* Parameter point taken from Mh125(light stau) scenario
* defined in Bahl et al. (arXiv:1808.07542)
*********************************************************************

       program HBwithFH
        use usefulbits, only : analysislist
        implicit none

        integer error
        integer,parameter :: fileid=78
c used by FHHiggsCorr
        double precision MHiggs(4)
        double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSelectUZ:
        integer uzint, uzext, mfeff

c used by FHCouplings:
#include "FHCouplings.h"
        double complex couplings(ncouplings), couplingsms(ncouplingsms)
        double precision gammas(ngammas), gammasms(ngammasms)
        integer fast

c used by FHHiggsProd:
        double precision sqrts, prodxs(nprodxs)

c used by FHGetPara:
        integer nmfv
c This is for FH 2.10.x
        double precision MSf(2,5,3), MASf(6,5), MCha(2), MNeu(4)
        double complex USf(2,2,5,3), UASf(6,6,5)
        double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
        double complex DeltaMB
        double precision MGl
        double precision MHtree(4), SAtree
        double precision AlfasMT
c used by FHRetrieveSMPara:
        double precision invAlfa0, invAlfaMZ, AlfasMZ, GF
        double precision ME, MU, MD, MM, MC, MS, ML, MB
        double precision MW, MZ
        double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
        double precision GammaW, GammaZ
c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
       	character(LEN=5) whichanalyses

c used by HiggsBounds_neutral_input_hadr
C       character(LEN=5) collider_str
        double precision Mh(3),GammaTotal_hj(3)
        integer CP_value(3), collider_s
        double precision  CS_lep_hjZ_ratio(3),
     &   CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &   CS_lep_hjhi_ratio(3,3),CS_hj_ratio(4,3),
     &   CS_gg_hj_ratio(4,3),CS_bb_hj_ratio(4,3),
     &   CS_hjW_ratio(4,3),CS_hjZ_ratio(4,3),
     &   CS_vbf_ratio(4,3),CS_tthj_ratio(4,3),
     &   CS_thj_tchan_ratio(4,3),CS_thj_schan_ratio(4,3),
     &   CS_qq_hjZ_ratio(4,3),CS_gg_hjZ_ratio(4,3),
     &   CS_hjhi(4,3,3),
     &   CS_tWhj_ratio(4,3),
     &   BR_hjss(3),BR_hjcc(3),
     &   BR_hjbb(3),BR_hjtt(3),BR_hjmumu(3),BR_hjtautau(3),
     &   BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),
     &   BR_hjgaga(3),BR_hjgg(3),
     &   BR_hjinvisible(3),BR_hkhjhi(3,3,3),
     &   BR_hjhiZ(3,3),BR_hjemu(3),BR_hjetau(3),BR_hjmutau(3),
     &   BR_hjHpiW(3,1)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &   CS_lep_HpjHmj_ratio(1),
     &   BR_tWpb,BR_tHpjb(1),BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1),
     &   BR_Hpjtb(1),BR_HpjWZ(1),BR_HpjhiW(1,3),
     &   CS_Hpjtb(4,1), CS_Hpjbjet(4,1),CS_Hpjcb(4,1), CS_Hpjcjet(4,1),
     &   CS_Hpjjetjet(4,1), CS_HpjW(4,1),
     &   CS_HpjZ(4,1), CS_vbf_Hpj(4,1), CS_HpjHmj(4,1), CS_Hpjhi(4,1,3)

c used by run_HiggsBounds
        integer HBresult,chan,ncombined
        double precision obsratio

c misc:
        integer i,j,k,as,t,collider
        double precision norm,CW2,Pi
        double precision
     &   g2hjbb(3),g2hjWW(3),g2hjZZ(3),
     &   g2hjgg(3),g2hjhiZ(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass

        Pi = 3.1415926535897932384626433832795029D0

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c Set number of neutral and charged Higgs bosons in the MSSM:
        nHiggsneut=3
        nHiggsplus=1

c The string 'whichanalyses' determines which subset of experimental
c results are used. In this example, we've used the option 'LandH',
c which instructs HiggsBounds to use tables of results
c from LEP, Tevatron and LHC (i.e. the full set of
c results supplied with HiggsBounds).
       whichanalyses='LandH'
C       whichanalyses='list '
C       allocate(analysislist(2))
C       analysislist(1) = 170271
C       analysislist(2) = 170272

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before any of the other HiggsBounds subroutines.
        call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c calls to FH subroutines:

c the subroutines setFlags, setPara, setSLHA
c are also contained in this file
	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
	call setPara
c	call setSLHA("tmp.slha")

c       // User this line for FH < 2.9.5
c	call FHGetPara(error, nmfv, MASf, UASf,

	call FHGetPara(error, nmfv, MSf, USf,MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &    MHtree, SAtree, AlfasMT)
	if( error .ne. 0 ) stop

	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
	if( error .ne. 0 ) stop

c NOTE: we are setting uzint=uzext
C       mfeff=1
C       uzint=2
C       uzext=2
C       call FHSelectUZ(error, uzint, uzext, mfeff)
C    if( error .ne. 0 ) stop

        fast=1
        call FHCouplings(error,
     &    couplings, couplingsms, gammas, gammasms, fast)
        if( error .ne. 0 ) stop

        call FHRetrieveSMPara(error,
     &   invAlfa0, invAlfaMZ, AlfasMZ, GF,
     &   ME, MU, MD, MM, MC, MS, ML, MB,
     &   MW, MZ, GammaW, GammaZ,
     &   CKMlambda, CKMA, CKMrhobar, CKMetabar)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c Set variables needed by HiggsBounds (using results from FeynHiggs).
c See HiggsBounds documentation for definition of variables used
c as arguments to HiggsBounds_neutral_input_part and run_HiggsBounds
c and FeynHiggs documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine HiggsBounds_neutral_input_part
c rather than the subroutine HiggsBounds_neutral_input_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in FeynHiggs

        do i=1,3
         Mh(i)=MHiggs(i)
         GammaTotal_hj(i) = GammaTot(i)

         BR_hjss(i)       = BR(H0FF(i,4,2,2))
         BR_hjcc(i)       = BR(H0FF(i,3,2,2))
         BR_hjbb(i)       = BR(H0FF(i,4,3,3))
         BR_hjtt(i)       = BR(H0FF(i,3,3,3))
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3))
         BR_hjemu(i)      = BR(H0FF(i,2,1,2))
         BR_hjetau(i)     = BR(H0FF(i,2,1,3))
         BR_hjmutau(i)    = BR(H0FF(i,2,2,3))

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         BR_hjHpiW(i,1) = 0.0D0 ! NOT CALCULATED BY FH!

        do collider=1,4
c We would like FH to calculate LHC cross sections
         select case(collider)
          case(1)
           sqrts=2.0D0
          case(2)
           sqrts=7.0D0
          case(3)
           sqrts=8.0D0
          case(4)
           sqrts=13.0D0
         end select

         call FHHiggsProd(error, sqrts, prodxs)
         if( error .ne. 0 ) stop

         CS_gg_hj_ratio(collider,i)=ggh(i)/gghSM(i)
         CS_bb_hj_ratio(collider,i)=bbh(i)/bbhSM(i)
         CS_hj_ratio(collider,i)=(ggh(i)+bbh(i))/(gghSM(i)+bbhSM(i))
         CS_hjW_ratio(collider,i)=Wh(i)/WhSM(i)
         CS_hjZ_ratio(collider,i)=Zh(i)/ZhSM(i)
         CS_qq_hjZ_ratio(collider,i) = CS_hjZ_ratio(collider,i) ! NOT CALCULATED EXPLICITLY BY FH
         CS_gg_hjZ_ratio(collider,i) = CS_hjZ_ratio(collider,i) ! NOT CALCULATED EXPLICITLY BY FH
         CS_vbf_ratio(collider,i)=qqh(i)/qqhSM(i)
         CS_tthj_ratio(collider,i)=tth(i)/tthSM(i)
         CS_thj_tchan_ratio(collider,i)=0.0D0! NOT CALCULATED YET BY FH
         CS_thj_schan_ratio(collider,i)=0.0D0! NOT CALCULATED YET BY FH
         CS_tWhj_ratio(collider,i)=0.0D0
         do j=1,3
          CS_hjhi(collider,i,j) = 0.0D0 ! NOT CALCULATED YET BY FH
         enddo
		! Charged Higgs cross section (HB input requires the sum of process and its charge conjugate!)
         CS_Hpjtb(collider,1)=2.0D0 * tHm2/1000. ! n.b.: Rescale from fb to pb! tHm2 == sigma(gb->tH-)
         CS_Hpjcb(collider,1)=0.0D0  ! NOT CALCULATED YET BY FH
         CS_Hpjbjet(collider,1)=0.0D0  ! NOT CALCULATED YET BY FH
         CS_Hpjcjet(collider,1)=0.0D0  ! NOT CALCULATED YET BY FH
         CS_Hpjjetjet(collider,1)=0.0D0  ! NOT CALCULATED YET BY FH
         enddo

         BR_HpjhiW(1,i)=BR(HpHV(i))

         g2hjbb_s(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))+
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0
         g2hjbb_p(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))-
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0

         g2hjtautau_s(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))+
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0
         g2hjtautau_p(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))-
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0

         if(       g2hjbb_p(i).lt.1.0D-10)then
            CP_value(i) = 1
         elseif(   g2hjbb_s(i).lt.1.0D-10)then
            CP_value(i) = -1
         else
            CP_value(i) = 0
         endif

         CS_lep_bbhj_ratio(i)     = g2hjbb_s(i)+g2hjbb_p(i)
         CS_lep_tautauhj_ratio(i) = g2hjtautau_s(i)+g2hjtautau_p(i)

         g2hjWW(i)= dble(  Coupling(H0VV(i,4))
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4))
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3))
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3))
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0

         CS_lep_hjZ_ratio(i)        = g2hjZZ(i)

         enddo
         norm=GF*sqrt(2.0D0)*MZ**2.0D0

         do j=1,3
           do i=1,3
            g2hjhiZ(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm

            CS_lep_hjhi_ratio(j,i) = g2hjhiZ(j,i)

            BR_hjhiZ(j,i)=BR(H0HV(j,i))
            do k=1,3
             BR_hkhjhi(k,j,i)=BR(H0HH(k,j,i))
            enddo
           enddo
          enddo

c higgs->neutralino1 neutralino1 contributes the invisible Higgs decay width
c when neutralino1 or sneutrino is the LSP

        do i=1,3
         sneutrino_lspcandidate_number=0
         invisible_lsp=.True.

c first determine whether lightest sneutrino is lighter than the lightest neutralino
c
c sneutrino_lspcandidate_number=0 indicates that lightest neutralino is
c lighter than all the sneutrinos
         lspcandidate_mass=MNeu(1)
         do as=1,3
            if( MASf(as,1) .lt. lspcandidate_mass )then
               lspcandidate_mass=MASf(as,1)
               sneutrino_lspcandidate_number=as
            endif
         enddo

         if(     MCha(1) .lt. lspcandidate_mass )then
           invisible_lsp=.False.
         elseif( MGl     .lt. lspcandidate_mass )then
           invisible_lsp=.False.
         else
          do as=1,6
           do t=2,4
            if( MASf(as,t) .lt. lspcandidate_mass )then
               invisible_lsp=.False.
            endif
           enddo
          enddo
         endif

         if(invisible_lsp)then
          if(sneutrino_lspcandidate_number.eq.0)then
            BR_hjinvisible(i) = BR(H0NeuNeu(i,1,1))
          else
            BR_hjinvisible(i) = BR(H0SfSf(i,1,1,1,as))
          endif
         else
          BR_hjinvisible(i) = 0.0D0
         endif

        enddo

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c	    Charged Higgs input

        Mhplus(1)              = MHiggs(4)
        GammaTotal_Hpj(1)      = GammaTot(4)
        CS_lep_HpjHmj_ratio(1) = 1.0D0
        BR_tWpb                = BR( tBF(1) )
        BR_tHpjb(1)            = BR( tBF(2) )
        BR_Hpjcs(1)            = BR( HpFF(2,2,2) )
        BR_Hpjcb(1)            = BR( HpFF(2,2,3) )
        BR_Hpjtaunu(1)         = BR( HpFF(1,3,3) )
        BR_Hpjtb(1)            = BR( HpFF(2,3,3) )
        BR_HpjWZ(1)            = 0.0D0 ! NOT CALCULATED BY FH (0 in the MSSM at tree-level)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

       write(*,*) "Higgs masses: ", Mh

        call HiggsBounds_set_mass_uncertainties((/2.0D0,0.0D0,0.0D0/), 0.0D0)

C       print*,'calling HiggsBounds_neutral_input_properties'

        call HiggsBounds_neutral_input_properties(Mh, GammaTotal_hj, CP_value)

C       print*,'calling HiggsBounds_charged_input'

        call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj,
     &          CS_lep_HpjHmj_ratio,
     &          BR_tWpb,BR_tHpjb,
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu,BR_Hpjtb,
     &          BR_HpjWZ,BR_HpjhiW)

C       print*,'calling HiggsBounds_neutral_input_SMBR'

        call HiggsBounds_neutral_input_SMBR(BR_hjss,BR_hjcc,BR_hjbb,
     &                           BR_hjtt,BR_hjmumu,
     &                           BR_hjtautau,BR_hjWW,
     &                           BR_hjZZ,BR_hjZga,BR_hjgaga,
     &                           BR_hjgg)

C       print*,'calling HiggsBounds_neutral_input_nonSMBR'

        call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible,BR_hkhjhi,BR_hjhiZ,
     &                                 BR_hjemu,BR_hjetau,BR_hjmutau,BR_hjHpiW)

C       print*,'calling HiggsBounds_neutral_input_LEP'

        call HiggsBounds_neutral_input_LEP(CS_lep_hjZ_ratio,CS_lep_bbhj_ratio,
     &                                 CS_lep_tautauhj_ratio,CS_lep_hjhi_ratio)

C       print*,'calling HiggsBounds_neutral_input_hadr'
C       print*,'    and HiggsBounds_charged_input_hadr'

        do collider=1,4
         select case(collider)
          case(1)
           collider_s = 2
          case(2)
           collider_s = 7
          case(3)
           collider_s = 8
          case(4)
           collider_s = 13
         end select

        call HiggsBounds_neutral_input_hadr(collider_s,CS_hj_ratio(collider,:),
     &          CS_gg_hj_ratio(collider,:),CS_bb_hj_ratio(collider,:),
     &          CS_hjW_ratio(collider,:),CS_hjZ_ratio(collider,:),
     &          CS_vbf_ratio(collider,:),CS_tthj_ratio(collider,:),
     &          CS_thj_tchan_ratio(collider,:),CS_thj_schan_ratio(collider,:),
     &          CS_qq_hjZ_ratio(collider,:),CS_gg_hjZ_ratio(collider,:),
     &          CS_tWhj_ratio(collider,:),
     &          CS_hjhi(collider,:,:))

        call HiggsBounds_charged_input_hadr(collider_s,CS_Hpjtb(collider,:),
     &   CS_Hpjcb(collider,:), CS_Hpjbjet(collider,:), CS_Hpjcjet(collider,:),
     &   CS_Hpjjetjet(collider,:), CS_HpjW(collider,:),
     &   CS_HpjZ(collider,:), CS_vbf_Hpj(collider,:),
     &   CS_HpjHmj(collider,:), CS_Hpjhi(collider,:,:))

        enddo


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
C       print*,' running HiggsBounds...'
        call run_HiggsBounds( HBresult,chan, obsratio, ncombined )

        print*,' '
        print*,'*************    HiggsBounds Results  **************'
        print*,' '
        print*,'Is this parameter point excluded by LEP, Tevatron'
        print*,'or LHC data?'
        print*, HBresult, ',  where'
        print*,'               0 = yes, it is excluded'
        print*,'               1 = no, it has not been excluded'
        print*,'              -1 = invalid parameter set'
        print*,' '
        print*,'The process with the highest statistical sensitivity'
        print*,'is'
        print*, chan,'(see Key.dat)'
        print*,'This process has a theoretical rate vs. limit of'
        print*, obsratio
        print*,' '
        print*,'The number of Higgs which have contributed to the'
        print*,'theoretical rate of this process was'
        print*, ncombined
        print*,' '
        print*,'See HiggsBounds documentation for more information.'
        print*,'****************************************************'
        print*,' '

        open(fileid, file="HBwithFH-output.dat")
        write(fileid,'(3F10.4,I2,I5,F8.4,I2)')
     &      Mh,HBresult,chan,obsratio,ncombined
        close(fileid)


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds

	end


************************************************************************

	subroutine setFlags
	implicit none

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, loglevel, runningMT, botResum, tlCplxApprox

c   Using default (recommended) values of all FH flags
	parameter (mssmpart = 4)
Cparameter (fieldren = 0)
Cparameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 4)
	parameter (looplevel = 2)
	parameter (loglevel = 3)
	parameter (runningMT = 1)
	parameter (botResum = 1)
	parameter (tlCplxApprox = 0)

	integer error

Ccall FHSetFlags(error, mssmpart, fieldren, tanbren,
C    &    higgsmix, p2approx, looplevel, loglevel,
C    &    runningMT, botResum, tlCplxApprox)

	call FHSetFlags(error,mssmpart, higgsmix, p2approx,
     &   looplevel, loglevel, runningMT, botResum, tlCplxApprox)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setPara
	implicit none

	double precision invAlfa0, invAlfaMZ,  AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB, MZ, MW
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
	double precision GammaW, GammaZ
	parameter (invAlfa0 = -1)
	parameter (invAlfaMZ = -1)
	parameter (AlfasMZ = -1)
	parameter (GF = -1)
	parameter (ME = -1)
	parameter (MU = -1)
	parameter (MD = -1)
	parameter (MM = -1)
	parameter (MC = -1)
	parameter (MS = -1)
	parameter (ML = -1)
	parameter (MB = -1)
	parameter (MW = -1)
	parameter (MZ = -1)
	parameter (GammaW = -1)
	parameter (GammaZ = -1)
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrhobar = -1)
	parameter (CKMetabar = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 172.5)
	parameter (TB = 10.)
	parameter (MA0 = 1000)
Cparameter (TB = 1.25)
Cparameter (MA0 = 500)

	parameter (MHp = -1)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 2000)
	parameter (M3SL = 350)
	parameter (M2SL = MSusy)
	parameter (M1SL = MSusy)
	parameter (M3SE = M3SL)
	parameter (M2SE = MSusy)
	parameter (M1SE = MSusy)
	parameter (M3SQ = 1500)
	parameter (M2SQ = MSusy)
	parameter (M1SQ = MSusy)
	parameter (M3SU = M3SQ)
	parameter (M2SU = MSusy)
	parameter (M1SU = MSusy)
	parameter (M3SD = M3SQ)
	parameter (M2SD = MSusy)
	parameter (M1SD = MSusy)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 2900)
	parameter (Ab = At)
	parameter (Atau = 800)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 1000)
	parameter (M_1 = 180)
	parameter (M_2 = 300)
	parameter (M_3 = 2500)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	write(*,*) "Xt/M3SQ = ", At/M3SQ - MUE/(M3SQ*TB)

	call FHSetSMPara(error,
     &   invAlfa0, invAlfaMZ, AlfasMZ, GF,
     &   ME, MU, MD, MM, MC, MS, ML, MB,
     &   MW, MZ, GammaW, GammaZ,
     &   CKMlambda, CKMA, CKMrhobar, CKMetabar)

	if( error .ne. 0 ) stop

	call FHSetPara(error, scalefactor,
     &    MT, TB, MA0, MHp,
     &    M3SL, M3SE, M3SQ, M3SU, M3SD,
     &    M2SL, M2SE, M2SQ, M2SU, M2SD,
     &    M1SL, M1SE, M1SQ, M1SU, M1SD,
     &    MUE,
     &    Atau, At, Ab,
     &    Amu, Ac, As,
     &    Ae, Au, Ad,
     &    M_1, M_2, M_3,
     &    Qtau, Qt, Qb)
        if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setSLHA(filename)
	implicit none
	character*(*) filename

#include "SLHA.h"

	integer error
	double complex slhadata(nslhadata)

	call SLHARead(error, slhadata, filename, 1)
	if( error .ne. 0 ) stop

	call FHSetSLHA(error, slhadata)
	if( error .ne. 0 ) stop
	end



