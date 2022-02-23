!--------------------------------------------------------------------
	program HBandHSwithFH
!--------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 06/03/2013).
!
! The program runs FeynHiggs (FH version > 2.9.4 required) to evaluate
! the partonic input of HiggsBounds/HiggsSignals for a example MSSM point.
! Then, we run HiggsBounds and HiggsSignals.
!
! Note: You have to set the correct FeynHiggs path in the configure file.
!       Then, this example program can be compiled by running
!
!       make HBandHSwithFH
!
!       in the HiggsSignals main directory.
!--------------------------------------------------------------------
	use usefulbits_hs, only : HSres
	use STXS, only : print_STXS
	implicit none

		integer error

c used by FHHiggsCorr
		double precision MHiggs(4)
		double complex SAeff, DeltaSAeff, UHiggs(3,3), ZHiggs(3,3)

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
c This is for FH 2.14.x
	double precision MSf(2,5,3), MASf(6,5), MCha(2), MNeu(4)
	double complex USf(2,2,5,3), UASf(6,6,5)
	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex DeltaMB
	double precision MGl, MHp
	double precision MHtree(4), SAtree
	double precision AlfasMT

Cdouble precision MHiggs(4)
Cdouble complex SAeff, DeltaSAeff, UHiggs(3,3), ZHiggs(3,3)


c used by FHRetrieveSMPara:
		double precision invAlfa0,invAlfaMZ, AlfasMZ, GF
		double precision ME, MU, MD, MM, MC, MS, ML, MB
		double precision MW, MZ, GammaW, GammaZ
		double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
		double precision RMUE, RAt, RM2

c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
        parameter (nHiggsneut = 3)
        parameter (nHiggsplus = 1)
       	character(LEN=5) whichanalyses

c used by HiggsBounds_neutral_input_part
        double precision Mh(3),GammaTotal_hj(3)
        integer CP_value(3)

        double precision CS_hj_ratio(4,3),
     &    CS_gg_hj_ratio(4,3),CS_bb_hj_ratio(4,3),
     &    CS_hjW_ratio(4,3),CS_hjZ_ratio(4,3),
     &    CS_vbf_ratio(4,3),CS_tthj_ratio(4,3),
     &    CS_thj_tchan_ratio(4,3),CS_thj_schan_ratio(4,3),
     &    CS_qq_hjZ_ratio(4,3),CS_gg_hjZ_ratio(4,3),
     &    CS_hjhi(4,3,3)

         double precision CS_ee_hjZ_ratio(3),CS_ee_bbhj_ratio(3),
     &    CS_ee_tautauhj_ratio(3),CS_ee_hjhi_ratio(3,3)

        double precision BR_hjss(3),BR_hjcc(3),
     &   BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),
     &   BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),
     &   BR_hjgaga(3),BR_hjgg(3),BR_hjtt(3),
     &   BR_hjinvisible(3),BR_hkhjhi(3,3,3),BR_hjhiZ(3,3),
     &   BR_hjemu(3),BR_hjetau(3),BR_hjmutau(3),BR_hjHpiW(3,1)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &   CS_ee_HpjHmj_ratio(1),
     &   BR_tWpb(1),BR_tHpjb(1),
     &   BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1),
     &   BR_Hpjtb(1),BR_HpjWZ(1),BR_HpjhiW(1,3),CS_Hpjtb(4)


c used in HiggsBounds
        double precision dmhneut_hb(nHiggsneut)
        double precision dmhch_hb(nHiggsplus)
c used in HiggsSignals (can be different)
        double precision dmhneut_hs(nHiggsneut)

        double precision obsratio_tmp,predratio_tmp
        integer HBresult_tmp,chan_tmp,ncombined_tmp

c return values of run_HiggsBounds
        integer HBresult,chan,ncombined
        double precision obsratio


c run options in HiggsSignals
        integer pdf, output_level

		double precision :: dmhneut(3)
		double precision :: dmhch(1)

c HiggsSignals output
        double precision Chisq_mu, Chisq_mh, Chisq, Pvalue
        integer nobs
      double precision Chisq_LHCRun1_mu,Chisq_LHCRun1_mh
      double precision Chisq_LHCRun1, Pvalue_LHCRun1
      integer nobs_LHCRun1

      double precision Chisq_STXS_rates, Chisq_STXS_mh
      double precision Chisq_STXS, Pvalue_STXS
      integer nobs_STXS

      integer ntotal, npeakmu, npeakmh, nmpred, nanalyses
c misc:
        integer i,j,k,as,t, collider
        double precision norm,CW2,Pi
        double precision
     &   g2hjbb(3),g2hjWW(3),g2hjZZ(3),
     &   g2hjgg(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass



        Pi = 3.1415926535897932384626433832795029D0

!--------------------------------------- Run Options ------------------------------------!
c The string 'whichanalyses' determines which subset of experimental
c results are used. In this example, we've used the option 'LandH',
c which instructs HiggsBounds to use tables of results
c from LEP, Tevatron and LHC (i.e. the full set of
c results supplied with HiggsBounds).
        whichanalyses='LandH'
c Choose a Higgs mass parametrization for HiggsSignals
c (1: box, 2:gaussian, 3:box+gaussian)
        pdf = 2
c Choose the output level of HiggsSignals
c (0: silent, 1: screen output, 2: even more output,...)
        output_level = 0


		dmhneut = (/3.0D0, 3.0D0, 0.0D0/)
		dmhch   = (/0.D0/)

!---------------------------- HiggsBounds and HiggsSignals ------------------------------!
c
c If you would like to perform scans over variables, the subroutines
c initialize_HiggsBounds, initialize_HiggsSignals (and finish_HiggsBounds,
c finish HiggsSignals) should be called outside the do-loops in order to save time.
c
!---- Initialize HiggsBounds and specify the dataset it should use 				     ----!
        call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
		call initialize_HiggsSignals(nHiggsneut,nHiggsplus,"LHC13")

		call HiggsBounds_set_mass_uncertainties(dmhneut,dmhch)

		call HiggsSignals_neutral_input_MassUncertainty(dmhneut)

!------------------------------ HiggsSignals options ------------------------------------!
!---- Set the output level                                                      	 ----!
		call setup_output_level(output_level)
!---- Set the Higgs mass parametrization                                     	 	 ----!
		call setup_pdf(pdf)

C	call setup_assignmentrange_massobservables(3.0D0)

!---------------------------------- Run FeynHiggs ---------------------------------------!
c the subroutines setFlags, setPara, setSLHA
c are also contained in this file
	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
!	call setPara
	call setPara_HH
c	call setSLHA("tmp.slha")

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
     &    invAlfa0, invAlfaMZ, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ, GammaW, GammaZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)

!--------------- Turn FeynHiggs results into HiggsBounds/Signals input ------------------!

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
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3))
         BR_hjtt(i)       = BR(H0FF(i,3,3,3))

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         BR_hjHpiW(i,1) = 0.0D0

         BR_hjemu(i) = BR(H0FF(i,2,1,2))
         BR_hjetau(i) = BR(H0FF(i,2,1,3))
         BR_hjmutau(i) = BR(H0FF(i,2,2,3))

         if(GammaSM(H0FF(i,4,3,3)).le.0.0D0)then
          g2hjbb(i)=0.0D0
         else
          g2hjbb(i)=Gamma(H0FF(i,4,3,3))
     &      /GammaSM(H0FF(i,4,3,3))
         endif

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

         CS_ee_bbhj_ratio(i)     = g2hjbb_s(i)+g2hjbb_p(i)
         CS_ee_tautauhj_ratio(i) = g2hjtautau_s(i)+g2hjtautau_p(i)

         g2hjWW(i)= dble(  Coupling(H0VV(i,4))
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4))
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3))
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3))
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0

         CS_ee_hjZ_ratio(i)        = g2hjZZ(i)

        enddo

        norm=GF*sqrt(2.0D0)*MZ**2.0D0

        do j=1,3
          do i=1,3
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm

            CS_ee_hjhi_ratio(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhiZ(j,i) = BR(H0HV(j,i))
            do k=1,3
             BR_hkhjhi(k,j,i) = BR(H0HH(k,j,i))
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


      do j=1,4
      select case(j)
       case(1)
        collider = 2
        sqrts = 2.0D0
       case(2)
        collider = 7
        sqrts = 7.0D0
       case(3)
        collider = 8
        sqrts = 8.0D0
       case(4)
        collider = 13
        sqrts = 13.0D0
      end select

      call FHHiggsProd(error, sqrts, prodxs)
      if( error .ne. 0 ) stop

      do i=1,3

      CS_gg_hj_ratio(j,i) = ggh(i)/gghSM(i)
      CS_bb_hj_ratio(j,i) = bbh(i)/bbhSM(i)
      CS_hj_ratio(j,i) = (ggh(i)+bbh(i))/
     & (gghSM(i)+bbhSM(i))
      CS_vbf_ratio(j,i) = qqh(i)/qqhSM(i)
      CS_hjW_ratio(j,i) = Wh(i)/WhSM(i)
      CS_hjZ_ratio(j,i) = Zh(i)/ZhSM(i)
      CS_qq_hjZ_ratio(j,i) = CS_hjZ_ratio(j,i) ! NOT CALCULATED EXPLICITLY BY FH
      CS_gg_hjZ_ratio(j,i) = CS_hjZ_ratio(j,i) ! NOT CALCULATED EXPLICITLY BY FH
      CS_tthj_ratio(j,i) = tth(i)/tthSM(i)
      CS_thj_tchan_ratio(j,i) = tth(i)/tthSM(i)  ! Use tth cross section here
      CS_thj_schan_ratio(j,i) = tth(i)/tthSM(i)  ! Use tth cross section here
      enddo

      CS_Hpjtb(j) = 2.0D0 * tHm2 / 1000.0D0 ! convert to pb.
      enddo

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c	    Charged Higgs input

        Mhplus(1)              = MHiggs(4)
        GammaTotal_Hpj(1)      = GammaTot(4)
        CS_ee_HpjHmj_ratio(1) = 1.0D0
        BR_tWpb(1)             = BR( tBF(1) )
        BR_tHpjb(1)            = BR( tBF(2) )
        BR_Hpjcs(1)            = BR( HpFF(2,2,2) )
        BR_Hpjcb(1)            = BR( HpFF(2,2,3) )
        BR_Hpjtaunu(1)         = BR( HpFF(1,3,3) )
        BR_Hpjtb(1)            = BR( HpFF(2,3,3) )
		BR_HpjWZ(1)            = 0.0D0
		do i=1,3
		BR_HpjhiW(1,i)         = BR( HpHV(i) )
		enddo

		print*,' '
		print*,'The Higgs boson masses are (h, H, A, H+):'
		print*,Mh,Mhplus
!--------------------------------- Run HiggsBounds --------------------------------------!

      call HiggsBounds_neutral_input_properties(Mh,
     &  GammaTotal_hj,CP_value)

      call HiggsBounds_neutral_input_SMBR(BR_hjss,BR_hjcc,
     & BR_hjbb,BR_hjtt,BR_hjmumu,BR_hjtautau,BR_hjWW,
     & BR_hjZZ,BR_hjZga,BR_hjgaga,BR_hjgg)

      call HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible,
     & BR_hkhjhi,BR_hjhiZ, BR_hjemu,BR_hjetau,BR_hjmutau,BR_hjHpiW)

      call HiggsBounds_neutral_input_LEP(CS_ee_hjZ_ratio,
     & CS_ee_bbhj_ratio, CS_ee_tautauhj_ratio,CS_ee_hjhi_ratio)

      do i=1,4
      select case(i)
       case(1)
       collider = 2
       case(2)
       collider = 7
       case(3)
       collider = 8
       case(4)
       collider = 13
      end select

      call HiggsBounds_neutral_input_hadr(collider,CS_hj_ratio(i,:),
     &  CS_gg_hj_ratio(i,:),CS_bb_hj_ratio(i,:),CS_hjW_ratio(i,:),
     &  CS_hjZ_ratio(i,:),CS_vbf_ratio(i,:),CS_tthj_ratio(i,:),
     &  CS_thj_tchan_ratio(i,:),CS_thj_schan_ratio(i,:),
     &  CS_qq_hjZ_ratio(i,:), CS_gg_hjZ_ratio(i,:), CS_hjhi(i,:,:))

! Only H+tb production can be calculated by FeynHiggs at the moment, set
! the other cross sections to zero (no experimental limits anyways, except on
! VBF->H+, which is zero in the MSSM.)

      call HiggsBounds_charged_input_hadr(collider, CS_Hpjtb(i), 0.0D0,
     & 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0)

      enddo

      call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj,
     & CS_ee_HpjHmj_ratio, BR_tWpb,BR_tHpjb,BR_Hpjcs,BR_Hpjcb,
     & BR_Hpjtaunu,BR_Hpjtb,BR_HpjWZ,BR_HpjhiW)


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
        call run_HiggsBounds( HBresult,chan,obsratio,ncombined)

      write(*,*) ' '
      write(*,*) '*************    HiggsBounds Results  **************'
      write(*,*) ' '
      write(*,*) 'Is this parameter point excluded by LEP, Tevatron'
      write(*,*) 'or LHC data?'
      write(*,*)  HBresult, ',  where'
      write(*,*) '               0 = yes, it is excluded'
      write(*,*) '               1 = no, it has not been excluded'
      write(*,*) '              -1 = invalid parameter set'
      write(*,*) ' '
      write(*,*) 'The process with the highest statistical sensitivity'
      write(*,*) 'is'
      write(*,*)  chan,'(see Key.dat)'
      write(*,*) 'This process has a theoretical rate vs. limit of'
      write(*,*)  obsratio
      write(*,*) ' '
      write(*,*) 'The number of Higgs which have contributed to the'
      write(*,*) 'theoretical rate of this process was'
      write(*,*)  ncombined
      write(*,*) ' '

      call HiggsBounds_get_most_sensitive_channels_per_Higgs(1,
     &     1,HBresult_tmp,chan_tmp,obsratio_tmp,predratio_tmp,
     &     ncombined_tmp)

      write(*,*) 'Detailed information about the most sensitive ',
     &'channels (per Higgs):'

      write(*,*) ' '
      write(*,*) ' Higgs boson HBresult channel     obsratio ',
     & '    predratio    ncombined'
        write(*,'(A10,2I10,2F15.5,1I5)') ' h ',HBresult_tmp,chan_tmp,
     &   obsratio_tmp,predratio_tmp, ncombined_tmp

      call HiggsBounds_get_most_sensitive_channels_per_Higgs(2,
     &     1,HBresult_tmp,chan_tmp,obsratio_tmp,predratio_tmp,
     &     ncombined_tmp)

        write(*,'(A10,2I10,2F15.5,1I5)') ' H ',HBresult_tmp,chan_tmp,
     &   obsratio_tmp, predratio_tmp, ncombined_tmp


      call HiggsBounds_get_most_sensitive_channels_per_Higgs(3,
     &     1,HBresult_tmp,chan_tmp,obsratio_tmp,predratio_tmp,
     &     ncombined_tmp)
        write(*,'(A10,2I10,2F15.5,1I5)') ' A ',HBresult_tmp,chan_tmp,
     &   obsratio_tmp, predratio_tmp, ncombined_tmp

      write(*,*) 'See HiggsBounds documentation for more information.'
      write(*,*) '****************************************************'
      write(*,*) ' '

!--------------------------------- Run HiggsSignals -------------------------------------!

c Set the theory mass uncertainties. This will affect the Chi^2 contribution
c from the mass measurements as well as the Higgs-to-peaks assignment in
c HiggsSignals .
C	dmhneut_hs = (/2.D0, 0.D0, 0.D0/)
C	call HiggsSignals_neutral_input_MassUncertainty(dmhneut_hs)
!---- Run HiggsSignals																 ----!
      call run_HiggsSignals(Chisq_mu,Chisq_mh,
     &  Chisq,nobs,Pvalue)

      call run_HiggsSignals_LHC_Run1_combination(Chisq_LHCRun1_mu,
     &  Chisq_LHCRun1_mh, Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1)

      call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh,
     &  Chisq_STXS, nobs_STXS, Pvalue_STXS)

      call complete_HS_results()

C     call print_STXS()

      write(*,*) ' '
      write(*,*) '*************  HiggsSignals Results  **************'
      write(*,*) ' '
      write(*,'(A40,1F10.5,A, 1I4)')
     &  "Chi^2 (signal strength obs.)/nobs = ",
     &  Chisq_mu, " / ", HSres(1)%nobs_peak_mu
      write(*,'(A40,1F10.5,A, 1I4)')
     &  "Chi^2 (STXS obs.)/nobs = ",
     &  Chisq_STXS, " / ", HSres(1)%nobs_STXS_rates
      write(*,'(A40,1F10.5,A, 1I4)')
     &  "Chi^2 (LHC-Run1 rate obs.)/nobs = ",
     &  Chisq_LHCRun1_mu, " / ", HSres(1)%nobs_LHCRun1_mu
      write(*,'(A40,1F10.5,A, 1I4)')
     &  "Chi^2 (mass obs.)/nobs = ",
     &  Chisq_mh + Chisq_LHCRun1_mh + Chisq_STXS_mh, " / ",
     &  HSres(1)%nobs_mh
      write(*,*) ' '
      write(*,*) 'The total chi^2 is ', HSres(1)%Chisq
      write(*,*) 'The number of observables is ',HSres(1)%nobs
C     write(*,*) ' '
C     write(*,*) 'In the case that the number of observables is equal'
C     write(*,*) 'to the number of degrees of freedom, this corresponds'
C     write(*,*) 'to a probability of'
C     write(*,*) Pvalue
      write(*,*) ' '
      write(*,*) 'See HiggsSignals documentation for more information.'
      write(*,*) '****************************************************'
      write(*,*) ' '

!---------------------- Finish HiggsBounds and HiggsSignals -----------------------------!

c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds
        call finish_HiggsSignals

	end


************************************************************************

	subroutine setFlags
	implicit none

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, loglevel, runningMT, botResum, tlCplxApprox

c   Using default (recommended) values of all FH flags
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 4)
	parameter (looplevel = 2)
	parameter (loglevel = 3)
	parameter (runningMT = 1)
	parameter (botResum = 1)
Cparameter (tlCplxApprox = 0) ! Use this if MA is input
	parameter (tlCplxApprox = 1) ! Use this if MHp is input
	integer error
c Since FH 2.14:
	call FHSetFlags(error, mssmpart, higgsmix, p2approx, looplevel,
     &   loglevel, runningMT, botResum, tlCplxApprox)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setPara
	implicit none

	double precision invAlfa0, invAlfaMZ, AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB, MZ, MW
	double precision GammaW, GammaZ
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

	parameter (invAlfa0 = -1)
	parameter (invAlfaMZ = -1)
	parameter (AlfasMZ = -1)
	parameter (GammaW = -1)
	parameter (GammaZ = -1)
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
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrhobar = -1)
	parameter (CKMetabar = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 173.2)
	parameter (TB = 10.)
	parameter (MA0 = 500)
	parameter (MHp = -1)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 2000)
	parameter (M3SL = MSusy)
	parameter (M2SL = M3SL)
	parameter (M1SL = M2SL)
	parameter (M3SE = MSusy)
	parameter (M2SE = M3SE)
	parameter (M1SE = M2SE)
	parameter (M3SQ = MSusy)
	parameter (M2SQ = M3SQ)
	parameter (M1SQ = M2SQ)
	parameter (M3SU = MSusy)
	parameter (M2SU = M3SU)
	parameter (M1SU = M2SU)
	parameter (M3SD = MSusy)
	parameter (M2SD = M3SD)
	parameter (M1SD = M2SD)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 2500)
	parameter (Ab = At)
	parameter (Atau = At)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 200)
	parameter (M_1 = 100)
	parameter (M_2 = 200)
	parameter (M_3 = 800)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	call FHSetSMPara(error,
     & invAlfa0, invAlfaMZ, AlfasMZ, GF,
     & ME, MU, MD, MM, MC, MS, ML, MB,
     & MW, MZ, GammaW, GammaZ,
     & CKMlambda, CKMA, CKMrhobar, CKMetabar)

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


	subroutine setPara_HH
	implicit none

	double precision invAlfa0, invAlfaMZ, AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB, MZ, MW
	double precision GammaW, GammaZ
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

	parameter (invAlfa0 = -1)
	parameter (invAlfaMZ = -1)
	parameter (AlfasMZ = -1)
	parameter (GammaW = -1)
	parameter (GammaZ = -1)
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
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrhobar = -1)
	parameter (CKMetabar = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 172.5)
	parameter (TB = 5.3)
	parameter (MA0 = -1)
	parameter (MHp = 170)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 2000)
	parameter (M3SL = MSusy)
	parameter (M2SL = M3SL)
	parameter (M1SL = M2SL)
	parameter (M3SE = MSusy)
	parameter (M2SE = M3SE)
	parameter (M1SE = M2SE)
	parameter (M3SQ = 710.)
	parameter (M2SQ = MSusy)
	parameter (M1SQ = M2SQ)
	parameter (M3SU = M3SQ)
	parameter (M2SU = MSusy)
	parameter (M1SU = M2SU)
	parameter (M3SD = M3SQ)
	parameter (M2SD = MSusy)
	parameter (M1SD = M2SD)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 461.5)
	parameter (Ab = At)
	parameter (Atau = At)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 5870)
	parameter (M_1 = 635)
	parameter (M_2 = 1000)
	parameter (M_3 = 2500)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	call FHSetSMPara(error,
     & invAlfa0, invAlfaMZ, AlfasMZ, GF,
     & ME, MU, MD, MM, MC, MS, ML, MB,
     & MW, MZ, GammaW, GammaZ,
     & CKMlambda, CKMA, CKMrhobar, CKMetabar)

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


!**************************************************************
	subroutine get_singleH_uncertainty(dCS,dggh,dbbh,g2hgg,g2hbb,mh)
	use theory_XS_SM_functions

	implicit none
	double precision dCS
 	double precision dggh, dbbh, g2hgg, g2hbb, mh
 	double precision vsmall

 	vsmall=1.0D-10

	if(g2hgg.le.vsmall.and.g2hbb.le.vsmall) then
		dCS = 0.0D0
	else
		dCS=(g2hgg*XS_lhc8_gg_H_SM(mh)*dggh+
     & 		g2hbb*XS_lhc8_bb_H_SM(mh)*dbbh)/
     &      (g2hgg*XS_lhc8_gg_H_SM(mh)+g2hbb*XS_lhc8_bb_H_SM(mh))
	endif

	end
!**************************************************************

