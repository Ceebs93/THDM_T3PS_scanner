! This file is part of SusHi.
!
! It includes the calls of HiggsBounds and HiggsSignals.
!

!--------------------------------------------------------------------
!	Code taken from the program HBandHSwithFH,
!       which is part of HiggsSignals (TS 06/03/2013).
!--------------------------------------------------------------------

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
c Set the HiggsSignals runmode (1: peak-centered, 2: mass-centered, 3: both) Chi^2 method(s)
      runmode = 1
c If the mass-centered chi^2 method is used, can specify the dm_theory treatment
c (1: mass variation, 2: smearing of mu-plot with mass pdf)
      mass_centered_method = 1
!---------------------------- HiggsBounds and HiggsSignals ------------------------------!
c
c If you would like to perform scans over variables, the subroutines
c initialize_HiggsBounds, initialize_HiggsSignals (and finish_HiggsBounds,
c finish HiggsSignals) should be called outside the do-loops in order to save time.
c
!---- Initialize HiggsBounds and specify the dataset it should use 				     ----!
      call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!

#ifdef HIGGSSIGNALS
      citations(34) = 1
      call initialize_HiggsSignals(nHiggsneut,nHiggsplus,
     &     "latestresults")
!------------------------------ HiggsSignals options ------------------------------------!
!---- Set the output level                                                      	 ----!
      call setup_output_level(output_level)
!---- Set the Higgs mass parametrization                                     	 	 ----!
      call setup_pdf(pdf)
!---- If the mass-centered chi^2 method is used, can specify the dm_theory treatment ----!
!     1: mass variation, 2: smearing of mu-plot with mass pdf
      call setup_mcmethod_dm_theory(mass_centered_method)
#endif

!---------------------------------- Run FeynHiggs ---------------------------------------!
c the subroutines setFlags, setPara, setSLHA
c are also contained in this file
c	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
c	call setPara
c	call setSLHA("tmp.slha")

c       Use this for FH < 2.9.5
c	call FHGetPara(error, nmfv, MASf, UASf,
c	call FHGetPara(error, nmfv, MSf, USf,MASf, UASf,
c     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
c     &    MHtree, SAtree)
c	if( error .ne. 0 ) stop

c	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
c	if( error .ne. 0 ) stop

c NOTE: we are setting uzint=uzext
      mfeff=1
      uzint=2
      uzext=2
      call FHSelectUZ(error, uzint, uzext, mfeff)
      if( error .ne. 0 ) stop

      fast=1
      call FHCouplings(error,
     &     couplings, couplingsms, gammas, gammasms, fast)
      if( error .ne. 0 ) stop

c We would like FH to calculate LHC cross sections
      sqrts_HB = 8.d0
      call FHHiggsProd(error, sqrts_HB, prodxs)
      if( error .ne. 0 ) stop 

      call FHRetrieveSMPara(error,
     &     invAlfa, AlfasMZ, GF,
     &     ME, MU, MD, MM, MC, MS, ML, MB,
     &     MW, MZ,
     &     CKMlambda, CKMA, CKMrhobar, CKMetabar)

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
         Mh_HB(i)=MHiggs(i) 
         GammaTotal_hj(i) = GammaTot(i)

         BR_hjss(i)       = BR(H0FF(i,4,2,2))
         BR_hjcc(i)       = BR(H0FF(i,3,2,2))
         BR_hjbb(i)       = BR(H0FF(i,4,3,3))
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3)) 

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         if(GammaSM(H0FF(i,4,3,3)).le.0.0D0)then
            g2hjbb(i)=0.0D0
         else
            g2hjbb(i)=Gamma(H0FF(i,4,3,3)) 
     &           /GammaSM(H0FF(i,4,3,3))
         endif

c Note that this is currently equivalent to
c         g2hjbb(i)= bbh(i)/bbhSM(i)
c         g2hjbb(i)= btagbh(i)/btagbhSM(i)
c as long as MH>80 GeV

         CS_bg_hjb_ratio(i) = g2hjbb(i)
         CS_bb_hj_ratio(i)  = g2hjbb(i)
   
         g2hjbb_s(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &        /RCouplingSM(H0FF(i,4,3,3))+ 
     &        LCoupling(H0FF(i,4,3,3))
     &        /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0
         g2hjbb_p(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &        /RCouplingSM(H0FF(i,4,3,3))- 
     &        LCoupling(H0FF(i,4,3,3))
     &        /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0

         g2hjtautau_s(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &        /RCouplingSM(H0FF(i,2,3,3))+ 
     &        LCoupling(H0FF(i,2,3,3))
     &        /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0
         g2hjtautau_p(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &        /RCouplingSM(H0FF(i,2,3,3))- 
     &        LCoupling(H0FF(i,2,3,3))
     &        /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0

         if(    g2hjbb_p(i).lt.1.0D-10)then
            CP_value(i) = 1
         elseif(g2hjbb_s(i).lt.1.0D-10)then
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
c Note that this is currently equivalent to
c         g2hjWW(i)= WhTev(i)/WhTevSM(i
c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
c as long as MH>80 GeV and uzint=uzext

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
c as long as MH>80 GeV and uzint=uzext
c It is also equivalent to g2hjWW(i)

         CS_lep_hjZ_ratio(i)    = g2hjZZ(i)

         CS_gg_hjZ_ratio(i)     = 0.0D0
         CS_dd_hjZ_ratio(i)     = g2hjZZ(i)
         CS_uu_hjZ_ratio(i)     = g2hjZZ(i)
         CS_ss_hjZ_ratio(i)     = g2hjZZ(i)
         CS_cc_hjZ_ratio(i)     = g2hjZZ(i)
         CS_bb_hjZ_ratio(i)     = g2hjZZ(i)

         CS_ud_hjWp_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWp_ratio(i)    = g2hjZZ(i)
         CS_ud_hjWm_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWm_ratio(i)    = g2hjZZ(i)

         CS_tev_vbf_ratio(i)    = g2hjZZ(i)
         CS_lhc7_vbf_ratio(i)   = g2hjZZ(i)
         CS_lhc8_vbf_ratio(i)   = g2hjZZ(i)


         if(tthSM(i).gt.0.0D0)then
            CS_tev_tthj_ratio(i)    = tth(i)/tthSM(i)  
         else
            CS_tev_tthj_ratio(i)    = 0.0D0
         endif

         CS_lhc7_tthj_ratio(i) = CS_tev_tthj_ratio(i)
         CS_lhc8_tthj_ratio(i) = CS_tev_tthj_ratio(i)

c tevatron gluon fusion XS is not calculated in FH is MH<90 geV
         if(Mh_HB(i).gt.90.0001D0)then      
            if(gghSM(i).gt.0.0D0)then
               CS_gg_hj_ratio(i) = ggh(i)/gghSM(i)
            else
               CS_gg_hj_ratio(i) = 0.0D0
            endif 
         else
            if(GammaSM(H0VV(i,5)).le.0.0D0)then
               CS_gg_hj_ratio(i)=0.0D0
            else
               CS_gg_hj_ratio(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
            endif
         endif

      enddo
  
      norm=GF*sqrt(2.0D0)*MZ**2.0D0

      do j=1,3
         do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &             dble(  Coupling(H0HV(j,i)) )**2.0D0
     &           + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &           )
     &           /norm

            CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhihi_nHbynH(j,i)=BR(H0HH(j,i,i))
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

c      print*,' '
c      print*,'The Higgs boson masses are (h, H, A, H+):'
c      print*,Mh_HB,Mhplus
!--------------------------------- Run HiggsBounds --------------------------------------!

c calls to HiggsBounds_neutral_input_part,HiggsBounds_charged_input,set_mass_uncertainties
c which give input to HiggsBounds

c        print*,'calling HiggsBounds_neutral_input_part'  
      call HiggsBounds_neutral_input_part(Mh_HB,GammaTotal_hj,CP_value, 
     &     CS_lep_hjZ_ratio,                            
     &     CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     
     &     CS_lep_hjhi_ratio_nHbynH,                    
     &     CS_gg_hj_ratio,CS_bb_hj_ratio,       
     &     CS_bg_hjb_ratio,                         
     &     CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,    
     &     CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,    
     &     CS_gg_hjZ_ratio,     
     &     CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     
     &     CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     
     &     CS_bb_hjZ_ratio,                         
     &     CS_tev_vbf_ratio,CS_tev_tthj_ratio,    
     &     CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    
     &     CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,    
     &     BR_hjss,BR_hjcc,                             
     &     BR_hjbb,BR_hjmumu,BR_hjtautau,               
     &     BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg,  
     &     BR_hjinvisible,BR_hjhihi_nHbynH              )
     
c        print*,'calling HiggsBounds_charged_input'
      call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, 
     &     CS_lep_HpjHmj_ratio,                        
     &     BR_tWpb,BR_tHpjb,                           
     &     BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)

c Set the theory mass uncertainties. HiggsBounds will vary the masses accordingly.
      dmhneut_hb = (/2.D0, 0.D0, 0.D0/)
      dmhch_hb   = (/0.D0/)
      call HiggsBounds_set_mass_uncertainties(dmhneut_hb,dmhch_hb)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
      call run_HiggsBounds( HB_result,HB_chan,HB_obsratio,HB_ncombined)

      if(verbose.gt.0) then
         print*,' '
         print*,'*************    HiggsBounds Results  **************'
         print*,' '
         print*,'Is this parameter point excluded by LEP, Tevatron'
         print*,'or LHC data?'
         print*, HB_result, ',  where'
         print*,'               0 = yes, it is excluded'
         print*,'               1 = no, it has not been excluded'
         print*,'              -1 = invalid parameter set'
         print*,' '
         print*,'The process with the highest statistical sensitivity'
         print*,'is'
         print*, HB_chan,'(see Key.dat)'
         print*,'This process has a theoretical rate vs. limit of'
         print*, HB_obsratio
         print*,' '
         print*,'The number of Higgs which have contributed to the'
         print*,'theoretical rate of this process was'
         print*, HB_ncombined
         print*,' '
         print*,'See HiggsBounds documentation for more information.'
         print*,'****************************************************'
         print*,' '
      endif

#ifdef HIGGSSIGNALS
!--------------------------------- Run HiggsSignals -------------------------------------!

c Set the theory mass uncertainties. This will affect the Chi^2 contribution
c from the mass measurements as well as the Higgs-to-peaks assignment in
c HiggsSignals .
      dmhneut_hs = (/2.D0, 0.D0, 0.D0/)
      call HiggsSignals_neutral_input_MassUncertainty(dmhneut_hs)
c Set the MSSM rate uncertainties:
!---dCS and dBR hold the model's rate uncertainties. Can be changed by user
!---with subroutine setup_rate_uncertainties. First, set to SM values.
      dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)  
      dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
c Now, we set relative uncertainties for the processes gg->H and bb->H
c and interpolate from the effective couplings/CS-ratios (of the lightest Higgs)
c the uncertainty of single Higgs production, dCS(1).
      dggh = 0.147D0
      dbbh = 0.200D0		
      call get_singleH_uncertainty(dCS(1),dggh,dbbh,CS_gg_hj_ratio(1),
     &     g2hjbb_s(1)+g2hjbb_p(1),Mh_HB(1))
      print*,"Interpolated singleH uncertainty (in %): ",dCS(1)
      call setup_rate_uncertainties(dCS, dBR)
!---- Run HiggsSignals																 ----!
      call run_HiggsSignals(runmode,HS_Chisq_mu,HS_Chisq_mh,
     &     HS_Chisq,HS_nobs,HS_Pvalue)
		
      if(verbose.gt.0) then
         print*,' '
         print*,'*************  HiggsSignals Results  **************'
         print*,' '
         print*,'The chi^2 contr. from signal strength observables is'
         print*,HS_Chisq_mu
         print*,'The chi^2 contr. from Higgs mass observables is'
         print*,HS_Chisq_mh
         print*,'The total chi^2 is'
         print*,HS_Chisq
         print*,'The number of observables is'
         print*,HS_nobs
         print*,' '
         print*,'In the case that the number of observables is equal'
         print*,'to the number of degrees of freedom, this corresponds'
         print*,'to a probability of'
         print*,HS_Pvalue
         print*,' '
         print*,'See HiggsSignals documentation for more information.'
         print*,'****************************************************'
         print*,' '
      endif

!---------------------- Finish HiggsBounds and HiggsSignals -----------------------------!

c deallocates arrays used by HiggsBounds:

      call finish_HiggsBounds
      call finish_HiggsSignals

#else

      call finish_HiggsBounds

#endif

c	end
