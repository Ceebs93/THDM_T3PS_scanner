!******************************************************
program HBwithLHClikelihood
!
! In this example we evaluate the likelihood from the CMS Higgs search with tautau
! final states (arXiv:1408.3316) for the mhmax scenario. The input is provided from
! a datafile created with SusHi obtained here:
! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMNeutral
!
! Note, that we deactivate the latest ATLAS and CMS tau tau 95% CL limits (the latter
! is reconstructed from the provided likelihood information internally in HiggsBounds)
! from the standard HiggsBounds procedure, since we want to use the likelihood instead.
! This is optional.
!
! (TS 27/05/2015)
!******************************************************
 use theory_XS_SM_functions
 use theory_BRfunctions
 use channels, only : HiggsBounds_deactivate_analyses,HiggsBounds_activate_all_analyses
 use output, only : createKey
 implicit none

 integer :: nHzero,nHplus,ndf, i, j, k, ii, jj, CP_value(3)
 double precision :: Mh(3),GammaTotal(3),CS_lep_hjZ_ratio(3),   &
     &          CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),    &   
     &          CS_lep_hjhi_ratio_nHbynH(3),                   &
     &          CS_tev_hj_ratio(3) ,CS_tev_hjb_ratio(3),    		&
     &          CS_tev_hjW_ratio(3),CS_tev_hjZ_ratio(3),    		&
     &          CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),   		&
     &          CS_lhc7_hj_ratio(3) ,CS_lhc7_hjb_ratio(3),    	&
     &          CS_lhc7_hjW_ratio(3),CS_lhc7_hjZ_ratio(3),    	&
     &          CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),   	&
     &          CS_lhc8_hj_ratio(3) ,CS_lhc8_hjb_ratio(3),    	&
     &          CS_lhc8_hjW_ratio(3),CS_lhc8_hjZ_ratio(3),    	&
     &          CS_lhc8_vbf_ratio(3),CS_lhc8_tthj_ratio(3),   	&
     &          BR_hjss(3),BR_hjcc(3),                            &
     &          BR_hjbb(3),                                    &
     &          BR_hjmumu(3),                                  &
     &          BR_hjtautau(3),                                &
     &          BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),BR_hjgaga(3),         &
     &          BR_hjgg(3), BR_hjinvisible(3),                    &
     &          BR_hjhihi_nHbynH(3,3)
 double precision :: ggH(3), bbH(3), tb, llh_h(3)                                 
 character(len=100)::filename_in, filename_out
 integer :: HBresult,chan,ncombined,HBresult_all,chan_all,ncombined_all
 double precision :: obsratio,obsratio_all
 integer :: error, status, n,  Hindex2, nc2,cbin
 double precision :: M_av, llh, mu, M_av2, llh2, mass, combllh,  M_av_exp, llh_exp
 integer ::  Hindex, nc, Hindex_exp, nc_exp
 
 nHzero = 3
 nHplus = 0
 
 call initialize_HiggsBounds(nHzero, nHplus, 'onlyH')
 
! Optionally, deactivate the 95%CL limit extraction for CMS MSSM h/H/A->tautau, since
! we want to use the likelihood information instead, as well as all relevant previous
! results from non-standard Higgs to tau tau searches from ATLAS/CMS:

! call HiggsBounds_deactivate_analyses((/3316, 2014049, 20140492/))

! Note: Deactivation of analyses can be changed before every HiggsBounds run 
! (as currently done, see below). If all analyses need to be activated again, just call
! call HiggsBounds_activate_all_analyses
 
 
 filename_in ="../example_data/mhmax_SusHi.dat"
 filename_out = "mhmax_HBwithLHClikelihood.dat"
 
 call system('rm -f mhmax_HBwithLHClikelihood.dat')
 
 open(432,file=trim(adjustl(filename_in)),action='read',status='old',iostat=status)
 open(433,file=trim(adjustl(filename_out)),action='write',status='new') 
  
  if(status.ne.0) then
   write(*,*) 'Bad status', status, 'with the following file:'
   write(*,*) trim(adjustl(filename_in))
   stop
  endif 
  n = 0
   do
! Read in the relevant cross section / BR predictions from the data grid.   
    read(432,*, iostat = error) Mh, tb, ggH, bbH, BR_hjtautau
    if (error == -1) exit
    n = n + 1
!----
! QUICK STOP FOR TESTING
!    if (n.le.9999) cycle
!    if (n.ge.10001) exit
!----    
    if(mod(n,100).eq.0) write(*,*) "number of processed points: ",n, " MA,TB = ", Mh(3),tb

    CP_value(1) = 1
    CP_value(2) = 1 
    CP_value(3) = -1 
    
      do i=1,3
       GammaTotal(i)=BRSM_GammaTot(Mh(i))
! Normalize the cross sections to the SM predictions
! (using HiggsBounds SM cross section routines)       
 	   CS_lhc8_hj_ratio(i)=ggH(i)/XS_lhc8_gg_H_SM(Mh(i))
	   CS_lhc8_hjb_ratio(i)=bbH(i)/XS_lhc8_bb_H_SM(Mh(i))
      enddo
! Set the other predictions to zero here.       
	  CS_lep_hjZ_ratio=0.0D0
	  CS_lep_bbhj_ratio=0.0D0
	  CS_lep_tautauhj_ratio=0.0D0
	  CS_lep_hjhi_ratio_nHbynH=0.0D0
	  CS_tev_hj_ratio=0.0D0
	  CS_tev_hjb_ratio=0.0D0
	  CS_tev_hjW_ratio=0.0D0
	  CS_tev_hjZ_ratio=0.0D0
	  CS_tev_vbf_ratio=0.0D0
	  CS_tev_tthj_ratio=0.0D0
 	  CS_lhc7_hj_ratio=0.0D0
 	  CS_lhc7_hjb_ratio=0.0D0
      CS_lhc7_hjW_ratio=0.0D0
      CS_lhc7_hjZ_ratio=0.0D0
	  CS_lhc7_vbf_ratio=0.0D0
	  CS_lhc7_tthj_ratio=0.0D0
	  CS_lhc8_hjW_ratio=0.0D0
	  CS_lhc8_hjZ_ratio=0.0D0
	  CS_lhc8_vbf_ratio=0.0D0
	  CS_lhc8_tthj_ratio=0.0D0
	  BR_hjss=0.0D0
	  BR_hjcc=0.0D0
	  BR_hjbb=0.0D0
	  BR_hjmumu=0.0D0
      BR_hjWW=0.0D0
      BR_hjZZ=0.0D0
      BR_hjZga=0.0D0
      BR_hjgaga=0.0D0
	  BR_hjgg=0.0D0
	  BR_hjinvisible=0.0D0
      BR_hjhihi_nHbynH=0.0D0              	  


	  call HiggsBounds_neutral_input_hadr(Mh,GammaTotal,CP_value,      &
     &          CS_lep_hjZ_ratio,                           &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
     &          CS_lep_hjhi_ratio_nHbynH,                   &
     &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    		&
     &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    		&
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   		&
     &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    	&
     &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    	&
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   	&
     &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    	&
     &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    	&
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   	&
     &          BR_hjss,BR_hjcc,                            &
     &          BR_hjbb,                                    &
     &          BR_hjmumu,                                  &
     &          BR_hjtautau,                                &
     &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
     &          BR_hjgg, BR_hjinvisible,                    &
     &          BR_hjhihi_nHbynH                            )


! Activate all analyses (in case some of them have been deactivated before)
      call HiggsBounds_activate_all_analyses

! Run the standard HiggsBounds routine considering all analyses
      call run_HiggsBounds( HBresult_all,chan_all, obsratio_all, ncombined_all )

	  call HiggsBounds_neutral_input_hadr(Mh,GammaTotal,CP_value,      &
     &          CS_lep_hjZ_ratio,                           &
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
     &          CS_lep_hjhi_ratio_nHbynH,                   &
     &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    		&
     &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    		&
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   		&
     &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    	&
     &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    	&
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   	&
     &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    	&
     &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    	&
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   	&
     &          BR_hjss,BR_hjcc,                            &
     &          BR_hjbb,                                    &
     &          BR_hjmumu,                                  &
     &          BR_hjtautau,                                &
     &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
     &          BR_hjgg, BR_hjinvisible,                    &
     &          BR_hjhihi_nHbynH                            )


! Deactivate present CMS and ATLAS searches for non-standard Higgs to tautau
      call HiggsBounds_deactivate_analyses((/14029, 2014049, 20140492/))! ,) 
! Standard HiggsBounds run (gives 95% CL limit):
      call run_HiggsBounds( HBresult,chan, obsratio, ncombined )

! Obtain exclusion-likelihood from CMS MSSM H->tautau search (ID 3316):
!
!  The arguments of the following subroutines mean the following:
!	3316 -> Analysis-ID of the CMS h/H/A->tautau search [int input]
!   Hindex -> Index of the Higgs boson that was selected as most sensitive [int output]
!   M_av -> mass position where limit is extracted (a signal strength average in case
!                                                   of combined Higgs bosons) [dbl output]
!   nc -> number of combined Higgs bosons [int output]
!   cbin -> binary code of the combined Higgs bosons [int output]
!   llh -> -2ln L value [dbl output]
!   obspred -> 'obs' or 'pred' to chose whether the observed or expected likelihood should be
!              extracted. [char input]
!
! Get expected/predicted likelihood 
      call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh_exp, 'pred')
! Get observed likelihood       
      call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh, 'obs')
      
      write(433,*) n, Mh, tb, HBresult, chan, obsratio, ncombined, llh, &
      &            HBresult_all,chan_all, obsratio_all, ncombined_all, Hindex,&
      &            M_av, nc, cbin, llh_exp

   enddo

!--------------------------------------------------------
! Example for the auxiliary chi^2 functions (run on the last point of the scan)
!--------------------------------------------------------
! Get observed likelihood for lightest Higgs

      write(*,*) "#------------------------------------------------------------------#"            
      call HiggsBounds_get_likelihood_for_Higgs(14029, 0, 1, nc, cbin, M_av, llh, 'obs')
      write(*,*) "The observed likelihood value for the light Higgs h is ",llh
      write(*,*) "The binary code and number of Higgses of the formed combination is ", cbin, nc
      write(*,*) "The likelihood has been evaluated at an average mass value of ", M_av
      write(*,*) "#------------------------------------------------------------------#"            
      
! Get observed likelihood for a subset of available Higgs bosons
! (here, e.g., exclude h and H from possible combination)

      call HiggsBounds_get_likelihood_for_comb(14029, 3, Hindex, nc, cbin, M_av, llh, 'obs')      
      write(*,*) "The observed likelihood value is ",llh
      write(*,*) "The binary code, number of Higgses and Higgs index of the formed combination is ", cbin, nc, Hindex
      write(*,*) "The likelihood has been evaluated at an average mass value of ", M_av
      write(*,*) "#------------------------------------------------------------------#"                 
!--------------------------------------------------------

! Write out the key with the used analyses (indicating possibly deactivated analyses)   
   call createKey("HB_with_deactivated_analyses_")
   close(432)
   close(433)
   
end program HBwithLHClikelihood