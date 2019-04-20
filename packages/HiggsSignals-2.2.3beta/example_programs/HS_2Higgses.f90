!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 30/10/2014).
!
! In this example we scan over the masses of two Higgs bosons, for fixed values of
! universal signal strength scale factors (scale1, scale2) and theoretical mass un-
! certainties (dm(1), dm(2)). We loop over the three different pdf choices.
! The output (total chi^2) can be plotted by using the gnuplot script plot_2Higgses.gnu.
!--------------------------------------------------------------------------------------
program HS_2Higgses
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall
 use usefulbits_hs, only : HSres
 use pc_chisq, only : print_inverse_cov_mh_to_file, get_masschi2_from_separation
 implicit none

 integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
 double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
 double precision :: Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, Chisq_LHCRun1, Pvalue_LHCRun1
 double precision :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, Pvalue_STXS
 integer ::  nobs_LHCRun1,  nobs_STXS
 double precision :: SMGammaTotal(2)
double precision :: scale_bbh, scale_ggh, dggh, dbbh
 double precision :: Mh(2),GammaTotal(2),ghjss_s(2),ghjss_p(2),ghjcc_s(2),ghjcc_p(2), &
&                    ghjbb_s(2),ghjbb_p(2),ghjtt_s(2),ghjtt_p(2), &
&                    ghjmumu_s(2),ghjmumu_p(2),ghjtautau_s(2),ghjtautau_p(2), &
&                    ghjWW(2),ghjZZ(2),ghjZga(2),ghjgaga(2),ghjgg(2),	&
&                    ghjhiZ(2,2)
 character(len=100)::filename
 character(len=1)::pdfchar
 double precision :: dm(2)
 integer		  :: pdf
!-HiggsBounds internal functions to obtain SM branching ratios 
 double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau,&
 &                   SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg,SMGamma_h
 double precision :: scale1, scale2, csq_sep
 nHzero=2
 nHplus=0

!--Enter the Higgs mass and its theory uncertainty here: 
 Mh = (/ 125.0D0, 125.0D0 /)
 dm = (/ 0.5D0, 0.5D0 /)
! dm = (/ 0.0D0, 0.0D0 /)

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
 call initialize_HiggsSignals(nHzero,nHplus,"LHC13")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
 call setup_output_level(0) 
! call setup_anticorrelations_in_mh(1)
! call setup_assignmentrange_massobservables(2.0D0)
! Always normalize rate prediction w.r.t. to predicted Higgs mass
 call setup_rate_normalization(.False.,.False.)
!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
 call HiggsSignals_neutral_input_MassUncertainty(dm)
! Set number of free parameters (for p-value evaluation)
 call setup_nparam(2)

 do pdf = 2,2
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
  call setup_pdf(pdf) 
  write(pdfchar,'(I1)') pdf
  filename='results/2Higgses_pdf'//pdfchar//'.dat'
  open(21,file=trim(adjustl(filename)))
  write(21,*) '# Mh(1) Mh(2) dm(1) dm(2) scale(1) scale(2) chisq_mu chisq_mh chisq ndf Pvalue chisq_sep'
  write(21,*) '#--------------------------------------------------------------------------'

   scale1 = sqrt(0.75D0)
   scale2 = sqrt(1.0D0-scale1**2.0D0)

    ghjss_s(1) = scale1
    ghjss_s(2) = scale2
    ghjss_p=0.0d0
    ghjcc_s(1) = scale1
    ghjcc_s(2) = scale2    
    ghjcc_p=0.0d0
    ghjbb_s(1) = scale1
    ghjbb_s(2) = scale2    
    ghjbb_p=0.0d0
    ghjtt_s(1) = scale1
    ghjtt_s(2) = scale2    
    ghjtt_p=0.0d0         
    ghjmumu_s(1) = scale1
    ghjmumu_s(2) = scale2    
    ghjmumu_p=0.0d0  
    ghjtautau_s(1) = scale1
    ghjtautau_s(2) = scale2
    ghjtautau_p=0.0d0
    ghjWW(1) = scale1
    ghjWW(2) = scale2    
    ghjZZ(1) = scale1
    ghjZZ(2) = scale2    
    ghjZga(1) = scale1
    ghjZga(2) = scale2    
    ghjgg(1) = scale1
    ghjgg(2) = scale2    
!     ghjggZ(1) = scale1
!     ghjggZ(2) = scale2    
    ghjhiZ=0d0
    ghjgaga(1) = scale1
    ghjgaga(2) = scale2 

 do i=1,41!,16
  do j=1,41!16
   Mh(1) = 124D0 + (i-1)*0.05D0
   Mh(2) = 124D0 + (j-1)*0.05D0

!    Mh(1) = 123.5D0 + (i-1)*0.2D0
!    Mh(2) = 123.5D0 + (j-1)*0.2D0
             
   SMGammaTotal(1)=SMGamma_h(Mh(1))
   SMGammaTotal(2)=SMGamma_h(Mh(2))   
      
!----Calculate the new total decay width:
    GammaTotal(1) = SMGammaTotal(1)*(1 + &
	&  	            (ghjWW(1)**2.0D0 - 1)*SMBR_HWW(Mh(1))+(ghjZZ(1)**2.0D0 - 1)*SMBR_HZZ(Mh(1)) + &
	&               (ghjgg(1)**2.0D0 - 1)*SMBR_Hgg(Mh(1))+(ghjtt_s(1)**2.0D0 - 1)*SMBR_Htoptop(Mh(1))+ &
	&               (ghjbb_s(1)**2.0D0 - 1)*SMBR_Hbb(Mh(1))+(ghjtautau_s(1)**2.0D0 - 1)*SMBR_Htautau(Mh(1))+ &
	&               (ghjss_s(1)**2.0D0 - 1)*SMBR_Hss(Mh(1))+(ghjcc_s(1)**2.0D0 - 1)*SMBR_Hcc(Mh(1))+ &
	&               (ghjZga(1)**2.0D0 - 1)*SMBR_HZgam(Mh(1))+(ghjmumu_s(1)**2.0D0 - 1)*SMBR_Hmumu(Mh(1))+ &		
	&               (ghjgaga(1)**2.0D0 - 1)*SMBR_Hgamgam(Mh(1))	)

    GammaTotal(2) = SMGammaTotal(2)*(1 + &
	&  	            (ghjWW(2)**2.0D0 - 1)*SMBR_HWW(Mh(2))+(ghjZZ(2)**2.0D0 - 1)*SMBR_HZZ(Mh(2)) + &
	&               (ghjgg(2)**2.0D0 - 1)*SMBR_Hgg(Mh(2))+(ghjtt_s(2)**2.0D0 - 1)*SMBR_Htoptop(Mh(2))+ &
	&               (ghjbb_s(2)**2.0D0 - 1)*SMBR_Hbb(Mh(2))+(ghjtautau_s(2)**2.0D0 - 1)*SMBR_Htautau(Mh(2))+ &
	&               (ghjss_s(2)**2.0D0 - 1)*SMBR_Hss(Mh(2))+(ghjcc_s(2)**2.0D0 - 1)*SMBR_Hcc(Mh(2))+ &
	&               (ghjZga(2)**2.0D0 - 1)*SMBR_HZgam(Mh(2))+(ghjmumu_s(2)**2.0D0 - 1)*SMBR_Hmumu(Mh(2))+ &		
	&               (ghjgaga(2)**2.0D0 - 1)*SMBR_Hgamgam(Mh(2))	)

	call HiggsBounds_neutral_input_properties(Mh,GammaTotal,(/1, 1/))

	call HiggsBounds_neutral_input_effC(                     &  
     &          ghjss_s,ghjss_p,ghjcc_s,ghjcc_p,       &
     &          ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p,       &
     &          ghjmumu_s,ghjmumu_p,                   &
     &          ghjtautau_s,ghjtautau_p,               &
     &          ghjWW,ghjZZ,ghjZga,                    &
     &          ghjgaga,ghjgg,ghjhiZ)

    call run_HiggsSignals( 1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

	call run_HiggsSignals_LHC_Run1_combination(Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, Chisq_LHCRun1, nobs_LHCRun1, Pvalue_LHCRun1) 

	call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)

    call complete_HS_results()
    
    call get_masschi2_from_separation(csq_sep)

	 write(21,*) Mh,dm,scale1,scale2,HSres(1)%Chisq_mu,HSres(1)%Chisq_mh,HSres(1)%Chisq,HSres(1)%nobs-2,HSres(1)%Pvalue, csq_sep
   enddo
  enddo
  close(21)
 enddo

 write(*,*) "Finishing HiggsSignals..."
 call finish_HiggsSignals

 end program HS_2Higgses
