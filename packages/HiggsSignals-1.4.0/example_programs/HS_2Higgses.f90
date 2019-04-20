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
 use pc_chisq, only : print_inverse_cov_mh_to_file, get_masschi2_from_separation
 implicit none

 integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
 double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
 double precision :: SMGammaTotal(2)
double precision :: scale_bbh, scale_ggh, dggh, dbbh
 double precision :: Mh(2),GammaTotal(2),g2hjss_s(2),g2hjss_p(2),g2hjcc_s(2),g2hjcc_p(2), &
&                    g2hjbb_s(2),g2hjbb_p(2),g2hjtt_s(2),g2hjtt_p(2), &
&                    g2hjmumu_s(2),g2hjmumu_p(2),g2hjtautau_s(2),g2hjtautau_p(2), &
&                    g2hjWW(2),g2hjZZ(2),g2hjZga(2),g2hjgaga(2),g2hjgg(2),g2hjggZ(2),	&
&                    g2hjhiZ(2,2),BR_hjhihi(2,2),BR_hjinvisible(2,2)
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
 Mh = (/ 125.0D0, 127.0D0 /)
 dm = (/ 0.5D0, 1.0D0 /)
! dm = (/ 0.0D0, 0.0D0 /)

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
 call initialize_HiggsSignals(nHzero,nHplus,"latestresults-1.3.0-LHCinclusive")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
 call setup_output_level(0) 
! call setup_anticorrelations_in_mh(1)
! call setup_assignmentrange_massobservables(2.0D0)
!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
 call HiggsSignals_neutral_input_MassUncertainty(dm)

 do pdf = 1,3
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
  call setup_pdf(pdf) 
  write(pdfchar,'(I1)') pdf
  filename='results/pdf'//pdfchar//'.dat'
  open(21,file=trim(adjustl(filename)))
  write(21,*) '# Mh(1) Mh(2) dm(1) dm(2) mu(1) mu(2) chisq_mu chisq_mh chisq ndf Pvalue chisq_sep'
  write(21,*) '#--------------------------------------------------------------------------'

 do i=1,31
  do j=1,31
   Mh(1) = 122.0D0 + (i-1)*0.2D0
   Mh(2) = 122.0D0 + (j-1)*0.2D0
             
   scale1 = 0.75D0
   scale2 = 0.25D0

   SMGammaTotal(1)=SMGamma_h(Mh(1))
   SMGammaTotal(2)=SMGamma_h(Mh(2))   

    g2hjss_s(1) = scale1
    g2hjss_s(2) = scale2
    g2hjss_p=0.0d0
    g2hjcc_s(1) = scale1
    g2hjcc_s(2) = scale2    
    g2hjcc_p=0.0d0
    g2hjbb_s(1) = scale1
    g2hjbb_s(2) = scale2    
    g2hjbb_p=0.0d0
    g2hjtt_s(1) = scale1
    g2hjtt_s(2) = scale2    
    g2hjtt_p=0.0d0         
    g2hjmumu_s(1) = scale1
    g2hjmumu_s(2) = scale2    
    g2hjmumu_p=0.0d0  
    g2hjtautau_s(1) = scale1
    g2hjtautau_s(2) = scale2
    g2hjtautau_p=0.0d0
    g2hjWW(1) = scale1
    g2hjWW(2) = scale2    
    g2hjZZ(1) = scale1
    g2hjZZ(2) = scale2    
    g2hjZga(1) = scale1
    g2hjZga(2) = scale2    
    g2hjgg(1) = scale1
    g2hjgg(2) = scale2    
    g2hjggZ(1) = scale1
    g2hjggZ(2) = scale2    
    g2hjhiZ=0d0
    g2hjgaga(1) = scale1
    g2hjgaga(2) = scale2    
    BR_hjhihi=0d0
    BR_hjinvisible=0d0
      
!----Calculate the new total decay width:
    GammaTotal(1) = SMGammaTotal(1)*(1 + &
	&  	            (g2hjWW(1) - 1)*SMBR_HWW(Mh(1))+(g2hjZZ(1) - 1)*SMBR_HZZ(Mh(1)) + &
	&               (g2hjgg(1) - 1)*SMBR_Hgg(Mh(1))+(g2hjtt_s(1) - 1)*SMBR_Htoptop(Mh(1))+ &
	&               (g2hjbb_s(1) - 1)*SMBR_Hbb(Mh(1))+(g2hjtautau_s(1) - 1)*SMBR_Htautau(Mh(1))+ &
	&               (g2hjss_s(1) - 1)*SMBR_Hss(Mh(1))+(g2hjcc_s(1) - 1)*SMBR_Hcc(Mh(1))+ &
	&               (g2hjZga(1) - 1)*SMBR_HZgam(Mh(1))+(g2hjmumu_s(1) - 1)*SMBR_Hmumu(Mh(1))+ &		
	&               (g2hjgaga(1) - 1)*SMBR_Hgamgam(Mh(1))	)

    GammaTotal(2) = SMGammaTotal(2)*(1 + &
	&  	            (g2hjWW(2) - 1)*SMBR_HWW(Mh(2))+(g2hjZZ(2) - 1)*SMBR_HZZ(Mh(2)) + &
	&               (g2hjgg(2) - 1)*SMBR_Hgg(Mh(2))+(g2hjtt_s(2) - 1)*SMBR_Htoptop(Mh(2))+ &
	&               (g2hjbb_s(2) - 1)*SMBR_Hbb(Mh(2))+(g2hjtautau_s(2) - 1)*SMBR_Htautau(Mh(2))+ &
	&               (g2hjss_s(2) - 1)*SMBR_Hss(Mh(2))+(g2hjcc_s(2) - 1)*SMBR_Hcc(Mh(2))+ &
	&               (g2hjZga(2) - 1)*SMBR_HZgam(Mh(2))+(g2hjmumu_s(2) - 1)*SMBR_Hmumu(Mh(2))+ &		
	&               (g2hjgaga(2) - 1)*SMBR_Hgamgam(Mh(2))	)

    call HiggsBounds_neutral_input_effC(Mh,GammaTotal, &
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p, &
     &    g2hjtt_s,g2hjtt_p, &
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ, &
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

    call run_HiggsSignals( 1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

    call get_masschi2_from_separation(csq_sep)

	 write(21,*) Mh,dm,scale1,scale2,Chisq_mu,Chisq_mh,Chisq,ndf,Pvalue, csq_sep
   enddo
  enddo
  close(21)
 enddo

 write(*,*) "Finishing HiggsSignals..."
 call finish_HiggsSignals

 end program HS_2Higgses
