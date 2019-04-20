!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals-2 (TS 29/03/2017).
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
program HSeffC
!--------------------------------------------------------------------------------------
! In this example we scan over two kappa scale factors, kappaF and kappaV, of the
! 125 GeV Higgs boson, using the effective couplings input. The kappa
! factor for the H to gamma-gamma partial width is derived from the two kappa factors.
! Experimental input are the LHC13 measurements, both peak and STXS observables.
! The output is written into /results/HSeffC.dat, which can be plotted with the
! python script plot_HSeffC.py in the results folder.
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall
 implicit none

 integer :: nHzero, nHplus, nobs_peak, nobs_STXS, i, j, k, ii, jj
 double precision :: obsratio, mass, mu
 double precision :: Pvalue_peak, Chisq_peak,  Chisq_peak_mu, Chisq_peak_mh
 double precision :: Pvalue_STXS, Chisq_STXS,  Chisq_STXS_rates, Chisq_STXS_mh
 double precision :: SMGammaTotal
 double precision :: kappaF, kappaV
 double precision :: Mh,GammaTotal,ghjss_s,ghjss_p,ghjcc_s,ghjcc_p, &
&                    ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p, &
&                    ghjmumu_s,ghjmumu_p,ghjtautau_s,ghjtautau_p, &
&                    ghjWW,ghjZZ,ghjZga,ghjgaga,ghjgg,	&
&                    ghjhiZ
 character(len=100)::filename
 double precision :: dm
 integer		  :: pdf
!-HiggsBounds internal functions to obtain SM branching ratios 
 double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau,&
 &                   SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg,SMGamma_h
 double precision :: Htogaga_rate, HtoVV_rate, HtoFF_rate
 nHzero=1
 nHplus=0

!--Setting up the output
 filename='results/HSeffC.dat'
 open(21,file=filename)
 write(21,*) '# mh kappaF kappaV Chisq_mu Chisq ndf Htogaga_rate HtoVV_rate HtoFF_rate' 
 write(21,*) '#--------------------------------------------------------------------#'

!--Enter the Higgs mass and its theory uncertainty here: 
 Mh = 125.09D0
 dm = 0.0D0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
!  call initialize_HiggsSignals(nHzero,nHplus,"LHC13_CMS_H-gaga")
 call initialize_HiggsSignals(nHzero,nHplus,"LHC13")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...)    ----!
 call setup_output_level(0)
 !---- Set the assignment range for the peak-centered method (optional)              ----! 
!  call setup_assignmentrange_massobservables(4.0D0)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)        ----!
 pdf = 2
 call setup_pdf(pdf) 
!---- Pass the Higgs mass uncertainty to HiggsSignals (if relevant)                  ----!
! call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set number of free model parameters ----!
 call setup_Nparam(2)

 do i=1,26!81
  do j=1,16!81
   kappaF = 0.80D0+(i-1)*0.02D0
   kappaV = 0.90D0+(j-1)*0.02D0

   SMGammaTotal=SMGamma_h(Mh)

   if(.not. (SMGammaTotal .lt. 0)) then
    ghjss_s=kappaF
    ghjss_p=0.0d0
    ghjcc_s=kappaF
    ghjcc_p=0.0d0
    ghjbb_s=kappaF
    ghjbb_p=0.0d0
    ghjtt_s=kappaF
    ghjtt_p=0.0d0         
    ghjmumu_s=kappaF
    ghjmumu_p=0.0d0  
    ghjtautau_s=kappaF
    ghjtautau_p=0.0d0
    ghjWW=kappaV
    ghjZZ=kappaV
    ghjZga=kappaV
    ghjgg=kappaF
    ghjhiZ=0d0
    ghjgaga=sqrt(get_g2hgaga(kappaF,kappaF,kappaF,kappaV,kappaV))
      
!----Calculate the new total decay width (not 100% accurate due to independent BR fit functions):
    GammaTotal = SMGammaTotal*(1 + &
	&  	(ghjWW**2.0 - 1)*SMBR_HWW(Mh)+(ghjZZ**2.0 - 1)*SMBR_HZZ(Mh) + &
	&   (ghjgg**2.0 - 1)*SMBR_Hgg(Mh)+(ghjtt_s**2.0 - 1)*SMBR_Htoptop(Mh)+ &
	&   (ghjbb_s**2.0 - 1)*SMBR_Hbb(Mh)+(ghjtautau_s**2.0 - 1)*SMBR_Htautau(Mh)+ &
	&   (ghjss_s**2.0 - 1)*SMBR_Hss(Mh)+(ghjcc_s**2.0 - 1)*SMBR_Hcc(Mh)+ &
	&   (ghjZga**2.0 - 1)*SMBR_HZgam(Mh)+(ghjmumu_s**2.0 - 1)*SMBR_Hmumu(Mh)+ &		
	&   (ghjgaga**2.0 - 1)*SMBR_Hgamgam(Mh)	)

    call HiggsBounds_neutral_input_properties(Mh,GammaTotal)

    call HiggsBounds_neutral_input_effC(               &  
     &          ghjss_s,ghjss_p,ghjcc_s,ghjcc_p,       &
     &          ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p,       &
     &          ghjmumu_s,ghjmumu_p,                   &
     &          ghjtautau_s,ghjtautau_p,               &
     &          ghjWW,ghjZZ,ghjZga,                    &
     &          ghjgaga,ghjgg,ghjhiZ)

    call run_HiggsSignals( 1, Chisq_peak_mu, Chisq_peak_mh, Chisq_peak, nobs_peak, Pvalue_peak)

	call run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)

! This will get the SM normalized rates for inclusive Higgs production,
! with H-> gamma gamma, VV and FF decays:
	call get_rates(1,4,5,(/11,21,31,41,51/),Htogaga_rate)
	call get_rates(1,4,5,(/12,22,32,42,52/),HtoVV_rate)
	call get_rates(1,4,5,(/14,24,34,44,54/),HtoFF_rate)

! 	write(*,*) "Htogaga_rate, Chi2 = ", Htogaga_rate, Chisq_mu
! Try new interface:
!     call HiggsBounds_neutral_input_hadr_channelrates_single(13,1,1,1,2.0d0*SMBR_Hgamgam(Mh))
!     call run_HiggsSignals( 1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
! 	call get_rates(1,4,5,(/11,21,31,41,51/),Htogaga_rate)
! 	write(*,*) "Htogaga_rate (mod), Chi2 = ", Htogaga_rate, Chisq_mu

! This will collect the main HiggsSignals results together into one file
    write(21,*) mh,kappaF,kappaV,Chisq_peak_mu+Chisq_STXS_rates &
    &          ,Chisq_peak+Chisq_STXS,nobs_peak+nobs_STXS,Htogaga_rate,&
    &           HtoVV_rate,HtoFF_rate

   endif
  enddo
 enddo
 
 close(21)

 write(*,*) "Finishing HiggsSignals..."
 call finish_HiggsSignals

contains

!************************************************************** 
 function get_g2hgaga(ghbb,ghtt,ghtautau,ghWW,ghZZ)
! Evaluates g2hgaga from other effective couplings, using partial widths informations
! at a Higgs mass of 126 GeV (calculated with HDECAY and taken from
! http://people.web.psi.ch/spira/higgscoup/ ).
!**************************************************************
 double precision, intent(in) :: ghbb,ghtt,ghtautau,ghWW,ghZZ
 double precision :: get_g2hgaga
  
 get_g2hgaga = (ghtt**2)*0.70904D-01 + (ghbb**2)*0.18760D-04 + (ghWW**2)*1.5863 + &
 & ghtt*ghbb*(-0.17319D-02) + ghtt*ghWW*(-0.67074) + &
 & ghbb*ghWW*0.82093D-02 + (ghtautau**2)*0.22663E-04 + &
 & ghtt*ghtautau*(-0.18696E-02) + ghbb*ghtautau*0.41239E-04 +&
 & ghtautau*ghWW*0.88634E-02
 
 end function get_g2hgaga
!**************************************************************
 end program HSeffC
