!------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013).
!------------------------------------------------------------
subroutine initialize_HiggsSignals_latestresults(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals with the experimental
! dataset "latestresults", avoiding to specify this via a string argument.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=13) :: Expt_string
 
 Expt_string = "latestresults"

 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_latestresults
!------------------------------------------------------------
subroutine initialize_HiggsSignals_LHC13(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals with the experimental
! dataset "latestresults", avoiding to specify this via a string argument.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=13) :: Expt_string
 
 Expt_string = "LHC13"

 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_LHC13
!------------------------------------------------------------
subroutine initialize_HiggsSignals_empty(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals without dataset.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=13) :: Expt_string
 
 Expt_string = "none"

 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_empty
!------------------------------------------------------------
subroutine initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
!------------------------------------------------------------
! This the first HiggsSignals subroutine that should be called
! by the user.
! It calls subroutines to read in the tables of Standard Model 
! decay and production rates from HiggsBounds, sets up the 
! experimental data from Tevatron and LHC, allocate arrays, etc.
! Arguments (input):
!   * nHiggs = number of neutral Higgs in the model 
!   * nHiggsplus = number of singly, positively charged Higgs in the model
!   * Expt_string = name of experimental dataset to be used
!------------------------------------------------------------
 use usefulbits, only : np,Hneut,Hplus,Chineut,Chiplus,debug,inputmethod,&
  &   theo,whichanalyses,just_after_run,&
  &   file_id_debug1,file_id_debug2,allocate_if_stats_required
 use usefulbits_HS, only : HiggsSignals_info, nanalys, eps, Exptdir, obs
 use datatables, only: setup_observables, setup_LHC_Run1_combination
 use STXS, only : load_STXS
 use input, only : check_number_of_particles,check_whichanalyses
 use io, only : setup_input_for_hs, setup_output_for_hs
 use theory_BRfunctions, only : setup_BRSM, BRSM
 use theory_XS_SM_functions, only : setup_XSSM, XSSM 

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=*), intent(in) :: Expt_string
 !-----------------------------------internal
 integer :: i
 logical :: exptdirpresent = .False.
 !----------------------------------parameter
 eps=5.0D0
 np(Hneut)=nHiggsneut

 np(Hplus)=nHiggsplus
 
 if(Expt_string.ne.'none') then
  Exptdir = Expt_string
  exptdirpresent = .True.
 endif
  
 
 np(Chineut)=0! not considering bounds on neutralinos here
 np(Chiplus)=0! not considering bounds on charginos here
        
 debug=.False.

 select case(whichanalyses)
  case('onlyL')
   whichanalyses='LandH'
  case('onlyH','onlyP','list ','LandH')
  case default
   whichanalyses='onlyH'
 end select 

 call HiggsSignals_info 
 if(inputmethod=='subrout') then 
  if(allocated(theo))then
   if(debug) write(*,*) "HiggsBounds/HiggsSignals internal structure already initialized!"
  else
   if(debug)write(*,*)'doing other preliminary tasks...'      ; call flush(6)
   call setup_input_for_hs

!    allocate(inputsub( 2 )) !(1)np(Hneut)>0 (2)np(Hplus)>0
!    inputsub(1)%desc='HiggsBounds_neutral_input_*'    ; inputsub(1)%req=req(   0,   1)
!    inputsub(2)%desc='HiggsBounds_charged_input'      ; inputsub(2)%req=req(   1,   0)
!  
!    do i=1,ubound(inputsub,dim=1)
!     inputsub(i)%stat=0
!    enddo
  endif
 endif 
 
 if(debug)write(*,*)'reading in Standard Model tables...'   ; call flush(6) 
 if(.not.allocated(BRSM)) call setup_BRSM 
 if(.not.allocated(XSSM)) call setup_XSSM
 call setup_uncertainties
             
 if(debug)write(*,*)'reading in experimental data...'       ; call flush(6)
 if(exptdirpresent) call setup_observables
 if(exptdirpresent) call load_STXS(Expt_string)
  
 call setup_LHC_Run1_combination

 if(debug)write(*,*)'sorting out processes to be checked...'; call flush(6)
 nanalys = size(obs)

 if(debug)write(*,*)'preparing output arrays...'            ; call flush(6)
 call setup_output_for_hs

 if(debug)write(*,*)'HiggsSignals has been initialized...'  ; call flush(6)

 just_after_run=.False.
  
!  contains 
   !         |   np
   !         |Hneu Hcha 
   !         | ==0  ==0 
!  function req(Hneu,Hcha)
!   integer, intent(in) ::Hneu,Hcha
!   integer :: req
!   
!   req=1
!   if(np(Hneut)==0)  req= Hneu  * req
!   if(np(Hplus)==0)  req= Hcha  * req
! 
!  end function req 

end subroutine initialize_HiggsSignals
!------------------------------------------------------------
subroutine HiggsSignals_neutral_input_MassUncertainty(dMh)
! Sets the theoretical mass uncertainty of the Higgs bosons.
!------------------------------------------------------------
 use usefulbits, only: theo,np,Hneut
 
 implicit none
 double precision,intent(in) :: dMh(np(Hneut))

 if(.not.allocated(theo))then
  stop 'subroutine HiggsSignals_initialize must be called first'
 endif
 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsSignal_neutral_input_MassUncertainty should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsSignal_neutral_input_MassUncertainty'
 endif
 
 theo(1)%particle(Hneut)%dM = dMh
 
end subroutine HiggsSignals_neutral_input_MassUncertainty
!------------------------------------------------------------
subroutine setup_uncertainties
!------------------------------------------------------------
 use usefulbits, only : file_id_common3
 use store_pathname_hs, only : pathname_HS
 use usefulbits_hs, only : delta_rate
 use io, only : read_matrix_from_file
 
 logical :: BRmodel, BRSM, XSmodel, XSSM

 call read_matrix_from_file(9,pathname_HS//"BRcov.in",delta_rate%BRcov, BRmodel)
 call read_matrix_from_file(9,pathname_HS//"BRcovSM.in",delta_rate%BRcovSM, BRSM)
 call read_matrix_from_file(11,pathname_HS//"XScov.in",delta_rate%CScov, XSmodel)
 call read_matrix_from_file(11,pathname_HS//"XScovSM.in",delta_rate%CScovSM, XSSM)
 call read_matrix_from_file(11,pathname_HS//"XScov_13TeV.in",delta_rate%CS13cov, XSmodel)
 call read_matrix_from_file(11,pathname_HS//"XScovSM_13TeV.in",delta_rate%CS13covSM, XSSM)

 if(BRmodel.and.BRSM) then
  delta_rate%BRcov_ok=.True.
! write(*,*) "Covariance matrix for relative branching ratio uncertainties read in successfully."
 else
! write(*,*) "Covariance matrix for relative branching ratio uncertainties not provided. Using default values."
 endif
 if(XSmodel.and.XSSM) then
  delta_rate%CScov_ok=.True.
! write(*,*) "Covariance matrix for relative cross section uncertainties read in successfully."
 else
! write(*,*) "Covariance matrix for relative cross section uncertainties not provided. Using default values."
 endif
    
end subroutine setup_uncertainties
!------------------------------------------------------------
subroutine setup_rate_normalization(normalize_to_refmass, normalize_to_refmass_outside_dmtheo)
 use usefulbits_hs, only : normalize_rates_to_reference_position,&
 &                          normalize_rates_to_reference_position_outside_dmtheo
 implicit none

 logical, intent(in) :: normalize_to_refmass
 logical, intent(in) ::  normalize_to_refmass_outside_dmtheo

 if(normalize_to_refmass) then
  write(*,*) "Using SM rate prediction at observed mass for signal strength calculation."
 else
  write(*,*) "Using SM rate prediction at predicted mass for signal strength calculation."
 endif 
 if(normalize_to_refmass_outside_dmtheo) then
  write(*,*) "If predicted mass and observed mass do not agree within theory uncertainty:",&
  &          " SM rate prediction at observed mass is used for signal strength calculation."
 else
  write(*,*) "If predicted mass and observed mass do not agree within theory uncertainty:",&
  &          " SM rate prediction at predicted mass is used for signal strength calculation."
 endif
 normalize_rates_to_reference_position = normalize_to_refmass
 normalize_rates_to_reference_position_outside_dmtheo = normalize_to_refmass_outside_dmtheo
 
end subroutine setup_rate_normalization
!------------------------------------------------------------
subroutine setup_model_rate_uncertainties(filename_XS, filename_XS13, filename_BR)
!------------------------------------------------------------
 use usefulbits, only : file_id_common3
 use store_pathname_hs, only : pathname_HS
 use usefulbits_hs, only : delta_rate
 use io, only : read_matrix_from_file
 
 character(LEN=*),intent(in) :: filename_XS, filename_XS13, filename_BR
 logical :: BRmodel, XSmodel

 call read_matrix_from_file(9,filename_BR,delta_rate%BRcov, BRmodel)
 call read_matrix_from_file(11,filename_XS,delta_rate%CScov, XSmodel)
 call read_matrix_from_file(11,filename_XS13,delta_rate%CS13cov, XSmodel)
 
 if(BRmodel.and.XSmodel) then
  delta_rate%BRcov_ok=.True.
  delta_rate%CScov_ok=.True.
  write(*,*) "Covariance matrices for rate uncertainties read in successfully."
 else
  write(*,*) "Covariance matrix for rate uncertainties not provided. Using default values."
 endif
    
end subroutine setup_model_rate_uncertainties
!------------------------------------------------------------
subroutine setup_rate_uncertainties( dCS, dBR )
! OUTDATED !!!
!------------------------------------------------------------
! Sets (relative) systematic uncertainties of the model for:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!------------------------------------------------------------
 use usefulbits_hs, only : delta_rate
 implicit none
 
 double precision, intent(in) :: dCS(5)
 double precision, intent(in) :: dBR(5)
 integer :: i
   
 delta_rate%dCS = dCS

 do i=lbound(dBR,dim=1),ubound(dBR,dim=1) 
  call setup_dbr(i,dBR(i))
 enddo
 
end subroutine setup_rate_uncertainties
!------------------------------------------------------------
subroutine setup_dbr(BRid, value)
! OUTDATED !!!
!------------------------------------------------------------
 use usefulbits_hs, only : delta_rate

 integer,intent(in) :: BRid
 double precision, intent(in) :: value

 if(BRid.gt.0.and.BRid.lt.10) then
  delta_rate%dBR(BRid) = value
 else
  write(*,*) "Warning in setup_dbr: Unknown decay mode."
 endif
   
end subroutine setup_dbr
!------------------------------------------------------------
subroutine setup_correlations(corr_mu, corr_mh)
!------------------------------------------------------------
! With this subroutine the user may switch off/on correlations
! (default=on) by setting corr = 0/1.
!------------------------------------------------------------
 use usefulbits_hs, only : correlations_mu, correlations_mh
 implicit none
 
 integer, intent(in) :: corr_mu, corr_mh
 if(corr_mu.eq.0) then
  correlations_mu = .False.
  write(*,*) 'Correlations in signal strength observables are switched off.'
 elseif(corr_mu.eq.1) then
  correlations_mu = .True.
 else
  stop 'Error: Correlations must be switched on/off by an integer value of 0 or 1.' 
 endif
 if(corr_mh.eq.0) then
  correlations_mh = .False.
  write(*,*) 'Correlations in Higgs mass observables are switched off.'
 elseif(corr_mh.eq.1) then
  correlations_mh = .True.
 else
  stop 'Error: Correlations must be switched on/off by an integer value of 0 or 1.' 
 endif
end subroutine setup_correlations
!------------------------------------------------------------
subroutine setup_symmetricerrors(symm)
! Sets the measured rate uncertainties to either a symmetrical average
! of the upper and lower cyan band widths (symm==1) or else uses the 
! original (asymmetrical) errors.
!------------------------------------------------------------
 use usefulbits_hs, only : symmetricerrors
 implicit none

 integer, intent(in) :: symm
 if(symm.eq.1) then
  write(*,*) "Using averaged (symmetrical) experimental rate uncertainties."
  symmetricerrors = .True.
 else
  write(*,*) "Using original (asymmetrical) experimental rate uncertainties."
  symmetricerrors = .False.
 endif
  
end subroutine setup_symmetricerrors
!------------------------------------------------------------
subroutine setup_absolute_errors(absol)
! Treats the measured rate uncertainties as either absolute
! uncertainties (1) or relative (0). By default, they are
! treated as relative uncertainties.
!------------------------------------------------------------
 use usefulbits_hs, only : absolute_errors
 implicit none

 integer, intent(in) :: absol
 if(absol.eq.1) then
  write(*,*) "Using absolute experimental rate uncertainties."
  absolute_errors = .True.
 else
  write(*,*) "Using relative experimental rate uncertainties."
  absolute_errors = .False.
 endif
  
end subroutine setup_absolute_errors
!------------------------------------------------------------
subroutine setup_correlated_rate_uncertainties(corr)
! OUTDATED
!------------------------------------------------------------
use usefulbits_hs, only : delta_rate
integer, intent(in) :: corr

if(corr.eq.0) then
 delta_rate%usecov = .False.
 write(*,*) "Deactivated correlated CS and BR uncertainties. Using approximated maximum error."
elseif(corr.eq.1) then
 delta_rate%usecov = .True.
 write(*,*) "Activated correlated CS and BR uncertainties. Using them if covariance matrices are present."
else
 write(*,*) "Warning in subroutine setup_correlated_rate_uncertainties: Argument ",corr," is not equal to 0 or 1."
endif
end subroutine setup_correlated_rate_uncertainties
!------------------------------------------------------------
subroutine setup_SMweights(useweight)
! If set to 1 (true), HiggsSignals assumes the same signal decomposition
! (weights) as in the SM for the given model. This will enter the determination
! of the theoretical rate uncertainty.
!
! OUTDATED
!------------------------------------------------------------
 use usefulbits_hs, only : useSMweights
 implicit none

 integer, intent(in) :: useweight
 if(useweight.eq.1) then
  write(*,*) "Using SM weights for theoretical rate uncertainties of the model."
  useSMweights = .True.
 else
  write(*,*) "Using true model weights for theoretical rate uncertainties of the model."
  useSMweights = .False.
 endif
  
end subroutine setup_SMweights
!------------------------------------------------------------
subroutine setup_anticorrelations_in_mu(acorr)
! Allows for anti-correlations in the signal strength covariance
! matrix if there is a relative sign difference in two mu measurements
! (acorr==1) or else uses only correlations irrespective of the relative
! (acorr==0).
!
! OUTDATED
!------------------------------------------------------------
 use usefulbits_hs, only : anticorrmu
 implicit none

 integer, intent(in) :: acorr
 if(acorr.eq.1) then
  write(*,*) "Allow anti-correlated signal strength measurements."
  anticorrmu = .True.
 else
  write(*,*) "Prohibit anti-correlated signal strength measurements."
  anticorrmu = .False.
 endif
  
end subroutine setup_anticorrelations_in_mu
!------------------------------------------------------------
subroutine setup_anticorrelations_in_mh(acorr)
! Allows for anti-correlations in the mass covariance
! matrix if there is a relative sign difference in two mu measurements
! (acorr==1) or else uses only correlations irrespective of the relative
! (acorr==0).
!
! OUTDATED
!------------------------------------------------------------
 use usefulbits_hs, only : anticorrmh
 implicit none

 integer, intent(in) :: acorr
 if(acorr.eq.1) then
  write(*,*) "Allow anti-correlated mass measurements."
  anticorrmh = .True.
 else
  write(*,*) "Prohibit anti-correlated mass measurements."
  anticorrmh = .False.
 endif
  
end subroutine setup_anticorrelations_in_mh
!------------------------------------------------------------
subroutine setup_assignmentrange(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange,assignmentrange_massobs, pdf
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange  
 else
  assignmentrange = range
  assignmentrange_massobs = range
 endif 

 if(assignmentrange.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
 
end subroutine setup_assignmentrange
!------------------------------------------------------------
subroutine setup_assignmentrange_LHCrun1(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange_LHCrun1, pdf
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange_LHCrun1
 else
  assignmentrange_LHCrun1 = range
 endif 

!  if(assignmentrange_LHCrun1.ne.1.0D0.and.pdf.eq.1) then
!   write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
!  endif
 
end subroutine setup_assignmentrange_LHCrun1
!------------------------------------------------------------
subroutine setup_assignmentrange_massobservables(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange_massobs, pdf
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange_massobs
 else
  assignmentrange_massobs = range
 endif 

 if(assignmentrange_massobs.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
 
end subroutine setup_assignmentrange_massobservables
!------------------------------------------------------------
subroutine setup_assignmentrange_STXS(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange_STXS
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange_STXS
 else
  assignmentrange_STXS = range
 endif 
 
end subroutine setup_assignmentrange_STXS
!------------------------------------------------------------

subroutine setup_nparam(Np)
!------------------------------------------------------------
 use usefulbits_hs, only : Nparam
 implicit none 
 integer, intent(in) :: Np
 Nparam = Np
end subroutine setup_nparam
!------------------------------------------------------------
subroutine setup_Higgs_to_peaks_assignment_iterations(iter)
! Sets the number of iterations for the Higgs-to-peak-assignment.
!------------------------------------------------------------
 use usefulbits_hs, only : iterations
 implicit none
 integer, intent(in) :: iter 
 iterations = iter 
 
end subroutine setup_Higgs_to_peaks_assignment_iterations
!------------------------------------------------------------
subroutine setup_mcmethod_dm_theory(mode)
! OUTDATED
 use mc_chisq, only : mc_mode
 implicit none
 integer, intent(in) :: mode
 character(LEN=14) :: mode_desc(2) = (/'mass variation','convolution   '/)

 if(mode.eq.1.or.mode.eq.2) then
  mc_mode = mode
  write(*,'(1X,A,A)') 'The mass-centered chi^2 method will treat the Higgs',&
& ' boson mass theory uncertainty by '//trim(mode_desc(mode))//'.'
 else
  stop 'Error in subroutine setup_mcmethod_dm_theory: Unknown mode (1 or 2 possible)!'
 endif 
end subroutine setup_mcmethod_dm_theory
!------------------------------------------------------------
subroutine setup_sm_test(int_SMtest,epsilon)
! With this subroutine the user may switch off the SM likeness test
! (default=on) or change the maximal deviation epsilon (default=5.0D-2)
! OUTDATED / NOT USED
!------------------------------------------------------------
 use usefulbits_hs, only : useSMtest, eps
 implicit none
 
 integer, intent(in) :: int_SMtest
 double precision, intent(in) :: epsilon
 
 if(int_SMtest.eq.0) then
  useSMtest = .False.
  write(*,*) 'SM likeness test has been switched off.'
 elseif(int_SMtest.eq.1) then
  useSMtest = .True.
  write(*,*) 'SM likeness test has been switched on.'  
 else
  stop 'Error: SM test must be switched on/off by an integer value of 0 or 1.' 
 endif
 eps = epsilon
end subroutine setup_sm_test
!------------------------------------------------------------
subroutine setup_thu_observables(thuobs)
! OUTDATED
 use usefulbits_hs, only : THU_included
 integer, intent(in) :: thuobs

 if(thuobs.eq.0) then
  THU_included = .False.
  write(*,*) 'Observables are assumed to NOT include theory errors.'
 else
  THU_included = .True.
  write(*,*) 'Observables are assumed to include theory errors.'  
 endif
  
end subroutine setup_thu_observables
!------------------------------------------------------------
subroutine setup_output_level(level)
! Controls the level of information output:
! 0 : silent mode
! 1 : screen output for each analysis with its peak/mass-centered observables and
!     their respective values predicted by the model
! 2 : screen output of detailed information on each analysis with its 
!     peak/mass-centered observables
! 3 : creates the files peak_information.txt and peak_massesandrates.txt
!------------------------------------------------------------
 use usefulbits_hs, only : output_level, additional_output
 implicit none
 integer, intent(in) :: level

 if(level.eq.0.or.level.eq.1.or.level.eq.2.or.level.eq.3) then
  output_level = level
 else
  stop 'Error in subroutine setup_output_level: level not equal to 0,1,2 or 3.' 
 endif
 if(level.eq.3) additional_output = .True.
 
end subroutine setup_output_level 
!------------------------------------------------------------
subroutine setup_pdf(pdf_in)
! Sets the probability density function for the Higgs mass uncertainty parametrization:
! 1 : box-shaped pdf
! 2 : Gaussian pdf
! 3 : box-shaped theory error + Gaussian experimental pdf
!------------------------------------------------------------
 use usefulbits_hs, only : pdf, assignmentrange

 implicit none
 integer, intent(in) :: pdf_in
 character(LEN=13) :: pdf_desc(3) = (/'box         ','Gaussian    ','box+Gaussian'/)
 
 pdf=pdf_in
 if((pdf.eq.1).or.(pdf.eq.2).or.(pdf.eq.3)) then
!   write(*,'(1X,A,A,1I1,A)') 'Use a '//trim(pdf_desc(pdf))//' probability density function ',&
! &  'for the Higgs mass(es) (pdf=',pdf,')'
 endif
   
 if(assignmentrange.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
   
end subroutine setup_pdf
!------------------------------------------------------------
!subroutine assign_toyvalues_to_observables(ii, peakindex, npeaks, mu_obs, mh_obs)
!! Assigns toy values to the peak's mass and mu value for analysis ii.
!! ii           :: analysis number (entry in mutables)
!! peakindex    :: index of the peak of analysis ii
!! npeaks 	   :: number of peaks found in analysis ii
!! mu_obs       :: toy value for mu to be given to the peak with peakindex
!! mh_obs       :: toy value for mh to be given to the peak with peakindex
!------------------------------------------------------------
! use usefulbits_hs, only: obs, usetoys
! 
! integer, intent(in) :: ii, peakindex, npeaks
! double precision, intent(in) :: mh_obs, mu_obs
!  
! if(peakindex.gt.npeaks) then
!  stop 'Error in subroutine assign_toyvalues_to_observables: Observable does not exist!'
! endif
!   
! obs(ii)%table%npeaks = npeaks
! if(.not.allocated(obs(ii)%table%Toys_muobs)) allocate(obs(ii)%table%Toys_muobs(npeaks))
! if(.not.allocated(obs(ii)%table%Toys_mhobs)) allocate(obs(ii)%table%Toys_mhobs(npeaks)) 
!
! obs(ii)%table%Toys_muobs(peakindex) = mu_obs
! obs(ii)%table%Toys_mhobs(peakindex) = mh_obs
! 
! usetoys = .True.
! 
!end subroutine assign_toyvalues_to_observables
!------------------------------------------------------------
subroutine assign_toyvalues_to_peak(ID, mu_obs, mh_obs)
! Assigns toy values to the peak's mass and mu value to a peak observable.
! ID           :: observable ID
! mu_obs       :: toy value for mu to be given to the peak
! mh_obs       :: toy value for mh to be given to the peak
!
! n.B.: Do we also want to set mu uncertainties here?
!------------------------------------------------------------
 use usefulbits_hs, only: obs, usetoys
 implicit none

 integer, intent(in) :: ID
 double precision, intent(in) :: mh_obs, mu_obs
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then    
  obs(pos)%peak%mpeak = mh_obs
  obs(pos)%peak%mu = mu_obs
  usetoys = .True.
 else
  write(*,*) "WARNING in assign_toyvalues_to_peak: ID unknown."
 endif
  
end subroutine assign_toyvalues_to_peak
!------------------------------------------------------------
subroutine assign_modelefficiencies_to_peak(ID, Nc, eff_ratios)
! Assigns to each channel of the observable the efficiency in the model 
! w.r.t the SM efficiency (as a ratio!)
! 
! ID           :: observable ID
! Nc			:: number of channels
! eff_ratios    :: array of length (Number of channels) giving the efficiency ratios
!
! Note: You can first employ the subroutine get_peak_channels (io module) to obtain
!       the relevant channel information of the observable.
!------------------------------------------------------------
 use usefulbits_hs, only: obs
 implicit none

 integer, intent(in) :: ID, Nc
 double precision, dimension(Nc), intent(in) :: eff_ratios
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then    
  if(size(eff_ratios,dim=1).ne.obs(pos)%table%Nc) then
   write(*,*) "WARNING in assign modelefficiencies_to_peak: Number of channels (",&
&  size(eff_ratios,dim=1),"!=",obs(pos)%table%Nc,"does not match for observable ID = ",ID
  else
   obs(pos)%table%channel_eff_ratios = eff_ratios
  endif 
 else
  write(*,*) "WARNING in assign_modelefficiencies_to_peak: ID unknown."
 endif
  
end subroutine assign_modelefficiencies_to_peak
!------------------------------------------------------------
subroutine assign_rate_uncertainty_scalefactor_to_peak(ID, scale_mu)
! Assigns a rate uncertainty scalefactor to the peak specified by ID.
! This scalefactor will only scale the experimental rate uncertainties.
! The theory rate uncertainties must be given manually via setup_rate_uncertainties.
!
! ID       :: observable ID of the peak observable
! scale_mu :: scale_mu by which the mu uncertainty is scaled
!
! OUTDATED / UNUSED
!------------------------------------------------------------
 use usefulbits_hs, only: obs, usescalefactor
 implicit none
 
 integer, intent(in) :: ID
 double precision, intent(in) :: scale_mu
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then
  obs(pos)%peak%scale_mu = scale_mu
 else
  write(*,*) "WARNING in assign_uncertainty_scalefactors_to_peak: ID unknown."
 endif
 usescalefactor = .True.
 
end subroutine assign_rate_uncertainty_scalefactor_to_peak
!------------------------------------------------------------
subroutine run_HiggsSignals_LHC_Run1_combination(Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)

 use usefulbits, only : theo,just_after_run, ndat
 use theo_manip, only : HB5_complete_theo 
 use usefulbits_HS, only : HSres, output_level, Nparam

 implicit none
 !----------------------------------------output
 integer,intent(out) ::           nobs
 double precision,intent(out) ::  Pvalue, Chisq, Chisq_mu, Chisq_mh
 !-------------------------------------internal
 integer :: n,i, nobs_mu, nobs_mh
 logical :: debug=.False.
 !---------------------------------------------

 if(.not.allocated(theo))then
  stop 'subroutine HiggsSignals_initialize must be called first'
 endif

 call HB5_complete_theo       

 do n=1,ndat

  call evaluate_LHC_Run1_combination(theo(n),n)       
      
  Pvalue  = HSres(n)%Pvalue_LHCRun1
  Chisq   = HSres(n)%Chisq_LHCRun1       
  Chisq_mu   = HSres(n)%Chisq_LHCRun1_mu
  Chisq_mh   = HSres(n)%Chisq_LHCRun1_mh
  nobs_mu = HSres(n)%nobs_LHCRun1_mu
  nobs_mh = HSres(n)%nobs_LHCRun1_mh
  nobs = nobs_mu+nobs_mh
  
 if(output_level.ne.0) then    
  write(*,*)
  write(*,*) '#*************************************************************************#' 
  write(*,*) '#         HIGGSSIGNALS RESULTS (LHC ATLAS + CMS Run1 combination)         #'
  write(*,*) '#*************************************************************************#' 
  write(*,'(A55,F21.8)') 'chi^2 from signal rate observables = ',Chisq_mu
  write(*,'(A55,F21.8)') 'chi^2 from Higgs mass observables = ',Chisq_mh
  write(*,'(A55,F21.8)') 'chi^2 (total) = ',Chisq
  write(*,'(A55,I21)') 'Number of rate observables = ', nobs_mu
  write(*,'(A55,I21)') 'Number of mass observables = ', nobs_mh
  write(*,'(A55,I21)') 'Number of observables (total) = ', nobs
  write(*,'(A48,I3,A4,F21.8)') 'Probability (ndf =',nobs-Nparam,') = ', Pvalue
  write(*,*) '#*************************************************************************#' 
  write(*,*)
 endif

 enddo

 just_after_run=.True.

end subroutine run_HiggsSignals_LHC_Run1_combination
!------------------------------------------------------------
subroutine setup_LHC_combination_run1_SMXS_from_paper(useSMXS_from_paper)
!------------------------------------------------------------
 use usefulbits_hs, only : LHC_combination_run1_SMXS_from_paper
 implicit none
 
 logical, intent(in) :: useSMXS_from_paper
 
 if(useSMXS_from_paper) then
  write(*,*) "Using SM cross sections from Tab.1 of arXiv:1606.02266 for LHC Run 1 combination chi^2 test."
 else
  write(*,*) "Using internal SM cross sections for LHC Run 1 combination chi^2 test."  
 endif
  
 LHC_combination_run1_SMXS_from_paper = useSMXS_from_paper
 
end subroutine setup_LHC_combination_run1_SMXS_from_paper
!------------------------------------------------------------
subroutine evaluate_LHC_Run1_combination( t , n )
!------------------------------------------------------------
! internal routine
!------------------------------------------------------------
 use usefulbits, only : np,Hneut,Hplus,dataset,results, vsmall
 use usefulbits_hs, only : HSresults, output_level, Nparam, pdf, &
& LHCrun1_rates, LHCrun1_correlationmatrix, useaveragemass, &
& assignmentrange_LHCrun1, HSres, normalize_rates_to_reference_position, &
& normalize_rates_to_reference_position_outside_dmtheo
 use pc_chisq, only : csq_mh
 use numerics, only : invmatrix, matmult, gammp
 
 implicit none
 !--------------------------------------input      
 type(dataset), intent(in) :: t
 integer, intent(in) :: n      
 !--------------------------------------output
!  type(HSresults), intent(inout) :: r
 !--------------------------------------internal
 integer :: p, d, id, i, j, k, ncomb
 double precision, allocatable :: covmat(:,:), invcovmat(:,:)
 double precision, allocatable :: covmatzero(:,:), invcovmatzero(:,:)
 double precision, dimension(20) :: v, v2, csq_mu, vzero, vzero2, csq_mu_max
 double precision, dimension(20,1) :: vmat, vzeromat 
 double precision :: mobs = 125.09D0
 double precision :: dmobs = 0.24D0
 double precision :: dmbbtautau = 20.0D0 
 double precision :: dmWW = 5.0D0
 double precision :: expmassrange, allowed_massrange
 double precision :: Higgs_signal_k
 double precision :: num1, num2, dnum1, dnum2, denom1, denom2, mav, dmav
   
 allocate(covmat(20,20),invcovmat(20,20))
 allocate(covmatzero(20,20),invcovmatzero(20,20)) 
 
 
 mav =0.0D0
 dmav = 0.0D0

  denom1 = 0.0D0
  denom2 = 0.0D0
  num1 = 0.0D0
  num2 = 0.0D0
  dnum1 = 0.0D0
  dnum2 = 0.0D0
 
 do i=lbound(LHCrun1_rates,dim=1),ubound(LHCrun1_rates,dim=1)
  id = LHCrun1_rates(i)%channel_id
  p = int((id-modulo(id,10))/dble(10))
  d = modulo(id,10)
  
  if(d.eq.4.or.d.eq.5) then
    expmassrange = dmbbtautau
  elseif(d.eq.2) then  
    expmassrange = dmWW
  else
    if(pdf.eq.1) then
     expmassrange = dmobs  
    else
     expmassrange = assignmentrange_LHCrun1*dmobs  
    endif
  endif
  
  LHCrun1_rates(i)%r_pred = 0.0D0

  ncomb = 0
  do k=1,np(Hneut)
  if(pdf.eq.1) then
   allowed_massrange = expmassrange + t%particle(Hneut)%dM(k)
  else
   allowed_massrange = sqrt(expmassrange**2.0D0 + t%particle(Hneut)%dM(k)**2.0D0)   
  endif
   if(abs(t%particle(Hneut)%M(k)-mobs).le.allowed_massrange ) then
    Higgs_signal_k = signalrate(k,p,d,mobs,t%particle(Hneut)%M(k),t%particle(Hneut)%dM(k))
    LHCrun1_rates(i)%r_pred = LHCrun1_rates(i)%r_pred + Higgs_signal_k
    if(id.eq.11) then ! gg -> h_k -> gaga weighted mass average
     num1 = num1 + Higgs_signal_k * t%particle(Hneut)%M(k)
     dnum1 = dnum1 + Higgs_signal_k * t%particle(Hneut)%dM(k)     
    else if(id.eq.13) then ! gg -> h_k -> ZZ -> 4l weighted mass average
     num2 = num2 + Higgs_signal_k * t%particle(Hneut)%M(k)    
     dnum2 = dnum2 + Higgs_signal_k * t%particle(Hneut)%dM(k)    
    endif       
    ncomb = ncomb+1
   endif
  enddo
  
  if(id.eq.11) then
   denom1 = LHCrun1_rates(i)%r_pred
  else if(id.eq.13) then
   denom2 = LHCrun1_rates(i)%r_pred   
  endif
   
  if(LHCrun1_rates(i)%r_pred.gt.LHCrun1_rates(i)%r) then
   LHCrun1_rates(i)%dr = LHCrun1_rates(i)%dr_up
  else 
   LHCrun1_rates(i)%dr = LHCrun1_rates(i)%dr_low
  endif
  
  if(LHCrun1_rates(i)%r.lt.0.0D0) then
   LHCrun1_rates(i)%dr0 = LHCrun1_rates(i)%dr_up
  else 
   LHCrun1_rates(i)%dr0 = LHCrun1_rates(i)%dr_low
  endif
  
  
   v(i) = LHCrun1_rates(i)%r_pred - LHCrun1_rates(i)%r
   vmat(i,1) = v(i)
   
   vzero(i) = LHCrun1_rates(i)%r
   vzeromat(i,1) = vzero(i)
!   write(*,'(2I3,3F10.5)') p, d, LHCrun1_rates(i)%r_pred, LHCrun1_rates(i)%r, LHCrun1_rates(i)%r/LHCrun1_rates(i)%r_pred 
 enddo
 
 if(denom1.gt.vsmall.and.denom2.gt.vsmall) then
  mav = 0.5D0 * (num1/denom1 + num2/denom2) 
  dmav = 0.5D0 * (dnum1/denom1 + dnum2/denom2) 
!   write(*,*) "Averaged mass is ",mav, " +- ",dmav
!  else
!    write(*,*) "denom1 and denom2 are ",denom1, denom2
 endif 
  
  do i=lbound(LHCrun1_rates,dim=1),ubound(LHCrun1_rates,dim=1)
   do j=lbound(LHCrun1_rates,dim=1),ubound(LHCrun1_rates,dim=1)
    covmat(i,j) = LHCrun1_correlationmatrix(i,j) * &
&                 LHCrun1_rates(i)%dr * LHCrun1_rates(j)%dr    
    covmatzero(i,j) = LHCrun1_correlationmatrix(i,j) * &
&                 LHCrun1_rates(i)%dr0 * LHCrun1_rates(j)%dr0
   enddo
  enddo  

 call invmatrix(covmat, invcovmat)
 call matmult(invcovmat,vmat,v2,20,1)

 call invmatrix(covmatzero, invcovmatzero)
 call matmult(invcovmatzero,vzeromat,vzero2,20,1)


 do i=1, 20
  csq_mu(i) = v(i)*v2(i)
 enddo

 do i=1, 20
  csq_mu_max(i) = vzero(i)*vzero2(i)
 enddo
  
 if(mav.lt.vsmall) then
  HSres(n)%Chisq_LHCRun1_mh=0.0D0
 else 
  HSres(n)%Chisq_LHCRun1_mh=csq_mh(mav,mobs,dmav,dmobs)
 endif
 
 if((HSres(n)%Chisq_LHCRun1_mh+sum(csq_mu)).gt.sum(csq_mu_max)) then
  HSres(n)%Chisq_LHCRun1_mu=sum(csq_mu_max)
  HSres(n)%Chisq_LHCRun1_mh=0.0D0  
 else  
  HSres(n)%Chisq_LHCRun1_mu=sum(csq_mu)
 endif
  
 HSres(n)%Chisq_LHCRun1= HSres(n)%Chisq_LHCRun1_mu +  HSres(n)%Chisq_LHCRun1_mh
 HSres(n)%nobs_LHCRun1_mu=20
 HSres(n)%nobs_LHCRun1_mh=1 
 if(HSres(n)%Chisq_LHCRun1.gt.vsmall.and.(HSres(n)%nobs_LHCRun1_mu+HSres(n)%nobs_LHCRun1_mh-Nparam).gt.0) then
  HSres(n)%Pvalue_LHCRun1=1 - gammp(dble(HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_LHCRun1_mh-Nparam)/2,HSres(n)%Chisq_LHCRun1/2)
 endif

 deallocate(covmat,invcovmat)
 deallocate(covmatzero,invcovmatzero)
 
contains 

!------------------------------------------------------------
 function signalrate(k,p,d,mobs,m,dmtheo)
!------------------------------------------------------------ 
 use usefulbits_hs, only : LHC_combination_run1_SMXS_from_paper
  !--------------------------------------external functions
 double precision :: SMCS_lhc8_gg_H,SMCS_lhc8_bb_H,SMCS_lhc8_vbf_H,		&
 &	SMCS_lhc8_HW, SMCS_lhc8_HZ, SMCS_lhc8_ttH, SMBR_Hgamgam,SMBR_HWW,   &
 &  SMBR_HZZ, SMBR_Htautau, SMBR_Hbb, SMBR_HZgam, SMBR_Hcc, SMBR_Hmumu, &
 &  SMBR_Hgg
  double precision, intent(in) :: mobs, m, dmtheo
  integer, intent(in) :: k,p,d
  double precision :: signalrate, production_rate, decay_rate, mass, refmass
  double precision :: production_rate_scalefactor, decay_rate_scalefactor
  
  mass = t%particle(Hneut)%M(k)

! TS (17/10/2018): Take reference mass for SM-normalization at mobs+dmtheo box boundary.
  if(mass.ge.(mobs+dmtheo)) then
   refmass = mobs + dmtheo
  else if(mass.le.(mobs-dmtheo)) then
   refmass = mobs - dmtheo
  else
   refmass = mass
  endif  
!---
  
  if(p.eq.1) then 
   if(LHC_combination_run1_SMXS_from_paper) then
   production_rate=  t%lhc8%XS_gg_hj_ratio(k) * 19.2D0 &
   &               + t%lhc8%XS_bb_hj_ratio(k) * 0.203D0
   else
   production_rate=  t%lhc8%XS_gg_hj_ratio(k) * SMCS_lhc8_gg_H(mass) &
   &               + t%lhc8%XS_bb_hj_ratio(k) * SMCS_lhc8_bb_H(mass)   
   endif
! NOTE: Here we make a small error in the scalefactor. Correct would be to rescale
!       the gg and bb contributions separately.   
   production_rate_scalefactor = (SMCS_lhc8_gg_H(mobs)+SMCS_lhc8_bb_H(mobs))/&
   &                             (SMCS_lhc8_gg_H(refmass)+SMCS_lhc8_bb_H(refmass))
  else if(p.eq.2) then
   if(LHC_combination_run1_SMXS_from_paper) then
   production_rate=  t%lhc8%XS_vbf_ratio(k) * 1.58D0
   else
   production_rate=  t%lhc8%XS_vbf_ratio(k) * SMCS_lhc8_vbf_H(mass)   
   endif      
   production_rate_scalefactor = SMCS_lhc8_vbf_H(mobs)/SMCS_lhc8_vbf_H(refmass)   
  else if(p.eq.3) then
   if(LHC_combination_run1_SMXS_from_paper) then  
   production_rate=  t%lhc8%XS_hjW_ratio(k) * 0.703D0
   else
   production_rate=  t%lhc8%XS_hjW_ratio(k) * SMCS_lhc8_HW(mass)     
   endif
   production_rate_scalefactor = SMCS_lhc8_HW(mobs)/SMCS_lhc8_HW(refmass)      
  else if(p.eq.4) then
   if(LHC_combination_run1_SMXS_from_paper) then
   production_rate=  t%lhc8%XS_hjZ_ratio(k) * 0.446D0
   else
   production_rate=  t%lhc8%XS_hjZ_ratio(k) * SMCS_lhc8_HZ(mass)
   endif 
   production_rate_scalefactor = SMCS_lhc8_HZ(mobs)/SMCS_lhc8_HZ(refmass)         
  else if(p.eq.5) then
   if(LHC_combination_run1_SMXS_from_paper) then  
   production_rate=  t%lhc8%XS_tthj_ratio(k) * 0.129D0
   else
   production_rate=  t%lhc8%XS_tthj_ratio(k) * SMCS_lhc8_ttH(mass)
   endif
   production_rate_scalefactor = SMCS_lhc8_ttH(mobs)/SMCS_lhc8_ttH(refmass)         
  endif 
  if(d.eq.1) then
   decay_rate = t%BR_hjgaga(k)
   decay_rate_scalefactor = SMBR_Hgamgam(mobs)/SMBR_Hgamgam(refmass)
  else if(d.eq.2) then
   decay_rate = t%BR_hjWW(k)
   decay_rate_scalefactor = SMBR_HWW(mobs)/SMBR_HWW(refmass)   
  else if(d.eq.3) then
   decay_rate = t%BR_hjZZ(k)
   decay_rate_scalefactor = SMBR_HZZ(mobs)/SMBR_HZZ(refmass)      
  else if(d.eq.4) then
   decay_rate = t%BR_hjtautau(k)
   decay_rate_scalefactor = SMBR_Htautau(mobs)/SMBR_Htautau(refmass)      
  else if(d.eq.5) then
   decay_rate = t%BR_hjbb(k)
   decay_rate_scalefactor = SMBR_Hbb(mobs)/SMBR_Hbb(refmass)      
  endif 

  if(normalize_rates_to_reference_position) then
   signalrate = production_rate * decay_rate
  else
! This is the default:  
   signalrate = production_rate * production_rate_scalefactor * &
&               decay_rate * decay_rate_scalefactor
  endif 
  
   if(normalize_rates_to_reference_position_outside_dmtheo) then
    if(abs(mobs-m).ge.dmtheo) then
     signalrate = production_rate * decay_rate
    endif
   endif  
  
 
 end function signalrate

!------------------------------------------------------------
end subroutine evaluate_LHC_Run1_combination
!------------------------------------------------------------
subroutine run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, nobs_STXS, Pvalue_STXS)
!------------------------------------------------------------
 use STXS, only : evaluate_model_for_STXS, get_chisq_from_STXS, print_STXS, &
&                 get_number_of_STXS_observables, STXSlist, print_STXS_to_file
 use usefulbits, only : theo,just_after_run, ndat, vsmall
 use usefulbits_hs, only : HSres, output_level
 use theo_manip, only : HB5_complete_theo
 use numerics, only : gammp

 double precision, intent(out) :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, Pvalue_STXS
 integer, intent(out) :: nobs_STXS
 double precision :: Pvalue
 integer :: nobs_STXS_rates, nobs_STXS_mh, i, n

 call HB5_complete_theo

 Chisq_STXS_mh = 0.0D0

 do n=1, ndat
  
 
 do i=lbound(STXSlist,dim=1), ubound(STXSlist,dim=1)
  call evaluate_model_for_STXS(STXSlist(i),theo(n))
 enddo
 call get_chisq_from_STXS(Chisq_STXS_rates, Pvalue_STXS)
 call get_number_of_STXS_observables(nobs_STXS_rates, nobs_STXS_mh)
 
 nobs_STXS = nobs_STXS_rates + nobs_STXS_mh
 ! Add routine for possible mh-observable in STXS here!

 Chisq_STXS = Chisq_STXS_rates + Chisq_STXS_mh
 
 HSres(n)%Chisq_STXS_rates = Chisq_STXS_rates
 HSres(n)%Chisq_STXS_mh = Chisq_STXS_mh
 HSres(n)%Chisq_STXS = Chisq_STXS
 HSres(n)%nobs_STXS_rates = nobs_STXS_rates
 HSres(n)%nobs_STXS_mh = nobs_STXS_mh
 HSres(n)%nobs_STXS = nobs_STXS
 
  Pvalue = 1.0D0
  if(Chisq_STXS.gt.vsmall.and.(nobs_STXS-Nparam).gt.0) then
   Pvalue = 1 - gammp(dble(nobs_STXS-Nparam)/2,Chisq_STXS/2)
  endif  

 HSres(n)%Pvalue_STXS = Pvalue
 
 enddo

  if(output_level.eq.1) call print_STXS
  if(output_level.eq.3) then
   call print_STXS_to_file
  endif 


 if(output_level.ne.0) then    
  write(*,*)
  write(*,*) '#*************************************************************************#' 
  write(*,*) '#                HIGGSSIGNALS RESULTS (STXS observables)                  #'
  write(*,*) '#*************************************************************************#' 
  write(*,'(A55,F21.8)') 'chi^2 (signal rate) from STXS observables = ',Chisq_STXS_rates
  write(*,'(A55,F21.8)') 'chi^2 (Higgs mass) from STXS observables = ',Chisq_STXS_mh
  write(*,'(A55,F21.8)') 'chi^2 (total) = ',Chisq_STXS
  write(*,'(A55,I21)') 'Number of STXS rate observables = ', nobs_STXS_rates
  write(*,'(A55,I21)') 'Number of STXS mass observables = ', nobs_STXS_mh
  write(*,'(A55,I21)') 'Number of STXS observables (total) = ', nobs_STXS
  write(*,'(A48,I3,A4,F21.8)') 'Probability (ndf =',nobs-Nparam,') = ', Pvalue
  write(*,*) '#*************************************************************************#' 
  write(*,*)
 endif



end subroutine run_HiggsSignals_STXS
!------------------------------------------------------------------------------------   
subroutine run_HiggsSignals(mode, Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
!------------------------------------------------------------
! This subroutine can be called by the user after HiggsSignals_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsSignals.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental results.
! Arguments (output):
!   * mode = 1,2 or 3 for peak-centered, mass-centered chi^2 method or both, respectively. 
!   * Chisq_mu = total chi^2 contribution from signal strength measurements
!   * Chisq_mh = total chi^2 contribution from Higgs mass measurements
!   * Chisq = total chi^2 value for the combination of the considered Higgs signals
!   * nobs = total number of observables
!   * Pvalue = total chi^2 probability for the agreement between model and data,
!              assuming number of observables == number of degrees of freedom
!    (see manual for more precise definitions))
!------------------------------------------------------------
 use usefulbits, only : theo,just_after_run, inputmethod, ndat!inputsub,
 use usefulbits_HS, only : HSres, runmode, output_level, usescalefactor, Nparam,Exptdir
 use channels, only : check_channels
 use theo_manip, only : HB5_complete_theo!, HB5_recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none               
 integer,intent(in) :: mode 
 !----------------------------------------output
 integer,intent(out) ::           nobs
 double precision,intent(out) ::  Pvalue, Chisq, Chisq_mu, Chisq_mh
 !-------------------------------------internal
 integer :: n,i
 logical :: debug=.False.
 !---------------------------------------------
 
if(mode.eq.1) then
 runmode="peak"
else if(mode.eq.2) then
!  runmode="mass"
 write(*,*) "Warning: The 'mass' method (runmode = 2) is no longer maintained."
 write(*,*) "         The peak-centered chi^2 method will be used instead."
 runmode="peak"
else if(mode.eq.3) then
! runmode="both"
 write(*,*) "Warning: The 'both' method (runmode = 3) is no longer maintained."
 write(*,*) "         The peak-centered chi^2 method will be used instead."
 runmode="peak"
else 
 stop 'Error in subroutine run_HiggsSignals: mode unknown'
endif

 if(.not.allocated(theo))then
  stop 'subroutine HiggsSignals_initialize must be called first'
 endif
!  if(inputmethod.eq.'subrout') then
!   do i=1,ubound(inputsub,dim=1)
!    if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
!    write(*,*) inputsub(i)%req, inputsub(i)%stat
!    write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
!    write(*,*)'should be called once and only once before each call to'
!    write(*,*)'subroutine run_HiggsSignals.'
!    stop 'error in subroutine run_HiggsSignals'
!    endif
! TS: Have to work on this bit to make it run simultaneously with HiggsBounds. Now,
!     commented out the =0 statement. HS thus has to be run before HB.   
!   inputsub(i)%stat=0!now we have used this input, set back to zero   
!   enddo
!  endif 

 if(debug)write(*,*)'manipulating input...'                 ; call flush(6)

 call HB5_complete_theo       
	
 if(debug)write(*,*)'compare each model to the experimental data...' ; call flush(6)                  

 do n=1,ndat

!  call recalculate_theo_for_datapoint(n)

  call evaluate_model(theo(n),n)       
      
  Pvalue  = HSres(n)%Pvalue_peak
  Chisq   = HSres(n)%Chisq_peak
  Chisq_mu   = HSres(n)%Chisq_peak_mu       
  Chisq_mh   = HSres(n)%Chisq_peak_mh         
  nobs = HSres(n)%nobs_peak

 if(output_level.ne.0) then    
  write(*,*)
  write(*,*) '#*************************************************************************#' 
  write(*,*) '#   HIGGSSIGNALS RESULTS (',trim(adjustl(Exptdir)),') -- peak observables                 #'
  write(*,*) '#*************************************************************************#' 
  write(*,'(A55,F21.8)') 'chi^2 (signal strength) from peak observables = ',&
 &  HSres(n)%Chisq_peak_mu
  write(*,'(A55,F21.8)') 'chi^2 (Higgs mass) from peak observables = ',HSres(n)%Chisq_peak_mh
!   write(*,'(A55,F21.8)') 'chi^2 from mass-centered observables = ',HSres(n)%Chisq_mpred
!   write(*,'(A55,F21.8)') 'chi^2 from signal strength peak observables (total) = ',HSres(n)%Chisq_mu
  write(*,'(A55,F21.8)') 'chi^2 (total) from peak observables = ',HSres(n)%Chisq
  write(*,'(A55,I21)') 'Number of signal strength peak observables = ',&
 &  HSres(n)%nobs_peak_mu
  write(*,'(A55,I21)') 'Number of Higgs mass peak observables = ',HSres(n)%nobs_peak_mh
!   write(*,'(A55,I21)') 'Number of mass-centered observables = ',HSres(n)%nobs_mpred
  write(*,'(A55,I21)') 'Number of peak observables (total) = ',HSres(n)%nobs_peak
  write(*,'(A48,I3,A4,F21.8)') 'Probability (ndf =',HSres(n)%nobs-Nparam,') using peak observables = ',HSres(n)%Pvalue_peak
  write(*,*) '#*************************************************************************#' 
  write(*,*)
 endif

 enddo

 just_after_run=.True.
 usescalefactor=.False.
 

end subroutine run_HiggsSignals
!------------------------------------------------------------
subroutine evaluate_model( t , n )
! internal routine
!------------------------------------------------------------
! This subroutine evaluates the signal strength modifier for every Higgs boson and
! considered analysis. It fills a matrix neutHiggs(:,:) of type neutHiggs with dimensions
! (number(considered analyses),nH).
!------------------------------------------------------------
 use usefulbits, only : np,Hneut,Hplus,dataset,results, vsmall
 use usefulbits_hs, only : neutHiggses, nanalys, runmode, HSresults, cov, obs, analyses,&
 &						   cov_mhneut, iterations, deallocate_covariance_matrices, &
 &						   output_level, Nparam, nanalys
 use datatables, only : setup_tablelist, check_available_Higgses
 use pc_chisq
 use mc_chisq
 use all_chisq
 use numerics
 implicit none
 !--------------------------------------input      
 type(dataset), intent(in) :: t   
 integer, intent(in) :: n   
 !-------------------------------------output
!  type(HSresults), intent(out) :: r

 integer :: ii, jj, iii, jjj
 
 double precision :: totchisq, muchisq, mhchisq, mpchisq, mpredchisq
 integer :: nobs, Nmu, Nmh, Nmpred
 character(LEN=100), allocatable :: assignmentgroups(:)
 integer, allocatable :: assignmentgroups_domH(:)
 integer, allocatable :: assignmentgroups_Higgs_comb(:,:)

 allocate(assignmentgroups(nanalys),assignmentgroups_domH(nanalys))
 allocate(assignmentgroups_Higgs_comb(nanalys,np(Hneut))) 

assignmentgroups = ''

!---Initialize assignmentgroups arrays with default values
do ii=lbound(assignmentgroups_domH,dim=1),ubound(assignmentgroups_domH,dim=1)
 assignmentgroups_domH(ii) = 0
 assignmentgroups_Higgs_comb(ii,:) = 0
enddo    
  
!---First, evaluate the model predictions  
 allocate(neutHiggses(nanalys,np(Hneut)))
!-Loop over considered analyses  
 do ii=lbound(neutHiggses,dim=1),ubound(neutHiggses,dim=1)
!-Loop over the neutral Higgs bosons of the model 
  do jj=lbound(neutHiggses,dim=2),ubound(neutHiggses,dim=2)
!!   write(*,*) "hello evaluate model:", ii, jj
   call calc_mupred(jj, t, obs(ii)%table, neutHiggses(ii,jj))
  enddo
  if(.not.allocated(obs(ii)%Higgses)) allocate(obs(ii)%Higgses(np(Hneut)))
  obs(ii)%Higgses(:) = neutHiggses(ii,:)  
 enddo

!-Pass the observables and their predicted Higgs properties (obs%Higgses)
!-to the tablelist "analyses"
 call setup_tablelist
!  select case(runmode)
!  
!  case('peak')
!-Peak-centered chisq method 
  jjj=0
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   call deallocate_covariance_matrices
   call assign_Higgs_to_peaks(analyses(ii)%table, analyses(ii)%peaks,0)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.1.and.&
&      len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))).ne.0) then
     jjj=jjj+1
     assignmentgroups(jjj)=analyses(ii)%peaks(iii)%assignmentgroup
     assignmentgroups_Higgs_comb(jjj,:)=analyses(ii)%peaks(iii)%Higgs_comb
     assignmentgroups_domH(jjj)=analyses(ii)%peaks(iii)%domH     
!     write(*,*) "Found leader of group ",assignmentgroups(jjj)
!     write(*,*) "ID ",analyses(ii)%peaks(iii)%id
!     write(*,*) "with Higgs combination ",assignmentgroups_Higgs_comb(jjj,:)  
!     write(*,*) "and dominant Higgs boson ",assignmentgroups_domH(jjj)          
    endif   
   enddo 
  enddo
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.0.and.&
&     len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))).ne.0) then
      !SELECT ASSIGNMENT GROUP FOLLOWERS
      do jjj=lbound(assignmentgroups,dim=1),ubound(assignmentgroups,dim=1)
       if(analyses(ii)%peaks(iii)%assignmentgroup.eq.assignmentgroups(jjj)) then
        !TAKE OVER THE HIGGS ASSIGNMENT OF THE LEADING PEAK
        analyses(ii)%peaks(iii)%Higgs_comb=assignmentgroups_Higgs_comb(jjj,:)
        analyses(ii)%peaks(iii)%domH=assignmentgroups_domH(jjj)
        if(assignmentgroups_domH(jjj).ne.0) then
         analyses(ii)%peaks(iii)%Higgs_assignment_forced=1
        endif 
        call evaluate_peak(analyses(ii)%peaks(iii),analyses(ii)%table)
       endif
      enddo
    endif
   enddo
  enddo     

!   write(*,*) "Starting assignment procedure..."

! Do the iterative Higgs-to-peak-assignment here:
  call assign_Higgs_to_peaks_with_correlations(iterations)
  
!   write(*,*) "Calculating chi^2..."  
  call calculate_total_pc_chisq(totchisq, muchisq, mhchisq, nobs, Nmu, Nmh)
!   write(*,*) "...done." 


  if(output_level.eq.1) call print_peakinformation
  if(output_level.eq.2) call print_peakinformation_essentials
  if(output_level.eq.3) then
   call print_peaks_to_file
   call print_peaks_signal_rates_to_file
  endif 

  call add_peaks_to_HSresults(HSres(n))
  
  
    
  HSres(n)%Chisq_peak=totchisq
  HSres(n)%Chisq_peak_mu = muchisq
  HSres(n)%Chisq_mpred = 0.0D0
  HSres(n)%Chisq_peak_mu=muchisq
  HSres(n)%Chisq_peak_mh=mhchisq  
  HSres(n)%nobs_mpred=0
  HSres(n)%nobs_peak_mu=Nmu
  HSres(n)%nobs_peak_mh=Nmh
  HSres(n)%nanalysis=size(analyses)
  HSres(n)%nobs_peak=nobs
!   
!   if(HSres(n)%Chisq.gt.vsmall.and.(HSres(n)%nobs-Nparam).gt.0) then
!    HSres(n)%Pvalue_peak = 1 - gammp(dble(HSres(n)%nobs-Nparam)/2,HSres(n)%Chisq/2)
!   endif

  if(HSres(n)%Chisq_peak.gt.vsmall.and.(HSres(n)%nobs_peak-Nparam).gt.0) then
   HSres(n)%Pvalue_peak = 1 - gammp(dble(HSres(n)%nobs_peak-Nparam)/2,HSres(n)%Chisq_peak/2)
  endif
  
!  case('mass')  
!   do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
!    call fill_mp_obs(ii)
!   enddo
!   if(mc_mode.eq.1) call mass_variation_by_theory_uncertainty
!   call create_covariance_matrix_mp
!   call calculate_mpred_chisq(mpchisq, nobs)
! 
!   if(output_level.eq.1) call print_mc_observables
!   if(output_level.eq.2) call print_mc_observables_essentials
!   if(output_level.eq.3) then
!    call print_mc_tables_to_file
!    call print_mc_observables_to_file   
!   endif
!   
!   HSres(n)%Chisq=mpchisq
!   HSres(n)%Chisq_peak_mu = 0.0D0
!   HSres(n)%Chisq_mpred = mpchisq 
!   HSres(n)%Chisq_mu=mpchisq
!   HSres(n)%Chisq_mh=0.0D0
!   HSres(n)%nobs_mpred=nobs
!   HSres(n)%nobs_peak_mu=0
!   HSres(n)%nobs_peak_mh=0
!   HSres(n)%nanalysis=size(analyses)  
!   HSres(n)%nobs=nobs    
!   if(HSres(n)%Chisq.gt.vsmall.and.(HSres(n)%nobs-Nparam).gt.0) then
!    HSres(n)%Pvalue=1 - gammp(dble(HSres(n)%nobs-Nparam)/2,HSres(n)%Chisq/2)
!   endif
! 
!  case('both')
!  jjj=0
!   do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
!    call deallocate_covariance_matrices
!    call assign_Higgs_to_peaks(analyses(ii)%table, analyses(ii)%peaks,0)
!    do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
!     if(analyses(ii)%table%mhchisq.eq.1.and.&
! &      len(trim(analyses(ii)%peaks(iii)%assignmentgroup)).ne.0) then
!      jjj=jjj+1
!      assignmentgroups(jjj)=analyses(ii)%peaks(iii)%assignmentgroup
!      assignmentgroups_Higgs_comb(jjj,:)=analyses(ii)%peaks(iii)%Higgs_comb
!      assignmentgroups_domH(jjj)=analyses(ii)%peaks(iii)%domH     
!     endif   
!    enddo    
!   enddo
!   do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
!    do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
!     if(analyses(ii)%table%mhchisq.eq.0.and.&
! &     len(trim(analyses(ii)%peaks(iii)%assignmentgroup)).ne.0) then
!       do jjj=lbound(assignmentgroups,dim=1),ubound(assignmentgroups,dim=1)
!        if(analyses(ii)%peaks(iii)%assignmentgroup.eq.assignmentgroups(jjj)) then
!         !TAKE OVER THE HIGGS ASSIGNMENT OF THE LEADING PEAK
!         analyses(ii)%peaks(iii)%Higgs_comb=assignmentgroups_Higgs_comb(jjj,:)
!         analyses(ii)%peaks(iii)%domH=assignmentgroups_domH(jjj)
!         if(assignmentgroups_domH(jjj).ne.0) then
!          analyses(ii)%peaks(iii)%Higgs_assignment_forced=1
!         endif 
!         ! TODO: Need to evaluate everything else here!
!         call evaluate_peak(analyses(ii)%peaks(iii),analyses(ii)%table)
!        endif
!       enddo
!     endif
!    enddo
!   enddo     
!   
!   call assign_Higgs_to_peaks_with_correlations(iterations) 
!   
!   do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
!    call check_available_Higgses(ii)
!    call fill_mp_obs(ii)
!   enddo  
!   if(mc_mode.eq.1) call mass_variation_by_theory_uncertainty
!  
!   call calculate_total_chisq(totchisq, muchisq, mhchisq, mpredchisq, nobs, Nmu, Nmh, Nmpred)
!  
!  !Have to write a new print method
!   if(output_level.eq.1)  call print_all_observables
!   if(output_level.eq.2) call print_peakinformation_essentials
!   if(output_level.eq.3) then
!    call print_peaks_to_file
!    call print_peaks_signal_rates_to_file
!   endif 
! 
!   call add_peaks_to_HSresults(r)
!  
!   HSres(n)%Chisq=totchisq
!   HSres(n)%Chisq_peak_mu = muchisq
!   HSres(n)%Chisq_mpred = mpredchisq 
!   HSres(n)%Chisq_mu=muchisq + mpredchisq
!   HSres(n)%Chisq_mh=mhchisq  
!   HSres(n)%nobs_mpred=Nmpred
!   HSres(n)%nobs_peak_mu=Nmu
!   HSres(n)%nobs_peak_mh=Nmh
!   HSres(n)%nanalysis=size(analyses)
!   HSres(n)%nobs=nobs
!   if(HSres(n)%Chisq.gt.vsmall.and.(HSres(n)%nobs-Nparam).gt.0) then
!    HSres(n)%Pvalue=1 - gammp(dble(HSres(n)%nobs-Nparam)/2,HSres(n)%Chisq/2)
!   endif
!   
!  case default
!   stop "Error in subroutine evaluate_model: Please specify runmode!"
!   
!  end select

 deallocate(neutHiggses)
 deallocate(assignmentgroups, assignmentgroups_domH, assignmentgroups_Higgs_comb)   
end subroutine evaluate_model
!------------------------------------------------------------
subroutine calc_mupred( j, t, mutab, Higgs )
! internal routine
! Calculates the model-predicted signal strength modifier
!------------------------------------------------------------
 use usefulbits, only : dataset, div, vsmall
 use usefulbits_HS, only : neutHiggs, mutable, useSMtest, eps
 implicit none
 
 integer, intent(in) :: j					! Higgs index
 type(dataset), intent(in) :: t
 type(mutable), intent(inout) :: mutab
 type(neutHiggs), intent(inout) :: Higgs

 integer :: i
 double precision :: c, dcbyc
 integer :: testSMratios 
 logical :: correct_properties

 Higgs%m = t%particle(mutab%particle_x)%M(j)
 Higgs%dm = t%particle(mutab%particle_x)%dM(j)
 Higgs%id = j
  
 call get_channelrates( j, t, mutab )

 correct_properties=.True.

!--Evaluate the predicted signal strength modifier c of the model
 c=0. 
 do i=1,mutab%Nc
!----use a weighted average of the channel rate ratios     
  c=c+mutab%channel_w(i,j)*mutab%channel_mu(i,j)
 enddo

!--Evaluate the deviation of each channel rate ratio to the signal
!--strength modifier c and test SM likeness criterium, if this is
!--activated.
 testSMratios= 1  !passes the SM-like ratios test 
 do i=1,mutab%Nc
  dcbyc=div((mutab%channel_mu(i,j)-c),c,0.0D0,1.0D9)
  if(dcbyc*mutab%channel_w(i,j).gt.eps.and.useSMtest) then
   testSMratios= -1  !fails the SM-like ratios test
  endif     
 enddo

 if(testSMratios.lt.0) correct_properties=.False.
  
 if(correct_properties) then
  Higgs%mu=c
 else
  Higgs%mu=0.0D0
 endif
  
end subroutine calc_mupred
!------------------------------------------------------------
subroutine get_channelrates( j, t, mutab )
! internal routine
!
! This subroutine assignes the rates, weights and systematic rate uncertainty of
! the Higgs boson (j) for the channels considered by the analysis (mutab).
!
! WARNING: if normalize_rates_to_reference_position is true
! The rates are normalized w.r.t. a reference rate at the (peak) mass position.
! This does not work with the mass-centered chi^2 method.
! Also, theoretical mass uncertainties are problematic!
!------------------------------------------------------------
 use usefulbits, only : dataset, div, small
 use usefulbits_HS, only : neutHiggs, mutable, delta_rate, normalize_rates_to_reference_position,&
 & normalize_rates_to_reference_position_outside_dmtheo
 use theory_XS_SM_functions
 use theory_BRfunctions
 
 integer, intent(in) :: j
 type(dataset), intent(in) :: t
 type(mutable), intent(inout) :: mutab


 integer :: i, p, d ! id
 integer :: ii, p1, p2, d1, d2 !id1, id2
 double precision :: rate, SMrate, modelrate, drsq_SM, drsq, dBR, dBRSM,drcov,drcovSM
 
 double precision :: rate_SMref,refmass,BR_SMref!,BR_SMref_mpeak
 double precision :: dynamicalmass, rate_SMdyn

! TS (17/10/2018: dynamicalmass is the default reference mass position for the SM normalization)

 if(size(mutab%mass,dim=1).eq.1) then
  refmass = mutab%mass(1)

! TS (17/10/2018): Take dynamical reference mass for SM-normalization at mobs+dmtheo box boundary.
  if(t%particle(mutab%particle_x)%M(j).ge.(mutab%mass(1)+t%particle(mutab%particle_x)%dM(j))) then
   dynamicalmass = mutab%mass(1) + t%particle(mutab%particle_x)%dM(j)
  else if(t%particle(mutab%particle_x)%M(j).le.(mutab%mass(1)-t%particle(mutab%particle_x)%dM(j))) then
   dynamicalmass = mutab%mass(1) - t%particle(mutab%particle_x)%dM(j)
  else
   dynamicalmass = t%particle(mutab%particle_x)%M(j)
  endif  
  
!   write(*,*) "HS debug, dynamicalmass, refmass = ",dynamicalmass, refmass
!---
  
 else
 ! Only relevant for the mass-centered chi^2 method
  refmass = t%particle(mutab%particle_x)%M(j)
 endif

!write(*,*) 'DEBUG HS: id = ', mutab%id
!write(*,*) 'DEBUG HS, m = ', t%particle(mutab%particle_x)%M(j)

 do i=1,mutab%Nc
!   id = mutab%channel_id(i)
!   p = int((id-modulo(id,10))/dble(10))
!   d = modulo(id,10)
  p = mutab%channel_p_id(i)
  d = mutab%channel_d_id(i)  
!--Do the production rate for the relevant experiment and cms-energy 
  if(mutab%collider.eq.'LHC') then
   if(abs(mutab%energy-7.0D0).le.small) then
    if(p.eq.1) then 
     rate=t%lhc7%XS_hj_ratio(j)
     SMrate=t%lhc7%XS_H_SM(j)
     rate_SMdyn=XS_lhc7_gg_H_SM(dynamicalmass)+XS_lhc7_bb_H_SM(dynamicalmass)
     rate_SMref=XS_lhc7_gg_H_SM(refmass)+XS_lhc7_bb_H_SM(refmass)
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%lhc7%XS_vbf_ratio(j)
     SMrate=t%lhc7%XS_vbf_SM(j)
     rate_SMdyn=XS_lhc7_vbf_SM(dynamicalmass)     
     rate_SMref=XS_lhc7_vbf_SM(refmass)
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%lhc7%XS_hjW_ratio(j)
     SMrate=t%lhc7%XS_HW_SM(j) 
     rate_SMdyn=WH_nnlo(dynamicalmass,'LHC7 ',1.0D0,1.0D0,1.0D0,.True.,.True.)
     rate_SMref=WH_nnlo(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,.True.,.True.)
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%lhc7%XS_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_ggqqbb(dynamicalmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc7%XS_tthj_ratio(j)
     SMrate=t%lhc7%XS_ttH_SM(j)
     rate_SMdyn=XS_lhc7_ttH_SM(dynamicalmass)
     rate_SMref=XS_lhc7_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH'     
    else if(p.eq.6) then
     rate=t%lhc7%XS_gg_hj_ratio(j)
     SMrate=t%lhc7%XS_gg_H_SM(j)
     rate_SMdyn=XS_lhc7_gg_H_SM(dynamicalmass)
     rate_SMref=XS_lhc7_gg_H_SM(refmass)
     mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     rate=t%lhc7%XS_bb_hj_ratio(j)
     SMrate=t%lhc7%XS_bb_H_SM(j)
     rate_SMdyn=XS_lhc7_bb_H_SM(dynamicalmass)     
     rate_SMref=XS_lhc7_bb_H_SM(refmass)
     mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     rate=t%lhc7%XS_thj_tchan_ratio(j)
     SMrate=t%lhc7%XS_tH_tchan_SM(j)
     rate_SMdyn=XS_lhc7_tH_tchan_SM(dynamicalmass)
     rate_SMref=XS_lhc7_tH_tchan_SM(refmass)
     mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     rate=t%lhc7%XS_thj_schan_ratio(j)
     SMrate=t%lhc7%XS_tH_schan_SM(j)
     rate_SMdyn=XS_lhc7_tH_schan_SM(dynamicalmass)
     rate_SMref=XS_lhc7_tH_schan_SM(refmass)
     mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     rate=t%lhc7%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_qq_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_qqbb(dynamicalmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_qqbb(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     rate=t%lhc7%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_gg_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_gg(dynamicalmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_gg(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='gg-HZ' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMdyn=1.0D0     
     rate_SMref=1.0D0
     mutab%channel_description(i,1)='none'
    endif 
   else if(abs(mutab%energy-8.0D0).le.small) then
    if(p.eq.1) then 
     rate=t%lhc8%XS_hj_ratio(j)
     SMrate=t%lhc8%XS_H_SM(j)
     rate_SMdyn=XS_lhc8_gg_H_SM(dynamicalmass)+XS_lhc8_bb_H_SM(dynamicalmass)
     rate_SMref=XS_lhc8_gg_H_SM(refmass)+XS_lhc8_bb_H_SM(refmass)
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%lhc8%XS_vbf_ratio(j)
     SMrate=t%lhc8%XS_vbf_SM(j)
     rate_SMdyn=XS_lhc8_vbf_SM(dynamicalmass)
     rate_SMref=XS_lhc8_vbf_SM(refmass)     
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j) 
     rate_SMdyn=WH_nnlo(dynamicalmass,'LHC8 ',1.0D0,1.0D0,1.0D0,.True.,.True.)
     rate_SMref=WH_nnlo(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,.True.,.True.)     
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_ggqqbb(dynamicalmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     rate_SMdyn=XS_lhc8_ttH_SM(dynamicalmass)
     rate_SMref=XS_lhc8_ttH_SM(refmass)     
     mutab%channel_description(i,1)='ttH'
    else if(p.eq.6) then
     rate=t%lhc8%XS_gg_hj_ratio(j)
     SMrate=t%lhc8%XS_gg_H_SM(j)
     rate_SMdyn=XS_lhc8_gg_H_SM(dynamicalmass)
     rate_SMref=XS_lhc8_gg_H_SM(refmass)
     mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     rate=t%lhc8%XS_bb_hj_ratio(j)
     SMrate=t%lhc8%XS_bb_H_SM(j)
     rate_SMdyn=XS_lhc8_bb_H_SM(dynamicalmass)     
     rate_SMref=XS_lhc8_bb_H_SM(refmass)
     mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     rate=t%lhc8%XS_thj_tchan_ratio(j)
     SMrate=t%lhc8%XS_tH_tchan_SM(j)
     rate_SMdyn=XS_lhc8_tH_tchan_SM(dynamicalmass)
     rate_SMref=XS_lhc8_tH_tchan_SM(refmass)
     mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     rate=t%lhc8%XS_thj_schan_ratio(j)
     SMrate=t%lhc8%XS_tH_schan_SM(j)
     rate_SMdyn=XS_lhc8_tH_schan_SM(dynamicalmass)
     rate_SMref=XS_lhc8_tH_schan_SM(refmass)
     mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     rate=t%lhc8%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_qq_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_qqbb(dynamicalmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_qqbb(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     rate=t%lhc8%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_gg_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_gg(dynamicalmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_gg(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='gg-HZ'       
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMdyn=1.0D0     
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif  
   else if(abs(mutab%energy-13.0D0).le.small) then
    if(p.eq.1) then 
     rate=t%lhc13%XS_hj_ratio(j)
     SMrate=t%lhc13%XS_H_SM(j)
     rate_SMdyn=XS_lhc13_gg_H_SM(dynamicalmass)+XS_lhc13_bb_H_SM(dynamicalmass)
     rate_SMref=XS_lhc13_gg_H_SM(refmass) + XS_lhc13_bb_H_SM(refmass)
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%lhc13%XS_vbf_ratio(j)
     SMrate=t%lhc13%XS_vbf_SM(j)
     rate_SMdyn=XS_lhc13_vbf_SM(dynamicalmass)
     rate_SMref=XS_lhc13_vbf_SM(refmass)     
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%lhc13%XS_hjW_ratio(j)
     SMrate=t%lhc13%XS_HW_SM(j) 
     rate_SMdyn=WH_nnlo(dynamicalmass,'LHC13',1.0D0,1.0D0,1.0D0,.True.,.True.)
     rate_SMref=WH_nnlo(refmass,'LHC13',1.0D0,1.0D0,1.0D0,.True.,.True.)     
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%lhc13%XS_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_ggqqbb(dynamicalmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc13%XS_tthj_ratio(j)
     SMrate=t%lhc13%XS_ttH_SM(j)
     rate_SMdyn=XS_lhc13_ttH_SM(dynamicalmass)
     rate_SMref=XS_lhc13_ttH_SM(refmass)     
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.6) then
     rate=t%lhc13%XS_gg_hj_ratio(j)
     SMrate=t%lhc13%XS_gg_H_SM(j)
     rate_SMdyn=XS_lhc13_gg_H_SM(dynamicalmass)
     rate_SMref=XS_lhc13_gg_H_SM(refmass)
     mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     rate=t%lhc13%XS_bb_hj_ratio(j)
     SMrate=t%lhc13%XS_bb_H_SM(j)
     rate_SMdyn=XS_lhc13_bb_H_SM(dynamicalmass)     
     rate_SMref=XS_lhc13_bb_H_SM(refmass)
     mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     rate=t%lhc13%XS_thj_tchan_ratio(j)
     SMrate=t%lhc13%XS_tH_tchan_SM(j)
     rate_SMdyn=XS_lhc13_tH_tchan_SM(dynamicalmass)
     rate_SMref=XS_lhc13_tH_tchan_SM(refmass)
     mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     rate=t%lhc13%XS_thj_schan_ratio(j)
     SMrate=t%lhc13%XS_tH_schan_SM(j)
     rate_SMdyn=XS_lhc13_tH_schan_SM(dynamicalmass)
     rate_SMref=XS_lhc13_tH_schan_SM(refmass)
     mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     rate=t%lhc13%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_qq_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_qqbb(dynamicalmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_qqbb(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     rate=t%lhc13%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_gg_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_gg(dynamicalmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_gg(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='gg-HZ'     
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMdyn=1.0D0     
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif  
   endif 
  else if(mutab%collider.eq.'TEV') then
    if(p.eq.1) then 
     rate=t%tev%XS_hj_ratio(j)
     SMrate=t%tev%XS_H_SM(j)
     rate_SMdyn=XS_tev_gg_H_SM(dynamicalmass)+XS_tev_bb_H_SM(dynamicalmass)
     rate_SMref=XS_tev_gg_H_SM(refmass)+XS_tev_bb_H_SM(refmass)
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%tev%XS_vbf_ratio(j)
     SMrate=t%tev%XS_vbf_SM(j)
     rate_SMdyn=XS_tev_vbf_SM(dynamicalmass)     
     rate_SMref=XS_tev_vbf_SM(refmass)     
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%tev%XS_hjW_ratio(j)
     SMrate=t%tev%XS_HW_SM(j) 
     rate_SMdyn=WH_nnlo(dynamicalmass,'TEV  ',1.0D0,1.0D0,1.0D0,.True.,.True.)     
     rate_SMref=WH_nnlo(refmass,'TEV  ',1.0D0,1.0D0,1.0D0,.True.,.True.)
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%tev%XS_hjZ_ratio(j)  
     SMrate=t%tev%XS_HZ_SM(j)
     rate_SMdyn=ZH_cpmix_nnlo_ggqqbb(dynamicalmass,'TEV  ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     rate_SMref=ZH_cpmix_nnlo_ggqqbb(refmass,'TEV  ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%tev%XS_tthj_ratio(j)
     SMrate=t%tev%XS_ttH_SM(j)
     rate_SMdyn=XS_tev_ttH_SM(dynamicalmass)          
     rate_SMref=XS_tev_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMdyn=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif       
  else if(mutab%collider.eq.'ILC') then
!--n.B.: As a first attempt, we use the LHC8 normalized cross sections for ZH, VBF, ttH.
!        In order to do this properly, a separate input for the ILC cross sections
!        has to be provided! It works only for single production mode observables (no
!        correct weighting of channels included!)Then, at least in the effective coupling
!        approximation, there is no difference to a full implementation.
!        The theoretical uncertainty of the ILC production modes will are defined in
!        usefulbits_HS.f90.
    if(p.eq.1.or.p.eq.2) then 
     write(*,*) 'Warning: Unknown ILC production mode (',p,') in table ',mutab%id
     rate=0.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='unknown'
    else if(p.eq.3) then
     rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j)
     rate_SMdyn=WH_nnlo(dynamicalmass,'LHC8 ',1.0D0,1.0D0,1.0D0,.True.,.True.)     
     rate_SMref=WH_nnlo(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,.True.,.True.)
     mutab%channel_description(i,1)='WBF'     
    else if(p.eq.4) then
     rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     rate_SMdyn=XS_lhc8_HZ_SM(dynamicalmass)          
     rate_SMref=XS_lhc8_HZ_SM(refmass)     
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     rate_SMdyn=XS_lhc8_ttH_SM(dynamicalmass)          
     rate_SMref=XS_lhc8_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMdyn=1.0D0     
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif           
   endif
   
!    write(*,*) "DEBUG, after production mode, rate_SMdyn = ",rate_SMdyn
!--Multiply now by the decay rate
  if(d.eq.1) then
   rate=rate*div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgaga_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hgaga(dynamicalmass)
   rate_SMref = rate_SMref*BRSM_Hgaga(refmass)
   mutab%channel_description(i,2)='gammagamma'   
  else if(d.eq.2) then
   rate=rate*div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0)   
   SMrate=SMrate*t%BR_HWW_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_HWW(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_HWW(refmass)
   mutab%channel_description(i,2)='WW'   
  else if(d.eq.3) then
   rate=rate*div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZZ_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_HZZ(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_HZZ(refmass)   
   mutab%channel_description(i,2)='ZZ'
  else if(d.eq.4) then
   rate=rate*div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Htautau_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Htautau(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Htautau(refmass)
   mutab%channel_description(i,2)='tautau'
  else if(d.eq.5) then
   rate=rate*div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hbb_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hbb(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Hbb(refmass)
   mutab%channel_description(i,2)='bb'  
  else if(d.eq.6) then
   rate=rate*div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZga_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_HZga(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_HZga(refmass)   
   mutab%channel_description(i,2)='Zgamma'
  else if(d.eq.7) then
   rate=rate*div(t%BR_hjcc(j),t%BR_Hcc_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hcc_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hcc(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Hcc(refmass)
   mutab%channel_description(i,2)='cc'
  else if(d.eq.8) then
   rate=rate*div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hmumu_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hmumu(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Hmumu(refmass)
   mutab%channel_description(i,2)='mumu'
  else if(d.eq.9) then
   rate=rate*div(t%BR_hjgg(j),t%BR_Hgg_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgg_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hgg(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Hgg(refmass)
   mutab%channel_description(i,2)='gg'
  else if(d.eq.10) then
   rate=rate*div(t%BR_hjss(j),t%BR_Hss_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hss_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Hss(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Hss(refmass)
   mutab%channel_description(i,2)='ss'
  else if(d.eq.11) then
   rate=rate*div(t%BR_hjtt(j),t%BR_Htt_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Htt_SM(j)
   rate_SMdyn = rate_SMdyn*BRSM_Htoptop(dynamicalmass)   
   rate_SMref = rate_SMref*BRSM_Htoptop(refmass)
   mutab%channel_description(i,2)='tt'
  else if(d.eq.0) then
   rate=rate*1.0D0
   SMrate=SMrate*1.0D0
   rate_SMdyn = rate_SMdyn*1.0D0
   rate_SMref = rate_SMref*1.0D0
   mutab%channel_description(i,2)='none'
  endif
  
!-------------------------
! NEW FEATURE (since HB-5.2): Enable to set channelrates directly.
  if(p.ne.0.and.d.ne.0) then
   select case(d)
    case(1)
     BR_SMref = t%BR_Hgaga_SM(j)
!      BR_SMref_mpeak = BRSM_Hgaga(refmass)
    case(2)
     BR_SMref = t%BR_HWW_SM(j)
!      BR_SMref_mpeak = BRSM_HWW(refmass)
    case(3)
     BR_SMref = t%BR_HZZ_SM(j)
!      BR_SMref_mpeak = BRSM_HZZ(refmass)
    case(4)
     BR_SMref = t%BR_Htautau_SM(j)
!      BR_SMref_mpeak = BRSM_Htautau(refmass)
    case(5)
     BR_SMref = t%BR_Hbb_SM(j)
!      BR_SMref_mpeak = BRSM_Hbb(refmass)
    case(6)
     BR_SMref = t%BR_HZga_SM(j)
!      BR_SMref_mpeak = BRSM_HZga(refmass)
    case(7)
     BR_SMref = t%BR_Hcc_SM(j)
!      BR_SMref_mpeak = BRSM_Hcc(refmass)
    case(8)
     BR_SMref = t%BR_Hmumu_SM(j)
!      BR_SMref_mpeak = BRSM_Hmumu(refmass)
    case(9)
     BR_SMref = t%BR_Hgg_SM(j)
!      BR_SMref_mpeak = BRSM_Hgg(refmass)
    case(10)
     BR_SMref = t%BR_Hss_SM(j)
    case(11)
     BR_SMref = t%BR_Htt_SM(j)
   end select
   if(mutab%collider.eq.'LHC') then
    if(abs(mutab%energy-7.0D0).le.small) then
     if(t%lhc7%channelrates(j,p,d).ge.0.0d0) then
      rate=div(t%lhc7%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    else if(abs(mutab%energy-8.0D0).le.small) then  
     if(t%lhc8%channelrates(j,p,d).ge.0.0d0) then
      rate=div(t%lhc8%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    else if(abs(mutab%energy-13.0D0).le.small) then  
     if(t%lhc13%channelrates(j,p,d).ge.0.0d0) then
      rate=div(t%lhc13%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    endif 
   else if(mutab%collider.eq.'TEV') then
     if(t%tev%channelrates(j,p,d).ge.0.0d0) then
      rate=div(t%tev%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
   endif
  endif  
!-------------------------
  
! write(*,*) 'DEBUG HS: SM BRs = ', t%BR_HWW_SM(j), t%BR_HZZ_SM(j), t%BR_Hgaga_SM(j)  
! write(*,*) 'DEBUG HS: rate, SMrate(i), rate_SMdyn = ', rate, SMrate, rate_SMdyn
! write(*,*) 'DEBUG HS: eff(i) = ', mutab%channel_eff(i) 
  
  
 if(normalize_rates_to_reference_position) then 
!! THIS IS STILL IN TESTING PHASE !!  
  mutab%channel_mu(i,j)=rate*SMrate/(rate_SMref)
!   write(*,*) "DEBUG HS: normalize_rates_to_reference_position, mu = ",mutab%channel_mu(i,j)
 else
!   mutab%channel_mu(i,j)=rate  !! OLD WAY 
  mutab%channel_mu(i,j)=rate*SMrate/(rate_SMdyn)  !! NEW WAY  (TS 17/10/2018)
!   write(*,*) "DEBUG HS: not normalize_rates_to_reference_position, mu = ",mutab%channel_mu(i,j)  
 endif
 
 
 if(normalize_rates_to_reference_position_outside_dmtheo) then
  if(abs(refmass-t%particle(mutab%particle_x)%M(j)).ge.t%particle(mutab%particle_x)%dM(j)) then
   mutab%channel_mu(i,j)=rate*SMrate/(rate_SMref)
   endif
  endif  
 
  
  mutab%channel_w(i,j)=mutab%channel_eff(i)*SMrate 
!  mutab%channel_w_corrected_eff(i,j)=mutab%channel_eff_ratios(i)*mutab%channel_eff(i)*SMrate   
 enddo


! write(*,*) 'DEBUG HS: BRs = ', t%BR_hjWW, t%BR_hjZZ, t%BR_hjgaga  
! write(*,*) 'DEBUG HS: LHC13 = ', t%lhc13%XS_hj_ratio, t%lhc13%XS_vbf_ratio, t%lhc13%XS_hjW_ratio,&
!                                 t%lhc13%XS_hjZ_ratio, t%lhc13%XS_tthj_ratio


 SMrate=sum(mutab%channel_w(:,j))
! write(*,*) 'DEBUG HS: SMrate = ', SMrate
! modelrate=sum(mutab%channel_w_corrected_eff(:,j))
 
 do i=1,mutab%Nc
  mutab%channel_w(i,j)=div(mutab%channel_w(i,j),SMrate,0.0D0,1.0D9)
!  mutab%channel_w_corrected_eff(i,j)=div(mutab%channel_w_corrected_eff(i,j),modelrate,0.0D0,1.0D9)
 enddo
 
! (TS 30/10/2013):
! write(*,*) "get_channelrates (mu, w, weff):" 
! write(*,*) mutab%channel_mu
! write(*,*)  mutab%channel_w
! write(*,*) mutab%channel_eff_ratios
 do i=1,mutab%Nc
  mutab%channel_w_corrected_eff(i,j)=mutab%channel_eff_ratios(i)*mutab%channel_w(i,j)
! n.b.: model weights are not normalized to 1!
 enddo
 
! write(*,*) j,mutab%id, "SM         = ", mutab%channel_w(:,j)
! write(*,*) j,mutab%id, "SM effcorr = ",mutab%channel_w_corrected_eff(:,j) 
 
 do i=1,mutab%Nc
  drsq_SM = 0.0D0
  drsq = 0.0D0

!   id1 = mutab%channel_id(i)
!   p1 = int((id1-modulo(id1,10))/dble(10))
!   d1 = modulo(id1,10)
  p1 = mutab%channel_p_id(i)
  d1 = mutab%channel_d_id(i)
  if(mutab%collider.ne.'ILC') then
   do ii=1,mutab%Nc 
    p2 = mutab%channel_p_id(ii)
    d2 = mutab%channel_d_id(ii)
!     id2 = mutab%channel_id(ii)
!     p2 = int((id2-modulo(id2,10))/dble(10))
!     d2 = modulo(id2,10)
    if(p1.eq.p2.and.p1.ne.0) then
     if(delta_rate%CScov_ok.and.delta_rate%usecov) then
!-- TS 29/03/2017: Add 13 TeV XS covariance matrix here     
      if(abs(mutab%energy-13.0D0).le.small) then
       drcov=delta_rate%CS13cov(p1,p1)
       drcovSM=delta_rate%CS13covSM(p1,p1)
      else
        drcov=delta_rate%CScov(p1,p1)
        drcovSM=delta_rate%CScovSM(p1,p1)
      endif 
      drsq=drsq+drcov*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
      drsq_SM=drsq_SM+drcovSM*mutab%channel_w(i,j)*mutab%channel_w(ii,j)     
     else
      drsq=drsq+delta_rate%dCS(p1)**2*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
      drsq_SM=drsq_SM+delta_rate%dCS_SM(p1)**2*mutab%channel_w(i,j)*mutab%channel_w(ii,j)
     endif
    endif 
    if(d1.eq.d2.and.d1.ne.0) then
     if(delta_rate%BRcov_ok.and.delta_rate%usecov) then
      dBRSM = delta_rate%BRcovSM(d1,d1)
      dBR = delta_rate%BRcov(d1,d1)
     else
      dBRSM = delta_rate%dBR_SM(d1)**2
      dBR = delta_rate%dBR(d1)**2
     endif 
     drsq=drsq+dBR*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
     drsq_SM=drsq_SM+dBRSM*mutab%channel_w(i,j)*mutab%channel_w(ii,j)   
    endif 
   enddo 
  endif 
  mutab%channel_syst(i,j)=sqrt(drsq)
  mutab%channel_systSM(i,j)=sqrt(drsq_SM)
 enddo
   

! write(*,*) 'DEBUG HS: mu = ', mutab%channel_mu  
! write(*,*) 'DEBUG HS: w = ', mutab%channel_w
! write(*,*) 'DEBUG HS: eff = ', mutab%channel_eff
   
end subroutine get_channelrates
!------------------------------------------------------------
subroutine get_Rvalues(ii,collider,R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
! Returns SM normalized signal rates of some relevant channels (w/o efficiencies) 
! for Higgs boson "ii" for a specific collider (see subroutine get_rates).
!------------------------------------------------------------
! use usefulbits, only : theo, np,Hneut
! use usefulbits_HS, only : mutable

 integer, intent(in) :: ii, collider
 double precision, intent(out) :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb
! type(mutable) :: dummytable
! integer :: i
 
 call get_rates(ii,collider,5,(/ 12, 22, 32, 42, 52 /),R_H_WW)
 call get_rates(ii,collider,5,(/ 13, 23, 33, 43, 53 /),R_H_ZZ)
 call get_rates(ii,collider,5,(/ 11, 21, 31, 41, 51 /),R_H_gaga)
 call get_rates(ii,collider,5,(/ 14, 24, 34, 44, 54 /),R_H_tautau)
 call get_rates(ii,collider,5,(/ 15, 25, 35, 45, 55 /),R_H_bb)
 call get_rates(ii,collider,2,(/ 35, 45 /),R_VH_bb)

end subroutine get_Rvalues
!************************************************************
subroutine get_rates(ii,collider,Nchannels,IDchannels,rate)
! Returns SM normalized signal rates (w/o efficiencies) for Higgs boson "ii" and collider
! experiment "collider"(=1,2,3 for TEV, LHC7, LHC8). "Nchannels" gives the total number
! and IDchannels the two-digit ID of the subchannels, which should be included in the rates.
! IDchannels is an array of size(Nchannels).
!------------------------------------------------------------
 use usefulbits, only : theo, np,Hneut
 use usefulbits_HS, only : mutable

 integer, intent(in) :: ii, collider, Nchannels
 integer, dimension(Nchannels), intent(in) :: IDchannels
 double precision, intent(out) :: rate
!-Internal 
 type(mutable) :: dummytable
 integer :: i
 
!-Initialize a dummy mutable in order to run get_channelrates for the channels we want.  
 if(collider.eq.1) then
  dummytable%collider = 'TEV'
 else if(collider.eq.2) then
  dummytable%collider = 'LHC'
  dummytable%energy = 7.0D0
 else if(collider.eq.3) then
  dummytable%collider = 'LHC'
  dummytable%energy = 8.0D0
 else if(collider.eq.4) then
  dummytable%collider = 'LHC'
  dummytable%energy = 13.0D0
 else
  write(*,*) 'WARNING: collider experiment for get_rates unknown.'
  continue
 endif 

 dummytable%id = 999999
 dummytable%particle_x = 1
 dummytable%Nc=Nchannels
 allocate(dummytable%mass(1))
 dummytable%mass(1) = theo(1)%particle(dummytable%particle_x)%M(ii)
!  allocate(dummytable%channel_id(Nchannels))
 allocate(dummytable%channel_p_id(Nchannels))
 allocate(dummytable%channel_d_id(Nchannels)) 
 allocate(dummytable%channel_eff(Nchannels))
 allocate(dummytable%channel_eff_ratios(Nchannels))
!-Set all efficiencies equal: 
 dummytable%channel_eff = 1.0D0
 dummytable%channel_eff_ratios = 1.0D0
 allocate(dummytable%channel_description(Nchannels,2))
 allocate(dummytable%channel_w(Nchannels,np(Hneut)))
 allocate(dummytable%channel_w_corrected_eff(Nchannels,np(Hneut)))
 allocate(dummytable%channel_systSM(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_syst(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_mu(Nchannels,np(Hneut)))
  
 do i=1,Nchannels
  if(IDchannels(i).le.99) then 
   dummytable%channel_p_id(i) = int((IDchannels(i)-modulo(id,10))/dble(10))
   dummytable%channel_d_id(i) = modulo(IDchannels(i),10)
  else
   write(*,*) "Error in get_rates: channel-ID not supported. Use get_rates_str instead!" 
  endif
 enddo
 
 call get_channelrates(ii, theo(1), dummytable)
 rate=0.0D0
 do i=lbound(dummytable%channel_mu,dim=1),ubound(dummytable%channel_mu,dim=1)
  rate = rate + dummytable%channel_mu(i,ii)*dummytable%channel_w(i,ii)
 enddo

 deallocate(dummytable%channel_p_id,dummytable%channel_d_id,dummytable%channel_eff,&
&           dummytable%channel_w,dummytable%channel_systSM,dummytable%channel_syst,		&
&           dummytable%channel_mu,dummytable%channel_eff_ratios,dummytable%channel_description, &
&           dummytable%channel_w_corrected_eff,dummytable%mass)

end subroutine get_rates
!************************************************************
subroutine get_rates_str(ii,collider,Nchannels,IDchannels_str,rate)
! Returns SM normalized signal rates (w/o efficiencies) for Higgs boson "ii" and collider
! experiment "collider"(=1,2,3 for TEV, LHC7, LHC8). "Nchannels" gives the total number
! and IDchannels_str the channel ID string of the subchannels, which should be included in the rates.
! IDchannels_str is an array of size(Nchannels).
!------------------------------------------------------------
 use usefulbits, only : theo, np,Hneut
 use usefulbits_HS, only : mutable

 integer, intent(in) :: ii, collider, Nchannels
 character(LEN=*), dimension(Nchannels), intent(in) :: IDchannels_str
 double precision, intent(out) :: rate
!-Internal 
 type(mutable) :: dummytable
 integer :: i,id,posperiod
 
!-Initialize a dummy mutable in order to run get_channelrates for the channels we want.  
 if(collider.eq.1) then
  dummytable%collider = 'TEV'
 else if(collider.eq.2) then
  dummytable%collider = 'LHC'
  dummytable%energy = 7.0D0
 else if(collider.eq.3) then
  dummytable%collider = 'LHC'
  dummytable%energy = 8.0D0
 else if(collider.eq.4) then
  dummytable%collider = 'LHC'
  dummytable%energy = 13.0D0
 else
  write(*,*) 'WARNING: collider experiment for get_rates unknown.'
  continue
 endif 

 dummytable%id = 999999
 dummytable%particle_x = 1
 dummytable%Nc=Nchannels
 allocate(dummytable%mass(1))
 dummytable%mass(1) = theo(1)%particle(dummytable%particle_x)%M(ii)
!  allocate(dummytable%channel_id(Nchannels))
 allocate(dummytable%channel_p_id(Nchannels))
 allocate(dummytable%channel_d_id(Nchannels)) 
 allocate(dummytable%channel_eff(Nchannels))
 allocate(dummytable%channel_eff_ratios(Nchannels))
!-Set all efficiencies equal: 
 dummytable%channel_eff = 1.0D0
 dummytable%channel_eff_ratios = 1.0D0
 allocate(dummytable%channel_description(Nchannels,2))
 allocate(dummytable%channel_w(Nchannels,np(Hneut)))
 allocate(dummytable%channel_w_corrected_eff(Nchannels,np(Hneut)))
 allocate(dummytable%channel_systSM(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_syst(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_mu(Nchannels,np(Hneut)))


! do i = 1,Nchannels
!  write(*,*) i, IDchannels_str(i)
! enddo
 
 do i=1,Nchannels
  posperiod = index(IDchannels_str(i),'.')
!   write(*,*) IDchannels_str(i)
  if(posperiod.eq.0) then
   if(len(trim(adjustl(IDchannels_str(i)))).eq.2) then
    read(IDchannels_str(i),*) id
    dummytable%channel_p_id(i) = int((id-modulo(id,10))/dble(10))
    dummytable%channel_d_id(i) = modulo(id,10)
   else
    stop " Error in get_rates_str: Cannot handle channel IDs!"
   endif
  else
!    write(*,*) dummytable%channel_p_id(i), dummytable%channel_d_id(i)
   read(IDchannels_str(i)(:posperiod-1),*) dummytable%channel_p_id(i)
   read(IDchannels_str(i)(posperiod+1:),*) dummytable%channel_d_id(i)
  endif
 enddo
 
!  write(*,*) dummytable%channel_p_id, dummytable%channel_d_id
 call get_channelrates(ii, theo(1), dummytable)
 rate=0.0D0
 do i=lbound(dummytable%channel_mu,dim=1),ubound(dummytable%channel_mu,dim=1)
  rate = rate + dummytable%channel_mu(i,ii)*dummytable%channel_w(i,ii)
 enddo

 deallocate(dummytable%channel_p_id,dummytable%channel_d_id,dummytable%channel_eff,&
&           dummytable%channel_w,dummytable%channel_systSM,dummytable%channel_syst,		&
&           dummytable%channel_mu,dummytable%channel_eff_ratios,dummytable%channel_description, &
&           dummytable%channel_w_corrected_eff,dummytable%mass)

end subroutine get_rates_str
!------------------------------------------------------------
subroutine get_Pvalue(nparam, Pvalue)
! Calculates the Chi^2 probability for the total Chi^2 value
! and the number of degrees of freedom given by the
! number of observables - nparam
!------------------------------------------------------------
 use usefulbits, only : vsmall
 use usefulbits_hs, only: HSres
 use numerics 
 implicit none
 integer, intent(in) :: nparam
 double precision, intent(out) :: Pvalue
 
 if(allocated(HSres)) then
  if(HSres(1)%Chisq.gt.vsmall.and.(HSres(1)%nobs-nparam).gt.0) then
   HSres(1)%Pvalue = 1 - gammp(dble(HSres(1)%nobs-nparam)/2,HSres(1)%Chisq/2)
  endif
 else
  write(*,*) "Warning: subroutine get_Pvalue should be called after run_HiggsSignals." 
 endif
 
 Pvalue = HSres(1)%Pvalue

end subroutine get_Pvalue
!------------------------------------------------------------
subroutine get_neutral_Higgs_masses(Mh, dMh)
! debugging/internal routine
! returns the mass and  theoretical mass uncertainty of the Higgs bosons.
!------------------------------------------------------------
 use usefulbits, only: theo,np,Hneut
 
 implicit none
 double precision,intent(out) :: Mh(np(Hneut)), dMh(np(Hneut))

 if(.not.allocated(theo))then
  stop 'No model information given!'
 endif
 if(np(Hneut).eq.0)then
  write(*,*)'Cannot access the neutral Higgs boson masses'
  write(*,*)'because np(Hneut) == 0.'
  stop 'error in subroutine get_neutral_Higgs_masses'
 endif

 Mh = theo(1)%particle(Hneut)%M 
 dMh = theo(1)%particle(Hneut)%dM
 
end subroutine get_neutral_Higgs_masses
!------------------------------------------------------------
subroutine complete_HS_results()
!------------------------------------------------------------
 use usefulbits, only : just_after_run, ndat
 use usefulbits_HS, only : HSres, Nparam
 use numerics, only : gammp
  
integer :: n

if(just_after_run) then
 do n=1,ndat
 

  HSres(n)%Chisq_mu = HSres(n)%Chisq_peak_mu + & !HSres(n)%Chisq_mpred + &
& HSres(n)%Chisq_STXS_rates + HSres(n)%Chisq_LHCRun1_mu

  HSres(n)%Chisq_mh = HSres(n)%Chisq_peak_mh + HSres(n)%Chisq_LHCRun1_mh + &
& HSres(n)%Chisq_STXS_mh

  HSres(n)%Chisq_STXS = HSres(n)%Chisq_STXS_rates + HSres(n)%Chisq_STXS_mh

  HSres(n)%Chisq_peak = HSres(n)%Chisq_peak_mu + HSres(n)%Chisq_peak_mh

  HSres(n)%Chisq_LHCRun1 = HSres(n)%Chisq_LHCRun1_mu + HSres(n)%Chisq_LHCRun1_mh

  HSres(n)%Chisq =   HSres(n)%Chisq_mu + HSres(n)%Chisq_mh

  HSres(n)%nobs_mu = HSres(n)%nobs_peak_mu + &!HSres(n)%nobs_mpred + &
& HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_STXS_rates

  HSres(n)%nobs_mh = HSres(n)%nobs_peak_mh + HSres(n)%nobs_LHCRun1_mh + &
& HSres(n)%nobs_STXS_mh

  HSres(n)%nobs_peak = HSres(n)%nobs_peak_mu + HSres(n)%nobs_peak_mh

  HSres(n)%nobs_STXS = HSres(n)%nobs_STXS_rates + HSres(n)%nobs_STXS_mh

  HSres(n)%nobs_LHCRun1 = HSres(n)%nobs_LHCRun1_mu + HSres(n)%nobs_LHCRun1_mh

  HSres(n)%nobs = HSres(n)%nobs_mu + HSres(n)%nobs_mh

  if(HSres(n)%Chisq.gt.vsmall.and.(HSres(n)%nobs-Nparam).gt.0) then
   HSres(n)%Pvalue=1 - gammp(dble(HSres(n)%nobs-Nparam)/2.0D0,HSres(n)%Chisq/2.0D0)
  endif

  if(HSres(n)%Chisq_peak.gt.vsmall.and.(HSres(n)%nobs_peak-Nparam).gt.0) then
   HSres(n)%Pvalue_peak=1 - gammp(dble(HSres(n)%nobs_peak-Nparam)/2.0D0,HSres(n)%Chisq_peak/2.0D0)
  endif

  if(HSres(n)%Chisq_LHCRun1.gt.vsmall.and.(HSres(n)%nobs_LHCRun1-Nparam).gt.0) then
   HSres(n)%Pvalue_LHCRun1=1 - gammp(dble(HSres(n)%nobs_LHCRun1-Nparam)/2.0D0,HSres(n)%Chisq_LHCRun1/2.0D0)
  endif

  if(HSres(n)%Chisq_STXS.gt.vsmall.and.(HSres(n)%nobs_STXS-Nparam).gt.0) then
   HSres(n)%Pvalue_STXS=1 - gammp(dble(HSres(n)%nobs_STXS-Nparam)/2.0D0,HSres(n)%Chisq_STXS/2.0D0)
  endif

 enddo
 else
  write(*,*) "Warning: complete_HS_results was called but just_after_run is", just_after_run
 endif
 
 
!------------------------------------------------------------
end subroutine complete_HS_results
!------------------------------------------------------------
subroutine finish_HiggsSignals
! This subroutine needs to be called right at the end, to close files
! and deallocate arrays
!------------------------------------------------------------
 use usefulbits, only : deallocate_usefulbits,debug,theo,debug, &!,inputsub
   & file_id_debug1,file_id_debug2
 use S95tables, only : deallocate_Exptranges
 use theory_BRfunctions, only : deallocate_BRSM
 use datatables, only : deallocate_observables
 use usefulbits_HS, only : deallocate_usefulbits_HS, analyses
 use mc_chisq, only : deallocate_mc_observables
 use store_pathname_HS
 
!#if defined(NAGf90Fortran)
! use F90_UNIX_IO, only : flush
!#endif
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(debug) write(*,*)'finishing off...'                      ; call flush(6)
 if(.not.allocated(theo))then
!  stop 'HiggsBounds_initialize  should be called first'
 if(debug) write(*,*) "HiggsBounds/HiggsSignals internal structure already deallocated!"
 else
  call deallocate_BRSM
  call deallocate_Exptranges
  call deallocate_usefulbits
!   if (allocated(inputsub)) deallocate(inputsub)
 endif
!  write(*,*) "before deallocate mc observables." 
 call deallocate_mc_observables
!  write(*,*) "after deallocate mc observables."  
 call deallocate_observables
 if(allocated(analyses)) deallocate(analyses)
 call deallocate_usefulbits_HS
! call system('rm -f '//trim(adjustl(pathname_HS))//'Expt_tables/analyses.txt')
 call system('rm -f HS_analyses.txt')
  
 if(debug) write(*,*)'finished'                              ; call flush(6)
 
end subroutine finish_HiggsSignals
!------------------------------------------------------------
subroutine finish_HiggsSignals_only
!------------------------------------------------------------
 use datatables, only : deallocate_observables
 use usefulbits_HS, only : deallocate_usefulbits_HS, analyses
 use mc_chisq, only : deallocate_mc_observables
 use store_pathname_HS
 
 call deallocate_mc_observables
 call deallocate_observables
 if(allocated(analyses)) deallocate(analyses)
 call deallocate_usefulbits_HS
 call system('rm -f HS_analyses.txt')
   
end subroutine finish_HiggsSignals_only
!------------------------------------------------------------
! SOME HANDY WRAPPER SUBROUTINES
!------------------------------------------------------------
subroutine initialize_HiggsSignals_for_Fittino(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals with the experimental
! dataset "latestresults", avoiding to specify this via a string argument.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
! character(LEN=19) :: Expt_string
 character(LEN=33) :: Expt_string 
! Expt_string = "Moriond2013_Fittino"
 Expt_string = "latestresults_April2013_inclusive"
 
 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_for_Fittino
!------------------------------------------------------------
subroutine get_number_of_observables_wrapper(ntotal, npeakmu, npeakmh, nmpred, nanalyses)
!------------------------------------------------------------
 use io, only : get_number_of_observables
 
 implicit none
 integer, intent(out) :: ntotal, npeakmu, npeakmh, nmpred, nanalyses
 
 call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses) 
end subroutine get_number_of_observables_wrapper
!------------------------------------------------------------
subroutine get_ID_of_peakobservable_wrapper(ii, ID)
!------------------------------------------------------------
 use io, only : get_ID_of_peakobservable 
 implicit none
 integer, intent(in) :: ii
 integer, intent(out) :: ID

 call get_ID_of_peakobservable(ii, ID)
end subroutine get_ID_of_peakobservable_wrapper
!------------------------------------------------------------
subroutine get_peakinfo_from_HSresults_wrapper(obsID, mupred, domH, nHcomb)
!--------------------------------------------------------------------
 use io, only : get_peakinfo_from_HSresults
 
 implicit none
 integer, intent(in) :: obsID
 double precision, intent(out) :: mupred
 integer, intent(out) :: domH, nHcomb

 call get_peakinfo_from_HSresults(obsID, mupred, domH, nHcomb)
end subroutine get_peakinfo_from_HSresults_wrapper
!------------------------------------------------------------
subroutine print_cov_mh_to_file_wrapper(Hindex)
!------------------------------------------------------------
 use pc_chisq, only : print_cov_mh_to_file

 implicit none
 integer, intent(in) :: Hindex
 
 call print_cov_mh_to_file(Hindex)
end subroutine print_cov_mh_to_file_wrapper
!------------------------------------------------------------
subroutine print_cov_mu_to_file_wrapper
!------------------------------------------------------------
 use pc_chisq, only : print_cov_mu_to_file

 implicit none
 call print_cov_mu_to_file
end subroutine print_cov_mu_to_file_wrapper
!------------------------------------------------------------
subroutine print_corr_mu_to_file_wrapper
!------------------------------------------------------------
 use pc_chisq, only : print_corr_mu_to_file

 implicit none
 call print_corr_mu_to_file
end subroutine print_corr_mu_to_file_wrapper
!------------------------------------------------------------
