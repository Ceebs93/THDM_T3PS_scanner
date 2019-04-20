! This file is part of HiggsBounds
!  -KW
!************************************************************      
subroutine initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses_in)
! This the first Higgsbounds subroutine that should be called
! by the user.
! It calls subroutines to read in the tables of Standard Model data,
! read in the tables of LEP, Tevatron and LHC data,
! set up lists of processes which should be checked against 
! the experimental results, allocate arrays etc
! Arguments (input):
!   * nHiggs= number of neutral Higgs in the model 
!     (see subroutine check_nH_nHplus in input.f90 for more details)
!   * nHiggsplus= number of singly,positively charged Higgs in the model 
!     (see subroutine check_nH_nHplus in input.f90 for more details)
!   * whichanalyses_in= which combination of experimental results to use
!     (see subroutine check_whichanalyses in input.f90 for more details)
!************************************************************
 use usefulbits, only : np,Hneut,Hplus,Chineut,Chiplus,debug,inputmethod,  &
  &   theo,whichanalyses,HiggsBounds_info,just_after_run,BRdirectinput,&
  &   file_id_debug1,file_id_debug2,allocate_if_stats_required,run_HB_classic! ,inputsub
 use input, only : setup_input,check_number_of_particles,check_whichanalyses
 use S95tables, only : setup_S95tables,S95_t2
 use likelihoods, only : setup_likelihoods
 use theory_BRfunctions, only : setup_BRSM
 use theory_XS_SM_functions, only : setup_XSSM 
 use channels, only : setup_channels
 use output, only : setup_output
#ifdef enableCHISQ
 use S95tables_type3, only : clsb_t3,fillt3needs_M2_gt_2M1
#endif

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

!#define FORFITTINO

 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
! integer,intent(in),optional :: nHiggsplus
! character(LEN=5),intent(in),optional :: whichanalyses_in
 integer,intent(in) :: nHiggsplus
 character(LEN=5),intent(in) :: whichanalyses_in
 !-----------------------------------internal
 integer :: i
 logical :: messages 
 !-------------------------------------------

! if((.not.present(nHiggsplus)).or.(.not.present(whichanalyses_in)))then 
   !Actually, this doesn't work as I wanted it to
   !because if initialize_HiggsBounds is called in the old way, the program
   !usually just crashes..... but leaving it in for now, in case
   !some compilers accept it
!   call attempting_to_use_an_old_HB_version('init')
! endif

#ifdef FORFITTINO
 write(*,*)'The arguments passed to initialize_HiggsBounds are:'
 write(*,*)'nHiggsneut=',nHiggsneut
 write(*,*)'nHiggsplus=',nHiggsplus
 write(*,*)'whichanalyses_in=','~'//trim(adjustl(whichanalyses_in))//'~'
#endif

#ifdef DEBUGGING
  debug=.True.
#else
  debug=.False.
#endif

 messages=debug.or.(inputmethod=='datfile')

! inputmethod='subrout' !('datfile' or 'website' are also possible, but not here)

 np(Hneut)=nHiggsneut
 np(Hplus)=nHiggsplus
 
 np(Chineut)=0! do not change this without contacting us first!
 np(Chiplus)=0! do not change this without contacting us first!

 whichanalyses=whichanalyses_in
 
 if(inputmethod=='subrout') then
  if(allocated(theo))then
   stop 'subroutine HiggsBounds_initialize has already been called once'
  endif

  if(messages)write(*,*)'doing other preliminary tasks...'      ; call flush(6)
  call setup_input

!   allocate(inputsub( 4 )) !(1)np(Hneut)>0 (2)np(Hplus)>0 (3)np(Chineut)>0 (4)np(Chineut)>0 and np(Chiplus)>0
   !                                                                           |        np
   !                                                                           |Hneu Hcha Chineut Chiplus
   !                                                                           | ==0  ==0  ==0     ==0 
!   inputsub(1)%desc='HiggsBounds_neutral_input_*'
!   inputsub(1)%req=req(   0,   1,   1,   1)
!   inputsub(2)%desc='HiggsBounds_charged_input'
!   inputsub(2)%req=req(   1,   0,   1,   1)
!   inputsub(3)%desc='SUSYBounds_neutralinoonly_input'
!   inputsub(3)%req=req(   1,   1,   0,   1)
!   inputsub(4)%desc='SUSYBounds_neutralinochargino_input'
!   inputsub(4)%req=req(   1,   1,   0,   0)
!   do i=1,ubound(inputsub,dim=1)
!    inputsub(i)%stat=0
!   enddo
  endif 

#ifndef WEBVERSION
 if(inputmethod.ne.'datfile') call HiggsBounds_info
 if (run_HB_classic.EQV..True.) then
  PRINT *, "run_HB_classic=True - HiggsBounds is running in classic mode"
 endif
  
#endif

 if(messages)write(*,*)'reading in Standard Model tables...'   ; call flush(6) 
 call setup_BRSM 

 call setup_XSSM
             
 if(messages)write(*,*)'reading in S95tables...'               ; call flush(6)
 call setup_S95tables                                    

 if(messages)write(*,*)'reading in likelihoods...'   ; call flush(6)
 call setup_likelihoods

 if(messages)then
  open(file_id_debug2,file='debug_predratio.txt')
  open(file_id_debug1,file='debug_channels.txt')
 endif

 if(messages)write(*,*)'sorting out processes to be checked...'; call flush(6)
 call setup_channels      
      
 if(messages)write(*,*)'preparing output arrays...'            ; call flush(6)
 call setup_output

#ifdef enableCHISQ
 if(allocated(allocate_if_stats_required))then
   call fillt3needs_M2_gt_2M1(clsb_t3,S95_t2)
 endif 
#endif

 just_after_run=.False.
 BRdirectinput=.False.
  
!  contains 
!    !         |        np
!    !         |Hneu Hcha Chineut Chiplus
!    !         | ==0  ==0  ==0     ==0 
!  function req(Hneu,Hcha, Chneu,  Chcha)
!   integer, intent(in) ::Hneu,Hcha, Chneu,  Chcha
!   integer :: req
!   
!   req=1
!   if(np(Hneut)==0)  req= Hneu  * req
!   if(np(Hplus)==0)  req= Hcha  * req
!   if(np(Chineut)==0)req= Chneu * req
!   if(np(Chiplus)==0)req= Chcha * req
! 
!  end function req 

end subroutine initialize_HiggsBounds
!************************************************************      


!************************************************************      
! Version of initialize_HiggsBounds which takes an integer as
! the third argument. More useful for library linking to 
! non-Fortran codes.
subroutine initialize_HiggsBounds_int(nHn,nHp,flag)
 
       implicit none

       integer nHn,nHp,flag

    interface
     subroutine initialize_HiggsBounds(nHiggsneut, nHiggsplus, whichanalyses_in)
     integer,intent(in) :: nHiggsneut
     integer,intent(in) :: nHiggsplus
     character(LEN=5),intent(in) :: whichanalyses_in
!      integer,intent(in),optional :: nHiggsplus
!      character(LEN=5),intent(in),optional :: whichanalyses_in

        end subroutine initialize_HiggsBounds
    end interface 

      IF (flag.EQ.1) then
        call initialize_HiggsBounds(nHn,nHp, "onlyL")
      elseif (flag.EQ.2) then
        call initialize_HiggsBounds(nHn,nHp, "onlyH")
      elseif (flag.EQ.3) then
        call initialize_HiggsBounds(nHn,nHp, "LandH")
      elseif (flag.EQ.4) then
        call initialize_HiggsBounds(nHn,nHp, "onlyP")
      else
        stop "Illegal value for flag in call to initialize_HB"
      endif
      

end subroutine
!************************************************************      

!************************************************************      
subroutine attempting_to_use_an_old_HB_version(subroutineid)
 use usefulbits, only : vers
 character(len=4),intent(in) :: subroutineid

 select case(subroutineid)
 case('init')
  write(*,*)'The subroutine initialize_HiggsBounds has been called with the'
  write(*,*)'wrong number of arguments. It should be called as:'
  write(*,*)'initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)'
  write(*,*)
  write(*,*)'Note that in early versions of HiggsBounds (HB 1.*.*)'
  write(*,*)'this subroutine was called as:'
  write(*,*)'initialize_HiggsBounds(nHiggsneut,whichanalyses)'
  write(*,*)
 case('effC','part','hadr')
  write(*,*)'The subroutine run_HiggsBounds_'//subroutineid//' has been discontinued in this'
  write(*,*)'version of HiggsBounds.'
 case default
  stop 'wrong input to subroutine attempting_to_use_an_old_HB_version'
 end select

 write(*,*)'If you have code written for use with HB 1.*.*, you have two choices:'
 write(*,*)
 write(*,*)' (1) You can edit your code, such that it works with this' 
 write(*,*)'     version of HiggsBounds (HB'//trim(adjustl(vers))//').'
 write(*,*)'     This has the advantage that you can test your model against many, many'
 write(*,*)'     more Higgs search limits , including charged Higgs search limits.'
 write(*,*)'     See the updated manual for more information.'
 write(*,*)
 write(*,*)' (2) You can download the most recent HB 1.*.* from the HiggsBounds'
 write(*,*)'     website. This contains the LEP Higgs search limits which are'
 write(*,*)'     generally the most useful when constraining new physics models.'
 write(*,*)'     We will continue to support this code.'

 stop 'Incorrect call to a HiggsBounds subroutine.'

end subroutine attempting_to_use_an_old_HB_version
!************************************************************      
subroutine HiggsBounds_input_SLHA(infile)
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): SLHA filename
!************************************************************
 use usefulbits, only : whichinput,infile1,theo,g2,effC,just_after_run, &
   & np,Hneut,Hplus!     ,inputsub
 use extra_bits_for_SLHA

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 character(len=300),intent(in) :: infile
 !--------------------------------------internal
 integer :: n
 !----------------------------------------------
 
 whichinput='SLHA'

!  if(np(Hneut).gt.0)inputsub(Hneut)%stat=inputsub(Hneut)%stat+1
!  if(np(Hplus).gt.0)inputsub(Hplus)%stat=inputsub(Hplus)%stat+1
 ! note: can't be used for charginos or neutralinos yet 

 n=1  

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif
 
 infile1=infile   

 call getSLHAdata(theo(n),effC(n),infile1)    

 just_after_run=.False.

end subroutine HiggsBounds_input_SLHA
!************************************************************      
!
!   HB5 GENERAL INPUT ROUTINES
!
!************************************************************      
subroutine HiggsBounds_neutral_input_properties(Mh,GammaTotal_hj,CP_value)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,just_after_run!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !----------------------------------------input
 double precision,intent(in) :: Mh(np(Hneut)),GammaTotal_hj(np(Hneut)),CP_value(np(Hneut))
 !--------------------------------------internal
 integer :: n
!  integer :: subtype
 !----------------------------------------------
!  subtype=1
  n=1  
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_mass_width should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_mass_width'
 endif

 theo(n)%particle(Hneut)%M       = Mh
 theo(n)%particle(Hneut)%Mc      = Mh    
 theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
 theo(n)%CP_value                = CP_value
 
 just_after_run=.False. 

end subroutine HiggsBounds_neutral_input_properties
!************************************************************   
subroutine HiggsBounds_neutral_input_effC(                     &  
     &          ghjss_s,ghjss_p,ghjcc_s,ghjcc_p,       &
     &          ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p,       &
     &          ghjmumu_s,ghjmumu_p,                   &
     &          ghjtautau_s,ghjtautau_p,               &
     &          ghjWW,ghjZZ,ghjZga,                    &
     &          ghjgaga,ghjgg,ghjhiZ)!,           &
!      &          BR_hjinvisible,BR_hjhihi_nHbynH)
! New neutral Higgs effective coupling input routine. 
! BR's are set separately.
!************************************************************
 use usefulbits, only : theo,np,Hneut,effC,whichinput,just_after_run!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: &!Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) ),&  
     &          ghjss_s( np(Hneut) ),ghjss_p( np(Hneut) ), &
     &          ghjcc_s( np(Hneut) ),ghjcc_p( np(Hneut) ), &
     &          ghjbb_s( np(Hneut) ),ghjbb_p( np(Hneut) ), &
     &          ghjtt_s( np(Hneut) ),ghjtt_p( np(Hneut) ), &
     &          ghjmumu_s( np(Hneut) ),ghjmumu_p( np(Hneut) ), &
     &          ghjtautau_s( np(Hneut) ),ghjtautau_p( np(Hneut) ), &
     &          ghjWW( np(Hneut) ),ghjZZ( np(Hneut) ),ghjZga( np(Hneut) ), &
     &          ghjgaga( np(Hneut) ),ghjgg( np(Hneut) ), &
     &          ghjhiZ(np(Hneut),np(Hneut))
!     &          BR_hjinvisible( np(Hneut) ),BR_hjhihi_nHbynH(np(Hneut),np(Hneut)) 
 !--------------------------------------internal
 integer :: n
!  integer :: subtype
 !----------------------------------------------
 
 whichinput='effC'
!  subtype=1
 n=1  
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_effC should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_effC'
 endif

!  theo(n)%particle(Hneut)%M       = Mh
!  theo(n)%particle(Hneut)%Mc      = Mh    
!  theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj

 effC(n)%hjss_s              = ghjss_s
 effC(n)%hjss_p              = ghjss_p
 effC(n)%hjcc_s              = ghjcc_s
 effC(n)%hjcc_p              = ghjcc_p
 effC(n)%hjbb_s              = ghjbb_s
 effC(n)%hjbb_p              = ghjbb_p
 effC(n)%hjtt_s              = ghjtt_s
 effC(n)%hjtt_p              = ghjtt_p
 effC(n)%hjmumu_s                = ghjmumu_s
 effC(n)%hjmumu_p                = ghjmumu_p
 effC(n)%hjtautau_s              = ghjtautau_s
 effC(n)%hjtautau_p              = ghjtautau_p
        
 effC(n)%hjWW              = ghjWW
 effC(n)%hjZZ              = ghjZZ  
 effC(n)%hjZga             = ghjZga      
 effC(n)%hjgaga            = ghjgaga
 effC(n)%hjgg              = ghjgg
!  g2(n)%hjggZ             = g2hjggZ

 effC(n)%hjhiZ = ghjhiZ

!  theo(n)%BR_hjinvisible  = BR_hjinvisible 
!  theo(n)%BR_hjhihi       = BR_hjhihi_nHbynH  

!  write(*,*) "HiggsBounds_neutral_input_effC hWW coupling = ",effC(n)%hjWW

 just_after_run=.False. 

end subroutine HiggsBounds_neutral_input_effC
!************************************************************      
subroutine HiggsBounds_neutral_input_SMBR(BR_hjss,BR_hjcc,BR_hjbb,    &
      &                           BR_hjtt,BR_hjmumu,          &
      &                           BR_hjtautau,BR_hjWW,        &
      &                           BR_hjZZ,BR_hjZga,BR_hjgaga, &
      &                           BR_hjgg)
! Input for the SM branching ratios      
!************************************************************
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run,BRdirectinput

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !----------------------------------------input
 double precision,intent(in) :: BR_hjss( np(Hneut) ),BR_hjcc( np(Hneut) ),       &
      &                         BR_hjbb( np(Hneut) ),BR_hjtt( np(Hneut) ),       &
      &                         BR_hjmumu( np(Hneut) ),BR_hjtautau( np(Hneut) ), &
      &                         BR_hjWW( np(Hneut) ),BR_hjZZ( np(Hneut) ),       &
      &                         BR_hjZga( np(Hneut) ),BR_hjgaga( np(Hneut) ),    &
      &                         BR_hjgg( np(Hneut) )
 !-------------------------------------internal
 integer :: n
 !---------------------------------------------

  n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_SMBR should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_SMBR'
 endif
! 
  theo(n)%BR_hjss           = BR_hjss
  theo(n)%BR_hjcc           = BR_hjcc 
  theo(n)%BR_hjbb           = BR_hjbb
  theo(n)%BR_hjtt           = BR_hjtt  
  theo(n)%BR_hjmumu         = BR_hjmumu
  theo(n)%BR_hjtautau       = BR_hjtautau                 
  theo(n)%BR_hjWW           = BR_hjWW
  theo(n)%BR_hjZZ           = BR_hjZZ
  theo(n)%BR_hjZga          = BR_hjZga 
  theo(n)%BR_hjgaga         = BR_hjgaga
  theo(n)%BR_hjgg           = BR_hjgg

 just_after_run=.False.
 BRdirectinput=.True. 
    
end subroutine HiggsBounds_neutral_input_SMBR
!************************************************************   
subroutine HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible,BR_hkhjhi,BR_hjhiZ,&
&                                    BR_hjemu,BR_hjetau,BR_hjmutau,BR_hjHpiW)
! Input for the non-SM branching ratios
!************************************************************
 use usefulbits, only : theo,np,Hneut,Hplus,whichinput,just_after_run!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: BR_hjinvisible( np(Hneut) ), &
 &                              BR_hkhjhi(np(Hneut),np(Hneut),np(Hneut)), &
 &                              BR_hjhiZ(np(Hneut),np(Hneut)), &
 &                              BR_hjemu(np(Hneut)),&
 &                              BR_hjetau(np(Hneut)),&
 &                              BR_hjmutau(np(Hneut))
 double precision,intent(in) :: BR_hjHpiW(np(Hneut),np(Hplus)) 
 !--------------------------------------internal
 integer :: n
 n=1  

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_nonSMBR should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_nonSMBR'
 endif

  theo(n)%BR_hjinvisible  = BR_hjinvisible 
  theo(n)%BR_hkhjhi       = BR_hkhjhi
  theo(n)%BR_hjhiZ        = BR_hjhiZ
  theo(n)%BR_hjemu        = BR_hjemu
  theo(n)%BR_hjetau       = BR_hjetau
  theo(n)%BR_hjmutau      = BR_hjmutau

!  write(*,*) "HiggsBounds_neutral_input_nonSMBR"
!  write(*,*) theo(n)%BR_hjHpiW

!  if(present(BR_hjHpiW)) then
  theo(n)%BR_hjHpiW     = BR_hjHpiW
!  endif

  
 just_after_run=.False. 

end subroutine HiggsBounds_neutral_input_nonSMBR
!************************************************************      
subroutine HiggsBounds_neutral_input_LEP(XS_ee_hjZ_ratio,XS_ee_bbhj_ratio, &
                                         XS_ee_tautauhj_ratio,XS_ee_hjhi_ratio)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run
 implicit none
 !---------------------------------------------
 double precision, intent(in) :: XS_ee_hjZ_ratio(np(Hneut)),&
 XS_ee_bbhj_ratio(np(Hneut)),XS_ee_tautauhj_ratio(np(Hneut)),&
 XS_ee_hjhi_ratio(np(Hneut),np(Hneut))
 !-------------------------------------internal
 integer :: n
 !--------------------------------------------- 
 whichinput='hadr' ! What if effC otherwise used?
 n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_LEP should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_LEP'
 endif

 theo(n)%lep%XS_hjZ_ratio = XS_ee_hjZ_ratio
 theo(n)%lep%XS_bbhj_ratio = XS_ee_bbhj_ratio
 theo(n)%lep%XS_tautauhj_ratio = XS_ee_tautauhj_ratio
 theo(n)%lep%XS_hjhi_ratio = XS_ee_hjhi_ratio
  
 just_after_run=.False.

end subroutine HiggsBounds_neutral_input_LEP
!************************************************************      
subroutine HiggsBounds_neutral_input_hadr(collider,CS_hj_ratio,      &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)
!************************************************************
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run,hadroncolliderdataset

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 double precision,intent(in) :: CS_hj_ratio( np(Hneut) ), &
     &          CS_gg_hj_ratio( np(Hneut) ),CS_bb_hj_ratio( np(Hneut) ), &
     &          CS_hjW_ratio( np(Hneut) ) ,CS_hjZ_ratio( np(Hneut) ),    &
     &          CS_vbf_ratio( np(Hneut) ) ,CS_tthj_ratio( np(Hneut) ),    &
     &          CS_thj_tchan_ratio( np(Hneut) ),CS_thj_schan_ratio( np(Hneut) ), &
     &          CS_hjhi( np(Hneut), np(Hneut) )
 integer, intent(in) :: collider
 !-------------------------------------internal
 integer :: n
!  type(hadroncolliderdataset) :: dataset
 !---------------------------------------------
 whichinput='hadr'
  n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr'
 endif

 select case(collider)
  case(2)
  call set_input(theo(n)%tev,CS_hj_ratio,                    &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)
  case(7)
  call set_input(theo(n)%lhc7,CS_hj_ratio,                   &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)
  case(8) 
  call set_input(theo(n)%lhc8,CS_hj_ratio,                   &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)
  case(13) 
  call set_input(theo(n)%lhc13,CS_hj_ratio,                  &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr'
 end select

 just_after_run=.False. 
      
 contains
       
  subroutine set_input(dataset,CS_hj_ratio,                  &
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,               &
     &          CS_hjW_ratio,CS_hjZ_ratio,                   &
     &          CS_vbf_ratio,CS_tthj_ratio,                  &
     &          CS_thj_tchan_ratio,CS_thj_schan_ratio,       &
     &          CS_hjhi)

 implicit none
 double precision,intent(in) :: CS_hj_ratio( np(Hneut) ), &
     &          CS_gg_hj_ratio( np(Hneut) ),CS_bb_hj_ratio( np(Hneut) ), &
     &          CS_hjW_ratio( np(Hneut) ) ,CS_hjZ_ratio( np(Hneut) ),    &
     &          CS_vbf_ratio( np(Hneut) ) ,CS_tthj_ratio( np(Hneut) ),    &
     &          CS_thj_tchan_ratio( np(Hneut) ),CS_thj_schan_ratio( np(Hneut) ), &
     &          CS_hjhi( np(Hneut), np(Hneut) )
   type(hadroncolliderdataset) :: dataset
 
   dataset%XS_hj_ratio = CS_hj_ratio
   dataset%XS_gg_hj_ratio = CS_gg_hj_ratio
   dataset%XS_bb_hj_ratio = CS_bb_hj_ratio
   dataset%XS_hjW_ratio = CS_hjW_ratio
   dataset%XS_hjZ_ratio = CS_hjZ_ratio
   dataset%XS_gg_hjZ_ratio = CS_hjZ_ratio ! assume here that the SM-normalized ratio is equal!
   dataset%XS_qq_hjZ_ratio = CS_hjZ_ratio ! assume here that the SM-normalized ratio is equal!  
   dataset%XS_vbf_ratio = CS_vbf_ratio
   dataset%XS_tthj_ratio = CS_tthj_ratio
   dataset%XS_thj_tchan_ratio = CS_thj_tchan_ratio
   dataset%XS_thj_schan_ratio = CS_thj_schan_ratio 
   dataset%XS_hjhi = CS_hjhi 

  end subroutine set_input      
      
end subroutine HiggsBounds_neutral_input_hadr
!************************************************************      
! subroutine HiggsBounds_neutral_input_ZHprod(collider,CS_qq_hjZ_ratio,CS_gg_hjZ_ratio)
!************************************************************

!************************************************************      
subroutine HiggsBounds_neutral_input_hadr_channelrates(collider,channelrates)
! n.b.: Elements of the matrix channelrates with values < 0 will be overwritten
!       by XS times BR using the narrow width approximation.
!************************************************************
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run,hadroncolliderdataset,&
 &                      Nprod,Ndecay

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 double precision,intent(in) :: channelrates(np(Hneut),Nprod,Ndecay)
 integer, intent(in) :: collider
 !-------------------------------------internal
 integer :: n
 !---------------------------------------------
 whichinput='hadr'
  n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr_channelrates should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr'
 endif

 select case(collider)
  case(2)
   theo(n)%tev%channelrates_tmp=channelrates
  case(7)
   theo(n)%lhc7%channelrates_tmp=channelrates  
  case(8) 
   theo(n)%lhc8%channelrates_tmp=channelrates    
  case(13) 
   theo(n)%lhc13%channelrates_tmp=channelrates  
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_channelrates'
 end select

 just_after_run=.False.       
      
end subroutine HiggsBounds_neutral_input_hadr_channelrates
!************************************************************      
subroutine HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, &
     &          CS_ee_HpjHmj_ratio,                        &
     &          BR_tWpb,BR_tHpjb,                           &
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu,BR_Hpjtb,     &
     &          BR_HpjWZ,BR_HpjhiW)
! This subroutine can be called by the user after subroutine 
! initialize_HiggsBounds has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
! HB-5: Extended input by charged Higgs decays to tb, WZ, hiW
!************************************************************
 use usefulbits, only : theo,np,Hplus,Hneut,just_after_run!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: Mhplus( np(Hplus) ),GammaTotal_Hpj( np(Hplus) ), &
     &          CS_ee_HpjHmj_ratio( np(Hplus) ),                              &
     &          BR_tWpb,BR_tHpjb( np(Hplus) ),                                 &
     &          BR_Hpjcs( np(Hplus) ),BR_Hpjcb( np(Hplus) ),BR_Hpjtaunu( np(Hplus) ), &
     &          BR_Hpjtb( np(Hplus) ),BR_HpjWZ( np(Hplus) ) 
 double precision,intent(in) :: BR_HpjhiW(np(Hplus),np(Hneut))
 !--------------------------------------internal
 integer :: n
!  integer :: j
!  integer :: subtype
 !----------------------------------------------
 
 n=1  
!  subtype=2
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hplus).eq.0)then
  write(*,*)'subroutine HiggsBounds_charged_input should'
  write(*,*)'only be called if np(Hplus)>0'
  stop 'error in subroutine HiggsBounds_charged_input'
 endif

 theo(n)%particle(Hplus)%M       = Mhplus 
 theo(n)%particle(Hplus)%Mc      = Mhplus 
 theo(n)%particle(Hplus)%GammaTot= GammaTotal_Hpj

 theo(n)%lep%XS_HpjHmj_ratio = CS_ee_HpjHmj_ratio

 theo(n)%BR_tWpb  = BR_tWpb
 theo(n)%BR_tHpjb = BR_tHpjb 

 theo(n)%BR_Hpjcs     = BR_Hpjcs
 theo(n)%BR_Hpjcb     = BR_Hpjcb
 theo(n)%BR_Hpjtaunu  = BR_Hpjtaunu
 theo(n)%BR_Hpjtb     = BR_Hpjtb
 theo(n)%BR_HpjWZ     = BR_HpjWZ
 theo(n)%BR_HpjhiW     = BR_HpjhiW
!  write(*,*) 'HiggsBounds_charged_input'
!  write(*,*) theo(n)%BR_HpjhiW 
!  if(present(BR_HpjhiW_in)) then
!   write(*,*) "BR_HpjhiW given: ", BR_HpjhiW_in
!  theo(n)%BR_HpjhiW     = BR_HpjhiW_in
!  else
!   if(np(Hneut).gt.0) then
!    theo(n)%BR_HpjhiW = 0.0D0
!   endif 
!  endif
!  write(*,*) theo(n)%BR_HpjhiW

 just_after_run=.False. 

end subroutine HiggsBounds_charged_input
!************************************************************      
subroutine HiggsBounds_charged_input_hadr(collider, CS_Hpjtb, CS_Hpjcb, &
&           CS_Hpjbjet, CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&           CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Hplus,Hneut,just_after_run,hadroncolliderdataset!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: CS_Hpjtb( np(Hplus) ), CS_Hpjcb( np(Hplus) ),&
&                               CS_Hpjbjet( np(Hplus) ), CS_Hpjcjet( np(Hplus) ),& 
&                               CS_Hpjjetjet( np(Hplus) ), &
&                               CS_HpjW( np(Hplus) ), CS_HpjZ( np(Hplus) ),&
&                               CS_vbf_Hpj( np(Hplus) ), CS_HpjHmj( np(Hplus) )
 integer, intent(in) :: collider
 double precision,intent(in) :: CS_Hpjhi( np(Hplus),np(Hneut) )
 !--------------------------------------internal
 integer :: n
 !----------------------------------------------
  n=1  
  
 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hplus).eq.0)then
  write(*,*)'subroutine HiggsBounds_charged_input should'
  write(*,*)'only be called if np(Hplus)>0'
  stop 'error in subroutine HiggsBounds_charged_input'
 endif

 select case(collider)
  case(2)
!    if(present(CS_Hpjhi)) then   
   call set_input(theo(n)%tev,CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, &
&                      CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)
!    else
!     call set_input(theo(n)%tev,CS_Hpjtb, CS_Hpjbjet, CS_HpjW, &
! &                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj)
!    endif
  case(7)
!    if(present(CS_Hpjhi)) then   
   call set_input(theo(n)%lhc7,CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, &
&                      CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)
!    else
!     call set_input(theo(n)%lhc7,CS_Hpjtb, CS_Hpjbjet, CS_HpjW, &
! &                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj)
!    endif
  case(8) 
!    if(present(CS_Hpjhi)) then   
   call set_input(theo(n)%lhc8,CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, &
&                      CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)
!    else
!     call set_input(theo(n)%lhc8,CS_Hpjtb, CS_Hpjbjet, CS_HpjW, &
! &                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj)
!    endif
  case(13) 
!    if(present(CS_Hpjhi)) then   
   call set_input(theo(n)%lhc13,CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, &
&                      CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)
!    else
!     call set_input(theo(n)%lhc13,CS_Hpjtb, CS_Hpjbjet, CS_HpjW, &
! &                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj)
!    endif
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_charged_input_hadr'
 end select

 just_after_run=.False. 

 contains
 
  subroutine set_input(dataset,CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet, &
&                      CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, &
&                      CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi)

 double precision,intent(in) :: CS_Hpjtb( np(Hplus) ), CS_Hpjcb( np(Hplus) ),&
&                               CS_Hpjbjet( np(Hplus) ), CS_Hpjcjet( np(Hplus) ),& 
&                               CS_Hpjjetjet( np(Hplus) ), &
&                               CS_HpjW( np(Hplus) ), CS_HpjZ( np(Hplus) ),&
&                               CS_vbf_Hpj( np(Hplus) ), CS_HpjHmj( np(Hplus) )
   double precision,intent(in) :: CS_Hpjhi( np(Hplus),np(Hneut) )
   type(hadroncolliderdataset) :: dataset

   dataset%XS_Hpjtb = CS_Hpjtb
   dataset%XS_Hpjcb = CS_Hpjcb   
   dataset%XS_Hpjbjet = CS_Hpjbjet 
   dataset%XS_Hpjcjet = CS_Hpjcjet 
   dataset%XS_Hpjjetjet = CS_Hpjjetjet 
   dataset%XS_vbf_Hpj = CS_vbf_Hpj
   dataset%XS_HpjW = CS_HpjW
   dataset%XS_HpjZ = CS_HpjZ 
   dataset%XS_HpjHmj = CS_HpjHmj
!    if(present(CS_Hpjhi)) then
   dataset%XS_Hpjhi = CS_Hpjhi 
!    endif 

  end subroutine set_input
  
end subroutine HiggsBounds_charged_input_hadr

!************************************************************      
subroutine HiggsBounds_get_neutral_hadr_CS(i,collider,&
&           singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)

 use usefulbits, only : theo, np, Hneut, hadroncolliderdataset
 implicit none
 
 integer, intent(in) :: i, collider
 double precision, intent(out) :: singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(i.gt.np(Hneut)) then
  write(*,"(A,I2,A)") 'WARNING: Requested neutral Higgs h',i,' not part of the model!'  
 else
  select case(collider)
   case(2)
    call get_cross_section(theo(1)%tev,i, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)
  case(7)
    call get_cross_section(theo(1)%lhc7,i, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)
  case(8)
    call get_cross_section(theo(1)%lhc8,i, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)
  case(13) 
    call get_cross_section(theo(1)%lhc13,i, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_get_neutral_SMnormalizedCS'
 end select
 endif

 contains
 
  subroutine get_cross_section(dataset,i, singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan)

   integer, intent(in) :: i
   double precision, intent(inout) :: singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan
   type(hadroncolliderdataset) :: dataset

   singleH = dataset%XS_hj_ratio(i)
   ggH = dataset%XS_gg_hj_ratio(i)
   bbH = dataset%XS_bb_hj_ratio(i)
   VBF = dataset%XS_vbf_ratio(i)
   WH = dataset%XS_hjW_ratio(i)
   ZH = dataset%XS_hjZ_ratio(i)
   ttH = dataset%XS_tthj_ratio(i)
   tH_tchan = dataset%XS_thj_tchan_ratio(i)
   tH_schan = dataset%XS_thj_schan_ratio(i)

  end subroutine get_cross_section
!************************************************************
end subroutine HiggsBounds_get_neutral_hadr_CS
!************************************************************
subroutine HiggsBounds_get_neutral_BR(i,BR_hjss,BR_hjcc,BR_hjbb,&
&           BR_hjtt,BR_hjmumu,BR_hjtautau,BR_hjWW,BR_hjZZ,BR_hjZga,&
&           BR_hjgaga,BR_hjgg)

 use usefulbits, only : theo, np, Hneut
 implicit none
 
 integer, intent(in) :: i
 double precision, intent(out) :: BR_hjss,BR_hjcc,BR_hjbb,&
&           BR_hjtt,BR_hjmumu,BR_hjtautau,BR_hjWW,BR_hjZZ,BR_hjZga,&
&           BR_hjgaga,BR_hjgg

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(i.gt.np(Hneut)) then
  write(*,"(A,I2,A)") 'WARNING: Requested neutral Higgs h',i,' not part of the model!'  
 else
  BR_hjss = theo(1)%BR_hjss(i)
  BR_hjcc = theo(1)%BR_hjcc(i)
  BR_hjbb = theo(1)%BR_hjbb(i)
  BR_hjtt = theo(1)%BR_hjtt(i)
  BR_hjmumu = theo(1)%BR_hjmumu(i)  
  BR_hjtautau = theo(1)%BR_hjtautau(i)
  BR_hjWW = theo(1)%BR_hjWW(i)
  BR_hjZZ = theo(1)%BR_hjZZ(i)  
  BR_hjZga = theo(1)%BR_hjZga(i)  
  BR_hjgaga = theo(1)%BR_hjgaga(i)
  BR_hjgg = theo(1)%BR_hjgg(i)  
 endif

end subroutine HiggsBounds_get_neutral_BR
!************************************************************
subroutine HiggsBounds_set_mass_uncertainties(dMhneut, dMhch)
!************************************************************      
! Assigns the mass uncertainties in the subroutine version.
!
   use usefulbits, only : theo,np,Hneut,Hplus
   
  implicit none
  
  double precision, intent(in) :: dMhneut(np(Hneut))
  double precision, intent(in) :: dMhch(np(Hplus))
    
  theo(1)%particle(Hneut)%dMh = dMhneut
  theo(1)%particle(Hplus)%dMh = dMhch 
  
end subroutine HiggsBounds_set_mass_uncertainties
!************************************************************      
subroutine get_mass_variation_param(n)

  use usefulbits, only : theo,np,Hneut,Hplus,diffMhneut,diffMhch,ndmh,dmhsteps,small_mh
  implicit none

  integer, intent(in) :: n

  double precision :: dMhneut(np(Hneut))
  double precision :: dMhch(np(Hplus))      
  integer :: km(np(Hneut)+np(Hplus))
  integer :: dm(dmhsteps**(np(Hneut)+np(Hplus)),np(Hneut)+np(Hplus))

  integer i,j,k,kp

  if(np(Hneut).gt.0) dMhneut = theo(n)%particle(Hneut)%dMh
  if(np(Hplus).gt.0) dMhch = theo(n)%particle(Hplus)%dMh

  if (modulo(dmhsteps,2).NE.1) then
     stop 'Wrong number of steps in set_mass_uncertainty: must be odd (>=3)'
  endif  

  ndmh = 0
  do i=1,np(Hneut)
   
   IF (dMhneut(i).GT.small_mh) THEN
    ndmh = ndmh + 1
   ENDIF
   km(i)=-(dmhsteps-1)/2
  enddo
   
  do i=1,np(Hplus)
   IF (dMhch(i).GT.small_mh) ndmh = ndmh + 1
    km(i+np(Hneut))=-(dmhsteps-1)/2
  enddo
  
  IF (ndmh.EQ.0) THEN
    RETURN
  ENDIF
    
!  print *, "Number of mass uncertainties: ", ndmh

  if(allocated(diffMhneut)) deallocate(diffMhneut)
  if(allocated(diffMhch)) deallocate(diffMhch)
    
  allocate(diffMhneut(dmhsteps**(np(Hneut)+np(Hplus)),np(Hneut)))
  allocate(diffMhch(dmhsteps**(np(Hneut)+np(Hplus)),np(Hplus)))

  k = 1
  do i=1,dmhsteps**ndmh

    do j=1,ndmh
     dm(i,j) = km(j)
    enddo

    km(k) = km(k)+1

    do j=2,ndmh
     IF (modulo(i,dmhsteps**(j-1)).EQ.0) THEN
      km(j) = km(j)+1
      km(j-1) = -1
     ENDIF
    ENDDO

  enddo    
  
  
  do i=1,dmhsteps**ndmh
   k=1
 
    do j=1,np(Hneut)
     IF (dMhneut(j).GT.small_mh) THEN
      diffMhneut(i,j)=theo(n)%particle(Hneut)%M(j)+dm(i,k)*dMhneut(k)/((dmhsteps-1)/2)
      k = k +1
     ELSE
      diffMhneut(i,j)=theo(n)%particle(Hneut)%M(j)
     ENDIF
    enddo
    
    kp = k
    do j=1,np(Hplus)
     IF (dMhch(j).GT.small_mh) THEN
      diffMhch(i,j)=theo(n)%particle(Hplus)%M(j)+dm(i,k)*dMhch(k-(kp-1))/((dmhsteps-1)/2)
      k = k +1
     ELSE
      diffMhch(i,j)=theo(n)%particle(Hplus)%M(j)
     ENDIF
    enddo
    
  
!    print *, i, (diffMhneut(i,j),j=1,np(Hneut)),(diffMhch(i,j),j=1,np(Hplus))
  
  enddo
  
end subroutine get_mass_variation_param


subroutine SUSYBounds_neutralinoonly_input(MN,GammaTotal_N, &
     &          CS_NjNi,                        &
     &          BR_NjqqNi,BR_NjZNi                           &
     &          )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Chineut,just_after_run!,inputsub,

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: MN( np(Chineut) ),GammaTotal_N( np(Chineut) ) , &
     &          CS_NjNi( np(Chineut),np(Chineut) ),                        &
     &          BR_NjqqNi( np(Chineut),np(Chineut) ),BR_NjZNi( np(Chineut),np(Chineut) )                        
 !--------------------------------------internal
 integer :: n
!  integer :: subtype
 !----------------------------------------------
 
 n=1  
!  subtype=3
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Chineut).eq.0)then
   write(*,*)'subroutine SUSYBounds_neutralinoonly_input should'
   write(*,*)'only be called if np(Chineut)>0'
   stop 'error in SUSYBounds_neutralinoonly_input'
 endif

 theo(n)%particle(Chineut)%M       = MN
 theo(n)%particle(Chineut)%GammaTot= GammaTotal_N

 theo(n)%lep%XS_NjNi = CS_NjNi

 theo(n)%BR_NjqqNi  = BR_NjqqNi
 theo(n)%BR_NjZNi   = BR_NjZNi  

 just_after_run=.False. 

end subroutine SUSYBounds_neutralinoonly_input
!************************************************************      
subroutine SUSYBounds_neutralinochargino_input(MC,GammaTotal_C, &
     &          CS_CpjCmj,                                      &
     &          BR_CjqqNi,                                      &
     &          BR_CjlnuNi,                                     &
     &          BR_CjWNi                                        &
     &          )
! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! has been called.
! Arguments (input): theoretical predictions (see manual for definitions)
!************************************************************
 use usefulbits, only : theo,np,Chineut,Chiplus,just_after_run!,inputsub

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 !----------------------------------------input
 double precision,intent(in) :: MC( np(Chiplus) ),GammaTotal_C( np(Chiplus) ),  &
     &          CS_CpjCmj(  np(Chiplus) ),                                     &
     &          BR_CjqqNi(  np(Chiplus),np(Chineut) ),                         &
     &          BR_CjlnuNi( np(Chiplus),np(Chineut) ),                         &
     &          BR_CjWNi(   np(Chiplus),np(Chineut) )                                        
 !--------------------------------------internal
 integer :: n
 integer :: subtype
 !----------------------------------------------
 
 n=1  
!  subtype=4
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if((np(Chineut).eq.0).or.(np(Chiplus).eq.0))then
   write(*,*)'subroutine SUSYBounds_neutralinochargino_input should'
   write(*,*)'only be called if np(Chineut)>0 and np(Chiplus)>0'
   stop 'error in subroutine SUSYBounds_neutralinochargino_input'
 endif

 theo(n)%particle(Chineut)%M       = MC
 theo(n)%particle(Chineut)%GammaTot= GammaTotal_C

 theo(n)%lep%XS_CpjCmj = CS_CpjCmj

 theo(n)%BR_CjqqNi  = BR_CjqqNi
 theo(n)%BR_CjlnuNi = BR_CjlnuNi
 theo(n)%BR_CjWNi   = BR_CjWNi  

 just_after_run=.False. 

end subroutine SUSYBounds_neutralinochargino_input
!************************************************************


subroutine run_HiggsBounds(HBresult, chan, obsratio, ncombined)
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
! (TS 30/01/2012): Note, that if many data points are tested at the same time (as for
!                  inputmethod==datfiles), this subroutine only returns the results of
!                  the last datapoint. The full results are saved in fullHBres.
 use usefulbits, only : np, Hneut, Hplus, run_HB_classic

 implicit none
 integer HBresult, chan, ncombined
 double precision obsratio
 
 integer hbres(0:np(Hneut)+np(Hplus)), hbchan(0:np(Hneut)+np(Hplus)), hbcomb(0:np(Hneut)+np(Hplus))
 double precision hbobs(0:np(Hneut)+np(Hplus))

! Check if we are using the old 'classic' method
 if (run_HB_classic.EQV..True.) then
   call run_HiggsBounds_classic(HBresult,chan,obsratio,ncombined)
   return
 endif
 
! Call the new ('full') method
 call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)
 
! Combined results are contained in the zero elements of result arrays
 HBresult = hbres(0)
 chan = hbchan(0)
 obsratio = hbobs(0)
 ncombined = hbcomb(0)
    
end subroutine run_HiggsBounds
!************************************************************



subroutine run_HiggsBounds_single(h, HBresult, chan, obsratio, ncombined)
! This subroutine can be used to get the exclusion results
! for a single Higgs boson (specified by the index h).
!
! To obtain individual results from more than one Higgs boson, it
! is more efficient to use run_HiggsBounds_full rather than this method.
 use usefulbits, only : np, Hneut, Hplus

 implicit none
 integer, intent(in) :: h
 integer, intent(out) :: HBresult, chan, ncombined
 double precision, intent(out) :: obsratio
 
 integer hbres(0:np(Hneut)+np(Hplus)), hbchan(0:np(Hneut)+np(Hplus)), hbcomb(0:np(Hneut)+np(Hplus))
 double precision hbobs(0:np(Hneut)+np(Hplus))
 
 IF (h.LT.0) stop "Illegal number of Higgs boson: h < 0"
 if (h.GT.np(Hneut)+np(Hplus)) stop "Illegal number of Higgs boson"
 
 call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)

 HBresult = hbres(h)
 chan = hbchan(h)
 obsratio = hbobs(h)
 ncombined = hbcomb(h)
end subroutine run_HiggsBounds_single
!************************************************************


subroutine run_HiggsBounds_full( HBresult,chan, &
     &                      obsratio, ncombined )
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits.
! 
! The results are given as (n+1)-component arrays (starting from 0), 
! where n is the total number of Higgs bosons in the model (neutral+charged).
! The zeroth component gives the combined results (equivalent to run_HiggsBounds).
!
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
 use usefulbits, only : theo,res,just_after_run,ndmh,debug,numres, &
&                       np,Hneut,Hplus,dmhsteps,ndat,fullHBres,small_mh,&
                        HBresult_all,ncombined_all,chan_all,obsratio_all,predratio_all
 use channels, only : check_channels
 !use input, only : test_input
 use theo_manip, only : HB5_complete_theo, HB5_recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none                      
 !----------------------------------------output
 integer,intent(out)::  HBresult(0:np(Hneut)+np(Hplus))
 integer,intent(out)::  chan(0:np(Hneut)+np(Hplus))
 integer,intent(out)::  ncombined(0:np(Hneut)+np(Hplus))
 double precision,intent(out) :: obsratio(0:np(Hneut)+np(Hplus))
 
 double precision :: Mhneut(np(Hneut))
 double precision :: Mhch(np(Hplus))
 !-------------------------------------internal
 integer :: n,i,j,ind,part,k
 !---------------------------------------------

! print *, "Running HiggsBounds in Normal Mode (most sensitive limit considered for each Higgs boson)"
 
 if (lbound(HBresult,dim=1).NE.0) stop "run_HiggsBounds_full: Array HBresult must begin with element 0"
 if (ubound(HBresult,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array HBresult must be equal to number of Higgses"
 endif 
 if (lbound(chan,dim=1).NE.0) stop "run_HiggsBounds_full: Array chan must begin with element 0"
 if (ubound(chan,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array chan must be equal to number of Higgses"
 endif 
 if (lbound(obsratio,dim=1).NE.0) stop "run_HiggsBounds_full: Array obsratio must begin with element 0"
 if (ubound(obsratio,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array obsratio must be equal to number of Higgses"
 endif 
 if (lbound(ncombined,dim=1).NE.0) stop "run_HiggsBounds_full: Array ncombined must begin with element 0"
 if (ubound(ncombined,dim=1).NE.(np(Hneut)+np(Hplus))) then
  stop "run_HiggsBounds_full: Upper limit of array ncombined must be equal to number of Higgses"
 endif 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(.not.allocated(HBresult_all)) allocate(HBresult_all(0:np(Hneut)+np(Hplus),numres))
 if(.not.allocated(chan_all)) allocate(chan_all(0:np(Hneut)+np(Hplus),numres))
 if(.not.allocated(ncombined_all)) allocate(ncombined_all(0:np(Hneut)+np(Hplus),numres)) 
 if(.not.allocated(obsratio_all)) allocate(obsratio_all(0:np(Hneut)+np(Hplus),numres))  
 if(.not.allocated(predratio_all)) allocate(predratio_all(0:np(Hneut)+np(Hplus),numres))   
!  do i=1,ubound(inputsub,dim=1)
!    if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
!        write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
!        write(*,*)'should be called once and only once before each call to'
!        write(*,*)'subroutine run_HiggsBounds.'
!        stop 'error in subroutine run_HiggsBounds'
!    endif
!    inputsub(i)%stat=0!now we have used this input, set back to zero
!  enddo


 call HB5_complete_theo    
 
  do n=1,ndat

!   if(debug) then
!    write(*,*) "DEBUG BRs: ", theo(n)%BR_hjWW, theo(n)%BR_hjZZ, theo(n)%BR_hjgaga
   
!   endif

        
   theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M
   theo(n)%particle(Hplus)%Mc = theo(n)%particle(Hplus)%M

   call get_mass_variation_param(n)

    do i=0,ubound(Hbresult,dim=1)
     obsratio_all(i,:) = -999d0
     predratio_all(i,:) = -999d0     
     HBresult_all(i,:) = 1
     chan_all(i,:) = -999
     ncombined_all(i,:) = -999
     obsratio(i) = -999d0
     HBresult(i) = 1
     chan(i) = -999
     ncombined(i) = -999
    enddo

!  Do we have mass uncertainties to take care off
    IF (ndmh.GT.0) THEN
   
!    print *, "Running HiggsBounds with Higgs mass uncertainties"
!    write(*,*) theo(n)%particle(Hplus)%dM
      
     if(np(Hneut).ne.0) Mhneut =  theo(n)%particle(Hneut)%M 
     if(np(Hplus).ne.0) Mhch = theo(n)%particle(Hplus)%M

!   Loop over all Higgses
     do i=1,np(Hneut)+np(Hplus)
      obsratio_all(i,:) = 1.D23

      IF (i.LE.np(Hneut)) THEN
       ind = i
       part = Hneut
      ELSE
       ind = i-np(Hneut)
       part = Hplus
      ENDIF
          
!     Check for mass steps for this particular Higgs boson
      IF(theo(n)%particle(part)%dMh(ind).GT.small_mh) THEN

!       theo(n)%particle(part)%M(ind)=theo(n)%particle(part)%M(ind)  &
!       &                             -(dmhsteps-1)/2*theo(n)%particle(part)%dMh(ind)

       theo(n)%particle(part)%M(ind)=theo(n)%particle(part)%M(ind)  &
       &                             -theo(n)%particle(part)%dMh(ind)

       do j=1,dmhsteps 
        
!        print *, theo(n)%particle(Hneut)%M, theo(n)%particle(Hplus)%M

        call HB5_recalculate_theo_for_datapoint(n)       

        call check_channels(theo(n),res(n),i)
     	
     	do k=1,size(res(n)%obsratio)
!       	write(*,*) "i,k,res(n)%obsratio(k),res(n)%chan(k) = ",i,k,res(n)%obsratio(k),res(n)%chan(k)
     	
        IF (res(n)%obsratio(k).LT.obsratio_all(i,k)) THEN
!          write(*,*) "i,k,res(n)%obsratio(k),res(n)%chan(k) = ",i,k,res(n)%obsratio(k),res(n)%chan(k)        
         HBresult_all(i,k)    = res(n)%allowed95(k)
         chan_all(i,k)        = res(n)%chan(k)       
         obsratio_all(i,k)    = res(n)%obsratio(k)
         predratio_all(i,k)    = res(n)%predratio(k)
         ncombined_all(i,k)   = res(n)%ncombined(k)        
        ENDIF
     	enddo
        
!        print *, i,theo(n)%particle(part)%M(ind),HBresult(i),chan(i),obsratio(i),ncombined(i)
        
        theo(n)%particle(part)%M(ind)= theo(n)%particle(part)%M(ind)  &
        &                             +theo(n)%particle(part)%dMh(ind)/(dmhsteps-1)*2
        
      enddo
     else
      call HB5_recalculate_theo_for_datapoint(n)
      call check_channels(theo(n),res(n),i)

      do k=1,size(res(n)%obsratio)
!        write(*,*) "i,k,res(n)%obsratio(k),res(n)%chan(k) = ",i,k,res(n)%obsratio(k),res(n)%chan(k)     	
       HBresult_all(i,k)    = res(n)%allowed95(k)
       chan_all(i,k)        = res(n)%chan(k)       
       obsratio_all(i,k)    = res(n)%obsratio(k)
       predratio_all(i,k)    = res(n)%predratio(k)       
       ncombined_all(i,k)   = res(n)%ncombined(k)

      enddo

!       HBresult(i)    = res(n)%allowed95(1)       
!       chan(i)        = res(n)%chan(1)       
!       obsratio(i)    = res(n)%obsratio(1)
!       ncombined(i)   = res(n)%ncombined(1)            
     endif 
     
     
     HBresult(i) = HBresult_all(i,1)
     chan(i) = chan_all(i,1)
     obsratio(i) = obsratio_all(i,1)
     ncombined(i) = ncombined_all(i,1)
     
!	 Logical OR between exclusions (one Higgs excluded = combined exclusion)
     HBresult(0) = HBresult(0) * HBresult(i)
    
!    Save the data for the Higgs that has the highest ratio of theory/obs 
     IF (obsratio(i).GT.obsratio(0)) THEN
      chan(0)        = chan(i)       
      obsratio(0)    = obsratio(i)
      ncombined(0)   = ncombined(i)        
     ENDIF
     
     theo(n)%particle(Hneut)%M = Mhneut
     theo(n)%particle(Hplus)%M = Mhch
      
    enddo
 !    return 
     
   ELSE
 
 
!    print *, "Running HiggsBounds without Higgs mass uncertainties"

    call HB5_recalculate_theo_for_datapoint(n)       
!    write(*,*) "Higgses = " , np(Hneut)+np(Hplus)

   do i=1,np(Hneut)+np(Hplus)
    call check_channels(theo(n),res(n),i)
    
!    	do k=1,size(res(n)%obsratio)
!    	write(*,*) "i,k,res(n)%obsratio(k),res(n)%chan(k) = ",i,k,res(n)%obsratio(k),res(n)%chan(k)
!    	enddo
    
     do k=1,size(res(n)%obsratio)
!        write(*,*) "i,k,res(n)%obsratio(k),res(n)%chan(k) = ",i,k,res(n)%obsratio(k),res(n)%chan(k)     	
       HBresult_all(i,k)    = res(n)%allowed95(k)
       chan_all(i,k)        = res(n)%chan(k)       
       obsratio_all(i,k)    = res(n)%obsratio(k)
       predratio_all(i,k)   = res(n)%predratio(k)       
       ncombined_all(i,k)   = res(n)%ncombined(k)
     enddo

     HBresult(i) = HBresult_all(i,1)
     chan(i) = chan_all(i,1)
     obsratio(i) = obsratio_all(i,1)
     ncombined(i) = ncombined_all(i,1)

!     HBresult(i)    = res(n)%allowed95(1)
!     chan(i)        = res(n)%chan(1)       
!     obsratio(i)    = res(n)%obsratio(1)
!     ncombined(i)   = res(n)%ncombined(1)        
!     
!    write(*,*) "hello: i=",i," HBres, chan, obsratio = ", HBresult(i), chan(i), obsratio(i)
    
    HBresult(0) = HBresult(0) * res(n)%allowed95(1)
    
    IF (obsratio(i).GT.obsratio(0)) THEN
!     write(*,*) "hello: ", n, i
     chan(0)        = res(n)%chan(1)       
     obsratio(0)    = res(n)%obsratio(1)
     ncombined(0)   = res(n)%ncombined(1)        
    ENDIF
       
!    IF (i.LE.np(Hneut)) THEN
!     print *, i,theo(n)%particle(Hneut)%M(i),HBresult(i),chan(i),obsratio(i),ncombined(i),HBresult(0), obsratio(0)   
!    ELSE
!     print *, i,theo(n)%particle(Hplus)%M(i-np(Hneut)),HBresult(i),chan(i),obsratio(i),ncombined(i),HBresult(0), obsratio(0)   
!    endif
   enddo
  ENDIF

!   write(*,*) "run_HB_full, obsratio: ", obsratio
!   write(*,*) "run_HB_full, chan    : ", chan  

  fullHBres(n)%allowed95=HBresult(0)
  fullHBres(n)%chan=chan(0)
  fullHBres(n)%obsratio=obsratio(0)
  fullHBres(n)%ncombined=ncombined(0)
 
 enddo


 just_after_run=.True.

 
! print *, "HB: run done"

end subroutine run_HiggsBounds_full
!************************************************************
subroutine HiggsBounds_get_most_sensitive_channels_per_Higgs(nH,pos,HBresult,chan,obsratio,predratio,ncombined)
!************************************************************
 use usefulbits, only : HBresult_all,obsratio_all,chan_all,ncombined_all,predratio_all,&
 &                      just_after_run,np,Hneut,Hplus,numres

 integer, intent(in) :: nH, pos
 integer, intent(out) :: HBresult, chan, ncombined
 double precision, intent(out) :: obsratio, predratio

 HBresult = 0
 chan = 0
 obsratio = 0
 predratio = 0
 ncombined = 0
 
 if(just_after_run.and.allocated(HBresult_all)) then
  if(nH.le.np(Hneut)+np(Hplus)) then
   if(pos.le.numres) then
    HBresult  = HBresult_all(nH,pos)
    chan      = chan_all(nH,pos)       
    obsratio  = obsratio_all(nH,pos)
    predratio  = predratio_all(nH,pos)    
    ncombined = ncombined_all(nH,pos)
   else
    write(*,*) 'WARNING: request exceeds the number of stored most sensitive channels (',numres,')'
   endif
  else
   write(*,*) 'WARNING: requested Higgs boson is invalid (choose between 1 and ',np(Hneut)+np(Hplus),'!)'
  endif
 else
  write(*,*) 'WARNING: Please call run_HiggsBounds or run_HiggsBounds_full before calling',&
  &          ' HiggsBounds_get_most_sensitive_channels!'
 endif 

end subroutine HiggsBounds_get_most_sensitive_channels_per_Higgs
!************************************************************
subroutine HiggsBounds_get_most_sensitive_channels(pos,HBresult,chan,obsratio,predratio,ncombined)
!************************************************************
 use usefulbits, only : HBresult_all,obsratio_all,predratio_all,chan_all,ncombined_all,&
 &                      just_after_run,np,Hneut,Hplus,numres

 integer, intent(in) :: pos
 integer, intent(out) :: HBresult, chan, ncombined
 double precision, intent(out) :: obsratio,predratio
 integer :: i,j,count

 integer,allocatable :: nH_rank(:),pos_rank(:), posflat(:) 
 double precision, allocatable :: predratio_tmp(:)
 
 allocate(nH_rank(numres),pos_rank(numres),posflat(numres),predratio_tmp(numres*(np(Hneut)+np(Hplus))))

 HBresult = 0
 chan = 0
 obsratio = 0
 ncombined = 0
 
 predratio_tmp = 0
 
 count=0
 if(just_after_run.and.allocated(HBresult_all)) then
   if(pos.le.numres) then
    do j=1,np(Hneut)+np(Hplus)
     do i=1,numres
      count=count+1
      predratio_tmp(count)=predratio_all(j,i)
     enddo
    enddo
    
    do i=1,numres
     posflat(i) = maxloc(predratio_tmp,1)
     predratio_tmp(posflat(i)) = -1.0D0
    enddo  

    count=0
    
    do j=1,np(Hneut)+np(Hplus)
     do i=1,numres
      count=count+1     
      do k=1,numres
       if(count.eq.posflat(k)) then
        nH_rank(k) = j
        pos_rank(k) = i
       endif
      enddo 
     enddo
    enddo
     
    HBresult  = HBresult_all(nH_rank(pos),pos_rank(pos))
    chan      = chan_all(nH_rank(pos),pos_rank(pos))       
    obsratio  = obsratio_all(nH_rank(pos),pos_rank(pos))
    predratio  = predratio_all(nH_rank(pos),pos_rank(pos))    
    ncombined = ncombined_all(nH_rank(pos),pos_rank(pos))
    
   else
    write(*,*) 'WARNING: request exceeds the number of stored most sensitive channels (',numres,')'
   endif
 else
  write(*,*) 'WARNING: Please call run_HiggsBounds or run_HiggsBounds_full before calling',&
  &          ' HiggsBounds_get_most_sensitive_channels!'
 endif 

 deallocate(nH_rank,pos_rank,posflat,predratio_tmp)

end subroutine HiggsBounds_get_most_sensitive_channels
!************************************************************
subroutine run_HiggsBounds_classic( HBresult,chan,obsratio,ncombined)
! This subroutine can be called by the user after HiggsBounds_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsBounds.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental limits
! Arguments (output): 
!   * HBresult = 1 if point is unexcluded, 0 if excluded, -1 if parameter point is invalid
!   * chan = number of channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!   * obsratio = ratio of the theoretical rate to the observed limit for this channel
!   * ncombined = number of Higgs combined in order to calculate this obsratio
!    (see manual for more precise definitions))
 use usefulbits, only : theo,res,debug,just_after_run,ndmh,diffmhneut,diffmhch, &
 np,Hneut,Hplus,full_dmth_variation,dmhsteps, ndat,fullHBres!,inputsub
 use channels, only : check_channels
 !use input, only : test_input
 use theo_manip, only : HB5_complete_theo, HB5_recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none                      
 !----------------------------------------output
 integer,intent(out)::     HBresult,chan,ncombined
 double precision,intent(out) :: obsratio
 
 double precision :: Mhneut(np(Hneut))
 double precision :: Mhch(np(Hplus))
 !-------------------------------------internal
 integer :: n,i
 integer :: HBresult_tmp,chan_tmp,ncombined_tmp
 double precision :: obsratio_tmp
 !---------------------------------------------
! n=1
 
!  print *, "Running HiggsBounds in Classic Mode (globally most sensitive limit only)"
    
  
 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

!  do i=1,ubound(inputsub,dim=1)
!    if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
!        write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
!        write(*,*)'should be called once and only once before each call to'
!        write(*,*)'subroutine run_HiggsBounds.'
!        stop 'error in subroutine run_HiggsBounds'
!    endif
!    inputsub(i)%stat=0!now we have used this input, set back to zero
!  enddo

 call HB5_complete_theo    

 do n=1,ndat
        
  theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M

  call get_mass_variation_param(n)

  IF (ndmh.GT.0) THEN
   
    if(np(Hneut).ne.0) Mhneut =  theo(n)%particle(Hneut)%M     
    if(np(Hplus).ne.0) Mhch = theo(n)%particle(Hplus)%M
  
    obsratio_tmp = 10.0E6 ! Set to very large initial value
    do i=1,dmhsteps**ndmh 
     theo(n)%particle(Hneut)%M = diffMhneut(i,:)
     theo(n)%particle(Hplus)%M = diffMhch(i,:)
		 
     if(debug)write(*,*)'manipulating input...'                 ; call flush(6)
     call HB5_recalculate_theo_for_datapoint(n)
	
     if(debug)write(*,*)'compare each data point to the experimental bounds...' ; call flush(6)                  
     call check_channels(theo(n),res(n),0)
		
     HBresult    = res(n)%allowed95(1)
     chan        = res(n)%chan(1)       
     obsratio    = res(n)%obsratio(1)
     ncombined   = res(n)%ncombined(1)        
!     print *, HBresult, chan, obsratio, ncombined
		
     IF (.NOT.full_dmth_variation) THEN
      IF (HBresult.EQ.1) THEN
!		 	 theo(n)%particle(Hneut)%M = Mhneut
! 			 theo(n)%particle(Hplus)%M = Mhch
       just_after_run=.True.
       exit
      ENDIF
     ELSE
      IF (obsratio.lt.obsratio_tmp) THEN
       HBresult_tmp  = HBresult
       chan_tmp      = chan
       obsratio_tmp  = obsratio
       ncombined_tmp = ncombined
      ENDIF
     ENDIF  
   enddo
   
   IF (full_dmth_variation) THEN
    HBresult  = HBresult_tmp
    chan      = chan_tmp
    obsratio  = obsratio_tmp
    ncombined = ncombined
!    theo(n)%particle(Hneut)%M = Mhneut
!    theo(n)%particle(Hplus)%M = Mhch
    just_after_run=.True.
!	return
   ENDIF
   
   theo(n)%particle(Hneut)%M = Mhneut
   theo(n)%particle(Hplus)%M = Mhch
   call HB5_recalculate_theo_for_datapoint(n)
   call check_channels(theo(n),res(n),0)

 ELSE
 
   if(debug)write(*,*)'manipulating input...'                 ; call flush(6)
   call HB5_recalculate_theo_for_datapoint(n)

   if(debug)write(*,*)'compare each data point to the experimental bounds...' ; call flush(6)                  
   
   call check_channels(theo(n),res(n),0)
   
   HBresult    = res(n)%allowed95(1)
   chan        = res(n)%chan(1)       
   obsratio    = res(n)%obsratio(1)
   ncombined   = res(n)%ncombined(1)        

   just_after_run=.True.
 ENDIF

 fullHBres(n)%allowed95=HBresult
 fullHBres(n)%chan=chan
 fullHBres(n)%obsratio=obsratio
 fullHBres(n)%ncombined=ncombined

 enddo 

 just_after_run=.True.    

end subroutine run_HiggsBounds_classic
!************************************************************
subroutine HiggsBounds_get_likelihood(analysisID, Hindex, nc, cbin, M, llh, obspred)
!************************************************************
 use usefulbits, only : theo,np,Hneut,Hplus
 use theo_manip, only : HB5_complete_theo
 use likelihoods, only : get_likelihood, calcpredratio_llh
 implicit none
 
 integer, intent(in) :: analysisID
 integer, intent(out) :: Hindex, nc, cbin
 double precision, intent(out) :: llh, M
 character(LEN=*), intent(in) :: obspred  

 integer :: i
 double precision, allocatable :: expllh(:) 

! double precision :: fact
 double precision, allocatable :: mass(:) ! predratio(:)
 integer, allocatable :: nclist(:) 
! call complete_theo
! allocate(predratio(np(Hneut)))
! predratio = 0.0D0

! write(*,*) "Calling HiggsBounds_get_likelihood..." 

 allocate(expllh(np(Hneut)),mass(np(Hneut)),nclist(np(Hneut)))
 expllh = 0.0D0
 
!  select case(analysisID)
!   case(14029)
!    c=1
!   case(16037)
!    c=2
!   case(170907242)
!    c=3
!   case default
!    stop 'Unknown analysisID in subroutine HiggsBounds_get_likelihood!' 
!  end select  
  
 call HB5_complete_theo   
  
! Determine most sensitive combination
 do i=1,np(Hneut)
  call get_likelihood(analysisID, i, theo(1), expllh(i), mass(i), nclist(i), cbin, 'pred')
 enddo

 Hindex = maxloc(expllh,dim=1)
 
 call get_likelihood(analysisID, Hindex, theo(1), llh, M, nc, cbin, obspred)
 
 deallocate(mass,nclist,expllh) !predratio

end subroutine HiggsBounds_get_likelihood
!************************************************************
subroutine HiggsBounds_get_combined_likelihood(analysisID, llh, obspred)
!************************************************************
 use usefulbits, only : theo,np,Hneut,Hplus, vsmall

 integer, intent(in) :: analysisID
 character(LEN=*), intent(in), optional :: obspred  
 double precision, intent(out) :: llh

 double precision :: M, llh_tmp
 integer :: i, j, nc, cbin, Hindex, cbin_end, cbin_in
 
 write(*,*) 'WARNING: The subroutine HiggsBounds_get_combined_likelihood is NOT '
 write(*,*) '         officially validated and approved. Use it on your own risk!' 
 
 cbin_end = 0
 do i= 1,np(Hneut)
  cbin_end = cbin_end + 2**(i-1)
 enddo 
 
 llh = -1.0D0
 cbin_in = 0
 llh_tmp = 0.0D0
 
 do while(cbin_in.lt.cbin_end)
  if(present(obspred)) then
   call HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, obspred)
  else
   call HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, 'obs')
  endif 
  if(llh.ge.0.0D0) then
   llh_tmp = llh_tmp + llh
  else
   exit
  endif  
  cbin_in = cbin_in + cbin
 enddo  

 if(llh_tmp.gt.0.0D0) then 
  llh = llh_tmp
 endif 

 
end subroutine HiggsBounds_get_combined_likelihood
!************************************************************
subroutine HiggsBounds_get_likelihood_for_Higgs(analysisID, cbin_in, Hindex, nc, cbin, M, llh, obspred)
!************************************************************
 use usefulbits, only : theo,np,Hneut,Hplus
 use theo_manip, only : HB5_complete_theo
 use likelihoods, only : get_likelihood, calcpredratio_llh
 implicit none
 
 integer, intent(in) :: analysisID,Hindex
 integer, intent(out) ::  nc, cbin
 double precision, intent(out) :: llh, M
 integer, intent(in) :: cbin_in
 character(LEN=*), intent(in) :: obspred   
 integer :: i
 
!  select case(analysisID)
!   case(3316,14029)
!    c=1
!   case default
!    stop 'Unknown analysisID in subroutine HiggsBounds_get_likelihood_for_Higgs!' 
!  end select

 call HB5_complete_theo 
 
 call get_likelihood(analysisID, Hindex, theo(1), llh, M, nc, cbin, obspred, cbin_in)   
 
end subroutine HiggsBounds_get_likelihood_for_Higgs
!************************************************************
subroutine HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, obspred)
!************************************************************
 use usefulbits, only : theo,np,Hneut,Hplus
 use theo_manip, only : HB5_complete_theo
 use likelihoods, only : get_likelihood, calcpredratio_llh
 implicit none
 
 integer, intent(in) :: analysisID, cbin_in
 integer, intent(out) :: Hindex, nc, cbin
 double precision, intent(out) :: llh, M
 character(LEN=*), intent(in) :: obspred  

 integer :: i
 double precision, allocatable :: obsllh(:) 
 double precision, allocatable :: mass(:)
 integer, allocatable :: nclist(:), cbinlist(:)

 allocate(obsllh(np(Hneut)),mass(np(Hneut)),nclist(np(Hneut)),cbinlist(np(Hneut)))
 obsllh = 0.0D0
 
!  select case(analysisID)
!   case(3316,14029)
!    c=1
!   case default
!    stop 'Unknown analysisID in subroutine HiggsBounds_get_likelihood_for_comb!' 
!  end select

 call HB5_complete_theo 
  
! Determine most sensitive combination
 do i=1,np(Hneut)
  call get_likelihood(analysisID, i, theo(1), obsllh(i), mass(i), nclist(i),cbinlist(i),obspred, cbin_in)
 enddo

 Hindex = maxloc(obsllh,dim=1) 
 llh = obsllh(Hindex)
 M = mass(Hindex)
 nc = nclist(Hindex)
 cbin = cbinlist(Hindex)
 
 deallocate(mass,nclist,obsllh,cbinlist)

end subroutine HiggsBounds_get_likelihood_for_comb
!************************************************************ 
subroutine HiggsBounds_SLHA_output
!**** ******************************************************** 
use usefulbits, only : whichinput,just_after_run
use output, only : do_output

  if(.not.just_after_run)then
   stop 'subroutine run_HiggsBounds should be called before subroutine HiggsBounds_SLHA_output' 
  endif  

  select case(whichinput)
  case('SLHA') 
    call do_output
  case default
    stop 'The subroutine HiggsBounds_SLHA_output should only be used when whichinput=SLHA'
  end select

end subroutine HiggsBounds_SLHA_output

#ifdef enableCHISQ
!************************************************************      
subroutine initialize_HiggsBounds_chisqtables
! use S95tables, only : S95_t2
 use S95tables_type3
 use usefulbits, only : allocate_if_stats_required,theo
 implicit none            
 
 if(allocated(theo))then
  stop 'subroutine initialize_HiggsBounds_chisqtables should be called before subroutine HiggsBounds_initialize'
 elseif(allocated(clsb_t3))then
  stop 'subroutine initialize_HiggsBounds_chisqtables has already been called once'
 endif 

 allocate(clsb_t3(ntable3))

 call initializetables_type3_blank(clsb_t3)
 call initializetables3(clsb_t3)

 call readclsbfiles_binary

 if(allocated(allocate_if_stats_required))then
   stop 'error in subroutine initialize_HiggsBounds_chisqtables'
 else
   allocate(allocate_if_stats_required(1))
 endif


end subroutine initialize_HiggsBounds_chisqtables

!************************************************************
subroutine finish_HiggsBounds_chisqtables
!************************************************************
 use S95tables_type3
 use usefulbits, only : allocate_if_stats_required
 implicit none      
 integer :: x      

 if(.not.allocated(clsb_t3))then
  stop 'initialize_HiggsBounds_chisqtables should be called first'
 endif

 do x=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)
  deallocate(clsb_t3(x)%dat)
 enddo
 deallocate(filename)
 deallocate(clsb_t3)
 deallocate(allocate_if_stats_required)

end subroutine finish_HiggsBounds_chisqtables

!************************************************************
subroutine HB_calc_stats(theory_uncertainty_1s,chisq_withouttheory,chisq_withtheory,chan2)
!************************************************************
! this is in the middle of development! DO NOT USE!
 use usefulbits, only : res,theo,pr,just_after_run,vsmall
 use interpolate
 use S95tables_type1
 use S95tables_type3
 use S95tables
 use extra_bits_for_chisquared
 implicit none      

  integer,intent(out)::chan2
  integer :: x,c,z,y
  integer :: id
  double precision, intent(in) :: theory_uncertainty_1s
  double precision :: chisq_withouttheory,chisq_withtheory
  double precision :: low_chisq,sigma

  x=1
  low_chisq=1.0D-2

  if(.not.allocated(theo))then
   stop 'subroutine HiggsBounds_initialize must be called first'
  elseif(.not.allocated(clsb_t3))then
   stop 'subroutine initialize_HiggsBounds_chisqtables must be called first'
  elseif(.not.just_after_run)then
   stop 'subroutine run_HiggsBounds must be called first'
  endif

  sigma=theory_uncertainty_1s
  if(sigma.lt.vsmall)then 
   write(*,*)'Warning: will not calculate chi^2 with theory uncertainty'
  endif

  chisq_withtheory              = -2.0D0
  chisq_withouttheory           = -2.0D0

  z=2;
  c= res(x)%chan(z)
  chan2=c

  if(res(x)%allowed95(z).eq.-1)then! labels an unphysical parameter point
   chisq_withtheory             =-1.0D0
   chisq_withouttheory          =-1.0D0
  elseif( c.gt.0 )then ! labels a physical parameter point and a real channel                      
  
   id=S95_t1_or_S95_t2_idfromelementnumber(pr(c)%ttype,pr(c)%tlist)
   y=clsb_t3elementnumber_from_S95table(pr(c)%ttype,id)

   if(y.gt.0)then

    !------------------------------

    call get_chisq(sigma,res(x)%axis_i(z),res(x)%axis_j(z),res(x)%sfactor(z), &
      &  y,chisq_withouttheory,chisq_withtheory)

    !-------------------------------

   else
     write(*,*)'hello y=',y
     stop 'problem here with y'
   endif

  else
   chisq_withtheory             =0.0D0
   chisq_withouttheory          =0.0D0
  endif  

end subroutine HB_calc_stats
#endif
!************************************************************
subroutine finish_HiggsBounds
! This subroutine needs to be called right at the end, to close files
! and deallocate arrays
!************************************************************
 use usefulbits, only : deallocate_usefulbits,debug,theo,debug, &
   & file_id_debug1,file_id_debug2!,inputsub
 use S95tables, only : deallocate_S95tables
 use theory_BRfunctions, only : deallocate_BRSM
 use theory_XS_SM_functions, only: deallocate_XSSM

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(.not.allocated(theo))then
  stop 'HiggsBounds_initialize  should be called first'
 endif

 if(debug)write(*,*)'finishing off...'                      ; call flush(6)
 call deallocate_BRSM
 call deallocate_XSSM
 call deallocate_S95tables
 call deallocate_usefulbits

 if(debug)write(*,*)'finished'                              ; call flush(6)
 
!  if(allocated(inputsub)) deallocate(inputsub)
end subroutine finish_HiggsBounds
!

! HB-5 additions

! Do we need control functions to guarantee all theory input is up-to-date and reset?

!subroutine HB5_reset_input
!end subroutine HB5_reset_input

  

!************************************************************      
!
!   SIMPLIFIED EFFC INPUT ROUTINES
!
!************************************************************      
subroutine HiggsBounds_neutral_input_effC_single(quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,effC,whichinput,just_after_run!,inputsub
 implicit none
 !---------------------------------------------
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut))
 !-------------------------------------internal
 integer :: n
!  integer :: subtype
 !--------------------------------------------- 
 whichinput='effC'
!  subtype=1
 n=1 
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1 ! WHAT IS THIS DOING?

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_effC_single should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_effC_single'
 endif

 select case(trim(adjustl(quantity)))
  case("ghjcc_s")
   effC(n)%hjcc_s=val
  case("ghjcc_p")
   effC(n)%hjcc_p=val
  case("ghjss_s")
   effC(n)%hjss_s=val
  case("ghjss_p")
   effC(n)%hjss_p=val
  case("ghjbb_s")
   effC(n)%hjbb_s=val
  case("ghjbb_p")
   effC(n)%hjbb_p=val
  case("ghjtt_s")
   effC(n)%hjtt_s=val
  case("ghjtt_p")
   effC(n)%hjtt_p=val
  case("ghjmumu_s")
   effC(n)%hjmumu_s=val
  case("ghjmumu_p")
   effC(n)%hjmumu_p=val
  case("ghjtautau_s")
   effC(n)%hjtautau_s=val
  case("ghjtautau_p")
   effC(n)%hjtautau_p=val
  case("ghjWW")
   effC(n)%hjWW=val
  case("ghjZZ")
   effC(n)%hjZZ=val
  case("ghjZga")
   effC(n)%hjZga=val
  case("ghjgaga")
   effC(n)%hjgaga=val
  case("ghjgg")
   effC(n)%hjgg=val


  case default
   stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_effC_single' 
 end select
  
 just_after_run=.False.

end subroutine HiggsBounds_neutral_input_effC_single
!************************************************************      
subroutine HiggsBounds_neutral_input_effC_double(quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,effC,whichinput,just_after_run!,inputsub
 implicit none
 !---------------------------------------------
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut),np(Hneut))
 !-------------------------------------internal
 integer :: n
!  integer :: subtype
 !--------------------------------------------- 
 whichinput='effC'
!  subtype=1
 n=1 
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1 ! WHAT IS THIS DOING?

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_effC_double should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_effC_double'
 endif

 select case(trim(adjustl(quantity)))
  case("ghjhiZ") 
   effC(n)%hjhiZ = val 
  case default
   stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_effC_double' 
 end select

 just_after_run=.False.
end subroutine HiggsBounds_neutral_input_effC_double
!************************************************************      
!
!   SIMPLIFIED LEP/HADRONIC XS INPUT ROUTINES
!
!************************************************************      
subroutine HiggsBounds_neutral_input_LEP_single(quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run!,inputsub
 implicit none
 !---------------------------------------------
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut))
 !-------------------------------------internal
 integer :: n
!  integer :: subtype
 !--------------------------------------------- 
 whichinput='hadr'
!  subtype=1
 n=1 
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1 ! WHAT IS THIS DOING?

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_LEP_single should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_LEP_single'
 endif

 select case(trim(adjustl(quantity)))
  case("XS_hjZ_ratio")
   theo(n)%lep%XS_hjZ_ratio = val
  case("XS_bbhj_ratio")
   theo(n)%lep%XS_bbhj_ratio = val
  case("XS_tautauhj_ratio")
   theo(n)%lep%XS_tautauhj_ratio = val
  case default
   stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_LEP_single' 
 end select
  
 just_after_run=.False.
end subroutine HiggsBounds_neutral_input_LEP_single
!************************************************************      
subroutine HiggsBounds_neutral_input_LEP_double(quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run!,inputsub
 implicit none
 !---------------------------------------------
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut),np(Hneut))
 !-------------------------------------internal
 integer :: n
!  integer :: subtype
 !--------------------------------------------- 
 whichinput='hadr'
!  subtype=1
 n=1 
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1 ! WHAT IS THIS DOING?

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_LEP_double should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_LEP_double'
 endif

 select case(trim(adjustl(quantity)))
  case("XS_hjhi_ratio") 
   theo(n)%lep%XS_hjhi_ratio = val 
  case default
   stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_LEP_double' 
 end select

 just_after_run=.False.
end subroutine HiggsBounds_neutral_input_LEP_double
!************************************************************      
subroutine HiggsBounds_neutral_input_hadr_single(collider,quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run, &
&                       hadroncolliderdataset !,inputsub
 implicit none
 !---------------------------------------------
 integer, intent(in) :: collider
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut))
 !-------------------------------------internal
 integer :: n
!  integer :: subtype
!  type(hadroncolliderdataset) :: dataset
 !--------------------------------------------- 
 whichinput='hadr'
!  subtype=1
 n=1 
!  inputsub(subtype)%stat=inputsub(subtype)%stat+1 ! WHAT IS THIS DOING?

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr_single should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr_single'
 endif

 select case(collider)
  case(2)
   call set_input(theo(n)%tev,quantity,val)
  case(7)
   call set_input(theo(n)%lhc7,quantity,val)
  case(8)
   call set_input(theo(n)%lhc8,quantity,val)
  case(13) 
   call set_input(theo(n)%lhc13,quantity,val)
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_single'
 end select

  
 just_after_run=.False.

 contains
 
  subroutine set_input(dataset,quantity,val)

   character(LEN=*), intent(in) :: quantity 
   double precision, intent(in) :: val(np(Hneut))
   type(hadroncolliderdataset) :: dataset

   select case(trim(adjustl(quantity)))
    case("XS_hj_ratio")
     dataset%XS_hj_ratio=val
    case("XS_gg_hj_ratio")
	 dataset%XS_gg_hj_ratio=val
	case("XS_bb_hj_ratio")
 	 dataset%XS_bb_hj_ratio=val 
	 dataset%XS_hjb_ratio=val
	case("XS_vbf_ratio")
	 dataset%XS_vbf_ratio=val
	case("XS_hjZ_ratio")
	 dataset%XS_hjZ_ratio=val
	case("XS_gg_hjZ_ratio")
	 dataset%XS_gg_hjZ_ratio=val	 
	case("XS_qq_hjZ_ratio")
	 dataset%XS_qq_hjZ_ratio=val
	case("XS_hjW_ratio")
	 dataset%XS_hjW_ratio=val
	case("XS_tthj_ratio")
	 dataset%XS_tthj_ratio=val
	case("XS_thj_tchan_ratio")
	 dataset%XS_thj_tchan_ratio=val
	case("XS_thj_schan_ratio")
	 dataset%XS_thj_schan_ratio=val
	case default
	 stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_hadr_single' 
   end select
  end subroutine set_input

end subroutine HiggsBounds_neutral_input_hadr_single
!************************************************************      
subroutine HiggsBounds_neutral_input_hadr_double(collider,quantity,val)
!************************************************************      
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run, &
&                       hadroncolliderdataset! ,inputsub
 implicit none
 !---------------------------------------------
 integer, intent(in) :: collider
 character(LEN=*), intent(in) :: quantity 
 double precision, intent(in) :: val(np(Hneut),np(Hneut))
 !-------------------------------------internal
 integer :: n
 !--------------------------------------------- 
 whichinput='hadr'
 n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr_double should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr_double'
 endif

 select case(trim(adjustl(quantity)))
  case("XS_hjhi")
   select case(collider)
    case(2)
     theo(n)%tev%XS_hjhi=val
    case(7)
     theo(n)%lhc7%XS_hjhi=val
    case(8) 
     theo(n)%lhc8%XS_hjhi=val
    case(13) 
     theo(n)%lhc13%XS_hjhi=val
    case default
     stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_double'
   end select
  case default
   stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_hadr_double' 
 end select
  
 just_after_run=.False.

end subroutine HiggsBounds_neutral_input_hadr_double
!************************************************************
subroutine HiggsBounds_neutral_input_hadr_channelrates_single(collider,nHiggs,p,d,val)
! n.b.: Elements of the matrix channelrates with values < 0 will be overwritten
!       by XS times BR using the narrow width approximation.
!************************************************************
 use usefulbits, only : theo,np,Hneut,whichinput,just_after_run,hadroncolliderdataset,&
 &                      Nprod,Ndecay

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 double precision,intent(in) :: val
 integer, intent(in) :: collider,p,d,nHiggs
 !-------------------------------------internal
 integer :: n
 !---------------------------------------------
 whichinput='hadr'
  n=1 

 if(.not.allocated(theo))then
  stop 'subroutine HiggsBounds_initialize must be called first'
 endif

 if(nHiggs.gt.np(Hneut))then
  write(*,*)'subroutine HiggsBounds_neutral_input_hadr_channelrates_single should'
  write(*,*)'only be called with nHiggs <= np(Hneut)'
  stop 'error in subroutine HiggsBounds_neutral_input_hadr_channelrates_single'
 endif

 select case(collider)
  case(2)
   theo(n)%tev%channelrates_tmp(nHiggs,p,d)=val
  case(7)
   theo(n)%lhc7%channelrates_tmp(nHiggs,p,d)=val
  case(8) 
   theo(n)%lhc8%channelrates_tmp(nHiggs,p,d)=val 
  case(13) 
   theo(n)%lhc13%channelrates_tmp(nHiggs,p,d)=val
  case default
   stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_channelrates_single'
 end select

 just_after_run=.False.       
      
end subroutine HiggsBounds_neutral_input_hadr_channelrates_single
!************************************************************
subroutine HiggsBounds_neutral_input_hadr_channelrates_clean
!************************************************************
 use theo_manip,only : clean_channelrates
 implicit none
  call clean_channelrates
  
end subroutine HiggsBounds_neutral_input_hadr_channelrates_clean
!************************************************************
! HB-4 legacy routines
!************************************************************     
! subroutine HiggsBounds_neutral_input_effC(Mh,GammaTotal_hj,        &  
!      &          g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,               &
!      &          g2hjbb_s,g2hjbb_p,g2hjtoptop_s,g2hjtoptop_p,       &
!      &          g2hjmumu_s,g2hjmumu_p,                             &
!      &          g2hjtautau_s,g2hjtautau_p,                         &
!      &          g2hjWW,g2hjZZ,g2hjZga,                             &
!      &          g2hjgaga,g2hjgg,g2hjggZ,g2hjhiZ_nHbynH,            &
!      &          BR_hjinvisible,BR_hjhihi_nHbynH                    )
! ! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! ! has been called.
! ! Arguments (input): theoretical predictions (see manual for definitions)
! !************************************************************
!  use usefulbits, only : theo,np,Hneut,g2,whichinput,just_after_run!,inputsub
! 
! #if defined(NAGf90Fortran)
!  use F90_UNIX_IO, only : flush
! #endif
! 
!  implicit none
! 
!  !----------------------------------------input
!  double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) ),                                 &  
!      &          g2hjss_s( np(Hneut) ),g2hjss_p( np(Hneut) ),g2hjcc_s( np(Hneut) ),g2hjcc_p( np(Hneut) ),        &
!      &          g2hjbb_s( np(Hneut) ),g2hjbb_p( np(Hneut) ),g2hjtoptop_s( np(Hneut) ),g2hjtoptop_p( np(Hneut) ),&
!      &          g2hjmumu_s( np(Hneut) ),g2hjmumu_p( np(Hneut) ),                                        &
!      &          g2hjtautau_s( np(Hneut) ),g2hjtautau_p( np(Hneut) ),                                        &
!      &          g2hjWW( np(Hneut) ),g2hjZZ( np(Hneut) ),g2hjZga( np(Hneut) ),                                 &
!      &          g2hjgaga( np(Hneut) ),g2hjgg( np(Hneut) ),g2hjggZ( np(Hneut) ),g2hjhiZ_nHbynH(np(Hneut),np(Hneut)),&
!      &          BR_hjinvisible( np(Hneut) ),BR_hjhihi_nHbynH(np(Hneut),np(Hneut)) 
!  !--------------------------------------internal
!  integer :: n
! !  integer :: subtype
!  !----------------------------------------------
!  
!  whichinput='effC'
!  ! subtype=1
!   n=1  
! !  inputsub(subtype)%stat=inputsub(subtype)%stat+1
! 
!  if(.not.allocated(theo))then
!   stop 'subroutine HiggsBounds_initialize must be called first'
!  endif
! 
!  if(np(Hneut).eq.0)then
!   write(*,*)'subroutine HiggsBounds_neutral_input_effC should'
!   write(*,*)'only be called if np(Hneut)>0'
!   stop 'error in subroutine HiggsBounds_neutral_input_effC'
!  endif
! 
!  theo(n)%particle(Hneut)%M       = Mh
!  theo(n)%particle(Hneut)%Mc      = Mh    
!  theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
! 
!  g2(n)%hjss_s              = g2hjss_s
!  g2(n)%hjss_p              = g2hjss_p
!  g2(n)%hjcc_s              = g2hjcc_s
!  g2(n)%hjcc_p              = g2hjcc_p
!  g2(n)%hjbb_s              = g2hjbb_s
!  g2(n)%hjbb_p              = g2hjbb_p
!  g2(n)%hjtoptop_s              = g2hjtoptop_s
!  g2(n)%hjtoptop_p              = g2hjtoptop_p
!  g2(n)%hjmumu_s                = g2hjmumu_s
!  g2(n)%hjmumu_p                = g2hjmumu_p
!  g2(n)%hjtautau_s              = g2hjtautau_s
!  g2(n)%hjtautau_p              = g2hjtautau_p
!         
!  g2(n)%hjWW              = g2hjWW
!  g2(n)%hjZZ              = g2hjZZ  
!  g2(n)%hjZga             = g2hjZga      
!  g2(n)%hjgaga            = g2hjgaga
!  g2(n)%hjgg              = g2hjgg
!  g2(n)%hjggZ             = g2hjggZ
! 
!  g2(n)%hjhiZ = g2hjhiZ_nHbynH
! 
!  theo(n)%BR_hjinvisible  = BR_hjinvisible 
!  theo(n)%BR_hjhihi       = BR_hjhihi_nHbynH  
! 
!  just_after_run=.False. 
! 
! end subroutine HiggsBounds_neutral_input_effC
! !************************************************************      
! subroutine HiggsBounds_neutral_input_part(Mh,GammaTotal_hj,CP_value,       &
!      &          CS_lep_hjZ_ratio,                            &
!      &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     &
!      &          CS_lep_hjhi_ratio_nHbynH,                    &
!      &          CS_gg_hj_ratio,CS_bb_hj_ratio,       &
!      &          CS_bg_hjb_ratio,                         &
!      &          CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,   & 
!      &          CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,   & 
!      &          CS_gg_hjZ_ratio,     &
!      &          CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     &
!      &          CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     &
!      &          CS_bb_hjZ_ratio,                         &
!      &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,    &
!      &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    &
!      &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,    &
!      &          BR_hjss,BR_hjcc,                             &
!      &          BR_hjbb,BR_hjmumu,BR_hjtautau,               &
!      &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg, & 
!      &          BR_hjinvisible,BR_hjhihi_nHbynH              )
! ! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! ! has been called.
! ! (see manual for full description)
! !************************************************************
!  use usefulbits, only : theo,np,Hneut,partR,whichinput,just_after_run!,inputsub
! 
! #if defined(NAGf90Fortran)
!  use F90_UNIX_IO, only : flush
! #endif
! 
!  implicit none
!  !----------------------------------------input
!  double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) )
!  integer,intent(in) ::CP_value( np(Hneut) )
!  double precision,intent(in) :: CS_lep_hjZ_ratio( np(Hneut) ),                         &
!      &          CS_lep_bbhj_ratio( np(Hneut) ),CS_lep_tautauhj_ratio( np(Hneut) ),     &
!      &          CS_lep_hjhi_ratio_nHbynH(np(Hneut),np(Hneut)),                         &
!      &          CS_gg_hj_ratio( np(Hneut) ),CS_bb_hj_ratio( np(Hneut) ),       &
!      &          CS_bg_hjb_ratio( np(Hneut) ),                                      &
!      &          CS_ud_hjWp_ratio( np(Hneut) ),CS_cs_hjWp_ratio( np(Hneut) ),   & 
!      &          CS_ud_hjWm_ratio( np(Hneut) ),CS_cs_hjWm_ratio( np(Hneut) ),   & 
!      &          CS_gg_hjZ_ratio( np(Hneut) ),    &
!      &          CS_dd_hjZ_ratio( np(Hneut) ),CS_uu_hjZ_ratio( np(Hneut) ),     &
!      &          CS_ss_hjZ_ratio( np(Hneut) ),CS_cc_hjZ_ratio( np(Hneut) ),     &
!      &          CS_bb_hjZ_ratio( np(Hneut) ),                                      &
!      &          CS_tev_vbf_ratio( np(Hneut) ),CS_tev_tthj_ratio( np(Hneut) ),    &
!      &          CS_lhc7_vbf_ratio( np(Hneut) ),CS_lhc7_tthj_ratio( np(Hneut) ),    &
!      &          CS_lhc8_vbf_ratio( np(Hneut) ),CS_lhc8_tthj_ratio( np(Hneut) ),    &
!      &          BR_hjss( np(Hneut) ),BR_hjcc( np(Hneut) ),                             &
!      &          BR_hjbb( np(Hneut) ),BR_hjmumu( np(Hneut) ),BR_hjtautau( np(Hneut) ),  &
!      &          BR_hjWW( np(Hneut) ),BR_hjZZ( np(Hneut) ),BR_hjZga( np(Hneut) ),       &
!      &          BR_hjgaga( np(Hneut) ),BR_hjgg( np(Hneut) ),                           & 
!      &          BR_hjinvisible( np(Hneut) ),BR_hjhihi_nHbynH(np(Hneut),np(Hneut)) 
!  !---------------------------------------internal
!  integer :: n
! !  integer :: subtype
!  !-----------------------------------------------
! 
!  whichinput='part'
! !  subtype=1
!   n=1 
! !  inputsub(subtype)%stat=inputsub(subtype)%stat+1
!       
!  if(.not.allocated(theo))then
!   stop 'subroutine HiggsBounds_initialize must be called first'
!  endif
! 
!  if(np(Hneut).eq.0)then
!   write(*,*)'subroutine HiggsBounds_neutral_input_part should'
!   write(*,*)'only be called if np(Hneut)>0'
!   stop 'error in subroutine HiggsBounds_neutral_input_part'
!  endif
! 
!  theo(n)%particle(Hneut)%M = Mh 
!  theo(n)%particle(Hneut)%Mc      = Mh        
!  theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
!  theo(n)%CP_value          = CP_value  
!  theo(n)%lep%XS_hjZ_ratio       = CS_lep_hjZ_ratio
!  theo(n)%lep%XS_bbhj_ratio      = CS_lep_bbhj_ratio
!  theo(n)%lep%XS_tautauhj_ratio  = CS_lep_tautauhj_ratio
!  theo(n)%lep%XS_hjhi_ratio      = CS_lep_hjhi_ratio_nHbynH 
!  partR(n)%gg_hj             = CS_gg_hj_ratio 
!  partR(n)%qq_hj(5,:)        = CS_bb_hj_ratio
!  partR(n)%bg_hjb            = CS_bg_hjb_ratio                    
!  partR(n)%qq_hjWp(1,:)      = CS_ud_hjWp_ratio
!  partR(n)%qq_hjWp(2,:)      = CS_cs_hjWp_ratio         
!  partR(n)%qq_hjWm(1,:)      = CS_ud_hjWm_ratio
!  partR(n)%qq_hjWm(2,:)      = CS_cs_hjWm_ratio  
!  partR(n)%gg_hjZ(:)         = CS_gg_hjZ_ratio      
!  partR(n)%qq_hjZ(1,:)       = CS_dd_hjZ_ratio
!  partR(n)%qq_hjZ(2,:)       = CS_uu_hjZ_ratio           
!  partR(n)%qq_hjZ(3,:)       = CS_ss_hjZ_ratio
!  partR(n)%qq_hjZ(4,:)       = CS_cc_hjZ_ratio             
!  partR(n)%qq_hjZ(5,:)       = CS_bb_hjZ_ratio                         
!  theo(n)%tev%XS_vbf_ratio  = CS_tev_vbf_ratio   
!  theo(n)%tev%XS_tthj_ratio = CS_tev_tthj_ratio  
!  theo(n)%lhc7%XS_vbf_ratio = CS_lhc7_vbf_ratio   
!  theo(n)%lhc7%XS_tthj_ratio= CS_lhc7_tthj_ratio
!  theo(n)%lhc8%XS_vbf_ratio = CS_lhc8_vbf_ratio   
!  theo(n)%lhc8%XS_tthj_ratio= CS_lhc8_tthj_ratio
!  theo(n)%BR_hjss           = BR_hjss  
!  theo(n)%BR_hjcc           = BR_hjcc                
!  theo(n)%BR_hjbb           = BR_hjbb
!  theo(n)%BR_hjmumu         = BR_hjmumu
!  theo(n)%BR_hjtautau       = BR_hjtautau              
!  theo(n)%BR_hjWW           = BR_hjWW
!  theo(n)%BR_hjZZ           = BR_hjZZ
!  theo(n)%BR_hjZga          = BR_hjZga  
!  theo(n)%BR_hjgaga         = BR_hjgaga
!  theo(n)%BR_hjgg           = BR_hjgg  
!  theo(n)%BR_hjinvisible    = BR_hjinvisible             
!  theo(n)%BR_hjhihi         = BR_hjhihi_nHbynH  
! 
!  just_after_run=.False. 
! 
! end subroutine HiggsBounds_neutral_input_part
! !************************************************************      
! subroutine HiggsBounds_neutral_input_hadr(Mh,GammaTotal_hj,CP_value,      &
!      &          CS_lep_hjZ_ratio,                           &
!      &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,    &   
!      &          CS_lep_hjhi_ratio_nHbynH,                   &
!      &          CS_tev_hj_ratio ,CS_tev_hjb_ratio,    &
!      &          CS_tev_hjW_ratio,CS_tev_hjZ_ratio,    &
!      &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,   &
!      &          CS_lhc7_hj_ratio ,CS_lhc7_hjb_ratio,    &
!      &          CS_lhc7_hjW_ratio,CS_lhc7_hjZ_ratio,    &
!      &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,   &
!      &          CS_lhc8_hj_ratio ,CS_lhc8_hjb_ratio,    &
!      &          CS_lhc8_hjW_ratio,CS_lhc8_hjZ_ratio,    &
!      &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,   &
!      &          BR_hjss,BR_hjcc,                            &
!      &          BR_hjbb,                                    &
!      &          BR_hjmumu,                                  &
!      &          BR_hjtautau,                                &
!      &          BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,         &
!      &          BR_hjgg, BR_hjinvisible,                    &
!      &          BR_hjhihi_nHbynH                            )
! ! This subroutine can be called by the user after subroutine initialize_HiggsBounds
! ! has been called.
! ! (see manual for full description)
! !************************************************************
!  use usefulbits, only : theo,np,Hneut,whichinput,just_after_run!,inputsub
! 
! #if defined(NAGf90Fortran)
!  use F90_UNIX_IO, only : flush
! #endif
! 
!  implicit none
!  !----------------------------------------input
!  double precision,intent(in) :: Mh( np(Hneut) ),GammaTotal_hj( np(Hneut) )
!  integer,intent(in) :: CP_value( np(Hneut) )
!  double precision,intent(in) :: CS_lep_hjZ_ratio( np(Hneut) ),                        &
!      &          CS_lep_bbhj_ratio( np(Hneut) ),CS_lep_tautauhj_ratio( np(Hneut) ),    &   
!      &          CS_lep_hjhi_ratio_nHbynH(np(Hneut),np(Hneut)),                        &
!      &          CS_tev_hj_ratio( np(Hneut)  ) ,CS_tev_hjb_ratio( np(Hneut) ),    &
!      &          CS_tev_hjW_ratio( np(Hneut) ) ,CS_tev_hjZ_ratio( np(Hneut) ),    &
!      &          CS_tev_vbf_ratio( np(Hneut) ) ,CS_tev_tthj_ratio( np(Hneut)),    &
!      &          CS_lhc7_hj_ratio( np(Hneut)  ),CS_lhc7_hjb_ratio( np(Hneut) ),    &
!      &          CS_lhc7_hjW_ratio( np(Hneut) ),CS_lhc7_hjZ_ratio( np(Hneut) ),    &
!      &          CS_lhc7_vbf_ratio( np(Hneut) ),CS_lhc7_tthj_ratio( np(Hneut)),    &
!      &          CS_lhc8_hj_ratio( np(Hneut)  ),CS_lhc8_hjb_ratio( np(Hneut) ),    &
!      &          CS_lhc8_hjW_ratio( np(Hneut) ),CS_lhc8_hjZ_ratio( np(Hneut) ),    &
!      &          CS_lhc8_vbf_ratio( np(Hneut) ),CS_lhc8_tthj_ratio( np(Hneut)),    &
!      &          BR_hjss( np(Hneut) ),BR_hjcc( np(Hneut) ),                            &
!      &          BR_hjbb( np(Hneut) ),                                                 &
!      &          BR_hjmumu( np(Hneut) ),BR_hjtautau( np(Hneut) ),                      &
!      &          BR_hjWW( np(Hneut) ),BR_hjZZ( np(Hneut) ),                            &
!      &          BR_hjZga( np(Hneut) ),BR_hjgaga( np(Hneut) ),                         &
!      &          BR_hjgg( np(Hneut) ), BR_hjinvisible( np(Hneut) ),                    &
!      &          BR_hjhihi_nHbynH(np(Hneut),np(Hneut))
!  !-------------------------------------internal
!  integer :: n
! !  integer :: subtype
!  !---------------------------------------------
!    
!  whichinput='hadr'
! !  subtype=1
!   n=1 
! !  inputsub(subtype)%stat=inputsub(subtype)%stat+1
! 
!  if(.not.allocated(theo))then
!   stop 'subroutine HiggsBounds_initialize must be called first'
!  endif
! 
!  if(np(Hneut).eq.0)then
!   write(*,*)'subroutine HiggsBounds_neutral_input_hadr should'
!   write(*,*)'only be called if np(Hneut)>0'
!   stop 'error in subroutine HiggsBounds_neutral_input_hadr'
!  endif
! 
! ! write(*,*) "DEBUG HB: before hadronic input. Mass is ",theo(n)%particle(Hneut)%M 
! 
! 
!  theo(n)%particle(Hneut)%M       = Mh
!  theo(n)%particle(Hneut)%Mc      = Mh  
!  theo(n)%particle(Hneut)%GammaTot= GammaTotal_hj
!  theo(n)%CP_value                = CP_value             
!  theo(n)%lep%XS_hjZ_ratio        = CS_lep_hjZ_ratio
!  theo(n)%lep%XS_bbhj_ratio       = CS_lep_bbhj_ratio
!  theo(n)%lep%XS_tautauhj_ratio   = CS_lep_tautauhj_ratio
!  theo(n)%lep%XS_hjhi_ratio       = CS_lep_hjhi_ratio_nHbynH 
!  theo(n)%tev%XS_hj_ratio         = CS_tev_hj_ratio
!  theo(n)%tev%XS_hjb_ratio        = CS_tev_hjb_ratio   
!  theo(n)%tev%XS_hjW_ratio        = CS_tev_hjW_ratio
!  theo(n)%tev%XS_hjZ_ratio        = CS_tev_hjZ_ratio    
!  theo(n)%tev%XS_vbf_ratio        = CS_tev_vbf_ratio      
!  theo(n)%tev%XS_tthj_ratio       = CS_tev_tthj_ratio
!  theo(n)%lhc7%XS_hj_ratio        = CS_lhc7_hj_ratio
!  theo(n)%lhc7%XS_hjb_ratio       = CS_lhc7_hjb_ratio   
!  theo(n)%lhc7%XS_hjW_ratio       = CS_lhc7_hjW_ratio
!  theo(n)%lhc7%XS_hjZ_ratio       = CS_lhc7_hjZ_ratio    
!  theo(n)%lhc7%XS_vbf_ratio       = CS_lhc7_vbf_ratio      
!  theo(n)%lhc7%XS_tthj_ratio      = CS_lhc7_tthj_ratio
!  theo(n)%lhc8%XS_hj_ratio        = CS_lhc8_hj_ratio
!  theo(n)%lhc8%XS_hjb_ratio       = CS_lhc8_hjb_ratio   
!  theo(n)%lhc8%XS_hjW_ratio       = CS_lhc8_hjW_ratio
!  theo(n)%lhc8%XS_hjZ_ratio       = CS_lhc8_hjZ_ratio    
!  theo(n)%lhc8%XS_vbf_ratio       = CS_lhc8_vbf_ratio      
!  theo(n)%lhc8%XS_tthj_ratio      = CS_lhc8_tthj_ratio
!  theo(n)%BR_hjss           = BR_hjss
!  theo(n)%BR_hjcc           = BR_hjcc 
!  theo(n)%BR_hjbb           = BR_hjbb
!  theo(n)%BR_hjmumu         = BR_hjmumu
!  theo(n)%BR_hjtautau       = BR_hjtautau                 
!  theo(n)%BR_hjWW           = BR_hjWW
!  theo(n)%BR_hjZZ           = BR_hjZZ
!  theo(n)%BR_hjZga          = BR_hjZga 
!  theo(n)%BR_hjgaga         = BR_hjgaga
!  theo(n)%BR_hjgg           = BR_hjgg
!  theo(n)%BR_hjinvisible    = BR_hjinvisible                  
!  theo(n)%BR_hjhihi         = BR_hjhihi_nHbynH  
! 
!  just_after_run=.False. 
!   
! ! write(*,*) "DEBUG HB: filled hadronic input. Mass is ",theo(n)%particle(Hneut)%M 
!     
! end subroutine HiggsBounds_neutral_input_hadr
!************************************************************     