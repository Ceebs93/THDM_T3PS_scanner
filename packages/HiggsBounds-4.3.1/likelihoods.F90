! This file is part of HiggsBounds
!  -TS (23/09/2014)
!******************************************************************
module likelihoods
!******************************************************************
! Some options:
 logical :: rescale = .False.

!--
 type likelihood2D
  integer :: analysisID
  integer :: particle_x 
  character(LEN=45) :: label
  character(LEN=3) :: expt 
  double precision :: lumi, energy
  character(LEN=100) :: description
  double precision :: mass_sep
  double precision :: mass_range_tolerance
  double precision :: mass
  double precision :: xstep
  double precision :: ystep
  double precision :: xmin
  double precision :: xmax
  double precision :: ymin
  double precision :: ymax
  integer :: size
  integer :: xsize
  integer :: ysize  
  double precision, dimension(200,200) :: llhgrid_pred, llhgrid_obs
  double precision, dimension(1:3) :: bestfit_pred, bestfit_obs
 end type
 
 type threetuple
  double precision, dimension(1:3) :: val
 end type 

 integer, parameter :: nllhs = 1

! If more likelihoods become available (nllhs > 1), have to generalize this
! to a list of likelihoods

! CMS related stuff 
 integer, parameter :: n2Dslices = 38
 type(likelihood2D), dimension(1:n2Dslices) :: CMS_llhdata ! should be renamed to something more general at some point!
 
 contains
!--------------------------------------------------------------
 subroutine setup_likelihoods
!--------------------------------------------------------------
  use usefulbits, only : file_id_common4, Hneut
  use store_pathname, only : pathname  
!  implicit none
  
  integer :: ios,i,j,status
  character(LEN=4) :: charstatus
  
! for CMS-3316 likelihoods
!   integer :: analysisID = 3316   
!   character(LEN=100) :: label = '[hep-ex] arXiv:1408.3316 (CMS)'
!   character(LEN=100) :: stem = 'Expt_tables/CMStables/CMS_tautau_llh_1408.3316/'
!  character(LEN=100) :: stem = 'Expt_tables/CMStables/CMS_tautau_llh_1408.3316_withSMHinBG/'
  integer :: analysisID = 14029
  character(LEN=100) :: label = 'CMS-HIG-PAS 14-029'
  character(LEN=100) :: stem = 'Expt_tables/CMStables/CMS_tautau_llh_14029/'  
  character(LEN=100) :: description = '(pp) -> h -> tautau, using -2ln(L) reconstruction'
  character(LEN=3) :: expt = 'CMS'
  integer :: particle_type = Hneut
  double precision :: lumi = 19.7D0
  double precision :: energy = 8.0D0
  double precision :: mass_sep = 0.20D0  ! in percentage
  double precision :: mass_range_tolerance = 10.0D0 ! If the predicted mass lies outside the mass range within the tolerance,
                                                    ! the mass will be rounded to the closest value within the grid.
  double precision, dimension(1:n2Dslices) :: masses, xmin, xmax, ymin, ymax, xstep, ystep
  integer, dimension(1:n2Dslices) :: xsize, ysize
  

  masses = (/ 90.0D0, 100.0D0, 110.0D0, 120.0D0, 125.0D0, 130.0D0, 140.0D0, &
&            150.0D0, 160.0D0, 170.0D0, 180.0D0, 190.0D0, 200.0D0, 210.0D0,&
&            220.0D0, 230.0D0, 240.0D0, 250.0D0, 275.0D0, 300.0D0, 325.0D0,&
&            350.0D0, 375.0D0, 400.0D0, 425.0D0, 450.0D0, 475.0D0, 500.0D0,&
&            550.0D0, 600.0D0, 650.0D0, 700.0D0, 750.0D0, 800.0D0, 850.0D0,&
&            900.0D0, 950.0D0, 1000.0D0 /)
!-------- Settings for the 1408.3316 analysis ---------!
!   xmin = (/ 0.1750000, 0.1000000, 0.0300000, 0.0250000, 0.0250000, 0.0150000, 0.0075000,&
! &         0.0050000, 0.0037500, 0.0025000, 0.0025000, 0.0025000, 0.0020000, 0.0020000,&
! &         0.0020000, 0.0020000, 0.0015000, 0.0012500, 0.0010000, 0.0005000, 0.0004500,&
! &         0.0003750, 0.0003750, 0.0003750, 0.0003250, 0.0002500, 0.0002000, 0.0001250,&
! &         0.0001000, 0.0000750, 0.0000700, 0.0000625, 0.0000575, 0.0000500, 0.0000500,&
! &         0.0000375, 0.0000375, 0.0000375/)
!   xmax = (/69.8250000,39.9000000,11.9700000, 9.9750000, 9.9750000, 5.9850000, 2.9925000,&
! &         1.9950000, 1.4962500, 0.9975000, 0.9975000, 0.9975000, 0.7980000, 0.7980000,&
! &         0.7980000, 0.7980000, 0.5985000, 0.4987500, 0.3990000, 0.1995000, 0.1795500,&
! &         0.1496250, 0.1496250, 0.1496250, 0.1296750, 0.0997500, 0.0798000, 0.0498750,&
! &         0.0399000, 0.0299250, 0.0279300, 0.0249375, 0.0229425, 0.0199500, 0.0199500,&
! &         0.0149625, 0.0149625, 0.0149625/)
!   ymin = (/ 0.0325000, 0.0200000, 0.0175000, 0.0087500, 0.0087500, 0.0062500, 0.0062500,&
! &         0.0050000, 0.0037500, 0.0030000, 0.0025000, 0.0025000, 0.0022500, 0.0020000,&
! &         0.0017500, 0.0015000, 0.0015000, 0.0010000, 0.0008750, 0.0005000, 0.0004500,&
! &         0.0003000, 0.0003000, 0.0003000, 0.0002500, 0.0002000, 0.0001750, 0.0001500,&
! &         0.0001250, 0.0000750, 0.0000700, 0.0000625, 0.0000625, 0.0000625, 0.0000625,&
! &         0.0000625, 0.0000625, 0.0000625/)
!   ymax = (/12.9675000, 7.9800000, 6.9825000, 3.4912500, 3.4912500, 2.4937500, 2.4937500,&
! &         1.9950000, 1.4962500, 1.1970000, 0.9975000, 0.9975000, 0.8977500, 0.7980000,&
! &         0.6982500, 0.5985000, 0.5985000, 0.3990000, 0.3491250, 0.1995000, 0.1795500,&
! &         0.1197000, 0.1197000, 0.1197000, 0.0997500, 0.0798000, 0.0698250, 0.0598500,&
! &         0.0498750, 0.0299250, 0.0279300, 0.0249375, 0.0249375, 0.0249375, 0.0249375,&
! &         0.0249375, 0.0249375, 0.0249375/)
!   xstep = (/0.350, 0.200, 0.060, 0.050, 0.050, 0.030, 0.015,&
! &         0.010, 0.0075, 0.0050, 0.0050, 0.0050, 0.0040, 0.0040,&
! &         0.0040, 0.0040, 0.0030, 0.0025, 0.0020, 0.0010, 0.00090,&
! &         0.00075, 0.00075, 0.00075, 0.00065, 0.00050, 0.00040, 0.00025,&
! &         0.00020, 0.00015, 0.00014, 0.000125, 0.000115, 0.00010, 0.00010,&
! &         0.000075, 0.000075, 0.000075/)
!   ystep = (/0.065, 0.040, 0.035, 0.0175, 0.0175, 0.0125, 0.0125,&
! &         0.010, 0.0075, 0.0060, 0.0050, 0.0050, 0.0045, 0.0040,&
! &         0.0035, 0.0030, 0.0030, 0.0020, 0.00175, 0.00100, 0.00090,&
! &         0.00060, 0.00060, 0.00060, 0.00050, 0.00040, 0.00035, 0.00030,&
! &         0.00025, 0.00015, 0.00014, 0.000125, 0.000125, 0.000125, 0.000125,&
! &         0.000125, 0.000125, 0.000125/)
!-------- Settings for the CMS-14-029 analysis ---------!
xmin = (/ 0.0875000, 0.0750000, 0.0500000, 0.0150000, 0.0125000, 0.0087500, 0.0050000,&
&         0.0050000, 0.0025000, 0.0025000, 0.0015000, 0.0015000, 0.0010000, 0.0010000,&
&         0.0010000, 0.0010000, 0.0010000, 0.0003750, 0.0003750, 0.0003000, 0.0003000,&
&         0.0003000, 0.0003000, 0.0003000, 0.0003000, 0.0001750, 0.0001750, 0.0001750,&
&         0.0001750, 0.0001000, 0.0001000, 0.0000750, 0.0000750, 0.0000625, 0.0000625,&
&         0.0000375, 0.0000375, 0.0000300 /)
xmax = (/34.9125000,29.9250000,19.9500000, 5.9850000, 4.9875000, 3.4912500, 1.9950000,&
&         1.9950000, 0.9975000, 0.9975000, 0.5985000, 0.5985000, 0.3990000, 0.3990000,&
&         0.3990000, 0.3990000, 0.3990000, 0.1496250, 0.1496250, 0.1197000, 0.1197000,&
&         0.1197000, 0.1197000, 0.1197000, 0.1197000, 0.0698250, 0.0698250, 0.0698250,&
&         0.0698250, 0.0399000, 0.0399000, 0.0299250, 0.0299250, 0.0249375, 0.0249375,&
&         0.0149625, 0.0149625, 0.0119700/)
ymin = (/ 0.0375000, 0.0250000, 0.0200000, 0.0100000, 0.0075000, 0.0050000, 0.0025000,&
&         0.0025000, 0.0015000, 0.0015000, 0.0010000, 0.0010000, 0.0007500, 0.0007500,&
&         0.0007500, 0.0007500, 0.0007500, 0.0005000, 0.0005000, 0.0003000, 0.0003000,&
&         0.0003000, 0.0003000, 0.0003000, 0.0003000, 0.0001750, 0.0001750, 0.0001000,&
&         0.0001000, 0.0000750, 0.0000750, 0.0000625, 0.0000625, 0.0000500, 0.0000500,&
&         0.0000375, 0.0000375, 0.0000300/)
ymax = (/14.9625000, 9.9750000, 7.9800000, 3.9900000, 2.9925000, 1.9950000, 0.9975000,&
&         0.9975000, 0.5985000, 0.5985000, 0.3990000, 0.3990000, 0.2992500, 0.2992500,&
&         0.2992500, 0.2992500, 0.2992500, 0.1995000, 0.1995000, 0.1197000, 0.1197000,&
&         0.1197000, 0.1197000, 0.1197000, 0.1197000, 0.0698250, 0.0698250, 0.0399000,&
&         0.0399000, 0.0299250, 0.0299250, 0.0249375, 0.0249375, 0.0199500, 0.0199500,&
&         0.0149625, 0.0149625, 0.0119700/)
xstep = (/0.1750000, 0.1500000, 0.1000000, 0.0300000, 0.0250000, 0.0175000, 0.0100000,&
&         0.0100000, 0.0050000, 0.0050000, 0.0030000, 0.0030000, 0.0020000, 0.0020000,&
&         0.0020000, 0.0020000, 0.0020000, 0.0007500, 0.0007500, 0.0006000, 0.0006000,&
&         0.0006000, 0.0006000, 0.0006000, 0.0006000, 0.0003500, 0.0003500, 0.0003500,&
&         0.0003500, 0.0002000, 0.0002000, 0.0001500, 0.0001500, 0.0001250, 0.0001250,&
&         0.0000750, 0.0000750, 0.0000600 /)
ystep = (/0.0750000, 0.0500000, 0.0400000, 0.0200000, 0.0150000, 0.0100000, 0.0050000,&
&         0.0050000, 0.0030000, 0.0030000, 0.0020000, 0.0020000, 0.0015000, 0.0015000,&
&         0.0015000, 0.0015000, 0.0015000, 0.0010000, 0.0010000, 0.0006000, 0.0006000,&
&         0.0006000, 0.0006000, 0.0006000, 0.0006000, 0.0003500, 0.0003500, 0.0002000,&
&         0.0002000, 0.0001500, 0.0001500, 0.0001250, 0.0001250, 0.0001000, 0.0001000,&
&         0.0000750, 0.0000750, 0.0000600 /)
!------------------
  xsize = (/ 200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200/)
  ysize = (/ 200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200,  200,  200,  200,  200,&
&          200,  200,  200/)
!------------------

  
 open(file_id_common4,file = trim(adjustl(pathname))// &
 &    '/Expt_tables/CMS_tautau_llh.binary',form='unformatted')
 read(file_id_common4,iostat=ios) CMS_llhdata

 if(ios.ne.0)then
  do j=1,n2Dslices
   call write_metadata(analysisID,label,description,expt,lumi,energy,particle_type,&
& mass_sep,mass_range_tolerance,masses(j),xmin(j),xmax(j),xstep(j),xsize(j),ymin(j),&
& ymax(j),ystep(j),ysize(j),CMS_llhdata(j))
   call read2Ddata(masses(j),CMS_llhdata(j),stem,'obs',status)
   if(status.ne.0) stop 'Error in setup_likelihoods: Data files not found!'
   call read2Ddata(masses(j),CMS_llhdata(j),stem,'pred',status)
   if(status.ne.0) stop 'Error in setup_likelihoods: Data files not found!'
  enddo
!  write(*,*) 'Creating binaries...'
#ifndef WEBVERSION
  rewind(file_id_common4)
  write(file_id_common4) CMS_llhdata
#endif  
 endif

 close(file_id_common4)
!--------------------------------------------------------------
  end subroutine setup_likelihoods
!--------------------------------------------------------------
 subroutine write_metadata(analysisID,label,description,expt,lumi,energy,particle_type,&
 &                         mass_sep,mass_range_tolerance,mass,xmin,xmax,xstep,xsize,ymin,&
 &                         ymax,ystep,ysize,llh2D)
!--------------------------------------------------------------  
  use usefulbits, only : small
  implicit none
  
  integer, intent(in) :: analysisID,particle_type
  character(LEN=*),intent(in) :: label, description,expt  
  double precision, intent(in) :: mass_sep,mass,xmin,xmax,ymin,ymax,xstep,ystep,lumi,&
 &                                energy,mass_range_tolerance
  integer, intent(in) :: xsize, ysize
  type(likelihood2D), intent(inout) :: llh2D


  llh2D%analysisID = analysisID
  llh2D%label = trim(adjustl(label))
  llh2D%description = trim(adjustl(description))
  llh2D%expt = expt
  llh2D%energy = energy
  llh2D%lumi = lumi
  llh2D%particle_x = particle_type
  llh2D%mass_sep = mass_sep
  llh2D%mass_range_tolerance = mass_range_tolerance
  llh2D%mass = mass
  llh2D%xmin = xmin
  llh2D%xmax = xmax
  llh2D%xstep = xstep
  llh2D%xsize = xsize  
  llh2D%ymin = ymin
  llh2D%ymax = ymax
  llh2D%ystep = ystep
  llh2D%ysize = ysize  

!  if(abs((xmax-xmin)/xstep+1.-xsize).gt.small) then
!   write(*,*) 'Warning: xsize does not match: ',(xmax-xmin)/xstep+1,' vs. ', xsize
!  endif
!  if(abs((ymax-ymin)/ystep+1.-ysize).gt.small) then
!   write(*,*) 'Warning: ysize does not match: ',(ymax-ymin)/ystep+1,' vs. ', ysize
!  endif
!-------------------------------------------------------------- 
 end subroutine write_metadata 
!-------------------------------------------------------------- 
 subroutine read2Ddata(mass,llh2D,stem,obspred,status)
!--------------------------------------------------------------
  use usefulbits, only : small,file_id_common2
  use store_pathname, only : pathname
  implicit none
  
  double precision, intent(in) :: mass
  type(likelihood2D), intent(inout) :: llh2D
  character(LEN=*), intent(in) :: stem,obspred
  integer, intent(out) :: status
  
  integer :: size,i,j,k,remove,posx,posy
  double precision, dimension(1:3) :: values
  type(threetuple), allocatable :: dummylikelihood(:), likelihood(:)
  double precision, dimension(200,200) :: llhgrid  
  double precision, dimension(1:3) :: bestfit
  character(LEN=4) :: intstring
  character(LEN=100) ::  filename
  double precision :: modx, mody  

  write(intstring,"(I4)") int(mass)  
  if(trim(adjustl(obspred)).eq.'obs') then
!   filename = 'L_data_SMHb_'//trim(adjustl(intstring))//'.out'
   filename = 'L_data_b_'//trim(adjustl(intstring))//'.out'   
  else if(trim(adjustl(obspred)).eq.'pred') then
!   filename = 'L_asimovSMH_SMHb_'//trim(adjustl(intstring))//'.out'
   filename = 'L_asimov_b_'//trim(adjustl(intstring))//'.out'   
  else
   stop 'Error in read2Ddata: obspred unknown.'   
  endif 

  llh2D%mass = mass

  call get_length_of_file(trim(adjustl(pathname))//trim(adjustl(stem))// &
  &                       trim(adjustl(filename)),size,status)
  if(status.ne.0) return
  
  open(file_id_common2,file=trim(adjustl(pathname))//trim(adjustl(stem))// &
  &                         trim(adjustl(filename)))
     
  allocate(dummylikelihood(size))  
  
  k=0
  remove=0
  do i=1, size
   read(file_id_common2,*) values
   if(values(3).eq.0.0D0) then
    bestfit = values
    remove = remove + 1
   else 
    k=k+1
    dummylikelihood(k)%val = values
   endif 
  enddo

  close(file_id_common2)

  llh2D%size = k
!  write(*,*) 'size = ',llh2D%size
    
  allocate(likelihood(k))
  
  do i=1,k
   likelihood(i)%val = dummylikelihood(i)%val
  enddo 
  
  deallocate(dummylikelihood)


!  write(*,*) "Reading ", obspred, "grid with mass = ", mass  
  llhgrid = 0.0D0
  
  do i=lbound(likelihood,dim=1),ubound(likelihood,dim=1)

   posx = nint((likelihood(i)%val(1)-llh2D%xmin)/llh2D%xstep)+1
   posy = nint((likelihood(i)%val(2)-llh2D%ymin)/llh2D%ystep)+1
!  write(*,*) likelihood(i)%val(1), llh2D%xmin, posx
   llhgrid(posx,posy) = likelihood(i)%val(3)     
  enddo
 
!  k=0
!  do i=1,llh2D%xsize
!   do j=1, llh2D%ysize
!    if(llhgrid(i,j).lt.small) then
!     k=k+1
!!     write(*,*) i,j,llh2D%llhgrid(i,j)
!    endif
!   enddo
!  enddo   
!  write(*,*) k, 'grid points are unfilled!'    
 
  if(trim(adjustl(obspred)).eq.'obs') then
   llh2D%llhgrid_obs = llhgrid
   llh2D%bestfit_obs = bestfit
  else if(trim(adjustl(obspred)).eq.'pred') then
   llh2D%llhgrid_pred = llhgrid
   llh2D%bestfit_pred = bestfit
  endif 
 
  deallocate(likelihood)
!--------------------------------------------------------------     
 end subroutine read2Ddata
 !--------------------------------------------------------------
 subroutine get_likelihood(analysisID, jj, t, llh, M_av, nc, cbin, obspred, cbin_in)
 !--------------------------------------------------------------
  use usefulbits, only : dataset
  implicit none
 
  type(dataset), intent(in) :: t
  integer, intent(in) :: analysisID, jj
  character(LEN=*), intent(in) :: obspred  
  double precision, intent(out) :: llh,M_av    
  integer, intent(out) :: nc, cbin
  integer, optional, intent(in) :: cbin_in   
  double precision :: cfact1,cfact2
  
  llh = 0.0D0
    
   select case(analysisID)
    case(3316,14029)
     if(present(cbin_in)) then
!      write(*,*) "Calling get_likelihood with cbin_in = ", cbin_in, ", obspred = ", obspred
      call calcfact_llh(1,jj,t,cfact1,cfact2,M_av,nc,cbin,cbin_in)
     else
      call calcfact_llh(1,jj,t,cfact1,cfact2,M_av,nc,cbin)     
     endif
     if(cbin.ne.0) then
      llh = get_llh_CMS_tautau(M_av, cfact1, cfact2, obspred)
     endif 
!     write(*,*) jj, cfact1, cfact2, M_av, nc, cbin, "llh = ", llh
!     if(nc.eq.3) then
!      write(*,*) "nc=3: ", jj, cfact1, cfact2, M_av, llh
!     endif      
    case default
     stop 'Error: Unknown analysis ID to retrieve chi^2 value.'
   end select   
   
 end subroutine get_likelihood
  
  !--------------------------------------------------------------
 subroutine calcpredratio_llh(c,jj,t,M_av,fact,nc,predratio)
 !--------------------------------------------------------------
  use usefulbits, only : dataset,small,vsmall
  implicit none
    
 !prsep(h,n)%tlist,prsep(h,n)%findj,t,axis_i(n),ncomb(n),predratio(n)
  type(dataset), intent(in) :: t
  integer, intent(in) :: c,jj 
  double precision, intent(out) :: M_av, predratio, fact
  integer, intent(out) ::  nc
 
  double precision :: cfact1,cfact2
  double precision :: sf, expllh, interval,expllh0
  double precision :: diff(2)
  integer :: kk, cbin

  double precision sf_low,sf_high,llh_low,llh_high, sf_max


  call calcfact_llh(c,jj,t,cfact1,cfact2,M_av,nc,cbin)
  fact = cfact1 + cfact2
 
  expllh = get_llh_CMS_tautau(M_av, cfact1, cfact2, 'pred')
  expllh0 = expllh
!  write(*,*) 'c = ',c,' jj = ',jj,'nc = ',nc,' M_av = ', M_av, 'ggH = ',cfact1, 'bbH = ',cfact2
  
  
  if(cfact1.lt.small.and.cfact2.lt.small) then
!   write(*,*) 'Very small cross sections. CMS tautau is not competitive, expected llh =',expllh
   predratio = -1.0D0   
   return
  endif
     
!  write(*,*) "Find predicted ratio..."
  if(expllh.lt.0.0D0) then
   predratio = -1.0D0
   return
! NEW
  endif
 
!  sf = 1.0D0
!  kk = 0
!  interval = 0.0005D0  
!  interval = 0.50D0  
!  diff = 0.0D0  
  
! TS implementation  
!  do while(abs(expllh-5.99D0).gt.0.01D0)
!   kk = kk + 1
!   if(abs(diff(2)-diff(1)).le.vsmall) then
!    interval = 0.5D0*sf
!   else    
!    interval = min(abs(rand(0)*interval/(diff(2)-diff(1))),0.5D0*sf)
!   endif 
!   
!   sf = sf + sign(1.0D0,5.99D0-expllh)*interval
!   if (sf.le.0.0D0) exit    
!   
!   diff(2)=abs(expllh-5.99D0)
!   
!   expllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'pred')
!   
!   diff(1)=abs(expllh-5.99D0)
!   
!!   write(*,*) kk, sf, interval, diff(1), diff(2), expllh
!  enddo
  

! OS implementation of bisection
! Maximum scale factor corresponds to min_predratio=1/sf_max
  sf_max = 1.d0
  llh_high = get_llh_CMS_tautau(M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
  kk=0
  do while(llh_high-5.99d0.lt.0d0.and.kk.lt.100)
   sf_max = sf_max*10.d0
   llh_high = get_llh_CMS_tautau(M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
   kk=kk+1
  enddo


  sf_low = 0.d0
  sf_high = sf_max
  kk=0  

  sf = (sf_high+sf_low)/2.

  llh_low = get_llh_CMS_tautau(M_av, sf_low*cfact1, sf_low*cfact2, 'pred')
  llh_high = get_llh_CMS_tautau(M_av, sf_high*cfact1, sf_high*cfact2, 'pred')

  ! Check that root exists in starting interval
  if ((llh_low-5.99)*(llh_high-5.99).ge.0) then
    sf = sf_max
  else
   expllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'pred')
   
!   write(*,*) kk, sf_low, sf_high, sf, expllh  
   do while(abs(expllh-5.99D0).gt.0.01D0.and.kk.lt.100)
   
    ! Five digit precision on scale factor in case of oscillating solutions
    if ((sf_high-sf_low).LT.1D-5) exit     
   
    kk = kk + 1
    if((expllh-5.99d0).lt.0) sf_low = sf
    if((expllh-5.99d0).gt.0) sf_high= sf
    sf = (sf_high+sf_low)/2.
   
    expllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'pred')
      
!    write(*,*) kk, sf_low, sf_high, sf, expllh
   enddo
  endif
  
  if( kk.ge.100) write(*,*) 'Warning: Maximum number of iterations reached with no convergence: predratio might be unreliable.'
!      write(*,*) kk, sf_low, sf_high, sf, expllh

    
!  write(*,*) "Ended predratio finder for h",jj," (",nc," combined, average mass = "&
!  &   ,M_av,") with ",kk,"steps: scalefactor = ", sf,", original -2ln(L) = ",expllh0    
! --

!  elseif(expllh.gt.5.99D0) then
!   do while(expllh.gt.5.99D0)
!
!! Adjust scanning intervals to distance
!    kk = kk + 1 
!    if(abs(expllh-5.99D0) > 2.5D0) then
!     interval = 0.05D0
!    elseif(abs(expllh-5.99D0) > 0.5D0) then
!     interval = 0.01D0
!    elseif(abs(expllh-5.99D0) > 0.25D0) then
!     interval = 0.001D0 
!    elseif(abs(expllh-5.99D0) > 0.05D0) then
!     interval = 0.0005D0 
!    endif
!     
!    sf = sf - interval
!    if (sf.le.0.0D0) exit
!    
!    expllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'pred')
!    
!    write(*,*) kk, sf, interval, expllh
!   enddo
!  elseif(expllh.lt.5.99D0) then
!   do while(expllh.lt.5.99D0)
!
!! Adjust scanning intervals to distance
!    kk = kk + 1 
!    if(abs(expllh-5.99D0) > 2.5D0) then
!     interval = 0.05D0
!    elseif(abs(expllh-5.99D0) > 0.5D0) then
!     interval = 0.01D0
!    elseif(abs(expllh-5.99D0) > 0.25D0) then
!     interval = 0.001D0 
!    elseif(abs(expllh-5.99D0) > 0.05D0) then
!     interval = 0.0005D0 
!    endif
!     
!    sf = sf + interval
!    expllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'pred')
!    
!    write(*,*) kk, sf, interval, expllh
!   enddo    
!  endif
  
  if(sf.le.0.0D0) sf = 1.0D-10

  predratio = 1./sf 
 
 end subroutine calcpredratio_llh 
 !--------------------------------------------------------------
 subroutine check_against_bound_llh(c,jj,t,M_av,nc,obsratio) 
 !-------------------------------------------------------------- 
  use usefulbits, only : dataset,vsmall
  implicit none
  
  type(dataset), intent(in) :: t
  integer, intent(in) :: c,jj , nc
  double precision, intent(in) :: M_av
  double precision, intent(out) :: obsratio
 
  double precision :: cfact1,cfact2
  double precision :: sf, obsllh,M_av2, interval,obsllh0
  double precision :: diff(2)  
  integer :: kk,nc2,cbin
  
  double precision sf_low,sf_high,llh_low,llh_high,sf_max

  call calcfact_llh(c,jj,t,cfact1,cfact2,M_av2,nc2,cbin)
 
  if(abs(M_av2-M_av).gt.vsmall.or.nc.ne.nc2) then
   write(*,*) M_av2, M_av, nc, nc2
   stop 'Error in subroutine check_against_bound_llh !'
   
  endif 
 
  obsllh = get_llh_CMS_tautau(M_av, cfact1, cfact2, 'obs')
  obsllh0 = obsllh
  
!  write(*,*) "Find observed ratio..."
  if(obsllh.lt.0.0D0) then
   obsratio = -1.0D0
   return
  endif 
   
  sf = 1.0D0
  kk = 0
  interval = 0.50D0  
  diff = 0.0D0  

!! TS implementation  
!  do while(abs(obsllh-5.99D0).gt.0.05D0)
!   kk = kk + 1
!   if(abs(diff(2)-diff(1)).le.vsmall) then
!    interval = 0.5D0*sf
!   else    
!    interval = min(abs(rand(0)*interval/(diff(2)-diff(1))),0.5D0*sf)
!   endif 
!   
!   sf = sf + sign(1.0D0,5.99D0-obsllh)*interval
!   if (sf.le.0.0D0) exit    
!   
!   diff(2)=abs(obsllh-5.99D0)
!   
!   obsllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'obs')
!   
!   diff(1)=abs(obsllh-5.99D0)
!
!   write(*,*) kk, sf, interval, diff(1), diff(2), obsllh
!   
!  enddo


!! OS implementation of bisection
  sf_max = 1.d0
  llh_high = get_llh_CMS_tautau(M_av, sf_max*cfact1, sf_max*cfact2, 'obs')
  kk=0
  do while(llh_high-5.99d0.lt.0d0.and.kk.lt.100)
   sf_max = sf_max*10.d0
   llh_high = get_llh_CMS_tautau(M_av, sf_max*cfact1, sf_max*cfact2, 'obs')
   kk=kk+1
  enddo

  sf_low = 0.d0
  sf_high = sf_max
  kk=0

  sf = (sf_high+sf_low)/2.

  llh_low = get_llh_CMS_tautau(M_av, sf_low*cfact1, sf_low*cfact2, 'obs')
  llh_high = get_llh_CMS_tautau(M_av, sf_high*cfact1, sf_high*cfact2, 'obs')

  ! Check that root exists in starting interval
  if ((llh_low-5.99d0)*(llh_high-5.99d0).ge.0) then
    sf = sf_max
  else
   obsllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'obs')
!   write(*,*) kk, sf_low, sf_high, sf, obsllh
     
   do while(abs(obsllh-5.99D0).gt.0.01D0.and.kk.lt.100)
   
    ! Five digit precision on scale factor in case of oscillating solutions
    if ((sf_high-sf_low).LT.1D-5) exit     
   
    kk = kk + 1
    if((obsllh-5.99d0).lt.0) sf_low = sf
    if((obsllh-5.99d0).gt.0) sf_high= sf
    sf = (sf_high+sf_low)/2.
   
    obsllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'obs')
      
!    write(*,*) kk, sf_low, sf_high, sf, obsllh
   enddo
  endif

  if( kk.ge.100) write(*,*) 'Warning: Maximum number of iterations reached with no convergence: obsratio might be unreliable.'
  
!  write(*,*) "Ended obsratio finder with ",kk,"steps: ", sf,obsllh, ", original -2ln(L) = ",obsllh0   
  
     
!  elseif(obsllh.gt.5.99D0) then
!   do while(obsllh.gt.5.99D0)
!
!! Adjust scanning intervals to distance
!    kk = kk + 1 
!    if(abs(obsllh-5.99D0) > 3.0D0) then
!     interval = 0.05D0
!    elseif(abs(obsllh-5.99D0) > 0.5D0) then
!     interval = 0.01D0
!    elseif(abs(obsllh-5.99D0) > 0.25D0) then
!     interval = 0.001D0 
!    elseif(abs(obsllh-5.99D0) > 0.05D0) then
!     interval = 0.0005D0 
!    endif
!     
!    sf = sf - interval
!    if(sf.le.0.0D0) exit
!    obsllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'obs')
!    
!    write(*,*) kk, sf, interval, obsllh
!   enddo
!  elseif(obsllh.lt.5.99D0) then
!   do while(obsllh.lt.5.99D0)
!
!! Adjust scanning intervals to distance
!    kk = kk + 1 
!    if(abs(obsllh-5.99D0) > 3.0D0) then
!     interval = 0.05D0
!    elseif(abs(obsllh-5.99D0) > 0.5D0) then
!     interval = 0.01D0
!    elseif(abs(obsllh-5.99D0) > 0.25D0) then
!     interval = 0.001D0 
!    elseif(abs(obsllh-5.99D0) > 0.05D0) then
!     interval = 0.0005D0 
!    endif
!     
!    sf = sf + interval
!    obsllh = get_llh_CMS_tautau(M_av, sf*cfact1, sf*cfact2, 'obs')
!    
!    write(*,*) kk, sf, interval, obsllh
!   enddo    
!  endif
  
  if(sf.le.0.0D0) sf = 1.0D-10  
  obsratio = 1./sf 
   
 end subroutine check_against_bound_llh
!-------------------------------------------------------------- 
 subroutine calcfact_llh(c,jj,t,cfact1,cfact2,M_av,nc,cbin,cbin_in)
 !-------------------------------------------------------------- 
 ! This routine calculates the predicted rates for each process
 ! Need to generalize this if more likelihoods become available.
  use usefulbits, only : dataset, np, div, vvsmall
  use theory_BRfunctions
  use theory_XS_SM_functions 
    
  implicit none      
  !--------------------------------------input
  type(dataset), intent(in) :: t
  integer, intent(in) :: c,jj
  integer, optional, intent(in) :: cbin_in
  !-----------------------------------output
  double precision, intent(out) :: cfact1, cfact2 ,M_av
  integer, intent(out) :: nc,cbin
  !-----------------------------------internal
  double precision,allocatable :: mass(:),fact1(:),fact2(:)
  integer :: npart,f,j !number of particles  
  double precision :: numerator, denominator
  logical, allocatable :: Havail(:)
  
  cfact1=0.0D0    
  cfact2=0.0D0    
  
  if(c.eq.1) then ! this is the CMS tautau likelihood
   npart=np( CMS_llhdata(1)%particle_x )
   allocate(mass(npart),fact1(npart),fact2(npart))

   mass(:)=t%particle( CMS_llhdata(1)%particle_x )%M(:) 
   fact1= 0.0D0
   fact2= 0.0D0      

  allocate(Havail(npart))
  
  if(present(cbin_in)) then
   call available_Higgses(Havail,npart,cbin_in)
  else
   Havail = .True.
  endif 

   
!   write(*,*) 'calcfact_llh: mass = ', mass
   cbin = 0
   
   !write(*,*) "Calling calcfact_llh..."
   
   do j=1,npart
    if(Havail(j)) then    
!   if(abs(mass(jj)-mass(j)).le.CMS_llhdata(1)%mass_sep*mass(jj))then
!   if((abs(mass(jj)-mass(j)).le.CMS_llhdata(1)%mass_sep*100.0D0)) then 
   if(abs(mass(jj)-mass(j)).le.CMS_llhdata(1)%mass_sep*max(mass(jj),mass(j)))then     
      if((t%gg_hj_ratio(j).ge.0.0D0).and.(t%bb_hj_ratio(j).ge.0.0D0)) then
!      write(*,*) "Using partonic input for CMS tautau analysis..."
       fact1(j)=t%gg_hj_ratio(j)*XS_lhc8_gg_H_SM(mass(j))*t%BR_hjtautau(j)
       fact2(j)=t%bb_hj_ratio(j)*XS_lhc8_bb_H_SM(mass(j))*t%BR_hjtautau(j)  
      else
!      write(*,*) "Using hadronic input for CMS tautau analysis..."
       fact1(j)=t%lhc8%XS_hj_ratio(j)*XS_lhc8_gg_H_SM(mass(j))*t%BR_hjtautau(j)
       fact2(j)=t%lhc8%XS_hjb_ratio(j)*XS_lhc8_bb_H_SM(mass(j))*t%BR_hjtautau(j)  
      endif
      cbin = cbin + 2**(j-1)
!     write(*,*) cbin
     endif
    endif 
   enddo
   
!   write(*,*) 'calcfact_llh: fact1 = ', fact1
!   write(*,*) 'calcfact_llh: fact2 = ', fact2
   
  if(fact1(jj).le.vvsmall.and.fact2(jj).le.vvsmall)then
!Higgs jj doesn't contribute - wait until another call of this subroutine before
!looking at nearby masses
   M_av = mass(jj)
   nc=0
   cfact1=0.0D0   
   cfact2=0.0D0      
  else
!find M_av weighted by the rates (only using higgs which have non-zero fact): 
   f=0  
   numerator=0.0D0
   denominator=0.0D0
   do j=1,npart    
    if(   fact1(j).gt.vvsmall.or.fact2(j).gt.vvsmall  )then       
     f=f+1
     numerator=numerator+mass(j)*(fact1(j)+fact2(j))
     denominator=denominator+(fact1(j)+fact2(j))     
    endif
   enddo 

  if(denominator.ne.0.0D0) then
   M_av = numerator/denominator
  else
   M_av = 0.0D0
   write(*,*) 'Warning: could not find average mass for CMS tautau analysis!'
  endif  

  nc=f !f will always be > 0 because we've already made sure that fact(jj)>0.0D0
  cfact1=sum(fact1)
  cfact2=sum(fact2)
  endif

  deallocate(mass,fact1,fact2)
  deallocate(Havail)
     
  endif 
 end subroutine calcfact_llh 
!--------------------------------------------------------------
 function get_llh_from_interpolation(llh2D, xval0, yval0,obspred)
 ! This routine returns the - 2 Delta ln(llh) from grid interpolation
!-------------------------------------------------------------- 
 use usefulbits, only : small
  type(likelihood2D), intent(in) :: llh2D  
  double precision, intent(in) :: xval0, yval0
  character(LEN=*), intent(in) :: obspred  
  double precision :: get_llh_from_interpolation
  integer :: posx1, posx2, posy1, posy2, posxlow, posxup
  double precision :: xval, yval, xval1, xval2, yval1, yval2, fvaly1, fvaly2, fval
  integer :: i, max_expansion
  logical :: lowerxfound, upperxfound
  double precision, dimension(200,200) :: llhgrid
! Josh Sayre extension --
  double precision :: truex, truey
! --

  if(trim(adjustl(obspred)).eq.'obs') then
   llhgrid = llh2D%llhgrid_obs
  else if(trim(adjustl(obspred)).eq.'pred') then
   llhgrid = llh2D%llhgrid_pred
  else
   stop 'Error in get_llh_from_interpolation: obspred unknown.'   
  endif 


  ! This number is the maximal number of tries, that the algorithm steps further up or
  ! downwards the grid in case of non-existent grid points.
  max_expansion = 3

  if(isnan(xval0)) then
   xval = 0.0D0
  else 
   xval = xval0
  endif
   
  if(isnan(yval0)) then
   yval = 0.0D0
  else 
   yval = yval0
  endif

! write(*,*) 'input x,y vals = ', xval0, xval, yval0,  yval

! Josh Sayre extension --
  truex=xval
  truey=yval
! --

 if(xval.lt.llh2D%xmin) then
!  write(*,*) 'Warning: ggH rate outside grid -- value too small.'
  xval = llh2D%xmin
 else if (xval.gt.llh2D%xmax) then
!  write(*,*) 'Warning: ggH rate outside grid -- value too large.'
  xval = llh2D%xmax-small 
 endif 
 if(yval.lt.llh2D%ymin) then
!  write(*,*) 'Warning: bbH rate outside grid -- value too small.'
  yval = llh2D%ymin
 else if (yval.gt.llh2D%ymax) then
!  write(*,*) 'Warning: ggH rate outside grid -- value too large.'
  yval = llh2D%ymax -small
 endif 

! Get coordinates of neighboring points
 posx1 = (xval-llh2D%xmin) / llh2D%xstep + 1
 if(posx1.lt.llh2D%xsize) then
  posx2 = (xval-llh2D%xmin) / llh2D%xstep + 2
 else
  write(*,*) "Warning: on the edge of the grid in x-position!"
  posx1 = (xval-llh2D%xmin) / llh2D%xstep
  posx2 = (xval-llh2D%xmin) / llh2D%xstep + 1  
 endif

 posy1 = (yval-llh2D%ymin) / llh2D%ystep + 1 
 if(posy1.lt.llh2D%ysize) then 
  posy2 = (yval-llh2D%ymin) / llh2D%ystep + 2
 else
  write(*,*) "Warning: on the edge of the grid in y-position!" 
  posy1 = (yval-llh2D%ymin) / llh2D%ystep
  posy2 = (yval-llh2D%ymin) / llh2D%ystep + 1
 endif
  
 ! Get x,y values at these points
 xval1 = (posx1 - 1)*llh2D%xstep + llh2D%xmin
 xval2 = (posx2 - 1)*llh2D%xstep + llh2D%xmin
 yval1 = (posy1 - 1)*llh2D%ystep + llh2D%ymin
 yval2 = (posy2 - 1)*llh2D%ystep + llh2D%ymin
 
 ! Do bilinear interpolation and retrieve likelihood value
 if(llhgrid(posx1,posy1).ne.0.0D0.and.llhgrid(posx2,posy1).ne.0.0D0) then
! Josh Sayre extension --  
  fvaly1 = (xval2 - truex) / (xval2 - xval1) * llhgrid(posx1,posy1) + &
  &        (truex - xval1) / (xval2 - xval1) * llhgrid(posx2,posy1)
! --
!   fvaly1 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posx1,posy1) + &
!   &        (xval - xval1) / (xval2 - xval1) * llhgrid(posx2,posy1)
 else if (llhgrid(posx1,posy1).ne.0.0D0) then
  fvaly1 = llhgrid(posx1,posy1)
 else if (llhgrid(posx2,posy1).ne.0.0D0) then
  fvaly1 = llhgrid(posx2,posy1)
 else
!  write(*,*) 'Warning in interpolating grid: Both gridpoints in x-direction missing!'
!  write(*,*) '(', xval1, yval1, '), (', xval2, yval1,')'
  lowerxfound = .False.
  upperxfound = .False.  
  do i=1,max_expansion
   if(posx1-i.gt.0) then
    if (llhgrid(posx1-i,posy1).ne.0.0D0) then
      xval1 = (posx1 - i - 1)*llh2D%xstep + llh2D%xmin
      lowerxfound = .True.
      posxlow = posx1-i
     endif 
    endif
   if(posx2+i.le.ubound(llhgrid,dim=1)) then 
    if (llhgrid(posx2+i,posy1).ne.0.0D0) then
      xval2 = (posx2 + i - 1)*llh2D%xstep + llh2D%xmin
      upperxfound = .True.
      posxup = posx2+i
     endif 
    endif
    if(lowerxfound.or.upperxfound) exit
   enddo 
   if(lowerxfound.and.upperxfound) then
! Josh Sayre extension --  
    fvaly1 = (xval2 - truex) / (xval2 - xval1) * llhgrid(posxlow,posy1) + &
    &        (truex - xval1) / (xval2 - xval1) * llhgrid(posxup,posy1)
! --    
!     fvaly1 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posxlow,posy1) + &
!     &        (xval - xval1) / (xval2 - xval1) * llhgrid(posxup,posy1)
   else if(lowerxfound) then
    fvaly1 = llhgrid(posxlow,posy1)     
   else if(upperxfound) then
    fvaly1 = llhgrid(posxup,posy1)     
   else
    write(*,*) 'Found no surrounding points!'
    fvaly1 = 99999.0D0
   endif   
 endif  
   
 if(llhgrid(posx1,posy2).ne.0.0D0.and.llhgrid(posx2,posy2).ne.0.0D0) then 
! Josh Sayre extension --  
  fvaly2 = (xval2 - truex) / (xval2 - xval1) * llhgrid(posx1,posy2) + &
  &        (truex - xval1) / (xval2 - xval1) * llhgrid(posx2,posy2)
! -- 
!   fvaly2 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posx1,posy2) + &
!   &        (xval - xval1) / (xval2 - xval1) * llhgrid(posx2,posy2)
 else if (llhgrid(posx1,posy2).ne.0.0D0) then
  fvaly2 = llhgrid(posx1,posy2)
 else if (llhgrid(posx2,posy2).ne.0.0D0) then
  fvaly2 = llhgrid(posx2,posy2)
 else
!  write(*,*) 'Warning in interpolating grid: Both gridpoints in x-direction missing!'
!  write(*,*) '(', xval1, yval2, '), (', xval2, yval2,')' 
  lowerxfound = .False.
  upperxfound = .False.  
  do i=1,max_expansion
   if(posx1-i.gt.0) then
    if (llhgrid(posx1-i,posy2).ne.0.0D0) then
      xval1 = (posx1 - i - 1)*llh2D%xstep + llh2D%xmin
      lowerxfound = .True.
      posxlow = posx1-i
     endif 
    endif
   if(posx2+i.le.ubound(llhgrid,dim=1)) then 
    if (llhgrid(posx2+i,posy2).ne.0.0D0) then
      xval2 = (posx2 + i - 1)*llh2D%xstep + llh2D%xmin
      upperxfound = .True.
      posxup = posx2+i
     endif 
    endif
!    write(*,*) i, lowerxfound, upperxfound
    if(lowerxfound.or.upperxfound) exit
   enddo 
   if(lowerxfound.and.upperxfound) then
! Josh Sayre extension --     
    fvaly2 = (xval2 - truex) / (xval2 - xval1) * llhgrid(posxlow,posy2) + &
    &        (truex - xval1) / (xval2 - xval1) * llhgrid(posxup,posy2)
! --   
!     fvaly2 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posxlow,posy2) + &
!     &        (xval - xval1) / (xval2 - xval1) * llhgrid(posxup,posy2)
   else if(lowerxfound) then
    fvaly2 = llhgrid(posxlow,posy2)     
   else if(upperxfound) then
    fvaly2 = llhgrid(posxup,posy2)     
   else
!    write(*,*) 'Found no surrounding points!'
    fvaly2 = 99999.0D0
   endif   
 endif  
 
 fval =  (yval2 - yval) / (yval2 - yval1) * fvaly1 + &
 &       (yval - yval1) / (yval2 - yval1) * fvaly2
  
 ! we want to obtain - 2 * Delta log(likelihood) \simeq chi^2
 get_llh_from_interpolation = 2*fval
 return
 
 end function get_llh_from_interpolation
!--------------------------------------------------------------  
 function get_llh_CMS_tautau(mpred, rate1, rate2, obspred)
!-------------------------------------------------------------- 
  use theory_XS_SM_functions, only : XS_lhc8_gg_H_SM, XS_lhc8_bb_H_SM
  use theory_BRfunctions, only : BRSM_Htautau
  
  implicit none
 
  double precision, intent(in) :: mpred, rate1, rate2
  character(LEN=*), intent(in) :: obspred

  double precision :: get_llh_CMS_tautau
  double precision :: ggH1, bbH1, llh1, ggH2, bbH2, llh2, mpred_tmp
  integer :: i,j

  mpred_tmp = mpred

  if (mpred_tmp.lt.(CMS_llhdata(1)%mass-CMS_llhdata(1)%mass_range_tolerance)) then
!   write(*,*) 'Warning: predicted mass below lowest mass value of CMS analysis.'
   get_llh_CMS_tautau = 0.0D0
  else if (mpred_tmp.gt.(CMS_llhdata(n2Dslices)%mass+CMS_llhdata(1)%mass_range_tolerance)) then
!   write(*,*) 'Warning: predicted mass above highest mass value of CMS analysis.'
   get_llh_CMS_tautau = 0.0D0
  else 
   do i=lbound(CMS_llhdata,dim=1),ubound(CMS_llhdata,dim=1)-1
    if(mpred_tmp.lt.CMS_llhdata(1)%mass) then
     mpred_tmp = CMS_llhdata(1)%mass
    else if(mpred_tmp.gt.CMS_llhdata(n2Dslices)%mass) then
     mpred_tmp = CMS_llhdata(n2Dslices)%mass
    endif      
    if(((mpred_tmp - CMS_llhdata(i)%mass).ge.0.0D0).and.(mpred_tmp-CMS_llhdata(i+1)%mass.lt.0.0D0)) then
     ! Rescale to grid mass values by SM rate ratios
     if(rescale) then
      ggH1 = rate1 * XS_lhc8_gg_H_SM(CMS_llhdata(i)%mass) * BRSM_Htautau(CMS_llhdata(i)%mass) / &
&                  ( XS_lhc8_gg_H_SM(mpred_tmp)           * BRSM_Htautau(mpred_tmp) )
      bbH1 = rate2 * XS_lhc8_bb_H_SM(CMS_llhdata(i)%mass) * BRSM_Htautau(CMS_llhdata(i)%mass) / &
&                  ( XS_lhc8_bb_H_SM(mpred_tmp)           * BRSM_Htautau(mpred_tmp) )     
     else
	  ggH1 = rate1      
      bbH1 = rate2
     endif 
     llh1 = get_llh_from_interpolation(CMS_llhdata(i),ggH1,bbH1,obspred)

	 if(rescale) then     
      ggH2 = rate1*XS_lhc8_gg_H_SM(CMS_llhdata(i+1)%mass)*BRSM_Htautau(CMS_llhdata(i+1)%mass) / &
&                ( XS_lhc8_gg_H_SM(mpred_tmp)           * BRSM_Htautau(mpred_tmp) )
      bbH2 = rate2*XS_lhc8_bb_H_SM(CMS_llhdata(i+1)%mass)*BRSM_Htautau(CMS_llhdata(i+1)%mass) / &
&                ( XS_lhc8_bb_H_SM(mpred_tmp)           * BRSM_Htautau(mpred_tmp) )     
     else
      ggH2 = rate1
      bbH2 = rate2
     endif
      
     llh2 = get_llh_from_interpolation(CMS_llhdata(i+1),ggH2,bbH2,obspred)

      get_llh_CMS_tautau = llh1 + (mpred_tmp - CMS_llhdata(i)%mass)/(CMS_llhdata(i+1)%mass-CMS_llhdata(i)%mass) * &
&                      (llh2 - llh1)
     exit
     else if(mpred_tmp-CMS_llhdata(i+1)%mass.eq.0.0D0) then
      get_llh_CMS_tautau = get_llh_from_interpolation(CMS_llhdata(i+1),rate1,rate2,obspred)
     endif 
   enddo
  endif
  return
!-------------------------------------------------------------- 
 end function get_llh_CMS_tautau 
!-------------------------------------------------------------- 
 subroutine outputproc_llh(proc,file_id_common,descrip) 
!--------------------------------------------------------------  
 use usefulbits, only : listprocesses
 
 type(listprocesses),intent(in) :: proc
 integer, intent(in) :: file_id_common
 character(LEN=200), intent(out) :: descrip

 character(LEN=45) :: label
 integer :: i,j
 character(LEN=1) :: jj

 i=proc%findi
 j=proc%findj

 write(jj,'(I1)')j
 
 label='('//trim(CMS_llhdata(1)%label)//')'

 descrip = '(pp)->h'//jj//'->tautau, using -2ln(L) reconstruction  '//label
 
 end subroutine outputproc_llh 
!--------------------------------------------------------------
 subroutine get_length_of_file(filename, n, status)
!-------------------------------------------------------------- 
  use usefulbits, only : file_id_common3
  implicit none

  character(LEN=*),intent(in) :: filename
  integer, intent(out) :: n, status
  integer :: error
 
  open(file_id_common3,file=trim(adjustl(filename)),action='read',status='old',iostat=status)
  if(status.ne.0) then
   write(*,*) 'Bad status', status, 'with the following file:'
   write(*,*) trim(adjustl(filename))
  endif 
  n = 0
   do
    read(file_id_common3,*, iostat = error)
    if (error == -1) exit
    n = n + 1
   enddo
  close(file_id_common3)
 end subroutine get_length_of_file
!-------------------------------------------------------------- 
 subroutine available_Higgses(Havail,nH,cbin_in)
  implicit none
  
  integer, intent(in) :: nH, cbin_in
  logical, intent(inout) :: Havail(:)
  
  integer :: c,i,r
  
  if(size(Havail,dim=1).ne.nH) then
   write(*,*) "Warning in subroutine available_Higgses: sizes do not match!"
   Havail = .True.
  else
   c = cbin_in
   do i=1,nH
    r = c / (2**(nH-i))
    if(r.eq.1) then
     Havail(nH-i+1) = .False.
    elseif (r.eq.0) then
     Havail(nH-i+1) = .True.
    else
     stop "Error in subroutine available_Higgses: not valid binary!"      
    endif 
    c = mod(c,2**(nH-i))
   enddo
!  write(*,*) "cbin = ",cbin_in, " Havail = ",Havail
  endif 
 
 end subroutine available_Higgses
!******************************************************************
end module likelihoods
!******************************************************************