! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type2
!******************************************************************
 use S95tables_type1
 implicit none
!#define fast
 
 !table type 2------------------------------
 type table2
  logical :: needs_M2_gt_2M1
  integer :: id,nx1,nx2,particle_x1,particle_x2 !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
  character(LEN=45) :: label
  character(LEN=3) :: expt      
  double precision :: lumi, energy  
  double precision :: xmax1,xmin1,xmax2,xmin2,sep1,sep2,deltax
  double precision, allocatable :: dat(:,:,:) !in dat(a,b,1:2)...obs=1,pred=2. 1st component of dat is y, 2nd is x 
  double precision :: maxdatval    
  double precision :: z !only used in slices_t2
 end type

 integer,parameter :: file_id_2_exp=10  !same as file_id_common in usefulbits.f90
 integer,parameter :: file_id_2_obs=11
                                    
 !------------------------------------------ 

 contains


 !************************************************************ 
 subroutine initializetables_type2_blank(tablet2)
 !***********************************************************  
 ! still leaves dat unallocated
  integer:: i
  type(table2) :: tablet2(:)

  do i=lbound(tablet2,dim=1),ubound(tablet2,dim=1)
   tablet2(i)%id          = -1
   tablet2(i)%nx1         = -1
   tablet2(i)%nx2         = -1
   tablet2(i)%particle_x1 = -1
   tablet2(i)%particle_x2 = -1
   tablet2(i)%label       = ''
   tablet2(i)%expt        = ''
   tablet2(i)%lumi       = -1.0D0
   tablet2(i)%energy     = -1.0D0   
   tablet2(i)%xmax1       = -1.0D0
   tablet2(i)%xmax2       = -1.0D0
   tablet2(i)%xmin1       = -1.0D0
   tablet2(i)%xmin2       = -1.0D0
   tablet2(i)%sep1        = -1.0D0
   tablet2(i)%sep2        = -1.0D0
   tablet2(i)%deltax      = -1.0D0
   tablet2(i)%maxdatval   = -1.0D0

   tablet2(i)%z   = -1.0D9 !only used in slices_t2

   tablet2(i)%needs_M2_gt_2M1     = .False.
  enddo
 
 end subroutine initializetables_type2_blank

 !*********************************************************** 
 subroutine initializetables2(S95_t2)
 !***********************************************************  
 ! fills S95_t2
 !***********************************************************  
  use store_pathname      
  use usefulbits, only: Hneut,Hplus,Chineut,Chiplus,small,file_id_common2,not_a_particle
  implicit none

  !--------------------------------------input
  type(table2) :: S95_t2(:)
  !-----------------------------------internal
  integer :: i,tno,j,x,xbeg,xend,k,ios
  character(LEN=2) :: tableno            
  character(len=100),allocatable :: filename(:) 
  double precision :: dummy
  double precision, allocatable :: testrow(:)
  integer :: file_id_arr(2)
  double precision :: maxdatval
  !-------------------------------------------       
  file_id_arr(1)=file_id_2_exp
  file_id_arr(2)=file_id_2_obs

  xbeg=lbound(S95_t2,dim=1)
  xend=ubound(S95_t2,dim=1)  
 
  allocate(filename(xbeg:xend))
  x=xbeg-1
   
  tno=14
  do i=1,8
   x=x+1
   tno=tno+1
   if((x.eq.3).or.(x.eq.7))tno=tno+1 
   write(tableno,'(I2)')tno      
   
   S95_t2(x)%id=tno*10
   S95_t2(x)%expt='LEP'
   S95_t2(x)%energy=0.208D0        
   S95_t2(x)%deltax=0.0D0 
   S95_t2(x)%particle_x1=Hneut
   S95_t2(x)%particle_x2=Hneut

   select case(S95_t2(x)%id)
   case(220,230,240)
    S95_t2(x)%label='hep-ex/0602042 (LEP)'
   case default
    S95_t2(x)%label='hep-ex/0602042, table '//tableno//' (LEP)'
   end select

   S95_t2(x)%sep1=1.0D0
   S95_t2(x)%sep2=1.0D0
   S95_t2(x)%maxdatval   = 1.0D2   
   !S95_t2(x)%OBid=x+2
       
   select case(S95_t2(x)%id)
   case(150,160,220)
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=60.0D0
    S95_t2(x)%xmin2=2.0D0
    S95_t2(x)%xmax2=120.0D0
    S95_t2(x)%needs_M2_gt_2M1=.True.
   case(180,190,230,240)   
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=180.0D0
    S95_t2(x)%xmin2=1.0D0
    S95_t2(x)%xmax2=180.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.False.
   case(200,210) 
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=90.0D0
    S95_t2(x)%xmin2=2.0D0
    S95_t2(x)%xmax2=180.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.True.                           
   case default
   write(*,*)'error in initializetables2 (a)' 
    stop
   end select
   
   filename(x)='table'//tableno//'full' 
  enddo
   
  do i=5,10
   x=x+1
   tno=i
   write(tableno,'(I2)')tno      
   
   S95_t2(x)%id=900+tno
   S95_t2(x)%expt='LEP' 
   S95_t2(x)%energy=0.208D0     
   S95_t2(x)%deltax=0.0D0    
   S95_t2(x)%label='hep-ex/0401026, fig '//trim(adjustl(tableno))//' (OPAL)'
   S95_t2(x)%sep1=1.0D0
   S95_t2(x)%sep2=1.0D0
   S95_t2(x)%maxdatval   = 1.0D6 !these tables are in fb   
       
   select case(tno)
   case(5,6,7,8)
    S95_t2(x)%xmin1=0.0D0
    S95_t2(x)%xmax1=100.0D0
    S95_t2(x)%xmin2=75.0D0
    S95_t2(x)%xmax2=120.0D0
    S95_t2(x)%needs_M2_gt_2M1=.False. 
    S95_t2(x)%particle_x1=Chineut
    S95_t2(x)%particle_x2=Chiplus
   case(9,10)   
    S95_t2(x)%xmin1=0.0D0
    S95_t2(x)%xmax1=100.0D0
    S95_t2(x)%xmin2=50.0D0
    S95_t2(x)%xmax2=200.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.False. 
    S95_t2(x)%particle_x1=Chineut
    S95_t2(x)%particle_x2=Chineut                    
   case default
   write(*,*)'error in initializetables2 (b)' 
    stop
   end select
   
   filename(x)='1026_fig'//trim(adjustl(tableno)) 
  enddo  

  x=x+1    
  S95_t2(x)%id=6065
  S95_t2(x)%expt='LEP'
  S95_t2(x)%energy=0.208D0        
  S95_t2(x)%deltax=0.0D0 
  S95_t2(x)%label='[hep-ex] arXiv:1301.6065 (LEP)'
  S95_t2(x)%sep1=0.1D0
  S95_t2(x)%sep2=0.5D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.0D0
  S95_t2(x)%xmax1=1.0D0
  S95_t2(x)%xmin2=43.0D0
  S95_t2(x)%xmax2=95.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hplus
  filename(x)="6065_LEP_HpHm_fig4" 

  x=x+1    
  S95_t2(x)%id=02671
  S95_t2(x)%expt='LEP'
  S95_t2(x)%energy=0.208D0        
  S95_t2(x)%deltax=0.0D0 
  S95_t2(x)%label='[hep-ex] arXiv:0812.0267, Fig.10a (OPAL,LEP)'
  S95_t2(x)%sep1=0.5D0
  S95_t2(x)%sep2=0.5D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=12.0D0
  S95_t2(x)%xmax1=90.0D0
  S95_t2(x)%xmin2=40.0D0
  S95_t2(x)%xmax2=93.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hplus
  filename(x)="OPAL_H+H-_AWAW_0267" 

  x=x+1    
  S95_t2(x)%id=02672
  S95_t2(x)%expt='LEP'
  S95_t2(x)%energy=0.208D0        
  S95_t2(x)%deltax=0.0D0 
  S95_t2(x)%label='[hep-ex] arXiv:0812.0267, Fig.10b (OPAL,LEP)'
  S95_t2(x)%sep1=0.5D0
  S95_t2(x)%sep2=0.5D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=12.0D0
  S95_t2(x)%xmax1=77.0D0
  S95_t2(x)%xmin2=40.0D0
  S95_t2(x)%xmax2=80.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hplus
  filename(x)="OPAL_H+H-_AWtaunu_0267" 


  x=x+1    
  S95_t2(x)%id=3381
  S95_t2(x)%expt=' D0'    
  S95_t2(x)%energy=1.96D0  
  S95_t2(x)%deltax=0.0D0   
  S95_t2(x)%label='[hep-ex] arXiv:0905.3381, table I (D0)'
  S95_t2(x)%sep1=0.1D0
  S95_t2(x)%sep2=5.0D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.2D0
  S95_t2(x)%xmax1=3.0D0
  S95_t2(x)%xmin2=80.0D0
  S95_t2(x)%xmax2=200.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)="D0_h-aa-mumumumu_3381" 

  x=x+1    
  S95_t2(x)%id=3382
  S95_t2(x)%expt=' D0'   
  S95_t2(x)%energy=1.96D0    
  S95_t2(x)%deltax=0.0D0      
  S95_t2(x)%label='[hep-ex] arXiv:0905.3381, table II (D0)'
  S95_t2(x)%sep1=0.2D0
  S95_t2(x)%sep2=5.0D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=3.6D0
  S95_t2(x)%xmax1=19.0D0
  S95_t2(x)%xmin2=85.0D0
  S95_t2(x)%xmax2=200.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)="D0_h-aa-tautaumumu_3381"

  x=x+1  
  S95_t2(x)%id=6227
  S95_t2(x)%expt=' D0'
  S95_t2(x)%label='D0 Note 6227'
  S95_t2(x)%energy=1.96D0    
  S95_t2(x)%sep1=0.04D0 
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.06D0
  S95_t2(x)%xmax1=0.18D0
  S95_t2(x)%xmin2=90.0D0
  S95_t2(x)%xmax2=300.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='D0_h-bb_h-tautau_comb_5.2-7.3fb_6227'  

  x=x+1  
  S95_t2(x)%id=5053
  S95_t2(x)%expt='ATL'
  S95_t2(x)%label='[hep-ex] arXiv:1406.5053 (ATLAS)'
  S95_t2(x)%energy=8.0D0    
  S95_t2(x)%sep1=5.0D0 
  S95_t2(x)%sep2=5.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=120.0D0
  S95_t2(x)%xmax1=130.0D0
  S95_t2(x)%xmin2=260.0D0
  S95_t2(x)%xmax2=500.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='5053_H-hh-gagabb_20fb-1'  

!   x=x+1  
!   S95_t2(x)%id=13032
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='CMS-PAS-HIG-13-032'
!   S95_t2(x)%energy=8.0D0    
!   S95_t2(x)%sep1=5.0D0 
!   S95_t2(x)%sep2=1.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=120.0D0
!   S95_t2(x)%xmax1=130.0D0
!   S95_t2(x)%xmin2=260.0D0
!   S95_t2(x)%xmax2=1100.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='13032_H-hh-gagabb_19.7fb-1'  

  x=x+1  
  S95_t2(x)%id=06896
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='[hep-ex] arXiv:1603.06896 (CMS)'
  S95_t2(x)%energy=8.0D0    
  S95_t2(x)%sep1=5.0D0 
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=120.0D0
  S95_t2(x)%xmax1=130.0D0
  S95_t2(x)%xmin2=260.0D0
  S95_t2(x)%xmax2=1100.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='06896_H-hh-gagabb_19.7fb-1'  

!   x=x+1  
!   S95_t2(x)%id=011811
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='[hep-ex] arXiv:1510.01181 (CMS)'
!   S95_t2(x)%energy=8.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=260.0D0
!   S95_t2(x)%xmax2=350.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='140341_CMS_H-hh-bbtautau_19.7fb-1'  

!   x=x+1  
!   S95_t2(x)%id=044781
!   S95_t2(x)%expt='ATLAS'
!   S95_t2(x)%label='[hep-ex] arXiv:1502.04478 (ATLAS)'
!   S95_t2(x)%energy=8.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=220.0D0
!   S95_t2(x)%xmax2=1000.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='044781_ATLAS_gg-A-hZ-tautaull_20.3fb-1'
! 
!   x=x+1  
!   S95_t2(x)%id=044782
!   S95_t2(x)%expt='ATLAS'
!   S95_t2(x)%label='[hep-ex] arXiv:1502.04478 (ATLAS)'
!   S95_t2(x)%energy=8.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=220.0D0
!   S95_t2(x)%xmax2=1000.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='044782_ATLAS_gg-A-hZ-bbll_20.3fb-1'
  
!   x=x+1  
!   S95_t2(x)%id=011812
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='[hep-ex] arXiv:1510.01181 (CMS)'
!   S95_t2(x)%energy=8.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=220.0D0
!   S95_t2(x)%xmax2=350.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='140342_CMS_A-hZ-tautaull_19.7fb-1'  

!   x=x+1  
!   S95_t2(x)%id=20160151
!   S95_t2(x)%expt='ATLAS'
!   S95_t2(x)%label='ATLAS-CONF-2016-015'
!   S95_t2(x)%energy=13.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=220.0D0
!   S95_t2(x)%xmax2=2000.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='20160151_ATLAS_gg-A-hZ-bbll_3.2fb-1'  
!   
!   x=x+1  
!   S95_t2(x)%id=20160152
!   S95_t2(x)%expt='ATLAS'
!   S95_t2(x)%label='ATLAS-CONF-2016-015'
!   S95_t2(x)%energy=13.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=220.0D0
!   S95_t2(x)%xmax2=2000.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='20160152_ATLAS_bb-A-hZ-bbll_3.2fb-1'   
!   
!   x=x+1  
!   S95_t2(x)%id=16002
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='CMS-PAS-HIG-16-002'
!   S95_t2(x)%energy=13.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=5.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=260.0D0
!   S95_t2(x)%xmax2=1200.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='16002_CMS_H-hh-bbbb_2.3fb-1'  
! 
!   x=x+1  
!   S95_t2(x)%id=16029
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='CMS-PAS-HIG-16-029'
!   S95_t2(x)%energy=13.0D0
!   S95_t2(x)%sep1=10.0D0 
!   S95_t2(x)%sep2=10.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=115.0D0
!   S95_t2(x)%xmax1=135.0D0
!   S95_t2(x)%xmin2=250.0D0
!   S95_t2(x)%xmax2=900.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='16029_CMS_H-hh-bbtautau_12.9fb-1'  
! 
!   x=x+1  
!   S95_t2(x)%id=14013
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='CMS-PAS-HIG-14-013,arXiv:1503.04114 (CMS)'
!   S95_t2(x)%energy=8.0D0    
!   S95_t2(x)%sep1=5.0D0 
!   S95_t2(x)%sep2=1.0D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=120.0D0
!   S95_t2(x)%xmax1=130.0D0
!   S95_t2(x)%xmin2=270.0D0
!   S95_t2(x)%xmax2=1097.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.True. 
!   S95_t2(x)%particle_x1=Hneut
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='14013_H-hh-bbbb_17.9fb-1'  

  x=x+1  
  S95_t2(x)%id=16032
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-16-032'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=2.5D0 
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=122.5D0
  S95_t2(x)%xmax1=127.5D0
  S95_t2(x)%xmin2=250.0D0
  S95_t2(x)%xmax2=900.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='16032_CMS_H-hh-gagabb_2.70fb-1'  

  x=x+1  
  S95_t2(x)%id=14022
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-14-022'
  S95_t2(x)%energy=8.0D0
  S95_t2(x)%sep1=1.0D0 
  S95_t2(x)%sep2=2.5D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=5.0D0
  S95_t2(x)%xmax1=15.0D0
  S95_t2(x)%xmin2=122.5D0
  S95_t2(x)%xmax2=127.5D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='14022_H-hh-tautau_19.7fb-1'  

  x=x+1  
  S95_t2(x)%id=150600424
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='arXiv:1506.00424 (CMS)'
  S95_t2(x)%energy=8.0D0
  S95_t2(x)%sep1=0.05D0 
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.25D0
  S95_t2(x)%xmax1=3.55D0
  S95_t2(x)%xmin2=86.0D0
  S95_t2(x)%xmax2=150.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='150600424_CMS_H-aa-mumu'  

  x=x+1  
  S95_t2(x)%id=16035
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-16-035'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=0.05D0 
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.25D0
  S95_t2(x)%xmax1=3.55D0
  S95_t2(x)%xmin2=86.0D0
  S95_t2(x)%xmax2=150.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='16035_CMS_H-aa-mumumu_2.8fb-1'

  x=x+1  
  S95_t2(x)%id=02301
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='[hep-ex] arXiv:1506.02301 (CMS)'
  S95_t2(x)%energy=8.0D0    
  S95_t2(x)%sep1=0.01D0 
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.00D0
  S95_t2(x)%xmax1=0.10D0
  S95_t2(x)%xmin2=150.0D0
  S95_t2(x)%xmax2=840.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='14006_CMS_h-gaga_2D'  

  x=x+1  
  S95_t2(x)%id=003892
  S95_t2(x)%expt='ATL'
  S95_t2(x)%label='[hep-ex] arXiv:1509.00389 (ATLAS)'
  S95_t2(x)%energy=8.0D0    
  S95_t2(x)%sep1=0.05D0
  S95_t2(x)%sep2=100.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.20D0
  S95_t2(x)%xmax1=0.80D0
  S95_t2(x)%xmin2=300.0D0
  S95_t2(x)%xmax2=1000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='00389_H-WW_20.3fb-1_2D'  

  x=x+1  
  S95_t2(x)%id=20160791
  S95_t2(x)%expt='ATL'
  S95_t2(x)%label='ATLAS-CONF-2016-079'
  S95_t2(x)%energy=13.0D0    
  S95_t2(x)%sep1=0.01D0
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.0D0
  S95_t2(x)%xmax1=0.10D0
  S95_t2(x)%xmin2=400.0D0
  S95_t2(x)%xmax2=1000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='20160791_gg-H-ZZ-4l_14.8fb-1' 

  x=x+1  
  S95_t2(x)%id=01123
  S95_t2(x)%expt='ATL'
  S95_t2(x)%label='arXiv:1710.01123 (ATLAS)'
  S95_t2(x)%energy=13.0D0    
  S95_t2(x)%sep1=5.0D0
  S95_t2(x)%sep2=100.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=5.0D0
  S95_t2(x)%xmax1=15.0D0
  S95_t2(x)%xmin2=200.0D0
  S95_t2(x)%xmax2=4000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='01123_ATLAS_gg-H-WW-e_nu_mu_nu_36.1fb-1' 

!   x=x+1  
!   S95_t2(x)%id=160331
!   S95_t2(x)%expt='CMS'
!   S95_t2(x)%label='CMS-PAS-HIG 16-033'
!   S95_t2(x)%energy=13.0D0    
!   S95_t2(x)%sep1=2.0D0
!   S95_t2(x)%sep2=0.1D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=0.0D0
!   S95_t2(x)%xmax1=40.0D0
!   S95_t2(x)%xmin2=130.0D0
!   S95_t2(x)%xmax2=2520.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.False. 
!   S95_t2(x)%particle_x1=not_a_particle
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='160331_CMS_H-ZZ-4l_ggF_12.9fb-1'  

  x=x+1  
  S95_t2(x)%id=06386
  S95_t2(x)%expt='ATLAS'
  S95_t2(x)%label='arXiv:1712.06386 (ATLAS)'
  S95_t2(x)%energy=13.0D0    
  S95_t2(x)%sep1=0.01D0
  S95_t2(x)%sep2=1.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.0D0
  S95_t2(x)%xmax1=.10D0
  S95_t2(x)%xmin2=400.0D0
  S95_t2(x)%xmax2=1000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='06386_ATLAS_gg-H-ZZ-4l+2l2nu_36.1fb-1'  


  x=x+1  
  S95_t2(x)%id=170121
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG 17-012'
  S95_t2(x)%energy=13.0D0    
  S95_t2(x)%sep1=10.0D0
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.0D0
  S95_t2(x)%xmax1=100.0D0
  S95_t2(x)%xmin2=130.0D0
  S95_t2(x)%xmax2=3000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='17012_CMS_gg-H-ZZ_35.9fb-1'

  x=x+1  
  S95_t2(x)%id=160332
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG 16-033'
  S95_t2(x)%energy=13.0D0    
  S95_t2(x)%sep1=2.0D0
  S95_t2(x)%sep2=0.1D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.0D0
  S95_t2(x)%xmax1=40.0D0
  S95_t2(x)%xmin2=130.0D0
  S95_t2(x)%xmax2=2520.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='160332_CMS_H-ZZ-4l_VBF_12.9fb-1'  

  x=x+1  
  S95_t2(x)%id=150011
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-15-001,arXiv:1603.02991 (CMS)'
  S95_t2(x)%energy=8.0D0
  S95_t2(x)%sep1=5.0D0 
  S95_t2(x)%sep2=5.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=15.0D0
  S95_t2(x)%xmax1=1000.0D0
  S95_t2(x)%xmin2=15.0D0
  S95_t2(x)%xmax2=1000.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='150011_H2-ZH1-lltautau_19.8fb-1'  

  x=x+1  
  S95_t2(x)%id=150012
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-15-001,arXiv:1603.02991 (CMS)'
  S95_t2(x)%energy=8.0D0
  S95_t2(x)%sep1=2.4D0 
  S95_t2(x)%sep2=2.4D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=24.0D0
  S95_t2(x)%xmax1=1200.0D0
  S95_t2(x)%xmin2=24.0D0
  S95_t2(x)%xmax2=1200.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='150012_H2-ZH1-llbb_19.8fb-1'  

  x=x+1  
  S95_t2(x)%id=16010
  S95_t2(x)%expt='CMS'
  S95_t2(x)%label='CMS-PAS-HIG-16-010'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=5.0D0 
  S95_t2(x)%sep2=5.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=50.0D0
  S95_t2(x)%xmax1=600.0D0
  S95_t2(x)%xmin2=300.0D0
  S95_t2(x)%xmax2=700.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='16010_H2-ZH1-llbb_2.3fb-1'

!   x=x+1  
!   S95_t2(x)%id=2016059
!   S95_t2(x)%expt='ATLAS'
!   S95_t2(x)%label='ATLAS-CONF-2016-059'
!   S95_t2(x)%energy=13.0D0
!   S95_t2(x)%sep1=0.02D0 
!   S95_t2(x)%sep2=0.5D0 
!   S95_t2(x)%maxdatval=1.0D6
!   S95_t2(x)%xmin1=0.00D0
!   S95_t2(x)%xmax1=0.10D0
!   S95_t2(x)%xmin2=200.0D0
!   S95_t2(x)%xmax2=2400.0D0
!   S95_t2(x)%needs_M2_gt_2M1=.False. 
!   S95_t2(x)%particle_x1=not_a_particle
!   S95_t2(x)%particle_x2=Hneut
!   filename(x)='2016059_Atlas_pp-H-gaga_15.4fb-1'

  x=x+1  
  S95_t2(x)%id=4147
  S95_t2(x)%expt='ATLAS'
  S95_t2(x)%label='[hep-ex] arXiv:1707.04147 (ATLAS)'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=0.02D0 
  S95_t2(x)%sep2=0.5D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.00D0
  S95_t2(x)%xmax1=0.10D0
  S95_t2(x)%xmin2=200.0D0
  S95_t2(x)%xmax2=2700.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='4147_Atlas_pp-H-gaga_36.7fb-1'

! By DDD, Nov 8 2018
  x=x+1  
  S95_t2(x)%id=20160341
  S95_t2(x)%expt='ATLAS'
  S95_t2(x)%label='ATLAS-EXOT-2016-034, bbh'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=10.0D0 
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=130.0D0
  S95_t2(x)%xmax1=700.0D0
  S95_t2(x)%xmin2=230.0D0
  S95_t2(x)%xmax2=800.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='20160341_ATLAS_A_HZ_llbb_36.1fb-1_bbh'

  
! By DDD, Nov 8 2018
  x=x+1  
  S95_t2(x)%id=20160342
  S95_t2(x)%expt='ATLAS'
  S95_t2(x)%label='ATLAS-EXOT-2016-034, ggf'
  S95_t2(x)%energy=13.0D0
  S95_t2(x)%sep1=10.0D0 
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=130.0D0
  S95_t2(x)%xmax1=700.0D0
  S95_t2(x)%xmin2=230.0D0
  S95_t2(x)%xmax2=800.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)='20160342_ATLAS_A_HZ_llbb_36.1fb-1_ggf'  

  ! checks we've filled the whole array  
  if(x.ne.xend)then
   write(*,*)'error in initializetables2 (c)',x,xend
   stop
  endif  
  
  ! read in the tables
  do x=xbeg,xend
   
   S95_t2(x)%nx2=nint((S95_t2(x)%xmax2-S95_t2(x)%xmin2)/S95_t2(x)%sep2)+1
   S95_t2(x)%nx1=nint((S95_t2(x)%xmax1-S95_t2(x)%xmin1)/S95_t2(x)%sep1)+1             

   allocate(S95_t2(x)%dat(S95_t2(x)%nx2,S95_t2(x)%nx1,2)) 
  enddo

  ! read in the tables
  open(file_id_common2,file = trim(adjustl(pathname))//'Expt_tables/' // &
      &  'S95_t2.binary',form='unformatted')

  read(file_id_common2,iostat=ios)S95_t2(xbeg)%dat
  if(ios.eq.0)then
    do x=xbeg+1,xend
     read(file_id_common2)S95_t2(x)%dat
    enddo

  else

    do x=xbeg,xend
     open(file_id_2_exp,file=trim(adjustl(pathname))//('Expt_tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                //trim(adjustl(filename(x)))//'_pred.txt'))            
     open(file_id_2_obs,file=trim(adjustl(pathname))//('Expt_tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                //trim(adjustl(filename(x)))//'_obs.txt'))                                 

     ! fill S95 from file
     ! row 0 and column 0 in LEP file contain higgs masses 
     ! and (0,0) ie top left set to -100
     ! so avoid them         
     allocate(testrow(0:S95_t2(x)%nx1))

     do k=lbound(file_id_arr,dim=1),ubound(file_id_arr,dim=1)
       read(file_id_arr(k),*)( testrow(i), i=0,S95_t2(x)%nx1 )
       if((testrow(0)+100.0D0).gt.small)stop 'error in initializetables2 (d)' !top left number should be -100
       do i=1,S95_t2(x)%nx1      
         if( abs(testrow(i)- (S95_t2(x)%xmin1 + dble(i-1)*S95_t2(x)%sep1) ).gt.small*S95_t2(x)%sep1 )then
            write(*,*)S95_t2(x)%id,testrow(i),(S95_t2(x)%xmin1 + dble(i-1)*S95_t2(x)%sep1)
            stop 'error in initializetables2 (e)'  
         endif
       enddo
     enddo
     deallocate(testrow)
  
     do j=1,S95_t2(x)%nx2  
      read(file_id_2_exp,*) dummy, ( S95_t2(x)%dat(j,i,2), i=1,S95_t2(x)%nx1 )
      if( abs(dummy- (S95_t2(x)%xmin2 + dble(j-1)*S95_t2(x)%sep2) ).gt.small*S95_t2(x)%sep2 ) then
!       write(*,*) S95_t2(x)%nx1, dummy, S95_t2(x)%dat(j,:,2)
!       write(*,*) j,S95_t2(x)%nx2  ,S95_t2(x)%xmin2 + dble(j-1)*S95_t2(x)%sep2   
       stop 'error in initializetables2 (f)'
      endif
      read(file_id_2_obs,*) dummy, ( S95_t2(x)%dat(j,i,1), i=1,S95_t2(x)%nx1 )   
      if( abs(dummy- (S95_t2(x)%xmin2 + dble(j-1)*S95_t2(x)%sep2) ).gt.small*S95_t2(x)%sep2 ) then
       write(*,*) "Problematic analysis: ",S95_t2(x)%id
       stop 'error in initializetables2 (g)'                 
      endif
     end do
   
     maxdatval=S95_t2(x)%maxdatval
     if( maxdatval .gt. 0.0D0 )then
       ! set entries .ge. S95_t2(x)%maxdatval to (-4): they will not be relevent
       where(  S95_t2(x)%dat .ge. maxdatval )  S95_t2(x)%dat= - 4.0D0                        
     endif

     close(file_id_2_exp)
     close(file_id_2_obs)       
    enddo   

    rewind(file_id_common2)
#ifndef WEBVERSION
    do x=xbeg,xend
     write(file_id_common2)S95_t2(x)%dat
    enddo
#endif

  endif

  close(file_id_common2)
  deallocate(filename)      
                    
 end subroutine initializetables2
 !***********************************************************
 function t2elementnumberfromid(t2,id)
  !--------------------------------------input
  type(table2), intent(in) :: t2(:)
  integer, intent(in) :: id
  !-----------------------------------function
  integer :: t2elementnumberfromid
  !-----------------------------------internal
  integer :: n,x
  !-------------------------------------------

  n=0
  do x=lbound(t2,dim=1),ubound(t2,dim=1)
   if(t2(x)%id.eq.id)then
    n=n+1
    t2elementnumberfromid=x
   endif
  enddo

  if(n.ne.1)stop 'problem in function t2elementnumberfromid 1'

 end function t2elementnumberfromid
 !*********************************************************** 
 subroutine fill_slices_t1_from_slices_of_t2(t2,v1orv2,xy_selection,ftype_selection,slices_t1)
 ! if this subroutine is used,  
 ! don't forget to deallocate slices_t1(x)%dat at some point
 !*********************************************************** 
 implicit none 
  !--------------------------------------input
  type(table2), intent(in) :: t2
  integer, intent(in) :: v1orv2
  integer, intent(in) :: xy_selection(:)
  integer, intent(in) :: ftype_selection(:)
  !-------------------------------------output
  type(table1) :: slices_t1(:)  !i.e. 2 slices
  !-----------------------------------internal
  integer :: i,j,k,n
  integer :: n_ftype_selection
  !-------------------------------------------
  n_ftype_selection=ubound(ftype_selection,dim=1)

  do n=lbound(ftype_selection,dim=1),n_ftype_selection
     if(ftype_selection(n).lt.lbound(t2%dat,dim=3))stop 'problem in fill_slices_t1_from_slices_of_t2 3a'
     if(ftype_selection(n).gt.ubound(t2%dat,dim=3))stop 'problem in fill_slices_t1_from_slices_of_t2 3b'
  enddo

  if(lbound(xy_selection,dim=1).ne.lbound(slices_t1,dim=1))then
     stop 'problem in fill_slices_t1_from_slices_of_t2 1a'
  endif
  if(ubound(xy_selection,dim=1).ne.ubound(slices_t1,dim=1))then
     stop 'problem in fill_slices_t1_from_slices_of_t2 1b'
  endif

  select case(v1orv2)
  case(1)

    do n=lbound(slices_t1,dim=1),ubound(slices_t1,dim=1)

     if(xy_selection(n).lt.lbound(t2%dat,dim=1))stop 'problem in fill_slices_t1_from_slices_of_t2 4a'
     if(xy_selection(n).gt.ubound(t2%dat,dim=1))stop 'problem in fill_slices_t1_from_slices_of_t2 4b'

     slices_t1(n)%id          =  t2%id    
     slices_t1(n)%nx          =  t2%nx1  
     slices_t1(n)%xmax        =  t2%xmax1 
     slices_t1(n)%xmin        =  t2%xmin1
     slices_t1(n)%sep         =  t2%sep1 
     slices_t1(n)%deltax      =  t2%deltax
  
     allocate( slices_t1(n)%dat(slices_t1(n)%nx,n_ftype_selection) ) 
     slices_t1(n)%dat = -1.0D0

     do i=1,slices_t1(n)%nx
       do k=1,n_ftype_selection
        slices_t1(n)%dat(i,k)=t2%dat(xy_selection(n),i,ftype_selection(k))
       enddo
     enddo

    enddo
  case(2)

    do n=lbound(slices_t1,dim=1),ubound(slices_t1,dim=1)

     if(xy_selection(n).lt.lbound(t2%dat,dim=2))stop 'problem in fill_slices_t1_from_slices_of_t2 4aa'
     if(xy_selection(n).gt.ubound(t2%dat,dim=2))stop 'problem in fill_slices_t1_from_slices_of_t2 4bb'

     slices_t1(n)%id          =  t2%id    
     slices_t1(n)%nx          =  t2%nx2  
     slices_t1(n)%xmax        =  t2%xmax2 
     slices_t1(n)%xmin        =  t2%xmin2
     slices_t1(n)%sep         =  t2%sep2 
     slices_t1(n)%deltax      =  t2%deltax
  
     allocate( slices_t1(n)%dat(slices_t1(n)%nx,n_ftype_selection) ) 
     slices_t1(n)%dat = -1.0D0

     do j=1,slices_t1(n)%nx
       do k=1,n_ftype_selection
        slices_t1(n)%dat(j,k)=t2%dat(j,xy_selection(n),ftype_selection(k))
       enddo
     enddo

    enddo
  case default
   stop 'problem in fill_slices_t1_from_slices_of_t2 5'
  end select

 end subroutine fill_slices_t1_from_slices_of_t2
 !***********************************************************

 !*********************************************************** 
 subroutine fill_t1_from_t2(t2,v1orv2,xy_selection,ftype_selection,t1)
 ! if this subroutine is used,  
 ! don't forget to deallocate slices_t1(x)%dat at some point
 !*********************************************************** 
 implicit none 
  !--------------------------------------input
  type(table2), intent(in) :: t2
  integer, intent(in) :: v1orv2
  integer, intent(in) :: xy_selection
  integer, intent(in) :: ftype_selection(:)
  !-------------------------------------output
  type(table1) :: t1  
  !-----------------------------------internal
  integer :: i,j,k,n
  integer :: n_ftype_selection
  !-------------------------------------------
  n_ftype_selection=ubound(ftype_selection,dim=1)

  do n=lbound(ftype_selection,dim=1),n_ftype_selection
     if(ftype_selection(n).lt.lbound(t2%dat,dim=3))stop 'problem in fill_t1_from_t2 3a'
     if(ftype_selection(n).gt.ubound(t2%dat,dim=3))stop 'problem in fill_t1_from_t2 3b'
  enddo

  t1%id          =  t2%id 
  t1%deltax      =  t2%deltax

  select case(v1orv2)
  case(1)

     if(xy_selection.lt.lbound(t2%dat,dim=1))stop 'problem in fill_t1_from_t2 4a'
     if(xy_selection.gt.ubound(t2%dat,dim=1))stop 'problem in fill_t1_from_t2 4b'
  
     t1%nx          =  t2%nx1  
     t1%xmax        =  t2%xmax1 
     t1%xmin        =  t2%xmin1
     t1%sep         =  t2%sep1 
  
     allocate( t1%dat(t1%nx,n_ftype_selection) ) 
     t1%dat = -1.0D0

     do i=1,t1%nx
       do k=1,n_ftype_selection
        t1%dat(i,k)=t2%dat(xy_selection,i,ftype_selection(k))
       enddo
     enddo

  case(2)

     if(xy_selection.lt.lbound(t2%dat,dim=2))stop 'problem in fill_t1_from_t2 4aa'
     if(xy_selection.gt.ubound(t2%dat,dim=2))stop 'problem in fill_t1_from_t2 4bb'

     t1%nx          =  t2%nx2  
     t1%xmax        =  t2%xmax2 
     t1%xmin        =  t2%xmin2
     t1%sep         =  t2%sep2 
  
     allocate( t1%dat(t1%nx,n_ftype_selection) ) 
     t1%dat = -1.0D0

     do j=1,t1%nx
       do k=1,n_ftype_selection
        t1%dat(j,k)=t2%dat(j,xy_selection,ftype_selection(k))
       enddo
     enddo

  case default
   stop 'problem in fill_t1_from_t2 5'
  end select

 end subroutine fill_t1_from_t2
 !*********************************************************** 
end module S95tables_type2
!************************************************************
