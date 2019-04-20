!******************************************************************
module store_pathname
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 79
 character(len=pathname_length),parameter :: pathname= &
     &     "/scratch/de3u14/2HDM/Magellan_THDM_T3PS_scanner/packages/HiggsBounds-5.3.2beta" // &
     &     "/"

end module store_pathname
!******************************************************************
