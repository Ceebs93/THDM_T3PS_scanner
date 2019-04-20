!******************************************************************
module store_pathname
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 57
 character(len=pathname_length),parameter :: pathname= &
     &     "/home/de3u14/lib/build/hep/HiggsBounds/HiggsBounds-4.3.1" // &
     &     "/"

end module store_pathname
!******************************************************************
