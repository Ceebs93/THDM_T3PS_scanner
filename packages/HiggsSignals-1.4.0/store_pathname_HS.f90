!******************************************************************
module store_pathname_HS
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 59
 character(len=pathname_length),parameter :: pathname_HS= &
     &     "/home/de3u14/lib/build/hep/HiggsSignals/HiggsSignals-1.4.0" // &
     &     "/"

end module store_pathname_HS
!******************************************************************
