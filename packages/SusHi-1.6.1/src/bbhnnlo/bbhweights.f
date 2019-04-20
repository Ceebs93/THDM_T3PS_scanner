C-{{{ subroutine bbhweights

      subroutine bbhweights(modin,thdmv,Hmx,Amx,lpseudo
     & ,Sindin,tb,muSUSY,lam,GF,deltamb,weH)
!--------------------------------------------------------
!This routine gives the weights between SM and MSSM
!cross sections for bbh in accordance to hep-ph/0305101
!(Guasch, Haefliger, Spira), formula 23.
!February 2012 - S. Liebler
!--------------------------------------------------------
      implicit none
      double precision weH(2)
      double precision deltamb, Hmx(3,3), Amx(3,3), tb
      double precision delmb_in, Sqrt2, vev, vevS, GF, lam, muSUSY
      integer i1, modin, thdmv, lpseudo, Sindin

      Sqrt2 = dsqrt(2.d0)
      vev = 1.d0/dsqrt(Sqrt2*GF)
      vevS = muSUSY*Sqrt2/lam
      
      delmb_in = deltamb

!MSSM and NMSSM
      if ((modin.eq.1).or.(modin.eq.3)) then
       if (lpseudo.eq.1) then
         weH(1) = Amx(Sindin,2)*tb
         weH(2) = weH(1) * (1.d0/(1.d0+delmb_in))
     &        * (1.d0 - (1.d0/tb**2)*delmb_in
     &           - Amx(Sindin,3)*vev/(Amx(Sindin,2)*vevS*tb)*delmb_in)   
       else
         weH(1) = Hmx(Sindin,1)/Cos(Atan(tb))
         weH(2) = weH(1) * (1.d0/(1.d0+delmb_in))
     &     * (1.d0 + delmb_in * Hmx(Sindin,2)/(tb*Hmx(Sindin,1))
     & + delmb_in*Hmx(Sindin,3)*vev*Cos(Atan(tb))/(Hmx(Sindin,1)*vevS))
       endif
      endif

!2HDM
      if (modin.eq.2) then
         if ((thdmv.eq.1).or.(thdmv.eq.4)) then
          if (lpseudo.eq.1) then
            weH(1) = -1.d0/tb
          else
            weH(1) = Hmx(Sindin,2)/Sin(Atan(tb))
          endif
         else if ((thdmv.eq.2).or.(thdmv.eq.3)) then
          if (lpseudo.eq.1) then
            weH(1) = tb
          else
            weH(1) = Hmx(Sindin,1)/Cos(Atan(tb))
          endif
         else
            write(*,*) "2HDM version unknown!"
            stop
         endif
         weH(2) = weH(1)
      endif

      do i1=1,2
         weH(i1) = weH(i1) * weH(i1)
      end do

      end

C-}}}
