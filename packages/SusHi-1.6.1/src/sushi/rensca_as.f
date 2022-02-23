! This file is part of SusHi.
! 
! It includes the routine for an analytic calculation of the muR dependence
! 
      subroutine rensca_as(api0,mu0,mu,nf,uu,nord,xsec,sigmu)

      implicit none
      integer i,j,k
      integer uu,nord,nas
      real*8 xsec(0:3),xscoef(0:3),cc(0:4,0:4),lm0,pi,api0,api,mu0,mu,nf
      real*8 sigmu(0:3)
      real*8 beta0,beta1,beta2,beta3
      common/bfunc/beta0,beta1,beta2,beta3
      include '../commons/common-keys.f'
      
      pi=4*datan(1.d0)

      xscoef(0) = xsec(0)
      xscoef(1) = (xsec(1)- xsec(0))/api0
      xscoef(2) = (xsec(2) - xsec(1))/api0**2
      xscoef(3) = (xsec(3) - xsec(2))/api0**3

      nas=nord
      if (nord.eq.3) nas=runasn3lo-1
      call runalpha(api0,mu0,mu,nf,nas+1,0,api)

      lm0 = 2*dlog(mu/mu0)
      cc=0.d0

      do i=0,nord
         cc(i,0) = (api/api0)**uu*xscoef(i)
      enddo

      cc(1,1) = uu*beta0*cc(0,0)
      cc(2,2) = (1+uu)*beta0*cc(1,1)/2.d0
      cc(2,1) = uu*beta1*cc(0,0) + (1+uu)*beta0*cc(1,0)
      cc(3,3) = ((2+uu)*beta0*cc(2,2))/3.d0
      cc(3,2) = ((1+uu)*beta1*cc(1,1) + (2+uu)*beta0*cc(2,1))/2.d0
      cc(3,1) = uu*beta2*cc(0,0) + (1+uu)*beta1*cc(1,0) + (2+uu)*beta0
     &     *cc(2,0)
!      cc(4,4) =  (beta0*(3 + uu)*cc(3,3))/4
!      cc(4,3) =  (beta1*(2 + uu)*cc(2,2) + beta0*(3 + uu)*cc(3,2))
!     &     /3.d0
!      cc(4,2) = (beta2*(1 + uu)*cc(1,1) + beta1*(2 + uu)*cc(2,1) +
!     &     beta0*(3 + uu)*cc(3,1))/2.d0
!      cc(4,1) = beta4*uu*cc(0,0) + beta2*(1 + uu)*cc(1,0) + 2*beta1
!     &     *cc(2,0) +  beta1*uu*cc(2,0) + 3*beta0*cc(3,0) + beta0
!     &     *uu*cc(3,0)

      sigmu=0.d0
      do i=0,nord
         do j=0,i
            do k=0,j
               sigmu(i) = sigmu(i) + api**j*cc(j,k)*lm0**k
            enddo
         enddo
      enddo

      return
      end
