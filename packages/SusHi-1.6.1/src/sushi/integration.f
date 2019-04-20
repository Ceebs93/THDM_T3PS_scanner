! This file is part of SusHi.
! 
! It includes the relevant routines for phase-space
! integration and the convolution with PDFs.
! 
C-{{{ integration

      subroutine integrate(ndim,integrand,integral,error,chi2)
      implicit none
      double precision integral,error,intl,errl,
     &xl(10),xu(10),
     &acc1,chi2,chi2l
      integer ndim,ndm,ncall,itmx,nprn,iv
      common /bveg1/ xl,xu,acc1,ndm,ncall,itmx,nprn
      include '../commons/common-vegpar.f'
      external integrand

      acc1 = 1.d-8
      xl = 0.d0
      xu = 1.d0
      ndm = ndim
      itmx = itmx1sushi
      ncall = ncall1sushi
      nprn = nprnvsushi

      call vegas(integrand,intl,errl,chi2l)
      if (lveg1sushi) then
         itmx = itmx2sushi
         ncall = ncall2sushi
         call vegas1(integrand,intl,errl,chi2l)
      endif
      integral = intl
      error = errl
      chi2 = chi2l
      end

C-}}}
C-{{{ variable transformations

c     All variable transformations below lead to integrals with bounds 0 and 1

C-{{{ transformations

      function trans(integrand,xx)
c     performs the integral
c     int_0^1 dv int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) integrand(z/x1/x2,s,t,u,x1,x2)
c     with z = mh**2 / cme**2, cme = hadronic center of mass energy
c     using the following variable transformation:
c     x1 -> w1 / w2
c     x2 -> z / w1
c     with
c     w1 = tau1 * tau2 * (1 - z) + z
c     w2 = tau1 * (1 - z) + z
c     and with the Jacobian determinant
c     | d(x1, x2)/d(tau1, tau2) | = (1 - z)**2 * tau1 * x2 / w2**2
      implicit none
      double precision trans,integrand,xx(*),
     &totxs,tau1,tau2,tau3,w1,w2,x1,x2,jdet,s,t,u
     &,totxs2
      include '../commons/common-int.f'

      if(pty) then
c     integrate over pt and y
         trans = totxs(integrand,xx)
c         trans = totxs2(integrand,xx)
         return
      endif

      tau1 = xx(1)
      tau2 = xx(2)
      tau3 = xx(3)
      w1 = tau1 * tau2 * (1 - z) + z
      w2 = tau1 * (1 - z) + z
      x1 = w1 / w2
      x2 = z / w1
      jdet = (1 - z)**2 * tau1 * x2 / w2**2

c     for the phase space integration:
      s = mh2 / w2
      t = mh2 * (1 - 1 / w2) * tau3
      u = mh2 * (1 - 1 / w2) * (1 - tau3)

c     calculate pt and y
      if(pseudorap) then
         y = dlog(x1*t/x2/u)/2.d0
      else
         y = dlog(x1*(mh2-t)/x2/(mh2-u))/2.d0
      endif
      pt = dsqrt(t*u/s)

      trans = integrand(w2,s,t,u,x1,x2) * jdet
      end

      function transplus(integrand1,integrand2,xx)
c     performs the integral
c     int_0^1 dv int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) integrand1(z/x1/x2,s,t,u,x1,x2)
c     where integrand1 is a Catani Seymour subtraction term (integrand2 contains only the plus distributions from integrand1)
c     using the following variable transformation:
c     x1 -> exp( y) * sqrt(z)
c     x2 -> exp(-y) * sqrt(z) / w2
c     with
c     w2 = tau1 * (1 - z) + z
c     y = ymin + (ymax - ymin) * tau2
c     and with the Jacobian determinant
c     | d(x1, x2)/d(tau1, tau2) | = z / w2**2 * (1 - z) (ymax - ymin)
      implicit none
      double precision transplus,integrand1,integrand2,xx(*),
     &tau1,tau2,w2,x1,x2,jdet1,jdet2
      include '../commons/common-int.f'

      tau1 = xx(1)
      tau2 = xx(2)
      jdet1 = 1 - z
      jdet2 = maxrap - minrap
      y = minrap + jdet2 * tau2
      x1 = dexp(y)*dsqrt(z)
      x2 = z/x1
      w2 = tau1 * (1 - z) + z
      transplus = -integrand2(w2,x1,x2) * z * jdet1 * jdet2
      transplus = transplus * 2.d0 ! add subtraction terms for t -> 0 and u -> 0
      if(x1.le.w2) then         ! subtraction term for t -> 0
         transplus = transplus
     &        + integrand1(w2,x1/w2,x2)
     &        * z/w2**2 * jdet1 * jdet2
      endif
      if(x1.ge.z/w2) then       ! subtraction term for u -> 0
         transplus = transplus
     &        + integrand1(w2,x2/w2,x1)
     &        * z/w2**2 * jdet1 * jdet2
      endif
      transplus = transplus * 2.d0 ! to include negative y
      end

      function transdelta(integrand,xx)
c     performs the integral
c     int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) pdf(x1)/x1 pdf(x2)/x2 * Dirac_delta(1 - z/x1/x2)
c   = int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) pdf(x1)/x1 pdf(x2)/x2 * x2 * Dirac_delta(x2 - z/x1)
c     using the following variable transformation:
c     x1 -> exp(y) * sqrt(z)
c     and
c     y ->  ymin + (ymax - ymin) * tau
c     with the Jacobian determinant
c     | dx1/dy | * | dy/dtau | = x1 * (ymax - ymin)
      implicit none
      double precision transdelta,integrand,xx(*),tau,x1,x2,jdet1,jdet2
      include '../commons/common-int.f'
      tau = xx(1)
      jdet1 = maxrap - minrap
      y = minrap + jdet1 * tau
      x1 = dexp(y)*sqrt(z)
      jdet2 = x1
      x2 = z / x1
      transdelta = x2 * jdet1 * jdet2 * integrand(x1,x2)
      transdelta = transdelta * 2.d0 ! to include negative y
      end

C-}}}
C-{{{ transformation {x1, x2, t} -> {x1, y, pt}

      function totxs(integrand,xx)
c     total cross section
      implicit none
      double precision totxs,integrand,xx(*),ydistrib
      external integrand
      totxs = ydistrib(integrand,xx,3)
      end

      function totxs2(integrand,xx)
c     total cross section (alternative version, just for cross checks)
      implicit none
      double precision totxs2,integrand,xx(*),ptdistrib
      external integrand
      totxs2 = ptdistrib(integrand,xx,3)
      end

      function ptdistrib(integrand,xx,i)
c     performs the integral
c     int dy Heaviside_theta(|y| - ymin) Heaviside_theta(ymax - |y|) integrand(y)
c     using the following variable transformation:
c     y ->   ymin + (ymax - ymin) * tau
c     with the Jacobian determinant
c     | dy/dtau | = ymax - ymin
c
c     i = 3: total cross section
c     i = 2: pt distribution
      implicit none
      double precision ptdistrib,integrand,xx(*),tau,ymin,ymax,jdet,
     &ydistrib,ptydistrib
      integer i
      include '../commons/common-int.f'
      external integrand

      tau = xx(i)
      ymin = 0.d0
      if(i.eq.3) then
         if(pseudorap) then
            call printerrorsushi(1,6,"Error in ptdistrib. Stop")            
         endif
         ymax = -dlog(z)/2.d0
      else if(i.eq.2) then
         if(pseudorap) then
            ymax = dlog( (1-z+dsqrt((1-z)**2-4*pt**2/mh2*z))
     &           /2.d0/(pt*dsqrt(z/mh2)) )
         else
            ymax = dlog( (1+z+dsqrt((1-z)**2-4*pt**2/mh2*z))
     &           /2.d0/dsqrt((1+pt**2/mh2)*z) )
         endif
      else
         call printerrorsushi(1,6,"Error in ptdistrib. Stop")
      endif
c     if subtraction term is active, integration over
c     the whole y-range is necessary
      if(subtr) then
         jdet = ymax
         y = jdet * tau
      else
         ymin = minrap
         if(ymin.gt.ymax) then
            ptdistrib = 0.d0
            return
         endif
         if((.not.pseudorap).and.(ymax.gt.maxrap)) then
            ymax = maxrap
         endif
         jdet = ymax - ymin
         y =   ymin + jdet * tau
      endif

      if(i.eq.3) then
         ptdistrib = jdet * ydistrib(integrand,xx,2)
      else
         ptdistrib = jdet * ptydistrib(integrand,xx)
      endif
      ptdistrib = ptdistrib * 2.d0 ! to include negative y
      end

      function ydistrib(integrand,xx,i)
c     Performs the integral
c     int_ptmin^ptmax dpt integrand(pt)
c     using the following variable transformation:
c     pt -> ptmin + (ptmax - ptmin) * tau
c     with the Jacobian determinant
c     | dpt/dtau | = ptmax - ptmin
c
c     i = 3: total cross section
c     i = 2: y distribution
      implicit none
      double precision ydistrib,integrand,xx(*),tau,ptmin,ptmax,jdet,
     &ptdistrib,ptydistrib
      integer i
      include '../commons/common-int.f'
      include '../commons/common-vars.f'
      include '../commons/common-inputoutput.f'

      external integrand

      tau = xx(i)
      ptmin = minpt
      if(i.eq.3) then
         ptmax = maxpt
      else if(i.eq.2) then
         if(pseudorap) then
            ptmax = (1-z)*dsqrt(mh2/z)/(dexp(y)+dexp(-y))
         else
            ptmax = dsqrt(mh2*(dexp(2*y)*(1/z+z)-1-dexp(4*y)))
     &           /(1+dexp(2*y))
         endif
         if(ptmin.ge.ptmax) then
            ydistrib = 0.d0
            return
         endif
      else
         call printerrorsushi(1,6,"Error in ydistrib. Stop")
      endif
c     if subtraction term is active, integration over
c     the whole pt-range is necessary
      if(.not.subtr) then
         if(ptmax.gt.maxptc) then
            ptmax = maxptc
         endif
      endif
      jdet = ptmax - ptmin
      pt = ptmin + jdet * tau
      if(i.eq.3) then
         ydistrib = jdet * ptdistrib(integrand,xx,2)
      else
         ydistrib = jdet * ptydistrib(integrand,xx)
      endif
      end

      function ptydistrib(integrand,xx)
c     performs the integral
c     int_0^1 dv int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) pdf(x1)/x1 pdf(x2)/x2 integrand(z/x1/x2,s,t,u,x1,x2)
c     using one of the following variable transformations:
c     v  -> (sqrt(z)*x1*(p1**2 - 1)) / (z*exp(y)*p1 - 2*sqrt(z)*x1 + x1**2*exp(-y)*p1)
c     x2 -> (z - sqrt(z)*x1*exp(-y)*p1) / (sqrt(z)*exp(y)*p1 - x1)
c     with p1 = sqrt(1 + pt**2/mh**2)
c     and with the Jacobian determinant
c     | d(x2, v)/d(pt,  y ) | = 2 * x2 * sqrt(s * t / u) / (s - mh**2)
c     for rapidity or
c     v  -> (z*exp(y) + sqrt(z)*x1*pt/mh) / (z*exp(y) + exp(-y)*x1**2)
c     x2 -> (z*exp(y) + sqrt(z)*x1*pt/mh) / exp(y) / (x1 - sqrt(z)*exp(y)*pt/mh)
c     with the Jacobian determinant
c     | d(x2, v)/d(pt, y) | = 2 * x2 * sqrt(s * t * u) / (s + t) / (s - mh**2)
c     for pseudorapidity
c     note: this function performs only the integration over x1 with
c     x1 = x1min + (x1max - x1min) * tau
c     with the Jacobian determinant
c     | dx1/dtau | = x1max - x1min
      implicit none
      double precision ptydistrib,integrand,xx(*),
     &tau1,x1max,x1min,mhpt,x1,x2,w2,s,t,u,jdet1,jdet2,scale,alphaspdfm
      include '../commons/common-int.f'
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'

      tau1 = xx(1)
      x1max = 1.d0
      if(pseudorap) then
         mhpt = dsqrt(mh2 * pt**2 * z)
         x1min = (dexp(y) * mhpt + mh2 * z) / (mh2 - dexp(-y) * mhpt)
         jdet2 = x1max - x1min
         x1 = x1min + jdet2 * tau1
         u = - x1 * dexp(-y) * mhpt / z
      else
         mhpt = dsqrt(mh2 * (mh2 + pt**2) * z)
         x1min = (dexp(y) * mhpt - mh2 * z) / (mh2 - dexp(-y) * mhpt)
         jdet2 = x1max - x1min
         x1 = x1min + jdet2 * tau1
         u = mh2 - x1 * dexp(-y) * mhpt / z
      endif
      x2 = (mh2 - u) * u / (pt**2 + u) / cme**2 / x1
      w2 = z / x1 / x2
      s = mh2 / w2
      t = mh2 - s - u

c     running scale
      if(((dist.eq.1).or.(dist.eq.3)).and.(scalespt)) then
         muFggh = dsqrt(mh2+pt**2)*muFfacggh
         apimz = alphasPDFm(2,mz)/Pi
         call runalpha(apimz,mz,dsqrt(mh2+pt**2)*muRfacggh
     &        ,nf,2,0,apimuR)
         scale = pi*apimur
      else
         scale = 1.d0
      endif

      if(pseudorap) then
         jdet1 = 2 * x2 * dsqrt(s * t * u) / (s + t) / (s - mh2)
      else
         jdet1 = 2 * x2 * dsqrt(s * t / u) / (s - mh2)
      endif

      ptydistrib = integrand(w2,s,t,u,x1,x2) * jdet1 * jdet2 * scale**3
      end

C-}}}

C-}}}

C-{{{ integrand for delta function

      function deltaggint(x1,x2)
      implicit none
      double precision deltaggint,x1,x2,PDFgg
      deltaggint = PDFgg(x1,x2) / ( x1 * x2 )
      end

      function deltagg(xx)
      implicit none
      double precision deltagg,xx(*),transdelta,deltaggint
      external deltaggint
      deltagg = transdelta(deltaggint,xx)
      end

C-}}}
C-{{{ integrand for Catani-Seymour subtraction term (gg channel)
      
      function intsubgg(xx)
      implicit none
      double precision intsubgg,xx(*),transplus,subggint1,subggint2
      external subggint1,subggint2
      intsubgg = transplus(subggint1,subggint2,xx)
      end

      function subggint1(x,x1,x2)
c     integrated Catani Seymour subtraction term
      implicit none
      double precision subggint1,x,x1,x2,ggterms,dterms,PDFgg
      subggint1 = ( ggterms(x) + dterms(x) ) / 2.d0
     &     * PDFgg(x1,x2) / (x1 * x2)
      end

      function subggint2(x,x1,x2)
c     integrated Catani Seymour subtraction term (only plus distributions)
      implicit none
      double precision subggint2,x,x1,x2,ggterms,dterms,PDFgg
      subggint2 = dterms(x) / 2.d0 * PDFgg(x1,x2) / (x1 * x2)
      end

      function ggterms(x)
      implicit none
      double precision ggterms,x
      include '../commons/common-vars.f'
      ggterms = (4*dlog(1 - x) - 2*lfh)*(-2*x + x**2 - x**3)
     &     -2*dlog(x)*(1 - x + x**2)**2/(1.d0 - x)
      end

      function dterms(x)
      implicit none
      double precision dterms,x,dd0,dd1
      include '../commons/common-consts.f'
      include '../commons/common-vars.f'
      dd0=1/(1.d0-x)
      dd1=dlog(1.d0-x)*dd0
      dterms = - 2 * lfh * dd0 + 4 * dd1
      end

      function dtermsz()
      implicit none
      double precision dtermsz,dd0,dd1
      include '../commons/common-consts.f'
      include '../commons/common-int.f'
      include '../commons/common-vars.f'
      dd0=dlog(1-z)/(1-z)
      dd1=dlog(1-z)/2.d0*dd0
      dtermsz = -2 * lfh * dd0 + 4 * dd1
      end

C-}}}
C-{{{ integrand for Catani-Seymour subtraction term (qg and gq channel)

      function intsubqg(xx)
      implicit none
      double precision intsubqg,xx(*),transplus,subqgint1,subqgint2
      external subqgint1,subqgint2
      intsubqg = transplus(subqgint1,subqgint2,xx)
      end

      function subqgint1(x,x1,x2)
c     integrated Catani Seymour subtraction term
      implicit none
      double precision subqgint1,x,x1,x2,qgterms,PDFqg
      subqgint1 = qgterms(x) * PDFqg(x1,x2) / (x1 * x2)
      end

      function subqgint2(x,x1,x2)
c     integrated Catani Seymour subtraction term (only plus distributions)
      implicit none
      double precision subqgint2,x,x1,x2
      subqgint2 = 0.d0
      end

      function qgterms(x)
      implicit none
      double precision qgterms,x
      include '../commons/common-vars.f'
      qgterms = ((1 + (1 - x)**2)*(2*dlog(1 - x)-dlog(x) - lfh) + x**2)
      end

C-}}}
C-{{{ integrand for real corrections (gg channel)

      function realggint(x,s,t,u,x1,x2)
      implicit none
      double precision realggint,x,s,t,u,x1,x2,tmp,AMPgg,
     &jet,jetsub,subggt,subggu,PDFgg
      include '../commons/common-int.f'

      if(t.ge.0.d0) then
         t = -s*1.d-16
      endif
      if(u.ge.0.d0) then
         u = -s*1.d-16
      endif

      tmp = AMPgg(s,t,u)
     &     * jet(pt,y)
      if(subtr) then
         if((-s/t.gt.1.d5).or.(-s/u.gt.1.d5)) then
            realggint = 0.d0
            return
         endif
         tmp = tmp
     &        - subggt(s,t,u) * jetsub(x1*x,x2)
     &        - subggu(s,t,u) * jetsub(x1,x2*x)
      endif
      realggint = x * (1 - x) * tmp * PDFgg(x1,x2) / ( x1 * x2 )
      end

c     total XS
      function realgg(xx)
      implicit none
      double precision realgg,xx(*),trans,realggint
      external realggint
      realgg = trans(realggint,xx)
      end

c     pt dist
      function realggpt(xx)
      implicit none
      double precision realggpt,xx(*),ptdistrib,realggint
      external realggint
      realggpt = ptdistrib(realggint,xx,2)
      end

c     y dist
      function realggy(xx)
      implicit none
      double precision realggy,xx(*),ydistrib,realggint
      external realggint
      realggy = ydistrib(realggint,xx,2)
      end

c     pt-y dist
      function realggpty(xx)
      implicit none
      double precision realggpty,xx(*),realggint,ptydistrib
      external realggint
      realggpty = ptydistrib(realggint,xx)
      end

C-}}}
C-{{{ integrand for real corrections (qg channel)

      function realqgint(x,s,t,u,x1,x2)
      implicit none
      double precision realqgint,x,s,t,u,x1,x2,tmpqg,tmpgq,AMPqg,
     &jet,subqg,jetsub,PDFqg
      include '../commons/common-int.f'
      tmpqg = AMPqg(s,t,u)
     &     * jet(pt,y)
      tmpgq = AMPqg(s,u,t)
     &     * jet(pt,y)
      if(subtr) then
         tmpqg = tmpqg
     &        - subqg(s,t,u) * jetsub(x1*x,x2)
         tmpgq = tmpgq      
     &        - subqg(s,u,t) * jetsub(x1,x2*x)
      endif
      realqgint = x * (1 - x)
     &     * ( tmpqg * PDFqg(x1,x2) / ( x1 * x2 )    ! qg channel
     &       + tmpgq * PDFqg(x2,x1) / ( x2 * x1 ) )  ! gq channel
      end

c     total XS
      function realqg(xx)
      implicit none
      double precision realqg,xx(*),trans,realqgint
      external realqgint
      realqg = trans(realqgint,xx)
      end

c     pt dist
      function realqgpt(xx)
      implicit none
      double precision realqgpt,xx(*),ptdistrib,realqgint
      external realqgint
      realqgpt = ptdistrib(realqgint,xx,2)
      end

c     y dist
      function realqgy(xx)
      implicit none
      double precision realqgy,xx(*),ydistrib,realqgint
      external realqgint
      realqgy = ydistrib(realqgint,xx,2)
      end

c     pt-y dist
      function realqgpty(xx)
      double precision realqgpty,xx(*),realqgint,ptydistrib
      external realqgint
      realqgpty = ptydistrib(realqgint,xx)
      end

C-}}}
C-{{{ integrand for real corrections (qq channel)

      function AMPqqjet(s,t,u)
      implicit none
      double precision AMPqqjet,s,t,u,AMPqq
c     amplitude without analytic phase space integration
      AMPqqjet = (t**2 + u**2)/(t + u)**2*3/2.d0*AMPqq(s,t + u)
      end

      function realqqint(x,s,t,u,x1,x2)
      implicit none
      double precision realqqint,x,s,t,u,x1,x2,AMPqqjet,jet,
     &PDFqq
      include '../commons/common-int.f'

      realqqint = (1 - x)**3 * AMPqqjet(s,t,u)
     &     * PDFqq(x1,x2) / (x1 * x2) * jet(pt,y)
      end

c     total XS
      function realqq(xx)
      implicit none
      double precision realqq,xx(*),trans,realqqint
      external realqqint
      realqq = trans(realqqint,xx)
      end

c     pt dist
      function realqqpt(xx)
      implicit none
      double precision realqqpt,xx(*),ptdistrib,realqqint
      external realqqint
      realqqpt = ptdistrib(realqqint,xx,2)
      end

c     y dist
      function realqqy(xx)
      implicit none
      double precision realqqy,xx(*),ydistrib,realqqint
      external realqqint
      realqqy = ydistrib(realqqint,xx,2)
      end

c     pt-y dist
      function realqqpty(xx)
      double precision realqqpty,xx(*),realqqint,ptydistrib
      external realqqint
      realqqpty = ptydistrib(realqqint,xx)
      end

C-}}}

C-{{{ PDFs
      
C-{{{ gg

      function PDFgg(x1,x2)
      implicit none
      double precision PDFgg,x1,x2,
     &upv1,dnv1,usea1,dsea1,str1,sbar1,chm1,cbar1,bot1,bbar1,glu1,
     &upv2,dnv2,usea2,dsea2,str2,sbar2,chm2,cbar2,bot2,bbar2,glu2

      call GetPDFs(x1,upv1,dnv1,usea1,dsea1,
     &str1,sbar1,chm1,cbar1,bot1,bbar1,glu1)
      call GetPDFs(x2,upv2,dnv2,usea2,dsea2,
     &str2,sbar2,chm2,cbar2,bot2,bbar2,glu2)

      PDFgg = glu1 * glu2
      end

C-}}}
C-{{{ qg

      function PDFqg(x1,x2)
      implicit none
      double precision PDFqg,x1,x2,q1,
     &upv1,dnv1,usea1,dsea1,str1,sbar1,chm1,cbar1,bot1,bbar1,glu1,
     &upv2,dnv2,usea2,dsea2,str2,sbar2,chm2,cbar2,bot2,bbar2,glu2

      call GetPDFs(x1,upv1,dnv1,usea1,dsea1,
     &str1,sbar1,chm1,cbar1,bot1,bbar1,glu1)
      call GetPDFs(x2,upv2,dnv2,usea2,dsea2,
     &str2,sbar2,chm2,cbar2,bot2,bbar2,glu2)

      q1 = upv1+dnv1+usea1+dsea1+str1+sbar1+chm1+cbar1+bot1+bbar1

      PDFqg = q1 * glu2
      end

C-}}}
C-{{{ qqbar

      function PDFqq(x1,x2)
      implicit none
      double precision PDFqq,x1,x2,
     &upv1,dnv1,usea1,dsea1,str1,sbar1,chm1,cbar1,bot1,bbar1,glu1,
     &upv2,dnv2,usea2,dsea2,str2,sbar2,chm2,cbar2,bot2,bbar2,glu2
      logical ppcoll
      common /coll/ ppcoll

      call GetPDFs(x1,upv1,dnv1,usea1,dsea1,
     &str1,sbar1,chm1,cbar1,bot1,bbar1,glu1)
      call GetPDFs(x2,upv2,dnv2,usea2,dsea2,
     &str2,sbar2,chm2,cbar2,bot2,bbar2,glu2)

      if(ppcoll) then
         PDFqq=upv1*usea2+upv2*usea1
     &        +dnv1*dsea2+dnv2*dsea1
     &        +str1*sbar2+str2*sbar1
     &        +chm1*cbar2+chm2*cbar1
     &        +bot1*bbar2+bot2*bbar1
      else
         PDFqq=upv1*upv2+usea1*usea2
     &        +dnv1*dnv2+dsea1*dsea2
     &        + 2 * (str1*str2
     &        +chm1*chm2+bot1*bot2)
      endif      
      end

C-}}}

C-}}}
C-{{{ subtraction terms

      function subggt(s,t,u)
      implicit none
      double precision subggt,s,t,u
      subggt = ( s**2 + s*(t + u) + (t + u)**2 )**2
     &     / ( s*t*(t + u)*(s + t + u) )
      end

      function subggu(s,t,u)
      implicit none
      double precision subggu,s,t,u
      subggu = ( s**2 + s*(t + u) + (t + u)**2 )**2
     &     / ( s*u*(t + u)*(s + t + u) )
      end

      function subqg(s,t,u)
      implicit none
      double precision subqg,s,t,u
      subqg = -( s**2 + (t + u)**2 ) / ( t*(s + t + u) )
      end

C-}}}
C-{{{ jet-functions

      function jetsub(x1,x2)
      implicit none
      double precision jetsub,x1,x2,jetuser,jet

      include '../commons/common-int.f'

      if(pseudorap) then
         jetsub = jet(0.d0,1.d100)
      else
         jetsub = jet(0.d0,dlog(x1/x2)/2.d0)
      endif
      end

      function jet(pt_in,y_in)
      implicit none
      double precision jet,pt_in,y_in,jetuser

      include '../commons/common-int.f'
      include '../commons/common-inputoutput.f'

      if(juser) then
         jet = jetuser(pt_in,y_in)
         return
      endif

      jet = 1.d0
      if((rapcut).and.(
     &     (dabs(y_in).gt.maxrap).or.(dabs(y_in).lt.minrap)) ) then
         jet = 0.d0
      endif
      if((ptcut).and.((pt_in.gt.maxptc).or.(pt_in.lt.minptc))) then
         jet = 0.d0
      endif
      end

c----------------------------------------------------
      function jetuser(pt,y)
      implicit none
      double precision jetuser,pt,y
      include '../commons/common-vars.f'
c     In this function, the user can apply arbitrary
c     cuts using the variables "pt" and "y", e.g.
c
c     jetuser = 0.d0
c     if((pt.gt.50.d0).and.(y.lt.3.d0)) then
c       jetuser = 1.d0
c     endif
c
c     applies the cuts pt > 50 GeV and y < 3.
c     NOTE that one needs to activate this function
c     by setting entry 5 in Block DISTRIB to 1
c     in the input file.
c
c     Please make changes only between
c     the two lines below:
c----------------------------------------------------

      jetuser = 0.d0
      if(.true.) then
         jetuser = 1.d0
      endif

c----------------------------------------------------
      end

C-}}}
