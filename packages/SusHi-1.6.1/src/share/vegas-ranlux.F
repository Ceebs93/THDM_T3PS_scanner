
c VEGAS

c --------------------------------------------------

      subroutine vegas(fxn,avgi,sd,chi2a)
c
c
c  subroutine performs n-dimensional monte carlo integ*n
c    - by g.p.  lepage  sept 1976/ (rev) apr 1978
c
c    -ftn5 version 21-8-1984
c    -hbook/hplot interface 6-1-1985
c
c  author                                       : g. p. lepage
c
      implicit double precision (a-h,o-z)
      external fxn
      real qran4(10)
      dimension xin(50),r(50),dx(10),ia(10),kg(10),dt(10)
      dimension qran(10),x(10)
      common/vgasio/ninp,noutp
      common/vgb2/ndo,it,si,si2,swgt,schi,xi(50,10),scalls
     +  ,d(50,10),di(50,10)
      COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
c
c
      data ndmx/50/,alph/1.5/,one/1./,mds/1/
c
c     call vgdat

      if(itmx.le.0)then
        write(noutp,199)'0vegas called with at max less equal zero'//
     +                   ' iterations. no execution.'
        return
      endif
      if(nprn.gt.0)then
        ipr=0
      else
        ipr=1
      endif
      ndo=1
      do 1 j=1,ndim
        xi(1,j)=one
    1 continue
c
      entry vegas1(fxn,avgi,sd,chi2a)
c...initializes cummulative variables,but not grid
c     call vgdat
      if(itmx.le.0)then
        write(noutp,199)'0vegas1 called with at max less equal zero'//
     +                   ' iterations. no execution.'
        return
      endif
c
      it=0
      si=0.
      si2=si
      swgt=si
      schi=si
      scalls=si
c
      entry vegas2(fxn,avgi,sd,chi2a)
c...no initialization
c     call vgdat
      if(itmx.le.0)then
        write(noutp,199)'0vegas2 called with at max less equal zero'//
     +                   ' iterations. no execution.'
        return
      endif
      nd=ndmx
      ng=1
      if(mds.ne.0)then
        ng=(ncall/2.)**(1./ndim)
        mds=1
        if((2*ng-ndmx).ge.0)then
          mds=-1
          npg=ng/ndmx+1
          nd=ng/npg
          ng=npg*nd
        endif
      endif
c
      k=ng**ndim
      npg=ncall/k
      if(npg.lt.2) npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=dxg**(2*ndim)/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one
      do 3 j=1,ndim
        dx(j)=xu(j)-xl(j)
        xjac=xjac*dx(j)
    3 continue
c
c  rebin preserving bin density
c
      if(nd.ne.ndo)then
        rc=ndo/xnd
        do 7 j=1,ndim
          k=0
          xn=0.
          dr=xn
          i=k
    4     k=k+1
          dr=dr+one
          xo=xn
          xn=xi(k,j)
    5     if(rc.gt.dr) go to 4
          i=i+1
          dr=dr-rc
          xin(i)=xn-(xn-xo)*dr
          if(i.lt.ndm) go to 5
          do 6  i=1,ndm
            xi(i,j)=xin(i)
    6     continue
          xi(nd,j)=one
    7   continue
        ndo=nd
      endif
c
      if(nprn.ne.0.and.nprn.ne.10)write(noutp,200)ndim,calls,it,itmx
     +  ,acc,mds,nd
      if(nprn.eq.10)write(noutp,290)ndim,calls,itmx,acc,mds,nd
c
      entry vegas3(fxn,avgi,sd,chi2a)
c     - main integration loop
      if(itmx.le.0)then
        write(noutp,199)'0vegas3 called with at max less equal zero'//
     +                   ' iterations. no execution.'
        return
      endif
    9 continue
      it=it+1
      ti=0.
      tsi=ti
c
      do 10 j=1,ndim
        kg(j)=1
        do 10 i=1,nd
          d(i,j)=ti
          di(i,j)=ti
   10 continue
c
   11 fb=0.
      f2b=fb
      k=0
c
   12 continue
      k=k+1
      CALL RanLux(QRAN4,NDIM)
      do i=1,10
         qran(i) = qran4(i)
      enddo
      wgt=xjac
      do 15 j=1,ndim
        xn=(kg(j)-qran(j))*dxg+one
        ia(j)=xn
        iaj=ia(j)
        iaj1=iaj-1
        if(iaj.le.1)then
          xo=xi(iaj,j)
          rc=(xn-iaj)*xo
        else
          xo=xi(iaj,j)-xi(iaj1,j)
          rc=xi(iaj1,j)+(xn-iaj)*xo
        endif
        x(j)=xl(j)+rc*dx(j)
        wgt=wgt*xo*xnd
   15 continue
c
      f=fxn(x)*wgt
      f1=f/calls
      w=wgt/calls
c
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do 16 j=1,ndim
        iaj=ia(j)
        di(iaj,j)=di(iaj,j)+f/calls
        if(mds.ge.0)  d(iaj,j)=d(iaj,j)+f2
   16 continue
      if(k.lt.npg) go to 12
c
      f2b=f2b*npg
      f2b=dsqrt(f2b)
      f2b=dabs((f2b-fb)*(f2b+fb))
      ti=ti+fb
      tsi=tsi+f2b
      if(mds.lt.0)then
        do 17 j=1,ndim
          iaj=ia(j)
          d(iaj,j)=d(iaj,j)+f2b
   17   continue
      endif
      k=ndim
   19 kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1) go to 11
      k=k-1
      if(k.gt.0) go to 19
c
c  final results for this iteration
c
      ti=ti/calls
      tsi=tsi*dv2g
      ti2=ti*ti
      if(tsi .eq. 0.)then
        wgt = 0.
      else
        wgt=ti2/tsi
      endif
      si=si+ti*wgt
      si2=si2+ti2
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      if(swgt .eq. 0.)then
        avgi=ti
      else
        avgi=si/swgt
      endif
      if(si2 .eq. 0.)then
        sd=tsi
      else
        sd=swgt*it/si2
      endif
      scalls=scalls+calls
      chi2a=0.
      if(it.gt.1)chi2a=sd*(schi/swgt-avgi*avgi)/(it-1)
      if(sd .ne. 0.)then
        sd=one/sd
        sd=dsqrt(sd)
      else
        sd=tsi
      endif
      if(nprn.ne.0)then
        tsi=dsqrt(tsi)
        if(nprn.ne.10)write(noutp,201)ipr,it,ti,tsi,avgi,sd,chi2a
        if(nprn.eq.10)write(noutp,203)it,ti,tsi,avgi,sd,chi2a
        if(nprn.lt.0)then
          do 20 j=1,ndim
            write(noutp,202)j
            write(noutp,204)(xi(i,j),di(i,j),d(i,j),i=1,nd)
   20     continue
        endif
      endif
c
c   refine grid
c
      if(sd .ne. 0.)then
        rel = dabs(sd/avgi)
      else
        rel = 0.
      endif
      if(rel.le.dabs(acc).or.it.ge.itmx)now=2
c
      do 23 j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2.
        dt(j)=d(1,j)
        do 22 i=2,ndm
          d(i,j)=xo+xn
          xo=xn
          xn=d(i+1,j)
          d(i,j)=(d(i,j)+xn)/3.
          dt(j)=dt(j)+d(i,j)
   22   continue
        d(nd,j)=(xn+xo)/2.
        dt(j)=dt(j)+d(nd,j)
   23 continue
c
      do 28 j=1,ndim
        rc=0.
        do 24 i=1,nd
          r(i)=0.
          if(d(i,j).gt.0.)then
            xo=dt(j)/d(i,j)
            r(i)=((xo-one)/xo/dlog(xo))**alph
          endif
          rc=rc+r(i)
   24   continue
        rc=rc/xnd
        k=0
        xn=0.
        dr=xn
        i=k
   25   k=k+1
        dr=dr+r(k)
        xo=xn
        xn=xi(k,j)
   26   if(rc.gt.dr) go to 25
        i=i+1
        dr=dr-rc
        if(dr .eq. 0.)then
          xin(i)=xn
        else
          xin(i)=xn-(xn-xo)*dr/r(k)
        endif
        if(i.lt.ndm) go to 26
        do 27 i=1,ndm
          xi(i,j)=xin(i)
   27   continue
        xi(nd,j)=one
   28 continue
c
      if(it.lt.itmx.and.dabs(acc).lt.rel)go to 9
c
      return
c
  199 format(a)
  200 format('0input parameters for vegas   ndim=',i3
     +,'  ncall=',f8.0/28x,'  it=',i5,'  itmx =',i5/28x
     +,'  acc=',d9.3/28x,'  mds=',i3,'   nd=',i4//)
  290 format('0vegas  ndim=',i3,'  ncall=',f8.0,'  itmx =',i5
     +  ,'  acc=',d9.3,'  mds=',i3,'   nd=',i4)
  201 format(/i1,'integration by vegas'/'0iteration no',i3,
     +'.   integral =',d14.8/20x,'std dev  =',d10.4/
     +' accumulated results.   integral =',d14.8/
     +24x,'std dev  =',d10.4 / 24x,'chi**2 per itn   =',d10.4)
  202 format('0data for axis',i2 / 7x,'x',7x,'  delt i  ',
     +2x,' convce    ',11x,'x',7x,'  delt i  ',2x,' convce     '
     +,11x,'x',7x,'  delt i  ',2x,' convce     '/)
  204 format(1x,3d12.4,5x,3d12.4,5x,3d12.4)
  203 format(1x,i3,d20.8,d12.4,d20.8,d12.4,d12.4)
c
      end

c --------------------------------------------------
