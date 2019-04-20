      character*3 qqhstring
      double precision sqrts,tauh,elwfac
      double precision apimuR,mbpole,ataut2,nf,nfrun
      double precision apimzlist(4),apiggh(4)
      double precision mbMSbarmuR,mbMSbarmuRbbh,mbMSbarmt
      double precision mbMSbarmuB,mbMSbarmuD
      double precision c1sm0,c1sm1,c1sm2,c1sm3,c2sm1
      double precision c1eff0,c1eff1,c1eff2,c1eff3,c2eff1
      complex*16 comc1eff0
      double precision gth,gth11,gth12,gth21,gth22
      double precision gbh,gbh11,gbh12,gbh21,gbh22
      double precision msbp1,msbp2,mstp1,mstp2

      double precision rho0,apimz,myapimz,gfermi,mt,mbmb,qqhyuk,qqhscale
      double precision gt,gb,dgb
      double precision c5sp(0:1,0:3),c5spin(0:1,0:3)
      integer dim5run

      double precision alphaggh,muBfactor,AlfasMZ,betanull
      double precision sigmanull,cme
      double precision Mw,Mh,muB,muD,Mz
      double precision tparam(15),thdmcres(15),thdmcm12
      double precision murfacggh,murggh,rmur,rmuf,lfh,lrh,lft,lfr,lrt
     &     ,lth
      double precision murfacbbh,muffacbbh,murbbh,mufbbh,apimurbbh
      double precision facmurgghmin,facmurgghmax
      double precision facmurint1,facmurint2

      integer norderggh,n3lo,norderbbh,renscheme,tanbresum,ew
      integer facmurgghstep

      common/strings/qqhstring

      common /bbhnnlo/ murfacbbh,muffacbbh,murbbh,mufbbh,mbMSbarmuRbbh
     &     ,apimurbbh

      common /ggh/ murfacggh,murggh,rmur,rmuf,facmurgghmin,facmurgghmax
     &     ,facmurint1,facmurint2,dim5run
      common /gghint/ facmurgghstep

      common/params/ sqrts,tauh,ataut2,lfh,lrh,lft,lfr,lrt,lth,elwfac
     &     ,apimuR,mbMSbarmuR,mbMSbarmt ,mbMSbarmuB,mbMSbarmuD,mbpole,nf
     &     ,nfrun,c1sm0,c1sm1,c1sm2,c1sm3,c2sm1,c5sp,c5spin,comc1eff0
     &     ,c1eff0,c1eff1,c1eff2,c1eff3,c2eff1 ,gth,gth11,gth12,gth21
     &     ,gth22,gbh ,gbh11,gbh12,gbh21,gbh22

      common/vars/rho0,apimz,myapimz,apimzlist,apiggh,gfermi,mt,mbmb
     &     ,qqhyuk,qqhscale,gt,gb,dgb

      common /inp/ alphaggh,muBfactor,AlfasMZ,cme,
     &norderggh,n3lo,norderbbh,renscheme,ew,tanbresum

      common /constants/ betanull,sigmanull,Mw,Mz,muB,muD

      common /mh_sushi/ Mh

      common /fourthgen/ msbp1,msbp2,mstp1,mstp2

      common /thdmc/ tparam,thdmcres,thdmcm12

