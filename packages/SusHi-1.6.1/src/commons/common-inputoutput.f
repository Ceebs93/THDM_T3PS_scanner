!input specific variables
      character SUversion*10,SUdate*10
      character jfilein*60
      double precision SUinvAlfa,GF
      double precision ME, MU, MD, MM, MS, ML, tanb
      double precision MA0, MHp
      double precision M3SL, M3SE, M3SQ, M3SU, M3SD
      double precision M2SL, M2SE, M2SQ, M2SU, M2SD
      double precision M1SL, M1SE, M1SQ, M1SU, M1SD
      double precision Xt,MSUSY
      double complex MUE,M_1,M_2,M_3
      double complex cAt,cAtau,cAc,cAs,cAmu,cAu,cAd,cAe,cAb
      double precision scalefactor
      double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
      double precision At
      double precision Qtau, Qt, Qb
      double precision yukfac(9),Mcharm,Mbp,Mtp,q0,q0c,q01,q02
      double precision maxptc,maxrapc,minptc,minrapc,ptref
      double precision lam, kap, Alam, Akap
      integer ptc,rapc,mbrunyuk,mbrunloop,delmbfh,mbsbch,Abch,thetabch
      integer model, twohdmver, thdmckey, thdmcbreak
      integer mssmpart, fieldren, tanbren, higgsmix, p2approx
      integer looplevel, loglevel, runningMT, botResum, tlCplxApprox

!output specific variables
      character jfnameact*60
      double precision c1sm13,c1sm23,c1sm33,c1susy1t3,c1susy2t3
      double precision c1susy1b3,c1susy2b3,gb3,gt3,gb13,gbh123
      double precision gbh213,gb23,gt13,gth123,gth213,gt23
      double precision c1sm14,c1sm24,c1sm34,c1susy1t4,c1susy2t4
      double precision c1susy1b4,c1susy2b4,gt4,gb4,gt14,gt24
      double precision gb14,gb24,gth124,gth214,gbh124,gbh214
      double precision xsec(2),xsecerr(2),xeff,xefferr
      double precision xeffvar(2),gghXSvar(2)
      double precision ggchan,qgchan,qqchan,ggchanerr,qgchanerr
     &     ,qqchanerr
      double precision elw
      double precision SUggh(4,4,-1:5),SUggherr(4,4,-1:5)
      double precision SUbbh(3),SUbbherr(3)
      double precision awithpt,awopt
      integer HB_result,HB_chan,HB_ncombined
      double precision HB_obsratio
      integer HS_nobs
      double precision HS_Chisq_mu,HS_Chisq_mh,HS_Chisq,HS_Pvalue

      common /FHflags/ mssmpart,fieldren,tanbren,higgsmix,p2approx,
     &looplevel,loglevel,runningMT,botResum,tlCplxApprox

      common /HBandHSoutput/ HB_obsratio,HS_Chisq_mu,HS_Chisq_mh,
     &HS_Chisq,HS_Pvalue,HB_result,HB_chan,HB_ncombined,HS_nobs

      common /inputvarint/ twohdmver,thdmckey,thdmcbreak,model,
     &     mbrunloop,thetabch,ptc,rapc,mbrunyuk,delmbfh,mbsbch,Abch

      common /inputvarcomplex/ cAt,cAtau,cAc,cAs,cAmu,cAu,
     &cAd,cAe,cAb,MUE,M_1,M_2,M_3 

      common /nmssmpar/ lam,kap,Alam,Akap

      common /inputvar/ SUinvAlfa,GF,ME,MU,MD,
     &MM,MS,ML,tanb,MA0,MHp,
     &M3SL,M3SE,M3SQ,M3SU,M3SD,
     &M2SL,M2SE,M2SQ,M2SU,M2SD,
     &M1SL,M1SE,M1SQ,M1SU,M1SD,
     &Xt,MSUSY,scalefactor,
     &CKMlambda,CKMA,CKMrhobar,CKMetabar,At,
     &Qtau,Qt,Qb,Mcharm,Mbp,Mtp,q0,q0c,q01,q02,
     &maxptc,maxrapc,minptc,minrapc,ptref,yukfac

      common /outputvar/ c1sm13,c1sm23,c1sm33,c1susy1t3,c1susy2t3,
     &     c1susy1b3,c1susy2b3,gb3,gt3,gb13,gbh123, gbh213,gb23,gt13
     &     ,gth123,gth213,gt23, c1sm14,c1sm24,c1sm34,c1susy1t4,c1susy2t4
     &     , c1susy1b4,c1susy2b4,gt4,gb4,gt14,gt24, gb14,gb24,gth124
     &     ,gth214,gbh124,gbh214, xsec,xsecerr,SUggh,SUggherr,SUbbh
     &     ,SUbbherr,ggchan,qgchan,qqchan
     &     ,xefferr,ggchanerr,qgchanerr ,qqchanerr,xeff,awithpt
     &     ,awopt,xeffvar,gghXSvar

      common /ew/ elw
      common /verdate/ SUversion,SUdate

      common /fileinfileout/ jfilein,jfnameact
