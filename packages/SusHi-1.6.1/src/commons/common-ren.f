      double precision Mgl,
     &     cthetab,sthetab,cthetat,sthetat,muSUSY,alpha,beta,dmb(3),
     &     dmbsbos(3),dmbsb(3),dAbds(3),
     &     dthetab(3),dmsb1(2),dmsb2(2),
     &     c2thetat,s2thetat,c2thetab,s2thetab,
     &     t2thetab,lb,lb1,lb2,lt,lt1,lt2,lg,
     &     mcmc,mcos,mc,mbos,mb,mbDRSMmuR,mbDRSMmt,mbDRSMmb,mbDRSMmuB,
     &     mbDRSMmuD,mbDRMSmuR,mbDRMSmb,mbDRMSmuB,mbDRMSmt,mbDRMSmuD,
     &     mbsb,mbossave,mbyuk,mbossbot,Abdrbar,
     &     Ab,thetab,dMLb,delta_mb,delmb
      
      double precision Hmx(3,3), Amx(3,3)
      
      integer Hind, Aind, Sind, htype
      
      common /evcsnew/ Mgl,
     &     cthetab,sthetab,cthetat,sthetat,muSUSY,alpha,beta,
     &     dmb,dmbsbos,dmbsb,dAbds,dthetab,dmsb1,dmsb2
      
      common /asdfg/ c2thetat,s2thetat,c2thetab,s2thetab,
     &     t2thetab,lb,lb1,lb2,lt,lt1,lt2,lg
      
      common /mixing/ Hmx, Amx, Hind, Aind, Sind, htype
      
      common /renorm/ mcmc,mcos,mc,mbos,mb,mbDRSMmuR,mbDRSMmt,
     &     mbDRSMmuB,mbDRSMmuD,mbDRMSmuR,mbDRMSmt,mbDRMSmuB,mbDRMSmuD,
     &     mbsb,mbossave,mbyuk,mbossbot,Abdrbar,Ab,
     &     thetab,dMLb,delta_mb,delmb
