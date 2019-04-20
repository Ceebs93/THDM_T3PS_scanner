
!von bbh:
      real*8 atauh,prefac,sigbbh(6),sigmabbh(0:4)
      real*8 sigerrbbh(0:4),sallbbh(0:4),serrbbh(0:4)

      real*8 sig(6),sigma(0:3,-1:5),sigerr(0:3,-1:5),sall(0:3,
     &     -1:5)
!      real*8 BStint1
      common/sigmac/ prefac,sigma,sigerr,sall,sig
!      common/tbstint/ BStint1,TStint1,Ttint1
!von bbh:
      common/sigmacbbh/ atauh,sigmabbh,sigerrbbh,sallbbh,serrbbh,sigbbh

