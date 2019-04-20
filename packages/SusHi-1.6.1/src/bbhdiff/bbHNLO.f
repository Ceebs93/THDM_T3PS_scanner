      SUBROUTINE Sigma_bbHNLO(SigLO,errLO,Sig_bb,err_bb,
     &Sig_bg,err_bg,Sig_total,err_total)
      !Routine by Marius Wiesemann, rewritten for SUSHI by Stefan Liebler
      implicit none
      DOUBLE PRECISION MBorn,err_born,M_mbb,err_mbb,M_mbg,err_mbg
     f ,M_virt,err_virt,M_Insbb,err_insbb,M_KPbb,err_kpbb,M_bbgh,
     f M_gbbH,err_gbbH,SigLO,errLO,Sig_bb,err_bb,Sig_bg,err_bg,err_bbgh
     f ,M_KPbg,err_kpbg,Sig_total,err_total
     f ,XL,XU,acc1,chi2
      INTEGER ndim,ncall,ncall1,itmx,nprn,itmx1,orderSAVE
      COMMON/BVEG1/XL(10),XU(10),acc1,ndim,ncall,itmx,nprn
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-cuts.f'
      include '../commons/common-vegpar.f'
      include '../commons/common-bbhdiff-ptdist.f'
      external MbbgH,MgbbH,Born_bbH,Mvirt_bbH,Ins_bbgH,KP_bbgH,KP_gbbH

      !pT_Hmin=0d0
      acc1 = 1.d-8
      acc = acc1
      ndim = 1       ! number of integration variables

      XL(1) = 0.00000000001d0
      XL(2) = 0.00000000001d0
      XL(3) = -0.99999999999d0
      XL(4) = 0.00000000001d0    ! these are the integration limits
      XU(1) = 0.99999999999d0    ! XL = X(lower), XU = X(upper)
      XU(2) = 0.99999999999d0
      XU(3) = 0.99999999999d0
      XU(4) = 0.99999999999d0

      itmx=10 !itmx1sushi
      nprn=nprnvsushi
      ncall = ncall1sushi*3

      if ((jets == 1).or.(jets == 2).or.(dist == 1)) then

      M_virt = 0.d0
      err_virt = 0.d0
      M_Insbb = 0.d0
      err_insbb = 0.d0
      M_KPbb = 0.d0
      err_kpbb = 0.d0
      M_KPbg = 0.d0
      err_kpbg = 0.d0
      M_mbb = 0.d0
      err_mbb = 0.d0
      M_mbg = 0.d0
      err_mbg = 0.d0
      M_bbgH = 0.d0
      err_bbgH = 0.d0
      M_gbbH = 0.d0
      err_gbbH = 0.d0
      M_bbgH = 0.d0
      err_bbgH = 0.d0
      M_gbbH = 0.d0
      err_gbbH = 0.d0

      if (bbhdifforder.eq.0) then
      Sig_total = 0.d0
      err_total = 0.d0
      return
      endif      

      else

      call vegas(Born_bbH,MBorn,err_born,chi2)

!      write(*,*) 'Born (bb->H)   :   ',MBorn,' +/-',err_born

      if (bbhdifforder.eq.0) then
      Sig_total = MBorn
      err_total = err_born
      return
      endif      

      ncall = ncall1sushi*5    ! in each iteration
      call vegas(Mvirt_bbH,M_virt,err_virt,chi2)
      
234   ncall = ncall1sushi*10
      call vegas(Ins_bbgH,M_Insbb,err_insbb,chi2)

      ndim = 2
      ncall = ncall1sushi*10
      call vegas(KP_bbgH,M_KPbb,err_kpbb,chi2)

      ndim = 2
      ncall = ncall1sushi*10
      call vegas(KP_gbbH,M_KPbg,err_kpbg,chi2)

      M_mbb=M_Virt+M_Insbb+M_KPbb
      err_mbb=dsqrt(err_Virt**2+err_Insbb**2+err_KPbb**2)
      
      M_mbg=M_KPbg
      err_mbg=err_KPbg**2

      end if

      ndim=3
      if (dist.eq.1) ndim = 2
      ncall = ncall1sushi*100
      call vegas(MbbgH,M_bbgH,err_bbgH,chi2)

      ndim = 3
      if (dist.eq.1) ndim = 2
      ncall = ncall1sushi*10
      call vegas(MgbbH,M_gbbH,err_gbbH,chi2)

!      write(*,*) 'Mvirt (bb->H)  :   ',M_virt,' +/-',err_virt
!      write(*,*) 'Born*I (bb->gH):   ',M_Insbb,' +/-',err_insbb
!      write(*,*) 'K+P (bb->gH)   :   ',M_KPbb,' +/-',err_kpbb
!      write(*,*) 'K+P (gb->bH)   :   ',M_KPbg,' +/-',err_kpbg
!      write(*,*) 'M^{m} (bb->gH) :   ',M_mbb,' +/-',err_mbb
!      write(*,*) 'M^{m} (gb->bH) :   ',M_mbg,' +/-',err_mbg
!      write(*,*) 'MbbgH (bb->H)  :   ',M_bbgH,' +/-',err_bbgH
!      write(*,*) 'MgbbH (bb->H)  :   ',M_gbbH,' +/-',err_gbbH
!      write(*,*) 'M^{m+1}(bb->gH):   ',M_bbgH,' +/-',err_bbgH
!      write(*,*) 'M^{m+1}(gb->bH):   ',M_gbbH,' +/-',err_gbbH

      Sig_bb=Mborn+M_mbb+M_bbgH
      err_bb=dsqrt(err_born**2+err_mbb**2+err_bbgH**2)

      Sig_bg=M_mbg+M_gbbH
      err_bg=dsqrt(err_mbg**2+err_gbbH**2)

      IF(jets == 1) THEN
      SigLO=0d0
      errLO=0d0
      Sig_bb=M_bbgH
      err_bb=err_bbgH
      Sig_bg=M_gbbH
      err_bg=err_gbbH
      ENDIF

      IF(jets == 2) THEN
      SigLO=0d0
      errLO=0d0
      Sig_bb=0d0
      err_bb=0d0
      Sig_bg=0d0
      err_bg=0d0
      ENDIF

      Sig_total = Sig_bb + Sig_bg
      err_total = err_bb + err_bg

!      write(*,105)
!      write(*,*) "SM differential cross section for bb->h:"
!      write(*,*) "NLO (bb->H)   :    ", Sig_bb," +/-", err_bb
!      write(*,*) "NLO (gb->bH)  :    ", Sig_bg," +/-", err_bg
!      write(*,*) "NLO sum       :    ", Sig_total," +/-", err_total
!      write(*,105)
! 105  format('#--------------------------------------------------#')

      END


      DOUBLE PRECISION FUNCTION Born_bbH(x)
      implicit none
      DOUBLE PRECISION x(10),Pf(4,3),PiA(4,3),PDFsbbbar,Sigma_0,Hbb
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      Born_bbH=0d0      
      
      !write(*,*) "Shad",shad,Hmass,eta1,mT,Hmin,muF

      eta1=x(1)
      eta2=Hmass**2/Shad/eta1
C-----------------***Variablen Reset***---------------------------------
      CALL eta1toP(eta1,PiA,Pf)

      IF(Fjet1(Pf,Hmin)==0) RETURN
      
      Sigma_0=Pi/6d0*Hbb(bbhdifforder)/Hmass**2 !col-,spin-av.,Flux, ME2

      Born_bbH=Sigma_0*PDFsbbbar(eta1,eta2,muF,bbhdifforder)
*** delta(x-1)=Hmass**2/Shad/eta1 ***
      Born_bbH=Born_bbH*Hmass**2/Shad/eta1
* cross section factor
      Born_bbH = Born_bbH * GeV2pb

      END
      
      SUBROUTINE eta1toP(eta1,PiA,Pf)
      DOUBLE PRECISION eta1,eta2,Pf(4,3),PiA(4,3)
      include '../commons/common-bbhdiff-scalesintvars.f'
      eta2=Hmass**2/Shad/eta1
      E = dsqrt(eta1*eta2*Shad)/2d0
** initial state momenta
      PiA(4,1) = E
      PiA(1,1) = 0d0
      PiA(2,1) = 0d0
      PiA(3,1) = E
      PiA(4,2) = E
      PiA(1,2) = 0d0
      PiA(2,2) = 0d0
      PiA(3,2) = -E
** final state momenta
      Pf(4,2)=PiA(4,1)+PiA(4,2)
      Pf(1,2)=PiA(1,1)+PiA(1,2)
      Pf(2,2)=PiA(2,1)+PiA(2,2)
      Pf(3,2)=PiA(3,1)+PiA(3,2)
      END

      DOUBLE PRECISION FUNCTION Mvirt_bbH(x)
      implicit none
      DOUBLE PRECISION x(10),Pf(4,3),PiA(4,3),PDFsbbbar,Sigma_0,Hbb
     f ,alpha_s
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      Mvirt_bbH=0d0
      
      eta1=x(1)
      eta2=Hmass**2/Shad/eta1
C-----------------***Variablen Reset***---------------------------------
      CALL eta1toP(eta1,PiA,Pf)
      IF(Fjet1(Pf,Hmin)==0) RETURN
      
      Sigma_0=Pi/6d0*Hbb(bbhdifforder)/Hmass**2 !col-,spin-av.,Flux, ME2

      Mvirt_bbH=Sigma_0
     f *PDFsbbbar(eta1,eta2,muF,bbhdifforder)*
     f (-4d0/3d0+ 14d0/3d0 *z2- 2d0/3d0 * dlog(muRen**2/Hmass**2)**2)
     f *alpha_s(bbhdifforder)/Pi
*** delta(x-1)=Hmass**2/Shad/eta1 ***
      Mvirt_bbH=Mvirt_bbH*Hmass**2/Shad/eta1
* cross section factor
      Mvirt_bbH = Mvirt_bbH * GeV2pb
      END

      DOUBLE PRECISION FUNCTION Ins_bbgH(x)
      implicit none
      DOUBLE PRECISION x(10),Pf(4,3),PiA(4,3),PDFsbbbar,Sigma_0,Hbb,dot
     f ,s,alpha_s 
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      Ins_bbgH=0d0
      
      eta1=x(1)
      eta2=Hmass**2/Shad/eta1
      
      CALL eta1toP(eta1,PiA,Pf)
      IF(Fjet1(Pf,Hmin)==0) RETURN
      
      s=2d0*dot(PiA,1,PiA,2)
      
      Sigma_0=Pi/6d0*Hbb(bbhdifforder)/Hmass**2 !col-,spin-av.,Flux, ME2

      Ins_bbgH=Sigma_0
     f *PDFsbbbar(eta1,eta2,muF,bbhdifforder) * 
     f ( 20d0/3d0 - 7d0*Pi**2/9d0 + 2d0*dlog(muRen**2/s) 
     f + 2d0/3d0 * dlog(muRen**2/s)**2) * alpha_s(bbhdifforder)/Pi
*** delta(x-1)=Hmass**2/Shad/eta1 ***
      Ins_bbgH=Ins_bbgH*Hmass**2/Shad/eta1
* cross section factor
      Ins_bbgH = Ins_bbgH * GeV2pb
      END

      DOUBLE PRECISION FUNCTION KP_bbgH(x)
      implicit none
      DOUBLE PRECISION x(10),PiA(4,3),Pf(4,3),xran
     f ,KOPbb_bbH,POPbb_bbH,s,Hbb,dot,Sigma_0,alpha_s
      Integer Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      KP_bbgH=0d0
      
      eta1=x(1)
      eta2=Hmass**2/Shad/eta1
      xran=x(2)
      
      CALL eta1toP(eta1,PiA,Pf)

      s=2*dot(PiA,1,PiA,2)

      IF(Fjet1(Pf,Hmin)==0) RETURN
      Sigma_0=Pi/6d0*Hbb(bbhdifforder)/Hmass**2 !col-,spin-av.,Flux, ME2
      
      KP_bbgH=KOPbb_bbH(xran)
      KP_bbgH=KP_bbgH+POPbb_bbH(xran,s)
      KP_bbgH=KP_bbgH*Sigma_0 * alpha_s(bbhdifforder)/Pi 
*** delta(x-1)=Hmass**2/Shad/eta1 ***
      KP_bbgH = KP_bbgH * Hmass**2/Shad/eta1
* cross section factor
      KP_bbgH = KP_bbgH * GeV2pb
      END
      
      DOUBLE PRECISION FUNCTION KOPbb_bbH(x)
      implicit none
      DOUBLE PRECISION x,Kbargq,Ktildegq,theta,PDFsbbbar
     f ,PDFx1,PDFx2,PDF1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      KOPbb_bbH=0d0
      
      PDFx1=PDFsbbbar(eta1/x,eta2,muF,bbhdifforder)
     f  *theta(x,eta1)/x
      PDFx2=PDFsbbbar(eta1,eta2/x,muF,bbhdifforder)
     f  *theta(x,eta2)/x
      PDF1=PDFsbbbar(eta1,eta2,muF,bbhdifforder)
      
      KOPbb_bbH=
     f  (-10d0/3d0+4d0*Pi**2/9d0)*(PDF1+PDF1)

      KOPbb_bbH= KOPbb_bbH +
     f  (2d0/3d0-2d0/3d0*x-2d0/3d0*dlog(1d0-x)-2d0/3d0*x*dlog(1d0-x)
     f   -2d0/3d0*dlog((1d0-x)/x)-2d0/3d0*x*dlog((1d0-x)/x)) *
     f ( PDFx1 + PDFx2 )
     
      KOPbb_bbH= KOPbb_bbH +
     f ( 2d0/3d0 * ((2d0 * dlog(1 - x))/(1 - x)) 
     f + 2d0/3d0 * ((2d0 * dlog((1 - x)/x))/(1 - x))) *
     f ( PDFx1-PDF1 + PDFx2-PDF1)
      END
      
      DOUBLE PRECISION FUNCTION POPbb_bbH(x,s)
      implicit none
      DOUBLE PRECISION x,s,P_qq,Farbe,theta,PDFsbbbar
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      POPbb_bbH=0d0

      POPbb_bbH= 1d0/2d0*P_qq(x)*(-dlog(muF**2/s)) *
     f ( PDFsbbbar(eta1/x,eta2,muF,bbhdifforder)
     f   *theta(x,eta1)/x
     f - PDFsbbbar(eta1,eta2,muF,bbhdifforder)
     f + PDFsbbbar(eta1,eta2/x,muF,bbhdifforder)
     f   *theta(x,eta2)/x
     f - PDFsbbbar(eta1,eta2,muF,bbhdifforder))
      END

      DOUBLE PRECISION FUNCTION KP_gbbH(x)
      implicit none
      DOUBLE PRECISION x(10),PiA(4,3),Pf(4,3),theta,PDFsbg,PDFsgb,xran
     f ,KOPgb_bbH,POPgb_bbH,s,Hbb,dot,Sigma_0,alpha_s
      Integer Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      KP_gbbH=0d0
      
      eta1=x(1)
      eta2=Hmass**2/Shad/eta1
      xran=x(2)
      
      CALL eta1toP(eta1,PiA,Pf)

      s=2*dot(PiA,1,PiA,2)

      IF(Fjet1(Pf,Hmin)==0) RETURN
      Sigma_0=Pi/6d0*Hbb(bbhdifforder)/Hmass**2 !col-,spin-av.,Flux, ME2
      
      KP_gbbH=KOPgb_bbH(xran)
      KP_gbbH=KP_gbbH+POPgb_bbH(xran,s)
      KP_gbbH=KP_gbbH*Sigma_0 * alpha_s(bbhdifforder)/Pi *
     f (PDFsgb(eta1/xran,eta2,muF,bbhdifforder)
     f   *theta(xran,eta1)/xran
     f + PDFsbg(eta1,eta2/xran,muF,bbhdifforder)
     f   *theta(xran,eta2)/xran)
*** delta(x-1)=Hmass**2/Shad/eta1 ***
      KP_gbbH = KP_gbbH * Hmass**2/Shad/eta1
* cross section factor
      KP_gbbH = KP_gbbH * GeV2pb
      END
********************!!!!! gb->bH (bb->H) !!!!!*************************
      DOUBLE PRECISION FUNCTION KOPgb_bbH(x)
      implicit none
      DOUBLE PRECISION x,Kbargq,Ktildegq
      KOPgb_bbH=0d0

      KOPgb_bbH=1d0/2d0 *(Kbargq(x) + Ktildegq(x))
      END
      
      DOUBLE PRECISION FUNCTION POPgb_bbH(x,s)
      implicit none
      DOUBLE PRECISION x,s,Psplitgq,Farbe
      include '../commons/common-bbhdiff-scalesintvars.f'
      POPgb_bbH=0d0

      POPgb_bbH= 1d0/2d0*Psplitgq(x)*(-dlog(muF**2/s))
      END
********************!!!!! gb->bH (bb->H) !!!!**************************

      SUBROUTINE bbhdiffinit(ordin,cme,higgsmass,higgsmin,topmass,
     &jetn,btagn,bnrn,pseudorapin,
     &pTmin,pTmax,yHmin,yHmax,distin,ptfixin,mu_R,mu_F)
      implicit none
      INTEGER jetn,btagn,bnrn,ordin,distin
      LOGICAL pseudorapin
      DOUBLE PRECISION cme,higgsmass,higgsmin,topmass
      DOUBLE PRECISION pTmin,pTmax,yHmin,YHmax
      DOUBLE PRECISION mu_R,mu_F,ptfixin
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-cuts.f'
      include '../commons/common-bbhdiff-btag.f'
      include '../commons/common-bbhdiff-ptdist.f'
      shad = cme**2
      Hmass = higgsmass
      pseudorapflag = pseudorapin
      Hmin=higgsmin
      mT = topmass
      bbhdifforder=ordin
      jets=jetn
      btag=btagn
      bnr=bnrn
      pT_Hmin = pTmin
      pT_Hmax = pTmax
      y_Hmin = yHmin
      y_Hmax = yHmax
      dist = distin
      ptfix = ptfixin
      muF = mu_F
      muRen = mu_R
      !not initialized: parameters of Fjet2!
      END
      
      
************************ dipoles bb-H NLO **************************

      DOUBLE PRECISION FUNCTION D2q1g3_bbH(PiA,Pf)
      implicit none
      DOUBLE PRECISION Pf(4,3),PiA(4,3),Pit(4,3),Pft(4,3)
      DOUBLE PRECISION Pftr(4,3),Pitr(4,3),Pittr(4,3),Pfttr(4,3)
      DOUBLE PRECISION x,dot,s,t,u,Mtree_bbH
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-ptdist.f'
      D2q1g3_bbH=0d0

      s=2*dot(PiA,1,PiA,2)
      t=-2*dot(PiA,1,Pf,1)
      u=-2*dot(PiA,2,Pf,1)

      x=(s+t+u)/s
      CALL xp(Pit,1,2,PiA,x)
      CALL kkk(Pit,Pft,1,2,PiA,1,2,Pf)
      IF(Fjet1(Pft,Hmin)==0) RETURN
      IF(dist.eq.1) RETURN
      
      D2q1g3_bbH= -2*cF/(t*x) *(2d0/(1d0-x)-(1d0+x)) *Mtree_bbH(Pit,Pft)
      END


      DOUBLE PRECISION FUNCTION D1q2g3_bbH(PiA,Pf)
      implicit none
      DOUBLE PRECISION Pf(4,3),PiA(4,3),Pit(4,3),Pft(4,3)
      DOUBLE PRECISION x,dot,s,t,u,Mtree_bbH
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-ptdist.f'
      D1q2g3_bbH=0d0

      s=2*dot(PiA,1,PiA,2)
      t=-2*dot(PiA,1,Pf,1)
      u=-2*dot(PiA,2,Pf,1)

      x=(s+t+u)/s
      CALL xp(Pit,2,1,PiA,x)
      CALL kkk(Pit,Pft,2,1,PiA,1,2,Pf)
      IF(Fjet1(Pft,Hmin)==0) RETURN
      IF(dist.eq.1) RETURN

      D1q2g3_bbH= -2*cF/(u*x) *(2d0/(1d0-x)-(1d0+x)) *Mtree_bbH(Pit,Pft)
      END


      DOUBLE PRECISION FUNCTION D1g2q3_bbH(PiA,Pf)
      implicit none
      DOUBLE PRECISION Pf(4,3),PiA(4,3),Pit(4,3),Pft(4,3)
      DOUBLE PRECISION x,dot,s,t,u,Mtree_bbH
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-ptdist.f'
      D1g2q3_bbH=0d0

      s=2*dot(PiA,1,PiA,2)
      t=-2*dot(PiA,1,Pf,1)
      u=-2*dot(PiA,2,Pf,1)

      x=(s+t+u)/s
      CALL xp(Pit,2,1,PiA,x)
      CALL kkk(Pit,Pft,2,1,PiA,1,2,Pf)
      IF(Fjet1(Pft,Hmin)==0) RETURN
      IF(dist.eq.1) RETURN

      D1g2q3_bbH= -2*tR /(u*x) * (1d0-2*x*(1-x)) * Mtree_bbH(Pit,Pft)
      END
      
      DOUBLE PRECISION FUNCTION D2g1q3_bbH(PiA,Pf)
      implicit none
      DOUBLE PRECISION Pf(4,3),PiA(4,3),Pit(4,3),Pft(4,3)
      DOUBLE PRECISION x,dot,s,t,u,Mtree_bbH
      INTEGER Fjet1
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-bbhdiff-ptdist.f'
      D2g1q3_bbH=0d0

      s=2*dot(PiA,1,PiA,2)
      t=-2*dot(PiA,1,Pf,1)
      u=-2*dot(PiA,2,Pf,1)

      x=(s+t+u)/s
      CALL xp(Pit,1,2,PiA,x)
      CALL kkk(Pit,Pft,1,2,PiA,1,2,Pf)
      IF(Fjet1(Pft,Hmin)==0) RETURN
      IF(dist.eq.1) RETURN

      D2g1q3_bbH= -2*tR /(t*x) * (1d0-2*x*(1-x)) * Mtree_bbH(Pit,Pft)
      END


      DOUBLE PRECISION FUNCTION Mtree_bbH(p12,p34)
      implicit none
      DOUBLE PRECISION p12(4,3),p34(4,3)
      DOUBLE PRECISION s,Mh,dot
      include '../commons/common-consts.f'
      s=2.d0*dot(p12,1,p12,2)
      Mh=dot(p34,2,p34,2)
      Mtree_bbH=2d0*cA*s
      RETURN
      END


      DOUBLE PRECISION FUNCTION MbbgH(x)
      implicit none
      DOUBLE PRECISION x(10),PiA(4,3), Pf(4,3),PDFsbbbar,Mtree
     f ,dot,Hbb,alpha_s,D2q1g3_bbH,D1q2g3_bbH,BORN,x3angle
      Integer Fjet2
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      include '../commons/common-bbhdiff-ptdist.f'
C-----------------***Center of MASS(ET), Particle Mass(XM), (Pi)***-----
      MbbgH=0d0
      eta1=x(1)
      eta2=x(2)
      x3angle=x(3)

C--obtain numbers for a fixed value of pt:
C change from x1,x2,x3 to x1,x3 and fix x2 from pt
      if (dist.eq.1) then
      eta1=x(1)
      x3angle=x(2)
      eta2= 2.d0*dsqrt(ptfix**2*(Hmass**2
     & + ptfix**2 - Hmass**2*x3angle**2))
      eta2 =(Hmass**2 + 2.d0 * ptfix**2 - Hmass**2 * x3angle**2 + eta2)
     & / (Shad * eta1 - Shad  * eta1 * x3angle**2) !formula to get eta2
      end if

C-----------------***Variablen Reset***---------------------------------
      CALL xtoP_12(eta1,eta2,x3angle,PiA,Pf)

      IF(dsqrt(eta1*eta2*Shad)<Hmass) RETURN

      Born=(Mtree(PiA,Pf,cA,cF)-D2q1g3_bbH(PiA,Pf)-D1q2g3_bbH(PiA,Pf))
      IF(Fjet2(Pf,Hmin)==0)
     f Born=(-D2q1g3_bbH(PiA,Pf)-D1q2g3_bbH(PiA,Pf))

      if (dist.eq.1) then
      Born = - Born * ((2.d0 * ptfix
     &* (2.d0 + (Hmass**2 + 2*ptfix**2 - Hmass**2 * x3angle**2)
     &/dsqrt( ptfix**2 *(Hmass**2 + ptfix**2 - Hmass**2 * x3angle**2)))
     &)  /(Shad * eta1 * (-1.d0 + x3angle**2)))
      Born = 2.d0 * Born
      end if

***Flux-Factor***
      Born=Born/(2d0*eta1*eta2*Shad)
***alpha_s-Energieabhaengig***
      Born=Born*4d0*Pi*alpha_s(bbhdifforder)
***PDFs***
      Born=Born*PDFsbbbar(eta1,eta2,muF,bbhdifforder)
***Phase-space-factors***
      BORN=BORN*(1-(Hmass**2/(eta1*eta2*Shad)))/(16d0*Pi)
* Hbb coupling
      BORN = BORN * Hbb(bbhdifforder)
* cross section factor
      BORN = BORN * GeV2pb
* cols and spins
      BORN = BORN /3d0/3d0/2d0/2d0
      MbbgH=Born
      END
      
      DOUBLE PRECISION FUNCTION MgbbH(x)
      implicit none
      DOUBLE PRECISION x(10),PiA(4,3),Pf(4,3),Born,PiTrafo(4,3),dot,Hbb
     f ,PfTrafo(4,3),Mtreebg,PDFsbg,PDFsgb,Mtreegb,alpha_s,D2g1q3_bbH
     f ,D1g2q3_bbH
      double precision x3angle,ptcheck,ptcheck2,calcpt
      Integer Fjet2
      include '../commons/common-bbhdiff-scalesintvars.f'
      include '../commons/common-consts.f'
      include '../commons/common-bbhdiff-order.f'
      include '../commons/common-bbhdiff-fraction.f'
      include '../commons/common-bbhdiff-ptdist.f'
      eta1=x(1)
      eta2=x(2)
      x3angle=x(3)

C--obtain numbers for a fixed value of pt:
C change from x1,x2,x3 to x1,x3 and fix x2 from pt
      if (dist.eq.1) then
      eta1=x(1)
      x3angle=x(2)
      eta2= 2.d0*dsqrt(ptfix**2*(Hmass**2
     & + ptfix**2 - Hmass**2*x3angle**2))
      eta2 =(Hmass**2 + 2.d0 * ptfix**2 - Hmass**2 * x3angle**2 + eta2)
     & / (Shad * eta1 - Shad  * eta1 * x3angle**2) !formula to get eta2
      ptcheck = dsqrt(1-x3angle**2)*(eta1*eta2*Shad-Hmass**2)
     & /(2*dsqrt(eta1*eta2*Shad))
      ptcheck2 = CALCpT(Pf,1)
      if (((ptcheck-ptfix)/ptfix).gt.0.00001d0) write(*,*) 'Error!'
      end if

C-----------------***Center of MASS(ET), Particle Mass(XM), (Pi)***-----
      MgbbH=0d0
      Born=0d0
C-----------------***Variablen Reset***---------------------------------

C      CALL xtoP_1(x,PiA,Pf)
      CALL xtoP_12(eta1,eta2,x3angle,PiA,Pf)

      IF(dsqrt(eta1*eta2*Shad)<Hmass) RETURN

      Born=(Mtreegb(PiA,Pf)/24d0 - D2g1q3_bbH(PiA,Pf)/9d0)
     f                 * PDFsgb(eta1,eta2,muF,bbhdifforder)
      Born=Born+(Mtreebg(PiA,Pf)/24d0 - D1g2q3_bbH(PiA,Pf)/9d0)
     f                 * PDFsbg(eta1,eta2,muF,bbhdifforder)

      IF(Fjet2(Pf,Hmin)==0) THEN
      Born=0d0
      Born=(-D2g1q3_bbH(PiA,Pf)/9d0)
     f                 * PDFsgb(eta1,eta2,muF,bbhdifforder)
      Born=Born+(- D1g2q3_bbH(PiA,Pf)/9d0)
     f                 * PDFsbg(eta1,eta2,muF,bbhdifforder)
      ENDIF

      if (dist.eq.1) then
      Born = - Born * ((2.d0 * ptfix
     &* (2.d0 + (Hmass**2 + 2*ptfix**2 - Hmass**2 * x3angle**2)
     &/dsqrt( ptfix**2 *(Hmass**2 + ptfix**2 - Hmass**2 * x3angle**2)))
     &)  /(Shad * eta1 * (-1.d0 + x3angle**2)))
      Born = 2.d0 * Born
      end if
  
***Flux-Factor***
      Born=Born/(2d0*eta1*eta2*Shad)
***alpha_s-Energieabhaengig***
      Born=Born*4*Pi*alpha_s(bbhdifforder)
***Phase-space-factors***
      BORN=BORN*(1-(Hmass**2/(eta1*eta2*Shad)))/(16d0*Pi)
* Hbb coupling
      BORN = BORN * Hbb(bbhdifforder)
* cross section factor
      BORN = BORN * GeV2pb
* cols and spins
      BORN = BORN /2d0/2d0

      MgbbH=Born
      END


