/*******************************************************************************
 2HDMC - two-Higgs-doublet model calculator
 SM class
 
 http://2hdmc.hepforge.org
*******************************************************************************/

#include "SM.h"
#include <cmath>
#include <iostream>

using namespace std;


SM::SM() {
  m_alpha_s = SM::alpha_s; 

  m_alpha   = SM::alpha;
  m_alpha0  = SM::alpha0;
  m_GF      = SM::GF;
  m_MZ      = SM::MZ;
  m_MW      = SM::MW;
  
  m_GammaW  = SM::GammaW;
  m_GammaZ  = SM::GammaZ;

  m_md_p[0] = SM::md_p;
  m_md_p[1] = SM::ms_p;
  m_md_p[2] = SM::mb_p;
  m_mu_p[0] = SM::mu_p;
  m_mu_p[1] = SM::mc_p;
  m_mu_p[2] = SM::mt_p;  
  m_ml_p[0] = SM::me_p;
  m_ml_p[1] = SM::mmu_p;
  m_ml_p[2] = SM::mtau_p;

  m_CKM[0][0]     = SM::Vud;
  m_CKM[0][1]     = SM::Vus;
  m_CKM[0][2]     = SM::Vub;
  m_CKM[1][0]     = SM::Vcd;
  m_CKM[1][1]     = SM::Vcs;
  m_CKM[1][2]     = SM::Vcb;
  m_CKM[2][0]     = SM::Vtd;
  m_CKM[2][1]     = SM::Vts;
  m_CKM[2][2]     = SM::Vtb;
  
  clear_lookup();
}


gsl_matrix* SM::get_MD() {

  gsl_matrix *m = gsl_matrix_alloc(3,3);

  gsl_matrix_set(m,0,0,m_md_p[0]);
  gsl_matrix_set(m,0,1,0.0);
  gsl_matrix_set(m,0,2,0.0);
  gsl_matrix_set(m,1,0,0.0);
  gsl_matrix_set(m,1,1,m_md_p[1]);
  gsl_matrix_set(m,1,2,0.0);
  gsl_matrix_set(m,2,0,0.0);
  gsl_matrix_set(m,2,1,0.0);
  gsl_matrix_set(m,2,2,m_md_p[2]);

  return m;
}


gsl_matrix* SM::get_MU() {

  gsl_matrix *m = gsl_matrix_alloc(3,3);

  gsl_matrix_set(m,0,0,m_mu_p[0]);
  gsl_matrix_set(m,0,1,0.0);
  gsl_matrix_set(m,0,2,0.0);
  gsl_matrix_set(m,1,0,0.0);
  gsl_matrix_set(m,1,1,m_mu_p[1]);
  gsl_matrix_set(m,1,2,0.0);
  gsl_matrix_set(m,2,0,0.0);
  gsl_matrix_set(m,2,1,0.0);
  gsl_matrix_set(m,2,2,m_mu_p[2]);

  return m;
}


gsl_matrix* SM::get_ML() {

  gsl_matrix *m = gsl_matrix_alloc(3,3);

  gsl_matrix_set(m,0,0,m_ml_p[0]);
  gsl_matrix_set(m,0,1,0.0);
  gsl_matrix_set(m,0,2,0.0);
  gsl_matrix_set(m,1,0,0.0);
  gsl_matrix_set(m,1,1,m_ml_p[1]);
  gsl_matrix_set(m,1,2,0.0);
  gsl_matrix_set(m,2,0,0.0);
  gsl_matrix_set(m,2,1,0.0);
  gsl_matrix_set(m,2,2,m_ml_p[2]);

  return m;
}


double SM::get_CKM_element(int i, int j) {

  if ((i>0)&&(i<=3)&&(j>0)&&(j<=3))
    return m_CKM[i-1][j-1];
  else
    return 0.;
}


gsl_matrix* SM::get_CKM_matrix() {
  gsl_matrix *ckm = gsl_matrix_alloc(3,3);
  for (int i=0;i<3;i++) {
    for (int j=0; j<3;j++) {
      gsl_matrix_set(ckm,i,j,m_CKM[i][j]);
    }
  }
  return ckm;
}


int SM::get_Nactivef(double M) {
  int nf = 0;
  for (int i=1;i<=6;i++) {
    if (get_qmass_pole(i)<M) nf++;
  }

  return nf;

}


double SM::get_v() {
  return 1./sqrt(sqrt(2)*m_GF);
}


double SM::get_v2() {
  return 1./(sqrt(2)*m_GF);
}


double SM::get_MW() {
  return m_MW;
}

  
double SM::get_MZ() {
  return m_MZ;
}


double SM::get_GF() {
  return m_GF;
}
  

double SM::get_gamma_W() {
  return m_GammaW;
}


double SM::get_gamma_Z() {
  return m_GammaZ;
}


double SM::get_gamma_V(int v) {
  if (v==2) return m_GammaZ;
  if (v==3) return m_GammaW;
  return 0.;  
}


double SM::get_gamma_top() {
  
  double tWb = get_gamma_tWd(3);
  double tWs = get_gamma_tWd(2);
  double tWd = get_gamma_tWd(1);
  
  double gtot = tWb+tWs+tWd;
  
/*  
  
  double Vtb = get_CKM_element(3,3);
  double GF = get_GF();
  double mt = get_umass_pole(3);
  double mb = get_dmass_pole(3);
  double mw = get_MW();

  double x = mw/mt;
  double y = mb/mt;

  double l = (1.-pow(x+y,2))*(1.-pow(x-y,2));
  double I = ((1.-pow(x,2))*(1.+2.*pow(x,2))-pow(y,2)*(2.-pow(x,2)-pow(y,2)))*sqrt(l);

  double gtot = GF*pow(mt,3)/(8.*M_PI*sqrt(2))*pow(Vtb,2)*I;

  printf("gtot top: %16.8E %16.8E %16.8E %16.8E\n", gtot, get_gamma_tWd(3), get_gamma_tWd(2), get_gamma_tWd(1));
*/
  return gtot;
}

double SM::get_gamma_tWd(int d) {

  if ((d<1)||(d>3)) return 0.;

  double VCKM = get_CKM_element(3,d);
  double GF = get_GF();
  double mt = get_umass_pole(3);
  double md = get_dmass_pole(d);
  double mb = get_dmass_pole(3);
  double mw = get_MW();

  double x = mw/mt;
  double y = md/mt;

  double l = (1.-pow(x+y,2))*(1.-pow(x-y,2));
  double I = ((1.-pow(x,2))*(1.+2.*pow(x,2))-pow(y,2)*(2.-pow(x,2)-pow(y,2)))*sqrt(l);

  double gtot = GF*pow(mt,3)/(8.*M_PI*sqrt(2))*pow(VCKM,2)*I;

  double as = run_alphas_MSbar(mt,mt,mb);
  double CF = 4./3.;
  double K = 1.+CF*as/M_PI*(5./4.-M_PI*M_PI/3.);

  gtot = gtot*K;

  return gtot;

}


double SM::get_vmass(int v) {
  if (v==2) return get_MZ();
  if (v==3) return get_MW();
  return 0.;
}


double SM::get_g() {
  double v = get_v();

  return 2.*get_MW()/v;
}


double SM::get_gprime() {
  double tantW = get_sintw()/get_costw();
  double g = get_g();

  double gp = g*tantW;

  return gp;
}


double SM::get_e() {
  double e = get_g()*get_sintw();
  return e;
}


double SM::get_costw() {
  double s = get_sintw();
  double c = sqrt(1.-s*s);
  return c;
}


double SM::get_sintw() {
  double s = sqrt(1.-get_MW()*get_MW()/(m_MZ*m_MZ));
  return s;
}


double SM::get_alpha() {

  return m_alpha;
}


double SM::get_alpha0() {

  return m_alpha0;
}



double SM::get_alpha_s() {
  return m_alpha_s;
}


void SM::set_diagonal_CKM() {
  m_CKM[0][0]     = 1.0;
  m_CKM[0][1]     = 0.0;
  m_CKM[0][2]     = 0.0;
  m_CKM[1][0]     = 0.0;
  m_CKM[1][1]     = 1.0;
  m_CKM[1][2]     = 0.0;
  m_CKM[2][0]     = 0.0;
  m_CKM[2][1]     = 0.0;
  m_CKM[2][2]     = 1.0;
  clear_lookup();
}


double SM::get_qmass_pole(int flav) {
  if ((flav>=1)&&(flav<=6)) {
    if ((flav%2)==1) {
      return m_md_p[(flav-1)/2];
    } else {
      return m_mu_p[flav/2-1];
    }
  } else {
    return 0.;
  }
}


double SM::get_lmass_pole(int l) {
  if ((l>=1)&&(l<=3))
    return m_ml_p[l-1];
  else
    return 0.;
}


double SM::get_dmass_pole(int dq) {
  if ((dq>=1)&&(dq<=3))
    return m_md_p[dq-1];
  else
    return 0.;
}


double SM::get_dmass_MSbar(int dq) {
  if ((dq>=1)&&(dq<=3))
    return get_qmass_MSbar(2*dq-1);
  else
    return 0.;
}


double SM::get_umass_pole(int dq) {
  if ((dq>=1)&&(dq<=3))
    return m_mu_p[dq-1];
  else
    return 0.;
}


double SM::get_umass_MSbar(int dq) {
  if ((dq>=1)&&(dq<=3))
    return get_qmass_MSbar(2*dq);
  else
    return 0.;
}


double SM::get_qmass_MSbar(int flav) {
  // Evaluation of running mass from pole mass
  // arXiv:0712.1419

  int loops = 1;

  if ((flav>6)||(flav<1)) return 0;

  if ((!b_HD)&&(flav<=3)) return get_qmass_pole(flav);  
  if ((b_HD) &&(flav<=2)) return get_qmass_pole(flav);  
  // Lookup table
  if (m_qmass_ms[flav]>=0.) return m_qmass_ms[flav];

  double Kq2 = 0.0;
  double Kq3 = 0.0;
  double mq_p = get_qmass_pole(flav);
  double mb = get_qmass_pole(5);
  double mt = get_qmass_pole(6);

  if (flav==4) {Kq2 = 11.21; Kq3=123.8;}
  if (flav==5) {Kq2 = 10.17; Kq3=101.5;}
  if (flav==6) {Kq2 = 9.13;  Kq3=80.4;}

  double as = run_alphas_MSbar(mq_p, get_qmass_pole(6), get_qmass_pole(5))/M_PI;

  // 1-loop conversion
  double mq_msbar = mq_p/(1.+4./3.*as);
  
   // Change to higher loops
   if (loops==2) {
     mq_msbar = mq_p/(1.+4./3.*as+Kq2*pow(as,2));
   } else if(loops==3) {
     mq_msbar = mq_p/(1.+4./3.*as+Kq2*pow(as,2)+Kq3*pow(as,3));
   }
 
  // Run to mq(mq)
  mq_msbar=run_qmass_MSbar(mq_msbar,mq_p,mq_msbar,mt,mb);


  //  HD comparison
  if ((b_HD)&&(flav==3)) {mq_msbar = ms_5;}
  if ((b_HD)&&(flav==4)) {mq_msbar = mc_5;}
  if ((b_HD)&&(flav==5)) {mq_msbar = mb_5;}
  if ((b_HD)&&(flav==6)) {mq_msbar = mt_5;}

  m_qmass_ms[flav] = mq_msbar;
  
  return mq_msbar;
}

void SM::set_qmass_msbar(int flav, double qmass_in) {
 
    double mbar = get_qmass_MSbar(flav);
    double mpole = get_qmass_pole(flav);
 
    while (abs(qmass_in-mbar)/qmass_in > 1E-5) {    
      set_qmass_pole(flav,mpole+(qmass_in-mbar));
      mbar = get_qmass_MSbar(flav);
      mpole = get_qmass_pole(flav);
    }    
}


void SM::set_alpha(double alpha_in) {
  
  if (alpha_in>=0) {
    m_alpha=alpha_in;
    clear_lookup();
  }
}


void SM::set_alpha0(double alpha_in) {
  
  if (alpha_in>=0) {
    m_alpha0=alpha_in;
    clear_lookup();
  }
}


void SM::set_alpha_s(double alpha_s_in) {
  if (alpha_s_in>=0) {
    m_alpha_s=alpha_s_in;
    clear_lookup();
  }
}


void SM::set_GF(double GF_in) {
  if (GF_in>=0) {
    m_GF=GF_in;
    clear_lookup();
  }
    
}


void SM::set_MZ(double MZ_in) {
  if (MZ_in>=0) {
    m_MZ=MZ_in;
    clear_lookup();
  }
}


void SM::set_MW(double MW_in) {
  if (MW_in>=0) {
    m_MW=MW_in;
    clear_lookup();
  }
}


void SM::set_gamma_Z(double GammaZ_in) {
  if (GammaZ_in>=0) {
    m_GammaZ=GammaZ_in;
    clear_lookup();
  }
}


void SM::set_gamma_W(double GammaW_in) {
  if (GammaW_in>=0) {
    m_GammaW=GammaW_in;
    clear_lookup();
  }
}


void SM::set_lmass_pole(int l, double lmass_in) {
  if ((l>=1)&&(l<=3)&&(lmass_in>=0)) {
    m_ml_p[l-1]=lmass_in;
    clear_lookup();
  }
}


void SM::set_qmass_pole(int flav, double qmass_in) {
  if ((flav>=1)&&(flav<=6)&&(qmass_in>=0)) {
    if ((flav%2)==1) {
      m_md_p[(flav-1)/2]=qmass_in;
    } else {
      m_mu_p[flav/2-1]=qmass_in;
    }
    clear_lookup();
  }
}


void SM::set_dmass_pole(int dq, double dmass_in) {
  if ((dq>=1)&&(dq<=3)&&(dmass_in>=0)) {
    m_md_p[dq-1]=dmass_in;
    clear_lookup();
  }
}


void SM::set_umass_pole(int uq, double umass_in) {
  if ((uq>=1)&&(uq<=3)&&(umass_in>=0)) {
    m_mu_p[uq-1]=umass_in;
    clear_lookup();
  }
}


// Quark mass evolution function R.
// Adapted from arXiv:0712.1419
double SM::RQ(double as, int nf, int nloops=3) {

  const double zeta3=1.202056903;
  const double zeta4=1.082323234;
  const double zeta5=1.036927755;
  
  double beta[4] = {0.,0.,0.,0.};
  double gamma[4] = {0.,0.,0.,0.};

  if (nloops>4) nloops = 3;

  if (nloops>=1) {
    beta[0]=1./4.*(11.-2./3.*nf);
    gamma[0]=1.;
  }
  if (nloops>=2) {
    beta[1] = 1./16.*(102.-38./3.*nf);
    gamma[1] = 1./16.*(202./3.-20./9.*nf);
  }
  if (nloops>=3) {
    beta[2] = 1./64.*(2857./2.-5033./18.*nf+325./54.*nf*nf);
    gamma[2] = 1./64.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
  } 
  if (nloops>=4) {
    beta[3] = 1./256.*(149753./6.+3564.*zeta3-(1078361./162.+6508./27.*zeta3)*nf+(50065./162.+6472./81.*zeta3)*nf*nf+1093./729.*nf*nf*nf);
    gamma[3]=1./256.*(4603055./162.+135680./27.*zeta3-8800.*zeta5-(91723./27.+34192./9.-880.*zeta4-18400./9.*zeta5)*nf+(5242./243.+800./9.*zeta3-160./3.*zeta4)*nf*nf-(332./243.-64./27.*zeta3)*nf*nf);
  }

  double ras = as/M_PI;

  double C1 = gamma[1]/beta[0]-beta[1]*gamma[0]/pow(beta[0],2);
  double C2 = gamma[2]/beta[0]-beta[1]*gamma[1]/pow(beta[0],2)-beta[2]*gamma[0]/pow(beta[0],2)+pow(beta[1],2)*gamma[0]/pow(beta[0],3);
//  double C3 = gamma[3]/beta[0]-beta[1]*gamma[2]/pow(beta[0],2)+pow(beta[1],2)*gamma[1]/pow(beta[0],3)-beta[2]*gamma[1]/pow(beta[0],2)-pow(beta[1],3)*gamma[0]/pow(beta[0],4)+2.*beta[1]*beta[2]*gamma[0]/pow(beta[0],3)-beta[3]*gamma[0]/pow(beta[0],2);   


  double R = pow(ras,gamma[0]/beta[0])*(1.+ras*C1+0.5*pow(ras,2)*(pow(C1,2)+C2));
//  double R = pow(ras,gamma[0]/beta[0])*(1.+ras*C1+pow(ras,2)*(pow(C1,2)+C2)+pow(ras,3)*(1./6.*pow(C1,3)+1./2.*C1*C2+1./3.*C3));

  return R;
}

// Matching function for the quark mass evolution valid at threshold,
// i.e. when the matching scale for the mass evolution coincides with the 
// scale at which the number of flavors change in the effective theory.
// Adopted from arXiv:0712.1419
double SM::m_threshold(double as) {
  //  const double zeta3  = 1.202056903;
  //  const double zeta4  = 1.082323234;
  //  const double B4     = -1.762800; 

  double MR = 1.;

//  if (!b_HD) {
//    double ras = as/M_PI;
//    MR = 1.+1./12.*pow(ras,2)*89./36.;
//  }
  
  return MR;

}


double SM::run_qmass_MSbar(double quark_mass, double Qinit, double Qfin, double mtop, double mbot) {
  // Running of quark masses

  double R_Qinit,R_Qfin,running_mass;

  double alphas_Qinit=run_alphas_MSbar(Qinit,mtop,mbot);
  double alphas_Qfin =run_alphas_MSbar(Qfin,mtop,mbot);
  double mcharm = get_qmass_pole(4);
  
  // Are we asking for the impossible?
  if ((alphas_Qinit!=alphas_Qinit)||(alphas_Qfin!=alphas_Qfin)) return quark_mass;

  int nf = 0;

  if (Qinit < mcharm) {
    nf = 3;    
    R_Qinit = RQ(alphas_Qinit,nf);

    if (Qfin <= mcharm) {
      R_Qfin = RQ(alphas_Qfin,nf);
      return R_Qfin/R_Qinit*quark_mass;
    } else {
      alphas_Qfin = run_alphas_MSbar(mcharm,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = 1./m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mcharm,Qfin,mtop,mbot);
    }
  } else if ((Qinit>=mcharm)&&(Qinit<mbot)) {
    nf = 4;
    R_Qinit = RQ(alphas_Qinit,nf);
    
    if ((Qfin >= mcharm)&&(Qfin<mbot)) {
      R_Qfin = RQ(alphas_Qfin,nf);
      return R_Qfin/R_Qinit*quark_mass;
    } else if (Qfin >= mbot) {
      alphas_Qfin = run_alphas_MSbar(mbot,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = 1./m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mbot,Qfin,mtop,mbot);
    } else if (Qfin < mcharm) {
      alphas_Qfin = run_alphas_MSbar(mcharm,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mcharm-1.E-12,Qfin,mtop,mbot);
    }
 } else if ((Qinit>=mbot)&&(Qinit<mtop)) {
    nf = 5;
    R_Qinit = RQ(alphas_Qinit,nf);
    
    if ((Qfin >= mbot)&&(Qfin<mtop)) {
      R_Qfin = RQ(alphas_Qfin,nf);
      return R_Qfin/R_Qinit*quark_mass;
    } else if (Qfin >= mtop) {
      alphas_Qfin = run_alphas_MSbar(mtop,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = 1./m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mtop,Qfin,mtop,mbot);
    } else if (Qfin < mbot) {
      alphas_Qfin = run_alphas_MSbar(mbot,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mbot-1.E-12,Qfin,mtop,mbot);
    }
  } else if (Qinit >= mtop) {
    nf = 6;
    R_Qinit = RQ(alphas_Qinit,nf);
    
    if (Qfin >= mtop) {
      R_Qfin = RQ(alphas_Qfin,nf);
      return R_Qfin/R_Qinit*quark_mass;
    } else {
      alphas_Qfin = run_alphas_MSbar(mtop,mtop,mbot);
      R_Qfin = RQ(alphas_Qfin,nf);
      running_mass = R_Qfin/R_Qinit*quark_mass;
      double K = m_threshold(alphas_Qfin);
      return run_qmass_MSbar(K*running_mass,mtop-1.E-12,Qfin,mtop,mbot);
    }
  }

  return 0.;
}
      

double SM::QCD_beta(int c, int nf) {

  if ((c < 0)||(c > 2)) return 0.;
  
  double beta = -1;
  
  switch(c) {
  
    case 0:
     beta=11.-2.*nf/3.;
     break;
    case 1:
     beta=51.-19.*nf/3.;
     break;
    case 2:
     beta=2857.-5033.*nf/9.+325.*nf*nf/27.;
     break;
  }
  
  return beta;
}

// Running QCD coupling constant alpha_s at scale Q
double SM::run_alphas_MSbar(double Q, double mtop, double mbot) {

  double MZ 		= m_MZ;
  double alphas_MZ 	= m_alpha_s;

// Only five flavours is used for the running of alpha_s in agreement with HD
  double beta0, beta1, beta2, alpha_s_running, Lambda4, Lambda5,  Lambda_min, Lambda_max, Lambda_moy, alphas_min, alphas_moy;
// double beta0, beta1, beta2, alpha_s_running, Lambda4, Lambda5, Lambda6, Lambda_min, Lambda_max, Lambda_moy, alphas_min, alphas_moy, alpha_s_match;
  int nf;
  
  double pi = M_PI;

  nf=5;
  beta0 = QCD_beta(0,nf);
  beta1 = QCD_beta(1,nf);
  beta2 = QCD_beta(2,nf);
  
  Lambda_min=1.e-3;
  Lambda_max=1.;
  alphas_min=0.;

  while(fabs(1.-alphas_min/alphas_MZ)>=1.e-6) {
    alphas_min=4.*pi/beta0/log(pow(MZ/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_min,2.)))/log(pow(MZ/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_min,2.)),2.)*(pow(log(log(pow(MZ/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
    Lambda_moy=(Lambda_min+Lambda_max)/2.;
    alphas_moy=4.*pi/beta0/log(pow(MZ/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_moy,2.)))/log(pow(MZ/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_moy,2.)),2.)*(pow(log(log(pow(MZ/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
    if((alphas_MZ>=alphas_min)&&(alphas_MZ<=alphas_moy))
      Lambda_max=Lambda_moy;
    else Lambda_min=Lambda_moy;
  }

  Lambda5=Lambda_min;

// only use five flavours in alpha_s running for agreement with Hdecay
//  if((Q<=mtop)&&(Q>=mbot))
  if(Q>=mbot) {
    /* 5 active flavors */
      alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda5,2.)))/log(pow(Q/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda5,2.)),2.)*(pow(log(log(pow(Q/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
      return alpha_s_running;

/*** These aren't the flavours you are looking for - please move along ***
  } if((Q>mtop)) {
    // 6 active flavors 
    alpha_s_match=4.*pi/beta0/log(pow(mtop/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda5,2.)))/log(pow(mtop/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda5,2.)),2.)*(pow(log(log(pow(mtop/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
    nf=6;
    beta0 = QCD_beta(0,nf);
    beta1 = QCD_beta(1,nf);
    beta2 = QCD_beta(2,nf);
  
    Lambda_min=1.e-3;
    Lambda_max=1.;
    alphas_min=0.;
  
    while(fabs(1.-alphas_min/alpha_s_match)>=1.e-5) {
      alphas_min=4.*pi/beta0/log(pow(mtop/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_min,2.)))/log(pow(mtop/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_min,2.)),2.)*(pow(log(log(pow(mtop/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
    
	  Lambda_moy=(Lambda_min+Lambda_max)/2.;
      alphas_moy=4.*pi/beta0/log(pow(mtop/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_moy,2.)))/log(pow(mtop/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_moy,2.)),2.)*(pow(log(log(pow(mtop/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
    
      if((alpha_s_match>=alphas_min)&&(alpha_s_match<=alphas_moy))
        Lambda_max=Lambda_moy;
	  else Lambda_min=Lambda_moy;
	}

    Lambda6=Lambda_min;

    alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda6,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda6,2.)))/log(pow(Q/Lambda6,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda6,2.)),2.)*(pow(log(log(pow(Q/Lambda6,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
    return alpha_s_running;
*/
  } else {
    /* 4 active flavors */
  
    alpha_s_running=4.*pi/beta0/log(pow(mbot/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda5,2.)))/log(pow(mbot/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda5,2.)),2.)*(pow(log(log(pow(mbot/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

//  1-loop ratio between pole mass and MSbar mass at pole mass
    double mratio=1+4./3.*alpha_s_running/pi;   
    
    double alpha_s_msbar=4.*pi/beta0/log(pow(mbot/mratio/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/mratio/Lambda5,2.)))/log(pow(mbot/mratio/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/mratio/Lambda5,2.)),2.)*(pow(log(log(pow(mbot/mratio/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

    double R_Qinit = RQ(alpha_s_running,4);
    double R_Qfin = RQ(alpha_s_msbar,4);
    mratio = mratio/R_Qfin*R_Qinit;
    
    double xi_factor = 1.-alpha_s_running/pi/6.*log(pow(mratio,2))
    +pow(alpha_s_running/pi,2)*(11./12.-11./24.*log(pow(mratio,2))+1./36.*pow(log(pow(mratio,2)),2));

    double alpha_s_match = alpha_s_running*pow(xi_factor,2);

    nf=4;
    beta0 = QCD_beta(0,nf);
    beta1 = QCD_beta(1,nf);
    beta2 = QCD_beta(2,nf);

    Lambda_min=1.e-3;
    Lambda_max=1.;
    alphas_min=0.;

    while(fabs(1.-alphas_min/alpha_s_match)>=1.e-5) {
	  alphas_min=4.*pi/beta0/log(pow(mbot/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_min,2.)))/log(pow(mbot/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_min,2.)),2.)*(pow(log(log(pow(mbot/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
   	  Lambda_moy=(Lambda_min+Lambda_max)/2.;
      alphas_moy=4.*pi/beta0/log(pow(mbot/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_moy,2.)))/log(pow(mbot/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_moy,2.)),2.)*(pow(log(log(pow(mbot/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
	  if((alpha_s_match>=alphas_min)&&(alpha_s_match<=alphas_moy))
	     Lambda_max=Lambda_moy;
	  else Lambda_min=Lambda_moy;
    }
    
    Lambda4=Lambda_min;

    alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda4,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda4,2.)))/log(pow(Q/Lambda4,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda4,2.)),2.)*(pow(log(log(pow(Q/Lambda4,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

    if (Q<=get_qmass_pole(4)) {
       /* 3 active flavors */
      double mcharm = get_qmass_pole(4);
  
      alpha_s_running=4.*pi/beta0/log(pow(mcharm/Lambda4,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mcharm/Lambda4,2.)))/log(pow(mcharm/Lambda4,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mcharm/Lambda4,2.)),2.)*(pow(log(log(pow(mcharm/Lambda4,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

//  1-loop ratio between pole mass and MSbar mass at pole mass
      double mratio=1+4./3.*alpha_s_running/pi;   
    
      double alpha_s_msbar=4.*pi/beta0/log(pow(mcharm/mratio/Lambda4,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mcharm/mratio/Lambda4,2.)))/log(pow(mcharm/mratio/Lambda4,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mcharm/mratio/Lambda4,2.)),2.)*(pow(log(log(pow(mcharm/mratio/Lambda4,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

//  run MSbar mass to itself i.e. (m(m(M))
      double R_Qinit = RQ(alpha_s_running,3);
      double R_Qfin = RQ(alpha_s_msbar,3);
      mratio = mratio/R_Qfin*R_Qinit;
    
      double xi_factor = 1.-alpha_s_running/pi/6.*log(pow(mratio,2))+pow(alpha_s_running/pi,2)*(11./12.-11./24.*log(pow(mratio,2))+1./36.*pow(log(pow(mratio,2)),2));

      double alpha_s_match = alpha_s_running*pow(xi_factor,2);

      nf=3;
      beta0 = QCD_beta(0,nf);
      beta1 = QCD_beta(1,nf);
      beta2 = QCD_beta(2,nf);

      Lambda_min=1.e-3;
      Lambda_max=1.;
      alphas_min=0.;

      while(fabs(1.-alphas_min/alpha_s_match)>=1.e-5) {
   	    alphas_min=4.*pi/beta0/log(pow(mcharm/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mcharm/Lambda_min,2.)))/log(pow(mcharm/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mcharm/Lambda_min,2.)),2.)*(pow(log(log(pow(mcharm/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
	    Lambda_moy=(Lambda_min+Lambda_max)/2.;
	    alphas_moy=4.*pi/beta0/log(pow(mcharm/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mcharm/Lambda_moy,2.)))/log(pow(mcharm/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mcharm/Lambda_moy,2.)),2.)*(pow(log(log(pow(mcharm/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
  
	    if((alpha_s_match>=alphas_min)&&(alpha_s_match<=alphas_moy))
	      Lambda_max=Lambda_moy;
	    else Lambda_min=Lambda_moy;
      }
    
      double Lambda3=Lambda_min;
      alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda3,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda3,2.)))/log(pow(Q/Lambda3,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda3,2.)),2.)*(pow(log(log(pow(Q/Lambda3,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
       
       
    }
    
    return alpha_s_running;
  }
}

void SM::clear_lookup() {

  for (int i=0;i<=6;i++) {
    m_qmass_ms[i] = -1.;
  }

  m_qmass_ms[4] = get_qmass_MSbar(4);
  m_qmass_ms[5] = get_qmass_MSbar(5);
  m_qmass_ms[6] = get_qmass_MSbar(6);

}
