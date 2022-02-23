/*******************************************************************************
 2HDMC - two-Higgs-doublet model calculator
 THDM class

 http://2hdmc.hepforge.org
*******************************************************************************/

#include "THDM.h"
#include "DecayTable.h"
#include "Constraints.h"
#include "HBHS.h"
#include "Util.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>

#define npara 15
#define nsmpara 6
#define nres 15


bool THDM::first_run = true;
const char *THDM::version = "1.8.0";

using namespace std;

/***************************************************************************
* Method used by external codes (SusHi) to set the model parameters
*
* key = 1: Generic basis  (lam1,lam2,lam3,lam4,lam5,lam6,lam7,m12_2,tb)
* key = 2: Physical basis (mh,mH,mA,mC,sba,l6,l7,m12_2,tb)
* key = 3: Hybrid basis with sba (mh, mH, sba, Z4, Z5, Z7, tb)
* key = 4: Hybrid basis with cba (mh, mH, cba, Z4, Z5, Z7, tb)
*
***************************************************************************/
int thdmc_set_param(int key, double smpara[nsmpara], double para[npara], double res[nres], int slha) {

/*
   cout << "Key="<< key << endl;
   for (int i=0; i<nsmpara; i++) {
     cout << i << " " << smpara[i] << endl;
  }
   for (int i=0; i<npara; i++) {
     cout << i << " " << para[i] << endl;
  }
*/

   for (int i=0; i<nres; i++) {
     res[i] = 0.;
  }

   THDM *model = new THDM();

   SM sm = model->get_SM();
   sm.set_alpha(1./smpara[0]);
   sm.set_GF(smpara[1]);
   sm.set_alpha_s(smpara[2]);
   sm.set_MZ(smpara[3]);
   sm.set_qmass_msbar(5,smpara[4]);
   sm.set_qmass_pole(6,smpara[5]);

   model->set_SM(sm);

   int type = 0;

 if (key==1) {

   double lam1 = para[0];
   double lam2 = para[1];
   double lam3 = para[2];
   double lam4 = para[3];
   double lam5 = para[4];
   double lam6 = para[5];
   double lam7 = para[6];
   double m12_2 = para[7];
   double tb = para[8];
   type = (int)para[9];

   bool ok = model->set_param_gen(lam1,lam2,lam3,lam4,lam5,lam6,lam7,m12_2,tb);

   if (!ok) return -1;
   model->set_yukawas_type(type);
 } else if (key==2) {

   double mh = para[0];
   double mH = para[1];
   double mA = para[2];
   double mC = para[3];
   double sba = para[4];
   double lam6 = para[5];
   double lam7 = para[6];
   double m12_2 = para[7];
   double tb = para[8];
   type = (int)para[9];

   bool ok = model->set_param_phys(mh,mH,mA,mC,sba,lam6,lam7,m12_2,tb);
   if (!ok) return -1;
   model->set_yukawas_type(type);

 } else if (key==3) {

   double mh = para[0];
   double mH = para[1];
   double sba = para[2];
   double tb = para[3];
   double Z4 = para[4];
   double Z5 = para[5];
   double Z7 = para[6];
   type = (int)para[9];

   bool ok = model->set_param_hybrid_sba(mh,mH,sba,Z4,Z5,Z7,tb);
   if (!ok) return -1;
   model->set_yukawas_type(type);

 } else if (key==4) {

   double mh = para[0];
   double mH = para[1];
   double cba = para[2];
   double tb = para[3];
   double Z4 = para[4];
   double Z5 = para[5];
   double Z7 = para[6];
   type = (int)para[9];

   bool ok = model->set_param_hybrid(mh,mH,cba,Z4,Z5,Z7,tb);
   if (!ok) return -1;
   model->set_yukawas_type(type);

 } else {
  return -1;
 }

//   model->print_param_gen();
//   model->print_param_phys();

   double mh,mH,mA,mC,sba,l6,l7,m12_2,tb;

   model->get_param_phys(mh,mH,mA,mC,sba,l6,l7,m12_2,tb);

   bool P = model->check_perturbativity(4.*M_PI);
   bool U = model->check_unitarity(16.*M_PI);
   bool S = model->check_stability();

    DecayTable table(*model);

    double gtoth = table.get_gammatot_h(1);
    double gtotH = table.get_gammatot_h(2);
    double gtotA = table.get_gammatot_h(3);

    res[0] = mh;
    res[1] = mH;
    res[2] = mA;
    res[3] = -asin(sba)+atan(tb);
    res[4] = tb;
    res[5] = (double)type;
    res[6] = S ? 1.0 : 0.0;
    res[7] = P ? 1.0 : 0.0;
    res[8] = U ? 1.0 : 0.0;
    res[9] = 0.;
    res[10] = gtoth;
    res[11] = gtotH;
    res[12] = gtotA;

    if (slha) {
     model->write_LesHouches("2HDMC.out",1,0,1);
    }

    if ((mh<0)||(mH<0)||(mA<0)||(mC<0)) return -1;

 return 0;

}

THDM::THDM() {
  init();
}


void THDM::init() {
  if (first_run) {
  	print_info();
  }

  params_set = false;
  yukawas_type = -1;

  v2 = sm.get_v2();

  for (int i=1;i<8;i++) {
    lambda[i]=0;
  }

  beta	=	0;
  m22_2	=	0;

  kappa_D = gsl_matrix_alloc(3,3);
  kappa_U = gsl_matrix_alloc(3,3);
  kappa_L = gsl_matrix_alloc(3,3);
  rho_D   = gsl_matrix_alloc(3,3);
  rho_U   = gsl_matrix_alloc(3,3);
  rho_L   = gsl_matrix_alloc(3,3);
  rho_N   = gsl_matrix_alloc(3,3);

  gsl_matrix_set_zero(kappa_D);
  gsl_matrix_set_zero(kappa_U);
  gsl_matrix_set_zero(kappa_L);
  gsl_matrix_set_zero(rho_D);
  gsl_matrix_set_zero(rho_U);
  gsl_matrix_set_zero(rho_L);
  gsl_matrix_set_zero(rho_N);
}


bool THDM::set_param_gen(double lambda1, double lambda2, double lambda3,
			 double lambda4, double lambda5, double lambda6,
			 double lambda7, double m12_2, double tan_beta) {

  // tan(beta) must be greater than 0 for valid basis
  if (tan_beta<=0) {
    params_set = false;
    return params_set;
  }

  double v2  = sm.get_v2();

  lambda[1] = lambda1;
  lambda[2] = lambda2;
  lambda[3] = lambda3;
  lambda[4] = lambda4;
  lambda[5] = lambda5;
  lambda[6] = lambda6;
  lambda[7] = lambda7;
  beta = atan(tan_beta);

  double tb  = tan_beta;
  double ctb = 1./tb;
  double cb  = 1./sqrt(1.+tb*tb);
  double sb  = tb*cb;
  double sb2 = sb*sb;
  double cb2 = cb*cb;
  double s2b = 2.*sb*cb;
  double c2b = cb2-sb2;
  double s3b = 3.*sb-4.*sb*sb2;
  double c3b = 4.*cb*cb2-3.*cb;

  m22_2 = m12_2*ctb-0.5*v2*(lambda[2]*sb2+(lambda[3]+lambda[4]+lambda[5])*cb2+lambda[6]*cb2*ctb+3.*lambda[7]*sb*cb);

  // Mass of CP-odd Higgs
  double mA_2 = m12_2/(sb*cb)-0.5*v2*(2.*lambda5+lambda6*ctb+lambda7*tb);

  double lambda345 = lambda3+lambda4+lambda5;
  double lhat = 0.5*s2b*(lambda1*cb2-lambda2*sb2-lambda345*c2b)-lambda6*cb*c3b-lambda7*sb*s3b;
  double lA = c2b*(lambda1*cb2-lambda2*sb2)+lambda345*s2b*s2b-lambda5+2.*lambda6*cb*s3b-2.*lambda7*sb*c3b;
  double s2ba = 2.*lhat*v2;
  double c2ba = -(mA_2 - lA*v2);
  // Check case with degenerate h, H
  if ((abs(s2ba)>sqrt(DBL_EPSILON))||(abs(c2ba)>sqrt(DBL_EPSILON))) {
    double bma = 0.5*atan2(s2ba,c2ba);
    sinba = sin(bma);
    params_set=true;
  } else {
    cerr << "WARNING: Requested masses m_h and m_H too close, precision problems expected\n";
    params_set=false;
  }

  return params_set;
}


bool THDM::set_param_HHG(double Lambda1, double Lambda2, double Lambda3,
			 double Lambda4, double Lambda5, double Lambda6,
			 double tan_beta) {

  if (tan_beta<=0) {
    params_set=false;
    return params_set;
  }

  double l1,l2,l3,l4,l5,l6,l7,m12_2;

  l1 = 2.*(Lambda1+Lambda3);
  l2 = 2.*(Lambda2+Lambda3);
  l3 = 2.*Lambda3+Lambda4;
  l4 = -Lambda4+0.5*(Lambda5+Lambda6);
  l5 = 0.5*(Lambda5-Lambda6);
  l6 = 0.;
  l7 = 0.;

  double beta = atan(tan_beta);
  m12_2 = sm.get_v2()*sin(beta)*cos(beta)*Lambda5;

  return set_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tan_beta);
}


bool THDM::set_param_phys(double m_h,double m_H, double m_A, double m_Hp,
			  double sba, double lambda6, double lambda7,
			  double m12_2,double tan_beta) {

  if (m_h>m_H) {
    cerr << "WARNING: Cannot set physical masses such that m_H < m_h\n";
    params_set = false;
    return params_set;
  }

  // Problematic parameter choices
  if ((tan_beta<=0)||(abs(sba)>1)||(m_h<0)||(m_H<0)||(m_A<0)||(m_Hp<0)) {
    params_set=false;
    return params_set;
  }

  lambda[6]=lambda6;
  lambda[7]=lambda7;
  beta=atan(tan_beta);

  double sb	 = sin(beta);
  double sb2 = sb*sb;
  double cb	 = cos(beta);
  double cb2 = cb*cb;
  double tb	 = tan(beta);
  double ctb = 1./tb;

  double alpha = -asin(sba)+beta;
  double sa  = sin(alpha);
  double sa2 = sa*sa;
  double ca  = cos(alpha);
  double ca2 = ca*ca;

  double cba = sqrt(1.-sba*sba);

  lambda[1]=(m_H*m_H*ca2+m_h*m_h*sa2-m12_2*tb)/v2/cb2-1.5*lambda6*tb+0.5*lambda7*tb*tb*tb;
  lambda[2]=(m_H*m_H*sa2+m_h*m_h*ca2-m12_2*ctb)/v2/sb2+0.5*lambda6*ctb*ctb*ctb-1.5*lambda7*ctb;
  lambda[3]=((m_H*m_H-m_h*m_h)*ca*sa+2.*m_Hp*m_Hp*sb*cb-m12_2)/v2/sb/cb-0.5*lambda6*ctb-0.5*lambda7*tb;
  lambda[4]=((m_A*m_A-2.*m_Hp*m_Hp)*cb*sb+m12_2)/v2/sb/cb-0.5*lambda6*ctb-0.5*lambda7*tb;
  lambda[5]=(m12_2-m_A*m_A*sb*cb)/v2/sb/cb-0.5*lambda6*ctb-0.5*lambda7*tb;
  m22_2 = -0.5/sb*(pow(m_h,2)*ca*sba+pow(m_H,2)*sa*cba)+m12_2*ctb;
  sinba = sba;

  params_set=true;

  return params_set;
}


bool THDM::set_param_higgs(double Lambda1, double Lambda2, double Lambda3,
			   double Lambda4, double Lambda5, double Lambda6,
			   double Lambda7, double m_Hp) {

  if (m_Hp<0) {
    params_set = false;
    return params_set;
  }

  lambda[1]=Lambda1;
  lambda[2]=Lambda2;
  lambda[3]=Lambda3;
  lambda[4]=Lambda4;
  lambda[5]=Lambda5;
  lambda[6]=Lambda6;
  lambda[7]=Lambda7;
  beta=0;
  m22_2=m_Hp*m_Hp-0.5*v2*Lambda3;

  double m_A2 = m_Hp*m_Hp-0.5*v2*(Lambda5-Lambda4);

  double s2ba = -2.*Lambda6*v2;
  double c2ba = -(m_A2+(Lambda5-Lambda1)*v2);
  // Handle special case with degenerate masses
  if ((abs(s2ba)>sqrt(DBL_EPSILON))||(abs(c2ba)>sqrt(DBL_EPSILON))) {
    double bma  =  0.5*atan2(s2ba,c2ba);
    sinba = sin(bma);
    params_set = true;
  } else {
    cerr << "WARNING: Requested masses m_h and m_H too close (" << m_A2+(Lambda5-Lambda1)*v2 << ", " << Lambda6 << "), precision problems expected\n";
    params_set = false;
  }

  return params_set;
}


// Hybrid basis from 1507.04281
bool THDM::set_param_hybrid(double mh, double mH, double cba, double Z4,
			   double Z5, double Z7, double tanb) {

  bool params_set = true;

  double vsmall = 1E-10;

  if ((mh<0)||(abs(cba) > 1)||(abs(mh-mH)<vsmall)) {
    params_set = false;
    return params_set;
  }

  // 2HDMC uses -1 < sin(beta-alpha) < 1 internally, so we convert to this convention
  // by transferring the sign from cba -> sba if necessary (outside first quadrant).
  double sba = sign(cba)*sqrt(1.-cba*cba);

  params_set = set_param_hybrid_sba(mh,mH,sba,Z4,Z5,Z7,tanb);

  return params_set;
}



// Hybrid basis from 1507.04281 (but with sba as the input)
bool THDM::set_param_hybrid_sba(double mh, double mH, double sba, double Z4,
			   double Z5, double Z7, double tanb) {

  double Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,mHp;

  bool params_set = true;
  double vsmall = 1E-10;

  if ((mh<0)||(abs(sba) > 1)||(abs(mh-mH)<vsmall)||(abs(tanb-1.)<vsmall)) {
    params_set = false;
    return params_set;
  }

  double cba = sqrt(1.-sba*sba);
  double v2 = sm.get_v2();

  double Z6 = (pow(mh,2)-pow(mH,2))*cba*sba/v2;

  Lambda1 = 0.;
  Lambda2 = 0.;
  Lambda3 = 0.;
  Lambda4 = Z4;
  Lambda5 = Z5;
  Lambda6 = Z6;
  Lambda7 = Z7;

  // Lambda1 fixed from Higgs masses
  Lambda1 = (mh*mh*sba*sba+mH*mH*cba*cba)/v2;

  double Lambda345;

  double t2b = 2.*tanb/(1.-tanb*tanb);

  if (abs(tanb-1.0)>vsmall) {

	// Lambda2 fixed from requested tan(beta) value
   	Lambda2 = Lambda1 + 2*(Z6+Z7)/t2b;

  	// Only one value of Lambda345 compatible with softly broken Z2 symmetry
	Lambda345 = Lambda1+2.*Z6/t2b-(Z6-Z7)*t2b;
	Lambda3 = Lambda345-Lambda4-Lambda5;
  } else {
    cout << "The point tanb=1 is problematic." << endl;
  }

  double mHp2 = mH*mH*sba*sba + mh*mh*cba*cba-0.5*(Z4+Z5)*v2;

  if (mHp2>0) {
	mHp = sqrt(mHp2);
  } else {
   	params_set = false;
   	return params_set;
  }

  // First set the potential of 2HDM using input in the Higgs basis
  params_set = set_param_higgs(Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, Lambda6, Lambda7, mHp);

  // Recalculate internal state from Higgs basis to desired tanb value
  // (where lambda6=lambda7=0 is manifest).
  if (params_set) recalc_tan_beta(tanb);

  return params_set;
}


bool THDM::set_inert(double m_h,double m_H, double m_A, double m_Hp, double lambda2, double lambda3) {

  if ((m_h<=0)||(m_H<=0)||(m_A<=0)||(m_Hp<=0)) {
    params_set = false;
    return params_set;
  }

  double v2 = sm.get_v2()/2;

  lambda[1]=0.5*pow(m_h,2)/v2;
  lambda[2]=lambda2;
  lambda[3]=lambda3;
  lambda[4]=(0.5*(m_H*m_H+m_A*m_A)-m_Hp*m_Hp)/v2;
  lambda[5]=0.5*(m_H*m_H-m_A*m_A)/v2;
  lambda[6]=0.0;
  lambda[7]=0.0;
  beta=0;
  m22_2=m_Hp*m_Hp-lambda3*v2;

  if (m_h>m_H) {
    sinba = 0.;
  } else {
    sinba = 1.;
  }

  set_yukawas_inert();

  params_set = true;
  return params_set;
}


bool THDM::set_param_sm(double mh) {

   double mA = 1000.;
   if (mA <= mh) mA = mh+1000;

   double m12_2_sm = (pow(mA,2)+pow(mh,2)/2.)/2.;

   bool pset = set_param_phys(mh,mA,mA,mA,1., 0., 0., m12_2_sm, 1.);

   set_yukawas_type(1);

   return pset;

}


bool THDM::set_MSSM(double mA, double tan_beta) {

  if ((mA<=0)||(tan_beta<=0)) {
    params_set = false;
    return params_set;
  }

  double l1,l2,l3,l4,l5,l6,l7,m12_2;

  double beta = atan(tan_beta);
  double cb = cos(beta);
  double sb = sin(beta);

  double g = sm.get_g();
  double gp = sm.get_gprime();

  l1 = (g*g+gp*gp)/4.;
  l2 = l1;
  l3 = (g*g-gp*gp)/4.;
  l4 = -0.5*g*g;
  l5 = 0.;
  l6 = 0.;
  l7 = 0.;
  m12_2=mA*mA*cb*sb;

  double pset = set_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tan_beta);

  if (pset) set_yukawas_type(2);

  return pset;
}


bool THDM::set_hMSSM(double mh, double mA, double tan_beta) {

  if ((mA<=0)||(tan_beta<=0)) {
    params_set = false;
    return params_set;
  }

  double l1,l2,l3,l4,l5,l6,l7,m12_2;

  double beta = atan(tan_beta);
  double cb = cos(beta);
  double sb = sin(beta);
  double cb2 = cb*cb;
  double sb2 = sb*sb;
  double m_A2 = mA*mA;

  double g = sm.get_g();
  double gp = sm.get_gprime();

  l1 = (g*g+gp*gp)/4.;
  l2 = l1;
  l3 = (g*g-gp*gp)/4.;
  l4 = -0.5*g*g;
  l5 = 0.;
  l6 = 0.;
  l7 = 0.;
  m12_2=mA*mA*cb*sb;

  double m_h2 = 0.;
  double dm = 1.E-5;


  while ((abs(sqrt(m_h2)-mh)>dm)&&(l2 < 4.*M_PI)) {

  	double M112   =  m_A2*sb2+v2*(l1*cb2);
   	double M122   = -m_A2*sb*cb+v2*(l3+l4)*sb*cb;
   	double M222   =  m_A2*cb2+v2*l2*sb2;
   	m_h2   =  0.5*(M112+M222-sqrt((M112-M222)*(M112-M222)+4.*M122*M122));

//   printf("%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n", l2, M112, M122, M222, m_h2, sqrt(m_h2));

   	double delta_l2 = (pow(mh,2)-m_h2)/(v2*sb2);
   	l2 = l2+delta_l2;

  }

  double pset = set_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tan_beta);
  if (pset) set_yukawas_type(2);

  return pset;
}


void THDM::set_kappa() {
  set_kappa_D();
  set_kappa_U();
  set_kappa_L();
}


void THDM::set_kappa_D() {
  gsl_matrix *md = sm.get_MD();
  gsl_matrix_memcpy(kappa_D,md);
  gsl_matrix_scale(kappa_D,sqrt(2)/sm.get_v());
  gsl_matrix_free(md);
}


void THDM::set_kappa_U() {
  gsl_matrix *mu = sm.get_MU();
  gsl_matrix_memcpy(kappa_U,mu);
  gsl_matrix_scale(kappa_U,sqrt(2)/sm.get_v());
  gsl_matrix_free(mu);
}


void THDM::set_kappa_L() {
  gsl_matrix *ml = sm.get_ML();
  gsl_matrix_memcpy(kappa_L,ml);
  gsl_matrix_scale(kappa_L,sqrt(2)/sm.get_v());
  gsl_matrix_free(ml);

}


void THDM::set_yukawas_type(int type) {
// Yukawa couplings in convention of hep-ph/0504050

  if (type==0) {
    set_yukawas_inert();
    return;
  }

  if ((tan(beta)==0)||(type<0)||(type>4)) return;

  yukawas_type=type;

  // Warn but do not stop the user from breaking Z2
  if ((abs(lambda[6])>1e-9)||(abs(lambda[7])>1e-9)) {
    cerr << "WARNING: Requested Yukawa type respects Z2-symmetry but lambda6 or lambda7 is not zero\n";
  }

  set_kappa();

  gsl_matrix_memcpy(rho_D,kappa_D);
  gsl_matrix_memcpy(rho_U,kappa_U);
  gsl_matrix_memcpy(rho_L,kappa_L);

  if (type==1) {
    gsl_matrix_scale(rho_D,1./tan(beta));
    gsl_matrix_scale(rho_U,1./tan(beta));
    gsl_matrix_scale(rho_L,1./tan(beta));
  } else if (type==2) {
    gsl_matrix_scale(rho_D,-tan(beta));
    gsl_matrix_scale(rho_U,1./tan(beta));
    gsl_matrix_scale(rho_L,-tan(beta));
  } else if (type==3) {
    gsl_matrix_scale(rho_D,-tan(beta));
    gsl_matrix_scale(rho_U,1./tan(beta));
    gsl_matrix_scale(rho_L,1./tan(beta));
  } else if (type==4) {
    gsl_matrix_scale(rho_D,1./tan(beta));
    gsl_matrix_scale(rho_U,1./tan(beta));
    gsl_matrix_scale(rho_L,-tan(beta));
  }

}

int THDM::get_yukawas_type() {
	return yukawas_type;
}

void THDM::set_yukawas_inert() {

  yukawas_type = 0;

  if ((abs(lambda[6])>1e-9)||(abs(lambda[7])>1e-9)) {
    cout << "WARNING: Requested Yukawa type respects Z2-symmetry but lambda6 or lambda7 is not zero\n";
  }

  set_kappa();

  // Set all rhos to 0
  gsl_matrix_set_zero(rho_D);
  gsl_matrix_set_zero(rho_U);
  gsl_matrix_set_zero(rho_L);
}


void THDM::set_yukawas_down(double rhod, double rhos, double rhob) {

  yukawas_type=-1;

  set_kappa_D();

  gsl_matrix_set_zero(rho_D);

  gsl_matrix_set(rho_D,0,0,rhod);
  gsl_matrix_set(rho_D,1,1,rhos);
  gsl_matrix_set(rho_D,2,2,rhob);
}


void THDM::set_yukawas_up(double rhou, double rhoc, double rhot) {

  yukawas_type=-1;

  set_kappa_U();

  gsl_matrix_set_zero(rho_U);

  gsl_matrix_set(rho_U,0,0,rhou);
  gsl_matrix_set(rho_U,1,1,rhoc);
  gsl_matrix_set(rho_U,2,2,rhot);
}


void THDM::set_yukawas_lepton(double rhoe, double rhomu, double rhotau) {

  yukawas_type=-1;

  set_kappa_L();

  gsl_matrix_set_zero(rho_L);

  gsl_matrix_set(rho_L,0,0,rhoe);
  gsl_matrix_set(rho_L,1,1,rhomu);
  gsl_matrix_set(rho_L,2,2,rhotau);
}


void THDM::set_yukawas_down(double rho11, double rho22, double rho33, double rho12, double rho13, double rho23) {
  yukawas_type=-1;
  set_kappa_D();

  gsl_matrix_set_zero(rho_D);

  gsl_matrix_set(rho_D,0,0,rho11);
  gsl_matrix_set(rho_D,1,1,rho22);
  gsl_matrix_set(rho_D,2,2,rho33);

	// Off-diagonal elements are symmetric
  gsl_matrix_set(rho_D,0,1,rho12);
  gsl_matrix_set(rho_D,0,2,rho13);
  gsl_matrix_set(rho_D,1,2,rho23);
  gsl_matrix_set(rho_D,1,0,rho12);
  gsl_matrix_set(rho_D,2,0,rho13);
  gsl_matrix_set(rho_D,2,1,rho23);
}


void THDM::set_yukawas_up(double rho11, double rho22, double rho33, double rho12, double rho13,double rho23) {
  yukawas_type=-1;
  set_kappa_U();

  gsl_matrix_set_zero(rho_U);

  gsl_matrix_set(rho_U,0,0,rho11);
  gsl_matrix_set(rho_U,1,1,rho22);
  gsl_matrix_set(rho_U,2,2,rho33);

	// Off-diagonal elements are symmetric
  gsl_matrix_set(rho_U,0,1,rho12);
  gsl_matrix_set(rho_U,0,2,rho13);
  gsl_matrix_set(rho_U,1,2,rho23);
  gsl_matrix_set(rho_U,1,0,rho12);
  gsl_matrix_set(rho_U,2,0,rho13);
  gsl_matrix_set(rho_U,2,1,rho23);
}


void THDM::set_yukawas_lepton(double rho11, double rho22, double rho33, double rho12, double rho13, double rho23) {
  yukawas_type=-1;
  set_kappa_L();

  gsl_matrix_set_zero(rho_L);

  gsl_matrix_set(rho_L,0,0,rho11);
  gsl_matrix_set(rho_L,1,1,rho22);
  gsl_matrix_set(rho_L,2,2,rho33);

	// Off-diagonal elements are symmetric
  gsl_matrix_set(rho_L,0,1,rho12);
  gsl_matrix_set(rho_L,0,2,rho13);
  gsl_matrix_set(rho_L,1,2,rho23);
  gsl_matrix_set(rho_L,1,0,rho12);
  gsl_matrix_set(rho_L,2,0,rho13);
  gsl_matrix_set(rho_L,2,1,rho23);
}


void THDM::get_yukawas_down(gsl_matrix *rho_D_out) {
  gsl_matrix_memcpy(rho_D_out,rho_D);
}


void THDM::get_yukawas_up(gsl_matrix *rho_U_out) {
  gsl_matrix_memcpy(rho_U_out,rho_U);
}


void THDM::get_yukawas_lepton(gsl_matrix *rho_L_out) {
  gsl_matrix_memcpy(rho_L_out,rho_L);
}

void THDM::get_kappa_down(double &kd, double &ks, double &kb) {
  kd=gsl_matrix_get(kappa_D,0,0);
  ks=gsl_matrix_get(kappa_D,1,1);
  kb=gsl_matrix_get(kappa_D,2,2);
}

void THDM::get_kappa_up(double &ku, double &kc, double &kt) {
  ku=gsl_matrix_get(kappa_U,0,0);
  kc=gsl_matrix_get(kappa_U,1,1);
  kt=gsl_matrix_get(kappa_U,2,2);
}

void THDM::get_kappa_lepton(double &ke, double &kmu, double &ktau) {
  ke=gsl_matrix_get(kappa_L,0,0);
  kmu=gsl_matrix_get(kappa_L,1,1);
  ktau=gsl_matrix_get(kappa_L,2,2);
}

void THDM::get_kappa_down(double mu,double &kd, double &ks, double &kb) {
  double kd_pole=gsl_matrix_get(kappa_D,0,0);
  double ks_pole=gsl_matrix_get(kappa_D,1,1);
  double kb_pole=gsl_matrix_get(kappa_D,2,2);

   double mb_mb = sm.get_qmass_MSbar(5);
  double mb_pole = sm.get_qmass_pole(5);

  double mb_run = sm.run_qmass_MSbar(mb_mb,mb_mb,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

  kd = kd_pole;
  ks = ks_pole;
  kb = kb_pole/mb_pole*mb_run;
}

void THDM::get_kappa_up(double mu, double &ku, double &kc, double &kt) {
  double ku_pole=gsl_matrix_get(kappa_U,0,0);
  double kc_pole=gsl_matrix_get(kappa_U,1,1);
  double kt_pole=gsl_matrix_get(kappa_U,2,2);

  double mt_mt = sm.get_qmass_MSbar(6);
  double mc_mc = sm.get_qmass_MSbar(4);
  double mc_pole = sm.get_qmass_pole(4);
  double mt_pole = sm.get_qmass_pole(6);

  double mc_run = sm.run_qmass_MSbar(mc_mc,mc_mc,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  double mt_run = sm.run_qmass_MSbar(mt_mt,mt_mt,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

  ku = ku_pole;
  kc = kc_pole/mc_pole*mc_run;
  kt = kt_pole/mt_pole*mt_run;
}

void THDM::get_kappa_lepton(double mu, double &ke, double &kmu, double &ktau) {
  get_kappa_lepton(ke,kmu,ktau);
}


void THDM::get_rho_down(double mu,double &rd, double &rs, double &rb) {
  double rd_pole=gsl_matrix_get(rho_D,0,0);
  double rs_pole=gsl_matrix_get(rho_D,1,1);
  double rb_pole=gsl_matrix_get(rho_D,2,2);

  double mb_mb = sm.get_qmass_MSbar(5);
  double mb_pole = sm.get_qmass_pole(5);

  double mb_run = sm.run_qmass_MSbar(mb_mb,mb_mb,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

  rd = rd_pole;
  rs = rs_pole;
  rb = rb_pole/mb_pole*mb_run;
}

void THDM::get_rho_up(double mu, double &ru, double &rc, double &rt) {
  double ru_pole=gsl_matrix_get(rho_U,0,0);
  double rc_pole=gsl_matrix_get(rho_U,1,1);
  double rt_pole=gsl_matrix_get(rho_U,2,2);

  double mt_mt = sm.get_qmass_MSbar(6);
  double mc_mc = sm.get_qmass_MSbar(4);
  double mc_pole = sm.get_qmass_pole(4);
  double mt_pole = sm.get_qmass_pole(6);

  double mc_run = sm.run_qmass_MSbar(mc_mc,mc_mc,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  double mt_run = sm.run_qmass_MSbar(mt_mt,mt_mt,mu,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

  ru = ru_pole;
  rc = rc_pole/mc_pole*mc_run;
  rt = rt_pole/mt_pole*mt_run;
}

void THDM::get_rho_lepton(double mu, double &re, double &rmu, double &rtau) {
  get_yukawas_lepton(re,rmu,rtau);
}



void THDM::get_yukawas_down(double &rhod, double &rhos, double &rhob) {
  rhod=gsl_matrix_get(rho_D,0,0);
  rhos=gsl_matrix_get(rho_D,1,1);
  rhob=gsl_matrix_get(rho_D,2,2);
}


void THDM::get_yukawas_up(double &rhou, double &rhoc, double &rhot) {
  rhou=gsl_matrix_get(rho_U,0,0);
  rhoc=gsl_matrix_get(rho_U,1,1);
  rhot=gsl_matrix_get(rho_U,2,2);
}


void THDM::get_yukawas_lepton(double &rhoe, double &rhomu, double &rhotau) {
  rhoe=gsl_matrix_get(rho_L,0,0);
  rhomu=gsl_matrix_get(rho_L,1,1);
  rhotau=gsl_matrix_get(rho_L,2,2);
}


void THDM::get_param_gen(double &lambda1, double &lambda2, double &lambda3,
			 double &lambda4, double &lambda5, double &lambda6,
			 double &lambda7, double &m12_2,double &tan_beta) {
  lambda1=lambda[1];
  lambda2=lambda[2];
  lambda3=lambda[3];
  lambda4=lambda[4];
  lambda5=lambda[5];
  lambda6=lambda[6];
  lambda7=lambda[7];
  tan_beta=tan(beta);
  m12_2=get_m12_2();
}


void THDM::get_param_HHG(double &Lambda1, double &Lambda2, double &Lambda3,
			 double &Lambda4, double &Lambda5, double &Lambda6,
			 double &tan_beta) {

  double sb=sin(beta);
  double cb=cos(beta);
  double sbcb=sb*cb;
  double v2 = sm.get_v2();
  double m12_2=get_m12_2();

  double lambda345 = lambda[3]+lambda[4]+lambda[5];

  Lambda1 = 0.5*(lambda[1]-lambda345+2.*m12_2/(v2*sbcb));
  Lambda2 = 0.5*(lambda[2]-lambda345+2.*m12_2/(v2*sbcb));
  Lambda3 = 0.5*(lambda345-2.*m12_2/(v2*sbcb));
  Lambda4 = 2.*m12_2/(v2*sbcb)-lambda[4]-lambda[5];
  Lambda5 = 2.*m12_2/(v2*sbcb);
  Lambda6 = 2.*m12_2/(v2*sbcb)-2.*lambda[5];

  tan_beta = tan(beta);
}


void THDM::get_param_phys(double &m_h,double &m_H, double &m_A, double &m_Hp,
			  double &sba, double &lambda6, double &lambda7,
			  double &m12_2,double &tan_beta) {

  lambda6=lambda[6];
  lambda7=lambda[7];
  tan_beta=tan(beta);
  double sb  = sin(beta);
  double sb2 = sb*sb;
  double cb	 = cos(beta);
  double cb2 = cb*cb;
  double tb	 = tan(beta);
  double ctb = 1./tb;
  double m_A2;
  m12_2=get_m12_2();
  if (tan_beta>0) {
    m_A2=m12_2/sb/cb-0.5*v2*(2*lambda[5]+lambda[6]*ctb+lambda[7]*tb);
  } else {
    m_A2=m22_2+0.5*v2*(lambda[3]+lambda[4]-lambda[5]);
  }
  double m_Hp2  =  m_A2+0.5*v2*(lambda[5]-lambda[4]);
  double M112   =  m_A2*sb2+v2*(lambda[1]*cb2+2.*lambda[6]*sb*cb+lambda[5]*sb2);
  double M122   = -m_A2*sb*cb+v2*((lambda[3]+lambda[4])*sb*cb+lambda[6]*cb2+lambda[7]*sb2);
  double M222   =  m_A2*cb2+v2*(lambda[2]*sb2+2.*lambda[7]*sb*cb+lambda[5]*cb2);
  double m_h2   =  0.5*(M112+M222-sqrt((M112-M222)*(M112-M222)+4.*M122*M122));
  double m_H2   =  0.5*(M112+M222+sqrt((M112-M222)*(M112-M222)+4.*M122*M122));

  sba = sinba;

  // Sanity checks. Masses set negative in case of troubles
  if (m_h2>0)   m_h=sqrt(m_h2);   else m_h=-sqrt(-m_h2);
  if (m_H2>0)   m_H=sqrt(m_H2);   else m_H=-sqrt(-m_H2);
  if (m_A2>0)   m_A=sqrt(m_A2);   else m_A=-sqrt(-m_A2);
  if (m_Hp2>0)  m_Hp=sqrt(m_Hp2); else m_Hp=-sqrt(-m_Hp2);
}


void THDM::get_param_higgs(double &Lambda1, double &Lambda2, double &Lambda3,
			   double &Lambda4, double &Lambda5, double &Lambda6,
			   double &Lambda7, double &m_Hp) {

  double sb=sin(beta);
  double s2b=sin(2.*beta);
  double s3b=sin(3.*beta);
  double s2b2=s2b*s2b;
  double sb2=sb*sb;
  double sb4=sb2*sb2;
  double cb=cos(beta);
  double c2b=cos(2.*beta);
  double c3b=cos(3.*beta);
  double cb2=cb*cb;
  double cb4=cb2*cb2;
  double tb=tan(beta);

  double lambda345=lambda[3]+lambda[4]+lambda[5];

  // See hep-ph/0504050
  Lambda1 =  lambda[1]*cb4+lambda[2]*sb4+0.5*lambda345*s2b2+2.*s2b*(cb2*lambda[6]+sb2*lambda[7]);
  Lambda2 =  lambda[1]*sb4+lambda[2]*cb4+0.5*lambda345*s2b2-2.*s2b*(sb2*lambda[6]+cb2*lambda[7]);
  Lambda3 =  0.25*s2b2*(lambda[1]+lambda[2]-2*lambda345)+lambda[3]-s2b*c2b*(lambda[6]-lambda[7]);
  Lambda4 =  0.25*s2b2*(lambda[1]+lambda[2]-2*lambda345)+lambda[4]-s2b*c2b*(lambda[6]-lambda[7]);
  Lambda5 =  0.25*s2b2*(lambda[1]+lambda[2]-2*lambda345)+lambda[5]-s2b*c2b*(lambda[6]-lambda[7]);
  Lambda6 = -0.5*s2b*(lambda[1]*cb2-lambda[2]*sb2-lambda345*c2b)+cb*c3b*lambda[6]+sb*s3b*lambda[7];
  Lambda7 = -0.5*s2b*(lambda[1]*sb2-lambda[2]*cb2+lambda345*c2b)+sb*s3b*lambda[6]+cb*c3b*lambda[7];

  double m12_2 = get_m12_2();
  double m11_2 = m12_2*tb-0.5*v2*(lambda[1]*cb2+(lambda[3]+lambda[4]+lambda[5])*sb2+3.*lambda[6]*sb*cb+lambda[7]*sb2*tb);
  double M22_2 = m11_2*sb2+m22_2*cb2+m12_2*s2b;

  m_Hp=sqrt(M22_2+0.5*v2*Lambda3);
}


void THDM::get_param_hybrid(double &m_h, double &m_H, double &cba,
			   double &Z4, double &Z5, double &Z7, double &tan_beta) {

   double Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,mHp;

   get_param_higgs(Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,mHp);

   double mh,mH,mA,sinba,lam6,lam7,m12_2,tb;

   get_param_phys(mh,mH,mA,mHp,sinba,lam6,lam7,m12_2,tb);
   if ((abs(lam6)>EPS)||(abs(lam7)>EPS)) {
    printf("\nWARNING: Model has hard Z_2-violation\n");
    printf("Output in H2 basis is inconsistent\n");
   }

   m_h = mh;
   m_H = mH;
   cba = sign(sinba)*sqrt(1.-sinba*sinba);
   tan_beta = tb;
   Z4 = Lambda4;
   Z5 = Lambda5;
   Z7 = Lambda7;
}


void THDM::recalc_tan_beta(double tan_beta) {

	// Only positive tan(beta) allowed
  if (tan_beta < 0) return;

  double l1,l2,l3,l4,l5,l6,l7,m_Hp;
  this->get_param_higgs(l1,l2,l3,l4,l5,l6,l7,m_Hp);
  lambda[1]=l1;
  lambda[2]=l2;
  lambda[3]=l3;
  lambda[4]=l4;
  lambda[5]=l5;
  lambda[6]=l6;
  lambda[7]=l7;
  beta=-atan(tan_beta);
  m22_2=m_Hp*m_Hp-0.5*v2*l3;
  double la1,la2,la3,la4,la5,la6,la7,ma_Hp;
  this->get_param_higgs(la1,la2,la3,la4,la5,la6,la7,ma_Hp);
  lambda[1]=la1;
  lambda[2]=la2;
  lambda[3]=la3;
  lambda[4]=la4;
  lambda[5]=la5;
  lambda[6]=la6;
  lambda[7]=la7;
  double sb=sin(beta);
  double s2b=sin(2.*beta);
  double sb2=sb*sb;
  double cb=cos(beta);
  double cb2=cb*cb;
  m22_2=-0.5*v2*l1*sb2+(m_Hp*m_Hp-0.5*v2*l3)*cb2+0.5*v2*l6*s2b;
  beta=atan(tan_beta);
}


complex <double> THDM::get_qki(int k, int i) {

  if ((k<1)||(k>4)||(i<1)||(i>2)) {
    complex <double> zero(0.0,0.0);
    return zero;
  }

  complex <double> I(0.0,1.0);

  double sba = get_sba();
  double cba = get_cba();

  complex <double> qki[4][2] = {{ sba, -cba,},
				{ cba, sba, },
				{ 0  , I },
				{ I  , 0 }};

  return qki[k-1][i-1];

}


void THDM::get_coupling_hdd(int h,int f1,int f2,complex <double> &cs, complex <double> &cp) {
  complex <double> I(0.0,1.0);

  complex <double> x(0.,0.);
  complex <double> y(0.,0.);
  cs = 0.;
  cp = 0.;

  if ((h<1)||(h>3)) return;

  double sba = get_sba();
  double cba = get_cba();

  if ((f1<=3)&&(f1>=1)&&(f2<=3)&&(f2>=1)) {
    double kd = gsl_matrix_get(kappa_D,f1-1,f2-1);
    double rd = gsl_matrix_get(rho_D,f1-1,f2-1);

    if (f1==f2) {
      double mms  = sm.get_dmass_MSbar(f1);
      double mp   = sm.get_dmass_pole(f1);

      // Starting scale for MSbar mass evolution
      double Qinit = mms;
      if (f1==2) {
       // Special case of strange mass where scale is not ms(ms) but ms(Q_ms=2 GeV)
       Qinit = SM::Q_ms;
      }

      // Starting scale for MSbar mass evolution for comp with HD
      if(sm.b_HD) Qinit = SM::Q_HD;

      double mrun = sm.run_qmass_MSbar(mms,Qinit,get_hmass(h),sm.get_qmass_pole(6),sm.get_qmass_pole(5));

      if (mrun>0) {
        kd = kd/mp*mrun;
        rd = rd/mp*mrun;
      }
    }

    switch(h) {
    case 1:
      x = -I/sqrt(2.)*(kd*sba+rd*cba);
      y = 0.0;
      break;
    case 2:
      x = -I/sqrt(2.)*(kd*cba-rd*sba);
      y = 0.0;
      break;
    case 3:
      x = 0.0;
      y = -I*I/sqrt(2.)*rd;
      break;
    }

    cs = x;
    cp = y;
  }
}

void THDM::get_coupling_huu(int h,int f1,int f2,complex <double> &cs, complex <double> &cp) {
  complex <double> I(0.0,1.0);

  complex <double> x(0.,0.);
  complex <double> y(0.,0.);
  cs = 0.;
  cp = 0.;


  if ((h<1)||(h>3)) return;

  double sba = get_sba();
  double cba = get_cba();

  if ((f1<=3)&&(f1>=1)&&(f2<=3)&&(f2>=1)) {
    double ku = gsl_matrix_get(kappa_U,f1-1,f2-1);
    double ru = gsl_matrix_get(rho_U,f1-1,f2-1);

    double mms  = sm.get_umass_MSbar(f1);
    double mp   = sm.get_umass_pole(f1);
    double Qinit = mms;
    // Starting scale for MSbar mass evolution for comp with HD
    if(sm.b_HD) Qinit = SM::Q_HD;

    double mrun = sm.run_qmass_MSbar(mms,Qinit,get_hmass(h),sm.get_qmass_pole(6),sm.get_qmass_pole(5));

    if (f1==f2) {
      if (mrun>0) {
        ku = ku/mp*mrun;
        ru = ru/mp*mrun;
      }
    }

    switch(h) {
    case 1:
      x = -I/sqrt(2.)*(ku*sba+ru*cba);
      y = 0.0;
      break;
    case 2:
      x = -I/sqrt(2.)*(ku*cba-ru*sba);
      y = 0.0;
      break;
    case 3:
      x = 0.0;
      y = -I*(-I/sqrt(2.))*ru;
      break;
    }

    cs = x;
    cp = y;

  }
}

void THDM::get_coupling_hll(int h,int f1,int f2,complex <double> &cs, complex <double> &cp) {
  complex <double> I(0.0,1.0);

  complex <double> x(0.,0.);
  complex <double> y(0.,0.);
  cs = 0.;
  cp = 0.;

  if ((h<1)||(h>3)) return;

  double sba = get_sba();
  double cba = get_cba();

  if ((f1<=3)&&(f1>=1)&&(f2<=3)&&(f2>=1)) {
    double kl = gsl_matrix_get(kappa_L,f1-1,f2-1);
    double rl = gsl_matrix_get(rho_L,f1-1,f2-1);

    switch(h) {
    case 1:
      x = -I/sqrt(2.)*(kl*sba+rl*cba);
      y = 0.0;
      break;
    case 2:
      x = -I/sqrt(2.)*(kl*cba-rl*sba);
      y = 0.0;
      break;
    case 3:
      x = 0.0;
      y = -I*I/sqrt(2.)*rl;
      break;
    }

    cs = x;
    cp = y;
  }
}

void THDM::get_coupling_hdu(int h,int d,int u,complex <double> &cs, complex <double> &cp) {
  complex <double> I(0.0,1.0);

  complex <double> x(0.,0.);
  complex <double> y(0.,0.);
  cs = 0.;
  cp = 0.;

  gsl_matrix* ckm = sm.get_CKM_matrix();

  if (h!=4) return;

  double mHp = get_hmass(4);

  if ((u<=3)&&(u>=1)&&(d<=3)&&(d>=1)) {
    gsl_matrix *RD = gsl_matrix_alloc(3,3);
    gsl_matrix *RU = gsl_matrix_alloc(3,3);

    gsl_matrix_view A = gsl_matrix_submatrix(ckm,0,0,3,3);
    gsl_matrix_view B = gsl_matrix_submatrix(rho_D,0,0,3,3);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0, &A.matrix, &B.matrix,0.0,RD);
    B = gsl_matrix_submatrix(rho_U,0,0,3,3);
    gsl_blas_dgemm(CblasConjTrans,CblasNoTrans,1.0, &B.matrix, &A.matrix,0.0,RU);

    double rd = gsl_matrix_get(RD,u-1,d-1);
    double ru = gsl_matrix_get(RU,u-1,d-1);
    gsl_matrix_free(RD);
    gsl_matrix_free(RU);
    gsl_matrix_free(ckm);


    double mdms  = sm.get_dmass_MSbar(d);
    double Qinit  = mdms;
    double mdp   = sm.get_dmass_pole(d);
    double mums  = sm.get_umass_MSbar(u);
    double mup   = sm.get_umass_pole(u);

    // Special case for strange quark not defined as ms(ms)
    if (d==2) {
     Qinit = SM::Q_ms;
    }

     // Starting scale for MSbar mass evolution for comp with HD
    if(sm.b_HD) Qinit = SM::Q_HD;
    double mdrun = sm.run_qmass_MSbar(mdms,Qinit,mHp,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

    Qinit  = mums;
    // Starting scale for MSbar mass evolution for comp with HD
    if(sm.b_HD) Qinit = SM::Q_HD;
    double murun = sm.run_qmass_MSbar(mums,Qinit,mHp,sm.get_qmass_pole(6),sm.get_qmass_pole(5));


    if (mdrun>0) {
      rd = rd/mdp*mdrun;
    }
    if (murun>0) {
      ru = ru/mup*murun;
    }

    cs = -I*0.5*(rd-ru);
    cp = -I*0.5*(rd+ru);
  }


}

void THDM::get_coupling_hln(int h,int l,int n,complex <double> &cs, complex <double> &cp) {
  complex <double> I(0.0,1.0);

  complex <double> x(0.,0.);
  complex <double> y(0.,0.);
  cs = 0.;
  cp = 0.;

  if (h!=4) return;

  if ((n<=3)&&(n>=1)&&(l<=3)&&(l>=1)) {
    double rl = gsl_matrix_get(rho_L,n-1,l-1);
    double rn = gsl_matrix_get(rho_N,n-1,l-1);

    cs = -I*0.5*(rl-rn);
    cp = -I*0.5*(rl+rn);
  }

}

void THDM::get_coupling_vhh(int v,int h1,int h2,complex <double> &c) {
  // Convention of hep-ph/0602242
  complex <double> I(0.0,1.0);

  double e=sm.get_e();
  double g=sm.get_g();
  double costw=sm.get_costw();
  double cos2tw=cos(2.*acos(costw));

  c=0;
  if (h1>h2) {
    int tmp=h1;
    h1=h2;
    h2=tmp;
  }

  // gamma
  if (v==1) {
    if ((h1==4)&&(h2==4))
      c=-I*e;
  }
  // Z
  if (v==2) {
    if ((h1<4)&&(h2<4))
      c=g/2./costw*imag(get_qki(h1,1)*conj(get_qki(h2,1))+get_qki(h1,2)*conj(get_qki(h2,2)));
    if ((h1==4)&&(h2==4))
      c=-I*g*cos2tw/2./costw;
  }
  // W+  (to get coupling for W- take complex conjugate)
  if ((v==3)&&(h1<4)&&(h2==4)) {
    c=conj(-I*g/2.*get_qki(h1,2));
  }

}


void THDM::get_coupling_vvh(int v1,int v2,int h,complex <double> &c) {
  // Convention of hep-ph/0602242
  complex <double> I(0.0,1.0);

  double g=sm.get_g();
  double costw=sm.get_costw();
  double MW=sm.get_MW();
  double MZ=sm.get_MZ();

  c=0;

  // Z Z
  if ((v1==2)&&(v2==2)&&(h<4)) {
    c=I*g/costw*MZ*real(get_qki(h,1));
  }
  // W+ W-
  if ((v1==3)&&(v2==3)&&(h<4)) {
    c=I*g*MW*real(get_qki(h,1));
  }

}


void THDM::get_coupling_vvhh(int v1,int v2,int h1,int h2,complex <double> &c) {
  // Convention of hep-ph/0602242
  complex <double> I(0.0,1.0);

  double g=sm.get_g();
  double e=sm.get_e();
  double costw=sm.get_costw();
  double cos2tw=cos(2.*acos(costw));

  c=0;
  if (v1>v2) {
    int tmp=v1;
    v1=v2;
    v2=tmp;
  }
  if (h1>h2) {
    int tmp=h1;
    h1=h2;
    h2=tmp;
  }

  // Z Z h/H/A h/H/A
  if ((v1==2)&&(v2==2)&&(h1<4)&&(h2<4)){
    c=I/2.*g*g/costw/costw*real(conj(get_qki(h1,1))*get_qki(h2,1)+conj(get_qki(h1,2))*get_qki(h2,2));
  }
  // W+ W- h/H/A h/H/A
  if ((v1==3)&&(v2==3)&&(h1<4)&&(h2<4)){
    c=I/2.*g*g*real(conj(get_qki(h1,1))*get_qki(h2,1)+conj(get_qki(h1,2))*get_qki(h2,2));
  }
  // gamma gamma H+ H-
  if ((v1==1)&&(v2==1)&&(h1==4)&&(h2==4)){
    c=I*2.*e*e;
  }
  // gamma Z H+ H-
  if ((v1==1)&&(v2==2)&&(h1==4)&&(h2==4)){
    c=I*g*e/costw*cos2tw;
  }
  // Z Z H+ H-
  if ((v1==2)&&(v2==2)&&(h1==4)&&(h2==4)){
    c=I/2.*g*g/costw/costw*cos2tw*cos2tw;
  }
  // W+ W- H+ H-
  if ((v1==3)&&(v2==3)&&(h1==4)&&(h2==4)){
    c=I/2.*g*g;
  }
  // gamma W+ h/H/A H-
  if ((v1==1)&&(v2==3)&&(h1<4)&&(h2==4)){
    c=conj(I/2.*g*e*get_qki(h1,2));
  }
  // Z W+ h/H/A H-
  if ((v1==2)&&(v2==3)&&(h1<4)&&(h2==4)){
    c=conj(-I/2.*g*g/costw*(1.-costw*costw)*get_qki(h1,2));
  }
}


void THDM::get_coupling_hhh(int h1,int h2,int h3,complex <double> &c) {
  // Conventions are according to hep-ph/0602242
  complex <double> I(0.0,1.0);

  double Z1,Z2,Z3,Z4,Z5,l6,l7,m_Hp;
  get_param_higgs(Z1,Z2,Z3,Z4,Z5,l6,l7,m_Hp);
  double Z6=-l6,Z7=-l7;
  double v=sm.get_v();

  c=0;
  if (h1>h2) {
    int tmp=h1;
    h1=h2;
    h2=tmp;
  }
  if (h1>h3) {
    int tmp=h1;
    h1=h3;
    h3=h2;
    h2=tmp;
  }
  if (h2>h3) {
    int tmp=h2;
    h2=h3;
    h3=tmp;
  }

  //Check for odd number of A
  if (((h1==3)&&(h2!=3)&&(h3!=3))||
      ((h1!=3)&&(h2==3)&&(h3!=3))||
      ((h1!=3)&&(h2!=3)&&(h3==3))||
      ((h1==3)&&(h2==3)&&(h3==3))) return;

  gsl_permutation *per=gsl_permutation_alloc(3);
  gsl_permutation_init(per);

  // Only neutral higgses
  int i;
  if ((h1<4)&&(h2<4)&&(h3<4))
    for (i=0;i<6;i++) {
      int h[3]={h1,h2,h3};
      gsl_permute_int(gsl_permutation_data(per),h,1,3);
      gsl_permutation_next(per);
      c+=-I*v/2.*(get_qki(h[0],1)*conj(get_qki(h[1],1))*real(get_qki(h[2],1))*Z1+
		  get_qki(h[0],2)*conj(get_qki(h[1],2))*real(get_qki(h[2],1))*(Z3+Z4)+
		  real(conj(get_qki(h[0],1))*get_qki(h[1],2)*get_qki(h[2],2)*Z5)+
		  real((2.*get_qki(h[0],1)+conj(get_qki(h[0],1)))*conj(get_qki(h[1],1))*get_qki(h[2],2)*Z6)+
		  real(conj(get_qki(h[0],2))*get_qki(h[1],2)*get_qki(h[2],2)*Z7));
    }
  // Neutral and charged higgses
  if ((h1<4)&&(h2==4)&&(h3==4))
    c=-I*v*(real(get_qki(h1,1)*Z3+real(get_qki(h1,2)*Z7)));

  gsl_permutation_free (per);

}

void THDM::get_coupling_hhhh(int h1,int h2,int h3,int h4,complex <double> &c) {
  // Conventions are according to hep-ph/0602242
  complex <double> I(0.0,1.0);

  double Z1,Z2,Z3,Z4,Z5,l6,l7,m_Hp;
  get_param_higgs(Z1,Z2,Z3,Z4,Z5,l6,l7,m_Hp);
  double Z6=-l6,Z7=-l7;

  c=0;
  
//     printf("before sort h1 h2 h3 h4 %d %d %d %d %16.8E %16.8E \n", h1, h2, h3, h4, real(c), imag(c));

  if (h1>h2) {
    int tmp=h1;
    h1=h2;
    h2=tmp;
  }
  if (h1>h3) {
    int tmp=h1;
    h1=h3;
    h3=h2;
    h2=tmp;
  }
  if (h2>h3) {
    int tmp=h2;
    h2=h3;
    h3=tmp;
  }
  if (h1>h4) {
    int tmp=h1;
    h1=h4;
    h4=h3;
    h3=h2;
    h2=tmp;
  }
  if (h2>h4) {
    int tmp=h2;
    h2=h4;
    h4=h3;
    h3=tmp;
  }
  if (h3>h4) {
    int tmp=h3;
    h3=h4;
    h4=tmp;
  }

//      printf("after sort h1 h2 h3 h4 %d %d %d %d %16.8E %16.8E \n", h1, h2, h3, h4, real(c), imag(c));

  //Check for odd number of A
  if (((h1==3)&&(h2!=3)&&(h3!=3)&&(h4!=3))||
      ((h1!=3)&&(h2==3)&&(h3!=3)&&(h4!=3))||
      ((h1!=3)&&(h2!=3)&&(h3==3)&&(h4!=3))||
      ((h1!=3)&&(h2!=3)&&(h3!=3)&&(h4==3))||
      ((h1!=3)&&(h2==3)&&(h3==3)&&(h4==3))||
      ((h1==3)&&(h2!=3)&&(h3==3)&&(h4==3))||
      ((h1==3)&&(h2==3)&&(h3!=3)&&(h4==3))||
      ((h1==3)&&(h2==3)&&(h3==3)&&(h4!=3))) return;

  gsl_permutation *per=gsl_permutation_alloc(4);
  gsl_permutation_init(per);

  // Only neutral higgses
  int i;
  if ((h1<4)&&(h2<4)&&(h3<4)&&(h4<4)) {
    for (i=0;i<24;i++) {
      int h[4]={h1,h2,h3,h4};
      gsl_permute_int(gsl_permutation_data(per),h,1,4);
      gsl_permutation_next(per);
      c+=-I/8.*(get_qki(h[0],1)*get_qki(h[1],1)*conj(get_qki(h[2],1))*conj(get_qki(h[3],1))*Z1+
		get_qki(h[0],2)*get_qki(h[1],2)*conj(get_qki(h[2],2))*conj(get_qki(h[3],2))*Z2+
		2.*get_qki(h[0],1)*conj(get_qki(h[1],1))*get_qki(h[2],2)*conj(get_qki(h[3],2))*(Z3+Z4)+
		2.*real(conj(get_qki(h[0],1))*conj(get_qki(h[1],1))*get_qki(h[2],2)*get_qki(h[3],2)*Z5)+
		4.*real(get_qki(h[0],1)*conj(get_qki(h[1],1))*conj(get_qki(h[2],1))*get_qki(h[3],2)*Z6)+
		4.*real(conj(get_qki(h[0],1))*get_qki(h[1],2)*get_qki(h[2],2)*conj(get_qki(h[3],2))*Z7));
    }
  }

  // Neutral and charged higgses
  if ((h1<4)&&(h2<4)&&(h3==4)&&(h4==4)) {
    c=-I/2.*(get_qki(h1,2)*conj(get_qki(h2,2))*Z2+get_qki(h1,1)*conj(get_qki(h2,1))*Z3+
	     2.*real(get_qki(h1,1)*get_qki(h2,2)*Z7)+
	     get_qki(h2,2)*conj(get_qki(h1,2))*Z2+get_qki(h2,1)*conj(get_qki(h1,1))*Z3+
	     2.*real(get_qki(h2,1)*get_qki(h1,2)*Z7));
  }

  // Only charged higgses
  if ((h1==4)&&(h2==4)&&(h3==4)&&(h4==4))  {
    c=-I*2.*Z2;
  }

  gsl_permutation_free (per);
}


double THDM::calc_unitarity() {

  double egmax=0;

  // S-Matrices from Ginzburg and Ivanov, hep-ph/0508020
  double S20 = lambda[3]-lambda[4];
  if (abs(S20) > abs(egmax)) egmax = S20;

  double s2 = sqrt(2.);

  double S21_data[] = { lambda[1],    lambda[5],    s2*lambda[6],
                        lambda[5],    lambda[2],    s2*lambda[7],
                        s2*lambda[6], s2*lambda[7], lambda[3]+lambda[4] };

  double S01_data[] = { lambda[1],  lambda[4],  lambda[6],  lambda[6],
                        lambda[4],  lambda[2],  lambda[7],  lambda[7],
                        lambda[6],  lambda[7],  lambda[3],  lambda[5],
                        lambda[6],  lambda[7],  lambda[5],  lambda[3] };

  double S00_data[] = { 3.*lambda[1],2.*lambda[3]+lambda[4],3.*lambda[6],3.*lambda[6],
                        2.*lambda[3]+lambda[4],3.*lambda[2],3.*lambda[7],3.*lambda[7],
                        3.*lambda[6],3.*lambda[7],lambda[3]+2.*lambda[4],3.*lambda[5],
                        3.*lambda[6],3.*lambda[7],3.*lambda[5],lambda[3]+2.*lambda[4] };

  gsl_matrix_view S21 = gsl_matrix_view_array(S21_data,3,3);
  gsl_matrix_view S01 = gsl_matrix_view_array(S01_data,4,4);
  gsl_matrix_view S00 = gsl_matrix_view_array(S00_data,4,4);

  gsl_eigen_symm_workspace *w3 = gsl_eigen_symm_alloc(3);
  gsl_eigen_symm_workspace *w4 = gsl_eigen_symm_alloc(4);
  gsl_vector *eval3 = gsl_vector_alloc(3);
  gsl_vector *eval4 = gsl_vector_alloc(4);

  double eg = 0.0;


  gsl_eigen_symm(&S21.matrix,eval3,w3);
  for (int i=0;i<3;i++) {
    eg = gsl_vector_get(eval3,i);
    if (abs(eg)>abs(egmax)) egmax = eg;
  }

  gsl_eigen_symm(&S01.matrix,eval4,w4);
  for (int i=0;i<4;i++) {
    eg = gsl_vector_get(eval4,i);
    if (abs(eg)>abs(egmax)) egmax = eg;
  }

  gsl_eigen_symm(&S00.matrix,eval4,w4);
  for (int i=0;i<4;i++) {
    eg = gsl_vector_get(eval4,i);
    if (abs(eg)>abs(egmax)) egmax = eg;
  }

  gsl_eigen_symm_free(w3);
  gsl_eigen_symm_free(w4);
  gsl_vector_free(eval3);
  gsl_vector_free(eval4);
  return egmax;
}


bool THDM::check_unitarity(double unitarity_limit) {

  double egmax;
  bool check=true;

  egmax = calc_unitarity();

  if (abs(egmax)>unitarity_limit) check=false;
  return check;
}

void THDM::calc_perturbativity(complex <double> &gmax,int &imax,int &jmax,int &kmax,int &lmax) {

  gmax=0.0;
  imax=0,jmax=0,kmax=0,lmax=0;

  complex <double> gval;

  for (int i=1;i<5;i++) {
     for (int j=1;j<5;j++) {
        for (int k=1;k<5;k++) {
           for (int l=1;l<5;l++) {
              get_coupling_hhhh(i,j,k,l,gval);
              if (abs(gval)>abs(gmax)) {
                 gmax = gval;
                 imax=i;
                 jmax=j;
                 kmax=k;
                 lmax=l;
              }
           }
        }
     }
  }
}

bool THDM::check_perturbativity(double perturbativity_limit) {

  complex <double> gmax(0.0,0.0);
  bool check=true;
  int imax=0,jmax=0,kmax=0,lmax=0;

  calc_perturbativity(gmax,imax,jmax,kmax,lmax);
   if (abs(gmax)>perturbativity_limit) check=false;
  return check;
}

bool THDM::check_stability() {
  stability_params p={{lambda[0],lambda[1],lambda[2],lambda[3],lambda[4],lambda[5],lambda[6],lambda[7]},0.,M_PI/2.,M_PI/4.,0};

  // Check for gamma=0
  if (lambda[1]<0) return false;
  // Check for gamma=Pi/2
  if (lambda[2]<0) return false;
  // Check for rho=0
  if (lambda[3]<-sqrt(lambda[1]*lambda[2])) return false;
  // Check for cos(theta)=0
  if (lambda[3]+lambda[4]-lambda[5]<-sqrt(lambda[1]*lambda[2])) return false;
  if ((abs(lambda[6])<EPS)&&(abs(lambda[7])<EPS)) {
    // Check if lambda6 and lambda7 are zero
    if (lambda[3]+lambda[4]-abs(lambda[5])<-sqrt(lambda[1]*lambda[2])) return false;
    // No more conditions exist if lambda6 and lambda7 are zero
    return true;
  }

  double rho,gamma,ct;
  // Check cos(theta)=+-1 with first gamma solution and rho!=1
  gamma=acos(sqrt((4*lambda[6]*lambda[7]-2*lambda[4]*lambda[3]+lambda[4]*lambda[2]+lambda[4]*lambda[1]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]+lambda[1]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7])*(-2*lambda[7]*lambda[7]+lambda[2]*lambda[5]+lambda[4]*lambda[2]-lambda[3]*lambda[5]-lambda[4]*lambda[3]+2*lambda[6]*lambda[7]))/(4*lambda[6]*lambda[7]-2*lambda[4]*lambda[3]+lambda[4]*lambda[2]+lambda[4]*lambda[1]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]+lambda[1]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7]));
  rho=sin(gamma)*(lambda[7]-cos(gamma)*cos(gamma)*lambda[7]+lambda[6]*cos(gamma)*cos(gamma))/(cos(gamma)*(-lambda[5]-lambda[4]+cos(gamma)*cos(gamma)*lambda[4]+cos(gamma)*cos(gamma)*lambda[5]));
  // Check abs(rho) [0,1]
  if ((abs(rho)<=1.) &&
      // Check gamma [0,Pi/2]
      ((gamma>=0.) && (gamma<=M_PI/2.)) &&
      ((1./2.)*(-lambda[4]*lambda[3]*lambda[3]+lambda[4]*lambda[2]*lambda[1]-lambda[3]*lambda[3]*lambda[5]+lambda[2]*lambda[1]*lambda[5]-2*lambda[7]*lambda[7]*lambda[1]+4*lambda[6]*lambda[7]*lambda[3]-2*lambda[2]*lambda[6]*lambda[6])/(lambda[4]*lambda[1]+lambda[4]*lambda[2]-2*lambda[4]*lambda[3]+lambda[1]*lambda[5]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7]+4*lambda[6]*lambda[7])<0))
    return false;
  // Check cos(theta)=+-1 with second gamma solution and rho!=1
  gamma=M_PI-acos(sqrt((4*lambda[6]*lambda[7]-2*lambda[4]*lambda[3]+lambda[4]*lambda[2]+lambda[4]*lambda[1]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]+lambda[1]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7])*(-2*lambda[7]*lambda[7]+lambda[2]*lambda[5]+lambda[4]*lambda[2]-lambda[3]*lambda[5]-lambda[4]*lambda[3]+2*lambda[6]*lambda[7]))/(4*lambda[6]*lambda[7]-2*lambda[4]*lambda[3]+lambda[4]*lambda[2]+lambda[4]*lambda[1]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]+lambda[1]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7]));
  rho=sin(gamma)*(lambda[7]-cos(gamma)*cos(gamma)*lambda[7]+lambda[6]*cos(gamma)*cos(gamma))/(cos(gamma)*(-lambda[5]-lambda[4]+cos(gamma)*cos(gamma)*lambda[4]+cos(gamma)*cos(gamma)*lambda[5]));
  // Check abs(rho) [0,1]
  if ((abs(rho)<=1.) &&
      // Check gamma [0,Pi/2]
      ((gamma>=0.) && (gamma<=M_PI/2.)) &&
      ((1./2.)*(-lambda[4]*lambda[3]*lambda[3]+lambda[4]*lambda[2]*lambda[1]-lambda[3]*lambda[3]*lambda[5]+lambda[2]*lambda[1]*lambda[5]-2*lambda[7]*lambda[7]*lambda[1]+4*lambda[6]*lambda[7]*lambda[3]-2*lambda[2]*lambda[6]*lambda[6])/(lambda[4]*lambda[1]+lambda[4]*lambda[2]-2*lambda[4]*lambda[3]+lambda[1]*lambda[5]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]-2*lambda[6]*lambda[6]-2*lambda[7]*lambda[7]+4*lambda[6]*lambda[7])<0))
    return false;

  // Check for rho=1 and ct!=+-1
  ct=(1./2.)*(-lambda[6]*lambda[3]-lambda[6]*lambda[4]+lambda[6]*lambda[2]+lambda[5]*lambda[6]+lambda[7]*lambda[1]-lambda[7]*lambda[3]-lambda[7]*lambda[4]+lambda[7]*lambda[5])/sqrt((-lambda[3]*lambda[5]-lambda[5]*lambda[4]+lambda[2]*lambda[5]+lambda[5]*lambda[5]+lambda[6]*lambda[7]-lambda[7]*lambda[7])*(lambda[1]*lambda[5]+lambda[6]*lambda[7]-lambda[3]*lambda[5]+lambda[5]*lambda[5]-lambda[5]*lambda[4]-lambda[6]*lambda[6]));
  gamma=atan(sqrt((-lambda[3]*lambda[5]-lambda[5]*lambda[4]+lambda[2]*lambda[5]+lambda[5]*lambda[5]+lambda[6]*lambda[7]-lambda[7]*lambda[7])*(lambda[1]*lambda[5]+lambda[6]*lambda[7]-lambda[3]*lambda[5]+lambda[5]*lambda[5]-lambda[5]*lambda[4]-lambda[6]*lambda[6]))/(-lambda[3]*lambda[5]-lambda[5]*lambda[4]+lambda[2]*lambda[5]+lambda[5]*lambda[5]+lambda[6]*lambda[7]-lambda[7]*lambda[7]));
  // Check ct [-1,1]
  if ((abs(ct)<=1.) &&
      // Check abs(gamma) [0,Pi/2]
      (abs(gamma)<=M_PI/2.) &&
      // Check stability
      ((1./2.)*(lambda[1]*lambda[2]*lambda[5]-lambda[1]*lambda[7]*lambda[7]-lambda[5]*lambda[5]*lambda[5]+2*lambda[5]*lambda[5]*lambda[4]+2*lambda[5]*lambda[5]*lambda[3]-lambda[5]*lambda[4]*lambda[4]-2*lambda[5]*lambda[6]*lambda[7]-lambda[5]*lambda[3]*lambda[3]-2*lambda[5]*lambda[4]*lambda[3]-lambda[6]*lambda[6]*lambda[2]+2*lambda[6]*lambda[7]*lambda[3]+2*lambda[6]*lambda[7]*lambda[4])/(lambda[1]*lambda[5]+2*lambda[6]*lambda[7]-2*lambda[3]*lambda[5]+lambda[2]*lambda[5]-2*lambda[5]*lambda[4]-lambda[6]*lambda[6]-lambda[7]*lambda[7]+2*lambda[5]*lambda[5])<0))
    return false;

  // Check rho=1 and ct=+1
  // Find first minimum
  if (stability_minimum(p)) {
    // Check first minimum
    if (p.v<0) return false;
    // Find second minimum
    stability_params p21=p;
    stability_params p22=p;
    p21.u=p.m;
    p21.m=(p21.l+p21.u)/2.;
    p22.l=p.m;
    p22.m=(p22.u+p22.l)/2.;
    if (stability_minimum(p21)) {
      // Check second minimum
      if (p21.v<0) return false;
      // Find third minimum
      stability_params p31=p21;
      stability_params p32=p21;
      stability_params p33=p;
      p31.u=p21.m;
      p31.m=(p31.l+p31.u)/2.;
      p32.l=p21.m;
      p32.m=(p32.l+p32.u)/2.;
      p33.l=p21.u;
      p33.m=(p33.l+p33.u)/2.;
      // Check third minimum
      if ((stability_minimum(p31))&&(p31.v<0)) return false;
      if ((stability_minimum(p32))&&(p32.v<0)) return false;
      if ((stability_minimum(p33))&&(p33.v<0)) return false;
    } else if (stability_minimum(p22)) {
      // Check second minimum
      if (p22.v<0) return false;
      // Find third minimum
      stability_params p31=p22;
      stability_params p32=p22;
      p31.u=p22.m;
      p31.m=(p31.l+p31.u)/2.;
      p32.l=p22.m;
      p32.m=(p32.l+p32.u)/2.;
      // Check third minimum
      if ((stability_minimum(p31))&&(p31.v<0)) return false;
      if ((stability_minimum(p32))&&(p32.v<0)) return false;
    }
  }
  // Check rho=1 and ct=-1
  p.lambda[6]=-lambda[6];
  p.lambda[7]=-lambda[7];
  p.m=(p.u+p.l)/2.;
  // Find first minimum
  if (stability_minimum(p)) {
    // Check first minimum
    if (p.v<0) return false;
    // Find second minimum
    stability_params p21=p;
    stability_params p22=p;
    p21.u=p.m;
    p21.m=(p21.l+p21.u)/2.;
    p22.l=p.m;
    p22.m=(p22.l+p22.u)/2.;
    if (stability_minimum(p21)) {
      // Check second minimum
      if (p21.v<0) return false;
      // Find third minimum
      stability_params p31=p21;
      stability_params p32=p21;
      stability_params p33=p;
      p31.u=p21.m;
      p31.m=(p31.l+p31.u)/2.;
      p32.l=p21.m;
      p32.m=(p32.l+p32.u)/2.;
      p33.l=p21.u;
      p33.m=(p33.l+p33.u)/2.;
      // Check third minimum
      if ((stability_minimum(p31))&&(p31.v<0)) return false;
      if ((stability_minimum(p32))&&(p32.v<0)) return false;
      if ((stability_minimum(p33))&&(p33.v<0)) return false;
    } else if (stability_minimum(p22)) {
      // Check second minimum
      if (p22.v<0) return false;
      // Find third minimum
      stability_params p31=p22;
      stability_params p32=p22;
      p31.u=p22.m;
      p31.m=(p31.l+p31.u)/2.;
      p32.l=p22.m;
      p32.m=(p32.l+p32.u)/2.;
      // Check third minimum
      if ((stability_minimum(p31))&&(p31.v<0)) return false;
      if ((stability_minimum(p32))&&(p32.v<0)) return false;
    }
  }

  return true;

}

void THDM::write_LesHouches(const char* file, bool fulldecay, bool couplings, bool qcd_on, const HBHSResult *hbhs) {
  complex <double> I(0.0,1.0);

  FILE* output;
  output=fopen(file,"w");

  if (output==NULL) {
    printf("Cannot open file \"%s\" for writing\n", file);
    return;
  }

  double mh,mH,mA,mHp,sba,lambda6,lambda7,tan_beta,m12_2;
  get_param_phys(mh,mH,mA,mHp,sba,lambda6,lambda7,m12_2,tan_beta);
  double cba = get_cba();
  double l1,l2,l3,l4,l5,l6,l7;
  get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_2,tan_beta);

  fprintf(output,"##################################################################\n");
  fprintf(output,"#                                                                #\n");
  fprintf(output,"#                 Two-Higgs Doublet Model Output                 #\n");
  fprintf(output,"#          Produced by 2HDMC: http://2hdmc.hepforge.org          #\n");
  fprintf(output,"#                                                                #\n");
  fprintf(output,"##################################################################\n");
  fprintf(output,"Block MODSEL # Select Model\n");
  fprintf(output,"    0   10    #  10 = 2HDM\n");
  fprintf(output,"Block FMODSEL # FLHA definitions, see [arXiv:1008.0762]\n");
  fprintf(output,"    1   %2d    #  Model\n", 30+yukawas_type);
  fprintf(output,"    5    0    #  No CP-violation\n");
  fprintf(output,"Block THDM\n");
  fprintf(output,"    1   %2d    #  Valid parameters (1=Yes, 0=no)\n", (mh>0)&&(mH>0)&&(mA>0)&&(mHp>0));
  fprintf(output,"    2   %2d    #  Tree-level unitarity (1=Yes, 0=no)\n", (int)this->check_unitarity());
  fprintf(output,"    3   %2d    #  Perturbativity (1=Yes, 0=no)\n", (int)this->check_perturbativity(4*M_PI));
  fprintf(output,"    4   %2d    #  Stability (1=Yes, 0=no)\n", (int)this->check_stability());
  fprintf(output,"Block SMINPUTS  # Standard Model inputs\n");
  fprintf(output,"    1      % 16.8e   # 1/alpha_em(MZ) SM MSbar\n",1./sm.get_alpha());
  fprintf(output,"    2      % 16.8e   # G Fermi\n",sm.get_GF());
  fprintf(output,"    3      % 16.8e   # alpha_s(MZ) SM MSbar\n",sm.get_alpha_s());
  fprintf(output,"    4      % 16.8e   # MZ\n",sm.get_MZ());
  fprintf(output,"    5      % 16.8e   # mb(mb)\n",sm.get_qmass_MSbar(5));
  fprintf(output,"    6      % 16.8e   # mt (pole)\n",sm.get_qmass_pole(6));
  fprintf(output,"    7      % 16.8e   # mtau(pole)\n",sm.get_lmass_pole(3));
  fprintf(output,"Block GAUGE  # SM Gauge couplings\n");
  fprintf(output,"    1      % 16.8e   # g'\n", sm.get_gprime());
  fprintf(output,"    2      % 16.8e   # g\n", sm.get_g());
  fprintf(output,"    3      % 16.8e   # g_3\n", 4.*M_PI*sm.get_alpha_s());
fprintf(output,"Block MINPAR    # Model parameters\n");
  fprintf(output,"# Parameters for general potential ");
  if (tan_beta!=0) {
    fprintf(output,"in generic basis\n");
  } else {
    fprintf(output,"in Higgs basis\n");
  }
  fprintf(output,"    3      % 16.8e   # tan(beta)\n",tan_beta);
  fprintf(output,"   11      % 16.8e   # lambda_1\n",l1);
  fprintf(output,"   12      % 16.8e   # lambda_2\n",l2);
  fprintf(output,"   13      % 16.8e   # lambda_3\n",l3);
  fprintf(output,"   14      % 16.8e   # lambda_4\n",l4);
  fprintf(output,"   15      % 16.8e   # lambda_5\n",l5);
  fprintf(output,"   16      % 16.8e   # lambda_6\n",l6);
  fprintf(output,"   17      % 16.8e   # lambda_7\n",l7);
  if (tan_beta!=0) {
    fprintf(output,"   18      % 16.8e   # m_12^2\n",m12_2);
  }
  fprintf(output,"   20      % 16.8e   # sin(beta-alpha)\n",sba);
  fprintf(output,"   21      % 16.8e   # cos(beta-alpha)\n",cba);

  if (yukawas_type>0) {
    fprintf(output,"   24     %1i                   # Yukawas Type\n",yukawas_type);
  }
  fprintf(output,"Block MGCKM     # CKM elements\n");
  fprintf(output,"    1     1      % 16.8e   # Vud\n",sm.get_CKM_element(1,1));
  fprintf(output,"    1     2      % 16.8e   # Vus\n",sm.get_CKM_element(1,2));
  fprintf(output,"    1     3      % 16.8e   # Vub\n",sm.get_CKM_element(1,3));
  fprintf(output,"    2     1      % 16.8e   # Vcd\n",sm.get_CKM_element(2,1));
  fprintf(output,"    2     2      % 16.8e   # Vcs\n",sm.get_CKM_element(2,2));
  fprintf(output,"    2     3      % 16.8e   # Vcb\n",sm.get_CKM_element(2,3));
  fprintf(output,"    3     1      % 16.8e   # Vtd\n",sm.get_CKM_element(3,1));
  fprintf(output,"    3     2      % 16.8e   # Vts\n",sm.get_CKM_element(3,2));
  fprintf(output,"    3     3      % 16.8e   # Vtb\n",sm.get_CKM_element(3,3));
  fprintf(output,"Block MASS      #  Mass spectrum (kinematic masses)\n");
  fprintf(output,"#  PDG      Mass\n");
  fprintf(output,"     1      % 16.8e   # Md\n",sm.get_qmass_pole(1));
  fprintf(output,"     2      % 16.8e   # Mu\n",sm.get_qmass_pole(2));
  fprintf(output,"     3      % 16.8e   # Ms\n",sm.get_qmass_pole(3));
  fprintf(output,"     4      % 16.8e   # Mc\n",sm.get_qmass_pole(4));
  fprintf(output,"     5      % 16.8e   # Mb\n",sm.get_qmass_pole(5));
  fprintf(output,"     6      % 16.8e   # Mt\n",sm.get_qmass_pole(6));
  fprintf(output,"    11      % 16.8e   # Me\n",sm.get_lmass_pole(1));
  fprintf(output,"    13      % 16.8e   # Mmu\n",sm.get_lmass_pole(2));
  fprintf(output,"    15      % 16.8e   # Mta\n",sm.get_lmass_pole(3));
  fprintf(output,"    23      % 16.8e   # MZ\n",sm.get_MZ());
  fprintf(output,"    24      % 16.8e   # MW\n",sm.get_MW());
  fprintf(output,"    25      % 16.8e   # Mh1, lightest CP-even Higgs\n",mh);
  fprintf(output,"    35      % 16.8e   # Mh2, heaviest CP-even Higgs\n",mH);
  fprintf(output,"    36      % 16.8e   # Mh3, CP-odd Higgs\n",mA);
  fprintf(output,"    37      % 16.8e   # Mhc\n",mHp);

  fprintf(output,"Block ALPHA     # Effective Higgs mixing parameter\n");
  fprintf(output,"            % 16.8e   # alpha\n",-asin(sba)+atan(tan_beta));

  double lu[4][4],ld[4][4],ll[4][4];
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      lu[i][j]=0.;
      ld[i][j]=0.;
      ll[i][j]=0.;
    }
  }
  double k1,k2,k3,r1,r2,r3;
  get_kappa_up(k1,k2,k3);
  get_yukawas_up(r1,r2,r3);
  (k1>0 ? lu[1][1] = r1/k1 : lu[1][1]=0.);
  (k2>0 ? lu[2][2] = r2/k2 : lu[2][2]=0.);
  (k3>0 ? lu[3][3] = r3/k3 : lu[3][3]=0.);

  get_kappa_down(k1,k2,k3);
  get_yukawas_down(r1,r2,r3);
  (k1>0 ? ld[1][1] = r1/k1 : ld[1][1]=0.);
  (k2>0 ? ld[2][2] = r2/k2 : ld[2][2]=0.);
  (k3>0 ? ld[3][3] = r3/k3 : ld[3][3]=0.);

  get_kappa_lepton(k1,k2,k3);
  get_yukawas_lepton(r1,r2,r3);
  (k1>0 ? ll[1][1] = r1/k1 : ll[1][1]=0.);
  (k2>0 ? ll[2][2] = r2/k2 : ll[2][2]=0.);
  (k3>0 ? ll[3][3] = r3/k3 : ll[3][3]=0.);


  fprintf(output,"Block UCOUPL\n");
  for (int i=1;i<=3;i++) {
    for (int j=1;j<=3;j++) {
       fprintf(output,"%5d%6d   % 16.8e   # LU_{%d,%d}\n",i,j,lu[i][j],i,j);
    }
  }

  fprintf(output,"Block DCOUPL\n");
  for (int i=1;i<=3;i++) {
    for (int j=1;j<=3;j++) {
       fprintf(output,"%5d%6d   % 16.8e   # LD_{%d,%d}\n",i,j,ld[i][j],i,j);
    }
  }

  fprintf(output,"Block LCOUPL\n");
  for (int i=1;i<=3;i++) {
    for (int j=1;j<=3;j++) {
       fprintf(output,"%5d%6d   % 16.8e   # LL_{%d,%d}\n",i,j,ll[i][j],i,j);
    }
  }

#if defined HiggsBounds

  if (hbhs) {

  fprintf(output, "Block HBRESULT\n");
  fprintf(output, "# Higgs  Result  Channel    Obsratio      Ncomb\n");
  fprintf(output, "%5d   %5d   %6d %16.8E %5d    # Full result\n", 0, hbhs->hb.result[0],hbhs->hb.channel[0],hbhs->hb.obsratio[0],hbhs->hb.ncombined[0]);
  fprintf(output, "%5d   %5d   %6d %16.8E %5d    # h\n",  1, hbhs->hb.result[1],hbhs->hb.channel[1],hbhs->hb.obsratio[1],hbhs->hb.ncombined[1]);
  fprintf(output, "%5d   %5d   %6d %16.8E %5d    # H\n",  2, hbhs->hb.result[2],hbhs->hb.channel[2],hbhs->hb.obsratio[2],hbhs->hb.ncombined[2]);
  fprintf(output, "%5d   %5d   %6d %16.8E %5d    # A\n",  3, hbhs->hb.result[3],hbhs->hb.channel[3],hbhs->hb.obsratio[3],hbhs->hb.ncombined[3]);
  fprintf(output, "%5d   %5d   %6d %16.8E %5d    # H+\n", 4, hbhs->hb.result[4],hbhs->hb.channel[4],hbhs->hb.obsratio[4],hbhs->hb.ncombined[4]);
  fprintf(output, "Block HSRESULT\n");
  fprintf(output, " %5d %16.8E # Total HS chi^2 \n", 0, hbhs->hs.chisq);
  fprintf(output, " %5d %16.8E # chi^2 from rates\n", 1, hbhs->hs.chisq_mu);
  fprintf(output, " %5d %16.8E # chi^2 from mass\n", 2, hbhs->hs.chisq_mh);
  fprintf(output, " %5d %16d # Number of observables\n", 3, hbhs->hs.nobs);

  }

#endif


  DecayTable table(*this);
  table.set_qcd(qcd_on);
  fprintf(output,"#     PDG   Width\n");
  table.print_top_decays_LesHouches(output,fulldecay);
  table.print_decays_LesHouches(output,1,fulldecay);
  table.print_decays_LesHouches(output,2,fulldecay);
  table.print_decays_LesHouches(output,3,fulldecay);
  table.print_decays_LesHouches(output,4,fulldecay);
  fprintf(output,"#\n");



  if (couplings) {

    complex <double> c1;
    complex <double> c2;

    fprintf(output,"BLOCK MGUSER\n");
    get_coupling_hll(1,1,1,c1,c2);
    fprintf(output,"         1     % 16.8e   # REPLGH1EE   , Real part of scalar part of h1ee coupling\n",real(c1));
    fprintf(output,"         2     % 16.8e   # IMPLGH1EE   , Imaginary part of scalar part of h1ee coupling\n",imag(c1));
    fprintf(output,"         3     % 16.8e   # REPRGH1EE   , Real part of pseudo scalar part of h1ee coupling\n",real(c2));
    fprintf(output,"         4     % 16.8e   # IMPRGH1EE   , Imaginary part of pseudo scalar part of h1ee coupling\n",imag(c2));
    get_coupling_hll(2,1,1,c1,c2);
    fprintf(output,"         5     % 16.8e   # REPLGH2EE   , Real part of scalar part of h2ee coupling\n",real(c1));
    fprintf(output,"         6     % 16.8e   # IMPLGH2EE   , Imaginary part of scalar part of h2ee coupling\n",imag(c1));
    fprintf(output,"         7     % 16.8e   # REPRGH2EE   , Real part of pseudo scalar part of h2ee coupling\n",real(c2));
    fprintf(output,"         8     % 16.8e   # IMPRGH2EE   , Imaginary part of pseudo scalar part of h2ee coupling\n",imag(c2));
    get_coupling_hll(3,1,1,c1,c2);
    fprintf(output,"         9     % 16.8e   # REPLGH3EE   , Real part of scalar part of h3ee coupling\n",real(c1));
    fprintf(output,"        10     % 16.8e   # IMPLGH3EE   , Imaginary part of scalar part of h3ee coupling\n",imag(c1));
    fprintf(output,"        11     % 16.8e   # REPRGH3EE   , Real part of pseudo scalar part of h3ee coupling\n",real(c2));
    fprintf(output,"        12     % 16.8e   # IMPRGH3EE   , Imaginary part of pseudo scalar part of h3ee coupling\n",imag(c2));

    get_coupling_hll(1,2,2,c1,c2);
    fprintf(output,"        13     % 16.8e   # REPLGH1MUMU , Real part of scalar part of h1mumu coupling\n",real(c1));
    fprintf(output,"        14     % 16.8e   # IMPLGH1MUMU , Imaginary part of scalar part of h1mumu coupling\n",imag(c1));
    fprintf(output,"        15     % 16.8e   # REPRGH1MUMU , Real part of pseudo scalar part of h1mumu coupling\n",real(c2));
    fprintf(output,"        16     % 16.8e   # IMPRGH1MUMU , Imaginary part of pseudo scalar part of h1mumu coupling\n",imag(c2));
    get_coupling_hll(2,2,2,c1,c2);
    fprintf(output,"        17     % 16.8e   # REPLGH2MUMU , Real part of scalar part of h2mumu coupling\n",real(c1));
    fprintf(output,"        18     % 16.8e   # IMPLGH2MUMU , Imaginary part of scalar part of h2mumu coupling\n",imag(c1));
    fprintf(output,"        19     % 16.8e   # REPRGH2MUMU , Real part of pseudo scalar part of h2mumu coupling\n",real(c2));
    fprintf(output,"        20     % 16.8e   # IMPRGH2MUMU , Imaginary part of pseudo scalar part of h2mumu coupling\n",imag(c2));
    get_coupling_hll(3,2,2,c1,c2);
    fprintf(output,"        21     % 16.8e   # REPLGH3MUMU , Real part of scalar part of h3mumu coupling\n",real(c1));
    fprintf(output,"        22     % 16.8e   # IMPLGH3MUMU , Imaginary part of scalar part of h3mumu coupling\n",imag(c1));
    fprintf(output,"        23     % 16.8e   # REPRGH3MUMU , Real part of pseudo scalar part of h3mumu coupling\n",real(c2));
    fprintf(output,"        24     % 16.8e   # IMPRGH3MUMU , Imaginary part of pseudo scalar part of h3mumu coupling\n",imag(c2));

    get_coupling_hll(1,3,3,c1,c2);
    fprintf(output,"        25     % 16.8e   # REPLGH1TATA , Real part of scalar part of h1tata coupling\n",real(c1));
    fprintf(output,"        26     % 16.8e   # IMPLGH1TATA , Imaginary part of scalar part of h1tata coupling\n",imag(c1));
    fprintf(output,"        27     % 16.8e   # REPRGH1TATA , Real part of pseudo scalar part of h1tata coupling\n",real(c2));
    fprintf(output,"        28     % 16.8e   # IMPRGH1TATA , Imaginary part of pseudo scalar part of h1tata coupling\n",imag(c2));
    get_coupling_hll(2,3,3,c1,c2);
    fprintf(output,"        29     % 16.8e   # REPLGH2TATA , Real part of scalar part of h2tata coupling\n",real(c1));
    fprintf(output,"        30     % 16.8e   # IMPLGH2TATA , Imaginary part of scalar part of h2tata coupling\n",imag(c1));
    fprintf(output,"        31     % 16.8e   # REPRGH2TATA , Real part of pseudo scalar part of h2tata coupling\n",real(c2));
    fprintf(output,"        32     % 16.8e   # IMPRGH2TATA , Imaginary part of pseudo scalar part of h2tata coupling\n",imag(c2));
    get_coupling_hll(3,3,3,c1,c2);
    fprintf(output,"        33     % 16.8e   # REPLGH3TATA , Real part of scalar part of h3tata coupling\n",real(c1));
    fprintf(output,"        34     % 16.8e   # IMPLGH3TATA , Imaginary part of scalar part of h3tata coupling\n",imag(c1));
    fprintf(output,"        35     % 16.8e   # REPRGH3TATA , Real part of pseudo scalar part of h3tata coupling\n",real(c2));
    fprintf(output,"        36     % 16.8e   # IMPRGH3TATA , Imaginary part of pseudo scalar part of h3tata coupling\n",imag(c2));

    get_coupling_huu(1,1,1,c1,c2);
    fprintf(output,"        37     % 16.8e   # REPLGH1UU   , Real part of scalar part of h1uu coupling\n",real(c1));
    fprintf(output,"        38     % 16.8e   # IMPLGH1UU   , Imaginary part of scalar part of h1uu coupling\n",imag(c1));
    fprintf(output,"        39     % 16.8e   # REPRGH1UU   , Real part of pseudo scalar part of h1uu coupling\n",real(c2));
    fprintf(output,"        40     % 16.8e   # IMPRGH1UU   , Imaginary part of pseudo scalar part of h1uu coupling\n",imag(c2));
    get_coupling_huu(2,1,1,c1,c2);
    fprintf(output,"        41     % 16.8e   # REPLGH2UU   , Real part of scalar part of h2uu coupling\n",real(c1));
    fprintf(output,"        42     % 16.8e   # IMPLGH2UU   , Imaginary part of scalar part of h2uu coupling\n",imag(c1));
    fprintf(output,"        43     % 16.8e   # REPRGH2UU   , Real part of pseudo scalar part of h2uu coupling\n",real(c2));
    fprintf(output,"        44     % 16.8e   # IMPRGH2UU   , Imaginary part of pseudo scalar part of h2uu coupling\n",imag(c2));
    get_coupling_huu(3,1,1,c1,c2);
    fprintf(output,"        45     % 16.8e   # REPLGH3UU   , Real part of scalar part of h3uu coupling\n",real(c1));
    fprintf(output,"        46     % 16.8e   # IMPLGH3UU   , Imaginary part of scalar part of h3uu coupling\n",imag(c1));
    fprintf(output,"        47     % 16.8e   # REPRGH3UU   , Real part of pseudo scalar part of h3uu coupling\n",real(c2));
    fprintf(output,"        48     % 16.8e   # IMPRGH3UU   , Imaginary part of pseudo scalar part of h3uu coupling\n",imag(c2));

    get_coupling_huu(1,2,2,c1,c2);
    fprintf(output,"        49     % 16.8e   # REPLGH1CC   , Real part of scalar part of h1cc coupling\n",real(c1));
    fprintf(output,"        50     % 16.8e   # IMPLGH1CC   , Imaginary part of scalar part of h1cc coupling\n",imag(c1));
    fprintf(output,"        51     % 16.8e   # REPRGH1CC   , Real part of pseudo scalar part of h1cc coupling\n",real(c2));
    fprintf(output,"        52     % 16.8e   # IMPRGH1CC   , Imaginary part of pseudo scalar part of h1cc coupling\n",imag(c2));
    get_coupling_huu(2,2,2,c1,c2);
    fprintf(output,"        53     % 16.8e   # REPLGH2CC   , Real part of scalar part of h2cc coupling\n",real(c1));
    fprintf(output,"        54     % 16.8e   # IMPLGH2CC   , Imaginary part of scalar part of h2cc coupling\n",imag(c1));
    fprintf(output,"        55     % 16.8e   # REPRGH2CC   , Real part of pseudo scalar part of h2cc coupling\n",real(c2));
    fprintf(output,"        56     % 16.8e   # IMPRGH2CC   , Imaginary part of pseudo scalar part of h2cc coupling\n",imag(c2));
    get_coupling_huu(3,2,2,c1,c2);
    fprintf(output,"        57     % 16.8e   # REPLGH3CC   , Real part of scalar part of h3cc coupling\n",real(c1));
    fprintf(output,"        58     % 16.8e   # IMPLGH3CC   , Imaginary part of scalar part of h3cc coupling\n",imag(c1));
    fprintf(output,"        59     % 16.8e   # REPRGH3CC   , Real part of pseudo scalar part of h3cc coupling\n",real(c2));
    fprintf(output,"        60     % 16.8e   # IMPRGH3CC   , Imaginary part of pseudo scalar part of h3cc coupling\n",imag(c2));

    get_coupling_huu(1,3,3,c1,c2);
    fprintf(output,"        61     % 16.8e   # REPLGH1TT   , Real part of scalar part of h1tt coupling\n",real(c1));
    fprintf(output,"        62     % 16.8e   # IMPLGH1TT   , Imaginary part of scalar part of h1tt coupling\n",imag(c1));
    fprintf(output,"        63     % 16.8e   # REPRGH1TT   , Real part of pseudo scalar part of h1tt coupling\n",real(c2));
    fprintf(output,"        64     % 16.8e   # IMPRGH1TT   , Imaginary part of pseudo scalar part of h1tt coupling\n",imag(c2));
    get_coupling_huu(2,3,3,c1,c2);
    fprintf(output,"        65     % 16.8e   # REPLGH2TT   , Real part of scalar part of h2tt coupling\n",real(c1));
    fprintf(output,"        66     % 16.8e   # IMPLGH2TT   , Imaginary part of scalar part of h2tt coupling\n",imag(c1));
    fprintf(output,"        67     % 16.8e   # REPRGH2TT   , Real part of pseudo scalar part of h2tt coupling\n",real(c2));
    fprintf(output,"        68     % 16.8e   # IMPRGH2TT   , Imaginary part of pseudo scalar part of h2tt coupling\n",imag(c2));
    get_coupling_huu(3,3,3,c1,c2);
    fprintf(output,"        69     % 16.8e   # REPLGH3TT   , Real part of scalar part of h3tt coupling\n",real(c1));
    fprintf(output,"        70     % 16.8e   # IMPLGH3TT   , Imaginary part of scalar part of h3tt coupling\n",imag(c1));
    fprintf(output,"        71     % 16.8e   # REPRGH3TT   , Real part of pseudo scalar part of h3tt coupling\n",real(c2));
    fprintf(output,"        72     % 16.8e   # IMPRGH3TT   , Imaginary part of pseudo scalar part of h3tt coupling\n",imag(c2));

    get_coupling_hdd(1,1,1,c1,c2);
    fprintf(output,"        73     % 16.8e   # REPLGH1DD   , Real part of scalar part of h1dd coupling\n",real(c1));
    fprintf(output,"        74     % 16.8e   # IMPLGH1DD   , Imaginary part of scalar part of h1dd coupling\n",imag(c1));
    fprintf(output,"        75     % 16.8e   # REPRGH1DD   , Real part of pseudo scalar part of h1dd coupling\n",real(c2));
    fprintf(output,"        76     % 16.8e   # IMPRGH1DD   , Imaginary part of pseudo scalar part of h1dd coupling\n",imag(c2));
    get_coupling_hdd(2,1,1,c1,c2);
    fprintf(output,"        77     % 16.8e   # REPLGH2DD   , Real part of scalar part of h2dd coupling\n",real(c1));
    fprintf(output,"        78     % 16.8e   # IMPLGH2DD   , Imaginary part of scalar part of h2dd coupling\n",imag(c1));
    fprintf(output,"        79     % 16.8e   # REPRGH2DD   , Real part of pseudo scalar part of h2dd coupling\n",real(c2));
    fprintf(output,"        80     % 16.8e   # IMPRGH2DD   , Imaginary part of pseudo scalar part of h2dd coupling\n",imag(c2));
    get_coupling_hdd(3,1,1,c1,c2);
    fprintf(output,"        81     % 16.8e   # REPLGH3DD   , Real part of scalar part of h3dd coupling\n",real(c1));
    fprintf(output,"        82     % 16.8e   # IMPLGH3DD   , Imaginary part of scalar part of h3dd coupling\n",imag(c1));
    fprintf(output,"        83     % 16.8e   # REPRGH3DD   , Real part of pseudo scalar part of h3dd coupling\n",real(c2));
    fprintf(output,"        84     % 16.8e   # IMPRGH3DD   , Imaginary part of pseudo scalar part of h3dd coupling\n",imag(c2));

    get_coupling_hdd(1,2,2,c1,c2);
    fprintf(output,"        85     % 16.8e   # REPLGH1SS   , Real part of scalar part of h1ss coupling\n",real(c1));
    fprintf(output,"        86     % 16.8e   # IMPLGH1SS   , Imaginary part of scalar part of h1ss coupling\n",imag(c1));
    fprintf(output,"        87     % 16.8e   # REPRGH1SS   , Real part of pseudo scalar part of h1ss coupling\n",real(c2));
    fprintf(output,"        88     % 16.8e   # IMPRGH1SS   , Imaginary part of pseudo scalar part of h1ss coupling\n",imag(c2));
    get_coupling_hdd(2,2,2,c1,c2);
    fprintf(output,"        89     % 16.8e   # REPLGH2SS   , Real part of scalar part of h2ss coupling\n",real(c1));
    fprintf(output,"        90     % 16.8e   # IMPLGH2SS   , Imaginary part of scalar part of h2ss coupling\n",imag(c1));
    fprintf(output,"        91     % 16.8e   # REPRGH2SS   , Real part of pseudo scalar part of h2ss coupling\n",real(c2));
    fprintf(output,"        92     % 16.8e   # IMPRGH2SS   , Imaginary part of pseudo scalar part of h2ss coupling\n",imag(c2));
    get_coupling_hdd(3,2,2,c1,c2);
    fprintf(output,"        93     % 16.8e   # REPLGH3SS   , Real part of scalar part of h3ss coupling\n",real(c1));
    fprintf(output,"        94     % 16.8e   # IMPLGH3SS   , Imaginary part of scalar part of h3ss coupling\n",imag(c1));
    fprintf(output,"        95     % 16.8e   # REPRGH3SS   , Real part of pseudo scalar part of h3ss coupling\n",real(c2));
    fprintf(output,"        96     % 16.8e   # IMPRGH3SS   , Imaginary part of pseudo scalar part of h3ss coupling\n",imag(c2));

    get_coupling_hdd(1,3,3,c1,c2);
    fprintf(output,"        97     % 16.8e   # REPLGH1BB   , Real part of scalar part of h1bb coupling\n",real(c1));
    fprintf(output,"        98     % 16.8e   # IMPLGH1BB   , Imaginary part of scalar part of h1bb coupling\n",imag(c1));
    fprintf(output,"        99     % 16.8e   # REPRGH1BB   , Real part of pseudo scalar part of h1bb coupling\n",real(c2));
    fprintf(output,"       100     % 16.8e   # IMPRGH1BB   , Imaginary part of pseudo scalar part of h1bb coupling\n",imag(c2));
    get_coupling_hdd(2,3,3,c1,c2);
    fprintf(output,"       101     % 16.8e   # REPLGH2BB   , Real part of scalar part of h2bb coupling\n",real(c1));
    fprintf(output,"       102     % 16.8e   # IMPLGH2BB   , Imaginary part of scalar part of h2bb coupling\n",imag(c1));
    fprintf(output,"       103     % 16.8e   # REPRGH2BB   , Real part of pseudo scalar part of h2bb coupling\n",real(c2));
    fprintf(output,"       104     % 16.8e   # IMPRGH2BB   , Imaginary part of pseudo scalar part of h2bb coupling\n",imag(c2));
    get_coupling_hdd(3,3,3,c1,c2);
    fprintf(output,"       105     % 16.8e   # REPLGH3BB   , Real part of scalar part of h3bb coupling\n",real(c1));
    fprintf(output,"       106     % 16.8e   # IMPLGH3BB   , Imaginary part of scalar part of h3bb coupling\n",imag(c1));
    fprintf(output,"       107     % 16.8e   # REPRGH3BB   , Real part of pseudo scalar part of h3bb coupling\n",real(c2));
    fprintf(output,"       108     % 16.8e   # IMPRGH3BB   , Imaginary part of pseudo scalar part of h3bb coupling\n",imag(c2));

    get_coupling_hdu(4,1,1,c1,c2);
    fprintf(output,"       109     % 16.8e   # REPLGHCUD   , Real part of scalar part of hcud coupling\n",real(c1));
    fprintf(output,"       110     % 16.8e   # IMPLGHCUD   , Imaginary part of scalar part of hcud coupling\n",imag(c1));
    fprintf(output,"       111     % 16.8e   # REPRGHCUD   , Real part of pseudo scalar part of hcud coupling\n",real(c2));
    fprintf(output,"       112     % 16.8e   # IMPRGHCUD   , Imaginary part of pseudo scalar part of hcud coupling\n",imag(c2));
    get_coupling_hdu(4,2,1,c1,c2);
    fprintf(output,"       113     % 16.8e   # REPLGHCUS   , Real part of scalar part of hcus coupling\n",real(c1));
    fprintf(output,"       114     % 16.8e   # IMPLGHCUS   , Imaginary part of scalar part of hcus coupling\n",imag(c1));
    fprintf(output,"       115     % 16.8e   # REPRGHCUS   , Real part of pseudo scalar part of hcus coupling\n",real(c2));
    fprintf(output,"       116     % 16.8e   # IMPRGHCUS   , Imaginary part of pseudo scalar part of hcus coupling\n",imag(c2));
    get_coupling_hdu(4,3,1,c1,c2);
    fprintf(output,"       117     % 16.8e   # REPLGHCUB   , Real part of scalar part of hcub coupling\n",real(c1));
    fprintf(output,"       118     % 16.8e   # IMPLGHCUB   , Imaginary part of scalar part of hcub coupling\n",imag(c1));
    fprintf(output,"       119     % 16.8e   # REPRGHCUB   , Real part of pseudo scalar part of hcub coupling\n",real(c2));
    fprintf(output,"       120     % 16.8e   # IMPRGHCUB   , Imaginary part of pseudo scalar part of hcub coupling\n",imag(c2));

    get_coupling_hdu(4,1,2,c1,c2);
    fprintf(output,"       121     % 16.8e   # REPLGHCCD   , Real part of scalar part of hccd coupling\n",real(c1));
    fprintf(output,"       122     % 16.8e   # IMPLGHCCD   , Imaginary part of scalar part of hccd coupling\n",imag(c1));
    fprintf(output,"       123     % 16.8e   # REPRGHCCD   , Real part of pseudo scalar part of hccd coupling\n",real(c2));
    fprintf(output,"       124     % 16.8e   # IMPRGHCCD   , Imaginary part of pseudo scalar part of hccd coupling\n",imag(c2));
    get_coupling_hdu(4,2,2,c1,c2);
    fprintf(output,"       125     % 16.8e   # REPLGHCCS   , Real part of scalar part of hccs coupling\n",real(c1));
    fprintf(output,"       126     % 16.8e   # IMPLGHCCS   , Imaginary part of scalar part of hccs coupling\n",imag(c1));
    fprintf(output,"       127     % 16.8e   # REPRGHCCS   , Real part of pseudo scalar part of hccs coupling\n",real(c2));
    fprintf(output,"       128     % 16.8e   # IMPRGHCCS   , Imaginary part of pseudo scalar part of hccs coupling\n",imag(c2));
    get_coupling_hdu(4,3,2,c1,c2);
    fprintf(output,"       129     % 16.8e   # REPLGHCCB   , Real part of scalar part of hccb coupling\n",real(c1));
    fprintf(output,"       130     % 16.8e   # IMPLGHCCB   , Imaginary part of scalar part of hccb coupling\n",imag(c1));
    fprintf(output,"       131     % 16.8e   # REPRGHCCB   , Real part of pseudo scalar part of hccb coupling\n",real(c2));
    fprintf(output,"       132     % 16.8e   # IMPRGHCCB   , Imaginary part of pseudo scalar part of hccb coupling\n",imag(c2));

    get_coupling_hdu(4,1,3,c1,c2);
    fprintf(output,"       133     % 16.8e   # REPLGHCTD   , Real part of scalar part of hctd coupling\n",real(c1));
    fprintf(output,"       134     % 16.8e   # IMPLGHCTD   , Imaginary part of scalar part of hctd coupling\n",imag(c1));
    fprintf(output,"       135     % 16.8e   # REPRGHCTD   , Real part of pseudo scalar part of hctd coupling\n",real(c2));
    fprintf(output,"       136     % 16.8e   # IMPRGHCTD   , Imaginary part of pseudo scalar part of hctd coupling\n",imag(c2));
    get_coupling_hdu(4,2,3,c1,c2);
    fprintf(output,"       137     % 16.8e   # REPLGHCTS   , Real part of scalar part of hcts coupling\n",real(c1));
    fprintf(output,"       138     % 16.8e   # IMPLGHCTS   , Imaginary part of scalar part of hcts coupling\n",imag(c1));
    fprintf(output,"       139     % 16.8e   # REPRGHCTS   , Real part of pseudo scalar part of hcts coupling\n",real(c2));
    fprintf(output,"       140     % 16.8e   # IMPRGHCTS   , Imaginary part of pseudo scalar part of hcts coupling\n",imag(c2));
    get_coupling_hdu(4,3,3,c1,c2);
    fprintf(output,"       141     % 16.8e   # REPLGHCTB   , Real part of scalar part of hctb coupling\n",real(c1));
    fprintf(output,"       142     % 16.8e   # IMPLGHCTB   , Imaginary part of scalar part of hctb coupling\n",imag(c1));
    fprintf(output,"       143     % 16.8e   # REPRGHCTB   , Real part of pseudo scalar part of hctb coupling\n",real(c2));
    fprintf(output,"       144     % 16.8e   # IMPRGHCTB   , Imaginary part of pseudo scalar part of hctb coupling\n",imag(c2));

    get_coupling_hln(4,1,1,c1,c2);
    fprintf(output,"       145     % 16.8e   # REPLGHCVEE  , Real part of scalar part of hcvee coupling\n",real(c1));
    fprintf(output,"       146     % 16.8e   # IMPLGHCVEE  , Imaginary part of scalar part of hcvee coupling\n",imag(c1));
    fprintf(output,"       147     % 16.8e   # REPRGHCVEE  , Real part of pseudo scalar part of hcvee coupling\n",real(c2));
    fprintf(output,"       148     % 16.8e   # IMPRGHCVEE  , Imaginary part of pseudo scalar part of hcvee coupling\n",imag(c2));
    get_coupling_hln(4,2,2,c1,c2);
    fprintf(output,"       149     % 16.8e   # REPLGHCVMMU , Real part of scalar part of hcvmmu coupling\n",real(c1));
    fprintf(output,"       150     % 16.8e   # IMPLGHCVMMU , Imaginary part of scalar part of hcvmmu coupling\n",imag(c1));
    fprintf(output,"       151     % 16.8e   # REPRGHCVMMU , Real part of pseudo scalar part of hcvmmu coupling\n",real(c2));
    fprintf(output,"       152     % 16.8e   # IMPRGHCVMMU , Imaginary part of pseudo scalar part of hcvmmu coupling\n",imag(c2));
    get_coupling_hln(4,3,3,c1,c2);
    fprintf(output,"       153     % 16.8e   # REPLGHCVTTA , Real part of scalar part of hcvtta coupling\n",real(c1));
    fprintf(output,"       154     % 16.8e   # IMPLGHCVTTA , Imaginary part of scalar part of hcvtta coupling\n",imag(c1));
    fprintf(output,"       155     % 16.8e   # REPRGHCVTTA , Real part of pseudo scalar part of hcvtta coupling\n",real(c2));
    fprintf(output,"       156     % 16.8e   # IMPRGHCVTTA , Imaginary part of pseudo scalar part of hcvtta coupling\n",imag(c2));

    get_coupling_vvh(3,3,1,c1);
    fprintf(output,"       157     % 16.8e   # REGWWH1     , Real part of wwh1 coupling\n",real(c1));
    fprintf(output,"       158     % 16.8e   # IMGWWH1     , Imaginary part of wwh1 coupling\n",imag(c1));
    get_coupling_vvh(3,3,2,c1);
    fprintf(output,"       159     % 16.8e   # REGWWH2     , Real part of wwh2 coupling\n",real(c1));
    fprintf(output,"       160     % 16.8e   # IMGWWH2     , Imaginary part of wwh2 coupling\n",imag(c1));

    get_coupling_vvh(2,2,1,c1);
    fprintf(output,"       161     % 16.8e   # REGZZH1     , Real part of zzh1 coupling\n",real(c1));
    fprintf(output,"       162     % 16.8e   # IMGZZH1     , Imaginary part of zzh1 coupling\n",imag(c1));
    get_coupling_vvh(2,2,2,c1);
    fprintf(output,"       163     % 16.8e   # REGZZH2     , Real part of zzh2 coupling\n",real(c1));
    fprintf(output,"       164     % 16.8e   # IMGZZH2     , Imaginary part of zzh2 coupling\n",imag(c1));

    get_coupling_vhh(1,4,4,c1);
    fprintf(output,"       165     % 16.8e   # REGAHCHC    , Real part of ahchc coupling\n",real(c1));
    fprintf(output,"       166     % 16.8e   # IMGAHCHC    , Imaginary part of ahchc coupling\n",imag(c1));

    get_coupling_vhh(2,1,2,c1);
    fprintf(output,"       167     % 16.8e   # REGZH1H2    , Real part of zh1h2 coupling\n",real(c1));
    fprintf(output,"       168     % 16.8e   # IMGZH1H2    , Imaginary part of zh1h2 coupling\n",imag(c1));
    get_coupling_vhh(2,1,3,c1);
    fprintf(output,"       169     % 16.8e   # REGZH1H3    , Real part of zh1h3 coupling\n",real(c1));
    fprintf(output,"       170     % 16.8e   # IMGZH1H3    , Imaginary part of zh1h3 coupling\n",imag(c1));
    get_coupling_vhh(2,2,3,c1);
    fprintf(output,"       171     % 16.8e   # REGZH2H3    , Real part of zh2h3 coupling\n",real(c1));
    fprintf(output,"       172     % 16.8e   # IMGZH2H3    , Imaginary part of zh2h3 coupling\n",imag(c1));
    get_coupling_vhh(2,4,4,c1);
    fprintf(output,"       173     % 16.8e   # REGZHCHC    , Real part of zhchc coupling\n",real(c1));
    fprintf(output,"       174     % 16.8e   # IMGZHCHC    , Imaginary part of zhchc coupling\n",imag(c1));

    get_coupling_vhh(3,1,4,c1);
    fprintf(output,"       175     % 16.8e   # REGWPHCH1   , Real part of wph1hc coupling\n",real(c1));
    fprintf(output,"       176     % 16.8e   # IMGWPHCH1   , Imaginary part of wph1hc coupling\n",imag(c1));
    get_coupling_vhh(3,2,4,c1);
    fprintf(output,"       177     % 16.8e   # REGWPHCH2   , Real part of wph2hc coupling\n",real(c1));
    fprintf(output,"       178     % 16.8e   # IMGWPHCH2   , Imaginary part of wph2hc coupling\n",imag(c1));
    get_coupling_vhh(3,3,4,c1);
    fprintf(output,"       179     % 16.8e   # REGWPHCH3   , Real part of wph3hc coupling\n",real(c1));
    fprintf(output,"       180     % 16.8e   # IMGWPHCH3   , Imaginary part of wph3hc coupling\n",imag(c1));

    get_coupling_hhh(1,1,1,c1);
    fprintf(output,"       181     % 16.8e   # REGH1H1H1   , Real part of h1h1h1 coupling\n",real(c1));
    fprintf(output,"       182     % 16.8e   # IMGH1H1H1   , Imaginary part of h1h1h1 coupling\n",imag(c1));
    get_coupling_hhh(1,1,2,c1);
    fprintf(output,"       183     % 16.8e   # REGH1H1H2   , Real part of h1h1h2 coupling\n",real(c1));
    fprintf(output,"       184     % 16.8e   # IMGH1H1H2   , Imaginary part of h1h1h2 coupling\n",imag(c1));
    get_coupling_hhh(1,2,2,c1);
    fprintf(output,"       185     % 16.8e   # REGH1H2H2   , Real part of h1h2h2 coupling\n",real(c1));
    fprintf(output,"       186     % 16.8e   # IMGH1H2H2   , Imaginary part of h1h2h2 coupling\n",imag(c1));
    get_coupling_hhh(1,3,3,c1);
    fprintf(output,"       187     % 16.8e   # REGH1H3H3   , Real part of h1h3h3 coupling\n",real(c1));
    fprintf(output,"       188     % 16.8e   # IMGH1H3H3   , Imaginary part of h1h3h3 coupling\n",imag(c1));
    get_coupling_hhh(2,2,2,c1);
    fprintf(output,"       189     % 16.8e   # REGH2H2H2   , Real part of h2h2h2 coupling\n",real(c1));
    fprintf(output,"       190     % 16.8e   # IMGH2H2H2   , Imaginary part of h2h2h2 coupling\n",imag(c1));
    get_coupling_hhh(2,3,3,c1);
    fprintf(output,"       191     % 16.8e   # REGH2H3H3   , Real part of h2h3h3 coupling\n",real(c1));
    fprintf(output,"       192     % 16.8e   # IMGH2H3H3   , Imaginary part of h2h3h3 coupling\n",imag(c1));
    get_coupling_hhh(1,4,4,c1);
    fprintf(output,"       193     % 16.8e   # REGH1HCHC   , Real part of h1hchc coupling\n",real(c1));
    fprintf(output,"       194     % 16.8e   # IMGH2HCHC   , Imaginary part of h1hchc coupling\n",imag(c1));
    get_coupling_hhh(2,4,4,c1);
    fprintf(output,"       195     % 16.8e   # REGH2HCHC   , Real part of h2hchc coupling\n",real(c1));
    fprintf(output,"       196     % 16.8e   # IMGH2HCHC   , Imaginary part of h2hchc coupling\n",imag(c1));

    get_coupling_vvhh(1,1,4,4,c1);
    fprintf(output,"       197     % 16.8e   # REGAAHCHC   ,  Real part of aahchc coupling\n",real(c1));
    fprintf(output,"       198     % 16.8e   # IMGAAHCHC   ,  Imaginary part of aahchc coupling\n",imag(c1));
    get_coupling_vvhh(1,2,4,4,c1);
    fprintf(output,"       199     % 16.8e   # REGAZHCHC   ,  Real part of azhchc coupling\n",real(c1));
    fprintf(output,"       200     % 16.8e   # IMGAZHCHC   ,  Imaginary part of azhchc coupling\n",imag(c1));
    get_coupling_vvhh(1,3,1,4,c1);
    fprintf(output,"       201     % 16.8e   # REGAWPHCH1  ,  Real part of awphch1 coupling\n",real(c1));
    fprintf(output,"       202     % 16.8e   # IMGAWPHCH1  ,  Imaginary part of awphch1 coupling\n",imag(c1));
    get_coupling_vvhh(1,3,2,4,c1);
    fprintf(output,"       203     % 16.8e   # REGAWPHCH2  ,  Real part of awphch2 coupling\n",real(c1));
    fprintf(output,"       204     % 16.8e   # IMGAWPHCH2  ,  Imaginary part of awphch2 coupling\n",imag(c1));
    get_coupling_vvhh(1,3,3,4,c1);
    fprintf(output,"       205     % 16.8e   # REGAWPHCH3  ,  Real part of awphch3 coupling\n",real(c1));
    fprintf(output,"       206     % 16.8e   # IMGAWPHCH3  ,  Imaginary part of awphch3 coupling\n",imag(c1));

    get_coupling_vvhh(2,2,1,1,c1);
    fprintf(output,"       207     % 16.8e   # REGZZH1H1   ,  Real part of azzh1h1 coupling\n",real(c1));
    fprintf(output,"       208     % 16.8e   # IMGZZH1H1   ,  Imaginary part of azzh1h1 coupling\n",imag(c1));
    get_coupling_vvhh(2,2,2,2,c1);
    fprintf(output,"       209     % 16.8e   # REGZZH2H2   ,  Real part of azzh2h2 coupling\n",real(c1));
    fprintf(output,"       210     % 16.8e   # IMGZZH2H2   ,  Imaginary part of azzh2h2 coupling\n",imag(c1));
    get_coupling_vvhh(2,2,3,3,c1);
    fprintf(output,"       211     % 16.8e   # REGZZH3H3   ,  Real part of azzh3h3 coupling\n",real(c1));
    fprintf(output,"       212     % 16.8e   # IMGZZH3H3   ,  Imaginary part of azzh3h3 coupling\n",imag(c1));
    get_coupling_vvhh(2,2,4,4,c1);
    fprintf(output,"       213     % 16.8e   # REGZZHCHC   ,  Real part of azzhchc coupling\n",real(c1));
    fprintf(output,"       214     % 16.8e   # IMGZZHCHC   ,  Imaginary part of azzhchc coupling\n",imag(c1));
    get_coupling_vvhh(2,3,1,4,c1);
    fprintf(output,"       215     % 16.8e   # REGZWPHCH1  ,  Real part of azwphch1 coupling\n",real(c1));
    fprintf(output,"       216     % 16.8e   # IMGZWPHCH1  ,  Imaginary part of azwphch1 coupling\n",imag(c1));
    get_coupling_vvhh(2,3,2,4,c1);
    fprintf(output,"       217     % 16.8e   # REGZWPHCH2  ,  Real part of azwphch2 coupling\n",real(c1));
    fprintf(output,"       218     % 16.8e   # IMGZWPHCH2  ,  Imaginary part of azwphch2 coupling\n",imag(c1));
    get_coupling_vvhh(2,3,3,4,c1);
    fprintf(output,"       219     % 16.8e   # REGZWPHCH3  ,  Real part of azwphch3 coupling\n",real(c1));
    fprintf(output,"       220     % 16.8e   # IMGZWPHCH3  ,  Imaginary part of azwphch3 coupling\n",imag(c1));

    get_coupling_vvhh(3,3,1,1,c1);
    fprintf(output,"       221     % 16.8e   # REGWWH1H1   ,  Real part of awwh1h1 coupling\n",real(c1));
    fprintf(output,"       222     % 16.8e   # IMGWWH1H1   ,  Imaginary part of awwh1h1 coupling\n",imag(c1));
    get_coupling_vvhh(3,3,2,2,c1);
    fprintf(output,"       223     % 16.8e   # REGWWH2H2   ,  Real part of awwh1h2 coupling\n",real(c1));
    fprintf(output,"       224     % 16.8e   # IMGWWH2H2   ,  Imaginary part of awwh1h2 coupling\n",imag(c1));
    get_coupling_vvhh(3,3,3,3,c1);
    fprintf(output,"       225     % 16.8e   # REGWWH3H3   ,  Real part of awwh1h3 coupling\n",real(c1));
    fprintf(output,"       226     % 16.8e   # IMGWWH3H3   ,  Imaginary part of awwh1h3 coupling\n",imag(c1));
    get_coupling_vvhh(3,3,4,4,c1);
    fprintf(output,"       227     % 16.8e   # REGWWHCHC   ,  Real part of awwhchc coupling\n",real(c1));
    fprintf(output,"       228     % 16.8e   # IMGWWHCHC   ,  Imaginary part of awwhchc coupling\n",imag(c1));
  }

  fclose(output);

  printf("LesHouches output written to file %s\n", file);
}


bool THDM::read_LesHouches(const char* file) {

  double  lambda[8];
  bool    par[32];
  double  masses[40];

  for (int i=0;i<8;i++)  lambda[i] = 0.;
  for (int i=0;i<40;i++) masses[i] = 0.;
  for (int i=0;i<32;i++) par[i] = false;

  int     model     = -1;
  double  tb        = 0.;
  double  m12_2     = 0.;
  double  sba       = 0.;
  int     type      = 0;

  int     a = 0;
  double  x = 0.;
  string  block = "EMPTY";
  string  s;

  ifstream data(file);

  if (!data) {
    cerr << "Error: Cannot read from file: " << file << "\n";
    return false;
  }

  while(getline(data,s)) {
    if (s.empty()) continue;

    // Check for comment lines and ignore them
    if(s[0] != '#') {
      istringstream ss(s);
      string cmd;
      ss >> cmd;

      // Ignore comments
      if (cmd.at(0)=='#') continue;

      // Transform to upper case letters
      transform(cmd.begin(),cmd.end(),cmd.begin(),::toupper);

      // Check for new block specifications
      if (!cmd.compare("BLOCK")) {
	block = "EMPTY";
	string blockname;
	ss >> blockname;
	transform(blockname.begin(),blockname.end(),blockname.begin(),::toupper);
	block = blockname;
      } else {
	int ncmd = atoi(cmd.c_str());

	if (block=="MODSEL") {
	  ss >> a;
	  model = a;
	}

	if (block=="SMINPUTS") {
	  ss >> x;
	  if (ncmd==1) sm.set_alpha(1./x);
	  if (ncmd==2) sm.set_GF(x);
	  if (ncmd==3) sm.set_alpha_s(x);
	  if (ncmd==4) sm.set_MZ(x);
	  if (ncmd==6) sm.set_umass_pole(3,x);
	  if (ncmd==7) sm.set_lmass_pole(3,x);
	}

	if (block=="MINPAR") {
	  ss >> x;
	  if (ncmd==3) tb = x;
	  if (ncmd==11) lambda[1] = x;
	  if (ncmd==12) lambda[2] = x;
	  if (ncmd==13) lambda[3] = x;
	  if (ncmd==14) lambda[4] = x;
	  if (ncmd==15) lambda[5] = x;
	  if (ncmd==16) lambda[6] = x;
	  if (ncmd==17) lambda[7] = x;
	  if (ncmd==18) m12_2 = x;
	  if (ncmd==20) sba = x;
	  if (ncmd==24) type = (int)x;
	  if (ncmd<32) par[ncmd] = true;
	}

	if (block=="MASS") {
	  ss >> x;
	  if (ncmd<40) masses[ncmd]=x;
	}

      }

    }
  }

  data.close();

  if (model!=10) return false;

  for (int i=1;i<=6;i++) {
    if (masses[i]>0) sm.set_qmass_pole(i,masses[i]);
  }

  if (masses[11]>0) sm.set_lmass_pole(1,masses[11]);
  if (masses[13]>0) sm.set_lmass_pole(2,masses[13]);
  if (masses[15]>0) sm.set_lmass_pole(3,masses[15]);

  if (masses[23]>0) sm.set_MZ(masses[23]);

  bool pset = false;

  if (par[3]&&(tb>0.)&&par[11]&&par[12]&&par[13]&&par[14]&&par[15]) {
    pset = set_param_gen(lambda[1],lambda[2],lambda[3],lambda[4],lambda[5],lambda[6],lambda[7],m12_2,tb);
  }

  if (!par[3]&&par[11]&&par[12]&&par[13]&&par[14]&&par[15]&&(masses[37]>0.)) {
    pset = set_param_higgs(lambda[1],lambda[2],lambda[3],lambda[4],lambda[5],lambda[6],lambda[7],masses[37]);
  }

  if (par[3]&&(tb==0.)&&par[11]&&par[12]&&par[13]&&par[14]&&par[15]&&(masses[37]>0.)) {
    pset = set_param_higgs(lambda[1],lambda[2],lambda[3],lambda[4],lambda[5],lambda[6],lambda[7],masses[37]);
  }

  if (par[3]&&(tb>0.)&&!par[11]&&!par[12]&&!par[13]&&!par[14]&&!par[15]&&par[20]&&(masses[25]>0.)&&(masses[35]>0.)&&(masses[36]>0.)&&(masses[37]>0.)) {
    pset = set_param_phys(masses[25],masses[35],masses[36],masses[37],sba,lambda[6],lambda[7],m12_2,tb);
  }

  if (par[24]) {
    set_yukawas_type(type);
  }

  return pset;

}


void THDM::print_yukawas() {
  printf("\nInvariant Yukawa matrices in convention of Haber");
  if (yukawas_type>0)
    printf(" for type %i",yukawas_type);
  printf("\n\nkappa_D\n");
  print_gsl_matrix(kappa_D,3,3);
  printf("kappa_U\n");
  print_gsl_matrix(kappa_U,3,3);
  printf("kappa_L\n");
  print_gsl_matrix(kappa_L,3,3);
  printf("rho_D\n");
  print_gsl_matrix(rho_D,3,3);
  printf("rho_U\n");
  print_gsl_matrix(rho_U,3,3);
  printf("rho_L\n");
  print_gsl_matrix(rho_L,3,3);
}


void THDM::print_param_gen() {
  double lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,tan_beta,m12_2;
  get_param_gen(lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,m12_2,tan_beta);

  printf("\n2HDM parameters in generic basis:\n");
  printf(" lambda_1: %12.5f\n",lambda1);
  printf(" lambda_2: %12.5f\n",lambda2);
  printf(" lambda_3: %12.5f\n",lambda3);
  printf(" lambda_4: %12.5f\n",lambda4);
  printf(" lambda_5: %12.5f\n",lambda5);
  printf(" lambda_6: %12.5f\n",lambda6);
  printf(" lambda_7: %12.5f\n",lambda7);
  printf("    m12^2: %12.5f\n",m12_2);
  printf("tan(beta): %12.5f\n",tan_beta);
}

void THDM::print_param_hybrid() {
  double mh,mH,cba,tb,Z4,Z5,Z7;
  get_param_hybrid(mh,mH,cba,Z4,Z5,Z7,tb);
  double lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,tan_beta,m12_2;
  get_param_gen(lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,m12_2,tan_beta);

  printf("\n2HDM parameters in Hybrid basis:\n");
  if ((abs(lambda6)>EPS)||(abs(lambda7)>EPS)) {
   printf("WARNING: Model has hard Z_2 violation\n");
   printf("Output in H2 basis is inconsistent\n");
  }
  printf("      m_h: %12.5f\n",mh);
  printf("      m_H: %12.5f\n",mH);
  printf(" cos(b-a): %12.5f\n",cba);
  printf("       Z4: %12.5f\n",Z4);
  printf("       Z5: %12.5f\n",Z5);
  printf("       Z7: %12.5f\n",Z7);
  printf("tan(beta): %12.5f\n",tb);
}

void THDM::print_hdecay() {
  double mh,mH,mA,mHp,sba,lambda6,lambda7,tan_beta,m12_2,alpha;
  get_param_phys(mh,mH,mA,mHp,sba,lambda6,lambda7,m12_2,tan_beta);

  alpha = -asin(sba)+atan(tan_beta);

  printf("\n2HDM output for HDECAY\n");
  printf(" *************************** 2 Higgs Doublet Model **************************\n");
  printf(" TYPE; 1 (I), 2 (II)\n");
  printf(" TYPE = %d\n", yukawas_type);
  printf(" TGBET2HDM = %.8fD0\n", tan_beta);
  printf(" ALPHA_H = %.8fD0\n", alpha);
  printf(" MHL = %.8fD0\n", mh);
  printf(" MHH = %.8fD0\n", mH);
  printf(" MHA = %.8fD0\n", mA);
  printf(" MH+- = %.8fD0\n", mHp);
  printf(" M12^2 = %.8fD0\n", m12_2);
  printf(" ****************************************************************************\n");

}


void THDM::print_param_HHG() {
  double Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,tan_beta;
  get_param_HHG(Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,tan_beta);

  printf("\n2HDM parameters in Higgs Hunter's basis:\n");
  printf(" Lambda_1: %12.5f\n",Lambda1);
  printf(" Lambda_2: %12.5f\n",Lambda2);
  printf(" Lambda_3: %12.5f\n",Lambda3);
  printf(" Lambda_4: %12.5f\n",Lambda4);
  printf(" Lambda_5: %12.5f\n",Lambda5);
  printf(" Lambda_6: %12.5f\n",Lambda6);
  printf("tan(beta): %12.5f\n",tan_beta);
}


void THDM::print_param_phys() {
  double mh,mH,mA,mHp,sba,lambda6,lambda7,tan_beta,m12_2;
  get_param_phys(mh,mH,mA,mHp,sba,lambda6,lambda7,m12_2,tan_beta);

  printf("\n2HDM parameters in physical mass basis:\n");
  printf("      m_h: %12.5f\n",mh);
  printf("      m_H: %12.5f\n",mH);
  printf("      m_A: %12.5f\n",mA);
  printf("     m_H+: %12.5f\n",mHp);
  printf(" sin(b-a): %12.5f\n",sba);
  printf(" lambda_6: %12.5f\n",lambda6);
  printf(" lambda_7: %12.5f\n",lambda7);
  printf("    m12^2: %12.5f\n",m12_2);
  printf("tan(beta): %12.5f\n",tan_beta);
}

void THDM::print_param_higgs() {
  double Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,mHp;
  get_param_higgs(Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,mHp);

  printf("\n2HDM parameters in Higgs basis:\n");
  printf(" Lambda_1: %12.5f\n",Lambda1);
  printf(" Lambda_2: %12.5f\n",Lambda2);
  printf(" Lambda_3: %12.5f\n",Lambda3);
  printf(" Lambda_4: %12.5f\n",Lambda4);
  printf(" Lambda_5: %12.5f\n",Lambda5);
  printf(" Lambda_6: %12.5f\n",Lambda6);
  printf(" Lambda_7: %12.5f\n",Lambda7);
  printf("     m_Hp: %12.5f\n",mHp);
}


double THDM::get_m12_2() {
  //hep-ph/0207010
  double sb=sin(beta);
  double sb2=sb*sb;
  double cb=cos(beta);
  double cb2=cb*cb;
  double tb=tan(beta);
  double ctb=1./tb;

  double m12_2 = 0.;

  if (tb>0) {
    m12_2=(m22_2+0.5*v2*(lambda[2]*sb2+(lambda[3]+lambda[4]+lambda[5])*cb2+lambda[6]*cb2*ctb+3.*lambda[7]*sb*cb))*tb;
  } else {
    m12_2=0.5*v2*lambda[6];
  }

  return m12_2;
}


SM THDM::get_SM() {
  return sm;
}


void THDM::set_SM(SM sm_in) {
  sm=sm_in;
  v2 = sm.get_v2();


}


double THDM::get_hmass(int h) {

  if ((h<1)||(h>4)) return 0.;

  double mh[5],a,l6,l7,tb,m12_2;
  get_param_phys(mh[1],mh[2],mh[3],mh[4],a,l6,l7,m12_2,tb);
  return mh[h];
}


double THDM::get_sba() {

  double mh,mH,mA,mHp,sba,l6,l7,tb,m12_2;
  get_param_phys(mh,mH,mA,mHp,sba,l6,l7,m12_2,tb);

  return sba;
}


double THDM::get_cba() {
  double c2 = sqrt(1.-pow(get_sba(),2));
  return c2;
}

double THDM::get_alpha() {

  double mh,mH,mA,mHp,sba,l6,l7,tb,m12_2;
  get_param_phys(mh,mH,mA,mHp,sba,l6,l7,m12_2,tb);

  double beta = atan(tb);
  double cba = get_cba();

  double ba = atan2(sba, cba);
  double a = -ba+beta;

  if (a>M_PI/2) a = a - M_PI;

  printf("Getting alpha: %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n", tb, sba, cba, ba, beta,  a);


  return a;

}


void THDM::print_info() {
  first_run=false;
  printf("****************************************************\n");
  printf("*                                                  *\n");
  printf("*    2HDMC - Two-Higgs-Doublet Model Calculator    *\n");
  printf("*             http://2hdmc.hepforge.org            *\n");
  printf("*                  Version %-10s              *\n",version);
  printf("*             Compiled on %-24s *\n",__DATE__);
  printf("*                                                  *\n");
  printf("****************************************************\n");
}
