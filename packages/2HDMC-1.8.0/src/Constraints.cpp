#include "Constraints.h"
#include "DecayTable.h"
#include "THDM.h"
#include "Util.h"
#include <fstream>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <cfloat>
#include <cstring>

using namespace std;


Constraints::Constraints() : table(model) {
//  init_externals();
  delta = exp(-log(10.)*DBL_DIG/5.);
}


Constraints::Constraints(THDM mod) : table(mod) {
//  init_externals();
  set_THDM(mod);
  delta = exp(-log(10.)*DBL_DIG/5.);
}

void Constraints::set_THDM(THDM mod) {
  model = mod;
  sm    = mod.get_SM();

//  table.set_model(mod);

}

void Constraints::init_externals() {

// This is deprecated. Use methods defined in HBHS.h instead.
	return;

#if defined NMSSMTools
  initialize_();
#endif
// #if defined HiggsBounds
//   int nH0=3;
//   int nHp=1;
//   int hbflag=3;
//   printf("\nInitializing HiggsBounds... ");
//
//
//    // Third argument is HB analysis setting: 1='onlyL', 2='onlyH' 3='LandH'
//     initialize_higgsbounds_int_(&nH0, &nHp, &hbflag);
// //  initialize_higgsbounds_(&nH0,&nHp,whichexpt);
//
//   printf("Please cite the relevant publications when using Higgs mass limits.\n");
//
// #endif
//  init_Hp();
}

void Constraints::init_Hp() {
  ifstream f;

  const char* Hpf_tau = "Hp_limit_tau";
  const char* Hpf_cs = "Hp_limit_cs";

  // tautau channel
  f.open(Hpf_tau);
  if (f.good()) {

    f >> nHp1;

    mHp1 = new double[nHp1];
    valHp1 = new double[nHp1];

    for (int i=0;(i<nHp1);i++) {
      f >> mHp1[i] >> valHp1[i];
    }

    f.close();
  } else {
    cerr << "Error: Cannot open \"" << Hpf_tau << "\" for reading LEP data\n";
  }

  // cs channel
  f.open(Hpf_cs);
  if (f.good()) {

    f >> nHp2;

    mHp2 = new double[nHp2];
    valHp2 = new double[nHp2];

    for (int i=0;(i<nHp2);i++) {
      f >> mHp2[i] >> valHp2[i];
    }

    f.close();
  } else {
    cerr << "Error: Cannot open \"" << Hpf_cs << "\" for reading LEP data\n";
  }


}

void Constraints::oblique_param(double mh, double &S, double &T, double &U, double &V, double &W, double &X) {

  double g = sm.get_g();
  double alpha = sm.get_alpha();
  double stw = sm.get_sintw();
  double ctw = sm.get_costw();
  double mw2 = pow(sm.get_MW(),2);
  double mz2 = pow(sm.get_MZ(),2);
  double mh2 = mh*mh;

  complex <double> I(0.0,1.0);

  double m2[3];
  double mu2[5];

  complex <double> q11 = model.get_qki(1,1);
  complex <double> q12 = model.get_qki(1,2);
  complex <double> q21 = model.get_qki(2,1);
  complex <double> q22 = model.get_qki(2,2);

  m2[1] = 0.;
  m2[2] = pow(model.get_hmass(4),2);

  mu2[1] = 0.;
  mu2[2] = pow(model.get_hmass(3),2);
  mu2[3] = pow(model.get_hmass(1),2);
  mu2[4] = pow(model.get_hmass(2),2);

  double UdU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  UdU[1][1]=1.0;
  UdU[2][2]=1.0;

  complex <double> UdV[3][5] = {{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.}};
  UdV[1][1]=I;
  UdV[2][2]=-I;
  UdV[1][3]=q11;
  UdV[1][4]=q21;
  UdV[2][3]=q12;
  UdV[2][4]=q22;

  complex <double> VdV[5][5] = {{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.}};
  VdV[1][1]=1.;
  VdV[1][3]=-I*q11;
  VdV[1][4]=-I*q21;
  VdV[2][2]=1.;
  VdV[2][3]=I*q12;
  VdV[2][4]=I*q22;
  VdV[3][1]=I*conj(q11);
  VdV[3][2]=-I*conj(q12);
  VdV[3][3]=conj(q11)*q11+conj(q12)*q12;
  VdV[3][4]=conj(q11)*q21+conj(q12)*q22;
  VdV[4][1]=I*conj(q21);
  VdV[4][2]=-I*conj(q22);
  VdV[4][3]=conj(q21)*q11+conj(q22)*q12;
  VdV[4][4]=conj(q21)*q21+conj(q22)*q22;

  double ImVdV[5][5] = {{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.}};
  for (int i=0; i<5;i++) {
    for (int j=0;j<5;j++) {
      ImVdV[i][j]=imag(VdV[i][j]);
    }
  }

  double S1=0.,S2=0.,S3=0.,S4=0.,S5=0.,S6=0.,S7=0.;
  double T1=0.,T2=0.,T3=0.,T4=0.;
  double U1=0.,U2=0.,U3=0.,U4=0.,U5=0.,U6=0.;
  double V1=0.,V2=0.,V3=0.,V4=0.;
  double CS = pow(g,2)/(384.*pow(M_PI,2)*pow(ctw,2));
  double CT = pow(g,2)/(64.*pow(M_PI,2)*mw2);
  double CU = pow(g,2)/(384.*pow(M_PI,2));
  double CV = pow(g,2)/(384.*pow(M_PI,2)*pow(ctw,2));
  double CW = pow(g,2)/(384.*pow(M_PI,2));
  double CX = -pow(g,2)*stw/(192.*pow(M_PI,2)*ctw);
  S = 0.;
  T = 0.;
  U = 0.;
  V = 0.;
  W = 0.;
  X = 0.;

// Oblique parameter S
  S1 += CS*pow((2.*pow(stw,2)-UdU[2][2]),2)*G_fcn(m2[2],m2[2],mz2);

  for (int b=2;b<=3;b++) {
    for (int bp=b+1;bp<=4;bp++) {
      S2 += CS*pow(ImVdV[b][bp],2)*G_fcn(mu2[b],mu2[bp],mz2);
    }
  }

  S3 += -2.*CS*UdU[2][2]*log(m2[2]/mw2);

  for (int b=2;b<=4;b++) {
    S4 += CS*real(VdV[b][b])*log(mu2[b]/mw2);
    S5 += CS*pow(ImVdV[1][b],2)*Ghat_fcn(mu2[b],mz2);
  }

     S6 += -CS*Ghat_fcn(mh2,mz2);
     S7 += -CS*log(mh2/mw2);


  S=S1+S2+S3+S4+S5+S6+S7;

  // Oblique parameter T
  for (int b=2;b<=4;b++) {
    T1 += CT*pow((abs(UdV[2][b])),2)*Fdrho(m2[2],mu2[b]);
  }

  for (int b=2;b<=3;b++) {
    for (int bp=b+1;bp<=4;bp++) {
      T2 += -CT*pow(ImVdV[b][bp],2)*Fdrho(mu2[b],mu2[bp]);
    }
  }

  for (int b=2;b<=4;b++) {
    T3 += CT*3.*pow(ImVdV[1][b],2)*(Fdrho(mz2,mu2[b])-Fdrho(mw2,mu2[b]));
  }

  T4 += -CT*3.*(Fdrho(mz2,mh2)-Fdrho(mw2,mh2));
  T = T1+T2+T3+T4;

  // Oblique parameter U
  for (int b=2;b<=4;b++) {
    U1 += CU*pow(abs(UdV[2][b]),2)*G_fcn(m2[2],mu2[b],mw2);
  }

  U2 = -CU*pow((2.*pow(stw,2)-UdU[2][2]),2)*G_fcn(m2[2],m2[2],mz2);

  for (int b=2;b<=3;b++) {
    for (int bp=b+1;bp<=4;bp++) {
      U3 += -CU*pow(ImVdV[b][bp],2)*G_fcn(mu2[b],mu2[bp],mz2);
    }
  }

  for (int b=2;b<=4;b++) {
    U4 += CU*pow(ImVdV[1][b],2)*(Ghat_fcn(mu2[b],mw2)-Ghat_fcn(mu2[b],mz2));
  }

  U5 = -CU*Ghat_fcn(mh2,mw2);
  U6 = CU*Ghat_fcn(mh2,mz2);


  U = U1+U2+U3+U4+U5+U6;

  // Oblique parameter V
  V1+=CV*pow((2.*pow(stw,2)-UdU[2][2]),2)*H_fcn(m2[2],m2[2],mz2);

  for (int b=2;b<=3;b++) {
    for (int bp=b+1;bp<=4;bp++) {
      V2+=CV*pow(ImVdV[b][bp],2)*H_fcn(mu2[b],mu2[bp],mz2);
    }
  }

  for (int b=2;b<=4;b++) {
    V3+=CV*pow(ImVdV[1][b],2)*Hhat_fcn(mu2[b],mz2);
  }

  V4+=-CV*Hhat_fcn(mh2,mz2);

  V=V1+V2+V3+V4;

  // Oblique parameter W
  for (int b=2;b<=4;b++) {
    W+=CW*pow(abs(UdV[2][b]),2)*H_fcn(m2[2],mu2[b],mw2);
  }

  for (int b=2;b<=4;b++) {
    W+=CW*pow(ImVdV[1][b],2)*Hhat_fcn(mu2[b],mw2);
  }

  W+=-CW*Hhat_fcn(mh2,mw2);

  // Oblique parameter X
  X = CX*(2.*pow(stw,2)-UdU[2][2])*G_fcn(m2[2],m2[2],mz2);

  // Rescaling
  S = S*4.*pow(stw,2)*pow(ctw,2)/alpha;
  T = T/alpha;
  U = U*4.*pow(stw,2)/alpha;
  V = V/alpha;
  W = W/alpha;
  X = X*stw*ctw/alpha;

}


double Constraints::G_fcn(double x, double y, double q) {

  double t = x+y-q;
  double r = q*q-2.*q*(x+y)+pow(x-y,2);
  double G=0.;

  if ((x>0)&&(y>0)&&(abs(x-y)>delta*x)) {
    G = -16./3.+5.*(x+y)/q-2.*pow(x-y,2)/pow(q,2)+3./q*((x*x+y*y)/(x-y)-(x*x-y*y)/q+pow(x-y,3)/(3.*q*q))*log(x/y)+r/pow(q,3)*f_fcn(t,r);
  } else {
    G = -16./3.+5.*(x+y)/q-2.*pow(x-y,2)/pow(q,2)+3./q*((x*x+y*y)-(x*x-y*y)*(x-y)/q+pow(x-y,4)/(3.*q*q))/y*(1.-(x-y)/y/2.+pow((x-y)/y,2)/3.-pow((x-y)/y,3)/4.+pow((x-y)/y,4)/5.)+r/pow(q,3)*f_fcn(t,r);
  }
  return G;
}


double Constraints::Gtilde_fcn(double x, double y, double q) {

  double t = x+y-q;
  double r = q*q-2.*q*(x+y)+pow(x-y,2);

  double G = 0.;

  if ((x>0)&&(y>0)&&(abs(x-y)>delta*x)) {
    G = -2.+((x-y)/q-(x+y)/(x-y))*log(x/y)+f_fcn(t,r)/q;
  } else {
    G = -2.+(pow(x-y,2)/q-(x+y))/y*(1.-(x-y)/y/2.+pow((x-y)/y,2)/3.-pow((x-y)/y,3)/4.+pow((x-y)/y,4)/5.)+f_fcn(t,r)/q;
  }
  return G;
}


double Constraints::Ghat_fcn(double x,double q) {

  double G = G_fcn(x,q,q)+12.*Gtilde_fcn(x,q,q);
  return G;
}


double Constraints::f_fcn(double t, double r) {

  if (r==0) return 0.;

  if (r>0) {
    return sqrt(r)*log(abs(t-sqrt(r))/abs(t+sqrt(r)));
  } else if (r<0) {
    return 2.*sqrt(-r)*atan(sqrt(-r)/t);
  }

  return 0.;
}


double Constraints::H_fcn(double x, double y, double q) {

  double t = x+y-q;
  double r = q*q-2.*q*(x+y)+pow(x-y,2);
  double H=0.;

  if ((x>0)&&(y>0)&&(abs(x-y)>delta*x)) {
    H = 2.-9.*(x+y)/q+6.*pow(x-y,2)/pow(q,2)+3./q*(-(x*x+y*y)/(x-y)+2.*(x*x-y*y)/q-pow(x-y,3)/pow(q,2))*log(x/y)+(x+y-pow(x-y,2)/q)*3.*f_fcn(t,r)/pow(q,2);
  } else {
    H= 2.-9.*(x+y)/q+6.*pow(x-y,2)/pow(q,2)+3./q*(-(x*x+y*y)+2.*(x*x-y*y)*(x-y)/q-pow(x-y,4)/pow(q,2))/y*(1.-(x-y)/y/2.+pow((x-y)/y,2)/3.-pow((x-y)/y,3)/4.+pow((x-y)/y,4)/5.)+(x+y-pow(x-y,2)/q)*3.*f_fcn(t,r)/pow(q,2);
  }
  return H;
}


double Constraints::Htilde_fcn(double x, double y, double q) {

  double t = x+y-q;
  double r = q*q-2.*q*(x+y)+pow(x-y,2);
  double H=0.;

  if ((x>0)&&(y>0)&&(abs(x-y)>delta*x)) {
    H=4.+((x+y)/(x-y)-2.*(x-y)/q)*log(x/y)+(-pow(q,2)+3.*q*(x+y)-2.*pow(x-y,2))/(r*q)*f_fcn(t,r);
  } else {
    H=4.+((x+y)-2.*pow(x-y,2)/q)/y*(1.-(x-y)/y/2.+pow((x-y)/y,2)/3.-pow((x-y)/y,3)/4.+pow((x-y)/y,4)/5.)+(-pow(q,2)+3.*q*(x+y)-2.*pow(x-y,2))/(r*q)*f_fcn(t,r);
  }
  return H;
}


double Constraints::Hhat_fcn(double x,double q) {

  double H = H_fcn(x,q,q)+12.*Htilde_fcn(x,q,q);
  return H;
}


double Constraints::delta_amu() {

  double alpha = sm.get_alpha();

  double Qd = -1./3.;
  double Qu =  2./3.;
  double Ql = -1.;
  int 	 Nc =  3;

  double m_mu = sm.get_lmass_pole(2);
  double mt = sm.get_qmass_pole(6);
  double mb = sm.get_qmass_pole(5);

  complex <double> cs, cp;

  double damutot 		= 0.;
  double damu_1loop = 0.;
  double damu_2loop = 0.;

  for (int h=1;h<=4;h++) {
    double mh = model.get_hmass(h);

    if (h<4) {
      model.get_coupling_hll(h,2,2,cs,cp);
    } else {
      model.get_coupling_hln(h,2,2,cs,cp);
    }

    double CS_mu=imag(cs);
    double CP_mu=real(cp);

    // CP conservation => Either S or P coupling for each Higgs boson
    double C = max(abs(CS_mu),abs(CP_mu));

    damu_1loop += pow(C,2)/(8.*pow(M_PI,2))*dmu_L(pow(mh/m_mu,2),h);
  }


  for (int h=1;h<=3;h++) {
    double mh = model.get_hmass(h);

    model.get_coupling_hll(h,2,2,cs,cp);
    double CS_mu=imag(cs);
    double CP_mu=real(cp);

    // Down-type quarks
    for (int i=1;i<=3;i++) {
      complex <double> gS(0.,0.), gP(0.,0.);
      model.get_coupling_hdd(h,i,i,gS,gP);
      double m = sm.get_dmass_MSbar(i);

      if (m > 1E-3) {
	double mrun = sm.run_qmass_MSbar(m,m,mh,mt,mb);
	double x = pow(mrun/mh,2);

	double damu = Nc*alpha/(4.*pow(M_PI,3))*pow(Qd,2)*m_mu*mrun/pow(mh,2)*(CP_mu*real(gP)*dmu_g(x)-CS_mu*imag(gS)*dmu_f(x));

	damu_2loop += damu;
      }
    }

    // Up-type quarks
    for (int i=1;i<=3;i++) {
      complex <double> gS(0.,0.), gP(0.,0.);
      model.get_coupling_huu(h,i,i,gS,gP);
      double m = sm.get_umass_MSbar(i);

      if (m > 1E-3) {
	double mrun = sm.run_qmass_MSbar(m,m,mh,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
	double x = pow(mrun/mh,2);

	double damu = Nc*alpha/(4.*pow(M_PI,3))*pow(Qu,2)*m_mu*mrun/pow(mh,2)*(CP_mu*real(gP)*dmu_g(x)-CS_mu*imag(gS)*dmu_f(x));

	damu_2loop += damu;
      }
    }

    // Leptons
    for (int i=1;i<=3;i++) {
      complex <double> gS(0.,0.), gP(0.,0.);
      model.get_coupling_hll(h,i,i,gS,gP);
      double m = sm.get_lmass_pole(i);

      if (m > 1E-3) {
	double x = pow(m/mh,2);

	double damu = alpha/(4.*pow(M_PI,3))*pow(Ql,2)*m_mu*m/pow(mh,2)*(CP_mu*real(gP)*dmu_g(x)-CS_mu*imag(gS)*dmu_f(x));

	damu_2loop += damu;
      }
    }

  }

  damutot = damu_1loop+damu_2loop;

  return damutot;
}


double Constraints::dmu_L(double z, int h) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;

  if (h==1) F.function = Lint_s;
  if (h==2) F.function = Lint_s;
  if (h==3) F.function = Lint_p;
  if (h==4) F.function = Lint_c;
  F.params = &z;

  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F, 0, 1, 0, 1e-10, 1000,
			w, &result, &error);
  if (status)
    if (status!=GSL_EROUND) {
      printf ("gsl error: %s\n", gsl_strerror (status));
      exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);

  return result;

}


double Constraints::dmu_f(double z) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  F.function = fint;
  F.params = &z;

  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F, 0, 1, 0, 1e-10, 1000,
			w, &result, &error);
  if (status)
    if (status!=GSL_EROUND) {
      printf ("gsl error: %s\n", gsl_strerror (status));
      exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);

  double f = 0.5*result;
  return f;

}


double Constraints::dmu_g(double z) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  F.function = gint;
  F.params = &z;

  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F, 0, 1, 0, 1e-10, 1000,
			w, &result, &error);
  if (status)
    if (status!=GSL_EROUND) {
      printf ("gsl error: %s\n", gsl_strerror (status));
      exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);

  double f = 0.5*result;
  return f;
}


bool Constraints::check_unitarity(double unitarity_limit) {
  return model.check_unitarity(unitarity_limit);
}

bool Constraints::check_perturbativity(double perturbativity_limit) {
  return model.check_perturbativity(perturbativity_limit);
}


bool Constraints::check_positivity() {
  return model.check_stability();
}


bool Constraints::check_stability() {
  return model.check_stability();
}


bool Constraints::check_lep() {
  return check_masses();
}

bool Constraints::check_masses() {

	cout << "Constraints::check_masses() has been obsoleted" << endl;
	cout << "Please use new HBHS interface" << endl;
	return true;
}

bool Constraints::check_charged(bool &HpHp, bool &HpHptau, bool &HpHpcs) {

  // Limit from Z width on charged Higgs mass, hep-ex/0404012
  const double Z_LIMIT_MCH = 39.6;

  double mcH = model.get_hmass(4);

  HpHp = (mcH > Z_LIMIT_MCH);

  double GF = sm.get_GF();
  double mw = sm.get_MW();
  double mz = sm.get_MZ();
  double sw = sm.get_sintw();
  double cw = sm.get_costw();

  // sqrt(s) for LEP-II run
  double E_LEP = 209.;
  double s = pow(E_LEP,2);

  double ae = -1./(4.*cw*sw);
  double ve = (-1.+4.*pow(sw,2))/(4.*cw*sw);
  double vH = (-1.+2.*pow(sw,2))/(2.*cw*sw);

  double beta = sqrt(1.-4.*pow(mcH,2)/s);

  // Convert (GeV)^-2 to picobarns
  double C = 0.3894*1E9;

  // e+e- --> H+H- cross-section in picobarn
  double s_eeHH = 2.*pow(GF,2)*pow(mw,4)*pow(sw,4)/(3.*M_PI*s)*(1.+2*ve*vH/(1.-pow(mz,2)/s)+(pow(ae,2)+pow(ve,2))*pow(vH,2)/pow(1-pow(mz,2)/s,2))*pow(beta,3)*C;

   double wHp_tot = table.get_gammatot_h(4);
   double wHp_tau = table.get_gamma_hln(4,3,3);
   double wHp_cs  = table.get_gamma_hdu(4,2,2);

  double lim_brHp_tau = 1E+6;
  double lim_brHp_cs  = 1E+6;

  HpHptau = true;
  HpHpcs = true;

  if (wHp_tot>0) {
    double brHp_tau = wHp_tau/wHp_tot;
    double brHp_cs  = wHp_cs/wHp_tot;

    int n = -1;
    for (int i=0;(i<nHp1);i++) {
      if (mHp1[i]>mcH) {
        n=i;
        break;
      }
    }

    if (n>0) {
      double valmin = valHp1[n-1];
      double valmax = valHp1[n];
      double mmin = mHp1[n-1];
      double mmax = mHp1[n];

      double myval = valmin + (mcH-mmin)/(mmax-mmin)*(valmax-valmin);
      lim_brHp_tau = myval;
    }

    n = -1;
    for (int i=0;(i<nHp2);i++) {
      if (mHp2[i]>mcH) {
      n = i;
      break;
      }
    }

    if (n>0) {
      double valmin = valHp2[n-1];
      double valmax = valHp2[n];
      double mmin = mHp2[n-1];
      double mmax = mHp2[n];

      double myval = valmin + (mcH-mmin)/(mmax-mmin)*(valmax-valmin);
      lim_brHp_cs = myval;
    }

    if (s_eeHH*brHp_tau>lim_brHp_tau) HpHptau = false;
    if (s_eeHH*brHp_cs>lim_brHp_cs) HpHpcs = false;
  }

  return HpHp&&HpHptau&&HpHpcs;
}

#if defined NMSSMTools
bool Constraints::check_NMSSMTools(bool &hZ, bool &hZ2b, bool &hZ2tau, bool &hZinv, bool &hZ2j,
				   bool &hZ2gamma, bool &hZ4b, bool &hZ4tau, bool &hZ2b2tau,
				   bool &hA, bool &hA4b, bool &hA4tau, bool &hA2b2tau,
				   bool &hA6b, bool &hA6tau, bool &ZhZjj) {

  // The constraints on the neutral Higgses are collected from NMSSMTools,
  // if available
  int i;
  gauge_.alsmz=sm.get_alpha_s();
  gauge_.alemz=sm.get_alpha();
  gauge_.gf=sm.get_GF();
  gauge_.g1=sm.get_gprime()*sm.get_gprime();
  gauge_.g2=sm.get_g()*sm.get_g();
  gauge_.s2tw=sm.get_sintw()*sm.get_sintw();

  smspec_.ms=sm.get_dmass_pole(2);
  smspec_.mc=sm.get_umass_pole(2);
  smspec_.mb=sm.get_dmass_MSbar(3);
  smspec_.mbp=sm.get_dmass_pole(3);
  smspec_.mt=sm.get_umass_pole(3);
  smspec_.mtau=sm.get_lmass_pole(3);
  smspec_.mmuon=sm.get_lmass_pole(2);
  smspec_.mz=sm.get_MZ();
  smspec_.mw=sm.get_MW();

  susyspec_.mgl=0.;
  int j,k;
  for (i=0;i<2;i++)
    susyspec_.mch[i]=0.;
  for (i=0;i<2;i++)
    for (j=0;j<2;j++) {
      susyspec_.u[i][j]=0.;
      susyspec_.v[i][j]=0.;
    }
  for (i=0;i<5;i++)
    susyspec_.mneu[i]=0.;
  for (i=0;i<5;i++)
    for (j=0;j<5;j++)
      susyspec_.neu[i][j]=0.;

  sfspec_.mur=0.;
  sfspec_.mul=0.;
  sfspec_.mdr=0.;
  sfspec_.mdl=0.;
  sfspec_.mlr=0.;
  sfspec_.mll=0.;
  sfspec_.mnl=0.;
  sfspec_.mst1=0.;
  sfspec_.mst2=0.;
  sfspec_.msb1=0.;
  sfspec_.msb2=0.;
  sfspec_.msl1=0.;
  sfspec_.msl2=0.;
  sfspec_.msnt=0.;
  sfspec_.cst=0.;
  sfspec_.csb=0.;
  sfspec_.csl=0.;
  sfspec_.msmu1=0.;
  sfspec_.msmu2=0.;
  sfspec_.msmunt=0.;
  sfspec_.csmu=0.;

  double mh,mH,mA,mHp,sba,m12_2,lam6,lam7,tb;
  model.get_param_phys(mh,mH,mA,mHp,sba,lam6,lam7,m12_2,tb);
  double alpha=-asin(sba)+atan(tb);
  higgspec_.smass[0]=mh;
  higgspec_.smass[1]=mH;
  higgspec_.smass[2]=1e10;
  higgspec_.pmass[0]=mA;
  higgspec_.pmass[1]=1e10;
  higgspec_.scomp[0][0]=cos(alpha);
  higgspec_.scomp[1][0]=-sin(alpha);
  higgspec_.scomp[2][0]=0.;
  higgspec_.scomp[0][1]=sin(alpha);
  higgspec_.scomp[1][1]=cos(alpha);
  higgspec_.scomp[2][1]=0.;
  higgspec_.scomp[0][2]=0.;
  higgspec_.scomp[1][2]=0.;
  higgspec_.scomp[2][2]=1.;
  higgspec_.pcomp[0][0]=1.;
  higgspec_.pcomp[1][0]=0.;
  higgspec_.pcomp[0][1]=0.;
  higgspec_.pcomp[1][1]=1.;
  higgspec_.cmass=mHp;

  // h1,h2,h3,a1,a2 are Higgs bosons in NMSSMTools so for us h3 and a2 does not exist
  // h1,h2,h3,a1,a2 total width
  brn_.width[0]=table.get_gammatot_h(1);
  brn_.width[1]=table.get_gammatot_h(2);
  brn_.width[3]=table.get_gammatot_h(3);
  brn_.width[2]=0.;
  brn_.width[4]=0.;
  // h1,h2,h3,a1,a2 -> gluon gluon
  for (i=0;i<5;i++)
    brn_.brjj[i]=0.;
  // h1,h2,h3,a1,a2 -> mu mu
  brn_.brmm[0]=table.get_gamma_hll(1,2,2)/brn_.width[0];
  brn_.brmm[1]=table.get_gamma_hll(2,2,2)/brn_.width[1];
  brn_.brmm[3]=table.get_gamma_hll(3,2,2)/brn_.width[3];
  brn_.brmm[2]=0.;
  brn_.brmm[4]=0.;
  // h1,h2,h3,a1,a2 -> tau tau
  brn_.brll[0]=table.get_gamma_hll(1,3,3)/brn_.width[0];
  brn_.brll[1]=table.get_gamma_hll(2,3,3)/brn_.width[1];
  brn_.brll[3]=table.get_gamma_hll(3,3,3)/brn_.width[3];
  brn_.brll[2]=0.;
  brn_.brll[4]=0.;
  // h1,h2,h3,a1,a2 -> s s
  brn_.brss[0]=table.get_gamma_hdd(1,2,2)/brn_.width[0];
  brn_.brss[1]=table.get_gamma_hdd(2,2,2)/brn_.width[1];
  brn_.brss[3]=table.get_gamma_hdd(3,2,2)/brn_.width[3];
  brn_.brss[2]=0.;
  brn_.brss[4]=0.;
  // h1,h2,h3,a1,a2 -> c c
  brn_.brcc[0]=table.get_gamma_huu(1,2,2)/brn_.width[0];
  brn_.brcc[1]=table.get_gamma_huu(2,2,2)/brn_.width[1];
  brn_.brcc[3]=table.get_gamma_huu(3,2,2)/brn_.width[3];
  brn_.brcc[2]=0.;
  brn_.brcc[4]=0.;
  // h1,h2,h3,a1,a2 -> b b
  brn_.brbb[0]=table.get_gamma_hdd(1,3,3)/brn_.width[0];
  brn_.brbb[1]=table.get_gamma_hdd(2,3,3)/brn_.width[1];
  brn_.brbb[3]=table.get_gamma_hdd(3,3,3)/brn_.width[3];
  brn_.brbb[2]=0.;
  brn_.brbb[4]=0.;
  // h1,h2,h3,a1,a2 -> t t
  brn_.brtt[0]=table.get_gamma_huu(1,3,3)/brn_.width[0];
  brn_.brtt[1]=table.get_gamma_huu(2,3,3)/brn_.width[1];
  brn_.brtt[3]=table.get_gamma_huu(3,3,3)/brn_.width[3];
  brn_.brtt[2]=0.;
  brn_.brtt[4]=0.;
  // h1,h2,h3,a1,a2 -> w w
  brn_.brww[0]=table.get_gamma_hvv(1,3)/brn_.width[0];
  brn_.brww[1]=table.get_gamma_hvv(2,3)/brn_.width[1];
  //brn_.brww[3]=table.get_gamma_hvv(3,3)/brn_.width[3];
  brn_.brww[2]=0.;
  //brn_.brww[4]=0.;
  // h1,h2,h3,a1,a2 -> z z
  brn_.brzz[0]=table.get_gamma_hvv(1,2)/brn_.width[0];
  brn_.brzz[1]=table.get_gamma_hvv(2,2)/brn_.width[1];
  //brn_.brzz[3]=table.get_gamma_hvv(3,2)/brn_.width[3];
  brn_.brzz[2]=0.;
  //brn_.brzz[4]=0.;
  // h1,h2,h3,a1,a2 -> gamma gamma
  brn_.brgg[0]=table.get_gamma_hvv(1,1)/brn_.width[0];
  brn_.brgg[1]=table.get_gamma_hvv(2,1)/brn_.width[1];
  brn_.brgg[3]=table.get_gamma_hvv(3,1)/brn_.width[3];
  brn_.brgg[2]=0.;
  brn_.brgg[4]=0.;
  // h1,h2,h3,a1,a2 -> Z gamma
  brn_.brzg[0]=0.;
  brn_.brzg[1]=0.;
  brn_.brzg[3]=0.;
  brn_.brzg[2]=0.;
  brn_.brzg[4]=0.;
  // h2 -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
  brn_.brhhh[0]=table.get_gamma_hhh(2,1,1)/brn_.width[1];
  brn_.brhhh[1]=0.;
  brn_.brhhh[2]=0.;
  brn_.brhhh[3]=0.;
  // hi -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
  for (i=0;i<2;i++) {
    brn_.brhaa[0][i]=table.get_gamma_hhh(i+1,3,3)/brn_.width[i];
    brn_.brhaa[1][i]=0.;
    brn_.brhaa[2][i]=0.;
  }
  brn_.brhaa[0][2]=0.;
  brn_.brhaa[1][2]=0.;
  brn_.brhaa[2][2]=0.;
  // hi -> h+h- (i=1..3)
  brn_.brhchc[0]=table.get_gamma_hhh(1,4,4)/brn_.width[0];
  brn_.brhchc[1]=table.get_gamma_hhh(2,4,4)/brn_.width[1];
  brn_.brhchc[2]=0.;
  // hi -> Zaj  (i=1..3, j=1..2)
  for (i=0;i<2;i++) {
    brn_.brhaz[0][i]=table.get_gamma_hvh(i+1,2,3)/brn_.width[i];
    brn_.brhaz[1][i]=0.;
  }
  brn_.brhaz[0][2]=0.;
  brn_.brhaz[1][2]=0.;
  // a2 -> a1hi (i=1..3)
  for (i=0;i<3;i++)
    brn_.braha[i]=0.;
  // ai -> Zhj  (i=1,2, j=1..3)
  brn_.brahz[0][0]=table.get_gamma_hvh(3,2,1)/brn_.width[3];
  brn_.brahz[1][0]=table.get_gamma_hvh(3,2,2)/brn_.width[3];
  brn_.brahz[2][0]=0.;
  brn_.brahz[0][1]=0.;
  brn_.brahz[1][1]=0.;
  brn_.brahz[2][1]=0.;
  // h1,h2,h3,a1,a2 -> h+ w-
  brn_.brhcw[0]=table.get_gamma_hvh(1,3,4)/brn_.width[0];
  brn_.brhcw[1]=table.get_gamma_hvh(2,3,4)/brn_.width[1];
  brn_.brhcw[3]=table.get_gamma_hvh(3,3,4)/brn_.width[3];
  brn_.brhcw[2]=0.;
  brn_.brhcw[4]=0.;
  // h1,h2,h3,a1,a2 -> other Higgses
  for (i=0;i<5;i++)
    brn_.brhiggs[i]=0.;
  brn_.brhiggs[1]+=brn_.brhhh[0];
  for (i=0;i<2;i++) {
    brn_.brhiggs[i]+=brn_.brhaa[0][i];
    brn_.brhiggs[i]+=brn_.brhchc[i];
    brn_.brhiggs[i]+=brn_.brhaz[0][i];
    brn_.brhiggs[i]+=brn_.brhcw[i];
  }
  brn_.brhiggs[3]+=brn_.brahz[0][0];
  brn_.brhiggs[3]+=brn_.brahz[1][0];
  brn_.brhiggs[3]+=brn_.brhcw[3];
  for (i=0;i<5;i++)
    for (j=0;j<5;j++)
      for (k=0;k<5;k++)
	brn_.brneu[i][j][k]=0.;
  // A stable Higgs is the same as a Higgs only decaying to neutralion_1
  if (brn_.width[0]<THDM::EPS) brn_.brneu[0][0][0]=1.;
  if (brn_.width[1]<THDM::EPS) brn_.brneu[0][0][1]=1.;
  if (brn_.width[3]<THDM::EPS) brn_.brneu[0][0][3]=1.;
  for (i=0;i<3;i++)
    for (j=0;j<5;j++)
      brn_.brcha[i][j]=0.;
  for (i=0;i<10;i++)
    for (j=0;j<3;j++)
      brn_.brhsq[i][j]=0.;
  for (i=0;i<7;i++)
    for (j=0;j<3;j++)
      brn_.brhsl[i][j]=0.;
  for (i=0;i<6;i++)
    for (j=0;j<2;j++)
      brn_.brasq[i][j]=0.;
  for (i=0;i<3;i++)
    for (j=0;j<2;j++)
      brn_.brasl[i][j]=0.;
  for (i=0;i<5;i++)
    brn_.brsusy[i]=0.;

  double par[24];
  for (i=0;i<24;i++)
    par[i]=0.;
  par[2]=tb;
  double prob[37];
  for (i=0;i<37;i++)
    prob[i]=0.;
  subexp_(par, prob);
  hZ=(prob[3]==0);
  hZ2b=(prob[4]==0);
  hZ2tau=(prob[5]==0);
  hZinv=(prob[6]==0);
  hZ2j=(prob[7]==0);
  hZ2gamma=(prob[8]==0);
  hZ4b=(prob[9]==0);
  hZ4tau=(prob[10]==0);
  hZ2b2tau=(prob[11]==0);
  hA=(prob[12]==0);
  hA4b=(prob[13]==0);
  hA4tau=(prob[14]==0);
  hA2b2tau=(prob[15]==0);
  hA6b=(prob[16]==0);
  hA6tau=(prob[17]==0);
  ZhZjj=(prob[18]==0);

  return hZ&&hZ2b&&hZ2tau&&hZinv&&hZ2j&&hZ2gamma&&hZ4b&&hZ4tau&&hZ2b2tau&&
    hA&&hA4b&&hA4tau&&hA2b2tau&&hA6b&&hA6tau&&ZhZjj;
}
#endif


double Constraints::delta_rho(double mh) {

  double S,T,U,V,W,X;
  double alpha = sm.get_alpha();
  oblique_param(mh, S, T, U, V, W, X);
  return alpha*T;

}


double Constraints::Fdrho(double x, double y) {

  double F = 0.0;

  if ((x>0)&&(y>0)&&(abs(x-y)>delta*x)) {
    F = 0.5*(x+y)-x*y/(x-y)*log(x/y);
  } else {
    F = 0.5*(x+y)-x*(1.-(x-y)/y/2.+pow((x-y)/y,2)/3.-pow((x-y)/y,3)/4.+pow((x-y)/y,4)/5.);
  }

  return F;

}

void Constraints::print_all(double mh_ref) {
  double S,T,U,V,W,X;
  oblique_param(mh_ref,S,T,U,V,W,X);

/*
  bool HpHp,HpHptau,HpHpcs;
  bool test_charged = check_charged(HpHp,HpHptau,HpHpcs);
  bool test_NMSSMTools = true;
#if defined NMSSMTools
  bool hZ,hZ2b,hZ2tau,hZinv,hZ2j,hZ2gamma,hZ4b,hZ4tau,hZ2b2tau,
    hA,hA4b,hA4tau,hA2b2tau,hA6b,hA6tau,ZhZjj;
  test_NMSSMTools = check_NMSSMTools(hZ,hZ2b,hZ2tau,hZinv,hZ2j,hZ2gamma,hZ4b,hZ4tau,hZ2b2tau,
				     hA,hA4b,hA4tau,hA2b2tau,hA6b,hA6tau,ZhZjj);
#endif

#if defined HiggsBounds

  bool test_HiggsBounds = true;
  int HBresult; int chan; double obsratio;
  int ncombined;
  test_HiggsBounds = check_HiggsBounds(HBresult, chan, obsratio, ncombined);
#endif
*/
//  bool test_masses=test_charged&&test_NMSSMTools&&test_HiggsBounds;

  printf("\n");
	printf("Constraints:\n");
  printf(" Tree-level unitarity: %u\n", model.check_unitarity());
  printf(" Perturbativity:       %u\n", model.check_perturbativity());
  printf(" Stability:            %u\n", model.check_stability());
//#if defined HiggsBounds
//  printf(" Mass constraints:     %u\n",test_HiggsBounds);
//  printf("  Charged Higgs        %u (HpHp:%u HpHptau:%u HpHpcs:%u)\n",test_charged,HpHp,HpHptau,HpHpcs);
//#if defined NMSSMTools
//  printf("  NMSSMTools           %u (hZ:%u hZ2b:%u hZ2tau:%u hZinv:%u hZ2j:%u hZ2gamma:%u hZ4b:%u hZ4tau:%u\n                        hZ2b2tau:%u hA:%u hA4b:%u hA4tau:%u hA2b2tau:%u hA6b:%u hA6tau:%u ZhZjj:%u)\n",test_NMSSMTools,hZ,hZ2b,hZ2tau,hZinv,hZ2j,hZ2gamma,hZ4b,hZ4tau,hZ2b2tau,hA,hA4b,hA4tau,hA2b2tau,hA6b,hA6tau,ZhZjj);
//#endif
//  printf("  HiggsBounds          %u (Process %i, see Key.dat, has theory/limit %8.2e combining %i Higgs)\n",test_HiggsBounds,chan,obsratio,ncombined);
//#endif
  printf("\nOblique parameters:\n");
  printf(" S            %12.5e\n", S);
  printf(" T            %12.5e\n", T);
  printf(" U            %12.5e\n", U);
  printf(" V            %12.5e\n", V);
  printf(" W            %12.5e\n", W);
  printf(" X            %12.5e\n", X);
  printf(" Delta_rho    %12.5e\n", delta_rho(91));
  printf(" Delta_amu    %12.5e\n\n", delta_amu());
}
