#include "DecayTable.h"
#include "Util.h"
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>

using namespace std;

// Particle naming
const char *dnames[4] = {" ","d ", "s ", "b "};
const char *unames[4] = {" ","u ", "c ", "t "};
const char *lnames[4] = {" ","e ", "mu", "ta"};
const char *nnames[4] = {" ","ve", "vm", "vt"};
const char *hnames[6] = {" ","h ", "H ", "A ", "H+", "H-"};
const char *vnames[5] = {" ","ga", "Z ", "W+", "W-"};
const int dPDG[4] = {0 ,1 , 3 , 5 };
const int uPDG[4] = {0 ,2 , 4 , 6 };
const int lPDG[4] = {0 ,11, 13, 15};
const int nPDG[4] = {0 ,12, 14, 16};
const int hPDG[5] = {0 ,25, 35, 36, 37};
const int vPDG[4] = {0 ,22, 23, 24};

// Limit on off-shell tails (# GeV above 3-body threshold)
const static double dmtt = 10.;
const static double dmtb = 10.;


DecayTable::DecayTable(THDM mod) {
  set_model(mod);
  qcd_on=true;
  

}

void DecayTable::set_model(THDM mod) {
  model = mod;
  sm = mod.get_SM();
  int i,j,k;
  for (i=1;i<5;i++) {
    gammatot_h[i]=-1.;
    gamma_hgg[i]=-1.;
    gamma_hgaga[i]=-1.;
    gamma_hZga[i]=-1.;
    for (j=1;j<5;j++) {
      gamma_hvv[i][j]=-1.;
      for (k=1;k<5;k++) {
        gamma_uhd[i][j][k]=-1.;
        gamma_hdd[i][j][k]=-1.;
        gamma_huu[i][j][k]=-1.;
        gamma_uhu[i][j][k]=-1.;
        gamma_hdu[i][j][k]=-1.;
        gamma_hll[i][j][k]=-1.;
        gamma_hln[i][j][k]=-1.;
        gamma_hvh[i][j][k]=-1.;
        gamma_hhh[i][j][k]=-1.;
      }
    }
  }
  
 
}


THDM DecayTable::get_model() {
  return model;
}


double DecayTable::get_gamma_uhd(int u, int h, int d) {

  if ((u<1)||(u>3)) return 0.;
  if (h!=4) return 0.;
  if ((d<1)||(d>3)) return 0.;

  if (gamma_uhd[u][h][d]>=0) return gamma_uhd[u][h][d];

  double M =  sm.get_umass_pole(u);
  double m1 = model.get_hmass(h);
  double m2 = sm.get_dmass_pole(d);
  
  if (M<(m1+m2)) {
    gamma_uhd[u][h][d]=0.;
    return gamma_uhd[u][h][d];
  }

  complex <double> cs,cp;
  model.get_coupling_hdu(h,d,u,cs,cp);
  cp = -cp;

  gamma_uhd[u][h][d] =  1./(16.*M_PI)*M*pow(1.-m1*m1/(M*M),2)*(pow(abs(cs),2)+pow(abs(cp),2));
  
  if (qcd_on) {
    double mt = sm.get_qmass_pole(6);
    double mb = sm.get_qmass_pole(5);
    double as = sm.run_alphas_MSbar(mt,mt,mb);

    double qH = pow(m1/M,2);
    
    double K = 1.+as/M_PI*(7.-8.*pow(M_PI,2)/9.-2.*log(1.-qH)+2.*(1.-qH)+(4./9.+2./3.*log(1.-qH))*pow(1.-qH,2));
    
    gamma_uhd[u][h][d] = gamma_uhd[u][h][d]*K;
  }
  

  return gamma_uhd[u][h][d];
}


double DecayTable::get_gamma_hdd(int h, int d1, int d2) {

  if ((h<1)||(h>3)) return 0.;
  if ((d1<1)||(d1>3)) return 0.;
  if ((d2<1)||(d2>3)) return 0.;

  if (gamma_hdd[h][d1][d2]>=0) return gamma_hdd[h][d1][d2];

  double M = model.get_hmass(h);
  double m1p = sm.get_dmass_pole(d1);
  double m2p = sm.get_dmass_pole(d2);
  
  double m1 = sm.get_dmass_MSbar(d1);
  double m2 = sm.get_dmass_MSbar(d2);
  double m1run = sm.get_dmass_MSbar(d1);
  double m2run = sm.get_dmass_MSbar(d2);
  if (m1 > 0) {
      double Qinit = m1;
      if (d1==2) {
       // Special case of strange mass where scale is not ms(ms) but ms(Q_ms=2 GeV)
       Qinit = SM::Q_ms;
      }
      if(sm.b_HD) Qinit = SM::Q_HD;
      m1run = sm.run_qmass_MSbar(m1,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }
  if (m2 > 0) {
      double Qinit = m2;
      if (d2==2) {
       // Special case of strange mass where scale is not ms(ms) but ms(Q_ms=2 GeV)
       Qinit = SM::Q_ms;
      }
      if(sm.b_HD) Qinit = SM::Q_HD;
      m2run = sm.run_qmass_MSbar(m2,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }

  if (M<(m1p+m2p)) {
    gamma_hdd[h][d1][d2] = 0.;
    return gamma_hdd[h][d1][d2];
  }

  complex <double> cs,cp;
  model.get_coupling_hdd(h,d1,d2,cs,cp);

  gamma_hdd[h][d1][d2] = hff_onshell(M,m1p,m1run,m2p,m2run,cs,cp,3,h,false);
  return gamma_hdd[h][d1][d2];
}


double DecayTable::get_gamma_huu(int h, int u1, int u2) {

  if ((h<1)||(h>3)) return 0.;
  if ((u1<1)||(u1>3)) return 0.;
  if ((u2<1)||(u2>3)) return 0.;

  if (gamma_huu[h][u1][u2]>=0) return gamma_huu[h][u1][u2];

  double M = model.get_hmass(h);
  double m1 = sm.get_umass_pole(u1);
  double m2 = sm.get_umass_pole(u2);

  double m1m1 = sm.get_umass_MSbar(u1);
  double m2m2 = sm.get_umass_MSbar(u2);
  double m1run = sm.get_umass_MSbar(u1);
  double m2run = sm.get_umass_MSbar(u2);
  
  if (m1m1 > 0) {
      double Qinit = m1m1;
      if(sm.b_HD) Qinit = SM::Q_HD;
      m1run = sm.run_qmass_MSbar(m1m1,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }
  if (m2m2 > 0) {
      double Qinit = m2m2;
      if(sm.b_HD) Qinit = SM::Q_HD;
      m2run = sm.run_qmass_MSbar(m2m2,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }


  if ((u1<3)&&(u2<3)&&(M<(m1+m2))) {
    gamma_huu[h][u1][u2] = 0.;
    return gamma_huu[h][u1][u2];
  }

  complex <double> cs,cp;
  model.get_coupling_huu(h,u1,u2,cs,cp);
 
 
 
  if ((u1 < 3)||(u2<3)) {
    if (M>(m1+m2)) { 
	    gamma_huu[h][u1][u2] = hff_onshell(M,m1,m1run,m2,m2run,cs,cp,3,h,false);
  	} else {
		gamma_huu[h][u1][u2] = 0.;
	}
    return gamma_huu[h][u1][u2];
  
  } else  {

      double mt_mt = sm.get_qmass_MSbar(6);
      double mb_mb = sm.get_qmass_MSbar(5);
      
      double mb = sm.get_qmass_pole(5);
      double gtop = sm.get_gamma_top();
      
      double mtr = sm.run_qmass_MSbar(mt_mt,mt_mt,M,mt_mt,mb_mb);
      double mW = sm.get_MW(); 

      int dl = 3;
      int du = 3; 

      if (M < m1+mb+mW+dmtt) {
        gamma_huu[h][u1][u2] = 0.;
      } else if (M < m1+m2-dl*gtop) {
        double htt = htt_offshell(M, m1,m2,3,h);
      	gamma_huu[h][u1][u2] = htt;
      } else if (M < m1+m2+du*gtop) {
        double x[4], y[4];
      
        x[0] = m1+m2-(dl+1)*gtop;
        x[1] = m1+m2-dl*gtop;
        x[2] = m1+m2+du*gtop;
        x[3] = m1+m2+(du+1)*gtop;
      
        y[0] = htt_offshell(x[0],m1,m2,3,h);    
        y[1] = htt_offshell(x[1],m1,m2,3,h);    
        y[2] = htt_onshell(x[2],m1,m2,3,h);    
        y[3] = htt_onshell(x[3],m1,m2,3,h);    
 
        double hinter = cubic(M,x,y);
      	gamma_huu[h][u1][u2] = hinter;
      
      
      } else {
        double htt = htt_onshell(M,m1,m2,3,h);
        double hff = hff_onshell(M,mtr,mtr,mtr,mtr,cs,cp,3,h,true);
      
        double R = 2.*m1/M;
        double hinter = interp(R, htt, hff, 0.5);
  
    	gamma_huu[h][u1][u2] = hinter;
      }
  }

  return gamma_huu[h][u1][u2];
}


double DecayTable::get_gamma_hdu(int h, int d, int u) {

  if ((h<4)||(h>4)) return 0.;
  if ((d<1)||(d>3)) return 0.;
  if ((u<1)||(u>3)) return 0.;

  if (gamma_hdu[h][d][u]>=0) return gamma_hdu[h][d][u];

  double M = model.get_hmass(h);
  double m1p = sm.get_dmass_pole(d);
  double m2p = sm.get_umass_pole(u);
  
  double m1 = sm.get_dmass_MSbar(d);
  double m2 = sm.get_umass_MSbar(u);
  double m1run = sm.get_dmass_MSbar(d);
  double m2run = sm.get_umass_MSbar(u);
  
  double mW = sm.get_MW();
  
  if (m1 > 0) {
      double Qinit = m1;
      if (d==2) {
       // Special case of strange mass where scale is not ms(ms) but ms(Q_ms=2 GeV)
       Qinit = SM::Q_ms;
      }
      if(sm.b_HD) Qinit = SM::Q_HD;
      m1run = sm.run_qmass_MSbar(m1,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }
  if (m2 > 0) {
      double Qinit = m2;
      if(sm.b_HD) Qinit = SM::Q_HD;
      m2run = sm.run_qmass_MSbar(m2,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  }

  complex <double> cs,cp;
  model.get_coupling_hdu(h,d,u,cs,cp);
  
  if (h<4) {
	gamma_hdu[h][d][u] = hff_onshell(M,m1p,m1run,m2p,m2run,cs,cp,3,h,false);
	return gamma_hdu[h][d][u];
  }

  if (u<3) {
	gamma_hdu[h][d][u] = hff_onshell(M,m1p,m1run,m2p,m2run,cs,cp,3,h,false);
	return gamma_hdu[h][d][u];
  }
 
	double mb = sm.get_qmass_pole(5);
	double gtop = sm.get_gamma_top();

	int dl = 0;
	int du = 1; 

	if (M < m1p+mb+mW+dmtb) {
		gamma_hdu[h][d][u] = 0.;
		return gamma_hdu[h][d][u];
	}
	
	if (M < m1p+m2p-dl*gtop) {
		gamma_hdu[h][d][u] = htb_offshell(M,m1p,m1run,m2p,m2run,cs,cp,3);
		return gamma_hdu[h][d][u];
	}
	
	if (M < m1p+m2p+du*gtop) {
		double x[4], y[4];

		x[0] = m1p+m2p-(dl+1)*gtop;
		x[1] = m1p+m2p-dl*gtop;
		x[2] = m1p+m2p+du*gtop;
		x[3] = m1p+m2p+(du+1)*gtop;

		y[0] = htb_offshell(x[0],m1p,m1run,m2p,m2run,cs,cp,3);    
		y[1] = htb_offshell(x[1],m1p,m1run,m2p,m2run,cs,cp,3);    
		y[2] = hpff_onshell(x[2],m1p,m1run,m2p,m2run,cs,cp,3,h);    
		y[3] = hpff_onshell(x[3],m1p,m1run,m2p,m2run,cs,cp,3,h);    

		double hinter = cubic(M,x,y);
		gamma_hdu[h][d][u] = hinter;
		return gamma_hdu[h][d][u];
	}

	double htd1 = hpff_onshell(M,m1p,m1run,m2p,m2run,cs,cp,3,h);
	double htd2 = hff_onshell(M,m1p,m1run,m2p,m2run,cs,cp,3,h,false);

	double hinter = htd1;

	double R = (m1p+m2p)/M;
	hinter = interp(R, htd1, htd2, 2.);  	
	gamma_hdu[h][d][u]=hinter;
	    
	return gamma_hdu[h][d][u];
	
}


double DecayTable::get_gamma_hll(int h, int l1, int l2) {

  if ((h<1)||(h>3)) return 0.;
  if ((l1<1)||(l1>3)) return 0.;
  if ((l2<1)||(l2>3)) return 0.;

  if (gamma_hll[h][l1][l2]>=0) return gamma_hll[h][l1][l2];

  double M = model.get_hmass(h);
  double m1 = sm.get_lmass_pole(l1);
  double m2 = sm.get_lmass_pole(l2);

  if (M<(m1+m2)) {
    gamma_hll[h][l1][l2] = 0.;
    return gamma_hll[h][l1][l2];
  }

  complex <double> cs,cp;
  model.get_coupling_hll(h,l1,l2,cs,cp);

  gamma_hll[h][l1][l2] = hff_onshell(M,m1,m1,m2,m2,cs,cp,1,h,false);

  return gamma_hll[h][l1][l2];
}


double DecayTable::get_gamma_hln(int h, int l, int n) {

  if ((h<1)||(h>4)) return 0.;
  if ((l<1)||(l>3)) return 0.;
  if ((n<1)||(n>3)) return 0.;

  if (gamma_hln[h][l][n]>=0) return gamma_hln[h][l][n];

  double M = model.get_hmass(h);
  double m1 = sm.get_lmass_pole(l);
  double m2 = 0.;

  if (M<(m1+m2)) {
    gamma_hln[h][l][n] = 0.;
    return gamma_hln[h][l][n];
  }

  complex <double> cs,cp;
  model.get_coupling_hln(h,l,n,cs,cp);

  
  gamma_hln[h][l][n] = hff_onshell(M,m1,m1,m2,m2,cs,cp,1,h,false);
  return gamma_hln[h][l][n];
}


double DecayTable::get_gamma_hgg(int h) {

  if (!qcd_on) return 0.;

  if ((h<1)||(h>3)) return 0.;

  if (gamma_hgg[h]>=0) return gamma_hgg[h];
  
  gamma_hgg[h] = hgg(h);
  return gamma_hgg[h];
}


double DecayTable::get_gamma_hgaga(int h) {

  if (!qcd_on) return 0.;

  if ((h<1)||(h>3)) return 0.;

  if (gamma_hgaga[h]>=0) return gamma_hgaga[h];
  
  gamma_hgaga[h] = hgaga(h);
  return gamma_hgaga[h];
}


double DecayTable::get_gamma_hZga(int h) {

  if ((h<1)||(h>3)) return 0.;

  if (gamma_hZga[h]>=0) return gamma_hZga[h];
  
  gamma_hZga[h] = hZga(h);
  return gamma_hZga[h];
}


double DecayTable::get_gamma_hvv(int h, int V) {

  if ((h<1)||(h>=4)) return 0.;
  if ((V<1)||(V>3)) return 0.;

  if ((h<=3)&&(V==1)) {
    return get_gamma_hgaga(h);
  }

  if (gamma_hvv[h][V]>=0) return gamma_hvv[h][V];
  
  double M = model.get_hmass(h);
  
  gamma_hvv[h][V] = hvv_all(h,V,M);
  return gamma_hvv[h][V];

}


double DecayTable::get_gamma_hvh(int H, int V, int h) {

  if ((H<1)||(H>4)||(h==H)) return 0.;
  if ((V<1)||(V>3)) return 0.;
  if ((h<1)||(h>4)) return 0.;
  
  if ((H==4)&&(V!=3)) return 0;
  if ((H!=4)&&(h==4)&&(V!=3)) return 0;
  
  if (V==1) return 0;

  if (gamma_hvh[H][V][h]>=0) return gamma_hvh[H][V][h];

  double M  = model.get_hmass(H);
  double m1  = model.get_hmass(h);
  double m2 = sm.get_vmass(V);
  
  double GV = sm.get_gamma_V(V);
  gamma_hvh[H][V][h] = 0.;
  
  if (M<=m1) {
    gamma_hvh[H][V][h] = 0.;
    return gamma_hvh[H][V][h];
  }

  if (M>(m1+m2+5*GV)) {
    gamma_hvh[H][V][h] = hvh_onshell(H,V,h,M);
    return gamma_hvh[H][V][h];
  } 
    
  if (M>m1+2*GV) {
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
        
    double result, error;
    
    integration_params ip;
    ip.M = M;
    ip.m1 = m1;
    ip.m2 = m2;
    ip.gamma = GV;
  
    gsl_function F;
    F.function = &hvh_fcn;
    F.params = &ip;
      
    double k = pow(m1/M,2);
    double imin = 0.;
    double imax = 1.-k;
  
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    int status = gsl_integration_qags (&F,imin,imax,0,1e-5,1000,
			  w, &result, &error); 
    if (status) 
      if (status!=GSL_EROUND) {
        printf("GSL integration warning in H%1d -> VH%1d (off-shell). Please check result.\n", H, h);
        if (EXIT_ON_GSL_ERROR) exit(-1);
      }
    gsl_set_error_handler(old_handler);
    old_handler = NULL;
    gsl_integration_workspace_free (w);
    w = NULL;
  
    complex <double> c;
    model.get_coupling_vhh(V,H,h,c);
  
    double dV = 0.;
    double stw = sm.get_sintw();
    double ctw = sm.get_costw();
    double GF = sm.get_GF();
    double MW = sm.get_MW();

    if (V==1) {
      gamma_hvh[H][V][h] = 0.;
      return gamma_hvh[H][V][h];
    } else if(V==2) {
      dV = 3./pow(ctw,2)*(7./12.-10./9.*pow(stw,2)+40./27.*pow(stw,4));
    } else if(V==3) {
      dV = 3.;
    } else {
      gamma_hvh[H][V][h] = 0.;
      return gamma_hvh[H][V][h];
    }

    double KHHV = 3.*GF/(16.*sqrt(2)*pow(M_PI,3))*pow(MW,2)*pow(abs(c),2)*M*dV;
    gamma_hvh[H][V][h] = KHHV*result;
    
    // Count both H+W- and H-W+ final states in partial width to charged states
    // also done in hvh_onshell
    if ((h==4)&&(V==3)) gamma_hvh[H][V][h] = 2*gamma_hvh[H][V][h];

    if (M>(m1+m2)) {
      double G2 = hvh_onshell(H,V,h,M); 
      gamma_hvh[H][V][h] = max(gamma_hvh[H][V][h],G2);
    }

    return gamma_hvh[H][V][h];
  }


  return gamma_hvh[H][V][h];
}


double DecayTable::get_gamma_hhh(int h, int h1, int h2) {

  if ((h<1)||(h>4)) return 0.;
  if ((h1<1)||(h1>4)) return 0.;
  if ((h2<1)||(h2>4)) return 0.;

  if (gamma_hhh[h][h1][h2]>=0) return gamma_hhh[h][h1][h2];
  double Sf = 1.;
  if ((h1==h2)&&(h1!=4)) Sf = 0.5;

  double M = model.get_hmass(h);
  double m1 = model.get_hmass(h1);
  double m2 = model.get_hmass(h2);

  complex <double> c;
  model.get_coupling_hhh(h,h1,h2,c);
  
  gamma_hhh[h][h1][h2] = 0.;

  if (M>(m1+m2)) {
    double M2    = pow(abs(c),2);
    gamma_hhh[h][h1][h2] = Sf/(8.*M_PI)*M2*PS2(M,m1,m2);
    return gamma_hhh[h][h1][h2];
  }

  return gamma_hhh[h][h1][h2];
}


double DecayTable::get_gammatot_top() {

  double gtot = 0.;
  
  gtot += sm.get_gamma_top();
  
  gtot += get_gamma_uhd(3,4,1);
  gtot += get_gamma_uhd(3,4,2);
  gtot += get_gamma_uhd(3,4,3);

  gtot += get_gamma_uhu(3,1,1);
  gtot += get_gamma_uhu(3,1,2);
  gtot += get_gamma_uhu(3,2,1);
  gtot += get_gamma_uhu(3,2,2);
  gtot += get_gamma_uhu(3,3,1);
  gtot += get_gamma_uhu(3,3,2);

  return gtot;

}


double DecayTable::get_gammatot_h(int h) {

  if (h>4) return 0.;

  if (gammatot_h[h]>=0) return gammatot_h[h];

  gammatot_h[h] = 0.;
  
  // Fermionic modes
  for (int i=1;i<4;i++) {
    for (int j=1;j<4;j++) {
       gammatot_h[h]+=get_gamma_hdd(h,i,j);
       gammatot_h[h]+=get_gamma_huu(h,i,j);
       gammatot_h[h]+=get_gamma_hdu(h,i,j);
       gammatot_h[h]+=get_gamma_hll(h,i,j);
       gammatot_h[h]+=get_gamma_hln(h,i,j);
    }
  }
  
  // Vector bosons
  for (int i=1;i<4;i++) {
     gammatot_h[h]+=get_gamma_hvv(h,i);
  }

  // Z gamma
   gammatot_h[h] += get_gamma_hZga(h);

  // Gluons
   gammatot_h[h] += get_gamma_hgg(h);

  // H -> VH
  for (int i=1;i<4;i++) {
    for (int j=1;j<=4;j++) {
       gammatot_h[h]+=get_gamma_hvh(h,i,j);
    }
  }

  // Higgses
  for (int i=1;i<5;i++) {
    for (int j=1;j<5;j++) {
       gammatot_h[h]+=get_gamma_hhh(h,i,j);
    }
  }

  return gammatot_h[h];
}


double  DecayTable::get_gammatot_v(int v) {
  return sm.get_gamma_V(v);
}



void DecayTable::print_decay_LesHouches(FILE* output, double br, int id1, int id2) {
  if (br>0) 
    fprintf(output,"     % 16.8e     %1i     %3i   %3i\n",br,2,id1,id2);
}


void DecayTable::print_decay(const char *h, const char *id1, const char *id2, double g, double br) { 
  if (br>THDM::EPS)
    printf("%2s -> %2s %2s %12.3e   %12.3e\n",h,id1,id2,g,br);
}


void DecayTable::print_top_decays() {
  double gtot = get_gammatot_top();
  double gt[12],br[12];

  gt[0] = sm.get_gamma_tWd(1);
  gt[1] = sm.get_gamma_tWd(2);
  gt[2] = sm.get_gamma_tWd(3);
  gt[3] = get_gamma_uhd(3,4,1);
  gt[4] = get_gamma_uhd(3,4,2);
  gt[5] = get_gamma_uhd(3,4,3);
  gt[6] = get_gamma_uhu(3,1,1);
  gt[7] = get_gamma_uhu(3,1,2);
  gt[8] = get_gamma_uhu(3,2,1);
  gt[9] = get_gamma_uhu(3,2,2);
  gt[10] = get_gamma_uhu(3,3,1);
  gt[11] = get_gamma_uhu(3,3,2);

  for (int i=0;i<12;i++) {
    br[i]=gt[i]/gtot;
  }

  printf("\nDecay table for %s\n", unames[3]);
  printf("Total width:%12.3e GeV      BR\n", gtot);
  print_decay(unames[3],vnames[3],dnames[1],gt[0],br[0]);
  print_decay(unames[3],vnames[3],dnames[2],gt[1],br[1]);
  print_decay(unames[3],vnames[3],dnames[3],gt[2],br[2]);
  print_decay(unames[3],hnames[4],dnames[1],gt[3],br[3]);
  print_decay(unames[3],hnames[4],dnames[2],gt[4],br[4]);
  print_decay(unames[3],hnames[4],dnames[3],gt[5],br[5]);
  print_decay(unames[3],hnames[1],unames[1],gt[6],br[6]);
  print_decay(unames[3],hnames[1],unames[2],gt[7],br[7]);
  print_decay(unames[3],hnames[2],unames[1],gt[8],br[8]);
  print_decay(unames[3],hnames[2],unames[2],gt[9],br[9]);
  print_decay(unames[3],hnames[3],unames[1],gt[10],br[10]);
  print_decay(unames[3],hnames[3],unames[2],gt[11],br[11]);

   printf("---------------------------------------\n");
}


void DecayTable::print_top_decays_LesHouches(FILE* output, bool full) {
  double gtot = get_gammatot_top();
  double gt[12],br[12];

  gt[0] = sm.get_gamma_tWd(1);
  gt[1] = sm.get_gamma_tWd(2);
  gt[2] = sm.get_gamma_tWd(3);
  gt[3] = get_gamma_uhd(3,4,1);
  gt[4] = get_gamma_uhd(3,4,2);
  gt[5] = get_gamma_uhd(3,4,3);
  gt[6] = get_gamma_uhu(3,1,1);
  gt[7] = get_gamma_uhu(3,1,2);
  gt[8] = get_gamma_uhu(3,2,1);
  gt[9] = get_gamma_uhu(3,2,2);
  gt[10] = get_gamma_uhu(3,3,1);
  gt[11] = get_gamma_uhu(3,3,2);

  for (int i=0;i<12;i++) {
    br[i]=gt[i]/gtot;
  }

  fprintf(output,"DECAY  6   % 16.8e   # top decays\n",gtot);
  if (full) {
    fprintf(output,"#            BR          NDA    ID1   ID2\n");
    print_decay_LesHouches(output,br[0],vPDG[3],dPDG[1]);
    print_decay_LesHouches(output,br[1],vPDG[3],dPDG[2]);
    print_decay_LesHouches(output,br[2],vPDG[3],dPDG[3]);
    print_decay_LesHouches(output,br[3],hPDG[4],dPDG[1]);
    print_decay_LesHouches(output,br[4],hPDG[4],dPDG[2]);
    print_decay_LesHouches(output,br[5],hPDG[4],dPDG[3]);
    print_decay_LesHouches(output,br[6],hPDG[1],uPDG[1]);
    print_decay_LesHouches(output,br[7],hPDG[1],uPDG[2]);
    print_decay_LesHouches(output,br[8],hPDG[2],uPDG[1]);
    print_decay_LesHouches(output,br[9],hPDG[2],uPDG[2]);
    print_decay_LesHouches(output,br[10],hPDG[3],uPDG[1]);
    print_decay_LesHouches(output,br[11],hPDG[3],uPDG[2]);
  }
}


void DecayTable::print_decays(int h) {
  print_decays(0,h,true,false);
}


void DecayTable::print_width(int h) {
  printf(" Total width for %s: %10.3e GeV\n", hnames[h], get_gammatot_h(h));
}


void DecayTable::print_decays_LesHouches(FILE* output, int h, bool full) {
  print_decays(output,h,full,true);
}


double DecayTable::br(double dG, double G) {

  double BR = 0.;

  if (G>0.) {
    BR = dG/G;
    if (BR<THDM::EPS) BR = 0.;
  }
  
  return BR;
}


void DecayTable::print_decays(FILE* output,int h, bool full, bool les) {

  if ((h<1)||(h>4)) return;

  double gtot = get_gammatot_h(h);

  if (les) {
    if (h==1) 
      fprintf(output,"DECAY  25   % 16.8e   # h1 decays, lightest CP-even Higgs\n",gtot);
    else if (h==2)
      fprintf(output,"DECAY  35   % 16.8e   # h2 decays, heaviest CP-even Higgs\n",gtot);
    else if (h==3)
      fprintf(output,"DECAY  36   % 16.8e   # h3 decays, CP-odd Higgs\n",gtot);
    else if (h==4)
      fprintf(output,"DECAY  37   % 16.8e   # Charged Higgs decays\n",gtot);
    fprintf(output,"#            BR          NDA    ID1   ID2\n");
  } else {
    printf("\nDecay table for %s\n", hnames[h]);
    printf("Total width:%12.3e GeV      BR\n", gtot);
  }

  if (!full) return;

  double gdd[4][4];
  double guu[4][4];
  double gdu[4][4];
  double gll[4][4];
  double gln[4][4];
  double gvv[4];
  double gvh[4][5];
  double ghh[5];
  double ghZga;
  double ghgg;
  double brdd[4][4];
  double bruu[4][4];
  double brdu[4][4];
  double brll[4][4];
  double brln[4][4];
  double brvv[4];
  double brvh[4][5];
  double brhh[5];
  double brhZga = 0.;
  double brhgg = 0.;
  
  // Fermion decay modes
  for (int i=1;i<4;i++) {
    for (int j=1;j<4;j++) {
      gdd[i][j]=get_gamma_hdd(h,i,j);
      guu[i][j]=get_gamma_huu(h,i,j);
      gdu[i][j]=get_gamma_hdu(h,i,j);
      gll[i][j]=get_gamma_hll(h,i,j);
      gln[i][j]=get_gamma_hln(h,i,j);
      brdd[i][j] = br(gdd[i][j],gtot);
      bruu[i][j] = br(guu[i][j],gtot);
      brdu[i][j] = br(gdu[i][j],gtot);
      brll[i][j] = br(gll[i][j],gtot);
      brln[i][j] = br(gln[i][j],gtot);
    }
  }

  // Vector bosons
  for (int i=1;i<4;i++) {
    gvv[i]=get_gamma_hvv(h,i);
    brvv[i]=br(gvv[i],gtot);
    for (int j=1;j<5;j++) {
      gvh[i][j]=get_gamma_hvh(h,i,j);
      brvh[i][j]=br(gvh[i][j],gtot);
    }
  }

  // Z gamma
  ghZga = get_gamma_hZga(h);
  brhZga = br(ghZga,gtot);

  // Gluons
  ghgg = get_gamma_hgg(h);
  brhgg = br(ghgg,gtot);

  for (int i=1;i<=4;i++) {
    ghh[i]=get_gamma_hhh(h,i,i);
    brhh[i]=br(ghh[i],gtot);
  }

  if (h==4) {
    for (int j=1;j<4;j++) {
      for (int i=1;i<4;i++) {
        if (les)
	  print_decay_LesHouches(output,brdu[i][j],uPDG[j],-dPDG[i]);
	else 
	  print_decay(hnames[h],unames[j],dnames[i],gdu[i][j],brdu[i][j]);
      }
    }
    for (int i=1;i<4;i++) {
      for (int j=1;j<4;j++) {
        if (les)
	  print_decay_LesHouches(output,brln[i][j],-lPDG[i],nPDG[j]);
	else 
	  print_decay(hnames[h],lnames[i],nnames[j],gln[i][j],brln[i][j]);
      }
    }
    for (int i=2;i<4;i++) {
      for (int j=1;j<=4;j++) {
        int sgn = 1;
        if ((i==3)&&(h!=4)) sgn = -1;
        if (les)
	  print_decay_LesHouches(output,brvh[i][j],vPDG[i],sgn*hPDG[j]);
	else 
	  print_decay(hnames[h],vnames[i],hnames[j],gvh[i][j],brvh[i][j]);
      }
    }
  } else {
    for (int i=1;i<4;i++) {
      for (int j=1;j<4;j++) {
	if (les)
	  print_decay_LesHouches(output,brdd[i][j],dPDG[i],-dPDG[j]);
	else 
	  print_decay(hnames[h],dnames[i],dnames[j],gdd[i][j],brdd[i][j]);
	if (les)
	  print_decay_LesHouches(output,bruu[i][j],uPDG[i],-uPDG[j]);
	else 
	  print_decay(hnames[h],unames[i],unames[j],guu[i][j],bruu[i][j]);
      }
    }
    for (int i=1;i<4;i++) {
      for (int j=1;j<4;j++) {
	if (les)
	  print_decay_LesHouches(output,brll[i][j],lPDG[i],-lPDG[j]);
	else 
	  print_decay(hnames[h],lnames[i],lnames[j],gll[i][j],brll[i][j]);
      }    
    }
    for (int i=1;i<4;i++) {
      int sgn = 1;
      if (i==3) sgn = -1;
      if (les) {
	print_decay_LesHouches(output,brvv[i],vPDG[i],sgn*vPDG[i]);
      } else {
	print_decay(hnames[h],vnames[i],vnames[i+(sgn==-1)],gvv[i],brvv[i]);
      }
    }

    if (les) 
      print_decay_LesHouches(output,brhZga,23,22);
    else
      print_decay(hnames[h],"Z ","ga",ghZga,brhZga);

    if (les) 
      print_decay_LesHouches(output,brhgg,21,21);
    else
      print_decay(hnames[h],"g ","g ",ghgg,brhgg);

    for (int i=1;i<=4;i++) {
      int sgn = 1;
      if (i==4) sgn = -1;
      if (les) {
        print_decay_LesHouches(output,brhh[i],hPDG[i],sgn*hPDG[i]);
      } else {
	print_decay(hnames[h],hnames[i],hnames[i+(sgn==-1)],ghh[i],brhh[i]);
      }
    }
    for (int i=1;i<4;i++) {
      for (int j=1;j<=3;j++) {
        if (les) {
	  print_decay_LesHouches(output,brvh[i][j],vPDG[i],hPDG[j]);
        } else 
	  print_decay(hnames[h],vnames[i],hnames[j],gvh[i][j],brvh[i][j]);
      }
    }

    // Special case with H+ W- charged conjugate final state
    if (les) {
      print_decay_LesHouches(output,0.5*brvh[3][4],vPDG[3],-hPDG[4]);
      print_decay_LesHouches(output,0.5*brvh[3][4],hPDG[4],-vPDG[3]);
    } else {
      print_decay(hnames[h],vnames[3],hnames[5],0.5*gvh[3][4],0.5*brvh[3][4]);
      print_decay(hnames[h],hnames[4],vnames[4],0.5*gvh[3][4],0.5*brvh[3][4]);
    }

  }
  
  if (!les) printf("---------------------------------------\n");
}

void DecayTable::set_qcd(bool set) {
  qcd_on=set;
}

double DecayTable::hvv_offshell(int h, int V,double M) {

  double m = sm.get_vmass(V);

  double GV = sm.get_gamma_V(V);
  double G = 0.;
 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
  double result, error;
    
  integration_params ip;
  ip.M = M;
  ip.m1 = m;
  ip.m2 = m;
  ip.gamma = GV;
  
  gsl_function F;
  F.function = &hvv_fcn;
  F.params = &ip;
      
  double k = pow(m/M,2);
  double imin = 0.;
  double imax = 1.-k;
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F,imin,imax,0,1e-5,1000,
			w, &result, &error); 
    
  if (status) 
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H%1d -> VV (off-shell). Please check result.\n", h);
      if (EXIT_ON_GSL_ERROR) exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;

  complex <double> c;
  model.get_coupling_vvh(V,V,h,c);
  
  double dV = 0.;
  double stw = sm.get_sintw();
  double GF = sm.get_GF();

  if (V==1) {
    return 0.;
  } else if(V==2) {
    dV = 3.*(7./12.-10./9.*pow(stw,2)+40./27.*pow(stw,4));
// Note that there is a typo in Djouadis review for the last term! (40/9 instead of 40/27)
  } else if(V==3) {
    dV = 3.;
  } else return 0.;

  double KHVV = 3.*pow(abs(c),2)*GF/(64.*sqrt(2.)*pow(M_PI,3))*M*dV;
  G = KHVV*result;


  return G;

}

double DecayTable::hvv_all(int h, int V,double M) {

  double m = sm.get_vmass(V);

  double GV = sm.get_gamma_V(V);
  double G = 0.;
 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
  double result, error;
    
  integration_params2 ip;
  ip.M = M;
  ip.m = m;
  ip.gamma = GV;
  
  gsl_function F;
  F.function = &hvv_fcn1;
  F.params = &ip;
      
  double imin = 0.;
  double imax = 1.;
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F,imin,imax,0,1e-5,1000,
			w, &result, &error); 
    
  if (status) 
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H%1d -> VV (off-shell). Please check the result.\n", h);
      if (EXIT_ON_GSL_ERROR) exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;

  complex <double> c;
  model.get_coupling_vvh(V,V,h,c);
  
  double dV = 0.;

  if (V==1) {
    return 0.;
  } else if(V==2) {
    dV = 0.5;
  } else if(V==3) {
    dV = 1.;
  } else return 0.;

  double KHVV = dV*pow(abs(c),2)/(64.*M_PI)*pow(M,3)/pow(m,4);
  G = KHVV*result;

  return G;

}


double DecayTable::hvv_onshell(int h, int V, double M) {
  double m = sm.get_vmass(V);
  double Sf = 0.;

  if (V==1) {
    return 0.;
  } else if(V==2) {
    Sf = 0.5;
  } else if(V==3) {
    Sf = 1.;
  } else return 0.;
  
  if (M<2.*m) return 0.0;

  complex <double> c;
  model.get_coupling_vvh(V,V,h,c);

  double M2 = 0.;
  
  if (m>0) {
    double x = pow(M/m,2);
    M2    = 4.*pow(abs(c),2)*(3.-x+x*x/4.);
  } else {
    M2    = 16.*pow(abs(c),2);
  }

  double G = Sf/(32.*M_PI)*M2*PS2(M,m,m);

  return G;
}


double DecayTable::hvh_onshell(int H, int V, int h, double M) {
  double m1 = sm.get_vmass(V);
  double m2  = model.get_hmass(h);

  if (M<m1+m2) return 0.;

  complex <double> c;
  model.get_coupling_vhh(V,H,h,c);
  
  int dV = 1;

  if (V==1) {
    return 0.;
  } 
  if ((V==3) && (h==4)) {
    dV = 2;
  }
  

  double G = pow(abs(c),2)/(16.*M_PI)*pow(m1,2)/M*sqrt(L(pow(m2,2),pow(m1,2),pow(M,2)))*L(pow(m2,2),pow(M,2),pow(m1,2))*dV;

  return G;
}


double DecayTable::PS2(double M, double m1, double m2) {
  double absp1 = sqrt((pow(M,2)-pow(m1+m2,2))*(pow(M,2)-pow(m1-m2,2)))/(2.*M);
  double PS = absp1/pow(M,2);

  return PS;
}



double DecayTable::hpff_onshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc, int h) {

  double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5))/M_PI;

  complex <double> cd = cs+cp;
  complex <double> cu = cs-cp;
  
  double cdcu = real(cd*cu);
  double cd2 = pow(abs(cd),2);
  double cu2 = pow(abs(cu),2);
  
  if (m1run > 0) {
  	cd2 = cd2*pow(m1/m1run,2);
  }
  
  if (m2run > 0) {
  	cu2 = cu2*pow(m2/m2run,2);
  	
  }
  
  if ((m1run>0)&&(m2run>0)) {
	cdcu=cdcu*m1*m2/(m1run*m2run);
  }
  
  double qd = pow(m1/M,2);
  double qu = pow(m2/M,2);
  double Q = 1.-qu-qd;
  double q = sqrt(qu*qd);
  double lam = L(pow(m1,2), pow(m2,2), pow(M,2));
  double sqL = sqrt(lam);

  double D_udp = 0.;
  double D_dup = 0.;
  double D_udm = 0.;
  
  if ((Nc==3)&&qcd_on) {

     double xd = 2.*qd/(Q+sqL);       
     double xu = 2.*qu/(Q+sqL);

     if(m1 > 0) {
       D_dup = DHp(qd,qu,xd,xu,sqL);      
     }
     
     if(m2 > 0) {
       D_udp = DHp(qu,qd,xu,xd,sqL);           
     }   

     if(q > 0) {
        D_udm = DHm(qu,qd,xu,xd,sqL);
     }

  }

  double G = Nc*M/(16.*M_PI)*sqL*(Q*(cu2*(1.+4./3.*as*D_udp)+cd2*(1.+4./3.*as*D_dup))+4*q*cdcu*(1.+4./3.*as*D_udm));

  return G;
}


double DecayTable::hff_onshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc, int h, bool tt=false) {

  if (M<(m1+m2)) return 0.;

  double M2    = 2.*(pow(M,2)-pow(m1+m2,2))*pow(abs(cs),2)+
    2.*(pow(M,2)-pow(m1-m2,2))*pow(abs(cp),2);

  double G = Nc*1./(8.*M_PI)*M2*PS2(M,m1,m2);

  complex <double> cst,cpt;
  model.get_coupling_huu(h,3,3,cst,cpt);
  double mt = sm.get_umass_MSbar(3);
  double Qinit = mt;
  if(sm.b_HD) Qinit = SM::Q_HD;
  double mtrun = sm.run_qmass_MSbar(mt,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
  
  // Apply QCD corrections for decay to quarks
  if ((Nc==3)&&qcd_on) {
    if (h<3) {
      int Nf = sm.get_Nactivef(M);
      double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5))/M_PI;
      double K = 1.+5.67*as+(35.94-1.36*Nf)*pow(as,2);      
      // massive corrections from hep-ph/9505358 only for q1=q2
      if((abs(cs)>0.)&&(m1run-m2run<THDM::EPS)) K=K+m1run/mtrun*real(cst*conj(cs))/real(cs*conj(cs))*(1.57-2./3.*log(pow(M,2)/pow(sm.get_qmass_pole(6),2)))*pow(as,2);
      if((abs(cs)>0.)&&(m1run-m2run<THDM::EPS)&&(m1run>0.)&&(m2run>0.)&&(!tt)) K=K+sqrt(m1run*m2run)/mtrun*real(cst*conj(cs))/real(cs*conj(cs))*(1./9.*pow(log(pow(m1run/M,2)),2))*pow(as,2);
      if (K<0) K=0;
      G = G*K;   
    } else if (h==3) {
      int Nf = sm.get_Nactivef(M);
      double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5))/M_PI;
      double K = 1.+5.67*as+(35.94-1.36*Nf)*pow(as,2);      
      // massive corrections from hep-ph/9505358 only for q1=q2
      if((abs(cp)>0.)&&(m1run-m2run<THDM::EPS)) K=K+sqrt(m1run*m2run)/mtrun*real(cpt*conj(cp))/real(cp*conj(cp))*(3.83-log(pow(M,2)/pow(sm.get_qmass_pole(6),2)))*pow(as,2);
      if((abs(cp)>0.)&&(m1run-m2run<THDM::EPS)&&(m1run>0.)&&(m2run>0)&&(!tt)) K=K+m1run/mtrun*real(cpt*conj(cp))/real(cp*conj(cp))*(1./6.*pow(log(pow(m1run/M,2)),2))*pow(as,2);     
      if (K<0) K=0;
      G = G*K;
    } else if (h==4) {
      double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5))/M_PI;
      double K = 1.+5.67*as;
      if (K<0) K=0;
      G = G*K;
    }
  }

  return G;
}


double DecayTable::get_gamma_uhu(int u1, int h, int u2) {

  if (u1!=3) return 0.;
  if ((h<1)||(h>3)) return 0.;
  if ((u2<1)||(u2>=u1)) return 0.;

  if (gamma_uhu[u1][h][u2]>=0) return gamma_uhu[u1][h][u2];

  double M =  sm.get_umass_pole(u1);
  double m1 = model.get_hmass(h);
  double m2 = sm.get_umass_pole(u2);
  
  if (M<(m1+m2)) {
    gamma_uhu[u1][h][u2]=0.;
    return gamma_uhu[u1][h][u2];
  }

  complex <double> cs,cp;
  model.get_coupling_huu(h,u1,u2,cs,cp);
 
  double x1 = m1/M;
  double x2 = m2/M; 	


  gamma_uhu[u1][h][u2] = 0;

  if (h<3) {
  	gamma_uhu[u1][h][u2] = real(cs*conj(cs))*M/(16*M_PI)*((1+x2)*(1+x2)-x1*x1)*sqrt(1-(x1+x2)*(x1+x2))*sqrt(1-(x1-x2)*(x1-x2));
  } else if (h==3) {
 	gamma_uhu[u1][h][u2] = real(cp*conj(cp))*M/(16*M_PI)*((1-x2)*(1-x2)-x1*x1)*sqrt(1-(x1+x2)*(x1+x2))*sqrt(1-(x1-x2)*(x1-x2));
  }
  
//  printf("uhu %d %d %d %16.8E %16.8E %16.8E %16.8E %16.8E\n", u1, h, u2, gamma_uhu[u1][h][u2], pow(1-x2,2), pow(1+x2,2),pow(1+x2,2)/(pow(1-x2,2)),x2);
         
  return gamma_uhu[u1][h][u2];
}


double DecayTable::htt_onshell(double M, double m1, double m2,  int Nc, int h) {

  complex <double> cs, cp;
  complex <double> cs_mtp, cp_mtp;
  model.get_coupling_huu(h,3, 3,cs,cp);

  double mtp = sm.get_qmass_pole(6);
  double mt_mt = sm.get_qmass_MSbar(6);
  double mbp = sm.get_qmass_pole(5);
  double mtr = sm.run_qmass_MSbar(mt_mt,mt_mt,M,mtp,mbp);

  cs_mtp = cs*mtp/mtr;
  cp_mtp = cp*mtp/mtr;
  
  double M2    = 2.*(pow(M,2)-pow(m1+m2,2))*pow(abs(cs_mtp),2)+
    2.*(pow(M,2)-pow(m1-m2,2))*pow(abs(cp_mtp),2);

  double G = Nc*1./(8.*M_PI)*M2*PS2(M,m1,m2);
  
  // Threshold QCD corrections for tt decay, with full mass dependence
  if (Nc==3&&qcd_on) {
    double K = 1;
    if (h!=4) {
      double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      
      double mf = m1;
      double b = sqrt(1.-4.*pow(mf/M,2));
      double L = log((1.+b)/(1.-b));
      double x = (1.-b)/(1.+b);
      double dL1 =  gsl_sf_dilog(x);
      double dL2 =  gsl_sf_dilog(-x);

      double b2 = pow(b, 2);
      
      // Threshold corrections for tt decay [hep-ph/0503172]
      double A = (1.+b2)*(4.*dL1+2.*dL2-3.*L*log(2./(1.+b))-2.*L*log(b))-3.*b*log(4./(1.-b2))-4.*b*log(b);
            
      double DH = 0;
      if (h<3) {
      	DH = A/b + 1./(16.*pow(b,3))*(3.+34.*pow(b,2)-13.*pow(b,4))*L+3./(8.*pow(b,2))*(7.*pow(b,2)-1);      
      } else if (h==3) {
      	DH = A/b + 1./(16.*b)*(19.+2*b2+3*pow(b,4))*L+3./8.*(7.-b2);
      }
      
      K = 1.+4./3.*as/M_PI*DH;
    } else if (h==4) {
      K = 1.;
    }
    
    G = G*K;
  }

  return G;
}


double DecayTable::htt_offshell(double M, double m1, double m2, int Nc, int h) {

  complex <double> cs, cp;
  complex <double> cs_mtp, cp_mtp;
  model.get_coupling_huu(h,3, 3,cs,cp);

  double mtp = sm.get_qmass_pole(6);
  double mt_mt = sm.get_qmass_MSbar(6);
  double mbp = sm.get_qmass_pole(5);
  double mtr = sm.run_qmass_MSbar(mt_mt,mt_mt,M,mtp,mbp);

  cs_mtp = cs*mtp/mtr;
  cp_mtp = cp*mtp/mtr;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
  double result, error;
    
  integration_params_tt ip;
  ip.M = M;
  ip.mt = m1;
  ip.mb = mbp;
  ip.mW = sm.get_MW();
  ip.gtop = sm.get_gamma_top();
  ip.h = h;
  gsl_function F;
  F.function = &htt_fcn1;
  F.params = &ip;
      
  double imin = 2.*m1/M;
  double imax = 2.;
    
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F,imin,imax,0,1e-3,1000,
			w, &result, &error); 
    
  if (status) {
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H%1d -> tt (offshell). Please check result.\n", h);
      if (EXIT_ON_GSL_ERROR) exit(-1);
  	}
  }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;


  double v2 = sm.get_v2();
  double G = Nc*pow(M,3)/(64.*pow(M_PI,3)*v2)*(pow(abs(cs_mtp),2)+pow(abs(cp_mtp),2))*result;

  return G;
}


double DecayTable::htb_offshell(double M, double m1, double m1run, double m2, double m2run, complex<double> cs, complex<double> cp, int Nc) {

  complex <double> I(0.0,1.0);

  double md = m1;
  double mb  = sm.get_dmass_pole(3);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		
  double result, error;
  
  integration_params_tb ip;
  ip.M = M;
  ip.mt = m2;
  ip.md = md;
  ip.mb = mb;
  ip.mW = sm.get_MW();
  ip.gtop = sm.get_gamma_top();  
  
  ip.Z_u = real(I*(-cs+cp))*sm.get_v()/(sqrt(2)*m2run);
  ip.Z_d = 0;
  if (m1run>0){
     ip.Z_d = real(I*( cs+cp))*sm.get_v()/(sqrt(2)*m1run);
  }
  gsl_function F;
  F.function = &htb_fcn1;
  F.params = &ip;
	
  double imin = pow((md+mb),2);
  double imax = pow((M - ip.mW),2);
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F,imin,imax,0.,1e-3,1000,
			w, &result, &error);
	
  if (status) {
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H+ -> tb (offshell). Please check the result.\n");
      if (EXIT_ON_GSL_ERROR) exit(-1);
	}
  }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;
  
  double GF = sm.get_GF();
  
  double KHtb = 3*pow(GF,2)/(16*pow(M*M_PI,3));

  double G = KHtb*result;
 
  return G;
}

double DecayTable::hgaga(int h) {

  complex <double> I(0.,1);

  double alpha 	= sm.get_alpha0();
  double v  		= sm.get_v();
  double v2 		= sm.get_v2();
  double mW 		= sm.get_MW();
  
  double M 		= model.get_hmass(h);
  double M2 		= M*M;  

  double mHp 		= model.get_hmass(4);

  double tau_W   = M2/(4.*mW*mW);
  double tau_Hp  = M2/(4.*mHp*mHp);

  int Nc = 3;

  complex <double> gS, gP;
  complex <double> S_sum(0.,0.), P_sum(0.,0.);

  // Down-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_hdd(h,i,i,gS,gP);
    double m = sm.get_dmass_MSbar(i);
    double mp = sm.get_dmass_pole(i);    

    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;

      double m_HD = sm.run_qmass_MSbar(m,Qinit,M/2.,sm.get_qmass_pole(6),sm.get_qmass_pole(5))*mp/sm.run_qmass_MSbar(m,Qinit,mp,sm.get_qmass_pole(6),sm.get_qmass_pole(5));      
      double tau = M2/(4.*m_HD*m_HD);
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      complex <double> Sd = 2.*Nc*pow(-1./3,2)*gS*v/mrun*F_sf(tau);
      complex <double> Pd = 2.*Nc*pow(-1./3,2)*gP*v/mrun*F_pf(tau);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }
  
  // Up-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_huu(h,i,i,gS,gP);
    double m = sm.get_umass_MSbar(i);
    double mp = sm.get_umass_pole(i);    

    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;
      double m_HD = sm.run_qmass_MSbar(m,Qinit,M/2.,sm.get_qmass_pole(6),sm.get_qmass_pole(5))*mp/sm.run_qmass_MSbar(m,Qinit,mp,sm.get_qmass_pole(6),sm.get_qmass_pole(5));      
      double tau = M2/(4.*m_HD*m_HD);
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      complex <double> Sd = 2.*Nc*pow(2./3,2)*gS*v/mrun*F_sf(tau);
      complex <double> Pd = 2.*Nc*pow(2./3,2)*gP*v/mrun*F_pf(tau);

// QCD corrections in heavy quark limit for CP-even Higgs
      if ((!sm.b_HD) && (qcd_on) && (i==3)) {
        double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
        Sd=Sd*(1.-as/M_PI);
      }
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Leptons
  for (int i=1;i<=3;i++) {
    model.get_coupling_hll(h,i,i,gS,gP);
    double m = sm.get_lmass_pole(i);
    
    if (m > 0) {
      double tau = M2/(4.*m*m);
      complex <double> Sd = 2.*pow(-1.,2)*gS*v/m*F_sf(tau);
      complex <double> Pd = 2.*pow(-1.,2)*gP*v/m*F_pf(tau);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Charged Higgs and W contribution to scalar operator
  complex <double> g_hww,g_hhchc;
  model.get_coupling_vvh(3,3,h,g_hww);
  model.get_coupling_hhh(h,4,4,g_hhchc);

  g_hww = g_hww*v/(2.*pow(mW,2));
  
  S_sum = S_sum + g_hww*F_1(tau_W) + g_hhchc/v*v2/(2.*pow(mHp,2))*F_0(tau_Hp);
    
  double G = pow(M,3)*pow(alpha,2)/(256.*pow(M_PI,3)*v2)*(pow(abs(S_sum),2)+pow(abs(P_sum),2));

  return G;
}



double DecayTable::hZga(int h) {

  complex <double> I(0.,1);

  double alpha 	= sm.get_alpha0();
  double GF     = sm.get_GF();
  double v  		= sm.get_v();
  double v2 		= sm.get_v2();
  double mW 		= sm.get_MW();
  double mW2 		= mW*mW;
  double mZ 		= sm.get_MZ();
  double mZ2 		= mZ*mZ;
  double stw            = sm.get_sintw();
  double ctw            = sm.get_costw();
  
  double M 		= model.get_hmass(h);
  double M2 		= M*M;  

  double mHp 		= model.get_hmass(4);

  double tau_W   = (4.*mW2)/M2;
  double tau_Hp  = (4.*mHp*mHp)/M2;

  double lambda_W   = (4.*mW2)/mZ2;
  double lambda_Hp  = (4.*mHp*mHp)/mZ2;

  int Nc = 3;

  if (M<mZ) return 0.;

  complex <double> gS, gP;
  complex <double> S_sum(0.,0.), P_sum(0.,0.);

  // Down-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_hdd(h,i,i,gS,gP);
    double m = sm.get_dmass_MSbar(i);
    double mp = sm.get_dmass_pole(i);    

    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      double tau = (4.*mp*mp)/M2;
      double lambda = (4.*mp*mp)/mZ2;
      double Qf = -1./3.;
      double I3f = -1./2.;
      complex <double> Sd =  2.*Nc*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gS*v/mrun*FF_s(tau,lambda);
      complex <double> Pd =  2.*Nc*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gP*v/mrun*FF_p(tau,lambda);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Up-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_huu(h,i,i,gS,gP);
    double m = sm.get_umass_MSbar(i);
    double mp = sm.get_umass_pole(i);        

    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      double tau = (4.*mp*mp)/M2;
      double lambda = (4.*mp*mp)/mZ2;
      double Qf = 2./3.;
      double I3f = 1./2.;
      complex <double> Sd =  2.*Nc*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gS*v/mrun*FF_s(tau,lambda);
      complex <double> Pd =  2.*Nc*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gP*v/mrun*FF_p(tau,lambda);

      // QCD corrections in heavy top quark limit
      if ((!sm.b_HD)&&(qcd_on)&&(i==3)) {
//        QCD corrections disabled (as in HDECAY)  
//        double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
//        Sd=Sd*(1.-as/M_PI);
      }
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Leptons
  for (int i=1;i<=3;i++) {
    model.get_coupling_hll(h,i,i,gS,gP);
    double m = sm.get_lmass_pole(i);
    
    if (m > 0) {
      double tau = (4.*m*m)/M2;
      double lambda = (4.*m*m)/mZ2;
      double Qf = -1.;
      double I3f = -1./2.;
      complex <double> Sd =  2.*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gS*v/m*FF_s(tau,lambda);
      complex <double> Pd =  2.*Qf*(I3f-2.*pow(stw,2)*Qf)/ctw*gP*v/m*FF_p(tau,lambda);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Charged Higgs and W contribution to scalar operator
  complex <double> g_hww,g_hhchc;
  model.get_coupling_vvh(3,3,h,g_hww);
  model.get_coupling_hhh(h,4,4,g_hhchc);

  g_hww = g_hww*v/(2.*pow(mW,2));

//  S_sum = S_sum - g_hww*FW(tau_W,lambda_W) - (2.*ctw*ctw-1.)*g_hhchc/v*v2/(2.*pow(mHp,2))*FHp(tau_Hp,lambda_Hp);

// Charged Higgs contribution above differs with a factor ctw compared to HDECAY. 
// The normalisation below gives same result but is not consistent with formulas 2.23 and 2.33 in Anatomy II
  S_sum = S_sum - g_hww*FW(tau_W,lambda_W) - (2.*ctw-1./ctw)*g_hhchc/v*v2/(2.*pow(mHp,2))*FHp(tau_Hp,lambda_Hp);

  double G = alpha*pow(GF,2)*mW2*pow(M,3)*pow((1.-mZ2/M2),3)/(64.*pow(M_PI,4))*(pow(abs(S_sum),2)+pow(abs(P_sum),2));

  return G;
}


double DecayTable::hgg(int h) {

  double v  = sm.get_v();
  double v2 = sm.get_v2();
  
  double M = model.get_hmass(h);
  double M2 = M*M;  

  double mt   = sm.get_qmass_MSbar(6);

  complex <double> S_sum(0.,0.);
  complex <double> P_sum(0.,0.);
  complex <double> gS,gP;

  // Down-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_hdd(h,i,i,gS,gP);
    double m = sm.get_dmass_MSbar(i);
    double mp = sm.get_dmass_pole(i);    

    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      double tau = M2/(4.*mp*mp);

      complex <double> Sd = gS*v/mrun*F_sf(tau);
      complex <double> Pd = gP*v/mrun*F_pf(tau);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  // Up-type quarks
  for (int i=1;i<=3;i++) {
    model.get_coupling_huu(h,i,i,gS,gP);
    double m = sm.get_umass_MSbar(i);
    
    if (m > 0) {
      double Qinit = m;
      if(sm.b_HD) Qinit = SM::Q_HD;
      double mrun = sm.run_qmass_MSbar(m,Qinit,M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));
      double mp = sm.get_umass_pole(i);    
      double tau = M2/(4.*mp*mp);

      complex <double> Sd = gS*v/mrun*F_sf(tau);
      complex <double> Pd = gP*v/mrun*F_pf(tau);
      S_sum = S_sum + Sd;
      P_sum = P_sum + Pd;
    }
  }

  int Nf = sm.get_Nactivef(M);
  double as = sm.run_alphas_MSbar(M,sm.get_qmass_pole(6),sm.get_qmass_pole(5));

  double KS = 1;
  double KP = 1;
  // NNLO QCD corrections in heavy top quark limit
  // hep-ph/9705240 (CP-even Higgses)
  // hep-ph/9807241 (CP-odd)
  if (qcd_on) {
    KS = 1.+as/M_PI*(95./4.-7./6.*Nf)+pow(as/M_PI,2)*(156.808-5.708*log(mt*mt/M2));
    KP = 1.+as/M_PI*(97./4.-7./6.*Nf)+pow(as/M_PI,2)*(171.544-5.*log(mt*mt/M2));
  }
  
  double G = pow(M,3)*pow(as,2)/(32.*pow(M_PI,3)*v2)*(KS*pow(abs(S_sum),2)+KP*pow(abs(P_sum),2));

  return G;
}


complex <double> DecayTable::F_sf(double t) {
  double ti = 1./t;
  complex <double> c;
  c = ti*(1.+(1.-ti)*ftau(t));
  return c;
}


complex <double> DecayTable::F_pf(double t) {
  double ti = 1./t;
  complex <double> c;
  c = ti*ftau(t);
  return c;
}


complex <double> DecayTable::F_0(double t) {
  double ti = 1./t;
  complex <double> c;
  c = ti*(-1.+ti*ftau(t));
// same sign as in Djouadi
  return c;
}


complex <double> DecayTable::F_1(double t) {
  double ti = 1./t;
  complex <double> c;
  c = 2.+3.*ti+3.*ti*(2.-ti)*ftau(t);
// opposite sign to Djouadi
 return c;
}


complex <double> DecayTable::ftau(double t) {

  complex <double> c;
  complex <double> I(0.,1.);

  if (t<=1.) {
    double x = asin(sqrt(t));
    c = x*x;
  }

  if (t>1.) {
    complex <double> x = log((sqrt(t)+sqrt(t-1))/(sqrt(t)-sqrt(t-1)))-I*M_PI;
    c = -1./4.*x*x;
  }

  return c;
}


complex <double> DecayTable::gtau(double t) {

  complex <double> c;
  complex <double> I(0.,1.);

  if (t<=1.) {
    double x = asin(sqrt(t));
    c = x*sqrt(1./t-1.);
  }

  if (t>1.) {
    complex <double> x = log((sqrt(t)+sqrt(t-1))/(sqrt(t)-sqrt(t-1)))-I*M_PI;
    c = 1./2.*x*sqrt(1.-1./t);
  }

  return c;
}



complex <double> DecayTable::I_2(double tau, double lambda) {
  complex <double> c;
  c = -tau*lambda/2./(tau-lambda)*(ftau(1./tau)-ftau(1./lambda));
  return c;
}

complex <double> DecayTable::I_1(double tau, double lambda) {
  complex <double> c;
  c = tau*lambda/2./(tau-lambda) + tau*tau*lambda*lambda/2./(tau-lambda)/(tau-lambda)*(ftau(1./tau)-ftau(1./lambda)) + 
                                   tau*tau*lambda/(tau-lambda)/(tau-lambda)*(gtau(1./tau)-gtau(1./lambda));
  return c;
}

complex <double> DecayTable::FF_s(double tau, double lambda) {
  complex <double> c;
  c = I_1(tau,lambda) - I_2(tau,lambda);
  return c;
}

complex <double> DecayTable::FF_p(double tau, double lambda) {
  complex <double> c;
  c = I_2(tau,lambda);
  return c;
}

complex <double> DecayTable::FW(double tau, double lambda) {
  complex <double> c;
  double ctw            = sm.get_costw();
  double stw            = sm.get_sintw();
  c = ctw*(4.*(3.-stw*stw/ctw/ctw)*I_2(tau,lambda) + ((1.+2./tau)*stw*stw/ctw/ctw-(5.+2/tau))*I_1(tau,lambda));
  return c;
}

complex <double> DecayTable::FHp(double tau, double lambda) {
  complex <double> c;
  c = I_1(tau,lambda);
  return c;
}

// Interpolation with power c
double DecayTable::interp(double R, double x, double y, double c) {

  double ival = pow(R,c)*x + (1.-pow(R,c))*y;
   return ival;
}

double DecayTable::DHp(double ui, double uj, double xi, double xj, double sqL) {
 
 double eps = 1.E-12;
 
 if (ui < eps) ui = eps;
 if (uj < eps) uj = eps;
 if (xi < eps) xi = eps;
 if (xj < eps) xj = eps;
 
 double D = 9./4. + (3.-2*ui+2.*uj)/4.*log(ui/uj) + ((3./2.-ui-uj)*pow(sqL,2)+5.*ui*uj)/(2.*sqL*(1.-ui-uj))*log(xi*xj)+BHp(ui,uj,xi,xj,sqL);

 return D;
 
}      

double DecayTable::DHm(double ui, double uj, double xi, double xj, double sqL) {
 
 double eps = 1.E-12;
 
 if (ui < eps) ui = eps;
 if (uj < eps) uj = eps;
 if (xi < eps) xi = eps;
 if (xj < eps) xj = eps;
 
 double D = 3. + (uj-ui)/2.*log(ui/uj)+(pow(sqL,2)+2.*(1.-ui-uj))/(2.*sqL)*log(xi*xj)+BHp(ui,uj,xi,xj,sqL);

 return D;
 
}     

double DecayTable::BHp(double ui, double uj, double xi, double xj, double sqL) {

	double B=(1.-ui-uj)/sqL*(4.*gsl_sf_dilog(xi*xj)-2.*gsl_sf_dilog(-xi)-2.*gsl_sf_dilog(-xj)+2.*log(xi*xj)*log(1.-xi*xj)-log(xi)*log(1.+xi)-log(xj)*log(1.+xj));
    B = B-4.*(log(1.-xi*xj)+xi*xj/(1.-xi*xj)*log(xi*xj));
    B = B+(sqL+ui-uj)/sqL*(log(1.+xi)-xi/(1.+xi)*log(xi))+(sqL+uj-ui)/sqL*(log(1.+xj)-xj/(1.+xj)*log(xj));

	return B;
}

