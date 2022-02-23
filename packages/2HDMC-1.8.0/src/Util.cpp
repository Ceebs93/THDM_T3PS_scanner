#include "Util.h"
#include <gsl/gsl_min.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


int sign(double x) {
	if (x>=0) return 1;
	if (x<0) return -1;
	
	return 0;
}

double stability_fcn(double x, void *params) {
  stability_params p=*(stability_params*)params;
  return (1./2.)*p.lambda[1]*pow(cos(x),4.)+(1./2.)*p.lambda[2]*pow(sin(x),4.)+(p.lambda[3]+p.lambda[4]+p.lambda[5])*cos(x)*cos(x)*sin(x)*sin(x)+(2*(p.lambda[6]*cos(x)*cos(x)+p.lambda[7]*sin(x)*sin(x)))*sin(x)*cos(x);
}

bool stability_minimum(stability_params &p) {
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = p.m;
  double a = p.l, b = p.u;
  gsl_function F;
  double abs_err=0.00001;

  F.function = &stability_fcn;
  F.params = &p;
  
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

  status=gsl_min_fminimizer_set (s, &F, m, a, b);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);
      
      status = gsl_min_test_interval (a, b, abs_err, 0.0);
     
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  p.m=m;
  p.v=gsl_min_fminimizer_f_minimum(s);
  gsl_min_fminimizer_free (s);
  gsl_set_error_handler(old_handler);
  if (((m-p.l)<=abs_err)||((p.u-m)<=abs_err)) return false;
  return (status == GSL_SUCCESS);
}


void print_gsl_matrix(gsl_matrix *mat,int m,int n) {
  printf("\n");
  for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) {
      double x = gsl_matrix_get(mat,i,j);
      printf("%16.8g",x);
    }
    printf("\n");
  }

}


double Lint_s(double x, void *param) {

  double z = *(double*)param;
  double f = x*x*(2.-x)/(x*x+(1.-x)*z);
  return f;
}

double Lint_p(double x, void *param) {

  double z = *(double*)param;
  double f = -x*x*x/(x*x+(1.-x)*z);
  return f;
}

double Lint_c(double x, void *param) {

  double z = *(double*)param;
  double f = -x*(1.-x)/(x+z-1);
  return f;
}


double fint(double x, void *param) {

  double z = *(double*)param;
  double f = (1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z);
  return f;
}

double gint(double x, void *param) {

  double z = *(double*)param;
  double f = 1./(x*(1.-x)-z)*log(x*(1.-x)/z);
  return f;
}

double hvv_fcn(double y, void *params) {

  integration_params ip = *((integration_params*)params);
  double mH = ip.M;
  double mV = ip.m1;
  double gV = ip.gamma;

  double k = pow(mV/mH,2);
  double rgV = pow(gV/mH,2);

  double cg0 =  k*rgV;
  double cg2 = k;
  double cg4 = y;

  double A[3];
  //   A[0] = cg1/sqrt(cg);
  //   A[1] = cg + pow(cg1,2);
  //   A[2] = 1./sqrt(cg);
  //   
  A[0] = k/sqrt(cg0);
  A[1] = cg0+pow(k,2);
  A[2] = 1./sqrt(cg0);

  double cg = -(-log(A[1]) * sqrt(cg0) + log(A[1]) * cg4 * sqrt(cg0) + 0.2e1 * log(A[1]) * cg2 * sqrt(cg0) - 0.2e1 * atan(A[0]) * cg4 + 0.2e1 * atan(A[0]) * cg2 - 0.4e1 * atan(A[0]) * pow(cg2, 0.2e1) + 0.2e1 * atan(A[0]) * pow(cg4, 0.2e1) + log(0.2e1 * cg2 * pow(cg4, 0.2e1) + pow(cg4, 0.4e1) + cg0 - 0.2e1 * cg2 * cg4 + pow(cg2, 0.2e1) + pow(cg4, 0.2e1) + cg0 * pow(cg4, 0.2e1) - 0.2e1 * cg0 * cg4 - 0.2e1 * pow(cg4, 0.3e1)) * sqrt(cg0) - 0.2e1 * log(0.1e1 - cg4) * sqrt(cg0) - cg4 * log(0.2e1 * cg2 * pow(cg4, 0.2e1) + pow(cg4, 0.4e1) + cg0 - 0.2e1 * cg2 * cg4 + pow(cg2, 0.2e1) + pow(cg4, 0.2e1) + cg0 * pow(cg4, 0.2e1) - 0.2e1 * cg0 * cg4 - 0.2e1 * pow(cg4, 0.3e1)) * sqrt(cg0) + 0.2e1 * cg4 * log(0.1e1 - cg4) * sqrt(cg0) - 0.2e1 * cg2 * log(0.2e1 * cg2 * pow(cg4, 0.2e1) + pow(cg4, 0.4e1) + cg0 - 0.2e1 * cg2 * cg4 + pow(cg2, 0.2e1) + pow(cg4, 0.2e1) + cg0 * pow(cg4, 0.2e1) - 0.2e1 * cg0 * cg4 - 0.2e1 * pow(cg4, 0.3e1)) * sqrt(cg0) + 0.4e1 * cg2 * log(0.1e1 - cg4) * sqrt(cg0) - 0.2e1 * atan((cg2 + pow(cg4, 0.2e1) - cg4) / (-0.1e1 + cg4) * A[2]) * cg4 + 0.2e1 * atan((cg2 + pow(cg4, 0.2e1) - cg4) / (-0.1e1 + cg4) * A[2]) * cg2 - 0.4e1 * atan((cg2 + pow(cg4, 0.2e1) - cg4) / (-0.1e1 + cg4) * A[2]) * pow(cg2, 0.2e1) + 0.2e1 * atan((cg2 + pow(cg4, 0.2e1) - cg4) / (-0.1e1 + cg4) * A[2]) * pow(cg4, 0.2e1)) * A[2] / 0.2e1;

  return cg;

}


double hvv_fcn1(double y, void *params) {

  integration_params2 ip = *((integration_params2*)params);


  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
  double result, error;
    
  ip.x1 = y;
   
  gsl_function F;
  F.function = &hvv_fcn2;
  F.params = &ip;
      
  double imin = 0.;
  double imax = pow(1.-sqrt(y),2);
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags (&F,imin,imax,0,1e-6,1000,
			w, &result, &error); 
    
  if (status) 
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H -> VV (off-shell). Please check result.\n");
      if (EXIT_ON_GSL_ERROR) exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;

  return result;

}


double hvv_fcn2(double y, void *params) {

  integration_params2 ip = *((integration_params2*)params);

  double mH = ip.M;
  double x1 = ip.x1;
  double x2 = y;
  double mv = ip.m;
  double gammav = ip.gamma;
  double mh = mH;
  double mh2 = pow(mh,2);
  double q12 = x1*mh2;
  double q22 = x2*mh2;
  double mv2 = pow(mv,2);
  double lambda = 1+pow(x1,2)+pow(x2,2)-2*x1*x2-2*x1-2*x2;
  double gamma0 = sqrt(lambda)*(1+pow(x1,2)+pow(x2,2)+10*x1*x2-2*x1-2*x2);
  double bw1 = mv*gammav/(pow((q12-mv2),2)+pow(mv*gammav,2))/M_PI;
  double bw2 = mv*gammav/(pow((q22-mv2),2)+pow(mv*gammav,2))/M_PI;
  double prod;
  if(sqrt(x1)+sqrt(x2)>=1) {
     prod = 0;
  }
  else {
     prod = bw1*bw2*gamma0*mh2*mh2;
  }
  return prod;
     
}

double htt_fcn1(double y, void *params) {

  integration_params_tt ip = *((integration_params_tt*)params);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
  double result, error;
    
  ip.x1 = y;
   
  gsl_function F;
  F.function = &htt_fcn2;
  F.params = &ip;
      
  double imin = 0.;
  double imax = 2.-ip.x1;
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
//  int status = gsl_integration_qags (&F,imin,imax,0,1e-4,1000,
//			w, &result, &error); 
  int status = gsl_integration_qag (&F,imin,imax,0,1e-4,1000,GSL_INTEG_GAUSS61,
			w, &result, &error); 



    
  if (status) 
    if (status!=GSL_EROUND) {
      printf("GSL integration warning in H/A -> tt (offshell). Please check the result., %d\n", status);
      if (EXIT_ON_GSL_ERROR) exit(-1);
    }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free (w);
  old_handler = NULL;
  w = NULL;

  return result;

}


double htt_fcn2(double y, void *params) {

  integration_params_tt ip = *((integration_params_tt*)params);

  double mh = ip.M;
  double x1 = ip.x1;
  double x2 = y;
  double mt = ip.mt;
  double mW = ip.mW;
  double mb = ip.mb;
  double gammat = ip.gtop;

  double y1 = 1.-x1;
  double y2 = 1.-x2;
  
  double kt = pow(mt/mh,2);
  double kb = pow(mb/mh,2);
  double kW = pow(mW/mh,2);
  
  double gt = pow(gammat/mh,2);
  
  double prod;
  
  if ((x1 < sqrt(4.*kt))||(x2 < sqrt(4.*kb))) return 0.;
  if (x1+x2>2.) return 0.;
  
  double cond = fabs((2.*(1.-x1-x2+kt+kb-kW)+x1*x2)) - (sqrt(pow(x1,2)-4.*kt)*sqrt(pow(x2,2)-4.*kb));
   
  if(cond>0) {
     prod = 0;
  }
  else {

     double GA = 0;
     
     if (ip.h<3) {
		GA = pow(y1,2)*(1.-y1-y2+kW-5.*kt)+2.*kW*(y1*y2-kW-2.*kt*y1+4.*kt*kW)-kt*y1*y2+kt*(1.-4.*kt)*(2.*y1+kW+kt);
     } else if (ip.h==3) {
        GA = pow(y1,2)*(1.-y1-y2+kW-kt)+2*kW*(y1*y2-kW)-kt*(y1*y2-2.*y1-kW-kt);
     }
     
     prod = GA/(pow(y1,2)+gt*kt);
//     printf("xx %16.8E %16.8E %16.8E %16.8E\n", x1, x2, cond, prod);

  }
  
  return prod;
     
}


double hvh_fcn(double x, void *params) {

  integration_params ip = *((integration_params*)params);

  double M = ip.M;
  double m1 = ip.m1;
  double m2 = ip.m2;
  double gV = ip.gamma;

  double k1 = pow(m1/M,2);
  double k2 = pow(m2/M,2);
  double g  = k2*pow(gV/M,2);

  double A[3];
  A[0] = k2/sqrt(g);
  A[1] = g+pow(k2,2); 
  A[2] = 1./sqrt(g);

  double cg0 = k2;
  double cg2 = g;
  double cg4 = k1;
  double cg6 = x;
  

  // Maple code for differential width
  double cg = (-log(A[1]) * cg6 * sqrt(cg2) + log(A[1]) * sqrt(cg2) + 0.2e1 * atan(A[0]) * cg6 - 0.2e1 * atan(A[0]) * pow(cg6, 0.2e1) - 0.2e1 * atan(A[0]) * cg6 * cg4 + 0.2e1 * atan(A[0]) * cg6 * cg0 - 0.2e1 * atan(A[0]) * cg0 + cg6 * log(cg2 - 0.2e1 * pow(cg6, 0.3e1) + cg2 * pow(cg6, 0.2e1) + pow(cg0, 0.2e1) + pow(cg6, 0.2e1) - 0.2e1 * cg6 * cg0 + pow(cg6, 0.4e1) - 0.2e1 * pow(cg0, 0.2e1) * cg6 + 0.4e1 * cg0 * pow(cg6, 0.2e1) + 0.2e1 * cg0 * cg6 * cg4 - 0.2e1 * cg4 * cg0 * pow(cg6, 0.2e1) + pow(cg0, 0.2e1) * pow(cg6, 0.2e1) + 0.2e1 * pow(cg6, 0.3e1) * cg4 - 0.2e1 * pow(cg6, 0.3e1) * cg0 + pow(cg4, 0.2e1) * pow(cg6, 0.2e1) - 0.2e1 * cg2 * cg6 - 0.2e1 * cg4 * pow(cg6, 0.2e1)) * sqrt(cg2) - 0.2e1 * cg6 * log(0.1e1 - cg6) * sqrt(cg2) - log(cg2 - 0.2e1 * pow(cg6, 0.3e1) + cg2 * pow(cg6, 0.2e1) + pow(cg0, 0.2e1) + pow(cg6, 0.2e1) - 0.2e1 * cg6 * cg0 + pow(cg6, 0.4e1) - 0.2e1 * pow(cg0, 0.2e1) * cg6 + 0.4e1 * cg0 * pow(cg6, 0.2e1) + 0.2e1 * cg0 * cg6 * cg4 - 0.2e1 * cg4 * cg0 * pow(cg6, 0.2e1) + pow(cg0, 0.2e1) * pow(cg6, 0.2e1) + 0.2e1 * pow(cg6, 0.3e1) * cg4 - 0.2e1 * pow(cg6, 0.3e1) * cg0 + pow(cg4, 0.2e1) * pow(cg6, 0.2e1) - 0.2e1 * cg2 * cg6 - 0.2e1 * cg4 * pow(cg6, 0.2e1)) * sqrt(cg2) + 0.2e1 * log(0.1e1 - cg6) * sqrt(cg2) + 0.2e1 * atan((cg4 * cg6 - cg6 * cg0 + pow(cg6, 0.2e1) - cg6 + cg0) / (-0.1e1 + cg6) * A[2]) * cg6 - 0.2e1 * atan((cg4 * cg6 - cg6 * cg0 + pow(cg6, 0.2e1) - cg6 + cg0) / (-0.1e1 + cg6) * A[2]) * pow(cg6, 0.2e1) - 0.2e1 * atan((cg4 * cg6 - cg6 * cg0 + pow(cg6, 0.2e1) - cg6 + cg0) / (-0.1e1 + cg6) * A[2]) * cg6 * cg4 + 0.2e1 * atan((cg4 * cg6 - cg6 * cg0 + pow(cg6, 0.2e1) - cg6 + cg0) / (-0.1e1 + cg6) * A[2]) * cg6 * cg0 - 0.2e1 * atan((cg4 * cg6 - cg6 * cg0 + pow(cg6, 0.2e1) - cg6 + cg0) / (-0.1e1 + cg6) * A[2]) * cg0) * A[2] / 0.2e1;


  return cg;

}

double cubic(double v, double X[4], double Y[4]) {

  if (v<0) return 0;
  
  double x = log(v);
   
  double x0 = log(X[0]);
  double x1 = log(X[1]);
  double x2 = log(X[2]);
  double x3 = log(X[3]);
  double y0 = log(Y[0]);
  double y1 = log(Y[1]);
  double y2 = log(Y[2]);
  double y3 = log(Y[3]);
   
  double A0 = (x-x1)*(x-x2)*(x-x3)/(x0-x1)/(x0-x2)/(x0-x3);
  double A1 = (x-x0)*(x-x2)*(x-x3)/(x1-x0)/(x1-x2)/(x1-x3);
  double A2 = (x-x0)*(x-x1)*(x-x3)/(x2-x0)/(x2-x1)/(x2-x3);
  double A3 = (x-x0)*(x-x1)*(x-x2)/(x3-x0)/(x3-x1)/(x3-x2);
  
  double c = exp(A0*y0+A1*y1+A2*y2+A3*y3);

  return c;
}


double htb_fcn1(double y, void *params) {

	integration_params_tb ip = *((integration_params_tb*)params);
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	
	double result, error;
	double msq = ip.M*ip.M - ip.mW*ip.mW;
	double lsq1 = sqrtlambda( ip.M*ip.M , y , ip.mW*ip.mW );
	double lsq2 = sqrtlambda( y , ip.mb*ip.mb , ip.md*ip.md );
	ip.x1 = y;
	 
	gsl_function F;
	F.function = &htb_fcn2;
	F.params = &ip;
	
	double imin = (msq*msq-pow((lsq1 + lsq2),2))/(4*y);
	double imax = (msq*msq-pow((lsq1 - lsq2),2))/(4*y);
	
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
	int status = gsl_integration_qag (&F,imin,imax,0,1e-4,1000,GSL_INTEG_GAUSS61,w, &result, &error); 
	
	if(status)
	if(status!=GSL_EROUND) {
		printf("GSL integration warning in H+ -> tb (off-shell). Please check the result.\n");
		if (EXIT_ON_GSL_ERROR) exit(-1);
	}
	gsl_set_error_handler(old_handler);
	gsl_integration_workspace_free (w);
	old_handler = NULL;
	w = NULL;

	return result;

}


double htb_fcn2(double y, void *params) 
{
	integration_params_tb ip = *((integration_params_tb*)params);
	
	double mh = ip.M;
	double x1 = ip.x1;
	double x2 = y;
	double mt = ip.mt;
	double mW = ip.mW;
	double md = ip.md;
	double mb = ip.mb;
	double Zu = ip.Z_u;
	double Zd = ip.Z_d;
	double gtop = ip.gtop;
		
	double p1p2 = .5*(mh*mh + mb*mb - x1 - x2);
	double p1p3 = .5*(x1 - md*md - mb*mb);
	double p2p3 = .5*(x2 - mW*mW - md*md);
	
	double num = (mW*mW*p1p3 + 2*p2p3*p1p2)*(Zu*Zu*pow(mt,4) - Zd*Zd*md*md*x2) + 
                     (mW*mW*md*md*(p2p3 + mb*mb) + 2*md*md*p2p3*(mW*mW + p2p3))*(2*Zd*Zd*(p1p2+p1p3) + 2*mt*mt*Zu*Zd);
	
	double prod = num/(pow((x2 - mt*mt),2)+mt*mt*gtop*gtop);
	
	return prod;
}


double L(double x, double y, double z) {
  double lam = pow(1.-x/z-y/z,2)-4.*x*y/pow(z,2);
  return lam;
}

double sqrtlambda(double a1, double a2, double a3)
{
	double lam = sqrt(pow(a1,2) + pow(a2,2)+pow(a3,2) - 2*(a1*a2 + a1*a3 + a2*a3));
	return lam;
}


