#if !defined(UTIL_H)
#define UTIL_H

#include <gsl/gsl_matrix.h>

// Flag to set if the program should exit on GSL error
const static bool EXIT_ON_GSL_ERROR = false;



struct stability_params {
  double lambda[8];
  double l;
  double u;
  double m;
  double v;
};

struct integration_params {
  double M;
  double m1;
  double m2;
  double gamma;
};

struct integration_params2 {
  double M;
  double x1;
  double m;
  double gamma;
};

struct integration_params_tt {
  double M;
  double mW;
  double mb;
  double mt;
  double x1;
  double gtop;
  int h;
};

struct integration_params_tb {
  double M;
  double mW;
  double mb;
  double mt;
  double x1;
  double gtop;
  double md;
  double Z_u;
  double Z_d;
};

int sign(double x);

double stability_fcn(double x, void *params);

bool stability_minimum(stability_params &p);

void print_gsl_matrix(gsl_matrix *mat,int m,int n);

double fint(double x, void* par);

double gint(double x, void* par);

double Lint_s(double x, void *param);

double Lint_p(double x, void *param);

double Lint_c(double x, void *param);

double hvv_fcn(double y, void *params);

double hvv_fcn1(double y, void *params);

double hvv_fcn2(double y, void *params);

double htt_fcn1(double y, void *params);

double htt_fcn2(double y, void *params);

double htb_fcn1(double y, void *params);

double htb_fcn2(double y, void *params);

double cubic(double x, double X[4], double Y[4]);

double hvh_fcn(double x, void *params);

double  L(double x, double y, double z);
double  DHp(double ui, double uj, double xi, double xj, double sqL);
double  DHm(double ui, double uj, double xi, double xj, double sqL);
double  BHp(double ui, double uj, double xi, double xj, double sqL);
double sqrtlambda(double a1, double a2, double a3);



#endif