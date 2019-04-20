#include <cmath>

void get_chi2_ST(double &S, double &T, double &chi2_ST_hepfit, double &chi2_ST_gfitter)
{
	 double mu_S, mu_T, sig_S, sig_T, rho_ST;
	 double chi2_S, chi2_corr, chi2_T;

    // - HEPfit (arXiv:1608.01509), U=0 fixed
    mu_S   = 0.10;
    mu_T   = 0.12;
    sig_S  = 0.08;
    sig_T  = 0.07;
    rho_ST = 0.86;

	 chi2_S    = ( pow( (S - mu_S), 2.0) )/( ( pow(sig_S, 2.0)  )*(1 - pow(rho_ST, 2.0)) );
    chi2_corr = - 2.0*rho_ST*(S - mu_S)*( T - mu_T)/(sig_S*sig_T*(1 - pow(rho_ST,2.0)));
	 chi2_T    = ( pow( (T - mu_T), 2.0) )/( ( pow(sig_T, 2.0)  )*(1 - pow(rho_ST, 2.0)) );

    chi2_ST_hepfit = chi2_S + chi2_corr + chi2_T;

    // - Gfitter (arXiv:1407.3792), U=0 fixed
    mu_S   = 0.06;
    mu_T   = 0.10;
    sig_S  = 0.09;
    sig_T  = 0.07;
    rho_ST = 0.91;

	 chi2_S    = ( pow( (S - mu_S), 2.0) )/( ( pow(sig_S, 2.0)  )*(1 - pow(rho_ST, 2.0)) );
    chi2_corr = - 2.0*rho_ST*(S - mu_S)*( T - mu_T)/(sig_S*sig_T*(1 - pow(rho_ST,2.0)));
	 chi2_T    = ( pow( (T - mu_T), 2.0) )/( ( pow(sig_T, 2.0)  )*(1 - pow(rho_ST, 2.0)) );

    chi2_ST_gfitter = chi2_S + chi2_corr + chi2_T;
}
