#if !defined(HBHS_H)
#define HBHS_H

#include "SM.h"
#include "THDM.h"
#include "DecayTable.h"

using namespace std;

// --------------- HiggsBounds -----------------
void HB_init();

void HB_set_input(THDM model);
void HB_set_input_effC(THDM model);

void HB_run_full(int hbres[],int hbchan[], double hbobs[], int hbcomb[]);

double HB_get_gammah(double m);

void HB_finish();


// --------------- HiggsSignals ----------------
void HS_init();

void HS_set_pdf(int pdf);
void HS_set_nparam(int nparam);
void HS_set_assignment_range(double range);
void HS_setup_assignment_range_massobservables(double range);
void HS_set_mass_uncertainties(double dM[]);
void HS_set_rate_uncertainties(double dCS[], double dBR[]);
void HS_set_output_level(int level);

void HS_run(double *csqmu, double *csqmh, double *csqtot, int *nobs, double *pval);

void HS_get_Rvalues(int i, int collider, double *R_H_WW, double *R_H_ZZ, double *R_H_gaga, double *R_H_tautau, double *R_H_bb, double *R_VH_bb);

void HS_finish();



// Structs used by the Wrapper to HiggsBounds Fortran subroutines and common blocks
extern "C"
{

	extern void initialize_higgsbounds_int_(int *nH0, int *nHp, int *flag);
    
//     extern void higgsbounds_neutral_input_effc_(
//     double Mh[3],
//     double GammaTotal[3], 
//     double g2hjss_s[3],
//     double g2hjss_p[3],
//     double g2hjcc_s[3],
//     double g2hjcc_p[3],
//     double g2hjbb_s[3],
//     double g2hjbb_p[3], 
//     double g2hjtt_s[3],
//     double g2hjtt_p[3],
//     double g2hjmumu_s[3],
//     double g2hjmumu_p[3],
//     double g2hjtautau_s[3],
//     double g2hjtautau_p[3], 
//     double g2hjWW[3],
//     double g2hjZZ[3],
//     double g2hjZga[3],
//     double g2hjgaga[3],
//     double g2hjgg[3],
//     double g2hjggZ[3], 
//     double g2hjhiZ[3][3], 
//     double BR_hjinvisible[3],
//     double BR_hjhihi[3][3]);

    extern void higgsbounds_neutral_input_properties_(
    double Mh[3],
    double GammaTotal[3],
    int CP_value[3]);

    extern void higgsbounds_neutral_input_effc_(
    double ghjss_s[3],
    double ghjss_p[3],
    double ghjcc_s[3],
    double ghjcc_p[3],
    double ghjbb_s[3],
    double ghjbb_p[3],
    double ghjtt_s[3],
    double ghjtt_p[3],
    double ghjmumu_s[3],
    double ghjmumu_p[3],
    double ghjtautau_s[3],
    double ghjtautau_p[3],
    double ghjWW[3],
    double ghjZZ[3],
    double ghjZga[3],
    double ghjgaga[3],
    double ghjgg[3],
    double ghjhiZ[3][3]);

	extern void higgsbounds_neutral_input_smbr_(
	double BR_hjss[3],
	double BR_hjcc[3],
	double BR_hjbb[3],
	double BR_hjtt[3],
    double BR_hjmumu[3],
    double BR_hjtautau[3],
    double BR_hjWW[3],
    double BR_hjZZ[3],
    double BR_hjZga[3],
    double BR_hjgaga[3],
    double BR_hjgg[3]);

    extern void higgsbounds_neutral_input_nonsmbr_(
    double BR_hjinvisible[3],
    double BR_hkhjhi[3][3][3],
    double BR_hjhiZ[3][3],
    double BR_hjemu[3],
    double BR_hjetau[3],
    double BR_hjmutau[3],
    double BR_hjHpW[1][3]);

    
//    extern void higgsbounds_neutral_input_part_(
// 		double Mh[3], 
// 		double MhGammaTot[3],
// 		int 	 CP[3],
// 		double CS_lep_hjZ_ratio[3], 
// 		double CS_lep_bbhj_ratio[3],
// 		double CS_lep_tautauhj_ratio[3],
// 		double CS_lep_hjhi_ratio[3][3], 
// 		double CS_tev_gg_hj_ratio[3],
// 		double CS_tev_bb_hj_ratio[3], 
// 		double CS_tev_bg_hjb_ratio[3], 
// 		double CS_tev_ud_hjWp_ratio[3],
// 		double CS_tev_cs_hjWp_ratio[3], 
// 		double CS_tev_ud_hjWm_ratio[3],
// 		double CS_tev_cs_hjWm_ratio[3], 
// 		double CS_tev_gg_hjZ_ratio[3],
// 		double CS_tev_dd_hjZ_ratio[3],
// 		double CS_tev_uu_hjZ_ratio[3], 
// 		double CS_tev_ss_hjZ_ratio[3],
// 		double CS_tev_cc_hjZ_ratio[3], 
// 		double CS_tev_bb_hjZ_ratio[3], 
// 		double CS_tev_pp_vbf_ratio[3], 
// 		double CS_tev_pp_tthj_ratio[3],
// 		double CS_lhc7_pp_vbf_ratio[3], 
// 		double CS_lhc7_pp_tthj_ratio[3],
// 		double CS_lhc8_pp_vbf_ratio[3], 
// 		double CS_lhc8_pp_tthj_ratio[3],
// 		double BR_hjss[3],
// 		double BR_hjcc[3],
// 		double BR_hjbb[3],
// 		double BR_hjtautau[3],
// 		double BR_hjmumu[3],
// 		double BR_hjWW[3],
// 		double BR_hjZZ[3],
// 		double BR_hjZga[3],
// 		double BR_hjgaga[3],
// 		double BR_hjgg[3],
// 		double BR_hjinvisible[3],
// 		double BR_hjhihi[3][3]);

	extern void higgsbounds_charged_input_(
		double MHplus[1],
		double MHplusGammaTot[1],
		double CS_lep_HpjHmi_ratio[1],
		double BR_tWpb[1],
		double BR_tHpjb[1],
		double BR_Hpjcs[1],
		double BR_Hpjcb[1],
		double BR_Hptaunu[1],
		double BR_Hptb[1],
		double BR_HpWZ[1],
		double BR_HphiW[3][1]);

	extern void higgsbounds_charged_input_hadr_(
		int *collider,
		double CS_Hpjtb[1],
		double CS_Hpjcb[1],
		double CS_Hpjbjet[1],
		double CS_Hpjcjet[1],
		double CS_Hpjjetjet[1],
		double CS_HpjW[1],
		double CS_HpjZ[1],
		double CS_vbf_Hpj[1],
		double CS_HpjHmj[1],
		double CS_Hpjhi[3][1]);
		 
	extern void run_higgsbounds_(int *HBresult, int *chan, 
				    double *obsratio, int *ncombined);

	extern void run_higgsbounds_full_(int HBresult[6], int chan[6], 
				    double obsratio[6], int ncombined[6]);

	extern double smgamma_h_(double *Mh);
	extern double smbr_hww_(double *Mh);
	extern double smbr_hzz_(double *Mh);
	extern double smbr_hgg_(double *Mh);
	extern double smbr_htoptop_(double *Mh);
	extern double smbr_hbb_(double *Mh);
	extern double smbr_hcc_(double *Mh);
	extern double smbr_hss_(double *Mh);
	extern double smbr_htautau_(double *Mh);
	extern double smbr_hmumu_(double *Mh);
	extern double smbr_hzgam_(double *Mh);
	extern double smbr_hgamgam_(double *Mh);

	extern void finish_higgsbounds_();        

    extern void initialize_higgssignals_latestresults_(int *nHzero, int *nHplus);
    
    extern void initialize_higgssignals_lhc13_(int *nHzero, int *nHplus);

    extern void setup_pdf_(int *pdf);
    
    extern void setup_output_level_(int *level);
    
    extern void setup_nparam_(int *npara);

    extern void higgssignals_neutral_input_massuncertainty_(double dMh[3]);
    extern void setup_rate_uncertainties_(double dCS[5], double dBR[5]);
     
    extern void setup_assignmentrange_(double *range);
	extern void setup_assignmentrange_massobservables_(double *range);
	
    extern void run_higgssignals_(int *mode, double *csqmu, double *csqmh, double *csqtot, int *nobs, double *pval);
    
    extern void run_higgssignals_lhc_run1_combination_(double *csqmu, double *csqmh, double *csqtot, int *nobs, double *pval);
	
	extern void finish_higgssignals_();        

    extern void get_rvalues_(int *i, int *collider, double *R_H_WW, double *R_H_ZZ, double *R_H_gaga, double *R_H_tautau, double *R_H_bb, double *R_VH_bb);
	
}

#endif
