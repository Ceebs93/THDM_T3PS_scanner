#include "HBHS.h"
#include "THDM.h"
#include "DecayTable.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>


static bool HB_initialized = false;
static bool HS_initialized = false;

#if defined HiggsBounds


void HB_init() {
  int nH0=3;
  int nHp=1;
  int hbflag=3;	
  
  if (HB_initialized) return;
  
//  printf("\nInitializing HiggsBounds... ");
   
   // Third argument is HB analysis setting: 1='onlyL', 2='onlyH' 3='LandH'
    initialize_higgsbounds_int_(&nH0, &nHp, &hbflag);
//  initialize_higgsbounds_(&nH0,&nHp,whichexpt);

 // printf("Please cite the relevant publications when using Higgs mass limits.\n");

  HB_initialized = true;

}


void HB_set_input(THDM model) {

	HB_set_input_effC(model);

}


void HB_set_input_effC(THDM model) {

	if (!HB_initialized) {
		cout << "WARNING: HiggsBounds must be initialized with HB_init() before usage" << endl;
		return;
	}

    bool debug =false;

    double Mh[3];
    double GammaTotal[3]; 
    double ghjss_s[3];
    double ghjss_p[3];
    double ghjcc_s[3];
    double ghjcc_p[3];
    double ghjbb_s[3];
    double ghjbb_p[3]; 
    double ghjtt_s[3];
    double ghjtt_p[3];
    double ghjmumu_s[3];
    double ghjmumu_p[3];
    double ghjtautau_s[3];
    double ghjtautau_p[3]; 
    double ghjWW[3];
    double ghjZZ[3];
    double ghjZga[3];
    double ghjgaga[3];
    double ghjgg[3];
//     double ghjggZ[3]; 
    double ghjhiZ[3][3];

// non-SM decay branching ratios     
    double BR_hjinvisible[3];
    double BR_hkhjhi[3][3][3];
    double BR_hjhiZ[3][3];
	double BR_hjHpW[1][3];
	double BR_hjemu[3];
	double BR_hjetau[3];
	double BR_hjmutau[3];

// SM decay branching ratios    
	double BR_hjss[3];
	double BR_hjcc[3];
	double BR_hjbb[3];
	double BR_hjtt[3];
	double BR_hjmumu[3];
	double BR_hjtautau[3];
	double BR_hjWW[3];
	double BR_hjZZ[3];
	double BR_hjZga[3];
	double BR_hjgaga[3];
	double BR_hjgg[3];    
    
    double MHp[1];
	double MHplusGammaTot[1];
	double CS_lep_HpjHmi_ratio[1];
	double BR_tWpb[1];
	double BR_tHpjb[1];
	double BR_Hpjcs[1];
	double BR_Hpjcb[1];
	double BR_Hptaunu[1];
	double BR_Hptb[1];
	double BR_HpWZ[1];
	double BR_HphiW[3][1];
	double CS_Hpjtb[1];
	double CS_Hpjbjet[1];
	double CS_Hpjcb[1];
	double CS_Hpjcjet[1];
	double CS_Hpjjetjet[1];
	double CS_HpjW[1];
	double CS_HpjZ[1];
	double CS_vbf_Hpj[1];
	double CS_HpjHmj[1];
	double CS_Hpjhi[3][1];
    int CP_value[3];
    
    int type;
    char str[100];
    
    double sba, tanb, lam6, lam7, m122;
    model.get_param_phys(Mh[0], Mh[1], Mh[2], MHp[0], sba, lam6, lam7, m122, tanb); 
	type = model.get_yukawas_type();
    
#if defined debug
    printf("Masses for HB/HS: %16.8E %16.8E %16.8E\n", Mh[0], Mh[1], Mh[2]);
#endif

    THDM sm_like;
    DecayTable table(model), sm_table(sm_like);
 
    SM sm = model.get_SM();

  	double g=sm.get_g();
  	double costw=sm.get_costw();
  	double mt = sm.get_umass_pole(3);
 
    complex <double> c,cs,cp,c_sm,cs_sm,cp_sm;
    
    for (int h=1;h<=3;h++) {
      double mh = Mh[h-1];
      CP_value[h-1] = 0;
//      sm_like.set_param_phys(Mh[h-1], Mh[h-1]*10., Mh[h-1]*10., Mh[h-1]*10., 1.0, 0., 0., 0., 1.0);
      sm_like.set_param_sm(Mh[h-1]);
      sm_like.set_yukawas_type(1);
      
      table.set_model(model);
      sm_table.set_model(sm_like);

//  NEUTRAL HIGGS INPUT
// n.b.: In general the g's should also contain the sign of the SM normalized coupling!

      model.get_coupling_hdd(h,2,2,cs,cp);
	  sm_like.get_coupling_hdd(1,2,2,cs_sm,cp_sm);
      ghjss_s[h-1] = abs(cs/cs_sm);
      ghjss_p[h-1] = abs(cp/cs_sm);

	  if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "ss", ghjss_s[h-1], ghjss_p[h-1]);	  

      model.get_coupling_hdd(h,3,3,cs,cp);
	  sm_like.get_coupling_hdd(1,3,3,cs_sm,cp_sm);
      ghjbb_s[h-1] = abs(cs/cs_sm);
      ghjbb_p[h-1] = abs(cp/cs_sm);

      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "bb", ghjbb_s[h-1], ghjbb_p[h-1]);
      
      model.get_coupling_huu(h,2,2,cs,cp);
	  sm_like.get_coupling_huu(1,2,2,cs_sm,cp_sm);
      ghjcc_s[h-1] = abs(cs/cs_sm);
      ghjcc_p[h-1] = abs(cp/cs_sm);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "cc", ghjcc_s[h-1], ghjcc_p[h-1]);
      
      model.get_coupling_huu(h,3,3,cs,cp);
	  sm_like.get_coupling_huu(1,3,3,cs_sm,cp_sm);
      ghjtt_s[h-1] = abs(cs/cs_sm);
      ghjtt_p[h-1] = abs(cp/cs_sm);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "tt", ghjtt_s[h-1], ghjtt_p[h-1]);

      model.get_coupling_hll(h,2,2,cs,cp);
	  sm_like.get_coupling_hll(1,2,2,cs_sm,cp_sm);
      ghjmumu_s[h-1] = abs(cs/cs_sm);
      ghjmumu_p[h-1] = abs(cp/cs_sm);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "mumu", ghjmumu_s[h-1], ghjmumu_p[h-1]);
      
      model.get_coupling_hll(h,3,3,cs,cp);
	  sm_like.get_coupling_hll(1,3,3,cs_sm,cp_sm);
      ghjtautau_s[h-1] = abs(cs/cs_sm);
      ghjtautau_p[h-1] = abs(cp/cs_sm);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "tata", ghjtautau_s[h-1], ghjtautau_p[h-1]);
      
      model.get_coupling_vvh(2, 2, h, c);
      sm_like.get_coupling_vvh(2, 2, 1, c_sm);
      ghjZZ[h-1] = abs(c/c_sm);
      if (debug) printf("%2d %5s %16.8E\n", h, "ZZ", ghjZZ[h-1]);      

      model.get_coupling_vvh(3, 3, h, c);
      sm_like.get_coupling_vvh(3, 3, 1, c_sm);
      ghjWW[h-1] = abs(c/c_sm);
      if (debug) printf("%2d %5s %16.8E\n", h, "WW", ghjWW[h-1]);
          
      double hgaga = table.get_gamma_hgaga(h);
      double hgaga_sm = sm_table.get_gamma_hgaga(1);
      ghjgaga[h-1] = sqrt(hgaga/hgaga_sm);
      if (debug) printf("%2d %5s %16.8E\n", h, "gaga",ghjgaga[h-1]);

      double hZga = table.get_gamma_hZga(h);
      double hZga_sm = sm_table.get_gamma_hZga(1);
      ghjZga[h-1] = sqrt(hZga/hZga_sm);
      if (debug) printf("%2d %5s %16.8E\n", h, "Zga", ghjZga[h-1]);
      
      double hgg = table.get_gamma_hgg(h);
      double hgg_sm = sm_table.get_gamma_hgg(1);
      ghjgg[h-1] = sqrt(hgg/hgg_sm);
      if (debug) printf("%2d %5s %16.8E\n", h, "gg", ghjgg[h-1]);

//  Set branching ratios directly:
      BR_hjss[h-1] = table.get_gamma_hdd(h,2,2)/table.get_gammatot_h(h);
      BR_hjcc[h-1] = table.get_gamma_huu(h,2,2)/table.get_gammatot_h(h);
      BR_hjbb[h-1] = table.get_gamma_hdd(h,3,3)/table.get_gammatot_h(h);
      BR_hjtt[h-1] = table.get_gamma_huu(h,3,3)/table.get_gammatot_h(h);
      BR_hjmumu[h-1] = table.get_gamma_hll(h,2,2)/table.get_gammatot_h(h);
      BR_hjtautau[h-1] = table.get_gamma_hll(h,3,3)/table.get_gammatot_h(h);
      BR_hjWW[h-1] = table.get_gamma_hvv(h,3)/table.get_gammatot_h(h);
      BR_hjZZ[h-1] = table.get_gamma_hvv(h,2)/table.get_gammatot_h(h);
      BR_hjZga[h-1] = table.get_gamma_hZga(h)/table.get_gammatot_h(h);
      BR_hjgaga[h-1] = table.get_gamma_hgaga(h)/table.get_gammatot_h(h);
      BR_hjgg[h-1] = table.get_gamma_hgg(h)/table.get_gammatot_h(h);

      BR_hjemu[h-1]=table.get_gamma_hll(h,1,2)/table.get_gammatot_h(h);
      BR_hjetau[h-1]=table.get_gamma_hll(h,1,3)/table.get_gammatot_h(h);
      BR_hjmutau[h-1]=table.get_gamma_hll(h,2,3)/table.get_gammatot_h(h);
      BR_hjHpW[1][h-1]=table.get_gamma_hvh(h,3,4)/table.get_gammatot_h(h);
      BR_hjinvisible[h-1]=0.0;
     	  
	  GammaTotal[h-1] = table.get_gammatot_h(h);
      
      if (debug) printf("gtot %16.8E %16.8E %16.8E %16.8E\n",  GammaTotal[h-1], table.get_gammatot_h(h), HB_get_gammah(Mh[h-1]), sm_table.get_gammatot_h(1));
    }
     
  	 for (int j=1;j<=3;j++) {  	 
      for (int i=1;i<=3;i++) {
       for (int k=1;k<=3;k++) {
        BR_hkhjhi[i-1][j-1][k-1]=table.get_gamma_hhh(k,j,i)/GammaTotal[k-1];
        }
       BR_hjhiZ[i-1][j-1]=table.get_gamma_hvh(j,2,i)/GammaTotal[j-1];        
       model.get_coupling_vhh(2,j,i,c);
       ghjhiZ[i-1][j-1]=abs(c)/(g/2./costw);
       if (debug) printf("%2d %2d hj->hihi %16.8E\n", j, i, BR_hkhjhi[j-1][i-1][i-1]);
      }
     }

    higgsbounds_neutral_input_properties_(Mh,GammaTotal,CP_value);

    higgsbounds_neutral_input_effc_(ghjss_s,ghjss_p,ghjcc_s,ghjcc_p,
    ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p,
    ghjmumu_s,ghjmumu_p,
    ghjtautau_s,ghjtautau_p,
    ghjWW,ghjZZ,ghjZga,
    ghjgaga,ghjgg,ghjhiZ);

    higgsbounds_neutral_input_smbr_(BR_hjss,BR_hjcc,BR_hjbb,BR_hjtt,
    BR_hjmumu,BR_hjtautau,BR_hjWW,BR_hjZZ,BR_hjZga,BR_hjgaga,BR_hjgg);

    higgsbounds_neutral_input_nonsmbr_(BR_hjinvisible,
    BR_hkhjhi,BR_hjhiZ,BR_hjemu,BR_hjetau,BR_hjmutau,BR_hjHpW);

//  CHARGED HIGGS INPUT
    
    CS_lep_HpjHmi_ratio[0] = 1.;
	BR_tWpb[0] = sm.get_gamma_top()/table.get_gammatot_top();
	BR_tHpjb[0]=table.get_gamma_uhd(3,4,3)/table.get_gammatot_top();

	BR_Hpjcs[0] = table.get_gamma_hdu(4,2,2)/table.get_gammatot_h(4);
	BR_Hpjcb[0] = table.get_gamma_hdu(4,3,2)/table.get_gammatot_h(4);
	BR_Hptaunu[0] = table.get_gamma_hln(4,3,3)/table.get_gammatot_h(4);
	BR_Hptb[0] = table.get_gamma_hdu(4,3,3)/table.get_gammatot_h(4);
	BR_HpWZ[0] = 0.0;
    for (int i=1;i<=3;i++){
     BR_HphiW[i-1][0]= table.get_gamma_hvh(4,3,i)/table.get_gammatot_h(4);
    }

    higgsbounds_charged_input_(MHp,MHplusGammaTot,
    CS_lep_HpjHmi_ratio,BR_tWpb,BR_tHpjb,
    BR_Hpjcs,BR_Hpjcb,BR_Hptaunu,BR_Hptb,
    BR_HpWZ,BR_HphiW);

// IF CALCULATED, WE CAN SET THE 13 TEV LHC CROSS SECTIONS FOR DIRECT
// CHARGED HIGGS PRODUCTION HERE:
	CS_Hpjtb[0]=0;
	CS_Hpjcb[0]=0;
	CS_Hpjbjet[0]=0;
	CS_Hpjcjet[0]=0;
	CS_Hpjjetjet[0]=0;
	CS_HpjW[0]=0;
	CS_HpjZ[0]=0;
	CS_vbf_Hpj[0]=0;
	CS_HpjHmj[0]=0;
	CS_Hpjhi[0][0]=0;
	CS_Hpjhi[1][0]=0;	
	CS_Hpjhi[2][0]=0;	
	
	

	int collider = 13;
    higgsbounds_charged_input_hadr_(&collider, CS_Hpjtb, CS_Hpjcb, CS_Hpjbjet,
    CS_Hpjcjet, CS_Hpjjetjet, CS_HpjW, CS_HpjZ, CS_vbf_Hpj, CS_HpjHmj, CS_Hpjhi);

}


double HB_get_gammah(double m) {
 
   return smgamma_h_(&m);
}


void HB_run_full(int hbres[], int hbchan[], double hbobs[], int hbcomb[]) {

//    printf("Running HB full\n");
	if (!HB_initialized) {
		cout << "WARNING: HiggsBounds must be initialized with HB_init() before usage" << endl;
		return;
	}
	run_higgsbounds_full_(hbres,hbchan,hbobs,hbcomb);
}

void HS_init() {

 int nH0 = 3;
 int nHp = 1;
 double range = 3.;
 
// printf("\nInitializing HiggsSignals... ");

//  initialize_higgssignals_latestresults_(&nH0,&nHp);
 initialize_higgssignals_lhc13_(&nH0,&nHp);
 
 HS_initialized = true;
}

void HS_finish() {

  finish_higgssignals_();
}

void HB_finish() {

  finish_higgsbounds_();
}


void HS_run(double *csqmu, double *csqmh, double *csqtot, int *nobs, double *pval) {
	if (!HS_initialized) {
		cout << "WARNING: HiggsSignals must be initialized with HS_init() before usage" << endl;
		return;
	}

	int mode = 1;
	int nobs1;
	int nobs2;
	double csqmu1;
	double csqmu2;	
	double csqmh1;
	double csqmh2;
	double csqtot1;
	double csqtot2;			
//   run_higgssignals_(&mode, csqmu, csqmh, csqtot, nobs, pval);
  run_higgssignals_(&mode, &csqmu1, &csqmh1, &csqtot1, &nobs1, pval);
  run_higgssignals_lhc_run1_combination_(&csqmu2, &csqmh2, &csqtot2, &nobs2, pval);
  
//   cout << csqmu1 << " " << csqmh1 << " " << csqtot1 << endl;
//   cout << csqmu2 << " " << csqmh2 << " " << csqtot2 << endl;  
    
  *csqmu = csqmu1 + csqmu2;
  *csqmh = csqmh1 + csqmh2;  
  *csqtot = *csqmu + *csqmh;
  *nobs = nobs1 + nobs2;
//   cout << *csqmu << " " << *csqmh << " " << *csqtot << endl;
}


void HS_set_pdf(int pdf) {

  setup_pdf_(&pdf);

}

void HS_set_assignment_range(double range) {
  setup_assignmentrange_(&range);
}

void HS_setup_assignment_range_massobservables(double range) {
  setup_assignmentrange_massobservables_(&range);
}


void HS_set_rate_uncertainties(double dCS[], double dBR[]) {
  setup_rate_uncertainties_(dCS, dBR);

}

void HS_set_mass_uncertainties(double dMh[]) {
	higgssignals_neutral_input_massuncertainty_(dMh);
}

void HS_get_Rvalues(int i, int collider, double *R_H_WW, double *R_H_ZZ, double *R_H_gaga, double *R_H_tautau, double *R_H_bb, double *R_VH_bb) {
  get_rvalues_(&i, &collider, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb);

}

void HS_set_output_level(int level) {

  setup_output_level_(&level);

}

#endif

