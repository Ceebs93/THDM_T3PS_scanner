#include "HiggsBounds.h"

#include <stdio.h>
int main() {
    initialize_HiggsBounds(3, 2, LandH);

    double Mh[3] = {100, 200, 300};
    double GammaTotal_hj[3] = {0.1, 0.2, 0.3};
    int CP_value[3] = {-1, 0, 1};
    HiggsBounds_neutral_input_properties(Mh, GammaTotal_hj, CP_value);

    const double ghjss_s[3] = {0.1, 1, 10};
    const double ghjss_p[3] = {0.2, 2, 20};
    const double ghjcc_s[3] = {0.3, 3, 30};
    const double ghjcc_p[3] = {0.4, 4, 40};
    const double ghjbb_s[3] = {0.5, 5, 50};
    const double ghjbb_p[3] = {0.6, 6, 60};
    const double ghjtt_s[3] = {0.7, 7, 70};
    const double ghjtt_p[3] = {0.8, 8, 80};
    const double ghjmumu_s[3] = {0.9, 9, 90};
    const double ghjmumu_p[3] = {1, 10, 100};
    const double ghjtautau_s[3] = {2, 20, 200};
    const double ghjtautau_p[3] = {3, 30, 300};
    const double ghjWW[3] = {4, 40, 400};
    const double ghjZZ[3] = {5, 50, 500};
    const double ghjZga[3] = {6, 60, 600};
    const double ghjgaga[3] = {7, 70, 700};
    const double ghjgg[3] = {8, 80, 800};
    double ghjhiZ[3][3]; // symmetric
    ghjhiZ[0][0] = 11;
    ghjhiZ[0][1] = 12;
    ghjhiZ[0][2] = 13;
    ghjhiZ[1][0] = 12;
    ghjhiZ[1][1] = 22;
    ghjhiZ[1][2] = 23;
    ghjhiZ[2][0] = 13;
    ghjhiZ[2][1] = 23;
    ghjhiZ[2][2] = 33;

    HiggsBounds_neutral_input_effC(ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, ghjbb_s,
                                   ghjbb_p, ghjtt_s, ghjtt_p, ghjmumu_s,
                                   ghjmumu_p, ghjtautau_s, ghjtautau_p, ghjWW,
                                   ghjZZ, ghjZga, ghjgaga, ghjgg, ghjhiZ[0]);

    const double BR_hjss[3] = {0.1, 1, 10};
    const double BR_hjcc[3] = {0.2, 2, 20};
    const double BR_hjbb[3] = {0.3, 3, 30};
    const double BR_hjtt[3] = {0.4, 4, 40};
    const double BR_hjmumu[3] = {0.5, 5, 50};
    const double BR_hjtautau[3] = {0.6, 6, 60};
    const double BR_hjWW[3] = {0.7, 7, 70};
    const double BR_hjZZ[3] = {0.8, 8, 80};
    const double BR_hjZga[3] = {0.9, 9, 90};
    const double BR_hjgaga[3] = {1, 10, 100};
    const double BR_hjgg[3] = {2, 20, 200};

    HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, BR_hjtt,
                                   BR_hjmumu, BR_hjtautau, BR_hjWW, BR_hjZZ,
                                   BR_hjZga, BR_hjgaga, BR_hjgg);

    const double BR_hjinvisible[3] = {0.1, 1, 10};
    const double BR_hjemu[3] = {0.2, 2, 20};
    const double BR_hjetau[3] = {0.3, 3, 30};
    const double BR_hjmutau[3] = {0.4, 4, 40};
    double BR_hkhjhi[3][3][3];
    BR_hkhjhi[0][0][0] = 111;
    BR_hkhjhi[1][0][0] = 211;
    BR_hkhjhi[2][0][0] = 311;
    BR_hkhjhi[0][1][0] = 121;
    BR_hkhjhi[1][1][0] = 221;
    BR_hkhjhi[2][1][0] = 321;
    BR_hkhjhi[0][2][0] = 131;
    BR_hkhjhi[1][2][0] = 231;
    BR_hkhjhi[2][2][0] = 331;

    BR_hkhjhi[0][0][1] = 112;
    BR_hkhjhi[1][0][1] = 212;
    BR_hkhjhi[2][0][1] = 312;
    BR_hkhjhi[0][1][1] = 122;
    BR_hkhjhi[1][1][1] = 222;
    BR_hkhjhi[2][1][1] = 322;
    BR_hkhjhi[0][2][1] = 132;
    BR_hkhjhi[1][2][1] = 232;
    BR_hkhjhi[2][2][1] = 332;

    BR_hkhjhi[0][0][2] = 113;
    BR_hkhjhi[1][0][2] = 213;
    BR_hkhjhi[2][0][2] = 313;
    BR_hkhjhi[0][1][2] = 123;
    BR_hkhjhi[1][1][2] = 223;
    BR_hkhjhi[2][1][2] = 323;
    BR_hkhjhi[0][2][2] = 133;
    BR_hkhjhi[1][2][2] = 233;
    BR_hkhjhi[2][2][2] = 333;

    double BR_hjhiZ[3][3];
    BR_hjhiZ[0][0] = 11;
    BR_hjhiZ[0][1] = 12;
    BR_hjhiZ[0][2] = 13;
    BR_hjhiZ[1][0] = 21;
    BR_hjhiZ[1][1] = 22;
    BR_hjhiZ[1][2] = 23;
    BR_hjhiZ[2][0] = 31;
    BR_hjhiZ[2][1] = 32;
    BR_hjhiZ[2][2] = 33;

    double BR_hjHpiW[3][2];
    BR_hjHpiW[0][0] = 11;
    BR_hjHpiW[0][1] = 12;
    BR_hjHpiW[1][0] = 21;
    BR_hjHpiW[1][1] = 22;
    BR_hjHpiW[2][0] = 31;
    BR_hjHpiW[2][1] = 32;

    HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi[0][0],
                                      BR_hjhiZ[0], BR_hjemu, BR_hjetau,
                                      BR_hjmutau, BR_hjHpiW[0]);

    const double XS_ee_hjZ_ratio[3] = {0.1, 1, 10};
    const double XS_ee_bbhj_ratio[3] = {0.2, 2, 20};
    const double XS_ee_tautauhj_ratio[3] = {0.3, 3, 30};
    double XS_ee_hjhi_ratio[3][3]; // symmetric
    XS_ee_hjhi_ratio[0][0] = 11;
    XS_ee_hjhi_ratio[0][1] = 12;
    XS_ee_hjhi_ratio[0][2] = 13;
    XS_ee_hjhi_ratio[1][0] = 12;
    XS_ee_hjhi_ratio[1][1] = 22;
    XS_ee_hjhi_ratio[1][2] = 23;
    XS_ee_hjhi_ratio[2][0] = 13;
    XS_ee_hjhi_ratio[2][1] = 23;
    XS_ee_hjhi_ratio[2][2] = 33;
    HiggsBounds_neutral_input_LEP(XS_ee_hjZ_ratio, XS_ee_bbhj_ratio,
                                  XS_ee_tautauhj_ratio, XS_ee_hjhi_ratio[0]);

    const double CS_hj_ratio[3] = {0.1, 1, 10};
    const double CS_gg_hj_ratio[3] = {0.2, 2, 20};
    const double CS_bb_hj_ratio[3] = {0.3, 3, 30};
    const double CS_hjW_ratio[3] = {0.4, 4, 40};
    const double CS_hjZ_ratio[3] = {0.5, 5, 50};
    const double CS_vbf_ratio[3] = {0.6, 6, 60};
    const double CS_tthj_ratio[3] = {0.7, 7, 70};
    const double CS_thj_tchan_ratio[3] = {0.8, 8, 80};
    const double CS_thj_schan_ratio[3] = {0.9, 9, 90};
    const double CS_qq_hjZ_ratio[3] = {1.0, 10, 100};
    const double CS_gg_hjZ_ratio[3] = {1.1, 11, 110};
    const double CS_tWhj_ratio[3] = {1.2, 12, 120};
    double CS_hjhi[3][3]; // symmetric
    CS_hjhi[0][0] = 11;
    CS_hjhi[0][1] = 12;
    CS_hjhi[0][2] = 13;
    CS_hjhi[1][0] = 12;
    CS_hjhi[1][1] = 22;
    CS_hjhi[1][2] = 23;
    CS_hjhi[2][0] = 13;
    CS_hjhi[2][1] = 23;
    CS_hjhi[2][2] = 33;
    HiggsBounds_neutral_input_hadr(
        LHC7, CS_hj_ratio, CS_gg_hj_ratio, CS_bb_hj_ratio, CS_hjW_ratio,
        CS_hjZ_ratio, CS_vbf_ratio, CS_tthj_ratio, CS_thj_tchan_ratio,
        CS_thj_schan_ratio, CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, CS_tWhj_ratio,
        CS_hjhi[0]);

    const double Mhplus[2] = {10, 100};
    const double GammaTotal_Hpj[2] = {0.2, 2};
    const double CS_ee_HpjHmj_ratio[2] = {0.3, 3};
    const double BR_tWpb = 0.4;
    const double BR_tHpjb[2] = {0.5, 5};
    const double BR_Hpjcs[2] = {0.6, 6};
    const double BR_Hpjcb[2] = {0.7, 7};
    const double BR_Hpjtaunu[2] = {0.8, 8};
    const double BR_Hpjtb[2] = {0.9, 9};
    const double BR_HpjWZ[2] = {1, 10};
    double BR_HpjhiW[2][3];
    BR_HpjhiW[0][0] = 11;
    BR_HpjhiW[0][1] = 12;
    BR_HpjhiW[0][2] = 13;
    BR_HpjhiW[1][0] = 21;
    BR_HpjhiW[1][1] = 22;
    BR_HpjhiW[1][2] = 23;
    HiggsBounds_charged_input(Mhplus, GammaTotal_Hpj, CS_ee_HpjHmj_ratio,
                              BR_tWpb, BR_tHpjb, BR_Hpjcs, BR_Hpjcb,
                              BR_Hpjtaunu, BR_Hpjtb, BR_HpjWZ, BR_HpjhiW[0]);

    const double CS_Hpmjtb[2] = {0.1, 1};
    const double CS_Hpmjcb[2] = {0.2, 2};

    const double CS_Hpmjbjet[2] = {0.3, 3};
    const double CS_Hpmjcjet[2] = {0.4, 4};
    const double CS_Hpmjjetjet[2] = {0.5, 5};
    const double CS_HpmjW[2] = {0.6, 6};
    const double CS_HpmjZ[2] = {0.7, 7};
    const double CS_vbf_Hpmj[2] = {0.8, 8};
    const double CS_HpjHmj[2] = {0.9, 9};
    double CS_Hpmjhi[2][3];
    CS_Hpmjhi[0][0] = 11;
    CS_Hpmjhi[0][1] = 12;
    CS_Hpmjhi[0][2] = 13;
    CS_Hpmjhi[1][0] = 21;
    CS_Hpmjhi[1][1] = 22;
    CS_Hpmjhi[1][2] = 23;

    HiggsBounds_charged_input_hadr(
        LHC13, CS_Hpmjtb, CS_Hpmjcb, CS_Hpmjbjet, CS_Hpmjcjet, CS_Hpmjjetjet,
        CS_HpmjW, CS_HpmjZ, CS_vbf_Hpmj, CS_HpjHmj, CS_Hpmjhi[0]);

    double singleH;
    double ggH;
    double bbH;
    double VBF;
    double WH;
    double ZH;
    double ttH;
    double tH_tchan;
    double tH_schan;
    double qqZH;
    double ggZH;

    HiggsBounds_get_neutral_hadr_CS(2, LHC7, &singleH, &ggH, &bbH, &VBF, &WH,
                                    &ZH, &ttH, &tH_tchan, &tH_schan, &qqZH,
                                    &ggZH);
    // printf("singleH %f\n",singleH);
    // printf("ggH %f\n",ggH);
    // printf("bbH %f\n",bbH);
    // printf("VBF %f\n",VBF);
    // printf("WH %f\n",WH);
    // printf("ZH %f\n",ZH);
    // printf("ttH %f\n",ttH);
    // printf("tH_tchan %f\n",tH_tchan);
    // printf("tH_schan %f\n",tH_schan);

    double out_BR_hjss;
    double out_BR_hjcc;
    double out_BR_hjbb;
    double out_BR_hjtt;
    double out_BR_hjmumu;
    double out_BR_hjtautau;
    double out_BR_hjWW;
    double out_BR_hjZZ;
    double out_BR_hjZga;
    double out_BR_hjgaga;
    double out_BR_hjgg;
    HiggsBounds_get_neutral_BR(3, &out_BR_hjss, &out_BR_hjcc, &out_BR_hjbb,
                               &out_BR_hjtt, &out_BR_hjmumu, &out_BR_hjtautau,
                               &out_BR_hjWW, &out_BR_hjZZ, &out_BR_hjZga,
                               &out_BR_hjgaga, &out_BR_hjgg);
    // printf("BR_hjss %f\n",out_BR_hjss);
    // printf("BR_hjcc %f\n",out_BR_hjcc);
    // printf("BR_hjbb %f\n",out_BR_hjbb);
    // printf("BR_hjtt %f\n",out_BR_hjtt);
    // printf("BR_hjmumu %f\n",out_BR_hjmumu);
    // printf("BR_hjtautau %f\n",out_BR_hjtautau);
    // printf("BR_hjWW %f\n",out_BR_hjWW);
    // printf("BR_hjZZ %f\n",out_BR_hjZZ);
    // printf("BR_hjZga %f\n",out_BR_hjZga);
    // printf("BR_hjgaga %f\n",out_BR_hjgaga);
    // printf("BR_hjgg %f\n",out_BR_hjgg);

    const double dMhneut[3] = {1, 2, 3};
    const double dMhch[2] = {0.1, 0.2};
    HiggsBounds_set_mass_uncertainties(dMhneut, dMhch);

    int HBresult;
    int chan;
    double obsratio;
    int ncombined;
    run_HiggsBounds(&HBresult, &chan, &obsratio, &ncombined);
    // printf("HBresult %d\n", HBresult);
    // printf("chan %d\n", chan);
    // printf("obsratio %f\n", obsratio);
    // printf("ncombined %d\n", ncombined);

    run_HiggsBounds_single(5, &HBresult, &chan, &obsratio, &ncombined);
    // printf("HBresult %d\n", HBresult);
    // printf("chan %d\n", chan);
    // printf("obsratio %f\n", obsratio);
    // printf("ncombined %d\n", ncombined);

    int HBresult_full[6];
    int chan_full[6];
    double obsratio_full[6];
    int ncombined_full[6];
    run_HiggsBounds_full(HBresult_full, chan_full, obsratio_full,
                         ncombined_full);
    //   printf("HBresult %d %d %d %d %d %d\n", HBresult_full[0],
    //   HBresult_full[1],
    //          HBresult_full[2], HBresult_full[3], HBresult_full[4],
    //          HBresult_full[5]);
    //   printf("chan %d %d %d %d %d %d\n", chan_full[0], chan_full[1],
    //   chan_full[2],
    //          chan_full[3], chan_full[4], chan_full[5]);
    //   printf("obsratio %f %f %f %f %f %f\n", obsratio_full[0],
    //   obsratio_full[1],
    //          obsratio_full[2], obsratio_full[3], obsratio_full[4],
    //          obsratio_full[5]);
    //   printf("ncombined %d %d %d %d %d %d\n", ncombined_full[0],
    //   ncombined_full[1],
    //          ncombined_full[2], ncombined_full[3], ncombined_full[4],
    //          ncombined_full[5]);

    double predratio;
    HiggsBounds_get_most_sensitive_channels_per_Higgs(
        5, 3, &HBresult, &chan, &obsratio, &predratio, &ncombined);

    //   printf("HBresult %d\n", HBresult);
    //   printf("chan %d\n", chan);
    //   printf("obsratio %f\n", obsratio);
    //   printf("predratio %f\n", predratio);
    //   printf("ncombined %d\n", ncombined);

    HiggsBounds_get_most_sensitive_channels(3, &HBresult, &chan, &obsratio,
                                            &predratio, &ncombined);
    //       printf("HBresult %d\n", HBresult);
    //       printf("chan %d\n", chan);
    //       printf("obsratio %f\n", obsratio);
    //       printf("predratio %f\n", predratio);
    //       printf("ncombined %d\n", ncombined);

    // printf("SMBR_Htautau %f\n", SMBR_Htautau(200));
    // printf("SMGamma_H %f\n", SMGamma_H(100));
    // printf("SMGamma_tWpb %f\n", SMGamma_tWpb(180));

    // printf("SMCS_tev_bg_Hb %f\n", SMCS_tev_bg_Hb(70));
    // printf("SMCS_lhc8_gg_H %f\n", SMCS_lhc8_gg_H(200));

    // printf("SMCS_effC_bb_HZ %f\n", SMCS_effC_bb_HZ(140, LHC13, 0.9, -1.1));

    finish_HiggsBounds();
    return 0;
}
