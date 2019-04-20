c		integer error

c used by FHHiggsCorr
c		double precision MHiggs(4)
c		double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSelectUZ:
		integer uzint, uzext, mfeff

c used by FHCouplings:
#include "FHCouplings.h"
		double complex couplings(ncouplings), couplingsms(ncouplingsms)
		double precision gammas(ngammas), gammasms(ngammasms)
		integer fast

c used by FHHiggsProd:
		double precision sqrts_HB, prodxs(nprodxs)

c used by FHGetPara:
c		integer nmfv
c		double precision MSf(2,4,3),MASf(6,4), MCha(2), MNeu(4)
c		double complex USf(2,2,4,3),UASf(6,6,4)
c		double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
c		double complex DeltaMB
c		double precision MGl
c		double precision MHtree(4), SAtree

c used by FHRetrieveSMPara:	
		double precision invAlfa
c		double precision SUinvAlfa, AlfasMZ, GF
c		double precision ME, MU, MD, MM, MC, MS, ML, MB
c		double precision MW, MZ
c		double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
        parameter (nHiggsneut = 3)
        parameter (nHiggsplus = 1)
       	character(LEN=5) whichanalyses

c used by HiggsBounds_neutral_input_part
        double precision Mh_HB(3),GammaTotal_hj(3)
        integer CP_value(3)
        double precision  CS_lep_hjZ_ratio(3),          
     &   CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &   CS_lep_hjhi_ratio_nHbynH(3,3),               
     &   CS_gg_hj_ratio(3),CS_bb_hj_ratio(3),  
     &   CS_bg_hjb_ratio(3),                       
     &   CS_ud_hjWp_ratio(3),CS_cs_hjWp_ratio(3),
     &   CS_ud_hjWm_ratio(3),CS_cs_hjWm_ratio(3), 
     &   CS_gg_hjZ_ratio(3),
     &   CS_dd_hjZ_ratio(3),CS_uu_hjZ_ratio(3),
     &   CS_ss_hjZ_ratio(3),CS_cc_hjZ_ratio(3), 
     &   CS_bb_hjZ_ratio(3),                        
     &   CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),
     &   CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),
     &   CS_lhc8_vbf_ratio(3),CS_lhc8_tthj_ratio(3),
     &   BR_hjss(3),BR_hjcc(3),                         
     &   BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),                     
     &   BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),                     
     &   BR_hjgaga(3),BR_hjgg(3),
     &   BR_hjinvisible(3),BR_hjhihi_nHbynH(3,3)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &   CS_lep_HpjHmj_ratio(1),                    
     &   BR_tWpb,BR_tHpjb(1),                     
     &   BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1) 

c used in HiggsBounds
        double precision dmhneut_hb(nHiggsneut)
        double precision dmhch_hb(nHiggsplus)
c used in HiggsSignals (can be different)
        double precision dmhneut_hs(nHiggsneut)
        
c return values of run_HiggsBounds
c        integer HB_result,HB_chan,HB_ncombined
c        double precision HB_obsratio

c used by set_rate_uncertainties in HiggsSignals
        double precision dCS(5),dBR(5)
        double precision dggh, dbbh

c return values of run_HiggsSignals
c        double precision HS_Chisq_mu,HS_Chisq_mh,HS_Chisq,HS_Pvalue
c        integer HS_nobs
        
c run options in HiggsSignals
        integer pdf, output_level, mass_centered_method,runmode
             
c misc:
        integer i,j,as,t
        double precision norm,CW2!,Pi
        double precision
     &   g2hjbb(3),g2hjWW(3),g2hjZZ(3),                
     &   g2hjgg(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass 

c        Pi = 3.1415926535897932384626433832795029D0
