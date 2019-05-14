#!/bin/sh

program="${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/scripts/calc_A_and_H_xsec.sh"
SLHA_TEMPLATE="${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"

${program} ${SLHA_TEMPLATE}
