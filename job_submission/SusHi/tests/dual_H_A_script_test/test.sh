#!/bin/sh

#program="${THDM_T3PS_SCANNER_DIR}/packages/SusHi-1.6.1/bin/sushi.2HDMC"
program="${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/scripts/calc_A_and_H_xsec.sh"
SLHA_TEMPLATE="./dual_A_and_H.slha"

${program} ${SLHA_TEMPLATE}
