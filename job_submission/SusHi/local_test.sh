#!/bin/sh

JOB_DIR=./tests/local_test
CONFIG_TEMPLATE_FILE="./config/default.conf"
FUNCTION_TEMPLATE="template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
nCores=8
program="${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/scripts/calc_A_and_H_xsec.sh"


### ----

mkdir -p ${JOB_DIR}
cp ./${CONFIG_TEMPLATE_FILE} ${JOB_DIR}/t3ps.conf
cp ./${FUNCTION_TEMPLATE} ${JOB_DIR}/job.template
sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf

echo -ne '\n\n' | t3ps --debug -o ${JOB_DIR} ${JOB_DIR}/t3ps.conf
