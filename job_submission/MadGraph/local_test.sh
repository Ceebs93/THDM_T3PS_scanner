#!/bin/sh

JOB_DIR=top_dir_/job_submission/SusHi/tests/local_test
CONFIG_TEMPLATE_FILE="top_dir_/job_submission/SusHi/config/default.conf"
FUNCTION_TEMPLATE="top_dir_/job_submission/SusHi/template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
nCores=8
program="top_dir_/job_submission/SusHi/scripts/calc_A_and_H_xsec.sh"
T3PSDIR=top_dir_/packages/T3PS


### ----

mkdir -p ${JOB_DIR}
cp ${CONFIG_TEMPLATE_FILE} ${JOB_DIR}/t3ps.conf
cp ${FUNCTION_TEMPLATE} ${JOB_DIR}/template.slha
sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
sed -i "s;T3PSDIR_;${T3PSDIR};g" ${JOB_DIR}/t3ps.conf

cd ${JOB_DIR}
echo -ne '\n\n' | t3ps --debug -o ${JOB_DIR} ${JOB_DIR}/t3ps.conf
echo $(pwd)
