#!/bin/sh

JOB_DIR=${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/tests/local_test
program=${THDM_T3PS_SCANNER_DIR}/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_FAST
MCMC_CONFIG_TEMPLATE_FILE="${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/config/mcmc_scan.conf"
FUNCTION_TEMPLATE="${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/template/mcmc_scan.func"
nCores=8
chain_length=100


### ----

mkdir -p ${JOB_DIR}
cp ${MCMC_CONFIG_TEMPLATE_FILE} ${JOB_DIR}/t3ps.conf
cp ${FUNCTION_TEMPLATE} ${JOB_DIR}/job.template
sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
sed -i "s/chain_length_/${chain_length}/g" ${JOB_DIR}/t3ps.conf

echo -ne '\n\n' | t3ps --debug -o ${JOB_DIR} ${JOB_DIR}/t3ps.conf
