#!/bin/sh

# Setting up which files should be used for the job
JOB_DIR=top_dir_/job_submission/MCMC/tests/local_test
program=top_dir/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_FAST
MCMC_CONFIG_TEMPLATE_FILE="top_dir_/job_submission/MCMC/config/mcmc_scan_default.conf"
FUNCTION_TEMPLATE="top_dir_/job_submission/MCMC/template/mcmc_scan.func"
nCores=8
chain_length=99


### ----

# Creating job directory and copying in templates of important files before
# replacing placeholders inside with actual paths/variables
mkdir -p ${JOB_DIR}
cp ${MCMC_CONFIG_TEMPLATE_FILE} ${JOB_DIR}/t3ps.conf
cp ${FUNCTION_TEMPLATE} ${JOB_DIR}/job.template
sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
sed -i "s/chain_length_/${chain_length}/g" ${JOB_DIR}/t3ps.conf

# Running T3PS with the configuration file supplied 
echo -ne '\n\n' | t3ps --debug -o ${JOB_DIR} ${JOB_DIR}/t3ps.conf
