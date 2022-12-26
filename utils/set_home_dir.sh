#!/usr/bin/env/bash

echo ${THDM_T3PS_SCANNER_DIR}

# Replacing placeholder "top_dir_" in key files with the path ${THDM_T3PS_SCANNER_DIR}. This avoids a requirement to source the env.sh from the top directory every time the user wishes to use the scanner. Running this means that all the files that need to know the top directory will have that saved to them.

echo -e "setting up env.sh"
cp ${THDM_T3PS_SCANNER_DIR}/utils/I_env.sh ${THDM_T3PS_SCANNER_DIR}/env.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${THDM_T3PS_SCANNER_DIR}/env.sh

# This section copies the Madgraph job section template files and replaces the placeholder term with the user path
echo -e "setting up Madgraph files"
MadGraph_dir=${THDM_T3PS_SCANNER_DIR}/job_submission/MadGraph_jobs/
MG_sub_dir=${MadGraph_dir}MG_utils/submission_script/
cp ${MG_sub_dir}I_slurm_submission.sh ${MG_sub_dir}slurm_submission.sh
cp ${MG_sub_dir}I_torque_submission.pbs ${MG_sub_dir}torque_submission.pbs
cp ${MadGraph_dir}utils/I_Makefile ${MadGraph_dir}Makefile
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MadGraph_dir}Makefile
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MG_sub_dir}slurm_submission.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MG_sub_dir}torque_submission.pbs

# This section copies the MCMC scanner template files and replaces the placeholder term with the user path
echo -e "setting up MCMC files"
MCMC_dir=${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/
cp ${MCMC_dir}/utils/I_local_test.sh ${MCMC_dir}local_test.sh
cp ${MCMC_dir}/utils/I_job.sh ${MCMC_dir}job.sh
cp ${MCMC_dir}/utils/I_Makefile ${MCMC_dir}Makefile
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}job.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}local_test.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}Makefile

# This section copies more scanner template files and replaces the placeholder term with the user path
cp ${MCMC_dir}config/I_hybrid_scan.conf ${MCMC_dir}config/hybrid_scan.conf
cp ${MCMC_dir}config/I_generic_scan.conf ${MCMC_dir}config/generic_scan.conf
cp ${MCMC_dir}job_task/I_job_torque.sh ${MCMC_dir}job_task/job_torque.sh
cp ${MCMC_dir}job_task/I_job_slurm.sh ${MCMC_dir}job_task/job_slurm.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}config/hybrid_scan.conf
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}config/generic_scan.conf
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}job_task/job_torque.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${MCMC_dir}job_task/job_slurm.sh

# This section copies the template files for external packages and replaces the placeholder term with the user path
Initial_pkg_dir=${THDM_T3PS_SCANNER_DIR}/packages/startup_versions/
cp ${Initial_pkg_dir}I_install.sh ${THDM_T3PS_SCANNER_DIR}/packages/install.sh
cp ${Initial_pkg_dir}I_HBHS.cpp ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.8.0/src/HBHS.cpp
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${THDM_T3PS_SCANNER_DIR}/packages/install.sh
sed -i "s;top_dir_;${THDM_T3PS_SCANNER_DIR};g" ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.8.0/src/HBHS.cpp

