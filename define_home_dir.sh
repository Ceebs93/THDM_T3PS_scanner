#!/usr/bin/env/bash

echo ${THDM_T3PS_SCANNER_DIR}

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/job.sh

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/config/mcmc_scan_default.conf

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/jobs_task/job_slurm.sh

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/jobs_task/job_torque.sh

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/local_test.sh

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/Makefile

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MadGraph/Makefile

sed -i "%s/top_dir_/${THDM_T3PS_SCANNER_DIR}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MadGraph/Makefile

#Reverse process

#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_}/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/job.sh
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/config/mcmc_scan_default.conf
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/jobs_task/job_slurm.sh
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/jobs_task/job_torque.sh
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/local_test.sh
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC/Makefile
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MadGraph/Makefile
#
#sed -i "%s/scratch/cb27g11/THDM_T3PS_SCANNER/top_dir_/g" ${THDM_T3PS_SCANNER_DIR}/job_submission/MadGraph/Makefile
