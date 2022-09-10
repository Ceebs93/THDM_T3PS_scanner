#!/bin/bash

#cd /scratch/cb27g11/THDM_T3PS_scanner

#source /scratch/cb27g11/THDM_T3PS_scanner/env.sh
echo "job_task/job.sh thinks scanner dir is: ${THDM_T3PS_SCANNER_DIR}"

PROGRAM=top_dir_/packages/T3PS/t3ps

############################################
echo "Program: ${PROGRAM}"

cd ${DIR}

CWD=$(pwd)
echo "Current dir: ${CWD}"

echo "Job starting.."
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf 
