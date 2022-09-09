#!/bin/bash

# Need to create some kind of first time set up bash script that will go through and replace all these absolute paths with the relevant ones for the user. Perhaps it is worth going through the whole package and jotting down where all these paths occur.

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
