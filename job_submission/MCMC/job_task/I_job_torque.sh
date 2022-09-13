#!/bin/bash

cd top_dir_

source top_dir_/env.sh
echo "job_task/job.sh thinks scanner dir is: ${THDM_T3PS_SCANNER_DIR}"

PROGRAM=top_dir_/packages/T3PS/t3ps

############################################
echo "Program: ${PROGRAM}"

cd ${DIR}

echo "Job starting.."
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf 
