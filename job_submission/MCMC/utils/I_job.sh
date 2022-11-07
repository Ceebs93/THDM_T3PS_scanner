#!/bin/bash

cd top_dir_

# Ensuring all nessecary paths/variables are in the $PATH
source top_dir_env.sh
echo "job_task/job.sh thinks scanner dir is: ${THDM_T3PS_SCANNER_DIR}"

# Script that will run the MCMCs
PROGRAM=top_dir_/packages/T3PS/t3ps

############################################
echo "Program: ${PROGRAM}"

# Moving to work directory
cd ${DIR}

CWD=$(pwd)
echo "Current dir: ${CWD}"

echo "Job starting.."
# Running T3PS with the configuration file supplied (recall that in setting up
# the chosen config file in the config directory was copied to each sub-job
# directory as 't3ps.conf'
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf
