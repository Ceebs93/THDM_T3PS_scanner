#!/bin/bash


PROGRAM=${THDM_T3PS_SCANNER_DIR}/packages/T3PS/t3ps

# - Local setup
module load gsl
source /home/de3u14/lib/build/miniconda3/bin/activate py27


############################################
echo "Program: ${PROGRAM}"

cd ${DIR}

CWD=$(pwd)
echo "Current dir: ${CWD}"

echo "Job starting.."
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf
