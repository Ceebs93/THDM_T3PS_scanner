#!/bin/bash

PROGRAM=${THDM_T3PS_SCANNER_DIR}/packages/T3PS/t3ps


# - Local setup
module load gsl
module load gcc/6.1.0
source /home/de3u14/lib/build/miniconda3/bin/activate py27

############################################

cd ${DIR}
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf
