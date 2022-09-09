#!/bin/bash

#cd /scratch/cb27g11/THDM_T3PS_scanner

source ${THDM_T3PS_SCANNER_DIR}/env.sh

cd $SLURM_SUBMIT_DIR

bash LOOPER_

