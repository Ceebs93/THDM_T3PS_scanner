#!/usr/bin/env bash

CWD=$(pwd)

export PATH=${CWD}/packages/T3PS:$PATH

export THDM_T3PS_SCANNER_DIR=/scratch/cb27g11/THDM_T3PS_scanner
echo ${THDM_T3PS_SCANNER_DIR}

cd /scratch/cb27g11/THDM_T3PS_scanner

if [ -f env_local.sh ]; then
	 echo "Found env_local.sh, sourcing it..."
    source env_local.sh
fi

