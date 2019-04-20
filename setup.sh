#!/usr/bin/env bash

CWD=$(pwd)

export PATH=${CWD}/packages/T3PS:$PATH

export THDM_T3PS_SCANNER_DIR=${CWD}

if [ -f setup_local.sh ]; then
	 echo "Found setup_local.sh, sourcing it..."
    source setup_local.sh
fi

