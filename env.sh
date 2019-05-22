#!/usr/bin/env bash

CWD=$(pwd)

export PATH=${CWD}/packages/T3PS:$PATH

export THDM_T3PS_SCANNER_DIR=${CWD}

if [ -f env_local.sh ]; then
	 echo "Found env_local.sh, sourcing it..."
    source env_local.sh
fi

