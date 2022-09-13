#!/usr/bin/env bash

CWD=$(pwd)

export PATH=${CWD}/packages/T3PS:$PATH

export THDM_T3PS_SCANNER_DIR=top_dir_
echo ${THDM_T3PS_SCANNER_DIR}

cd top_dir_

if [ -f env_local.sh ]; then
	 echo "Found env_local.sh, sourcing it..."
    source env_local.sh
fi

