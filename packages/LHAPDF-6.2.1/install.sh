#!/usr/bin/env bash

if [ -z ${THDM_T3PS_SCANNER_DIR} ]; then
	echo "Source setup.sh in the root dir first!"
	exit
fi

./configure --prefix=${THDM_T3PS_SCANNER_DIR}/packages/LHAPDF-6.2.1-install/
make
make install
