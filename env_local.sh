#!/usr/bin/env/bash

# Load Python environment

conda activate THDM

# Load libraries and compiler
echo "python version in use is"
PV=$(python --version)
echo "${PV}"
# Load THDM and T3PS
echo "Scanner Dir: ${THDM_T3PS_SCANNER_DIR}"                                             
