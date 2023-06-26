#!/usr/bin/env/bash

# Load Python environment
source activate THDM

# Check Python version and Scanner Dir
echo "python version in use is"
PV=$(python --version)
echo "${PV}"
echo "Scanner Dir: ${THDM_T3PS_SCANNER_DIR}"                                             
