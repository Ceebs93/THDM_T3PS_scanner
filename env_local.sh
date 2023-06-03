#!/usr/bin/env/bash

# Load Python environment
module purge
module load conda
source activate THDM

# Load libraries and compiler
module load gcc/12.1.0
module load gsl/2.7
echo "python version in use is"
PV=$(python --version)
echo "${PV}"
# Load THDM and T3PS
echo "Scanner Dir: ${THDM_T3PS_SCANNER_DIR}"                                             
