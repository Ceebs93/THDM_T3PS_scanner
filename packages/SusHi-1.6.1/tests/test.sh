#!/bin/bash


SUSHI_DIR="${THDM_T3PS_SCANNER_DIR}/packages/SusHi-1.6.1"
SUSHI_BIN="${SUSHI_DIR}/bin/sushi"




# - Physical basis
#${SUSHI_BIN} ./2HDMC_physicalbasis.in ./2HDMC_physicalbasis.out

# - Hybrid basis
${SUSHI_BIN} ./2HDMC_hybrid.in ./2HDMC_hybrid.out
