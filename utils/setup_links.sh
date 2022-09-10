#!/usr/bin/env bash

ln -f -s ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.8.0/src/ ${THDM_T3PS_SCANNER_DIR}/links/inc/2HDMC
ln -f -s ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.8.0/lib/lib2HDMC.a ${THDM_T3PS_SCANNER_DIR}/links/lib/lib2HDMC.a
ln -f -s ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.8.0/lib_HBHS/lib2HDMC_HBHS.a ${THDM_T3PS_SCANNER_DIR}/links/lib/lib2HDMC_HBHS.a
ln -f -s ${THDM_T3PS_SCANNER_DIR}/packages/HiggsBounds-5.10.1/build/lib/libHB.a ${THDM_T3PS_SCANNER_DIR}/links/lib/libHB-5.10.1.a
ln -f -s ${THDM_T3PS_SCANNER_DIR}/packages/HiggsSignals-2.6.2/build/lib/libHS.a ${THDM_T3PS_SCANNER_DIR}/links/lib/libHS-2.6.2.a
