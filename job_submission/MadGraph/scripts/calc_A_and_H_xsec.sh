#!/bin/bash

# -- Input:
#  - template.slha: an .slha config file with particle_ token


SUSHI_BIN=${THDM_T3PS_SCANNER_DIR}/packages/SusHi-1.6.1/bin/sushi.2HDMC

TEMPLATE=$1

TEMPLATE_A=template_A.slha
TEMPLATE_H=template_H.slha

sed "s/particle_/21/g" ${TEMPLATE} > ${TEMPLATE_A}
sed "s/particle_/12/g" ${TEMPLATE} > ${TEMPLATE_H}

${SUSHI_BIN} ./${TEMPLATE_A} sushi_A.out
${SUSHI_BIN} ./${TEMPLATE_H} sushi_H.out
