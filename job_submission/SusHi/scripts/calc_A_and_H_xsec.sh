#!/bin/bash

SUSHI_BIN=/home/de3u14/lib/build/hep/SusHi/SusHi-1.6.1/bin/sushi.2HDMC

TEMPLATE=$1

TEMPLATE_A=template_A.slha
TEMPLATE_H=template_H.slha

sed "s/particle_/21/g" ${TEMPLATE} > ${TEMPLATE_A}
sed "s/particle_/12/g" ${TEMPLATE} > ${TEMPLATE_H}

${SUSHI_BIN} ./${TEMPLATE_A} sushi_A.out
${SUSHI_BIN} ./${TEMPLATE_H} sushi_H.out
