#!/bin/bash

PROGRAM=/home/de3u14/lib/build/hep/T3PS/T3PS-1.0.2/src/t3ps

module load gsl
module load gcc/6.1.0
source /home/de3u14/lib/build/miniconda/envs/py27/bin/activate py27

############################################

cd ${DIR}
echo -ne '\n\n' | ${PROGRAM} -o ./ t3ps.conf
