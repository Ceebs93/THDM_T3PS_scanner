#!/bin/sh

source /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/MG_utils/MGSource.sh
echo -e "looking for env.sh ${THDM_T3PS_DIR}/env.sh"
#conda activate MGenv

#module load gcc/11.1.0

path=pwd

echo -e "I am here: ${PATH}"

cd $SLURM_SUBMIT_DIR

location=pwd
echo -e "I am here: ${location}"

bash LOOPER_

