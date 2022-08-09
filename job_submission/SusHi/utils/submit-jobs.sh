#!/bin/bash

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}

cd ${JOB_PROJECT_DIR}
cp ${ROOT_DIR}/${TASK} ./

echo 
echo "##################"
echo "JOB_PROJECT_DIR: ${JOB_PROJECT_DIR}"
echo "TASK: ${TASK}"
echo "NODES: ${NODES}"
echo "PPN: ${PPN}"
echo "TIME: ${TIME}"
echo "##################"
echo 

while read -r DIR; do
	 echo "DIR: ${DIR}"

	 echo "command being run is: sbatch ${NODES} ${PPN} ${TIME} --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${ROOT_DIR}/${TASK}"

	 sbatch ${NODES} ${PPN} ${TIME} --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}, ${ROOT_DIR}/${TASK}
done < ${LIST}

echo -e "\nJobs should be submitted now."
echo -e "Please check the output of the squeue command below.\n"
squeue --me
