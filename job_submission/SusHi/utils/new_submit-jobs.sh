#!/bin/bash

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}

cd ${JOB_PROJECT_DIR}
#cp ${ROOT_DIR}/${TASK} ./

echo 
echo "##################"
echo "JOB_PROJECT_DIR: ${JOB_PROJECT_DIR}"
echo "TASK:       ${TASK}"
echo "NODES:      ${NODES}"
echo "PPN:        ${PPN}"
echo "TIME:       ${TIME}"
echo "SCRIPT: ${SCRIPT}"
echo "USER_EMAIL: ${USER_EMAIL}"
echo "##################"
echo 



while read -r DIR; do
	 echo "DIR: ${DIR}"
	 echo "DIR/SCRIPT: ${DIR}/${SCRIPT}"
	 echo "command being run is: sbatch ${NODES} ${PPN} ${TIME} --mail-type=Begin,End --mail-user${USER_EMAIL}  --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${DIR}/${SCRIPT}"

	 sbatch ${NODES} ${PPN} ${TIME} --mail-type=None --mail-user=${USER_EMAIL} --export=DIR=${DIR} --export=THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${DIR}/${SCRIPT}
done < ${LIST}

echo -e "\nJobs should be submitted now."
echo -e "Please check the output of the squeue command below.\n"
squeue --me
