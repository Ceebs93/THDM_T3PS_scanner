#!/bin/bash

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${TAG}

cd ${JOB_PROJECT_DIR}
cp ${ROOT_DIR}/${TASK} ./

echo 
echo "##################"
echo "TASK: ${TASK}"
echo "RESOURCES: ${RESOURCES}"
echo "##################"
echo 

while read -r DIR; do
	 echo "DIR: ${DIR}"
	 qsub -l ${RESOURCES} -v "DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}" ${ROOT_DIR}/${TASK}
done < ${LIST}

echo -e "\nJobs should be submitted now."
echo -e "Please check the output of qstat command below.\n"
qstat
