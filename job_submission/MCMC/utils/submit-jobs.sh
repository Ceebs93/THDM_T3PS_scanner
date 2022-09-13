#!/bin/bash

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}

echo "JOB_PROJECT_DIR is ${JOB_PROJECT_DIR}"
echo "ROOT_DIR/TASK IS ${ROOT_DIR/${TASK}}"

cd ${JOB_PROJECT_DIR}
cp ${ROOT_DIR}/${TASK} ./
echo "ROOT_DIR is : ${ROOT_DIR}"
echo "TASK is : ${TASK}"


echo 
echo "#####################"
echo "PWD:      $(pwd)"
echo "JOB_DIR:  ${JOB_DIR}"
echo "JOB_PROJECT_DIR: ${JOB_PROJECT_DIR}"
echo "NODES:    ${NODES}"
echo "PPN:      ${PPN}"
echo "TIME:     ${TIME}"
echo "TASK:     ${TASK}"
echo "CLUSTERTYPE: ${CLUSTERTYPE}"
echo "#####################"
echo 

if [ -f "${LIST}" ]; then
	echo "I have actually entered the correct bit"
	if [ ${CLUSTERTYPE} == "SLURM" ]; then
		while read -r DIR; do
			 echo "DIR: ${DIR}"

			 echo "command being run is: sbatch ${NODES} ${PPN} ${TIME} --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${ROOT_DIR}/${TASK}"

			 sbatch ${NODES} ${PPN} ${TIME} --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}, ${ROOT_DIR}/${TASK}
		done < ${LIST}
		echo -e "\nJobs should be submitted now."
		echo -e "Please check the output of the squeue command below.\n"
		squeue --me
		
	elif [ ${CLUSTERTYPE} == "TORQUE" ]; then
		while read -r DIR; do
			 echo "DIR: ${DIR}"
			 
			 echo "command being run is: qsub -l ${RESOURCES} -v "DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}" ${ROOT_DIR}/${TASK}"
			 qsub -l ${RESOURCES} -v "DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}" ${ROOT_DIR}/${TASK}

		done < ${LIST}
		echo -e "\nJobs should be submitted now."
		echo -e "Please check the output of qstat command below.\n"
		qstat
	fi

else
	DIR="${JOB_PROJECT_DIR}"
	echo "DIR: ${DIR}"
	sed -i "s;JOBDIR_;${DIR};g" ${DIR}/job_local.sh



	bash ${DIR}/job_local.sh --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}


fi
echo -e "\nJobs should be submitted now."
echo -e "Please check the output of the squeue command below.\n"
squeue --me
