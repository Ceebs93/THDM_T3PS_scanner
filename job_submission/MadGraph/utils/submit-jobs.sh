#!/bin/bash

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}

cd ${JOB_PROJECT_DIR}

echo 
echo "##################"
echo "JOB_PROJECT_DIR: ${JOB_PROJECT_DIR}"
echo "NODES:      ${NODES}"
echo "PPN:        ${PPN}"
echo "TIME:       ${TIME}"
echo "SCRIPT: ${SCRIPT}"
echo "USER_EMAIL: ${USER_EMAIL}"
echo "CLUSTERTYPE: ${CLUSTERTYPE}"
echo "##################"
echo 

if [ ${LOCAL} == "yes" ]; then
	
	# Can simply run the Looper from the job folder in this case
	cd ${JOB_PROJECT_DIR}
	bash csv_looper.sh  


else
	if [ ${CLUSTERTYPE}=="SLURM" ]; then
		while read -r DIR; do
	
			# Uncomment below for debugging
			# echo "command being run is: sbatch ${NODES} ${PPN} ${TIME} --mail-type=Begin,End --mail-user${USER_EMAIL}  --export=DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${DIR}/${SCRIPT}"
		
			 sbatch ${NODES} ${PPN} ${TIME} --mail-type=None --mail-user=${USER_EMAIL} --export=DIR=${DIR} --export=THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR} ${DIR}/${SCRIPT}

		done < ${LIST}
		echo -e "\nJobs should be submitted now."
		echo -e "Please check the output of the squeue command below.\n"
		squeue --me



	elif [ ${CLUSTERTYPE}=="TORQUE" ]; then
		while read -r DIR; do

		#	 echo "command being run is: qsub -l ${RESOURCES} -v "DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}" ${ROOT_DIR}/${TASK}"
			qsub -l ${RESOURCES} -v "DIR=${DIR},THDM_T3PS_SCANNER_DIR=${THDM_T3PS_SCANNER_DIR}" ${ROOT_DIR}/${TASK}

		done < ${LIST}
		echo -e "\nJobs should be submitted now."
		echo -e "Please check the output of qstat command below.\n"
		qstat
	fi
fi
