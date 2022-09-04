#!/bin/bash

# - Display information
echo -e "./utils/create-jobs.sh called\n"
echo -e "####################"
echo -e "NAME:        ${NAME}"
echo -e "CONFIG:     ${CONFIG}"
echo -e "nJobs:      ${nJobs}"
echo -e "TEMPLATE:   ${TEMPLATE}"
echo -e "ROOT_DIR:   ${ROOT_DIR}"
echo -e "####################"

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}

# - Create job project directory
mkdir -p ${JOB_PROJECT_DIR}

# - Section creates setup for cluster jobs
if [ ${CLUSTER} == "yes" ]; then

	# - Job creation loop
	for ((i=0;i<${nJobs};i++));
	do

		id=$(printf "%03d" ${i})
		JOB_DIR=${ROOT_DIR}/jobs/${NAME}/job_${id}

		mkdir -p ${JOB_DIR}
		cp ${ROOT_DIR}/${CONFIG} ${JOB_DIR}/t3ps.conf
   	   sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
   	   echo "program is set as" ${program}
   	   sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
   	   sed -i "s/chain_length_/${chain_length}/g" ${JOB_DIR}/t3ps.conf
		cp ${ROOT_DIR}/${TEMPLATE} ${JOB_DIR}/job.template

	   echo "${JOB_DIR}" >> ${ROOT_DIR}/jobs/${NAME}/all.jobs

done
fi

# - Section creates setup for local job
if [ ${CLUSTER} == "no" ]; then

	if [ ${nJobs} != 1 ]; then
		echo "Can only run one job at a time locally. Please edit the Makefile to reflect this."

	else	
		echo "Running job/s..."

		for ((i=0;i<${nJobs};i++));
		do
			JOB_DIR=${ROOT_DIR}/jobs/${NAME}/job_local
			mkdir -p ${JOB_DIR}
			echo "JOB_DIR is : ${JOB_DIR}"
			cp ${ROOT_DIR}/${CONFIG} ${JOB_DIR}/t3ps.conf
   		   sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
  		   echo "program is set as" ${program}
  		   sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
   		   sed -i "s/chain_length_/${chain_length}/g" ${JOB_DIR}/t3ps.conf
        		cp ${ROOT_DIR}/${TEMPLATE} ${JOB_DIR}/job.template
		done
	fi
fi
 
