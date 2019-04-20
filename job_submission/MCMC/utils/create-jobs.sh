#!/bin/bash

# - Display information
echo -e "./utils/create-jobs.sh called\n"
echo -e "####################"
echo -e "TAG:        ${TAG}"
echo -e "CONFIG:     ${CONFIG}"
echo -e "nJobs:      ${nJobs}"
echo -e "TEMPLATE:   ${TEMPLATE}"
echo -e "ROOT_DIR:   ${ROOT_DIR}"
echo -e "####################"

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${TAG}

# - Create job project directory
mkdir -p ${JOB_PROJECT_DIR}

# - Job creation loop
for ((i=0;i<${nJobs};i++));
do

	id=$(printf "%03d" ${i})
	JOB_DIR=${ROOT_DIR}/jobs/${TAG}/job_${id}

	mkdir -p ${JOB_DIR}
	cp ${ROOT_DIR}/${CONFIG} ${JOB_DIR}/t3ps.conf
   sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
   sed -i "s/chain_length_/${chain_length}/g" ${JOB_DIR}/t3ps.conf
	cp ${ROOT_DIR}/${TEMPLATE} ${JOB_DIR}/job.template

	echo "${JOB_DIR}" >> ${ROOT_DIR}/jobs/${TAG}/all.jobs

done
