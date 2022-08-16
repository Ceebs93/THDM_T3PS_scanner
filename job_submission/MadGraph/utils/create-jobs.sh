#!/bin/bash

# - Display information
echo -e "./utils/create-jobs.sh called\n"
echo -e "####################"
echo -e "NAME:        ${NAME}"
echo -e "INPUT_DATA: ${INPUT_DATA}"
echo -e "CONFIG:     ${CONFIG}"
echo -e "nJobs:      ${nJobs}"
echo -e "TEMPLATE:   ${TEMPLATE}"
echo -e "ROOT_DIR:   ${ROOT_DIR}"
echo -e "####################"

JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}
JOB_PROJECT_DIR=${ROOT_DIR}/jobs/${NAME}


# - Create job project directory
mkdir -p ${JOB_PROJECT_DIR}

# - Calculate number of points per job
nPts=$(tail --lines=+2 ${INPUT_DATA} | wc -l)
nPtsPerJob=$(( ${nPts} / ${nJobs} ))
nPtsPerJob=$(( nPtsPerJob + 1 ))

# - Display information
echo -e "nPts:       ${nPts}"
echo -e "nPtsPerJob: ${nPtsPerJob}"

# - Split datafile into `nJobs` bits
tail --lines=+2 ${INPUT_DATA} | split --suffix-length=3 --lines=${nPtsPerJob} --numeric-suffixes - ${ROOT_DIR}/jobs/${TAG}/split_

# - Job creation loop
for ((i=0;i<${nJobs};i++));
do

	id=$(printf "%03d" ${i})
	JOB_DIR=${ROOT_DIR}/jobs/${NAME}/job_${id}

	mkdir -p ${JOB_DIR}
	cp ${ROOT_DIR}/${CONFIG} ${JOB_DIR}/t3ps.conf
   sed -i "s;program_;${program};g" ${JOB_DIR}/t3ps.conf
   sed -i "s/nCores_/${nCores}/g" ${JOB_DIR}/t3ps.conf
	cp ${ROOT_DIR}/${TEMPLATE} ${JOB_DIR}/template.slha
#	cp ${ROOT_DIR}/${TASK}     ${JOB_DIR}/job.sh
	mv ${JOB_PROJECT_DIR}/split_${id} ${JOB_DIR}/to_be_processed.dat

	echo "${JOB_DIR}" >> ${ROOT_DIR}/jobs/${NAME}/all.jobs

done
