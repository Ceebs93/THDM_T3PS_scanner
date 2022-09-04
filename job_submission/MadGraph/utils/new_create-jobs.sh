#!/bin/bash

# - Display information
echo -e "./utils/create-jobs.sh called\n"
echo -e "####################"
echo -e "NAME:        ${NAME}"
echo -e "INPUT_DATA:  ${INPUT_DATA}"
echo -e "LOCAL:       ${LOCAL}"
echo -e "SPLIT_NAME:  ${SPLIT_NAME}"
echo -e "PROCESS:     ${PROCESS}"
#echo -e "CONFIG:     ${CONFIG}"
echo -e "nJobs:       ${nJobs}"
echo -e "ROOT_DIR:    ${ROOT_DIR}"
echo -e "LOOPER:      ${LOOPER}"
echo -e "Card_Editor: ${Card_Editor}"
echo -e "BASECARD:    ${BASECARD}"
echo -e "TEMPLATE: ${TEMPLATE}"
echo -e "RIPPER:     ${RIPPER}"
echo -e "####################"

JOB_PROJECT_DIR=${ROOT_DIR}jobs/${NAME}

# - Create job project directory
mkdir -p ${JOB_PROJECT_DIR}

if [ ${LOCAL} == "yes" ]; then

	#Creating the job and corresponding folders
	JOB_DIR=${ROOT_DIR}jobs/${NAME}
	mkdir -p ${JOB_DIR}
	
	CSV_NAME="${JOB_PROJECT_DIR}/${SPLIT_NAME}"
	
	cp ${INPUT_DATA} ${CSV_NAME}.csv	

	echo -e ${CSV_NAME}

	# - Creating and editing 'input_editor' script which will input the different values of variables into the MadGraph runcards.
	cp ${ROOT_DIR}MG_utils/inputcard_editor/${Card_Editor} ${JOB_DIR}/${Card_Editor}
   sed -i "s;CARD_EDITOR_;${ROOT_DIR};g" ${JOB_DIR}/${Card_Editor}
   sed -i "s;JOB_PROJECT_DIR;${JOB_PROJECT_DIR};g" ${JOB_DIR}/${Card_Editor}

	# - Creating 'Data_coallator' for job.
	cp ${ROOT_DIR}MG_utils/Data_Ripper/Data_coallator.py ${JOB_DIR}/Data_coallator.py 


	# - Creating and editing the 'ripper' file which will extract desired values from MadGraph output
	cp ${ROOT_DIR}/MG_utils/Data_Ripper/${RIPPER} ${JOB_DIR}/Data_Ripper.py
   sed -i "s;JOB_DIR_;${JOB_DIR};g" ${JOB_DIR}/Data_Ripper.py
   sed -i "s;RESULTS_;${RESULTS};g" ${JOB_DIR}/Data_Ripper.py


	# - Creating and editing 'basecard' for MadGraph runs
	cp ${ROOT_DIR}MG_utils/mg_runcards/basecards/${BASECARD} ${JOB_DIR}/${BASECARD}
	BASECARDPATH=${JOB_DIR}/${BASECARD}
   sed -i "s;SCANNER_DIR_;${THDM_T3PS_SCANNER_DIR}/packages/;g" ${BASECARDPATH}


	# - Creating and editing the 'looper' file which will control the reading and feeding of a data csv to Madgraph
	cp ${ROOT_DIR}MG_utils/Looper/${LOOPER} ${JOB_DIR}/csv_looper.sh
	LOOPERPATH=${JOB_DIR}/csv_looper.sh

   sed -i "s;PROCESS_;${PROCESS};g" ${LOOPERPATH}
   sed -i "s;SCANNER_DIR_;${THDM_T3PS_SCANNER_DIR}/packages/;g" ${LOOPERPATH}
   sed -i "s;CSV_NAME_;${CSV_NAME};g" ${LOOPERPATH}
   sed -i "s;BASECARD_;${BASECARDPATH};g" ${LOOPERPATH}
   sed -i "s;CARD_EDITOR_;${JOB_DIR}/${Card_Editor};g" ${LOOPERPATH}
   sed -i "s;JOB_DIR_;${JOB_DIR};g" ${LOOPERPATH}
   sed -i "s;JOB_PROJECT_DIR_;${JOB_PROJECT_DIR};g" ${LOOPERPATH}
   sed -i "s;ROOT_DIR_;${ROOT_DIR};g" ${LOOPERPATH}
   sed -i "s;RUNCARD_;${JOB_DIR}/runcard.txt;g" ${LOOPERPATH}
   sed -i "s;DATA_RIPPER_VERSION_;${JOB_DIR}/Data_Ripper.py;g" ${LOOPERPATH}


else

	# - Calculate number of points per job
	nPts=$(tail --lines=+2 ${INPUT_DATA} | wc -l)
	nPtsPerJob=$(( ${nPts} / ${nJobs} ))
	nPtsPerJob=$(( nPtsPerJob + 1 ))


	# - Display information
	echo -e "nPts:       ${nPts}"
	echo -e "nPtsPerJob: ${nPtsPerJob}"


	# - Creating and editing python file to break csv into multiple parts for required number of jobs.
		cp ${ROOT_DIR}MG_utils/row_split.py ${JOB_PROJECT_DIR}/row_split.py
	   sed -i "s;INPUT_DATA_;${INPUT_DATA};g" ${JOB_PROJECT_DIR}/row_split.py
	   sed -i "s;nJobs_;${nJobs};g" ${JOB_PROJECT_DIR}/row_split.py
	   sed -i "s;datafile_;${JOB_PROJECT_DIR}/${SPLIT_NAME};g" ${JOB_PROJECT_DIR}/row_split.py

	split_no=$((nJobs - 1))
	echo ${split_no}
	tail -n +2 ${INPUT_DATA} | split -l ${nPtsPerJob} - split_
	declare -a csv_array=()
	
	for file in split_*
	do
	    
	    csv_array+=($file)
	    head -n 1 ${INPUT_DATA} > tmp_file
	    cat $file >> tmp_file
	    mv -f tmp_file ${JOB_PROJECT_DIR}/${file}.csv
	done

	rm ${ROOT_DIR}/split_*

	echo ${csv_array[@]}

	# - Running python to perform the splitting
#	python ${JOB_PROJECT_DIR}/row_split.py


	# - Job creation loop
	for ((i=0;i<${nJobs};i++));
	do

		id=$(printf "%03d" ${i})
		JOB_DIR=${ROOT_DIR}jobs/${NAME}/job_${id}
		mkdir -p ${JOB_DIR}
	
	#	csv_id=$((id + 1))
	#	CSV_NAME="${JOB_PROJECT_DIR}/${SPLIT_NAME}_${csv_id}"
		CSV_NAME="${csv_array[${id}]}"
		echo ${CSV_NAME}


		# - Creating and editing 'input_editor' script which will input the different values of variables into the MadGraph runcards.
		cp ${ROOT_DIR}MG_utils/inputcard_editor/${Card_Editor} ${JOB_DIR}/${Card_Editor}
   	   sed -i "s;CARD_EDITOR_;${ROOT_DIR};g" ${JOB_DIR}/${Card_Editor}
   	   sed -i "s;JOB_PROJECT_DIR;${JOB_PROJECT_DIR};g" ${JOB_DIR}/${Card_Editor}

		# - Creating 'Data_coallator' for job.
		cp ${ROOT_DIR}MG_utils/Data_Ripper/Data_coallator.py ${JOB_DIR}/Data_coallator.py 


		# - Creating and editing the 'ripper' file which will extract desired values from MadGraph output
		cp ${ROOT_DIR}/MG_utils/Data_Ripper/${RIPPER} ${JOB_DIR}/Data_Ripper.py
   	   sed -i "s;JOB_DIR_;${JOB_DIR};g" ${JOB_DIR}/Data_Ripper.py
   	   sed -i "s;RESULTS_;${RESULTS};g" ${JOB_DIR}/Data_Ripper.py


		# - Creating and editing 'basecard' for MadGraph runs
		cp ${ROOT_DIR}MG_utils/mg_runcards/basecards/${BASECARD} ${JOB_DIR}/${BASECARD}
		BASECARDPATH=${JOB_DIR}/${BASECARD}
   	   sed -i "s;SCANNER_DIR_;${THDM_T3PS_SCANNER_DIR}/packages/;g" ${BASECARDPATH}


		# - Creating and editing the 'looper' file which will control the reading and feeding of a data csv to Madgraph
		cp ${ROOT_DIR}MG_utils/Looper/${LOOPER} ${JOB_DIR}/csv_looper.sh
		LOOPERPATH=${JOB_DIR}/csv_looper.sh

   	   sed -i "s;PROCESS_;${PROCESS};g" ${LOOPERPATH}
   	   sed -i "s;SCANNER_DIR_;${THDM_T3PS_SCANNER_DIR}/packages/;g" ${LOOPERPATH}
   	   sed -i "s;CSV_NAME_;${CSV_NAME};g" ${LOOPERPATH}
  # sed -i "s;higgs_;${HIGGS};g" ${LOOPERPATH}
   	   sed -i "s;BASECARD_;${BASECARDPATH};g" ${LOOPERPATH}
   	   sed -i "s;CARD_EDITOR_;${JOB_DIR}/${Card_Editor};g" ${LOOPERPATH}
   	   sed -i "s;JOB_DIR_;${JOB_DIR};g" ${LOOPERPATH}
   	   sed -i "s;JOB_PROJECT_DIR_;${JOB_PROJECT_DIR};g" ${LOOPERPATH}
   	   sed -i "s;ROOT_DIR_;${ROOT_DIR};g" ${LOOPERPATH}
   	   sed -i "s;RUNCARD_;${JOB_DIR}/runcard.txt;g" ${LOOPERPATH}
   	   sed -i "s;DATA_RIPPER_VERSION_;${JOB_DIR}/Data_Ripper.py;g" ${LOOPERPATH}

	# - Creating and editing the submission script 
		cp ${ROOT_DIR}MG_utils/submission_script/${TEMPLATE} ${JOB_DIR}/${TEMPLATE}
   	   sed -i "s;LOOPER_;${LOOPERPATH};g" ${JOB_DIR}/${TEMPLATE}

		echo "${JOB_DIR}" >> ${ROOT_DIR}/jobs/${NAME}/all.jobs

done
fi
