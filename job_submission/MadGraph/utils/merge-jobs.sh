#!/bin/bash

echo "NAME: ${NAME}"
echo "Merging *.csv"

echo "${ROOT_DIR}"
readarray -t array < <( find "jobs/${NAME}" -name "job_*")

PROJECT_DIR=${ROOT_DIR}jobs/${NAME}/

# These two lines are to explicitly name the cross-section column with the
#process involved to improve understanding later or when different processes
#are being added to the same file
cp ${ROOT_DIR}utils/merge-csv.py ${PROJECT_DIR}merge-csv.py
   sed -i "s;X_SECT_COL_;${PROCESS}_X_sects;g" ${PROJECT_DIR}merge-csv.py


for ((i=0; i<${#array[*]}; i++)); do
	
	CSV_PATH=${array[$i]}/Data_Files/${PROCESS}.csv
	echo "CSV_PATH is ${CSV_PATH}"
	echo " full path is ${ROOT_DIR}${array[$i]}/Data_Files/${PROCESS}.csv"  


	nPts=$(wc -l ${ROOT_DIR}${CSV_PATH} | awk '{print $1}')
	echo "nPts is : ${nPts}"

	echo $(head -n ${nPts} ${ROOT_DIR}${CSV_PATH} && tail -n+2 -q ${ROOT_DIR}${CSV_PATH} ) 

	head -n 1 ${ROOT_DIR}${CSV_PATH} > combined.out && tail -n+2 -q ${ROOT_DIR}${CSV_PATH} >> combined.csv

done

touch ${PROJECT_DIR}${PROCESS}_combined.csv

cat combined.out >> ${PROJECT_DIR}${PROCESS}_combined.csv
cat combined.csv >> ${PROJECT_DIR}${PROCESS}_combined.csv

rm ${ROOT_DIR}/combined*
	
python2 ${PROJECT_DIR}merge-csv.py "${PROJECT_DIR}${PROCESS}_combined.csv" "${INPUT_DATA}" "${PROCESS}" "${PROJECT_DIR}${NAME}_final.csv" "${MG_sin}" "${OG_sin}"
