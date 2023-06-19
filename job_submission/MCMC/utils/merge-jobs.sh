#!/bin/bash

# Setting paths for later
OUTPUT_MERGED_FILE=./jobs/${NAME}/${NAME}_all_data_merged.dat
OUTPUT_MERGED_CLEANED_FILE=./jobs/${NAME}/${NAME}_all_data_merged_and_cleaned.dat
OUTPUT_FINAL_FILE=./jobs/${NAME}/${NAME}_all_data_final.dat
OUTPUT_HDF_FILE=./jobs/${NAME}/${NAME}_allowed_only_merged.h5f

echo "output file is called: ${OUTPUT_MERGED_FILE}"

# Checks to see if user wants to just convert .dat file to a hdf5 file
if [ ${CONVERT_ONLY} == "n" ]; then
	
	 # Gathers the .dat files from all sub-jobs in current job and merges them
	 echo "NAME: ${NAME}"
	 echo "Merging data..."
	 
	 #Starts from end of all *chain.* files and reads lines to write in ${OUTPUT_MERGED_FILE}
	 tail -n 1000 --quiet $(find ./jobs/${NAME} -name "*.chain.*") > ${OUTPUT_MERGED_FILE}
	 
	 echo -e "Removing duplicate spaces/tabs..."
	 sed -i "s/[[:space:]]\+/ /g" ${OUTPUT_MERGED_FILE}
	 
	 echo -e "Removing junk"
	 awk -f ${ROOT_DIR}/utils/remove_junk.awk ${OUTPUT_MERGED_FILE} > ${OUTPUT_MERGED_CLEANED_FILE}
	 
	 echo -e "Adding header"
	 cat ${ROOT_DIR}/${HEADER} ${OUTPUT_MERGED_CLEANED_FILE} > ${OUTPUT_FINAL_FILE}
	 
	 echo -e "Removing leading whitespace"
	 sed -i "s/^[ \t]*//" ${OUTPUT_FINAL_FILE}


fi

# Creates a csv file from the .dat file if user has indicated to do this
if [ ${MAKE_CSV} == "y" ]; then

	echo "NAME: ${NAME}"
	echo "Running data_to_csv.py"

	# We must pass the variables in an array as bash will only pass one variable to python after the name of the script we are running.
	python_input=("${NAME}" "${OUTPUT_MERGED_FILE}" "${Y}" "${BASIS}")
	echo ${python_input[@]}

	python utils/data_to_csv.py ${python_input[@]}

fi


# - Convert ASCII to HDF5
echo -e "Convering ASCII to HDF"
# Creates a hdf5 file from combined .dat files
if [ ${CONVERT} == "y" ]; then
	 echo -e "Converting to HDF5..."
    ${ROOT_DIR}/utils/convert_to_hdf.py ${OUTPUT_FINAL_FILE} ${OUTPUT_HDF_FILE} ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
#    ${ROOT_DIR}/utils/convert_to_hdf.py ./jobs/${TAG}/temp.dat ./jobs/${TAG}/${TAG}_allowed_only_merged.h5f ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
fi
