#!/bin/bash


OUTPUT_MERGED_FILE=./jobs/${TAG}/${TAG}_all_data_merged.dat
OUTPUT_MERGED_CLEANED_FILE=./jobs/${TAG}/${TAG}_all_data_merged_and_cleaned.dat
OUTPUT_FINAL_FILE=./jobs/${TAG}/${TAG}_all_data_final.dat
OUTPUT_HDF_FILE=./jobs/${TAG}/${TAG}_allowed_only_merged.h5f

if [ ${CONVERT_ONLY} == "no" ]; then

	 echo "TAG: ${TAG}"
	 echo "Merging data..."
	 
	 tail --lines=+1000 --quiet $(find ./jobs/${TAG} -name "*.chain.*") > ${OUTPUT_MERGED_FILE}
	 
	 echo -e "Removing duplicate spaces/tabs..."
	 sed -i "s/[[:space:]]\+/ /g" ${OUTPUT_MERGED_FILE}
	 
	 echo -e "Removing junk"
	 awk -f ${ROOT_DIR}/utils/remove_junk.awk ${OUTPUT_MERGED_FILE} > ${OUTPUT_MERGED_CLEANED_FILE}
	 
	 echo -e "Adding header"
	 cat ${ROOT_DIR}/${HEADER} ${OUTPUT_MERGED_CLEANED_FILE} > ${OUTPUT_FINAL_FILE}
	 
	 echo -e "Removing leading whitespace"
	 sed -i "s/^[ \t]*//" ${OUTPUT_FINAL_FILE}

#	 echo -e "Replacing spaces with tabs"
#	 sed -i "s/^//" ${OUTPUT_FINAL_FILE}

fi

# - Convert ASCII to HDF5
echo -e "Convering ASCII to HDF"

if [ ${CONVERT} == "yes" ]; then
	 echo -e "Converting to HDF5..."
    ${ROOT_DIR}/utils/convert_to_hdf.py ${OUTPUT_FINAL_FILE} ${OUTPUT_HDF_FILE} ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
#    ${ROOT_DIR}/utils/convert_to_hdf.py ./jobs/${TAG}/temp.dat ./jobs/${TAG}/${TAG}_allowed_only_merged.h5f ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
fi
