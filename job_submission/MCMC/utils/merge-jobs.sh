#!/bin/bash

if [ ${CONVERT_ONLY} == "no" ]; then

	 echo "TAG: ${TAG}"
	 echo "Merging data"
	 
	 tail --lines=+100 --quiet $(find ./jobs/${TAG} -name "*.chain.*") > ./jobs/${TAG}/${TAG}_all_data_merged.dat
	 
	 echo -e "Removing duplicate spaces/tabs..."
	 sed -i "s/[[:space:]]\+/ /g" ./jobs/${TAG}/${TAG}_all_data_merged.dat 
	 
	 echo -e "Removing junk"
	 awk -f ${ROOT_DIR}/utils/remove_junk.awk ./jobs/${TAG}/${TAG}_all_data_merged.dat > ./jobs/${TAG}/${TAG}_all_data_merged_and_polished.dat
	 
	 echo -e "Adding header"
	 cat ${ROOT_DIR}/${HEADER} ./jobs/${TAG}/${TAG}_all_data_merged_and_polished.dat > ./jobs/${TAG}/${TAG}_all_data_final.dat
	 
	 echo -e "Removing leading whitespace"
	 sed -i "s/^[ \t]*//" ./jobs/${TAG}/${TAG}_all_data_final.dat

fi

# - Convert ASCII to HDF5
echo -e "Convering ASCII to HDF"

if [ ${CONVERT} == "yes" ]; then
	 echo -e "Converting to HDF5..."
    ${ROOT_DIR}/utils/convert_to_hdf.py ./jobs/${TAG}/${TAG}_all_data_final.dat ./jobs/${TAG}/${TAG}_allowed_only_merged.h5f ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
#    ${ROOT_DIR}/utils/convert_to_hdf.py ./jobs/${TAG}/temp.dat ./jobs/${TAG}/${TAG}_allowed_only_merged.h5f ${DATASET_NAME} ${FORMAT} ${COMPRESSION}
fi
