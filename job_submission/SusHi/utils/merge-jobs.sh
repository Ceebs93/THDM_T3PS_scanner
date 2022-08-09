#!/bin/bash

echo "TAG: ${NAME}"
echo "Merging *.data"

cat ${ROOT_DIR}/${HEADER} $(find ./jobs/${NAME} -name "*.data") > ./jobs/${NAME}/${NAME}_all_data_with_sushi_merged.dat
