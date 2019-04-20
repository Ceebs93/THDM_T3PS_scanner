#!/bin/bash

echo "TAG: ${TAG}"
echo "Merging *.data"

cat ${ROOT_DIR}/${HEADER} $(find ./jobs/${TAG} -name "*.data") > ./jobs/${TAG}/${TAG}_all_data_with_sushi_merged.dat
