#!/bin/sh
# Creates working directory

# Protection for unset TAG, in order not to wipe out ./results directory.
if [ -z ${TAG+x} ]; then
	echo "Variable TAG is unset. Please specify TAG variable. e.g. make new TAG=test";
else
	if [ -d ./results/${TAG} ]; then
	echo "There is already a ./results/${TAG}. Making copy to ./backup";
	 	cp -rf ./results/${TAG}/ ./backup/${TAG}/;
		rm -rf results/${TAG};
	fi
	echo "Creating new working directory ./results/${TAG}.";
	mkdir -p results/${TAG}/figures/cross-section;
	mkdir -p results/${TAG}/figures/spectrum;
	mkdir -p results/${TAG}/figures/paramspace;
	mkdir -p results/${TAG}/LHA/;
fi
