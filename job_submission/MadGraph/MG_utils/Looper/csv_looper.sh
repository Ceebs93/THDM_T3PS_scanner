#!/usr/bin/env bash

#I inherit the following variables from create-jobs.sh: CSV_NAME_ PROCESS_ BASECARD_ RUNCARD_ CARD_EDITOR_ ROOT_DIR_ JOB_PROJECT_DIR_ JOB_DIR_ DATA_RIPPER_VERSION_ RESULTS_

COUNT=0
extractions=0

{	#Double read to skip the header line of data csv
	read
	#IFS=internal field separator. Variables being run on are read in from the data csv. Note that the position of these variables is given at the end of this file.
	while IFS="," read -r mH mHc mA Itb Isinba
	do
		echo ${Isinba} ${Itb}
		printf -v sinba "%.6f \n" $Isinba
		printf -v tb "%.6f \n" $Itb
		echo ${sinba} ${tb}

		COUNT=$(( $COUNT + 1 ))
		echo "Line number is: " ${COUNT}

                python CARD_EDITOR_ $mH $mHc $mA "BASECARD_" "RUNCARD_" $tb $sinba "PROCESS_" "RESULTS_" # Here we run, and pass variables to, the inputcard editor. This edits the inputcard for MadGraph

                python ${THDM_T3PS_SCANNER_DIR}/packages/MG5_aMC_v3_1_0/bin/mg5_aMC RUNCARD_ # Here we run MG with the edited inputcard

		mov_d=$(echo PROCESS__${tb})
		mov_dirs=$(echo ${mov_d}_${sinba})
		#Had to do this in two steps or bash added erroneous spaces

		echo "Move directory is: ${mov_dirs}"
		mv RESULTS_/"$mov_dirs" JOB_DIR_

		echo "About to use 'Data_Ripper_iridis.py'"
                python DATA_RIPPER_VERSION_ $mov_dirs $Higgs # Here the needed data is extracted from the output
		
		extractions=$(( $extractions + 1 ))
		echo "Data_Ripper_iridis.py used ${extractions} time/s"
                
		cd RESULTS_
	        rm -r "$mov_dirs"
		echo "Removed ${mov_dirs}"

	done
 
	echo "data_coallator starting..."
	python JOB_DIR_/Data_coallator.py "RESULTS_" "PROCESS_"
echo
#Here -d allows us to specify our delimiter, then -f indicates we want to cut by field as opposed to bytes. The numbers indicate which column from the csv we want, and these are named above in the line 'while IFS="," read -r etc

} < <(cut -d "," -f2,3,4,6,7 CSV_NAME_.csv)

