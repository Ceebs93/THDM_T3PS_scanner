#!/usr/bin/env bash

#I inherit the following variables from create-jobs.sh: split_ab bq_tqh1 /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/basecard_type1.txt /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/runcard.txt /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/inputcard_editor1.py  /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/ /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001 /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/Data_Ripper.py RESULTS_

#Higgs = "h1"

COUNT=0
extractions=0

{	#Double read to skip the header line of data csv
	read
	#IFS=internal field separator. Variables being run on are read in from the data csv. Note that the position of these variables is given at the end of this file.
	while IFS="," read -r mH mA mHc Isinba Itb
	do
		echo ${Isinba} ${Itb}
		printf -v sinba "%.6f \n" $Isinba
		printf -v tb "%.6f \n" $Itb
		echo ${sinba} ${tb}

		COUNT=$(( $COUNT + 1 ))
		echo "Line number is: " ${COUNT}

                python /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/inputcard_editor1.py  $mH $mHc $mA "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/basecard_type1.txt" "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/runcard.txt" $tb $sinba "bq_tqh1" "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/" # Here we run, and pass variables to, the inputcard editor. This edits the inputcard for MadGraph

                python /scratch/cb27g11/THDM_T3PS_scanner/packages/MG5_aMC_v3_1_0/bin/mg5_aMC /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/runcard.txt # Here we run MG with the edited inputcard

		mov_d=$(echo bq_tqh1_${tb})
		mov_dirs=$(echo ${mov_d}_${sinba})
		#Had to do this in two steps or bash added erroneous spaces

		echo "Move directory is: ${mov_dirs}"
		mv /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/"$mov_dirs" /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/

		echo "About to use 'Data_Ripper_iridis.py'"
                python /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/Data_Ripper.py $mov_dirs $Higgs # Here the needed data is extracted from the output
		
		extractions=$(( $extractions + 1 ))
		echo "Data_Ripper_iridis.py used ${extractions} time/s"
                
		cd /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001
	        rm -r "$mov_dirs"
		echo "Removed ${mov_dirs}"

	done
 
	echo "data_coallator starting..."
	python /scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/Data_coallator.py "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/job_001/Data_Files/" "bq_tqh1"
echo
#Here -d allows us to specify our delimiter, then -f indicates we want to cut by field as opposed to bytes. The numbers indicate which column from the csv we want, and these are named above in the line 'while IFS="," read -r etc

} < <(cut -d "," -f1,2,3,4,6 split_ab.csv)

