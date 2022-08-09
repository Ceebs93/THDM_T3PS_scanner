#!/bin/bash

if [ $# -eq 0 ]; then
    echo "No input provided from control file, running with nwmag_looper values"

    csv_name="Type1_reruns"

    Process="bq_tqhi_type1_again"
    Higgs="tot"
    basefilepath="/scratch/cb27g11/mg_run_basic/nw_mg_runT1.txt"
    testfilepath="/scratch/cb27g11/mg_run_iridis/nwmag_runcard.txt"

elif [ $# -eq 5 ]; then
    echo "Input provided from control file, running with these values"

    csv_name=$1
    
    Process=$2
    Higgs=$3
    basefilepath=$4
    testfilepath=$5

else
  echo "Please provide exactly 5 inputs, csv_name, process name, higgs combination, filepath to basic MG runcard and filepath to functional MG runcard; or run directly from nwmag_looper.sh"
  exit
fi

COUNT=0
ripruns=0

#python /scratch/cb27g11/dat_to_DF.py $file1 $file2 $csv_name
{	#Double read means we skip the header line
	read
	#IFS=internal field separator
	while IFS="," read -r mH mA mHc Isinba Itb
	do
		echo ${Isinba} ${Itb}
		printf -v sinba "%.6f \n" $Isinba
		printf -v tb "%.6f \n" $Itb
		echo ${sinba} ${tb}

		COUNT=$(( $COUNT + 1 ))
		echo "Line number is: " ${COUNT}
                 python /scratch/cb27g11/inputcard_editor/nwinputcard_editor.py $mH $mHc $mA $testfilepath $tb $sinba $Process $basefilepath # Here we run, and pass variables to, inputcard_editor. This edits the inputcard for MG

                python /scratch/cb27g11/MG5_aMC_v3_1_0/bin/mg5_aMC /scratch/cb27g11/mg_run_iridis/nwmag_runcard.txt # Here we run MG with the edited inputcard

		mov_d=$(echo ${Process}_${tb})
		mov_dirs=$(echo ${mov_d}_${sinba})
		#Had to do this in two steps or bash added erroneous spaces

		echo "Move directory is: ${mov_dirs}"
		mv /scratch/cb27g11/"$mov_dirs" /scratch/cb27g11/Data_Storage

              # mv /scratch/cb27g11/"$Process"* /scratch/cb27g11/Data_Storage # Here we move the output folder into Data_Storage
		echo "About to use 'Data_Ripper_iridis.py'"
                python /scratch/cb27g11/Data_Ripper_iridis/Data_Ripper_iridis.py $mov_dirs $Higgs # Here the needed data is extracted from the output
		
		ripruns=$(( $ripruns + 1 ))
		echo "Data_Ripper_iridis.py used ${ripruns} time/s"
                
		cd RESULTS_
	        rm -r "$mov_dirs"
		echo "Removed ${mov_dirs}"
               #rm -r "$Process"*
               #bash /scratch/cb27g11/Data_Storage/Compressor.sh # Original output is compressed to save storage space
               # cd ~/
	done

	results="/scratch/cb27g11/Data_Ripper_files/" 
	echo "data_coallator starting..."
	python /scratch/cb27g11/Data_Ripper_iridis/Data_coallator.py $results/results $Process
echo
#Here -d allows us to specify our delimiter, then -f indicates we want to cut by field as opposed to bytes etc

} < <(cut -d "," -f1,2,3,4,6 /scratch/cb27g11/Looping_Bash/${csv_name}.csv)

