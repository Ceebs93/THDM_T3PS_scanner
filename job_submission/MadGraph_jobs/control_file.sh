#!/bin/bash

#### --- Default variable values --- ####
run_choice="${run_choice:-'yes'}"

#### --- Gathering basic input from user --- ####
echo 'Starting setup for MadGraph jobs...'

echo 'Enter 1 to create a job and submit/run it, enter 2 to create ONLY, enter 3 to submit ONLY or enter 4 to merge previous job output.'
read run_type

echo 'Enter job_name'
read job_name

echo 'Enter the name of your process, e.g. "bg_twh3" as shorthand for bg -> t w A, or any other scheme that suits.'
read process_name

echo "### --- Makefile Settings --- ###" > 'utils/settings.mk'

# Here we define 3 variables as true or false, will_create, will_submit and will_merge, these determine what sections are written to the Makefile so that it can be automatically tailored to the user's needs without them needing to edit the Makefile themselves

if [ "${run_type}" == 1 ] ; then
	echo 'Create and submit job selected'		
	will_create=true 
	will_submit=true
	will_merge=false
fi

if [ "${run_type}" == 2 ] ; then
	echo 'Create job ONLY selected'		
	will_create=true
	will_submit=false
	will_merge=false
fi	

if [ "${run_type}" == 3 ] ; then
	echo 'Submit job ONLY selected'		
	will_create=false
	will_submit=true
	will_merge=false
fi

if [ "${run_type}" == 4 ] ; then
	echo 'Merge jobs selected'		
	will_create=false
	will_submit=false
	will_merge=true
fi

#### --- Create-jobs Section --- ###
if [ ${will_create} == true ] ; then
	echo 'Setting up create-jobs'

	echo 'Enter relative path to the input data you wish to use'
	read input_data

	echo 'Do you want to run a cluster job (enter yes for cluster, no for local)?'
	read run_choice

	echo 'Enter the ingoing particles, making sure to add a "\" before spaces e.g. "B\ g"'
	read proc_in

	echo 'Enter the outgoing particles, making sure to add a "\" before spaces e.g. T\ W"'
	read proc_out

	echo 'Enter number of cores to use'
	read ncores

	echo 'Enter number of jobs to run'
	read nJobs

	echo 'Enter name of looper to use, e.g. "csv_looper_1.sh"'
	read looper

	echo 'Enter name of runcard editor to use, e.g. "inputcard_editor_1.py"'
	read card_editor
	
	echo 'Enter name of basecard to use, e.g. "basecard_type2.txt"'
	read basecard

	echo 'Enter name of data_ripper to use, e.g. "Data_Ripper_typeII.py"'
	read ripper

	echo "CREATE_JOB_NAME            = ${job_name}" >> utils/settings.mk
	echo "CREATE_JOB_INPUT_DATA      = ${input_data}" >> utils/settings.mk
	echo "CREATE_JOB_CLUSTER         = ${run_choice}" >> utils/settings.mk
	echo "CREATE_JOB_PROCESS_NAME    = ${process_name}" >> utils/settings.mk
	echo "CREATE_JOB_PROC_IN         = ${proc_in}" >> utils/settings.mk
	echo "CREATE_JOB_PROC_OUT        = ${proc_out}" >> utils/settings.mk
	echo "CREATE_JOB_nCores          = ${ncores}" >> utils/settings.mk
	echo "CREATE_JOB_nJobs           = ${nJobs}" >> utils/settings.mk
	echo "CREATE_JOB_LOOPER          = ${looper}" >> utils/settings.mk
	echo "CREATE_JOB_Card_Editor     = ${card_editor}" >> utils/settings.mk
	echo "CREATE_JOB_BASECARD        = ${basecard}" >> utils/settings.mk
	echo "CREATE_JOB_RIPPER          = ${ripper}" >> utils/settings.mk

	make create-jobs

fi
######################################

#### --- Submit-jobs Section --- ####
if [ ${will_submit} == true ] ; then
	echo 'Setting up submit-jobs'

	if [ -z ${var+x} ]; then
		echo 'About to run a cluster job, is this okay (enter "no" to switch to local, or "yes" to continue)?'
		read run_choice
	fi

	if [ "$run_choice" == 'yes' ] ; then
		echo 'Please choose a cluster job manager option, Slurm or Torque. Enter "s" for Slurm or "t" for Torque'
		read cluster_choice_
	

		if [ "${cluster_choice_}" == "s" ] ; then	

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node (or per core)'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time
	

			echo "SUBMIT_JOB_NAME        = "${job_name}"" >> utils/settings.mk
			echo "SUBMIT_JOB_LIST        = "all.jobs"" >> utils/settings.mk
			echo "SUBMIT_JOB_NODES       = "--nodes=${node_no}"" >> utils/settings.mk
			echo "SUBMIT_JOB_PPN         = "--ntasks-per-node=${node_tasks}"" >> utils/settings.mk
			echo "SUBMIT_JOB_TIME        = "--time=${run_time}"" >> utils/settings.mk
			echo "SUBMIT_JOB_TASK        = "job_task/job_slurm.sh"" >> utils/settings.mk
			echo "SUBMIT_JOB_CLUSTERTYPE = "SLURM"" >> utils/settings.mk

		elif [ "${cluster_choice_}" == "t" ] ; then

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node (or per core)'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time

			echo "SUBMIT_JOB_LIST        = "all.jobs"" >> utils/settings.mk
			echo "SUBMIT_JOB_TAG         = "${job_name}"" >> utils/settings.mk
			echo "SUBMIT_JOB_RESOURCES   = "nodes=${node_no}:ppn=${node_tasks},walltime=${run_time}"" >> utils/settings.mk
			echo "SUBMIT_JOB_TASK        = "job_task/job_torque.sh"" >> utils/settings.mk
			echo "SUBMIT_JOB_CLUSTERTYPE = "TORQUE"" >> utils/settings.mk
	
		fi

	elif [ "${run_choice}" == "no" ] ; then

		echo 'Enter number of nodes'
		read node_no

		echo 'Enter the number of tasks to be assigned per node (or per core)'
		read node_tasks

		echo 'Enter run time for job/s (in format of hh:mm:ss)'
		read run_time
	
		echo "SUBMIT_JOB_LIST        = "all.jobs"" >> utils/settings.mk
		echo "SUBMIT_JOB_NAME        = "${job_name}"" >> utils/settings.mk
		echo "SUBMIT_JOB_NODES       = "--nodes=${node_no}"" >> utils/settings.mk
		echo "SUBMIT_JOB_PPN         = "--ntasks-per-node=${node_tasks}"" >> utils/settings.mk
		echo "SUBMIT_JOB_TIME        = "--time=${run_time}"" >> utils/settings.mk
		echo "SUBMIT_JOB_TASK        = "job.sh"" >> utils/settings.mk
	fi

	make submit-jobs
fi
#####################################


##### --- Merge-jobs Section --- #####
if [ "${will_merge}" == true ] ; then

	echo 'Setting up Merge-jobs'

	echo 'Enter path to your input data'
	read input_data

	echo 'Please enter the label of the first MadGraph variable for CSV combination'
	read mg_var1_label

	echo 'Now enter the equivalent label in the original, input CSV'
	read og_var1_label

	echo 'Please enter the label of the second MadGraph variable for CSV combination'
	read mg_var2_label

	echo 'Now enter the equivalent label in the original, input CSV'
	read og_var2_label
	
	echo 'Please enter the label of the third MadGraph variable for CSV combination'
	read mg_var3_label

	echo 'Now enter the equivalent label in the original, input CSV'
	read og_var3_label
	
	echo "MERGE_JOB_NAME         = "${job_name}"" >> utils/settings.mk
	echo "MERGE_JOB_PROCESS      = ${process_name}" >> utils/settings.mk
	echo "MERGE_JOB_INPUT_DATA   = ${input_data}" >> utils/settings.mk
	echo "MERGE_JOB_MGVAR1_LABEL = ${mg_var1_label}" >> utils/settings.mk
	echo "MERGE_JOB_OGVAR1_LABEL = ${og_var1_label}" >> utils/settings.mk
	echo "MERGE_JOB_MGVAR2_LABEL = ${mg_var2_label}" >> utils/settings.mk
	echo "MERGE_JOB_OGVAR2_LABEL = ${og_var2_label}" >> utils/settings.mk
	echo "MERGE_JOB_MGVAR3_LABEL = ${mg_var3_label}" >> utils/settings.mk
	echo "MERGE_JOB_OGVAR3_LABEL = ${og_var3_label}" >> utils/settings.mk

	make merge-jobs

fi
