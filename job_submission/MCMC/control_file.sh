#!/bin/bash

#### --- Default variable values --- ####
run_choice="${run_choice:-'yes'}"

#### --- Gathering basic input from user --- ####
echo 'Starting setup for MCMC jobs...'

echo 'Enter 1 to create a job and submit/run it, enter 2 to create ONLY, enter 3 to submit ONLY or enter 4 to merge previous job output.'
read run_type

echo 'Enter job_name'
read job_name

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

	echo 'Enter chosen parameterisation basis: higgs, hhg (higgs hunters guide), generic, phys (physical) or hybrid'
	read basis

	echo 'Do you want to run a cluster job (enter yes for cluster, no for local)?'
	read run_choice

	echo 'Enter number of cores to use'
	read ncores

	echo 'Enter number of jobs to run'
	read nJobs
	
	echo 'What 2HDM model type do you want to use, 1, 2, 3 or 4?'
	read Y

	echo 'Enter chain length'
	read chain_len

	echo "CREATE_JOB_NAME         = ${job_name}" >> utils/settings.mk
	echo "CREATE_JOB_CONFIG       = config/${basis}_scan.conf" >> utils/settings.mk
	echo "CREATE_JOB_CLUSTER      = ${run_choice}" >> utils/settings.mk
	echo "CREATE_JOB_nCores       = ${ncores}" >> utils/settings.mk
	echo "CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_${basis}_FAST'" >> utils/settings.mk
	echo "CREATE_JOB_nJobs        = ${nJobs}" >> utils/settings.mk
	echo "CREATE_JOB_chain_length = ${chain_len}" >> utils/settings.mk
	echo "CREATE_JOB_TEMPLATE     = template/${basis}_basis_scan.func" >> utils/settings.mk
	echo "CREATE_JOB_Y            = ${Y}" >> utils/settings.mk
	echo "             "

	make create-jobs

fi
######################################

#### --- Submit-jobs Section --- ####
if [ ${will_submit} == true ] ; then
	echo '                      '
	echo 'Setting up submit-jobs'

	if [ -z ${var+x} ]; then
		echo 'About to run a cluster job, is this correct (enter "no" to switch to local, or "yes" to continue)?'
		read run_choice
	fi

	if [ "$run_choice" == 'yes' ] ; then
		echo 'Please choose a cluster job manager option, Slurm or Torque. Enter "s" for Slurm or "t" for Torque'
		read cluster_choice_
	

		if [ "${cluster_choice_}" == "s" ] ; then	

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time
	

			echo "SUBMIT_JOB_LIST        = "all.jobs"" >> utils/settings.mk
			echo "SUBMIT_JOB_NAME        = "${job_name}"" >> utils/settings.mk
			echo "SUBMIT_JOB_NODES       = "--nodes=${node_no}"" >> utils/settings.mk
			echo "SUBMIT_JOB_PPN         = "--ntasks-per-node=${node_tasks}"" >> utils/settings.mk
			echo "SUBMIT_JOB_TIME        = "--time=${run_time}"" >> utils/settings.mk
			echo "SUBMIT_JOB_TASK        = "job_task/job_slurm.sh"" >> utils/settings.mk
			echo "SUBMIT_JOB_CLUSTERTYPE = "SLURM"" >> utils/settings.mk

			echo "                       "

		elif [ "${cluster_choice_}" == "t" ] ; then

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node'
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

		echo '                     '

		echo 'Enter number of nodes'
		read node_no

		echo 'Enter the number of tasks to be assigned per node'
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

	echo "                       "

	make submit-jobs
fi
#####################################


##### --- Merge-jobs Section --- #####
if [ "${will_merge}" == true ] ; then

	echo '                     '

	echo 'Setting up Merge-jobs'

	echo 'Do you want to convert your output into hd5f format? Answer "yes" or "no"'
	read to_convert

	if [ "${to_convert}" == "no" ] ; then
		only_convert="no"
	
	else
		echo 'Do you ONLY want to convert your output into hd5f (i.e. not perform merging of jobs)?Answer "yes" or "no"'
		read only_convert
	fi

	echo 'Do you want to create a CSV of your merged job output? Answer "yes" or "no"'
	read make_csv

	echo "MERGE_JOB_NAME         = "${job_name}"" >> utils/settings.mk
	echo "MERGE_JOB_HEADER       = "header/default.header"" >> utils/settings.mk
	echo "MERGE_JOB_DATASET_NAME = "${job_name}"" >> utils/settings.mk
	echo "MERGE_JOB_CONVERT      = "${to_convert}"" >> utils/settings.mk
	echo "MERGE_JOB_CONVERT_ONLY = "${only_convert}"" >> utils/settings.mk
	echo "MERGE_JOB_MAKE_CSV     = "${make_csv}"" >> utils/settings.mk
	echo "MERGE_JOB_FORMAT       = "table"" >> utils/settings.mk
	echo "MERGE_JOB_COMPRESSION  = "blosc"" >> utils/settings.mk

	echo "                      "

	make merge-jobs

fi
