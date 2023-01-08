
##### --- Create and Submit --- #####
if [ "${run_type}" ==  1 ] ; then
	echo 'Create and submit job selected'		

	echo 'Enter chosen parameterisation basis: higgs, hhg (higgs hunters guide), generic, phys (physical) or hybrid'
	read basis

	echo 'Enter number of cores to use'
	read ncores

	echo 'Enter chain length'
	read chain_len

	echo "CREATE_JOB_NAME         = '${job_name}' " >> utils/settings.mk
	echo "CREATE_JOB_CONFIG       = 'config/${basis}_scan.conf' " >> utils/settings.mk
	echo "CREATE_JOB_CLUSTER      = '${run_choice}' " >> utils/settings.mk
	echo "CREATE_JOB_nCores       = '${ncores}' " >> utils/settings.mk
	echo "CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_${basis}_FAST' " >> utils/settings.mk
	echo "CREATE_JOB_nJobs        = '${njobs}' " >> utils/settings.mk
	echo "CREATE_JOB_chain_length = '${chain_len}' " >> utils/settings.mk
	echo "CREATE_JOB_nJobs        = 'template/${basis}_basis_scan.func' " >> utils/settings.mk

	echo 'Do you want to run a cluster job (enter yes for cluster, no for local)?'
	read run_choice

	if [ "${run_choice}" == 'yes' ] ; then
		echo 'Please choose a cluster job manager option, Slurm or Torque. Enter "s" for Slurm or "t" for Torque'
		read cluster_choice_
	
		if [ "${cluster_choice_}" == "s" ] ; then		

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node (or per core)'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time
	
			echo "CREATE_JOB_chain_length = '${chain_len}' " >> utils/settings.mk

			echo "SUBMIT_JOB_LIST        = "all.jobs""
			echo "SUBMIT_JOB_NAME        = "${job_name}""
			echo "SUBMIT_JOB_NODES       = "--nodes=${node_no}""
			echo "SUBMIT_JOB_PPN         = "--ntasks-per-node=${node_tasks}""
			echo "SUBMIT_JOB_TIME        = "--time=${run_time}""
			echo "SUBMIT_JOB_TASK        = "job_task/job_slurm.sh""
			echo "SUBMIT_JOB_CLUSTERTYPE = "SLURM""

		elif [ "${cluster_choice_}" == "t" ] ; then

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node (or per core)'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time

			echo "SUBMIT_JOB_LIST        = "all.jobs""
			echo "SUBMIT_JOB_TAG         = "${job_name}""
			echo "SUBMIT_JOB_RESOURCES   = "nodes=${node_no}:ppn=${node_tasks},walltime=${run_time}""
			echo "SUBMIT_JOB_TASK        = "job_task/job_torque.sh""
			echo "SUBMIT_JOB_CLUSTERTYPE = "TORQUE""
	
		fi

	elif [ "${run_choice}" == "no" ] ; then

		echo 'Enter number of nodes'
		read node_no

		echo 'Enter the number of tasks to be assigned per node (or per core)'
		read node_tasks

		echo 'Enter run time for job/s (in format of hh:mm:ss)'
		read run_time
	
		echo "SUBMIT_JOB_LIST        = "all.jobs""
		echo "SUBMIT_JOB_NAME        = "${job_name}""
		echo "SUBMIT_JOB_NODES       = "--nodes=${node_no}""
		echo "SUBMIT_JOB_PPN         = "--ntasks-per-node=${node_tasks}""
		echo "SUBMIT_JOB_TIME        = "--time=${run_time}""
		echo "SUBMIT_JOB_TASK        = "job.sh""
	fi

	make create-jobs 
	make submit-jobs

  
##### --- Create Only --- #####
elif [ "${run_type}" == 2 ] ; then
	echo 'Create ONLY selected'  

	echo 'Enter chosen parameterisation basis: higgs, hhg (higgs hunters guide), generic, phys (physical) or hybrid'
	read basis

	echo 'Enter number of cores to use'
	read ncores

	echo 'Enter chain length'
	read chain_len

	echo "CREATE_JOB_NAME         = '${job_name}' " >> utils/settings.mk
	echo "CREATE_JOB_CONFIG       = 'config/${basis}_scan.conf' " >> utils/settings.mk
	echo "CREATE_JOB_CLUSTER      = '${run_choice}' " >> utils/settings.mk
	echo "CREATE_JOB_nCores       = '${ncores}' " >> utils/settings.mk
	echo "CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_${basis}_FAST' " >> utils/settings.mk
	echo "CREATE_JOB_nJobs        = '${njobs}' " >> utils/settings.mk
	echo "CREATE_JOB_chain_length = '${chain_len}' " >> utils/settings.mk
	echo "CREATE_JOB_nJobs        = 'template/${basis}_basis_scan.func' " >> utils/settings.mk

	make create-jobs


##### --- Submit Only --- #####
elif [ "${run_type}" == 3 ] ; then
	echo 'Submit ONLY selected'  

	sed -i "s;job_name_;${job_name};g" Makefile
	sed -i "s;Basis_;${basis};g" Makefile
	sed -i "s;nCores_;${ncores};g" Makefile

	echo 'Do you want to run a cluster job (enter yes for cluster, no for local)?'
	read run_choice
	
	if [ "${run_choice}" == 'yes' ] ; then
		echo 'Please choose a cluster job manager option, Slurm or Torque. Enter "s" for Slurm or "t" for Torque'
		read cluster_choice_
	
		if [ "${cluster_choice_}" == "s" ] ; then		

			echo 'Enter number of nodes'
			read node_no

			echo 'Enter the number of tasks to be assigned per node (or per core)'
			read node_tasks

			echo 'Enter run time for job/s (in format of hh:mm:ss)'
			read run_time
	
			echo "CREATE_JOB_chain_length = '${chain_len}' " >> utils/settings.mk

			echo "SUBMIT_JOB_LIST        = "all.jobs"" >> utils/settings.mk
			echo "SUBMIT_JOB_NAME        = "${job_name}"" >> utils/settings.mk
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


##### --- Merge --- #####
elif [ "${run_type}" == 4 ] ; then

	echo 'Merge jobs selected'
	
	echo 'Do you want to convert your output into hd5f format? Answer "yes" or "no"'
	read to_convert

	echo 'Do you ONLY want to convert your output into hd5f (i.e. not perform merging of jobs)?Answer "yes" or "no"'
	read only_convert

	echo 'Do you want to create a CSV of your merged job output? Answer "yes" or "no"'
	read make_csv

	echo "MERGE_JOB_NAME         = "job_name_"" >> utils/settings.mk
	echo "MERGE_JOB_HEADER       = "header/default.header"" >> utils/settings.mk
	echo "MERGE_JOB_DATASET_NAME = ${MERGE_JOB_NAME}" >> utils/settings.mk
	echo "MERGE_JOB_CONVERT      = "convert_"" >> utils/settings.mk
	echo "MERGE_JOB_CONVERT_ONLY = "just_convert_"" >> utils/settings.mk
	echo "MERGE_JOB_MAKE_CSV     = "to_csv_"" >> utils/settings.mk
	echo "MERGE_JOB_FORMAT       = "table"" >> utils/settings.mk
	echo "MERGE_JOB_COMPRESSION  = "blosc"" >> utils/settings.mk

	make merge-jobs

fi  


