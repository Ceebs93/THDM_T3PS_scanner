#!/bin/bash

# -
echo 'Starting setup for MCMC jobs...'
echo 'Enter 1 to create a job and submit/run it, enter 2 to create ONLY, enter 3 to submit ONLY or enter 4 to merge previous job output.'
read run_type

echo 'Enter job_name'
read job_name

echo 'Enter chosen parameterisation basis: higgs, hhg (higgs hunters guide), generic, phys (physical) or hybrid'
read basis

echo 'Do you want to run a cluster job (enter yes for cluster, no for local)?'
read run_choice

echo 'Enter number of cores to use'
read ncores

echo 'Enter chain length'
read chain_len

echo 'Enter number of nodes

echo 'Enter the number of tasks to be assigned per node (or per core)'
read node_tasks

cp utils/makefile_draft Makefile

if run_type == 1; then
	echo 'Create and submit job selected'		
	sed -i "s;job_name_;${job_name};g" Makefile
	sed -i "s;basis_;${basis};g" Makefile
	sed -i "s;run_choice_;${run_choice};g" Makefile
	sed -i "s;nCores__;${ncores};g" Makefile
	if run_choice == 'yes'; then
		echo ''
	sed -i "s;njobs_;${njobs};g" Makefile
	sed -i "s;chain_length_;${chain_len};g" Makefile
	
	make create-jobs 
	make submit-jobs
  
elif run_type == 2; then
	echo 'Create ONLY selected'  
	sed -i "s;job_name_;${job_name};g" Makefile
	sed -i "s;basis_;${basis};g" Makefile
	sed -i "s;run_choice_;${run_choice};g" Makefile
	sed -i "s;nCores__;${ncores};g" Makefile


elif run_type == 3; then
	echo 'Submit ONLY selected'  
	sed -i "s;job_name_;${job_name};g" Makefile
	sed -i "s;basis_;${basis};g" Makefile
	sed -i "s;nCores__;${ncores};g" Makefile


elif run_type == 4; then
	echo 'Merge jobs selected'
	sed -i "s;job_name_;${job_name};g" Makefile
	sed -i "s;basis_;${basis};g" Makefile
	sed -i "s;nCores__;${ncores};g" Makefile


fi  


