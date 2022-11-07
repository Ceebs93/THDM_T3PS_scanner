############################
### --- Create jobs  --- ###
############################

# - Specify the following settings for preparing the jobs

CREATE_JOB_NAME         = "0711_test"
CREATE_JOB_CONFIG       = "config/hybrid_scan.conf"
CREATE_JOB_CLUSTER      = "yes"
CREATE_JOB_nCores       = 1
CREATE_JOB_program      = /scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_FAST
CREATE_JOB_nJobs        = 2
CREATE_JOB_chain_length = 9999999
CREATE_JOB_TEMPLATE     = "template/hybrid_basis_scan.func"

##############################
### --- Job Submission --- ###
##############################

# - Specify the following settings for submitting the jobs

#For Local and Slurm jobs
SUBMIT_JOB_LIST        = "all.jobs"
SUBMIT_JOB_NAME        = "0711_test"
SUBMIT_JOB_NODES       = "--nodes=1"
SUBMIT_JOB_PPN         = "--ntasks-per-node=1"
SUBMIT_JOB_TIME        = "--time=10:00:00"
SUBMIT_JOB_TASK        = "job_task/job_slurm.sh"
SUBMIT_JOB_CLUSTERTYPE = "SLURM"

#For Torque jobs
#SUBMIT_JOB_LIST        = "all.jobs"
#SUBMIT_JOB_TAG         = "job_batch1"
#SUBMIT_JOB_RESOURCES   = "nodes=1:ppn=8,walltime=00:30:00"
#SUBMIT_JOB_TASK        = "job_task/job.sh"
#SUBMIT_JOB_CLUSTERTYPE = "TORQUE"

########################
### --- Merge jobs - ###
########################

# - Specify the following settings for merging the jobs

MERGE_JOB_NAME         = "1707_allnarrow"
MERGE_JOB_HEADER       = "header/default.header"
MERGE_JOB_DATASET_NAME = ${MERGE_JOB_NAME}
MERGE_JOB_CONVERT      = "no"
MERGE_JOB_CONVERT_ONLY = "no"
MERGE_JOB_MAKE_CSV     = "yes"
MERGE_JOB_FORMAT       = "table"
MERGE_JOB_COMPRESSION  = "blosc"

#######################

DECLARE_ROOT_DIR = ROOT_DIR="/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MCMC"

VAR_CREATE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/CREATE_JOB_/' | sed 's/CREATE_JOB_//g' )
EXPORT_CREATE_JOB = $(foreach v,$(VAR_CREATE_JOB),$(v)="$(CREATE_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_SUBMIT_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/SUBMIT_JOB_/' | sed 's/SUBMIT_JOB_//g' )
EXPORT_SUBMIT_JOB = $(foreach v,$(VAR_SUBMIT_JOB),$(v)="$(SUBMIT_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_MERGE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/MERGE_JOB_/' | sed 's/MERGE_JOB_//g' )
EXPORT_MERGE_JOB = $(foreach v,$(VAR_MERGE_JOB),$(v)="$(MERGE_JOB_$(v))") $(DECLARE_ROOT_DIR)

################################################################

help:
	@echo -e "\nAvailable commands: "
	@echo -e "-  create-jobs"
	@echo -e "-  submit-jobs"
	@echo -e "-  merge-jobs"

create-jobs:
	 @$(EXPORT_CREATE_JOB) ./utils/create-jobs.sh

submit-jobs:
	 @$(EXPORT_SUBMIT_JOB) ./utils/submit-jobs.sh

merge-jobs:
	 @$(EXPORT_MERGE_JOB) ./utils/merge-jobs.sh
