############################
### --- Create jobs  --- ###
############################

CREATE_JOB_TAG          = "2019_05_09_Diadem2"
CREATE_JOB_CONFIG       = "config/mcmc_scan_diadem.conf"
CREATE_JOB_nCores       = 8
CREATE_JOB_program      = ${THDM_T3PS_SCANNER_DIR}/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_FAST
CREATE_JOB_nJobs        = 200
CREATE_JOB_chain_length = 9999999
CREATE_JOB_TEMPLATE     = "template/mcmc_scan.func"

##############################
### --- Job Submission --- ###
##############################

#SUBMIT_JOB_LIST       = "test.jobs"
#SUBMIT_JOB_LIST       = "batch1.jobs"
#SUBMIT_JOB_LIST       = "batch2.jobs"
SUBMIT_JOB_LIST       = "batch3.jobs"
SUBMIT_JOB_TAG        = "2019_05_09_Diadem2"
SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=60:00:00"
SUBMIT_JOB_TASK       = "job_task/job.sh"

########################
### --- Merge jobs - ###
########################

MERGE_JOB_TAG          = "2019_04_03_Diadem"
#MERGE_JOB_HEADER      = "header/default.header"
MERGE_JOB_HEADER       = "header/default.header"
MERGE_JOB_DATASET_NAME = ${MERGE_JOB_TAG}
MERGE_JOB_CONVERT      = "yes"
MERGE_JOB_CONVERT_ONLY = "yes"
MERGE_JOB_FORMAT       = "table"
MERGE_JOB_COMPRESSION  = "blosc"

#######################

DECLARE_ROOT_DIR = ROOT_DIR="${THDM_T3PS_SCANNER_DIR}/job_submission/MCMC"

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
