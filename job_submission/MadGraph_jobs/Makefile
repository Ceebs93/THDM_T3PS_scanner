SHELL = /bin/sh

##################################
###### --- Create jobs  --- ######
##################################

CREATE_JOB_NAME         = "bgh_0206"
CREATE_JOB_INPUT_DATA   = "somehowstillmissing.csv"
CREATE_JOB_LOCAL        = "no"
CREATE_JOB_SPLIT_NAME   = "split_csv_name"
CREATE_JOB_PROCESS_NAME = "bg_twh1"
CREATE_JOB_PROCESS      = "Bs\ g\ \>\ Ts\ Ws\ h1"
CREATE_JOB_nCores       = 1
CREATE_JOB_nJobs        = 10
CREATE_JOB_LOOPER       = "csv_looper_NC.sh"
CREATE_JOB_Card_Editor  = "proccard_editor2.py" 
CREATE_JOB_BASECARD     = "basecard_type2.txt"
CREATE_JOB_RIPPER       = "Data_Ripper_typeII_template_new.py"
CREATE_JOB_LIST         = "all.jobs"
CREATE_JOB_MGVAR1_LABEL  = "tb"
CREATE_JOB_MGVAR2_LABEL  = "sinba"  
CREATE_JOB_MGVAR3_LABEL  = "mH"  
CREATE_JOB_MGVAR4_LABEL  = "mA"  


##################################
##### --- Job Submission --- #####
##################################

SUBMIT_JOB_NAME       = ${CREATE_JOB_NAME}
SUBMIT_JOB_LOCAL      = "no" 
SUBMIT_JOB_LIST       = ${CREATE_JOB_LIST}
SUBMIT_JOB_NODES      = "--nodes=1"
SUBMIT_JOB_PPN        = "--ntasks-per-node=4"
SUBMIT_JOB_TIME       = "--time=60:00:00"
SUBMIT_JOB_SCRIPT     = "slurm_submission.sh"
SUBMIT_JOB_CLUSTERTYPE = "SLURM"


#For Torque jobs
#SUBMIT_JOB_LIST        = ${CREATE_JOB_NAME}
#SUBMIT_JOB_TAG         = "job_batch1"
#SUBMIT_JOB_RESOURCES   = "nodes=1:ppn=8,walltime=00:30:00"
#SUBMIT_JOB_TASK        = "job_task/job.sh"
#SUBMIT_JOB_CLUSTERTYPE = "TORQUE"


#########################################
########### --- Merge jobs - ############
#########################################
MERGE_JOB_NAME         = "bgh_0206"
MERGE_JOB_PROCESS      = "bg_twh1"
MERGE_JOB_INPUT_DATA   = "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph_jobs/somehowstillmissing.csv"
MERGE_JOB_XSECT_COL    = "X_sect_bg_h"
MERGE_JOB_MGVAR1_LABEL = "Tb"
MERGE_JOB_OGVAR1_LABEL = "tb"
MERGE_JOB_MGVAR2_LABEL = "sinba"
MERGE_JOB_OGVAR2_LABEL = "sinba"
MERGE_JOB_MGVAR3_LABEL = "mH"
MERGE_JOB_OGVAR3_LABEL = "mH"
MERGE_JOB_MGVAR4_LABEL = "mA"
MERGE_JOB_OGVAR4_LABEL = "mA"
#########################################


DECLARE_ROOT_DIR = ROOT_DIR="/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph_jobs/"

VAR_CREATE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/CREATE_JOB_/' | sed 's/CREATE_JOB_//g' )
EXPORT_CREATE_JOB = $(foreach v,$(VAR_CREATE_JOB),$(v)="$(CREATE_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_SUBMIT_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/SUBMIT_JOB_/' | sed 's/SUBMIT_JOB_//g' )
EXPORT_SUBMIT_JOB = $(foreach v,$(VAR_SUBMIT_JOB),$(v)="$(SUBMIT_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_MERGE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/MERGE_JOB_/' | sed 's/MERGE_JOB_//g' )
EXPORT_MERGE_JOB = $(foreach v,$(VAR_MERGE_JOB),$(v)="$(MERGE_JOB_$(v))") $(DECLARE_ROOT_DIR)

#########################################################

create-jobs:
	@$(EXPORT_CREATE_JOB) ./utils/create-jobs.sh

submit-jobs:
	@$(EXPORT_SUBMIT_JOB) ./utils/submit-jobs.sh

merge-jobs:
	@$(EXPORT_MERGE_JOB) ./utils/merge-jobs.sh

