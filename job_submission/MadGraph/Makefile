############################
### --- Create jobs  --- ###
############################

CREATE_JOB_NAME         = "new_constraints_cont"
CREATE_JOB_INPUT_DATA   = "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/new_constraints_cont.csv"
#CREATE_JOB_CONVERT_DATA = "yes"
CREATE_JOB_LOCAL        = "no"
CREATE_JOB_SPLIT_NAME   = "split_csv_name"
CREATE_JOB_PROCESS_NAME = "bg_twh1"
CREATE_JOB_PROC_IN      = "Bs\ g"
CREATE_JOB_PROC_OUT     = "Ts\ Ws\ h1"
CREATE_JOB_nCores       = 1
CREATE_JOB_nJobs        = 8
CREATE_JOB_LOOPER       = "csv_looper_extra.sh"
CREATE_JOB_Card_Editor  = "inputcard_editor2.py" 
CREATE_JOB_BASECARD     = "basecard_type2.txt"
CREATE_JOB_TEMPLATE     = "slurm_submission.sh"
CREATE_JOB_RIPPER       = "Data_Ripper_typeII_template.py"
CREATE_JOB_LIST         = "all.jobs"


##############################
### --- Job Submission --- ###
##############################

SUBMIT_JOB_NAME       = "new_constraints_cont"
SUBMIT_JOB_LOCAL      = "no" 
SUBMIT_JOB_LIST       = "all.jobs"
SUBMIT_JOB_NODES      = "--nodes=1"
SUBMIT_JOB_PPN        = "--ntasks-per-node=8"
SUBMIT_JOB_TIME       = "--time=30:00:00"
SUBMIT_JOB_SCRIPT     = "slurm_submission.sh"
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

# - 2019_04_03_Diadem - #
MERGE_JOB_NAME        = "test_job"
MERGE_JOB_PROCESS     = "bq_tqh2"
MERGE_JOB_INPUT_DATA  = "/scratch/cb27g11/Looping_Bash/missing.csv"
#######################


DECLARE_ROOT_DIR = ROOT_DIR="/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/"

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

