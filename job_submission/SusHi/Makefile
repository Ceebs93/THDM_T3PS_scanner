############################
### --- Create jobs  --- ###
############################

CREATE_JOB_NAME        = "Name_for_job"
CREATE_JOB_INPUT_DATA  = "/scratch/cb27g11/Looping_Bash/missing.csv"
CREATE_JOB_SPLIT_NAME  = "split_csv_name"
CREATE_JOB_PROCESS     = "bq_tqh1"
#CREATE_JOB_CONFIG      = "config/default.conf"
CREATE_JOB_nCores      = 8
CREATE_JOB_nJobs       = 3
CREATE_JOB_LOOPER      = "csv_looper_draft.sh"
CREATE_JOB_Card_Editor = "inputcard_editor1.py" 
CREATE_JOB_BASECARD    = "basecard_type1.txt"
CREATE_JOB_TEMPLATE    = "slurm_submission.sh"
CREATE_JOB_RIPPER      = "Data_Ripper_typeI_template.py"


##############################
### --- Job Submission --- ###
##############################

SUBMIT_JOB_NAME       = "Name_for_job"
SUBMIT_JOB_LIST       = "all.jobs"
SUBMIT_JOB_NODES      = "--nodes=1"
SUBMIT_JOB_PPN        = "--ntasks-per-node=8"
SUBMIT_JOB_TIME       = "--time=00:05:00"
SUBMIT_JOB_SCRIPT     = "slurm_submission.sh"


########################
### --- Merge jobs - ###
########################

# - 2019_04_03_Diadem - #
MERGE_JOB_NAME         = "2019_04_03_Diadem"
MERGE_JOB_HEADER      = "header/default.header"

#######################

DECLARE_ROOT_DIR = ROOT_DIR="${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/"

VAR_CREATE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/CREATE_JOB_/' | sed 's/CREATE_JOB_//g' )
EXPORT_CREATE_JOB = $(foreach v,$(VAR_CREATE_JOB),$(v)="$(CREATE_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_SUBMIT_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/SUBMIT_JOB_/' | sed 's/SUBMIT_JOB_//g' )
EXPORT_SUBMIT_JOB = $(foreach v,$(VAR_SUBMIT_JOB),$(v)="$(SUBMIT_JOB_$(v))") $(DECLARE_ROOT_DIR)

VAR_MERGE_JOB    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/MERGE_JOB_/' | sed 's/MERGE_JOB_//g' )
EXPORT_MERGE_JOB = $(foreach v,$(VAR_MERGE_JOB),$(v)="$(MERGE_JOB_$(v))") $(DECLARE_ROOT_DIR)

#########################################################

create-jobs:
	 @$(EXPORT_CREATE_JOB) ./utils/new_create-jobs.sh

submit-jobs:
	 @$(EXPORT_SUBMIT_JOB) ./utils/new_submit-jobs.sh

merge-jobs:
	 @$(EXPORT_MERGE_JOB) ./utils/merge-jobs.sh
