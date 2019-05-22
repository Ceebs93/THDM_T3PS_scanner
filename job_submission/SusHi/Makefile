############################
### --- Create jobs  --- ###
############################

# - 2019_04_03_Diadem  - #
CREATE_JOB_TAG        = "2019_04_03_Diadem"
CREATE_JOB_INPUT_DATA = "/scratch/de3u14/2HDM/Magellan_THDM_T3PS_scanner/job_submission/MCMC/jobs/2019_04_03_Diadem/2019_04_03_Diadem_allowed_only_tabsep_noheader.dat"
CREATE_JOB_CONFIG     = "config/default.conf"
CREATE_JOB_nCores     = 8
CREATE_JOB_nJobs      = 200
CREATE_JOB_program    = "${THDM_T3PS_SCANNER_DIR}/job_submission/SusHi/scripts/calc_A_and_H_xsec.sh"
CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
CREATE_JOB_TASK       = "job_task/job.sh"


##############################
### --- Job Submission --- ###
##############################

# - 2019_04_03_Diadem - #
SUBMIT_JOB_TAG        = "2019_04_03_Diadem"
#SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_LIST       = "test.jobs"
SUBMIT_JOB_LIST       = "rest.jobs"
SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=16:00:00"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=1,walltime=00:10:00"
SUBMIT_JOB_TASK       = "job_task/job.sh"


########################
### --- Merge jobs - ###
########################

# - 2019_04_03_Diadem - #
MERGE_JOB_TAG         = "2019_04_03_Diadem"
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
	 @$(EXPORT_CREATE_JOB) ./utils/create-jobs.sh

submit-jobs:
	 @$(EXPORT_SUBMIT_JOB) ./utils/submit-jobs.sh

merge-jobs:
	 @$(EXPORT_MERGE_JOB) ./utils/merge-jobs.sh
