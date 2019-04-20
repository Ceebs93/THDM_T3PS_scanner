############################
### --- Create jobs  --- ###
############################

# - Beta - 89m - #
#CREATE_JOB_TAG        = "Beta_89m_pts"
#CREATE_JOB_INPUT_DATA = "data/Beta_89m_pts/allowed_pts_only_all_tab_separeted.dat"
#CREATE_JOB_CONFIG     = "config/job_hybrid_cba_A_and_H.conf"
#CREATE_JOB_nCores     = 8
#CREATE_JOB_nJobs      = 100
#CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
#CREATE_JOB_TASK       = "job_task/job.sh"

# - Gamma  - #
#CREATE_JOB_TAG        = "Gamma"
#CREATE_JOB_INPUT_DATA = "data/Gamma/sparsed.dat"
#CREATE_JOB_CONFIG     = "config/job_hybrid_cba_A_and_H_extended.conf"
#CREATE_JOB_nCores     = 8
#CREATE_JOB_nJobs      = 100
#CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
#CREATE_JOB_TASK       = "job_task/job.sh"

# - Delta  - #
#CREATE_JOB_TAG        = "Delta"
#CREATE_JOB_INPUT_DATA = "data/Delta/Delta_MCMC_allowed_only_merged_undersampled_by_factor_8.dat"
#CREATE_JOB_CONFIG     = "config/job_hybrid_cba_A_and_H_extended.conf"
#CREATE_JOB_nCores     = 8
#CREATE_JOB_nJobs      = 100
#CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
#CREATE_JOB_TASK       = "job_task/job.sh"

# - Epsilon  - #
#CREATE_JOB_TAG        = "Epsilon"
#CREATE_JOB_INPUT_DATA = "data/Epsilon/Epsilon_MCMC_allowed_only_tab_sep.dat"
#CREATE_JOB_CONFIG     = "config/job_hybrid_cba_A_and_H_extended.conf"
#CREATE_JOB_nCores     = 8
#CREATE_JOB_nJobs      = 100
#CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
#CREATE_JOB_TASK       = "job_task/job.sh"

# - Falcon  - #
CREATE_JOB_TAG        = "Falcon"
CREATE_JOB_INPUT_DATA = "data/Falcon/Falcon_MCMC_allowed_only_merged_no_xsec_tab_sep.dat"
CREATE_JOB_CONFIG     = "config/job_hybrid_cba_A_and_H_extended.conf"
CREATE_JOB_nCores     = 8
CREATE_JOB_nJobs      = 100
CREATE_JOB_TEMPLATE   = "template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha"
CREATE_JOB_TASK       = "job_task/job.sh"


##############################
### --- Job Submission --- ###
##############################

#SUBMIT_JOB_TAG        = "Beta_89m_pts"
##SUBMIT_JOB_LIST       = "test_submission.jobs"
#SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=01:30:00"
#SUBMIT_JOB_TASK       = "job_task/job.sh"

#SUBMIT_JOB_TAG        = "Gamma"
##SUBMIT_JOB_LIST       = "test_submission.job"
#SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=12:00:00"
##SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=1,walltime=00:10:00"
#SUBMIT_JOB_TASK       = "job_task/job.sh"

#SUBMIT_JOB_TAG        = "Delta"
##SUBMIT_JOB_LIST       = "test_submission.job"
#SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=6:00:00"
##SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=1,walltime=00:10:00"
#SUBMIT_JOB_TASK       = "job_task/job.sh"

# - Epsilon
#SUBMIT_JOB_TAG        = "Epsilon"
##SUBMIT_JOB_LIST       = "test_submission.job"
#SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=5:00:00"
##SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=1,walltime=00:10:00"
#SUBMIT_JOB_TASK       = "job_task/job.sh"


# - Falcon - #
SUBMIT_JOB_TAG        = "Falcon"
#SUBMIT_JOB_LIST       = "test_submission.job"
SUBMIT_JOB_LIST       = "all.jobs"
#SUBMIT_JOB_LIST       = "try.job"
SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=8,walltime=6:00:00"
#SUBMIT_JOB_RESOURCES  = "nodes=1:ppn=1,walltime=00:10:00"
SUBMIT_JOB_TASK       = "job_task/job.sh"


########################
### --- Merge jobs - ###
########################

#MERGE_JOB_TAG         = "Beta_89m_pts"
#MERGE_JOB_HEADER      = "header/default.header"

# - Epsilon
#MERGE_JOB_TAG         = "Epsilon"
#MERGE_JOB_HEADER      = "header/default_extended.header"

# - Falcon
MERGE_JOB_TAG         = "Falcon"
MERGE_JOB_HEADER      = "header/default_extended.header"

#######################

DECLARE_ROOT_DIR = ROOT_DIR="${THDM_T3PS_SCANNER_DIR}/sushi-submission"

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
