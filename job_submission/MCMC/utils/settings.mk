### --- Makefile Settings --- ###
CREATE_JOB_NAME         = 2901_hhgtest
CREATE_JOB_CONFIG       = config/hhg_scan.conf
CREATE_JOB_CLUSTER      = yes
CREATE_JOB_nCores       = 1
CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_hhg_FAST'
CREATE_JOB_nJobs        = 1
CREATE_JOB_chain_length = 10000
CREATE_JOB_TEMPLATE        = template/hhg_basis_scan.func
SUBMIT_JOB_LIST        = all.jobs
SUBMIT_JOB_NAME        = 2901_hhgtest
SUBMIT_JOB_NODES       = --nodes=1
SUBMIT_JOB_PPN         = --ntasks-per-node=2
SUBMIT_JOB_TIME        = --time=00:15:00
SUBMIT_JOB_TASK        = job_task/job_slurm.sh
SUBMIT_JOB_CLUSTERTYPE = SLURM
