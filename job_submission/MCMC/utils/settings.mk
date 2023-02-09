### --- Makefile Settings --- ###
CREATE_JOB_NAME         = extrastats
CREATE_JOB_CONFIG       = config/hybrid_scan.conf
CREATE_JOB_CLUSTER      = yes
CREATE_JOB_nCores       = 4
CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_hybrid_FAST'
CREATE_JOB_nJobs        = 8
CREATE_JOB_chain_length = 9999999
CREATE_JOB_TEMPLATE        = template/hybrid_basis_scan.func
SUBMIT_JOB_LIST        = all.jobs
SUBMIT_JOB_NAME        = extrastats
SUBMIT_JOB_NODES       = --nodes=4
SUBMIT_JOB_PPN         = --ntasks-per-node=2
SUBMIT_JOB_TIME        = --time=60:00:00
SUBMIT_JOB_TASK        = job_task/job_slurm.sh
SUBMIT_JOB_CLUSTERTYPE = SLURM
