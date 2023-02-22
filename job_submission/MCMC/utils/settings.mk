### --- Makefile Settings --- ###
CREATE_JOB_NAME         = demo_job
CREATE_JOB_CONFIG       = config/phys_scan.conf
CREATE_JOB_CLUSTER      = yes
CREATE_JOB_nCores       = 1
CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_phys_FAST'
CREATE_JOB_nJobs        = 8
CREATE_JOB_chain_length = 999999
CREATE_JOB_TEMPLATE     = template/phys_basis_scan.func
CREATE_JOB_Y            = 1
SUBMIT_JOB_LIST        = all.jobs
SUBMIT_JOB_NAME        = demo_job
SUBMIT_JOB_NODES       = --nodes=2
SUBMIT_JOB_PPN         = --ntasks-per-node=4
SUBMIT_JOB_TIME        = --time=20:00:00
SUBMIT_JOB_TASK        = job_task/job_slurm.sh
SUBMIT_JOB_CLUSTERTYPE = SLURM
