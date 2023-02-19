### --- Makefile Settings --- ###
CREATE_JOB_NAME         = Test_Job
CREATE_JOB_CONFIG       = config/generic_scan.conf
CREATE_JOB_CLUSTER      = yes
CREATE_JOB_nCores       = 1
CREATE_JOB_program      = '/scratch/cb27g11/THDM_T3PS_scanner/ParameterPointProcessor/bin/ParameterScan_T3PS_with_HB_HS_generic_FAST'
CREATE_JOB_nJobs        = 2
CREATE_JOB_chain_length = 100000
CREATE_JOB_TEMPLATE        = template/generic_basis_scan.func
