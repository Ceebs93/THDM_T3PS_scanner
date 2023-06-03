# - Test
#run_TASK     = "task_ParamScan_Hybrid_modified_mA_mHc_EWPO_survey.sh"
run_TASK     = "task_ParamScan_Hybrid_modified_mA_mHc.sh"
run_TAG      = "test123"
run_CONFIG   = "ParamSpace_Hybrid_modified_mA_mHc_EWPO_survey.config"
run_WRITELHA = 0
run_ytype    = 2

VAR_RUN    = $(shell echo '$(.VARIABLES)' |  awk -v RS=' ' '/run_/' | sed 's/run_//g' )
EXPORT_RUN = $(foreach v,$(VAR_RUN),$(v)="$(run_$(v))")
###################################################################################


all : build

yt = 2
cba =   0.4
tb  =   5.0
mh  = 125.0
mH  = 500.0
mHc = 500.0
mA  = 200.0
Z7  = 1.5

#testfile = testfile
#testfile = test_still_working
#testfile = test_supposed_to_trigger_bug
#testfile = test_bug_test
#testfile = test_sm_limit
testfile = test_EWPO_T

test :
#	./bin/ParameterScan_T3PS $(yt) $(Z7) $(mH) $(mHc) $(mA) $(cba) $(tb)
#	./bin/ParameterScan_T3PS $(testfile)
#	./bin/ParameterScan_T3PS_with_HB $(yt) $(Z7) $(mH) $(mHc) $(mA) $(cba) $(tb)
#	./bin/ParameterScan_T3PS_with_HB $(testfile)
#	./bin/ParameterScan_T3PS_with_HB_HS $(yt) $(Z7) $(mH) $(mHc) $(mA) $(cba) $(tb)
#  ./bin/ParameterScan_T3PS_with_HB_HS $(testfile)
#	./bin/ParameterScan_T3PS_with_HB_HS_Physical_FAST physical.test
	./bin/ParameterScan_T3PS_with_HB_HS_generic_FAST $(testfile)
	./bin/ParameterScan_T3PS_with_HB_HS_hybrid_FAST $(testfile)
	./bin/ParameterScan_T3PS_with_HB_HS_phys_FAST $(testfile)
	./bin/ParameterScan_T3PS_with_extra_data

test-strings : .PHONY
	@echo "This is VAR_RUN:"
	@echo "$(VAR_RUN)"
	@echo "This is EXPORT_JOB:"
	@echo "$(EXPORT_JOB)"
	@echo "This is EXPORT_FORM_DAT:"
	@echo "$(EXPORT_FORM_DAT)"
	@echo "This is SrV_format_data_eval:"
	@echo "$(SrV_format_data_eval)"

run : build
	@$(EXPORT_RUN) ./tasks/$(run_TASK); 

new : 
	@TAG=$(TAG); ./scripts/createWD.sh

build : build_allbinaries

build_allbinaries : 
	@ cd src; make binaries;

prep_xsec_sub :
	./scripts/prep_xsec_sub.sh $(xsec_job_tag) $(xsec_job_split)


clean :
	rm -f ./lib/*.o
	find ./bin/ -type f -not -name 'dummy' | xargs rm -f
	rm -f ./src/.depend_cpp
	touch ./src/.depend_cpp


.PHONY: figures
