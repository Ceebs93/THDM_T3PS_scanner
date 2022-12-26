
binaries : ../bin/ParameterScan_T3PS_with_HB_HS_FAST 
			#  ../bin/test

CMN =  ../lib/EWPO.cpp.o

OBJ_ParameterScan_T3PS_with_HB_HS               = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS.cpp.o
OBJ_ParameterScan_T3PS_with_HB_HS_FAST          = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS_FAST.cpp.o
#OBJ_ParameterScan_T3PS_with_HB_HS_Physical_FAST = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS_Physical_FAST.cpp.o
OBJ_test = ../lib/test.cpp.o

#######

LINKER = g++

####################

# - ParameterScan_T3PS_with_HB_HS
../bin/ParameterScan_T3PS_with_HB_HS: $(OBJ_ParameterScan_T3PS_with_HB_HS)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - ParameterScan_T3PS_with_HB_HS
../bin/ParameterScan_T3PS_with_HB_HS_FAST: $(OBJ_ParameterScan_T3PS_with_HB_HS_FAST)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - ParameterScan_T3PS_with_HB_HS
#../bin/ParameterScan_T3PS_with_HB_HS_Physical_FAST: $(OBJ_ParameterScan_T3PS_with_HB_HS_Physical_FAST)
#	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - test
#../bin/test: $(OBJ_test)
#	$(CC) $^ $(CLFLAGS_) -o $@

####################
## -- Linking -- ###
####################

LIBDIRLINK = ${THDM_T3PS_SCANNER_DIR}/links/lib
INCDIRLINK = ${THDM_T3PS_SCANNER_DIR}/links/inc

#################
### -- C++ -- ###
#################

# - Compiler
CC = g++

# - CFLAGS
CCFLAGS = -I../inc/ -I$(INCDIRLINK) $(2HDMC_CFLAGS) -MMD -MF .depend_cpp -fpermissive
2HDMC_CFLAGS = -I/local/software/gsl/2.6/include -DHiggsBounds

# - LFLAGS
#CLFLAGS_2HDMC_WITH_HB             = -L$(LIBDIRLINK)
#CLFLAGS_2HDMC_WITH_HB             = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS} 
CLFLAGS_2HDMC_WITH_HB_HS          = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS}
#CLFLAGS_2HDMC_STANDALONE          = -L$(LIBDIRLINK) ${2HDMC_STANDALONE_LFLAGS}
#CLFLAGS_2HDMC_WITH_HB_HS_SuperIso = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS}

############

# - 2HDMC flags
2HDMC_WITH_HB_HS_LFLAGS = -L$(LIBDIRLINK) -l2HDMC_HBHS -L/local/software/gsl/2.6/lib -lgsl -lgslcblas -lgfortran
#2HDMC_STANDALONE_LFLAGS = -L$(LIBDIRLINK) -l2HDMC -L/local/software/gsl/2.6/lib -lgsl -lgslcblas -lgfortran

# - HiggsBounds/HiggsSignals flags
HB_HS_LFLAGS = ${HB_LFLAGS} ${HS_LFLAGS}
HB_LFLAGS = -L$(LIBDIRLINK) -lHB-5.10.1
HS_LFLAGS = -L$(LIBDIRLINK) -lHS-2.6.2

#####################################################
### - General compilation rule for .cpp sources - ###
#####################################################

../lib/ParameterScan_T3PS_with_HB_HS_FAST.cpp.o : ParameterScan_T3PS_with_HB_HS.cpp
	$(CC) $(CCFLAGS) -c $< -DFAST -o $@

#../lib/ParameterScan_T3PS_with_HB_HS_Physical_FAST.cpp.o : ParameterScan_T3PS_with_HB_HS_Physical.cpp
#	$(CC) $(CCFLAGS) -c $< -o $@
#	$(CC) $(CCFLAGS) -c $< -DFAST -o $@

../lib/%.cpp.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

include .depend_cpp

clean:
	rm ../lib/*
