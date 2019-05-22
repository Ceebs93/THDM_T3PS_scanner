
binaries : ../bin/ParameterScan_T3PS_with_HB_HS_FAST \
			  ../bin/ParameterScan_T3PS_with_HB_HS_Physical_FAST \
			  ../bin/test

CMN =  ../lib/EWPO.cpp.o

OBJ_ParameterScan_T3PS                          = ../lib/ParameterScan_T3PS.cpp.o
OBJ_ParameterScan_T3PS_with_HB                  = ../lib/ParameterScan_T3PS_with_HB.cpp.o
OBJ_ParameterScan_T3PS_with_HB_HS               = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS.cpp.o
OBJ_ParameterScan_T3PS_with_HB_HS_FAST          = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS_FAST.cpp.o
OBJ_ParameterScan_T3PS_with_HB_HS_SuperIso_FAST = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS_SuperIso_FAST.cpp.o
OBJ_ParameterScan_T3PS_with_HB_HS_Physical_FAST = $(CMN) ../lib/ParameterScan_T3PS_with_HB_HS_Physical_FAST.cpp.o
OBJ_test = ../lib/test.cpp.o

#######

LINKER = g++

####################

# - ParameterScan_T3PS
../bin/ParameterScan_T3PS: $(OBJ_ParameterScan_T3PS)
	$(CC) $^ $(CLFLAGS_2HDMC_STANDALONE) -o $@

# - ParameterScan_T3PS_with_HB
../bin/ParameterScan_T3PS_with_HB: $(OBJ_ParameterScan_T3PS_with_HB)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB) -o $@

# - ParameterScan_T3PS_with_HB_HS
../bin/ParameterScan_T3PS_with_HB_HS: $(OBJ_ParameterScan_T3PS_with_HB_HS)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - ParameterScan_T3PS_with_HB_HS
../bin/ParameterScan_T3PS_with_HB_HS_FAST: $(OBJ_ParameterScan_T3PS_with_HB_HS_FAST)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - ParameterScan_T3PS_with_HB_HS
../bin/ParameterScan_T3PS_with_HB_HS_Physical_FAST: $(OBJ_ParameterScan_T3PS_with_HB_HS_Physical_FAST)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS) -o $@

# - ParameterScan_T3PS_with_HB_HS_SuperIso
../bin/ParameterScan_T3PS_with_HB_HS_SuperIso_FAST: $(OBJ_ParameterScan_T3PS_with_HB_HS_SuperIso_FAST)
	$(CC) $^ $(CLFLAGS_2HDMC_WITH_HB_HS_SuperIso) -o $@


# - test
../bin/test: $(OBJ_test)
	$(CC) $^ $(CLFLAGS_) -o $@

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
2HDMC_CFLAGS = -I/local/software/gsl/1.15/include -DHiggsBounds

# - LFLAGS
CLFLAGS_2HDMC_WITH_HB             = -L$(LIBDIRLINK)
CLFLAGS_2HDMC_WITH_HB             = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS} 
CLFLAGS_2HDMC_WITH_HB_HS          = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS}
CLFLAGS_2HDMC_STANDALONE          = -L$(LIBDIRLINK) ${2HDMC_STANDALONE_LFLAGS}
CLFLAGS_2HDMC_WITH_HB_HS_SuperIso = -L$(LIBDIRLINK) ${2HDMC_WITH_HB_HS_LFLAGS} ${HS_LFLAGS} ${HB_LFLAGS}
#CLFLAGS_2HDMCSTANDALONE = -L$(LIBDIRLINK) -l2HDMC -L/local/software/gsl/1.15/lib -lgsl -lgslcblas -lgfortran

############

# - 2HDMC flags
2HDMC_WITH_HB_HS_LFLAGS = -L$(LIBDIRLINK) -l2HDMC_HBHS -L/local/software/gsl/1.15/lib -lgsl -lgslcblas -lgfortran
2HDMC_STANDALONE_LFLAGS = -L$(LIBDIRLINK) -l2HDMC -L/local/software/gsl/1.15/lib -lgsl -lgslcblas -lgfortran

# - HiggsBounds/HiggsSignals flags
HB_HS_LFLAGS = ${HB_LFLAGS} ${HS_LFLAGS}
HB_LFLAGS = -L$(LIBDIRLINK) -lHB-5.3.2beta
HS_LFLAGS = -L$(LIBDIRLINK) -lHS-2.2.3beta

#HB_HS_LFLAGS = -L$(LIBDIRLINK) -lHB -lHS
#HB_HS_LFLAGS = -L$(LIBDIRLINK) -lHB-4.2.1 -lHS

# - SuperIso flags
SuperIso_LFLAGS = -L$(LIBDIRLINK) -lisospin


#####################################################
### - General compilation rule for .cpp sources - ###
#####################################################

../lib/ParameterScan_T3PS_with_HB_HS_FAST.cpp.o : ParameterScan_T3PS_with_HB_HS.cpp
	$(CC) $(CCFLAGS) -c $< -DFAST -o $@

../lib/ParameterScan_T3PS_with_HB_HS_SuperIso_FAST.cpp.o : ParameterScan_T3PS_with_HB_HS_SuperIso.cpp
	$(CC) $(CCFLAGS) -c $< -DFAST -o $@

../lib/ParameterScan_T3PS_with_HB_HS_Physical_FAST.cpp.o : ParameterScan_T3PS_with_HB_HS_Physical.cpp
	$(CC) $(CCFLAGS) -c $< -o $@
#	$(CC) $(CCFLAGS) -c $< -DFAST -o $@

../lib/%.cpp.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

include .depend_cpp

clean:
	rm ../lib/*
