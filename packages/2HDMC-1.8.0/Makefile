# Makefile used for iridis4 installation of 2HDMC-1.8.0

# Choose your C++ compiler here (in general g++ on Linux systems):
CC = g++
LDFLAGS= -I/local/software/gsl/2.6/include -lgsl -lgslcblas

#Optimisation level, eg: -O3
OPT=-O3
#OR debug level: -g(n=1,2,3)
DEBUG=

VPATH=src

CFLAGS= -std=c++11 -Wall $(DEBUG) $(OPT) -I/local/software/gsl/2.6/include -lgsl -lgslcblas
CFLAGS_HBHS= $(CFLAGS) -DHiggsBounds -fpermissive

OBJDIR=lib
OBJDIR_HBHS=lib_HBHS

LIBDIR=$(OBJDIR)
LIBDIR_HBHS=$(OBJDIR_HBHS)
SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp
SOURCES_HBHS= THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp HBHS.cpp

OBJECTS=$(SOURCES:.cpp=.o)
LIB=lib2HDMC.a
LIB_HBHS=lib2HDMC_HBHS.a
LDFLAGS+=-L$(LIBDIR) -l2HDMC  -L/local/software/gsl/2.6/lib -lgsl -lgslcblas -lm

LIBS=
PROG=CalcPhys CalcGen CalcHiggs CalcHybrid CalcHMSSM CalcMSSM CalcInert CalcLH Demo #check_point check_point_2hdmc
INCLUDE=

# To use HiggsBounds/HiggsSignals for Higgs constraints, set both of the
# following path variables to the corresponding build directories.
# Requires HiggsBounds>=5.7.0 and HiggsSignals>=2.4.0
HiggsBounds_DIR=/${THDM_T3PS_SCANNER_DIR}/packages/HiggsBounds-5.10.1/build
HiggsSignals_DIR=/${THDM_T3PS_SCANNER_DIR}/packages/HiggsSignals-2.6.2/build

ifdef HiggsBounds_DIR
ifdef HiggsSignals_DIR
CFLAGS+=-DHiggsBounds
LDFLAGS+=-L$(HiggsBounds_DIR)/lib -L$(HiggsSignals_DIR)/lib -lHS -lHB -lgfortran
LDFLAGS_HBHS=-L$(LDFLAGS) -lHS -HB -lgfortran -L$(LIBDIR) -l2HDMC
INCLUDE+=-I$(HiggsBounds_DIR)/../include -I$(HiggsSignals_DIR)/../include
SOURCES+=HBHS.cpp
endif
endif

.PHONY: lib clean distclean

all: lib lib_HBHS $(PROG)

$(OBJDIR)/%.o : %.cpp %.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@
$(OBJDIR_HBHS)/%.o : %.cpp %.h
	$(CC) $(CFLAGS_HBHS) -c $< -o $@

lib: $(addprefix $(LIBDIR)/, $(LIB))
lib_HBHS: $(addprefix $(LIBDIR_HBHS)/, $(LIB_HBHS))

$(addprefix $(LIBDIR)/, $(LIB)): $(addprefix $(OBJDIR)/, $(OBJECTS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS))

$(addprefix $(LIBDIR_HBHS)/, $(LIB_HBHS)): $(addprefix $(OBJDIR_HBHS)/, $(OBJECTS_HBHS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR_HBHS)/, $(OBJECTS_HBHS))

%: %.cpp $(addprefix $(LIBDIR)/, $(LIB))
	@ echo $(CC) $< -Isrc $(CFLAGS) $(LDFLAGS) $(addprefix $(LIBDIR)/, $(LIBS)) -o $@
	@ $(CC) $< -Isrc $(CFLAGS) $(LDFLAGS) $(addprefix $(LIBDIR)/, $(LIBS)) -o $@

clean:
	@ echo "Cleaning library"
	@ rm -f $(addprefix $(OBJDIR)/, *.o)
	@ rm -f $(addprefix $(OBJDIR_HBHS)/, *.o) 
	@ rm -f $(addprefix $(LIBDIR)/, $(LIB))
	@ rm -f $(addprefix $(LIBDIR_HBHS)/, $(LIB))

distclean:
	@ make -s clean
	@ echo "Cleaning programs"
	@ rm -f $(PROG)

