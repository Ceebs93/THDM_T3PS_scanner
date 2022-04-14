# Makefile for 2HDMC 

# Choose your C++ compiler here (in general g++ on Linux systems):
CC = g++

#Optimisation level, eg: -O3
OPT=-O3

#OR debug level: -g(n=1,2,3)
DEBUG= 

VPATH=src 

# *x* Note to self: Added gsl/1.15/include here
CFLAGS= -Wall $(DEBUG) $(OPT) -I/local/software/gsl/1.15/include -fpermissive
CFLAGS_HBHS= $(CFLAGS) -DHiggsBounds -fpermissive

OBJDIR=lib
LIBDIR=$(OBJDIR)

OBJDIR_HBHS=lib_HBHS
LIBDIR_HBHS=$(OBJDIR_HBHS)


# - HiggsBounds and HiggsSignals paths
HIGGSBOUNDS_PATH=${THDM_T3PS_SCANNER_DIR}/packages/HiggsBounds-5.3.2beta
HIGGSSIGNALS_PATH=${THDM_T3PS_SCANNER_DIR}/packages/HiggsSignals-2.2.3beta


SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES_HBHS=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp HBHS.cpp
OBJECTS_HBHS=$(SOURCES_HBHS:.cpp=.o)


LDFLAGS = -lgsl -lgslcblas
LDFLAGS += -L$(LIBDIR) -L$(HIGGSBOUNDS_PATH) -L$(HIGGSSIGNALS_PATH) -L/local/software/gsl/1.15/lib -l2HDMC -lgsl -lgslcblas -lm
LDFLAGS += -lgfortran

LDFLAGS_HBHS = $(LDFLAGS) -lHS -lHB -lgfortran -L$(LIBDIR) -l2HDMC

LIB=lib2HDMC.a
LIB_HBHS=lib2HDMC_HBHS.a

LIBS=lib/Constraints.o lib/DecayTable.o lib/SM.o lib/THDM.o lib/runTHDM.o
LIBS_HS=lib_HS/Constraints.o lib_HS/DecayTable.o lib_HS/SM.o lib_HS/THDM.o lib_HS/runTHDM.o
PROG=CalcPhys CalcGen CalcHiggs CalcHybrid CalcHMSSM CalcMSSM CalcInert CalcLH Demo

# To use HiggsBounds/HiggsSignals for Higgs constraints, put the corresponding libraries
# in the "lib" subdirectory and uncomment the following lines.
# Note that both programs: HiggsBounds (> 4.2.0) and HiggsSignals (> 1.4.0) must be present


.PHONY: lib clean distclean

all: lib lib_HBHS $(PROG)

$(OBJDIR)/%.o : %.cpp %.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR_HBHS)/%.o : %.cpp %.h
	$(CC) $(CFLAGS_HBHS) -c $< -o $@


lib: $(addprefix $(LIBDIR)/, $(LIB))
lib_HBHS: $(addprefix $(LIBDIR_HBHS)/, $(LIB_HBHS))

$(addprefix $(LIBDIR)/, $(LIB)): $(addprefix $(OBJDIR)/, $(OBJECTS))
	echo "Making library $@"
	ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS))

$(addprefix $(LIBDIR_HBHS)/, $(LIB_HBHS)): $(addprefix $(OBJDIR_HBHS)/, $(OBJECTS_HBHS))
	echo "Making library $@"
	ar rcs $@ $(addprefix $(OBJDIR_HBHS)/, $(OBJECTS_HBHS))

%: %.cpp $(addprefix $(LIBDIR)/, $(LIB))
	$(CC) $< -Isrc $(CFLAGS) $(LDFLAGS) -o $@

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
