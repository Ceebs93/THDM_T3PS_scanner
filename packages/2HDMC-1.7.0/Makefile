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

OBJDIR=lib
LIBDIR=$(OBJDIR)
SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp HBHS.cpp
OBJECTS=$(SOURCES:.cpp=.o)

LIB=lib2HDMC.a
LIBSTANDALONE=lib2HDMCstandalone.a
# *x* Note to self: Added gsl/1.15/lib here
#
HIGGSBOUNDS_PATH=${THDM_T3PS_SCANNER_DIR}/packages/HiggsBounds-5.3.2beta
HIGGSSIGNALS_PATH=${THDM_T3PS_SCANNER_DIR}/packages/HiggsSignals-2.2.3beta

LDFLAGS=-lgsl -lgslcblas
LDFLAGS+=-L$(LIBDIR) -L$(HIGGSBOUNDS_PATH) -L$(HIGGSSIGNALS_PATH) -L/local/software/gsl/1.15/lib -l2HDMC -lgsl -lgslcblas -lm
LDFLAGS+=-lHS -lHB -lgfortran

#LDFLAGS+=-L$(LIBDIR) -l2HDMC
LIBS=
PROG=CalcPhys CalcGen CalcHiggs CalcHybrid CalcHMSSM CalcMSSM CalcInert CalcLH Demo

# To use HiggsBounds/HiggsSignals for Higgs constraints, put the corresponding libraries
# in the "lib" subdirectory and uncomment the following lines.
# Note that both programs: HiggsBounds (> 4.2.0) and HiggsSignals (> 1.4.0) must be present
CFLAGS+=-DHiggsBounds
LDFLAGS+=-L$(LIBDIR) -lHS-2.2.3beta -lHB-4.3.1 -lgfortran
#LDFLAGS+=-L$(LIBDIR) -lHS -lHB-4.3.1 -lgfortran
#LDFLAGS+=-L$(LIBDIR) -lHS -lHB-4.2.1 -lgfortran
SOURCES+=HBHS.cpp

.PHONY: lib clean distclean

all: lib $(PROG)

$(OBJDIR)/%.o : %.cpp %.h
	$(CC) $(CFLAGS) -c $< -o $@

lib: $(addprefix $(LIBDIR)/, $(LIB))
libstandalone: $(addprefix $(LIBDIR)/, $(LIBSTANDALONE))


$(addprefix $(LIBDIR)/, $(LIB)): $(addprefix $(OBJDIR)/, $(OBJECTS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS))

$(addprefix $(LIBDIR)/, $(LIBSTANDALONE)): $(addprefix $(OBJDIR)/, $(OBJECTS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS))

%: %.cpp $(addprefix $(LIBDIR)/, $(LIB))
	$(CC) $< -Isrc $(CFLAGS) $(LDFLAGS) $(addprefix $(LIBDIR)/, $(LIBS)) -o $@

clean:
	@ echo "Cleaning library"
	@ rm -f $(addprefix $(OBJDIR)/, *.o)
	@ rm -f $(addprefix $(LIBDIR)/, $(LIB))

distclean:
	@ make -s clean
	@ echo "Cleaning programs"
	@ rm -f $(PROG)
