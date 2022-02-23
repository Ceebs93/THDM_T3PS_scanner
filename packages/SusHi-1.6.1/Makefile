#
# SusHi Makefile
#
# Before running make, type "./configure"
# and adjust the paths in Makefile below.
#
# possible options: make [option]
#	predef=FH     : link to FeynHiggs
#       predef=2HDMC  : link to 2HDMC
#       predef=HB     : link to HiggsBounds+FeynHiggs
#	predef=HS     : link to HiggsSignals+FeynHiggs
#	clean      : delete all executables
# The compiler specified in compilerissues - written
# by the ./configure script - is used.
#
SHELL=/bin/bash
# Specify the path to a compiled FeynHiggs version (for "make predef=FH, HB or HS"):
FHPATH = ../FeynHiggs-2.11.2
FHVERSION = 2.11.2
# Specify the path to a compiled HiggsBound version (for "make predef=HB or HS"):
#HBPATH = ../HiggsBounds-4.1.3
HBPATH = ${THDM_T3PS_SCANNER_DIR}/packages/HiggsBounds-5.3.2beta
HBVERSION = 5.3.2beta
# Specify the path to a compiled HiggsSignals version (for "make predef=HS"):
HSPATH = ${THDM_T3PS_SCANNER_DIR}/packages/HiggsSignals-2.2.3beta
HSVERSION = 2.2.3beta
# Specify the path to a compiled 2HDMC version (for "make predef=2HDMC"):
2HDMCPATH = ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.7.0/
2HDMCVERSION = 1.7.0
# Specify the path to the compiled LHAPDF libraries (might be found automatically):
LHAPATH = ${THDM_T3PS_SCANNER_DIR}/packages/LHAPDF-6.2.1_build/lib
# Specify the path to standard C++ libraries (libstdc) 
# (only needed when linked to 2HDMC, might be found automatically):
STDLIB = -L/local/software/gsl/2.6/lib
# =============================================================================
include compilerissues
-include compilerissues.local

#Splitting the version number of FeynHiggs 1.2.3
FH_VER_1 := $(shell echo $(FHVERSION) | cut -f1 -d.)
FH_VER_2 := $(shell echo $(FHVERSION) | cut -f2 -d.)
FH_VER_3 := $(shell echo $(FHVERSION) | cut -f3 -d.)

FHFLAGS=
ifeq (${FH_VER_2},11)
  FH2113 := $(shell [ $(FH_VER_3) -le 2 ] && echo true)
  ifeq (${FH2113},true)
    FH_VER_2 = 10
  endif
endif
ifeq (${FH_VER_2},9)
  ifeq (${FH_VER_3},5)
    FHFLAGS = -DFEYNHIGGS295
  endif
endif
ifeq (${FH_VER_2},10)
 FHFLAGS = -DFEYNHIGGS295 -DFEYNHIGGS2100
endif
ifeq (${FH_VER_2},11)
 FHFLAGS = -DFEYNHIGGS295 -DFEYNHIGGS2100 -DFEYNHIGGS2113
endif
ifeq (${FH_VER_2},12)
 FHFLAGS = -DFEYNHIGGS295 -DFEYNHIGGS2100 -DFEYNHIGGS2113 -DFEYNHIGGS2120
endif

ifeq ($(shell [ -e ${FHPATH}/build/version.h ] && echo true),true)
  $(shell cp ${FHPATH}/build/version.h src/sushi/tmp.h)
else
  $(shell echo "#define FHVERSION \"${FHVERSION}\"" > src/sushi/tmp.h)
endif

FLAGS=
ADDRUNCPPSRC=
ifeq (${predef},FH)
  PREDEFLIBS := ${PREDEFLIBS} -L${FHPATH}/build -lFH
  FLAGS=-DFEYNHIGGS ${FHFLAGS}
  ADDRUNCPPSRC:=src/sushi/runfeynhiggs.f
else
  ifeq (${predef},2HDMC)
  PREDEFLIBS := ${PREDEFLIBS} ${2HDMCPATH}/lib/lib2HDMC.a
  FLAGS := ${FLAGS} -DTHDMC
  LDFLAGS := -lgsl -lgslcblas -lm $(STDLIB) -lstdc++
  ADDRUNCPPSRC:=src/sushi/runthdmc.f
  else
    ifeq (${predef},HB)
      PREDEFLIBS := ${PREDEFLIBS} -L${FHPATH}/build -L${HBPATH} -lFH -lHB
      FLAGS=-DFEYNHIGGS -DHIGGSBOUNDS ${FHFLAGS}
      ADDRUNCPPSRC:=src/sushi/runfeynhiggs.f
    else
      ifeq (${predef},HS)
        PREDEFLIBS := ${PREDEFLIBS} -L${FHPATH}/build -L${HBPATH} -L${HSPATH} -lFH -lHS -lHB
        FLAGS=-DFEYNHIGGS -DHIGGSBOUNDS -DHIGGSSIGNALS ${FHFLAGS}
        ADDRUNCPPSRC:=src/sushi/runfeynhiggs.f
      else
        predef=PLAIN
      endif
    endif
  endif
endif

ifeq (${LHAPATH},)
 LHAPDF_CONFIG=lhapdf-config
 LHAPATH = $(shell $(LHAPDF_CONFIG) --libdir) -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)
else
 LHAPATH := $(LHAPATH) -Wl,-rpath,$(LHAPATH)
endif

ifeq (${F77},ifort)
	F77 = ifort
	OPT = -O3
else
	F77 = gfortran
	OPT = -O3
endif

fort=${F77} ${F77FLAGS} ${OPT}
RM=/bin/rm -f
CP=/bin/cp -f
LN = /bin/ln -sf

GGHNNLOSRC=controlggh.f	gghcall.f gghnnlo.f gghn3lo.f nnloexpand.f \
	n3loexpand.f ggh.f gghmatch.f nnlomtexpand.f
BBHNNLOSRC=bbh.f bbhcall.f bbhnnloexact.f bbhweights.f controlbbh.f
BBHDIFFSRC=bbHNLO.f functions_for_bbHNLO.f
SHARESRC=controlshare.F lumis.F elw.F ranlux.F vegas-ranlux.F
SUSYNLOSRC=ggamix_bare.f gghmix_CT3.f gghmix_CT6.f gghmix_bare3.f \
         gghmix_bare6.f gghmix_CT1.f gghmix_CT4.f gghmix_bare1.f \
         gghmix_bare4.f gghmix_CT2.f gghmix_CT5.f gghmix_bare2.f \
         gghmix_bare5.f
SUSHISRC=libsushi.f slhablocks.f cite.f inputoutput.f readslha.f real.f \
         evalcsusy.f integration.f renscheme.f sigma.f virtual.f runalpha.f \
         elw_lf.f functions.f tostring.f rensca_as.f
RUNCPPSRC=${ADDRUNCPPSRC} src/sushi/integrals.f src/sushi/init.f \
          src/sushi/virtsqqgluino.f


GGHNNLOOBJ=${GGHNNLOSRC:%.f=src/gghnnlo/%.o}
BBHNNLOOBJ=${BBHNNLOSRC:%.f=src/bbhnnlo/%.o}
BBHDIFFOBJ=${BBHDIFFSRC:%.f=src/bbhdiff/%.o}
SHAREOBJ=${SHARESRC:%.F=src/share/%.o}
SUSYNLOOBJ=${SUSYNLOSRC:%.f=src/susynlo/%.o}
TMP=${SUSHISRC:%.f=src/sushi/%.o}
SUSHIOBJ=${TMP:%.F=src/sushi/%.o}
RUNCPPOBJ=${RUNCPPSRC:%.f=%.o}

all:
	${RM} src/pdfspec/lhapdf.o
	${RM} src/sushi/init.o
	${RM} src/sushi/runfeynhiggs.o
	${RM} src/sushi/versions.h
	make bin/sushi.${predef}

bin/sushi.${predef}: src/sushi.F lib/libsushi${predef}.a lib/libshare.a
	${fort} -o bin/sushi.${predef} src/sushi.F \
	lib/libsushi${predef}.a lib/libshare.a ${PREDEFLIBS} \
        -lgfortran -L${LHAPATH} -lLHAPDF ${LDFLAGS} ${STDLIB}
	cd bin
	${LN} sushi.${predef} bin/sushi
	${LN} libsushi${predef}.a lib/libsushi.a

lib/libsushi${predef}.a: ${RUNCPPOBJ} ${SUSHIOBJ} ${SUSYNLOOBJ}
	ar -ruc lib/libsushi${predef}.a ${RUNCPPOBJ} ${SUSHIOBJ} ${SUSYNLOOBJ}

lib/libshare.a:  ${GGHNNLOOBJ} ${BBHNNLOOBJ} ${BBHDIFFOBJ} ${SHAREOBJ} \
           src/pdfspec/lhapdf.o
	ar -ruc lib/libshare.a ${GGHNNLOOBJ} ${BBHNNLOOBJ} \
             ${BBHDIFFOBJ} ${SHAREOBJ} src/pdfspec/lhapdf.o

gghnnlo: ${GGHNNLOOBJ}
bbhnnlo: ${BBHNNLOOBJ}
bbhdiff: ${BBHDIFFOBJ}
share: ${SHAREOBJ}
susynlo: ${SUSYNLOOBJ}
src: ${RUNCPPOBJ} ${SUSHIOBJ}	

src/sushi/virtsqqgluino.o: src/sushi/virtsqqgluino.f
	${fort} -c -cpp -o $@ $<
src/sushi/integrals.o: src/sushi/integrals.f
	${fort} -c -cpp -o $@ $<
src/sushi/init.o: src/sushi/init.f
	${fort} ${FLAGS} -c -cpp -o $@ $<
src/sushi/versions.h:
	mv src/sushi/tmp.h $@
	echo "#define HBVERSION \"${HBVERSION}\"" >> src/sushi/versions.h
	echo "#define HSVERSION \"${HSVERSION}\"" >> src/sushi/versions.h
	echo "#define THDMCVERSION \"${2HDMCVERSION}\"" >> src/sushi/versions.h
src/sushi/runfeynhiggs.o: src/sushi/runfeynhiggs.f src/sushi/HBandHSwithFH.h src/sushi/HBandHSwithFH.F src/sushi/versions.h
	${fort} ${FLAGS} -c -cpp -I${FHPATH}/src/include -I${HBPATH} -o $@ $<
src/sushi/runthdmc.o: src/sushi/runthdmc.f src/sushi/versions.h
	${fort} ${FLAGS} -c -cpp -o $@ $<
src/pdfspec/lhapdf.o: src/pdfspec/lhapdf.f
	${fort} ${FLAGS} -c -cpp -o $@ $<
src/%.o: src/%.f
	${fort} -c -o $@ $<
src/%.o: src/%.F
	${fort} -c -o $@ $<
src/gghnnlo/%.o: src/gghnnlo/%.f
	${fort} -c -o $@ $<
src/bbhnnlo/%.o: src/bbhnnlo/%.f
	${fort} -c -o $@ $<
src/bbhdiff/%.o: src/bbhdiff/%.f
	${fort} -c -o $@ $<
src/share/%.o: src/share/%.F
	${fort} -c -o $@ $<
src/susynlo/%.o: src/susynlo/%.f
	${fort} -c -o $@ $<
src/pdfspec/%.o: src/pdfspec/%.f
	${fort} -c -o $@ $<

someclean:
	${RM} lib/libsushi*.a
	${RM} lib/libshare.a
	${RM} src/sushi.o
	${RM} src/sushi/*.o
	${RM} src/share/*.o
	${RM} bin/sushi*
	${RM} src/pdfspec/lhapdf.o
	${RM} ${GGHNNLOOBJ}	
	${RM} ${BBHNNLOOBJ}	
	${RM} ${BBHDIFFOBJ}

clean: someclean
	${RM} ${SUSYNLOOBJ}
