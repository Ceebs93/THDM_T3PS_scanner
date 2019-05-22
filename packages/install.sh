#!/usr/bin/env bash

## -- This script attempts to automatically install all the following packages:
#	 - HiggsBounds
#	 - HiggsSignals
#	 - 2HDMC
#	 - LHAPDF
#	 - SusHi

# - Exit if any error is found
set -e

if [ -z "${THDM_T3PS_SCANNER_DIR}" ]; then
    echo "Variable THDM_T3PS_SCANNER_DIR is not defined."
    echo "Please source setup.sh first in the root directory of the package."
	 exit 1
fi

###########################
### --- HiggsBounds --- ###
###########################

HiggsBounds_folder_name=HiggsBounds-5.3.2beta
HiggsBounds_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${HiggsBounds_folder_name}

echo "############################################################"
echo "### --- Attempting to install ${HiggsBounds_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${HiggsBounds_pkg_path}"
echo "./configure"
echo "make"
echo ""

cd ${HiggsBounds_pkg_path}
make clean
./configure
make


############################
### --- HiggsSignals --- ###
############################

HiggsSignals_folder_name=HiggsSignals-2.2.3beta
HiggsSignals_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${HiggsSignals_folder_name}

echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${HiggsSignals_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${HiggsSignals_pkg_path}"
echo "./configure"
echo "make"
echo ""

cd ${HiggsSignals_pkg_path}
make clean
./configure
make


#####################
### --- 2HDMC --- ###
#####################

THDMC_folder_name=2HDMC-1.7.0
THDMC_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${THDMC_folder_name}


echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${THDMC_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${THDMC_pkg_path}"
echo "make clean"
echo "make"
echo ""


cd ${THDMC_pkg_path}
make clean
make
make lib


######################
### --- LHAPDF --- ###
######################

LHAPDF_folder_name=LHAPDF-6.2.1
LHAPDF_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${LHAPDF_folder_name}
LHAPDF_pkg_build_path=${THDM_T3PS_SCANNER_DIR}/packages/${LHAPDF_folder_name}_build

echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${LHAPDF_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${LHAPDF_pkg_path}"
echo "make clean"
echo "./configure --prefix=${LHAPDF_pkg_path}"
echo "make"
echo "make install"

mkdir -p ${LHAPDF_pkg_build_path}
cd ${LHAPDF_pkg_path}
make clean
./configure --prefix=${LHAPDF_pkg_build_path}
make
make install


#####################
### --- SusHi --- ###
#####################

SusHi_folder_name=SusHi-1.6.1
SusHi_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${SusHi_folder_name}

echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${SusHi_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${SusHi_pkg_path}"
echo "./configure"
echo "make"
echo ""

cd ${SusHi_pkg_path}
make clean
./configure 
make predef=2HDMC


echo -e "install.sh script has finished."
echo -e "All packages should be installed."
