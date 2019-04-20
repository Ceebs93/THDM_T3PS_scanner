#!/usr/bin/env bash

## -- This script attempts to automatically install all the following packages:
#	 - HiggsBounds
#	 - HiggsSignals
#	 - LHAPDF
#	 - SusHi
#	 - 2HDMC


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
./configure
make

#####################
### --- 2HDMC --- ###
#####################

THDMC_folder_name=2HDMC-1.7.0
THMDC_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${THDMC_folder_name}


echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${THDMC_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${THDMC_pkg_path}"
echo "make clear"
echo "make"
echo ""


cd ${THDMC_pkg_path}
make clear
make


######################
### --- LHAPDF --- ###
######################

LHAPDF_folder_name=LHAPDF-6.2.1
LHAPDF_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${LHAPDF_folder_name}

echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${LHAPDF_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd {LHAPDF_pkg_path}"
echo ""
echo ""
echo ""

cd ${LHAPDF_pkg_path}
./configure --prefix=${LHAPDF_pkg_path}
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
echo "cd {SusHi_pkg_path}"
echo ""
echo ""
echo ""

cd ${SusHi_pkg_path}
./configure 
make
