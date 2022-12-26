#!/usr/bin/env bash

## -- This script attempts to automatically install all the following packages:
#	 - HiggsBounds
#	 - HiggsSignals
#	 - 2HDMC
#	 - LHAPDF
#	 - MG5

# - Exit if any error is found
set -e

if [ -z "${THDM_T3PS_SCANNER_DIR}" ]; then
    echo "Variable THDM_T3PS_SCANNER_DIR is not defined."
    echo "Please source setup.sh first in the root directory of the package."
	 exit 1
fi

# Want to set up compilers explicitly and export them to ensure that compliations are using the correct ones



###########################
### --- HiggsBounds --- ###
###########################

HiggsBounds_folder_name=HiggsBounds-5.10.1
HiggsBounds_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${HiggsBounds_folder_name}

echo "############################################################"
echo "### --- Attempting to install ${HiggsBounds_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${HiggsBounds_pkg_path}"
echo "mkdir build && cd build"
echo "cmake .."
echo "make"
echo ""

cd ${HiggsBounds_pkg_path}
mkdir build && cd build
cmake ..
make


############################
### --- HiggsSignals --- ###
############################

HiggsSignals_folder_name=HiggsSignals-2.6.2
HiggsSignals_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${HiggsSignals_folder_name}

echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${HiggsSignals_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${HiggsSignals_pkg_path}"
echo "mkdir build && cd build"
echo "cmake .."
echo "make"
echo ""

cd ${HiggsSignals_pkg_path}
mkdir build && cd  build
cmake ..
make


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
echo "./configure --prefix=${LHAPDF_pkg_path}"
echo "make"
echo "make install"

mkdir -p ${LHAPDF_pkg_build_path}
cd ${LHAPDF_pkg_path}
./configure --prefix=${LHAPDF_pkg_build_path}
make
make install


#####################
### --- 2HDMC --- ###
#####################

THDMC_folder_name=2HDMC-1.8.0
THDMC_pkg_path=${THDM_T3PS_SCANNER_DIR}/packages/${THDMC_folder_name}


echo -e "\n\n\n"
echo "############################################################"
echo "### --- Attempting to install ${THDMC_folder_name} --- ####"
echo "############################################################"
echo ""
echo "cd ${THDMC_pkg_path}"
echo "cp ../HiggsBounds-5.10.1/build/lib/libHB.a lib/libHB.a"
echo "cp ../HiggsSignals-2.6.2/build/lib/libHS.a lib/libHS.a"
echo "make"
echo ""


cd ${THDMC_pkg_path}
# copying HiggsBounds and HiggsSignals libraries in order to run 2HDMC with them
cp ../HiggsBounds-5.10.1/build/lib/libHB.a lib/libHB.a
cp ../HiggsSignals-2.6.2/build/lib/libHS.a lib/libHS.a
make
make lib


########################
### --- MadGraph --- ###
########################

MG5_dir=MG5_aMC_3_1_0

echo -e "\n\n\n"
echo "######################################################"
echo "### --- Downloading and setting up ${MG5_dir} --- ####"
echo "######################################################"

echo "cd ${THDM_T3PS_SCANNER_DIR}/packages"
echo "wget https://launchpad.net/mg5amcnlo/3.0/3.1.x/+download/MG5_aMC_v3.1.0.tar.gz"
echo "tar -zxf MG5_v3.1.0.tar.gz"
echo "cd ${MG5_dir}"
echo "mv ../THDM_type1_UFO ./models"
echo "mv ../2HDMtII_NLO ./models"

cd ${THDM_T3PS_SCANNER_DIR}/packages
wget https://launchpad.net/mg5amcnlo/3.0/3.1.x/+download/MG5_aMC_v3.1.0.tar.gz
tar -zxf MG5_v3.1.0.tar.gz
cd ${MG5_dir}
mv ../THDM_type1_UFO ./models
mv ../2HDMtII_NLO ./models

echo "You should now run MadGraph in order to install LHAPDF and other software you may wish to use within it."

echo -e "install.sh script has finished."
echo -e "All packages should be installed."
