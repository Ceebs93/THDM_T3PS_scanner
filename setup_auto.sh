#!/usr/bin/env/bash

echo -e "Defining THDM_T3PS_SCANNER_DIR variable and exporting through scanner"
export THDM_T3PS_SCANNER_DIR=${PWD}

echo ${THDM_T3PS_SCANNER_DIR}
#export ${THDM_T3PS_SCANNER_DIR}

bash utils/set_home_dir.sh

#echo -e "Loading modules needed for compilation"
#echo -e "module load gcc/11.1.0"
#echo -e "module load gsl/2.6"
#echo -e "module load cmake/3.22.0"
#echo -e "module load conda"
#module load gcc/11.1.0
#module load gsl/2.6
#module load cmake/3.22.0
#module load conda

echo -e "Sourcing THDM_T3PS environment variables"
echo -e "source env.sh"
source env.sh

echo -e ""
echo -e "Installing external packages"
echo -e "./packages/install.sh"
./packages/install.sh


echo -e ""
echo -e "Setting up symbolic links"
echo -e "./utils/setup_links.sh"
./utils/setup_links.sh

#echo -e ""
#echo -e "Compiling ParameterPointProcessor binary"
#echo -e "cd ParameterPointProcessor"
#echo -e "make"
#cd ParameterPointProcessor
#make
