#!/usr/bin/env bash

echo -e "Sourcing THTDM_T3PS environment variables"
echo -e "source env.sh"
source env.sh

echo -e ""
echo -e "Installing external packages"
echo -e "./packages/install.sh"
./packages/install.sh


echo -e ""
echo -e "Setting up symbolic links"
echo -e "./setup_links.sh"
./setup_links.sh

echo -e ""
echo -e "Compiling ParameterPointProcessor binary"
echo -e "cd ParameterPointProcessor"
echo -e "make"
cd ParameterPointProcessor
make
