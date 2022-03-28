# Tasks for Magellan 

- [ ] Confirm upgrade of 2HDMC/HB/HS has been successful
- [ ] Create job running interface for local jobs
- [ ] Create job running interface for Torque jobs
- [ ] Install LHAPDF successfully
- [ ] Ensure LHAPDF is integrated into 2HDMC etc
- [ ] Install SuShi and test running of pipeline including SuSHi
- [ ] Update auto-install file to reflect the package upgrades
- [ ] Create a "initial-start-up.sh" to identify the top directory of installation and set this equal to $THDM_T3PS_scanner and write this directory to normal-start-up so that we can source it from any dir and still have the correct path, to identify current environment variables such as compilers, should set variables equal to these so that in auto_install.sh we can ask the user if they are happy to use these
- [ ] Edit auto-install file to run intial-start-up.sh and request user confirmation/input for compilers and environments
- [ ] Create a "normal-start-up.sh" to be run at the begining of each use, use code from source env.sh so that it will set $THDM_T3PS_scanner and look for an env_local.sh to source
- [ ] LHAPDF needs to be edited so that the config file will find the correct parts of python in conda
- [ ] Create job running interfaces for the SuSHi pipeline
- [ ] Work on installation guide for main directory
- [ ] Create docker and install ubuntu 20, Make 3.82, GCC 11.1.0, GSL 2.6, python2.7
- [ ] Install necessecary python modules to docker, pandas, numpy, scipy, hd5f etc
- [ ] Install Magellan to the docker
- [ ] Test Magellan in docker works
- [ ] Add MCMC theory to manual
- [ ] Add description of T3PS to manual
- [ ] Add 2HDM theory to manual
- [ ] Add usage of parameter-point generation pipeline, on cluster and locally to manual
- [ ] Add usage of full SuSHi pipline, on cluster and locally to manual
- [ ] Add description of parameterprocessor set up and compilation to manual
- [ ] Add description of how to change model type
- [ ] Add brief description of SuSHi usage
- [ ] Perform checks on points produced by final installation
- [ ] Add performance to manual
- [ ] Add explanation of config file for MCMC runs and how to change it to manual
- [ ] Add current python processing modules into main part of Magellan, get working in python2.7
- [ ] Add usage of processing modules to get csv or hd5f output.
- [ ] Add future plans for Magellan to manual
- [ ] Create some example config files

# Stretch Tasks for Magellan

- [ ] Look into uncoupling 2HDMC from parameter-processor in order to allow for people to use other models. This will require editing the config file for MCMC runs
- [ ] Look into integrating the plotting aspects of Magellan back into the toolbox as a whole
- [ ] Look into creating a processor to read SLHA files, and a corresponding config file for running this
