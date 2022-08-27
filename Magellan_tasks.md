# Tasks for Magellan 

- [X] Confirm upgrade of 2HDMC/HB/HS has been successful
- [X] Create job running interface for local MCMC jobs
- [ ] Create job running interface for local madgraph jobs
- [ ] Create job running interface for Torque jobs
- [X] Install LHAPDF successfully
- [ ] ~~Ensure LHAPDF is integrated into 2HDMC etc~~
- [ ] Upload conversion tools for dat to csv files from MCMC output and edit appropriately
- [ ] Edit merge-jobs for MadGraph to ensure it is functional for csves
- [X] Rename SusHi directory to MadGraph
- [X] Install MadGraph and test running of pipeline including MadGraph
- [ ] Update auto-install file to reflect the package upgrades
- [ ] Create a "initial-start-up.sh" to identify the top directory of installation and set this equal to $THDM_T3PS_scanner and write this directory to normal-start-up so that we can source it from any dir and still have the correct path, to identify current environment variables such as compilers, should set variables equal to these so that in auto_install.sh we can ask the user if they are happy to use these
- [ ] Edit auto-install file to run intial-start-up.sh and request user confirmation/input for compilers and environments
- [ ] Create a "normal-start-up.sh" to be run at the begining of each use, use code from source env.sh so that it will set $THDM_T3PS_scanner and look for an env_local.sh to source
- [X] LHAPDF needs to be edited so that the config file will find the correct parts of python in conda - Ended up using the LHAPDF installation in MadGraph instead. Should still be reachable by 2HDMC
- [X] Create job running interfaces for the ~~SusHi~~ Madgraph pipeline
- [ ] Work on installation guide for main directory
- [ ] Create docker and install ubuntu 20, Make 3.22, GCC 11.1.0, GSL 2.6, python2.7
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
1
# Tasks for Magellan 
2
​
3
- [X] Confirm upgrade of 2HDMC/HB/HS has been successful
4
- [X] Create job running interface for local jobs
5
- [ ] Create job running interface for Torque jobs
6
- [X] Install LHAPDF successfully
7
- [ ] ~~Ensure LHAPDF is integrated into 2HDMC etc~~
8
- [ ] Upload conversion tools for dat to csv files from MCMC output and edit appropriately
9
- [ ] Edit merge-jobs for MadGraph to ensure it is functional for csves
10
- [X] Rename SusHi directory to MadGraph
11
- [X] Install MadGraph and test running of pipeline including MadGraph
12
- [ ] Update auto-install file to reflect the package upgrades
13
- [ ] Create a "initial-start-up.sh" to identify the top directory of installation and set this equal to $THDM_T3PS_scanner and write this directory to normal-start-up so that we can source it from any dir and still have the correct path, to identify current environment variables such as compilers, should set variables equal to these so that in auto_install.sh we can ask the user if they are happy to use these
14
- [ ] Edit auto-install file to run intial-start-up.sh and request user confirmation/input for compilers and environments
15
- [ ] Create a "normal-start-up.sh" to be run at the begining of each use, use code from source env.sh so that it will set $THDM_T3PS_scanner and look for an env_local.sh to source
16
- [X] LHAPDF needs to be edited so that the config file will find the correct parts of python in conda - Ended up using the LHAPDF installation in MadGraph instead. Should still be reachable by 2HDMC
17
- [X] Create job running interfaces for the ~~SusHi~~ Madgraph pipeline
18
- [ ] Work on installation guide for main directory
19
- [ ] Create docker and install ubuntu 20, Make 3.82, GCC 11.1.0, GSL 2.6, python2.7
20
- [ ] Install necessecary python modules to docker, pandas, numpy, scipy, hd5f etc
21
- [ ] Install Magellan to the docker
22
- [ ] Test Magellan in docker works
23
- [ ] Add MCMC theory to manual
24
- [ ] Add description of T3PS to manual
25
- [ ] Add 2HDM theory to manual
26
- [ ] Add usage of parameter-point generation pipeline, on cluster and locally to manual
27
- [ ] Add usage of full SuSHi pipline, on cluster and locally to manual
28
- [ ] Add description of parameterprocessor set up and compilation to manual
29
- [ ] Add description of how to change model type
30
- [ ] Add brief description of SuSHi usage
31
- [ ] Perform checks on points produced by final installation
32
- [ ] Add performance to manual
33
- [ ] Add explanation of config file for MCMC runs and how to change it to manual
34
- [ ] Add current python processing modules into main part of Magellan, get working in python2.7
35
- [ ] Add usage of processing modules to get csv or hd5f output.
36
- [ ] Add future plans for Magellan to manual
37
- [ ] Create some example config files
38
​
39
# Stretch Tasks for Magellan
40
​
41
- [ ] Look into uncoupling 2HDMC from parameter-processor in order to allow for people to use other models. This will require editing the config file for MCMC runs

- [ ] Add future plans for Magellan to manual
- [ ] Create some example config files

# Stretch Tasks for Magellan

- [ ] Look into uncoupling 2HDMC from parameter-processor in order to allow for people to use other models. This will require editing the config file for MCMC runs
- [ ] Look into integrating the plotting aspects of Magellan back into the toolbox as a whole
- [ ] Look into creating a processor to read SLHA files, and a corresponding config file for running this
