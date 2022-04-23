# Setup

----------Current Instructions-----------------

env.sh sets the path for ${THDM_T3PS_scanner} that gets used throughout everything. Unfortunately it currently sets it by using pwd so env.sh needs to always be sourced from that top directory or everything will break. When doing the install comment out the if statement, some programs have tests that can be run only using python3 so its just easier to ignore env_local.sh which loads a python2 env.
load gcc/11.1.0, gsl/2.6 and cmake/3.22.0

cd into packages/HiggsBounds-5.10.1 and do mkdir build && cd build , cmake .. , make. There should be a tests dir now, and with python>3.5 pytest, numpy and pandas it should be possible to run, ctest --output-on-failure to check the install. If in doubt, there is a github for this https://gitlab.com/higgsbounds/higgsbounds which should cover anything I`ve forgotten.

Set the HiggsBound_DIR environment variable to the build folder of HiggsBounds to make sure the next program can find it.
cd into packages/HiggsSignals-2.6.2 and do mkdir build && cd build , cmake .. , make . Again there should be a tests dir now and ctest --output-on-failure can be run with the same python + modules as before. This shares a github with HiggsBounds https://gitlab.com/higgsbounds/higgssignals

Next is LHAPDF, there`s currently two versions just due to issues installing the later one, so go with the older. This isn`t hugely important though and if theres problems it can be skipped. Do ./configure --prefix/path/for/installation , make , make install - in case I`ve forgotten something https://lhapdf.hepforge.org/install.html#install_quickstart .

Then 2HDMC, a library from HiggsBounds and HiggsSignals needs to be copied into the lib directory of 2HDMC, then do make . After that if ./Demo runs successfully it should be okay.

SusHi isn`t currently in use so can be ignored.

There is a THDM-env.txt which contains the information on the conda environment used and conda env create --file THDM-env.txt can be used to recreate it.
Then env_local.sh can be sourced
run setup_links.sh
cd into ParameterPointProcessor and do make
That should be the installation done.


-----------Old Instructions---------------------
## Instructions

There is a main installation script `setup_auto.sh`, which aims to setup everything necessary for
the parameters scanner.
It is unlikely that it will work out-of-the-box due to mismatching package versions, missing
environment variables etc. When an error message is prompted, the user expected to debug.
Below are the steps that `setup_auto.sh` aims to do.

In case you run into errors, see the [setup troubleshooting docs][setup_troubleshoot] which might be
able to help with your issue.

1. **Setting up local directives, environment variables (`env_local.sh`).**

	If you have any local setup directives you may want to store them in `env_local.sh` which is
   automatically sourced by `env.sh` if it exists.

    E.g. my `env_local.sh` contains:

    ~~~~
    source /home/de3u14/lib/build/miniconda3/bin/activate py27
    module load gsl
    module load gcc/6.1.0
    ~~~~

    Which loads:
    - `python 2.7` environment
    - the `gsl` libraries (needed by 2HDMC)
    - `gcc` 6.1.0 version

2. **Source `env.sh`**

	Source the `env.sh` script standing in the root directory of this package
	to set up some environment variables.
	
	~~~~
	source env.sh
	~~~~
	
	The environment variables which are set up are:
	- `THDM_T3PS_SCANNER_DIR`: The path to the root of the scanner package
	- `PATH` environment variables is prepended by `packages/T3PS` in order the give access to the `T3PS` executable.


3. **Installing the external packages.**

    - [`HiggsBounds`][HiggsBounds-url]
    - [`HiggsSignals`][HiggsSignals-url]
    - [`2HDMC`][2HDMC-url]
    - [`LHAPDF`][LHAPDF-url]
	- [`SusHi`][SusHi-url]

    Your options:

    - a) This repo already ships with the source code of the above packages, locate in the [`./packages`](../packages/) directory.
		There is a script provided which attempts to automatically install these packages, located at:
		
		~~~~
        ./packages/install.sh
		~~~~
	
		Before executing this script.

        1. Check if the path variables are correctly set within the package `Makefile`
			- `2HDMC`: Check if `HIGGSBOUNDS_PATH` and `HIGGSSIGNALS_PATH` are correctly specified in the `Makefile`.
			- `SusHi`: Check if `LHAPATH` and `2HDMCPATH` path variable are correctly specified in the `Makefile`

        2. Run `./packages/install.sh`

	- b) You can also opt to manually install these packages one-by-one. In
	 	this case please refer to the documentation of the packages.


4. **Setup symbolis links.**

    This is done by the `setup_links.sh` script.
    It creates symbolic links of the headers and libraries of the various packages in
    the `./links` directory.

5. **Compile `ParameterPointProcessor`.**

    ~~~~
    cd ParameterPointProcessor
    make
    ~~~~

6. **Installation complete. Start using the scanner.**
	
	At this point hopefully you have a working setup. For instructions on how
	to submit jobs please refer to [usage.md](./usage.md).


------------------------------------------------------------

[setup_troubleshoot]: ./setup_troubleshoot.md
[2HDMC-url]: https://2hdmc.hepforge.org/
[HiggsBounds-url]: https://higgsbounds.hepforge.org/
[HiggsSignals-url]: https://higgsbounds.hepforge.org/
[SusHi-url]: https://sushi.hepforge.org/
[LHAPDF-url]: https://lhapdf.hepforge.org/
