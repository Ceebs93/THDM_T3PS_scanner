# Setup

## Instructions



1. If you have any local setup directives you might want to store them in `setup_local.sh` which is
   automatically sourced by `setup.sh` if it exists.

    E.g. my `setup_local.sh` contains:

    ~~~~
    source /home/de3u14/lib/build/miniconda3/bin/activate py27
    module load gsl
    ~~~~

    to setup `python 2.7` and load `gsl` libraries.


2. Source the `setup.sh` script standing in the root directory of this package
	to set up some environment variables.

    ~~~~
    source setup.sh
    ~~~~


3. Install the external packages

    - [`HiggsBounds`][HiggsBounds-url]
    - [`HiggsSignals`][HiggsSignals-url]
    - [`2HDMC`][2HDMC-url]
    - [`LHAPDF`][LHAPDF-url]

    Options:

    - a) There is 
        `./packages/install.sh`

        1. Check if the path variables are correctly set within the package `Makefile`s
			- `2HDMC`: Check if `HIGGSBOUNDS_PATH` and `HIGGSSIGNALS_PATH` are correctly specified in the `Makefile`.
			- `SusHi`: Check if `LHAPATH` and `2HDMCPATH` path variable are correctly specified in the `Makefile`
        2. Run `./packages/install.sh`

    - b) Manual install. Please refer to the packages manuals


------------------------------------------------------------

## `setup_local.sh``

~~~~
source <miniconda>/bin/activate py27
module load gsl
~~~~


## Paths


In `./MCMC/job_submission/config/mcmc_scan.conf`

~~~~
program = <path-of-THDM-T3PS-scanner>/2HDM-Processor/bin/ParameterScan_T3PS_with_HB_HS_FAST
~~~~


[2HDMC-url]: https://2hdmc.hepforge.org/
[HiggsBounds-url]: https://higgsbounds.hepforge.org/
[HiggsSignals-url]: https://higgsbounds.hepforge.org/
[SusHi-url]: https://sushi.hepforge.org/
[LHAPDF-url]: https://lhapdf.hepforge.org/
