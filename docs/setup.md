# Setup

## Instructions

Automatic installation `setup_auto.sh`.



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



4. **Installation complete. Start using the scanner.**
	
	At this point hopefully you have a working setup. For instructions on how
	to submit jobs please refer to [usage.md](./usage.md).


------------------------------------------------------------


[2HDMC-url]: https://2hdmc.hepforge.org/
[HiggsBounds-url]: https://higgsbounds.hepforge.org/
[HiggsSignals-url]: https://higgsbounds.hepforge.org/
[SusHi-url]: https://sushi.hepforge.org/
[LHAPDF-url]: https://lhapdf.hepforge.org/
