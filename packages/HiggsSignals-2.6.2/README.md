# HiggsSignals-2

HiggsSignals performs a statistical test of the Higgs sector predictions of
arbitrary models with the measurements of Higgs boson signal rates and masses
from the Tevatron and the LHC.

**Looking for the newest [release] or the [documentation]?**

The input routines are documented in the [HiggsBounds] [documentation][hbdoc] and [manual][hbmanual].

### The HiggsBounds Collaboration
HiggsSignals is developed by the HiggsBounds Collaboration. Current members are
Philip Bechtle, Sven Heinemeyer, Tobias Klingl, Tim Stefaniak, Georg Weiglein,
and Jonas Wittbrodt. 

In case of questions regarding HiggsSignals please contact [Jonas].
If you you want to report a bug please open an [issue].

Former members are Oliver Brein, Daniel Dercks, Oscar Stål and Karina E.
Williams.

### Journal References

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein    <br/>
    "HiggsSignals: Confronting arbitrary Higgs sectors with measurements at the
    Tevatron and the LHC"                                                         <br/>
    Eur.Phys.J. C74 no.2, 2711, 2014, e-Print: [arxiv:1305.1933] [hep-ph]

  - Oscar Stål, Tim Stefaniak  <br/>
    "Constraining extended Higgs sectors with HiggsSignals"                       <br/>
    PoS EPS-HEP2013 314, 2013, e-Print: [arxiv:1310.4039] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein    <br/>
    "Probing the Standard Model with Higgs signal rates from the Tevatron, 
    the LHC and a future ILC"                                                     <br/>
    JHEP 1411 039, 2014, e-Print: [arxiv:1403.1582] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Tobias Klingl, Tim Stefaniak, 
    Georg Weiglein, Jonas Wittbrodt                                                <br/>
    "HiggsSignals-2: Probing new physics with precision Higgs measurements 
    in the LHC 13 TeV era"                                                         <br/>
    e-Print: [arxiv:2012.09197] [hep-ph]


HiggsSignals relies on [HiggsBounds] for input, so please also cite the
references listed there.

## Installation

HiggsBounds requires a Fortran compiler supporting at least Fortran 95 (like
`gfortran`), `cmake` as well as [HiggsBounds]. HiggsBounds has to be available
on the system in compiled form. Its location will be automatically found by
CMake. If you have several versions of HiggsBounds installed, set the
`HiggsBounds_DIR` environment variable to the build folder of the version you
want to use before calling `cmake`.

The code is compiled by 
  
    mkdir build && cd build
    cmake ..
    make

If your distribution does not provide a recent version of CMake (e.g. Scientific
Linux or CentOS) simply grab the binary version from [here][cmake-download].

If you wish to use a different Fortran compiler set the `FC` environment
variable before running `cmake`.

### Tests
A test suite that checks the results of the example codes against reference
values is included and can be executed through `make check` or (for more
detailed results) `ctest --output-on-failure`. Running the tests requires
`python>3.5` with the `pytest`, `numpy` and `pandas` packages. The test suite
takes a few minutes to complete.


### FeynHiggs examples
Some of the `example_programs` illustrate the use of HiggsSignals together with
[FeynHiggs]. To use these a compiled version of FeynHiggs is required and the
path to the FeynHiggs folder has to specified by adding

    -DFeynHiggs_ROOT=/path/to/feynhiggs

to the `cmake ..` line. See e.g. `example_programs/HBandHSwithFH.F90` for an
example.

## Usage

Compilation generates a binary which can be invoked from the command line as 
	
	./HiggsSignals <options>

and the library

	libHS.a

which can be linked to other codes to use the subroutine interface (see
`src/HiggsSignals_subroutines.F90`). If your project builds using CMake you can simply use
```cmake
find_package(HiggsBounds)
find_package(HiggsSignals)
# ...
target_link_libraries(YourTarget PRIVATE HiggsBounds::HB HiggsSignals::HS)
```
and CMake will take care of everything else.

### C Interface
HiggsSignals provides a C interface in `include/HiggsSignals.h` that wraps the
most important Fortran subroutines and takes care of type conversions. When
using the interface, make sure to also link any libraries required by the
Fortran C interface of your compiler setup (e.g. `-lgfortran -lm -lquadmath -lm`
for gcc-gfortran). If you use CMake this will be done automatically.

## License

[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

<!-- links -->
[Jonas]: mailto:jonas.wittbrodt@desy.de
[issue]: https://gitlab.com/higgsbounds/higgssignals/issues
[arxiv:1305.1933]: https://inspirehep.net/record/1232501
[arxiv:1310.4039]: https://inspirehep.net/record/1258615
[arxiv:1403.1582]: https://inspirehep.net/record/1284976
[arxiv:2012.09197]: https://inspirehep.net/literature/1837082
[HiggsBounds]: https://gitlab.com/higgsbounds/higgsbounds
[FeynHiggs]: http://www.feynhiggs.de/
[cmake-download]: https://cmake.org/download/
[release]: https://gitlab.com/higgsbounds/higgssignals/-/releases
[hbdoc]: https://higgsbounds.gitlab.io/higgsbounds/
[hbmanual]: https://arxiv.org/abs/2006.06007
[documentation]: https://higgsbounds.gitlab.io/higgssignals/
