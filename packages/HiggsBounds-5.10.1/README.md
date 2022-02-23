# HiggsBounds-5

HiggsBounds takes a selection of Higgs sector predictions for any particular
model as input and then uses the experimental topological cross section limits
from Higgs searches at LEP, the Tevatron and the LHC to determine if this
parameter point has been excluded at 95% C.L..

**Download the newest [release] or take a look at the [documentation].**

**[The HiggsBounds-5 manual is available here][manual].**

### The HiggsBounds Collaboration
Current members of the HiggsBounds team are Philip Bechtle, Sven Heinemeyer,
Tobias Klingl, Tim Stefaniak, Georg Weiglein and Jonas Wittbrodt.

In case of questions regarding HiggsBounds please contact [Tim] and/or [Jonas].
If you you want to report a bug please open an [issue].

Former members are Oliver Brein, Daniel Dercks, Oscar Stål and Karina E. Williams.

### Journal References

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Georg Weiglein, 
    Karina E. Williams                                                                <br/>
    *HiggsBounds: Confronting Arbitrary Higgs Sectors with Exclusion Bounds from
    LEP and the Tevatron*                                                             <br/>
    Comput.Phys.Commun.181:138-167, 2010, e-Print: [arXiv:0811.4169] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Georg Weiglein, 
    Karina E. Williams                                                                <br/>
    *HiggsBounds 2.0.0: Confronting Neutral and Charged Higgs Sector Predictions 
    with Exclusion Bounds from LEP and the Tevatron*                                  <br/>
    Comput.Phys.Commun. 182:2605-2631, 2011, e-Print: [arXiv:1102.1898] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, 
    Georg Weiglein, Karina Williams                                                   <br/>
    *Recent Developments in HiggsBounds and a Preview of HiggsSignals*                <br/>
    PoS CHARGED2012 (2012) 024, e-Print: [arXiv:1301.2345] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Oscar Stål, Tim Stefaniak,
    Georg Weiglein, Karina Williams                                                   <br/>
    *HiggsBounds-4: Improved Tests of Extended Higgs Sectors against Exclusion
    Bounds from LEP, the Tevatron and the LHC*                                        <br/>
    Eur.Phys.J C74 (2014) 2693, e-Print: [arXiv:1311.0055] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein        <br/>
    *Applying Exclusion Likelihoods from LHC Searches to Extended Higgs Sectors*      <br/>
    Eur.Phys.J. C75 (2015) no.9, 421, e-Print: [arXiv:1507.06706] [hep-ph]

  - Philip Bechtle, Daniel Dercks, Sven Heinemeyer, Tobias Klingl,  Tim Stefaniak,
    Georg Weiglein, Jonas Wittbrodt                                                   <br/>
    *HiggsBounds-5: Testing Higgs Sectors in the LHC 13 TeV Era*                      <br/>
    Eur.Phys.J.C 80 (2020) 12, 1211, e-Print: [arxiv:2006.06007] [hep-ph]

  - Henning Bahl, Victor Martin Lozanoa, Tim Stefaniak, Jonas Wittbrodt              <br/>
    *Testing Exotic Scalars with HiggsBounds*                                         <br/>
    e-Print [arxiv:2109.?????] [hep-ph]

### Experimental Results

HiggsBounds incorporates a large number of experimental results. The list of
implemented analyses can be accessed locally through the `AllAnalyses`
executable (see below) including InspireHEP cite keys and arXiv or experimental
report numbers for easy citing. An up-to date version of that list is also
available [here][AllAnalyses] (as `HBAnalyses.txt`) together with a
`HBAnalyses.bib` file and an example tex file citing all HiggsBounds analyses. 

## Compilation

HiggsBounds requires a Fortran compiler supporting at least Fortran 95
(like `gfortran`) and `cmake`.

The code is compiled by

    mkdir build && cd build
    cmake ..
    make

If your distribution does not provide a recent version of CMake (e.g. Scientific
Linux or CentOS) simply grab the binary version from [here][cmake-download].

If you wish to use a different Fortran compiler set the `FC` environment
variable before running `cmake`.

### LEP chi-squared extension
In order to use the chi-squared values from LEP simply add

    -DLEP_CHISQ=ON

to the `cmake ..` line. This automatically downloads the required data files
(see e.g. `example_programs/HBwithLEPlikelihood.F90` for an example using this
extension).


### FeynHiggs examples
Some of the example_programs illustrate the use of HiggsBounds together with
[FeynHiggs]. To use these a compiled version of FeynHiggs is required and

    -DFeynHiggs_ROOT=/path/to/feynhiggs

has to be set in the to the `cmake ..` line. See e.g.
`example_programs/HBwithFH.F90` for an example.

### Tests
A test suite that checks the results of the example codes against reference
values is included and can be executed through `make check` or (for more
detailed results) `ctest --output-on-failure`. Running the tests requires
`python>3.5` with the `pytest`, `numpy` and `pandas` packages. The test suite
takes a few minutes to complete.

## Library
HiggsBounds can be used as a library of Fortran subroutines (see
`src/HiggsBounds_subroutines.F90` and the [subroutine documentation](doc/subroutine.md)) by linking the `build/lib/libHB.a` library.
If your project builds using CMake you can simply use
```cmake
find_package(HiggsBounds)
# ...
target_link_libraries(YourTarget PRIVATE HiggsBounds::HB)
```
and CMake will take care of everything else.

### C Interface
HiggsBounds provides a C interface in `include/HiggsBounds.h` that wraps the
Fortran subroutines and takes care of type conversions. When using the C
interface, make sure to link any libraries required by the Fortran C interface
of your compiler setup (e.g. `-lgfortran -lm -lquadmath -lm` for gcc-gfortran).
If you use CMake this will be handled automatically. 


## Executables
Compilation generates a binary which can be invoked from the command line as

	./HiggsBounds <options>

In addition an `AllAnalysis` executable which prints all implemented
experimental results, and all `example_programs` whose dependencies are met are
compiled.

## License

[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

<!-- links -->
[Tim]: mailto:tim.stefaniak@desy.de
[Jonas]: mailto:jonas.wittbrodt@desy.de
[issue]: https://gitlab.com/higgsbounds/higgsbounds/issues
[arXiv:0811.4169]: https://inspirehep.net/record/803530
[arXiv:1102.1898]: https://inspirehep.net/record/889030
[arXiv:1301.2345]: https://inspirehep.net/record/1210431
[arXiv:1311.0055]: https://inspirehep.net/record/1263076
[arXiv:1507.06706]: https://inspirehep.net/record/1384775
[arxiv:2006.06007]: https://inspirehep.net/literature/1800733
[AllAnalyses]: https://gitlab.com/higgsbounds/higgsbounds/-/jobs/artifacts/master/browse/HBAnalyses?job=HBAnalyses
[FeynHiggs]: http://www.feynhiggs.de/
[cmake-download]: https://cmake.org/download/
[release]: https://gitlab.com/higgsbounds/higgsbounds/-/releases
[documentation]: https://higgsbounds.gitlab.io/higgsbounds/
[Doxygen]: http://www.doxygen.nl/
[manual]: https://arxiv.org/abs/2006.06007
