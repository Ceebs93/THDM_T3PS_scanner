# External packages

## `HiggsBounds`

**Website:** https://higgsbounds.hepforge.org/

> HiggsBounds takes a selection of Higgs sector predictions for any particular
> model as input and then uses the experimental topological cross section limits
> from Higgs searches at LEP, the Tevatron and the LHC to determine if this
> parameter point has been excluded at 95% C.L..


## `HiggsSignals`

**Website:** https://higgsbounds.hepforge.org/

> HiggsSignals performs a statistical test of the Higgs sector predictions of
> arbitrary models (using the HiggsBounds input routines) with the measurements
> of Higgs boson signal rates and masses from the Tevatron and the LHC.


## `2HDMC`

**Website:**  https://2hdmc.hepforge.org/

> 2HDMC is a general-purpose calculator for the two-Higgs doublet model. It allows parametrization of
> the Higgs potential in many different ways, convenient specification of generic Yukawa sectors, the
> evaluation of decay widths (including higher-order QCD corrections), theoretical constraints and
> much more.

`2HDMC` can be linked together with `HiggsBounds` and `HiggsSignals`, which
makes it possible to create a C++ binary that has access to all three packages.


## `LHAPDF`

**Website:** https://lhapdf.hepforge.org/

> LHAPDF is a general purpose C++ interpolator, used for evaluating PDFs from
> discretised data files.


## `SusHi`

**Website:** https://sushi.hepforge.org/

> SusHi (Supersymmetric Higgs) is a Fortran code, which calculates Higgs cross
> sections in gluon fusion and bottom-quark annihilation at hadron colliders in
> the SM, the 2HDM, the MSSM and the NMSSM.

**Dependencies:**
- 2HDMC
- LHAPDF


Inside the `Makefile`

~~~~
2HDMCPATH = ${THDM_T3PS_SCANNER_DIR}/packages/2HDMC-1.7.0
2HDMCVERSION = 1.7.0
# Specify the path to the compiled LHAPDF libraries (might be found automatically):
LHAPATH = ${THDM_T3PS_SCANNER_DIR}/packages/LHAPDF-6.2.1-install/lib
# Specify the path to standard C++ libraries (libstdc)
~~~~
