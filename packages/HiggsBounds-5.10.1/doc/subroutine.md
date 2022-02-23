# The HiggsBounds Subroutine Interface

## Using HiggsBounds through the Subroutines

The HiggsBounds subroutine interface is the recommended way of using HiggsBounds
in your code. All functions of this interface are located and documented in
`src/HiggsBounds_subroutines.F90`. The interface can be separated into four
parts.

### 1. Initialization

Before HiggsBounds can be used, the correct amount of storage needs to be
allocated and all required tabulated datafiles need to be read. This is done by
the subroutines initialize_higgsbounds() (or alternatively
initialize_higgsbounds_int()). These subroutines take the number of neutral and
charged Higgs bosons (needed to allocate the correct amount of storage) and a
flag specifying which analyses should be included (to read only the required
datafiles). It is recommended to use the full dataset `"LandH"`.

### 2. Input

In order to obtain the model predictions for all of the implemented analyses,
HiggsBounds requires a detailed input of the properties of the scalars in your
model. 

#### Neutral Input
For the neutral Higgs predictions you can either use the **effective coupling
input** or the **hadronic input** (the *partonic input* option from
HiggsBounds-4 is no longer available). It may also be useful to mix the two
approaches (see below).

1. For either input method use  higgsbounds_neutral_input_properties() to set
   the basic properties of the neutral Higgs bosons.

2. The remaining properties are set through the subroutines:
     - **effective coupling input**: higgsbounds_neutral_input_effc(),
       higgsbounds_neutral_input_nonsmbr()
     - **hadronic input**: higgsbounds_neutral_input_smbr(),
       higgsbounds_neutral_input_nonsmbr(), higgsbounds_neutral_input_lep(),
       higgsbounds_neutral_input_hadr()

It can be useful to combine the two methods by first using the effective
coupling input and then replacing some of the quantities by more precisely
calculated values. However, note that using any of the above input subroutines
will replace the values for *all* of their arguments. To set values for only a
specific quantity, use the higgsbounds_neutral_input_hadr_single(),
higgsbounds_neutral_input_hadr_double(), higgsbounds_neutral_input_lep_single()
and higgsbounds_neutral_input_lep_double() subroutines. 

#### Charged Input
There is no effective coupling input for the charged Higgs bosons, since there
is nothing in the SM to normalize their couplings to. The subroutines for the
charged Higgs input are higgsbounds_charged_input() and
higgsbounds_charged_input_hadr().


#### Uncertainties
In models where the Higgs masses are not input parameters you should also set
their theoretical uncertainties through higgsbounds_set_mass_uncertainties().


#### SLHA input
An alternative to all of the above is to use the HiggsBounds SLHA interface
through the higgsbounds_input_slha() subroutine. This interface is documented in
more detail in the [manual].


### 3. Running HiggsBounds

Once you have provided all the input you can actually run HiggsBounds though the
run_higgsbounds() subroutine and save the returned results. If you are
interested in the next-sensitive channels these can be obtained through
higgsbounds_get_most_sensitive_channels() and
higgsbounds_get_most_sensitive_channels_per_higgs().

If you used the SLHA input you may also want to use the SLHA output through
higgsbounds_slha_output().

### 4. Cleanup

When you are done with HiggsBounds you need to call finish_higgsbounds() to
clean up.

If you use HiggsBounds on many parameter points, make sure to only repeat steps
2 and 3 for each point.


## Phenomenological Quantities
HiggsBounds also provides access to many of the phenomenological quantities that
are used internally. This may be especially useful if you are using the
effective coupling input, in which case HiggsBounds can provide you with values
for many quantities derived from this input.

The higgsbounds_get_neutral_br() and higgsbounds_get_neutral_hadr_cs() provide
this information by returning the BRs into SM particles or the hadronic
production cross sections for the specified neutral Higgs boson. These can only
used after input has been provided, typically through
higgsbounds_neutral_input_effc().

The tabulated data used by HiggsBounds can also be accessed through a subroutine
interface found in `src/access_SM.f90` and `src/access_effC.f90`. Using these
subroutines only requires HiggsBounds to be initialized.


<!-- Links -->
[manual]: https://arxiv.org/????.?????
