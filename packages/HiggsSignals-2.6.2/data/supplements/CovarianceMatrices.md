# Theory covariance matrices for cross sections and BRs
In HiggsSignals there is one covariance matrix for the BSM model and one for the SM for each of the BRs (`data/BRcov.in` and `data/BRcovSM.in`), the 7+8TeV cross sections (`data/XScov.in` and `data/XScovSM.in`), and the 13TeV cross sections (`data/XScov_13TeV.in` and `data/XScovSM_13TeV.in`). By default the model specific matrices are always identical to the SM ones. If model specific theoretical correlations should be taken into account the corresponding files can be replaces by user generated ones.

## Covaiance matrices for the BRs
The covariance matrices can be evaluated by Toy MC, using the root scripts
`data/supplements/smearErrorsBR.cpp` and `data/supplements/smearErrorsBR_YR4.cpp`. Various parametric and
theoretical uncertainties on the partial widths and cross sections are smeared
according to the estimates from the YR.

### YR3 BR covariance matrix
The HiggsSignals-1 matrix can be obtained with the setting:

```c++
addPUandTHUlin = true;
deriveTHUsf = false; 
THUshape    = 1;
THUBoxErrorScaleFactor = 1.;
THU100percCorr = false;
THU100percCorrToPU = false;
```

This corresponds to

  BR    | relative error
--------|----------------
 gaga   | 0.0234686
 WW     | 0.0228754
 ZZ     | 0.0230075
 tautau | 0.0252033
 bb     | 0.0195539
 Zga    | 0.0368504
 cc     | 0.0880989
 mumu   | 0.0255934
 gg     | 0.0594597

The resulting covariance matrix on the BR's is contained in `BRcovSM.in.HS1`.

### LHCHXSWG covariance matrix

The BR uncertainty estimates of the LHC Higgs Cross Section WG are approximately reproduced with the following settings (in `smearErrors.cpp`):

```c++
addPUandTHUlin = true;
deriveTHUsf = false; 
THUshape    = 1;
THUBoxErrorScaleFactor = 5.;
THU100percCorr = false;
THU100percCorrToPU = false;
```

The resulting relative errors are:

  BR    | relative error
--------|----------------
 gaga   | 0.049463
 WW     | 0.0416186
 ZZ     | 0.0427421
 tautau | 0.06737
 bb     | 0.0329192
 Zga    | 0.150493
 cc     | 0.109137
 mumu   | 0.0704453
 gg     | 0.103211


The covariance matrix for the BR uncertainties is contained in `BRcovSM.in.LHCHXSWG`.
	
### YR4 covariance matrix

A covariance matrix that matches the YR4 estimates is the default in
HiggsSignals-2. It can be generated using
`data/supplements/smearErrorsBR_YR4.cpp` with the settings

```c++
bool addPUandTHUlin = true;
bool deriveTHUsf = false; 
int  THUshape    = 1; 
double THUBoxErrorScaleFactor = 1.;
bool THU100percCorr = false;
bool THU100percCorrToPU = false;
```

The resulting covariance matrix is in `BRcov.in.HS2`.

## Covariance matrix for the cross section uncertainties.

The toy MC evaluation for the cross sections uncertainties smears PDF+alpha_s, QCD scale and EW theoretical uncertainties. The scripts generate bot the correlation and the covariance matrices.

### 7+8 TEV covariance matrix

The 7+8TeV matrices can be obtained with the setting (in
`data/supplements/smearErrorsXS.cpp`):
```c++
shape            = 0; 
errorScaleFactor = 1.;
error100percCorr = false;
errorScaleZW100percCorr = true;
errorPdf100percCorr = false;
errorPdfOnlyProc100percCorr = true;
```

### 13 TeV covariance matrix

The 13TeV matrices use the same settings but in
`data/supplements/smearErrorsXS_LHC13.cpp`.
