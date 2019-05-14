# Usage instructions

The model scan conventionally consists of two separate parts, in order to solve the computational
bottleneck of calculating cross sections.
1. **MCMC scan.**
    1. Scan phase: Here only the Higgs signal strength and EWPOs are taken into account, and no check
    is made on theoretical constraints during the scan. 
    2. After all the chains have finished, they are trimmed (to discard points from the warm-up phase) and merged together.
    3. Checks on theoretical consistency are now enforced, creating a much smaller dataset in the end.
2. **Calculation of production cross sections (SusHi).**
    1. The calculation of the production cross section is only performed for points surviving the theoretical consistency checks.


## MCMC scan

Here the likelihood is composed of:
- Higgs signal strength measurements (HiggsSignals)
- EWPO (S,T)
Saved variables:
- Model input parameters
- Theoretical properties
    - potential shape parameters
    - coupling (modifiers)
    - validity (perturbativity, unitarity, stability)
- Experimental properties
    - Branching fractions

In the MCMC scan theoretical consistency constraints such as perturbativity, unitarity or stability
are not enforced


## Cross-section computation scan

Computation of cross sections with `SusHi`.
- gg -> A (NNLO QCD)
- gg -> H (NNLO QCD)


The computational speed, (processed points)/(second) is much slower in this step

