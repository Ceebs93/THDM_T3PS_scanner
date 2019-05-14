# Usage instructions


**Overview:**

The model scan conventionally consists of two separate parts, in order to solve the computational
bottleneck of calculating cross sections.

1. **MCMC scan.**
    1. Scan phase: Here only the Higgs signal strength and EWPOs are taken into account, and no check
    is made on theoretical constraints during the scan. 
    2. After all the chains have finished, they are trimmed (to discard points from the warm-up phase) and merged together.
    3. Checks on theoretical consistency are now enforced, creating a significantly smaller dataset in the end.
2. **Calculation of production cross sections (SusHi).**
    1. The calculation of the production cross section is only performed for points surviving the theoretical consistency checks.

Note that due to the need for adequeate number of points, both steps are performed on a cluster.
The framework in its current form is compatible only with PBS cluster.

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------


## MCMC scan

### Description 

Here the likelihood is composed of:
- Higgs signal strength measurements (HiggsSignals)
- EWPO (S,T)

The saved variables are:
- Model input parameters
- Theoretical properties
    - potential shape parameters
    - coupling (modifiers)
    - validity (perturbativity, unitarity, stability)
- Experimental properties
    - Branching fractions

The scan engine is `T3PS`.

During the course of MCMC scan theoretical consistency constraints such as perturbativity, unitarity
or stability are not enforced, these only happen at the end. 


### Parameter point processor

The parameter point processor is a `C++` binary, which is called by `T3PS`.
It is found under the `ParameterPointProcessor` directory.

**Input:**
- Model parameters. 

**Output:**
- Model parameters. 
- Theoretical properties
- Branching fractions
- &chi;<sup>2</sup> of HiggsSignals, EWPO
- ...

The &chi;<sup>2</sup> values are read off by T3PS and used to construct the likelihood.

Default parameter point processor:
`ParameterScan_T3PS_with_HB_HS.cpp`

### Job submission

Associated work area is found in [./job_submission/MCMC/](../job_submission/MCMC/).

`Makefile`

**Directories:**
- `utils`: contains utility `shell` scripts which received input from the `Makefile` and help with
    tasks, such as job submission preperations, job submission, and merging of the jobs.
    Contains:
    - `create-jobs.sh`: Prepares the jobs passed on to it by the `Makefile`.
    - `submit-jobs.sh`: Submits jobs to PBS cluster.
    - `merge-jobs.sh`: Once the jobs are finished helps to merge the output files.
    - `convert_to_hdf.py`: Converts ASCII to HDF5.
- `template`: contains `T3PS` template file.
- `headers`: contains the header for the merged datafile, specifying each columns in the final
    dataset. Needs to be consistent with Parameter Point Processor.
- `config`: contains the `T3PS` config templates. The token placeholders `program_`, `nCores_`,
    `chain_length_` are replaced by 
- `jobs`: folder contains the invididual jobs in format of `JOB_TAG/job_XXX`. This directory contains:
    - `t3ps.conf`: the `T3PS` configuration of the submitted jobs
    - `job.template`: the job template needed by `T3PS`.


**Instructions:**
1. Create cluster jobs.
    Specify the following options in the `Makefile`:
    - `CREATE_JOB_TAG`: Name/tag of the job batch.
    - `CREATE_JOB_CONFIG`: `T3PS` MCMC configuration file.
    - `CREATE_JOB_nCores`: Number of cores to use.
    - `CREATE_JOB_program`: The path to the Paramater Point Processor binary.
    - `CREATE_JOB_nJobs`: Total number of jobs to submit
    - `CREATE_JOB_chain_length`: Length of the chains.
    - `CREATE_JOB_TEMPLATE`: Path to the `T3PS` job template.
    
    Having done this issue the command:
    ~~~~
    make create-jobs
    ~~~~

2. Submit cluster jobs.
    


------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

## Cross-section computation scan


### Overview

Computation of cross sections with `SusHi`.
- gg -> A (NNLO QCD)
- gg -> H (NNLO QCD)

The computational speed, (processed points)/(second) is much slower in this step



### Job submission

