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

The work area of the MCMC scan jobs is found in [./job_submission/MCMC/](../job_submission/MCMC/).

The job preparation, submission and merging is handled by the [`Makefile`](../job_submission/MCMC/Makefile).

#### Instructions:

1. **Create cluster jobs.**

    Specify the following options in the `Makefile`:
    - `CREATE_JOB_TAG`: Name/tag of the job batch.
    - `CREATE_JOB_CONFIG`: `T3PS` MCMC configuration file.
    - `CREATE_JOB_nCores`: Number of cores to use.
    - `CREATE_JOB_program`: The path to the Paramater Point Processor binary.
    - `CREATE_JOB_nJobs`: Total number of jobs to submit
    - `CREATE_JOB_chain_length`: Length of the chains.
    - `CREATE_JOB_TEMPLATE`: Path to the `T3PS` job template.
    
    After you have cross-checked the specifications, you can prepare the jobs (job folders are created, and
    relevant config files are copied over) by issuing:

    ~~~~
    make create-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/MCMC/utils/create-jobs.sh`](../job_submission/MCMC/utils/create-jobs.sh).

2. **Submit cluster jobs.**

    Specify the following options in the `Makefile`:
    - `SUBMIT_JOB_TAG`: Name/tag of the to-be-submitted job batch.
    - `SUBMIT_JOB_LIST`: Name of the file in the `jobbs/$JOB_TAG` folder containing the list of job
    directories. The submission script loops over these.
    - `SUBMIT_JOB_RESOURCES`: PBS job resource specifications.
    - `SUBMIT_JOB_TASK`: Path to the job task script.

    After you have cross-checked the job submission specifications, you can submit the jobs by issuing:

    ~~~~
    make submit-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/MCMC/utils/submit-jobs.sh`](../job_submission/MCMC/utils/submit-jobs.sh).
    
3. **Merge cluster jobs.**

    After the cluster jobs have finished you can merge the individual chains together.

    Specify the following options in the `Makefile`:
    - `MERGE_JOB_TAG`: Name/tag of the job batch.
    - `MERGE_JOB_HEADER`: Header file by which the merged ASCII file is prepended.
    - `MERGE_JOB_CONVERT_ONLY`: If you would like to convert and already merged ASCII file, without re-merging the jobs use this option.
    - `MERGE_JOB_CONVERT`: Convert the ASCII file to HDF5.
    - `MERGE_JOB_COMPRESSION`: Compression library used. Recommended: `"blosc"`
    - `MERGE_JOB_DATASET_NAME`: Name of the dataset within HDF5. Do not start with numerical characters.
    - `MERGE_JOB_FORMAT`: Format of the HDF5 dataset. Recommended: `"table"`

    After you have cross-checked the merge specifications, you can merge the jobs by issuing:

    ~~~~
    make merge-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/MCMC/utils/merge-jobs.sh`](../job_submission/MCMC/utils/merge-jobs.sh).

#### Contents of the work area

**Directories:**
- `utils`: contains utility `shell` scripts which received input from the `Makefile` and help with
    tasks, such as job submission preparations, job submission, and merging of the jobs.
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
- `jobs`: folder contains the individual jobs in format of `JOB_TAG/job_XXX`. This directory contains:
    - `t3ps.conf`: the `T3PS` configuration of the submitted jobs
    - `job.template`: the job template needed by `T3PS`.


------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

## Cross-section computation scan


### Overview

Computation of cross sections with `SusHi`.
- gg -> A (NNLO QCD)
- gg -> H (NNLO QCD)

The computational speed i.e. (processed points)/(second) is much slower in this step, therefore this
is only applied to the points which survive the theoretical consistency conditions.


### Job submission

The work area of the MCMC scan jobs is found in [./job_submission/SusHi/](../job_submission/SusHi/).

The job preparation, submission and merging is handled by the [`Makefile`](../job_submission/SusHi/Makefile).


#### Instructions:

1. **Create cluster jobs.**

    Specify the following options in the `Makefile`:
    - `CREATE_JOB_TAG`: Name/tag of the job batch.
    - `CREATE_JOB_INPUT_DATA`: Path to the input parameter points which are to be appended by the cross section values.
        This is on default the merged output of the MCMC scan.
    - `CREATE_JOB_CONFIG`: `T3PS` MCMC configuration file.
    - `CREATE_JOB_nCores`: Number of cores to use.
    - `CREATE_JOB_nJobs`: Total number of jobs to submit
    - `CREATE_JOB_program`: The path to the program which is called by T3PS to process a single point.
        The current default program is
        [`/job_submission/SusHi/template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha`](../job_submission/SusHi/template/SusHi_2HDMC_hybrid_cba_template_A_and_H.slha),
        which calculates both &sigma;(pp -> H) and &sigma;(pp -> A) for the same point.
    - `CREATE_JOB_TEMPLATE`: Path to the `T3PS` job template.
    
    After you have cross-checked the specifications, you can prepare the jobs (job folders are created, and
    relevant config files are copied over) by issuing:

    ~~~~
    make submit-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/SusHi/utils/create-jobs.sh`](../job_submission/SusHi/utils/create-jobs.sh).


2. **Submit cluster jobs.**

    Specify the following options in the `Makefile`:
    - `SUBMIT_JOB_TAG`: Name/tag of the to-be-submitted job batch.
    - `SUBMIT_JOB_LIST`: Name of the file in the `jobbs/$JOB_TAG` folder containing the list of job
    directories. The submission script loops over these.
    - `SUBMIT_JOB_RESOURCES`: PBS job resource specifications.
    - `SUBMIT_JOB_TASK`: Path to the job task script.

    After you have cross-checked the job submission specifications, you can submit the jobs by issuing:

    ~~~~
    make submit-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/SusHi/utils/submit-jobs.sh`](../job_submission/SusHi/utils/submit-jobs.sh).

3. **Merge cluster jobs.**

    After the cluster jobs have finished you can merge the individual jobs together.

    Specify the following options in the `Makefile`:
    - `MERGE_JOB_TAG`: Name/tag of the job batch.
    - `MERGE_JOB_HEADER`: Header file by which the merged ASCII file is prepended.

    After you have cross-checked the merge specifications, you can merge the jobs by issuing:

    ~~~~
    make merge-jobs
    ~~~~

    which exports the relevant Makefile variables and calls
    [`job_submission/SusHi/utils/merge-jobs.sh`](../job_submission/SusHi/utils/merge-jobs.sh).

