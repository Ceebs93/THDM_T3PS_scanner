#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:03:32 2022

@author: cb27g11
"""

import numpy as np
import pandas as pd
import os
import sys

from scipy import stats

#I inherit the following variables from merge-jobs.sh:
# X_SECT_COL_

# This file will inherit values via 'sys.argv' from the corresponding
# looper.sh file

#location of combined madgraph output
madgraph_out = str(sys.argv[1])
#location of the original csv file prior to splitting
OG_csv = str(sys.argv[2])
FINAL_CSV = str(sys.argv[3])

##############################################################################
def add_xsect(Filename, Chkdfile):
    """
    Takes the name of an output CSV from the MG5 looper and the name of the
    file that was originally given to this and combines them to create a CSV
    with the original data as well as the MG5 calculated cross section

    Parameters
    ----------
    Filename : string
        This is the name (and path to) a csv file from running Madgraph. Note
        that one should not add '.csv' to the end of this as it is added by 
        the function.

    Chkdfile : string
        This is the name (and path to) the original csv file that was fed to
        MadGraph (prior to any splitting being done in setting up a MadGraph
        job). It should contain the same sin(beta-alpha), tan(beta) points as
        'Filename' so that they can be combined.

    proc : string
        This is used to name the final output csv, it is intended to be the 
        process that was fed to MadGraph to produce 'Filename' i.e.
        bq_tqh1 to indicate the process b q > t q h1 .

    Returns
    -------
    Chkd_df : pandas dataframe
        Combination of the Filename and Chkdfile dataframes.

    """

    # Reading in the two csv files
    MG_df = pd.read_csv(Filename)
    Chkd_df = pd.read_csv(Chkdfile)

    # Checking that the dataframes have the same number of rows before adding
    # cross-section column onto Chkd_df
    if len(MG_df['Sbma']) == len(Chkd_df.sinba):
        if 'X_sections' in MG_df.columns:
            
            bg_twh1 = MG_df['X_sections']


            Chkd_df['bg_twh1'] = bg_twh1
        else:
            print("Could not find a column 'X_sections' in provided file")
    else:
        print("Row lengths of CSVes do not match. Check for correct files")

    return Chkd_df
################################################################################


combi_df = add_xsect(madgraph_out, OG_csv)

# These lines remove duplicates and any rows with missing entries
combi_df.dropna()
final_df = combi_df.drop_duplicates()
combi_df.info()
final_df = final_df.astype(float)

# Saving the final output as a csv to FINAL_CSV
final_df.to_csv(FINAL_CSV)
