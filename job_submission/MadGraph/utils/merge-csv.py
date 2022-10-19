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
FINAL_CSV = str(sys.argv[3]) + ".csv"
input_sin = str(sys.argv[4])
OG_sin = str(sys.argv[5])

##############################################################################
def add_xsect(Filename, Chkdfile, sin_in, sin_OG):
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

    sin_in : string
        This is the column label for the sin(beta-alpha) value for the Madgraph
        output that is being merged onto the pre-existing csv.

    sin_OG : string
        This is the column label for the sin(beta-alpha) value for the pre-exi-
        sting csv file onto which the cross-section values are being added.

    Returns
    -------
    Chkd_df : pandas dataframe
        Combination of the Filename and Chkdfile dataframes.

    """

    # Reading in the two csv files
    MG_df = pd.read_csv(Filename)
    print('MG_DF', MG_df.columns.tolist())
    Chkd_df = pd.read_csv(Chkdfile)
    print('Chkd_df', Chkd_df.columns.tolist())

    # Checking that the dataframes have the same number of rows before adding
    # cross-section column onto Chkd_df
    if len(MG_df.Sbma) == len(Chkd_df.sinba):
        if 'X_sections' in MG_df.columns:
            
            X_SECT_COL_ = MG_df['X_sections']


            Chkd_df['X_SECT_COL_'] = X_SECT_COL_
        else:
            print("Could not find a column 'X_sections' in provided file")
    else:
        print("Row lengths of CSVes do not match. Check for correct files")

    return Chkd_df
################################################################################


combi_df = add_xsect(madgraph_out, OG_csv, input_sin, OG_sin)

# These lines remove duplicates and any rows with missing entries
combi_df.dropna()
final_df = combi_df.drop_duplicates()
combi_df.info()
final_df = final_df.astype(float)

# Saving the final output as a csv to FINAL_CSV
final_df.to_csv(FINAL_CSV, index=False)
