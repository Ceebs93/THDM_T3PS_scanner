#!/usr/bin/env python2
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
print('FINAL_CSV', FINAL_CSV)

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

    Returns
    -------
    OG_csv_df : pandas dataframe
        Combination of the Filename and Chkdfile dataframes.

    """

    # Reading in the two csv files
    MG_df = pd.read_csv(Filename)
    print('MG_DF', MG_df.columns.tolist())
    OG_csv_df = pd.read_csv(Chkdfile)
    print('OG_csv_df', OG_csv_df.columns.tolist())

    # Checking that the dataframes have the same number of rows before adding
    # cross-section column onto OG_csv_df
    if len(MG_df['MG_SIN_LABEL_']) == len(OG_csv_df['OG_SIN_LABEL_']):
        MG_df.sort_values(by=['MG_SIN_LABEL_'], inplace=True)
        OG_csv_df.sort_values(by=['OG_SIN_LABEL_'], inplace=True)

        MG_sin_list = MG_df['MG_SIN_LABEL_'].to_list()
        OG_sin_list = OG_csv_df['OG_SIN_LABEL_'].to_list()
        print(MG_df['X_sections'])
        
        # Checking to ensure sin values match up
        if round(MG_sin_list[0], 5) == round(OG_sin_list[0], 5) and\
        round(MG_sin_list[-1], 5) == round(OG_sin_list[-1], 5):
            print('sin match confirmed')
            
            MG_df.sort_values(by=['MG_TAN_LABEL_'], inplace=True)
            OG_csv_df.sort_values(by=['OG_TAN_LABEL_'], inplace=True)

            MG_tan_list = MG_df['MG_TAN_LABEL_'].to_list()
            OG_tan_list = OG_csv_df['OG_TAN_LABEL_'].to_list()
            
            midpoint = int(len(MG_tan_list)/2)
            # Checking to ensure tan values also match up, including at one
            # different point than used for sin
            if round(MG_tan_list[0], 5) == round(OG_tan_list[0], 5) and\
            round(MG_tan_list[midpoint], 5) == round(OG_tan_list[midpoint], 5):

                print('tan match confirmed')
                if 'X_sections' in MG_df.columns:
                    print('X_sections located')            
                    X_SECT_COL_ = MG_df['X_sections'].tolist()

                    print(X_SECT_COL_[0:5])
                    OG_csv_df['X_SECT_COL_'] = X_SECT_COL_
                else:
                    print("Could not find a column 'X_sections' in provided \
                          file")
    else:
        print("Row lengths of CSVes do not match. Check for correct files")

    return OG_csv_df
###############################################################################
###############################################################################
def round_no(x, sf):

    """Simple function to round a float to a given number of significant figures.
    Does not work for negative numbers.

    Parameters
    ----------
    x : float
        The float the user wishes to round

    sf : int
        The number of significant figures the user wishes to round to

    Returns
    -------
    rounded : float
        The rounded input number 
    """
    rounded = round(x, -int(floor(log10(abs(x)))))

    return rounded
###############################################################################
###############################################################################
def round_to_sf(num, sf):

    """Simple function to round a float to a given number of significant figures.
    Does not work for negative numbers.

    Parameters
    ----------
    num : float/list of floats
        The float/s the user wishes to round

    sf : int
        The number of significant figures the user wishes to round to

    Returns
    -------
    rounded_no : float/list of floats
        The rounded input numbers 
    """

    if type(num) == float:
        rounded_no = round_no(num, sf)

    else:
        rounded_no = []
        for i in range(0,len(num)):
            rounded_no.append(round_no(num[i]),sf)

    return rounded_no
###############################################################################

combi_df = add_xsect(madgraph_out, OG_csv)

# These lines remove duplicates and any rows with missing entries
combi_df.dropna()
final_df = combi_df.drop_duplicates()
combi_df.info()
final_df = final_df.astype(float)

# Saving the final output as a csv to FINAL_CSV
final_df.to_csv(FINAL_CSV, index=False)
