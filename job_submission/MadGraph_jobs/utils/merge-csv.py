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
# X_SECT_COL_ MG_VAR1_LABEL_  OG_VAR1_LABEL_  MG_VAR2_LABEL_  OG_VAR2_LABEL_ 
# MG_VAR3_LABEL_  OG_VAR3_LABEL_ 

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
    final_df : pandas dataframe
        Combination of the Filename and Chkdfile dataframes.

    """

    # Reading in the two csv files and rounding everything to 6 decimal places
    MG_df = pd.read_csv(Filename)
    MG_df = MG_df.round(6)  
    MG_df.describe 
    data_base = pd.read_csv(Chkdfile)
    data_base = data_base.round(6)
    print("                       ")
    
    
    # Dropping any duplicate rows
    data_base = data_base.drop_duplicates()
    MG_df = MG_df.drop_duplicates()
       
    # Reading in multiple columns from the original csv file to lists 
    if "X_SECT_COL_" not in data_base.columns:
        L = len(data_base)
        new_col = np.zeroes(L)
        data_base["X_SECT_COL_"] = new_col
        
    main_in = data_base["X_SECT_COL_"].tolist()
    print(main_in)
    print("OG_VAR1_LABEL_", "OG_VAR2_LABEL_", "OG_VAR3_LABEL_", "OG_VAR4_LABEL_")
    main_OG_VAR1_LABEL_ = data_base.OG_VAR1_LABEL_.tolist()
    print("main_OG_VAR1_LABEL_ = data_base.OG_VAR1_LABEL_.tolist()")
    main_OG_VAR1_LABEL_ = np.array("main_OG_VAR1_LABEL_")
    main_OG_VAR2_LABEL_ = data_base.OG_VAR2_LABEL_.tolist()
    main_OG_VAR2_LABEL_ = np.array("main_OG_VAR2_LABEL_")
    main_OG_VAR3_LABEL_ = data_base.OG_VAR3_LABEL_.tolist()
    main_OG_VAR3_LABEL_ = np.array("main_OG_VAR3_LABEL_")
    main_OG_VAR4_LABEL_ = data_base.OG_VAR4_LABEL_.tolist()
    main_OG_VAR4_LABEL_ = np.array("main_OG_VAR4_LABEL_")
        
    # Reading in multiple columns from the new csv file to lists
    new_in = MG_df.X_sections.tolist()
    print(new_in)
    print("MG_VAR1_LABEL_","MG_VAR2_LABEL_","MG_VAR3_LABEL_","MG_VAR4_LABEL_")
    new_MG_VAR1_LABEL_ = MG_df.MG_VAR1_LABEL_.tolist()
    new_MG_VAR1_LABEL_ = np.array(new_MG_VAR1_LABEL_ )
    new_MG_VAR2_LABEL_ = MG_df.MG_VAR2_LABEL_.tolist()
    new_MG_VAR2_LABEL_ = np.array(new_MG_VAR2_LABEL_)
    new_MG_VAR3_LABEL_ = MG_df.MG_VAR3_LABEL_.tolist()
    new_MG_VAR3_LABEL_ = np.array(new_MG_VAR3_LABEL_)
    new_MG_VAR4_LABEL_ = MG_df.MG_VAR4_LABEL_.tolist()
    new_MG_VAR4_LABEL_ = np.array(new_MG_VAR4_LABEL_)
    
    for i in range(0,len(new_in)):
        if new_MG_VAR1_LABEL_[i] in main_OG_VAR1_LABEL_:
            #print(new_MG_VAR1_LABEL[i])
            var1_index = np.where(main_OG_VAR1_LABEL_ == new_MG_VAR1_LABEL_[i])
            var1_index = var1_index[0][0]

            if new_MG_VAR3_LABEL_[i] in main_OG_VAR3_LABEL_:
                #print(new_MG_VAR3_LABEL_[i],"I found this mass!")
                var3_index = np.where(main_OG_VAR3_LABEL_ == new_MG_VAR3_LABEL_[i])
                var3_index = var3_index[0][0]

                if var3_index != var1_index:
                    print("mass was for a different point")
                    continue
                elif new_MG_VAR2_LABEL_[i] == main_OG_VAR2_LABEL_[var3_index]:

                    if new_MG_VAR4_LABEL_[i] == main_OG_VAR4_LABEL_[var3_index]:
                        main_in[var3_index] = new_in[i]
                        print("success!")
                    else:
                        print("var4 index didn't match")
                        continue
                else:
                    print("var3 was correct but var2 wasn't")
                    continue
    
    data_base["X_SECT_COL_"] = main_in
    
    # These lines remove duplicates and any rows with missing entries
    data_base.dropna()
    final_df = data_base.drop_duplicates()
    final_df.info()
    final_df.describe()
    final_df = final_df.astype(float)

    return final_df
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
combi_df.to_csv(FINAL_CSV, index=False)
