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
        job). It should contain the same sin(beta-alpha), VAR1_(beta) points as
        'Filename' so that they can be combined.

    Returns
    -------
    data_base : pandas dataframe
        Combination of the Filename and Chkdfile dataframes.

    """

    # Reading in the two csv files and rounding everything to 6 decimal places
    MG_df = pd.read_csv(Filename)
    MG_df = MG_df.round(6)   
    data_base = pd.read_csv(Chkdfile)
    data_base = data_base.round(6)
    
    # Dropping any duplicate rows
    data_base = data_base.drop_duplicates()
    MG_df = MG_df.drop_duplicates()
       
    # Reading in multiple columns from the original csv file to lists 
    if NEW_COL_STR_ not in data_base.columns:
        L = len(data_base)
        new_col = np.zeroes(L)
        data_base[NEW_COL_STR_] = new_col

    else:
        continue
        
    main_in = data_base.[NEW_COL_STR_].tolist()
    main_VAR1_ = data_base.asVAR1_.tolist()
    main_VAR1_ = np.array(main_VAR1_)
    main_VAR1_ = np.around(main_VAR1_, decimals=4)

    main_VAR2_ = data_base.asVAR2_.tolist()
    main_VAR2_ = np.array(main_VAR2_)
    main_VAR2_ = np.around(main_VAR2_, decimals=5)

    main_VAR3_ = data_base.asVAR3_.tolist()
    main_VAR3_ = np.array(main_VAR3_)
    main_VAR3_ = np.around(main_VAR3_, decimals=3)

    main_VAR4_ = data_base.asVAR4_.tolist()
    main_VAR4_ = np.array(main_VAR4_)
    main_VAR4_ = np.around(main_VAR4_, decimals=3)
        
    # Reading in multiple columns from the new csv file to lists
    new_in = MG_df.X_sections.tolist()
    new_VAR1_ = MG_df.asVAR_1.tolist()
    new_VAR1_ = np.array(new_VAR1_)

    new_VAR1_ = np.around(new_VAR1_, decimals=4)
    new_VAR2_ = MG_df.asVAR2_.tolist()
    new_VAR2_ = np.array(new_VAR2_)

    new_VAR2_ = np.around(new_VAR2_, decimals=5)
    new_VAR3_ = MG_df.asVAR3_.tolist()
    new_VAR3_ = np.array(new_VAR3_)
    new_VAR3_ = np.around(new_VAR3_, decimals=3)
  
    new_VAR4_ = MG_df.asVAR4_.tolist()
    new_VAR4_ = np.array(new_VAR4_)
    new_VAR4_ = np.around(new_VAR4_, decimals=3)
    
    for i in range(0,len(new_in)):
        if new_VAR1_[i] in main_VAR1_:
            VAR1__index = np.where(main_VAR1_ == new_VAR1_[i])
            VAR1__index = VAR1__index[0][0]

            if new_VAR3_[i] in main_VAR3_:
                VAR3__index = np.where(main_VAR3_ == new_VAR3_[i])
                VAR3__index = VAR3__index[0][0]

                if VAR3__index != VAR1__index:
                    print("VAR3_ didn't match")
                    continue
                elif new_VAR2_[i] == main_VAR2_[VAR3__index]:

                    if new_VAR4_[i] == main_VAR4_[VAR3__index]:
                        main_in[VAR3__index] = new_in[i]
                        print("success!")
                    else:
                        print("VAR4_ didn't match")
                        continue
                else:
                    print("VAR3_ was correct but VAR2_ value wasn't")
                    continue
    
    data_base[NEW_COL_STR_] = main_in
    
    # These lines remove duplicates and any rows with missing entries
    combi_df.dropna()
    final_df = combi_df.drop_duplicates()
    combi_df.info()
    final_df = final_df.astype(float)

    return data_base
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
