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
    main_tb = data_base.tb.tolist()
    main_tb = np.array(main_tb)
    main_tb = np.around(main_tb, decimals=4)
    main_sinba = data_base.sinba.tolist()
    main_sinba = np.array(main_sinba)
    main_sinba = np.around(main_sinba, decimals=5)
    main_mH = data_base.mH.tolist()
    main_mH = np.array(main_mH)
    main_mH = np.around(main_mH, decimals=3)
    main_mHc = data_base.mHc.tolist()
    main_mHc = np.array(main_mHc)
    main_mHc = np.around(main_mHc, decimals=3)
        
    # Reading in multiple columns from the new csv file to lists
    new_in = MG_df.X_sections.tolist()
    new_tb = MG_df.Tb.tolist()
    new_tb = np.array(new_tb)
    new_tb = np.around(new_tb, decimals=4)
    new_sinba = MG_df.sinba.tolist()
    new_sinba = np.array(new_sinba)
    new_sinba = np.around(new_sinba, decimals=5)
    new_mH = MG_df.mH.tolist()
    new_mH = np.array(new_mH)
    new_mH = np.around(new_mH, decimals=3)
    new_mHc = MG_df.mHp.tolist()
    new_mHc = np.array(new_mHc)
    new_mHc = np.around(new_mHc, decimals=3)
    
    for i in range(0,len(new_in)):
        if new_tb[i] in main_tb:
            #print(new_tb[i])
            tan_index = np.where(main_tb == new_tb[i])
            tan_index = tan_index[0][0]

            if new_mH[i] in main_mH:
                #print(new_mH[i],"I found this mass!")
                mH_index = np.where(main_mH == new_mH[i])
                mH_index = mH_index[0][0]

                if mH_index != tan_index:
                    print("mass was for a different point")
                    continue
                elif new_sinba[i] == main_sinba[mH_index]:

                    if new_mHc[i] == main_mHc[mH_index]:
                        main_in[mH_index] = new_in[i]
                        print("success!")
                    else:
                        print("mHc didn't match")
                        continue
                else:
                    print("mass was correct but sin value wasn't")
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
