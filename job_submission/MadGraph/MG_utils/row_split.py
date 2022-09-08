#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 16:27:12 2022

@author: cb27g11
"""

import pandas as pd
import math

#I inherit the following variables from merge-jobs.sh:
# INPUT_DATA_, nJobs_, datafile_

###############################################################################
def row_dropper(filepath, x, newfile):
    """
    Reads in a csv files and drops x rows from the top.
    
    Parameters
    ----------
    filepath : string
        the name and path to the csv file the user wishes to use.

    x : int
        the number of rows to be removed from the top of the csv file

    newfile : string
        the name and path the user wishes to save the new csv file to

    Returns
    -------
    None : saves new csv to 'newfile' but does not return anything.
    """

    old_df = pd.read_csv(filepath)
    print(old_df.tail(5))
    to_drop = list(range(0, x))
    new_df = old_df.drop(old_df.index[to_drop])
    print(new_df.head(5))

    new_df.to_csv(newfile, index=False)
###############################################################################

##############################################################################
def csv_splitter(filepath, x, newnametag):
    """
    Reads in a csv located at 'filepath', splits this into 'x' csves with sizes
    as close to even as is possible. Names are assigned as 'newnametag'_00,
    'newnametag'_01 and so on.
    
    Parameters
    ----------
    filepath : string
        the name and path to the csv file the user wishes to split.

    x : int
        the number of files to split the initial csv file into.

    newnametag : string
        the base name to give to each of the newly created csv files. These
        will have the form 'newnametag'_** where the asterisks are unique
        numbers.

    Returns
    -------
    None : function creates and saves multiple new csv files but does not 
        directly return these.

    """
    old_df = pd.read_csv(filepath)
    csv_len = float(len(old_df))
    x = float(x)
    print(csv_len/x)
    print(5.0/3.0)
    csv_size = int(math.ceil(csv_len/x))

    print("csv_size", csv_size)

    print("x", x)

    i = 1

    while i !=x:
        print("Loop number: " + str(i))
        if i != x:

            print(" creating csv number: " + str(i))
            selected_rows = list(range(0, csv_size))
            print("Selected rows: " + str(selected_rows))
            
            new_df = old_df.iloc[selected_rows]
            print("new_df", new_df)
            
            print("Saving new_df...")
            new_df.to_csv(str(newnametag) + "_" + str(i) + ".csv", index=False)
            
            print("dropping saved rows from old_df. Before:  ")
            print(old_df)
            old_df = old_df.drop(old_df.index[selected_rows])
            print("old_df", old_df)

            i += 1

    if i == x:
        old_df.to_csv(str(newnametag) + "_" + str(i) + ".csv", index=False)


csv_splitter("INPUT_DATA_", nJobs_, "datafile_")



