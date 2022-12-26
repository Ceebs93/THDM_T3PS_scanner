#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: cb27g11
"""
# Data collator

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
#from builtins import *

import os
import shutil
import numpy as np
import pandas as pd
import sys
import glob

# This file inherits variables from looper.sh file as well as create-jobs.


#Folder containing all the runs data
store_path = str(sys.argv[1]) 
print("store_path ", store_path)
process = str(sys.argv[2])

print("search term for file_paths", str(store_path) + "/" + str(process + "*.csv*"))
# Creates list of all the MadGraph runs to use
file_paths = glob.glob(str(store_path) + "/" + str(process) + "*.csv")
print("file_paths", file_paths)
#path for csv file to be created
csv_path = str(store_path)+ "/" + str(process) + ".csv"
print("csv_path ", csv_path)

###############################################################################
def combine_csvs(Scsvpaths, Savepath):
    """
    Retrieves all the csv files from file folders and creates a list.

    Parameters
    ----------
    Scsv_paths : list of strings
        each item in the list is a path to a csv file that will be combined
        with the rest in the list.

    Savepath : string
        the path (and name) for the final csv file containing the combined
        input csvs.

    Returns
    -------
    None : Instead function directly writes new csv file to the specified
        location.
    """
    
    df_list = []

    for file in Scsvpaths:
        csv_loc = file
        read_in = pd.read_csv(str(csv_loc))
        df_list.append(read_in)

    combi_df = pd.concat(df_list, ignore_index = True)
    final_df = combi_df.drop_duplicates()
    combi_df.info()
    final_df = final_df.astype(float)
    print(final_df.head(5))

    final_df.to_csv(str(Savepath), index=False)
                            
###############################################################################

combine_csvs(file_paths, csv_path)
