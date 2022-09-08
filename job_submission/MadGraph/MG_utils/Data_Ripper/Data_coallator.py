# Data collator

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: cb27g11
"""

import os
import shutil
import numpy as np
import pandas as pd
import sys
import glob

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

#I inherit the following variables from create-jobs.sh:
# X_SECT_COL_

# This file will inherit values via 'sys.argv' from the corresponding
# looper.sh file

#Folder containing all the runs data
store_path = str(sys.argv[1]) 

process = str(sys.argv[2])

# Creates list of all the MadGraph runs to use
file_paths = glob.glob(str(store_path) + "/" + str(process) + "*.csv")

#path for csv file to be created
csv_path = str(store_path)+ "/" + str(process) + ".csv"


###############################################################################
def combine_csvs(Scsvpaths, Savepath):
    """
    Retrieves all the csv files from the directorys and creates a list.

    Parameters
    ----------
    Scsvpaths : list of strings
        each string in the list is a path and name of a csv file

    Savepath : string
        a path and name for the final writing of the combined csv file

    Returns
    -------
    None : combine_csvs does not directly return anything, it simply writes
        the new csv file to memory.

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
