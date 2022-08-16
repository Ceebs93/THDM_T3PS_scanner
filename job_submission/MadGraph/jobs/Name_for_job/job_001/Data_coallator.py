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

# This file inherits variables from looper.sh file as well as create-jobs.s

store_path = str(sys.argv[1]) #Folder containing all the runs data
print("store_path ", store_path)
process = str(sys.argv[2])
print("process ", process)
print("search term for file_paths", str(store_path) + str(process + "*.csv*"))
file_paths = glob.glob(str(store_path) + "/" + str(process) + "*.csv") # Creates list of all the MadGraph runs to use
print("file_paths", file_paths)
csv_path = str(store_path)+ "/" + str(process) + ".csv" #path for csv file to be created
print("csv_path ", csv_path)
#########################################################################################################
def combine_csvs(Scsvpaths, Savepath):
    """ retrieves all the csv files from file folders and creates a list"""
    
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
                            
###########################################################################################################

combine_csvs(file_paths, csv_path)
