#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: cb27g11
"""
# Data Ripper
import os
import shutil
import numpy as np
import sys
import pandas as pd

# Note: This version of Data_Ripper_iridis is set to look for tan(beta) values 
# by scanning for\'# tanbeta' in the .lhe file. THis is the setup for the
# type-II model in use, it will NOT work for other models as they may use a
# different name for the variable.

#I inherit the following variables from merge-jobs.sh:
# JOB_DIR_

# This file will inherit values via 'sys.argv' from the corresponding
# looper.sh file

Process = str(sys.argv[1]) #Name of results folder
uncomp_dir = "JOB_DIR_/" + str(Process) #Folder containing all the run data
output_dir = "JOB_DIR_" + "/Data_Files"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

###############################################################################
def get_uncomp_file(Uncomp_dir):
    """
    Retrieves all the banner.txt/unweighted_events files from file folders and
creates a list.

    Parameters
    ----------
    Uncomp_dir : string
        path to the directory containing files user wishes to extract data from

    Returns
    -------
    uncomped_files : list
        a list of the paths to the various "*banner.txt" files new location 
        after being moved by this function.
    """

    uncomped_files = [] # list of zipped lhe files
    print("Data_Ripper is looking in the following places: ")
    print(uncomp_dir)
    print("                                                      ")
    path_to_file = uncomp_dir + "/Events/run_01/run_01_tag_1_banner.txt"
    path_to_lhe = uncomp_dir + "/Events/run_01/unweighted_events.lhe.gz"

    if os.path.exists(path_to_file) and os.path.exists(path_to_lhe) :
        print("Banner-path found")
        # name of new folder for storing
        file_storage = "JOB_DIR_" + "/Processed_results/" + str(Process)
        # storing the data for later use
        if not os.path.exists(file_storage):
            os.makedirs(file_storage)
        
        shutil.copy(path_to_file, file_storage)# copying the file
        shutil.copy(path_to_lhe, file_storage)# moving the lhe file
        txt_stored = str(file_storage) + "/run_01_tag_1_banner.txt"
        # adding the path for the copied file to our list
        uncomped_files.append(txt_stored)
    else:
        uncomped_files.append("empty")
    return uncomped_files                
###############################################################################    

###############################################################################
def variable_values(* args):
    """
    Extracts values of tan(beta) and sin(beta-alpha) from banner.txt files and 
    returns them in lists

    Parameters
    ----------
    * args : list
        a list of paths to files user wishes to extract data from

    Returns
    -------
    text_tan_data : list
        a list of tan values extracted from MadGraph output
    text_sin_data : list
        a list of sin values extracted from MadGraph output
    """
    text_tan_data = []
    text_sin_data = []

    for a_list in args:
        for a_file in a_list:
            if os.path.isfile(a_file):
                with open(a_file,'r') as data_file:
                    # open the current file in the list and read each line
                    lines = data_file.readlines()
                    for line in lines:
                        if '# tanbeta' in line:
                            print("Data_Ripper has found tanbeta")
                            # find line with tanbeta and add value to list
                            temp_line_holder = line.split()
                            print("Found tanbeta=" +  str(temp_line_holder[1]))
                            text_tan_data.append(temp_line_holder[1])
                    
                        elif '# sinbma' in line:
                            # works as previous part did but for sin(beta-alpha)
                            temp_line_holder = line.split()
			    print("Found sinbma=" + str(temp_line_holder[1]))
                            text_sin_data.append(temp_line_holder[1])
                        else:
                            continue
            else:
                continue

        return(text_tan_data, text_sin_data)	 
###############################################################################

###############################################################################
def Xsection_values(file_set):
    """
    Takes the cross-section and errors from MadGraph output files
    
    Parameters
    ----------
    file_set : list
        a list of paths to data files the user wishes to extract values from

    Returns
    -------
    X_section : list
        a list of cross-section values in pb (as this is MadGraph's output
        unit)
    """

    X_Section = list()

    for a_file in file_set:
        if os.path.isfile(a_file):
            with open(a_file,'r') as data_file:
                # open the current file in the list and read each line
                lines = data_file.readlines()
                for line in lines:
                    if '#  Integrated weight (pb)  :' in line:
                        # find line with cross-section and add value to list
                        temp_line_holder = line.split()
			print("Found x_section=" + str(temp_line_holder[5]))
                        X_Section.append(temp_line_holder[5])
                    else:
                        continue
        
    return(X_Section)
###############################################################################    
    
File_List = get_uncomp_file(uncomp_dir)
print("Data_Ripper will look for sin and tan values in: " + str(uncomp_dir))
print("                                     ")
Tb_list, Sbma_list = variable_values(File_List)
X_sections = Xsection_values(File_List)

data = [Tb_list, Sbma_list, [X_sections[0]]]
col_names = ["Tb", "Sbma", "X_section"]

nw_df = pd.DataFrame({'Tb':Tb_list, 'Sbma':Sbma_list, 'X_sections':X_sections})
nw_df.to_csv(str(output_dir) + "/" + str(Process) + ".csv", index=False)

print(len(Sbma_list), len(Tb_list), len(X_sections))
