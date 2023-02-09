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
proc_out_dir = "JOB_DIR_/" + str(Process) #Folder containing all the run data

output_dir = "JOB_DIR_" + "/Data_Files"
paths = os.listdir(proc_out_dir)# All the different runs we need to include
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
    print(proc_out_dir)
    print("                                                      ")
    path_to_file = proc_out_dir + "/Events/run_01/run_01_tag_1_banner.txt"
    path_to_lhe = proc_out_dir + "/Events/run_01/unweighted_events.lhe.gz"

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

##################################################################################
def Xtrctd_values(file_set, searchstring, pos, Multi=True):
    """
    Extracts values of tan(beta) and sin(beta-alpha) from banner.txt files and 
    returns them in lists

    Parameters
    ----------
    file_set : list
        a list of paths to files user wishes to extract data from.

    searchstring : string
        a string containing a key element from the line in the 'banner' file
        that the variable value should appear in.

    pos : integer
        the position that your variable value will be read from when the line
        is read in as a list. E.g. " This is 1 line" has 4 entries, 'this', 'is'
        , '1' and finally 'line'. If you wanted 'is' you should provide pos = 1

    Multi : Boolean
        Indicates whether Data_Ripper is acting on a single file or a list of
        them (equally it will work if your target file is within a directory and 
        you provide the path to the parent directory of this)

    Returns
    -------
    xtrctd_vals : list
        a list of values for your variable extracted from MadGraph output
    """

    xtrctd_vals = []

    if Multi == True:
        for a_list in file_set:
            # selects a list of files (i.e. from one process)
            xtrctd_vals = []
            for a_file in a_list:
                print(a_file)
                # selecting an individual file from the list
                if os.path.isfile(a_file):
                    with open(a_file,'r') as data_file:
                        # open the current file in the list and go on to read each line
                        lines = data_file.readlines()
                        for line in lines:
        
                            if searchstring in line:
                                # we find the line containing the cross-section and
                                #add the value to our list
                                temp_line_holder = line.split()
                                xtrctd_vals.append(temp_line_holder[pos])
                            else:
                                continue
                        #print(X_Section)
    
                else:
                    continue

            #print(X_Section)
            #print(len(X_Sections))
        xtrctd_vals.append(xtrctd_vals)

    else:
        for a_file in file_set:
            # selecting an individual file from the list
            if os.path.isfile(a_file):
                with open(a_file,'r') as data_file:
                    # open the current file in the list and go on to read each line
                    lines = data_file.readlines()
                    for line in lines:

                        if searchstring in line:
                            # we find the line containing the cross-section and
                            #add the value to our list
                            temp_line_holder = line.split()
                            xtrctd_vals.append(temp_line_holder[pos])
                        else:
                            continue


    return(xtrctd_vals)
##################################################################################
    
File_List = get_uncomp_file(proc_out_dir)

print("Data_Ripper will look for sin and tan values in: " + str(proc_out_dir))
print("                                     ")

# The strings here will be used to find the values for each variable, make sure
# to change them if you are looking for something else
Tb_list = Xtrctd_values(File_List, '# tanbeta', 1, Multi=False)
Sbma_list = Xtrctd_values(File_List, '# sinbma', 1, Multi=False)
H_mass = Xtrctd_values(File_List, "# mh2", 1, Multi=False)
A_mass = Xtrctd_values(File_List, "# mh3", 1, Multi=False)
Hp_mass = Xtrctd_values(File_List, "# mhc", 1, Multi=False)
X_sections = Xsection_values(File_List, '#  Integrated weight (pb)  :', 5,
             Multi=False)

# Changing our variable lists to arrays
tb = np.array(Tb_list)
sinba = np.array(Sbma_list)
X_sections = np.array(X_sections)
H_mass = np.array(H_mass)
A_mass = np.array(A_mass)
Hp_mass = np.array(Hp_mass)

# Ensuring all values have been correctly read in as floats not strings
tb = tb.astype(float)
sinba = sinba.astype(float)
X_sections = X_sections.astype(float)
H_mass = H_mass.astype(float)
A_mass = A_mass.astype(float)
Hp_mass = Hp_mass.astype(float)

# Create a pandas dataframe to contain data
nw_df = pd.DataFrame({'Tb':tb, 'sinba':sinba, 'X_sections':X_sections,
                      'mH':H_mass, 'mA':A_mass, 'mHp':Hp_mass})

# Writing dataframe to csv
nw_df.to_csv(str(output_dir) + "/" + str(Process) + ".csv", index=False)
print(nw_df)
