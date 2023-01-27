#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

# Data Ripper

@author: Ciara Byers
"""
#import glob
#import gzip
import os
import numpy as np
import pandas as pd
#import shutil
import tarfile
#from Benchmark_extractor import bench_compare
#from Data_manipulator import interpolate_gaps, imputer
from kappaconverter import transkappa
#from pandas import DataFrame as df

user_home = os.path.expanduser('~') #Self explanatory
#output_dir = user_home + "/outputs/NoCh_Allhi_hcm150_A250_hi50" #directory we start from
#paths = os.listdir(output_dir)# All the different runs we need to include
parent_dir = user_home + "/combo_wombo/PR_reruns_qqhi/"


##################################################################################
def get_uncomp_file(Parent_dir, subprocesses=True):
    """ retrieves all the uncompressed 'banner' files from file folders and
    creates a list"""

    print("inside get_uncomp_file")

    uncomprsd_files = os.listdir(Parent_dir)# either value folder, or process folder list

    if subprocesses == True:
        #Expect this to be a list of directories, e.g. 'total', 'process-h',
        #'process-H' etc
        process_dirs = os.listdir(Parent_dir)
        #print(process_dirs)
        uncomped_files = []
        #Now we go into each subprocess directory
        for process in process_dirs:
            #print("this is my process" + str(process))
            sub_proc_list = []
            #creating a list with the process as first entry
            sub_proc_list.append(process)
            sub_file = Parent_dir + "/" + process
            #print("this is my sub_file" + str(sub_file))
            path_to_runs = os.listdir(sub_file)
            # Finally we obtain the full path to the file
            for path in path_to_runs:
                path_to_file = sub_file + "/" + path + "/run_01_tag_1_banner.txt"
                if os.path.isfile(path_to_file):
                #Now we add the full paths to the process list
                    sub_proc_list.append(path_to_file)
                else:
                    sub_proc_list.append("empty")

            uncomped_files.append(sub_proc_list)


    else:
        uncomped_files = [] # list of zipped lhe files
        for folder in uncomprsd_files:
            path_to_file = Parent_dir + "/" + folder + "/run_01_tag_1_banner.txt"
            if os.path.isfile(path_to_file):
                uncomped_files.append(path_to_file)# adding the path to the copied
                #file to our list
            else:
                uncomped_files.append("empty")

    return uncomped_files
###################################################################################

###################################################################################
def untarrer(gzball, path, nw_path):
    """takes a list of compressed (tar.gz) files and outputs a list of
    uncompressed versions"""
    un_tarred = []
    #print(gzball)
    os.chdir(path)
    for file in gzball:
        if file == "empty":
            print('shit')
            un_tarred.append(file)
        else:
            untar_entry = tarfile.open(file)
            untar_entry.extractall(nw_path)
            untar_entry.close()

    return un_tarred
###################################################################################

###################################################################################
def variable_values(file_lists, Multi=True):
    """"takes values of tan(beta) and sin(beta-alpha) from data file_lists and
    returns them as lists"""
    tan_data_lists = []
    sin_data_lists = []

    if Multi == True:
        for a_list in file_lists:
            text_sin_data = []
            text_tan_data = []
            #print(text_sin_data)
            # selecting a list out of those provided
            #start from 2nd entry in list as first will be the name of the process
            #print(a_list[0])
            for a_file in a_list:
                # selecting an individual file from the list
                if os.path.isfile(a_file):
                    #checking the file exists
    
                    with open(a_file,'r') as data_file:
                        # open the current file in the list and goes on to read
                        #each line
                        lines = data_file.readlines()
                        for line in lines:
    
                            if '# tanbeta' in line:
                                # we find the line containing tanbeta and add the
                                #value to our list
                                temp_line_holder = line.split()
                                text_tan_data.append(temp_line_holder[1])
    
                            elif '# sinbma' in line:
                                # this works the same way as the previous part but
                                #for sin(beta-alpha)
                                temp_line_holder = line.split()
                                text_sin_data.append(temp_line_holder[1])
    
            tan_data_lists.append(text_tan_data)
            sin_data_lists.append(text_sin_data)

    else:
        for a_file in file_lists:
            # selecting an individual file from the list
            if os.path.isfile(a_file):
                #checking the file exists

                with open(a_file,'r') as data_file:
                    # open the current file in the list and goes on to read
                    #each line
                    lines = data_file.readlines()
                    for line in lines:

                        if '# tanbeta' in line:
                            # we find the line containing tanbeta and add the
                            #value to our list
                            temp_line_holder = line.split()
                            tan_data_lists.append(temp_line_holder[1])

                        elif '# sinbma' in line:
                            # this works the same way as the previous part but
                            #for sin(beta-alpha)
                            temp_line_holder = line.split()
                            sin_data_lists.append(temp_line_holder[1])


    return(tan_data_lists, sin_data_lists)
###################################################################################

##################################################################################
def Xtrctd_values(file_set, searchstring, pos, Multi=True):
    """"takes the cross-section and errors from files"""

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


##################################################################################
def Xsection_values(file_set, Multi=True):
    """"takes the cross-section and errors from files"""

    X_Sections = []

    if Multi == True:
        for a_list in file_set:
            # selects a list of files (i.e. from one process)
            X_Section = []
            for a_file in a_list:
                print(a_file)
                # selecting an individual file from the list
                if os.path.isfile(a_file):
                    with open(a_file,'r') as data_file:
                        # open the current file in the list and go on to read each line
                        lines = data_file.readlines()
                        for line in lines:
        
                            if '#  Integrated weight (pb)  :' in line:
                                # we find the line containing the cross-section and
                                #add the value to our list
                                temp_line_holder = line.split()
                                X_Section.append(temp_line_holder[5])
                            else:
                                continue
                        #print(X_Section)
    
                else:
                    continue

            #print(X_Section)
            #print(len(X_Sections))
        X_Sections.append(X_Section)

    else:
        for a_file in file_set:
            # selecting an individual file from the list
            if os.path.isfile(a_file):
                with open(a_file,'r') as data_file:
                    # open the current file in the list and go on to read each line
                    lines = data_file.readlines()
                    for line in lines:

                        if '#  Integrated weight (pb)  :' in line:
                            # we find the line containing the cross-section and
                            #add the value to our list
                            temp_line_holder = line.split()
                            X_Sections.append(temp_line_holder[5])
                        else:
                            continue


    return(X_Sections)
##################################################################################

#################################################################################
def dframe_maker(list_of_data, column_names):
    """Takes a list of data lists (i.e. tan values, cross-sections etc) as well
    as a list of column names and creates a pandas dataframe from these"""
    #Here we convert our data lists into a numpy array so that we can transpose
    #it into the arrangement we want for our dataframe.
    numpy_lists = np.array(list_of_data)
    transposed = numpy_lists.T
    transposed_lists = transposed.tolist()
    new_df = pd.DataFrame(transposed_lists, columns=column_names)

    print(new_df.head(5))

    return(new_df)
###################################################################################


# untarred_files = np.array(os.listdir(uncomp_dir))
# untarred_files = [uncomp_dir + "/" + s + "/run_01_tag_1_banner.txt" for s in
#                   untarred_files]

untarred_files = get_uncomp_file(parent_dir, subprocesses=False)

# creating lists of tan, sin and x-section values
#listylist =[ ['one',1, 2, 3], ['two', 4, 5, 6], ['three', 7, 8, 9]]
Tb_list, Sbma_list = variable_values(untarred_files, Multi=False)
H_mass = Xtrctd_values(untarred_files, "# mh2", 1, Multi=False)
A_mass = Xtrctd_values(untarred_files, "# mh3", 1, Multi=False)
Hp_mass = Xtrctd_values(untarred_files, "# mhc", 1, Multi=False)
X_sections = Xsection_values(untarred_files, Multi=False)
print(len(Tb_list))
print(len(X_sections))

Tb_list = np.array(Tb_list)
Sbma_list = np.array(Sbma_list)
X_sections = np.array(X_sections)
H_mass = np.array(H_mass)
A_mass = np.array(A_mass)
Hp_mass = np.array(Hp_mass)

Tb_list = Tb_list.astype(float)
Sbma_list = Sbma_list.astype(float)
X_sections = X_sections.astype(float)
H_mass = H_mass.astype(float)
A_mass = A_mass.astype(float)
Hp_mass = Hp_mass.astype(float)

data_lists = [X_sections, Tb_list, Sbma_list, H_mass, A_mass, Hp_mass]
columns = ["qq_tot_xsect", "Tb_list", "Sbma_list", "mH", "mA", "mHp"]

#Create a pandas dataframe to contain data

data_in_df = dframe_maker(data_lists, columns)

##############################################################################
def var_to_csv(tan_list, sin_list, xsections, higgs, name):
    """Takes numpy array of tan, sin and cross section values with a string
    indicating the desired higgs/es ( tot, h1, h2 or h2) to be used and 
    calculates the kappa and cos values for the data before saving them in a 
    csv file called name"""

    #As in the case that we consider the total cross-section we compare with
    # the SM by assuming we only have the SM-like higgs
    if str(higgs) == "tot":
        higgs = "h1"

  #  rdcd_xsects, k_tt, k_bb, k_ww, Cbma_list = transkappa(Tb_list, Sbma_list,
   #                                                       ("tt", "bb", "WW"), X_sections, (str(higgs)))

    data_in_df = dframe_maker(data_lists, columns)

    #Filter dataframe according to limits on k_tt1(or k_tot) (https://cds.cern.ch/record/2725523/files/HIG-19-008-pas.pdf)
  #  data_in_df = data_in_df[data_in_df['k_tt_tot']>=-0.9]

   # data_in_df = data_in_df[data_in_df['k_tt_tot']<=1.1]

    #lower_half = data_in_df[data_in_df['k_tt_tot']<=-0.7]

    #upper_half = data_in_df[data_in_df['k_tt_tot']>0.7]

    #if not lower_half.empty:
     #   if not upper_half.empty:
      #      data_in_df = pd.concat(upper_half, lower_half)
       # else:
        #    data_in_df = lower_half
   # elif not upper_half.empty:
    #    data_in_df = upper_half

    data_in_df.to_csv("/home/cb27g11/combo_wombo/" + str(name) + 'csv', index=False)

##############################################################################

#Sort dataframe
#data_in_df.sort_values(by='k_tt_tot', ascending=False, inplace=True)

#Filter dataframe according to limits on k_tt1(or k_tot) (https://cds.cern.ch/record/2725523/files/HIG-19-008-pas.pdf)
#data_in_df = data_in_df[data_in_df['k_tt_tot']>=-0.9]
#print(data_in_df)
#data_in_df = data_in_df[data_in_df['k_tt_tot']<=1.1]
#print(data_in_df)
#lower_half = data_in_df[data_in_df['k_tt_tot']<=-0.7]
#print(data_in_df)
#upper_half = data_in_df[data_in_df['k_tt_tot']>0.7]
#print(data_in_df)

# #check the two halves of the dataframe are non-empty and combime them accordingly
# #to create new version of data_in_df
#if not lower_half.empty:
#    if not upper_half.empty:
#        data_in_df = pd.concat(upper_half, lower_half)
 #   else:
  #      data_in_df = lower_half
#elif not upper_half.empty:
 #   data_in_df = upper_half

#print(data_in_df.head(10))

data_in_df.to_csv('/home/cb27g11/combo_wombo/PR_reruns_qqhi.csv', index=False)


