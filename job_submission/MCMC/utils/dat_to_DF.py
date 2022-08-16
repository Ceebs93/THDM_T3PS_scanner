#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 19:34:32 2021

@author: cb27g11
"""

import numpy as np
import pandas as pd

#from kappaextractor import transkappa

################################################################################
def tan_trnsfm(dataf):
    """"calculates beta values and converts them to the lowest equivalent value.
This ensures we don't accidentally get any massively inflated cross-sections
back from MadGraph."""

    #Extracts tanbeta values and puts them in a list
    tanbeta = dataf['tb'].tolist()
    beta_list = []
    print("Starting beta transformations")
    for i in range(0, len(tanbeta)):
        # Beta values are calculated
        beta_val = np.arctan(tanbeta[i])
        # Checks if beta_val is too high, if it is then 2pi is subtracted until
        # it is no longer too high, this value is added to a "beta_list"
        if beta_val > 2*np.pi:
            while beta_val > 2*np.pi:
                beta_val = beta_val - 2*np.pi
        beta_list.append(beta_val)

    # List of converted tanbetas is created
    nw_tans = np.tan(np.array(beta_list))
    print("Completed tan transformations")

    dataf["tb"]=nw_tans
    print("tan values added to df")
    return dataf
################################################################################

################################################################################
def dat_to_DF(Filename, Csvname):
    """Reads in the dat file from the t3ps scanner and converts it to a csv
    file"""
    file1 = open(Filename, 'r')
    Lines = file1.readlines()
    lens = []
    #print(len(Lines[0]), len(Lines[276]))
    print(len(Lines))
    for i in range(0, len(Lines)):
        if len(Lines[i]) <= 700 :
            #or len(Lines[i])>= 800
            print(len(Lines[i]))
            lens.append(i)
    file1.close()
    print(len(lens))
    df = pd.read_csv(Filename, skiprows= lens, sep="\s+")
    print("Dataframe columns will be: " + str(df.columns.tolist()))
    df.dropna(inplace=True)

    print(df.head(5))

    df.to_csv(Csvname, index=False)
################################################################################

################################################################################
def add_xsect(Filename, Chkdfile, higgs, proc):
    """"Takes the name of an output CSV from the MG5 looper and the name of the
    file that was originally given to this and combines them to create a CSV
    with the original data as well as the MG5 calculated cross section"""

    #Adding .csv to end of filenames
    Filename = Filename + ".csv"
    Chkdfile = Chkdfile + ".csv"

    # Reading in the two csv files
    MG_df = pd.read_csv(Filename)
    Chkd_df = pd.read_csv(Chkdfile)
    print(type(Chkd_df))
    # Checking that the dataframes have the same number of rows before adding
    # cross section column onto Chkd_df
    if len(MG_df['Sbma']) == len(Chkd_df.sinba):
        if 'X_sections' in MG_df.columns:
            if int(higgs) == 1:
                x_column = MG_df['X_sections']
                Chkd_df = Chkd_df.assign(X_section_bq=x_column)
            elif int(higgs)==2:
                x_column = MG_df['X_sections']
                Chkd_df = Chkd_df.assign(X_section_bg=x_column)
            elif int(higgs)==3:
                x_column = MG_df['X_sections']
                Chkd_df = Chkd_df.assign(X_section_qq=x_column)
        # X_section = pd.Series(x_column)
        # print(type(X_section))
        # Chkd_df.insert(loc=0, column="X_section", value=X_section)
        # print(type(Chkdfile))
        else:
            print("Error! Could not find a column 'X_sections' in provided file")
    else:
        print("Error, row lengths of CSVes do not match. It is likely data "
              "has been lost somewhere, or you have paired up the wrong files!")

    return Chkd_df
################################################################################

################################################################################
def make_useful(Filename, inc_Hc, THDMC_output=False):
    """ takes data file, Filename, reads data into a pandas dataframe then 
    extracts columns mH, mHc, mA, cba, tb, sinba, m12_2, k_huu, and k_hdd; and 
    creates a new dataframe (new_df) from these. This is then sorted and values
    are set to be floats. THDMC_output indicates if the data in question is
    output from 2HDMC or not, if not, it is assumed to come from Magellan's 
    scanner. It will be treated as having one of these two formats.
    """

    print("Starting make_useful ...")

    # Reads in data from Filename

    if THDMC_output == True:
        df = pd.read_csv(Filename, sep=",")

        # Output from 2HDMC has a column called pass that indicates whether a
        # point has passed checks or not. We do not need this as by the time
        # the data runs through here it has already been filtered using this.
        if 'pass' in df:
            del df['pass']
        df = df.rename(columns={'mHp': 'mHc', 'tanb': 'tb'})
        if 'mh' in df:
            del df['mh']

    else:
        df = pd.read_csv(Filename, sep=",", names= ['Z7', 'mH', 'mHc', 'mA',
                         'cba', 'tb', 'sinba', 'Z4', 'Z5', 'm12_2',
                         'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7',
                         'g_HpHmh', 'Gamma_h', 'Gamma_H', 'Gamma_Hc', 'Gamma_A',
                         'br_h_bb', 'br_h_tautau', 'br_h_gg', 'br_h_WW', 'br_h_ZZ', 'br_h_gaga',
                         'br_A_tt', 'br_A_bb', 'br_A_gg', 'br_A_mumu', 'br_A_tautau', 'br_A_Zga', 'br_A_Zh',
                         'br_A_ZH', 'br_A_gaga',
                         'br_H_tt', 'br_H_bb', 'br_H_gg', 'br_H_mumu',
                         'br_H_tautau', 'br_H_Zga', 'br_H_Zh', 'br_H_WW', 'br_H_ZZ', 'br_H_ZA',
                         'br_H_hh', 'br_H_AA', 'br_H_gaga',
                         'br_Hp_tb', 'br_Hp_taunu', 'br_Hp_Wh', 'br_Hp_WH', 'br_Hp_WA',
                         'sta', 'uni', 'per_4pi', 'per_8pi',
                         'S', 'T', 'U', 'V', 'W', 'X',
                         'delta_rho', 'delta_amu', 'tot_hbobs', 'sens_ch',
                         'chi2_HS', 'chi2_ST_hepfit', 'chi2_ST_gfitter', 'chi2_Tot_hepfit', 'chi2_Tot_gfitter',
                         'k_huu', 'k_hdd',
                         'likelihood', 'stay_count', 'ratio'], skiprows=[0])

    temp_df = df.astype(float)

    # Drops all unwanted columns
    print("Columns in dataframe will be, " + str(temp_df.columns.tolist()))

    print("Inital dataframe created")

# Below we apply limits to our csv, done in two seperate chunks, the positive
# and negative values of sinba. These are then recombined.
    upper_df = temp_df[temp_df['sinba'] >= 0.88]
    upper_df = upper_df[upper_df['sinba'] <= 1.2]

    lower_df = temp_df[temp_df['sinba']<-0.7]
    lower_df = lower_df[lower_df['sinba']>=-1.1]

    temp_df = combiner([upper_df, lower_df])

    temp_df = temp_df[temp_df['mH'] < 1000]
    temp_df = temp_df[temp_df['mA'] < 1000]
    temp_df.reindex()

    df = tan_trnsfm(temp_df)

    print("about to drop NaNs")
    temp_df.dropna(inplace=True)

    # When excluding the charged higgs we require masses to be above 500GeV
    if inc_Hc == False:
        temp_df = temp_df.drop(temp_df[temp_df.mHc < 500].index)

    # Creates new dataframe to hold all the columns we are interested in and
    # adds these columns in
    print("About to create final df")
    # Dataframe created with 0 rows and columns from temp_df
    new_df = pd.DataFrame(columns = temp_df.columns.tolist())
    for i in range(temp_df.shape[0]):
        new_row = {}
        for j in range(temp_df.shape[1]):
            value = temp_df.iat[i, j]
            colname = temp_df.columns[j]
            new_row[colname] = float(value)
        new_df = new_df.append(new_row, ignore_index=True)


    print("About to deal with NaNs")
    for clmn in new_df.columns.tolist():

        nans = np.isnan(new_df[str(clmn)].tolist())
        if True in nans:
            indices = np.where(nans == True)
            indices = indices[0]

            for i in range(0, len(indices)):
                new_df = new_df.drop(indices[i])
            new_df.reset_index(drop = True, inplace = True)
            new_df = new_df.reset_index(drop = True)

    #ensures values are floats and resets index
    new_df = new_df.astype(float)
    new_df = new_df.reset_index(drop = True)

    return new_df
###############################################################################

################################################################################
def combiner(df_list):
    """ Takes a list of pandas dataframes, combines them and removes any duplicate
    rows"""


    # Combines the list of dataframes
    combi_df = pd.concat(df_list, ignore_index = True)
    combi_df.dropna()
    final_df = combi_df.drop_duplicates()
    combi_df.info()
    final_df = final_df.astype(float)

    return final_df
################################################################################

################################################################################
def alpha_beta(tanbeta, cosba):
    """Takes lists of tan(beta) and cos(beta-alpha) values and uses them to
    calculate the values of alpha and beta then return them as lists"""

    if len(tanbeta) != len(cosba):
        raise Exception("Lists are not the same length, there are missing values!")

    for i in range(0, len(tanbeta)):
        print("Calculating beta values")
        beta_val = np.arctan(tanbeta[i])
        # Condition below removes unphysical cos(beta-alpha) values
        if -1<=(cosba[i])<=1 :
            print("Calculating alpha values")
            alpha_val = -1*np.arcsin(cosba[i])+beta_val

    return(alpha_val, beta_val)
###############################################################################

###############################################################################
def user_interface():

    print("Do you wish to: 1) Process a file  or  2) Combine csvs?")
    option = int(input("Processing includes: conversion to csv, fixing tanbeta"
                   + "values, applying limits on variables and adding MG5 "
                   "results onto CSVs. "))

    if option == 1:
        print("What kind of file do you want to process?")
        print("1) A .dat file    2) A .csv file pre 2HDMC check    3) A .csv file"
              +" post 2HDMC check .csv    4) A .csv post MG5)")

        file_type = int(input("Please enter 1, 2, 3 or 4 only: "))
    
        print("What is the path to, and name of the file(If you are adding "
              "cross sections to a datafile please enter the name of the file "
              "containing these cross sections)?")
        file_name = input("You do not need to add .dat/.csv: ")
    
        if file_type == 1:
            nther_file = file_name + ".dat"
            dat_to_DF(nther_file, file_name + '.csv')
            print(file_name + ".dat has been converted to " + file_name + ".csv")
    
        elif file_type == 2:
            tn_fxd = make_useful(file_name + '.csv', inc_Hc=False, THDMC_output=False)
            print("The tanbeta values in " + file_name + ".csv have been fixed")
            print("Please enter name for final csv ")
            nw_name = input("You do not need to add .csv ")
            tn_fxd.to_csv(nw_name +'.csv', index=False)
    
        elif file_type == 3:
            tn_fxd = make_useful(file_name + '.csv', inc_Hc=False, THDMC_output=True)
            tn_fxd.to_csv(file_name +'.csv', index=False)
            print("The tanbeta values in " + file_name + ".csv have been fixed")

        elif file_type == 4:
            print("Please enter the name of the post 2HDMC CSV file you wish"
                  "to use: ")
            chkdfile = input("You do not need to add .csv ")
            process = input( "Enter the process; 1) bq  2) bg  3) qq or 4) define your own process/tag")
            if process == 4:
                process == input("Enter process/tag")
            higgses = input("Which Higgs/es are in process? Answer only; hi , h , H or A ")
            x_sect_df = add_xsect(file_name, chkdfile, higgses, process)
            print("Please enter the name for the final CSV file you wish"
                  " to use: ")
            endname = input("You do not need to add .csv ")
            x_sect_df.to_csv(endname + '.csv', index=False)

    elif option == 2:
        print("Please enter a list, separated by spaces only, of csvs to combine")
        list_input = input("You do not need to add .csv to the end of files ")
        # creates list of csv names from user input
        csv_list = list_input.split()
        # adds .csv to the end of each csv name
        csv_list = ["{}{}".format(i,".csv") for i in csv_list]
        # for loop creates list of loaded csvs
        loaded_lists = []
        for filename in csv_list:
            df = pd.read_csv(filename, index_col=None)
            loaded_lists.append(df)
        #combining dataframes
        frame = pd.concat(loaded_lists, axis=0, ignore_index=True)
        combined_csv = tan_trnsfm(frame)

        print("Please enter name for combined csv ")
        nw_name = input("You do not need to add .csv ")
        combined_csv.to_csv(str(nw_name) + ".csv", index=False)
###############################################################################

user_interface()
# The following two lines will read in a dat file from Magellan and save it as
# a csv suitable for the 2HDMC point checker

#csv_clipped = make_useful('testbatch_0510.csv', inc_Hc=False, THDMC_output=False)

#csv_clipped.to_csv('testbatch_0510.csv', index=False)

