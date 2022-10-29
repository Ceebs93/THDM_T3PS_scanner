#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 19:34:32 2021

@author: cb27g11
"""

import numpy as np
import pandas as pd
import sys

from scipy import stats


################################################################################
def tan_trnsfm(dataf):
    """
    Calculates beta values and converts them to the lowest equivalent value.
This ensures we don't accidentally get any massively inflated cross-sections
back from MadGraph.
    
    Parameters
    ----------
    dataf : pandas dataframe
        a dataframe containing a column of 'tanbeta' values to be transformed

    Returns
    -------
    dataf : pandas dataframe
        the given dataframe but with the 'tanbeta' values transformed to the 
        lowest valid values. i.e. we wish for the values for the 'betas' to
        be < 2*pi 
    
    """

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
    """
    Reads in the dat file from the t3ps scanner and converts it to a csv
    file

    Parameters
    ----------
    Filename : string
        This is the name (and path to) a .dat file the user wishes to read into
        a csv file.

    Csvname : string
        This is the name (and path for) the output csv the user wishes to use.

    Returns
    -------
    None : writes the new csv file to the given path.

    """

    file1 = open(Filename, 'r')
    Lines = file1.readlines()
    lens = []

    for i in range(0, len(Lines)):
        if len(Lines[i]) <= 700 :
            #or len(Lines[i])>= 800
            print(len(Lines[i]))
            lens.append(i)
    file1.close()
    df = pd.read_csv(Filename, skiprows= lens, sep="\s+")

    print("Dataframe columns will be: " + str(df.columns.tolist()))
    df.dropna(inplace=True)

    print(df.head(5))

    df.to_csv(Csvname, index=False)
################################################################################


################################################################################
def make_useful(Filename, inc_Hc, col_names=names):
    """
    Takes data file, Filename, reads data into a pandas dataframe then extracts
     columns mH, mHc, mA, cba, tb, sinba, m12_2, k_huu, and k_hdd; and creates 
    a new dataframe (new_df) from these based on filtering out points which are
    improbable or violate EWPO constraints. The user can find the arxiv numbers
    of the papers from which these constraints are taken within the function 
    itself as well as the manual. This is then sorted and values are set to be 
    floats. Function assumes that the input file is in the form of a .dat
    file and that it has the same number of columns as the global variable, 
    "names", by default; however the user may set this to any list that suits
    their data better.

    Parameters
    ----------
    Filename : string
        the name of the .dat file the user wishes to filter and transform into 
        a .csv file.

    inc_Hc : boolean
        if False the function will remove all points where the mass of the 
        charged Higgs is below a certain value.

    col_names : list
        a list of column names to be given to the created dataframe. By default
        this will be set to the global variable "names".
        
    Returns
    -------
    new_df : pandas dataframe
        This is a new pandas dataframe containing the filtered data from the
        input file.
    """

    print("Starting make_useful ...")

    # Reads in data from Filename
    df = pd.read_csv(Filename, sep=",", names=col_names, skiprows=[0])

    temp_df = df.astype(float)

    print("Columns in dataframe will be: " + str(temp_df.columns.tolist()))

    # Below we apply limits to our csv, done in two seperate chunks, the 
    # positive and negative values of sinba. These are then recombined.
    # values used as limits here come from arxiv:2011.03652
    upper_df = temp_df[temp_df['khuu'] >= 0.88]
    upper_df = upper_df[upper_df['khuu'] <= 1.2]

    lower_df = temp_df[temp_df['khuu']<-0.7]
    lower_df = lower_df[lower_df['khuu']>=-1.1]

    temp_df = combiner([upper_df, lower_df])

    temp_df = temp_df[temp_df['mH'] < 1000]
    temp_df = temp_df[temp_df['mA'] < 1000]
    temp_df.reindex()

    df = tan_trnsfm(temp_df)

    df.dropna(inplace=True)

    # When excluding the charged higgs we require masses to be above 580GeV
    # arxiv:1702.04571 
    if inc_Hc == False:
        temp_df = df.drop(df[df.mHc < 580].index)

    print("Filtering points...")
  
    ### - Filtering
    
    sig1 = 0.6827
    sig2 = 0.9545
    sig3 = 0.9973
    
    chi2_1sig_temp_df2 = stats.chi2.ppf( q=sig1, df=2)
    chi2_2sig_temp_df2 = stats.chi2.ppf( q=sig2, df=2)
    chi2_3sig_temp_df2 = stats.chi2.ppf( q=sig3, df=2)
    chi2_1sig_temp_df5 = stats.chi2.ppf( q=sig1, df=5)
    chi2_2sig_temp_df5 = stats.chi2.ppf( q=sig2, df=5)
    chi2_3sig_temp_df5 = stats.chi2.ppf( q=sig3, df=5)
    
    temp_df_bad = temp_df.query('chi2_HS > 10.0 & chi2_HS < 100.0')
    
    chi2_HS_min          = temp_df_bad.chi2_HS.min()
    print("chi2_HS_min", chi2_HS_min)

    chi2_Tot_gfitter_min = temp_df_bad.chi2_Tot_gfitter.min()
    chi2_Tot_hepfit_min  = temp_df_bad.chi2_Tot_hepfit.min()
    
    chi2_Tot_sig1_lim = chi2_Tot_gfitter_min + chi2_1sig_temp_df5
    chi2_Tot_sig2_lim = chi2_Tot_gfitter_min + chi2_2sig_temp_df5
    chi2_Tot_sig3_lim = chi2_Tot_gfitter_min + chi2_3sig_temp_df5
    
    chi2_Tot_gfitter_upper_bound = 120.0
    
    chi2_Tot_mask = (temp_df["chi2_Tot_gfitter"] < chi2_Tot_gfitter_upper_bound)
    print("chi2_Tot_mask")
    print chi2_Tot_mask.head(10)


    # Using EWPO to check if points are allowed and discarding any that are not
    per_4pi = 1.0
    sta     = 1.0
    uni     = 1.0
        
    # Setting up pandas 'mask' to apply to data    
    theory_mask = (temp_df["per_8pi"] == per_4pi) & (temp_df["sta"] == sta) & (temp_df["uni"] == uni)
        
    mask = theory_mask & chi2_Tot_mask
    print("temp_df", temp_df.head(10)) 
    
    allowed_pts_df = temp_df[mask]
    print("allowed_pts", allowed_pts_df.head(10))

    allowed_pts_df.reset_index(drop=True, inplace=True)    
    
    print("Number of data points in the reduced DataFrame: {}".format(len(allowed_pts_df)) )

    #ensures values are floats and resets index
    new_df = allowed_pts_df.astype(float)
    new_df = new_df.reset_index(drop = True)

    return new_df
###############################################################################


################################################################################
def combiner(df_list):
    """
    Takes a list of pandas dataframes, combines them and removes any duplicate
    rows.

    Parameters
    ----------
    df_list : list
        This is a list of pandas dataframes that has already been read in from 
        a csv or other file type.

    Returns
    -------
    final_df : pandas dataframe
        Output pandas dataframe which will be the result of concatenating and
        'tidying' input dataframes.
    """


    combi_df = pd.concat(df_list, ignore_index = True)
    combi_df.dropna()
    final_df = combi_df.drop_duplicates()
    combi_df.info()
    final_df = final_df.astype(float)

    return final_df
################################################################################


################################################################################
def alpha_beta(tanbeta, cosba):
    """
    Takes lists of tan(beta) and cos(beta-alpha) values and uses them to
    calculate the values of alpha and beta then return them as lists

    Parameters
    ----------
    tanbeta : list of floats
        the values for tan(beta) as decimals

    cosba : list of floats
        the values for cos(beta-alpha) as decimals

    Returns
    -------
    alpha_val : list of floats
        the corresponding alpha values for the input data
   
    beta_val : list of floats
        the corresponding beta values for the input data

 """

    if len(tanbeta) != len(cosba):
        raise Exception("Lists are not the same length, there are missing\
 values!")

    for i in range(0, len(tanbeta)):
        print("Calculating beta values")
        beta_val = np.arctan(tanbeta[i])

        # Condition below removes unphysical cos(beta-alpha) values
        if -1<=(cosba[i])<=1 :
            print("Calculating alpha values")
            alpha_val = -1*np.arcsin(cosba[i])+beta_val
        else:
            Exception("Error! Unphysical cos(beta-alpha) value provided!")

    return(alpha_val, beta_val)
###############################################################################


###############################################################################
def user_interface():
    """
    Runs a basic interface to prompt the user for input about which functions
    they wish to run and on what files they wish to run them on. Can convert
    .dat files to csv files and combine csv files. Will run filtering on data
    to ensure it complies with current constraints.
    """


    try:
        option = input(" Enter 1 to convert a dat file to csv or enter 2 to \
 combine multiple csves.")
    
        if option == 1:
            print("What is the path to, and name of the file (in the form of a\
 string), for conversion?")
            file_name = input("You do not need to add .dat/.csv, this will be\
 added automatically: ")
        
            nther_file = file_name + ".dat"
            dat_to_DF(nther_file, file_name + '.csv')

            tn_fxd = make_useful(file_name + '.csv', inc_Hc=False)
            print(file_name + ".dat has been converted to " + file_name + 
".csv")
                
            print("Please enter name for final csv ")
            nw_name = input("You do not need to add .csv, this will be added\
 automatically ")

            tn_fxd.to_csv(nw_name +'.csv', index=False) 

        elif option == 2:
            print("Please enter a list, separated by spaces only, of csvs to\
 combine: ")

            list_input = input("You do not need to add .csv to the end of\
 files, this will be added automatically")

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
            nw_name = input("You do not need to add .csv, this will be added\
 automatically")
            combined_csv.to_csv(str(nw_name) + ".csv", index=False)
    except ValueError as e:
        print("Error! Please enter 1 or 2 only.")

###############################################################################

names= ['Z7', 'mH', 'mHc', 'mA', 'cba',
         'tb', 'sinba', 'Z4', 'Z5', 'm12_2', 'l1', 'l2', 'l3', 'l4', 'l5',
         'l6', 'l7', 'g_HpHmh', 'Gamma_h', 'Gamma_H', 'Gamma_Hc', 'Gamma_A',
         'br_h_bb', 'br_h_tautau', 'br_h_gg', 'br_h_WW', 'br_h_ZZ',
         'br_h_gaga', 'br_A_tt', 'br_A_bb', 'br_A_gg', 'br_A_mumu',
         'br_A_tautau', 'br_A_Zga', 'br_A_Zh', 'br_A_ZH', 'br_A_gaga',
         'br_H_tt', 'br_H_bb', 'br_H_gg', 'br_H_mumu', 'br_H_tautau',
         'br_H_Zga', 'br_H_Zh', 'br_H_WW', 'br_H_ZZ', 'br_H_ZA', 'br_H_hh',
         'br_H_AA', 'br_H_gaga', 'br_Hp_tb', 'br_Hp_taunu', 'br_Hp_Wh',
         'br_Hp_WH', 'br_Hp_WA', 'sta', 'uni', 'per_4pi', 'per_8pi', 'S', 'T',
         'U', 'V', 'W', 'X', 'delta_rho', 'delta_amu', 'tot_hbobs', 'sens_ch',
         'chi2_HS', 'chi2_ST_hepfit', 'chi2_ST_gfitter', 'chi2_Tot_hepfit',
         'chi2_Tot_gfitter', 'k_huu', 'k_hdd', 'likelihood', 'stay_count',
         'ratio'],

user_interface()

