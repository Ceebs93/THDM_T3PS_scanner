#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 19:34:32 2021

@author: cb27g11
"""

from __future__ import absolute_import, division, print_function
import fileinput
import numpy as np
import pandas as pd
import sys

from scipy import stats

# Filename
file_name = str(sys.argv[1])
# dat_to_DF will assume file it needs is jobs/file_name/file_name_all_final.dat
file_path = "jobs/" + file_name + "/" + file_name + "_all_data_merged.dat"

# Yukawa type for applying bounds
yukawa = str(sys.argv[2])

# Basis for column names
base = str(sys.argv[3])

bases = {
         "hybrid" :
         ['Z7_out', 'mH_out', 'mHc_out', 'mA_out', 'cba_out', 'tb_out', 'Z7',
         'mH', 'mHc', 'mA', 'cba', 'tb', 'sinba', 'Z4', 'Z5', 'm12_2', 'l1',
         'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'g_HpHmh', 'Gamma_h', 'Gamma_H',
         'Gamma_Hc', 'Gamma_A', 'br_h_bb', 'br_h_tautau', 'br_h_gg',
         'br_h_WW', 'br_h_ZZ', 'br_h_gaga', 'br_A_tt', 'br_A_bb', 'br_A_gg',
         'br_A_mumu', 'br_A_tautau', 'br_A_Zga', 'br_A_Zh', 'br_A_ZH',
         'br_A_gaga', 'br_H_tt', 'br_H_bb', 'br_H_gg', 'br_H_mumu',
         'br_H_tautau', 'br_H_Zga', 'br_H_Zh', 'br_H_WW', 'br_H_ZZ', 'br_H_ZA',
         'br_H_hh', 'br_H_AA', 'br_H_gaga', 'br_Hp_tb', 'br_Hp_taunu',
         'br_Hp_Wh', 'br_Hp_WH', 'br_Hp_WA', 'sta', 'uni', 'per_4pi',
         'per_8pi', 'S', 'T', 'U', 'V', 'W', 'X', 'delta_rho', 'delta_amu',
         'tot_hbobs', 'sens_ch', 'chi2_HS', 'chi2_ST_hepfit',
         'chi2_ST_gfitter', 'chi2_Tot_hepfit', 'chi2_Tot_gfitter', 'k_huu',
         'k_hdd', 'likelihood', 'stay_count'],
         "hhg" :
         [],
         "generic" :
         [],
         "higgs" :
         [],
         "phys" :
         []
     }

# Dictionary of limits based on yukawa coupling type. Values are given as:
#[khuu +ve lowerbound, khuu +ve upperbound, khuu -ve upperbound,
#khuu -ve lowerbound, min_mHp]
yukawas = { "1" : [],
            "2" : [0.88, 1.2, -0.7, -1.1, 580],
            "3" : [],
            "4" : []
     }
# Type-I khuu values used as limits here come from arxiv:2011.03652
# Type-II khuu values used as limits here come from arxiv:
# Type-III khuu values used as limits here come from arxiv:
# Type-IV khuu values used as limits here come from arxiv:
# Type-II lower charged higgs mass bound arxiv:1702.04571 

################################################################################
def dat_to_DF(Filename, Filepath, Yukawa, Base):
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

    # Select relevant dictionary entry for column names and coupling bounds etc
    col_names = bases[Base]
    bounds = yukawas[Yukawa]

    df = pd.read_csv(Filepath, sep="\s+", skiprows=range(0,201))
    print(len(df))
    df.set_axis(col_names, axis=1, inplace=True)

    print("Dataframe columns will be: " + str(df.columns.tolist()))

    df = df.dropna()
    df = df.drop_duplicates()
    print(df.head(5))

    temp_df = df.astype(float)

    # Below we apply limits to our csv, done in two seperate chunks, the 
    # positive and negative values of sinba. These are then recombined.
    # values used as limits here come from arxiv:2011.03652
    upper_df = temp_df[temp_df['k_huu'] >= bounds[0]]
    upper_df = upper_df[upper_df['k_huu'] <= bounds[1]]

    lower_df = temp_df[temp_df['k_huu']< bounds[2]]
    lower_df = lower_df[lower_df['k_huu']>= bounds[3]]

    temp_df = combiner([upper_df, lower_df])
    temp_df.reindex()

    temp_df.dropna(inplace=True)

    # When excluding the charged higgs we require masses to be above 580GeV
    # arxiv:1702.04571 
    temp_df = temp_df.drop(df[df.mHc < bounds[4]].index)

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
    print(chi2_Tot_mask.head(10))


    # Using EWPO to check if points are allowed and discarding any that are not
    per_8pi = 1.0
    sta     = 1.0
    uni     = 1.0
        
    # Setting up pandas 'mask' to apply to data    
    theory_mask = (temp_df["per_8pi"] == per_8pi) & (temp_df["sta"] == sta) & (temp_df["uni"] == uni)
        
    mask = theory_mask & chi2_Tot_mask
    print("temp_df", temp_df.head(10)) 
    
    allowed_pts_df = temp_df[mask]
    print("allowed_pts", allowed_pts_df.head(10))

    allowed_pts_df.reset_index(drop=True, inplace=True)    
    
    print("Number of data points in the reduced DataFrame: {}".format(len(allowed_pts_df)) )

    #ensures values are floats and resets index
    new_df = allowed_pts_df.astype(float)
    new_df = new_df.reset_index(drop = True)

    save_path = "../jobs/" + file_name + "/" + file_name + "_all_final.csv"
    
    return new_df
###############################################################################


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

dat_to_DF(file_name, file_path, yukawa, base)
