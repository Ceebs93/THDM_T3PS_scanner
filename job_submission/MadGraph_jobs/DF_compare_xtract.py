#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:50:19 2022

@author: cb27g11
"""

import pandas as pd
import numpy as np

#csvname = "Type2Data_mark2.csv"
OG_dat = "jobs/1209_midtanbq_h"
Chkd_dat = "Type2Data_xsecs_13_6TeV"

proc = 'T2_2607allnarrow_allcols.csv'

#df = pd.read_csv(Filename, sep=",", skiprows=[0])

#new_df = df.drop_duplicates()

#new_df.to_csv("tot_xsecs/bg_twhi_type1_OG_data.csv", index=False)

#column_list = ['Z7','mA','cba','Z4','Z5','l1','l2','l3','l4','l5','l6','l7','g_HpHmh',
#                'Gamma_h','Gamma_H','Gamma_Hc','Gamma_A','br_h_bb','br_h_tautau'
#                ,'br_h_gg','br_h_WW','br_h_ZZ','br_h_gaga','br_A_tt','br_A_bb',
#                'br_A_gg','br_A_mumu','br_A_tautau','br_A_Zga','br_A_Zh',
#                'br_A_ZH','br_A_gaga','br_H_tt','br_H_bb','br_H_gg','br_H_mumu',
#                'br_H_tautau','br_H_Zga','br_H_Zh','br_H_WW','br_H_ZZ','br_H_ZA'
#                ,'br_H_AA','br_H_hh','br_H_gaga','br_Hp_tb','br_Hp_taunu',
#                'br_Hp_Wh','br_Hp_WH','br_Hp_WA','sta','uni','per_4pi',
#                'per_8pi','S','T','U','V','W','X','delta_rho','delta_amu',
#                'tot_hbobs','sens_ch','chi2_HS','chi2_ST_hepfit',
#                'chi2_ST_gfitter','chi2_Tot_hepfit','chi2_Tot_gfitter','k_huu',
#                'k_hdd','likelihood']
#
#odf = pd.read_csv("2607_allnarrow.dat", sep=",", names =column_list)
#odf.to_csv("2607_allnarrow.csv", index=False)

column_list = [u'bq_tqhi_X_sects', u'bq_tqh1_X_sects', u'bq_tqh2_X_sects',
               u'bq_tqh3_X_sects', u'bg_twhi_X_sects', u'bg_twh1_X_sects',
               u'bg_twh2_X_sects', u'bg_twh3_X_sects', u'qq_tbhi_X_sects',
               u'qq_tbh1_X_sects', u'qq_tbh2_X_sects', u'qq_tbh3_X_sects']

def add_xsect(OG_data, Chkdfile, cols, proc):
    """Takes the name of an output CSV from the MG5 looper and the name of the
    file that was originally given to this and combines them to create a CSV
    with the original data as well as the MG5 calculated cross section"""

    #Adding .csv to end of filenames
    OG_data = OG_data + ".csv"
    Chkdfile = Chkdfile + ".csv"

    # Reading in the two csv files
    OG_df = pd.read_csv(OG_data)
    print(len(cols))
    print('len(OG_df)',len(OG_df['mA']))
    Chkd_df = pd.read_csv(Chkdfile)
    print('len(Chkd_df)',len(Chkd_df['mA']))
    # Checking that the dataframes have the same number of rows before adding
    # cross section column onto OG_df

    total_df = pd.DataFrame(columns=cols)

    for m in range(0,len(Chkd_df)):
        print('m is ', m)
        data_val = Chkd_df.loc[m,'mA']
        print('current data_val is ', data_val)
        if data_val in OG_df['mA'].values:
            #print('found ' + str(m) + ' in OG_df')
            idx = (OG_df[OG_df['mA']==data_val].index.values)
            fullrow_df = OG_df.loc[idx]
            new_row = fullrow_df[cols]
            print("new_row",new_row)
            print(type(new_row), type(total_df))
            total_df = pd.concat([total_df, new_row], ignore_index=True)
            print(total_df.head(5))

        #else:
            #print("couldn't find ", m)

    Both_DFs = pd.merge(OG_df.set_index('mA', drop=True),total_df.set_index('mA', drop=True), left_index=True, right_index=True).dropna().reset_index()
    Both_DFs.to_csv(proc, index=False)

add_xsect(OG_dat, Chkd_dat, column_list, proc)
