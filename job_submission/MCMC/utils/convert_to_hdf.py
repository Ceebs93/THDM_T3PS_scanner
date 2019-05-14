#!/usr/bin/env python

import tables
import pandas as pd
import dask.dataframe as dd
import sys
import time
import os
from tqdm import tqdm
from scipy import stats
from glob import glob



headers = {
            "default" : ['Z7_mcmc', 'mH_mcmc', 'mHc_mcmc', 'mA_mcmc', 'cba_mcmc', 'tb_mcmc',
                         'Z7_2HDMC', 'mH_2HDMC', 'mHc_2HDMC', 'mA_2HDMC', 'cba_2HDMC', 'tb_2HDMC',
                         'sinba', 'Z4_c', 'Z5_c', 'm12_2',
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
                         'likelihood', 'stay_count']

          }


def ASCII_to_pd_h5f( input, output, out_ds_name, form, compr ):

    start = time.time()
    all_files = glob( input )
    print('Input:\n', all_files)
    
    print('Reading in file(s) {}...'.format(input) )
#   df_from_each_file = (pd.read_csv(f, sep="\s+|\t+|\s+\t+|\t+\s+") for f in all_files)
#   df_from_each_file = (pd.read_csv(f, delim_whitespace=True) for f in all_files)
    
#   print( 'Concatenating...' )
#   df = pd.concat(df_from_each_file, ignore_index=True)
 
    #df = pd.read_table(input, index_col=None,  delim_whitespace=True,  error_bad_lines=False)
    df_list = []

    chunksize = 100000

    for df_chunk in tqdm(pd.read_csv(input, index_col=None,  delim_whitespace=True,  error_bad_lines=False, chunksize=chunksize)):

#        df_list.append(df_chunk)

        ### - Filtering
    
        sig1 = 0.6827
        sig2 = 0.9545
        sig3 = 0.9973
    
        chi2_1sig_df_chunk2 = stats.chi2.ppf( q=sig1, df=2)
        chi2_2sig_df_chunk2 = stats.chi2.ppf( q=sig2, df=2)
        chi2_3sig_df_chunk2 = stats.chi2.ppf( q=sig3, df=2)
        chi2_1sig_df_chunk5 = stats.chi2.ppf( q=sig1, df=5)
        chi2_2sig_df_chunk5 = stats.chi2.ppf( q=sig2, df=5)
        chi2_3sig_df_chunk5 = stats.chi2.ppf( q=sig3, df=5)
    
        df_chunk_excl_bad = df_chunk.query('chi2_HS > 10.0 & chi2_HS < 100.0')
    
        chi2_HS_min          = df_chunk_excl_bad.chi2_HS.min()
        chi2_Tot_gfitter_min = df_chunk_excl_bad.chi2_Tot_gfitter.min()
        chi2_Tot_hepfit_min  = df_chunk_excl_bad.chi2_Tot_hepfit.min()
    
        chi2_Tot_sig1_lim = chi2_Tot_gfitter_min + chi2_1sig_df_chunk5
        chi2_Tot_sig2_lim = chi2_Tot_gfitter_min + chi2_2sig_df_chunk5
        chi2_Tot_sig3_lim = chi2_Tot_gfitter_min + chi2_3sig_df_chunk5
    
        chi2_Tot_gfitter_upper_bound = 120.0
    
        mask_chi2_Tot = (df_chunk["chi2_Tot_gfitter"] < chi2_Tot_gfitter_upper_bound)
    
        per_4pi = 1.0
        sta     = 1.0
        uni     = 1.0
    
        mask_theory = (df_chunk["per_8pi"] == per_4pi) & (df_chunk["sta"] == sta) & (df_chunk["uni"] == uni)
    
    #    allowed_df_chunk = df_chunk_excl_bad.query('chi2_Tot_gfitter < {} & tot_hbobs < 1.0'.format(  chi2_Tot_sig3_lim ) )
    #    allowed_df_chunk = df_chunk_excl_bad.query('chi2_Tot_gfitter < {}'.format(  chi2_Tot_sig3_lim ) )
    #
        mask = mask_theory & mask_chi2_Tot
    
        allowed_df_chunk = df_chunk[mask]
    
        allowed_df_chunk.reset_index(drop=True, inplace=True)
    #    reduced = allowed_df_chunk[:2000000]

        df_list.append(allowed_df_chunk)

    reduced = pd.concat(df_list)

    print("Number of data points in the reduced DataFrame: {}".format(len(reduced)) )

    #########
    
    print('Creating file {}...'.format(output) )
    reduced.to_hdf( output, key=out_ds_name, format=form, complib=compr )
    
    print( 'Finished.' )
    
    #####################
    ### --- Stats --- ###
    #####################
    
    input_size  = os.path.getsize( input  )/1024/1024
    output_size = os.path.getsize( output )/1024/1024
    
    end = time.time()
    elapsed = end-start
    
    print('Conversion time:')
    print('{:.0f}s.'.format(elapsed))
    print('Change in size:')
    print('input -> output')
    print('{:.0f} MB -> {:.0f} MB'.format( input_size, output_size))

if __name__ == "__main__":

    input       = str(sys.argv[1])
    output      = str(sys.argv[2])
    out_ds_name = str(sys.argv[3])
    form        = str(sys.argv[4])
    compr       = str(sys.argv[5])
    
    ASCII_to_pd_h5f(input, output, out_ds_name, form, compr)
