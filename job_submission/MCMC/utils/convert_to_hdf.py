#!/usr/bin/env python

import tables
import pandas as pd
import sys
import time
import os
from scipy import stats
from glob import glob

def ASCII_to_pd_h5f( input, output, out_ds_name, form, compr ):

    import pandas as pd

    start = time.time()
    all_files = glob( input )
    print('Input:\n', all_files)
    
    print('Reading in file(s) {}...'.format(input) )
#   df_from_each_file = (pd.read_csv(f, sep="\s+|\t+|\s+\t+|\t+\s+") for f in all_files)
#   df_from_each_file = (pd.read_csv(f, delim_whitespace=True) for f in all_files)
    
#   print( 'Concatenating...' )
#   df = pd.concat(df_from_each_file, ignore_index=True)
 
    df = pd.read_csv(input, delim_whitespace=True)

    ### - Filtering

    sig1 = 0.6827
    sig2 = 0.9545
    sig3 = 0.9973

    chi2_1sig_df2 = stats.chi2.ppf( q=sig1, df=2)
    chi2_2sig_df2 = stats.chi2.ppf( q=sig2, df=2)
    chi2_3sig_df2 = stats.chi2.ppf( q=sig3, df=2)
    chi2_1sig_df5 = stats.chi2.ppf( q=sig1, df=5)
    chi2_2sig_df5 = stats.chi2.ppf( q=sig2, df=5)
    chi2_3sig_df5 = stats.chi2.ppf( q=sig3, df=5)

    df_excl_bad = df.query('chi2_HS > 10.0 & chi2_HS < 100.0')

    chi2_HS_min          = df_excl_bad.chi2_HS.min()
    chi2_Tot_gfitter_min = df_excl_bad.chi2_Tot_gfitter.min()
    chi2_Tot_hepfit_min  = df_excl_bad.chi2_Tot_hepfit.min()

    chi2_Tot_sig1_lim = chi2_Tot_gfitter_min + chi2_1sig_df5
    chi2_Tot_sig2_lim = chi2_Tot_gfitter_min + chi2_2sig_df5
    chi2_Tot_sig3_lim = chi2_Tot_gfitter_min + chi2_3sig_df5

    chi2_Tot_gfitter_upper_bound = 120.0

    mask_chi2_Tot = (df["chi2_Tot_gfitter"] < chi2_Tot_gfitter_upper_bound)

    per_4pi = 1.0
    sta     = 1.0
    uni     = 1.0

    mask_theory = (df["per_4pi"] == per_4pi) & (df["sta"] == sta) & (df["uni"] == uni)

#    allowed_df = df_excl_bad.query('chi2_Tot_gfitter < {} & tot_hbobs < 1.0'.format(  chi2_Tot_sig3_lim ) )
#    allowed_df = df_excl_bad.query('chi2_Tot_gfitter < {}'.format(  chi2_Tot_sig3_lim ) )
#
    mask = mask_theory & mask_chi2_Tot

    allowed_df = df[mask]

    allowed_df.reset_index(drop=True, inplace=True)
#    reduced = allowed_df[:2000000]
    reduced = allowed_df[:]

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
