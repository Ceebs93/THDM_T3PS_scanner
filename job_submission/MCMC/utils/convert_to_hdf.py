#!/usr/bin/env python

import os
import pandas as pd
import tables
import time
import sys

from data_to_csv import tan_trnsfm
from glob import glob
from tqdm import tqdm
from scipy import stats

# A list of different regularly used headers can be added here.
          
headers = {
            "default" : ['Z7_mcmc', 'mH_mcmc', 'mHc_mcmc', 'mA_mcmc', 'cba_mcmc',
                'tb_mcmc', 'Z7_2HDMC', 'mH_2HDMC', 'mHc_2HDMC', 'mA_2HDMC',
                'cba_2HDMC', 'tb_2HDMC', 'sinba', 'Z4_c', 'Z5_c', 'm12_2','l1',
                'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'g_HpHmh', 'Gamma_h',
                'Gamma_H', 'Gamma_Hc', 'Gamma_A', 'br_h_bb', 'br_h_tautau',
                'br_h_gg', 'br_h_WW', 'br_h_ZZ', 'br_h_gaga', 'br_A_tt',
                'br_A_bb', 'br_A_gg', 'br_A_mumu', 'br_A_tautau', 'br_A_Zga',
                'br_A_Zh', 'br_A_ZH', 'br_A_gaga', 'br_H_tt', 'br_H_bb',
                'br_H_gg', 'br_H_mumu', 'br_H_tautau', 'br_H_Zga', 'br_H_Zh',
                'br_H_WW', 'br_H_ZZ', 'br_H_ZA', 'br_H_hh', 'br_H_AA',
                'br_H_gaga', 'br_Hp_tb', 'br_Hp_taunu', 'br_Hp_Wh', 'br_Hp_WH',
                'br_Hp_WA', 'sta', 'uni', 'per_4pi', 'per_8pi', 'S', 'T', 'U',
                'V', 'W', 'X', 'delta_rho', 'delta_amu', 'tot_hbobs', 
                'sens_ch', 'chi2_HS', 'chi2_ST_hepfit', 'chi2_ST_gfitter',
                'chi2_Tot_hepfit', 'chi2_Tot_gfitter', 'k_huu', 'k_hdd',
                'likelihood', 'stay_count']
          }


constraints = { "2HDM_TII" : { 
# values here come from arxiv:2011.03652
                            'khuu' : [-1.1, -0.7, 0.88, 1.2],
# values here come from      
                            'mH'   : [150, 1000],
# value here comes from arxiv1702.04571
                            'mHc'  : [580, 'L'],
# Values here come from Chinese Physics C Vol. 44, No. 7 (2020) 073101 (upper)
# arXiv:1305.1649v2 (lower)
                            'tb'   : [5,30]
                             }
              }



def ASCII_to_pd_h5f( input, output, out_ds_name, form, compr ):
    """
    Takes a datafile and converts it to a h5df file according to given
    parameters.    

    Parameters
    ---------
    input : string
        path to and name of the datafile to convert to h5df

    output : string
        path to and name for the output h5df file

    out_ds_name : string
        name of the job the data has been taken from

    form : string
        indicates the format for the h5df file

    compr : string
        indicates the compression to be used in creating the h5df file

    Returns
    -------
    None : converts datafile to h5df form and writes to desired path

    """


    start = time.time()
    all_files = glob( input )
    print('Input:\n', all_files)
    
    print('Reading in file(s) {}...'.format(input) )
    df_list = []

    chunksize = 100000

    # Splitting datapoints from input into chunks to increase processing speed 
    for df_chunk in tqdm(pd.read_csv(input, index_col=None,\
              delim_whitespace=True, error_bad_lines=False,\
              chunksize=chunksize)):

        ### - Tan correction
        
        



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
	print("chi2_HS_min")
	print chi2_HS_min
        chi2_Tot_gfitter_min = df_chunk_excl_bad.chi2_Tot_gfitter.min()
        chi2_Tot_hepfit_min  = df_chunk_excl_bad.chi2_Tot_hepfit.min()
    
        # Calculating limits on chi2_Tot values for different sigma values
        chi2_Tot_sig1_lim = chi2_Tot_gfitter_min + chi2_1sig_df_chunk5
        chi2_Tot_sig2_lim = chi2_Tot_gfitter_min + chi2_2sig_df_chunk5
        chi2_Tot_sig3_lim = chi2_Tot_gfitter_min + chi2_3sig_df_chunk5
    
        chi2_Tot_gfitter_upper_bound = 120.0
    
        # Creating 'mask' to filter out datapoints with chi2_Tot_gfitter values
        # higher than the upper bound value
        mask_chi2_Tot = (df_chunk["chi2_Tot_gfitter"] < chi2_Tot_gfitter_upper_bound)
    	print("mask_chi2_Tot")
	print mask_chi2_Tot.head(10)
        per_4pi = 1.0
        sta     = 1.0
        uni     = 1.0
    
        # Creating 'mask' to filter out datapoints that do not pass theory
        # constraints (i.e. perturbativity, stability and unitarity)
        mask_theory = (df_chunk["per_8pi"] == per_4pi) & \
                      (df_chunk["sta"] == sta) & (df_chunk["uni"] == uni)
   
        mask = mask_theory & mask_chi2_Tot 
   	print(df_chunk.head(10)) 
        # Applying masks to datapoints to filter them
        allowed_df_chunk = df_chunk[mask]
    	print(allowed_df_chunk.head(10))
        allowed_df_chunk.reset_index(drop=True, inplace=True)

        df_list.append(allowed_df_chunk)

    # Combining all filtered chunks of data into one dataframe
    reduced = pd.concat(df_list)

    # Correcting beta and tan(beta) values so that beta<=2pi. This prevents
    # erroneously 

    print("Number of data points in the reduced DataFrame: {}".format(len(reduced)) )

    
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
