#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 19:25:07 2023

@author: cb27g11
"""
import numpy as np
import pandas as pd
import os
import sys

from scipy import stats

#I inherit the following variables from merge-jobs.sh:
# X_SECT_COL_

# This file will inherit values via 'sys.argv' from the corresponding
# looper.sh file

#location of combined madgraph output
madgraph_out = [str(sys.argv[1])]

#location of the original csv file prior to splitting
OG_csv = str(sys.argv[2])

#location to save the output to
FINAL_CSV = str(sys.argv[3])
print('FINAL_CSV', FINAL_CSV)




#Self explanatory (will make TOP_DIR_ to be replaced with user dir on install)
Top_dir = "/scratch/cb27g11/THDM_T3PS_scanner/"

#directory we start from
job_dir = Top_dir + "job_submission/Madgraph_jobs/" 



for data_file in madgraph_out:

    current_i = madgraph_out.index(data_file)
    print(data_file, "this was file number " + str(current_i) )
    data_base = pd.read_csv(OG_csv)
    print(data_base.columns.tolist())
    data_in = pd.read_csv(job_dir + data_file)
    
    data_base = data_base.drop_duplicates()
    data_in = data_in.drop_duplicates()

    var_matching(data_base, csv2, 

###############################################################################
def var_matching(csv1, csv2, var_list, search_list, rd_to):
"""

"""
 
    main_in = data_base.X_sect_bg_A.tolist()
    main_tb = data_base.tb.tolist()
    main_tb = np.array(main_tb)
    main_tb = np.around(main_tb, decimals=4)
    main_sinba = data_base.sinba.tolist()
    main_sinba = np.array(main_sinba)
    main_sinba = np.around(main_sinba, decimals=5)
    main_mH = data_base.mH.tolist()
    main_mH = np.array(main_mH)
    main_mH = np.around(main_mH, decimals=3)
    main_mHc = data_base.mHc.tolist()
    main_mHc = np.array(main_mHc)
    main_mHc = np.around(main_mHc, decimals=3)
    #print(len(main_tb), len(main_sinba), len(main_mH))
    
    
    new_in = data_in.bg_A_xsect.tolist()
    new_tb = data_in.Tb_list.tolist()
    new_tb = np.array(new_tb)
    new_tb = np.around(new_tb, decimals=4)
    new_sinba = data_in.Sbma_list.tolist()
    new_sinba = np.array(new_sinba)
    new_sinba = np.around(new_sinba, decimals=5)
    new_mH = data_in.mH.tolist()
    new_mH = np.array(new_mH)
    new_mH = np.around(new_mH, decimals=3)
    new_mHc = data_in.mHp.tolist()
    new_mHc = np.array(new_mHc)
    new_mHc = np.around(new_mHc, decimals=3)
    #print(len(new_tb), len(new_sinba), len(new_mH))
    
    
    def check_4_vals(main_list, new_list):
    
        for i in range(0,len(new_list)):
            if new_tb[i] in main_tb:
                tan_index = np.where(main_tb == new_tb[i])
                tan_index = tan_index[0][0]
                if new_mH[i] in main_mH:
                    mH_index = np.where(main_mH == new_mH[i])
                    mH_index = mH_index[0][0]
                    if mH_index != tan_index:
                        print("mass was for a different point")
                        continue
                    elif new_sinba[i] == main_sinba[mH_index]:
                        print(mH_index, main_list[mH_index], new_list[i])
                        if new_mHc[i] == main_mHc[mH_index]:
                            main_list[mH_index] = new_list[i]
                            print("success!")
                        else:
                            print("mHc didn't match")
                            continue
                    else:
                        print("mass was correct but sin value wasn't")
                        continue
    
        return main_list
    
    main_in = check_4_vals(main_in, new_in)
    data_base['X_sect_bg_A'] = main_in
    
    data_base.to_csv("Type2masscorrect_test.csv", index=False)

###############################################################################


