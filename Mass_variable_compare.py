#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 19:25:07 2023

@author: cb27g11
"""
import numpy as np
import pandas as pd
import os


user_home = os.path.expanduser('~') #Self explanatory
output_dir = user_home #+ "/Development/Ciara/Unharvested/" #directory we start from
#paths = os.listdir(output_dir)# All the different runs we need to include
paths =["/combo_wombo/PR_1209_midtanbg_h3.csv"]
#print(paths)

for data_file in paths:

    current_i = paths.index(data_file)
    print(data_file, "this was file number " + str(current_i) )
    data_base = pd.read_csv("/home/cb27g11/Development/Ciara/Type2masscorrect_test.csv")
    #data_base = pd.read_csv("/home/cb27g11/Development/Ciara/to_combine/new_constraints_datacombo.csv")
    print(data_base.columns.tolist())
    #data_in = pd.read_csv("Type2_fulldata.csv")
    #data_in = pd.read_csv("unclppd_h2aldata.csv")
    data_in = pd.read_csv(output_dir + data_file)
    
    data_base = data_base.drop_duplicates()
    data_in = data_in.drop_duplicates()
    
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
                #print(new_tb[i])
                tan_index = np.where(main_tb == new_tb[i])
                tan_index = tan_index[0][0]
#                print("tan values", main_tb[tan_index], new_tb[i])
                if new_mH[i] in main_mH:
                    #print(new_mH[i],"I found this mass!")
                    mH_index = np.where(main_mH == new_mH[i])
                    mH_index = mH_index[0][0]
#                    print("mH is ",main_mH[mH_index], new_mH[i])
                    if mH_index != tan_index:
                        print("mass was for a different point")
                        continue
                    elif new_sinba[i] == main_sinba[mH_index]:
                        print(mH_index, main_list[mH_index], new_list[i])
#                        print("mHc values are, ", main_mHc[mH_index], new_mHc[i])
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