#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 13:20:37 2020

@author: cb27g11
"""
import numpy as np

def transfmkappa(tanb, sinba, kappas, xsection, higgs):
    """
    Designed for 2HDMs. Takes an array of data points that are in terms of
    tan(beta)/sin(beta-alpha) and transforms it into terms of a given kappa.
    i.e. kappatt, the coupling of the SM higgs to a kappa term.

    Parameters
   
    """

    rdcd_xsects = []
    cos_vals = []
    k_tt1 = []
    k_bb1 = []
    k_ww1 = []
    k_tt2 = []
    k_bb2 = []
    k_ww2 = []
    k_tt3 = []
    k_bb3 = []

    beta_no = 0
    alpha_no = 0

    for i in range(0, len(xsection)):
        tanb = np.array(tanb)
        sinba = np.array(sinba)
        xsection = np.array(xsection)
        
        tanb = tanb.astype(float)
        sinba = sinba.astype(float)
        xsection = xsection.astype(float)
        beta_val = np.arctan(tanb[i])

        if -1<=(sinba[i])<=1 :
            beta_no += 1
            alpha_no += 1
            alpha_val = -1*np.arcsin(sinba[i]) + beta_val
            rdcd_xsects.append(xsection[i])

            cos_val = np.cos(beta_val - alpha_val)
            cos_vals.append(cos_val)
            #Checking values are usable
            
            if "h1" in higgs:

                if "tt" in kappas:
                    # now we can convert our values into ones in terms of kappatt
                    k_tt1.append(np.around(np.cos(alpha_val)/np.sin(beta_val), 10))

                if "bb" in kappas:
                    # now we can convert our values into ones in terms of kappabb
                    k_bb1.append(np.around(-1*np.sin(alpha_val)/np.cos(beta_val), 10))

                if "WW" in kappas:
                    # now we can convert our values into ones in terms of kappaWW
                    k_ww1.append(np.around(np.sin(beta_val-alpha_val), 10))


            if "h2" in higgs:

                if "tt" in kappas:
                    # now we can convert our values into ones in terms of kappatt
                    k_tt2.append(np.around(np.sin(alpha_val)/np.sin(beta_val), 10))

                if "bb" in kappas:
                    # now we can convert our values into ones in terms of kappabb
                    k_bb2.append(np.around(np.cos(alpha_val)/np.cos(beta_val), 10))

                if "WW" in kappas:
                    # now we can convert our values into ones in terms of kappaWW
                    k_ww2.append(np.around(np.cos(beta_val-alpha_val), 10))


            if "h3" in higgs:

                if "tt" in kappas:
                    # now we can convert our values into ones in terms of kappatt
                    k_tt3.append(np.around(1/np.tan(beta_val), 10))

            
                if "bb" in kappas:
                    # now we can convert our values into ones in terms of kappabb
                    k_bb3.append(np.around(np.tan(beta_val), 10))

        else:
            continue

    k_tt = [k_tt1, k_tt2, k_tt3]
    k_bb = [k_bb1, k_bb2, k_bb3]
    k_ww = [k_ww1, k_ww2]
    if len(k_tt) <= 1:
            k_tt = k_tt[0]
    if len(k_bb) <= 1:
            k_bb = k_bb[0]
    if len(k_ww) == 1:
            k_ww = k_ww[0]

    return(rdcd_xsects, k_tt, k_bb, k_ww, cos_vals)
