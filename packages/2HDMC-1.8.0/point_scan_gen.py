# Script to scan over parameters in 2HDMC, with HB and HS checks
import os
#import math
import numpy as np
import pandas as pd
import sys
#from multiprocessing import Pool
#from functools import partial 
#import random
#from datetime import datetime

csv_file = sys.argv[1]
new_csv = sys.argv[2]

# Default (lambda6 = lambda7 = 0)
def check_point(mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt, SLHA_filename):
	"""
	Runs 2HDMC on the given values, if 2HDMC returns Pass it then tests HB+HS on
	input points, writes to file and labels whether point passed checks.

	Input: mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt

	Output: "Pass" or "Fail", depending on result of 2HDMC+HB+HS
	"""
	# Compute alpha and beta for point
	beta = np.arctan(tanb)
	alpha = beta - np.arcsin(sinba)

	# Compute kappas
	# higgs h (h1)
	k1tt = np.cos(alpha)/np.sin(beta)
	k1bb = -1*np.sin(alpha)/np.cos(beta)
	k1ww = np.sin(beta-alpha)

	# higgs H (h2)
	k2tt = np.sin(alpha)/np.sin(beta)
	k2bb = np.cos(alpha)/np.cos(beta)
	k2ww = np.cos(beta-alpha)

	# higgs A (h3)
	k3tt = 1/np.tan(beta)
	k3bb = np.tan(beta)

	# Execute 2HDMC
	os.system("/scratch/cb27g11/THDM_T3PS_scanner/packages/2HDMC-1.8.0/check_point_2hdmc "\
		+str(mh)+" "+str(mH)+" "+str(mA)+" "+str(mHp)+" "+str(sinba)\
		+" "+str(lam6)+" "+str(lam7)+" "+str(m12)+" "+str(tanb)+" "\
		+str(yt)+" "+SLHA_filename+" ")

	######### LET 2HDMC RUN ##########

	# If point passes 2HDMC and HBHS checks, write to file
	if os.path.isfile("PlaceHolder") == True:
		
		with open(new_csv,'a') as f:
			f.write(str(1)+","+str(mh)+","+str(mH)+","+str(mA)\
				+","+str(mHp)+","+str(sinba)+","+str(m12)\
				+","+str(tanb)+","+str(alpha)+","+str(beta)\
				+","+str(k1tt)+","+str(k1bb)+","+str(k1ww)\
				+","+str(k2tt)+","+str(k2bb)+","+str(k2ww)\
				+","+str(k3tt)+","+str(k3bb))				
			f.write("\n")
			f.close()
		
		# remove placeholder otherwise failed points will be saved
				os.system("rm PlaceHolder")
		
			print("A pass!")
			return "Pass"

	if os.path.isfile("PlaceHolder") == False:

		with open("fail.csv",'a') as f:

			f.write(str(0)+","+str(mh)+","+str(mH)+","+str(mA)+","\
				+str(mHp)+","+str(sinba)+","+str(m12)+","\
				+str(tanb)+","+str(alpha)+","+str(beta)+","\
				+str(k1tt)+","+str(k1bb)+","+str(k1ww)+","\
				+str(k2tt)+","+str(k2bb)+","+str(k2ww)+","\
				+str(k3tt)+","+str(k3bb))
				
			f.write("\n")
			f.close()
		print("Fail")
		return "Fail"

def point_input(filename, yt, savefile):
	"""
	Reads variables (mH, mHp, mA, sinba, tanb, l6, l7) in from given csv
	file, calls 'check_point' on each point, creates a results csv file
	and writes passing points to it.

	Input: name of file and value of yt (yt is yukawa coupling type) 

	Output: True/False dependent on whether point passes checks
	"""
	# create results file
	if os.path.isfile(savefile) == False:

		os.system("touch " + savefile)

		with open(savefile,'a') as f:
			f.write("pass,mh,mH,mA,mHp,sinba,m12,tanb,alpha,beta,\
				kappa1tt,kappa1bb,kappa1ww,kappa2tt,kappa2bb,\
				kappa2ww,kappa3tt,kappa3bb")
			f.write("\n")
	
	# Reads input csv file into DataFrame
	dataf = pd.read_csv(filename)
	print(dataf.head(5))

	# Read off variable values
	mh_all, mH_all, mA_all, mHp_all = 125.09, dataf['mH'].tolist(),\
				dataf['mA'].tolist(), dataf['mHc'].tolist()

	sinba_all, lam6_all, lam7_all = dataf['sinba'].tolist(),\ 
				dataf['l6'].tolist(), dataf['l7'].tolist()

	m12_all, tanb_all = dataf['m12'].tolist(), dataf['tb'].tolist()

	# scan over points
	for i in range(len(mH_all)):
	# Assign variable values for each point

		mh, mH, mA, mHp = mh_all, mH_all[i], mA_all[i], mHp_all[i]

		sinba, lam6, lam7 =  sinba_all[i], lam6_all[i], lam7_all[i]

		m12_all[i], tanb_all[i] = m12, tanb

		file_name = "SLHA_savefile" + str(i)

		check_point(mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt,
			 file_name)
	return True


# Set parameters and execute scan 
point_input(csv_file, 2, new_csv)
