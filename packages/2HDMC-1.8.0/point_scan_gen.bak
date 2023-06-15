# Script to scan over parameters in 2HDMC, with HB and HS checks
import os
import math
import numpy as np
import pandas as pd
import sys
from multiprocessing import Pool
from functools import partial 
import random
from datetime import datetime

# Define function to check set of points (lambda6 = lambda7 = 0)
def check_point(mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt, SLHA_filename):
	"""
	Input: 2HDMC input parameters; mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt
	Output: "Pass" of "Fail", depending on result of 2HDMC+HB+HS
	Runs 2HDMC, then if 2HDMC passes, tests HB+HS on input points, writes to file and labels whether point passed checks.
	"""
	# Compute alpha and beta for point
	beta = math.atan(tanb)
	alpha = beta - math.asin(sinba)
	# compute kappas
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
	#os.chdir('/scratch/cb27g11/2HDMC-1.8.0') 
	os.system("/scratch/cb27g11/THDM_T3PS_scanner/packages/2HDMC-1.8.0/check_point_2hdmc "+str(mh)+" "+str(mH)+" "+str(mA)+" "+str(mHp)+" "+str(sinba)+" "+str(lam6)+" "+str(lam7)+" "+str(m12)+" "+str(tanb)+" "+str(yt)+" "+SLHA_filename+" ")

	############################################## LET 2HDMC RUN ##############################################

	# If 2HDMC and HBHS check passes write to file
	if os.path.isfile("bmark") == True:
			with open(new_csv,'a') as f:
				f.write(str(1)+","+str(mh)+","+str(mH)+","+str(mA)+","+str(mHp)+","+str(sinba)+","+str(m12)\
					+","+str(tanb)+","+str(alpha)+","+str(beta)\
					+","+str(k1tt)+","+str(k1bb)+","+str(k1ww)+","+str(k2tt)\
					+","+str(k2bb)+","+str(k2ww)+","+str(k3tt)+","+str(k3bb))				
				f.write("\n")
				f.close()
				# remove benchmark otherwise all other point will be saved
				os.system("rm bmark")
			print("Hooray! A pass!")
			return "Pass"

	if os.path.isfile("bmark") == False:
		with open("fail.csv",'a') as f:
			f.write(str(0)+","+str(mh)+","+str(mH)+","+str(mA)+","+str(mHp)+","+str(sinba)+","+str(m12)\
				+","+str(tanb)+","+str(alpha)+","+str(beta)\
				+","+str(k1tt)+","+str(k1bb)+","+str(k1ww)+","+str(k2tt)\
				+","+str(k2bb)+","+str(k2ww)+","+str(k3tt)+","+str(k3bb))				
			f.write("\n")
			f.close()
		print("Nope...")
		return "Fail"

def point_input(filename, yt, savefile):
	"""
		Reads required variables (mH, mHp, mA, sinba, tanb, l6, l7) in from given csv file.
	Input: name of file and value of yt (yt is determining the type of yukawa couplings and corresponds to the type of 2HDM in use) 
	Output: file_out, containing all points checked and whether they pass 2HDMC/HBHS
	"""
	# create results file
	if os.path.isfile(savefile) == False:
		os.system("touch " + savefile)
		with open(savefile,'a') as f:
			f.write("pass,mh,mH,mA,mHp,sinba,m12,tanb,alpha,beta,kappa1tt,kappa1bb,kappa1ww,kappa2tt,kappa2bb,kappa2ww,kappa3tt,kappa3bb")
			f.write("\n")
	
	dataf = pd.read_csv(filename)
	print(dataf.head(5))

	# scan over points
	mh_all, mH_all, mA_all, mHp_all, sinba_all, lam6_all, lam7_all, m12_all, tanb_all =125, dataf['mH'].tolist(), dataf['mA'].tolist(), dataf['mHc'].tolist(), dataf['sinba'].tolist(),dataf['l6'].tolist(), dataf['l7'].tolist(), dataf['m12'].tolist(), dataf['tb'].tolist()

	# scan over points
	for i in range(len(mH_all)):
		mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb = mh_all, mH_all[i], mA_all[i], mHp_all[i], sinba_all[i], lam6_all[i], lam7_all[i], m12_all[i], tanb_all[i]
		file_name = "SLHA_savefile" + str(i)

		check_point(mh, mH, mA, mHp, sinba, lam6, lam7, m12, tanb, yt, file_name)
	return True

########################################################### SET SCAN PARAMETERS #################################################################

csv_file = sys.argv[1]
new_csv = sys.argv[2]
# Set params and execute scan!
# 
point_input(csv_file, 2, new_csv)
