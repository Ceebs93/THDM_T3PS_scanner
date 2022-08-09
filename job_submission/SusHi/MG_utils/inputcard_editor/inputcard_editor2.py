from __future__ import absolute_import, division, print_function
import numpy as np
import fileinput
import sys

#sys.argv[0] - this is the name of this script
#sys.argv[1] - this is the first passed argument from bash
Hm = float(sys.argv[1]) #CP even higgs mass
Hcm = float(sys.argv[2]) #Charged higgs mass
Am = float(sys.argv[3]) #CP odd higs mass
textfilepath = str(sys.argv[4]) #path to txt file madgraph will use
print("inputcard_editor has textfilepath = " + str(textfilepath))
Tb = sys.argv[5] #Current tan value
print("inputcard_editor has Tb = " + str(Tb))
Sbma = sys.argv[6] #Current sin value
print("inputcard_editor has Sbma = " + str(Sbma))

Process = str(sys.argv[7]) #Name of folder
print("inputcard_editor has Process = " + str(Process))

def the_main_event():
	print('I AM ACTUALLY RUNNING')
	with open('/scratch/cb27g11/mg_run_basic/nw_mg_runT2.txt','r') as old_card:
		text = old_card.read() #stores string of old_card ready for editing
		make_input(Tb, Sbma, Hcm, text, Process)
                print('about to run madgraph ' + str(Tb) + ' ' + str(Sbma))

def make_input(Tanb, Sinbma, Hcpm, datatext, procname):
# inputs are the value of tan_beta, the value of sin(beta-alpha) values, the desired mass for the charged higgses and a string of text

        with open(textfilepath, 'w') as new_card:
        #simulation card, the .txt file that gets fed to madgraph
		print(textfilepath)
                sim_card = datatext.replace('TBV', str(Tanb))
		print('did tan')
                sim_card = sim_card.replace('SBMAV', str(Sinbma))
		sim_card = sim_card.replace('HMV', str(Hm))
		sim_card = sim_card.replace('AMV', str(Am))
                sim_card = sim_card.replace('HCMV', str(Hcm))
                new_name = '/scratch/cb27g11/' + str(procname) + '_' + str(Tanb) + '_' + str(Sinbma)
		print(new_name)
                sim_card = sim_card.replace('NAME', str(new_name))
		print('did the name!')
                print('/scratch/cb27g11/' + str(procname) + '_' + str(Tanb) + '_' + str(Sinbma))
                new_card.write(sim_card) # saves new txt file for madgraph

the_main_event()
