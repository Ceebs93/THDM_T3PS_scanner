from __future__ import absolute_import, division, print_function
import numpy as np
import fileinput
import sys

# This file will inherit values via 'sys.argv' from the corresponding
# looper.sh file
# For reference:
# sys.argv[0] - this is the name of this script
# sys.argv[1] - this is the first passed argument from bash

Hm = float(sys.argv[1]) #CP even higgs mass
Hcm = float(sys.argv[2]) #Charged higgs mass
Am = float(sys.argv[3]) #CP odd higs mass
basecardpath = str(sys.argv[4]) #path to 'basic' runcard to edit and create final
runcardpath = str(sys.argv[5]) #path to 'final' runcard madgraph will use
Tb = sys.argv[6] #Current tan value
Sbma = sys.argv[7] #Current sin value
Process = str(sys.argv[8]) #Process running and name of folder
results_folder = str(sys.argv[9])
print(Process)

def the_main_event():
    """
    Uses 'make_input' and global arguments to read the basecard in and edit 
    the text from it to create the runcard that will be used by MadGraph.

    Parameters
    ---------
    None : 'the_main_event' uses global arguments only. These are 'basecardpath
        ', 'Tb', 'Sbma', 'Hcm', 'old_card' and 'Process'

    Returns
    -------
    None : creates and writes the runcard
    """
    
    print('inputcard_editor running...')
    with open(basecardpath,'r') as old_card:
        #stores string of old_card ready for editing
        text = old_card.read()
        make_input(Tb, Sbma, Hcm, text, Process)
        print('about to run madgraph ' + str(Tb) + ' ' + str(Sbma))

def make_input(Tanb, Sinbma, Hcpm, datatext, procname):
    """
    Creates the "runcard" for MadGraph. Takes inputs, tan_beta, sin(beta-alpha), 
    the desired mass for the charged higgses and a string of text to name the new 
    runcard with.

    Parameters
    ----------
    Tanb : float
    tan(beta) value of point
    Sinbma : float
        sin(beta-alpha) value of point
    Hcpm : mass of the charged Higgs
    datatext : string
        contents of basecard in string form
    procname : string
        name of the process to be run
    
    Returns
    -------
    N/A : writes the newly created file to be used as the 'runcard'
    """
    
    with open(runcardpath, 'w') as new_card:
        #simulation card, the .txt file that gets fed to madgraph
        print(runcardpath)

        sim_card = datatext.replace('TBV', str(Tanb))
        sim_card = sim_card.replace('SBMAV', str(Sinbma))
        sim_card = sim_card.replace('HMV', str(Hm))
        sim_card = sim_card.replace('AMV', str(Am))
        sim_card = sim_card.replace('HCMV', str(Hcm))
        sim_card = sim_card.replace('NAME', results_folder + str(procname) + 
                                    '_' + str(Tanb) + '_' + str(Sinbma))
        print(results_folder + str(procname) + '_' + str(Tanb) + '_' + str(Sinbma))

        new_card.write(sim_card) # saves new txt file for madgraph

the_main_event()
