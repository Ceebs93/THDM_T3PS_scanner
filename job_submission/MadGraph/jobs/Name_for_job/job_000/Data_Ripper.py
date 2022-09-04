# Data Ripper
import os
import shutil
import numpy as np
import sys
import pandas as pd

#Note: This version of Data_Ripper_iridis is set to look for tan(beta) values by scanning for\'# tb' in the
# .lhe file. This is the setup for the type-I model in use, it will NOT work for type-II as this uses a
# different name for the variable.

Process = str(sys.argv[1]) #Name of results folder
uncomp_dir = "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/Name_for_job/job_000/" + str(Process) #Folder containing all the run data
output_dir = "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/Name_for_job/job_000" + "/Data_Files"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

######################################################################################################
def get_uncomp_file(Uncomp_dir):
    """ retrieves all the banner.txt/unweighted_events files from file folders and creates a list"""
    uncomped_files = [] # list of zipped lhe files
    print("Data_Ripper is looking in the following places: ")
    print(uncomp_dir)
    print("                                                      ")
    path_to_file = uncomp_dir + "/Events/run_01/run_01_tag_1_banner.txt"
    path_to_lhe = uncomp_dir + "/Events/run_01/unweighted_events.lhe.gz"

    if os.path.exists(path_to_file) and os.path.exists(path_to_lhe) :
        print("Banner-path found")
        file_storage = "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/Name_for_job/job_000" + "/Processed_results/" + str(Process) # name of new folder for storing
        # storing the data for later use
        if not os.path.exists(file_storage):
            os.makedirs(file_storage)
        
        shutil.copy(path_to_file, file_storage)# copying the file
        shutil.copy(path_to_lhe, file_storage)# moving the lhe file
        txt_stored = str(file_storage) + "/run_01_tag_1_banner.txt"
        uncomped_files.append(txt_stored)# adding the path to the copied file to our list
    else:
        uncomped_files.append("empty")
    return uncomped_files                
####################################################################################################    

####################################################################################################
def variable_values(* args):
    """"takes values of tan(beta) and sin(beta-alpha) from banner.txt files and returns them as lists"""
    text_tan_data = []
    text_sin_data = []

    for a_list in args:
        for a_file in a_list:
            if os.path.isfile(a_file):
                with open(a_file,'r') as data_file:
                    # open the current file in the list and go on to read each line
                    lines = data_file.readlines()
                    for line in lines:
                        if '# tb' in line:
                            print("Data_Ripper has found tanbeta")
                            # we find the line containing tanbeta and add the value to our list
                            temp_line_holder = line.split()
                            print("Found a tanbeta value of: " +  str(temp_line_holder[1]))
                            text_tan_data.append(temp_line_holder[1])
                    
                        elif '# sinbma' in line:
                            # this works the same way as the previous part but for sin(beta-alpha)
                            temp_line_holder = line.split()
			    print("Data_Ripper has found  a sinbma value of: " + str(temp_line_holder[1]))
                            text_sin_data.append(temp_line_holder[1])
                        else:
                            continue
            else:
                continue

        return(text_tan_data, text_sin_data)	 
###################################################################################################

##################################################################################################
def Xsection_values(file_set):
    """"takes the cross-section and errors from files"""
    X_Section = list()

    for a_file in file_set:
        if os.path.isfile(a_file):
            with open(a_file,'r') as data_file:
                # open the current file in the list and go on to read each line
                lines = data_file.readlines()
                for line in lines:
                    if '#  Integrated weight (pb)  :' in line:
                        # we find the line containing the cross-section and add the value to our list
                        temp_line_holder = line.split()
			print("Data_Ripper found a x_section value of: " + str(temp_line_holder[5]))
                        X_Section.append(temp_line_holder[5])
                    else:
                        continue
        
    return(X_Section)
###################################################################################################    
    
File_List = get_uncomp_file(uncomp_dir)
print("Data_Ripper will look for sin and tan values in: " + str(uncomp_dir))
print("                                     ")
Tb_list, Sbma_list = variable_values(File_List)
X_sections = Xsection_values(File_List)
print("X_sections", X_sections)

data = [Tb_list, Sbma_list, [X_sections[0]]]
col_names = ["Tb", "Sbma", "X_section"]

nw_df = pd.DataFrame({'Tb':Tb_list, 'Sbma':Sbma_list, 'X_sections':X_sections})
nw_df.to_csv(str(output_dir) + "/" + str(Process) + ".csv", index=False)

print(len(Sbma_list), len(Tb_list), len(X_sections))
