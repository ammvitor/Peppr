# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys
import subprocess
import os
import numpy as np
import getopt
from Bio.PDB import PDBParser, PDBIO, Chain, Residue
#you need to have biopython and autodock vina crankpep - export ADCPHOME where you installed it


#declaring the input variables
args=sys.argv[1:]
input_seq = ''
input_receptor =''
percentage = ''
cores = ''
input_ntrys = ''
input_Replicas = ''

#optarg for the input example-c 2 -i PLDXPAL -t done.trg -p 50 -N 6 -n 80000
try:
   opts, args = getopt.getopt(args,"h:i:t:p:N:n:c:",["help","input_seq =",
                                    "receptor =",
                                    "percentage ="
                                    "replicas =",
                                    "steps =",
                                    "cores ="])
except getopt.GetoptError:
   print ('test.py -i <inputfile> -o <outputfile>')
   sys.exit(2)
 
for opt, arg in opts:
   if opt == '-h':
      print ('args.py -i <inputfile> -o <outputfile>')
      sys.exit()
   elif opt in ("-i", "--input_seq"):
      input_seq = arg
   elif opt in ("-t", "--receptor"):
      input_receptor = arg
   elif opt in ("-p", "--percentage"):
      percentage = arg
   elif opt in ("-N", "--replicas"):
      input_Replicas = arg
   elif opt in ("-n", "--steps"):
      input_ntrys = arg
   elif opt in ("-c", "--cores"):
      cores = arg
      
adcphome = os.environ['ADCPHOME']+"/bin/adcp"
res = 0                                # Prepare this variable for a global using purpose

def set_sequence():   
    ##set up the sequences 
    sequences = []
    #finds the index of the X in the inpout sequence
    result = input_seq.find("X")       # Should be used to find the index number of "X", so here result is a number.
    global aminoacids                  # # Call aminoacids as a global variable for the using in next function
    aminoacids =["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

    #runs crankpep for each of the different aminoacids - first , it craets the 20 sequences to be tested
    for res in aminoacids:
        global new                     # Call new as a global variable for the using in next function
        new = list(input_seq)
        new[result]=res                # Point change the X into one residue in "aminoacids"
        sequences.append(''.join(new)) # Add the mutated "new" sequence into the "sequences" list as a single element
    global residues                    # Define residues as global variable, originally named sequences
    residues = sequences               # Give sequences content to the variable-"residues"

def structure_correct(outputname):
    #calls biopdb to fix the broken pdb results for the top1 poses per sequence
    io = PDBIO()
    #fix it by loading to biopdb and printing 
    total_number_of_output=["1","2","3","4","5","6","7","8","9","10"]
    for rank in total_number_of_output:
        try:
            pdb = PDBParser(QUIET=True).get_structure("UGLY", outputname+"_ranked_"+rank+".pdb")
            io.set_structure(pdb)        
            io.save(outputname+"_ranked_"+rank+"_corrected.pdb")
            #reads the output files with the energy to rank poses
        except:
            pass

    
def extract_score(pathtooutput,listofsequences,sequences):
    file = open(pathtooutput,'r')                                        # open the file.txt in a read mode
    Lines = file.readlines()                                                    # output the content in file.txt line by line into "Lines"
    count =0
    #skips the header of the file 
    for line in Lines:
        if line[0] != "-":
            count = count + 1                                                   # find the best score line's previous line number
        else:
            break
    #takes the first line - lowest energy and rank per residue
    stripped_line = Lines[count+1].strip()                                      # count+1 means to locate to the line with best scores; .strip() means only extract this line
    line_list = stripped_line.split()                                           # split() means to treat the extracted line into pices and save as elements in "line_list"
    listofsequences.append([])                                                    # vector <vector>
    listofsequences[len(listofsequences)-1].append(aminoacids[residues.index(sequences)]) # append the point changed residue name into the first element in list_of_lists
    listofsequences[len(listofsequences)-1].append(sequences)                             # append the changed peptide sequence into the second element in list_of_lists
    listofsequences[len(listofsequences)-1].append(line_list[1])                    # append the best scores into the third element in list_of_lists
        

def adcp_run():   
    list_of_lists=[]
    #run crankpep for all 20 sequences                                                                 # for global Structure Correct
    for seq in residues:

        print("Sequence tested now is: " + seq)
        #creates a varaible to output the files inside different folders
        # global outputname_ADCP                                                  # for global Structure Correct
        outputname_ADCP=seq+"/"+seq
        #creates a folder for each of the sequences
        subprocess.Popen(["mkdir",seq]).communicate()
                                                      # for global Structure Correct
        path_to_output_file =outputname_ADCP+".txt"
        myoutput = open(path_to_output_file,'w+') # equal to open seq.txt then write down
        # the path to adcp should be fixed and sorted with a enviromental variable such as ADCPHOME
        p = subprocess.Popen([adcphome,
                              "-t",input_receptor,"-s",
                              ''.join(new),  # Question: is this "new" = to the previous function's "new"?
                              "-N",
                              input_Replicas,
                              "-n",
                              input_ntrys,
                              "-p",
                              percentage,
                              "-o",
                              outputname_ADCP,
                              "-c",
                              cores,
                              "-O"],stdout=myoutput).communicate()
        structure_correct(outputname_ADCP)
        extract_score(path_to_output_file,list_of_lists,seq)
        sort_list(list_of_lists)
        



def sort_list(listofsequences): # positive numbers should be fixed
    #sort the list of 20 residues
    sorted_multi_list = sorted(listofsequences, reverse=True,key=lambda x: x[2])
    sorted_multi_array = np.array(sorted_multi_list)
    # print(sorted_multi_array)
    print("Results:")
    #prints in the screen the sorted list of 20 residues 
    for peptide in sorted_multi_list:
        print("Score ", peptide[0]," ", peptide[1], " ", peptide[2] + " kcal/mol")
# print("Congratulations! It's finished!")

# peppr_sup() This is for testing whether def is working.