# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import sys
import subprocess
import os
import numpy as np
import getopt
from Bio.PDB import PDBParser, PDBIO, Chain, Residue
import matplotlib.pyplot as plt
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
    #aminoacids =["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    aminoacids =["A","I","L","M","F","W","Y","V","C","G","P","S","T","N","Q","R","H","D","E","K"]
    #aminoacids =["A","I","L"]

    #runs crankpep for each of the different aminoacids - first , it craets the 20 sequences to be tested
    for res in aminoacids:
        global new                     # Call new as a global variable for the using in next function
        new = list(input_seq)
        new[result]=res                # Point change the X into one residue in "aminoacids"
        sequences.append(''.join(new)) # Add the mutated "new" sequence into the "sequences" list as a single element
    global residues                    # Define residues as global variable, originally named sequences
    residues = sequences               # Give sequences content to the variable-"residues"

def adcp_run():
    global list_of_lists 
    global CA_coordinates_final
    CA_coordinates_final=[]
    list_of_lists=[]
    #run crankpep for all 20 sequences
    for seq in residues:
        print("Sequence tested now is: " + seq)
        #creates a varaible to output the files inside different folders
        outputname_ADCP=seq+"/"+seq
        #creates a folder for each of the sequences
        subprocess.Popen(["mkdir",seq]).communicate()
        path_to_output_file =outputname_ADCP+".txt"
        myoutput = open(path_to_output_file,'w+') # equal to open seq.txt then write down
        # the path to adcp should be fixed and sorted with a enviromental variable such as ADCPHOME
        tmpname = "tmp_"+seq
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
        subprocess.Popen(["rm" , "-r" ,  tmpname]).communicate()    #delete tmp folders 
        def structure_correct():
            CA_XYZ=[]
            #calls biopdb to fix the broken pdb results for the top1 poses per sequence
            io = PDBIO()
            #fix it by loading to biopdb and printing
            pdb = PDBParser(QUIET=True).get_structure("UGLY", outputname_ADCP+"_ranked_1.pdb")
            io.set_structure(pdb)        
            io.save(outputname_ADCP+"_ranked_1_corrected.pdb")
            chains = list(pdb[0].get_chains())
            residue = list(chains[0].get_residues())
            for i in range(len(residue)):
                atoms = list(residue[i].get_atoms())
                for k in range(len(atoms)):
                    if atoms[k].get_name()=="CA":
                        XYZ=[]
                        XYZ.append(atoms[k].get_vector()[0])
                        XYZ.append(atoms[k].get_vector()[1])
                        XYZ.append(atoms[k].get_vector()[2])
                        CA_XYZ.append(XYZ)
            #reads the output files with the energy to rank poses
            return CA_XYZ
        
        CA_coordinates_final.append(structure_correct())
        
        def extract_score():
            file = open(path_to_output_file,'r')                                        # open the file.txt in a read mode
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
            list_of_lists.append([])                                                    # vector <vector>
            list_of_lists[len(list_of_lists)-1].append(aminoacids[residues.index(seq)]) # append the point changed residue name into the first element in list_of_lists
            list_of_lists[len(list_of_lists)-1].append(seq)                             # append the changed peptide sequence into the second element in list_of_lists
            list_of_lists[len(list_of_lists)-1].append(float(line_list[1]))                   # append the best scores into the third element in list_of_lists
        
        extract_score()

def cross_rmsdmatrix():
    CROSS_RMSD = []
    for peptide_i in CA_coordinates_final:
       RMSD_peptidei= []
       for peptide_k in CA_coordinates_final:
           RMSD = 0
           for CA in range(len(CA_coordinates_final[0])):
               RMSD = RMSD + sqrt((((peptide_i[CA][0] - peptide_k[CA][0])**2)+((peptide_i[CA][1] - peptide_k[CA][1])**2)+((peptide_i[CA][2] - peptide_k[CA][2])**2))/3)
           RMSD_peptidei.append(RMSD/len(CA_coordinates_final[0]))
       CROSS_RMSD.append(RMSD_peptidei)
    fig, ax = plt.subplots()
    plt.imshow(CROSS_RMSD, cmap='hot', interpolation='nearest')
    ax.set_xticklabels(aminoacids)
    ax.set_yticklabels(aminoacids)
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)
    ax.set_xticks(np.arange(len(aminoacids)))
    ax.set_yticks(np.arange(len(aminoacids)))

    plt.savefig('corssrmsd.png')
    
    
def sort_list():
    #sort the list of 20 residues
    sorted_multi_list = sorted(list_of_lists, reverse=False,key=lambda x: x[2])
    sorted_multi_array = np.array(sorted_multi_list)
    #prints in the screen the sorted list of 20 residues 
    for peptide in sorted_multi_list:
        print(peptide[0]," ", peptide[1], " ", peptide[2] , " kcal/mol")


# peppr_sup() This is for testing whether def is working.
