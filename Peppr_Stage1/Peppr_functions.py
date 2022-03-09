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
args=sys.argv[1:]                                                                   # get arguments from usr's input; filename is sys.arg[0], so args start from [1:]
input_seq = ''
input_receptor =''
percentage = ''
cores = ''
input_ntrys = ''
input_Replicas = ''

#optarg for the input example-c 2 -i PLDXPAL -t done.trg -p 50 -N 6 -n 80000
try:
   opts, args = getopt.getopt(args,"h:i:t:p:N:n:c:",["help","input_seq =",          # getopot.getopt(sys.arg, short_option‘-h,-i,-t,-p,etc’, long_option'--help,--input_seq,--receptor)
                                    "receptor =",                                   # with usr's input e.g -i PLDXPAL -c2, the 'getopt' function can grab them and save them seprately into 'opts' and 'args'
                                    "percentage ="
                                    "replicas =",
                                    "steps =",
                                    "cores ="])
except getopt.GetoptError:
   print ('test.py -i <inputfile> -o <outputfile>')
   sys.exit(2)                                                                      # Exiting the program raises a SystemExit exception, 0 means normal exit, and others are abnormal exits.
 
for opt, arg in opts:                                                               # Generate several pairs of value, e.g: opr,arg = -i,PLDXPAL
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
      
adcphome = os.environ['ADCPHOME']+"/bin/adcp"                                       # call the system enviroment variable'ADCPHOME' which is pre-defined by usr, then add /bin/adcp to point out the excutable file.
res = 0                                                                             # Prepare this variable for a global using purpose

def set_sequence(aminoacids_list):                                                  # The aminaacids_list="A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"
    ##set up the sequences 
    sequences = []
    #finds the index of the X in the inpout sequence
    result = input_seq.find("X")                                                    # Should be used to find the index number of "X", 'result' = index[X].

    #runs crankpep for each of the different aminoacids - first , it craets the 20 sequences to be tested
    for res in aminoacids_list:                                                     # Here the 'res' = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
       # Call new as a global variable for the using in next function
        new = list(input_seq)                                                       # function 'list' transfer the 'input_seq' into a list
        new[result]=res                                                             # Point change the X into one residue in "aminoacids"
        sequences.append(''.join(new))                                              # Add the mutated "new" sequence into the "sequences" list as a single element
    return sequences                                                                # 'sequence' is a list which includes 20 sequences with the X changed into 20 normal residues.

def structure_correct(outputname):                                                  # To fix the ugly structure and save the fixed structures, finally return a path; e.g 'outputname' = PLDYS/PLDYS
    #calls biopdb to fix the broken pdb results for the top1 poses per sequence
    io = PDBIO()                                                                    # call the class 'PDBIO()' within BIO.PDB
    #fix it by loading to biopdb and printing 
    total_number_of_output=["1","2","3","4","5","6","7","8","9","10"]
    for rank in total_number_of_output:
        try:
            pdb = PDBParser(QUIET=True).get_structure("UGLY", outputname+"_ranked_"+rank+".pdb")    # 'get_structure(self, ID, file)' is a method in class 'PDBParser', 
            io.set_structure(pdb)                                                                   # 'set_structure(self, pdb_object)' is a method in class 'PDBIO'    
            io.save(outputname+"_ranked_"+rank+"_corrected.pdb")                                    # 'save(self, file, select=_select, write_end=True, preserve_atom_numbering=False)' is a method in class 'PDBIO'
            #reads the output files with the energy to rank poses
            
        except:
            pass
        
    return(os.getcwd()+"/"+outputname+"_ranked_"+"1"+"_corrected.pdb")              # only return the path of 'outputname'+'ranked_1_corrected.pdb' file name

    
def extract_score(pathtooutput,listofsequences,sequences):                          # 'listofsequences' is an empty list while calling 'extract_score'; 'sequences' is one of the 20 sequences
    aminoacids =["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] #irts preferable to have a variable within a smaller scope or as an argument
    residues = set_sequence(aminoacids)                                             # 'residues' = [the 20 sequence list]
    file = open(pathtooutput,'r')                                                   # open the file.txt in a read mode, e.g PLDYS/PLDYS.txt
    Lines = file.readlines()                                                        # output the content in file.txt line by line into "Lines"
    count =0                                                                        # To get the previous line number of score value's line
    #skips the header of the file 
    for line in Lines:
        if line[0] != "-":
            count = count + 1                                                       # find the best score line's previous line number
        else:
            break
    #takes the first line - lowest energy and rank per residue
    stripped_line = Lines[count+1].strip()                                          # count+1 means to locate to the line with best scores; .strip() means only extract this line
    line_list = stripped_line.split()                                               # split() means to treat the extracted line into pices and save as elements in "line_list"
    listofsequences.append([])                                                      # vector <vector>
    listofsequences[len(listofsequences)-1].append(aminoacids[residues.index(sequences)]) # find the 'sequences' rank in 20 sequences, then pass this index value to aminoacids to extract the point changed residue name; append the point changed residue name into the first element in list_of_lists
    listofsequences[len(listofsequences)-1].append(sequences)                       # append the changed peptide sequence into the second element in list_of_lists
    listofsequences[len(listofsequences)-1].append(line_list[1])                    # append the best scores into the third element in list_of_lists
        

def adcp_run():   
    number = 0
    print(adcphome)
    aminoacids =["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    list_of_lists=[]
    #run crankpep for all 20 sequences                
    residues = set_sequence(aminoacids)
    topmols=[]
    topmols_seq=[]
    for seq in residues:
        number+=1
        print("Sequence tested now is: " + seq)
        #creates a varaible to output the files inside different folders 
        # global outputname_ADCP                                                    
        outputname_ADCP=seq+"/"+seq                                                 # 'seq+"/"+seq' the first 'seq' is the folder named 'seq', the second 'seq' means the output file's prefix
        # print(seq)
        #creates a folder for each of the sequences
        subprocess.Popen(["mkdir",seq]).communicate()                               # mkdir 20 folder with name is their sequence; 'Popen(args)' is a method in class 'subprocess'; execute operating system-level commands in Python code, 'communicate' will put the output in memory, not in the pipe
                                                      # for global Structure Correct
        path_to_output_file =outputname_ADCP+".txt"                                 # used in 'extract_score' e.g: PLADY/PLADY.txt
        myoutput = open(path_to_output_file,'w+') # equal to open seq.txt then write down
        # the path to adcp should be fixed and sorted with a enviromental variable such as ADCPHOME
        p = subprocess.Popen([adcphome,
                              "-t",input_receptor,"-s",
                              seq,  # Question: is this "new" = to the previous function's "new"?
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
       
        topmols.append(structure_correct(outputname_ADCP))                          # get the PATH for each $sequence_rank_1_corrected.pdb
        topmols_seq.append(seq)                                                     # get the sequence for each $sequence_rank_1_corrected.pdb
        extract_score(path_to_output_file,list_of_lists,seq)                        # 'list_of_lists' is 'listofsequences' in function 'extract_score', so after calling 'extract_score', 'list_of_lists' should have lots of information, e.g R PLDRYS  -5.4 
        sort_list(list_of_lists,number)                                             # if 'number' >= 20 then 'sort_list' will run.
        subprocess.Popen(["rm","-r","tmp_"+seq]).communicate()
    listname = open("outlist.dat","w")
    for ndx in  range(0,len(topmols)):
        listname.write(topmols_seq[ndx]+" "+topmols[ndx]+"\n")                      # this is to out put one sequence and one path of each sequence 
    listname.close()
()
        



def sort_list(listofsequences,loops): # positive numbers should be fixed
    if(loops >= 20):
        #sort the list of 20 residues
        sorted_multi_list = sorted(listofsequences, reverse=True,key=lambda x: x[2])    # 'key=lambda x: x[2]' is to tell 'sort' only sort the third element within 'list_of_lists'
        sorted_multi_array = np.array(sorted_multi_list)
        # print(sorted_multi_array)
        print("Results:")
        #prints in the screen the sorted list of 20 residues 
        for peptide in sorted_multi_list:
            print("Score ", peptide[0]," ", peptide[1], " ", peptide[2] + " kcal/mol")
# print("Congratulations! It's finished!")

# peppr_sup() This is for testing whether def is working.