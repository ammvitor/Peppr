"""
import peppr_sup_HelpMessageAdded_def_main_ 

peppr_sup_HelpMessageAdded_def_main_.peppr_sup()
"""
import Peppr_functions  #like the head file in c++
import sys
import subprocess
import os
import getopt
import numpy as np
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

#Peppr_functions.set_sequence()
Peppr_functions.adcp_run()             # updated, adcp_run() includes "set_sequence", "structure_correct", "extract_score", "sort_list"

print("Congratulations! It's finished!")

