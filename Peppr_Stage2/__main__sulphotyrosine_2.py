from pyrsistent import v
from PDB_parser_sulphotyrosine_2 import PDBfile
import os
import subprocess
from Bio.PDB import PDBParser, PDBIO, Chain, Residue


tyrosynetobemod = 5                                                                 # can we make it as a user input? point out the selected tyrosine which to be modified, it also can get the value from the usr's input

with open('outlist.dat',"r") as list_pdb:                                           # 'with open(arg1) as arg2' don't need to file.close(); format of outlist.dat: PDAYL /home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/test/test2/PDAYL/PDAYL_ranked_1_corrected.pdb
    Lines = list_pdb.readlines()
    count = 0
    # Strips the newline character
    pdbfilenames_tocat =[]                                                          # create library for pdb files(20*natural,20*PTM-so3, 20*PTM-po3)
    mol2filenames_tocat =[]                                                         # create library for mol2 files(20*natural,20*PTM-so3, 20*PTM-po3)        
    cycpdbfilenames_tocat=[]
    cycmol2filenames_tocat=[]                                                       # create library for cyclic mol2 files(20*natural,20*PTM-so3, 20*PTM-po3)
    for line in Lines:                                                              # read the seq and path of $sequence_rank_1_corrected.pdb one by one
        
        linesplit=line.split()                                                      # to split the sequence and path, here 'line' is actually one line in the 'Lines'
        x=PDBfile()                                                                 # creat an object 'x'
        x.PDBreader(linesplit[1])                                                   # call the method 'PDBreader' within class 'PDBfile', and the input is the path of $sequence_rank_1_corrected.pdb; finally can generate a series of 'charactors' of the object 'x'
        
        for residue_ndx in range(0,len(x.residue_name)):
            if(x.residue_name[residue_ndx] == "TYR" and x.residue_index[residue_ndx] == tyrosynetobemod ):
                count = 0 
                flag = "FALSE"
                while(count < 4 and flag == "FALSE"):
                # for i in range(4):
                    print(count, "  "+linesplit[1])  # for checking errors
                    try:                                                            # to deal with the erro ‘matrix is numerically singular’
                        x.PDBreader(linesplit[1])                                   # PDBreader method can clean up all previous contents, and re-generate the charactors for object x
                        x.addPO3_toTYR(count,tyrosynetobemod)                           # 'i' is used for setting the addition values which will be used for calculat PTM values, 'tyrosynetobemod' is used for locate the specific TYR in method 'add_PO3'                                              
                        x.PDBwriter(linesplit[0]+"_PO3.pdb")                        # add po3 to target residue
                        pdbfilenames_tocat.append(linesplit[0]+"_PO3.pdb")          # add file name to list which will be used to create the pdb library
                        x.obabel_mol2_em(linesplit[0]+"_PO3.pdb",linesplit[0]+"_PO3.mol2",tyrosynetobemod,"PTM") # Convert pdb file to mol2 file; "PTM" is to tell method "obabel_mol2_em", to call the checking progress, remove the wrong bond between side chains 
                        mol2filenames_tocat.append(linesplit[0]+"_PO3.mol2")        # add mol2 file name to list which will be used to create the mol2 library
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_PO3_N-C_cyc.pdb",tyrosynetobemod) # call cyclic method to add CONECT information at the bottom of pdb files, to guide the process of pdb to mol2
                        cycpdbfilenames_tocat.append(linesplit[0]+"_PO3_N-C_cyc.pdb")   # add file name to list which will be used to create the cyclic pdb library
                        x.obabel_mol2_cyc(linesplit[0]+"_PO3_N-C_cyc.pdb", linesplit[0]+"_PO3_N-C_cyc.mol2",tyrosynetobemod,"PTM") # Convert pdb file to mol2 file, remove wrong bond record
                        cycmol2filenames_tocat.append(linesplit[0]+"_PO3_N-C_cyc.mol2") # add mol2 file name to list which will be used to create the cyclic mol2 library
                        os.system('rm ' +linesplit[0]+"_PO3_N-C_cyc.pdb")           # if there is any problem just comment this line and check the pdb files afterwards; remove the preparation pdb file
                        
                        x.PDBreader(linesplit[1])                                   # Basicly the comments of these part are the same as above part, the only difference is this part is for adding SO3 ptm
                        x.addSO3_toTYR(count,tyrosynetobemod)                                          
                        x.PDBwriter(linesplit[0]+"_SO3.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_SO3.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_SO3.pdb",linesplit[0]+"_SO3.mol2",tyrosynetobemod,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3.mol2")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3.mol2")
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_SO3_N-C_cyc.pdb",tyrosynetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_SO3_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_SO3_N-C_cyc.pdb", linesplit[0]+"_SO3_N-C_cyc.mol2",tyrosynetobemod,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_SO3_N-C_cyc.mol2")
                        os.system('rm '+linesplit[0]+"_SO3_N-C_cyc.pdb")            # if there is any problem just comment this line and check the pdb files afterwards
                        flag = "TRUE"
                    except:
                        count += 1
                        pass
                    
                break
                
        x.PDBreader(linesplit[1])                                                   # re-creat the object x with all charactors un-modified
        x.PDBwriter(linesplit[0]+".pdb")                                            # creat the natural pdb file
        pdbfilenames_tocat.append(linesplit[0]+".pdb")                              # add file name to list which will be used to create the pdb library
        x.obabel_mol2_em(linesplit[0]+".pdb",linesplit[0]+".mol2",tyrosynetobemod,"NOTPTM") # Convert pdb file to mol2 file; "NOTPTM" tells the obabel_mol2_em method not to remove any bonds
        mol2filenames_tocat.append(linesplit[0]+".mol2")                            # add mol2 file name to list which will be used to create the mol2 library
        
        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_cyc.pdb",tyrosynetobemod)             # add all side chains' CONECT information to pdb file 
        cycpdbfilenames_tocat.append(linesplit[0]+"_cyc.pdb")                       # add file name to list which will be used to create the cyclic pdb library
        x.obabel_mol2_cyc(linesplit[0]+"_cyc.pdb", linesplit[0]+"_cyc.mol2",tyrosynetobemod,"NOTPTM") # "NOTPTM" for telling method don't need to find the atoms from PTM, only check the wrong bond between normal  side chains.
        cycmol2filenames_tocat.append(linesplit[0]+"_cyc.mol2")                     # add mol2 file name to list which will be used to create the cyclic mol2 library
        os.system('rm '+linesplit[0]+"_cyc.pdb")                                    # if there is any problem just comment this line and check the pdb files afterwards
        
pdbfilenames_tocat = list(set(pdbfilenames_tocat))                                  # remove duplicate values
mol2filenames_tocat = list(set(mol2filenames_tocat))                                # remove duplicate values
cycmol2filenames_tocat = list(set(cycmol2filenames_tocat))                          # remove duplicate values


with open('pdb_lib_done.pdb',"w") as pdb_lib_done:                                  # build the pdb library; with__open__as grammar don't need to use f.close
    count = 0
    # Strips the newline character
    for pdbname in pdbfilenames_tocat:                                              # iterate over the entire pdb name list    
        count += 1
        with open(pdbname,"r") as pdbnamefile:                                      # pdbname is the name of each pdb file, here open them one by one and read them as "pdbnamefile"
            pdb_lib_done.write(str("MODEL "+str(count)+"\n"))                       # "pdb_lib_done" is the new created library file, first line is "MODEL 1"
            Lines = pdbnamefile.readlines()                                         # pdbnamefile is each pdb file with atom informations, and make a list named "Lines" include line by line 
            for line in Lines:                                                      # iterate over the entire "Lines" list
                pdb_lib_done.write(line)                                            # with the iteration, write down line by line.
        os.system('rm '+pdbname)                                                    # empty the pdbname file

                
with open('mol2_lib_done.mol2',"w") as mol2_lib_done:
    count = 0
    # Strips the newline character
    for mol2name in mol2filenames_tocat:
        with open(mol2name,"r") as mol2namefile:
            Lines = mol2namefile.readlines()
            for line in Lines:
                mol2_lib_done.write(line)
        os.system('rm '+mol2name)
        # myoutput = open("ouput.dat",'w+')
        # p = subprocess.Popen(["/home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/smina.static",
        #                      "-l","mol2_lib_done.mol2",
        #                      "--autobox_ligand","center.mol2",
        #                      "-r","receptor.pdb",
        #                      "--autobox_add","6",
        #                      "-o","results.sdf"],stdout=myoutput).communicate()
        
with open('cyclic_lib_done.mol2',"w") as cyclic_lib_done:                       
    count=0
    for cycname in cycmol2filenames_tocat:
        with open(cycname,"r") as mol2namefile:
            Lines = mol2namefile.readlines()
            for line in Lines:
                cyclic_lib_done.write(line)
        os.system('rm '+cycname)
        # myoutput = open("ouput.dat",'w+')
        # p = subprocess.Popen(["/home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/smina.static",
        #                      "-l","cyclic_lib_done.mol2",
        #                      "--autobox_ligand","center.mol2",
        #                      "-r","receptor.pdb",
        #                      "--autobox_add","6",
        #                      "-o","results.sdf"],stdout=myoutput).communicate()
        

