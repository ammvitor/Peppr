from pyrsistent import v
from PDB_parser_sulphotyrosine_2 import PDBfile
import os
import subprocess
from Bio.PDB import PDBParser, PDBIO, Chain, Residue


tyrosynetobemod = 5                                                                 # point out the selected tyrosine which to be modified, it also can get the value from the usr's input

with open('outlist.dat',"r") as list_pdb:                                           # 'with open(arg1) as arg2' don't need to file.close(); format of outlist.dat: PDAYL /home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/test/test2/PDAYL/PDAYL_ranked_1_corrected.pdb
    Lines = list_pdb.readlines()
    count = 0
    # Strips the newline character
    pdbfilenames_tocat =[]
    mol2filenames_tocat =[]
    cycpdbfilenames_tocat=[]
    cycmol2filenames_tocat=[]
    for line in Lines:                                                              # read the seq and path of $sequence_rank_1_corrected.pdb one by one
        
        linesplit=line.split()                                                      # to split the sequence and path, here 'line' is actually one line in the 'Lines'
        x=PDBfile()                                                                 # creat an object 'x'
        x.PDBreader(linesplit[1])                                                   # call the method 'PDBreader' within class 'PDBfile', and the input is the path of $sequence_rank_1_corrected.pdb; finally can generate a series of 'charactors' of the object 'x'
        
        for residue_ndx in range(0,len(x.residue_name)):
            if(x.residue_name[residue_ndx] == "TYR" and x.residue_index[residue_ndx] == tyrosynetobemod ):
                for i in range(4):
                    print(i, "  "+linesplit[1])  # for checking errors
                    try:                                                            # to deal with the erro ‘matrix is numerically singular’
                        x.PDBreader(linesplit[1]) 
                        x.addPO3_toTYR(i,tyrosynetobemod)                           # 'i' is used for setting the addition values, 'tyrosynetobemod' is used for locate the specific TYR in method 'add_PO3'                                              
                        x.PDBwriter(linesplit[0]+"_PO3.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_PO3.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_PO3.pdb",linesplit[0]+"_PO3.mol2",tyrosynetobemod,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_PO3.mol2")
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_PO3_N-C_cyc.pdb",tyrosynetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_PO3_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_PO3_N-C_cyc.pdb", linesplit[0]+"_PO3_N-C_cyc.mol2",tyrosynetobemod,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_PO3_N-C_cyc.mol2")
                        os.system('rm ' +linesplit[0]+"_PO3_N-C_cyc.pdb")
                        
                        x.PDBreader(linesplit[1])
                        x.addSO3_toTYR(i,tyrosynetobemod)                                           # 'i' is used for setting the addition values in method 'add_SO3'
                        x.PDBwriter(linesplit[0]+"_SO3.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_SO3.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_SO3.pdb",linesplit[0]+"_SO3.mol2",tyrosynetobemod,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3.mol2")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3.mol2")
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_SO3_N-C_cyc.pdb",tyrosynetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_SO3_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_SO3_N-C_cyc.pdb", linesplit[0]+"_SO3_N-C_cyc.mol2",tyrosynetobemod,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_SO3_N-C_cyc.mol2")
                        os.system('rm '+linesplit[0]+"_SO3_N-C_cyc.pdb")
                    except:
                        pass
                break
                
        x.PDBreader(linesplit[1])
        x.PDBwriter(linesplit[0]+".pdb")
        pdbfilenames_tocat.append(linesplit[0]+".pdb")
        x.obabel_mol2_em(linesplit[0]+".pdb",linesplit[0]+".mol2",tyrosynetobemod,"NOTPTM")
        mol2filenames_tocat.append(linesplit[0]+".mol2")
        
        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_cyc.pdb",tyrosynetobemod)
        cycpdbfilenames_tocat.append(linesplit[0]+"_cyc.pdb")
        x.obabel_mol2_cyc(linesplit[0]+"_cyc.pdb", linesplit[0]+"_cyc.mol2",tyrosynetobemod,"NOTPTM")
        cycmol2filenames_tocat.append(linesplit[0]+"_cyc.mol2")
        os.system('rm '+linesplit[0]+"_cyc.pdb")
        
pdbfilenames_tocat = list(set(pdbfilenames_tocat))                                  # remove duplicate values
mol2filenames_tocat = list(set(mol2filenames_tocat))                                # remove duplicate values
cycmol2filenames_tocat = list(set(cycmol2filenames_tocat))                          # remove duplicate values


with open('pdb_lib_done.pdb',"w") as pdb_lib_done:
    count = 0
    # Strips the newline character
    for pdbname in pdbfilenames_tocat:
        count += 1
        with open(pdbname,"r") as pdbnamefile:
            pdb_lib_done.write(str("MODEL "+str(count)+"\n"))
            Lines = pdbnamefile.readlines()
            for line in Lines:
                pdb_lib_done.write(line)
        os.system('rm '+pdbname)

                
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
        
with open('cyclic_lib_done.mol2',"w") as cyclic_lib_done:                       # 我加的
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
        

