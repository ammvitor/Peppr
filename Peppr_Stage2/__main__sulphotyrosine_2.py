from pyrsistent import v
from PDB_parser_sulphotyrosine_2 import PDBfile
import os
import subprocess
import traceback
from Bio.PDB import PDBParser, PDBIO, Chain, Residue


tyrosynetobemod = 5                                                                 # can we make it as a user input? point out the selected tyrosine which to be modified, it also can get the value from the usr's input
serinetobemod = 3
threoninetobemod = 0
cystocyc = 6                                                                        # indicate the residue_index for the Cys we want to cyclized

with open('outlist.dat',"r") as list_pdb:                                           # 'with open(arg1) as arg2' don't need to file.close(); format of outlist.dat: PDAYL /home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/test/test2/PDAYL/PDAYL_ranked_1_corrected.pdb
    Lines = list_pdb.readlines()
    count = 0
    # Strips the newline character
    pdbfilenames_tocat = []                                                          # create library for pdb files(20*natural,20*PTM-so3, 20*PTM-po3)
    mol2filenames_tocat = []                                                         # create library for mol2 files(20*natural,20*PTM-so3, 20*PTM-po3)        
    cycpdbfilenames_tocat = []
    cycmol2filenames_tocat = []                                                       # create library for cyclic mol2 files(20*natural,20*PTM-so3, 20*PTM-po3)
    chloroacetylate_cyclic_tocat = []                                                       # Attention! Based on the paper: PMID: 22419118 DOI: 10.1039/c2ob25306b  Though called chloroacetylate cyclic, but there is no Cl. 
    for line in Lines:                                                              # read the seq and path of $sequence_rank_1_corrected.pdb one by one
        
        linesplit=line.split()                                                      # to split the sequence and path, here 'line' is actually one line in the 'Lines'
        x=PDBfile()                                                                 # creat an object 'x'
        x.PDBreader(linesplit[1])                                                   # call the method 'PDBreader' within class 'PDBfile', and the input is the path of $sequence_rank_1_corrected.pdb; finally can generate a series of 'charactors' of the object 'x'
        
        for residue_ndx in range(0,len(x.residue_name)):
            if(x.residue_name[residue_ndx] == "TYR" and x.residue_index[residue_ndx] == tyrosynetobemod ):
                count = 0 
                flag = "FALSE"
                while(count < 9 and flag == "FALSE"):
                    print(count, "PO3_TYR  "+linesplit[1])  # for checking errors
                    try:                                                            # to deal with the erro ‘matrix is numerically singular’
                        x.PDBreader(linesplit[1])                                   # PDBreader method can clean up all previous contents, and re-generate the charactors for object x
                        x.addPO3_toTYR(count,tyrosynetobemod)                           # 'i' is used for setting the addition values which will be used for calculat PTM values, 'tyrosynetobemod' is used for locate the specific TYR in method 'add_PO3'                                               
                        x.PDBwriter(linesplit[0]+"_PO3_TYR.pdb")                        # add po3 to target residue
                        pdbfilenames_tocat.append(linesplit[0]+"_PO3_TYR.pdb")          # add file name to list which will be used to create the pdb library
                        x.obabel_mol2_em(linesplit[0]+"_PO3_TYR.pdb",linesplit[0]+"_PO3_TYR.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM") # Convert pdb file to mol2 file; "PTM" is to tell method "obabel_mol2_em", to call the checking progress, remove the wrong bond between side chains 
                        mol2filenames_tocat.append(linesplit[0]+"_PO3_TYR.mol2")        # add mol2 file name to list which will be used to create the mol2 library
                        
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_PO3_TYR_N-C_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod) # call cyclic method to add CONECT information at the bottom of pdb files, to guide the process of pdb to mol2
                        cycpdbfilenames_tocat.append(linesplit[0]+"_PO3_TYR_N-C_cyc.pdb")   # add file name to list which will be used to create the cyclic pdb library
                        x.obabel_mol2_cyc(linesplit[0]+"_PO3_TYR_N-C_cyc.pdb", linesplit[0]+"_PO3_TYR_N-C_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM") # Convert pdb file to mol2 file, remove wrong bond record
                        cycmol2filenames_tocat.append(linesplit[0]+"_PO3_TYR_N-C_cyc.mol2") # add mol2 file name to list which will be used to create the cyclic mol2 library
                        # try:
                        if(cystocyc != 0):
                            x.addchloroacetyl_toNt(count)                          #
                            x.chloroacetylate_n_cys_Cyclic_PDBwriter(linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc)
                            chloroacetylate_cyclic_tocat.append(linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.mol2")
                            x.obabel_mol2_cyc(linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.pdb",linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                            # cycmol2filenames_tocat.append(linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.mol2")
                        # except:
                        #     pass
                        os.system('rm ' +linesplit[0]+"_PO3_TYR_N-C_cyc.pdb" +" "+linesplit[0]+"_PO3_TYR_Nacetyl-Cys_cyc.pdb")           # if there is any problem just comment this line and check the pdb files afterwards; remove the preparation pdb file
                        flag = "TRUE"
                    except:
                        count += 1
                        pass
                break
        for residue_ndx in range(0,len(x.residue_name)):
            if(x.residue_name[residue_ndx] == "TYR" and x.residue_index[residue_ndx] == tyrosynetobemod ):
                count = 0 
                flag = "FALSE"
                while(count < 9 and flag == "FALSE"):
                    print(count, "SO3_TYR  "+linesplit[1])  # for checking errors
                    try:                                                            # to deal with the erro ‘matrix is numerically singular’                        
                        x.PDBreader(linesplit[1])                                   # Basicly the comments of these part are the same as above part, the only difference is this part is for adding SO3 ptm
                        x.addSO3_toTYR(count,tyrosynetobemod)                                                                 
                        x.PDBwriter(linesplit[0]+"_SO3_TYR.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_SO3_TYR.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_SO3_TYR.pdb",linesplit[0]+"_SO3_TYR.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3_TYR.mol2")
                        mol2filenames_tocat.append(linesplit[0]+"_SO3_TYR.mol2")
                        
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_SO3_TYR_N-C_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_SO3_TYR_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_SO3_TYR_N-C_cyc.pdb", linesplit[0]+"_SO3_TYR_N-C_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_SO3_TYR_N-C_cyc.mol2")                        
                        # try:
                        if(cystocyc != 0):
                            x.addchloroacetyl_toNt(count)                           #
                            x.chloroacetylate_n_cys_Cyclic_PDBwriter(linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc)
                            chloroacetylate_cyclic_tocat.append(linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.mol2")
                            x.obabel_mol2_cyc(linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.pdb",linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                            # cycmol2filenames_tocat.append(linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.mol2")
                        # except:
                        #     pass
                        os.system('rm '+linesplit[0]+"_SO3_TYR_N-C_cyc.pdb"+" "+linesplit[0]+"_SO3_TYR_Nacetyl-Cys_cyc.pdb")            # if there is any problem just comment this line and check the pdb files afterwards
                        flag = "TRUE"
                    except:
                        count += 1
                        pass
                    
                break
       
        for residue_ndx in range(0,len(x.residue_name)):                            # Add po3 to Serine
            if(x.residue_name[residue_ndx] == "SER" and x.residue_index[residue_ndx] == serinetobemod ):
                count = 0 
                flag = "FALSE"
                while(count < 4 and flag == "FALSE"):
                # for i in range(4):
                    print(count, "_SER  "+linesplit[1])  # for checking errors
                    try:
                        x.PDBreader(linesplit[1])                                   # Basicly the comments of these part are the same as above part, the only difference is this part is for adding SO3 ptm
                        x.addPO3_toSER(count,serinetobemod)                                          
                        x.PDBwriter(linesplit[0]+"_PO3_SER.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_PO3_SER.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_PO3_SER.pdb",linesplit[0]+"_PO3_SER.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_PO3_SER.mol2")
                        mol2filenames_tocat.append(linesplit[0]+"_PO3_SER.mol2")
                        
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_PO3_SER_N-C_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_PO3_SER_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_PO3_SER_N-C_cyc.pdb", linesplit[0]+"_PO3_SER_N-C_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_PO3_SER_N-C_cyc.mol2")
                        # try:
                        if(cystocyc != 0):
                            x.addchloroacetyl_toNt(count)                           #
                            x.chloroacetylate_n_cys_Cyclic_PDBwriter(linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc)
                            chloroacetylate_cyclic_tocat.append(linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.mol2")
                            x.obabel_mol2_cyc(linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.pdb",linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                            # cycmol2filenames_tocat.append(linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.mol2")
                        # except:
                        #     pass
                        os.system('rm '+linesplit[0]+"_PO3_SER_N-C_cyc.pdb"+" "+linesplit[0]+"_PO3_SER_Nacetyl-Cys_cyc.pdb")            # if there is any problem just comment this line and check the pdb files afterwards
                        flag = "TRUE"
                    except:
                        count += 1
                        pass
                    
                break 
        for residue_ndx in range(0,len(x.residue_name)):                            # # Add po3 to Threonion 
            if(x.residue_name[residue_ndx] == "THR" and x.residue_index[residue_ndx] == threoninetobemod ):
                count = 0 
                flag = "FALSE"
                while(count < 4 and flag == "FALSE"):
                # for i in range(4):
                    print(count,"_THR  "+linesplit[1])  # for checking errors
                    try:
                        x.PDBreader(linesplit[1])                                   # Basicly the comments of these part are the same as above part, the only difference is this part is for adding SO3 ptm
                        x.addPO3_toTHR(count,threoninetobemod)                                          
                        x.PDBwriter(linesplit[0]+"_PO3_THR.pdb")
                        pdbfilenames_tocat.append(linesplit[0]+"_PO3_THR.pdb")
                        x.obabel_mol2_em(linesplit[0]+"_PO3_THR.pdb",linesplit[0]+"_PO3_THR.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        mol2filenames_tocat.append(linesplit[0]+"_PO3_THR.mol2")
                        mol2filenames_tocat.append(linesplit[0]+"_PO3_THR.mol2")
                        
                        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_PO3_THR_N-C_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod)
                        cycpdbfilenames_tocat.append(linesplit[0]+"_PO3_THR_N-C_cyc.pdb")
                        x.obabel_mol2_cyc(linesplit[0]+"_PO3_THR_N-C_cyc.pdb", linesplit[0]+"_PO3_THR_N-C_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                        cycmol2filenames_tocat.append(linesplit[0]+"_PO3_THR_N-C_cyc.mol2")
                        # try:
                        if(cystocyc != 0):
                            x.addchloroacetyl_toNt(count)                           #
                            x.chloroacetylate_n_cys_Cyclic_PDBwriter(linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc)
                            chloroacetylate_cyclic_tocat.append(linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.mol2")
                            x.obabel_mol2_cyc(linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.pdb",linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"PTM")
                            # cycmol2filenames_tocat.append(linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.mol2")
                        # except:
                        #     pass
                        os.system('rm '+linesplit[0]+"_PO3_THR_N-C_cyc.pdb"+" "+linesplit[0]+"_PO3_THR_Nacetyl-Cys_cyc.pdb")            # if there is any problem just comment this line and check the pdb files afterwards
                        flag = "TRUE"
                    except:
                        count += 1
                        pass
                        
                break 
                
        x.PDBreader(linesplit[1])                                                   # re-creat the object x with all charactors un-modified
        x.PDBwriter(linesplit[0]+".pdb")                                            # creat the natural pdb file
        pdbfilenames_tocat.append(linesplit[0]+".pdb")                              # add file name to list which will be used to create the pdb library
        x.obabel_mol2_em(linesplit[0]+".pdb",linesplit[0]+".mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"NOTPTM") # Convert pdb file to mol2 file; "NOTPTM" tells the obabel_mol2_em method not to remove any bonds
        mol2filenames_tocat.append(linesplit[0]+".mol2")                            # add mol2 file name to list which will be used to create the mol2 library
        
        x.n_c_Cyclic_PDBwriter(linesplit[0]+"_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod)             # add all side chains' CONECT information to pdb file 
        cycpdbfilenames_tocat.append(linesplit[0]+"_cyc.pdb")                       # add file name to list which will be used to create the cyclic pdb library
        x.obabel_mol2_cyc(linesplit[0]+"_cyc.pdb", linesplit[0]+"_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"NOTPTM") # "NOTPTM" for telling method don't need to find the atoms from PTM, only check the wrong bond between normal  side chains.
        cycmol2filenames_tocat.append(linesplit[0]+"_cyc.mol2")                     # add mol2 file name to list which will be used to create the cyclic mol2 library
        # try:
        if(cystocyc != 0):
            count = 0
            flag = "FALSE"
            while(count < 9 and flag == "FALSE"):
                try:
                    x.addchloroacetyl_toNt(count)                           #
                    x.chloroacetylate_n_cys_Cyclic_PDBwriter(linesplit[0]+"_Nacetyl_cys_cyc.pdb",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc)
                    x.obabel_mol2_cyc(linesplit[0]+"_Nacetyl_cys_cyc.pdb",linesplit[0]+"_Nacetyl_cys_cyc.mol2",tyrosynetobemod,serinetobemod,threoninetobemod,cystocyc,"NOPTM")
                    chloroacetylate_cyclic_tocat.append(linesplit[0]+"_Nacetyl_cys_cyc.mol2")
                    # cycmol2filenames_tocat.append(linesplit[0]+"_Nacetyl_cys_cyc.mol2")
                    flag = "TRUE"
                except:
                    count += 1
                    pass
        # except:
        #     pass
        os.system('rm '+linesplit[0]+"_cyc.pdb"+" "+linesplit[0]+"_Nacetyl_cys_cyc.pdb")                                    # if there is any problem just comment this line and check the pdb files afterwards
        
pdbfilenames_tocat = list(set(pdbfilenames_tocat))                                  # remove duplicate values
mol2filenames_tocat = list(set(mol2filenames_tocat))                                # remove duplicate values
cycmol2filenames_tocat = list(set(cycmol2filenames_tocat))                          # remove duplicate values
chloroacetylate_cyclic_tocat = list(set(chloroacetylate_cyclic_tocat))

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
        
with open('Nacetyl_Cys_cyclic_lib_done.mol2',"w") as Nacetyl_Cys_cyclic_lib_done:                       
    count=0
    for nacetyl_cys_cycname in chloroacetylate_cyclic_tocat:
        with open(nacetyl_cys_cycname,"r") as mol2namefile:
            Lines = mol2namefile.readlines()
            for line in Lines:
                Nacetyl_Cys_cyclic_lib_done.write(line)
        os.system('rm '+nacetyl_cys_cycname)
        # myoutput = open("ouput.dat",'w+')
        # p = subprocess.Popen(["/home/dozeduck/test/scrip_test/github/ammvitor-Peppr_0307/smina.static",
        #                      "-l","cyclic_lib_done.mol2",
        #                      "--autobox_ligand","center.mol2",
        #                      "-r","receptor.pdb",
        #                      "--autobox_add","6",
        #                      "-o","results.sdf"],stdout=myoutput).communicate()
  
