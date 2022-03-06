from pyrsistent import v
from PDB_parser_sulphotyrosine_2 import PDBfile
import os
import subprocess



tyrosynetobemod = 4

with open('outlist.dat',"r") as list_pdb:
    Lines = list_pdb.readlines()
    count = 0
    # Strips the newline character
    pdbfilenames_tocat =[]
    mol2filenames_tocat =[]
    for line in Lines:
        
        linesplit=line.split()
        x=PDBfile()
        x.PDBreader(linesplit[1])
        
        for residue_ndx in range(0,len(x.residue_name)):
            if(x.residue_name[residue_ndx] == "TYR" and x.residue_index[residue_ndx] == tyrosynetobemod ):
                x.addPO3_toTYR()
                x.PDBwriter(linesplit[0]+"_PO3.pdb")
                pdbfilenames_tocat.append(linesplit[0]+"_PO3.pdb")
                x.obabel_mol2_em(linesplit[0]+"_PO3.pdb",linesplit[0]+"_PO3.mol2")
                mol2filenames_tocat.append(linesplit[0]+"_PO3.mol2")
                x.PDBreader(linesplit[1])
                x.addSO3_toTYR()
                x.PDBwriter(linesplit[0]+"_SO3.pdb")
                pdbfilenames_tocat.append(linesplit[0]+"_SO3.pdb")
                x.obabel_mol2_em(linesplit[0]+"_SO3.pdb",linesplit[0]+"_SO3.mol2")
                mol2filenames_tocat.append(linesplit[0]+"_SO3.mol2")

                break
        x.PDBwriter(linesplit[0]+".pdb")
        pdbfilenames_tocat.append(linesplit[0]+".pdb")
        x.obabel_mol2_em(linesplit[0]+".pdb",linesplit[0]+".mol2")
        mol2filenames_tocat.append(linesplit[0]+".mol2")


with open('pdb_lib_done.pdb',"w") as pdb_lib_done:
    count = 0
    # Strips the newline character
    for pdbname in pdbfilenames_tocat:
        count = 1
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
        myoutput = open("ouput.dat",'w+')
        p = subprocess.Popen(["/home/jvscunha/Downloads/smina.static",
                             "-l","mol2_lib_done.mol2",
                             "--autobox_ligand","center.mol2",
                             "-r","receptor.pdb",
                             "--autobox_add","6",
                             "-o","results.sdf"],stdout=myoutput).communicate()

