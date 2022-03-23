from cmath import sqrt
from math import cos,sin 
from scipy.optimize import fsolve
import numpy
import math
import os
from sympy import *
from Bio.PDB import PDBParser, PDBIO, Chain, Residue

class PDBfile:

    atomic_index = []                                                               # each empty list here is to used as the charactors for object.
    atomic_name = []
    residue_name = []
    chain_name = []
    residue_index = []
    X_peratom = []
    Y_peratom = []
    Z_peratom = []
    bfactor_per_factor = []
    charge_per_factor = []
    Atomtype_per_atom = []
    
    #def __init__(self,filename,tyrosynetobemod):
       
    def DistanceCalculator(self,a,b):                                               # Calculate the distance between point_a and point_b, here a and b are two list include the 3D coordinations
        dist =math.sqrt(numpy.square(a[0]-b[0]) + numpy.square(a[1]-b[1])+numpy.square(a[2]-b[2])) 
        return dist 
    
    
    def define_area(self,point1,point2,point3):                                     # define a plane, for the purpose of calculating the distance between a point and a plane
        point1 = numpy.asarray(point1)
        point2 = numpy.asarray(point2)
        point3 = numpy.asarray(point3)
        p12    = numpy.asmatrix(point2 - point1)
        p13    = numpy.asmatrix(point3 - point1)
        n      = numpy.cross(p12,p13)                                               # set the normal vector
        px     = n[0,0]
        py     = n[0,1]
        pz     = n[0,2]
        d      = -(px * point1[0] + py * point2[1] + pz * point3[2])
        return px, py, pz, d   
        
    def point_to_area_distance(self,point1,point2,point3,point4):                   # calculate the distance between a point and a plane, here, point1,2,3 are the 3 points for define the plane; where point4 is the target point
        px, py, pz, d = self.define_area(point1, point2, point3)
        mod_d = px * point4[0] + py * point4[1] + pz * point4[2] + d
        mod_area = numpy.sqrt(numpy.sum(numpy.square([px,py,pz])))
        d = abs(mod_d) / mod_area
        return d
    

    
    def obabel_mol2_em(self, filename,outputname,tyrtomod,ptm):                     # For dealing with non_cyclic.pdb
        os.system('obabel -ipdb '+filename+' -O pep.mol2 -d')                       # call obabel to convert PLDYL_PO3.pdb fille to pep.mol2
        if(ptm == "PTM"):                                                           # if this is for PTM residues then we still need to correct the wrong bonds between sidechains
           self.ptm_MOL2readerwriter("pep.mol2",filename,tyrtomod)                  # mainly used for delete mis connected sidechain atoms and replace them with empty line, finally save the edited file as modified_PLDYL_PO3.mol2
           os.system('obabel modified_'+filename+'.mol2 -O '+outputname)            # re-arrange the modified mol2 file and output as PLDYL_PO3.mol2 or PLDYL_SO3.mol2
        else:
            os.system('obabel pep.mol2 -O '+outputname)                             # if the pdb file is not added ptms then save as PLDYL.mol2 directly

        os.system('rm pep.mol2 modified_'+filename+'.mol2')

    def obabel_mol2_cyc(self, filename,outputname,tyrtomod,ptm):                    # For dealing with cyclic.pdp
        os.system('obabel -ipdb '+filename+' -O pep.mol2 -d')                       # call obabel to convert PLDYL_PO3.pdb fille to pep.mol2
        # print(ptm)
        if(ptm == "PTM"):                                                           # if this is for PTM residues then we still need to correct the wrong bonds between sidechains
            # print(ptm)
            self.ptm_MOL2readerwriter("pep.mol2",filename,tyrtomod)                 # mainly used for delete mis connected sidechain atoms and replace them with empty line, finally save the edited file as modified_PLDYL_PO3.mol2
            os.system('obabel modified_'+filename+'.mol2 -O '+outputname+' --minimize --steps 1500 --sd')   # call obabel to minimize the modified cyclic structure, mainly to build cyclic structure based on bond information
        else:
            self.not_ptm_MOL2readerwriter("pep.mol2", filename)                     # the difference between "not_ptm_MOL2readerwriter" and "ptm_MOL2readerwriter" is this one not to need to find PTM atoms
            os.system('obabel modified_'+filename+'.mol2 -O '+outputname+' --minimize --steps 1500 --sd')

        os.system('rm pep.mol2 modified_'+filename+".mol2")
        

     
    def find_p_s_atom_index(self,tyrtomod):                                         # this is used in "ptm_MOL2readerwriter" and "n_c_Cyclic_PDBwriter" methods, to find the "s" atom in SO3 and "p" atom in PO3 
        for i in range(0,len(self.atomic_index)):
            if(self.atomic_name[i] == "P" or self.atomic_name[i] == "S" and self.residue_name[i] == "TYR" and self.residue_index[i] == tyrtomod):
                p_s_atom_index = self.atomic_index[i]
        return int(p_s_atom_index)

    def find_n_c_atom_index(self):                                                  # this is used in "ptm_MOL2readerwriter" and "not_ptm_MOL2readerwriter" methods, to find every "N" atom and "C" atom
        n_c_indexs = []
        for i in range(0,len(self.atomic_index)):
            if(self.atomic_name[i] == "N" or self.atomic_name[i] == "C"):
                n_c_indexs.append(self.atomic_index[i])
        return n_c_indexs
            
    def ptm_MOL2readerwriter(self,filename,outputname,tyrtomod):                    # For deleting wrongly linked bond, "filename" = pep.mol2 ; outputname = TLDYRL_SO3.pdb
        bond_index = []                                                             # list for collecting bond indexs
        bond_a1 = []                                                                # list for collecting the first atom in the bond line
        bond_a2 = []                                                                # list for collecting the second atom in the bond line
        bond_type = []                                                              # list for collecting the bond type information
        not_bond_lines = []
        bond_lines = []
        p_s = self.find_p_s_atom_index(tyrtomod)                                    # call find_p_s_atom_index method to grep the atom index of p in PO3 or s in SO3
        f = open(filename, "r")                                                     # f is the content wiithin filename(pep.mol2), operator is reading
        count_line=0                                                                # count_line is used for count the iteration times
        count_sign=0                                                                # count_sign is used for counting the times while "@" showed in pep.mol2
        for line in f:                                                              # here the f is a list which includes all lines within in pep.mol2
            count_line += 1                                                         # line number
            not_bond_lines.append(line)
            if(line.startswith('@')):                                               # "startwith" is a method to determine whether this line is started with the specified sign
                count_sign += 1                                                     # each time find the line is started with "@", count_sign + 1
                if(count_sign == 3):                                                # while count_sign == 3, which means current line is @< bond > line, from the next line on is the bond information
                    bond_start = count_line + 1                                     # because while count_sign == 3, the current line is @< bond > line, so the next line is what we want
                    
        f.close()                                                                   # close file f
        f = open(filename, "r")                                                     # re-open pep.mol2 as f for reading
        count = 0 
        for line in f:
            count += 1
            if(count >= bond_start):                                                # count = the current line number, as count >= bond_start which means current line is bond information line
                bond_lines.append(line)                                             # append full bond informations to bond_lines list
                bond_index.append(line.split()[0])                                  # append the first column which is the bond idexs of each bond to bond_index
                bond_a1.append(int(line.split()[1]))                                # append the second column which is the first atom in the bond to bond_a1
                bond_a2.append(int(line.split()[2]))                                # append the third column which is the second atom in the bond to bond_a2
                bond_type.append(line.split()[3])                                   # append the fourth column which is the bond type in the bond to bond_type 
        list_a1 = [p_s + 1, p_s + 2, p_s + 3, p_s]                                  # creat the list for PTM sidechains, O1 = p_s + 1; O2 = p_s + 2; O3 = p_s + 3
        list_a2 = [p_s + 1, p_s + 2, p_s + 3, p_s, p_s - 2]                         # OH = p_s - 2
        n_c_index = self.find_n_c_atom_index()                                      # call method "find_n_c_atom_index" to grep all start atom "n" and tail atom "c" for each residue
        
        del_line_index = []                                                         # creat the list for wrongly generated bond information while conver pdb to mol2.
        for i in range(0,len(bond_lines)):                                          # iterat the list of bond lines
            if(bond_a1[i] in list_a1 and bond_a2[i] not in list_a2):                # if first atom in bond is P, S, O1, O2 or O3; in the mean time the 2nd atom in the bond is not P, S, O1, O2 or O3 or OH  then add this line number to list
                del_line_index.append(bond_start + i - 1)                           # 'bond_start' = number of the 1st line of bond lines，i started from 0，so bond_start + i = current line number;  current line number - 1 = current line index
            elif(bond_a1[i] not in n_c_index and abs(bond_a1[i] - bond_a2[i]) > 6 ):    # As long as the 1st atom of bond isn't c or n terminal for each residue，the absolute value between bond_a1 and bond_a2 > 6(because in TRP-W CB-CG=6，which is the biggest difference in natural residues)
                del_line_index.append(bond_start + i - 1)                           # 'bond_start' = number of the 1st line of bond lines，i started from 0，so bond_start + i = current line number;  current line number - 1 = current line index


        with open(filename) as fp_in:                                               # this time started to edit pep.mol2, firstly we open pep.mol2 as fp_in
            with open("modified_"+outputname+".mol2", 'w') as fp_out:               # create the output file as fp_out
                fp_out.writelines(line for i, line in enumerate(fp_in) if i not in del_line_index)  # while the line index in pep.mol2 = the error bond line's index, we skip this line
        for del_line_indexes in del_line_index:                                     # in the previous step, we skipped the wrong bond line, but we still need to add empty line to replace their position, otherwise obabel won't work for us
            os.system('sed -i '+"'" +str(del_line_indexes)+"G' modified_"+outputname+".mol2")  # replace the skipped lines with empty line, in order to let obabel work sucessfully

    def not_ptm_MOL2readerwriter(self,filename,outputname):                         # basicly the same as previous method "ptm_MOL2_readerwriter", except this one don't need to concern PTM sidechain atoms
        bond_index = []
        bond_a1 = []
        bond_a2 = []
        bond_type = []
        not_bond_lines = []
        bond_lines = []
        f = open(filename, "r")                                             
        count_line=0
        count_sign=0
        for line in f:                                                      
            count_line += 1
            not_bond_lines.append(line)
            if(line.startswith('@')):                                       
                count_sign += 1
                if(count_sign == 3):
                    bond_start = count_line + 1                             
                    
        f.close()
        f = open(filename, "r")
        count = 0 
        for line in f:
            count += 1
            if(count >= bond_start):                                                
                bond_lines.append(line)
                bond_index.append(line.split()[0])
                bond_a1.append(int(line.split()[1]))
                bond_a2.append(int(line.split()[2]))
                bond_type.append(line.split()[3])
        n_c_index = self.find_n_c_atom_index()                               
        del_line_index = []                                                         
        for i in range(0,len(bond_lines)):
            if(bond_a1[i] not in n_c_index and abs(bond_a1[i] - bond_a2[i]) > 6 ):
                del_line_index.append(bond_start + i - 1)                           

        with open(filename) as fp_in:
            with open("modified_"+outputname+".mol2", 'w') as fp_out:
                fp_out.writelines(line for i, line in enumerate(fp_in) if i not in del_line_index)  
        for del_line_indexes in del_line_index:                                                     
            os.system('sed -i '+"'" +str(del_line_indexes)+"G' modified_"+outputname+".mol2")   
    
    def addSO3_toTYR(self,addition,tyrtobemod):                                     # 'addition' used to adjust the addition value of coordinations, which is "count" in __main__file
        a = 3
        refCE1=[round(14.180,a), round(5.372,a), round(37.509,a)]                   # reference coordination from crystal structure:1H8I, modified the bond length between S-OH; S-CZ; CZ-CE2; CZ-CE1 based on the output bond length of Crankpep   
        refCZ=[round(15.504,a), round(5.242,a), round(37.912,a)]
        refOH=[round(16.403,a), round(6.281,a), round(37.877,a)]
        refS=[round(17.374,a), round(6.411,a), round(36.652,a)]
        refO1=[round(18.162,a), round(5.216,a), round(36.410,a)]
        refO2=[round(16.359,a), round(6.985,a), round(35.794,a)]
        refO3=[round(18.347,a), round(7.484,a), round(37.142,a)]
        refO=[round(10.075,a), round(0.739,a), round(38.038,a)]
        refC=[round(10.871,a), round(0.787,a), round(38.988,a)]
        refN=[round(12.324,a), round(2.361,a), round(40.323,a)]
        refCE2=[round(15.970,a), round(4.014,a), round(38.340,a)]
        
        distS_OH=self.DistanceCalculator(refS,refOH)                                 # reference distance between atom S and atom OH
        distS_CZ=self.DistanceCalculator(refS,refCZ)                                 # reference distance between atom S and atom CZ
        distS_CE1=self.DistanceCalculator(refS,refCE1)                               # reference distance between atom S and atom CE1
        # distS_CE2=self.DistanceCalculator(refS,refCE2)                             # reference distance between atom S and atom CE2        
        distS_N=self.DistanceCalculator(refS,refN)
        distS_C=self.DistanceCalculator(refS,refC)
        
        distO1_OH=self.DistanceCalculator(refO1,refOH)                               # reference distance between atom O1 and atom OH
        distO1_CZ=self.DistanceCalculator(refO1,refCZ)                               # reference distance between atom O1 and atom CZ
        distO1_CE1=self.DistanceCalculator(refO1,refCE1)                             # reference distance between atom O1 and atomCE1
        # distO1_CE2=self.DistanceCalculator(refO1,refCE2)
        distO1_O=self.DistanceCalculator(refO1,refO)
        distO1_C=self.DistanceCalculator(refO1,refC)
        # distO1_N=self.DistanceCalculator(refO1,refN)
        # distO1_S=self.DistanceCalculator(refO1,refS)
        
        distO2_OH=self.DistanceCalculator(refO2,refOH)
        distO2_CZ=self.DistanceCalculator(refO2,refCZ)
        distO2_CE1=self.DistanceCalculator(refO2,refCE1)
        # distO2_CE2=self.DistanceCalculator(refO2,refCE2)
        distO2_O=self.DistanceCalculator(refO2,refO)
        distO2_C=self.DistanceCalculator(refO2,refC)
        # distO2_N=self.DistanceCalculator(refO2,refN)
        # distO2_S=self.DistanceCalculator(refO2,refS)
        
        distO3_OH=self.DistanceCalculator(refO3,refOH)
        distO3_CZ=self.DistanceCalculator(refO3,refCZ)
        distO3_CE1=self.DistanceCalculator(refO3,refCE1)
        # distO3_CE2=self.DistanceCalculator(refO3,refCE2)
        distO3_O=self.DistanceCalculator(refO3,refO)
        distO3_C=self.DistanceCalculator(refO3,refC)
        # distO3_N=self.DistanceCalculator(refO3,refN)
        # distO3_S=self.DistanceCalculator(refO3,refS)      

        CZ=[]
        OH=[]
        CE1=[]
        CE2=[]
        N=[]
        C=[]
        b=[0.001,0.002,0.003,0.004]                                                 # To adjust the addition value
        for i in range (0, len(self.X_peratom)):
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                CZ.append(self.X_peratom[i]+b[addition])
                CZ.append(self.Y_peratom[i]+b[addition])
                CZ.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "OH" and self.residue_index[i] == tyrtobemod):
                OH.append(self.X_peratom[i]+b[addition])
                OH.append(self.Y_peratom[i]+b[addition])
                OH.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE1" and self.residue_index[i] == tyrtobemod):
                CE1.append(self.X_peratom[i]+b[addition])
                CE1.append(self.Y_peratom[i]+b[addition])
                CE1.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE2" and self.residue_index[i] == tyrtobemod):
                CE2.append(self.X_peratom[i])
                CE2.append(self.Y_peratom[i])
                CE2.append(self.Z_peratom[i])


        # Calculate the coordination for atom "S"
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedS=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distS_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distS_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distS_CE1**2],
                            [x,y,z],[OH[0],OH[1],OH[2]])
        print(solvedS)
        for i in range (0, len(self.X_peratom)):                                # add S to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+1, 1+float(len(self.X_peratom)))
                self.atomic_name.insert(i+1, "S")
                self.residue_name.insert(i+1,"TYR")
                self.chain_name.insert(i+1,"A")
                self.residue_index.insert(i+1,self.residue_index[i])
                self.X_peratom.insert(i+1,float(solvedS[0]))
                self.Y_peratom.insert(i+1,float(solvedS[1]))
                self.Z_peratom.insert(i+1,float(solvedS[2]))
                self.bfactor_per_factor.insert(i+1,float(1))
                self.charge_per_factor.insert(i+1,float(1))
                self.Atomtype_per_atom.insert(i+1,"S")
        
        # Calculate the coordination for atom "O1"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO1=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO1_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO1_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO1_CE1**2],
                            [x,y,z],[solvedS[0],solvedS[1],solvedS[2]])
        print(solvedO1)
        for i in range (0, len(self.X_peratom)):                                # add O1 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+2, 2+float(len(self.X_peratom)))
                self.atomic_name.insert(i+2, "O1")
                self.residue_name.insert(i+2,"TYR")
                self.chain_name.insert(i+2,"A")
                self.residue_index.insert(i+2,self.residue_index[i])
                self.X_peratom.insert(i+2,float(solvedO1[0]))
                self.Y_peratom.insert(i+2,float(solvedO1[1]))
                self.Z_peratom.insert(i+2,float(solvedO1[2]))
                self.bfactor_per_factor.insert(i+2,float(1))
                self.charge_per_factor.insert(i+2,float(1))
                self.Atomtype_per_atom.insert(i+2,"O") 

        
        # Calculate the coordination for atom "O2"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO2=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO2_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO2_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO2_CE1**2],
                            [x,y,z],[solvedS[0],solvedS[1],solvedS[2]])
        #print(solvedO2)
        for i in range (0, len(self.X_peratom)):                                # add O2 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+3, 3+float(len(self.X_peratom)))
                self.atomic_name.insert(i+3, "O2")
                self.residue_name.insert(i+3,"TYR")                
                self.chain_name.insert(i+3,"A")
                self.residue_index.insert(i+3,self.residue_index[i])
                self.X_peratom.insert(i+3,float(solvedO2[0]))
                self.Y_peratom.insert(i+3,float(solvedO2[1]))
                self.Z_peratom.insert(i+3,float(solvedO2[2]))
                self.bfactor_per_factor.insert(i+3,float(1))
                self.charge_per_factor.insert(i+3,float(1))
                self.Atomtype_per_atom.insert(i+3,"O") 

        
        # Calculate the coordination for atom "O3"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO3=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO3_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO3_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO3_CE1**2],
                            [x,y,z],[solvedS[0],solvedS[1],solvedS[2]])
        #print(solvedO3)
        for i in range (0, len(self.X_peratom)):                                    # add O3 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+4, 3+float(len(self.X_peratom)))
                self.atomic_name.insert(i+4, "O3")
                self.residue_name.insert(i+4,"TYR")
                self.chain_name.insert(i+4,"A")
                self.residue_index.insert(i+4,self.residue_index[i])
                self.X_peratom.insert(i+4,float(solvedO3[0]))
                self.Y_peratom.insert(i+4,float(solvedO3[1]))
                self.Z_peratom.insert(i+4,float(solvedO3[2]))
                self.bfactor_per_factor.insert(i+4,float(1))
                self.charge_per_factor.insert(i+4,float(1))
                self.Atomtype_per_atom.insert(i+4,"O") 
        for i in range(0,len(self.X_peratom)):                                      # rearrange the atomic index numbers
            self.atomic_index[i] = float(i+1)
        
    def addPO3_toTYR(self,addition,tyrtobemod):                                     # 'addition' used to adjust the addition value of coordinations
        a = 3
        refCE1=[round(14.180,a), round(5.372,a), round(37.509,a)]                   # reference coordination from crystal structure:1H8I   
        refCZ=[round(15.504,a), round(5.242,a), round(37.912,a)]
        refOH=[round(16.403,a), round(6.281,a), round(37.877,a)]
        refP=[round(17.374,a), round(6.411,a), round(36.652,a)]
        refO1=[round(18.162,a), round(5.216,a), round(36.410,a)]
        refO2=[round(16.359,a), round(6.985,a), round(35.794,a)]
        refO3=[round(18.347,a), round(7.484,a), round(37.142,a)]
        refO=[round(10.075,a), round(0.739,a), round(38.038,a)]
        refC=[round(10.871,a), round(0.787,a), round(38.988,a)]
        refN=[round(12.324,a), round(2.361,a), round(40.323,a)]
        refCE2=[round(15.970,a), round(4.014,a), round(38.340,a)]
        

        
        distP_OH=self.DistanceCalculator(refP,refOH)                                 # reference distance between atom S and atom OH
        distP_CZ=self.DistanceCalculator(refP,refCZ)                                 # reference distance between atom S and atom CZ
        distP_CE1=self.DistanceCalculator(refP,refCE1)                               # reference distance between atom S and atom CE1
        # distP_CE2=self.DistanceCalculator(refP,refCE2)                             # reference distance between atom S and atom CE2        
        
        distO1_OH=self.DistanceCalculator(refO1,refOH)                               # reference distance between atom O1 and atom OH
        distO1_CZ=self.DistanceCalculator(refO1,refCZ)                               # reference distance between atom O1 and atom CZ
        distO1_CE1=self.DistanceCalculator(refO1,refCE1)                             # reference distance between atom O1 and atomCE1
        # distO1_CE2=self.DistanceCalculator(refO1,refCE2)
        # distO1_O=self.DistanceCalculator(refO1,refO)
        # distO1_N=self.DistanceCalculator(refO1,refN)
        # distO1_P=self.DistanceCalculator(refO1,refP)
        
        distO2_OH=self.DistanceCalculator(refO2,refOH)
        distO2_CZ=self.DistanceCalculator(refO2,refCZ)
        distO2_CE1=self.DistanceCalculator(refO2,refCE1)
        # distO2_CE2=self.DistanceCalculator(refO2,refCE2)
        # distO2_O=self.DistanceCalculator(refO2,refO)
        # distO2_N=self.DistanceCalculator(refO2,refN)
        # distO2_P=self.DistanceCalculator(refO2,refP)
        
        distO3_OH=self.DistanceCalculator(refO3,refOH)
        distO3_CZ=self.DistanceCalculator(refO3,refCZ)
        distO3_CE1=self.DistanceCalculator(refO3,refCE1)
        # distO3_CE2=self.DistanceCalculator(refO3,refCE2)
        # distO3_O=self.DistanceCalculator(refO3,refO)
        # distO3_N=self.DistanceCalculator(refO3,refN)
        # distO3_P=self.DistanceCalculator(refO3,refP)      

        CZ=[]
        OH=[]
        CE1=[]
        CE2=[]
        b=[0.001,0.002,0.003,0.004]
        for i in range (0, len(self.X_peratom)):
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                CZ.append(self.X_peratom[i]+b[addition])
                CZ.append(self.Y_peratom[i]+b[addition])
                CZ.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "OH" and self.residue_index[i] == tyrtobemod):
                OH.append(self.X_peratom[i]+b[addition])
                OH.append(self.Y_peratom[i]+b[addition])
                OH.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE1" and self.residue_index[i] == tyrtobemod):
                CE1.append(self.X_peratom[i]+b[addition])
                CE1.append(self.Y_peratom[i]+b[addition])
                CE1.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE2" and self.residue_index[i] == tyrtobemod):
                CE2.append(self.X_peratom[i])
                CE2.append(self.Y_peratom[i])
                CE2.append(self.Z_peratom[i])
       

         
        # Calculate the coordination for atom "P"        
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedP=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distP_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distP_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distP_CE1**2],
                            [x,y,z],[OH[0],OH[1],OH[2]])
        #print(solvedP)
        for i in range (0, len(self.X_peratom)):                                    # add S to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+1, 1+float(len(self.X_peratom)))
                self.atomic_name.insert(i+1, "P")
                self.residue_name.insert(i+1,"TYR")
                self.chain_name.insert(i+1,"A")
                self.residue_index.insert(i+1,self.residue_index[i])
                self.X_peratom.insert(i+1,float(solvedP[0]))
                self.Y_peratom.insert(i+1,float(solvedP[1]))
                self.Z_peratom.insert(i+1,float(solvedP[2]))
                self.bfactor_per_factor.insert(i+1,float(1))
                self.charge_per_factor.insert(i+1,float(1))
                self.Atomtype_per_atom.insert(i+1,"P")

        

        # Calculate the coordination for atom "O1"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO1=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO1_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO1_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO1_CE1**2],
                            [x,y,z],[solvedP[0],solvedP[1],solvedP[2]])
        #print(solvedO1)
        for i in range (0, len(self.X_peratom)):                                    # add O1 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+2, 2+float(len(self.X_peratom)))
                self.atomic_name.insert(i+2, "O1")
                self.residue_name.insert(i+2,"TYR")
                self.chain_name.insert(i+2,"A")
                self.residue_index.insert(i+2,self.residue_index[i])
                self.X_peratom.insert(i+2,float(solvedO1[0]))
                self.Y_peratom.insert(i+2,float(solvedO1[1]))
                self.Z_peratom.insert(i+2,float(solvedO1[2]))
                self.bfactor_per_factor.insert(i+2,float(1))
                self.charge_per_factor.insert(i+2,float(1))
                self.Atomtype_per_atom.insert(i+2,"O")

        

        # Calculate the coordination for atom "O2"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO2=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO2_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO2_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO2_CE1**2],
                            [x,y,z],[solvedP[0],solvedP[1],solvedP[2]])
        #print(solvedO2)
        for i in range (0, len(self.X_peratom)):                                    # add O2 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+3, 3+float(len(self.X_peratom)))
                self.atomic_name.insert(i+3, "O2")
                self.residue_name.insert(i+3,"TYR")
                self.chain_name.insert(i+3,"A")
                self.residue_index.insert(i+3,self.residue_index[i])
                self.X_peratom.insert(i+3,float(solvedO2[0]))
                self.Y_peratom.insert(i+3,float(solvedO2[1]))
                self.Z_peratom.insert(i+3,float(solvedO2[2]))
                self.bfactor_per_factor.insert(i+3,float(1))
                self.charge_per_factor.insert(i+3,float(1))
                self.Atomtype_per_atom.insert(i+3,"O") 

        

        # Calculate the coordination for atom "O3"     
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        solvedO3=nsolve([(x-OH[0])**2+(y-OH[1])**2+(z-OH[2])**2-distO3_OH**2,
                             (x-CZ[0])**2+(y-CZ[1])**2+(z-CZ[2])**2-distO3_CZ**2,
                             (x-CE1[0])**2+(y-CE1[1])**2+(z-CE1[2])**2-distO3_CE1**2],
                            [x,y,z],[solvedO2[0],solvedO2[1],solvedO2[2]])
        #print(solvedO3)
        for i in range (0, len(self.X_peratom)):                                    # add O3 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ" and self.residue_index[i] == tyrtobemod):
                self.atomic_index.insert(i+4, 4+float(len(self.X_peratom)))
                self.atomic_name.insert(i+4, "O3")
                self.residue_name.insert(i+4,"TYR")
                self.chain_name.insert(i+4,"A")
                self.residue_index.insert(i+4,self.residue_index[i])
                self.X_peratom.insert(i+4,float(solvedO3[0]))
                self.Y_peratom.insert(i+4,float(solvedO3[1]))
                self.Z_peratom.insert(i+4,float(solvedO3[2]))
                self.bfactor_per_factor.insert(i+4,float(1))
                self.charge_per_factor.insert(i+4,float(1))
                self.Atomtype_per_atom.insert(i+4,"O") 
        for i in range(0,len(self.X_peratom)):                                      # rearrange the atomic index numbers
            self.atomic_index[i] = float(i+1)

    
    
    def addO_toCYS(self):
        CB_SG_dist = 2
        SG_OG_dist = 2
        angle = 2*3.1415*60/360
        CBX =0.0
        CBY =0.0
        CBZ =0.0
        SGX =0.0
        SGY =0.0
        SGZ =0.0
        for i in range (0, len(self.X_peratom)):
            if(self.residue_name[i] == "CYS" and self.atomic_name[i] == "CB"):
                CBX = -self.X_peratom[i]
                CBY = -self.Y_peratom[i]
                CBZ = -self.Z_peratom[i]
            if(self.residue_name[i] == "CYS" and self.atomic_name[i] == "SG"):
                SGX = self.X_peratom[i]
                SGY = self.Y_peratom[i]
                SGZ = self.Z_peratom[i]
             
        alpha = (CBX - SGX)
        beta = (CBY - SGY)
        omega = (CBZ - SGZ)
        OZ = 1.2*(alpha/abs(alpha))*cos(angle)+SGZ
        OY = SGY
        OX = -1.2*(beta/abs(beta))*sin(angle)+SGX
        self.atomic_index.append(float(50))
        self.atomic_name.append("OG")
        self.residue_name.append("CYS")
        self.chain_name.append("A")
        self.residue_index.append(float(4))
        self.X_peratom.append(float(OX))
        self.Y_peratom.append(float(OY))
        self.Z_peratom.append(float(OZ))
        self.bfactor_per_factor.append(float(1))
        self.charge_per_factor.append(float(1))
        self.Atomtype_per_atom.append("O")

    

    
    def PDBreader(self,filename):                                                   # first method to be called in __main__, used for creating object and charactors.
        self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors
        self.atomic_index.clear()
        self.atomic_name.clear()
        self.residue_name.clear()
        self.chain_name.clear()
        self.residue_index.clear()
        self.X_peratom.clear()
        self.Y_peratom.clear()
        self.Z_peratom.clear()
        self.bfactor_per_factor.clear()
        self.charge_per_factor.clear()
        self.Atomtype_per_atom.clear()
        f = open(filename, "r")                                                     # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
        for line in f:                                                              # iterate each line in file "f"                     
                
                if(line.split()[0] == "ATOM" or line.split()[0] == "HETATM"):       # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
                    self.atomic_index.append(float(line.split()[1]))                # The second column is the atomic number
                    self.atomic_name.append(line.split()[2])                        # The 3rd column is the atom name C CA CD1 CD2 and so on
                    self.residue_name.append(line.split()[3])                       # Column 4 is the residue name TYR ALA etc.
                    self.chain_name.append(line.split()[4])                         # The 5th column is the name of the chain it is on
                    self.residue_index.append(float(line.split()[5]))               # The sixth column is the residue number
                    self.X_peratom.append(float(line.split()[6]))                   # Column 7 is the x-coordinate of the atom
                    self.Y_peratom.append(float(line.split()[7]))                   # The 8th column is the Y-coordinate of the atom
                    self.Z_peratom.append(float(line.split()[8]))                   # The ninth column is the Z-coordinate of the atom
                    self.bfactor_per_factor.append(float(line.split()[9]))          # The 10th column is the B-factor of the atom, which is used to judge the activity level
                    self.charge_per_factor.append(float(line.split()[10]))          # Column 11 is the charge of the residue
                    self.Atomtype_per_atom.append(line.split()[11])                 # Column 12 is the atomic type of the atom C, H, O, N, S, CL, etc.
                    #print(line)

    
    def PDBwriter(self,filename):
        f = open(filename, "w")                                                             # e.g: f = linesplit[0]+"_PO3.pdb"
        for i in range (0 ,len(self.atomic_index)):                                         # Create a loop, i is a sequence starting from 0, and the number of atoms is the length  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # Formatted output, %4s, right-aligned, the output occupies 4 columns in total. If the length is less than 4 columns, the left end will be filled with spaces. If it is greater than 4 columns, the actual length will be output as a string
                                             self.atomic_index[i],                          # %7d, right-aligned, the output occupies a total of 7 columns, if the length is less than 7 columns, the left end is filled with spaces, signed decimal certificate integer
                                             self.atomic_name[i],                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                             self.residue_name[i],                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                             self.chain_name[i],                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                             self.residue_index[i],                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                             self.X_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Y_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Z_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.bfactor_per_factor[i],                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.charge_per_factor[i],                     # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.Atomtype_per_atom[i]), file = f )         # %12s, right-aligned, the output occupies a total of 12 columns, if it is less than 12 columns, it will be filled with spaces from the left end
        print("END", file = f)
        f.close()
            

    def n_c_Cyclic_PDBwriter(self,filename,tyrtomod):
        find_c_atom_index = []                                                              # create the list for saving c-terminal atom for each residues
        find_n_atom_index = []                                                              # create the list for saving n-terminal atom for each residues
        find_ca_tom_index = []                                                              # create the list for saving "CA" atom for each residues
        for i in range(0,len(self.atomic_index)):
            if(self.atomic_name[i] == "C"):
                find_c_atom_index.append(self.atomic_index[i])
            if(self.atomic_name[i] == "N"):
                find_n_atom_index.append(self.atomic_index[i])
            if(self.atomic_name[i] == "CA"):
                find_ca_tom_index.append(self.atomic_index[i])
            if(self.residue_name[i]== "PRO" and self.residue_index[i]==1 and self.atomic_name[i] == "CD"):      # if the 1st residue is "PRO" then we need to add its special bond CONECT information individually.   
                cd_terminal_atomIndex = self.atomic_index[i]
        c_terminal_atomIndex = max(find_c_atom_index)                                       # max() to find the last residue's C-termial atom
        n_terminal_atomIndex = min(find_n_atom_index)                                       # min() to find the first residue's N-termial atom    
        ca_terminal_atomIndex= min(find_ca_tom_index)                                       # min() to find the first residues's CA atom
        # return n_terminal_atomIndex, c_terminal_atomIndex
        f = open(filename, "w")                                                             
        for i in range (0 ,len(self.atomic_index)):                                           
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     
                                             self.atomic_index[i],                          
                                             self.atomic_name[i],                           
                                             self.residue_name[i],                          
                                             self.chain_name[i],                            
                                             self.residue_index[i],                         
                                             self.X_peratom[i],                             
                                             self.Y_peratom[i],                             
                                             self.Z_peratom[i],                             
                                             self.bfactor_per_factor[i],                    
                                             self.charge_per_factor[i],                     
                                             self.Atomtype_per_atom[i]), file = f )         
        if(self.residue_name[1]=="PRO"):                                                    # If the first residue is PRO, then the 'N' in n-terminal should like this: N-C; N-CA; N-CD
            print("%-6s%5d%5d%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(c_terminal_atomIndex), int(ca_terminal_atomIndex), int(cd_terminal_atomIndex)), file = f)
        else:                                                                               # If the first residue is not PRO, then the 'N' in n-terminal
            print("%-6s%5d%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(c_terminal_atomIndex), int(ca_terminal_atomIndex)), file = f)
            print("%-6s%5d%5d%5d%5d" % ("CONECT", int(c_terminal_atomIndex), int(n_terminal_atomIndex), int(c_terminal_atomIndex)-1, int(c_terminal_atomIndex)+1), file = f)
            print("%-6s%5d%5d" % ("CONECT", int(c_terminal_atomIndex)-1, int(c_terminal_atomIndex)), file = f)
        try:                                                                                # if it's PTM added peptide we still need to add PTM sidechain atom CONECT information to guide obabel generate mol2 file
            ptm_index = self.find_p_s_atom_index(tyrtomod)                                  
            if ptm_index:                                                                   # if ptm_index is not null then we can add the CONECT informatioins below
                print("%-6s%5d%5d%5d%5d%5d" % ("CONECT", int(ptm_index), int(ptm_index)-2, int(ptm_index)+1, int(ptm_index)+2, int(ptm_index)+3), file = f)     # S OH O1 O2 O3 
                print("%-6s%5d%5d" % ("CONECT", int(ptm_index)-2, int(ptm_index)), file = f)                                                                    # OH S
            # print("END", file = f)
        except:
            pass
        f.close()
        self.pdb_residues_side_chain_conect(filename)                                       # call "pdb_residues_side_chain_conect" method to add sidechain CONECT record for each residue in the sequence
    
    def n_cysteine_Cyclic_PDBwriter(self,filename):
        find_s_atom_index = []
        find_n_atom_index = []
        for i in range(0,len(self.atomic_index)):
            if(self.residue_name[i] == "CYS" and self.atomic_name[i] == "SG"):
                find_s_atom_index.append(self.atomic_index[i])
            if(self.atomic_name[i] == "N"):
                find_n_atom_index.append(self.atomic_index[i])

        s_terminal_atomIndex = max(find_s_atom_index)
        n_terminal_atomIndex = min(find_n_atom_index)
        # return n_terminal_atomIndex, s_terminal_atomIndex
        f = open(filename, "w")                                                             
        for i in range (0 ,len(self.atomic_index)):                                           
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     
                                             self.atomic_index[i],                          
                                             self.atomic_name[i],                           
                                             self.residue_name[i],                          
                                             self.chain_name[i],                            
                                             self.residue_index[i],                         
                                             self.X_peratom[i],                             
                                             self.Y_peratom[i],                             
                                             self.Z_peratom[i],                             
                                             self.bfactor_per_factor[i],                    
                                             self.charge_per_factor[i],                     
                                             self.Atomtype_per_atom[i]), file = f )         
        # n, s = self.n_cysteine_Cyclic()
        print("%-6s%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(s_terminal_atomIndex)), file = f)
        print("END", file = f)
        f.close()
        
    

        
    def pdb_residues_side_chain_conect(self, filename):                                     # this method used in "n_c_Cyclic_PDBwriter" mainly for add all residues sidechain "CONECT" information to PDB files
        aa = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
        n_index = []                                                                        # creat the list of the n terminal atom for each residues
        pep_seq = []                                                                        # creat the list of the residues of this sequence    
        for i in range(len(self.atomic_index)):
            if(self.atomic_name[i] == "N"):
                n_index.append(self.atomic_index[i])
                pep_seq.append(self.residue_name[i])
        self.side_chain_conect(pep_seq,n_index,filename)
                
    def side_chain_conect(self,pepseq,nindex,filename):                                    # this method used in "pdb_residues_side_chain_conect"; pepseq is the sequence of peptide
        f = open(filename, "a")                                                            # 'a' means append new content from the bottom of the file
        for i in range(len(pepseq)):
            if(pepseq[i] == "ALA"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
            elif(pepseq[i] == "CYS"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+5), file = f) # CB-CG
            elif(pepseq[i] == "ASP"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # OD1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+7), file = f) # OD2-CG
            elif(pepseq[i] == "GLU"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+8), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+8), file = f) # CD-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD-OE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD-OE2
            elif(pepseq[i] == "PHE"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+9), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+9), file = f) # CD1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD1-CE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+9), file = f) # CD2-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+8), file = f) # CD2-CE2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+7, int(nindex[i])+10), file = f) # CE2-CZ
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+8, int(nindex[i])+10), file = f) # CE2-CZ
            # elif(pepseq[i] == "GLY"):
            elif(pepseq[i] == "HIS"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+9), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+9), file = f) # CD2-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+8), file = f) # CD2-NE2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+9), file = f) # ND1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+7), file = f) # ND1-CE1
            elif(pepseq[i] == "ILE"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+6), file = f) # CB-CG1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD1-CG1    
            elif(pepseq[i] == "LYS"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD-CE
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+8), file = f) # CE-NZ
            elif(pepseq[i] == "LEU"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+7), file = f) # CD2-CG
            elif(pepseq[i] == "MET"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # SD-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # SD-CE
            elif(pepseq[i] == "ASN"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # ND2-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+7), file = f) # OD1-CG
            elif(pepseq[i] == "PRO"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i]), int(nindex[i])+1), file = f) # N-CA
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+6), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD-CG                
            elif(pepseq[i] == "GLN"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+8), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+8), file = f) # CD-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD-OE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD-NE2
            elif(pepseq[i] == "ARG"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+7), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+6), file = f) # CD-NE
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+10), file = f) # NE-CZ
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+8, int(nindex[i])+10), file = f) # NH1-CZ
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+9, int(nindex[i])+10), file = f) # NH2-CZ
            elif(pepseq[i] == "SER"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+5), file = f) # CB-OG
            elif(pepseq[i] == "THR"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+5), file = f) # CB-CG2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+6), file = f) # CB-OG1
            elif(pepseq[i] == "VAL"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+5), file = f) # CB-CG1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+6), file = f) # CB-CG2
            elif(pepseq[i] == "TRP"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+10), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+10), file = f) # CD1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+9), file = f) # CD1-NE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+10), file = f) # CD2-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+7), file = f) # CD2-CE2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+8), file = f) # CD2-CE3
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+7, int(nindex[i])+9), file = f) # CE2-NE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+7, int(nindex[i])+12), file = f) # CE2-CZ2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+8, int(nindex[i])+13), file = f) # CE3-CZ3
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+11, int(nindex[i])+12), file = f) # CH2-CZ2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+11, int(nindex[i])+13), file = f) # CH2-CZ3
            elif(pepseq[i] == "TYR"):
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+1, int(nindex[i])+4), file = f) # CA-CB
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+4, int(nindex[i])+9), file = f) # CB-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+9), file = f) # CD1-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+5, int(nindex[i])+7), file = f) # CD1-CE1
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+9), file = f) # CD2-CG
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+6, int(nindex[i])+8), file = f) # CD2-CE2
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+7, int(nindex[i])+11), file = f) # CE2-CZ
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+8, int(nindex[i])+11), file = f) # CE2-CZ
                print("%-6s%5d%5d" % ("CONECT", int(nindex[i])+10, int(nindex[i])+11), file = f) # OH-CZ
                
        print("END", file = f)
        f.close()        
                
            
                       
            
        
