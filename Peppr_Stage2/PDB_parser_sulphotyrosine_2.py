from cmath import sqrt
from math import cos,sin 
from scipy.optimize import fsolve
import numpy
import math
import os
from sympy import *

class PDBfile:

    atomic_index = []
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
       
    def DistanceCalculator(self,a,b):                                            # Calculate the distance between point_a and point_b, here a and b are two list include the 3D coordinations
        dist =math.sqrt(numpy.square(a[0]-b[0]) + numpy.square(a[1]-b[1])+numpy.square(a[2]-b[2])) 
        return dist 
    
    
    def define_area(self,point1,point2,point3):                                  # define a plane, for the purpose of calculating the distance between a point and a plane
        point1 = numpy.asarray(point1)
        point2 = numpy.asarray(point2)
        point3 = numpy.asarray(point3)
        p12    = numpy.asmatrix(point2 - point1)
        p13    = numpy.asmatrix(point3 - point1)
        n      = numpy.cross(p12,p13)                                            # set the normal vector
        px     = n[0,0]
        py     = n[0,1]
        pz     = n[0,2]
        d      = -(px * point1[0] + py * point2[1] + pz * point3[2])
        return px, py, pz, d   
        
    def point_to_area_distance(self,point1,point2,point3,point4):                # calculate the distance between a point and a plane, here, point1,2,3 are the 3 points for define the plane; where point4 is the target point
        px, py, pz, d = self.define_area(point1, point2, point3)
        mod_d = px * point4[0] + py * point4[1] + pz * point4[2] + d
        mod_area = numpy.sqrt(numpy.sum(numpy.square([px,py,pz])))
        d = abs(mod_d) / mod_area
        return d
    

    
    def obabel_mol2_em(self, filename,outputname):
        os.system('obabel -ipdb '+filename+' -O pep.mol2 -d')
        cyclic=false
        if(cyclic):
            os.system('obabel pep.mol2 -O '+outputname+' --minimize --steps 1500 --sd')
        else:
            os.system('obabel pep.mol2 -O '+outputname+' ')

        os.system('rm pep.mol2')

    def obabel_mol2_cyc(self, filename,outputname):
        os.system('obabel -ipdb '+filename+' -O pep.mol2 -d')
        os.system('obabel pep.mol2 -O '+outputname+' --minimize --steps 1500 --sd')
        # os.system('obabel -ipdb '+filename+' -O '+outputname+' --minimize --steps 1500 --sd ')
        # os.system('obabel pep.mol2 -O '+outputname+' --minimize --steps 1500 --sd')
        os.system('rm pep.mol2')
            


    
    def addSO3_toTYR(self,addition):                                                # 'addition' used to adjust the addition value of coordinations
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
        # distS_CE2=self.DistanceCalculator(refS,refCE2)                               # reference distance between atom S and atom CE2        
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
        b=[0.001,0.002,0.003,0.004]                                                       # To adjust the addition value
        for i in range (0, len(self.X_peratom)):
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
                CZ.append(self.X_peratom[i]+b[addition])
                CZ.append(self.Y_peratom[i]+b[addition])
                CZ.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "OH"):
                OH.append(self.X_peratom[i]+b[addition])
                OH.append(self.Y_peratom[i]+b[addition])
                OH.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE1"):
                CE1.append(self.X_peratom[i]+b[addition])
                CE1.append(self.Y_peratom[i]+b[addition])
                CE1.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE2"):
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
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
        
    def addPO3_toTYR(self,addition):                                                # 'addition' used to adjust the addition value of coordinations
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
        # distP_CE2=self.DistanceCalculator(refP,refCE2)                               # reference distance between atom S and atom CE2        
        
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
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
                CZ.append(self.X_peratom[i]+b[addition])
                CZ.append(self.Y_peratom[i]+b[addition])
                CZ.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "OH"):
                OH.append(self.X_peratom[i]+b[addition])
                OH.append(self.Y_peratom[i]+b[addition])
                OH.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE1"):
                CE1.append(self.X_peratom[i]+b[addition])
                CE1.append(self.Y_peratom[i]+b[addition])
                CE1.append(self.Z_peratom[i]+b[addition])
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CE2"):
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
        for i in range (0, len(self.X_peratom)):                                # add S to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
        for i in range (0, len(self.X_peratom)):                                # add O1 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
        for i in range (0, len(self.X_peratom)):                                # add O2 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
        for i in range (0, len(self.X_peratom)):                                # add O3 to the end of TYR
            if(self.residue_name[i] == "TYR" and self.atomic_name[i] == "CZ"):
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
        for i in range(0,len(self.X_peratom)):                                  # rearrange the atomic index numbers
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

    

    
    def PDBreader(self,filename):
        self.atomic_index.clear()
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
        f = open(filename, "r")                                             # f 的内容为打开filename，操作为读取
        for line in f:                                                      # 创建循环，line = filename中每一行                     
                
                if(line.split()[0] == "ATOM" or line.split()[0] == "HETATM"):   # 判断句，用于将每一行split然后判断该行第一列是否==ATOM或HETATM
                    self.atomic_index.append(float(line.split()[1]))            # 第2列为原子序数
                    self.atomic_name.append(line.split()[2])                    # 第3列为原子名称C CA CD1 CD2等等
                    self.residue_name.append(line.split()[3])                   # 第4列为残基名称TYR ALA 等等
                    self.chain_name.append(line.split()[4])                     # 第5列为所在链的名称
                    self.residue_index.append(float(line.split()[5]))           # 第6列为残基序数
                    self.X_peratom.append(float(line.split()[6]))               # 第7列为原子的X轴坐标
                    self.Y_peratom.append(float(line.split()[7]))               # 第8列为原子的Y轴坐标
                    self.Z_peratom.append(float(line.split()[8]))               # 第9列为原子的Z轴坐标
                    self.bfactor_per_factor.append(float(line.split()[9]))      # 第10列为该原子的B因子，用于判断活跃程度
                    self.charge_per_factor.append(float(line.split()[10]))      # 第11列为该残基的电荷
                    self.Atomtype_per_atom.append(line.split()[11])             # 第12列为该原子的原子类型C，H，O，N，S，CL，等等
                    #print(line)

    
    def PDBwriter(self,filename):
        f = open(filename, "w")                                                             # f的内容为打开filename，操作为将print的写入原文件
        for i in range (0 ,len(self.atomic_index)):                                         # 创建循环，i 为从0开始的数列，原子数量相同  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # 格式化输出,%4s,右对齐,输出共占4列,若长度小于4列，则左端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.atomic_index[i],                          # %7d,右对齐,输出共占7列，若长度小于7列，则左端以空格补齐,有符号的十进制证整数
                                             self.atomic_name[i],                           # %-4s,左对齐,输出共占4列,若长度小于4列，则右端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.residue_name[i],                          # %1s,右对齐，输出共占1列，若小于1列，则从左端以空格补齐，若大于1列，则输出实际长度,字符串
                                             self.chain_name[i],                            # %2s,右对齐，输出共占2列，若小于2列，则从左端以空格补齐，若大于2列，则输出实际长度,字符串
                                             self.residue_index[i],                         # %4d,右对齐，输出共占4列，若长度小于4列，则左端以空格补齐，有符号的十进制证整数
                                             self.X_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Y_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Z_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.bfactor_per_factor[i],                    # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.charge_per_factor[i],                     # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.Atomtype_per_atom[i]), file = f )         # %12s,右对齐，输出共占12列，若小于12列，则从左端以空格补齐
        print("END", file = f)
        f.close()
            

    def n_c_Cyclic_PDBwriter(self,filename):
        find_c_atom_index = []
        find_n_atom_index = []
        find_ca_tom_index = []
        for i in range(0,len(self.atomic_index)):
            if(self.atomic_name[i] == "C"):
                find_c_atom_index.append(self.atomic_index[i])
            if(self.atomic_name[i] == "N"):
                find_n_atom_index.append(self.atomic_index[i])
            if(self.atomic_name[i] == "CA"):
                find_ca_tom_index.append(self.atomic_index[i])
            if(self.residue_name[i]== "PRO" and self.residue_index[i]==1 and self.atomic_name[i] == "CD"):
                cd_terminal_atomIndex = self.atomic_index[i]
        c_terminal_atomIndex = max(find_c_atom_index)
        n_terminal_atomIndex = min(find_n_atom_index)
        ca_terminal_atomIndex= min(find_ca_tom_index)
        # return n_terminal_atomIndex, c_terminal_atomIndex
        f = open(filename, "w")                                                             # f的内容为打开filename，操作为将print的写入原文件
        for i in range (0 ,len(self.atomic_index)):                                         # 创建循环，i 为从0开始的数列，原子数量相同  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # 格式化输出,%4s,右对齐,输出共占4列,若长度小于4列，则左端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.atomic_index[i],                          # %7d,右对齐,输出共占7列，若长度小于7列，则左端以空格补齐,有符号的十进制证整数
                                             self.atomic_name[i],                           # %-4s,左对齐,输出共占4列,若长度小于4列，则右端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.residue_name[i],                          # %1s,右对齐，输出共占1列，若小于1列，则从左端以空格补齐，若大于1列，则输出实际长度,字符串
                                             self.chain_name[i],                            # %2s,右对齐，输出共占2列，若小于2列，则从左端以空格补齐，若大于2列，则输出实际长度,字符串
                                             self.residue_index[i],                         # %4d,右对齐，输出共占4列，若长度小于4列，则左端以空格补齐，有符号的十进制证整数
                                             self.X_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Y_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Z_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.bfactor_per_factor[i],                    # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.charge_per_factor[i],                     # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.Atomtype_per_atom[i]), file = f )         # %12s,右对齐，输出共占12列，若小于12列，则从左端以空格补齐
        if(self.residue_name[1]=="PRO"):                                                    # If the first residue is PRO, then the 'N' in n-terminal should like this: N-C; N-CA; N-CD他
            print("%-6s%5d%5d%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(c_terminal_atomIndex), int(ca_terminal_atomIndex), int(cd_terminal_atomIndex)), file = f)
        else:                                                                               # If the first residue is not PRO, then the 'N' in n-terminal
            print("%-6s%5d%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(c_terminal_atomIndex), int(ca_terminal_atomIndex)), file = f)
            print("%-6s%5d%5d%5d%5d" % ("CONECT", int(c_terminal_atomIndex), int(n_terminal_atomIndex), int(c_terminal_atomIndex)-1, int(c_terminal_atomIndex)+1), file = f)
            print("%-6s%5d%5d" % ("CONECT", int(c_terminal_atomIndex)-1, int(c_terminal_atomIndex)), file = f)
            print("END", file = f)
            f.close()
    
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
        f = open(filename, "w")                                                             # f的内容为打开filename，操作为将print的写入原文件
        for i in range (0 ,len(self.atomic_index)):                                         # 创建循环，i 为从0开始的数列，原子数量相同  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # 格式化输出,%4s,右对齐,输出共占4列,若长度小于4列，则左端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.atomic_index[i],                          # %7d,右对齐,输出共占7列，若长度小于7列，则左端以空格补齐,有符号的十进制证整数
                                             self.atomic_name[i],                           # %-4s,左对齐,输出共占4列,若长度小于4列，则右端以空格补齐，若大于4列，则输出实际长度,字符串
                                             self.residue_name[i],                          # %1s,右对齐，输出共占1列，若小于1列，则从左端以空格补齐，若大于1列，则输出实际长度,字符串
                                             self.chain_name[i],                            # %2s,右对齐，输出共占2列，若小于2列，则从左端以空格补齐，若大于2列，则输出实际长度,字符串
                                             self.residue_index[i],                         # %4d,右对齐，输出共占4列，若长度小于4列，则左端以空格补齐，有符号的十进制证整数
                                             self.X_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Y_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.Z_peratom[i],                             # %8.3f,右对齐，输出共占8列，其中有3位小数，若数值宽度小于8左端补空格,小数
                                             self.bfactor_per_factor[i],                    # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.charge_per_factor[i],                     # %6.2f,右对齐，输出共占6列，其中有2位小数，若数值宽度小于6左端补空格,小数
                                             self.Atomtype_per_atom[i]), file = f )         # %12s,右对齐，输出共占12列，若小于12列，则从左端以空格补齐
        # n, s = self.n_cysteine_Cyclic()
        print("%-6s%5d%5d" % ("CONECT", int(n_terminal_atomIndex), int(s_terminal_atomIndex)), file = f)
        print("END", file = f)
        f.close()
        
    
        # if user ask to N-C Terminal then:  
        # n, c = self.n_c_Cyclic()
        # print("%-6s%5d%5d" % ("CONECT", int(n), int(c)), file = f)

        # # if user ask to N-Cysteine then:
        # n, s = self.n_cysteine_Cyclic()
        # print("%-6s%5d%5d" % ("CONECT", int(n), int(s)), file = f)
        # print("END", file = f)
        
        # f.close()
        
