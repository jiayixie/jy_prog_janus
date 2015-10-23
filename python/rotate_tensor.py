#ths is a subroutine contains some fundamental tensor operations: readin, rotation, transfer ...

import sys
import os
from numpy import *

set_printoptions(precision=4)
set_printoptions(suppress=True)

# rotate ET is necessary
def rotate(thetaz,thetay,thetax,flagreadin,Cm):
        #--rotating angles ---
        o3=thetaz*pi/180.;  #rotating coordinate CW through o3 about Z
        o2=thetay*pi/180.; # .............................o2 about transformmed Y
        o1=thetax*pi/180;   # .............................o1 about transformmed X
        #---coordinate transformation matrix --------
        a3m=matrix([[cos(o3),   sin(o3),        0],
                   [-sin(o3),   cos(o3),        0],
                   [0,           0,             1]])

        a2m=matrix([[cos(o2),   0,              -sin(o2)],
                  [0,           1,              0 ],
                  [sin(o2),     0,              cos(o2)]])

        a1m=matrix([[1,         0,              0],                  [0,           cos(o1),        sin(o1)],
                  [0,           -sin(o1),       cos(o1)]])
        am=a3m*a2m*a1m
        #---define coefficient transformation matrix ---------
        axx=am[0,0];axy=am[0,1];axz=am[0,2];
        ayx=am[1,0];ayy=am[1,1];ayz=am[1,2];
        azx=am[2,0];azy=am[2,1];azz=am[2,2];

        M1m=square(am)
        M2m=matrix([[2*axy*axz,         2*axz*axx,              2*axx*axy],
                    [2*ayy*ayz,         2*ayz*ayx,              2*ayx*ayy],
                    [2*azy*azz,         2*azz*azx,              2*azx*azy]])
        M3m=matrix([[ayx*azx,           ayy*azy,                ayz*azz],
                    [azx*axx,           azy*axy,                azz*axz],
                    [axx*ayx,           axy*ayy,                axz*ayz]])
        M4m=matrix([[ayy*azz+ayz*azy,   ayx*azz+ayz*azx,        ayy*azx+ayx*azy],
                    [axy*azz+axz*azy,   axz*azx+axx*azz,        axx*azy+axy*azx],
                    [axy*ayz+axz*ayy,   axz*ayx+axx*ayz,        axx*ayy+axy*ayx]])

        Mm=vstack((hstack((M1m,M2m)),hstack((M3m,M4m))))
        #print "\ncoefficient transformation matrix:\n",Mm
        Cm_=dot(dot(Mm,Cm),Mm.T)
        #print "\nrotated Elastic constant matrix:\n",Cm_
        return Cm_



# read in ET
def readinET(filein, C_unit): 
#form1 the ET starts from the 4th line; the rho is in the 2ed line
#C_unit = 1 # if the unit of the input ET is GPa. C11 in the order of 100
#C_unit = 100 # if the unit of the input ET is 100GPa. C11 in the order of 1
        C=[[0 for i in range(6)] for j in range(6)]
        nl=0
        i=0
        #global rho
        for line in open(filein):
                i=i+1
                if(i==2):
                        l=line.rstrip().split()
                        rho=float(l[3])
                        continue
                elif(i<4):
                        continue #skip the title and other info
                l=line.rstrip().split()
                for nc in range(len(l)):
                        #print nl,nc,l[nc]
                        C[nl][nc]=float(l[nc])*C_unit
                nl=nl+1
        ET= matrix(C)
        return ET,rho
def readinET_tatham(filein, C_unit): 
#form2 the ET starts from the 3rd line; the rho is 3.0340 g/cm3; this is for tatham's ET
#C_unit = 1 # if the unit of the input ET is GPa. C11 in the order of 100
#C_unit = 100 # if the unit of the input ET is 100GPa. C11 in the order of 1
        C=[[0 for i in range(6)] for j in range(6)]
        nl=0
        flag=0
        #global rho
        rho=3.0340
        for line in open(filein):
                if "average" in line:
                        flag=1
                        i=0
                        continue
                if(flag==0):
                        continue
                else:
                        i=i+1
                if(i<2):
                        continue #skip the title and other info
                l=line.rstrip().split()
                for nc in range(len(l)):
                        #print nl,nc,l[nc]
                        C[nl][nc]=float(l[nc])*C_unit
                nl=nl+1
        ET= matrix(C)
        return ET


def transfer_C_C (Cs): # this is used to transfer value from Cshort (6X6) to Clong (9X9)
        N=3
        Cl=[[[[0 for ii in xrange(N)] for ij in xrange(N)] for ik in xrange(N)] for il in xrange(N)]
        for i in range(N):
          for j in range(N):
            for k in range(N):
              for l in range(N):

                        if (i==j):
                                m=i;
                        else:
                                m=9-(i+j)-3
                        if(k==l):
                                n=k;
                        else:
                                n=9-(k+l)-3
                        """
                        print i+1,j+1,k+1,l+1,"==>", m+1,n+1
                        print i,j,k,l,"==>", m,n
                        print Cs[m,n]
                        """
                        Cl[i][j][k][l]=Cs[m,n]
        #print Cl
        return Cl

def compute_Tjl(Cl,n): # compute Tik = Cijkl*nj*nl
        T=[[0 for i1 in xrange(3)] for i2 in xrange(3)]

        for i in range(3):
          for k in range(3):
            #T[i][k]
            for j in range(3):
                for l in range(3):
                        T[i][k]=T[i][k]+Cl[i][j][k][l]*n[j]*n[l]
        return matrix(T)

def transferToHex (ET): # get the Hex part of the ET
        A=3./8.*(ET[0,0]+ET[1,1])+ET[0,1]/4.+ET[5,5]/2.
        C=ET[2,2]
        F=(ET[0,2]+ET[1,2])/2.
        L=(ET[3,3]+ET[4,4])/2.
        N=(ET[0,0]+ET[1,1])/8.-ET[0,1]/4.+ET[5,5]/2.
        Cout=matrix([[A,        A-2*N,  F,      0,      0,      0],
                   [A-2*N,      A,      F,      0,      0,      0],
                   [F,          F,      C,      0,      0,      0],
                   [0,          0,      0,      L,      0,      0],
                   [0,          0,      0,      0,      L,      0],
                   [0,          0,      0,      0,      0,      N]])
        return Cout


