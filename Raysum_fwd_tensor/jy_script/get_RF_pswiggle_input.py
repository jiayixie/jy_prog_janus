# this reads in a sac file, change it to asc file, and write in terms of the input of pswiggle

import os
import sys

inputsac=sys.argv[1]
outputasc=sys.argv[2]
fgeom = sys.argv[3]
idtrace = int(sys.argv[4])
azishift = float(sys.argv[5])

Ntrace=0
azilst=[];NSsftlst=[];EWsftlst=[]
for line in open(fgeom):
        if ("#" in line):
                continue
        l=line.rstrip().split()
        azi=float(l[0])
        while(azi<0):
                azi=azi+360
        azilst.append(azi)
        NSsftlst.append(float(l[2]))
        EWsftlst.append(float(l[3]))
        Ntrace+=1
if (idtrace>Ntrace):
	print "The input idtrace%d is larger than the number of trace%d"%(idtrace,Ntrace)
#--sac2asc
str =  "/home/jiayi/progs/cv/sac_ascII %s temp.asc.txt"%(inputsac)
os.system(str)
azi = azilst[idtrace-1]
azi=azi+azishift
while(azi<0):
	azi+=360.
while(azi>360):
	azi-=360.

fout=open(outputasc,"w")
for line in open("temp.asc.txt","r"):
	if ("#" in line):
		continue
	l=line.rstrip().split()
	fout.write("%8s %8f %8s\n"%(l[0],azi,l[1]))
fout.close()


