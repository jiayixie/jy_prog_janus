import os
import string
import sys
#write out the station.lst for the stations with existing SACs only.

ftfile=sys.argv[1]
fstaAll=sys.argv[2]
errMsg = sys.argv[3]
foutE=open(sys.argv[4],"w")
foutN=open(sys.argv[5],"w")
Esta=[]
Nsta=[]
ferr=open(errMsg,"a")
for line in open(ftfile,"r"):
	l=line.rstrip().split(".")
	ll=l[0].split("_")
	name=ll[1]+"."+l[1]
	cmp=l[2]
	if cmp == "LHE":
	   	if name not in Esta:
		   Esta.append(name)
	elif cmp == "LHN":
		if name not in Nsta:
		   Nsta.append(name)
	else:
		ferr.write(line)
ferr.close()
print "======ALL stations #=",len(Esta),len(Nsta)
for line in open(fstaAll,"r"):
	l=line.rstrip().split()
	name=l[0]+"."+l[1]
	if name in Esta:
		foutE.write(line)
	if name in Nsta:
		foutN.write(line)
foutE.close()
foutN.close()
