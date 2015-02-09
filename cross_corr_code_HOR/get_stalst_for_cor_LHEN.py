import os
import string
import sys
#write out the station.lst for the stations with existing SACs only.
#require both  E and N exist
ftfile=sys.argv[1]
errMsg = sys.argv[2]
foutEN=open(sys.argv[3],"w")
Esta=[]
Nsta=[]
ferr=open(errMsg,"a")
for line in open(ftfile,"r"):
	l=line.rstrip().split(".")
	ll=l[0].split("_")
	name=ll[1]+"."+l[1]+"\0"  # modified on Jan 17, 2013
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
for i in range(len(Esta)):
	if Esta[i] in Nsta:
		foutEN.write("%s\n"%(Esta[i]))
foutEN.close()
