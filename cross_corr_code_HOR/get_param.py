import sys
import os
import string

inftlst = sys.argv[1]
outfilter = sys.argv[2]
outwhiten = sys.argv[3]
fin=open(inftlst,"r")
nmE=[]
nmN=[]
nE=0
for line in fin:
	l=line.rstrip().split(".")
	l2=l[0].split("_")
	name=l2[1]+"."+l[1]
	if l[2]=="LHE":
		nmE.append(name)
		nE=nE+1
	elif l[2]=="LHN":
		nmN.append(name)
	else:
		print "!!!###Wrong sacnm, no comp info"		
fin.close()

fout1=open(outfilter,"w")
fout2=open(outwhiten,"w")
for i in range(nE):
    if nmE[i] in nmN:
	fout1.write("200 150 5 4 1 1 ft_%s.LHE.SAC\n200 150 5 4 1 1 ft_%s.LHN.SAC\n"%(nmE[i],nmE[i]))
	fout2.write("200 150 5 4 1 1 ft_%s.LHE.SAC  ft_%s.LHN.SAC\n"%(nmE[i],nmE[i]))		
fout1.close()
fout2.close()
