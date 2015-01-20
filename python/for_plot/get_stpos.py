#!/usr/bin/python
import sys

fnm = sys.argv[1]
stnm = [];stlo=[];stla=[]
for line in open("/home/jiayi/work/SC/Info/station.lst","r"):
	l=line.rstrip().split()
	stnm.append(l[0])
	stlo.append(l[1])
	stla.append(l[2])
if fnm in stnm:
	id = stnm.index(fnm)
	print stlo[id],stla[id]



