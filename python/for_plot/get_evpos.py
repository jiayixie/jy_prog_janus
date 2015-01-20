#!/usr/bin/python
import sys

fnm = sys.argv[1]
fevt = fnm[:14]
evnm = [];evlo=[];evla=[]
for line in open("/home/jiayi/work/SC/Info/evt.txt","r"):
	l=line.rstrip().split()
	evnm.append(l[0])
	evla.append(l[1])
	evlo.append(l[2])
if fevt in evnm:
	id = evnm.index(fevt)
	print evlo[id],evla[id]



