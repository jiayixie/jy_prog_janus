#!/usr/bin/python
import sys

f1=sys.argv[1]
f2=sys.argv[2]

fout=open("tomo.xyz34","w")
lon=[];lat=[];amp=[];azi=[];
for line in open(f1,"r"):
	l = line.rstrip().split()
	lon.append(l[0])
	lat.append(l[1])
	amp.append(l[2])
i=0;
for line in open(f2,"r"):
	l = line.rstrip().split()
	flon=l[0];flat=l[1];fazi=l[2]
	if flon == lon[i] and flat==lat[i] and fazi != "NaN" :
		fout.write("%s %s %s %s\n"%(flon,flat,amp[i],fazi))
#	else:
	#	print "%s %s not in order! or NaN"%(flon,flat)
	
	i=i+1	
fout.close()
