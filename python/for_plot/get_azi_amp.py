#this is used to get the lon lat azi amp file used as C_plot_2psi_topo2 input.
#previous algorithm is not good, since angel cannot be averaged, averaging gives inaccurate result
#here, we just select some points out of origional input file. no average at all
#for 2psi->2psi

import sys
import string
import os
if(len(sys.argv)!=4):
	print "Usage: python xx.py 1]infile 2]outfile 3]step"
	sys.exit()
N=500
file = sys.argv[1]
fout = open(sys.argv[2],"w")
step = int(sys.argv[3])
lonlst=[]
latlst=[]
ang = []
amp=[]
for i in range(N):
	ang.append([-999]*N)
	amp.append([-999]*N)

for line in open(file):
	l=line.rstrip().split()
	lon = l[0]
 	lat = l[1]
	if lon not in lonlst:
		lonlst.append(lon)
	if lat not in latlst:
		latlst.append(lat)
lonlst.sort()
latlst.sort()	
for line in open(file):
	l=line.rstrip().split()
	lon=l[0];lat=l[1]
	ilon = lonlst.index(lon)
	ilat = latlst.index(lat)
	ang[ilon][ilat]=45.0-float(l[4])/2.0
	amp[ilon][ilat]=float(l[3])/float(l[2])*200
Nlon=len(lonlst)
Nlat=len(latlst)
for i in range(Nlon):
	if (i%step!=0):
		continue
	for j in range(Nlat):
		if(j%step !=0):
			continue
		if (ang[i][j]<-990 or amp[i][j]<-990):
			continue
		fout.write("%8s %8s %8f %8f\n"%(lonlst[i],latlst[j],amp[i][j],ang[i][j]))

fout.close()
