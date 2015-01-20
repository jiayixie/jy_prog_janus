# this is used to smooth the angle 

import sys
import os
if ( len(sys.argv)!=4):
	print "Input:\n1] the grd file, lon lat tan(theta)\n2] the output file name (lon lat theta)\n3] the period of this psi (e.g., 1psi-360; 2psi-180; 4psi-90)"
	sys.exit()

from numpy import *

tangrd = sys.argv[1]
outxyz = open(sys.argv[2],"w")
T = float(sys.argv[3])
os.system("grd2xyz %s > temp.xyz"%(tangrd))

for line in open("temp.xyz"):
	l=line.rstrip().split()
	azi=arctan(float(l[2]))*180/pi
	while(azi<0):
		azi=azi+T
	while(azi>T):
		azi=azi-T
	outxyz.write("%s %s %f\n"%(l[0],l[1],azi))

os.system("rm -f temp.xyz")
outxyz.close()

