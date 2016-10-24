# this is used to compute the most likely value from the posterior distribution.
# the posterior distribution is presented at very 1 km (current setting, the 1st 3 para in get_posteriaDist(0.,100.,1....) )

# the posterior distributin file:  
# every row is one accepted model
# every row contains (maxdep-mindep)/step+1 groups of parameter 
# each group represent model value at one depth, each group contains 9 values (vsv, vsh,vpv,vph,eta,theta,phi,rho,vpv/vsv)
# posteria    1.278    1.278    2.556    2.556        1        0        0   2.2055        2 

import sys
import os
import numpy as np

filenm = "/lustre/janus_scratch/jixi7887/Tibet/inv_ET_BS/inv_Apr28data_v1/inv_v1.2.1.scale3.3laytheta.freeta.AZmat.constRA_Tibet_v1_laytheta.20disc.XYdip.aug11.vpvs.3.AZmat4.5perc/9633_96.0_33.0_inv_v1.2.1.scale3.3laytheta.freeta.AZmat.constRA_Tibet_v1_laytheta.20disc.XYdip.aug11.vpvs.3.AZmat4.5perc/bin_avg/post_96.0_33.0.txt_phigp0"


ivalueLst=[0,1] # vsv-0; vsh-1
depmin=0.
depmax=100.
step=1.;

#--- compute maximum likeli value
def getMaxMode(vlst):
	bins=np.linspace(min(vlst),max(vlst),30)
	histo=np.histogram(vlst,bins)
	maxMode=histo[1][np.argmax(histo[0])]
	return maxMode
#--- for each type of value, each depth, from all posterior records --> one list
vsvlst2d=[]
for line in open(filenm):
# loop over all posterior recods
	l=line.rstrip().split()
	# loop over all depth
	vsvlst=[]
	for depth in range(100):
		igp=(depth-depmin)/step # compute group based on depth
		icolvsv=int(1+igp*9)
		icolvsh=icolvsv+1
		vsv=float(l[icolvsv])
		vsh=float(l[icolvsv+1])
		vpv=float(l[icolvsv+2])
		vph=float(l[icolvsv+3])
		eta=float(l[icolvsv+4])
		theta=float(l[icolvsv+5])
		phi=float(l[icolvsv+6])
		rho=float(l[icolvsv+7])
		vpvs=float(l[icolvsv+8])

		ani=(vsh*vsh-vsv*vsv)/(2*vsv*vsv)
		vsvlst.append(ani)	
	vsvlst2d.append(vsvlst)
#-- now each row of vsvlst2d is the values of vsv at differetn depth from one posterior record, will transpose it to reverse row and column
vsvlst2d=np.array(vsvlst2d)
vsvlst2d=np.transpose(vsvlst2d)
#-- ok, each row of vsv2d is the values of vsv at one depth, from different posterior records
for idep in range(vsvlst2d.shape[0]): # loop over each depth
	maxmode=getMaxMode(vsvlst2d[idep])
	print idep,np.mean(vsvlst2d[idep]),np.median(vsvlst2d[idep]),maxmode


