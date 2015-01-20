#!/usr/bin/python

import sys
from numpy import *
import os
id = int(sys.argv[2])
val=[]
for line in open(sys.argv[1],"r"):
	l=line.rstrip().split()
	val.append(float(l[0]))
mean1 = average(val)	
med1 = median(val)
std1 = std(val)
#min1=max(med1-std1,0.0)
min1 = med1-std1
max1=med1+std1
step1 = ((max1-min1)/10)
#str = "makecpt -Chot -T%d/%d/%d -I > pycpt.cpt"%(min,max,step)
if id == 1:
	str = "makecpt -Cno_green -T%f/%f/%f -I > pycpt.cpt"%(min1,max1,step1)
	os.system(str);
	sys.exit()
if id == 3:
	str = "makecpt -Cno_green -T%d/%d/%d -I > pycpt.cpt"%(min1,max1,step1)	
	os.system(str);
	sys.exit()
if id == 4:
	step1 = ((std1)/5)
#	print mean1-std1,mean1+std1,step1,mean1,std1
	str = " makecpt -Cno_green -T%f/%f/%f -I > pycpt.cpt"%(mean1-std1,mean1+std1,step1)
	os.system(str)
	sys.exit()
valnew=[]
for i in range(len(val)):
	if fabs(val[i]-med1)>std1:
		continue
	valnew.append(val[i])	
mean2=average(valnew)
med2 = median(valnew)
std2 = std(valnew)
min2=int(max(mean2-std2,0))
max2 = int(mean2+std2)
step2 = int((max2-min2)/10)
if id == 2:
	str = "makecpt -Chot -T%d/%d/%d -I > pycpt.cpt"%(min2,max2,step2)
os.system(str)



