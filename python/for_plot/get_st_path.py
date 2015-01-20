# this is for the  C_plot_FTAN_AMP_v2
# output the path location for a given _s SAC name

import sys
import string

# input the SAC name and station list (nt.stnm lon lat stnm)
if (len(sys.argv)!=5):
	print "Usage: input \n1]SAC_name \n2]station_list \n3]format of name (1: ntnm.stnm  2:stnm) \n4]out_file_name"
	print len(sys.argv)
	sys.exit()
#
sacnm=sys.argv[1]
l=sacnm.split("/")[1].split("_")
stnm1=l[1]
ll=l[2].split(".SAC")
stnm2=ll[0]

if sys.argv[3]=="1":
	key=0
elif sys.argv[3]=="2":
	key=3
else:
	print "Wrong 3rd input!"
	sys.exit()

#
flag1=-2;flag2=-2;
for line in open(sys.argv[2]):
	l=line.rstrip().split()
	lon=float(l[1])
	lat=float(l[2])
	stnm=l[key]
	#print stnm,lon,lat
	if stnm==stnm1:
		loc1="%.1f %.1f"%(lon,lat)
		flag1=1
	elif stnm == stnm2:
		loc2="%.1f %.1f"%(lon,lat)
		flag2=1
	if (flag1+flag2>0):
		break;
#
if(flag1+flag2<0):
	print "######################\n###########\n##########either %s or %s is not in the station list! Incompelete station list??"%(stnm1,stnm2)
	sys.exit()
#
fout=open(sys.argv[4],"w")
fout.write("%s\n%s"%(loc1,loc2))
fout.close()


