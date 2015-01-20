# this is used to generate a contour for my inregular map. 
# then I can use this contour to grdclip
# the v2 comes from /media/JAPAN/jiayi/WC_cv/do_ani/NN_RL_iso1lay_joint_v4.lst which have a moon shape (?? not sure)
# the v1 comes from ?? which have a full shape without NW corner
# the v3 comes from /media/JAPAN/jiayi/WC_cv/do_ani/NN_RL_iso1lay_joint_Jun26_v1.lst which have a full shape with NW corner
 
import os
import sys
import string

lonmin=90.;lonmax=110.;
latmin=21.;latmax=40.;
step=0.5
# NN_RL_contour.lst the outer contour
# NN_RL_contour_hole.lst the inner contour
# NN_RL_contour_v2.lst the modified contour


file_out = "./NN_RL_contour_v3.lst"
#file_out2="NN_RL_contour_hole.lst"
#file_exist="/media/JAPAN/jiayi/WC_cv/do_ani/NN_RL_fill_hole.lst"
#file_exist="/media/JAPAN/jiayi/WC_cv/do_ani/NN_RL_iso1lay_joint_v4.lst"
file_exist = "/media/JAPAN/jiayi/WC_cv/do_ani/NN_RL_iso1lay_joint_Jun26_v1.lst"
fout=open(file_out,"w")
#fout2=open(file_out2,"w")
point=[];
for line in open(file_exist):
	l=line.rstrip().split()
	lon=float(l[1]);lat=float(l[2])
	id = int(l[4])
	if ((lon-lonmin)*(lon-lonmax)>0 ) or (lat-latmin)*(lat-latmax)>0:
		print "point outside region",lon,lat
#	if (lon>95 and lon<101 and lat>29 and lat<32 and id==0):
#		print lon,lat,id
#		continue
	point.append("%.1f_%.1f"%(lon,lat))
print len(point)
"""
Nlon2=int((101.-96.)/step+1)
Nlat2=int((35.-31.)/step+1)
n=0
for i in range(Nlon2):
	lon2=96.+step*i
	for j in range(Nlat2):
		lat2=31.+step*j
		p="%.1f_%.1f"%(lon2,lat2)
		if p not in point:
			fout2.write("%.1f %.1f 1\n"%(lon2,lat2))
			n=n+1
fout2.close()
print n, file_out2
"""
Nlon=int((lonmax-lonmin)/step+1)
Nlat=int((latmax-latmin)/step+1)
n=0
for i in range(Nlon):
	lon=lonmin+step*i
	for j in range(Nlat):
		lat=latmin+step*j
		p="%.1f_%.1f"%(lon,lat)
		if p in point:
			n=n+1
			fout.write("%.1f %.1f 1\n"%(lon,lat))
		else:
			fout.write("%.1f %.1f 0\n"%(lon,lat))	
#		if (lon>95 and lon<102 and lat>30 and lat<35 and p not in point):
#			fout2.write("%.1f %.1f 1\n"%(lon,lat))
#		else:
#			fout2.write("%.1f %.1f 0\n"%(lon,lat))
				

fout.close()
#fout2.close()
print n,len(point)
