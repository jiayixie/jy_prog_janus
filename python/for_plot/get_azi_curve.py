#!/urs/bin/python
#output: 1]period 2]theta 3]untheta 4]amp_of_psi 5]un_amp 6]v0 7]un_v0 8]chi2
import sys
import string

perlst =[]
perlst=[8,10,12,14,16,20,22,24,26,28,30, 32, 35 ,38,40]
"""
for i in range(8,31,2):
	perlst.append(i)
"""
n1=len(perlst)
perlst.append(25)
for i in range ( 30,75,10):
	perlst.append(i)
n3=len(perlst)
#stlolst=[101.50,106.00,116.00,117.00]
#stlalst=[27.00,30.00,31.00,25.50]
#stlolst = [114,113,110,108,105,103,102 ]
#stlalst = [ 33,30,24,25,26,29,32]
stlolst=[]
stlalst=[]
for line in open("points.txt"):
	l=line.rstrip().split()
	stlolst.append(float(l[0]))
	stlalst.append(float(l[1]))

val1={};val2={};valc={};valcun={};valv0={};valv0un={};valchi={}
Npoint = len(stlolst)
for i in range(Npoint):
	val1[i]=[1000]*n3
	val2[i]=[1000]*n3
	valc[i]=[1000]*n3
	valcun[i]=[1000]*n3
	valv0[i]=[1000]*n3
	valv0un[i]=[1000]*n3
	valchi[i]=[1000]*n3
#
for i in range(n3):
    per = perlst[i]
    if i < n1 :
 	   #fin =  "/home/jiayi/work/SC/from_lq/psi/fit_out_%d_v1_lon_lat_V0_c_psi_un"%(per)
	   fin = "/home/jiayi/work/SC/try_1psi2psi/TXTlq_aug17/fit_out_%d_v1_lon_lat_V0_c_psi_un"%(per)
	   finchi = "/home/jiayi/work/SC/try_1psi2psi/TXTlq_aug17/fit_out_%d_v1_lon_lat_chi2"%(per)
    else :
           #fin = "/home/jiayi/work/SC/try_1psi2psi/packAug14/TXT8.5.3_25p_0_4e3_2psi/fit_out_%d_v1_lon_lat_V0_c_psi_un"%(per)
	   fin = "/home/jiayi/work/SC/try_1psi2psi/TXT1_0_1e4_12psi_varT1_390/fit_out_%d_v1_lon_lat_V0_c_psi_un"%(per)
	   finchi = "/home/jiayi/work/SC/try_1psi2psi/TXT1_0_1e4_12psi_varT1_390/fit_out_%d_v1_lon_lat_chi2"%(per)
    for line in open(fin,"r"):
	l=line.rstrip().split()
	lon = float(l[0]);lat=float(l[1])
	theta = 45.0-float(l[4])/2.0; untheta = float(l[7])		#calculate the fast direction(CCW to East)
	c2 = float(l[3]);unc2=float(l[6]) #amplitude of psi signal
	v0 = float(l[2]);unv0=float(l[5]) # velocity of the point
	for np in range(Npoint):
		if (stlolst[np]-lon)**2<0.001 and (stlalst[np]-lat)**2<0.001:
			val1[np][i]=theta;val2[np][i]=untheta
			valc[np][i]=c2;valcun[np][i]=unc2;
			valv0[np][i]=v0;valv0un[np][i]=unv0;
    for line in open(finchi,"r"):
	l=line.rstrip().split()
	lon = float(l[0]);lat=float(l[1]);
	for np in range(Npoint):
                if (stlolst[np]-lon)**2<0.001 and (stlalst[np]-lat)**2<0.001:
			valchi[np][i]=float(l[2])
for i in range(Npoint):
	fout = open("sta%d_theta_un.txt"%(i),"w")
	var = val1[i][0]
	fout.write("%d %f %f %f %f %f %f %f\n"%(perlst[0],val1[i][0],val2[i][0],valc[i][0],valcun[i][0],valv0[i][0],valv0un[i][0],valchi[i][0]))
	range1=var-90;range2=var+90 #the theta is within [-90,90], it can also ba changed into [v1-90,v1+90]
	for j in range(1,n3):
		if val1[i][j]>range2:
			val1[i][j]=val1[i][j]-180
		elif val1[i][j]<range1:
			val1[i][j]=val1[i][j]+180
		fout.write("%d %f %f %f %f %f %f %f\n"%(perlst[j],val1[i][j],val2[i][j],valc[i][j],valcun[i][j],valv0[i][j],valv0un[i][j],valchi[i][j]))
		range1=val1[i][j]-90;range2=val1[i][j]+90
	fout.close()




