# based on do_fourier.py
# add another 3 inversion: 
# A0; A0 + A1; A0 + A2;
# output the match files:
# cvA0.dat, cvA0_A1.dat, cvA0_A2.dat. cvA0_A1_A2.dat
# add to compute the average misfit.

import sys;
import math;
import string;
import numpy as np;
import scipy as sp;

#### prediction ##########################################

def A0pre ( inbaz, A0 ):
	return A0;

def A1pre ( inbaz, A0, A1, SIG1):
	return A0 + A1*math.sin(inbaz+SIG1);

def A1pre1 ( inbaz, A1, SIG1):
	return A1*math.sin(inbaz+SIG1); 

def A2pre ( inbaz, A0, A2, SIG2):
	return A0 + A2*math.sin(2*inbaz+SIG2);

def A2pre1 ( inbaz, A2, SIG2):
	return A2*math.sin(2*inbaz+SIG2);

def A3pre ( inbaz, A0, A1, SIG1, A2, SIG2):
	return A0 + A1*math.sin(inbaz + SIG1) + A2*math.sin(2*inbaz + SIG2);

def A3pre1 ( inbaz, A1, SIG1):
	return A1*math.sin(inbaz + SIG1);

def A3pre2 ( inbaz, A2, SIG2):
	return A2*math.sin(2*inbaz + SIG2);

def A3pre3 ( inbaz, A1, SIG1, A2, SIG2 ):
	return A1*math.sin(inbaz + SIG1) + A2*math.sin(2*inbaz + SIG2);
##########################################################

def match1 ( data1, data2 ):
	nn = min(len(data1), len(data2));
	X1 = 0.
	di = [];
	tempdata2 = [];
	for i in range(nn):
		di.append(data1[i] - data2[i]);
		tempdata2.append(math.fabs(data2[i]));
	for i in range(nn):
		X1 = X1 + (di[i] - np.mean(di))**2;
	return math.sqrt(X1/nn);

# grid search::::
def invert ( inbaz, indat, inun):
#	print len(inbaz),len(indat),len(inun);
#	sys.exit();
	n = 3;
	m = len(inbaz);
	maA0 = max(indat);
	miA0 = min(indat);
	maA1 = math.fabs(max(indat));
	miA1 = 0;
	maPHI1  = 2*math.pi;
	miPHI1 = 0;
	misfit0 = 1000;
	ii = 0;
	iii = 0;
	for i in range (100):
		tA0 = miA0 + (maA0-miA0)/100*i;
		for j in range (1):
			tA1 = miA1 + (maA1-miA1)/100*j;
			for k in range (1):
				tPHI1 = miPHI1 + (maPHI1)/100*k;
				misfit = 0;
				iii = iii+1;
				for l in range (len(inbaz)):
#					value = tA0 + tA1 * math.cos(inbaz[l]+tPHI1);
					value = tA0;
					misfit = misfit + (value-indat[l])**2/(inun[l]**2);
				if misfit0 > math.sqrt(misfit):
					ii = ii+1;
					misfit0 = math.sqrt(misfit);
					outA0 = tA0;
					outA1 = tA1;
					outPHI1 = tPHI1;
					if (1):
						print misfit0,ii,tA0,tA1,tPHI1,iii;
	A0 = outA0;
	A1 = outA1;
	PHI1 = outPHI1;
	return (A0,A1,PHI1);


def invert_1 ( inbaz, indat, inun):
	m = len(inbaz); # data space;
	n = 5; # model space;
	G = [];
	A0 = 0;
	A1 = 0;
	A2 = 0;
	fi1 = 0;
	fi2 = 0;
#	print m,len(inun);
	U = np.zeros((m,m));
	for i in range (m):
		U[i][i] = 1/inun[i];
	for i in range (m):
		tG = [];
		tbaz = inbaz[i]/180*(math.pi);
		tG.append(1);
		tG.append(math.sin(tbaz));
		tG.append(math.cos(tbaz));
		tG.append(math.sin(tbaz*2));
		tG.append(math.cos(tbaz*2));
		G.append(tG);
	G1 = np.array(G);
	G1 = np.dot(U,G1);
#	print G1;
	d = np.array(indat);
	d = d.T;
	d = np.dot(U,d);
#	print d;
	model = np.linalg.lstsq(G1,d)[0];
	resid = np.linalg.lstsq(G1,d)[1];
#	print result,resid;
#	sys.exit();
	A0 = model[0];
	A1 = math.sqrt(model[1]**2 + model[2]**2);
	fi1 = math.atan2(model[2],model[1]);
	A2 = math.sqrt(model[3]**2 + model[4]**2); 
	fi2 = math.atan2(model[4],model[3])

# compute forward:::
	odat = [];
	ccdat = [];
	ccdat = np.dot(G1,model);
	for i in range (len(inbaz)):
		ccdat[i] = ccdat[i]*inun[i];
		odat.append(ccdat[i]);
	return A0,A1,fi1,A2,fi2,odat;


def invert_A0 ( inbaz, indat, inun ):   #only invert for A0 part
        m = len(inbaz); # data space;
        n = 1; # model space;
        G = [];
        cvA0 = 0;
        U = np.zeros((m,m));
        for i in range (m):
                U[i][i] = 1/inun[i];
        for i in range (m):
                tG = [];
                tG.append(1);
		G.append(tG);
        G1 = np.array(G);
        G1 = np.dot(U,G1);

        d = np.array(indat);
        d = d.T;
        d = np.dot(U,d);

	model = np.linalg.lstsq(G1,d)[0];
	cvA0 = model[0];

	odat = [];
	ccdat = [];
        ccdat = np.dot(G1,model);
        for i in range (len(inbaz)):
                ccdat[i] = ccdat[i]*inun[i];
                odat.append(ccdat[i]);
#		print indat[i],ccdat[i];
	return cvA0,odat;


def invert_A2 ( inbaz, indat, inun ):
        m = len(inbaz); # data space;
        n = 3; # model space;
        G = [];
        A0 = 0;
        A1 = 0;
        fi1 = 0;
        U = np.zeros((m,m));
        for i in range (m):
                U[i][i] = 1/inun[i];
        for i in range (m):
                tG = [];
                tbaz = inbaz[i]/180*(math.pi);
                tG.append(1);
                tG.append(math.sin(tbaz*2));
                tG.append(math.cos(tbaz*2));
                G.append(tG);
        G1 = np.array(G);
        G1 = np.dot(U,G1);
        d = np.array(indat);
        d = d.T;
        d = np.dot(U,d);
        model = np.linalg.lstsq(G1,d)[0];
        resid = np.linalg.lstsq(G1,d)[1];
        A0 = model[0];
        A2 = math.sqrt(model[1]**2 + model[2]**2);
        fi2 = math.atan2(model[2],model[1]);
	
	odat = [];
        ccdat = [];
        ccdat = np.dot(G1,model);
        for i in range (len(inbaz)):
                ccdat[i] = ccdat[i]*inun[i];
                odat.append(ccdat[i]);

        return A0,A2,fi2,odat;


def invert_A1 ( inbaz, indat, inun ):
        m = len(inbaz); # data space;
        n = 3; # model space;
        G = [];
        A0 = 0;
        A2 = 0;
        fi2 = 0;
        U = np.zeros((m,m));
        for i in range (m):
                U[i][i] = 1/inun[i];
        for i in range (m):
                tG = [];
                tbaz = inbaz[i]/180*(math.pi);
                tG.append(1);
                tG.append(math.sin(tbaz));
                tG.append(math.cos(tbaz));
                G.append(tG);
        G1 = np.array(G);
        G1 = np.dot(U,G1);
        d = np.array(indat);
        d = d.T;
        d = np.dot(U,d);
        model = np.linalg.lstsq(G1,d)[0];
        resid = np.linalg.lstsq(G1,d)[1];
        A0 = model[0];
        A1 = math.sqrt(model[1]**2 + model[2]**2);
        fi1 = math.atan2(model[2],model[1]);

	odat = [];
        ccdat = [];
        ccdat = np.dot(G1,model);
        for i in range (len(inbaz)):
                ccdat[i] = ccdat[i]/U[i][i];
                odat.append(ccdat[i]);

        return A0,A1,fi1,odat;




def group ( inbaz, indat, outbaz,outdat,outun):
#	outbaz = [];
#	outdat = [];
#	print indat;
#	sys.exit();
	binwidth = 30;
	nbin = int((360+1)/binwidth);
#	print binwidth,nbin;
#	sys.exit();
	for i in range(nbin):
		mi = i*binwidth;
		ma = (i+1)*binwidth;
		tbaz = i*binwidth + float(binwidth)/2;
		ttdat = [];
#		tun = 0;
		for j in range (len(inbaz)):
			if inbaz[j] >=mi and inbaz[j] < ma:
				ttdat.append(indat[j]);
		if (len(ttdat) > 0):
			outbaz.append(tbaz);
			outdat.append(np.mean(ttdat));	
			if (len(ttdat)>1):
				outun.append(np.std(ttdat)/(math.sqrt(len(ttdat))));
			if (len(ttdat)==1):
				outun.append(0.1);



################### Main ################################################

if (len(sys.argv) != 2):
	print "input [in.lst]";
	sys.exit();



baz = [];
names = [];
nameoutlst=[]

for line in open(sys.argv[1]):
	l = line.rstrip().split();
	names.append(l[0]);
	baz.append(int(l[1]));
	nameoutlst.append(l[2])

atime = [];
adata = [];
rdata = []; # this is 0+1+2
rdata0 = []; # only 0
cvrdata1 = []; # 0+1
cvrdata2 = []; # 0+2
print id(cvrdata1), id(cvrdata2);
lens = [];
A0 = [];
dt = 0;
for i in range (len(names)):
	time = [];
	data = [];
	c = [];
	c1 = [];
	c2 = [];
	c3 = [];
	for ll in open(names[i],"r"):
		ll = ll.rstrip();
		ll1 = ll.split();
		ttime = float(ll1[0])-0.;
		if (ttime < -2.):
			continue;
		if (ttime > 20.):
			break;
		time.append(ttime-0.);
		data.append(float(ll1[1]));
		c.append(0.);
		c1.append(0.);
		c2.append(0.);
		c3.append(0.);

	atime.append(time);
	adata.append(data);
	rdata.append(c);
	rdata0.append(c1);
	cvrdata1.append(c2);
	cvrdata2.append(c3);
	lens.append(len(time));

#print "read over!!!";
dt = atime[0][1] - atime[0][0];

ofi1 = [];
#########################################parameters in 3 different inversion
zA0 = [];

oA0 = []
oA1 = [];
oSIG1 = [];

tA0 =[];
tA2 = [];
tSIG2 = [];

A0 = [];
A1 = [];
A2 = [];
SIG1 = [];
SIG2 = [];
A0123 = [];
MF0 = [];  # misfit between A0 and R[i]
MF1 = [];  # misfit between A0+A1+A2 and R[i]
MF2 = [];  # misfit between A0+A1+A2 and binned data

############################################################################
#GROUPED DATA
############################################################################
gbaz = [];
gdata = [];
gun = [];
############################################################################
#odat1 = [];

#print "now do !"
print lens,min(lens);
for i in range (min(lens)):
	tdat = [];
	for j in range (len(atime)):
		tdat.append(adata[j][i]);
	baz1 = [];
	tdat1 = [];
	udat1 = [];

	group(baz,tdat,baz1,tdat1,udat1);
	if (i>1220):
		print "group over",i;
	gbaz.append(baz1);
	gdata.append(tdat1);
	gun.append(udat1);

	
# now do inversions ########################################################

	tempv0 = 0.
	tempv1 = 0.
	tempv2 = 0.
	tempv3 = 0.
	tempv4 = 0.


	odat1 = [];
	if (i>1220):
		print baz1,tdat1,udat1;
	(tempv0,odat1)= invert_A0 (baz1,tdat1,udat1);
	zA0.append(tempv0);
	print "tempv0=",tempv0
	print odat1
	sys.exit()
#	fft = open("test0.dat","w");
#	for j in range (len(baz1)):
#		tempstr = "%g %g %g %g\n" % (baz1[j],tdat1[j],udat1[j],odat1[j]);
#		fft.write(tempstr);
#	fft.close();

#	print tempv0;
	if ( i > 1220):
		print "invert_A0!"
	odat1 = [];

	(tempv0,tempv1, tempv2, odat1) = invert_A1 (baz1,tdat1,udat1);
	oA0.append(tempv0);
	oA1.append(tempv1);
	oSIG1.append(tempv2);
	if ( i > 1220):
                print "invert_A1!"

#        fft = open("test1.dat","w");
#        for j in range (len(baz1)):
#                tempstr = "%g %g %g %g\n" % (baz1[j],tdat1[j],udat1[j],odat1[j]);
#                fft.write(tempstr);
#        fft.close();
#	print tempv0;

	if ( i > 1220):
                print "invert_A2!"

	odat1 = [];
        (tempv0,tempv1, tempv2,odat1) = invert_A2 (baz1,tdat1,udat1);
        tA0.append(tempv0);
        tA2.append(tempv1);
        tSIG2.append(tempv2);

#        fft = open("test2.dat","w");
#        for j in range (len(baz1)):
#                tempstr = "%g %g %g %g\n" % (baz1[j],tdat1[j],udat1[j],odat1[j]);
#                fft.write(tempstr);
#        fft.close();	
#	print tempv0;
	if ( i > 1220):
                print "invert_1!"	

	odat1 = [];
	(tempv0,tempv1, tempv2,tempv3,tempv4,odat1) = invert_1 (baz1,tdat1,udat1);
#	print tempv0, tempv1, tempv2, tempv3, tempv4;
#	sys.exit();
        A0.append(tempv0);
	A1.append(tempv1);
	SIG1.append(tempv2);
        A2.append(tempv3);
        SIG2.append(tempv4);
	if ( i > 1220):
		print "all inversion over",i;


	mf = 0.;
	mf1 = 0.
	for j in range (len(baz)):
		mf = mf + (tempv0 - adata[j][i])**2;
		vv = A3pre(baz[j],tempv0,tempv1,tempv2,tempv3,tempv4);
		mf1 = mf1 + (vv - adata[j][i])**2;
	mf = math.sqrt(mf/len(baz));
	mf1 = math.sqrt(mf1/len(baz));
	MF0.append(mf-0.);
	MF1.append(mf1-0.);
	if ( i > 1220):
		print "put in over!";

	mf2 = 0.;
	for j in range (len(baz1)):
		vv = A3pre(baz1[j],tempv0,tempv1,tempv2,tempv3,tempv4);
		mf2 = mf2 + (vv - tdat1[j])**2;
	mf2 = math.sqrt(mf2/len(baz1));
	MF2.append(mf2-0.);
	if ( i > 1220):
		print "A3pre now!";
	if ( i > 1220):
		print tempv0;

#	sys.exit();
# out put #########################################################
############## GROUPED DATA #######################################
print "now begin to group!";



###################################################################

#pen("A0.dat","w");
ff0 = open("A0.dat","w");
ff1 = open("A1.dat","w");
ff2 = open("A2.dat","w");
ff3 = open("A0_A1_A2.dat","w");

for i in range(min(lens)):

	#--write the A0 component signal, time,A0(time)
	ttA = zA0[i];
	if (ttA>-2 and ttA<2):
		tempstr = "%g %g\n" % (atime[0][i],ttA); 
		ff0.write(tempstr);

	#--write the A0+A1 component signal, time,A0(time),A1[time],phi1
	ttA = oA0[i];
	ttA1 = oA1[i];
	PHI1 = oSIG1[i];
        if (ttA>-2 and ttA<2):
		if (PHI1<0):
			PHI1 = PHI1 + math.pi;
                tempstr = "%g %g %g %g\n" % (atime[0][i],ttA,ttA1,PHI1);
		ff1.write(tempstr);

	#--write the A0+A2 component signal, time,A0(time),A2[time],phi2
	ttA = tA0[i];
	ttA2 = tA2[i];
	PHI2 = tSIG2[i];
	if (ttA>-2 and ttA<2):
		if (PHI2 < 0):
                        PHI2 = PHI2 + math.pi;
		tempstr = "%g %g %g %g\n" % (atime[0][i],ttA,ttA2,PHI2);
                ff2.write(tempstr);

	ttA = A0[i];
	ttA1 = A1[i];
	ttA2 = A2[i];
	PHI1 = SIG1[i];
	PHI2 = SIG2[i];

	
	#--the syn data at each Baz. rdata=A0+A1+A2; cvdata1=A1; cvdata2=A2
        for j in range(len(atime)):
                temp1 = ttA1*math.sin(float(baz[j])/180.*math.pi + PHI1);
                temp2 = ttA2*math.sin(2*float(baz[j])/180.*math.pi + PHI2);
                temp3 = temp1 + temp2 + ttA;
                rdata[j][i] = float(temp3);
                cvrdata1[j][i] = float(temp1);
                cvrdata2[j][i] = float(temp2);
        ofi1.append(PHI1);
	
	#--write A0+A1+A2 signal. time, A0(time), A1(time),phi1,phi2, misfit between A0 and R[i],misfit between A0+A1+A2 and R[i], misfit between A0+A1+A2 and binned data
        if ( ttA>-200 and ttA<200):
#                if (PHI1 < 0):
#                        PHI1 = PHI1 + math.pi;
#                if (PHI2 < 0):
#                        PHI2 = PHI2 + math.pi;
                tempstr = "%g %g %g %g %g %g %g %g %g\n" % (atime[0][i],ttA, ttA1,PHI1*180/math.pi,ttA2,PHI2*180/math.pi,MF0[i],MF1[i],MF2[i]);
#		print 
#		sys.exit();
                ff3.write(tempstr);
ff0.close();
ff1.close();
ff2.close();
ff3.close();

print "over !";
#########################################################################
sys.exit(); ##########EXIT######


ff1 = open("variance_reduction.dat","w");
vr0 = [];
vr1 = [];
vr2 = [];
vr3 = [];

vr11 = [];
vr21 = [];
vr31 = [];
vr32 = [];
for i in range (len(names)):
	tempbaz = baz[i];
	tempbaz1 = float(baz[i])*math.pi/180.;
	outname = "pre" + names[i];
	ff0 = open(outname,"w");


	obs = [];
	cvA0 =[];
	cvA1n = [];
	cvA1 = [];
	cvA2n = [];
	cvA2 = [];
	cvA3 = [];
	cvA3n1 = [];
	cvA3n2 = [];	



	cvdiffA0 = 0.;
	cvdiffA1 = 0.;
	cvdiffA2 = 0.;
	cvdiffA3 = 0.;
	cvobs = 0.;

	for j in range (min(lens)):
		if atime[0][j] > 10:
			break;
		raw = adata[i][j];

		cvA0.append(A0pre(tempbaz1,zA0[j]));

		cvA1.append(A1pre(tempbaz1,oA0[j],oA1[j],oSIG1[j]));
		cvA1n.append(A1pre1(tempbaz1,oA1[j],oSIG1[j]));

		cvA2.append(A2pre(tempbaz1,tA0[j],tA2[j],tSIG2[j]));
		cvA2n.append(A2pre1(tempbaz1,tA2[j],tSIG2[j]));

		cvA3.append(A3pre(tempbaz1,A0[j],A1[j],SIG1[j],A2[j],SIG2[j]));
		cvA3n1.append(A3pre1(tempbaz1,A1[j],SIG1[j]));
		cvA3n2.append(A3pre2(tempbaz1,A2[j],SIG2[j]));

		obs.append(raw);
#		print A0[j],A1[j],SIG1[j],A2[j],SIG2[j];
#		print outname;
#		print atime[0][j],cvA3;
#		sys.exit();		
		tempstr = "%g %g %g %g %g %g %g %g %g %g\n" % (atime[0][j], raw,cvA0[j],cvA1[j],cvA2[j],cvA3[j],cvA1n[j],cvA2n[j],cvA3n1[j],cvA3n2[j]);
		ff0.write(tempstr);

	ff0.close();
	vr0.append(match1(cvA0,adata[i]));
	vr1.append(match1(cvA1,adata[i]));
	vr2.append(match1(cvA2,adata[i]));
	vr3.append(match1(cvA3,adata[i]));



#	print vr0, vr1,vr2, vr3;
	tempstr = "%d %g %g %g %g %s\n" %(baz[i],vr0[i],vr1[i],vr2[i],vr3[i],names[i]);
	ff1.write(tempstr);
#	sys.exit();
	print names[i], " write ok!";
ff1.close();


ff1 = open("average_vr.dat","w");
tempstr = "%g %g %g %g\n" %(np.mean(vr0), np.mean(vr1), np.mean(vr2), np.mean(vr3));
ff1.write(tempstr);
ff1.close();
###########################################################################
#print cvrdata1[0][0];
#print cvrdata1[1][0];

for i in range (len(names)):
	outname = "rep" + names[i];
	ff = open(outname,"w");
	for j in range (min(lens)):
		tempstr = "%g %g\n" % (atime[0][j],rdata[i][j] );
		ff.write(tempstr);
	ff.close();


for i in range (len(names)):
        outname = "0rep" + names[i];
        ff = open(outname,"w");
        for j in range (min(lens)):
                tempstr = "%g %g\n" % (atime[0][j],A0[j] );
                ff.write(tempstr);
        ff.close();

for i in range (len(names)):
        outname = "1rep" + names[i];
        ff = open(outname,"w");
        for j in range (min(lens)):
                tempstr = "%g %g\n" % (atime[0][j],cvrdata1[i][j] );
                ff.write(tempstr);
        ff.close();


for i in range (len(names)):
        outname = "2rep" + names[i];
        ff = open(outname,"w");
        for j in range (min(lens)):
                tempstr = "%g %g\n" % (atime[0][j],cvrdata2[i][j] );
                ff.write(tempstr);
        ff.close();
