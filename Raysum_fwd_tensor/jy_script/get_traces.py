# this read in the trace geometry file and the trace file, output each trace is separate files (time, azimuth, amplitude)
# the output trace from Raysum is rotated, R T Z
import os
import sys

fgeom = sys.argv[1]
ftrace = sys.argv[2]
outdir = sys.argv[3]

Ntrace=0
azilst=[];NSsftlst=[];EWsftlst=[]
for line in open(fgeom):
	if ("#" in line):
		#print line
		continue
	l=line.rstrip().split()
	azi=float(l[0])
	while(azi<0):
		azi=azi+360
	azilst.append(azi)
	NSsftlst.append(float(l[2]))
	EWsftlst.append(float(l[3]))
	Ntrace+=1

print "there are %d traces"%Ntrace

flagdata=0;i=0;Ndata=-100
for line in open(ftrace):
	i+=1
	if ( i==2):
		l=line.rstrip().split()
		Ntrace2=int(l[0]);Nsample=int(l[1]);dt=float(l[2]);shift=float(l[4])
		if(Ntrace2!=Ntrace):
			print "inconsistent trace numbers between %s and %s\n"%(fgeom,ftrace)
		continue
	elif ("Trace number" in line ):
		Ndata=i+2
		l=line.rstrip().split()
		idtrace=int(l[3])
		print "find trace %d"%idtrace
		fout1=open("%s/trace_%d.R.txt"%(outdir,idtrace),"w")
		fout2=open("%s/trace_%d.T.txt"%(outdir,idtrace),"w")
		fout3=open("%s/trace_%d.Z.txt"%(outdir,idtrace),"w")
		#continue

	elif(flagdata==0):
		if(i==Ndata-1):
			flagdata=1;	
			count=0
		continue
	else:
		count+=1
		if(count>Nsample):
			flagdata=0
			fout1.close()
			fout2.close()
			fout3.close()
			continue
		l=line.rstrip().split()
		t=(count-1)*dt	
		azi=azilst[idtrace-1]
		fout1.write("%5f %5f %5f\n"%(t,azi,float(l[0])))
		fout2.write("%5f %5f %5f\n"%(t,azi,float(l[1])))
		fout3.write("%5f %5f %5f\n"%(t,azi,float(l[2])))


