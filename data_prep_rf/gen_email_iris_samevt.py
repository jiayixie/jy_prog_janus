#only write events whose distance satisfeies some criteria
#write out two file: one is the email for requesting data; Another one is used to change the SAC header, by default, the data is named after its beginning time. 
# ALL STATIONS HAVE THE SAME NUMBER OF EVENTS
import sys
import string
import datetime
import os
sys.path.append('/home/jiayi/progs/jy/python')
import get_dist

if(len(sys.argv)!=4):
	print "Usage python xx.py 1]event_iris.lst 2]station.lst 3]title"
	sys.exit()
#====parameters========
deg_min=30.
deg_max=120.
tob=60 #time before t0
toe=90 #time after t0
cmp="BH?"
samplesac="/home/jiayi/Tool/sac1.sac"
c_stla = 38.86 # all stations have the same event list, events are selected according to the dist between evt and cst_location
c_stlo = -106.91
title = sys.argv[3]
fevt=sys.argv[1]
fsta=sys.argv[2]
stlo=[]
stla=[]
ntnm=[]
stnm=[]
Nst=0;
#======station.lst: stlo stla network stnm =====================
for line in open(fsta):  
	l=line.rstrip().split()
	stlo.append(float(l[0]))
	stla.append(float(l[1]))
	ntnm.append(l[2])
	stnm.append(l[3])
	Nst=Nst+1;

os.system("rm -f temp_taup.txt\n");
evttime=[];
evdeplst=[];
evlalst=[];
evlolst=[];
evmaglst=[];
iev=-1;
ftemp=open("check_taup.txt","w");
#======event file, the iris format ============================
for line in open(fevt,"r"):
    l1=line.rstrip().split()
    ya=l1[1][0:4]
    mon=l1[1][5:7]
    da=l1[1][8:10]
    hr=l1[2][0:2]
    min=l1[2][3:5]
    sec=l1[2][6:8]
    evla=l1[3][:-1]
    evlo=l1[4][:-1]
    evdep=l1[5][:-1]
    evmag=l1[9]
    temp_t=datetime.datetime(int(ya),int(mon),int(da),int(hr),int(min),int(sec))
    dist = get_dist.get_dist(float(evlo),float(evla),c_stlo,c_stla)
#====throw events according to its distance from c_stla c_stlo =========
    if (dist < deg_min*111 or dist > deg_max*111): 
        continue
    evttime.append(temp_t)    
    evdeplst.append(evdep);evlalst.append(evla);evlolst.append(evlo);evmaglst.append(evmag)
    iev=iev+1
#====uset Taup to get the first arrival time ==========================
    for ist in range(Nst):	
	st_dist=get_dist.get_dist(float(evlo),float(evla),stlo[ist],stla[ist])
	string = "taup_time -mod prem -h %f -ph P,Pdiff -km %f | awk '{if(NR==6){printf\"%d %d\";print $0}}' >> temp_taup.txt"%(float(evdep),st_dist,iev,ist)
	os.system(string);
	temp="%d %d : %s/%s/%s/%s/%s : %s\n"%(iev,ist,ya,mon,da,hr,min,string)
	ftemp.write(temp);
ftemp.close();
print "step1 finished!!"
print len(evttime),iev+1;

#====finish calculation begin to write ==========================================
fout=open("email_to_iris.txt","w")
fout2=open("SAC_ch.txt","w")
fout.write(".NAME Jiayi Xie\n.INST University of Colorado\n.MAIL Boulder,CO\n.EMAIL jiayi.seis@gmail.com\n.PHONE\n.FAX\n.MEDIA: Electronic(FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL %s\n.QUALITY B\n.END\n"%(title))
for line in open("temp_taup.txt"):
	l=line.rstrip().split()
	iev=int(l[0])
	ist=int(l[1])
	p=float(l[6])
	t_o=evttime[iev]
	tstnm=stnm[ist]
	tntnm=ntnm[ist]
	dt=datetime.timedelta(seconds=int(float(l[5]))) #ingore the time smaller than 1 sec
	t0=t_o+dt #the first p arrival time(GMT time)
	dt1=datetime.timedelta(seconds=tob)
	dt2=datetime.timedelta(seconds=toe)
	b=t0-dt1;
	e=t0+dt2;
#=======write out string for email, ntw sta begin_time end_time #cmp cmp==========================
	fout.write("%-5s %5s %4d %02d %02d %02d %02d %02d.00 %4d %02d %02d %02d %02d %02d.00 1 %s\n"%(tstnm,tntnm,b.year,b.month,b.day,b.hour,b.minute,b.second,e.year,e.month,e.day,e.hour,e.minute,e.second,cmp)) 
#=======write out change header macro, write in event info and o, t0, user0(ray param(s/deg))=====
#=======read a sample sac first, in case the target sac doesn't exist ==================================
	fout2.write("r %s\nr %4d.%03d.%02d.%02d.%02d*\n ch O GMT %4d %03d %02d %02d %02d 000 evlo %10s evla %10s evdp %s mag %s t0 %d user0 %f\n wh over\n"%(samplesac,b.year,b.timetuple().tm_yday,b.hour,b.minute,b.second,t_o.year,t_o.timetuple().tm_yday,t_o.hour,t_o.minute,t_o.second,evlolst[iev],evlalst[iev],evdeplst[iev],evmaglst[iev],tob,p))
fout.close()
fout2.close()

