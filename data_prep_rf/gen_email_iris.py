#get the event list from iris event_list format.
#only write those events whose distance satisfeies some criteria
#write out two file: one is the email for requesting data; Another one is used to change the SAC header, by default, the data is named by its beginning time. 
import sys
import string
import datetime
import os
sys.path.append('/home/jiayi/progs/jy/python')
import get_dist

if(len(sys.argv)!=3):
	print "Usage python xx.py 1]event_iris.lst 2]station.lst"
	sys.exit()
deg_min=30.
deg_max=120.
tob=60
toe=90
cmp="BH?"

fevt=sys.argv[1]
fsta=sys.argv[2]
stlo=[]
stla=[]
ntnm=[]
stnm=[]
Nst=0;
for line in open(fsta):
	l=line.rstrip().split()
	stlo.append(float(l[0]))
	stla.append(float(l[1]))
	ntnm.append(l[2])
	stnm.append(l[3])
	Nst=Nst+1;
c_stlo = 33.61 #II.PFO
c_stla = -116.46
#os.system("rm -f temp_taup.txt\n");
evttime=[];
evdeplst=[];
evlalst=[];
evlolst=[];
evmaglst=[];
iev=0;
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
    evttime.append(temp_t)    
    evdeplst.append(evdep);evlalst.append(evla);evlolst.append(evlo);evmaglst.append(evmag)
    dist = get_dist.get_dist(float(evlo),float(evla),c_stlo,c_stla)
    if (dist < deg_min*111 or dist > deg_max*111):
        continue
    iev=iev+1
#    for ist in range(Nst):	
#	string = "taup_time -mod prem -h %f -ph P,ttp -km %f | awk '{if(NR==6){printf\"%d %d\";print $0}}' >> temp_taup.txt"%(float(evdep),dist,iev,ist)
#	os.system(string);
print "step1 finished!!"

fout=open("email_to_iris.txt","w")
fout2=open("ch_SAC.txt","w")
fout.write(".NAME Jiayi Xie\n.INST University of Colorado\n.MAIL Boulder,CO\n.EMAIL jiayi.seis@gmail.com\n.PHONE\n.FAX\n.MEDIA: Electronic(FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL RF_Jiayi_Oct24\n.QUALITY B\n.END\n")
for line in open("temp_taup.txt"):
	l=line.rstrip().split()
	iev=int(l[0])
	ist=int(l[1])
	p=float(l[6])
	t_o=evttime[iev]
	tstnm=stnm[ist]
	tntnm=ntnm[ist]
	dt=datetime.timedelta(seconds=float(l[5]))
	t0=t_o+dt #the first p arrival time(GMT time)
	dt1=datetime.timedelta(seconds=tob)
	dt2=datetime.timedelta(seconds=toe)
	b=t0-dt1;
	e=t0+dt2;
	#write out string for email, ntw sta begin_time end_time #cmp cmp
	fout.write("%-5s %5s %4d %02d %02d %02d %02d %02d.00 %4d %02d %02d %02d %02d %02d.00 1 %s\n"%(tstnm,tntnm,b.year,b.month,b.day,b.hour,b.minute,b.second,e.year,e.month,e.day,e.hour,e.minute,e.second,cmp)) 
	#write out change header macro, write in event info and o, t0, user0(ray param(s/deg))
	fout2.write("r %4d.%03d.%02d.%02d.%02d*\n ch O GMT %4d %03d %02d %02d %02d 000 evlo %10s evla %10s evdp %s mag %s t0 %d user0 %f\n wh over\n"%(b.year,b.timetuple().tm_yday,b.hour,b.minute,b.second,t_o.year,t_o.timetuple().tm_yday,t_o.hour,t_o.minute,t_o.second,evlolst[iev],evlalst[iev],evdeplst[iev],evmaglst[iev],tob,p))
fout.close()
fout2.close()

