// this is a function that's similar to the para_avg, but this function could handle the case of multiple peak.
// decide if there is multiple peaks(by defulat <=2 peaks) in the phi(strike) parameter. if there is, then group the parameters into 2 groups. one group with phi~phi1, ther other group with phi~phi1+90, and Npk=2; if there isn't, Npk=1
// do normal para_avg for each para group (#of para gp is Npk), and return avg, std

// this version could handle multiple(<=2) phi group in both crust and mantle. The final output has Nphi_group_crust*Nphi_group_mantle  paraavgs
// the seperation of mantle group can be disabled by setting idphiM<0

// this version, get parabest for each phi group
// previous v3 version has bug for identifying multiple gorup for 0.90,180 peak case; now change it by 1) find peak; 2) move all points to [peak-130,peak+60[ range 3) then do average and compare average&peak
// this version, the pk will be passed to the main program, only pass the pkC (peak of phi in the crust) value, not the pkM value
// only the phi value in the 1st layer/grid of that group has the right value, since only that layer/grid has converted phi value

// this version, para_avg_multiple_gp_v6.C, considers the averaging for AZcos AZsin case. need to convert the fast-direction correlated with them before doing averaging. AZcos, AZsin --> AMP,PHI,FA (function of PHI) --> convert FA --> convert PHI --> AZcos, AZsin --> do avg
// debugged the serperate_gp function

double convert(double vin,double vref,double T){
  double v;
  v=vin;
  //while(vref>T)vref-=T;
  //while(vref<0)vref+=T;

  while(v-vref>0.5*T)v-=T;
  while(v-vref<-0.5*T)v+=T;
  return v;
}//convert
//--------------------------

//int seperate_gp(vector<double> &vlst,int &Ngp, double &pk, vector<int> &indexflaglst){
int seperate_gp(vector<double> vlst,int &Ngp, double &pk, vector<int> &indexflaglst, double distcri, double T){
  //this is used to group the parameters according to the phi(strke) value in the crust (HOW ABOUT MANTLE?)
  // input the para list (a list of phi) that will be judged during the grouping process. 
  //return two index lists, each list belongs to one phi group
  // the input T and angles are in degree
  // the output pk is within [0,T]
  //#############PARAMETER
  //double distcri=5;//if |peak-avg|>distcri, then we think there are multiple groups of phi value
  double dv,dist,tv,dist1,dist2;
  //double T=180.; //period of phi
  double avg1,avg2,tmpv; //,pk;
  int i,j,k,Nv,nmax,Nbin;  
  vector<int> binlst;//,indexflaglst;  
  binlst.reserve(100);

  Nv=vlst.size();
  //printf("Begin another seperate_group, list size=%d\n",Nv);//---test---
  if(Nv==0){
	printf("### seperate_gp, the input list has no value inside!\n");
	exit(0);

  }//if Nv
  //--compute two kind of average values; avg1=sum(v)/N; avg2=sum(v if v<T/2; v-T if v>T/2)/N
  /*avg1=0.;avg2=0.;
  for(i=0;i<Nv;i++){
	if(vlst[i]<0 or vlst[i]>T){
		printf("#### seperate_gp the value (%g) in input list is outside the given range[0,%g]",vlst[i],T);
		exit(0);
	}
	avg1+=vlst[i];
	if(vlst[i]<T/2)avg2+=vlst[i];
	else{avg2+=(vlst[i]-T);}
  }//for i<Nv
  avg1/=Nv;
  avg2/=Nv;  
  */
  //-- compute the peak value; bin the list and get the value correlated with most data
  dv=5.;//bin width
  Nbin=int(T/dv);
  if(int(T*10)%int(dv*10)!=0){
	printf("### seperate_gp, the dv(%g) is not good! mod(T,dv)=mod(%g,%g)!=0\n",dv,T,dv);
	exit(0);
  }
  
  for(i=0;i<Nbin;i++){
	binlst.push_back(0);
  }
  //FILE *ftemp2;
  //ftemp2=fopen("temp_v.txt","w");
  for(i=0;i<Nv;i++){
        while(vlst[i]<0)vlst[i]+=T;
        while(vlst[i]>T)vlst[i]-=T;
	k=int(floor(vlst[i]/dv));
	//fprintf(ftemp2,"%d %d %g\n",i,k,vlst[i]);
	if(k==Nbin){printf("#### seperate_gp, strange case, vlst[%d]=%g,k=%d==%d\n",i,vlst[i],k,Nbin);exit(0);}
	binlst[k]++;
  }
  //fclose(ftemp2);  
  //get the value with max binlst value
  k=0;
  nmax=-100;
  //FILE *ftemp;
  //ftemp=fopen("temp_bin.txt","w");
  for(i=0;i<Nbin;i++){
	if(binlst[i]>nmax){nmax=binlst[i];k=i;}
	//fprintf(ftemp,"%d %f %d\n",i,(i+0.5)*dv,binlst[i]);
  }
  //fclose(ftemp);
  pk=(k+0.5)*dv;//the mid value of that bin; normally, I think pk is correlated with the group with c<=0, but I might be wrong
  
  //--compute avg in the [peak-130,peak+60] range
  avg1=0.;
  for(i=0;i<Nv;i++){
	tmpv=vlst[i];
	while(tmpv<pk-130)tmpv+=T;
	while(tmpv>pk+(T-130))tmpv-=T;
	avg1+=tmpv;
  }
  avg1/=Nv;
  //printf("hey, peak value = %g avg=%g === ibin=%d\n",pk,avg1,k);//---check---

  //-- choose between avg1, avg2 based on their dist to pk, and get the |avg-pk|
  avg1=convert(avg1,pk,T);
  //avg2=convert(avg2,pk,T);
  dist=fabs(avg1-pk);
  //dist2=fabs(avg2-pk);
  //dist=dist1<dist2?dist1:dist2;
  
  //--decide if there is multiple peaks in the vlst based on the value of dist
  if(dist>distcri)Ngp=2;
  else Ngp=1;

  printf("dist=%g distcri=%g, Ngp=%d\n",dist,distcri,Ngp);//--check--
  //--group the index of vlst into Ngp group(s)
  indexflaglst.clear();
  if(Ngp==1){//only one group
	for(i=0;i<Nv;i++){
		indexflaglst.push_back(1);//this lst tells if this para belongs to gp1 or 2 
		tv=vlst[i];
		vlst[i]=convert(tv,pk,T);		
	}//for i<Nv
  }//if Ngp==1
  else{//two groups, one with center value around pk, the other around pk+T*0.5 (pk+90)
 	for(i=0;i<Nv;i++){
		tv=vlst[i];
		tv=convert(tv,pk,T);
		if(fabs(tv-pk)<0.23*T){//### 0.23 is an arbitrary number
			//this para belongs to the groups with center value ~pk
			vlst[i]=tv;
			indexflaglst.push_back(1);	
			continue;
		}
		tv=vlst[i];
		tv=convert(tv,pk+0.5*T,T);
		if(fabs(tv-pk-0.5*T)<0.23*T){//belongs to the gp2
			//printf("%g->%g\n",vlst[i],tv);//---check--
			vlst[i]=tv;
			indexflaglst.push_back(2);
			continue;
		}
		indexflaglst.push_back(0);
		printf("---value %g, no gp (peak~%g and %g) wants it\n",vlst[i],pk,pk+0.5*T);//---check---
	}//for i<Nv
  }//else two groups
  if(indexflaglst.size()!=Nv){printf("something wrong, size doesn't match! %d!=%d\n",indexflaglst.size(),Nv);exit(0);}//---check---
  
  printf("pause\n");
  return 1;//indexflaglst;
  //return vlst;
}//seperate_gp
//--------------------------
int cs2ap_v2(double Ac,double As, double &amp,double &phi, double &fa, int phiflag){
// this is slightly different from the cs2ap (in the CALforward*.C). This does not give the fast direction, but just the angle directely correlated with sin and cos
  //the returned phi & fa are in radius
  //Ac*cos(X*theta)+As*sin(X*theta)=A*sin(X*theta+phi)
  //where A=sqrt(Ac*Ac+As*As); phi=atan(Ac/As); fast_axis=T/4-phi/X, and T=2pi/X
  double T=M_PI*2/phiflag;
  amp=sqrt(Ac*Ac+As*As);
  phi=atan2(Ac,As); //rad
  fa=T/4.-phi/phiflag;
  return 1;
}//cs2ap_v2
//--------------------------
int ap2cs(double amp,double phi,double &Ac,double &As){
  //phi is in radius
  //for a givin phi and amp, compute the coefficient for cos and sin
  //A*sin[X*theta+phi]=Ac*cos(X*theta)+As*sin(X*theta)
  //==> Ac=A*sin(phi); As=A*cos(phi)
  if(fabs(phi)>10){
	printf("### ap2cs, Hey, the input phi should be in radius NOT degree! here phi=%g\n",phi);
	exit(0);
  }
  Ac=amp*sin(phi);
  As=amp*cos(phi);
  return 1;
}//ap2cs
//--------------------------
vector<int> convert_AZpara(vector<paradef> &paralst,vector<int> AZcosidlst, vector<int> idlst, int Ngood, int &Ngp, double sepangcri){
  // this is used to group and convert the AZcos, AZsin: AZcos, AZsin --> FA, AMP --> convert FA, move values close by taking the periodicity into account --> converted AZcos, AZsin based on AMP and converted_FA
  // return the indexflaglst of the last AZcos para, so by default, the input AZcos list should have the same angle
  // --the sepangcri define the criteria that defines a parameter has multiple group or not. generally, this value is 5-10. But, depend on the need, you can set it very large, to eliminate multiple-group seperation.
  
  int i,j,ip,k;
  //int size;
  int Rphi;
  //vector<double> AZcoslst,AZsinlst;
  double T,AZcos,AZsin,phi,fa,amp,pk,pk2;
  vector<double> anglst,amplst;
  vector<int> indexflaglst;
  double rad2deg,deg2rad;

 
  rad2deg=180./M_PI;
  deg2rad=M_PI/180.;
  indexflaglst.reserve(Ngood);
  anglst.reserve(Ngood);
  amplst.reserve(Ngood);

  //########## IMPORTANT PARAMETER ########### CONTROL WHICH RAYLEIGH-WAVE LOVE-WAVE AZI TO COMPUT
  Rphi=2;
  T=M_PI*2/Rphi*rad2deg; // in degree

  for(i=0;i<AZcosidlst.size();i++){
	ip=AZcosidlst[i];
	anglst.clear();
	amplst.clear();
	indexflaglst.clear();
	//---test--- record the cos, sin, and fa before and after the change
  	vector<double> Vcosb,Vsinb,Vcosa,Vsina,fab,faa;
	for(j=0;j<Ngood;j++){
		k=idlst[j];
		AZcos=paralst[k].parameter[ip];
		AZsin=paralst[k].parameter[ip+1];
		cs2ap_v2(AZcos,AZsin,amp,phi,fa,Rphi); // fa and phi are in rad
		anglst.push_back(fa*rad2deg);
		amplst.push_back(amp);

		//---test--- record the cos, sin, and fa before and after the change
		Vcosb.push_back(AZcos);	
		Vsinb.push_back(AZsin);	
		fab.push_back(fa*rad2deg);
	}//j<size
	seperate_gp(anglst,Ngp,pk,indexflaglst,sepangcri,T);
	printf("in convert_AZ, the %dth AZcospara, peak=%g\n",i,pk);
	//-----------
	if(Ngp==1){
		for(j=0;j<Ngood;j++){
			fa=anglst[j]=convert(anglst[j],pk,T);//deg
			amp=amplst[j];
			phi=(T/4-fa)*Rphi*deg2rad;//rad
			//anglstcvt.push_back(ang);//---test---
			ap2cs(amp,phi,AZcos,AZsin);
			paralst[idlst[j]].parameter[ip]=AZcos;
			paralst[idlst[j]].parameter[ip+1]=AZsin;

			//---test--- record the cos, sin, and fa before and after the change
			Vcosa.push_back(AZcos);	
			Vsina.push_back(AZsin);
			faa.push_back(fa);
		}//for j
	}//if Ngp
	else{
		pk2=pk+0.5*T;
		while(pk2>T)pk2-=T;
		while(pk2<0)pk2+=T;
		printf("#### convert_AZpara, dealing with the %dth AZcos parameter (%dth parameter), there are %d groups of angles peak=%g and %g\n",i,ip,Ngp,pk,pk2); 
//---to be modified----
        	for(j=0;j<Ngood;j++){
                	if(indexflaglst[j]==1){//belongs to gp1
				fa=anglst[j]=convert(anglst[j],pk,T);//deg
			}
                	else if (indexflaglst[j]==2){//belongs to gp2
				fa=anglst[j]=convert(anglst[j],pk2,T);//deg
                	}
			else if (indexflaglst[j]==0){//does not belong to any groups
				fa=anglst[j];
			}
			else{
				printf("### convert_AZ, hey, the indexflaglst has problem! value should be either 0, 1 or 2. but here it equals to %d\n",indexflaglst[j]);
				exit(0);
			}
			amp=amplst[j];
			phi=(T/4-fa)*Rphi*deg2rad;//rad
			ap2cs(amp,phi,AZcos,AZsin);
			//printf("cos: %g -> %g sin: %g-> %g\n",paralst[idlst[j]].parameter[ip],AZcos,paralst[idlst[j]].parameter[ip+1],AZsin);
                        paralst[idlst[j]].parameter[ip]=AZcos;
                        paralst[idlst[j]].parameter[ip+1]=AZsin;
			//---test--- record the cos, sin, and fa before and after the change
			Vcosa.push_back(AZcos);	
			Vsina.push_back(AZsin);

  			phi=atan2(AZcos,AZsin);
  			phi=phi/Rphi;
			double T2=M_PI*2/Rphi;
  			fa=T2/4.-phi; //fast axis direction
  			while(fa>T2)
          			fa-=T2;
  			while(fa<0)
          			fa+=T2;
			fa=fa*180./M_PI; //rad2deg
			faa.push_back(fa);
        	}//for i
	

		//exit(0); //---test---
	}//else Ngp>1
	printf("hey, finish converting paraCos%d para%d\n",i,ip);
  }//for i < AZcosidlst.size()
  
  
  return indexflaglst;
  //return 1;
}//convert_AZpara

//--------------------------
vector<double> compute_paraavg(vector<paradef> paralst,vector<int> idlst, vector<double> &stdlst, int id ){
        int i,j,k;
        int Npara,Ngood;
        vector<double> vavglst;

        Npara=paralst[0].npara;
        Ngood=idlst.size();

	//printf("begin compute_paraavg %d %d\n",paralst[0].parameter.size(),paralst[0].LoveRAparameter.size());//--check--
        if((id-1)*(id-2)!=0){printf("### compute_paraavg, wrong value for id! should be 1 or 2\n");exit(0);}

        stdlst.clear();
        for(i=0;i<Npara;i++){
                vavglst.push_back(0.);
                stdlst.push_back(0.);
        }

	vector<double> tempV;
        for(j=0;j<Npara;j++){
		tempV.clear();
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        if(id==1){vavglst[j]+=paralst[k].parameter[j];tempV.push_back(paralst[k].parameter[j]);}
                        else if(id==2){vavglst[j]+=paralst[k].LoveRAparameter[j];tempV.push_back(paralst[k].LoveRAparameter[j]);}
                }
                vavglst[j]/=Ngood;
        }//for j<Npara

        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        if(id==1)stdlst[j]+=pow(paralst[k].parameter[j]-vavglst[j],2);     
                        else {stdlst[j]+=pow(paralst[k].LoveRAparameter[j]-vavglst[j],2);}
                }
                stdlst[j]=sqrt(stdlst[j]/Ngood);
		//---check--
		if(stdlst[j]>50){
                        printf("###compute_paraavg, something wrong(?),para%d, avg=%g,large std %g\n",j,vavglst[j],stdlst[j]);
			//exit(0);
		}
		//
	}//for j<Npara                                                                              
  //printf("finish compute_paraavg\n");//--check--
  return vavglst;
}//compute_paraavg
//--------------------------
vector<vector<double> > compute_paraavgAZ(vector<paradef> paralst,vector<int> idlst, vector<vector<double> > &stdlst){
        int i,j,k;
        int Npara,Ngood;
        vector<double> LAZp0(2,0.);
	vector<vector<double> > vavglst;

        Npara=paralst[0].npara;
        Ngood=idlst.size();

	//printf("begin compute_paraavg %d %d\n",paralst[0].parameter.size(),paralst[0].LoveRAparameter.size());//--check--

        stdlst.clear();
        for(i=0;i<Npara;i++){
                stdlst.push_back(LAZp0);
		vavglst.push_back(LAZp0);
        }

	vector<double> tempV1,tempV2;
	tempV1.reserve(1000);tempV2.reserve(1000);
        for(j=0;j<Npara;j++){
		tempV1.clear();tempV2.clear();
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        vavglst[j][0]+=paralst[k].LoveAZparameter[j][0];
			vavglst[j][1]+=paralst[k].LoveAZparameter[j][1];                   
			tempV1.push_back(paralst[k].LoveAZparameter[j][0]); tempV2.push_back(paralst[k].LoveAZparameter[j][1]);
                }
                vavglst[j][0]/=Ngood;
		vavglst[j][0]/=Ngood;
        }//for j<Npara

        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
			stdlst[j][0]+=pow(paralst[k].LoveAZparameter[j][0]-vavglst[j][0],2);     
			stdlst[j][1]+=pow(paralst[k].LoveAZparameter[j][1]-vavglst[j][1],2);
			
                }
                stdlst[j][0]=sqrt(stdlst[j][0]/Ngood);
		stdlst[j][1]=sqrt(stdlst[j][1]/Ngood);
		
	}//for j<Npara       
                                                                   
  //printf("finish compute_paraavg AZ\n");//--check--
  return vavglst;
}//compute_paraavgAZ
//--------------------------
vector<int> para_avg_multiple_gp(vector<int> idphilst,vector<int> idphiMlst, vector<paradef> &paralst, vector<paradef> &parabestlst, vector<paradef> &paraavglst, vector<paradef> &parastdlst, vector<vector<int> > &idlstlst, int flag, double &pkC, vector<int> AZcosidlst, vector<int> AZcosidlstM,modeldef model, int flagupdaterho, paradef para0){
  // in this function, if idphiM<0 then, won't do mantle group seperation based on mantle phi
  // flag indicate if average is for para.parameter(flag=1) or for para.parameter/LoveRAparameter/LoveAZparameter (flag=3)
  // the newly added (May 14, 2015) para0 parameter is used to transfer the para0.para0 information (actually, this is useless. previously, i was worrying that the periodicity in phi(180deg) will affect the Gc,s Bc,s Hc,s Ec,s parameters, but actually they will not be affected, that is phi, phi+/-180 give the same  Gc,s Bc,s Hc,s Ec,s for both elliptical and non-elliptical Hex tensors)
  // on Jun 17, 2015, modified the idphi to idphilst. now, we will convert all the given phi parameters (based on idphilst). But there is a problem, right now, the Ngp(M) and indexflaglst[] is based purely on the last idphi (i.e., idphilst[size-1]). Actually, should consider Ngp=Ngp_node1*Ngp_node2*...*Ngp_nodeN, maybe do this later

  int i,j,k,n,size,idmin,Ngp,NgpM,igp,pflag,Ngood,idphi,idphiM;
  vector<int> indexflaglst,indexflaglstM,idlst,idlstnew,idminlst;
  //vector<int> idlstnewtest;
  vector<double> philst,philstM;  
  double mismin,pkM,pk,tv,T;  
  //vector<double> parastd;
  paradef paraavg,parastd;
  T=180.; //period of phi

  if((flag-1)*(flag-3)!=0){
	printf("### para_avg_multiple_gp, flag should be 1(compute avg for para.parameter), or 3 (compute avg for para.parameter, para.LoveRAparameter, and para.LoveAZparameter), but input is %d\n",flag);
	exit(0);
  }
  //printf("hey begin!\n");//--check--
  idlstlst.clear();
  size=paralst.size();
  if(size<1) return idminlst;
 
  paraavg=paralst[0];
  //parabest=paralst[0];
  parastd=paralst[0];
  //parastd.clear();
 
  //--initialize the paralst and parastd --> will be handled by the function compute_paraavg


  //--get the smallest misfit and get the index lst for model within misfit criteria-- this is just a groups of model with acceptable misfit; then after group seperation, we will re-select the model with higher criteria;(two criteria, one in seperate_gp[based on its dist to the avg phi value], one in the later part of this subroutine[based on the mismin of that phi gorup])
  mismin=1e10;
  for(i=0;i<size;i++){
	if(paralst[i].misfit<mismin){idmin=i;mismin=paralst[i].misfit;}
  }
  //parabest=paralst[idmin];
  printf("@@@ para_avg mismin=%g\n",mismin);
  //mismin=mismin+0.5;
  mismin=(mismin*2>(mismin+1.0))?mismin*2:(mismin+1.0); //arbitrary selection criteria
  //mismin=1e10;
  printf("@@@ selection criter, misfit<%g\n",mismin);
  Ngood=0;
  idlst.clear();
  for(i=0;i<size;i++){
	if(paralst[i].misfit<mismin){
		Ngood++;
		idlst.push_back(i);
	}
  }


  //--- by default, a group (e.g., crust) can only be either TTI/TI/iso or AZcos, cannot be both at the same time
  if(AZcosidlst.size()>0){// if this group is AZcos
  //---- convert the AZcos AZsin parameters; group the angle correlated with AZcos and AZsin, then recompute AZcos and AZsin. 
  indexflaglst=convert_AZpara(paralst, AZcosidlst,idlst,Ngood,Ngp,10.);
  }
  else{// else this group is TTI or TI/iso
  /*
  //---get philst(only from mod with small misfit), and call function seperate_gp to group the philst
  //philst.clear();
  if(idphilst.size()>0){
  for(n=0;n<idphilst.size();n++){
  idphi=idphilst[n];
  for(i=0;i<Ngood;i++){
	k=idlst[i];
	philst.push_back(paralst[k].parameter[idphi]);}  

  Ngp=0;
  seperate_gp(philst,Ngp,pkC,indexflaglst,5.,T);
  }// for n
  }//if idphi>=0
  else{
    Ngp=1; // no idphi informatino, then do not seperate the group
    for(i=0;i<Ngood;i++)indexflaglst.push_back(1);
  }
  */
  //--- put the modified (+/-T) phi back into the paralst; and compute the avg for all parameters for each gp(seperated based on philst)
  /*for(i=0;i<Ngood;i++){
	// since I disabled the &vlst in the seperate_gp function, the philst is not changed, so no need to transfer its value back to paralst; I disabled it only b.c. of the memory problem correlated with the '&'
	k=idlst[i];
  }*/
  //

  //--since philst is not transfered back, need to modify the phi values here
  if(idphilst.size()>0){
  for(n=0;n<idphilst.size();n++){
    philst.clear();
    idphi=idphilst[n];
    for(i=0;i<Ngood;i++){
        k=idlst[i];
        philst.push_back(paralst[k].parameter[idphi]);}
    Ngp=0;
    seperate_gp(philst,Ngp,pkC,indexflaglst,5.,T);
    pk=pkC;
    if(Ngp==1){//only one group
	for(i=0;i<Ngood;i++){
		k=idlst[i];
		tv=paralst[k].parameter[idphi];
		paralst[k].parameter[idphi]=convert(tv,pk,T);		
	}//for i<Nv
    }//if Ngp==1
    else{//two groups, one with center value around pk, the other around pk+T*0.5 (pk+90)
 	for(i=0;i<Ngood;i++){
		if(indexflaglst[i]==1){
			tv=paralst[idlst[i]].parameter[idphi];
			paralst[idlst[i]].parameter[idphi]=convert(tv,pk,T);
		}
		else if (indexflaglst[i]==2){
			tv=paralst[idlst[i]].parameter[idphi];
			paralst[idlst[i]].parameter[idphi]=convert(tv,pk+0.5*T,T);
			//paralst[idlst[i]].para0=para0.para0;
			//Vpara2Lovepara(paralst[idlst[i]],model,flagupdaterho);//---modified May 6, 2015
		}
	}//for i
    }//else two groups
  }// for n<idphilst.size()
  }//if idphi>=0
  else{
    indexflaglst.clear();
    Ngp=1; // no idphi informatino, then do not seperate the group
    for(i=0;i<Ngood;i++)indexflaglst.push_back(1);
  }//else
  //---
  }// else AZcosidlst.size()

  if(AZcosidlstM.size()>0){// if this group is AZcos
  //---- convert the AZcos AZsin parameters; group the angle correlated with AZcos and AZsin, then recompute AZcos and AZsin. Without this step, the averaged values. 
  indexflaglstM=convert_AZpara(paralst,AZcosidlstM,idlst,Ngood,NgpM,100.);// here I set the sepangcri=100. to eliminate group seperation in the mantle
  }
  else{// else this group is TTI or TI/iso
  //## if want to seperate the model also based on mantle phi,then
  if(idphiMlst.size()>0){ // do do group seperation for mantle
  for(n=0;n<idphiMlst.size();n++){
    idphiM=idphiMlst[n];
    philstM.clear();
    for(i=0;i<Ngood;i++){
	k=idlst[i];
	philstM.push_back(paralst[k].parameter[idphiM]);
    }
    NgpM=0;
    seperate_gp(philstM,NgpM,pkM,indexflaglstM,15.,T);
    /*for(i=0;i<Ngood;i++){
	k=idlst[i];
	paralst[k].parameter[idphiM]=philstM[i];
      }*/
    pk=pkM;
    if(Ngp==1){//only one group
	for(i=0;i<Ngood;i++){
		k=idlst[i];
		tv=paralst[k].parameter[idphiM];
		paralst[k].parameter[idphiM]=convert(tv,pk,T); 
	}//for i<Nv
    }//if Ngp==1
    else{//two groups, one with center value around pk, the other around pk+T*0.5 (pk+90)
 	for(i=0;i<Ngood;i++){
		if(indexflaglstM[i]==1){
			tv=paralst[idlst[i]].parameter[idphiM];
			paralst[idlst[i]].parameter[idphiM]=convert(tv,pk,T);
		}
		else if (indexflaglstM[i]==2){
			tv=paralst[idlst[i]].parameter[idphiM];
			paralst[idlst[i]].parameter[idphiM]=convert(tv,pk+0.5*T,T);
		}
	}//for i
    }//else two groups
  }// for n
  }//if idphiM>0 
  else{// do not do group seperation for mantle
	indexflaglstM.clear();
	NgpM=1;
	for(i=0;i<Ngood;i++)indexflaglstM.push_back(1);
  }
  //printf("hi 331\n");//--check---
  }// else AZcosidlstM.size()

//---to be modified----
  for(i=0;i<Ngp;i++){
	for(j=0;j<NgpM;j++){
  		idlstnew.clear();
		//---find the parabest for this phi group
  		//printf("hi 337\n");//--check---
		mismin=1e10;
		idmin=-1;
		for(k=0;k<Ngood;k++){
			if(paralst[idlst[k]].misfit<mismin and indexflaglst[k]==i+1 and indexflaglstM[k]==j+1){
				mismin=paralst[idlst[k]].misfit;
				idmin=idlst[k];
			}
		}
		idminlst.push_back(idmin);
		parabestlst.push_back(paralst[idmin]);
		//mismin=mismin+0.5;
		mismin=(mismin*2>(mismin+1.0))?mismin*2:(mismin+0.5); //arbitrary selection criteria
  		//printf("@@@ selection criter, misfit<%g\n",mismin);
		//mismin=1e10;//this is for prior distribution generation

  		//--prepare the new idlst
  		for(k=0;k<Ngood;k++){
			if(paralst[idlst[k]].misfit<mismin and indexflaglst[k]==i+1 and indexflaglstM[k]==j+1)idlstnew.push_back(idlst[k]);
			//idlstnew.push_back(idlst[k]);
		}//for k<Ngood
		printf("====the %dth group has %d models,mismin=%g selection criteria=====%g\n",i*Ngp+j,idlstnew.size(),paralst[idmin].misfit,mismin);//--check--
		if(flag==1)paraavg.parameter=compute_paraavg(paralst,idlstnew,parastd.parameter,1);
		else if(flag==3){
			paraavg.parameter=compute_paraavg(paralst,idlstnew,parastd.parameter,1);
			paraavg.LoveRAparameter=compute_paraavg(paralst,idlstnew,parastd.LoveRAparameter,2);
			paraavg.LoveAZparameter=compute_paraavgAZ(paralst,idlstnew,parastd.LoveAZparameter);
		}
		//printf("%d %d\n",idlstlst.size(),idlstnew.size());
		if(idlstnew.size()<2)continue;
		//idlstnewtest.clear();	
		//for(int nn=0;nn<100;nn++){
		//	idlstnewtest.push_back(idlstnew[nn]);
		//	printf("%d/%d idlstnew[%d]=%d\n",nn,idlstnew.size(),nn,idlstnew[nn]);
		//}

		idlstlst.push_back(idlstnew);
    		paraavglst.push_back(paraavg);
		parastdlst.push_back(parastd);
  	}//for j<NgpM
  }//for i<Ngp

  return idminlst;
  //return idmin;
}//para_avg_multiple_gp
