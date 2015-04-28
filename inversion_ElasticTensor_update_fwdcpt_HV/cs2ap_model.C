int model_cs2ap(vector<modeldef> &modlst, paradef para){
  //---this is used to convert the Acos, Asin values to amp and fast-axis for each model (at appropriate depth)
  int i,j,k,N,p6,p4;
  int N1,N2,phiflag,igp,Ngp;
  double Ac,As,amp,phi,fa;
  double dep,dep1,T;
  vector<int> gplst,iv;
  vector<double> tv;
  vector<vector<double> > Deplst2;
  vector<vector<int> > IDlst2;

  N=modlst.size();Ngp=modlst[0].ngroup;
  phiflag=2;
  T=360./phiflag;//deg
  Deplst2.reserve(N);
  IDlst2.reserve(N);
  
  //--check if there are AZ parameters, if not, then no need to run this function, return 0; if yes, find the group number that has AZ parameters
  for(i=0;i<para.npara;i++){
	p6=(int)para.para0[i][6];
	p4=(int)para.para0[i][4];
        if(p6==10){
		if(find(gplst.begin(),gplst.end(),p4)!=gplst.end())continue; //this gpid is already in gplst
		else gplst.push_back(p4);
	}
  }
  if(gplst.size()==0)return 0;
  //--find the layer number that corresponds to the group with AZ parameter (by default, the number of layer should be the same among models? so I only need the information from one model)
  for(i=0;i<N;i++){//store the depth of each group, for every model Deplst2[imod][igp]
	tv.clear();
	dep=0.;
	for(j=0;j<Ngp;j++){
		dep+=modlst[i].groups[j].thick;
		tv.push_back(dep);
    	}//for j
	Deplst2.push_back(tv);
  }// for i<N
  for(i=0;i<N;i++){//obtain the id for the depth of each group, for every model. IDlst2[imod][igp] 
	iv.clear();
	for(j=0;j<Ngp;j++){
		dep=Deplst2[i][j];
		dep1=0.;
		for(k=0;k<modlst[i].laym0.nlayer;k++){
			if(fabs(dep1-dep)<1e-4){iv.push_back(k);break;}
			dep1+=modlst[i].laym0.thick[k];
		}//for k
	}//for j
	IDlst2.push_back(iv);
  }//for i<N  

  for(i=0;i<N;i++){
	for(j=0;j<gplst.size();j++){
		igp=gplst[j];
		if(igp==0)N1=0;
		else N1=IDlst2[i][igp-1]+1;// at each interface, there are 2 values, one belong to the current group, one belong to the next group
		N2=IDlst2[i][igp];
		for(k=N1;k<=N2;k++){
			Ac=modlst[i].laym0.theta[k];
			As=modlst[i].laym0.phi[k];
			cs2ap_v2(Ac,As,amp,phi,fa,phiflag);
			modlst[i].laym0.theta[k]=amp;		
			fa=fa*180./M_PI;
			while(fa<0)fa+=T;
			while(fa>T)fa-=T;
			modlst[i].laym0.phi[k]=fa;// degree		
		}//for k
	}//for j<gplst.size	
  }//for i<N

  //--then loop around all the models, layers, perform cs2ap, and restore ap in the model.laym0.theta/phi

  return 1;
}// int model_cs2ap



