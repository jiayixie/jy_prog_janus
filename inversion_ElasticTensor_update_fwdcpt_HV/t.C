//-----------------------------------------------------	
	int model_avg_sub(vector<vector<float> > &vlst,vector<vector<float> > &stdlst,vector<float> &hlst,vector<modeldef> &modlst,int ng,float h0,float h,float dth)
	{
	  //h0=0;h=Hsed;hstd=hsedstd;dth=0.4;
	  //vlst[ithick]=[vsv,vsh,vpv,vph,eta,theta,phi,  vsvmin~phimin,  vsvmax~phimax];
	  //stdlst[ithick]=[vsvstd~phistd]

	  vector<float> tv,dep2lst,dep1lst;
	  vector<int> iddep1lst;
	  float th,tth,fm1,fm2,dep;
	  float vsv,vsh,dep1,dep2,vsvmin,vsvmax,vshmin,vshmax;
	  float tvsh,tvsv,tani,ani,animin,animax;
	  float fm3;
	  char str[500];
	  int i,j,Ncount,size,k,flag;
	  vector<vector<float> > tvlst;
	  float vpvmin,vpvmax,vphmin,vphmax,etamin,etamax,thetamin,thetamax,phimin,phimax;
	  float tvpv,tvph,teta,ttheta,tphi,vpv,vph,eta,theta,phi;

	  hlst.clear();vlst.clear();stdlst.clear();tvlst.clear();
		  
	  size=modlst.size();
	
	  for(i=0;i<size;i++){	// obtain the end depth for this mod_avg
	    dep2=0.;  
	    for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
	    dep2lst.push_back(dep2);
	    dep1lst.push_back(0.);
	    iddep1lst.push_back(0);
	  }
	  flag=0;
	  for(th=h0;th<=h+dth;th=th+dth){//thickness
		//printf("th=%g h=%g dep1[0]=%g dep2[0]=%g\n",th,h,dep1lst[0],dep2lst[0]);//---test---
		if(th>h){//modified on Nov 24, 2013, allow the th to go very close to h 
			if(flag==0){flag=1;th=h-0.001;}
			else{break;}
		}
		if(fabs(th-h+dth)<1e-4)th=h;
		vsv=0.;vsh=0.;Ncount=0;fm1=0.;fm2=0.;ani=0.;fm3=0.;tth=0.;
		vpv=0.;vph=0.;eta=0.;theta=0.;phi=0.;
		tvlst.clear();
		vsvmin=1e10;vsvmax=-1;
		vshmin=1e10;vshmax=-1;
		animin=1e10;animax=-1;
		vpvmin=1e10;vpvmax=-1.;
		vphmin=1e10;vphmax=-1.;
		etamin=1e10;etamax=-1.;
		thetamin=1e10;thetamin=-1.;
		phimin=1e10;phimax=-1.;
		for(i=0;i<size;i++){//iterate within model list
		  //dep2=0.;
		  //for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
		  dep2=dep2lst[i];// end depth
		  //dep1=dep2-modlst[i].groups[ng].thick;//dep1--upper of the layer; dep2--lower of the layer
		  if(th>dep2)continue;//tth=dep2; //tth = min(dep2,th) ####if th>dep2, continue, skip present mod, jumpt to next mod??, this means th is outside the present layer of present model, so skip.
		  else tth=th;
		  dep1=0;
		  for(j=0;j<modlst[i].laym0.nlayer-1;j++){ // modeified on Aug 27, 2012. For the smooth model, the thickness of the last layer is 0 
		    dep1=dep1+modlst[i].laym0.thick[j];
		  //for(j=iddep1lst[i];j<modlst[i].laym0.nlayer;j++){
		//	dep1=dep1+dep1lst[i]+modlst[i].laym0.thick[j];
			if(dep1>=tth ){//  the interpolation, there was a bug here (for smooth model, though it works for layered model), changed on Aug27, 2012. another bug modified on Oct 8, 2013
				tvsv=(modlst[i].laym0.vsv[j+1]-modlst[i].laym0.vsv[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vsv[j];
				tvsh=(modlst[i].laym0.vsh[j+1]-modlst[i].laym0.vsh[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vsh[j];
				tvpv=(modlst[i].laym0.vpv[j+1]-modlst[i].laym0.vpv[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vpv[j];
				tvph=(modlst[i].laym0.vph[j+1]-modlst[i].laym0.vph[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vph[j];
				teta=(modlst[i].laym0.eta[j+1]-modlst[i].laym0.eta[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.eta[j];
				ttheta=(modlst[i].laym0.theta[j+1]-modlst[i].laym0.theta[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.theta[j];
				tphi=(modlst[i].laym0.phi[j+1]-modlst[i].laym0.phi[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.phi[j];	
				
				if(tvsv<0 or tvsv>5){
					printf("test--- th=%g, imod=%d, ilay=%d, tvsv=%g=[(%g-%g)/(%g)+%g]\n",th,i,j,tvsv,modlst[i].laym0.vsv[j+1],modlst[i].laym0.vsv[j],modlst[i].laym0.thick[j],modlst[i].laym0.vsv[j]);}
				tani=100*(tvsh-tvsv)/(sqrt((2*tvsv*tvsv+tvsh*tvsh)/3.0));
				vsv=vsv+tvsv;// modified on Aug27, 2012
				vsh=vsh+tvsh;//
				ani=ani+tani;
				vpv=vpv+tvpv;
				vph=vph+tvph;
				eta=eta+teta;
				theta=theta+ttheta;
				phi=phi+tphi;
				Ncount++;	
				if (tvsv>vsvmax)vsvmax=tvsv;
				if(tvsv<vsvmin)vsvmin=tvsv;
				if(tvsh>vshmax)vshmax=tvsh;
				if(tvsh<vshmin)vshmin=tvsh;
				if(tani>animax)animax=tani;if(tani<animin)animin=tani;
				if(tvpv>vpvmax)vpvmax=tvpv;if(tvpv<vpvmin)vpvmin=tvpv;
				if(tvph>vphmax)vphmax=tvph;if(tvph<vphmin)vphmin=tvph;
				if(teta>etamax)etamax=teta;if(teta<etamin)etamin=teta;
				if(ttheta>thetamax)thetamax=ttheta;if(ttheta<thetamin)thetamin=ttheta;
				if(tphi>phimax)phimax=tphi;if(tphi<phimin)phimin=tphi;
				//---tvlst: vsv,vsh ..., phi
				tv.clear();tv.push_back(tvsv);tv.push_back(tvsh);//tv.push_back(tani);
				tv.push_back(tvpv);tv.push_back(tvph);tv.push_back(teta);tv.push_back(ttheta);tv.push_back(tphi);
				tvlst.push_back(tv);
			//	dep1lst[i]=dep1-modlst[i].laym0.thick[j];//
			//	iddep1lst[i]=j;//
				break;
			}//if dep>tth
		  }//for j Nlay
		}//for i Nmodel,size
		vsv=vsv/Ncount;
		vsh=vsh/Ncount;
		ani=ani/Ncount;
		vpv=vpv/Ncount;
		vph/=Ncount;
		eta/=Ncount;
		theta/=Ncount;
		phi/=Ncount;
		if (isnan(vsv) or isnan(vsh) or isnan(ani) or isnan(vpv) or isnan(vph) or isnan(eta) or isnan(theta) or isnan(phi)){
			printf("Hey, NaN happen!!! something wrong!!\n Ncount=%d th=%g tth=%g h0=%g h=%g\n",Ncount,th,tth,h0,h);
			sprintf(str,"echo WRONG Ncount = %d >> point_finished_Ani.txt ",Ncount);
			system(str);
			//exit(0);
		}
		//---vlst: vsv~phi, vsvmin~phimin, vsvmax~phimax---
		tv.clear();
		//tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vsvmin);tv.push_back(vsvmax);tv.push_back(vshmin);tv.push_back(vshmax);tv.push_back(ani);tv.push_back(animin);tv.push_back(animax);
		tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vpv);tv.push_back(vph);tv.push_back(eta);tv.push_back(theta);tv.push_back(phi);
		tv.push_back(vsvmin);tv.push_back(vshmin);tv.push_back(vpvmin);tv.push_back(vphmin);tv.push_back(etamin);tv.push_back(thetamin);tv.push_back(phimin);
		tv.push_back(vsvmax);tv.push_back(vshmax);tv.push_back(vpvmax);tv.push_back(vphmax);tv.push_back(etamax);tv.push_back(thetamax);tv.push_back(phimax);

		vlst.push_back(tv);
		hlst.push_back(th);
	
		//---stdlst: vsvstd~phistd
		tv.clear();
		fm1=0.;fm2=0.;	
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][0]-vsv,2);
			fm2=fm2+pow(tvlst[i][1]-vsh,2);//----there was a bug here, fixed on May 18, 2012
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);
		fm1=0.;fm2=0.;		
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][2]-vpv,2);
			fm2=fm2+pow(tvlst[i][3]-vph,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);
		fm1=0.;fm2=0.;fm3=0.;
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][4]-eta,2);
			fm2=fm2+pow(tvlst[i][5]-theta,2);
			fm3=fm3+pow(tvlst[i][6]-phi,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		fm3=sqrt(fm3/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);tv.push_back(fm3);
		
		stdlst.push_back(tv);
	  }//for th
        vector<vector<float> >().swap(tvlst);
	vector<float>().swap(tv);
	vector<float>().swap(dep2lst);
	vector<float>().swap(dep1lst);
	vector<int>().swap(iddep1lst);
	return 1;
	}//mod avg sub

//-----------------------------------------------------------------
