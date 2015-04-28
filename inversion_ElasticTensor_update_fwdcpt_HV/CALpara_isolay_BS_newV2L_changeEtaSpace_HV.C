/*========CONTENT============
int initpara(paradef &para)
int checkParaModel(paradef para, modeldef model)
int readpara(vector<vector<double> > &para0, const char* fname)
int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)
int gen_newpara(paradef inpara, modeldef model, paradef &outpara, int pflag)
int para_avg(vector<paradef> &paralst, paradef &paraavg,vector<double> &Rparastd,vector<double> &Lparastd, vector<int> &idlst)
// this version ,change the space of eta, require it to be <1.1
// modified Mar12, 2014;For eta parameter,  when the sigma <-5, use the constant c in that group; when sigma>-5and<-3 use constant eta/RAvs; (c measures the non-ellipticity of the tensor)
//============================
// in this version, around line 704, change the way dealing with newpara that's outside the model space range. previously, if the newpara is outside the range, then we won't accept it, and will re-generate a newpara until its within the range--> but this makes the prior non-uniform. now, we make a reflection of the para that's outside the model space, we think this will make the prior more uniform
// this version (HV), change checkParaModel (reorder it in a more reasobale way, add/rm some criteria), change the gen_newpara (take the scaling into account at this step), and para2mod (remove the scaling inside this step, make it a simple para2mod)
*/
//class paracal{
//public:
//---modify, Jan 17, 2014.
vector<int> get_index(vector<double> motherlst, vector<double> kidlst){
  vector<int> index;
  vector<double>::iterator id;
  for(int i=0;i<kidlst.size();i++){
        id=find(motherlst.begin(),motherlst.end(),kidlst[i]);
        if(id==motherlst.end()){printf("#### cannot find kid %g in motherlst\n",kidlst[i]);exit(0);}
        else{
                index.push_back(id-motherlst.begin());
        }
  }
  return index;
}
//---------------------------------------
	int initpara(paradef &para)
	{
	  //para.Rnpara=0;para.Lnpara=0;
	  para.npara=0;
     	  para.L=0.;
	  para.misfit=0.;
	  para.flag=0;
	  //para.Rparameter.clear();para.Rpara0.clear();para.Rspace1.clear();
	  //para.Lparameter.clear();para.Lpara0.clear();para.Lspace1.clear();
	  para.parameter.clear();para.LoveRAparameter.clear();
	  para.para0.clear();para.LoveAZparameter.clear();para.space1.clear();
 	  /* vector<double>().swap(para.Rparameter);vector<double>().swap(para.Lparameter);
	  vector<vector<double> >().swap(para.Rspace1);vector<vector<double> >().swap(para.Lspace1);
	  vector<vector<string> >().swap(para.Rpara0);vector<vector<string> >().swap(para.Lpara0);
	*/
	  return 1;
	}
//-----------------------------------------------------	 
	int get_para_numbers(paradef para, int ig, int &countHS, int &countCOS, int &countsigma, int &countISOsigma, int &countANIsigma){
	//this counts the number of different type of parameters.
	// the para that describes a hexagonally symmetric medium (vsv-eta,theta,phi); 
	// the para that described the azimuthal variation of Vsv (Vsvcos, Vsvsin);
	// the para that will need partial derivative computation.
	// the para that use isotropic scaling: vsh=vsv, vph=vpv=vpvs*vsv
	// the para that use anisotropic scaling: based on anisotropy, e.g., VsRA
	  
	  int i,intsigma;
	  int p6,p4;
 	  double sigma;
	  countHS=0;countCOS=0;countsigma=0;countISOsigma=0;countANIsigma=0;

	  for(i=0;i<para.npara;i++){
	    p6=(int)para.para0[i][6];//flag of parameters 1--vsv ...
	    sigma=para.space1[i][2];
	    p4=(int)para.para0[i][4];//the groups number this para belongs to
	    if(p4!=ig)continue;

	    if(p6<1 or p6>11){
		printf("#### checkParaModel/get_para_numbers, the %d th parameter (starts from 0)=%d, is outside the range, should be 1-11, modify the para.in file\n",i,p6);
		exit(0);
	    }
		
	    if(p6<8)countHS+=1;//
	    if((p6-10)*(p6-11)==0)countCOS+=1;
	    if(p6<6){
		if(sigma*(sigma+0.5)>0 )countsigma+=1;//if sigma>0 or sigma<-0.5 this para will need paratial derivative computation
		if(sigma<0){
		  intsigma=(int)sigma;
		  if(intsigma==-1 or intsigma==-3)countISOsigma+=1;
		  else countANIsigma+=1;
		}
	    }
 	  }//for npara	
	  return 1;
	}
//-----------------------------------------------------	 
	int checkParaModel(paradef para, modeldef model, vector<int> Viso){
	//check if the para satisfies certain criteria
	  int i,j;
	  int p0,flagcpt,gflag,p4,p5,p6,np,countHS,countCOS,countsigma,countISOsigma,countANIsigma;
	  int ppflagid,nvid,ngid;
	  int k2,k3,k4,k5,k6,k7;
	  vector<int>::iterator id;

	  if(para.flag<1){
	  	printf("Use checkParaModel after the para.space1 is filled, i.e., after mod2para\n");
	  	exit(0);
	  }//if
	
	  //---criterion 
	  //---modify, Jan 17, 2014 -- require the perlst of  AziampRdisp&AziphiRdisp to be a subset of the perlst of Rdisp, see the Jan 17, 2014 debug record on green:~/NOTE/code_debugging for detail
	  //just check if the get_index will go through, if yes, then it's fine
	  vector<int> index;
	  index=get_index(model.data.Rdisp.pper,model.data.AziampRdisp.pper);
	  index=get_index(model.data.Rdisp.pper,model.data.AziphiRdisp.pper);
	  index=get_index(model.data.Ldisp.pper,model.data.AziampLdisp.pper);
	  index=get_index(model.data.Ldisp.pper,model.data.AziphiLdisp.pper);

	  for(i=0;i<model.ngroup;i++){
	    flagcpt=model.groups[i].flagcpttype;// tells the forward computation type;
	    gflag = model.groups[i].flag; //tells the type of the paramters (e.g., lay, gradient, BSpline, grid ...)
	    np=model.groups[i].np;	
	    k2=2*np;k3=3*np;k4=4*np;k5=5*np;k6=6*np;k7=7*np;

	    //-----------------------------------------------------------------------------------------
	    //---criterion 
	    if(flagcpt==1){//do forward computation using Vkernel & Vpara. model is isotropic or TI.
		//---m.gp.flag !=3 and !=6
		if((gflag-3)*(gflag-6)==0){
			printf("#### checkParaModel, group%d, iso/TI model, do not suggest using Bspline that will becomed grids, or use grid model directly. Although this kind of parameterization is acceptable, at this moment, we don't prefer using it.\n",i);
			exit(0);
		}	
	
		//---#of parameters 
		get_para_numbers(para,i,countHS,countCOS,countsigma,countISOsigma,countANIsigma);		
		if(countHS+countCOS+countsigma+countISOsigma+countANIsigma==0){
			printf("group%d has no para to be purterbed.\n",i);
			continue;
		}
	    	id=find(Viso.begin(),Viso.end(),i);
		if(countCOS!=0){printf("#### checkParaModel, group%d, iso or TI model, # of VsvCOSpara =%d, should be %d, modify para.in\n",i,countCOS,0);exit(0);}
	    	if(id!=Viso.end()){//this is an isotropic group, HSpara=4 (vsv,vsh,vpv,vph); isotropic scaling vsh,vpv,vph or vsh, vph 
			if(countHS!=k4){printf("#### checkParaModel, group%d, iso model, # of HSpara =%d, should be %d, modify para.in\n",i,countHS,k4);exit(0);}
			if(countsigma!=k4  ){printf("#### checkParaModel, group%d, iso model, # of para that computes partial derivatives =%d, should be %d, modify para.in\n",i,countsigma,k4);exit(0);}
			if(countISOsigma!=k3 and countISOsigma!=k2){printf("#### checkParaModel, group%d, iso model, # of HSpara that uses isotropic scaling =%d, should be %d, modify para.in\n",i,countISOsigma,k3);exit(0);}

	    	}
	    	else{// this is a TI group. HSpara=5, 5 para computes partial derivatives
                        if(countHS!=k5){printf("#### checkParaModel, group%d, TI model, # of HSpara =%d, should be %d, modify para.in\n",i,countHS,k5);exit(0);}
                        if(countsigma!=k4 and countsigma!=k5  ){printf("#### checkParaModel, group%d, TI model, # of para that computes partial derivatives =%d, should be %d or %d, modify para.in\n",i,countsigma,k4,k5);exit(0);}
		}//else	Viso
	    }//if flagcpt==1

	    //-----------------------------------------------------------------------------------------
	    //---criterion 
	    else if (flagcpt==2 or flagcpt==4){
		//--if the type of parameters are used properly
		if((gflag-1)*(gflag-2)*(gflag-4)*(gflag-5)==0){
			printf("#### checkParaModel, group%d, Vkernel-->RAdisp Lkernel-->AZdisp (if flagcpt=2); OR Lkernel--> RAdisp&AZdisp (if flagcpt=4). Can only have m.gp.flag=3 or 6 (Bspline to grid, or grid), but here m.gp.flag=%d, change mod.in\n",i,gflag);
			exit(0);
		}
		//--if there are isotropic groups
		id=find(Viso.begin(),Viso.end(),i);
		if (id!=Viso.end()){printf("#### checkParaModel, group%d, anisotropic, whose m.gp.flagcpttype=%d, cannot be isotropic. for isotropic, set flagcpttype=1 or 3\n",i,flagcpt);exit(0);}//has isotropic group

		//---# of parameters
		get_para_numbers(para,i,countHS,countCOS,countsigma,countISOsigma,countANIsigma);
		if(countHS+countCOS+countsigma+countISOsigma+countANIsigma==0){
			printf("group%d has no para to be purterbed.\n",i);
			continue;
		}
		if(countCOS!=0){printf("#### checkParaModel, group%d, anisotropic, whose m.gp.flagcpttype=%d, # of VsvCOSpara =%d, should be %d, modify para.in\n",i,flagcpt,countCOS,0);exit(0);}
		if(countHS!=k7){printf("#### checkParaModel, group%d, m.gp.flagcpttype=%d, (uses Lkernel entirely(flagcpt=4) or partially(flagcpt=2)), # of HSpara =%d, should be %d, modify para.in\n",i,flagcpt,countHS,k7);exit(0);}
		if(countsigma!=k5  ){printf("#### checkParaModel, group%d, m.gp.flagcpttype=%d, (uses Lkernel entirely(flagcpt=4) or partially(flagcpt=2)), # of para that computes partial derivatives =%d, should be %d, modify para.in\n",i,flagcpt,countsigma,k5);exit(0);}		
	    }// else if flagcpt=2 or 4

	    //-----------------------------------------------------------------------------------------
	    //---criterion 
	    else if (flagcpt==3){
		//--if the type of parameters are used properly
		if((gflag-3)*(gflag-5)*(gflag-6)==0){printf("#### checkParaModel, group%d,flagcpt=3, Vkernel->RAdiap, Vsvkernel->AZdisp. m.gp.flag cannot be %d (water or grid model.). change mod.in\n",i,gflag);exit(0);}
		//---# of parameters
		get_para_numbers(para,i,countHS,countCOS,countsigma,countISOsigma,countANIsigma);
		if(countHS+countCOS+countsigma+countISOsigma+countANIsigma==0){
			printf("group%d has no para to be purterbed.\n",i);
			continue;
		}
		if(countCOS!=k2){printf("#### checkParaModel, group%d, m.gp.flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. # of cos/sin paramters should be %d, but here is %d, modify para.in)\n",i,k2,countCOS);exit(0);}
		id=find(Viso.begin(),Viso.end(),i);
		if(id!=Viso.end()){//isotropic group + Vsv azimuthal anisotropy
			if(countHS!=k4){printf("#### checkParaModel, group%d, flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. isotropic group, # HSpara should be %d, but not %d modify para.in\n",i,k4,countHS);exit(0);}
			if(countsigma!=k4  ){printf("#### checkParaModel, group%d, flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. isotropic group, # of non-zero sigma should be %d, but not %d modify para.in\n",i,k4,countISOsigma);exit(0);}
			if(countISOsigma!=k3 and countISOsigma!=k2){printf("#### checkParaModel, group%d, flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. isotropic group, # of ISOsigma(=-1) should be %d, but not %d modify para.in\n",i,k3,countISOsigma);exit(0);}
		}
		else{// anisotropic group + Vsv azimuthal anisotropy
			if(countHS!=k5){printf("#### checkParaModel, group%d, flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. anisotropic group, # HSpara should be %d, but not %d modify para.in\n",i,k5,countHS);exit(0);}
			if(countsigma!=k5  ){printf("#### checkParaModel, group%d, flagcpttype=3, (Vkernel->RAdiap, Vsvkernel->AZdisp. anisotropic group, # of non-zero sigma should be %d, but not %d modify para.in\n",i,k5,countISOsigma);exit(0);}
		}		

	    }//else if flagcpt=3
	    else{printf("#### checkParaModel, group%d, unconsidered flagcpttype! %d should be 1-4, modify mod.in\n",i,flagcpt);exit(0);}	

	    //-----------------------------------------------------------------------------------------
	  }//for ngroup 
	
	  //---check the order of parameters ------------
	  ppflagid=-1;nvid=-1;
	  ngid=-1;
	  for(j=0;j<para.npara;j++){
		
		p0=(int)para.para0[j][0];
		p5=(int)para.para0[j][5];//nv
		p6=(int)para.para0[j][6];//flag of para
		p4=(int)para.para0[j][4];//ng 
		

		if(p0!=0)continue;//if para is not value/Bspline

		if(p4!=ngid){nvid=-1;ngid=p4;}//come to next group
	    	//---criterion, check the order of the layer
		if(p5>nvid)ppflagid=-1;//reach the next layer in this group
		else if (p5<nvid){
			printf("### checkParaModel, in group %d, we require the layer number nv (teh 6th colunm) is increasing!\n",(int)para.para0[j][4]);
			exit(0);
		}//if	
		nvid=p5;
		//---- check the order of parameters within each layer
		if (p6<=ppflagid){
	      	  printf("### checkParaModel, in group %d, layer %d,  we require the para's ppflag (the 7th column) is ordered in an increasing order!\n",(int)para.para0[j][4],nvid);
		  exit(0);
	        }//if
	        ppflagid=p6;

	  }//for npara

	  //---criterion 
	  return 1;
	}//checkParaModel

//-----------------------------------------------------	 
//	0		1	2	3	4   5	6	7		8			9		10	
//	flag		dv_type	dv	sigma	ng  nv	ppflag	LVflag		RW_flag			LW_flag		AZ_flag
//	0--value	1--dv				1-vsv	1--use Vkernel	1--for RAdisp cpt	0--not for	0/2/4 	
//	1--gp_thick	else dv*100			2-vsh 	to do forward	2--for Azidisp cpt	   L wave	indicating cos(X*phi)
//	-1--vpvs					3-vpv	computation	0--not for Rayleigh	1/-1--is for
//							4-vph			   wave
//							5-eta	-1--use Lkernel
//							6-theta
//							7-phi
//							8-rho
//							9-h
//							10-dVsv_cos
//							11-dVsv_sin
//							12-vpvs
//
	//int readpara(vector<vector<string> > &para0, int &npara,const char* fname)
	int readpara(paradef &para, const char* fname)
	{  
	  fstream mff;
	  int k,i;
	  string line;
	  vector<string> v;
	  vector<double> vd;


	  //if (surflag>0){
	  mff.open(fname);
	  if ( not mff.is_open()){cout<<"########para file"<<fname<<"does not exist!\n exit!!\n";exit (0);}
	  k=0;
	  para.para0.clear();
	  while(getline(mff,line))
	  {
		v.clear();
		Split(line,v," ");
		//v.push_back("0"); //add a column indicating if this parameter belongs to vsv~eta and belongs to a groups whose flagcpttype==2 or ==4 (part or all of foward cpt is through the Love parameters); initial value is '0', modified in the checkParaModel function

		vd.clear();
		for(i=0;i<v.size();i++){
		  vd.push_back(atof(v[i].c_str()));
		}//for i
		if((int)vd[0]==0 and v.size()<7){
		  printf("#### readpara, wrong column #(should be >=7,1 flag,2 dv_type,3 dv,4 sigma,5 ng,6 nv,7 para_type,(8 LVflag)) in the line %d\n",k+1);
		  exit(0);
		}
		//--make the width of para0 to be 11; 11 flags for each para, see INITstructure.h for detail
		for(i=v.size();i<11;i++){
		  vd.push_back(0.);
		}
		if((int)vd[0]==1)vd[6]=9.;//thickness; vd[6] is the ppflag, indicating the type of parameter
		else if((int)vd[0]==-1)vd[6]=12.;//vpvs
		//printf("@@@ check, readpara, width of of para0 is %d, should be 11!\n",vd.size());
		k=k+1;
		para.para0.push_back(vd);	
	  }
	  mff.close();
	  para.npara=k;
	  //}//if surflag		
	  return 1;
	}//readpara
//----------------------------------------------------

//-----------------------------------------------------	 
/*
the inpara tells us what para need to be purterbed.
inpara.para0[i][j], i=npara.
outpara.para0=inpara.para0

j:	0		1		2		3		4		5
    0-value/Bcoeff	1-[i][2]=dv     dv or dv*100	sigma		group id	value id
    1-gp thickness	else...=dv*100  based on [i][1]
   -1- vpvs											  model
	==								==		==      ====>>>> outpara.parameter[] (value/thik/vpvs), length=npara=#of lines in file para.in(=>inpara.para0)
													 outpara.parameter.push_back(model.groups[ng].velue[nv]/thick/vpvs)
			-------------------------------------					====>>>> outpara.space1[min,max,sigma], length=npara
*/

	int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
	{
	  int i,ng,nv,p0,p1,p3,p4,p5,p6,p7;
	  double tv,tmin,tmax,sigma,p2;
	  vector<double> vt(3,0.);

	  outpara=inpara;
	  outpara.parameter.clear();
	  outpara.L=model.data.L;
	  outpara.misfit=model.data.misfit;


	  for(i=0;i<inpara.npara;i++){
	  	p0=(int)inpara.para0[i][0]; // flag; 0--value/Bcoeff 1--gp thickness -1--vpvs
	  	p1=(int)inpara.para0[i][1]; // type of perturbation; 1--abs dv; else--percentation
	  	p2=(float)inpara.para0[i][2]; // perturbatin; dv
		ng=(int)inpara.para0[i][4]; // grounp number; 
		p6=(int)inpara.para0[i][6]; // ppflag; 1--vsv 2--vsh 3--vpv 4--vph 5--eta 6--theta 7--phi  8-rho; 9-h; 10-dVsvcos; 11-dVsvsin; 12-vpvs;

		if(p0==0){
		  p5=(int)inpara.para0[i][5];
		  nv=p5;
		  if(p6==1)
		  	tv=model.groups[ng].vsvvalue[nv];
		  else if (p6==2)
		  	tv=model.groups[ng].vshvalue[nv];
		  else if (p6==3)
		  	tv=model.groups[ng].vpvvalue[nv];
		  else if (p6==4)
		  	tv=model.groups[ng].vphvalue[nv];
		  else if (p6==5)
		  	tv=model.groups[ng].etavalue[nv];
		  else if (p6==6 or p6==10)
		  	tv=model.groups[ng].thetavalue[nv];
		  else if (p6==7 or p6==11)
		  	tv=model.groups[ng].phivalue[nv];
		  else if (p6==8)
			tv=model.groups[ng].rhovalue[nv];
		  else {	
			printf("### mod2para, wrong number for the 7th column in para.in! ppflag=%d\n",p6);
			exit(0);}
		}// p0==0
		else if (p0==1){//and p6==9
			tv=model.groups[ng].thick;
		}
		else if (p0==-1){// and p6==12
			tv=model.groups[ng].vpvs;
		}
		else{
			printf("#### mod2para, wrong number for the 1st column in para.in!;");
			exit(0);
		}
		outpara.parameter.push_back(tv);

		if(inpara.flag<1){
			//---fill the para.space1
			sigma=outpara.para0[i][3]; // if don't want to perturb this parameter, then set sigma=0.
			if(p1==1){
				tmin=tv-p2;tmax=tv+p2;
				if(fabs(tmin-tmax)<1E-5 and sigma>1E-5){
					printf("### inproper para.in &mod set up for the %dth parameter(starts from 0), min,max=[%g,%g] RESET IT\n",i,tmin,tmax);
					exit(0);
				}
			}
			else{
				tmin=tv*(1-p2/100.);tmax=tv*(1+p2/100.);
				if(fabs(tmin-tmax)<1E-5 and sigma>1E-5){
					printf("### inproper para.in &mod set up for the %dth parameter(starts from 0), min,max=[%g,%g] RESET IT\n",i,tmin,tmax);
					exit(0);
				}
			}
			if(p6<=9 or p6==12){ // velocity or eta or theta or phi or rho or h or vpvs, should >0
				tmin=max(0.,tmin);
				tmax=max(0.,tmax);
				tmax=max(tmin+0.001,tmax);
				if (p6==5 ){//modified Mar11, 2014
					//if(ng==1)tmax=min(0.8,tmax);
					//else
					tmax=min(1.1,tmax);
				}
				// there was a constraint on the sedimental velocity, don't know why need it, so it's not added here;
			}
			vt[0]=tmin;vt[1]=tmax;vt[2]=sigma;
			if(tmin>tmax){ //modified Mar 12, 2014
				printf("### in mod2para, fill model space part. for para%d, the tmin(%g)<tmax(%g)!! something wrong, reset the para.in\n",i,tmin,tmax);
				exit(0);
			}
			//---test---
			//if(p6==7){vt[0]=90.;vt[1]=180.;}
			//if(p6==6){vt[0]=60.;vt[1]=90.;}
			outpara.space1.push_back(vt);
			printf("@@@ npara%d, v=%.3f [%.3f,%.3f] sigma=%.3f p0=%d p6=%d\n",i,tv,vt[0],vt[1],sigma,p0,p6);
		}//inpara.flag<1
	  }//for i<npara	  


	  if(inpara.flag<1){ //initial flag
		for(i=0;i<outpara.npara;i++){
			ng=(int)inpara.para0[i][4];
			p6=(int)inpara.para0[i][6];
			//---fill the flags in para0
			if(model.groups[ng].flagcpttype==2 or model.groups[ng].flagcpttype==4){
		  	  outpara.para0[i][7]=-1;//the LVflag, 1--use Vpara&Vkernel to do forward computation; -1--use Lovepara&Lovekernel to do forward computation
			}
			else{
			  outpara.para0[i][7]=1;
			}
			//printf("@@@ check, mod2para, the %dth LVflag=%g ng=%d,flagcpttype=%d \n",i,outpara.para0[i][7],ng,model.groups[ng].flagcpttype);
			if(p6==1){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=-1;
			  outpara.para0[i][10]=2;
			}
			else if(p6==2){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=-1;
			  outpara.para0[i][10]=4;
			}
			else if(p6==3){
			  outpara.para0[i][8]=1;
			}
			else if(p6==4  or p6==5){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][10]=2;
			}
			else if(p6==9){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=1;
			}
			else if(p6==10 or p6==11){
			  outpara.para0[i][8]=2;
			  outpara.para0[i][10]=2;
			}


	  	}// if inpara.flag<1
	  
	  }// for i
	  outpara.flag=1;
	
	  return 1;
	}// mod2para

//-----------------------------------------------------	 
/*
combine information in para.para0[] and para.parameter[], we can fill model.groups.value[]/thick/vpvs
for i<para.npara
        para.para0[i][j], var=para.parameter[i]
        j:      0                       4                       5
                flag                    group id
                0-value/Bcoeff  ----------------------------->  value id,nv     ===> model.groups[ng].value[nv]=para.parameter[i]
                1-gp thickness                                                  ===> model.groups[ng].thick=para.parameter[i]
                -1- vpvs                                                        ===> model.groups[ng].vpvs=para.parameter[i]


*/
	int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)
	{
	  //this is different from the para2mod defined below this function; it won't modify the para anymore, it just use the given para, and transfer it to model, regardless if space1[i][2]<-1 or not;
	  // this is used in computing Vkernel!	
	  int i, p0,p1,ng,p6,nv,pflag;
 	  double newv,dh;
	  int np,flagvpv,flagvph,flagvsv,flagvsh;
	  flagvpv=flagvph=flagvsv=flagvsh=-1;
	  //float factor;

	  outmodel=inmodel;
	  for(i=0;i<para.npara;i++){
	  	p0=(int)para.para0[i][0];//flag, explanations are in mod2para;
	  	p1=(int)para.para0[i][1];
	  	ng=(int)para.para0[i][4];
	  	p6=(int)para.para0[i][6];//ppflag, explanations are in mod2para;
		
		newv=para.parameter[i];

		if(p0==0){ //value
		  nv=(int)para.para0[i][5];
		  if(p6==1)
		  	  {outmodel.groups[ng].vsvvalue[nv]=newv;flagvsv=1;}
		  else if (p6==2)
		  	  {outmodel.groups[ng].vshvalue[nv]=newv;flagvsh=1;}
		  else if (p6==3){//vpv
			  outmodel.groups[ng].vpvvalue[nv]=newv;flagvpv=1;}
		  else if (p6==4){//vph
			  outmodel.groups[ng].vphvalue[nv]=newv;flagvph=1;}
		  else if (p6==5){//eta
			   outmodel.groups[ng].etavalue[nv]=newv;}
		  else if (p6==6 or p6==10){//theta
			  outmodel.groups[ng].thetavalue[nv]=newv;
		  	  }//p6=6 or 10	
		  else if (p6==7 or p6==11){//phi
			  outmodel.groups[ng].phivalue[nv]=newv;
		  	  }//if p6==7 or 11
		  else if (p6==8){//rho
			  outmodel.groups[ng].rhovalue[nv]=newv;
			  }
		  else{printf("### para2mod, the %dth para, from gp%d nv%d, has p0=0(indicating it's a value other than h or vpvs), p6=%d, which is unrecognized\n",i,ng,nv,p6);exit(0);}
		}//if p0==0

		else if(p0==1){// groups thickness
		//---should only change the depth of one gp, keep the DEPTH of others grups unchanged (so the thickness of the following group will chnage)
                  if(ng==outmodel.ngroup-1){// the last group, so only change the thickness of this group and the total thickness
                  // actually it is kind of meaningless to perturb the last group's thickness ...
                        dh=newv-outmodel.groups[ng].thick;
			//printf("group%d is the last\n, thcik %g->%g dh=%g tthick%g->%g",ng,outmodel.groups[ng].thick,newv,dh,outmodel.tthick,outmodel.tthick+dh);
                        outmodel.tthick+=dh;
                        outmodel.groups[ng].thick+=dh;
			//printf("==> thick[%d]=%g tthick=%g\n",ng,outmodel.groups[ng].thick,outmodel.tthick);
                  }
                  else{// not the last group, then change the thickness of this group, and that of the following group, total thickness will not change
                        dh=newv-outmodel.groups[ng].thick;
			//printf("group%d is not the last\n, thcik %g->%g dh=%g \n",ng,outmodel.groups[ng].thick,newv,dh);
                        outmodel.groups[ng].thick+=dh;
                        outmodel.groups[ng+1].thick-=dh;
			//printf("==> thick[%d]=%g thick[%d]=%g thick[%d]=%g, sum=%g, tthick=%g\n",ng,outmodel.groups[ng].thick,ng+1,outmodel.groups[ng+1].thick,ng+2,outmodel.groups[ng+2].thick,outmodel.groups[ng].thick+outmodel.groups[ng+1].thick+outmodel.groups[ng+2].thick,outmodel.tthick);
                  }
		  /*/------------------test----------
		  for(int k=0;k<outmodel.ngroup;k++){
			if(outmodel.groups[k].thick<0){
			printf("gp[%d].thick=%g<0!",k,outmodel.groups[k].thick);
			exit(0);
			}
		  }
		  //-----------
		  */

		}// else if p0==1
		
		else if (p0==-1){//vpvs
			outmodel.groups[ng].vpvs=newv;
		}
		else {printf("### para2mod, wrong value for p0, should be either -1,1, or0, not %d\n",p0);exit(0);}

	  }// for i


	  return 1;
	}// para2mod
	 


//-----------------------------------------------------	 

	int para2mod_not_in_use(paradef &para, modeldef inmodel, modeldef &outmodel)
	{
	  //for a parameteri, if its space[i][2] (i.e., the sigma) <-1,then, use the default RAvp, RAvs, eta relationship to scale vpv,vph,eta according to vsv,vsh;
	  //ALSO, since I add & in front of para, the para will be changed accordingly (change to the scaled value, so don't need mod2para after para2mod in order to make sure para and mod have the same value)
	  //RAvs=(vsh-vsv)/((vsv+vsh)/2)
	  //vpv=vsv*vpvs
	  //A=0.25*RAvs
	  //vph=(1+A)/(1-A)*vpv
	  //eta=1-4.2*RAvs 	
	  //
	  //another way to scale eta, make the media to have purely ellipsoidal anisotropy
	  //eta=0.5*(1+(vpv^2-2*vsv^2)/(vph^2-2*vsv^2))
	  //
	  //for vsh, if the sigma<-1; then scale it according to the vsv in the current layer & the RAvs value in the 1st layer of that group
	  //set them to have the same VsRA (now I use the simplified computation in which RAvs=(vsh-vsv)/(vsv+vsh)/2)
	  //
	  //
	  //if scaling relationship is not linear, then the scaling can only be applied to layered,gradient,point model, not Bsp model!; BS
	  //so, for Bsp model, I only apply the linear scaling relationship (vsv->vpv) while leave vph, eta(both have non-linear scaling relationship) unscaled.
	  int i, p0,p1,ng,p6,nv,pflag,gflag;
 	  double newv;
	  int np,flagvpv,flagvph,flagvsv,flagvsh;
	  flagvpv=flagvph=flagvsv=flagvsh=-1;
	  float factor,A,C,L,F,cc,eta;

	  outmodel=inmodel;
	  for(i=0;i<para.npara;i++){
	  	p0=(int)para.para0[i][0];//flag, explanations are in mod2para;
	  	p1=(int)para.para0[i][1];
	  	ng=(int)para.para0[i][4];
	  	p6=(int)para.para0[i][6];//ppflag, explanations are in mod2para;
		gflag=inmodel.groups[ng].flag; //indicating the type of values within this group; gradient, layer, Bsp, point etc.		

		newv=para.parameter[i];

		if(p0==0){ //value
		  nv=(int)para.para0[i][5];
		  if(p6==1){
				outmodel.groups[ng].vsvvalue[nv]=newv;flagvsv=1;}
		  else if (p6==2){//vsh
			  if(para.space1[i][2]<-1 ){// linear scale; scale it according to the VsRA in the layer 1 of this gorup 
				factor=outmodel.groups[ng].vshvalue[0]/outmodel.groups[ng].vsvvalue[0];
				outmodel.groups[ng].vshvalue[nv]=factor*outmodel.groups[ng].vsvvalue[nv];
				para.parameter[i]=outmodel.groups[ng].vshvalue[nv];	
			  }
	       		  else		  
			  	outmodel.groups[ng].vshvalue[nv]=newv;
			  flagvsh=1;}
		  else if (p6==3){//vpv ; 
			  if(para.space1[i][2]<-1.)//linear scale; use the RAvp, eta, RAvs relation 
				{outmodel.groups[ng].vpvvalue[nv]=outmodel.groups[ng].vpvs*outmodel.groups[ng].vsvvalue[nv];
				 para.parameter[i]=outmodel.groups[ng].vpvvalue[nv];
				}
			  else	  
			   	{outmodel.groups[ng].vpvvalue[nv]=newv;}
			  flagvpv=1;
			  }//p6==3

		  else if (p6==4){//vph ; BS
			  // if sigma<-1; then 1) if non-Bspline crust, then use RAvp-RAvs relation to get vph; 2) if mantle or Bspline crust, then vph=vpv ATTENTION:sometimes I allow scaling in mantle by rm the ng<2 criteria
			  if(para.space1[i][2]<-1.)// non-linear scale; use the RAvp, eta, RAvs relation ; BS
			  {
			      	//if((gflag-2)*(gflag-3)!=0 and ng<2){//non-linear scaling for non-Bsp crust(sed+crystalline crust); since mantle's scaling realtion is unknown, don;t apply this scaling to mantle(ng>=2)    ; BS
			      	if((gflag-2)*(gflag-3)!=0){ //---check---non-Bspline layer
					if(para.space1[i][2]<-3.){//use the RAvp-RAvs relation in the 1st layer of this group
						factor=(outmodel.groups[ng].vphvalue[0]-outmodel.groups[ng].vpvvalue[0])/(outmodel.groups[ng].vphvalue[0]+outmodel.groups[ng].vpvvalue[0])/(outmodel.groups[ng].vshvalue[0]-outmodel.groups[ng].vsvvalue[0])*(outmodel.groups[ng].vshvalue[0]+outmodel.groups[ng].vsvvalue[0]) ; //RAvp0/RAvs0
						if(isnan(factor))factor=0.;
						//printf("check, vph, sigma<-3.factor=%g====\n",factor);//--check--
						factor=factor*(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv]);
					}//if <-3
					else{//use a given RAvp-RAvs relationship
						factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv])*0.5; //0.25*RAvs (VpRA=0.5VsRA)
					}//else >-3
					//factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv]); //0.5*RAvs (VpRA=VsRA)
					factor=(1+factor)/(1-factor);	  	
					if(isnan(factor))factor=1.;
					outmodel.groups[ng].vphvalue[nv]=factor*outmodel.groups[ng].vpvvalue[nv];}//if (gflag-2)*(gflag-3)
				else{outmodel.groups[ng].vphvalue[nv]=outmodel.groups[ng].vpvvalue[nv];}
				para.parameter[i]=outmodel.groups[ng].vphvalue[nv];
			   }
			   else	
			  	{outmodel.groups[ng].vphvalue[nv]=newv;} 
			   flagvph=1;	
		  	   }//p6==4

		  else if (p6==5){//eta ; BS
			   if(para.space1[i][2]<-1. )//non-linear scale; use the RAvp, eta, RAvs relation ; BS
			   {
				  //if(ng<2 and (gflag-2)*(gflag-3)!=0){ //; since mantle's scaling realtion is unknown, don;t apply this scaling to mantle(ng>=2)
				  if( (gflag-2)*(gflag-3)!=0){//---check--- non-Bspline layer
				  //if(1){//-----check---- Apr 20, 2014.
				
					if(para.space1[i][2]<-5.){//constant c (non-ellipticity) modified Mar 12, 2014
						A=pow(outmodel.groups[ng].vphvalue[0],2);//A/rho
						C=pow(outmodel.groups[ng].vpvvalue[0],2);//C/rho
						L=pow(outmodel.groups[ng].vsvvalue[0],2);//L/rho
						F=outmodel.groups[ng].etavalue[0]*(A-2*L);//F/rho
						cc=0.125*(A+C-2*F-4*L); //c/rho=0.125*(A+C-2F-4L)/rho

	                                        A=pow(outmodel.groups[ng].vphvalue[nv],2);//A/rho
                                                C=pow(outmodel.groups[ng].vpvvalue[nv],2);//C/rho
                                                L=pow(outmodel.groups[ng].vsvvalue[nv],2);//L/rho

						eta=(A+C-4*L-8*cc)*0.5/(A-2*L);//
						//printf("ng=%d nv=%d c=%g,eta=%g\n",ng,nv,cc*outmodel.groups[ng].rhovalue[0],eta);//---check---
						outmodel.groups[ng].etavalue[nv]=eta;
	
					}
					else if(para.space1[i][2]<-3.){//constant eta/RAvs
						factor=(outmodel.groups[ng].vshvalue[0]-outmodel.groups[ng].vsvvalue[0])/(outmodel.groups[ng].vshvalue[0]+outmodel.groups[ng].vsvvalue[0]);//  (vsh-vsv)/(vsh+vsv)=vsra0
						factor=(1.-outmodel.groups[ng].etavalue[0])/factor; // (1-eta)/vsra0
						if(isnan(factor) or factor>1E20)factor=0.;
						outmodel.groups[ng].etavalue[nv]=1.-factor*(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv]);//1-c0*vsra
						//printf("ng=%d nv=%d factor=%g eta=%g\n",ng,nv,factor,outmodel.groups[ng].etavalue[nv]);//---check---
					}//if < -3
				   	else{
				   	  //factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv])*2;//RAvs
					  //outmodel.groups[ng].etavalue[nv]=1.0-4.2*factor;
					  //c==0, elliptical case
 				  		factor=(pow(outmodel.groups[ng].vpvvalue[nv],2)-2*pow(outmodel.groups[ng].vsvvalue[nv],2))/(pow(outmodel.groups[ng].vphvalue[nv],2)-2*pow(outmodel.groups[ng].vsvvalue[nv],2));
				  		outmodel.groups[ng].etavalue[nv]=0.5*(1+factor);}//else >=-3
				  }
				  else{outmodel.groups[ng].etavalue[nv]=1.0;}
   				  para.parameter[i]=outmodel.groups[ng].etavalue[nv];	
			   }
		  	   else
			  	{outmodel.groups[ng].etavalue[nv]=newv;
				
				}
			   }//p6==5
		  else if (p6==6 or p6==10){//theta or cos
			  if(para.space1[i][2]<-1)//set the theta to the value of the 1st layer in that group
			  {
			  	outmodel.groups[ng].thetavalue[nv]=outmodel.groups[ng].thetavalue[0];
				para.parameter[i]=outmodel.groups[ng].thetavalue[nv];
			  }
			  else
			  	{outmodel.groups[ng].thetavalue[nv]=newv;}
		  	  }//p6=6 or 10	
		  else if (p6==7 or p6==11){//phi or sin
			  if(para.space1[i][2]<-1)
			  {
			  	outmodel.groups[ng].phivalue[nv]=outmodel.groups[ng].phivalue[0];
				para.parameter[i]=outmodel.groups[ng].phivalue[nv];
			  }	  
			  else
			  	{outmodel.groups[ng].phivalue[nv]=newv;}
		  	  }//if p6==7 or 11
		}//if p0==0

		else if(p0==1){// groups thickness
			outmodel.tthick=outmodel.tthick-outmodel.groups[ng].thick+newv; // the tthick is changing
			//outmodel.groups[outmodel.ngroup-1].thick=outmodel.groups[outmodel.ngroup-1].thick+(outmodel.groups[ng].thick-newv);//change the thickness of the last group while keeping the tthcik the same // find a bug, if ngroup<3, this choice would cause problem!
			outmodel.groups[ng].thick=newv;
		
		}// else if p0==1
		
		else if (p0==-1){//vpvs
			outmodel.groups[ng].vpvs=newv;
		}
		else{cout<<"########wrong!!!!! in para2mode"<<endl;exit(0);}
	  }// for i


	  //--the following value filling process have been performed in the readmodAniso function.
	  //---depending on the pflag of each group, fill all the unfilled vectors among those vsv...phi 7 vectors with appropriate values derived from certain relationships between vp&vs;
	/*
	  for(ng=0;ng<outmodel.ngroup;ng++){
		np=outmodel.groups[ng].np;
		pflag=outmodel.groups[ng].pflag;
	  	if(pflag<7){
			for(i=0;i<np;i++)
				{outmodel.groups[ng].thetavalue[i]=outmodel.groups[ng].phivalue[i]=0.;}
			
			if(pflag<5){
			  for(i=0;i<np;i++)outmodel.groups[ng].etavalue[i]=1.0;
			  
			  if(pflag==3){
			    for(i=0;i<np;i++){
			      if(flagvpv>0)outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i];
			      else if(flagvph>0)outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vphvalue[i];
			    }
			  }//pflag==3
			  
			  else if(pflag==2){
			      for(i=0;i<np;i++){
				      outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vsvvalue[i]*outmodel.groups[ng].vpvs;

			      }  
			  }//pflag==2

			  else if(pflag==1){
			  	if(flagvsv>0)outmodel.groups[ng].vshvalue[i]=outmodel.groups[ng].vsvvalue[i];
				else if(flagvsh>0)outmodel.groups[ng].vsvvalue[i]=outmodel.groups[ng].vshvalue[i];
				outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vsvvalue[i]*outmodel.groups[ng].vpvs;
			  }//pflag==1 

			}//pflag<5		
		} // if pflag<7	  
	  } // for ng
	*/
	  return 1;
	}// para2mod
	 

/*
//-----------------------------------------------------	 
// randomly generating new parameters, uniform or normal distribution
	int gen_newpara_single (  vector<vector<double> > space1, vector<double> &parameter, int npara,int pflag)
	//vector<double> gen_newpara_single ( const vector<vector<double> > &space1, vector<double> parameter, int npara,int pflag)
	{
	  int i,flag,nx;
	  double newv,sigma,mean,perc;

	  //cout<<"gen_newpara! pflag="<<pflag<<endl;
	  
	  if(pflag==0)//uniform
	  {
		for(i=0;i<npara;i++)
			{
			  if(space1[i][2]<0.001)continue;
			  newv=gen_random_unif01()*(space1[i][1]-space1[i][0])+space1[i][0];//gen_random_unif01 => [0,1)
			  //printf("$$$$ gen_newpara unif:ipara=%d/%d min=%g max=%g newv=%g\n",i,npara,space1[i][0],space1[i][1],newv);
			  parameter[i]=newv;
			}
	  }//if		  	  
	  else if(pflag==1)//normal
	  {
		for(i=0;i<npara;i++)
		{
		  flag=0;
		  mean=parameter[i];
		  sigma=space1[i][2];
		  if(sigma<0.001)continue;
		  //printf("test-- gen_newpara i=%d (out of %d-1) mean=%g sigma=%g space=[%g~%g]\n",i,npara,mean,sigma,space1[i][0],space1[i][1]);
		  while(flag<1)//throw away newv that's outside para range
		  {
			newv=gen_random_normal(mean,sigma); // normal distribution
			if(newv>space1[i][1] or newv<space1[i][0])continue;
			else flag=2;
		  }//while flag
		  
		  newv=gen_random_normal(mean,sigma); // normal distribution
		  // mirror reflect newv that's outside para range
		  // in some rare cases, even after relection, the value is still outside the m[min,max] range, and that may cause trouble ==> new use multiple reflection, see below:
		  //if(newv>space1[i][1]){newv=space1[i][1]-(newv-space1[i][1]);}
		  //else if(newv<space1[i][0]){newv=space1[i][0]+(space1[i][0]-newv);}
		  if(newv>space1[i][1]){//x>c2
		    nx=1;
		    newv=2*space1[i][1]-newv; //x1=2*c2-x
		    //x2=2*c1-x1;...; x[n]=2*ck-x[n-1] where k=1 if n is even, k=2 if n is odd
		    while((newv-space1[i][0])*(newv-space1[i][1])>0){
			nx+=1;
			if(nx%2==0)//n is even
				newv=2*space1[i][0]-newv;
			else //n is odd
				newv=2*space1[i][1]-newv;
		    }
		  }//if newv>space[i][1]
		  else if (newv<space1[i][0]){ //x<c1
		    nx=1;
		    newv=2*space1[i][0]-newv;//x1=2*c1-x
                    //x2=2*c2-x1;...; x[n]=2*ck-x[n-1] where k=2 if n is even, k=1 if n is odd
		    while((newv-space1[i][0])*(newv-space1[i][1])>0){
			nx+=1;
			if(nx%2==0)//n is even
				newv=2*space1[i][1]-newv;
			else //n is odd
				newv=2*space1[i][0]-newv;
		    }//while	
		  }//else if newv<space1[i][0]

		  parameter[i]=newv;
		}//fori

	  }//else if 1
	  else if(pflag==2)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
          {
                perc=gen_random_unif01();
                for(i=0;i<npara;i++)                
		{
		//cout<<"test--- i="<<i<<endl;
		  if (space1[i][2]<0.001){continue;} // skip the para don't require any purterbation
		  newv=perc*(space1[i][1]-space1[i][0])+space1[i][0];
		  //printf("f2 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		  parameter[i]=newv; 
               }

          }//else if 2
	  else if(pflag==3)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
	  {
		//printf("pglag=%d,npara=%d\n",pflag,inpara.npara);
		flag=0;
		for (i=0;i<npara;i++)
		{
		 if (space1[i][2]<0.001)continue;
		 mean=parameter[i];
		 sigma=space1[i][2];
		 break;
		}
		if (i==npara){
		 // printf("gen_newPara, pflag==3, no para need to be purterbed! skip!\n");
		  return 1;
		  //return parameter;
		}
                while(flag<1)
                  {
                        newv=gen_random_normal(mean,sigma); // normal distribution
                        if(newv>space1[i][1] or newv<space1[i][0])continue;
                        else flag=2;
                  }//while flag
		perc=newv/mean;
		for(i=0;i<npara;i++)
		{
		  if (space1[i][2]<0.0001)continue; // skip the para don't require any purterbation
		  newv=parameter[i]*perc;
	//	  printf("f3 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		  parameter[i]=newv;
		}
	  }//else if 3
	  else{printf("### wrong pflag for gen_newpara!!\n");exit(0);}
	  return 1;
	  //return parameter;
	}//gen new para single
*/
//-----------------------------------------------------

//----------------------------------------------------- 
	double gen_newpara_single_v2(vector<vector<double> > space1, vector<double> parameter,int ip,int pflag){
	  double newv,min,max,mean,sigma;
	  int nx;

	  mean=parameter[ip];
	  newv=mean;
	  min=space1[ip][0];
	  max=space1[ip][1];
	  sigma=space1[ip][2];

	  if(sigma<0.001)return 1;
	  //if(sigma<0.001){printf("Hey, this para %d, (=%g) either does not need any perturbation, or it is perturbed through scaling relationship. Should not be send to this gen_newpara_single_v2 function.\n",ip,mean);exit(0);}

	  if(pflag==0){//uniform
		newv=gen_random_unif01()*(max-min)+min;
		}
	  else if (pflag==1){//normal
		//newv=gen_random_normal(mean,sigma);
		//if(newv>max){newv=max-(newv-max);}
		//else if(newv<min){newv=min+(min-newv);}

		  newv=gen_random_normal(mean,sigma); // normal distribution
		  // mirror reflect newv that's outside para range
		  // in some rare cases, even after relection, the value is still outside the m[min,max] range, and that may cause trouble ==> new use multiple reflection, see below:
		  //if(newv>space1[i][1]){newv=space1[i][1]-(newv-space1[i][1]);}
		  //else if(newv<space1[i][0]){newv=space1[i][0]+(space1[i][0]-newv);}
		  if(newv>max){//x>c2
		    nx=1;
		    newv=2*max-newv; //x1=2*c2-x
		    //x2=2*c1-x1;...; x[n]=2*ck-x[n-1] where k=1 if n is even, k=2 if n is odd
		    while((newv-min)*(newv-max)>0){
			nx+=1;
			if(nx%2==0)//n is even
				newv=2*min-newv;
			else //n is odd
				newv=2*max-newv;
		    }
		  }//if newv>space[i][1]
		  else if (newv<min){ //x<c1
		    nx=1;
		    newv=2*min-newv;//x1=2*c1-x
                    //x2=2*c2-x1;...; x[n]=2*ck-x[n-1] where k=2 if n is even, k=1 if n is odd
		    while((newv-min)*(newv-max)>0){
			nx+=1;
			if(nx%2==0)//n is even
				newv=2*max-newv;
			else //n is odd
				newv=2*min-newv;
		    }//while	
		  }//else if newv<space1[i][0]
		  //if(newv<0){system("pause");}
	  }
	  else {printf("### wrong pflag for gen_newpara!!\n");exit(0);}
	  //removed the section that "the perturbed parameter to have the same perturbation", this is enabled through the value of sigma, fulfilled in the main gen_newpara function
	  //printf("gen_newpara %g\n",newv);//---test---
	  return newv;
	}//gen_newpara_single_v2
//----------------------------------------------------- 
	double gen_newpara_single_scale(double flag, modeldef model,int ng,int nv,int p6){
	//generate parameter based on scaling relationships
	  double newv,c,ts;
	  int intflag,gpnum;
	  intflag=(int)flag;
	  if(intflag==-1){
	  //isotropic scaling, this group is isotropic
	  //vsh=vsv;vph=vpv=vsv*vpvs
	  //actually, if para was set properly, eta or theta or phi should not appear here, but i still allow setting eta=0 theta=phi=0 in this function at this moment
		if(p6==1){printf("###Hey, Vsv should not be scaled! reset para.in\n");exit(0);}
		else if(p6==2){//vsh
			newv=model.groups[ng].vsvvalue[nv];}
		else if(p6==3){//vpv
			newv=model.groups[ng].vsvvalue[nv]*model.groups[ng].vpvs;}
		else if (p6==4){//vph
			newv=model.groups[ng].vpvvalue[nv];}
		else if (p6==5){//eta
			newv=1.0;}
		else if (p6==6 or p6==7)newv=0.;//theta or phi
		else if (p6==8){//rho; keep this part consistent with readmodAniso function
			if(ng==2){//mantle default  assume this is mantle!!! THIS NEED TO BE MODIFIED IF GROUP2 IS NOT MANTLE
			  ts=0.5*(model.groups[ng].vsvvalue[nv]+model.groups[ng].vshvalue[nv]);
			  newv=3.42+0.01*100*(ts-4.5)/4.5;
			}
			else if (model.groups[ng].flag==5){//water layer
				newv=1.02;
			}
			else{
			  ts=0.5*(model.groups[ng].vsvvalue[nv]+model.groups[ng].vshvalue[nv]);
			  newv=1.22679 + 1.53201*ts -0.83668*ts*ts + 0.20673*ts*ts*ts -0.01656*ts*ts*ts*ts;
			}
		}
		else{printf("###inproper para.in, para with p6=%d should not apprear in the isotropic scaling\n",p6);exit(0);}
	  }
	  else if (intflag==-2){
	  //anisotropic scaling, this group is anisotropic
	  //scale based on gpnum layer's value, keep constant magnitude of anisotropy throughout this whole group
	  //if flag=-2.0 then scaled based on the layer0's value
	  //if flag=-2.2 then scale based on the layer2's vaue	
	  //for AZI case, Acos=Acos[x], Asin=Asin[x]; azi has the same direction and amplitude 
		gpnum=((int)(flag*10)%10);
		if(p6==1){printf("###Hey, Vsv should not be scaled! reset para.in\n");exit(0);}
		else if (p6==2){//vsh
			c=model.groups[ng].vshvalue[gpnum]/model.groups[ng].vsvvalue[gpnum];
			newv=c*model.groups[ng].vsvvalue[nv];}
		else if (p6==3){//vpv
			newv=model.groups[ng].vsvvalue[nv]*model.groups[ng].vpvs;}
		else if (p6==4){//vph
			c=model.groups[ng].vphvalue[gpnum]/model.groups[ng].vpvvalue[gpnum];
			newv=c*model.groups[ng].vpvvalue[nv];}
		else if (p6==5){//eta
			newv=model.groups[ng].etavalue[gpnum];}
		else if (p6==6 or p6==10 ){//theta or AZcos
			newv=model.groups[ng].thetavalue[gpnum];}
		else if (p6==7 or p6==11 ){//phi or AZsin
			newv=model.groups[ng].phivalue[gpnum];}
		else{printf("###inproper para.in, para with p6=%d should not apprear in the anisotropic scaling 1\n",p6);exit(0);}
	  }
	  else if (intflag==-3){
		if(p6==1){printf("###Hey, Vsv should not be scaled! reset para.in\n");exit(0);}
               else if(p6==2){//vsh
                       newv=model.groups[ng].vsvvalue[nv];}
               else if(p6==3){//vpv
		       ts=model.groups[ng].vsvvalue[nv];
                       newv=0.9409 + 2.0947*ts - 0.8206*ts*ts + 0.2683*ts*ts*ts -0.0251*ts*ts*ts*ts;}
                else if (p6==4){//vph
                       newv=model.groups[ng].vpvvalue[nv];}
                else if (p6==5){//eta
                        newv=1.0;}
                else if (p6==6 or p6==7)newv=0.;//theta or phi
                else{printf("###inproper para.in, para with p6=%d should not apprear in the isotropic scaling\n",p6);exit(0);}
          }
	
//tp = 0.9409 + 2.0947*ts - 0.8206*ts*ts + \
//                 0.2683*ts*ts*ts -0.0251*ts*ts*ts*ts;
//                 trho = 1.22679 + 1.53201*ts -0.83668*ts*ts + 0.20673*ts*ts*ts -0.01656*ts*ts*ts*ts;


	  else if (intflag==-4){
	  //anisotropic scaling, this group is anisotropic
	  //scale based on a given scaling relationship
	  //e.g., VpRA=c1*VsRA; eta=c2+c3*VsRA	  
	  // Acos/Asin=Acos[x]/Asin[x] --> this group has the same azi direction. 
	  //this kind of scaling need to be dealt with carefully. need to be careful about the updating order of the parameters. the parameters are ordered from vsv->eta, but when doing scaling, may need some later parameters to be updated 1st.
	 	if(p6==4){//vph --> [(vph-vpv)/vpv]/[(vsh-vsv)/vsv]=c
			if(ng==0){c=1.0;}//in sediment VpRA=VsRA;
			else{c=0.5;}//in cst & mat VpRA=0.5*VsRA
			c=c*model.groups[ng].vpvvalue[nv]/model.groups[ng].vsvvalue[nv];
			newv=c*(model.groups[ng].vshvalue[nv]-model.groups[ng].vsvvalue[nv])+model.groups[ng].vpvvalue[nv];
		}
		else if (p6==5){//eta --> eta=1.0-4.2*vsRA
			c=(model.groups[ng].vshvalue[nv]-model.groups[ng].vsvvalue[nv])/model.groups[ng].vsvvalue[nv];
			newv=1.0-4.2*c;
		}
		//-------under construction -----
		else if (p6==11){// Asin --> Asin/Acos=Asin[x]/Acos[x]; theta-Acos, phi-Asin
			gpnum=((int)(flag*10)%10);
			if(fabs(model.groups[ng].thetavalue[gpnum])<1e-4){// the Acos[x]==0; then should set Acos=0 and Asin=random();
				return -999.;}

			else{// Asin=Asin[x]/Acos[x]*Acos
				newv=model.groups[ng].phivalue[gpnum]/model.groups[ng].thetavalue[gpnum]*model.groups[ng].thetavalue[nv];
				//--check the sign, if sign is flipped, then return 1k+newv
				if(newv*model.groups[ng].phivalue[gpnum]<0)return 1000.0+newv;
			}
		}
		//-----------
		else{
	  	printf("### this anisotropic scaling is still under construction...\n");
	  	exit(0);
		}
  	  }
	  else if (intflag==-5){
	  // vpvs scaling, this group can be isotropic or anisotropic
	  // scale the the vpv based on 1st layer's vsv and vpv0/vsv0
	  	if(p6==3){//vpv
			c=model.groups[ng].vpvvalue[0]/model.groups[ng].vsvvalue[0];
			newv=model.groups[ng].vsvvalue[nv]*c;
		}
		else{printf("###inproper para.in, para with p6=%d should not apprear in the vpvs scaling\n",p6);exit(0);}
	  }
	  return newv;
	}//gen_newpara_single_scale

//----------------------------------------------------- 
	int gen_newpara(paradef inpara, modeldef model, paradef &outpara, int pflag)
	{
	  //before using this function, make sure that para&model are consistent with each other
	  //once a parameter is changed, need to update the corresponding value in model because I may use this updated value for another scaled parameter.
	  //within this function, I only generate parameter, the modified model values are not transferred back. should use para2mod after this function to keep para&model consistent with each other
	  //it is easy to transfer the updated model values back by just using '&model'. But at this points, I just want to keep this function simple
	  //half of this function actually take the role of para2mod... but that part CANNOT be removed, because in the gen_newpara_single_scale function, we need to use the UPDATED model (not para) values! 
	  int i;
	  int intsigma,p0,ng,p6,nv;
	  double newv,sigma;
	  double dh;
	  modeldef tmodel;

	  outpara=inpara;
	  tmodel=model;

	  for(i=0;i<inpara.npara;i++){
		p0=(int)inpara.para0[i][0];
		ng=(int)inpara.para0[i][4];
		nv=(int)inpara.para0[i][5];
		p6=(int)inpara.para0[i][6];
		sigma=inpara.space1[i][2];
	
		if(p0==0){//value
		  if(sigma>1e-5)
			outpara.parameter[i]=gen_newpara_single_v2(outpara.space1,outpara.parameter,i,pflag);
		  else if (sigma<-1e-5){
		    //intsigma=(int)sigma;
		    outpara.parameter[i]=gen_newpara_single_scale(sigma,model,ng,nv,p6);
		    if(outpara.parameter[i]<-900. and p6==11){// AZsin, and the Acos[x] used for scaling is 0; so should set Acos=0, and Asin=random
			if((int)inpara.para0[i-1][6]!=10){
				printf("### gen_newpara, wrong para.in setting!, the paramter before Asin should be Aco!\n");
				exit(0);
			}
			outpara.parameter[i-1]=0.;
			outpara.parameter[i]=gen_newpara_single_v2(outpara.space1,outpara.parameter,i,pflag);

		    }//if Asin<-900
		    else if (outpara.parameter[i]>900){//the sign of cos and sin is not right, should multiply both of them by -1
			outpara.parameter[i-1]*=-1;
			outpara.parameter[i]=(outpara.parameter[i]-1000.0)*(-1);
		    }// else if >900
		  }//else if	
		  // then else indicate sigma==0; keep parameter unchanged

		  newv=outpara.parameter[i];
                  if(p6==1)
                          {model.groups[ng].vsvvalue[nv]=newv;}
                  else if (p6==2)
                          {model.groups[ng].vshvalue[nv]=newv;}
                  else if (p6==3){//vpv
                           model.groups[ng].vpvvalue[nv]=newv;}
                  else if (p6==4){//vph
                           model.groups[ng].vphvalue[nv]=newv;}
                  else if (p6==5){//eta
                           model.groups[ng].etavalue[nv]=newv;}
                  else if (p6==6 or p6==10){//theta or dVScos
                           model.groups[ng].thetavalue[nv]=newv;
                          }//p6=6 or 10 
                  else if (p6==7 or p6==11){//phi or dVSsin
                           model.groups[ng].phivalue[nv]=newv;
                          }//if p6==7 or 11
		  else if (p6==8){//rho
			   model.groups[ng].rhovalue[nv]=newv;
			  }
		  else{
			printf("#### wrong flag for paramter%d from gp%d nv%d, p0==0 (indicate it is some value other than thickness or vpvs), p6=%d is unrecognized\n",i,ng,nv,p6);
			exit(0);
		       }

		}//if p0==0 value
		else if(p0==1){// thickness parameter
		  outpara.parameter[i]=gen_newpara_single_v2(outpara.space1,outpara.parameter,i,pflag);		  
		  //---should only change the depth of one gp, keep the DEPTH of others grups unchanged (so the thickness of the following group will chnage)
		  newv=outpara.parameter[i];
		  if(ng==model.ngroup-1){// the last group, so only change the thickness of this group and the total thickness
		  // actually it is kind of meaningless to perturb the last group's thickness ...
			model.tthick=model.tthick-model.groups[ng].thick+newv;
			model.groups[ng].thick=newv;
		  }
		  else{// not the last group, then change the thickness of this group, and that of the following group, total thickness will not change
			dh=newv-model.groups[ng].thick;
			model.groups[ng].thick+=dh;
			model.groups[ng+1].thick-=dh;
		  }

		}//if p0=1
		else if (p0==-1)//vpvs
		  {outpara.parameter[i]=gen_newpara_single_v2(outpara.space1,outpara.parameter,i,pflag);
  		   newv=outpara.parameter[i];
		   model.groups[ng].vpvs=newv;
		  }
		else {printf("## gen_newpara, wrong value for p0, should be either -1,1, or0, not %d\n",p0);exit(0);}
	  } //for i
	  return 1;
	}//gen_newpara
	

//----------------------------------------------------- 
	int para_avg(vector<paradef> &paralst, paradef &parabest, paradef &paraavg,vector<double> &parastd,vector<double> &LoveRAparastd, vector<vector<double> > &LoveAZparastd,vector<int> &idlst){
	//get the average parameters from a list of acceptable paras
	//also, get the best parameter (min misfit)
	    int i,j,k,size,Ngood,ibadp,idmin;
	    double mismin=1e10;
	    vector<double> LAZp0(2,0.);
	    vector<int> badv;

	    initpara(paraavg);
  	    initpara(parabest);
	    size=paralst.size();
	    if(size<1) return 0;
		
	    paraavg=paralst[0];
	    parabest=paralst[0];
	    parastd.clear();LoveRAparastd.clear();LoveAZparastd.clear();

//what's the size of love para
	    for(j=0;j<paralst[0].npara;j++){// the size of Love para is the same as V para
                parastd.push_back(0.);
                LoveRAparastd.push_back(0.);
		LoveAZparastd.push_back(LAZp0);
		paraavg.parameter[j]=paraavg.LoveRAparameter[j]=paraavg.LoveAZparameter[j][0]=paraavg.LoveAZparameter[j][1]=0.;
            }//for j

	    for(i=0;i<size;i++)
		{if(paralst[i].misfit<mismin) {idmin=i;mismin=paralst[i].misfit;}}
	    //min(2*Mis_min,Mis_min+0.5)
	    parabest=paralst[idmin];
	    printf("@@@ para_avg mismin=%g\n",mismin);
	    //if(mismin>1.5) mismin=mismin*2;
	    //else mismin=mismin+0.5;
	    mismin=mismin+0.5;
	    // mismin=1e10;
	    Ngood=0;
	    for(i=0;i<size;i++){
		if(paralst[i].misfit<mismin){
		  Ngood++;
		  idlst.push_back(i);
		  for(j=0;j<paralst[0].npara;j++){
			paraavg.parameter[j]+=paralst[i].parameter[j];
			paraavg.LoveRAparameter[j]+=paralst[i].LoveRAparameter[j];
			paraavg.LoveAZparameter[j][0]+=paralst[i].LoveAZparameter[j][0];
			paraavg.LoveAZparameter[j][1]+=paralst[i].LoveAZparameter[j][1];}
		}//if
	    }//for i
	   for(i=0;i<paralst[0].npara;i++){
		paraavg.parameter[i]/=Ngood;
		paraavg.LoveRAparameter[i]/=Ngood;
		paraavg.LoveAZparameter[i][0]/=Ngood;
		paraavg.LoveAZparameter[i][1]/=Ngood;
	   }//for i
	   //cout<<"Ngood="<<Ngood<<endl; 
	   printf("@@@ para_avg size =%d, Ngood=%d mismin=%g\n",size,Ngood,mismin);
	   //if(Ngood>size){printf("WRONG HERE! Ngood > size!!\n");exit(0);}
	   for(i=0;i<paralst[0].npara;i++){
		for(j=0;j<Ngood;j++){
		    k=idlst[j];
		    parastd[i]+=pow(paralst[k].parameter[i]-paraavg.parameter[i],2);
		    LoveRAparastd[i]+=pow(paralst[k].LoveRAparameter[i]-paraavg.LoveRAparameter[i],2);
		    LoveAZparastd[i][0]+=pow(paralst[k].LoveAZparameter[i][0]-paraavg.LoveAZparameter[i][0],2);
		    LoveAZparastd[i][1]+=pow(paralst[k].LoveAZparameter[i][1]-paraavg.LoveAZparameter[i][1],2);
		}
		parastd[i]=sqrt(parastd[i]/Ngood);
		LoveRAparastd[i]=sqrt(LoveRAparastd[i]/Ngood);
		LoveAZparastd[i][0]=sqrt(LoveAZparastd[i][0]/Ngood);
		LoveAZparastd[i][1]=sqrt(LoveAZparastd[i][1]/Ngood);
		if (parastd[i]>50){//indicating it's phi and has period problem that hasn't been taken into account
			badv.push_back(i);
		}
	    }//for i
	    //========added Nov23, 2013; check the para_std, if std too large (>50), indicating it's phi and has period problem that hasn't been taken into account
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		paraavg.parameter[ibadp]=0.;//paraavg.LoveRAparameter[ibadp]=paraavg.LoveAZparameter[ibadp][0]=paraavg.LoveAZparameter[ibadp][1]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			if(paralst[k].parameter[ibadp]>90){//the phi should be the same among parameter, LoveRApara.. LoveAZpara..
				paralst[k].parameter[ibadp]-=180;
				//paralst[k].LoveRAparameter[ibadp]-=180;
				//paralst[k].LoveAZparameter[ibadp][0]-=180;
				//paralst[k].LoveAZparameter[ibadp][1]-=180;
			}
			paraavg.parameter[ibadp]+=paralst[k].parameter[ibadp];
		}//for j
		paraavg.parameter[ibadp]/=Ngood;
		paraavg.LoveRAparameter[ibadp]=paraavg.parameter[ibadp];
		paraavg.LoveAZparameter[ibadp][0]=paraavg.parameter[ibadp];
		paraavg.LoveAZparameter[ibadp][1]=paraavg.parameter[ibadp];
	    }//for i<badv.size
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		parastd[ibadp]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			parastd[ibadp]+=pow(paralst[k].parameter[ibadp]-paraavg.parameter[ibadp],2);
		}//for j
		parastd[ibadp]=sqrt(parastd[ibadp]/Ngood);
		LoveRAparastd[ibadp]=parastd[ibadp];
		LoveAZparastd[ibadp][0]=parastd[ibadp];
		LoveAZparastd[ibadp][1]=parastd[ibadp];
	    }//for i
	    return 1;
	}//para_avg
//-----------------------------------------------------
//----------------------------------------------------- 
	int para_avg_VOnly(vector<paradef> &paralst, paradef &parabest, paradef &paraavg,vector<double> &parastd,vector<double> &LoveRAparastd, vector<vector<double> > &LoveAZparastd,vector<int> &idlst){
 	// this compute para_avg for Vpara only, ignore the Lovepara
	//get the average parameters from a list of acceptable paras
	//also, get the best parameter (min misfit)
	    int i,j,k,size,Ngood,ibadp,idmin;
	    double mismin=1e10;
	    vector<double> LAZp0(2,0.);
	    vector<int> badv;

	    initpara(paraavg);
  	    initpara(parabest);
	    size=paralst.size();
	    if(size<1) return 0;
		
	    paraavg=paralst[0];
	    parabest=paralst[0];
	    parastd.clear();LoveRAparastd.clear();LoveAZparastd.clear();

//what's the size of love para
	    for(j=0;j<paralst[0].npara;j++){// the size of Love para is the same as V para
                parastd.push_back(0.);
		paraavg.parameter[j]=0.;
            }//for j

	    for(i=0;i<size;i++)
		{if(paralst[i].misfit<mismin) {idmin=i;mismin=paralst[i].misfit;}}
	    //min(2*Mis_min,Mis_min+0.5)
	    parabest=paralst[idmin];
	    printf("@@@ para_avg mismin=%g\n",mismin);
	    //if(mismin>1.5) mismin=mismin*2;
	    //else mismin=mismin+0.5;
	    mismin=mismin+0.5;
	    // mismin=1e10;
	    Ngood=0;
	    for(i=0;i<size;i++){
		if(paralst[i].misfit<mismin){
		  Ngood++;
		  idlst.push_back(i);
		  for(j=0;j<paralst[0].npara;j++){
			paraavg.parameter[j]+=paralst[i].parameter[j];
		   }
		}//if
	    }//for i
	   for(i=0;i<paralst[0].npara;i++){
		paraavg.parameter[i]/=Ngood;
	   }//for i
	   //cout<<"Ngood="<<Ngood<<endl; 
	   printf("@@@ para_avg size =%d, Ngood=%d mismin=%g\n",size,Ngood,mismin);
	   //if(Ngood>size){printf("WRONG HERE! Ngood > size!!\n");exit(0);}
	   for(i=0;i<paralst[0].npara;i++){
		for(j=0;j<Ngood;j++){
		    k=idlst[j];
		    parastd[i]+=pow(paralst[k].parameter[i]-paraavg.parameter[i],2);
		}
		parastd[i]=sqrt(parastd[i]/Ngood);
		if (parastd[i]>50){//indicating it's phi and has period problem that hasn't been taken into account
			badv.push_back(i);
		}
	    }//for i
	    //========added Nov23, 2013; check the para_std, if std too large (>50), indicating it's phi and has period problem that hasn't been taken into account
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		paraavg.parameter[ibadp]=0.;//paraavg.LoveRAparameter[ibadp]=paraavg.LoveAZparameter[ibadp][0]=paraavg.LoveAZparameter[ibadp][1]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			if(paralst[k].parameter[ibadp]>90){//the phi should be the same among parameter, LoveRApara.. LoveAZpara..
				paralst[k].parameter[ibadp]-=180;
			}
			paraavg.parameter[ibadp]+=paralst[k].parameter[ibadp];
		}//for j
		paraavg.parameter[ibadp]/=Ngood;
	    }//for i<badv.size
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		parastd[ibadp]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			parastd[ibadp]+=pow(paralst[k].parameter[ibadp]-paraavg.parameter[ibadp],2);
		}//for j
		parastd[ibadp]=sqrt(parastd[ibadp]/Ngood);
	    }//for i
	    return 1;
	}//para_avg_VOnly
//-----------------------------------------------------
