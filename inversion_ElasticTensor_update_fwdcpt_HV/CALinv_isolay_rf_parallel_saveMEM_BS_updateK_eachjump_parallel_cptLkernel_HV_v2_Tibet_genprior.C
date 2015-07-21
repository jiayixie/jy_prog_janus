// this rf version, the get_misfitKernel is slightly different
// right now, the isoflag is useless, since it's not used in the gen_newpara function
// the cpu time counted is strange, not sure what the problem is, but the wall time is working properly
// this version, will update the kernel every N steps.
// this version, update K inside each jump; i.e., inside each jump, if iaccp>X then update the kernel. modified Mar 12, 2014
//  this version, parallel the inversin process (previous only paralell the kernel_computation process)
//  this version (HV), incooperate the HV version subroutines; use the new para2mod+gen_newpara and slightly modified Bsp2Point(rm the para2mod inside it) functions.
// this version(genprior) is used to generate the prior distributions. will remove the effect of likelihiood/misfit
//
//-------
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
//-------
int temp_writepara(FILE *f, paradef para,modeldef model, int iaccp, int ithread, int flag)
{
//
// (1)ithread (2)flag (3)iaccp  (4)'L' (5)m.d.L  (6)'RL' (7)m.d.R.L  (8)m.d.AZampR.L (9)m.d.AZphiR.L (10)'LL' (11)m.d.L.L  (12)m.d.AZampL.L  (13)m.d.AZphiL.L  (14)'misfit' (15)m.d.misfit  (16)'Rm' (17)m.d.R.misfit  (18)m.d.AZampR.misfit  (19)m.d.AZphiR.misfit  (20)'Lm' (21)m.d.L.m   (22)m.d.AZampL.m  (23)m.d.AZphiL.m	
	int i;
	fprintf(f,"%d %d %d L %8g  RL %8g %8g %8g %8g LL %8g %8g %8g  misfit %8g Rm %8g %8g %8g %8g Lm %8g %8g %8g\t\t",ithread,flag,iaccp,model.data.L,model.data.Rdisp.L,model.data.AziampRdisp.L,model.data.AziphiRdisp.L, model.data.Rdisp.hvL,model.data.Ldisp.L, model.data.AziampLdisp.L,model.data.AziphiLdisp.L,  model.data.misfit, model.data.Rdisp.misfit, model.data.AziampRdisp.misfit, model.data.AziphiRdisp.misfit,model.data.Rdisp.hvmisfit,model.data.Ldisp.misfit, model.data.AziampLdisp.misfit,model.data.AziphiLdisp.misfit);
	
	for(i=0;i<para.npara;i++)
	{ fprintf(f,"%8g ",para.parameter[i]);
	}
	fprintf(f," RApara: ");
	for(i=0;i<para.npara;i++)
	{ fprintf(f,"%8g ",para.LoveRAparameter[i]);
	}
	/*
	//---check---write out the Love wave 4psi signal amp, write peak and avg of the 1st 4 values
	double max=-1.;
	double sum=0.;
	int n=4;
	for(i=0;i<model.data.AziampLdisp.npper;i++){
		if(model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200>max)max=model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200;
		if(i<n)sum+=model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200;
	}
	fprintf(f,"\t\t %g %g",max,sum/n);
	*/
	//---check---write the c in the crust layers-------
	//fprintf(f,"c: ");
	//for(int m=10;m<15;m++){
	//float c=model.laym0.rho[m]*(0.125*pow(model.laym0.vph[m],2)+0.125*pow(model.laym0.vpv[m],2)-0.25*model.laym0.eta[m]*(pow(model.laym0.vph[m],2)-2*pow(model.laym0.vsv[m],2))-0.5*pow(model.laym0.vsv[m],2));
	//fprintf(f,"\t%5f %5f %5f %5f %5f %5f   %g",model.laym0.vsv[m],model.laym0.vsh[m],model.laym0.vpv[m],model.laym0.vph[m],model.laym0.eta[m],model.laym0.rho[m],c);}
	//-----
	fprintf(f,"\n");
	return 1;
}//temp_writepara


//int do_inv_BS(int id,double misfitcri,vector<paradef> &aralst,paradef refparaBS,modeldef refmodelBS,paradef refpara,modeldef refmodel,vector<int> Rvmono,vector<int> Lvmono,vector<int> Rvgrad, vector<int> Lvgrad,vector<vector<double> > PREM, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel,int k1,int k2,time_t start,int isoflag, int Rsurflag, int Lsurflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag,int Nprem,int Rmonoc,int Lmonoc, int PosAni, vector<int> Vposani, int iitercri,int ijumpcri,char *fbinnm,float inpamp, float inpphi, int flagupdaterho){
//int do_inv_BS(const int num_thread,const int id,const double misfitcri, vector<paradef> &paralst,const paradef refparaBS,const modeldef refmodelBS,const paradef refpara,const modeldef refmodel,const vector<int> Rvmono,const vector<int> Lvmono,const vector<int> Rvgrad,const vector<int> Lvgrad,const vector<vector<double> > PREM,const  vector<vector<vector<double> > > Vkernel, const vector<vector<vector<double> > > Lkernel,const int k1,const int k2,const time_t start,const int isoflag,const int Rsurflag, const int Lsurflag,const int AziampRflag,const int AziampLflag,const int AziphiRflag,const int AziphiLflag,const int Nprem,const int Rmonoc,const int Lmonoc,const int PosAni, const vector<int> Vposani, const int iitercri,const int ijumpcri,const char *fbinnm,const float inpamp, const float inpphi, const int flagupdaterho){
int do_inv_BS(const int num_thread,const int id,const double misfitcri, vector<paradef> &paralst,const paradef refparaBS,const modeldef refmodelBS, const paradef refpara0, const modeldef refmodel0,const vector<int> Rvmono,const vector<int> Lvmono,const vector<int> Rvgrad,const vector<int> Lvgrad,const vector<vector<double> > PREM, const vector<vector<vector<double> > > Vkernel0, const vector<vector<vector<double> > > Lkernel0,const int k1,const int k2,const time_t start,const int isoflag,const int Rsurflag, const int Lsurflag,const int AziampRflag,const int AziampLflag,const int AziphiRflag,const int AziphiLflag,const int Nprem,const int Rmonoc,const int Lmonoc,const int PosAni, const vector<int> Vposani, const int iitercri,const int ijumpcri,const char *fbinnm,const float inpamp, const float inpphi, const int flagupdaterho){
  int i,p6;
  int tflag,iiter,iaccp,ibad,lastiaccp,idphiC,flagidphiC;
  double oldL,oldRL,oldLL,newL,newRL,newLL,oldmisfit,newmisfit;
  double prob,prandom;
  char str[500];
  time_t now,dtime;
  const time_t start2=time(0);
  const double start_walltime0=time(NULL);
  const double start_cputime0=clock();
  const double start_walltime=omp_get_wtime();//get_wall_time();
  const double start_cputime=get_cpu_time();
  vector<vector<vector<double> > > Vkernel,Lkernel;
  modeldef refmodel;
  paradef refpara;  

  modeldef tmodel,tempmod;
  paradef para1,para2,temppara;
  ofstream outbin,outbinRA;
  if ( id >1 ){
  outbin.open(fbinnm,ios_base::out|ios_base::binary);
  sprintf(str,"%s_effTI",fbinnm);
  outbinRA.open(str,ios_base::out|ios_base::binary);
  }
  //---get the phi id for crust, the group1's phi_para,
  //---Here, I assume that group[1] is crust, should be adjusted for different situations----
 
  for( i=0;i<refpara0.npara;i++){
    if((int)refpara0.para0[i][6]==7 and (int)refpara0.para0[i][4]==1){idphiC=i;break;}
  }//for
  if(i==refpara0.npara){printf("This model has no crustal phi value\n");flagidphiC=0;}
  else
  {printf("inv, crust idphi=%d i=%d\n",idphiC,i);flagidphiC=1;}

  const int Nacc_updateK_init = 505;//505;//505;//500; //500;
  //const int accpcri=5000; //5000;//3000;//5000 ;//3000;
  //const int Nacc_updateK_init = 700;//505;//505;//500; //500;
  const int accpcri=10000; //5000;//3000;//5000 ;//3000;
  //const int accpcri=5000; //5000;//3000;//5000 ;//3000;
  const int Njump=ijumpcri; //set this to the times of the max_thread_num to max the efficiency of the code; Also, need to pay attetion if Njump is too small to be enough for the Monte-Carlo search. (usually, i require Njump>=5)
  int iloop=0;

  printf("num of jump=%d\n",Njump);
  iiter=iaccp=ibad=0;
  lastiaccp=0;

  int countitt=0,countacc=0,NKdid=0,Tacc;
  //-----
  if(num_thread>10){
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(int(num_thread/10)); // Use 10 threads for all consecutive parallel regions
  }
  //omp_set_num_threads(1);
  int flagbreak=0;
  while(iloop<600 and iaccp<accpcri){
  //maybe, just use Njump jumps all the time. eacho jump itterate iitercri times before it terminates. No break out of the code is necessary then.
  //num_threads(Njump)
  #pragma omp parallel default(none) shared(flagidphiC,flagbreak,idphiC,outbin,outbinRA,paralst,iloop,iiter,iaccp) private(para2,para1,tmodel,oldLL,oldRL,oldL,countitt,countacc,ibad,tflag,newL,newRL,newLL,newmisfit,now,prob,oldmisfit,prandom,refmodel,Vkernel,Lkernel,lastiaccp,refpara,Tacc,NKdid,p6) 
  {
  printf("jump threads=%d\n",omp_get_num_threads());//--check---
  #pragma omp for schedule(dynamic,1)
  for(int i=0;i<Njump;i++){
 	if(flagbreak==1)continue;
	refmodel=refmodel0;
	refpara=refpara0;
	Vkernel=Vkernel0;
	Lkernel=Lkernel0;
  	//--temp---
  	printf("do the jump %d/%d\n",i,Njump);
  	FILE *ftemp;
	char str[100];

	sprintf(str,"temp_LMisfit_id%d_loop%d_thread%d.txt",id,iloop,i);
	if((ftemp=fopen(str,"w"))==NULL){
		printf("Cannot open temp_LMsifit file to write! %d\n",i);
		exit(0);
	}
	paradef RApara,para1BS,para2BS;
	modeldef RAmodel,model1BS,model2BS;
   	modeldef model1,model2;
  	//-----

	int ijump=i;  

	para1=refpara;
	para2=refpara;
	model1=refmodel;
	model2=refmodel;

	para1BS=refparaBS;
	para2BS=refparaBS;
	model1BS=refmodelBS;
	model2BS=refmodelBS;

	oldLL=oldRL=oldL=0.;
	countacc=0;
	Tacc=Nacc_updateK_init;//200;
	NKdid=0;
	lastiaccp=0;

	//--for each jump, first get the starting para
	gen_newpara(refparaBS,refmodelBS, para1BS,k1);
	para2mod(para1BS,refmodelBS,model1BS);
	Bsp2Point(model1BS,para1BS,model1,para1,flagupdaterho); //para1BS->para1; there is mod2para(Pmodel,Ppara) inside the Bsp2Point function
	/*for(int n=0;n<para1BS.npara;n++){
		printf("para%d=%g->%g\n",n,refparaBS.parameter[n],para1BS.parameter[n]);
	}
	exit(0);
	*/
	
	if(PosAni>0){//require the L always faster than R
	    //cout<<"check posAni\n";//test---
            updatemodel(model1,flagupdaterho);
	    //the updatemodel, the values in the layer is the value when the ET is flat, i.e., the effect of theta is not taken into account, it's not the vel of the effective TI medium, but the vel of the un-rotated medium
            ibad=0;
	    while(positiveAni(model1,Vposani)==0){
		ibad++;
		gen_newpara(refparaBS,refmodelBS,para1BS,k1);
		para2mod(para1BS,refmodelBS,model1BS);//???
		Bsp2Point(model1BS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
		updatemodel(model1,flagupdaterho);
		if(ibad%10000==0){printf("######POSITIVE ANISO MODEL CANNOT BE SATISFIED !!!! ijump=%d\n",i);
		if(ibad>100000){break;}
		//break;//---test---
		}
	    }	
	}//if PosAni>0
	/*--check--
	*/
	
	//if(Rmonoc+Lmonoc>0){
	if(Rsurflag+Lsurflag>0){
	    //cout<<"check monoc\n";//test---
	    if(PosAni<=0){updatemodel(model1,flagupdaterho);}
	    ibad=0;
	    while(goodmodel(model1,Rvmono,Rvgrad,Lvmono,Lvgrad,Rsurflag,Lsurflag)==0)
		{ibad++;
		 gen_newpara(refparaBS,refmodelBS,para1BS,k1);
		 para2mod(para1BS,refmodelBS,model1BS);//???
		 Bsp2Point(model1BS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
		 updatemodel(model1,flagupdaterho);
		 if(ibad%10000==0){printf("##### ibad=%d GOOD MODEL CANNOT BE SATISFIED UNDER MONOC==1!! ijump=%d\n",ibad,ijump);
		 if(ibad>100000){break;}
		  //break; //---test---
		   }
		}//while	  
	}//if monoc  
	//printf("\n\n\ntest--begin get misfit Kernel\n");
	get_misfitKernel(model1,para1,refmodel,refpara,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi,flagupdaterho);
	para2=para1;
	para2BS=para1BS;
	model2=model1;
	model2BS=model1BS;

	//----then begin to sample the model space----------------INVERSION--------------------------------
	for(countitt=0;countitt<iitercri;countitt++){
		if(iiter%10000==0){
			printf("iiter=%d iaccp=%d dtime=%.1f ijump=%d iloop=%d;  ",iiter,iaccp,(float)(time(0)-start),i,iloop);
			//if((float)(time(0)-start)>3600.*2)
			if((float)(time(0)-start)>3600.*2)
			{	printf("exceed 2hr limit! break\n");
				flagbreak=1;
				break;} //break the search after 50min. modified Mar 13, 2014
			printf("inv time used wtime=%.2fs cputime=%.2fs\n",omp_get_wtime()-start_walltime,get_cpu_time()-start_cputime);
		}
	  	iiter++;
		if(countacc>5000)break; //added Jul 21, 2015
		//if (iaccp%500==0 and iaccp>0 and iaccp!=lastiaccp){// update the kernel after X step walk in the model space, what if this model is a strange model??
		/*
		 //----test prior------
		if(countacc>1500)break;//2500//modified on Mar 13, 2014
		if (countacc%Tacc==0 and countacc>0 and countacc!=lastiaccp){// update the kernel after X step walk in the model space, what if this model is a strange model??
			printf("re-compute Vkernel & Lkernel!! iloop=%d ijump=%d countacc=%d\n",iloop,ijump,countacc);
			lastiaccp=countacc;
			NKdid++;
			if(NKdid>=2)Tacc=Nacc_updateK_init*2;//500;//modified on Mar 13,2014. change the frequency of updating kernel
			if(NKdid>=4)Tacc=Nacc_updateK_init*4;//1000;
			if(NKdid>=10)Tacc=Nacc_updateK_init*10;//10000;
			if(NKdid>=20)Tacc=Nacc_updateK_init*20;//10000;
			
			if(compute_dispMineos(RAmodel,PREM,Nprem,Rsurflag,Lsurflag,ijump)!=0){// only if this model does give reasonable Mineos output, will i use it as new reference model
				
			compute_Vkernel(RApara,RAmodel,Vkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho,ijump);
			//Vkernel2Lkernel(RApara,RAmodel,Vkernel,Lkernel,flagupdaterho);
			compute_Lkernel(RApara,RAmodel,Lkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho,ijump);
			RApara.LoveAZparameter=refpara.LoveAZparameter; //set LoveAZparameter_ref to 0!
			for(int ip=0;ip<RApara.npara;ip++){// set AZcos AZsin parameters to 0! clear the AZcos AZsin parameters, modified Apr 17, 2015
				p6=(int)RApara.para0[ip][6];
				if((p6-10)*(p6-11)==0){RApara.parameter[ip]=0.;}
			}
			refmodel=RAmodel;
			refpara=RApara;
			}// if compute_dispMineos
			else{
				printf("update kernel skipped!\n");
			}
			//RAmodel is the RA part of model model2BS; RApara is the RApart & point version of para2BS
		}//update kernel
		*/
  		//---------
  		//tmodel=model2BS;
		//para2mod(para2BS,tmodel,model2BS);
		gen_newpara(para2BS,model2BS,para1BS,k2);
		para2mod(para1BS,refmodelBS,model1BS);//???
		Bsp2Point(model1BS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
  		//-----------------Criteria 1, goodmodel------------------------
         	 if(PosAni>0){//require the L always faster than R
	   		 //cout<<"check 2posani | ";//test---
            		updatemodel(model1,flagupdaterho);
            		if(positiveAni(model1,Vposani)==0){
		    		tflag=1;
				//temp_writepara(ftemp,para1,model1,countacc,i,-1);//---test
		    		continue;
            		}

          	}//if PosAni>0
	
		//printf("pass posAni ijump=%d\n",ijump);//---test---
	 	//if(Rmonoc+Lmonoc>0){
	 	if(Rsurflag+Lsurflag>0){
            		if(PosAni<=0){updatemodel(model1,flagupdaterho);}
			if(goodmodel(model1,Rvmono,Rvgrad,Lvmono,Lvgrad,Rsurflag,Lsurflag)==0){
	    		tflag=2;	
			//temp_writepara(ftemp,para1,model1,countacc,i,-2);//---test
    	    		continue;}
        	}//if monoc 
		//printf("pass goodmodel ijump=%d\n",ijump);//---test---

		get_misfitKernel(model1,para1,refmodel,refpara,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi,flagupdaterho);
		if (std::isnan(model1.data.L) or std::isnan(model1.data.Rdisp.L) or std::isnan(model1.data.Ldisp.L)){
		printf("get nan likelyhood!\n");
		continue;}//modified Mar 13, 2014
		newL=model1.data.L;//para1.L;
		newRL=model1.data.Rdisp.L;
		newLL=model1.data.Ldisp.L;
		newmisfit=para1.misfit; 

		
		//----------Criteria 2, acceptable L-----------------------------------
		/* ---test prior ---
  		if(misfitcri<-1){
	  	    if(newL<oldL){
	    		prob=(oldL-newL)/oldL;
	   		prandom=gen_random_unif01();
	    		if(prandom<prob) { 
		    		tflag=3;
				//temp_writepara(ftemp,para1,model1,countacc,i,-3);//---test
		    		continue; }
	  	    }	

		}
		else {
			if(newmisfit>misfitcri)continue;
		}
		*/
		//----------Pass criterio, accept para----------------------------------
		// para1 is the point model version of para1BS
		para2BS=para1BS;
		model2BS=model1BS;
		para2=para1;
		model2=model1;
		oldL=newL;
		oldRL=newRL;
		oldLL=newLL;
		oldmisfit=newmisfit;
		//SAVEMEM tparalst.push_back(para1);
		//printf("test-- accp Rnpara=%d,%d Lnpara=%d,%d\n",para1.Rnpara,para1.Rparameter.size(),para1.Lnpara,para1.Lparameter.size());
		//if(PosAni+Lmonoc+Rmonoc<1){// no PosAni or monoc requirement, so model hasn't been updated
		if(PosAni+Lsurflag+Rsurflag<1){// no PosAni or goodmodel requirement, so model hasn't been updated
			updatemodel(model1,flagupdaterho);
			
		}
		countacc++;		
		temp_writepara(ftemp,para1,model1,countacc,i,1);

		#pragma omp critical (updateflag_write_model)
		{
		  //printf("Hey, write model iaccp %d\n",iaccp);
		  //paralstBS.push_back(para1BS);
		  iaccp++;
		  if(id>1){	 
		  paralst.push_back(para1);
		  write_bin(model1,outbin,para1,1,i,iaccp);
		  RApara=para1;
		  tmodel=model1;
		  Lovepara2Vpara(RApara,tmodel);
		  para2mod(RApara,tmodel,RAmodel);
		  updatemodel(RAmodel,flagupdaterho);
		  if(flagidphiC>0)
		    RApara.parameter[idphiC]=para1.parameter[idphiC];//################# this is a test ###### this is added for the convience of later average process (this is the only way to seperate different phi groups)
		  write_bin(RAmodel,outbinRA,RApara,1,i,iaccp); 
		  if(flagidphiC>0)
		    RApara.parameter[idphiC]=0.; //################# this is a test ######
		  }
		  else{//--this part is just a test,---check---
			paralst.push_back(para1BS);
		  }
			//#pragma omp flush(key,iaccp,dtime)
		  //}//if
		}//critical (updateflag)
		
	}//for count

  //----temp---
  fclose(ftemp);
  }//for i	
  }//pragma
  iloop++;
  }//while
  printf("finish do_inv iiter=%d, iaccp=%d iloop=%d\n",iiter,iaccp,iloop);
  printf("TOTAL inv time used dtime=%.2fs wtime=%.2fs cputime=%.2fs\n",(float)(time(0)-start),omp_get_wtime()-start_walltime,get_cpu_time()-start_cputime);

  //if(iaccp<accpcri)return 0;
  if (id>1){
  outbin.close();
  outbinRA.close();
  }
  return 1;
}
