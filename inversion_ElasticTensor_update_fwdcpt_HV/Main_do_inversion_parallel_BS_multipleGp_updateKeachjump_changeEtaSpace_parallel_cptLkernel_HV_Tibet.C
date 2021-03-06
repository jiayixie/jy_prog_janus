// pay attention to the included CALMineos, the "computeDisp4IsoVs" version is used to compute disp correlated with Isotropic Vs, not anisotropic model in which Vsv&Vsh are different
// this version, could handle multiple-peak phi case (<=2 peaks); in such case, the modavg would be multiple
// this version, use the para_avg_multiple_gp_v3.C, which enables the computation of parabest for each phi group, so there are multiple parabest output
// this version, use CALinv_isolay_rf_parallel_saveMEM_BS_updateK.C, which updated the Vkernel&Lkernel during the inversion
// this version, use CALpara_isolay_BS_newV2L_changeEtaSpace.C, which require eta<=1.1
// this Tibet version is for the Tibet's Elastic tensor inversion

#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include <Eigen/Core>
#include<string>
#include <chrono>
#include <random>

using namespace std;
using namespace Eigen;

unsigned seed=chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

#include"./string_split.C"
#include"./generate_Bs.C"
#include"./gen_random_cpp.C"
#include"./INITstructure_BS_HV.h"
//#include"CALpara_isolay_BS.C"
#include "CALpara_isolay_BS_newV2L_changeEtaSpace_HV.C"
#include"./CALgroup_smooth_BS.C"
#include"CALmodel_LVZ_ET_BS_HV_vpconstr_Tibet.C"
#include"CALforward_Mineos_readK_parallel_BS_newV2L_parallel_cptLkernel_HV.C"
#include "./ASC_rw_HV.C"
#include "./BIN_rw_Love.C"
#include "CALinv_isolay_rf_parallel_saveMEM_BS_updateK_eachjump_parallel_cptLkernel_HV_v2.C"
//#include "para_avg_multiple_gp_v4.C" 
//#include "Test_fwd_cpt.C"
#define _USE_MATH_DEFINES



int main(int argc, char *argv[])
{
int npoint,Rsurflag,Lsurflag,Rmonoc,Lmonoc,PosAnic,iitercri1,ijumpcri1,ijumpcri2;
int i,j,k,isoflag,Nprem,k1,k2;
int AziampRsurflag,AziphiRsurflag,AziampLsurflag,AziphiLsurflag;
double bestmisfit, misfit, L,misfitcri;
float depcri1,depcri2,qpcri,qscri,lon,lat;
char inponm[200],Rgpindir[100],Rphindir[100],Lphindir[100],Lgpindir[100],kernelnmR[200],kernelnmL[200],kernelnmRHV[200];
vector<string> AziampRdispnm,AziphiRdispnm,AziampLdispnm,AziphiLdispnm;
vector<string> Rdispnm,Ldispnm;
char nodeid[5],str[150],modnm[100],Lparanm[100],Rparanm[100],PREMnm[100],fparanm[100];
time_t start;
modeldef model0,model1,ttmodel,modelref,modelnew;
modeldef modeltemp,modelavg1;
paradef para0,para1,para2,pararef,paranew,tpara,paraavg1,ttpara,parabest;
FILE *inpo,*fkernel,*fmisfit;
vector<vector<double> >  PREM;
vector<vector<vector<double> > > Vkernel,Lkernel;
vector<int> Lvmono,Lvgrad,Rvmono,Rvgrad,Vposani,Viso,idlst;
vector<paradef> paralst,paralstBS,parabestlst;
vector<double> parastd,LoveRAparastd;
vector<vector<double> > LoveAZparastd;
char modnm1[200],fRdispnm[200],fLdispnm[200],modnm2[200],fAZRdispnm[200],fAZLdispnm[200];
char outinitnm[200],lay[20],dirlay[100],tmpstr[500],fbinnm1[200],fbinnm2[200];
float inpamp,inpphi;
int flagreadVkernel,flagreadLkernel,flagupdaterho;

/*if(argc!=11){
printf("Usage: xx 1]input_point_file 2]output_dir_name 3]input_vsv_dir_name 4]Rphindir 5]Rgpindir 6]Lphindir 7]Lgindir 8]fparanm(para.in file) 9]flagreadVkernel(1-Y; 0-N) 10]num_of_thread\n");
printf("1] node_id node_lon node_lat\n3]dir/initmod/vsv_node_lon_lat.mod\n");
exit(0);
}
*/

  //----------------PARAMETERS-----------------------------------------
  isoflag=0; //isoflag==1: Vsv=Vsh, isoflag==0: Vsv!=Vsh
  Rsurflag=1; //surflag==1: open phase only. surfalg ==3 open phase and group, surflag==2: open group only; surflag=4: hv only; surflag=5:p+hv; surflag=6: g+hv; surflag=7: g+p+hv
  Lsurflag=1;
  AziampRsurflag=1;
  AziphiRsurflag=1;
  AziampLsurflag=0;
  AziphiLsurflag=0;
  //after I have changed the way the misfit is computed (compute_misfitDISP), the inpamp & inpphi becomes useless
  inpamp=0.25;//useless 0.25; //the weight of the azi_aniso disp curve, amp part (0~1)
  inpphi=0.25;//useless 0.25; //the weight of the azi_aniso disp curve, ang part (0-1)
  //the weight of iso dispersion curve is 1-inpamp-inpphi  
  iitercri1=50000;//100000;//12000 (mod1, 1cstlay)
  ijumpcri1=10; //atoi(argv[10]); // set it to be the same as number_of_thread
  depcri1=20.0;
  depcri2=80.0;
  qpcri=900.;//900.;
  qscri=250.;
  Rmonoc=1;//useless
  Lmonoc=1;//useless
  PosAnic=1;
  flagupdaterho=0;
  //Rvmono.push_back(0);
  Rvmono.push_back(1);
  Rvmono.push_back(2);
  Lvmono.push_back(1);
  Lvmono.push_back(2);
  //Rvgrad.push_back(0); // require the 1st two values in that group are increasing
  Rvgrad.push_back(1);
  //Rvgrad.push_back(2);
  //Lvgrad.push_back(0);
  Lvgrad.push_back(1);
  //Vposani.push_back(0);
  Vposani.push_back(1);
  Vposani.push_back(2);
  /*
  Viso.push_back(0);
  Viso.push_back(1);
  Viso.push_back(2);
  */
  //Vposani.push_back(1);Vposani.push_back(2);
  k1=0;k2=1;
  //----------------------------------------------------------------------

  //sprintf(PREMnm,"/home/jiayi/progs/jy/Mineos/Mineos-Linux64-1_0_2/DEMO/models/prem_noocean.txt");
  sprintf(PREMnm,"/home/jixi7887/progs/jy/Mineos/Mineos-Linux64-1_0_2/DEMO/models/ak135_iso_nowater.txt");
  ///*
  sprintf(inponm,argv[1]);
  sprintf(dirlay,argv[2]);
  sprintf(modnm,argv[3]);//starting model
  sprintf(Rphindir,argv[4]);
  sprintf(Rgpindir,argv[5]);
  sprintf(Lphindir,argv[6]);
  sprintf(Lgpindir,argv[7]);
  sprintf(fparanm,argv[8]);
  flagreadLkernel=flagreadVkernel=atoi(argv[9]);
  int num_thread=atoi(argv[10]);
  //*/
  /*
  sprintf(inponm,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/point1.txt");
  sprintf(dirlay,"/lustre/janus_scratch/jixi7887/code_test/inv_v1_testHV/P12A_-115.0_39.4_inv_v1_testHV");
  sprintf(modnm,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/Data/test/P12A_iso.mod");
  sprintf(Rphindir,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/Data/test");
  sprintf(Rgpindir,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/Data/test");
  sprintf(Lphindir,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/Data/test");
  sprintf(Lgpindir,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/Data/test");
  sprintf(fparanm,"/projects/jixi7887/work/code_test/test_HVratio_Mineos/para_v1.txt");
  flagreadVkernel=1;//atoi(argv[1]);
  int num_thread=atoi(argv[1]);
  */

  readPREM(PREMnm,PREM,Nprem);

  sprintf(tmpstr,"if [ ! -d %s ]; then mkdir %s; fi",dirlay,dirlay);
  system(tmpstr);
  sprintf(tmpstr,"if [ ! -d %s/initmod ]; then mkdir %s/initmod; fi",dirlay,dirlay);
  system(tmpstr);
sprintf(tmpstr,"if [ ! -d %s/binmod ]; then mkdir %s/binmod; fi",dirlay,dirlay);
  system(tmpstr);

  //---------------------------------------------------------
  if((inpo=fopen(inponm,"r"))==NULL){
	printf("Cannot open points file %s!\n",inponm);exit(0);
  }
  npoint=0;

  while(1){
    if(fscanf(inpo,"%s %f %f",&nodeid[0],&lon,&lat)==EOF)
	break;
    npoint++;
    printf("Begin to work on point %d: id=%s lon=%f lat=%f\n",npoint,nodeid,lon,lat);

    start=time(0);
    bestmisfit=1e10;

    Rdispnm.clear();
    sprintf(str,"%s/dispISO_%.1f_%.1f.txt",Rphindir,lon,lat);
    Rdispnm.push_back(str);
    if(Rsurflag==5)
    {sprintf(str,"%s/HV.Ray_%.1f_%.1f.txt",Rphindir,lon,lat);
    Rdispnm.push_back(str);}

    Ldispnm.clear();
    sprintf(str,"%s/dispISO_%.1f_%.1f.txt",Lphindir,lon,lat);
    Ldispnm.push_back(str);

    AziampRdispnm.clear();
    sprintf(str,"%s/dispAMP_%.1f_%.1f.txt",Rphindir,lon,lat);
    AziampRdispnm.push_back(str);

    AziphiRdispnm.clear();
    sprintf(str,"%s/dispPHI_%.1f_%.1f.txt",Rphindir,lon,lat);
    //sprintf(str,"%s/aziphi_restore_unc/aziphi_%.1f_%.1f.txt",Rphindir,lon,lat);//#############HEY TEMPERARY########## TEST
    AziphiRdispnm.push_back(str);

    AziampLdispnm.clear();
    sprintf(str,"%s/dispAMP_%.1f_%.1f.txt",Lphindir,lon,lat);
    AziampLdispnm.push_back(str);

    AziphiLdispnm.clear();
    sprintf(str,"%s/dispPHI_%.1f_%.1f.txt",Lphindir,lon,lat);
    AziphiLdispnm.push_back(str);



    //-----the final output model, that will serve as input model for the next step 
    sprintf(outinitnm,"%s/initmod/ani_%s_%.1f_%.1f.mod",dirlay,nodeid,lon,lat);
    //-----the outpur binary files, all the accepted models during inversion
    sprintf(fbinnm1,"%s/binmod/ani_%s_%.1f_%.1f.bin",dirlay,nodeid,lon,lat);
    sprintf(fbinnm2,fbinnm1);//overwrite the first binary output


    //---------------------------------------------------------
    initmodel(model0);
    initmodel(modelref);

    initpara(para0);
    initpara(para1);
    initpara(pararef);
  
    readdisp(model0,Rdispnm,Ldispnm,AziampRdispnm,AziphiRdispnm,AziampLdispnm,AziphiLdispnm,Rsurflag,Lsurflag,AziampRsurflag,AziphiRsurflag,AziampLsurflag,AziphiLsurflag); 
    /*==check===
    //printf("test-- the number of AZ data: AZRamp.npper=%d AZRamp.pvel.size()=%d AZRamp.pvelo.size=%d\n",model0.data.AziampRdisp.npper,model0.data.AziampRdisp.pvel.size(),model0.data.AziampRdisp.pvelo.size());
    printf("the Rayleigh wave data:\n");
    for(i=0;i<model0.data.Rdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.Rdisp.pper[i],model0.data.Rdisp.pvelo[i],model0.data.Rdisp.unpvelo[i]);
    }
    printf("the Rayleigh Azi wave data:\n");
    for(i=0;i<model0.data.AziampRdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.AziampRdisp.pper[i],model0.data.AziampRdisp.pvelo[i],model0.data.AziampRdisp.unpvelo[i]);
    }
    for(i=0;i<model0.data.AziphiRdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.AziphiRdisp.pper[i],model0.data.AziphiRdisp.pvelo[i],model0.data.AziphiRdisp.unpvelo[i]);
    }
    */
    
    readmodAniso(model0,modnm);// both m.g.LV/Rv are filled regardless of flags. (readin iso model)
    printf("finish read mod\n");
    printf("#########group[0].vsvvalue[0]=%g vshvalue[0]=%g\n",model0.groups[0].vsvvalue[0],model0.groups[0].vshvalue[0]);
    /*===check===
    for(i=0;i<model0.ngroup;i++){
    	printf("group %d thick=%f\n",i,model0.groups[i].thick);
	for(j=0;j<model0.groups[i].np;j++){
		printf("  np %d\n vsv=%g vsh=%g vpv=%g vph=%g eta=%g theta=%g phi=%g vpvs=%g\n",j,model0.groups[i].vsvvalue[j],model0.groups[i].vshvalue[j],model0.groups[i].vpvvalue[j],model0.groups[i].vphvalue[j],model0.groups[i].etavalue[j],model0.groups[i].thetavalue[j],model0.groups[i].phivalue[j],model0.groups[i].vpvs);
	}
    
    }
    //exit(0);
    */
    readpara(para0,fparanm);


    mod2para(model0,para0,para1);////fill both para.R/Lpara0 (they could be inequal if Rf*Lf>0, they are equal if Rf*Lf=0)

    checkParaModel(para1,model0,Viso);

    /*===check===
    printf("check readpara & mod2para----\n");
    printf("inpara.flag=%d\n",para0.flag);
    for(i=0;i<para1.npara;i++){
	printf("parameter %d\n",i);
	printf(" value=%g  dv=%g ng=%g nv=%g pflag=%g LVflag=%g RWflag=%g LWflag=%g AZflag=%g\n",para1.parameter[i],para1.para0[i][2],para1.para0[i][4],para1.para0[i][5],para1.para0[i][6],para1.para0[i][7],para1.para0[i][8],para1.para0[i][9],para1.para0[i][10]);
    }*/
    // ; BS
    //---check---
    ttmodel=model0;
    printf("Before Bsp2P\n");
    for(i=0;i<ttmodel.ngroup;i++){
        printf("BSmodel group%d\n",i);
        for(j=0;j<ttmodel.groups[i].np;j++){
                printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
        }
    }
    //
    
    modeldef modelP,modelBS;
    paradef paraP,paraBS;
    Bsp2Point(model0,para1,modelP,paraP,flagupdaterho);
    modelBS=model0;paraBS=para1;
    model0=modelP;para1=paraP;

    //---check-- 
    printf("\n\nafter Bsp2P\n");
    ttmodel=modelBS;
    for(i=0;i<ttmodel.ngroup;i++){
	printf("BSmodel group%d\n",i);
	for(j=0;j<ttmodel.groups[i].np;j++){
                printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
	}
    }
    ttmodel=modelP;
    for(i=0;i<ttmodel.ngroup;i++){
        printf("Pmodel group%d\n",i);
        for(j=0;j<ttmodel.groups[i].np;j++){
                printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
        }
    }
    //

    /*printf("\nafter Bsp2Point inpara.flag=%d\n",para0.flag);
    for(i=0;i<para1.npara;i++){
	printf("parameter %d\n",i);
	printf(" value=%g  dv=%g ng=%g nv=%g pflag=%g LVflag=%g RWflag=%g LWflag=%g AZflag=%g\n",para1.parameter[i],para1.para0[i][2],para1.para0[i][4],para1.para0[i][5],para1.para0[i][6],para1.para0[i][7],para1.para0[i][8],para1.para0[i][9],para1.para0[i][10]);
    }*/
    // ; BS
    Vpara2Lovepara(para1,model0,flagupdaterho);//there is para2mod inside the code



    if(model0.flag==0){updatemodel(model0,flagupdaterho);}
    compute_dispMineos(model0,PREM,Nprem,Rsurflag,Lsurflag,0);

    printf("@@@ check, main, from initial model\n");
    for(j=0;j<model0.data.Rdisp.npper;j++)printf("  @@@ check, T=%g,vin=%g vold=%g\n",model0.data.Rdisp.pper[j],model0.data.Rdisp.pvelo[j],model0.data.Rdisp.pvel[j]);

    printf("@@@ check, do V kernel ====\n");
    //---obtain V kernel ---
    sprintf(kernelnmR,"%s/VkernelRp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/VkernelLp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmRHV,"%s/VkernelRHV1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);

    if(flagreadVkernel==1){
    	if((read_kernel(para1,model0,Vkernel,kernelnmR,kernelnmL,kernelnmRHV,Rsurflag,Lsurflag,PREM,Nprem))==0){
		 printf ("#####!! read_kernel failed\n");
    	 	 sprintf(str,"echo point %d: id=%s lon=%f lat=%f >> point_rdKernel_failed_Ani.txt",npoint,nodeid,lon,lat);
         	 system(str);
		 continue;
     	} // if readkernel()==0    
    }//if flagreadkernel==1
    else {
    	compute_Vkernel(para1,model0,Vkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho,0);
	cout<<"check finish compute_Vkernel\n";
	write_kernel(Vkernel,model0,para1,kernelnmR,kernelnmL,kernelnmRHV,Rsurflag,Lsurflag);
	cout<<"check finish write_Vkernel\n";
    }//else flagreadkernel==1
    //====


    printf("@@@ check, do L kernel ====\n");
    //---obtain Love kernel ---
    sprintf(kernelnmR,"%s/LkernelRp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/LkernelLp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmRHV,"%s/LkernelRHV1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    if(flagreadLkernel==1){
    	if((read_kernel(para1,model0,Lkernel,kernelnmR,kernelnmL,kernelnmRHV,Rsurflag,Lsurflag,PREM,Nprem))==0){
		 printf ("#####!! read_kernel failed\n");
    	 	 sprintf(str,"echo point %d: id=%s lon=%f lat=%f >> point_rdKernel_failed_Ani.txt",npoint,nodeid,lon,lat);
         	 system(str);
		 continue;
     	} // if readkernel()==0    
    	
    }
    else{
	compute_Lkernel(para1,model0,Lkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho,0);
	cout<<"check finish compute Lkernel\n";
    	//Vkernel2Lkernel(para1,model0,Vkernel,Lkernel,flagupdaterho);
    	write_kernel(Lkernel,model0,para1,kernelnmR,kernelnmL,kernelnmRHV,Rsurflag,Lsurflag);
	cout<<"check finish write_Lkernel\n";
    }
/*  //----------------------------
  //test the accuracy of the forward computation
  paradef RApara;
  modeldef RAmodel;

  pararef=para1;
  modelref = model0;
  //---perturb the parameters
  for (int kk=4;kk<=18;kk=kk+7)
  	paraBS.parameter[kk]=paraBS.parameter[kk]*1.04;
  Bsp2Point(modelBS,paraBS,model0,para1,flagupdaterho);
  ttmodel=model0;
  para2mod(para1,ttmodel,model0);
  updatemodel(model0,flagupdaterho);
  Vpara2Lovepara(para1,model0,flagupdaterho);//--before this, only the para.parameters are updated, not the para.LoveRAparameter
  //---get the RAmodel, RApara from the Point para and model (para1 and model0)
  get_RAmodpara(para1,model0,RApara,RAmodel,flagupdaterho);
  
  //
  //--do compute RAdisp with
  //Mineos
  printf("compute RAdisp M\n");
  computeRAdisp_Mineos(RAmodel,PREM,Nprem);
  //Lkernel
  printf("compute RAdisp L\n");
  computeRAdisp_Lkernel(RApara,RAmodel,pararef,modelref,Vkernel,Lkernel);
  //Vkernel
  printf("compute RAdisp V\n");
  computeRAdisp_Vkernel(RApara,RAmodel,pararef,modelref,Vkernel,Lkernel);  
  //-----------------------------
printf("\n\nexit!!\n\n");
exit(0);
*/
/*    
 int d = Math.abs(a - b) % 360;
     int r = d > 180 ? 360 - d : d;
*/
    printf("test-- cpt misfit\n");
    //somethihng to do, 1st, in generating the ang disp curve, make it smooth
    //something to do, in computing the misfit for the angle phi, need to take the period of the angle into account. modify the compute_misfitDISP in CALmodel_LVZ_ET.C
    compute_misfitDISP(model0,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,inpamp,inpphi);
    //is the AZ disp filled? so are they taken into account in the misfit cpt? A: yes, both amp and angle are 0

    modelref=model0;
    pararef=para1;
   
    ttmodel=modelBS;
    para2mod(paraBS,ttmodel,modelBS);
    updatemodel(modelBS,flagupdaterho);

    float theta;
    printf("test-- do inv\n");
    if((do_inv_BS(num_thread,2,-2,paralst,paraBS,modelBS,pararef,modelref,Rvmono,Lvmono,Rvgrad,Lvgrad,PREM,Vkernel,Lkernel,k1,k2,start,isoflag,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,Nprem,Rmonoc,Lmonoc,PosAnic,Vposani,iitercri1,ijumpcri1,fbinnm1,inpamp,inpphi,flagupdaterho))==0){
    //if((do_inv(2,-2,paralst,pararef,modelref,Rvmono,Lvmono,Rvgrad,Lvgrad,PREM,Vkernel,Lkernel,k1,k2,start,isoflag,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,Nprem,Rmonoc,Lmonoc,PosAnic,Vposani,iitercri1,ijumpcri1,fbinnm1,inpamp,inpphi,flagupdaterho))==0){
    	sprintf(str,"echo %s %.1f %.1f >> point_do_inv_failed.txt",nodeid,lon,lat);
	system(str);
    	continue;
    }
   printf("test-- END of inversion\n");

     
  }//while1
  return 1;
}
//main
