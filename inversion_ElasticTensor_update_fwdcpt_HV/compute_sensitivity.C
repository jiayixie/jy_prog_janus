// pay attention to the included CALMineos, the "computeDisp4IsoVs" version is used to compute disp correlated with Isotropic Vs, not anisotropic model in which Vsv&Vsh are different
// this version, could handle multiple-peak phi case (<=2 peaks); in such case, the modavg would be multiple
// this version, use the para_avg_multiple_gp_v3.C, which enables the computation of parabest for each phi group, so there are multiple parabest output
// this version, use CALinv_isolay_rf_parallel_saveMEM_BS_updateK.C, which updated the Vkernel&Lkernel during the inversion
// this version, use CALpara_isolay_BS_newV2L_changeEtaSpace.C, which require eta<=1.1

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
#include"CALmodel_LVZ_ET_BS_HV.C"
#include"CALforward_Mineos_readK_parallel_BS_newV2L_parallel_cptLkernel_HV.C"
#include "./ASC_rw_HV.C"
#include "./BIN_rw_Love.C"
#include "CALinv_isolay_rf_parallel_saveMEM_BS_updateK_eachjump_parallel_cptLkernel_HV_v2.C"
//#include "para_avg_multiple_gp_v4.C" 
//#include "Test_fwd_cpt.C"
#define _USE_MATH_DEFINES

int write_sensitivity_one_dep(FILE *f,double depth, vector<vector<double> >  Ddisp, double ddp,int flag){
  int i;
  
  if(flag==0){
  fprintf(f,"%7.3f TRHV ",(float)depth);}

  for(i=0;i<Ddisp.size();i++){
	fprintf(f,"%8.4f ",(float)(Ddisp[i][3]/ddp)); }

  return 1;
}//int write_sensitivity_one_dep

int compute_write_sensitivity(const char* fname, modeldef modelref, float c, int flag,vector<vector<double> > PREM,int Nprem ) 
{
  modeldef modelnew;
  FILE *fout;
  vector<vector<double> > DRpvel,DRgvel,DLpvel,DLgvel,DRhvratio,DLdump;
  float ddp;
  double depth;
  int i;
  double h1,h2,s;

  if((fout=fopen(fname,"w"))==NULL){printf("Cannot open sensitivity file %s to write\n",fname);exit(0);}
  //--write the titile period information
  fprintf(fout,"depth   TRHV ");
  for(i=0;i<modelref.data.Rdisp.nhvper;i++){
 	fprintf(fout,"%8.1f ",(float)modelref.data.Rdisp.hvper[i]);}
  fprintf(fout,"TRph ");
  for(i=0;i<modelref.data.Rdisp.npper;i++){
 	fprintf(fout,"%8.1f ",(float)modelref.data.Rdisp.pper[i]);}
  fprintf(fout,"TLph ");
  for(i=0;i<modelref.data.Ldisp.npper;i++){
        fprintf(fout,"%8.1f ",(float)modelref.data.Ldisp.pper[i]);}
  fprintf(fout,"\n");
  
  
  //#pragma omp parallel default(none) shared(flag,c,fout,modelref,PREM,Nprem) private(modelnew,ddp,DRpvel,DRgvel,DLpvel,DLgvel,DRhvratio,DLdump,depth,i)
  {
  depth=0.;
  printf("model.nlayer=%d\n",modelref.laym0.nlayer);
  printf("jump threads=%d\n",omp_get_num_threads());
  //#pragma omp for schedule(dynamic,1)  
  for(i=0;i<modelref.laym0.nlayer;i++){
	modelnew=modelref;
	if(flag==1){
		ddp=modelref.laym0.vsv[i]*c;
		modelnew.laym0.vsv[i]+=ddp;}
	if(flag==2){
		ddp=modelref.laym0.vsh[i]*c;
		modelnew.laym0.vsh[i]+=ddp;}
	if(flag==3){
		ddp=modelref.laym0.vpv[i]*c;
		modelnew.laym0.vpv[i]+=ddp;}
	if(flag==4){
		ddp=modelref.laym0.vph[i]*c;
		modelnew.laym0.vph[i]+=ddp;}
	if(flag==5){
		ddp=modelref.laym0.eta[i]*c;
		modelnew.laym0.eta[i]+=ddp;}
	if(flag==6){
		ddp=modelref.laym0.rho[i]*c;
		modelnew.laym0.rho[i]+=ddp;}
	if(i==0){h1=0;h2=modelnew.laym0.thick[i];}
	else if (i==modelref.laym0.nlayer-1){h1=modelnew.laym0.thick[i-1];h2=0;}
	else{h1=h1=modelnew.laym0.thick[i-1];h2=modelnew.laym0.thick[i];}
	s=ddp*(h1+h2)/2.;
	compute_dispMineos(modelnew,PREM,Nprem,5,1,flag*100+i);
	compute_diff(modelnew.data.Rdisp,modelref.data.Rdisp,DRpvel,DRgvel,DRhvratio);
	compute_diff(modelnew.data.Ldisp,modelref.data.Ldisp,DLpvel,DLgvel,DLdump);
	write_sensitivity_one_dep(fout,depth,DRhvratio,s,0);	
	fprintf(fout,"TRph ");
	write_sensitivity_one_dep(fout,depth,DRpvel,s,1);	
	fprintf(fout,"TLph ");
	write_sensitivity_one_dep(fout,depth,DLpvel,s,2);	
	fprintf(fout,"\n");
	depth+=modelnew.laym0.thick[i];
	//if(depth>100.)break;
  }//for i
  }//pragma
  fclose(fout);
}

int main(int argc, char *argv[])
{
int npoint,Rsurflag,Lsurflag,Rmonoc,Lmonoc,PosAnic,iitercri1,iitercri2,ijumpcri1,ijumpcri2;
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
  isoflag=1; //isoflag==1: Vsv=Vsh, isoflag==0: Vsv!=Vsh
  Rsurflag=5; //surflag==1: open phase only. surfalg ==3 open phase and group, surflag==2: open group only; surflag=4: hv only; surflag=5:p+hv; surflag=6: g+hv; surflag=7: g+p+hv
  Lsurflag=1;
  AziampRsurflag=0;
  AziphiRsurflag=0;
  AziampLsurflag=0;
  AziphiLsurflag=0;
  //after I have changed the way the misfit is computed (compute_misfitDISP), the inpamp & inpphi becomes useless
  inpamp=0.25;//0.25; //the weight of the azi_aniso disp curve, amp part (0~1)
  inpphi=0.25;//0.25; //the weight of the azi_aniso disp curve, ang part (0-1)
  //the weight of iso dispersion curve is 1-inpamp-inpphi  
  iitercri1=100000;//100000;//12000 (mod1, 1cstlay)
  iitercri2=15000;
  ijumpcri1=10; //atoi(argv[10]); // set it to be the same as number_of_thread
  depcri1=20.0;
  depcri2=80.0;
  qpcri=900.;//900.;
  qscri=250.;
  Rmonoc=1;
  Lmonoc=1;
  PosAnic=1;
  flagreadLkernel=0;
  flagupdaterho=1;
  //Rvmono.push_back(0);
  Rvmono.push_back(1);
  //Rvmono.push_back(2);
  Lvmono.push_back(1);
  //Lvmono.push_back(2);
  //Rvgrad.push_back(0);
  Rvgrad.push_back(1);
  //Rvgrad.push_back(2);
  //Lvgrad.push_back(0);
  Lvgrad.push_back(1);
  Vposani.push_back(1);
  Vposani.push_back(2);
  Viso.push_back(0);
  Viso.push_back(1);
  Viso.push_back(2);
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
  flagreadVkernel=atoi(argv[9]);
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
    sprintf(str,"%s/disp.Ray_%.1f_%.1f.txt",Rphindir,lon,lat);
    Rdispnm.push_back(str);
    sprintf(str,"%s/HV.Ray_%.1f_%.1f.txt",Rphindir,lon,lat);
    Rdispnm.push_back(str);

    Ldispnm.clear();
    sprintf(str,"%s/disp.Lov_%.1f_%.1f.txt",Lphindir,lon,lat);
    Ldispnm.push_back(str);

    AziampRdispnm.clear();
    sprintf(str,"%s/aziamp.Ray_%.1f_%.1f.txt",Rphindir,lon,lat);
    AziampRdispnm.push_back(str);

    AziphiRdispnm.clear();
    sprintf(str,"%s/aziphi.Ray_%.1f_%.1f.txt",Rphindir,lon,lat);
    //sprintf(str,"%s/aziphi_restore_unc/aziphi_%.1f_%.1f.txt",Rphindir,lon,lat);//#############HEY TEMPERARY########## TEST
    AziphiRdispnm.push_back(str);

    AziampLdispnm.clear();
    sprintf(str,"%s/aziamp_%.1f_%.1f.txt",Lphindir,lon,lat);
    AziampLdispnm.push_back(str);

    AziphiLdispnm.clear();
    sprintf(str,"%s/aziphi_%.1f_%.1f.txt",Lphindir,lon,lat);
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
    
    readmodAniso(model0,modnm);// both m.g.LV/Rv are filled regardless of flags. (readin iso model)
    printf("finish read mod\n");
    printf("#########group[0].vsvvalue[0]=%g vshvalue[0]=%g\n",model0.groups[0].vsvvalue[0],model0.groups[0].vshvalue[0]);
    readpara(para0,fparanm);


    mod2para(model0,para0,para1);////fill both para.R/Lpara0 (they could be inequal if Rf*Lf>0, they are equal if Rf*Lf=0)

    checkParaModel(para1,model0,Viso);

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

    Vpara2Lovepara(para1,model0,flagupdaterho);//there is para2mod inside the code



    if(model0.flag==0){updatemodel(model0,flagupdaterho);}
    compute_dispMineos(model0,PREM,Nprem,Rsurflag,Lsurflag,0);

    printf("@@@ check, main, from initial model\n");
    for(j=0;j<model0.data.Rdisp.npper;j++)printf("  @@@ check, T=%g,vin=%g vold=%g\n",model0.data.Rdisp.pper[j],model0.data.Rdisp.pvelo[j],model0.data.Rdisp.pvel[j]);

  //----------------------------
  //---- compute the sensitivity of the HV disp to all five parameters (vsv,vsh,vpv,vph,eta)
  modelref = model0;
  float c = 0.01; //perturbation
  char fname[300];

  #pragma omp parallel default(none) shared(modelref,PREM,Nprem,c) private(fname,i)
  {
  printf("jump threads=%d\n",omp_get_num_threads());//--check---
  #pragma omp for schedule(dynamic,1)
  for(i=1;i<7;i++){
  printf("i=%d\n",i);
  if(i==1){sprintf(fname,"sens_vsv.txt");}
  else if(i==2){sprintf(fname,"sens_vsh.txt");}
  else if(i==3){sprintf(fname,"sens_vpv.txt");}
  else if(i==4){sprintf(fname,"sens_vph.txt");}
  else if(i==5){sprintf(fname,"sens_eta.txt");}
  else if(i==6){sprintf(fname,"sens_rho.txt");}
  
  compute_write_sensitivity(fname,modelref,c,i,PREM,Nprem);
  }//for i
  }//pragma

   printf("test-- END of code\n");

     
  }//while1
  return 1;
}
//main
