#define MAIN
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// this code is used to rotate the SAC from E-N to T-R.
// we are only checking the name looking like COR_sta1_sta2..., and not checking the existance of COR_sta2_sta1, SO the order of station name in sta.lst MUST be the same as that used in stack.c. Or, we'll miss some paths.
using namespace std;


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
        if((fsac=fopen(fname, "rb")) == NULL) return NULL;

        if ( !SHD ) SHD = &SAC_HEADER;

         fread(SHD,sizeof(SAC_HD),1,fsac);

         if ( SHD->npts > nmax ) {
           fprintf(stderr,
           "ATTENTION !!! in the file %s npts exceeds limit  %d",fname,nmax);
           SHD->npts = nmax;
         }

         fread(sig,sizeof(float),(int)(SHD->npts),fsac);

         fclose (fsac);

   /*-------------  calculate from t0  ----------------*/
   {
        int eh, em ,i;
        float fes;
        char koo[9];

        for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
        koo[8] = '\0';

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

        return SHD;
}




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];
 
   for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }
   fwrite(SHD,sizeof(SAC_HD),1,fsac);

   fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


   fclose (fsac);
}


SAC_HD shd_cor_EE;
SAC_HD shd_cor_EN;
SAC_HD shd_cor_NE;
SAC_HD shd_cor_NN;
SAC_HD shd_cor_temp;

int main(int na, char *arg[])
{ 
  //cout<<"GILL!!!"<<endl;
  if(na!=3)
    {
      cout<<"usage:lfroate stalist path"<<endl;
      return 0;
    }
  FILE *f1,*f2,*f3;
  int N=3000;
  double PI;  
  PI=atan(1.0)*4;
  double temp1,temp2;
  double cos1,cos2,sin1,sin2;
  char name[N][10];
  char name_EE[200],name_EN[200],name_NE[200],name_NN[200];
  char name_TT[200],name_RR[200],name_TR[200],name_RT[200];
  float sig_EE[100000],sig_EN[100000],sig_NE[100000],sig_NN[100000];
  float sig_RR[100000],sig_TT[100000],sig_TR[100000],sig_RT[100000];
  float sig_temp[100000],dump1,dump2;
  if(na!=3)
    {
      cout<<"usage:lfroate stalist path"<<endl;
	return 0;
    }
  f1=fopen(arg[1],"r");
  int i,j,k,ii,jj;
  for(i=0;i<N;i++)
    {
      if(fscanf(f1,"%s %g %g",name[i],&dump1,&dump2)==EOF) break;
      cout<<name[i]<<endl;
    }
  fclose(f1);

  cout<<"number of stations read "<<i<<endl;
  char buff[300];
  sprintf(buff,"if [ ! -d STACK_TT ]; then mkdir STACK_TT; fi");
  system(buff);
  sprintf(buff,"if [ ! -d STACK_RR ]; then mkdir STACK_RR; fi");
  system(buff);
  sprintf(buff,"if [ ! -d STACK_TR ]; then mkdir STACK_TR; fi");
  system(buff);
  sprintf(buff,"if [ ! -d STACK_RT ]; then mkdir STACK_RT; fi");
  system(buff);
  for(j=0;j<i-1;j++)
    {
      //cout<<"working on sta "<<j<<" "<<name[j]<<" "<<endl;
      continue;
      sprintf(buff,"if [ ! -d STACK_TT/%s ]; then mkdir STACK_TT/%s; fi",name[j],name[j]);
      system(buff);
      sprintf(buff,"if [ ! -d STACK_RR/%s ]; then mkdir STACK_RR/%s; fi",name[j],name[j]);
      system(buff);
      sprintf(buff,"if [ ! -d STACK_TR/%s ]; then mkdir STACK_TR/%s; fi",name[j],name[j]);
      system(buff);
      sprintf(buff,"if [ ! -d STACK_RT/%s ]; then mkdir STACK_RT/%s; fi",name[j],name[j]);
      system(buff);
      continue;
      for(k=j+1;k<i;k++)
	{
	  //cout<<"working on "<<name[j]<<" "<<name[k]<<endl;
	  sprintf(name_EE,"%s/%s/COR_%s_%s.SAC_EE",arg[2],name[j],name[j],name[k]);
	  sprintf(name_EN,"%s/%s/COR_%s_%s.SAC_EN",arg[2],name[j],name[j],name[k]);
	  sprintf(name_NE,"%s/%s/COR_%s_%s.SAC_NE",arg[2],name[j],name[j],name[k]);
	  sprintf(name_NN,"%s/%s/COR_%s_%s.SAC_NN",arg[2],name[j],name[j],name[k]);
	  sprintf(name_TT,"./STACK_TT/%s/COR_%s_%s.SAC_TT",name[j],name[j],name[k]);
	  sprintf(name_RR,"./STACK_RR/%s/COR_%s_%s.SAC_RR",name[j],name[j],name[k]);
	  sprintf(name_TR,"./STACK_TR/%s/COR_%s_%s.SAC_TR",name[j],name[j],name[k]);
	  sprintf(name_RT,"./STACK_RT/%s/COR_%s_%s.SAC_RT",name[j],name[j],name[k]);
	  if(read_sac (name_TT, sig_temp, &shd_cor_temp, 100000 ))
            {
        //      cout<<name_TT<<"  exist!!"<<endl;
              continue;
            }
	  if(!read_sac (name_EE, sig_temp, &shd_cor_temp, 100000 ))
            {
       //       cout<<name_EE<<" not exist!!"<<endl;
              continue;
            }
	  

	  f2=fopen("runsac.csh","w");
          fprintf(f2,"sac <<END\n");
          fprintf(f2,"r %s %s %s %s\n",name_EE,name_EN,name_NE,name_NN);
          fprintf(f2,"w %s %s %s %s\n",name_EE,name_EN,name_NE,name_NN);
          fprintf(f2,"q\nEND\n\n");
          fclose(f2);
          system("csh runsac.csh");
	  if(!read_sac (name_EE, sig_EE, &shd_cor_EE, 100000 ))
	    {
	      cout<<name_EE<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_EN, sig_EN, &shd_cor_EN, 100000 ))
	    {
	      cout<<name_EN<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_NE, sig_NE, &shd_cor_NE, 100000 ))
	    {
	      cout<<name_NE<<" not found!!"<<endl;
	      continue;
	    }	
	  if(!read_sac (name_NN, sig_NN, &shd_cor_NN, 100000 ))
	    {
	      cout<<name_NN<<" not found!!"<<endl;
	      continue;
	    }
	  temp1=fabs(cos((shd_cor_NN.cmpaz-shd_cor_EE.cmpaz)/180.0*PI));
	  temp2=fabs(cos((shd_cor_NN.user9-shd_cor_EE.user9)/180.0*PI)); //???? user9
//	  cout<<"user9EE,cmpazEE="<<shd_cor_EE.user9<<" "<<shd_cor_EE.cmpaz<<"\n user9 N,cmpaz N="<<shd_cor_NN.user9<<" "<<shd_cor_NN.cmpaz<<endl;
//	  cout<<temp1<<" "<<temp2<<endl;
//	  abort();
	  if(temp1>0.1||temp2>0.1)
	    {
	      cout<<"cmpaz incompatible!!"<<endl;
	      continue;
	    }
	  temp1=(shd_cor_NN.az-shd_cor_NN.cmpaz-180)/180.0*PI;
	  cos1=cos(temp1);
	  sin1=sin(temp1);
	  temp2=(shd_cor_NN.baz-shd_cor_NN.user9-180)/180.0*PI;
	  cos2=cos(temp2);
	  sin2=sin(temp2);
	  //cout<<temp1<<" "<<sin1<<" "<<cos1<<" "<<temp2<<" "<<sin2<<" "<<cos2<<endl;
	  for(ii=0;ii<shd_cor_NN.npts;ii++)
	    {
	      sig_TT[ii]=-1*cos1*cos2*sig_EE[ii]+cos1*sin2*sig_EN[ii]-sin1*sin2*sig_NN[ii]+sin1*cos2*sig_NE[ii];
	      sig_RR[ii]=-1*sin1*sin2*sig_EE[ii]-sin1*cos2*sig_EN[ii]-cos1*cos2*sig_NN[ii]-cos1*sin2*sig_NE[ii];
	      sig_TR[ii]=-1*cos1*sin2*sig_EE[ii]-cos1*cos2*sig_EN[ii]+sin1*cos2*sig_NN[ii]+sin1*sin2*sig_NE[ii];
              sig_RT[ii]=-1*sin1*cos2*sig_EE[ii]+sin1*sin2*sig_EN[ii]+cos1*sin2*sig_NN[ii]-cos1*cos2*sig_NE[ii];
	    }
	  write_sac (name_TT, sig_TT, &shd_cor_NN);
	  write_sac (name_RR, sig_RR, &shd_cor_NN);
	  write_sac (name_TR, sig_TR, &shd_cor_NN);
	  write_sac (name_RT, sig_RT, &shd_cor_NN);
	}
	  

    }

 
}
