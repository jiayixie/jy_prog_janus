#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"
/* Function prorotypes */


void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,
              double *dt,int *n, float seis_in_E[],float seis_in_N[],float seis_out_E[],float seis_out_N[],
              float seis_outamp_E[],float seis_outamp_N[],
              float seis_outph_E[],float seis_outph_N[],int *ns,double *dom);

void swapn(unsigned char *b, int N, int n);




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
         fsac = fopen(fname, "rb");
	 if ( !fsac )
	   {
	     //fclose (fsac);
	     return NULL;
	   }
	 if ( !SHD ) SHD = &SAC_HEADER;

	 fread(SHD,sizeof(SAC_HD),1,fsac);

	 if ( SHD->npts > nmax )
	 {
	  fprintf(stderr,
	   "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);

	  SHD->npts = nmax;
	 }

	 fread(sig,sizeof(float),(int)(SHD->npts),fsac);

	fclose (fsac);

   /*-------------  calcule de t0  ----------------*/
   {
	int eh, em ,i;
	float fes;
	char koo[9];

	for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
	koo[8] = NULL;

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
 
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

	 fwrite(SHD,sizeof(SAC_HD),1,fsac);

	 fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


	 fclose (fsac);
}



/*c/////////////////////////////////////////////////////////////////////////*/
float sig[1000000];
SAC_HD shd1;
SAC_HD shd2; 


/*c/////////////////////////////////////////////////////////////////////////*/

int main (int argc, char *argv[])
{
static int n, ns,npow;
static double f1, f2, f3, f4, dt,dom;
static float seis_in_E[400000],seis_out_E[400000];
static float seis_in_N[400000],seis_out_N[400000];
static float seis_outamp_E[400000],seis_outph_E[400000];
static float seis_outamp_N[400000],seis_outph_N[400000];
double t1,t2,t3,t4;
char  name_E[160],name_N[160];
char  name_amp_E[160],name_ph_E[160];
char  name_amp_N[160],name_ph_N[160];
FILE  *in, *ff;
int   i, j, nn,iii;


  if( argc != 2) {
      printf("Usage: whiten_phamp  parameter_file\n");
      exit(-1);
  }

// open and read parameter file param.dat
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }
//by default, both E and N component exist! so, we don't check their existence here
  while((nn = fscanf(in,"%lf %lf %lf %lf %lf %d %s %s",&t1,&t2,&t3,&t4,
            &dt,&npow,name_E,name_N)) != EOF) 
    { // start main loop
      if(nn == 0 || nn != 8) break;
//      printf("Corners periods. Low: %f - %f, High: %f - %f\n",t1, t2, t3, t4);
//      printf("Step: %f, Cosine power: %d\n",dt, npow);
      
      // remove quotes from name
/*      j = 0;
      for(i = 0; i < strlen(name_E); i++) {
	if(name_E[i] == '\'' || name_E[i] == '\"') continue;
	name_E[j] = name_E[i]; j++;
      }
      name_E[j] = '\0';    
*/
/*      printf("%s\n",name_E);
      if(name_E[9]=='E')
	{
	  for(iii=0;iii<=20;iii++)
	    {
	      if(&name_E[iii]==NULL)
		break;
	      if(iii==9)
		name_N[iii]='N';
	      else
		name_N[iii]=name_E[iii];
	    }
	  printf("%s\n",name_N);
	}
      else if(name_E[10]=='E')
	{	
	  for(iii=0;iii<=20;iii++)
	    {
	      if(&name_E[iii]==NULL)
		break;
	      if(iii==10)
		name_N[iii]='N';
	      else
		name_N[iii]=name_E[iii];
	    }
	  printf("%s\n",name_N);
	}
      else
	{
	  printf("input name incompatible!!\n");
	  return 0;
	}
*/      // do running average before whitening
      
      ff = fopen("sac_one_cor","w");
      fprintf(ff,"/bin/rm -f smooth_E.sac smooth_N.sac\n");
      fprintf(ff,"sac << END\n");
      fprintf(ff,"r %s %s\n",name_E,name_N);
      fprintf(ff,"abs\n");
 //     fprintf(ff,"smooth mean h 128\n");
      fprintf(ff,"smooth mean h 128\n");
      fprintf(ff,"w aaa bbb\n");
      fprintf(ff,"r aaa\n");
      fprintf(ff,"subf bbb\n");
      fprintf(ff,"abs\n");
      fprintf(ff,"addf aaa\n");
      fprintf(ff,"addf bbb\n");
      fprintf(ff,"div 2\n");
      fprintf(ff,"w a1.avg\n");
      fprintf(ff,"r %s %s\n",name_E,name_N);
      fprintf(ff,"divf a1.avg\n");
      fprintf(ff,"w smooth_E.sac smooth_N.sac\n");
      fprintf(ff,"quit\n");
      fprintf(ff,"END\n");
      fclose(ff);
      system("sh sac_one_cor");
      
      // end of running average
      
      
      if ( !read_sac("smooth_E.sac", sig, &shd1, 1000000 ) )
	{
	  fprintf(stderr,"!!!####file smooth_E.sac did not found\n" );
	  continue;
	}
      n  = shd1.npts;
      dt = shd1.delta;
      
      for( i =0; i< n; i++)
	{  
	  seis_in_E[i] = sig[i];  
	  //     printf(" seis_in1  %d %f\n", i,sig[i]);
	}
      if ( !read_sac("smooth_N.sac", sig, &shd2, 1000000 ) )
	{
	  fprintf(stderr,"!!!####file smooth_N.sac did not found\n" );
	  continue;
       }
      for( i =0; i< n; i++)
	{  
	  seis_in_N[i] = sig[i];  
	  //     printf(" seis_in1  %d %f\n", i,sig[i]);
	}
//      printf(" Dt1= %f, Nsamples1= %d\n",dt, n);
      
      
      f1 = 1.0/t1; f2 = 1.0/t2; f3 = 1.0/t3; f4 = 1.0/t4;
      filter4_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_in_E,seis_in_N,seis_out_E,seis_out_N,seis_outamp_E,seis_outamp_N,seis_outph_E,seis_outph_N,&ns,&dom);
      
      shd1.delta=dt;shd2.delta=dt;
      shd1.npts=n;shd2.npts=n;
      write_sac(name_E,seis_out_E, &shd1);
      write_sac(name_N,seis_out_N, &shd2);
      
      strcpy(name_amp_E,name_E);
      strcpy(name_ph_E,name_E);
      strcat(name_amp_E,".am");
      strcat(name_ph_E, ".ph");
      strcpy(name_amp_N,name_N);
      strcpy(name_ph_N,name_N);
      strcat(name_amp_N,".am");
      strcat(name_ph_N, ".ph");
      shd1.npts = ns/2 + 1;
      shd1.delta = dom;
      shd1.b = 0;
      shd1.iftype = IXY;
      shd2.npts = ns/2 + 1;
      shd2.delta = dom;
      shd2.b = 0;
      shd2.iftype = IXY;
      write_sac(name_amp_E,seis_outamp_E, &shd1 );
      write_sac(name_ph_E, seis_outph_E,  &shd1 );
      
      write_sac(name_amp_N,seis_outamp_N, &shd2 );
      write_sac(name_ph_N, seis_outph_N,  &shd2 );
      
    }
  
  return 0;
}
