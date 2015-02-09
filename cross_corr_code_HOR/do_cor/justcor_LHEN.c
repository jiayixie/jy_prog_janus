#define MAIN

#include <stdio.h>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the header file:
	cmpaz, the azi of the 1st component and 1st sta relative to N.
	user9, the azi of the 2ed component and 2ed sta relative to N.
	user8, number of days that have been stacked
*/
/* Finction prorotypes */

void dcommon_(int *len, float *amp,float *phase);
void dmultifft_(int *len,float *amp,float *phase, int *lag,float *seis_out, int *ns);

//void read_sac(char *name,char *stnam,float *stlat,float *stlon,int *n,float *sei);
void swapn(unsigned char *b, int N, int n);
//void write_cor(char *nam_f1,char *nam_corr,int *lag,float *corr,
//     char *stnam,float *stlat,float *stlon);


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




/*c/////////////////////////////////////////////////////////////////////////*/
 float sig_old[1000000];
 SAC_HD shdamp1, shdph1, shdamp2, shdph2, shd_cor;

/*--------------------------------------------------------------------------*/
int do_cor( int lag, int iii, int kkk )
/*--------------------------------------------------------------------------*/
{
  FILE *fst,*fev;
  if((fst=fopen("stalst","r"))==NULL)
    return 0;
  if((fev=fopen("evlst","r"))==NULL)
    return 0;
  char sta_name[2000][10];
  //  double sta_lat[2000],sta_lon[2000];
  char ev_name[40][20];
  int ist,iev,ine, jsta1, jsta2, k,count;
  ist=0;
  for(;;)
    {
      if(fscanf(fst,"%s",sta_name[ist])==EOF)
	break;
      ist++;
    }
  fclose(fst);
  iev=0;
  for(;;)
    {
      if(fscanf(fev,"%s",ev_name[iev])==EOF)
	break;
      iev++;
    }
  fclose(fev);

  int len,ns,i; 
  // int iii,jjj,kkk;
  char name1_N[200],name1_E[200];
  char name2_N[200],name2_E[200]; 
  char filename[200],tempfilename[200];
  char amp_sac[200], phase_sac[200];
  float amp[400000], phase[400000];
  float cor[400000];
  float seis_out[400000];
  double lat1,lat2,lon1,lon2,cmpaz;
  FILE *ff;
  char cp_am_ph[200];
  char buff[300];
  sprintf(buff,"if [ ! -d COR ]; then mkdir COR; fi");
  system(buff);

  for( jsta1 = 0; jsta1 < ist; jsta1++ ) 
    {  
      fprintf(stderr,"---jsta1 %d: %s--%d%d\n",jsta1,sta_name[jsta1],iii,kkk);
      sprintf(buff,"if [ ! -d COR/%s ]; then mkdir COR/%s; fi",sta_name[jsta1],sta_name[jsta1]);
      system(buff);
      for( jsta2 = (jsta1+1); jsta2 < ist; jsta2++ )
	{	
	  count=0;
	  if(iii==0)
	    return 0;
	  else if(iii==1)
	    {
	      sprintf(filename, "COR/%s/COR_%s_%s.SAC_N",
		      sta_name[jsta1], sta_name[jsta1], sta_name[jsta2]);
	      sprintf(tempfilename,"COR/%s/COR_%s_%s.SAC_N",
		      sta_name[jsta2], sta_name[jsta2], sta_name[jsta1]);
	    }
	  else if(iii==2)
	    {
	      sprintf(filename, "COR/%s/COR_%s_%s.SAC_E",
		      sta_name[jsta1], sta_name[jsta1], sta_name[jsta2]);
	      sprintf(tempfilename, "COR/%s/COR_%s_%s.SAC_E",
		      sta_name[jsta2], sta_name[jsta2], sta_name[jsta1]);
	    }
	  if(kkk==0)
	    return 0;
	  else if(kkk==1)
	    {strcat(filename,"N");strcat(tempfilename,"N");}
	  else
	    {strcat(filename,"E");strcat(tempfilename,"E");}
	  if(access(filename, F_OK) == 0)
	    {
	      fprintf(stderr,"%s exist!! skip\n",filename);
	      continue;
	    }
	  if(access(tempfilename,F_OK)==0)// added on Oct 24, 2012
	    {
		fprintf(stderr,"%s exist!! skip\n",tempfilename);
		continue;
	    }
//	  fprintf(stderr,"jsta1 %d jsta2 %d\n",jsta1,jsta2);
	  for( ine = 0; ine < iev; ine++ ) 
	    {
	      //	      if( sdb_N->rec[ine][jsta1].n > 0 && sdb_E->rec[ine][jsta1].n > 0)
	      //  
	      //  if(  sdb_N->rec[ine][jsta2].n > 0 && sdb_E->rec[ine][jsta2].n > 0)
	      if(iii==0)
		{
		  return 0;
		  //			  sprintf( amp_sac, "%s.am", sdb->rz[ine][jsta1].ft_fname );
		  //sprintf( phase_sac, "%s.ph", sdb->rz[ine][jsta1].ft_fname );
		}
	      else if(iii==1) 
		{
		  sprintf( amp_sac,   "%s/ft_%s.LHN.SAC.am", ev_name[ine],sta_name[jsta1] );
		  sprintf( phase_sac, "%s/ft_%s.LHN.SAC.ph", ev_name[ine],sta_name[jsta1] );
		}
	      else 
		{
		  sprintf( amp_sac,   "%s/ft_%s.LHE.SAC.am", ev_name[ine],sta_name[jsta1] );
		  sprintf( phase_sac, "%s/ft_%s.LHE.SAC.ph", ev_name[ine],sta_name[jsta1] );
		}
	      
	      
	      // read amp and phase files and read into common memory
	      if ( read_sac(amp_sac, amp, &shdamp1, 1000000 )==NULL )
		{
		  //		  fprintf( stderr,"file %s did not found\n", amp_sac );
		  continue;
		  //return 0;
		}
	      if (isnan(amp[0]))
		{fprintf(stderr,"%s is NaN file\n",amp_sac);
		 continue;}
	      if ( read_sac(phase_sac, phase, &shdph1, 1000000)== NULL )
		{
		  //		  fprintf( stderr,"file %s did not found\n", phase_sac );
		  continue;
		  //return 0;
		}
	      if (isnan(phase[0]))
		{fprintf(stderr,"%s is NaN file\n",phase_sac);
		 continue;}
	      len = shdamp1.npts;
	      
	      dcommon_( &len, amp, phase ); // reads amp and phase files into common memory
	      
	      // compute correlation
	      
	      if(kkk==0)
		{
		  return 0;
		  //			  sprintf( amp_sac, "%s.am", sdb->rz[ine][jsta2].ft_fname );
		  //sprintf( phase_sac, "%s.ph", sdb->rz[ine][jsta2].ft_fname );
		}
	      else if(kkk==1)
		{
		  sprintf( amp_sac,   "%s/ft_%s.LHN.SAC.am", ev_name[ine],sta_name[jsta2] );
		  sprintf( phase_sac, "%s/ft_%s.LHN.SAC.ph", ev_name[ine],sta_name[jsta2] );
		  // sprintf( amp_sac, "%s.am", sdb_N->rec[ine][jsta2].ft_fname );
		  //  sprintf( phase_sac, "%s.ph", sdb_N->rec[ine][jsta2].ft_fname  );
		}
	      else
		{
		  sprintf( amp_sac,   "%s/ft_%s.LHE.SAC.am", ev_name[ine],sta_name[jsta2] );
		  sprintf( phase_sac, "%s/ft_%s.LHE.SAC.ph", ev_name[ine],sta_name[jsta2] );
		  //	  sprintf( amp_sac, "%s.am", sdb_E->rec[ine][jsta2].ft_fname );
		  //  sprintf( phase_sac, "%s.ph", sdb_E->rec[ine][jsta2].ft_fname  );
		}
	      
	      
	      //fprintf(stderr,"file %s  %s\n", sdb->rec[ine][jsta1].ft_fname,sdb->rec[ine][jsta2].ft_fname );
	      // get array of floats for amp and phase of first signal
	      
	      if ( (read_sac(amp_sac, amp, &shdamp2, 100000) ==NULL)  ) 
		{
		  //		  fprintf(stderr,"file %s did not found\n", amp_sac );
		  continue;
		  //return 0;
		}
 	      if (isnan(amp[0]))
		{fprintf(stderr,"%s is NaN file\n",amp_sac);
		 continue;}
	      
	      if ( read_sac(phase_sac, phase, &shdph2, 100000)==NULL ) 
		{
		  //		  fprintf(stderr,"file %s did not found\n", phase_sac );
		  continue;
		  //		  return 0;
		}
 	      if (isnan(phase[0]))
		{fprintf(stderr,"%s is NaN file\n",phase_sac);
		 continue;}
	      
	      
	      len = shdamp2.npts;
		      
	      //if(!check_info(sdb_N, sdb_E, ine, jsta1, jsta2,iii,kkk )) 
	      //		{
	      //  fprintf(stderr,"files incompatible\n");
	      //  return 0;
	      //}
	      
	      //else
			  
	      dmultifft_(&len, amp, phase, &lag, seis_out,&ns);
	      cor[lag] = seis_out[0];
	      for( i = 1; i< (lag+1); i++)
		{ 
		  cor[lag-i] =  seis_out[i];
		  cor[lag+i] =  seis_out[ns-i];
		}
	      count++;
	      if(count!=1)
		{
                  
		  for(k = 0; k < (2*lag+1); k++) {

//			if (isnan(cor[k]))
				sig_old[k] += cor[k];
//			else
//				sig_old[k] += 0.;
		        }
		}
	      else
		{
		  //		  shdamp1.delta = 1;
		  lat1 = shdamp1.stla; 
		  lon1 = shdamp1.stlo;
		  lat2 = shdamp2.stla;
		  lon2 = shdamp2.stlo;
		  cmpaz= shdamp2.cmpaz;
		  for(k = 0; k < (2*lag+1); k++) sig_old[k] = cor[k];
		}
	      
	    }// for ine
	  if(count>0)
	    {
	      shdamp1.delta = 1;
	      shdamp1.user8=  count;
	      shdamp1.user9=  cmpaz;
	      //	      shdamp1.delta = sdb->rec[ine][jsta1].dt;
	      shdamp1.evla =  lat1;
	      shdamp1.evlo =  lon1;
	      shdamp1.stla =  lat2;
	      shdamp1.stlo =  lon2;
	      shdamp1.npts =  2*lag+1;
	      shdamp1.b    = -(lag)*shdamp1.delta;
	      //for(k = 0; k < (2*lag+1); k++) sig_old[k] = cor[k];
	      write_sac (filename, sig_old, &shdamp1);
	      //	      write_sac (filename, sig, &shd_cor );
		  
	    }
	      
	}  //loop over jsta2
    }  //loop over jsta1

  return 0;
}


//SAC_DB3 sdb;
//SAC_DB sdb_N,sdb_E;

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
int main (int na, char *arg[])
/*--------------------------------------------------------------------------*/
{
  FILE *ff;
  int iii,kkk;
  int ns1 = 0, ns2 = 0,lag;
  char str[600], filename[200], ch1[5],ch2[5];

  if ( na != 2 )
    {
      fprintf(stderr,"usage: corr_filtered_5000 lag \n");
      exit(1);
    }

  sscanf(arg[1],"%d", &lag );
  //  iii=atoi(arg[2]);
  //kkk=atoi(arg[3]);
  //  if((iii!=1 && iii!=2) || (kkk!=1 && kkk!=2) )
  // {
  //  fprintf(stderr,"iii or kkk wrong!!\n");
  //  return;
  //}
  if((ff = fopen("stalst","r"))==NULL)
    {
      fprintf(stderr,"stalst not exist!!\n");
      return;
    }
  
  fclose(ff);

  for(iii=1;iii<=2;iii++)
    {
	for(kkk=1;kkk<=2;kkk++)
	{
	printf("---doing componet[1--N 2--E] %d %d---------------\n",iii,kkk);
	do_cor(lag,iii,kkk);}  }


  fprintf(stderr, "finished correlations\n");

  // move COR/COR_STA1_STA2.SAC.prelim to COR/COR_STA1_STA2.SAC
  //  for ( ns2 = 1; ns2 < sdb_N.nst; ns2++ ) for ( ns1 = 0; ns1 < ns2; ns1++ ) 
  //{
  //  sprintf(filename, "COR/%s/COR_%s_%s.SAC.prelim_%s%s", sdb_N.st[ns1].name,sdb_N.st[ns1].name, sdb_N.st[ns2].name,ch1,ch2);
  //  sprintf(str, "mv COR/%s/COR_%s_%s.SAC.prelim_%s%s COR/%s/COR_%s_%s.SAC_%s%s",
  //      sdb_N.st[ns1].name, sdb_N.st[ns1].name, sdb_N.st[ns2].name,ch1,ch2, sdb_N.st[ns1].name,sdb_N.st[ns1].name, sdb_N.st[ns2].name,ch1,ch2);
  //  if(access(filename, F_OK) == 0) system(str);
  //}
  
  return 0;
}
