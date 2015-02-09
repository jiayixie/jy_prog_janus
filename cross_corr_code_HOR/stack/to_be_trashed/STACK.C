#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <iostream>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"

using namespace std;


SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig
, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    fprintf(stderr,"read_sac: Could not open %s\n", fname);
    return NULL;
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

//  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n", fname, (int)nmax);
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
     koo[8] = '\0';

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    fprintf(stderr,"write_sac: Could not open %s to write\n", fname);
  }
  else {

    if ( !SHD ) {
//      SHD = &SAC_HEADER;
    }

    SHD->iftype = (int)ITIME;
    SHD->leven = (int)TRUE;
    SHD->lovrok = (int)TRUE;
    SHD->internal4 = 6L;
    SHD->depmin = sig[0];
    SHD->depmax = sig[0];

    for ( i = 0; i < SHD->npts ; i++ ) {
      if ( SHD->depmin > sig[i] ) {
        SHD->depmin = sig[i];
      }
      if ( SHD->depmax < sig[i] ) {
        SHD->depmax = sig[i];
      }
    }

    fwrite(SHD,sizeof(SAC_HD),1,fsac);
    fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

    fclose (fsac);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////

int sac_add ( SAC_HD *shd, float *fsig, char *filename, int flag ) {
  int ns,i;
  float tsig[1000000];
  SAC_HD shd1;
  if (read_sac(filename, tsig, &shd1, 10001) == NULL) {
    fprintf(stderr,"cannot open file %s\n",filename);
    return 0;
    }
  ns = shd1.npts;
  if (isnan(tsig[0]))
  {
   fprintf(stderr,"file %s is NaN\n",filename);
   return 0;
  }
  if (flag == 1) {
    for (i=1;i<=3000;i++) {
	fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 + i];
	fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 - i];
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2];
    *shd = shd1;
    }
  else if (flag == 0 ) {
//    fprintf (stderr,"reverse!!!!!\n");
    for (i=1;i<=3000;i++) {
        fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 - i];
        fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 + i];
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2];    
    *shd=shd1; // ---- modified on Jan 17, 2013. here the stlo and evlo are not correct, but it doesn't matter, it will be changed right before we write the sac. The main purpose is transfer the user8,user9,and cmpaz to the header file.
    }
  return 1;
  }

////////////////////////////////////////////////////////////////////////////////////////////

#define NSTA 5000

int main (int argn, char *argv[]) {
 
  int i,j,k,nch,nsta,ndir,ii,jj,jjj,jjjj,flag;
  char staname[NSTA][10],outname[100],dirlist[200][300],tname1[300],tname2[300],outfname[300],outdir1[300],tstr[300];
  char tname3[300],tname4[300],ch[4][10];
  float stalat[NSTA],stalon[NSTA],fsig[6001];
  FILE *ff1;  

  SAC_HD shd,shd1;


  if (argn != 4) {
    fprintf (stderr,"input [station.lst] [dir.lst] [out.dir]\n");
    return 0;
    }

  if ((ff1 = fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[1]);
    return 0;
    }

  for (i=0;;i++) {
    if (fscanf(ff1,"%s %g %g", &staname[i], &stalon[i], &stalat[i]) != 3) break;
    }
  nsta = i;
  fclose(ff1);  
  sprintf(ch[0],"NN");
  sprintf(ch[1],"EE");
  sprintf(ch[2],"NE");
  sprintf(ch[3],"EN");
  if ((ff1 = fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[1]);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%s", &dirlist[i]) != 1) break;
    }
  ndir = i;
  fclose(ff1); 
  
  fprintf (stderr,"%d %d\n",nsta,ndir);

  jj = 0;
  k = 0;

  for (i=0;i<nsta;i++) {
     fprintf(stderr,"Begin--Nsta1= %d, %s\n",i,staname[i]);
     sprintf(outdir1,"%s/%s\0",argv[3],staname[i]);
     if (access(outdir1,F_OK) != 1) {
//	fprintf (stderr,"open directory %s\n",outdir1);
	sprintf(tstr,"mkdir %s\0" , outdir1);
	system(tstr);
	}

     for (nch=0;nch<4;nch++){
     jjj = 0;
     jjjj = 0; // #of stacked station pair for all evtdir
     for (j=i+1;j<nsta;j++) {
       sprintf(outname,"COR_%s_%s.SAC_%s",staname[i],staname[j],ch[nch]);
       sprintf(outfname,"%s/%s/COR_%s_%s.SAC_%s",argv[3],staname[i],staname[i],staname[j],ch[nch]);
       flag = 0;

       jj = jj +1 ;
       jjj = jjj + 1;
//       if (fmod(jj,50) == 0) fprintf(stderr,"%d %d %s %d\n",jj,jjj,outname,jjjj);


       if (access(outfname,F_OK) == 0) continue;

//       jjjj = 0;
       for (k=0;k<ndir;k++) {
         sprintf(tname1,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[i],staname[i],staname[j],ch[nch]);
         sprintf(tname2,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[j],staname[j],staname[i],ch[nch]);
	 sprintf(tname3,"%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[i],staname[j],ch[nch]);
	 sprintf(tname4,"%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[j],staname[i],ch[nch]);
	 if (access(tname1,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname1);
	    if (flag == 0) {
		shd = sac_null;
		for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
		if ( (sac_add (&shd,fsig,tname1,1))==1) 
		{flag = 1;
		jjjj ++;}
		continue;
		}
	    else {
		if((sac_add(&shd,fsig,tname1,1))==1)
		{jjjj ++;}
		continue;
		}
	    }
	 if (access(tname2,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname2);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                if((sac_add (&shd,fsig,tname2,0))==1)
                {flag = 1;
		jjjj ++;}
		continue;
                }
            else {
                if((sac_add(&shd,fsig,tname2,0))==1)
		{jjjj ++;}
		continue;
                }
            }

/*	 if (access(tname3,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname3);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                sac_add (&shd,fsig,tname3,1);
                flag = 1;
		jjjj ++;
		continue;
                }
            else {
                sac_add(&shd,fsig,tname3,1);
		jjjj ++;
		continue;
                }
            }

	  if (access(tname4,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname4);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                sac_add (&shd,fsig,tname4,0);
                flag = 1;
		jjjj ++;
                }
            else {
                sac_add(&shd,fsig,tname4,0);
		jjjj ++;
                }
            }
*/
         }//k ndir
       if (flag==1) {
	 shd.evlo = stalon[i];
	 shd.evla = stalat[i];
	 sprintf(shd.kevnm,"%s\0",staname[i]);
	 shd.stlo = stalon[j];
	 shd.stla = stalat[j];
	 sprintf(shd.kstnm,"%s\0",staname[j]);
	 shd.npts = 6001;
	 shd.b = -3000;
	 shd.delta = 1;
	 shd.e = shd.b+shd.npts*shd.delta-1;
	 shd.lcalda = 1;
	 shd.dist = 100.;
	 write_sac(outfname,fsig,&shd);
	 }
	}//j
	fprintf(stderr,"====Num of cor added= %d, out of %d cordir %d nsta, ch %s====\n",jjjj,ndir,nsta,ch[nch]);
     }//nch
     }//i

  return 1;
  }
