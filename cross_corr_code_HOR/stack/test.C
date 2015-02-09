#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <iostream>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"

// this version, we add weight to the stack_dir
// this version, count month stacked, and save the info in user4
// this version, handles the different station format used in the input SAC file names.
// this is used to make up a bug in previous version! the header of some files are missed, for those missed_HD(user8,user9,cmpaz), redo stack.
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

int sac_add ( SAC_HD *shd, float *fsig, char *filename, int flag, float wei ) {
  int ns,i;
  float tsig[1000000],Ncount;
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

  if(shd->user4<-12000.){Ncount=0;}
  else{Ncount=shd->user4;}
  if(shd1.user4<-12000){Ncount++;}
  else{Ncount=Ncount+shd1.user4;}

  if (flag == 1) {
    for (i=1;i<=3000;i++) {
	fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 + i]*wei;
	fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 - i]*wei;
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2]*wei;
    *shd = shd1;
    }
  else if (flag == 0 ) {
//    fprintf (stderr,"reverse!!!!!\n");
    for (i=1;i<=3000;i++) {
        fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 - i]*wei;
        fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 + i]*wei;
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2]*wei;    
    *shd = shd1;
    }

  shd->user4=Ncount;

  return 1;
  }

////////////////////////////////////////////////////////////////////////////////////////////

#define NSTA 5000

int main (int argn, char *argv[]) {
 
  int i,j,k,nch,nsta,ndir,ii,jj,jjj,jjjj,flag;
  char staname[NSTA][10],staname2[NSTA][10],outname[100],dirlist[200][300],tname1[300],tname2[300],outfname[300],outdir1[300],tstr[300];
  char tname3[300],tname4[300],ch[4][10];
  float stalat[NSTA],stalon[NSTA],fsig[6001],kflag[200],wei;
  FILE *ff1;  
  SAC_HD shd,shd1;

  float ftsig[1000000];
  SAC_HD ftshd;
  FILE *fwrong;
  char fwrongnm[300];

  if (argn != 4) {
    fprintf (stderr,"input [station.lst] [dir.lst] [out.dir]\n");
    return 0;
    }

  if ((ff1 = fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[1]);
    return 0;
    }

  for (i=0;;i++) {
    if (fscanf(ff1,"%s %g %g %s", &staname[i], &stalon[i], &stalat[i],&staname2[i]) != 4) {break;}
    }
  nsta = i;
  fclose(ff1);  
  sprintf(ch[0],"NN");
  sprintf(ch[1],"EE");
  sprintf(ch[2],"NE");
  sprintf(ch[3],"EN");
  if ((ff1 = fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[2]);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%s %g", &dirlist[i], &kflag[i]) != 2) break;
    }
  ndir = i;
  fclose(ff1); 
  
  fprintf (stderr,"%d %d\n",nsta,ndir);

  jj = 0;
  k = 0;

  sprintf(fwrongnm,"problematic_stack.txt");
  fwrong=fopen(fwrongnm,"w");
 
  for (i=0;i<nsta;i++) {
     fprintf(stderr,"Begin--Nsta1= %d, %s %s\n",i,staname[i],staname2[i]);
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

 


       if (access(outfname,F_OK) == 0 ) {
  //e/if (read_sac(filename, tsig, &shd1, 10001) == NULL)
		ftshd=sac_null;
		read_sac(outfname,ftsig,&ftshd,10001);
		//cout<<outfname<<"--- user9="<<ftshd.user9<<endl;
		if (ftshd.user9>-12340)continue;
		fprintf(fwrong,"rewriteHD i=%5d nch=%5d j=%5d  %s\n",i,nch,j,outfname);
	}
	else{
		//printf("############## strange!! i=%5d nch=%5d j=%5d  fname=%s ####################\n",i,nch,j,outfname);
		fprintf(fwrong,"not_exist i=%5d nch=%5d j=%5d  %s\n",i,nch,j,outfname);
	}
	

//       jjjj = 0;
       for (k=0;k<ndir;k++) { // just do some change here to make it fit for the stack_cv_jy case.
				// the station name cv and jy used are in different format; so when read in the staname, read it in two format, and if k > X, use the second format. This change should be simple
	 if (k>=1){
	  sprintf(tname1,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname2[i],staname2[i],staname2[j],ch[nch]);
         sprintf(tname2,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname2[j],staname2[j],staname2[i],ch[nch]);
	 }// if k>1
 	 else{
         sprintf(tname1,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[i],staname[i],staname[j],ch[nch]);
         sprintf(tname2,"%s/%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[j],staname[j],staname[i],ch[nch]);
	 } //else
	 //sprintf(tname3,"%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[i],staname[j],ch[nch]);
	 //sprintf(tname4,"%s/COR_%s_%s.SAC_%s\0",dirlist[k],staname[j],staname[i],ch[nch]);
	 wei = kflag[k];
	 //printf("i=%d j=%d k=%d --- \ntname1 %s\ntname2 %s\n",i,j,k,tname1,tname2);
	 if (access(tname1,F_OK) == 0) {
	    //fprintf (stderr,"yes!! file %s\n",tname1);
	    if (flag == 0) {
		shd = sac_null;
		for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
		if ( (sac_add (&shd,fsig,tname1,1,wei))==1) 
		{flag = 1;
		jjjj ++;}
		//printf("jjjj added = %d\n",jjjj);
		continue;
		}
	    else {
		if((sac_add(&shd,fsig,tname1,1,wei))==1)
		{jjjj ++;}
		//printf("jjjj added = %d\n",jjjj);
		continue;
		}
	    }
	 if (access(tname2,F_OK) == 0) {
	    //fprintf (stderr,"yes!! file %s\n",tname2);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                if((sac_add (&shd,fsig,tname2,0,wei))==1)
                {flag = 1;
		jjjj ++;}
		//printf("jjjj added = %d\n",jjjj);
		continue;
                }
            else {
                if((sac_add(&shd,fsig,tname2,0,wei))==1)
		{jjjj ++;}
		//printf("jjjj added = %d\n",jjjj);
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
	fprintf(stderr,"====Num of cor added= %d, out of %d cordir %d nsta, user4=%g, ch %s====\n",jjjj,ndir,nsta,shd.user4,ch[nch]);
     }//nch
     }//i
  fclose(fwrong);
  return 1;
  }
