#define MAIN

// this version will do 
// 1. check the station name and event name, and make sure that the location of each station is 
// compatitive with that in station.lst file
//
// this version also do
// 2. calculates the distance of 2 stations using Vincenty's formulea
//
// 3. cut it into 2 parts (both of which would be [nnpts] points long
// )using C and stack them together.
//
// this version is updated by cv


#include "/home/jiayi/progs/jy/HEAD/head_noise/64_koftan.h"
/*#include "../SAC_FROM_SEED/sac_db.h"*/
//#include "/home/jiayi/progs/jy/HEAD/head_noise/64_sac_db.h"
#include <math.h>
//#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"
//#include <iostream>
//#include <unistd.h>

//using namespace std;


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
        if((fsac = fopen(fname, "rb")) == NULL) {
          printf("could not open sac file %s to read\n",fname);
          return NULL;
        }

        if ( !fsac )
        {
          /*fprintf(stderr,"file %s not find\n", fname);*/
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
        koo[8] = 0;

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
        if((fsac = fopen(fname, "wb"))==NULL) {
           printf("could not open sac file to write\n");
           exit(1);
        }

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

int check_station_loc(float *lat1, float *lon1,char *name ){
FILE *fstalist,*RECORD;
char tempname[20];
double templat,templon;
float cvlat,cvlon;
cvlat = *lat1;
cvlon = *lon1;
char cvbuff[300];
templat=0.0;
templon=0.0;
int flag=0;
//fprintf(stderr,"cvlat: %g cvlon: %g\n",cvlat,cvlon);
//fclose(fstalist);
//
if ((fstalist = fopen("station.lst","r"))==NULL) {
   fprintf(stderr,"cannot open station.lst file \n");
   return 1;}


rewind(fstalist);
//fclose(fstalist);
//fstalist = fopen("station.lst","r");

RECORD = fopen("record_change_location","a");
//fseek(RECORD,0L,2);
//abort();
for(;;) {
    flag = 0;
//  if(fgets(&(cvbuff[0]),300,fstalist)==NULL) break;
//  fprintf(stderr,"cvbuff: %s",cvbuff);
//  sscanf(cvbuff, "%s %g %g",&(tempname[0]),&templon,&templat); 

    if(fscanf(fstalist,"%s %lf %lf",&(tempname[0]),&templon,&templat)==EOF) 
    {
//     fprintf(stderr,"end of stationlist file and flag is : %d\n",flag);
//     abort();
     break;
    }
//    fprintf(stderr, "%s\n",name);
//  fprintf(stderr,"%s %lf %lf\n",tempname,templon,templat);
//  if(templon < 0.1) abort();
//  fprintf(stderr,"%s %s",tempname,name);
  if (strcmp(tempname,name)==0) {
//    fprintf(stderr,"sta1: %s sta2: %s templat: %g cvlat: %g  templon: %g cvlon: %g\n",tempname,name,templat,cvlat,templon,cvlon);
    //abort();
    if (fabs(templat - cvlat)>0.001) {
       fprintf(RECORD,"sta: %s  original lat: %g  now lat: %g \n",name,cvlat,templat);
       printf("sta1: %s sta2: %s templat: %g cvlat: %g  templon: %g cvlon: %g\n",tempname,name,templat,cvlat,templon,cvlon);
       cvlat = templat;
       flag =1;
       //abort();
       }
    if (fabs(templon - cvlon)> 0.001 ) {
       fprintf(RECORD,"sta: %s  original lon: %g  now lon: %g \n",name,cvlon,templon);
       fprintf(stderr,"sta1: %s sta2: %s templat: %g cvlat: %g  templon: %g cvlon: %g\n",tempname,name,templat,cvlat,templon,cvlon);
       cvlon = templon;
       flag = 1;
       }
    break; 
    }
  }
fclose(RECORD);
*lat1 =cvlat;
*lon1 =cvlon;
fclose(fstalist);
if (flag == 1) fprintf(stderr,"changed lat and lon: %g %g for %s\n",*lat1,*lon1,name);
//if (flag == 0) abort();
return 1;

}





double get_dist2(double lat1, double lon1, double lat2, double lon2,double *az, double *baz) {

double pi;
pi = 4.0*atan(1.0);
double cva = 6378.137;
double cvb = 6356.7523142;
double f = 1/298.257223563;
double L = 0.00;
L = lon2-lon1;
//cout<<"L "<<L<<endl;
//if (L > 180.000)  L =360.000000 - L;
//if (L < -180.000) L =  360.000 - abs(L);
//L = fabs(L);
//cout<<"L: "<<L<<endl;
double U1 = 0;
U1 = atan((1-f)*tan(lat1/180*pi));
double U2 = 0;
U2 = atan((1-f)*tan(lat2/180*pi));
//cout<<"U1 "<<U1<<"U2 "<<U2<<" "<<sin(U1)<<" "<<sin(U2)<<" "<<cos(U1)<<" "<<cos(U2)<<endl;
double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
L = L*pi/180;


double numda = L;
numda1 = numda;
//cout<<"numda "<<numda<<"cos numda "<<cos(numda)<<endl;
do {
numda = numda1;
cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) ); // cv1 sin(quan)
cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
//cout<<"cv1 cv2: "<<cv1<<" "<<cv2<<endl;
cv = atan2(cv1,cv2);
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
cv4 = 1 - cv3*cv3;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
//cout<<"numda1 "<<numda1<<endl;
} while (fabs(numda - numda1) > 0.0000000001);

double mius, cvA, cvB, deltacv,s;
mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
//cout<<"cvA "<<cvA<<endl;
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
//cout<<"delatacv "<<deltacv<<"cvb "<<cvb<<"cv "<<cv<<endl;

s = cvb * cvA *(cv - deltacv);
//fprintf(stderr,"check : %g %g  %g %g \n", cos(U2)*sin(numda),(cos(U1)*sin(U1)-sin(U1)*cos(U2)*cos(numda)), cos(U1)*sin(numda), (-sin(U1)*cos(U2)+cos(U1)*sin(U2)*cos(numda)));
//abort();

//*az = atan2(cos(U2)*sin(numda),(cos(U1)*sin(U1)-sin(U1)*cos(U2)*cos(numda)) );
//*baz = atan2(cos(U1)*sin(numda),(-sin(U1)*cos(U2)+cos(U1)*sin(U2)*cos(numda)));
double temp1, temp2;
temp1= atan2(cos(U2)*sin(numda),(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(numda)) )*180/pi;
temp2= atan2(cos(U1)*sin(numda),(-sin(U1)*cos(U2) +cos(U1)*sin(U2)*cos(numda)))*180/pi;
//az=atan2(0.227546,0.0439109);
//baz=atan2(0.241236,0.0611887);
//if (temp1<0) temp1=temp1+360;
//if (temp2<0) temp2=temp2+360;
if (temp1<0) temp1=temp1+360;
temp2=temp2+180;
*az = temp1;
*baz= temp2;
//fprintf(stderr,"az:%g  baz:%g numda: %g, U1:%g U2:%g: s:%g %g %g %g\n",*az,*baz,numda,U1,U2,s,temp1,temp2,L);
//abort();
return s;
}



/*c/////////////////////////////////////////////////////////////*/
/* one_pair is the main routine to run the other subroutines.  */
#define SLEN 200000
float sig0[SLEN];
char fname[300];
/*--------------------------------------------------------------*/
int one_pair(char *name1, char *staname1, char *staname2, int nnpts )
/*--------------------------------------------------------------*/
{
  float b, dt, tmax, Umin = 1.;
  double dist;
  float Evla, Evlo, Stla, Stlo;
  char Statname[10], Evname[10];
  char ch1;
  FILE *fp1;
  float *pointer1, *pointer2; 

  char cvread[300],cvsac1[300],cvsac2[300],cvsacout[300];
  SAC_HD cvSAC_HEADER, cvSAC_HEADER1, cvSAC_HEADER2,cvSAC_HEADERout;
  float cvsig[SLEN], cvsig1[nnpts+1],cvsig2[nnpts+1],cvsigout[nnpts+1];
  double cvaz,cvbaz;

  sprintf(fname, "%s", name1);

  //printf("fname is %s\n\n", fname);

  /*---------------- reading sac file  -------------------*/
  if ( read_sac (fname, sig0, &SAC_HEADER, SLEN) == NULL )
    {
      fprintf(stderr,"file %s not found\n", fname);
      return 0;
    }

  b = SAC_HEADER.b;
  dt = SAC_HEADER.delta;
  dist = SAC_HEADER.dist;
  Evla = SAC_HEADER.evla;
  Evlo = SAC_HEADER.evlo;
  Stla = SAC_HEADER.stla;
  Stlo = SAC_HEADER.stlo;
  //fprintf(stderr,"Evla: %g Evlo: %g NAME: %s\n",Evla,Evlo,staname1);
  pointer1=&Evla;
  pointer2=&Evlo;
//  check_station_loc(pointer1,pointer2,staname1);
  Evla = *pointer1;
  Evlo = *pointer2;
  pointer1=&Stla;
  pointer2=&Stlo;
//  check_station_loc(pointer1,pointer2,staname2);
  Stla = *pointer1;
  Stlo = *pointer2;
  //fprintf(stderr,"Evla: %g Evlo: %g NAME: %s\n",Evla,Evlo,staname1);
  //abort();  


//  SAC_HEADER.evla=Evla;
//  SAC_HEADER.evlo=Evlo;
//  SAC_HEADER.stla=Stla;
//  SAC_HEADER.stlo=Stlo;
//  fprintf (stderr, "last check : %g %g %g %g\n",SAC_HEADER.evla, SAC_HEADER.evlo, SAC_HEADER.stla, SAC_HEADER.stlo);


  dist = get_dist2(Evla,Evlo,Stla,Stlo,&cvaz,&cvbaz);
//  fprintf(stderr,"dist: %g and cvaz: %g cvbaz: %g\n",dist,cvaz,cvbaz);
//  if((fp1 = fopen("change.csh", "w"))==NULL) {
//    printf("cannot open change.csh.\n");
//    exit(1);
//  }
  sprintf(Statname, staname1);
  strcpy(SAC_HEADER.kstnm, Statname);

  //printf("sta1 %s sta2 %s\n", staname1, staname2);
  //printf("sta1len %d sta2len %d\n", strlen(staname1), strlen(staname2));
  
  int cvleng,cvleng1;
  cvleng = SAC_HEADER.npts;
  if (fmod(cvleng, 2)!=0)
     cvleng1 = (cvleng-1)/2;
  else 
     cvleng1 = cvleng/2-1;
  
  read_sac(fname,cvsig,&cvSAC_HEADER,SLEN);

  int cvi=0;
  
  for (cvi=1;cvi<=nnpts;cvi++)  {
     cvsig1[cvi]=cvsig[nnpts+cvi];
     cvsig2[cvi]=cvsig[nnpts-cvi];
     cvsigout[cvi]=(cvsig1[cvi]+cvsig2[cvi]);
   }
  cvsig1[0]=cvsig[nnpts];
  cvsig2[0]=cvsig[nnpts];
  cvsigout[0]=cvsig[nnpts];

  strcpy(cvSAC_HEADER.kevnm,staname1);
  strcpy(cvSAC_HEADER.kstnm,staname2);
  cvSAC_HEADER.dist = dist;

  cvSAC_HEADER.evla=Evla;
  cvSAC_HEADER.evlo=Evlo;
  cvSAC_HEADER.stla=Stla;
  cvSAC_HEADER.stlo=Stlo;

  cvSAC_HEADER.az = (float)(cvaz);
  cvSAC_HEADER.baz = (float)(cvbaz);




  memcpy(&cvSAC_HEADERout,&cvSAC_HEADER,sizeof(struct sac));
  memcpy(&cvSAC_HEADER1,&cvSAC_HEADER,sizeof(struct sac));
  memcpy(&cvSAC_HEADER2,&cvSAC_HEADER,sizeof(struct sac));

  cvSAC_HEADER1.o = 0;
  cvSAC_HEADER1.b = 0;
  cvSAC_HEADER1.e = nnpts;
  cvSAC_HEADER1.npts = nnpts+1;

  cvSAC_HEADER2.o = 0;
  cvSAC_HEADER2.b = 0;
  cvSAC_HEADER2.e = nnpts;
  cvSAC_HEADER2.npts = nnpts+1;

  cvSAC_HEADERout.o = 0;
  cvSAC_HEADERout.b = 0;
  cvSAC_HEADERout.e = nnpts;
  cvSAC_HEADERout.npts = nnpts+1;


  strcpy(cvsac1,fname);
  strcpy(cvsac2,fname);
  strcpy(cvsacout,fname);
  strcat(cvsac1,"_p\0");
  strcat(cvsac2,"_n\0"); 
  strcat(cvsacout,"_s\0"); 

  
  write_sac(cvsac1,&cvsig1[0],&cvSAC_HEADER1);
  write_sac(cvsac2,&cvsig2[0],&cvSAC_HEADER2);  
  write_sac(cvsacout,&cvsigout[0],&cvSAC_HEADERout);
  /*fprintf(fp1, "sac << END\n");
  fprintf(fp1, "r %s\n", fname);
  fprintf(fp1, "ch kevnm %s\n", staname1);
  fprintf(fp1, "w %s\n", fname);
  fprintf(fp1, "cut 0 5000\n");
  fprintf(fp1, "r %s\n", fname);
  fprintf(fp1, "ch o 0\n");
  fprintf(fp1, "ch kevnm %s\n", staname1);
  fprintf(fp1, "w %s_p\n", fname);

  fprintf(fp1, "cut -5000 0\n");
  fprintf(fp1, "r %s\n", fname);
  fprintf(fp1, "reverse\n");
  fprintf(fp1, "ch b 0\n");
  fprintf(fp1, "ch e 5000\n");
  fprintf(fp1, "ch o 0\n");
  fprintf(fp1, "w %s_n\n", fname);

  fprintf(fp1, "cut off\n");
  fprintf(fp1, "addf %s_p\n", fname);
  fprintf(fp1, "ch o 0\n");
  fprintf(fp1, "div 2\nw %s_s\n", fname);
  fprintf(fp1, "END\n\n");

  fclose(fp1);
  system("csh change.csh");
  */
  return 1;
}


/*c/////////////////////////////////////////////////////////////*/
//char fname[300], str[300];
char str[300];
//SAC_DB sdb;
/*--------------------------------------------------------------*/
int main (int arg, char *argv[])
/*--------------------------------------------------------------*/
{
  float stla, stlo, evla, evlo, dist;
  FILE *ff;
//  int ns1 = 0, ns2 = 4, n = 0, len = 0, i;
//  char filename[29], fname[29], shortname[26];
//  char staname1[10], staname2[10], dump[500];
//  char ch, ch1;
  char dir[100],sacnm[200],nt1[10],sta1[10],nt2[10],sta2[10],ntst1[20],ntst2[20];
  char name[300];
  if(arg!=3) {
  printf("please enter : 1]dirsac_stack_rotate 2]the filelist which contains the names of sac file\n");
  exit(1);
  }
  sprintf(dir,"%s",argv[1]);
  if((ff = fopen(argv[2], "r"))==NULL) {
    printf("cannot open filelist.\n");
    exit(1);
  }
  for (;;)
  {
    if ((fscanf(ff,"%s",sacnm))==EOF)
	break;
    sprintf(name,"%s/%s",dir,sacnm);
    sscanf(sacnm, "COR_%[A-Z,a-z,0-9].%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].%[A-Z,a-z,0-9].", nt1,sta1,nt2,sta2);
    sprintf(ntst1,"%s.%s",nt1,sta1);
    sprintf(ntst2,"%s.%s",nt2,sta2);
   // printf("===%s,%s  %s\n",ntst1,ntst2,name);
    if ( !one_pair(name, ntst1, ntst2, 3000) ) continue;
  }
//  abort(); 
/*  if ((ch = fgetc(ff))==EOF) exit(1);
  do
   {
    fgets(dump, 4, ff);
    i = 0;
    ch = fgetc(ff);
    do{
      staname1[i] = ch;
      ch = fgetc(ff);
      i++;
    }while(ch!='_');
    staname1[i] = '\0';
    i = 0;
    ch = fgetc(ff);
    do{
      staname2[i] = ch;
      ch = fgetc(ff);
      i++;
    }while(ch!='.');
    staname2[i] = '\0';
    fgets(dump, 200, ff);
//    abort();
    sprintf(filename, "COR_%s_%s.SAC", staname1, staname2);

    if ( !one_pair(filename, staname1, staname2, 3000) ) continue;

    if (fmod(n,100) == 0) {
      fprintf (stderr,"finish %d %s\n",n,filename);
      }

    n++;
    ch = fgetc(ff);
   } while(!feof(ff));
*/
   fclose(ff);

}
