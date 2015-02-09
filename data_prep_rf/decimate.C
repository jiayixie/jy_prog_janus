// this is used to decimate the SAC files, since the RFiter cannot handle sac files with too many points
// npts<MAXPTS/2 (MAXPTS = 8192 in iterdeconfd.f)
// Also, mv SAC data with problems( 1]incompelete 3 comp, 2]with no header info, 3]cannot be decimated/ too long)
#define MAIN
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_mysac.h"
#include "/home/jiayi/progs/jy/HEAD/head_noise/64_sac_db.h"
using namespace std;

SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
	fsac = fopen(fname, "rb");
	if ( !fsac )
	{
	 fclose (fsac);
	 return NULL;
	}

	if ( !SHD ) SHD = &SAC_HEADER;

	 fread(SHD,sizeof(SAC_HD),1,fsac);

	 if ( SHD->npts > nmax )
	 {
	   /*fprintf(stderr,
	     "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);*/

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
/*c////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------*/
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
/////////////////////////////////////////////////////////////////////////

char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[4096];
//  char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;

  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';

  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));

  return buffer;
}

bool FileNotExists(const char* filename) 
{
  struct stat info;
  int ret = -1;
 
  //get the file attributes
  ret = stat(filename, &info);
  if(ret == 0) 
    {
    //stat() is able to get the file attributes,
    //so the file obviously exists
    return false;
    } 
  else 
    {
    //stat() is not able to get the file attributes,
    //so the file obviously does not exist or
    //more capabilities is required
    return true;
    }
}

#define NPTSMAX 15000
#define FILEMAX 5000
#define RFNPTS 4096 //  8192/2 the maximum npts could be handled by Rf program
#define DELTAMAX 0.1
int main(int argc, char *argv[] )
{
  if(argc!=6)
  {
	 cout<<"input 1]saclsist  2]out_too_long/short_sac_dir 3]out_no_header_sac_dir 4]out_in_comp_sac_dir 5]min_SACtime\n";
	return 0;
  }

  FILE *fin,*fout,*frm,*fcmp;
  SAC_HD sd;
  float sig[NPTSMAX],delta,maxtime,dt,tmin;
  int i,npts,factor;
  int nptsz,nptse,nptsn;
  char fnameZ[300],fnameN[300],fnameE[300],fname[300],str[300],dirlong[300],dirhead[300],dir_incomp[300];
  char str1[300],str2[300];
  if(!(fin = fopen(argv[1],"r")))
	{
	 cout<<"cannot read saclist "<<argv[1]<<endl;return 0;
	}
  fout=fopen("tmp_dec.csh","w"); //script used to decimate SAC
  frm=fopen("tmp_mv_longshort.csh","w"); // script used to mv bad SAC
  fcmp=fopen("tmp_mv_incmp.csh","w");
  sprintf(dirlong,argv[2]);
  sprintf(dirhead, argv[3]);
  sprintf(dir_incomp,argv[4]);
  tmin=atof(argv[5]);
  fprintf(fout,"sac << eof\n");
  for(i=0;;)
  {
	if(fscanf(fin,"%s",fnameZ)==EOF) break;
//	cout<<fnameZ<<"------\n";
        sprintf(fnameN,replace_str(fnameZ,"BHZ","BHN")) ;
        sprintf(fnameE,replace_str(fnameZ,"BHZ","BHE")) ;
        sprintf(fname,replace_str(fnameZ,"BHZ","BH?")) ;
	//------- check the availibility of 3 components.
	if (  FileNotExists(fnameN) or  FileNotExists(fnameE)) //normally, the fnameZ does exists, since it apperars in saclist 
	 {
	    cout<<"####file "<<fname<<" has no N or E comp\n";
 	   fprintf(fcmp,"mv %s %s/\n",fname,dir_incomp);
	   continue;
	  }
        if( !read_sac(fnameZ,sig,&sd,NPTSMAX))
	{
	    cout<<"######file "<<fnameZ<<" cannot be found!!\n";
	   fprintf(fcmp,"mv -f %s %s\n",fname,dir_incomp);
	   continue;	
	}
	//------ check the length of the SAC 
	if(sd.e < tmin)
	{
	 fprintf(stderr, "####### the input SAC %s it shorter than %f\n",fname,tmin);
	 fprintf(frm,"mv -f %s %s\n",fname,dirlong);
	 continue;
	}
	//factor = sd.npts/RFNPTS+1;
	factor = int((tmin+3)/(sd.delta*RFNPTS))+1; // check the number of points after cutting from 0 tmin
        dt = sd.delta*factor;
        if(dt>DELTAMAX+0.01)
	{fprintf(stderr,"###### the input SAC %s is toooo loooong\n",fname);
	 fprintf(frm,"mv -f %s %s\n",fname,dirlong);
	 continue;}
	// ---- check if header is written in
	if(sd.dist<0)
	{
         fprintf(stderr,"##### the input SAC %s does not have evla evlo info\n",fname);
	 fprintf(frm,"mv -f %s %s\n",fname,dirhead);
	 continue;}
	// --- check if header cmpaz cmpinc is right
	read_sac(fnameZ,sig,&sd,NPTSMAX);
	nptsz=sd.npts;
	if ( fabs(sd.cmpaz)>0.0001 or fabs(sd.cmpinc)>0.0001)
	{
	  sd.cmpaz=0.;
	  sd.cmpinc=0.;
	  write_sac(fnameZ,sig,&sd);
	}

	read_sac(fnameN,sig,&sd,NPTSMAX);
	nptsn=sd.npts;
	if ( fabs(sd.cmpaz)>0.0001 or fabs(sd.cmpinc-90.0)>0.0001)
	{
	  sd.cmpaz=0.;
	  sd.cmpinc=90.;
	  write_sac(fnameN,sig,&sd);
	}

	read_sac(fnameE,sig,&sd,NPTSMAX);
	nptse=sd.npts;
	if ( fabs(sd.cmpaz-90.)>0.0001 or fabs(sd.cmpinc-90)>0.0001)
	{
	  sd.cmpaz=90.;
	  sd.cmpinc=90.;
	  write_sac(fnameE,sig,&sd);
	}
	//---- sometimes, E and N have diff poits/length than Z
	if( nptse*sd.delta<tmin or nptsn*sd.delta<tmin )
	{
	 fprintf(stderr, "####### the input SAC %s N or E is shorter than %f\n",fname,tmin);
	 fprintf(frm,"mv -f %s %s\n",fname,dirlong);
	 continue;
	}
	//----- if points satified criteria, don't do decimate
	//if (sd.nptsz<=RFNPTS) continue;
	if(factor<2)continue;
	fprintf(fout,"r %s\ndecimate %d\nw over\n",fname,factor);	
  }//for i read saclst
  fprintf(fout,"q\neof\n");
  fclose(fin);  
  fclose(fout);
  fclose(frm);
/*  sprintf(str,"csh tmp_dec.csh\n");
  cout<< system(str)<<endl;
  sprintf(str1,"csh tmp_mv_bad.csh\n");
//   cout<<system(str1)<<endl;
  sprintf(str2,"csh tmp_mv_incmp.csh\n");
     cout<<"do!\n"<<str2<<endl;
   cout<< system(str2)<<endl;
     cout<<"end do\n";
  system("csh tmp_mv_incmp.csh");
//  cout<<system("ls")<<endl;
*/
  return 1;
}


