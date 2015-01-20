// based on /home/weisen/PROGS_64/EIKONAL/SCRIPT/slow_maps_to_iso_map_ani_data_v4_ANT_265_weight_robust_cv_v1_2
// reduce the weighting level for ultra high values
#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
using namespace std;

double get_dist(double lat1,double lon1,double lat2,double lon2)
{
  double theta,pi,temp;
  double radius=6371;
  pi=4.0*atan(1.0);

  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)*sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)*sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)+cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if(temp>1)
    {
      cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<endl;
      temp=1;
    }
  if(temp<-1)
    {
      cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
}

int main(int na, char *arg[])
{
  if(na!=14)
    {
      cout<<"usage:travel_time_to_velocity_map 1)file contains list of slowness dir 2)file contains list of station.lst 3)min 4)max 5)N_bin 6)out_name 7)dx 8)x0 9)nx 10)y0 11)ny 12)pflag 13)cridist"<<endl;
      return 0;
    }
  FILE *finslowlst,*finstalst,*finslow,*finsta,*file_anin,*file_iso,*file_ani;
  int i,j,k;
  int npts_x,npts_y;
  char buff1[300],sta1[10],name_iso[100],name_ani[100],name_ani_n[100],pflag[5];
  double lat,lon,t_lat,t_lon,radius,pi,sta1_lon,sta1_lat;
  int t_i,t_j,nsta;
  int ii,jj,kk,kkk,min_n,ipflag,Nslowlst,Nstalst,islow;
  double min,max,d_bin;
  int N_bin;
  
  vector<string> slowlst,stalst; // contains the list of input directory that contains slowness files, and corresponding station.lst list
  char filenmslowlst[200],filenmstalst[200],str[200],fileslowdir[200],filestalst[200];  
  char event_name[300],tstr[300];
  double dx,dy,x0,y0,x1,y1,temp,lat_temp,temp2,trash1,trash2,cridist;
  FILE *fall,*fall1;
  
  i=1;
  sprintf(filenmslowlst,arg[i++]);
  sprintf(filenmstalst,arg[i++]);
  min=atof(arg[i++]);
  max=atof(arg[i++]);
  N_bin=atoi(arg[i++]);
  sprintf(name_iso,"%s.iso",arg[i]);
  sprintf(name_ani,"%s.ani",arg[i]);
  sprintf(name_ani_n,"%s_ani_n",arg[i++]);  
  dx = atof(arg[i++]);
  dy = dx;
  x0 = atof(arg[i++]);
  npts_x = atoi(arg[i++]);
  y0 = atof(arg[i++]);
  npts_y = atoi(arg[i++]);
  ipflag=atoi(arg[i++]);
  cridist = atof(arg[i++]); // minimum distance to central station

  fprintf(stderr,"Memory check!!\n");
  double hist[N_bin];
  double slow_sum1[N_bin];
  double slow_un[N_bin];
  double slow[npts_x][npts_y][2000];
  double azi[npts_x][npts_y][2000];
  double flag[npts_x][npts_y][2000];
  double weight[2000];
  double nw[2000];
  double cdist1 = 250.;  // minimum station distance
  int idw[2000];
  double weight_sum;
  int n[npts_x][npts_y];
  double trash;
  double slow_sum[npts_x][npts_y],slow_std[npts_x][npts_y];
  double cvlon,cvlat,tdist;
  fprintf(stderr,"Memory enough!!\n");

  d_bin=(max-min)/N_bin;
  radius=6371.1391285;
  pi=4.0*atan(1.0);

  //---read in the slowness directory list and the corresponding station.lst list
  if((finslowlst=fopen(filenmslowlst,"r"))==NULL){
        printf("Cannot open points file %s!\n",filenmslowlst);exit(0);
  }
  i=0;
  while(1){
    if(fscanf(finslowlst,"%s",&tstr[0])==EOF)
    	break;
    //printf("%s,dir=%s\n",filenmslowlst,tstr);
    slowlst.push_back(tstr);
    i+=1;
  }//while 1 
  Nslowlst=i;
  
  if((finstalst=fopen(filenmstalst,"r"))==NULL){
        printf("Cannot open points file %s!\n",filenmstalst);exit(0);
  }
  j=0;
  while(1){
    if(fscanf(finstalst,"%s",&tstr[0])==EOF)
	break;
    stalst.push_back(tstr);
    j+=1;
  }
  Nstalst=j;
  printf("====read in %d slowness directories, and %d station.lst\n",Nslowlst,Nstalst);
  if(i!=j){
    printf("number of slowness_dir and number of station.lst does not match!\n");
    exit(0);
  }


  fall = fopen("all.measurement.dat","w");
  fall1 = fopen("all.measurement.used.dat","w");
  

  x1=x0+(npts_x-1)*dx;
  y1=y0+(npts_y-1)*dy;

  if (ipflag == 1) {
      sprintf(pflag,"phase.c");
      }
  else {
      sprintf(pflag,"group.c");
      }


  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  slow_sum[i][j]=0;
	  n[i][j]=0;
	}
    }


  int Nexist,Nne;
  Nexist=Nne=0;

  for (islow=0;islow<Nslowlst;islow++){
  sprintf(fileslowdir,"%s",slowlst[islow].c_str());
  sprintf(filestalst,"%s",stalst[islow].c_str());

  printf("===%d/%d, slowdir=%s; station.lst=%s\n",islow,Nslowlst,fileslowdir,filestalst);

  finsta=fopen(filestalst,"r");
  nsta=0;
  char *pch;
  char *pch1[20];
  for(;;)
    {
      if(fscanf(finsta,"%s %lf %lf",&event_name,&cvlon,&cvlat)==EOF)
        break;
      nsta++;
      sprintf(buff1,"%s/slow_azi_%s.%s.txt.HD.2.v2",fileslowdir,event_name,pflag);
      if((finslow=fopen(buff1,"r"))==NULL)
       {
	 //cout<<buff1<<" not exist!!"<<endl;
//	 return 1;
	 Nne+=1;
         continue;
       }
      Nexist+=1;
      //cout<<buff1<<endl;
      for(;;)
       {
	 //if(fscanf(finslow,"%lf %lf %lf %lf %lf",&lon,&lat,&temp,&temp2,&trash)!= 5) break;

         if (fgets(tstr,300,finslow) == NULL ) break;
         pch = strtok(tstr," \n");
         k = 0;
         while (pch!=NULL) {
           pch1[k] = pch;
           pch = strtok(NULL," \n");
           k = k + 1;
           }
         if (k == 5) continue;
         lon = atof(pch1[0]);
         lat = atof(pch1[1]);
         temp = atof(pch1[2]);
         temp2 = atof(pch1[3]);

         tdist = get_dist(cvlat,cvlon,lat,lon);
         if (tdist < cridist + 50.)
             continue;
	 //if (trash >= 2 || temp2 > 900) 
         //    continue;
	 if (trash > cdist1 || temp2 > 900)
             continue;

         ////////////// 2396 and 2398 ///
         /*
         if (fabs(lon-239.6)<0.1 && fabs(lat-44) < 0.1) {
            fprintf(f2396,"%g %g %s \n",temp,temp2,event_name);
            }
         if (fabs(lon-239.8)<0.1 && fabs(lat-44) < 0.1) {
            fprintf(f2398,"%g %g %s \n",temp,temp2,event_name);
            }
         */
         /////////////////////////////////
        

//         cout<<lon<<lat<<endl;
	 if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	   continue;
	 i=int((lon-x0)/dx+0.1);
	 j=int((lat-y0)/dy+0.1);
//         cout<<i<<" "<<lon<<" "<<j<<" "<<lat<<" "<<endl;
	 if(temp<0.6 &&temp>0.15)
	   {
             fprintf(fall,"%g %g %g %g 1 %s\n",lon,lat,1./temp,temp2,event_name);
	     slow[i][j][n[i][j]]=temp;
	     azi[i][j][n[i][j]]=temp2;
             flag[i][j][n[i][j]]=0;
	     n[i][j]++;
	     //	     slow_sum[i][j]+=temp;
	   }
       }
      //cout<<buff1<<endl;
      fclose(finslow);
    }// for ;;, read in each station.lst
  fclose(finsta);
  printf("%d slowness files exist, %d not exist\n",Nexist,Nne);
  }//for islow 
  fclose(fall);
  cout<<"ok read in slowness"<<endl;
  //-----------------------------------
  file_iso=fopen(name_iso,"w");
  nsta=50;
  double w2,tave,tstd;
  double temp_slow_sum;
  int temp_n;
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  if(n[i][j]<0.3*nsta)
	    {
	      fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	      continue;
	    }
	  w2=0;
	  weight_sum=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      weight[k]=0;
              nw[k] = 0;
              tave = 0.;
	      for(kk=0;kk<n[i][j];kk++)
		{
		  if(fabs(azi[i][j][kk]-azi[i][j][k])<20. || fabs(azi[i][j][kk]-azi[i][j][k])>(360-20.))
		    weight[k] = weight[k]++;
                    nw[k]++;
                    idw[k] = kk;
                    tave = tave + slow[i][j][k];
		}
              /*
              tave = tave/nw[k];
              tstd = 0.;
              for(kk=0;kk<nw[k];kk++) {
                   tstd = tstd + (tave - slow[i][j][idw[kk]])*(tave - slow[i][j][idw[kk]]);
                   }
              tstd = sqrt(tstd/(nw[k]-1));
              if (fabs(slow[i][j][k]-tave) > tstd && nw[k]< 30.) weight[k] = weight[k]*pow((fabs(slow[i][j][k]-tave)/tstd),2);              
              */
	      weight[k]=1/weight[k];
              //while (nw[k]<20) { weight[k] = weight[k]/2.; nw[k]=nw[k]*2; }
              //if (tstd > 0.05) weight[k] = 0.;
	      weight_sum+=weight[k];
	    }
          /// reduce the largest weight[k] to some value.
          tave = weight_sum/n[i][j];
          tstd = 0.;
          for(k=0;k<n[i][j];k++) {
            tstd = tstd + (weight[k] - tave)*(weight[k] - tave);
            }
          tstd = sqrt(tstd/n[i][j]);
          weight_sum = 0.;
          for(k=0;k<n[i][j];k++) {
            if (weight[k] > tave + 3*tstd) weight[k] = tave + 3*tstd;
            weight_sum+=weight[k];
            }
          ///
	  for(k=0;k<n[i][j];k++)
	    {
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
	    }
	  //	  if(n[i][j]>=0.5*nsta)
	  //{
	  //  
	  //  slow_sum[i][j]=slow_sum[i][j]/n[i][j];
	  temp_slow_sum=slow_sum[i][j];
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2));
	  w2=0;
	  weight_sum=0;
	  slow_sum[i][j]=0;
	  temp_n=0;

	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      weight_sum+=weight[k];
	      temp_n++;
	    }
	  //	   if(i==marker_i&&j==marker_j)
	  // fprintf(stderr,"%d %d\n",n[i][j],temp_n);
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j]) 
		continue;
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
              flag[i][j][k]=1;
              fprintf(fall1,"%g %g %g %g %g\n",x0+i*dx,y0+j*dy,1./slow[i][j][k],azi[i][j][k],weight[k]);
	    }
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2));
	  temp=slow_std[i][j]*sqrt(w2)/slow_sum[i][j]/slow_sum[i][j];
	  //	  temp=sqrt(temp/(n[i][j]-1)/n[i][j])/slow_sum[i][j]/slow_sum[i][j];
	  // cout<<x0+i*dx<<" "<<y0+j*dy<<" "<<1/slow_sum[i][j]<<" "<<temp<<" "<<n[i][j]<<endl;
	  fprintf(file_iso,"%lf %lf %lf %lf %d\n",x0+i*dx,y0+j*dy,1/slow_sum[i][j],temp,temp_n);
	  //}
	  //else
	  //{
	  //  //    cout<<x0+i*dx<<" "<<y0+j*dy<<" 0 "<<" 9999 "<<n[i][j]<<endl;
	  //  fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	  //}
	}
    }
  fclose(fall1);
  fclose(file_iso);

  //abort();

  file_ani=fopen(name_ani,"w");
  file_anin=fopen(name_ani_n,"w");  

  // 2396 , 2398
  //fclose(f2396);
  //fclose(f2398);
  //
  double test_lon,test_lat;
  //test_lon=240.6;
  //test_lat=36.6;
  //
  char t_name[300];
  double tslow;
  FILE *file_point;

  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
        {
          //sprintf(t_name,"%g_%g.raw\0",x0+i*dx,y0+j*dy);
          //file_point = fopen(t_name,"w");
	  //  i=int((test_lon-x0)/dx+0.1);
	  //j=int((test_lat-y0)/dy+0.1);
	  //cout<<i<<" "<<j<<endl;
	  //fprintf(file_ani,"\n>\n%lf %lf\n",x0+i*dx,y0+j*dy);
	  for(k=0;k<N_bin;k++)
	    {
	      hist[k]=0;
	      slow_sum1[k]=0;
	      slow_un[k]=0;
	    }
	  //if(i-3<0||i+3>=npts_x||j-3<0||j+3>=npts_y)
	  if(i-4<0||i+3>=npts_x||j-4<0||j+3>=npts_y)
	    continue;
	  kkk=0;
	  //	  min_n=999999999;
	  //for(ii=i-3;ii<=i+3;ii+=3)
	  for(ii=i-4;ii<=i+4;ii+=4)
	    {
	      //for(jj=j-3;jj<=j+3;jj+=3)
	      for(jj=j-4;jj<=j+4;jj+=4)
		{
		  if(n[ii][jj]<nsta*0.3)
		    continue;
		  kkk+=n[ii][jj];
		  //  if(n[ii][jj]<min_n)
		  //min_n=n[ii][jj];
		}
	    }
	  fprintf(file_anin,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kkk);
	  if(kkk<9*nsta*0.3||n[i][j]<nsta*0.3)
	    continue;
	  //for(ii=i-3;ii<=i+3;ii+=3)
	  //
	  //
          sprintf(t_name,"%g_%g.raw\0",x0+i*dx,y0+j*dy);
          file_point = fopen(t_name,"w");

	  //for(ii=i-4;ii<=i+4;ii+=4)
	  for(ii=i-3;ii<=i+3;ii+=3)
	    {
	      for(jj=j-3;jj<=j+3;jj+=3)
	      //for(jj=j-4;jj<=j+4;jj+=4)
		{
		  if(n[ii][jj]<nsta*0.3)
		    continue;
		  for(k=0;k<n[ii][jj];k++)
		    {
                      if (flag[ii][jj][k]<=0)
                        {
                          continue;
                        }

		      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
			{
			  fprintf(stderr,"out of range!!");
			  return 1;
			}
		      hist[int((azi[ii][jj][k]-min)/d_bin)]++;
		      slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]+=slow[ii][jj][k]-slow_sum[ii][jj];
                      tslow=slow_sum[i][j]+(slow[ii][jj][k]-slow_sum[ii][jj]);
                      fprintf(file_point,"%g %g %g %g %g\n",tslow,1./tslow,azi[ii][jj][k],slow[ii][jj][k],1./slow[ii][jj][k]);
		      //slow_un[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*(slow[ii][jj][k]-slow_sum[ii][jj]);
		    }
		}
	    }
          fclose(file_point);
	  kk=0;
	  
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  kk++;
		}
	    }
	  fprintf(file_ani,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kk);
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  slow_sum1[k]=slow_sum1[k]/hist[k];
		  //		  slow_un[k]=slow_un[k]-slow_sum1[k]*slow_sum1[k]*hist[k];
		  //slow_un[k]=sqrt(slow_un[k]/(hist[k]-1)/hist[k])/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]); //uncertainty of vel not slow
		  slow_un[k]=slow_std[i][j]/sqrt(double(hist[k]));
		  slow_un[k]=slow_un[k]/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]);//uncertainty of vel not slow
		  //slow_un[k]=slow_std[i][j];
		  //cout<<min+(0.5+k)*d_bin<<" "<<1/slow_sum1[k]<<" "<<slow_un[k]<<endl;
		  fprintf(file_ani,"%lf %lf %lf\n",min+(0.5+k)*d_bin,1/(slow_sum[i][j]+slow_sum1[k]),slow_un[k]);
		}
	    }
	  //	  fprintf(stderr,"done!!\n");
	  //return 0;
	}
    }
  fclose(file_ani);
  fclose(file_anin);
  return 0;
}
