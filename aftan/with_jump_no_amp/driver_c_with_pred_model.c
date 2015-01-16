/*
 * The sample of test driver for FTAN with phase match filter for
 * subroutines aftanpg and aftanipg
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aftan.h"


int main (int argc, char *argv[])
{
  static int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
  static double t0, dt, delta, vmin, vmax, tmin, tmax;
  static double snr, tresh, ffact, perc, taperl,fmatch,piover4;
  static float sei[32768];
  static double arr1[100][8],arr2[100][7];
  static double tamp, ampo[32][32768], pred[2][300];
  static int nrow, ncol, npred;
  static double prpvper[300],prpvvel[300]; // phase vel prediction files

  char  *p,name[160],name1[160],buf[200],str[160],phvelname[160],root[160];
  char sta1[6], sta2[6];
  FILE  *in, *fd, *inv;
  int   i, j;
  int   nn,sac = 1; // =1 - SAC, =0 - ftat files

// input command line arguments treatment
  if(argc != 2) {
      printf("Usage: aftan4_c_test parameter_file\n");
      exit(-1);
  }
// open and read contents of parameter file
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch,
             name)) != EOF) { // start main loop
  strcpy(root,name);
  p = strrchr(root,'.');
  *(p+1) = '\0';
  //strcpy(phvelname,root);
  //strcat(phvelname,"SAC_PHP");

      if(n == 0 || n != 11) break;

      printf("vmin= %lf, vmax= %lf, tmin= %lf, tmax= %lf\n",
              vmin,vmax,tmin,tmax);
// remove quotes from file names
      j = 0;
      for(i = 0; i < strlen(name); i++) {
          if(name[i] == '\'' || name[i] == '\"') continue;
          name[j] = name[i]; j++;
      }
      name[j] = '\0';
      printf("Tresh= %lf, Filter factor= %lf, SNR= %lf, Match = %lf\nData file name=%s\n",
             tresh,ffact,snr,fmatch,name);
// if presents, read phase velocity prediction file
// ---
nprpv = -1;
                                                                                
/*
 *  read phase velocity information
 */
 //   sscanf(root, "COR_%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].", sta1, sta2);
	  sscanf(root,"%[0-9].%[A-Z,a-z].",sta1,sta2);
	  printf("======evnm:%s  stanm:%s root:%s\n",sta1,sta2,root);
//    sprintf(phvelname, "/home/linf/California/scalifornia_avg_phvel.dat");
//sprintf(phvelname, "/Users/jiayixie/progs/prog_CV/FTA/aftani_c_pgl/FTA_with_pred_group/scalifornia_avg_phvel.dat");  
	  sprintf(phvelname,"/home/jiayi/work/SC/data/PREDICTION/evt_sta_pair/R/%s.%s.dat",sta1,sta2);
//	  sprintf(phvelname,"/Users/jiayixie/Tianshan/Info/PRED_MODEL/PRED_R/%s.%s.dat",sta1,sta2);
//	  sprintf(phvelname,"/Users/jiayixie/Tianshan/Info/PRED_MODEL/DISP_v1_glob_R/%s.%s.dat",sta1,sta2);
	  //	  fprintf(stderr, "predicted phase velocity %s \n",phvelname);

  if((inv = fopen(phvelname,"r")) == NULL) {
      printf("Can not find file %s.\n",phvelname);
      nprpv = 0;
  } else {
  while(fgets(buf,200,inv) != NULL) {
         if(nprpv == -1) { nprpv++; continue; }
         if((n = sscanf(buf,"%lf %lf",&prpvper[nprpv],&prpvvel[nprpv])) < 2) break;
         nprpv++;
	}
         fclose(inv);
         printf("Phase velocity prediction file name= %s\n",phvelname);
  }

/*
 *   read SAC or ascii data 
 */
      readdata(sac,name,&n,&dt,&delta,&t0,sei);
/*
 * Read group velocity prediction file
 */
//      strcpy(name1,root);
//      strcat(name1,"SAC_GRP");
     sscanf(root, "COR_%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].", sta1, sta2);
//	  sscanf(root,"%[0-9].%[A-Z,a-z].",sta1,sta2);
	  //sprintf(name1, "/home/linf/California/PRED_DISP/COR_%s_%s.SAC_PRED", sta1, sta2);
      //printf("Group velocity prediction curve: %s\n",name1);
      //if((fd = fopen(name1,"r")) == NULL) {
      //    printf("Can not find file %s.\n",name1);
      //    exit(1);
      //}
      //i = 0;
      //fgets(str,100,fd);
      //while((nn = fscanf(fd,"%lf %lf",&pred[0][i],&pred[1][i])) == 2) i++;
      //npred = i;
      //fclose(fd);
/* ====================================================================
 * Parameters for aftanipg function:
 * Input parameters:
 * piover4 - phase shift = pi/4*piover4, for cross-correlation
 *           piover4 should be -1.0 !!!!     (double)
 * n       - number of input samples, (int)
 * sei     - input array length of n, (float)
 * t0      - time shift of SAC file in seconds, (double)
 * dt      - sampling rate in seconds, (double)
 * delta   - distance, km (double)
 * vmin    - minimal group velocity, km/s (double)
 * vmax    - maximal value of the group velocity, km/s (double)
 * tmin    - minimal period, s (double)
 * tmax    - maximal period, s (double)
 * tresh   - treshold, usually = 10, (double)
 * ffact   - factor to automatic filter parameter, (double)
 * perc    - minimal length of of output segment vs freq. range, % (double)
 * npoints - max number points in jump, (int)
 * taperl  - factor for the left end seismogram tapering,
 *           taper = taperl*tmax,    (double)
 * nfin    - starting number of frequencies, nfin <= 32, (int)
 * snr     - phase match filter parameter, spectra ratio to
 *           determine cutting point    (double)
 * fmatch  - factor to length of phase matching window
 * npred   - length of prediction table
 * pred    - prediction table: pred[0][] - periods in sec,
 *                             pred[1][] - pedicted velocity, km/s
 * ==========================================================
 * Output parameters are placed in 2-D arrays arr1 and arr2,
 * arr1 contains preliminary results and arr2 - final.
 * ==========================================================
 * nfout1 - output number of frequencies for arr1, (int)
 * arr1   - the first nfout1 raws contain preliminary data,
 *          (double arr1[n][5], n >= nfout1)
 *          arr1[:,0] -  central periods, s (double)
 *          arr1[:,1] -  apparent periods, s (double)
 *          arr1[:,2] -  group velocities, km/s (double)
 *          arr1[:,3] -  phase velocities, km/s (double)
 *          arr1[:,4] -  amplitudes, Db (double)
 *          arr1[:,5] -  discrimination function, (double)
 *          arr1[:,6] -  signal/noise ratio, Db (double)
 *          arr1[:,7] -  maximum half width, s (double)
 * nfout2 - output number of frequencies for arr2, (int)
 *          If nfout2 == 0, no final result.
 * arr2   - the first nfout2 raws contains final data,
 *          (double arr2[n][5], n >= nfout2)
 *          arr2[:,0] -  central periods, s (double)
 *          arr2[:,1] -  apparent periods, s (double)
 *          arr2[:,2] -  group velocities, km/s (double)
 *          arr2[:,3] -  amplitudes, Db (double)
 *          arr2[:,4] -  signal/noise ratio, Db (double)
 *          arr2[:,5] -  maximum half width, s (double)
 *          tamp      -  time to the beginning of ampo table, s (double)
 *          nrow      -  number of rows in array ampo, (int)
 *          ncol      -  number of columns in array ampo, (int)
 *          ampo      -  Ftan amplitude array, Db, (double [32][32768])
 * ierr   - completion status, =0 - O.K.,           (int)
 *                             =1 - some problems occures
 *                             =2 - no final results
 */

//  t0      = 0.0;
  nfin    = 32;
  npoints = 10;        // only 3 points in jump
  perc    = 50.0;     // 50 % for output segment
//  taperl  = 2.0;      // factor to the left end tapering
  printf("pi/4 = %5.1lf, t0 = %9.3lf\n",piover4,t0);
  printf("#filters= %d, Perc= %6.2f %s, npoints= %d, Taper factor= %6.2f\n",
          nfin,perc,"%",npoints,taperl);
/* Call aftanipg function, FTAN + prediction         */

  // printf("FTAN + prediction curve\n");

  //ffact =2.0;
  //aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
  //      &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
  //      &nprpv,prpvper,prpvvel,
  //      &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  //printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_P");
  //if(nfout2 == 0) continue;   // break aftan sequence
  //printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* FTAN with phase with phase match filter. First Iteration. */

  printf("FTAN - the first ineration\n");
  ffact =1.0;
  aftanpg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_1");
  if(nfout2 == 0) continue;   // break aftan sequence
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Make prediction based on the first iteration               */

  npred = nfout2;
  tmin = arr2[0][1];
  tmax = arr2[nfout2-1][1];
  for(i = 0; i < nfout2; i++) {
      pred[0][i] = arr2[i][1];   // apparent period // central periods
      pred[1][i] = arr2[i][2];   // group velocities
  }

/* FTAN with phase with phase match filter. Second Iteration. */

  printf("FTAN - the second iteration (phase match filter)\n");
  ffact = 2.0;
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
        &nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_2");
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

  }
  fclose(in);
  return 0;
}
