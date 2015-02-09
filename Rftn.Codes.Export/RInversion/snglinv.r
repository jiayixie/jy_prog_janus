#
#  SINGLINV - a single receiver function inversion program
#
#  VERSION 2.1 George Randall and Chuck Ammon 1997
#
#     This version uses the Poisson's Ratio of the Initial Model
#
#
#
      program snglinv

      parameter(NLMAX = 45, NTMAX = 520, NSMAX = 2, NDAT = NTMAX*NSMAX+2*NLMAX)
      dimension alpha(NLMAX),beta(NLMAX),rho(NLMAX),thiki(NLMAX)
      character*32 modela,title
      integer inunit,ounit,oun2
      logical porsv(NSMAX), yesno
      common /seismo/ seis(NTMAX,NSMAX), dt(NSMAX), dura(NSMAX), dly(NSMAX),gauss(NSMAX),p(NSMAX),nt(NSMAX),porsv
      common /imodel/alpha,beta,thiki,rho,nlyrs
      common /innout/ inunit,ounit
      real tfraction
      character*24 todays_date
      real fmin
      integer npasses
      logical hpfilter
      common /filter/ fmin, npasses, hpfilter
    
#
      inunit = 5;  ounit = 6;  oun2 = 8
      fmin = 0.03; npasses = 2; hpfilter = .false.
#
#**********************************************************************************************
#     Where to place blame
      write(ounit,'(/)')
      write(ounit,*) '**********************************************************'
      write(ounit,*)'snglinv - Receiver function inversion program.'
      write(ounit,*)'          VERSION 2.1 July 1997'
      write(ounit,*)'    Charles J. Ammon and George Randall.'
      write(ounit,*)'Additional routines by George Zandt and Tom Owens.'
      write(ounit,*) '**********************************************************'
#**********************************************************************************************
      call fdate(todays_date)
      write(ounit,*) 'Inversion run on: ',todays_date
      write(ounit,*) '**********************************************************'
      write(ounit,*)'Maximum Number of points in each waveform = 512'
      write(ounit,*) '**********************************************************'
#**********************************************************************************************
#
      do i = 1,NLMAX {
	   alpha(i) = 0.
	   beta(i) = 0.
	   rho(i) = 0.
	   thiki(i) = 0.
	   }
#
#      p = true
#      sv = false
#
       do i = 1,NSMAX{
	   porsv(i) = .true.}
#
      write(ounit,*)'input velocity model:'
      read(inunit,'(a)')modela
      write(ounit,*)'Enter the max number of iterations per inversion'
      read(inunit,*) maxiter
      write(ounit,*)'Enter the smoothing trade-off parameter'
      read(inunit,*) sigjmpb
      write(ounit,*)'Enter inversion ID number (for output naming)'
      read(inunit,*) invnum
      write(ounit,*)'Enter Singular Value truncation fraction'
      read(inunit,*) tfraction
      
      hpfilter = yesno('Apply a high-pass filter to waveforms? ')
      if(hpfilter) {
         write(ounit,*)'Enter the corner frequency.'
         read(inunit,*) fmin
         write(ounit,*) 'Enter the number of filter passes (1 or 2).'
         read(inunit,*) npasses
      }
#
#
#************************************************************************************************
#
# - - read in the waveform for the inversions
#
      call getseis(ns,seis,NTMAX,NSMAX,dt,dura,dly,gauss,p,nt,porsv)
#

#************************************************************************************************
#     Echo the input parameters
#
      write(ounit,*) ' '
      write(ounit,*) '**********************************************************'
      write(ounit,*) 'Inversion Input Parameters:'
      write(ounit,*) '**********************************************************'
      write(ounit,*) ' '
      write(ounit,*)'Initial Velocity Model: ', modela
      write(ounit,*)'Max Number of Iterations: ',maxiter
      write(ounit,*)'Smoothing trade-off parameter: ',sigjmpb
      write(ounit,*)'Singular-Value Truncation Fraction: ', tfraction
      write(ounit,'(a22,1x,i2.2)')' Inversion ID Number: ', invnum
      write(ounit,*) ' '
      write(ounit,*) '**********************************************************'
      write(ounit,*) 'Inversion Diagnostics:'
      write(ounit,*) '**********************************************************'
#
#************************************************************************************************
#
# - - read in the initial velocity model
#
      open(unit=oun2,file=modela)
      rewind=oun2
      read(oun2,100)nlyrs,title
100   format(i3,1x,a32)
      do i1 = 1,nlyrs {
	 read(oun2,110)idum,alpha(i1),beta(i1),rho(i1),thiki(i1),dum1,dum2,dum3,dum4,dum5
                      }
110   format(i3,1x,9f8.4)
      close(unit=oun2)
      
#
#     invert the waveform
#
      call jinv(sigjmpb,maxiter,ns,invnum,tfraction)

      stop
      end
