c ==========================================================
c Function filter4. Broadband filreting.
c ==========================================================
c Parameters for filter4 function:
c Input parameters:
c f1,f2   - low corner frequences, f2 > f1, Hz, (double)
c f3,f4   - high corner frequences, f4 > f3, Hz, (double)
c npow    - power of cosine tapering,  (int)
c dt      - sampling rate of seismogram in seconds, (double)
c n       - number of input samples, (int)
c seis_in - input array length of n, (float)
c Output parameters:
c seis_out - output array length of n, (float)
c ==========================================================

      subroutine filter4(f1,f2,f3,f4,npow,dt,n,seis_in_E,seis_in_N,
     1     seis_out_E,seis_out_N,seis_outamp_E,seis_outamp_N,seis_outph_E,seis_outph_N,ns,dom)
      implicit none
      include 'fftw3.h'
      integer*4 npow,n
      real*8    f1,f2,f3,f4,dt
      real*4    seis_in_E(400000),seis_out_E(400000)
      real*4   seis_outamp_E(400000), seis_outph_E(400000)
      real*4    seis_in_N(400000),seis_out_N(400000)
      real*4   seis_outamp_N(400000), seis_outph_N(400000)
c ---
      integer*4 k,ns,nk,i
      real*8    plan1_E,plan2_E
      real*8    plan1_N,plan2_N
      real*8    dom
      double complex czero,s_E(400000),sf_E(400000)
      double complex s_N(400000),sf_N(400000)
c ---
      czero = (0.0d0,0.0d0)



c determin the power of FFT
      ns = 2**max0(int(dlog(dble(n))/dlog(2.0d0))+1,13)
      dom = 1.0d0/dt/ns

      do k = 1,ns
        s_E(k) = czero
        S_N(k) = czero
      enddo

      do k = 1,n
        s_E(k) = seis_in_E(k)
        s_N(k) = seis_in_N(k)
      enddo

c make backward FFT for seismogram: s ==> sf
      call dfftw_plan_dft_1d(plan1_E,ns,s_E,sf_E,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1_E)
      call dfftw_destroy_plan(plan1_E)
      call dfftw_plan_dft_1d(plan1_N,ns,s_N,sf_N,
     *     FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1_N)
      call dfftw_destroy_plan(plan1_N)
c kill half spectra and correct ends

      nk = ns/2+1
      do k = nk+1,ns
        sf_E(k) = czero
        sf_N(k) = czero
      enddo
 
      sf_E(1) = sf_E(1)/2.0d0
      sf_N(1) = sf_N(1)/2.0d0
      sf_E(nk) = dcmplx(dreal(sf_E(n)),0.0d0)
      sf_N(nk) = dcmplx(dreal(sf_N(n)),0.0d0)
c  
       do k = 1,nk
          seis_outamp_E(k) = 0.0
          seis_outph_E(k)  = 0.0
          seis_outamp_N(k) = 0.0
          seis_outph_N(k)  = 0.0
       enddo


c=============================================================
c     do smoothing on sf equivalent to do " smooth mean h 20" in SAC

      call smooth(f1,f2,f3,f4,dom,nk,sf_E,sf_N,20)
C=============================================================
c===============================================================
c   make tapering
      call flt4(f1,f2,f3,f4,dom,nk,npow,sf_E,sf_N)

       do i = 1,nk
        seis_outamp_E(i)= real(dsqrt(dreal(sf_E(i))**2 +
     1                        dimag(sf_E(i))**2))
        seis_outph_E(i) = real(datan2(dimag(sf_E(i)),dreal(sf_E(i))))
        seis_outamp_N(i)= real(dsqrt(dreal(sf_N(i))**2 +
     1                        dimag(sf_N(i))**2))
        seis_outph_N(i) = real(datan2(dimag(sf_N(i)),dreal(sf_N(i))))

      enddo


c make forward FFT for seismogram: sf ==> s
      call dfftw_plan_dft_1d(plan2_E,ns,sf_E,s_E,
     *                         FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan2_E)
      call dfftw_destroy_plan(plan2_E)
      call dfftw_plan_dft_1d(plan2_N,ns,sf_N,s_N,
     *     FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan2_N)
      call dfftw_destroy_plan(plan2_N)

c forming final result

      do k = 1,n
        seis_out_E(k) = 2.0*real(dreal(s_E(k)))/ns
        seis_out_N(k) = 2.0*real(dreal(s_N(k)))/ns
      enddo



      return
      end



c===============================================================
c Tapering subroutine itself
c===============================================================
      subroutine flt4(f1,f2,f3,f4,dom,nk,npow,sf_E,sf_N)
      real*8    f1,f2,f3,f4,dom
      integer*4 nk,npow
      double complex sf_E(400000),sf_N(400000)
      real*8    d1,d2,f,dpi,ss,s(400000)
      integer*4 i,j
c ---
      dpi = datan(1.0d0)*4.0d0
      do i = 1,nk
         s(i) = 0.0d0
      enddo
      do i = 1,nk
        f = (i-1)*dom
        if(f.le.f1) then
          goto 1
        else if(f.le.f2) then
          d1 = dpi/(f2-f1)
          ss = 1.0d0
          do j = 1,npow
            ss = ss*(1-dcos(d1*(f1-f)))/2.0d0
          enddo
          s(i) = ss
       else if(f.le.f3) then
           s(i) = 1.0d0
        else if(f.le.f4) then
          d2 = dpi/(f4-f3)
          ss = 1.0d0
          do j = 1,npow
            ss = ss*(1+dcos(d2*(f3-f)))/2.0d0
          enddo
          s(i) = ss
        endif
  1     continue
      enddo
      do i = 1,nk
        sf_E(i) = sf_E(i)*s(i)
        sf_N(i) = sf_N(i)*s(i)
      enddo
      return
      end


c===============================================================
c  smoothing routine      call smooth(f1,f2,f3,f4,dom,nk,sf_E,sf_N,20)
c=s==============================================================
      subroutine smooth(f1,f2,f3,f4,dom,nk,sf_E,sf_N,number)
      real*8    f1,f2,f3,f4
      integer*4 number,nk
      double complex sf_E(400000), sf_N(400000)
      real*8    sorig(400000), sout(400000),dom
      real*8   f,sum, avg
c ---
        do i = 1,nk
         sorig(i) = dsqrt(dreal(sf_E(i))**2+dimag(sf_E(i))**2)
        enddo
     
        do i = 1,nk

        f = (i-1)*dom

        if( f .ge. f1 .and. f .le. f4 ) then
            sum = 0. 
          do jk = -number,number
             ijk = i+jk
             sum = sum + sorig(ijk)
          enddo
            sout(i) = sum/(2.*number+1.)
        else
            sout(i) = sorig(i)
        endif

       enddo


       do i = 1,nk
         f = (i-1)*dom
       if( f .ge. f1 .and. f .le. f4 ) then
          sout(i) = 1.0d0/sout(i)
       else
          sout(i) = 0.0d0
       endif
          
       enddo



        do i = 1,nk
           sf_E(i) = sf_E(i)*sout(i)
           sf_N(i) = sf_N(i)*sout(i)
        enddo

       return
 
       end


