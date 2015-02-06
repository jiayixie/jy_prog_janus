c Calculates a spread of seismograms for a given model.
c Requires readmodel.f,...
c Usage: seis-spread modelfile geometryfile phasefile arrivalfile

c####&

      program seis_spread

c Anyone who fails to start a Fortran program with this line
c should be severely beaten:
        implicit none
    
c Include constants:
        include 'params.h'

c Scratch variables:
        integer i,j,iargc
        character modname*(namelen),geomname*(namelen),phname*(namelen)
        character arrname*(namelen),tracename*(namelen)
        character iphname*(namelen),handedness*(namelen)
        character tensortop*(namelen),tensormid*(namelen)
        character tensorbottom*(namelen)
        real amp_in
c Model parameters:
        integer nlay
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
c Geometry parameters
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer ntr
c Phase parameters
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph,iphase
c Arrivals
        real travel_time(maxph,maxtr),amplitude(3,maxph,maxtr)
c Traces
        real Tr_cart(3,maxsamp,maxtr),dt,width,shift
        real Tr_ph(3,maxsamp,maxtr)
        integer align,nsamp,mults,out_rot
        
        
c   aa is a list of rank-4 tensors (a_ijkl = c_ijkl/rho)
c   rot is a list of rotator matrices, used to rotate into the local
c   coordinate system of each interface.
        real aa(3,3,3,3,maxlay), rot(3,3,maxlay)
c   ar_list is aa, pre-rotated into reference frames for the interfaces
c   above (last index=2) and below (1) each respective layer.
        real ar_list(3,3,3,3,maxlay,2)
        
        
c Get filenames from command-line argument.
        i = iargc()
        if (i .lt. 10) then
          write(*,*) 'Usage: seis-tensor modelfile geometryfile ',
     &               'phasefile arrivalfile tracefile iphase',
     &               'tensortop tensormid tensorbottom handedness'
          write(*,*) '-- Model, geometry, and tensor files must exist.'
          write(*,*) '-- Other files will be overwritten.'
          write(*,*) '-- iphase may be P, SV, or SH.'
          write(*,*) '-- handedness L or R for L/R handed system'
          write(*,*) '-- righthanded hex tensor in Voigt/Love notation:'
          write(*,*) '-- A    A-2N F   0   0   0'
          write(*,*) '-- A-2N A    F   0   0   0'
          write(*,*) '-- F    F    C   0   0   0'
          write(*,*) '-- 0    0    0   L   0   0'
          write(*,*) '-- 0    0    0   0   L   0'
          write(*,*) '-- 0    0    0   0   0   N'
          write(*,*) '-- lefthanded hex tensor in Voigt/Love notation:'
          write(*,*) '-- C   F    F    0   0   0'
          write(*,*) '-- F   A    A-2N 0   0   0'
          write(*,*) '-- F   A-2N A    0   0   0'
          write(*,*) '-- 0   0    0    N   0   0'
          write(*,*) '-- 0   0    0    0   L   0'
          write(*,*) '-- 0   0    0    0   0   L'
          write(*,*) 'tensorfiles has rho in kg/m^3 in first line,'
          write(*,*) 'followed by C11,C12,..,C16 in next line,'
          write(*,*) 'C21,C22,...,C26 in next line etc., with listed'
          write(*,*) 'Cij assumed to be GPa not normalized by rho'
          stop
        end if
        call getarg(1,modname)
        write(*,*) 'Model is ',modname

        call getarg(6,iphname)
        if (iphname .eq. 'P') then
          iphase=1
         else if (iphname .eq. 'SV') then
           iphase=2
         else if (iphname .eq. 'SH') then
           iphase=3
         else
           write (*,*) iphname,' is not a valid phase type.'
           write (*,*) 'Valid phases are P, SV and SH.'
           stop
         end if
        
        write (*,*) 'Initial phase is ',iphname

        call getarg(7,tensortop)
        call getarg(8,tensormid)
        call getarg(9,tensorbottom)
        write(*,*) 'top layer tensor ',tensortop
        write(*,*) 'middle layer tensor ',tensormid
        write(*,*) 'bottom layer tensor ',tensorbottom

        call getarg(10,handedness)


c Read in model
        call readmodel(modname,thick,rho,alpha,beta,isoflag,
     &                 pct_a,pct_b,trend,plunge,strike,dip,nlay)
c print the read in model to screen
        call writemodel(6,thick,rho,alpha,beta,isoflag,
     &                  pct_a,pct_b,trend,plunge,strike,dip,nlay)
          
c Set up model for calculation, make rotators
        write(*,*) 'calling buildmodel'
        call buildmodel(aa,ar_list,rot,thick,rho,alpha,beta,isoflag,
     &                  pct_a,pct_b,trend,plunge,strike,dip,nlay,
     &                  tensortop,tensormid,tensorbottom,handedness)
     
c Read in geometry (desired traces)
        call getarg(2,geomname)
        write (*,*) 'Geometry is ',geomname
        call readgeom(geomname,baz,slow,sta_dx,sta_dy,ntr)
c        call writegeom(6,baz,slow,sta_dx,sta_dy,ntr)
        
c Read in parameters from file 'raysum-params', if it exists
        geomname='raysum-params'
        call readparams(geomname,mults,nsamp,dt,width,align,
     &                  shift,out_rot)
        
c Generate phase list
        call getarg(3,phname)
        numph=0
        if (mults .ne. 3) then
          call ph_direct(phaselist,nseg,numph,nlay,iphase)
        end if
        if (mults .eq. 1) then
          call ph_fsmults(phaselist,nseg,numph,nlay,1,iphase)
        end if
        if (mults .eq. 2) then
          do j=1,nlay-1
            call ph_fsmults(phaselist,nseg,numph,nlay,j,iphase)
          end do
        end if
        if (mults .eq. 3) then
          write(*,*) 'Reading phases from file ',phname
          call readphases(phname,phaselist,nseg,numph)
        end if
c        call printphases(phaselist,nseg,numph)
        if (mults .ne. 3) then
          open(unit=iounit1,file=phname,status='unknown')
          call writephases(iounit1,phaselist,nseg,numph)
          close(unit=iounit1)
          write(*,*) 'Phases written to ',phname
        end if

        call getarg(4,arrname)
        write(*,*) 'Arrivals will be written to ',arrname
        open(unit=iounit1,file=arrname,status='unknown')
        
        call getarg(5,tracename)
        write(*,*) 'Traces will be written to ',tracename
        open(unit=iounit2,file=tracename,status='unknown')
        
c        Perform calculation                   
        amp_in=1.
        write(*,*) 'getting arrivals'
        call get_arrivals(travel_time,amplitude,thick,rho,isoflag,
     &       strike,dip,aa,ar_list,rot,baz,slow,sta_dx,sta_dy,
     &       phaselist,ntr,nseg,numph,nlay,amp_in)
     
c        Normalize arrivals
        if (iphase .eq. 1) then
          WRITE(*,*) 'normalizing arrivals'
          call norm_arrivals(amplitude,baz,slow,alpha(1),beta(1),
     &                       rho(1),ntr,numph,1,1)
        end if
         
c        Write out arrivals
        call writearrivals(iounit1,travel_time,amplitude,ntr,numph)
        close(unit=iounit1)
        
c        Assemble traces
        write(*,*) 'assembling traces'
        call make_traces(travel_time,amplitude,ntr,numph,nsamp,
     &                   dt,width,align,shift,Tr_cart)
     
        if (out_rot .eq. 0) then
          call writetraces(iounit2,Tr_cart,ntr,nsamp,dt,align,shift)
        else
          if (out_rot .eq. 1) then
c            Rotate to RTZ
            call rot_traces(Tr_cart,baz,ntr,nsamp,Tr_ph)
          else
c            Rotate to wavevector coordinates
            call fs_traces(Tr_cart,baz,slow,alpha(1),beta(1),
     &                     rho(1),ntr,nsamp,Tr_ph)
          end if
c          Write results
          write(*,*) 'writing results'
          call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
        end if
        close(unit=iounit2)
                
      end
      
      
      subroutine readparams(filename,mults,nsamp,dt,width,align,shift,
     &                      out_rot)
      
        implicit none
        include 'params.h'
        
        character filename*(namelen),buffer*(buffsize)
        integer mults,nsamp,align,ios,eof,out_rot
        real dt,width,shift
        
c          Default values
        mults=1
        nsamp=600
        dt=0.05
        width=1.
        align=1
        shift=5.
        out_rot=2
        
        open(unit=iounit1,status='old',file=filename,iostat=ios)
        
        if (ios .eq. 0) then
          write(*,*) 'Reading parameters from ',filename
          call getline(iounit1,buffer,eof)
          read (buffer,*) mults
          call getline(iounit1,buffer,eof)
          read (buffer,*) nsamp
          call getline(iounit1,buffer,eof)
          read (buffer,*) dt
          call getline(iounit1,buffer,eof)
          read (buffer,*) width
          call getline(iounit1,buffer,eof)
          read (buffer,*) align
          call getline(iounit1,buffer,eof)
          read (buffer,*) shift
          call getline(iounit1,buffer,eof)
          read (buffer,*) out_rot
        else
          write (*,*) 'Parameter file ',filename,' does not exist.'
          write (*,*) 'Writing default values.'
          close(unit=iounit1)
          open(unit=iounit1,status='unknown',file=filename)
          write(iounit1,*) '# Multiples: 0 for none, 1 for Moho, 2 ',
     &                      'for all first-order, 3 to read file'
          write(iounit1,*) mults
          write(iounit1,*) '# Number of samples per trace'
          write(iounit1,*) nsamp
          write(iounit1,*) '# Sample rate (seconds)'
          write(iounit1,*) dt
          write(iounit1,*) '# Gaussian pulse width (seconds)'
          write(iounit1,*) width
          write(iounit1,*) '# Alignment: 0 is none, 1 aligns on P'
          write(iounit1,*) align
          write(iounit1,*) '# Shift of traces -- t=0 at this time (sec)'
          write(iounit1,*) shift
          write(iounit1,*) '# Rotation to output: 0 is NS/EW/Z, ',
     &                     '1 is R/T/Z, 2 is P/SV/SH'
          write(iounit1,*) out_rot
        end if
        close(unit=iounit1)
                  
      end







