subroutine sacout(ofil,x,NPTS,b,delta,user2)
{
	dimension x(1)

	character*8 kdummy
	character*8 kdummy2
	character*32 ofil
	real adummy
	real delta,depmin,depmen,depmax,b,e,user0,user1,user2
	integer idummy,n,NPTS,ounit
	integer itern1,itern2,itern3
	integer leven,lpsol,lovrok,lcalda

#       common /sheaderf/delta,depmin,depmax,b,e,user0,user1,user2
#	common/sheaderi/NPTS
#	include 'Sheader.inc'

	ounit = 8

#	set up the dummy values for unused header fields

	idummy = -12345
	adummy = -12345.
	kdummy = '-12345'
	kdummy2 = '-12345'

#	set upt some required SAC parameters describing the file

	iftype=1; leven=1; lpsol=1; lovrok=1; lcalda=1
	itern1=6; itern2=0; itern3=0



#	compute the min, max, and duration of the file

	call minmax(x,NPTS,depmin,depmax,depmen)

	e = b + float(NPTS - 1) * delta


#	open the SAC file for writing

	open(unit=ounit,file=ofil)
	rewind(ounit)

#	write the SAC header

	write(ounit,1000)  delta,depmin,depmax,adummy,adummy
	write(ounit,1000)      b,     e,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy, user2,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,depmen,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy
	write(ounit,1000) adummy,adummy,adummy,adummy,adummy

	write(ounit,1010) idummy,idummy,idummy,idummy,idummy
	write(ounit,1010) idummy,itern1,itern2,itern3,NPTS
	write(ounit,1010) idummy,idummy,idummy,idummy,idummy
	write(ounit,1010) iftype,idummy,idummy,idummy,idummy
	write(ounit,1010) idummy,idummy,idummy,idummy,idummy
	write(ounit,1010) idummy,idummy,idummy,idummy,idummy
	write(ounit,1010) idummy,idummy,idummy,idummy,idummy
	write(ounit,1010)  leven, lpsol,lovrok,lcalda,idummy

	write(ounit,1015) kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy
	write(ounit,1020) kdummy,kdummy,kdummy

#	write the data
	
		write(ounit,1000) (x(n), n = 1,NPTS)

	close(ounit)

1000	format(5g15.7)
1010	format(5i10)
1015	format(a8,a8)
1020	format(3a8)

	return
	end
}
      
      
      subroutine minmax(x,npts,min,max,mean)
{
	dimension x(1)
      real min,max,mean
      min=9.0e+19
      max=-9.0e+19
      mean=0.
      do i=1,npts {
           if(x(i).gt.max) max=x(i)
           if(x(i).lt.min) min=x(i)
           mean=mean + x(i)
      }
      mean=mean/float(npts)
      return
      end
}
