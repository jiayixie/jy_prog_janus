      subroutine sacin(infil,x,npts,b,delta,user2)
      dimension x(1)
      character*8 kdummy
      character*16 kdummy2
      character*32 infil
      real adummy,delta,depmin,depmax,b,e,user0,user1,user2
      integer idummy,n,npts,inunit
      integer inter1,itern2,itern3
      integer iftype,leven,lpsol,lovrok,lcalda
c       common /sheaderf/delta,depmin,depmax,b,e,user0,user1,user2
c	common/sheaderi/NPTS
c	include 'Sheader.inc'
      inunit = 8
c	open the SAC file for writing
      open(unit=inunit,file=infil)
      rewind(inunit)
c	read the SAC header
      read(inunit,1000) delta,depmin,depmax,adummy,adummy
      read(inunit,1000) b, e,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) user0, user1, user2,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1000) adummy,adummy,adummy,adummy,adummy
      read(inunit,1010) idummy,idummy,idummy,idummy,idummy
      read(inunit,1010) idummy,inter1,itern2,itern3, npts
      read(inunit,1010) idummy,idummy,idummy,idummy,idummy
      read(inunit,1010) iftype,idummy,idummy,idummy,idummy
      read(inunit,1010) idummy,idummy,idummy,idummy,idummy
      read(inunit,1010) idummy,idummy,idummy,idummy,idummy
      read(inunit,1010) idummy,idummy,idummy,idummy,idummy
      read(inunit,1010) leven, lpsol,lovrok,lcalda,idummy
      read(inunit,1015) kdummy,kdummy2
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
      read(inunit,1020) kdummy,kdummy,kdummy
c	read the data
      read(inunit,1000) (x(n), n = 1,npts)
      close(inunit)
1000  format(5g15.7)
1010  format(5i10)
1015  format(a8,a16)
1020  format(3a8)
      return
      end
