set yalst = ( 2011 ) #2006 2007 2008 )
set monlst = (01 02 03 04 05 06 07 08 09 10 11 12 )
set daylst1  = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
set here = `pwd`
set stationlst = /data/jiayi/work/WUS/irisInfo/station.lst #/media/CHINA/jiayi/WesternChina/station.lst
set transcutDir = /media/CHINA/jiayi/WUS
unlimit
foreach iya ( 1  )
  set ya = ${yalst[$iya]}
  echo "year ----- "$ya
  set flag = `pwd | awk '{if(int('$ya'/4)*4=='$ya'){print "2";}else{print "1"}}'` #this is not a precise way to check leap year     
  foreach im ( 1 2 3 4 5 6 7 8 9 10 11 12 )
      set mon = ${monlst[$im]}
      set allday = ${daylst1[$im]}
      if ( $flag == 2 && $im == 2 ) then
        set allday = 29
      endif
      mkdir ${ya}_${mon}
      cd ${ya}_${mon}
      set dd = 1
      rm -f tempAll.ft evlst
      while ( $dd <= $allday )
        set datedir = `echo $dd | awk '{printf "%04d_%02d_%02d_0_0_0",'$ya','$mon','$dd'}'`	
	mkdir $datedir
	cd $datedir
	echo $datedir
	ls ${transcutDir}/${ya}_${mon}/$datedir/ft_*SAC | awk '{split($0,a,"'$datedir'/");print a[2]}' > temp.ft
	more temp.ft >> ../tempAll.ft
	python /home/jiayi/progs/jy/pre_cor/get_param.py temp.ft param_filter.dat param_whiten.dat
	/home/jiayi/progs/jy/pre_cor/filter4/filter4 param_filter.dat ${transcutDir}/${ya}_${mon}/$datedir .
	/home/jiayi/progs/jy/pre_cor/whiten_outphamp_LOVE/whiten_phamp param_whiten.dat
	echo $datedir >> ../evlst
	rm -f ft_*SAC
	rm -f temp.ft param_filter.dat param_whiten.dat 
	cd ..
	@ dd = $dd + 1
      end
#############################
#      python	/home/jiayi/progs/jy/do_cor/get_stalst_for_cor_LHEN.py tempAll.ft $here/errCorMsg.txt stalst
#      /home/jiayi/progs/jy/do_cor/justCOR_LHEN 3000 #stalst and evlst are default input
      echo "end of Cor LHE-----------"
      rm -f evlst tempAll.ft #stalst
   	
      cd $here
  end #im
end #iya
