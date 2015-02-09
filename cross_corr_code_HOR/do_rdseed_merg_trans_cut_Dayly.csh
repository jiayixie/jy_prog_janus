set yalst = ( 2009 2010 2011 )
set monlst  = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
set daylst1 = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
#set daylst2 = ( 31 29 31 30 31 30 31 31 30 31 30 31 )
set here = `pwd`
set seeddir = $here/seeds_00_11
set stalst =  $here/irisInfo/station.lst
set channel1 = "HE"
set channel2 = "HN"
unlimit
foreach iya ( 1 2 3 )
  set ya = ${yalst[$iya]}
  echo "year -- "$ya
  set flag = `pwd | awk '{if(int('$ya'/4)*4=='$ya'){print "2";}else{print "1"}}'`
#  echo $flag
  foreach im ( 1 2 3 4 5 6 7 8 9 10 11 12 )
    set mon = ${monlst[$im]}
    set allday = ${daylst1[$im]}
    mkdir ${ya}_${mon}
    cd ${ya}_${mon}
    if ( $flag == 2 && $im == 2 ) then
	set allday = 29
    endif
    echo "month --" $mon,$allday
    set dd = 1
    ls $seeddir/Nov9_${ya}_${mon}_*.seed | awk '{split($0,a,".");print "mv "$0,a[1]"."a[3]}' | csh
    
    while ( $dd <= $allday )
        #---rdseed------------------------------------------------------------------------------------------------------------------    
	set date = `echo $dd | awk '{printf "%04d_%02d_%02d",'$ya','$mon','$dd'}'`
	echo $date
	if ( -e  ${date}_0_0_0 )then
		echo $date "exist! skip!"
  		@ dd = $dd + 1 
		continue
	endif
	set seednm =  $seeddir/Nov9_$date.seed 
	echo $seednm
        if !( -e $seednm  )then
                echo $date >> $here/unavailableSeed.lst
		@ dd = $dd + 1
		continue
        endif
        csh $here/do_rdseed_dayly.csh $seednm

        #--- do merge, trans, cut ---------------------------------------------------------------------------------------------------
        #===rm and change name
 #       if ( 1 == 0 )then
        echo "----rm and change name-----"
        #mkdir sac_10
        rm  -f *.*.*.*.*.*.*.*.10.*.SAC 
        rm -f RESP.*.*.10.* 
        rm -f SAC_PZs_*_*_*_10_* 
        ls RESP.*H? | awk '{split($1,a,".");if(a[5]~/BH/){sub("BH","LH",a[5]);};print "\\mv "$1,a[1]"."a[2]"."a[3]"."a[5]}' | csh
        ls SAC_PZs_*H?_* | awk '{split($1,a,"_");if(a[5]~/BH/){sub("BH","LH",a[5]);};print "\\mv "$1,"SAC_PZs_"a[3]"_"a[4]"_"a[5]}' | csh
  #      endif
        #======begin======
        echo "---- do merge trans cut ----"
        echo $ya $mon $dd > event.dat  
	mkdir ${date}_0_0_0 #this could be done in the c program as well
        #----- do LHE---------------------------------
        #get a shorter station.lst based on the data in present dir. Also,seperate station into two groups(SingleSac and MultiSac). And, generate a sac macro used to deresp MultiSac
        python    /home/jiayi/progs/jy/LOVE_noise/get_stalst_and_macro.py $channel1 $stalst station_single1.lst station_multi1.lst chResp1.m . .  ${date}_0_0_0
	#=== single sac decimate-->merge-->deresp-->cut 
        /home/jiayi/progs/jy/LOVE_noise/deci_merg_and_setdb_SingleSac station_single1.lst $channel1 . . $date
        /home/jiayi/progs/jy/LOVE_noise/cut_and_trans 1000 83000 $channel1
	#=== multiple sac deresp-->decimate-->merge-->cut
	csh chResp1.m
	/home/jiayi/progs/jy/LOVE_noise/deci_merg_and_setdb_MultiSac station_multi1.lst $channel1 . . $date
	/home/jiayi/progs/jy/LOVE_noise/cut_only 1000 83000 $channel1

	#------ do LHN--------------------------------
        python    /home/jiayi/progs/jy/LOVE_noise/get_stalst_and_macro.py $channel2 $stalst station_single2.lst station_multi2.lst chResp2.m . .  ${date}_0_0_0
	#=== single sac decimate-->merge-->deresp-->cut 
        /home/jiayi/progs/jy/LOVE_noise/deci_merg_and_setdb_SingleSac station_single2.lst $channel2 . . $date
        /home/jiayi/progs/jy/LOVE_noise/cut_and_trans 1000 83000 $channel2
	#=== multiple sac deresp-->decimate-->merge-->cut
	csh chResp2.m
	/home/jiayi/progs/jy/LOVE_noise/deci_merg_and_setdb_MultiSac station_multi2.lst $channel2 . . $date
	/home/jiayi/progs/jy/LOVE_noise/cut_only 1000 83000 $channel2

	#-------------rm trash-----------------------
        rm -f rdseed.err_log*  decimate_one_LH?.csh sac_bp_rmInstr_LH?.csh sac_list_LH? temp_LH?_list station_single?.lst station_multi?.lst chResp?.m  sac_db_LH?_*.out
	rm -f 20*.M.SAC SAC_PZs_* RESP*
        #rm -f 20*.?${channel}.M.SAC SAC_PZs_*.L${channel} RESP*.L${channel} 
        #rm -f -r sac_10
        @ dd = $dd + 1
    end
    cd $here
  end
  set tm = `date`
  echo  $ya $mon : $tm >> $here/Data_finished.txt
end
date
