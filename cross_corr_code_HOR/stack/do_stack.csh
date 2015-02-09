#!/csh
set stalst1 = /data/jiayi/work/WUS/irisInfo/station_jy_sort.lst  #   CI     ARV  -118.83  35.1269
#set stalst1 = ./t
#awk '{printf "%s.%s %f %f\n",$1,$2,$3,$4}' $stalst1  > temp_station.lst
awk '{printf "%s.%s %f %f\n",$1,$2,$3,$4}' $stalst1 | sort -d > temp_station.lst
#awk '{if($0~/R39A/||$0~/Z43A/)printf "%s.%s %f %f\n",$1,$2,$3,$4}' $stalst1 > temp_station.lst
set stalst = temp_station.lst
set stackdir = /data/jiayi/work/WUS/stack_green
mkdir $stackdir

set evlst = stack_info.lst
rm -f $evlst
#----

#--
set cordir = /media/disk-2/WUS/WUS_fil_whi_cor #/data/jiayi/work/WUS/do_cross/WUS_fil_whi_cor
echo "$cordir/2012_01/COR  1 " >> $evlst
echo "$cordir/2012_02/COR  1 " >> $evlst
echo "$cordir/2012_03/COR  1 " >> $evlst
echo "$cordir/2012_04/COR  1 " >> $evlst
echo "$cordir/2012_05/COR  1 " >> $evlst
echo "$cordir/2012_06/COR  1 " >> $evlst
echo "$cordir/2012_07/COR  1 " >> $evlst
echo "$cordir/2012_08/COR  1 " >> $evlst

#--
set cordir = /media/JAPAN/jiayi/WUS/WUS_fil_whi
echo "$cordir/2010_04/COR  1 " >> $evlst
echo "$cordir/2011_01/COR  1 " >> $evlst
echo "$cordir/2011_05/COR  1 " >> $evlst
echo "$cordir/2011_08/COR  1 " >> $evlst
echo "$cordir/2011_09/COR  1 " >> $evlst
echo "$cordir/2011_10/COR  1 " >> $evlst

#echo "$cordir/2011_08/COR  1 " >> $evlst
#echo "$cordir/2011_10/COR  1 " >> $evlst
#---
#set cordir = /media/disk/WUS/WUS_fil_whi_cor
#echo "$cordir/2012_01/COR  1 " >> $evlst
#echo "$cordir/2012_05/COR  1 " >> $evlst
#---
#set cordir = /media/JAPAN/jiayi/WUS/WUS_fil_whi
#echo "$cordir/2011_01/COR  1 " >> $evlst
#echo "$cordir/2011_05/COR  1 " >> $evlst
#echo "$cordir/2011_09/COR  1 " >> $evlst

/home/jiayi/progs/jy/stack/STACK_weight_count temp_station.lst $evlst $stackdir 

\cp $evlst $stackdir/

