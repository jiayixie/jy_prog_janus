#data preparation, write sac head and remove instrument response
#de resp can be done one by one only!!

set dir_resp = /home/jiayi/work/RF_test/data_prep/resp
set f1 = `echo 170. | awk '{print 1./$0}'`
set f2 = `echo 150. | awk '{print 1./$0}'`
set f3 = `echo 5. | awk '{print 1./$0}'`
set f4 = `echo 4. | awk '{print 1./$0}'`

#==========write SAC header by using the ouput of gen_email_iris.py==============
sac << eof1
m ch.txt
eof1

#==========de response===========================================================
mkdir SAC_deresp
foreach sacz(  `ls -1 20*.*BHZ.*SAC` )
set sacn = ` echo $sacz | awk '{sub("BHZ","BHN",$1);print $1}'`
set sace = ` echo $sacz | awk '{sub("BHZ","BHE",$1);print $1}'`
set sacnmz = `ls $sacz | awk '{split($0,a,".");print a[7]"."a[8]"."a[9]"."a[10]}'`
set sacnmn = ` echo $sacnmz | awk '{sub("BHZ","BHN",$1);print $1}'`
set sacnme = ` echo $sacnmz | awk '{sub("BHZ","BHE",$1);print $1}'`
set respz = $dir_resp/RESP.$sacnmz
set respn = $dir_resp/RESP.$sacnmn
set respe = $dir_resp/RESP.$sacnme
set sac
#if !  ( (-e  $respz) && ( -e $respn) && (-e $respe) ) then
if !   (-e  $respz)  then
	echo "####RESP for $sacnmz doesn't exist!!!\n"
	continue
endif

sac << eof
r $sacz $sacn $sace
rmean
rtr
w over
r $sacz
transfer from evalresp fname $respz to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sacz
r $sacn
transfer from evalresp fname $respn to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sacn
r $sace
transfer from evalresp fname $respe to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sace
q
eof

end
#===============
mkdir SAC_raw
mv 20*.SAC SAC_raw/
