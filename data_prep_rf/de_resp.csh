#data preparation, write sac head and remove instrument response
#de resp can be done one by one only!!

set dir_resp = /home/jiayi/work/RF_test/data_prep/resp
set f1 = `echo 120. | awk '{print 1./$0}'`
set f2 = `echo 100. | awk '{print 1./$0}'`
set f3 = `echo 0.2  | awk '{print 1./$0}'`
set f4 = `echo 0.1 | awk '{print 1./$0}'`
echo "*******corner frequency $f1 $f2 $f3 $f4"
#==========write SAC header by using the ouput of gen_email_iris.py==============
if !( -e SAC_ch.txt )then
echo "file SAC_ch.txt does not exist!"
exit
endif
sac << eof1
m SAC_ch.txt
q
eof1
echo "*******finish change header"
#==========de response===========================================================
mkdir SAC_deresp
foreach sacz(  `ls -1 20*.*BHZ.*SAC ` ) #find SAC files according to the BHZ files, if only BHN/E/1/2 exist, ignore
echo "**********begin to work on $sacz"
set sacnmz = `ls $sacz | awk '{split($0,a,".");print a[7]"."a[8]"."a[9]"."a[10]}'`
set respz = $dir_resp/RESP.$sacnmz
set sac
#-------- z n e -------------------------------
set sacn = ` echo $sacz | awk '{sub("BHZ","BHN",$1);print $1}'`
set sace = ` echo $sacz | awk '{sub("BHZ","BHE",$1);print $1}'`
if ( ( -e $sacn) && ( -e $sace )) then
set sacnmn = ` echo $sacnmz | awk '{sub("BHZ","BHN",$1);print $1}'`
set sacnme = ` echo $sacnmz | awk '{sub("BHZ","BHE",$1);print $1}'`
set respn = $dir_resp/RESP.$sacnmn
set respe = $dir_resp/RESP.$sacnme
if ( ( -e  $respz) && ( -e $respn) && (-e $respe) )then
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
echo "*******finish de resp z n e"
endif
endif
#--------z 1 2 -------------------
set sac1 = ` echo $sacz | awk '{sub("BHZ","BH1",$1);print $1}'`
set sac2 = ` echo $sacz | awk '{sub("BHZ","BH2",$1);print $1}'`
if ( ( -e $sac1) && ( -e $sac2 )) then
set sacnm1 = ` echo $sacnmz | awk '{sub("BHZ","BH1",$1);print $1}'`
set sacnm2 = ` echo $sacnmz | awk '{sub("BHZ","BH2",$1);print $1}'`
set resp1 = $dir_resp/RESP.$sacnm1
set resp2 = $dir_resp/RESP.$sacnm2
if ( ( -e  $respz) && ( -e $resp1) && (-e $resp2) )then
sac << eof
r $sacz $sac1 $sac2
rmean
rtr
w over
r $sacz
transfer from evalresp fname $respz to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sacz
r $sac1
transfer from evalresp fname $resp1 to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sac1
r $sac2
transfer from evalresp fname $resp2 to vel freqlimits $f1 $f2 $f3 $f4
w SAC_deresp/$sac2
q
eof
echo "*******finish de resp z 1 2"
endif
endif

end
#===============
mkdir SAC_raw
mv 20*.SAC SAC_raw/
#=====rotate the seismogram to r t =============
#cd SAC_deresp
#sac << eof
#m ../SAC_rot_ne.m 
#m ../SAC_rot_12.m 
#q
#eof
#mv SAC_rot ../
#==========================================



