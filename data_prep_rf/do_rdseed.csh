set title = RF_Tibet_Oct27
set pzfile = PZs
set respfile = RESP
mkdir $pzfile
mkdir $respfile
rm -f sta_resp.lst sta_pz.lst

foreach seednm ( `ls  *.seed   ` )
set sacfile = ` ls $seednm | cut -d_ -f4,5 | cut -d. -f1 `
mkdir $sacfile
set stnm = ` ls $seednm | awk '{split($0,a,"_");split(a[5],b,".");print a[4]"."b[1]}'`
rdseed << eof
$seednm


d







Y




Y
Quit
eof

mv 20*${stnm}*SAC $sacfile
#mv ${title}*${sacfile}*seed $sacfile
mv RESP.*  $respfile
mv SAC_PZs*  $pzfile
end #seednm

ls -1 $respfile/* | awk '{split($0,a,"'$respfile'/");split(a[2],b,".");print b[2]"."b[3]}' >> sta_resp.lst
ls -1 $pzfile/* | awk '{split($0,a,"'$pzfile'/");split(a[2],b,"_");print b[3]"."b[4]}' >> sta_pz.lst

rm -f rdseed.err_log*
