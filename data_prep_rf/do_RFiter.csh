set here = `pwd`
set RFdir = $here/RF_iter_deresp
set rotdirall = sacRot
mkdir $RFdir
foreach rotdir ( `ls -d $rotdirall/MagLE5.7/*_* | head -220` )
set nt_st = `echo $rotdir | cut -d"/" -f3`
echo $nt_st
mkdir $RFdir/$nt_st
cd $here/$rotdir
foreach sacr ( `ls *.r ` )
set sact = `echo $sacr | awk '{sub(".r",".t",$0);print $0}'`
set sacz = `echo $sacr | awk '{sub(".r",".z",$0);print $0}'`
set outr = `echo $sacr | awk '{sub(".r",".eqr",$0);print $0}'`
iterdecon <<eof
$sacr
$sacz 
100      
5.0     
0.001   
2.5      
1      
0      
eof
#
# rename the result (decon.out) to a more useful value
#
\mv decon.out $RFdir/$nt_st/$outr
#rm -f 
end
end
