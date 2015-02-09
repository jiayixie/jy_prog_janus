#this is used to rotate NE to rt and rename Z to z
#Also cut data into 0~te, the te should be the same value as that used in do_deciamte.csh.

set main = `pwd`
set dirot = sacRot
#mkdir $dirot
set tmin = 130
set dirsac = sacData_LE5.7

foreach sacdir ( `ls -d $dirsac/* `  )
#foreach sacdir ( `more rot_Input_errSta.txt` )
cd $main/$sacdir
set nt_st = ` echo $sacdir | cut -d"/"  -f2`
echo $nt_st
#rm -f -r $main/$dirot/$nt_st
mkdir $main/$dirot/$nt_st
echo $nt_st
sac << eof
m $main/rot.m dirot $main/$dirot/$nt_st te $tmin
q
eof

end
