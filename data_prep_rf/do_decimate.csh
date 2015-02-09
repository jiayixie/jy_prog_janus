#execute the decimate program. decimate SAC data to appropriate npts, and mv SAC data with incompelete 3 comp or incompelete header.
mkdir sacBad
set dirlong = ./sacBad/long_short
set dirhead = ./sacBad/no_head
set dirincmp = ./sacBad/incmp
set sacdirall = ./sacData
set saclst = saclist
set tmin = 130 #this would be the cut_end in rotate_step
mkdir $dirlong $dirhead $dirincmp
foreach dir ( `ls -d $sacdirall/*` )
ls $dir/*BHZ*   > $saclst
./decimate $saclst $dirlong $dirhead $dirincmp $tmin
end
