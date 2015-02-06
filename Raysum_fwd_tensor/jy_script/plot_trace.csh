# this is used to plot the traces with gmt
#
# step 1: read the outpu .tr file, output the azi, and Amp(azi,t) for each tace
# step 2: for each trace, plot it with pswiggle
#

set title = test
set ftr = sample.tr
set fgeom = sample.geom
set outdir =  test
set outps = test.ps
set comp = Z

mkdir -p $outdir
python get_traces.py $fgeom $ftr $outdir

echo plotting $outps
set Ntrace = `awk '{if(NR==2)print $1}' $ftr`
set N = 1

set REG = -R0/60/-10/380 
set SCA = -JX10/20
set Z = -Z0.2
set p = 0.5

while ( $N <= $Ntrace ) 
	echo trace $N
	set ftrace = $outdir/trace_${N}.${comp}.txt

	if ( $N == 1 ) then 
		pswiggle $ftrace $REG $SCA $Z -Ba20f10:"time (sec)":/a90f10:"azimuth (deg)"::."$title":WSen -Gred  -T0.5p/black -K -P -Y 5  -W${p}p/black > $outps
	else
		pswiggle $ftrace -R -J $Z -Gred -T0.5p/black -K -O -W${p}p/black >> $outps
	endif

	@ N = $N + 1
end


pwd | psxy -R -J -O >> $outps

 
