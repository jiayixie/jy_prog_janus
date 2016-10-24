#!/bin/csh -f
# this version, get parabest for each phi group
# this version, free the vpv, vph and eta
# 
limit datasize 5000000; 
limit stacksize 1000000
setenv OMP_NESTED TRUE
set finponm = point.txt

set fpoint = $argv[1] #point_info_v2.txt: nm lon lat
set num_thread = $argv[2] #10
set vv = $argv[3]
set here = `pwd`

foreach info ( `awk '{printf"%s_%s_%s\n",$1,$2,$3}' $fpoint ` )
echo $info
set nm = `echo $info | cut -d_ -f1`
set lon = `echo $info | cut -d_ -f2`
set lat = `echo $info | cut -d_ -f3`

set id = ${vv}_Tibet_v1_laytheta.20disc.XYdip.aug11.vpvs.3.AZmat4.5perc.negRA.prior #May 4 2016
set exe =  /home/jixi7887/progs/jy/inversion_ElasticTensor_update_fwdcpt_HV
set dirout = /lustre/janus_scratch/jixi7887/Tibet/inv_ET_BS/inv_Apr28data_v1/inv_${id}/${nm}_${lon}_${lat}_inv_${id}
mkdir -p /lustre/janus_scratch/jixi7887/Tibet/inv_ET_BS/inv_Apr28data_v1
mkdir -p /lustre/janus_scratch/jixi7887/Tibet/inv_ET_BS/inv_Apr28data_v1/inv_${id}
mkdir -p $dirout
#set fstartmod = /projects/jixi7887/work/YunNanTibet/Data/input_mod_iso_May4/mod_${lon}_${lat}.AZcm.txt
#set fstartmod = /projects/jixi7887/work/YunNanTibet/Data/input_mod_ani_May8/mod_${lon}_${lat}.AZcm.txt
set fstartmod = /projects/jixi7887/work/YunNanTibet/Data/input_mod_ani_May8/mod_${lon}_${lat}.AZmat4.5perc.txt
set fparanm = $here/para_${vv}.txt
set flagreadVkernel = 0
#--------------
#set dateL = Jun29 # Apr28
#set dateR = Jun23
set dateL = Aug11 #Jul10
set dateR = Jul9
set dirdata = /projects/jixi7887/work/YunNanTibet/Data
set dirRphvel = $dirdata/dispRay_mish_${dateR}_2015
set dirRgpvel = $dirdata/dispRay_mish_${dateR}_2015
set dirLphvel = $dirdata/dispLov_mish_${dateL}_2015
set dirLgpvel = $dirdata/dispLov_mish_${dateL}_2015



#--- do inversion
if ( 1 == 1 )then
	mkdir -p $dirout
	#rm -f $dirout/*
	rm -f $dirout/*disp*
	rm -f $dirout/temp_LMisfit_*
	rm -fr $dirout/bin_avg
	\cp ./run_Mineos_bran.csh $dirout 
	\cp ./run_Mineos_bran_HV.csh $dirout 
	cd $dirout
	pwd
	echo "$nm $lon $lat" > $finponm
	set dirout = .

	echo $dirout
	$exe/invTibet_posRALT10gp12SP_nosign_20disc_XYdip_vpvs_negAni_prior $finponm $dirout $fstartmod $dirRphvel $dirRgpvel $dirLphvel $dirLgpvel $fparanm $flagreadVkernel $num_thread
	#exit

	set dirbin = $dirout/binmod
	set diroutbin = $dirout/bin_avg
	mkdir -p $diroutbin
        $exe/CALavg_Tibet_20disc_XYdip $finponm $dirbin $diroutbin 1 $dirout $fstartmod $dirRphvel $dirRgpvel $dirLphvel $dirLgpvel $fparanm 1
        $exe/CALavg_Tibet_20disc_XYdip $finponm $dirbin $diroutbin 2 $dirout $fstartmod $dirRphvel $dirRgpvel $dirLphvel $dirLgpvel $fparanm 1

	rm -f  MineosInputMod*
	wc temp_LMisfit*.txt | awk '{if($1==0)printf"rm -f %s\n",$4}'  | csh
	cd $here
endif
end
