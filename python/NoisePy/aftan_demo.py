import obspy
import noisepy as npy
import matplotlib.pyplot as plt
import numpy as np
# f='/lustre/janus_scratch/life9360/ALASKA_COR/COR/SAW/COR_SAW_BHZ_SII_BHZ.SAC'
# f1='/lustre/janus_scratch/life9360/SES3D_WorkingDir/OUTPUT_SAC_ak135/LF.S005032..BXZ.SAC'
f1='/projects/life9360/SYN_TEST_ARTIE/DATA/TOMO_1D.ak135/SAC_Z/S017021.SS.BXZ.sem.sac'
f1='/lustre/janus_scratch/life9360/EA_ses3d_working_dir/OUTPUT_SAC_EA_10sec_1km_001/LF.EA123S32..BXZ.SAC'
# f2='/projects/life9360/SYN_TEST_ARTIE/DATA/TOMO_1D.ak135/SAC_Z/S005035.SS.BXZ.sem.sac'
# f1='INSTASEIS.LXZ.SAC'
# f2='/projects/life9360/SYN_TEST_ARTIE/DATA/TOMO_1D.ak135/SAC_Z/S019025.SS.BXZ.sem.sac'
# f2='/lustre/janus_scratch/life9360/SES3D_WorkingDir/OUTPUT_SAC/LF.S005032..BXZ.SAC'
# f2='/projects/life9360/SYN_TEST_ARTIE/DATA/TOMO_1D.ak135/SAC_Z/S006001.SS.BXZ.sem.sac'
# f2='/lustre/janus_scratch/life9360/SES3D_WorkingDir/OUTPUT_SAC_vel/LF.S002008..BXZ.SAC'
# f2='/lustre/janus_scratch/life9360/SES3D_WorkingDir/OUTPUT_SAC/LF.S005020..BXZ.SAC'

#prefname='/lustre/janus_scratch/life9360/EA_ses3d_working_dir/OUTPUT_DISP_EA_10sec_1km_001/PREPHASE_R/NKNT.EA123S32.pre'
prefname='/home/life9360/CLCO.N25K.pre'
#f='~/COR_TA.O11A_TA.Q21A.SAC_TT_s'
f='/home/jixi7887/COR_TA.O11A_US.WUAZ.SAC_TT_s'
st=obspy.read(f)
tr=st[0]
tr=npy.noisetrace(tr.data, tr.stats)
tr.aftan(piover4=-1., pmf=True, tmin=2.0, tmax=50.0, phvelname=prefname)
#tr1.aftan(piover4=0, pmf=True, tmin=2.0, tmax=50.0, phvelname=prefname)
#tr1.ftanparam.FTANcomp(tr.ftanparam)
# 
# st=obspy.read(f2)
# tr=st[0]
# tr2=npy.noisetrace(tr.data, tr.stats)
# tr2.aftan(piover4=0., pmf=True, tmin=2.0, tmax=100.0,phvelname='ak135.disp')
tr.ftanparam.writeDISP(f)
tr.getSNR (tmin=4.0, tmax=70.0)
tr.SNRParam.writeAMPSNR(f)
#tr.ftanparam.FTANcomp(tr2.ftanparam)
# 
# InAK135Arr=np.loadtxt('prem.disp')
# T=InAK135Arr[:,0];
# V=InAK135Arr[:,1];
# ax = plt.subplot()
# ax.plot(T, V, '--r', lw=3);
# # 
# InAK135Arr=np.loadtxt('ak135.disp')
# T=InAK135Arr[:,0];
# V=InAK135Arr[:,1];
# ax = plt.subplot()
# ax.plot(T, V, '--g', lw=3);
tr.plotftan()
plt.show()


