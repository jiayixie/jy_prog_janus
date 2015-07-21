set name = "_2040disc"
#set name = "_20disc_X0dip"
#set name = "_"
set sub = $1
if ( $name == "_" ) then
	set name1 = "_"
	set name2 = ""
else
	set name1 = ${name}_
	set name2 = $name
endif

echo invTibet_posRAgp12SP_nosign${name1}vpvs CALavg_Tibet${name2} ${sub}
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV_Tibet.C -std=c++11 -o  invTibet -g
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_Tibet_v1.C -std=c++11 -o  invTibet_NposRAgp_f2o4_prior_$argv[1] -g
g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_Tibet_v1.C -std=c++11 -o  invTibet_posRAgp12SP_nosign${name1}vpvs.prior${sub}   -g 
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_Tibet_v1.iso.C  -std=c++11 -o  invTibet_NposRAgp_nosign_TI  -g 

g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp   CALavg_getposteria_cptLkernel_HV_Tibet.C -std=c++11 -o  CALavg_Tibet${name2}.prior${sub} -g
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp   CALavg_getposteria_cptLkernel_HV_test.C -std=c++11 -o  CALavg_test -g
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp CALavg_getposteria_cptLkernel_HV_Tibet.iso.C  -std=c++11 -o  CALavg_Tibet_TI -g

#--- another fwd cpt code,
#g++ -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp  do_fwdsyn_Tibet.C -std=c++11 -o do_fwdsyn_Tibet -g

#--- yet another one for TTI model
#g++ -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp  do_fwdsyn_tti.C -std=c++11 -o do_fwdsyn_tti -g

