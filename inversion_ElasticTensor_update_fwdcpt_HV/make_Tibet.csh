#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV_Tibet.C -std=c++11 -o  invTibet -g
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_Tibet_v1.C -std=c++11 -o  invTibet_posRAgp2 -g

g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp   CALavg_getposteria_cptLkernel_HV_Tibet.C -std=c++11 -o  CALavg_Tibet -g

#--- another fwd cpt code,
#g++ -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp  do_fwdsyn_Tibet.C -std=c++11 -o do_fwdsyn_Tibet -g
