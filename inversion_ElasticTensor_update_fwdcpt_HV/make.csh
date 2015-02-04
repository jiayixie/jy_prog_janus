#----make the main code--
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV.C -std=c++11 -o  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV_RLHVnoposani -g

#-----make the CALavg code --
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_cptLkernel_HV.C -o CALavg_getposteria_cptLkernel_RL -std=c++11  -g

#--- make the fwd cpt test code --
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  test_fwdHV.C -std=c++11 -o test_fwdHV -g

#--- make the sensitivity computation code --
g++  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -std=c++11 -fopenmp  compute_sensitivity.C -o compute_sensitivity  -g

