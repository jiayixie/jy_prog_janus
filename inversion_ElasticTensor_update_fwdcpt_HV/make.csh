#----make the main code--
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV.C -std=c++11 -o  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV_RLposani -g
g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  Main_do_inversion_parallel_BS_multipleGp_updateKeachjump_changeEtaSpace_parallel_cptLkernel_HV.C -std=c++11 -o  test_wtvpconstr_smsign_mono01_grd1_vpvsGT1.7_BS5f5o4_newgenpara_L_posRA1 -g

#-----make the CALavg code --
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_cptLkernel_HV.C -o CALavg_getposteria_cptLkernel_RLHV_BS5f5o4_newgenpara_L -std=c++11  -g
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_cptLkernel_HV_test.C -o CALavg_getposteria_cptLkernel_RLHV_BS5f5o4_newgenpara_L_test -std=c++11  -g

#--- make the fwd cpt test code --
#g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp  test_fwdHV.C -std=c++11 -o test_fwdHV -g
#--- another fwd cpt code,
#g++ -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp  do_fwdsyn.C -std=c++11 -o do_fwdsyn 

#--- make the sensitivity computation code --
#g++  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -std=c++11 -fopenmp  compute_sensitivity.C -o compute_sensitivity  -g



#--------

