#include<stdio.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;
#include "./generate_Bs.C"

int main(int argc, char *argv[]){
vector<double> Bspline;
int factor = atoi(argv[1]);
int order = atoi(argv[2]);
//int nBs=4,order=4,thick=100;
int nBs=5;
printf("factor=%d order=%d\n",factor,order);
gen_B_spline(nBs,order,0.,100.,factor,100,Bspline);
return 1;
}

