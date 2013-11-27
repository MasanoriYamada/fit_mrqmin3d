#ifndef ANALYS_H_
#define ANALYS_H_

#include <sys/stat.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>


static const int T_in=7;
static const int T_fi=7;
static const int XnodeSites =16;
static const int YnodeSites =16;
static const int ZnodeSites =16;
static const int TnodeSites =32/2;
static const int binsize=50;
static const int Confsize=700;
//--------------------------------------------------//
static const int binnumber=(Confsize/binsize);
static const int XYZnodeSites =XnodeSites*XnodeSites*XnodeSites;
static std::string dir_path;
static std::string in_dir_path;
static std::string out_dir_path;




inline void root_mkdir(const char* str){
if(mkdir(str,0775)!=0){std::cout<<"Warning	cat't make dir	already exist::"<<str<<std::endl;}
}

#endif
