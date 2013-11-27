//--------------------------------------------------------------------------
/**
 * @File jack_sLength.cpp
 * @brief jack ave error
 * @ingroup YAMADA
 * @author  M.YAMADA 
 * @date    Sat Jun 13 22:09:45 2013
 */
//--------------------------------------------------------------------------


#include <new>
#include "../include/io.h"
#include "../include/analys.h"
#include "../include/jack.h"
#include "../include/fit.h"


using namespace std;

std::string inPath = "/home/sinyamada/results/set1/Spin0-0Bin50/fitPot";
std::string physInfo = "RC16x32_B1830Kud013760Ks013710C1761";
std::string outPath = "/home/sinyamada/results/set1/Spin0-0Bin50/fitPot/jack";
//modify
std::string inStaticsInfo = "1g1y_sph_p";
std::string outStaticsInfo = "1g1y_sph_jack_p";
const int dataSize = 6;
//end



bool inBinary = false;
bool outBinary = false; 

main(){


  IODATA Data;

  Data.setReadBinaryMode(inBinary);
  Data.setWriteBinaryMode(outBinary);
  Data.setConfSize(binnumber);
  Data.setD("1D");

  root_mkdir(outPath.c_str());
  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){

    JACK jackp;
    jackp.set(Confsize,binsize,dataSize);
    
    for (int iconf=0; iconf< binnumber; iconf++) {
      double * y = new double[dataSize]();
       Data.callData(y,2,inPath,inStaticsInfo,physInfo,iconf,iT);
       for(int id = 0;id<dataSize;id++){cout<<"Input Data it ="<<iT<<" "<<y[id]<<endl;}
       jackp.setBinData(y,iconf);
    }
    
    double* xdata = new double[dataSize]();
    double* err= new double[dataSize]();
    double* ave= new double[dataSize]();
    ave = jackp.calcAve();    
    err = jackp.calcErr();    
    for(int id = 0;id<dataSize;id++){cout<<"Output Data it= "<<iT<<" "<<ave[id]<<" "<<err[id]<<endl;}
    Data.outErr(xdata,ave,err,outPath,outStaticsInfo,physInfo,0,iT,dataSize);   
    
    delete [] xdata;
  }//It
  cout <<"@End all jobs"<<endl; 
}


