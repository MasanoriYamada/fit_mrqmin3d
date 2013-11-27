//--------------------------------------------------------------------------
/**
 * @Filex   fit_pot.1gy.cpp
 * @brief   fiting potential gauss + yukawa
 * @ingroup YAMADA
 * @author  M.YAMADA 
 * @date    Tue Nov 5 2013 
 */
//--------------------------------------------------------------------------


#include "../include/io.h"
#include "../include/analys.h"
#include "../include/jack.h"
#include "../include/func_1g1y3D.h"
#include "../include/func_1g1y1D.h"
#include "../include/func_1g1yy3D.h"
#include "../include/func_1g1yy1D.h"


using namespace std;

std::string		inPath	      = "/home/sinyamada/results/set1/Spin0-0Bin50/Potential/binPotential/xyz";
std::string		outPath	      = "/home/sinyamada/results/set1/Spin0-0Bin50/fitPot";
std::string             physInfo      = "RC16x32_B1830Kud013760Ks013710C1761";
std::string		inStaticsInfo = "Potential";
std::string		osi1g1yy1D    = "1g1yy1D";
std::string		osi1g1y1D     = "1g1y1D";
bool			inBinary      = true;
bool			outBinary     = false;
static const int	in_line	      = 1;	//which line read at inpufile ?

int main(){

  IODATA	io;
  io.setReadBinaryMode(inBinary);
  io.setWriteBinaryMode(outBinary);
  io.setConfSize(binnumber);


  for(int iT=T_in  ; iT< T_fi +1  ; iT++ )
    {
      //calc Jack error for using fit 
      JACK	jackPot;
      jackPot.set(Confsize,binsize,XYZnodeSites);
      COMPLEX*	binPot = new COMPLEX[XYZnodeSites*binnumber]();
      
    for (int iconf=0; iconf< binnumber; iconf++) 
      {
	COMPLEX*	pot = new COMPLEX[XYZnodeSites]();
	io.callData(pot,in_line,inPath,inStaticsInfo,physInfo,iconf,iT);
	reScale(pot);
	memcpy(binPot + iconf*XYZnodeSites, pot, sizeof(pot)*XYZnodeSites*2);    
	jackPot.setBinData(pot,iconf);
	
	delete [] pot;
      }
    
    double*	err;
    double*	ave;

    err	= jackPot.calcErr();
    ave	= jackPot.calcAve();

    //for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){cout<<"ave "<<ixyz<<" "<<ave[ixyz]<<" "<<err[ixyz]<<endl;}
    //inPot.outErr(ave,err,outPath,outStaticsInfo,physInfo,700,iT,XYZnodeSites);
    
    //fit part
    
    Fit* fit1g1y3D = new Fit(new Fit1g1y3D());
    Fit* fit1g1y1D = new Fit(new Fit1g1y1D());
    Fit* fit1g1yy3D = new Fit(new Fit1g1yy3D());
    Fit* fit1g1yy1D = new Fit(new Fit1g1yy1D());

    double	chisq;
    cout << "@Start fit t ="<<iT <<"fit function = "<< fit1g1y3D->getName()<<endl;
    
    for (int iconf=0; iconf< binnumber; iconf++) 
      {
	double*	binPotIn = new double[XYZnodeSites]();
	
	for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++)
	  {
	    binPotIn[ixyz] = binPot[ixyz + iconf*XYZnodeSites].real();
	  }
	
	double* a1 = new double[fit1g1y3D->getDofA()];
	double* a2 = new double[fit1g1y1D->getDofA()];
	double* a3 = new double[fit1g1yy3D->getDofA()];
	double* a4 = new double[fit1g1yy1D->getDofA()];
	/*
	try{	
	fit1g1y3D->setRange(0,7);
	fit1g1y3D->fit(binPotIn, err, a1,chisq);
	cout<<iconf<<"@param " <<a1[0]<<" "<<a1[1]<<" "<<a1[2]<<" "<<a1[3]<<" "<<a1[4]<<" chisq = "<<chisq<<endl;
	}
	catch (const char* error1)      {
	  cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
	  cerr<<endl;
	}
	*/
	try{
	fit1g1y1D->setRange(0,7);
	fit1g1y1D->fit(binPotIn, err, a2,chisq);
	cout<<iconf<<"@param " <<a2[0]<<" "<<a2[1]<<" "<<a2[2]<<" "<<a2[3]<<" "<<a2[4]<<" chisq = "<<chisq<<endl;
	}
	catch (const char* error1)      {
	  cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
	  cerr<<endl;
	}
	/*
	try{
	fit1g1yy3D->setRange(0,7);
	fit1g1yy3D->fit(binPotIn, err, a3,chisq);
	cout<<iconf<<"@param " <<a3[0]<<" "<<a3[1]<<" "<<a3[2]<<" "<<a3[3]<<" "<<a3[4]<<" chisq = "<<chisq<<endl;
	}
	catch (const char* error1)      {
	  cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
	  cerr<<endl;
	}
	*/
	try{
	fit1g1yy1D->setRange(0,7);
	fit1g1yy1D->fit(binPotIn, err, a4,chisq);
	cout<<iconf<<"@param " <<a4[0]<<" "<<a4[1]<<" "<<a4[2]<<" "<<a4[3]<<" "<<a4[4]<<" chisq = "<<chisq<<endl;
	}
	catch (const char* error1)      {
	  cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
	  cerr<<endl;
	}

	//double param[6] =	{a[0],a[1],a[2],a[3],a[4],chisq};


	delete [] binPotIn;
	
	io.outData(a2,outPath,osi1g1y1D,physInfo,iconf,iT,fit1g1y1D->getDofA()+1);
	io.outData(a4,outPath,osi1g1yy1D,physInfo,iconf,iT,fit1g1yy1D->getDofA()+1);
    }				//iconf  
    delete [] binPot;
    
  }				//It
  cout <<"@End all jobs"<<endl; 
}


