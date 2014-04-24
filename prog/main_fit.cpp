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

//boost
#include <boost/lexical_cast.hpp>

using namespace std;
using boost::lexical_cast;
//set1
/*
std::string		inPath	      = "/home/sinyamada/results/set1/Spin0-0Bin50/Potential/binPotential/xyz";
std::string		outPath	      = "/home/sinyamada/results/set1/Spin0-0Bin50/fitPot/bin";
std::string             physInfo      = "RC16x32_B1830Kud013760Ks013710C1761";
std::string		inStaticsInfo = "Potential";
std::string		osi1g1yy1D    = "1g1yy1D";
std::string		osi1g1y1D     = "1g1y1D";
bool			inBinary      = true;
bool			outBinary     = false;
static const int	in_line	      = 1;	//which line read at inpufile ?
*/

//kappa_013700.013640
std::string		inPath	      = "/home/sinyamada/results/kappa_013700.013640/ts32/spin00.bin"
  + lexical_cast<string>(binsize) + "/Potential/binPotential/xyz";
std::string		outPath	      = "/home/sinyamada/results/kappa_013700.013640/ts32/spin00.bin"
  + lexical_cast<string>(binsize) + "/fitPot/bin";
std::string		outAvePath    = "/home/sinyamada/results/kappa_013700.013640/ts32/spin00.bin"
  + lexical_cast<string>(binsize) + "/ave/fitPot/bin";
std::string             physInfo      = "RC32x64_B1900Kud01370000Ks01364000C1715";
std::string		inStaticsInfo = "Potential";
std::string		osi1g1yy1D    = "1g1yy1D";
std::string		osi1g1y1D     = "1g1y1D";
std::string		osi1g1yy3D    = "1g1yy3D";
std::string		osi1g1y3D     = "1g1y3D";
bool			inBinary      = true;
bool			outBinary     = false;
static const int	in_line	      = 1;	//which line read at inpufile ?


int main(){

  IODATA	io;
  io.setReadBinaryMode(inBinary);
  io.setWriteBinaryMode(outBinary);
  io.setConfSize(binnumber);

  root_mkdir(outPath.c_str());
  root_mkdir(outAvePath.c_str());

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

    for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){cout<<"ave "<<ixyz<<" "<<ave[ixyz]<<" "<<err[ixyz]<<endl;}
    //inPot.outErr(ave,err,outPath,outStaticsInfo,physInfo,700,iT,XYZnodeSites);


    //-------------------------------------------//
    //---------------fit part -------------------//
    //-------------------------------------------//

    cout << "@Start ave fit t ="<<iT <<endl;


    Fit* fit1g1y3D = new Fit(new Fit1g1y3D());
    Fit* fit1g1y1D = new Fit(new Fit1g1y1D());
    Fit* fit1g1yy3D = new Fit(new Fit1g1yy3D());
    Fit* fit1g1yy1D = new Fit(new Fit1g1yy1D());
    double	chisq = 0.0;


    //------------potential ave fit for test -----------------//
    double* paAvefit1g1y3D = new double[fit1g1y3D->getDofA()+1];
    double* paAvefit1g1y1D = new double[fit1g1y1D->getDofA()+1];
    double* paAvefit1g1yy3D = new double[fit1g1yy3D->getDofA()+1];
    double* paAvefit1g1yy1D = new double[fit1g1yy1D->getDofA()+1];
 
    try{
      fit1g1y3D->setRange(0,15);
      fit1g1y3D->fit(ave, err, paAvefit1g1y3D,chisq);
      cout<<"@ave param "<<fit1g1y3D->getName()<<" " <<paAvefit1g1y3D[0]<<" "<<paAvefit1g1y3D[1]<<" "<<paAvefit1g1y3D[2]<<" "<<paAvefit1g1y3D[3]<<" "<<paAvefit1g1y3D[4]<<" chisq = "<<chisq<<endl;
      paAvefit1g1y3D[5]=chisq;
    }
    catch (const char* error1)	    {
      cerr <<error1<<"ave ERR convergence is not achieved"<<endl;
      cerr<<endl;
    }
 
 
    try{
      fit1g1y1D->setRange(0,15);
      fit1g1y1D->fit(ave, err, paAvefit1g1y1D,chisq);
      cout<<"@ave param " <<fit1g1y1D->getName()<<" " <<paAvefit1g1y1D[0]<<" "<<paAvefit1g1y1D[1]<<" "<<paAvefit1g1y1D[2]<<" "<<paAvefit1g1y1D[3]<<" "<<paAvefit1g1y1D[4]<<" chisq = "<<chisq<<endl;
      paAvefit1g1y1D[5]=chisq;
    }
    catch (const char* error1)	    {
      cerr <<error1<<"ave ERR convergence is not achieved"<<endl;
      cerr<<endl;
	}
 
 
    try{
      fit1g1yy3D->setRange(0,15);
      fit1g1yy3D->fit(ave, err, paAvefit1g1yy3D,chisq);
      cout<<"@ave param " <<fit1g1yy3D->getName()<<" " <<paAvefit1g1yy3D[0]<<" "<<paAvefit1g1yy3D[1]<<" "<<paAvefit1g1yy3D[2]<<" "<<paAvefit1g1yy3D[3]<<" "<<paAvefit1g1yy3D[4]<<" chisq = "<<chisq<<endl;
      paAvefit1g1yy3D[5]=chisq;
    }
    catch (const char* error1)	    {
      cerr <<error1<<"ave ERR convergence is not achieved"<<endl;
      cerr<<endl;
    }
 
    try{
	fit1g1yy1D->setRange(0,15);
	fit1g1yy1D->fit(ave, err, paAvefit1g1yy1D,chisq);
	cout<<"@ave param " <<fit1g1yy1D->getName()<<" " <<paAvefit1g1yy1D[0]<<" "<<paAvefit1g1yy1D[1]<<" "<<paAvefit1g1yy1D[2]<<" "<<paAvefit1g1yy1D[3]<<" "<<paAvefit1g1yy1D[4]<<" chisq = "<<chisq<<endl;
	paAvefit1g1yy1D[5]=chisq;
    }
    catch (const char* error1)	    {
      cerr <<error1<<"ave ERR convergence is not achieved"<<endl;
      cerr<<endl;
    }
 
	io.outData(paAvefit1g1y3D ,outAvePath,osi1g1y3D,physInfo,Confsize,iT,fit1g1y3D->getDofA()+1);
	io.outData(paAvefit1g1y1D ,outAvePath,osi1g1y1D,physInfo,Confsize,iT,fit1g1y1D->getDofA()+1);
	io.outData(paAvefit1g1yy3D,outAvePath,osi1g1yy3D,physInfo,Confsize,iT,fit1g1yy3D->getDofA()+1);
	io.outData(paAvefit1g1yy1D,outAvePath,osi1g1yy1D,physInfo,Confsize,iT,fit1g1yy1D->getDofA()+1);
 
 
	delete []  paAvefit1g1y3D  ;
	delete []  paAvefit1g1y1D  ;
	delete []  paAvefit1g1yy3D ;
	delete []  paAvefit1g1yy1D ;
 

    //------------main fit (each conf)-----------------//

    cout << "@Start main fit t ="<<iT <<endl;

    for (int iconf=0; iconf< binnumber; iconf++)
      {
	double*	binPotIn = new double[XYZnodeSites]();

	for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++)
	  {
	    binPotIn[ixyz] = binPot[ixyz + iconf*XYZnodeSites].real();
	  }

	double* a1 = new double[fit1g1y3D->getDofA()+1];
	double* a2 = new double[fit1g1y1D->getDofA()+1];
	double* a3 = new double[fit1g1yy3D->getDofA()+1];
	double* a4 = new double[fit1g1yy1D->getDofA()+1];

// 	  try{
// 	  fit1g1y3D->setRange(0,15);
// 	  fit1g1y3D->fit(binPotIn, err, a1,chisq);
// 	  cout<<iconf<<"@param " <<a1[0]<<" "<<a1[1]<<" "<<a1[2]<<" "<<a1[3]<<" "<<a1[4]<<" chisq = "<<chisq<<endl;
// 	  a1[5] = chisq;
// 	  }
// 	  catch (const char* error1)	  {
// 	  cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
// 	  cerr<<endl;
// 	  }
// 	  io.outData(a1,outPath,osi1g1y3D,physInfo,iconf,iT,fit1g1y1D->getDofA()+1);
// 
// 
// 	  try{
// 	    fit1g1y1D->setRange(0,15);
// 	    fit1g1y1D->fit(binPotIn, err, a2,chisq);
// 	    cout<<iconf<<"@param " <<a2[0]<<" "<<a2[1]<<" "<<a2[2]<<" "<<a2[3]<<" "<<a2[4]<<" chisq = "<<chisq<<endl;
// 	    a2[5] = chisq;
// 	  }
// 	  catch (const char* error1)	{
// 	    cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
// 	    cerr<<endl;
// 	  }
// 	  io.outData(a2,outPath,osi1g1y1D,physInfo,iconf,iT,fit1g1y1D->getDofA()+1);	  
//  
// 	  try{
// 	    fit1g1yy3D->setRange(0,15);
// 	    fit1g1yy3D->fit(binPotIn, err, a3,chisq);
// 	    cout<<iconf<<"@param " <<a3[0]<<" "<<a3[1]<<" "<<a3[2]<<" "<<a3[3]<<" "<<a3[4]<<" chisq = "<<chisq<<endl;
// 	    a3[5] = chisq;
// 	  }
// 	  catch (const char* error1)	{
// 	    cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
// 	    cerr<<endl;
// 	  }
// 	io.outData(a3,outPath,osi1g1yy3D,physInfo,iconf,iT,fit1g1yy1D->getDofA()+1);
//  
// 	  try{
// 	    fit1g1yy1D->setRange(0,15);
// 	    fit1g1yy1D->fit(binPotIn, err, a4,chisq);
// 	    cout<<iconf<<"@param " <<a4[0]<<" "<<a4[1]<<" "<<a4[2]<<" "<<a4[3]<<" "<<a4[4]<<" chisq = "<<chisq<<endl;
// 	    a4[5] = chisq;
// 	  }
// 	  catch (const char* error1)	{
// 	    cerr <<error1<<"  CONF = " << iconf<<" ERR convergence is not achieved"<<endl;
// 	  cerr<<endl;
// 	  }
// 
// 	  
	  
	  delete [] binPotIn;

    }				//iconf
    delete [] binPot;

  }				//It
  cout <<"@End all jobs"<<endl;
}
