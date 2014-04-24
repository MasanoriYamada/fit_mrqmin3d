//--------------------------------------------------------------------------
/**
 * @Filex   fit_pot.1gy.cpp (periodic boundary condition)
 * @brief   fiting potential gauss + yukawa^2
 * @ingroup YAMADA
 * @author  M.YAMADA
 * @date    Tue Nov 5 2013
 */
//--------------------------------------------------------------------------

#include "../../include/analys.h"
#include "../../include/func_1g1yy1D.h"



using namespace std;


//---------------------fit for 3D---------------------------------------//
void Func1g1yy1D(double& x,double a[],double& yfit,double dyda[],int& ma)
{
  double        b1 = a[0];
  double        b2 = a[1];
  double        b3 = a[2];
  double        b4 = a[3];
  double        b5 = a[4];

  if ( x == 0 ){yfit = b1;}
  else if ( x > 0) {
    yfit = b1*exp(-b2*x*x)+ b3*((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x))/x)*(exp(-b5*x))*(exp(-b5*x)/x);
  }
  else {cout<<"ERR:: fit data x is warng"<<endl;}

  if ( x == 0 ){
    dyda[0] = 1;
    dyda[1] = 0;
    dyda[2] = 0;
    dyda[3] = 0;
    dyda[4] = 0;
  }
  else if ( x > 0) {
    dyda[0] = exp(-b2*x*x);
    dyda[1] = -b1*x*x*exp(-b2*x*x);
    dyda[2] = ((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x)))*(exp(-b5*x)/x)*(exp(-b5*x)/x);
    dyda[3] = 2*b3*(1-exp(-b4*x*x))*exp(-b4*x*x-2*b5*x);
    dyda[4] = -2*b3*(1 - exp(-b4*x*x))*(1 - exp(-b4*x*x))*(exp(-2*b5*x))/x;
  }
  else {cout<<"ERR:: fit data x is warng"<<endl;}

}

void  Fit1g1yy1D::fit(double* datain, double* datain_sigma, double* b, double& chisq)
{

#define datain(ix,iy,iz)        datain[ix+XnodeSites*(iy+YnodeSites*(iz))]
#define datain_sigma(ix,iy,iz)  datain_sigma[ix+Nx*(iy+Ny*(iz))]
  double*       x  = new double[fit_lengh*fit_lengh*fit_lengh]();
  double*       y  = new double[fit_lengh*fit_lengh*fit_lengh]();
  double*       dy = new double[fit_lengh*fit_lengh*fit_lengh]();


  int ixyz = 0;
  for(int iz = r_in_; iz < r_fi_ + 1; iz++){
    for(int iy = r_in_; iy < r_fi_ + 1; iy++){
      for(int ix = r_in_; ix < r_fi_ + 1; ix++){

	x[ixyz]  = ascale*sqrt(ix*ix + iy*iy + iz*iz);
	y[ixyz]  = datain[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	dy[ixyz] = datain_sigma[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	ixyz++;
      }
    }
  }
                           //ixyz
  int                 ndata      = fit_lengh*fit_lengh*fit_lengh;
  static double       a[Nparam];
  static int          is_initial = 1;
  int   ia[]                     = {1,1,1,1,1}, ma =Nparam;
  double              covar[Nparam*Nparam],alpha[Nparam*Nparam];
  int                 nca        = Nparam;
  double              alambda;

  if(is_initial==1){
    is_initial = 0;

    //set1 it07
   // a[0] = 1148.23;
   // a[1] = 53.1;
   // a[2] = -333.105;
   // a[3] = 2.51011;
   // a[4] = 2.23122;
    //kappa_013700.013640 it11
    a[0] = 1479.5 ;  
    a[1] = 47.1749 ;  
    a[2] = -470.246;  
    a[3] = 1.24333 ;  
    a[4] = 1.03994 ;  
  }
  

  // The initial call
  alambda = -1.0;
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy1D, alambda);

  // iterations
  for(int iter=0; iter<65536; iter++){
    mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy1D, alambda);
    //cout<<chisq<<" "<<alambda<<endl;
    //cout<<"@param " <<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" "<<a[4]<<" chisq = "<<chisq<<endl;
    if (alambda > 1.0e64) break;
  }
  if (alambda <= 1.0e64) {
    //cerr << "convergence is not achieved"<<endl;
    throw (fitname_.c_str());
  }

  // The final call

  alambda = 0.0;
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy1D, alambda);
  //debug
  //  for (int id=0 ;id< 5; id++) cout << "param"<< a[id];
  //cout<<" "<<chisq<<endl;
  //debug end
  for (int id = 0; id < Nparam ;id++ )
    {
      b[id] = a[id];
    }

  chisq	  = chisq /(ndata - Nparam - 1);

  delete [] x;
  delete [] y;
  delete [] dy;
}
