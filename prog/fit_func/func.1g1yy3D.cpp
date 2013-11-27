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
#include "../../include/func_1g1yy3D.h"



using namespace std;

inline double y_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){return b1;}
  else if (0 < x) {
   return b1*exp(-b2*x*x)+ b3*((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x))/x)*(exp(-b5*x))*(exp(-b5*x)/x);
  }
  else if (x < 0) {
   return b1*exp(-b2*x*x)+ b3*((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x))/x)*(exp(-b5*x))*(exp(-b5*x)/x);
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }

}

inline double dy1_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){
    return 1;
  }
  else if ( x != 0) {
    return exp(-b2*x*x);
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }
}
inline double dy2_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){
    return 0;
  }
  else if ( x != 0) {
    return -b1*x*x*exp(-b2*x*x);
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }
}
inline double dy3_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){
    return 0;
  }
  else if(0 < x){
    return ((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x)))*(exp(-b5*x)/x)*(exp(-b5*x)/x);
  }
  else  if (x < 0){
    return  ((1 - exp(-b4*x*x)))*((1 - exp(-b4*x*x)))*(exp(-b5*x)/x)*(exp(-b5*x)/x);
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }
}
inline double dy4_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){
    return 0;
  }
  else if(0 < x){
    return 2*b3*(1-exp(-b4*x*x))*exp(-b4*x*x-2*b5*x);
  }
  else  if (x < 0){
    return 2*b3*(1-exp(-b4*x*x))*exp(-b4*x*x-2*b5*x);
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }
}
inline double dy5_function(double x, double b1, double b2, double b3, double b4, double b5 ){
  if ( x == 0 ){
    return 0;
  }
  else if(0 < x){
    return -2*b3*(1 - exp(-b4*x*x))*(1 - exp(-b4*x*x))*(exp(-2*b5*x))/x;
  }
  else  if (x < 0){
    return -2*b3*(1 - exp(-b4*x*x))*(1 - exp(-b4*x*x))*(exp(-2*b5*x))/x;
  }
  else {
    cout<<"ERR:: fit data x is warng"<<endl;
    exit(0);
  }
}

//---------------------fit function---------------------------------------//

//---------------------fit for 3D---------------------------------------//

triple *xyz_1g1yy3D	   = NULL;
triple *origins_1g1yy3D	   = NULL;
const int Norigins_1g1yy3D = 3*3*3;

void Func1g1yy3D(double x,double y,double z, double a[],double& yfit, double dyda[],int& ma)
{
  double r = sqrt(x*x + y*y + z*z);
  double yfit_tmp;
  double dyda_tmp[ma];

      double	b1 = a[0];
      double	b2 = a[1];
      double	b3 = a[2];
      double	b4 = a[3];
      double	b5 = a[4];
  
      yfit_tmp    = y_function(r,b1,b2,b3,b4,b5);
      dyda_tmp[0] = dy1_function(r,b1,b2,b3,b4,b5);
      dyda_tmp[1] = dy2_function(r,b1,b2,b3,b4,b5);
      dyda_tmp[2] = dy3_function(r,b1,b2,b3,b4,b5);
      dyda_tmp[3] = dy4_function(r,b1,b2,b3,b4,b5);
      dyda_tmp[4] = dy5_function(r,b1,b2,b3,b4,b5);


  yfit += yfit_tmp;
  for(int i = 0; i < ma; i++){
    dyda[i] += dyda_tmp[i];
  }
} 

inline void Func1g1yy3D(const triple& x,const triple& x0,
		 double	 a[],
		 double& yfit,
		 double	 dyda[],
		 int&	 ma)
{
  Func1g1yy3D(x.x - x0.x,
       x.y - x0.y,
       x.z - x0.z,
       a, yfit, dyda, ma);
}
void Func1g1yy3D(double& index,double a[],double& yfit,double dyda[],int& ma)
{
  int n = ((int*)&index)[0];
  yfit = 0.0;
  for(int i = 0; i < ma; i++) dyda[i]=0.0;
  for(int i = 0; i < Norigins_1g1yy3D; i++){
    Func1g1yy3D(xyz_1g1yy3D[n], origins_1g1yy3D[i], a,yfit,dyda,ma);
  }
  
}

void  Fit1g1yy3D::fit(double datain[], double datain_sigma[], double* b, double& chisq)
{
  
#define datain(ix,iy,iz)	datain[ix+XnodeSites*(iy+YnodeSites*(iz))]
#define datain_sigma(ix,iy,iz)	datain_sigma[ix+Nx*(iy+Ny*(iz))]
  double*	x  = new double[fit_lengh*fit_lengh*fit_lengh]();
  double*	y  = new double[fit_lengh*fit_lengh*fit_lengh]();
  double*	dy = new double[fit_lengh*fit_lengh*fit_lengh]();
  xyz_1g1yy3D		   = new triple[fit_lengh*fit_lengh*fit_lengh]();
  int ixyz = 0;
  for(int iz = r_in_; iz < r_fi_ + 1; iz++){
    for(int iy = r_in_; iy < r_fi_ + 1; iy++){
      for(int ix = r_in_; ix < r_fi_ + 1; ix++){
	triple tmp(ascale*ix,ascale*iy,ascale*iz); 
	( (int*)& x[ixyz])[0]= ixyz;
	y[ixyz]  = datain[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	dy[ixyz] = datain_sigma[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	xyz_1g1yy3D[ixyz] =tmp;
	ixyz++;
      }
    }
  }
  
  
  // for boundary condition
  origins_1g1yy3D = new triple[Norigins_1g1yy3D];
  
  {
    double Lx = ascale*XnodeSites/2.0;
    double Ly = ascale*YnodeSites/2.0;
    double Lz = ascale*ZnodeSites/2.0;
    int count = 0;
    
    for(    double x = 0; x <1*Lx; x+=Lx){
      for(  double y = 0; y <1*Ly; y+=Ly){
	for(double z = 0; z <1*Lz; z+=Lz){
	  origins_1g1yy3D[count] = triple(x,y,z);
	  count++;
	}
      }
    }
  }
  
    
  
  int		ndata	   = fit_lengh*fit_lengh*fit_lengh;
  int	        ia[]	   = {1,1,1,1,1}; 
  int		ma	   = Nparam;
  int		nca	   = Nparam;
  double	covar[Nparam*Nparam];
  double	alpha[Nparam*Nparam];
  double	alambda;
  static double	a[Nparam];
  static int	is_initial = 1;
  
  
  if(is_initial==1){
        is_initial =0;

    //it07
	a[0] = 1084.59;
	a[1] = 54.265;
	a[2] = -342.863;
	a[3] = 2.35015;
	a[4] =  2.65749;
	/*
    a[0] = 1079.14;
    a[1] = 53.9753;
    a[2] = -361.082;
    a[3] = 2.3395;
    a[4] = 2.7287;
    */
  }


  // The initial call
  alambda = -1.0;
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy3D, alambda);

  // iterations
  for(int iter=0; iter<65536; iter++){
    mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy3D, alambda);
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
  mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, Func1g1yy3D, alambda);
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
