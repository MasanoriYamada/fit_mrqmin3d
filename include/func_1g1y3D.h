//sub class
#ifndef FIT_FUNC_1g1y3D_20131121_H_
#define FIT_FUNC_1g1y3D_20131121_H_

#include "fit.h"
#include <string>

class Fit1g1y3D : public FitInterface{
public:

  Fit1g1y3D(){
    r_fi_     = 0;
    r_in_     = 0;
    fit_lengh = 0;
    fitname_  =	"1g1y3D";
  }
  ~Fit1g1y3D(){}


  std::string getName()
{
  return fitname_;
}
  int getDofA()
{
  return Nparam;
}
protected:
  void fit(double* datain,double* datain_sigma, double* a, double& chisq);
  void setRange(int in, int fi)
{
  r_in_=in;
  r_fi_=fi;
  fit_lengh = r_fi_ - r_in_ + 1;
}

private:
  static const int     Nparam	  = 5;
  int r_fi_;
  int r_in_;
  int fit_lengh;
  std::string fitname_;
  
};

#endif
