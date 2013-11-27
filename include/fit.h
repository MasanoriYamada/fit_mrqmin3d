#ifndef FIT_H_
#define FIT_H_

#include <complex>
#include <string>
#include "../include/analys.h"


typedef std::complex<double> COMPLEX;



static const double ascale = 0.1209;
static const double hbarc = 197.327;
static const double ToMev = hbarc/ascale;

//set in out info


extern "C" {
  void mrqmin_(double x[],double y[],double sig[], int& ndata,double a[],int ia[],int& ma,
	       double covar[], double alpha[],int& nca, double &chisq,
	       void (*func)(double& x,double a[],double& yfit,double dyda[],int& ma),double& alambda);
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma);
}

//(Super class) Strategy pattern of Gof design
class FitInterface{
public:
  void fitDo(double* datain,double* datain_sigma, double* a, double& chisq)
  {
    fit(datain, datain_sigma, a, chisq);
  }
  void setRangeDo(int in, int fi)
  {
    setRange(in, fi);
  }
  virtual std::string getName() = 0;
  virtual int getDofA() = 0;

protected:
  virtual void fit(double* datain,double* datain_sigma, double* a, double& chisq) = 0;
  virtual void setRange(int in, int fi) = 0;
};




//Context class
class Fit{
public:

  Fit(FitInterface* fit_obj)
  {
    fit_obj_ = fit_obj;
  }

  ~ Fit()
  {
    free (fit_obj_);
  }

  void fit(double* datain,double* datain_sigma, double* a, double& chisq)
  {
    fit_obj_->fitDo(datain, datain_sigma, a, chisq);
  }
  void setRange(int in, int fi)
  {
    fit_obj_->setRangeDo(in, fi);
      }
  std::string getName()
  {
    return fit_obj_->getName();
  }
  int getDofA()
  {
    return  fit_obj_->getDofA();
  }

private:
  FitInterface* fit_obj_;

};

struct triple {
  double x;
  double y;
  double z;
  ~triple(){}
  triple(){}
  triple(double x0,double y0,double z0) : x(x0), y(y0), z(z0) {}

  double r() const { return sqrt(x*x + y*y + z*z); }
};

//#define m_po 0.597  // redused mass of PROTON_omega Mev
#define binPot(ixyz,iconf) binPot[ixyz + XYZnodeSites*iconf]

inline void reScale(COMPLEX* pot){
  for(int id = 0; id < XYZnodeSites; id++){
    pot[id] = ToMev*pot[id];
  }
}

#endif
