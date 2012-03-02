#ifndef __ErrorEstimate_SPME_Ana_SelfCorr_h_wanghan__
#define __ErrorEstimate_SPME_Ana_SelfCorr_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"

class ErrorEstimate_SPME_Ana_SelfCorr
{
  IntVectorType K;
  // IntVectorType Kmax;
  int kcut;
  int numk;
  double beta;
  int order;
  int nele;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;
private:
  void freeAll();
  void calV();
  void calAStar();
  void calKernel();
  inline int  index3to1 (const IntVectorType i,
			 const IntVectorType N) const;
  inline void index1to3 (const unsigned input,
			 const IntVectorType N,
			 IntVectorType * i) const;
  fftw_complex ** myMalloc (int numOfk, int nele);
  void myFree (int numOfk, fftw_complex ** a);
private:
  fftw_complex *rho1;
  fftw_complex *rho2;
  fftw_complex *k1mx, *k1my, *k1mz;
  fftw_complex *k1rx, *k1ry, *k1rz;
  fftw_complex *k2ma, *k2mb;
  fftw_complex *error1x, *error1y, *error1z;
  fftw_complex *error1;
  fftw_complex *error2;
  double self_error;
  double qself_error;

  fftw_plan p_forward_rho1, p_forward_rho2;
  fftw_plan p_backward_k1mx;
  fftw_plan p_backward_k1my;
  fftw_plan p_backward_k1mz;
  fftw_plan p_forward_k2a;
  fftw_plan p_forward_k2b;
  fftw_plan p_backward_error1x;
  fftw_plan p_backward_error1y;
  fftw_plan p_backward_error1z;
  fftw_plan p_backward_error2;

  fftw_complex **Fmx; // k
  fftw_complex **Fmy;
  fftw_complex **Fmz;
  fftw_complex **FConvRho1x;
  fftw_complex **FConvRho1y;
  fftw_complex **FConvRho1z;
  fftw_complex **Frx;
  fftw_complex **Fry;
  fftw_complex **Frz;
  fftw_complex *selfx;
  fftw_complex *selfy;
  fftw_complex *selfz;
  fftw_plan *p_backward_FConvRho1x; // k
  fftw_plan *p_backward_FConvRho1y;
  fftw_plan *p_backward_FConvRho1z;
  fftw_plan *p_backward_Fmx; // k
  fftw_plan *p_backward_Fmy;
  fftw_plan *p_backward_Fmz;
  
  bool malloced;
public:
  ErrorEstimate_SPME_Ana_SelfCorr();
  ~ErrorEstimate_SPME_Ana_SelfCorr();
  void reinit (const double & beta,
	       const int & order,
	       // const IntVectorType Kmax,
	       const DensityProfile_PiecewiseConst & dp);
  void calError (const DensityProfile_PiecewiseConst & dp,
		 const double charge = 1.);
  void print_error (const std::string & file) const;
  void print_meanf (const std::string & file,
		    const DensityProfile_PiecewiseConst & dp) const;
}
    ;


int ErrorEstimate_SPME_Ana_SelfCorr::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_SPME_Ana_SelfCorr::
index1to3 (const unsigned input,
	   const IntVectorType N,
	   IntVectorType * i) const
{
  unsigned tmp = input;
  i->z = tmp % (N.z);
  tmp = (tmp - i->z) / N.z;
  i->y = tmp % (N.y);
  i->x =  (tmp - i->y) / N.y;
}

#endif
