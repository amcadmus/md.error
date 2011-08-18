#ifndef __ErrorEstimate_Ewald_h_wanghan__
#define __ErrorEstimate_Ewald_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"

class ErrorEstimate_Ewald
{
  IntVectorType K;
  IntVectorType Kmax;
  int nele;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;
  double beta;
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
private:
  fftw_complex *rho1;
  fftw_complex *rho2;
  fftw_complex *k1mx, *k1my, *k1mz;
  fftw_complex *k1rx, *k1ry, *k1rz;
  fftw_complex *k2m;
  fftw_complex *error1x, *error1y, *error1z;
  fftw_complex *error1;
  fftw_complex *error2;
  fftw_plan p_forward_rho1, p_forward_rho2;
  fftw_plan p_backward_k1mx;
  fftw_plan p_backward_k1my;
  fftw_plan p_backward_k1mz;
  fftw_plan p_forward_k2;
  fftw_plan p_backward_error1x;
  fftw_plan p_backward_error1y;
  fftw_plan p_backward_error1z;
  fftw_plan p_backward_error2;
  bool malloced;
public:
  ErrorEstimate_Ewald();
  ~ErrorEstimate_Ewald();
  void reinit (const double & beta,
	       const IntVectorType & Kmax,
	       const DensityProfile_PiecewiseConst & dp);
  void calError (const DensityProfile_PiecewiseConst & dp);
  void print_error (const std::string & file) const ;
}
    ;

int ErrorEstimate_Ewald::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_Ewald::
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
