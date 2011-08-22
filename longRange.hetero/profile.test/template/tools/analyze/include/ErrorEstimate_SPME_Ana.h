#ifndef __ErrorEstimate_SPME_Ana_h_wanghan__
#define __ErrorEstimate_SPME_Ana_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"
#include "ErrorEstimate_SPME_Ik.h"

class ErrorEstimate_SPME_Ana
{
  IntVectorType refined_K;
  IntVectorType refine;
  int kcut;
  int numk;
  int numAlphak;

  IntVectorType K;
  int nele;
  double beta;
  int order;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;

  // ErrorEstimate_SPME_Ik spmeik;
  fftw_complex *refined_error1x, *refined_error1y, *refined_error1z;
  fftw_complex *refined_error1;

  fftw_complex *rho1;
  fftw_complex **k1mx, **k1my, **k1mz;
  fftw_complex **k1rx, **k1ry, **k1rz;
  fftw_complex **error1x, **error1y, **error1z;
  VectorType *selfx, *selfy, *selfz;
  fftw_plan p_forward_rho1;
  fftw_plan * p_backward_k1mx;
  fftw_plan * p_backward_k1my;
  fftw_plan * p_backward_k1mz;
  fftw_plan * p_backward_error1x;
  fftw_plan * p_backward_error1y;
  fftw_plan * p_backward_error1z;

  fftw_complex **  myMalloc (int numOfk, int nele);
  void myFree (int numOfk, fftw_complex ** a);
  bool malloced;
  
  ErrorEstimate_SPME_Ik spmeik;
private:
  void freeAll();
  void calV();
  void calAStar();
  void calKernel ();
  inline int  index3to1 (const IntVectorType i,
			 const IntVectorType N) const;
  inline void index1to3 (const unsigned input,
			 const IntVectorType N,
			 IntVectorType * i) const;
public:
  ErrorEstimate_SPME_Ana();
  ~ErrorEstimate_SPME_Ana();
  void reinit (const double & beta,
	       const int & order,
	       const DensityProfile_PiecewiseConst & dp,
      	       const IntVectorType refine);
  void calError (const DensityProfile_PiecewiseConst & dp);
  void print_meanf (const std::string & file,
		    const DensityProfile_PiecewiseConst & dp) const;
};


int ErrorEstimate_SPME_Ana::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_SPME_Ana::
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
