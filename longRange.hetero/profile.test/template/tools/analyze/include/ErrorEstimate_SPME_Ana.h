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
  inline unsigned kl2to1 (const int & k,
			  const int & l) const;
  inline void     kl1to2 (const unsigned & index,
			  int & k,
			  int & l);
  IntVectorType K;
  int nele;
  double beta;
  int order;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;

  // ErrorEstimate_SPME_Ik spmeik;
  fftw_complex *refined_error1x, *refined_error1y, *refined_error1z;
  fftw_complex *refined_termIx,  *refined_termIy,  *refined_termIz;
  fftw_complex *refined_termIIx, *refined_termIIy, *refined_termIIz;
  fftw_complex *refined_error1;
  fftw_complex *refined_error2;

  fftw_complex *rho1;
  fftw_complex *rho2;
  
  fftw_complex **Fmx; // k
  fftw_complex **Fmy;
  fftw_complex **Fmz;
  fftw_complex **FConvRho1x;
  fftw_complex **FConvRho1y;
  fftw_complex **FConvRho1z;
  fftw_complex **Frx;
  fftw_complex **Fry;
  fftw_complex **Frz;

  fftw_complex **FxFx; // k*l
  fftw_complex **FyFy;
  fftw_complex **FzFz;
  fftw_complex **FxFxConvRho2; // k*l
  fftw_complex **FyFyConvRho2;
  fftw_complex **FzFzConvRho2;

  fftw_complex **Gmxx; // k
  fftw_complex **Gmxy;
  fftw_complex **Gmxz;
  fftw_complex **Gmyx; // k
  fftw_complex **Gmyy;
  fftw_complex **Gmyz;
  fftw_complex **Gmzx; // k
  fftw_complex **Gmzy;
  fftw_complex **Gmzz;

  fftw_complex **GxxFx; // k*l
  fftw_complex **GxyFy;
  fftw_complex **GxzFz;
  fftw_complex **GyxFx;
  fftw_complex **GyyFy;
  fftw_complex **GyzFz;
  fftw_complex **GzxFx;
  fftw_complex **GzyFy;
  fftw_complex **GzzFz;
  fftw_complex **GxxFxConvRho2; // k*l
  fftw_complex **GxyFyConvRho2;
  fftw_complex **GxzFzConvRho2;
  fftw_complex **GyxFxConvRho2;
  fftw_complex **GyyFyConvRho2;
  fftw_complex **GyzFzConvRho2;
  fftw_complex **GzxFxConvRho2;
  fftw_complex **GzyFyConvRho2;
  fftw_complex **GzzFzConvRho2;
  
  fftw_complex *selfx;
  fftw_complex *selfy;
  fftw_complex *selfz;
  
  fftw_plan p_forward_rho1;
  fftw_plan p_forward_rho2;
  
  // fftw_plan *p_backward_F;
  fftw_plan *p_backward_FConvRho1x; // k
  fftw_plan *p_backward_FConvRho1y;
  fftw_plan *p_backward_FConvRho1z;

  fftw_plan *p_backward_Fmx; // k
  fftw_plan *p_backward_Fmy;
  fftw_plan *p_backward_Fmz;
  fftw_plan *p_backward_Gmxx;
  fftw_plan *p_backward_Gmxy;
  fftw_plan *p_backward_Gmxz;
  fftw_plan *p_backward_Gmyx;
  fftw_plan *p_backward_Gmyy;
  fftw_plan *p_backward_Gmyz;
  fftw_plan *p_backward_Gmzx;
  fftw_plan *p_backward_Gmzy;
  fftw_plan *p_backward_Gmzz;
  
  fftw_plan *p_forward_FxFx; // k*l
  fftw_plan *p_forward_FyFy;
  fftw_plan *p_forward_FzFz;
  fftw_plan *p_forward_GxxFx;
  fftw_plan *p_forward_GxyFy;
  fftw_plan *p_forward_GxzFz;
  fftw_plan *p_forward_GyxFx;
  fftw_plan *p_forward_GyyFy;
  fftw_plan *p_forward_GyzFz;
  fftw_plan *p_forward_GzxFx;
  fftw_plan *p_forward_GzyFy;
  fftw_plan *p_forward_GzzFz;

  fftw_plan *p_backward_FxFxConvRho2; // k*l
  fftw_plan *p_backward_FyFyConvRho2;
  fftw_plan *p_backward_FzFzConvRho2;  
  fftw_plan *p_backward_GxxFxConvRho2;
  fftw_plan *p_backward_GxyFyConvRho2;
  fftw_plan *p_backward_GxzFzConvRho2;
  fftw_plan *p_backward_GyxFxConvRho2;
  fftw_plan *p_backward_GyyFyConvRho2;
  fftw_plan *p_backward_GyzFzConvRho2;
  fftw_plan *p_backward_GzxFxConvRho2;
  fftw_plan *p_backward_GzyFyConvRho2;
  fftw_plan *p_backward_GzzFzConvRho2;

  bool malloced;
  
  ErrorEstimate_SPME_Ik spmeik;
private:
  void freeAll();
  void calKernel ();
  fftw_complex ** myMalloc (int numOfk, int nele);
  void myFree (int numOfk, fftw_complex ** a);

  inline int  index3to1 (const IntVectorType i,
			 const IntVectorType N) const;
  inline void index1to3 (const unsigned input,
			 const IntVectorType N,
			 IntVectorType * i) const;

  void interpolate (const IntVectorType ii,		// on the refined grid
		    const fftw_complex * value,		// on the coarse grid
		    fftw_complex & result);
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
  void print_error (const std::string & file) const;
};

unsigned ErrorEstimate_SPME_Ana::
kl2to1 (const int & k,
	const int & l) const
{
  return l + numk * k;
}

void ErrorEstimate_SPME_Ana::
kl1to2 (const unsigned & index,
	int & k,
	int & l)
{
  unsigned tmp = index;
  l = tmp % numk;
  k = (tmp - l) / numk;
}

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
