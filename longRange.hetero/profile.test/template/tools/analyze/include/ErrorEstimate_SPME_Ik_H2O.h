#ifndef __ErrorEstimate_SPME_Ik_H2O_h_wanghan__
#define __ErrorEstimate_SPME_Ik_H2O_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"

class ErrorEstimate_SPME_Ik_H2O
{
  friend class ErrorEstimate_SPME_Ana_old;
  
  IntVectorType K;
  // IntVectorType Kmax;
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
  void calTm (const std::vector<std::vector<double > > & dirs,
	      const std::vector<double > & charges,
	      fftw_complex * out);
  void calTOTH (const std::vector<std::vector<double > > & moles,
		const std::vector<double > & charges);
  inline int  index3to1 (const IntVectorType i,
			 const IntVectorType N) const;
  inline void index1to3 (const unsigned input,
			 const IntVectorType N,
			 IntVectorType * i) const;
private:
  fftw_complex *rho1;
  fftw_complex *rho2;
  fftw_complex *rhoH;
  fftw_complex *rhoO;
  fftw_complex *k1mx, *k1my, *k1mz;
  fftw_complex *k1rx, *k1ry, *k1rz;
  fftw_complex *coO_k1mx, *coO_k1my, *coO_k1mz;
  fftw_complex *coO_k1rx, *coO_k1ry, *coO_k1rz;
  fftw_complex *coH_k1mx, *coH_k1my, *coH_k1mz;
  fftw_complex *coH_k1rx, *coH_k1ry, *coH_k1rz;
  fftw_complex *k2ma, *k2mb;
  fftw_complex *coO_k2ma, *coO_k2mb;
  fftw_complex *coH_k2ma, *coH_k2mb;
  fftw_complex *error1x, *error1y, *error1z;
  fftw_complex *error1;
  fftw_complex *error2;
  fftw_complex *coO_error2;
  fftw_complex *coH_error2;
  fftw_complex *TO;
  fftw_complex *TH;

  fftw_plan p_forward_rho1, p_forward_rho2;
  fftw_plan p_forward_rhoO, p_forward_rhoH;
  fftw_plan p_backward_k1mx;
  fftw_plan p_backward_k1my;
  fftw_plan p_backward_k1mz;
  fftw_plan p_backward_coO_k1mx;
  fftw_plan p_backward_coO_k1my;
  fftw_plan p_backward_coO_k1mz;
  fftw_plan p_backward_coH_k1mx;
  fftw_plan p_backward_coH_k1my;
  fftw_plan p_backward_coH_k1mz;
  fftw_plan p_forward_k2a;
  fftw_plan p_forward_k2b;
  fftw_plan p_forward_coO_k2a;
  fftw_plan p_forward_coO_k2b;
  fftw_plan p_forward_coH_k2a;
  fftw_plan p_forward_coH_k2b;
  fftw_plan p_backward_error1x;
  fftw_plan p_backward_error1y;
  fftw_plan p_backward_error1z;
  fftw_plan p_backward_error2;
  fftw_plan p_backward_coO_error2;
  fftw_plan p_backward_coH_error2;
  bool malloced;
public:
  ErrorEstimate_SPME_Ik_H2O();
  ~ErrorEstimate_SPME_Ik_H2O();
  void reinit (const double & beta,
	       const int & order,
	       // const IntVectorType Kmax,
	       const DensityProfile_PiecewiseConst & dp,
	       const std::vector<std::vector<double > > & moles,
	       const std::vector<double > & charges);
  void calError (const DensityProfile_PiecewiseConst & dp,
		 const double charge = 1.);
  void print_error (const std::string & file) const;
  void print_meanf (const std::string & file,
		    const DensityProfile_PiecewiseConst & dp) const;
}
    ;


int ErrorEstimate_SPME_Ik_H2O::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_SPME_Ik_H2O::
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
