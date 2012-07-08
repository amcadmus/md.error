#ifndef __ErrorEstimate_SPME_St_H2O_gr_h_wanghan__
#define __ErrorEstimate_SPME_St_H2O_gr_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"
#include "Poly.h"
#include <string>

using namespace std;
// const double global_integral_zero_epsilon = 1e-3;

class ErrorEstimate_SPME_St_H2O_Function_gr 
{
  const PiecewisePoly * ps;
  double mm;
public:
  ErrorEstimate_SPME_St_H2O_Function_gr (const PiecewisePoly & ps_,
					 const double & mm_)
      : ps (&ps_), mm(mm_) {}
  void setm (const double & mm_);
  inline double operator () (const double & r) const;
}
    ;

inline void ErrorEstimate_SPME_St_H2O_Function_gr::
setm (const double & mm_)
{
  mm = mm_;
}

inline double ErrorEstimate_SPME_St_H2O_Function_gr::
operator () (const double & r) const
{
  return 2. * r / mm * (ps->value(r) - 1.) * sin(2. * M_PI * mm * r);
}

class ErrorEstimate_SPME_St_H2O_gr
{
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
  void calTm (const std::vector<double > & grx,
	      const std::vector<double > & grv,
	      fftw_complex * out);
  void calTm (const std::vector<std::vector<double > > & dirs,
	      const std::vector<double > & charges,
	      fftw_complex * out);
  void calTOTH (const std::vector<std::vector<double > > & moles,
		const std::vector<double > & charges,
		const double & rhoMol,
		const string & fileOO,
		const string & fileOH,
		const string & fileHH);
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
  fftw_complex *k2ma;
  fftw_complex *error1x, *error1y, *error1z;
  fftw_complex *error1;
  fftw_complex *error2;
  fftw_complex *TO;
  fftw_complex *TH;
  fftw_complex *TOKx, *TOKy, *TOKz;
  fftw_complex *THKx, *THKy, *THKz;
  fftw_complex *coKO, *coKH;
  fftw_complex *rhoO, *rhoH;
  fftw_complex *coerrorO, *coerrorH;
  
  fftw_plan p_forward_rho1, p_forward_rho2;
  fftw_plan p_backward_k1mx;
  fftw_plan p_backward_k1my;
  fftw_plan p_backward_k1mz;
  fftw_plan p_forward_k2a;
  fftw_plan p_backward_error1x;
  fftw_plan p_backward_error1y;
  fftw_plan p_backward_error1z;
  fftw_plan p_backward_error2;
  
  fftw_plan p_backward_TOKx;
  fftw_plan p_backward_TOKy;
  fftw_plan p_backward_TOKz;
  fftw_plan p_backward_THKx;
  fftw_plan p_backward_THKy;
  fftw_plan p_backward_THKz;
  fftw_plan p_forward_coKO;
  fftw_plan p_forward_coKH;
  fftw_plan p_backward_coerrorO;
  fftw_plan p_backward_coerrorH;
  fftw_plan p_forward_rhoO;
  fftw_plan p_forward_rhoH;
  
  bool malloced;
public:
  ErrorEstimate_SPME_St_H2O_gr();
  ~ErrorEstimate_SPME_St_H2O_gr();
  void reinit (const double & beta,
	       const int & order,
	       // const IntVectorType Kmax,
	       const DensityProfile_PiecewiseConst & dp,
	       const std::vector<std::vector<double > > & moles,
	       const std::vector<double > & charges,
	       const double & rhoMol,
	       const string & fileOO,
	       const string & fileOH,
	       const string & fileHH);
  void calError (const DensityProfile_PiecewiseConst & dp,
		 const double charge = 1.);
  void print_error (const std::string & file) const;
  void print_meanf (const std::string & file,
		    const DensityProfile_PiecewiseConst & dp) const;
}
    ;


int ErrorEstimate_SPME_St_H2O_gr::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_SPME_St_H2O_gr::
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
