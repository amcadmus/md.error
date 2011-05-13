#ifndef __ErrorEstimator_h_wanghan__
#define __ErrorEstimator_h_wanghan__

#include "DensityProfile.h"
#include "ForceKernel.h"
#include <vector>
#include <cmath>

class ErrorEstimatorConvN2 
{
  const ForceKernel & fk;
  std::vector<double > boxsize;
  unsigned nx, ny, nz;
  double   hx, hy, hz;
  std::vector<double > profile;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  ErrorEstimatorConvN2 (const ForceKernel & fk,
			const std::vector<double > & box,
			const double refh);
  void estimate (const DensityProfile_PiecewiseConst & dp);
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
};

unsigned ErrorEstimatorConvN2::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ErrorEstimatorConvN2::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}




//
// consider correlation
//
class ErrorEstimatorConvN2_Corr
{
  const ForceKernel & fk;
  std::vector<double > boxsize;
  int  nx, ny, nz;
  double   hx, hy, hz;
  std::vector<double > profile;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  ErrorEstimatorConvN2_Corr (const ForceKernel & fk,
			const std::vector<double > & box,
			const double refh);
  void estimate (const DensityProfile_Corr_PiecewiseConst & dp,
		 const int & corrBond);
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
};

unsigned ErrorEstimatorConvN2_Corr::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ErrorEstimatorConvN2_Corr::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}



//
// with FFT
//

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
// #include "specialfunctions.h"
#include "Integral1D.h"

class f13 
{
public:
  double k;
  double operator () (const double & x) const
      {
	double x2 = x*x;
	double x4 = x2*x2;
	return 1./(x4*x4*x4*x) * sin(2.*M_PI*k*x);
      }
};

class f5
{
public:
  double k;
  double operator () (const double & x) const
      {
	double x2 = x*x;
	return 1./(x2*x2*x) * sin(2.*M_PI*k*x);
      }
};

class f15
{
public:
  double k;
  double operator () (const double & x) const
      {
	return gsl_pow_int(x, -15) * sin(2.*M_PI*k*x);
      }
};


class F_I3 
{
public:
  double k;
  double operator () (const double & x) const
      {
	return (-4032) / (gsl_pow_int(x, 13) * sqrt(x)) *
	    gsl_sf_bessel_Jnu (1.5, 2 * M_PI * k * x);
      }
};

class FF_I1
{
public:
  double k;
  double operator () (const double & x) const
      {
	return (64. * 24. * 24.) / (gsl_pow_int(x, 14) * sqrt(x)) *
	    gsl_sf_bessel_Jnu (0.5, 2 * M_PI * k * x);
      }
};

class FF_I5
{
public:
  double k;
  double operator () (const double & x) const
      {
	return (64. * 24. * 24.) / (gsl_pow_int(x, 14) * sqrt(x)) *
	    gsl_sf_bessel_Jnu (2.5, 2 * M_PI * k * x);
      }
};


class ErrorEstimatorFFT_Corr
{
  const ForceKernel & fk;
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double   hx, hy, hz;
  int corrBond, corrDim;
  int numCorr;
  double corrCut;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  fftw_complex *f2k;
  fftw_complex *s1r, *s1k;
  fftw_complex *s2rx, *s2ry, *s2rz, *s2r;
  fftw_complex *s2kx, *s2ky, *s2kz;
  fftw_complex *s3k, *s3r;
  fftw_complex *rhor, *rhok;
  fftw_complex *corrr, *corrk;
  fftw_complex *corr0r, *corr1r, *corr2r;
  fftw_complex *corr0k, *corr1k, *corr2k;
  fftw_complex *corr00r, *corr01r, *corr02r, *corr11r, *corr12r, *corr22r;
  fftw_complex *corr00k, *corr01k, *corr02k, *corr11k, *corr12k, *corr22k;
  fftw_complex *b0k, *b1k, *b2k;
  fftw_complex *a00k, *a01k, *a02k, *a11k, *a12k, *a22k;
  // fftw_complex *a00r, *a01r, *a02r, *a11r, *a12r, *a22r;
  fftw_plan p_forward_rho, p_backward_rho;
  fftw_plan p_forward_corr;
  fftw_plan p_forward_corr0, p_forward_corr1, p_forward_corr2;
  fftw_plan p_forward_corr00, p_forward_corr01, p_forward_corr02;
  fftw_plan p_forward_corr11, p_forward_corr12, p_forward_corr22;
  fftw_plan p_backward_s1;
  fftw_plan p_backward_s2x;
  fftw_plan p_backward_s2y;
  fftw_plan p_backward_s2z;
  fftw_plan p_backward_s3;
  
  // fftw_plan p_backward_a00, p_backward_a01, p_backward_a02;
  // fftw_plan p_backward_a11, p_backward_a12, p_backward_a22;
private:
  f5  f_inte5;
  f13 f_inte13;
  f15 f_inte15;
  Integral1D<f5,  double> inte5;
  Integral1D<f13, double> inte13;
  Integral1D<f15, double> inte15;
  F_I3 f_i3;
  FF_I1 ff_i1;
  FF_I5 ff_i5;
  Integral1D<F_I3, double > inte_f_i3;
  Integral1D<FF_I1, double > inte_ff_i1;
  Integral1D<FF_I5, double > inte_ff_i5;
  double bond_f_i3;
  double prec_f_i3;
  double bond_ff_i1;
  double prec_ff_i1;
  double bond_ff_i5;
  double prec_ff_i5;
private:
  void formCorrMat (const DensityProfile_Corr_PiecewiseConst & dp);
  void makef2k (const double & sigma,
		const double & rc) ;
  void makeS2k (const double & epsilon,
		const double & sigma,
		const double & rc);
  void makeS3k (const double & epsilon,
		const double & sigma,
		const double & rc);
public:
  ErrorEstimatorFFT_Corr (const ForceKernel & fk,
			  const DensityProfile_Corr_PiecewiseConst & dp);
  ~ErrorEstimatorFFT_Corr ();
  void estimate (const DensityProfile_Corr_PiecewiseConst & dp,
		 const int & corrBond,
		 const double & rc);
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
public:
  double integral_an (const int & n_,
		      const double & k,
		      const double & rc);
  double integral_an1 (const int & n_,
		       const double & k,
		       const double & rc);
  double integral_an5_numerical  (const double & k,
				  const double & rc);
  double integral_an13_numerical (const double & k,
				  const double & rc);
  double integral_an15_numerical (const double & k,
				  const double & rc);
  double integral_f_i3_numerical (const double & k,
				  const double & rc);
  double integral_ff_i1_numerical (const double & k,
				   const double & rc);
  double integral_ff_i5_numerical (const double & k,
				   const double & rc);
};

unsigned ErrorEstimatorFFT_Corr::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ErrorEstimatorFFT_Corr::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}



#endif
