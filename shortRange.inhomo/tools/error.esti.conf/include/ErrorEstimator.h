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
#include "specialfunctions.h"
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

class ErrorEstimatorFFT_Corr
{
  const ForceKernel & fk;
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double   hx, hy, hz;
  int corrBond, corrDim;
  int numCorr;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  fftw_complex * f2k;
  fftw_complex * s1r, * s1k;
  fftw_complex * rhor, * rhok;
  fftw_plan p_forward_rho, p_backward_rho;
  fftw_plan p_forward_s1,  p_backward_s1;
private:
  f13 f_inte;
  Integral1D<f13, double> inte;
private:
  void makef2k (const double & sigma,
		const double & rc) ;
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
  double integral_an13_numerical (const double & k,
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
