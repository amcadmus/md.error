#ifndef __ForceCorr_h_wanghan__
#define __ForceCorr_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Integral1D.h"
#include "DensityProfile.h"
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

class F13 
{
public:
  double k;
  double operator () (const double & x) const
      {
	return gsl_pow_int(x, -13) * sin(2.*M_PI*k*x);
      }
};

class F5
{
public:
  double k;
  double operator () (const double & x) const
      {
	return gsl_pow_int(x, -5) * sin(2.*M_PI*k*x);
      }
};


class ForceCorr 
{
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double hx, hy, hz;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  inline int  index3to1 (int  ix, int  iy, int  iz) const;
  inline void index1to3 (int& input,
			 int& ix, int& iy, int& iz) const;
private:
  bool malloced;
  double rc;

  cufftComplex *copyBuff;
  cufftComplex *d_rhor;
  cufftComplex *d_rhok;
  cufftComplex *d_s2kx, *d_s2ky, *d_s2kz;
  cufftComplex *d_error2kx, *d_error2ky, *d_error2kz;
  cufftHandle plan;
public:
  cufftComplex *d_error2rx, *d_error2ry, *d_error2rz;
private:
  // fftw_complex *rhor, *rhok;
  fftw_complex *s2kx, *s2ky, *s2kz;
  // fftw_complex *error2rx, *error2ry, *error2rz;
  // fftw_complex *error2kx, *error2ky, *error2kz;
  // fftw_plan p_forward_rho;
  // fftw_plan p_backward_error2x;
  // fftw_plan p_backward_error2y;
  // fftw_plan p_backward_error2z;
  F5  f_inte5;
  F13 f_inte13;
  Integral1D<F5,  double> inte5;
  Integral1D<F13, double> inte13;
private:
  void freeAll ();
  double integral_an13_numerical (const double & k,
				  const double & r1,
				  const double & r2,
				  const double & prec);
  double integral_an5_numerical (const double & k,
				 const double & r1,
				 const double & r2,
				 const double & prec);
  void makeS2k (const double & epsilon,
		const double & sigma);
public:
  ForceCorr ();
  ~ForceCorr ();
public:
  unsigned getNx () const {return nx;}
  unsigned getNy () const {return ny;}
  unsigned getNz () const {return nz;}
  std::vector<double > getBox () const {return boxsize;}
  // fftw_complex * getForceCorrX () {return error2rx;}
  // fftw_complex * getForceCorrY () {return error2ry;}
  // fftw_complex * getForceCorrZ () {return error2rz;}
  // const fftw_complex * getForceCorrX () const {return error2rx;}
  // const fftw_complex * getForceCorrY () const {return error2ry;}
  // const fftw_complex * getForceCorrZ () const {return error2rz;}
public:
  void reinit (const double & rc,
	       DensityProfile_PiecewiseConst & dp);
  void calError (const DensityProfile_PiecewiseConst & dp);
public:
  void print_x  (const std::string & file) const;
}
    ;



unsigned ForceCorr::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ForceCorr::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}



int ForceCorr::
index3to1 (int  ix, int  iy, int  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ForceCorr::
index1to3 (int& input,
	   int& ix, int& iy, int& iz) const
{
  int tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}


#endif 
