#ifndef __ErrorEstimate_dir_h_wanghan__
#define __ErrorEstimate_dir_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Defines.h"
#include "DensityProfile.h"
#include "Integral1D.h"

class Fs1k
{
public:
  double k;
  double beta;
  double operator () (const double & r) const
      {
	double betar = beta * r;
	double ri = 1./r;
	double fr = ( 2. * beta / sqrt(M_PI) * exp (-betar * betar) +
		      erfc(betar) * ri ) * ri;
	// K_{homo} = fr * fr
	return  r * fr * fr * sin(2.*M_PI*k*r);
      }
};


class Fs1k1
{
public:
  double k;
  double beta;
  double operator () (const double & r) const
      {
	double betar = beta * r;
	double ri = 1./r;
	double fr = ( 2. * beta / sqrt(M_PI) * exp (-betar * betar) +
		      erfc(betar) * ri ) * ri;
	return r * r * fr * fr;
      }
};

class Fs2k
{
public:
  double k;
  double beta;
  double operator () (const double & r) const
      {
	double betar = beta * r;
	return -2. * erfc (betar) * sin(2.*M_PI*k*r);
      }
};
  

class ErrorEstimate_dir
{
  std::vector<double > boxsize;
  IntVectorType K;
  int nele;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;
  double beta;
  double rcut;
  Fs1k  f_s1k;
  Fs1k1 f_s1k1;
  Fs2k  f_s2k;
  Integral1D<Fs1k,  double> inte_s1k;
  Integral1D<Fs1k1, double> inte_s1k1;
  Integral1D<Fs2k,  double> inte_s2k;
private:
  void freeAll();
  void calV();
  void calAStar();
  void makeS1k ();
  void makeS2k ();
  inline int  index3to1 (const IntVectorType i,
			 const IntVectorType N) const;
  inline void index1to3 (const unsigned input,
			 const IntVectorType N,
			 IntVectorType * i) const;
private:
  fftw_complex *rho1;
  fftw_complex *rho2;
  fftw_complex *s2kx, *s2ky, *s2kz;
  // fftw_complex *k1rx, *k1ry, *k1rz;
  fftw_complex *s1k;
  fftw_complex *error1x, *error1y, *error1z;
  fftw_complex *error1;
  fftw_complex *error2;
  fftw_plan p_forward_rho1, p_forward_rho2;
  // fftw_plan p_backward_s2kx;
  // fftw_plan p_backward_s2ky;
  // fftw_plan p_backward_s2kz;
  // fftw_plan p_forward_s1k;
  fftw_plan p_backward_error1x;
  fftw_plan p_backward_error1y;
  fftw_plan p_backward_error1z;
  fftw_plan p_backward_error2;
  bool malloced;
private:
  double integral_s1k_numerical (const double & k,
				 const double & r1,
				 const double & r2,
				 const double & prec);
  double integral_s1k1_numerical (const double & k,
				  const double & r1,
				  const double & r2,
				  const double & prec);
  double integral_s2k_numerical (const double & k,
				 const double & r1,
				 const double & r2,
				 const double & prec);
public:
  ErrorEstimate_dir();
  ~ErrorEstimate_dir();
  void reinit (const double & beta,
	       const double & rcut,
	       const DensityProfile_PiecewiseConst & dp);
  void calError (const DensityProfile_PiecewiseConst & dp);
  void print_error (const std::string & file) const;
  void print_meanf (const std::string & file,
		    const DensityProfile_PiecewiseConst & dp) const;
}
    ;

int ErrorEstimate_dir::
index3to1 (const IntVectorType i,
	   const IntVectorType N) const
{
  return i.z + N.z * (i.y + N.y * i.x);
}

void ErrorEstimate_dir::
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
