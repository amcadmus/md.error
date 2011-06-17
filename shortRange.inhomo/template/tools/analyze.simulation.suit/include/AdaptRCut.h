#ifndef __AdaptRCut_h_wanghan__
#define __AdaptRCut_h_wanghan__

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "Integral1D.h"
#include "DensityProfile.h"
#include <vector>

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


class AdaptRCut 
{
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double hx, hy, hz;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
private:
  std::vector<double > rcList;
  int nrc;
  bool malloced;
  fftw_complex *rhor, *rhok;
  fftw_complex **s1k;
  fftw_complex **s2kx, **s2ky, **s2kz;
  fftw_complex **error1r, **error1k;
  fftw_complex **error2rx, **error2ry, **error2rz;
  fftw_complex **error2kx, **error2ky, **error2kz;
  fftw_complex **error;
  double *rcut;
  int *rcutIndex;
  double *result_error;
  fftw_plan p_forward_rho;
  fftw_plan *p_backward_error1;
  fftw_plan *p_backward_error2x;
  fftw_plan *p_backward_error2y;
  fftw_plan *p_backward_error2z;
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
  void makeS1k (const double & epsilon,
		const double & sigma);
  void makeS2k (const double & epsilon,
		const double & sigma);
  void calRCutOnePoint (const double & prec,
			const unsigned & idx);
public:
  AdaptRCut ();
  ~AdaptRCut ();
public:
  unsigned getNx () const {return nx;}
  unsigned getNy () const {return ny;}
  unsigned getNz () const {return nz;}
  std::vector<double > getBox () const {return boxsize;}
  const double * getRCut () const {return rcut;}
  const std::vector<double > & getRcList () const {return rcList;}
  const int * getProfileIndex () const {return rcutIndex;}
public:
  void reinit (const double & rcmin,
	       const double & rcmax,
	       const double & rcstep,
	       DensityProfile_PiecewiseConst & dp);
  void calError (const DensityProfile_PiecewiseConst & dp);
  void calRCut  (const double & prec);
  void uniformRCut (const double & rcut);
  void print_rc (const std::string & file) const;
  void print_error  (const std::string & file) const;
  void save_rc  (const std::string & file) const;
  void write_rc (const std::string & file) const;
}
    ;



unsigned AdaptRCut::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void AdaptRCut::
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
