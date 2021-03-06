#ifndef __PressureCorrection_h_wanghan__
#define __PressureCorrection_h_wanghan__

#include "DensityProfile.h"
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

class FF_I1
{
public:
  double k;
  double * bessel_table;
  int bessel_table_length;
  double bessel_table_h;
  double bessel_table_hi;
  double operator () (const double & x) const
      {
	double xx = 2 * M_PI * k * x;
	int posi = xx * bessel_table_hi;
	return (-24.) / (gsl_pow_4(x) * sqrt(x)) *
	    ((xx * bessel_table_hi - posi) *
	     (bessel_table[posi+1] - bessel_table[posi]) +
	     bessel_table[posi]);
	// return (-24.) / (gsl_pow_4(x) * sqrt(x)) *
	//     gsl_sf_bessel_Jnu (0.5, 2 * M_PI * k * x);
      }
};

class FF_I5
{
public:
  double k;
  double * bessel_table;
  int bessel_table_length;
  double bessel_table_h;
  double bessel_table_hi;
  double operator () (const double & x) const
      {
	double xx = 2 * M_PI * k * x;
	int posi = xx * bessel_table_hi;
	return (-24.) / (gsl_pow_4(x) * sqrt(x)) *
	    ((xx * bessel_table_hi - posi) *
	     (bessel_table[posi+1] - bessel_table[posi]) +
	     bessel_table[posi]);
	// return (-24.) / (gsl_pow_4(x) * sqrt(x)) *
	//     gsl_sf_bessel_Jnu (2.5, 2 * M_PI * k * x);
      }
};

class PresureCorrection
{
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double hx, hy, hz;
  double rcut;
  inline int  index3to1 (int  ix, int  iy, int  iz) const;
  inline void index1to3 (int& input,
			 int& ix, int& iy, int& iz) const;
  fftw_complex *rhor, *rhok;
  fftw_complex *kxxk, *kyyk, *kzzk;
  fftw_complex *convxx, *convyy, *convzz;
  fftw_plan p_forward_rho;
  fftw_plan p_backward_kxx;
  fftw_plan p_backward_kyy;
  fftw_plan p_backward_kzz;
  bool malloced;
private:
  double * bessel_table_1;
  double * bessel_table_5;
  int bessel_table_length;
  double integral_upper;
  double bessel_table_h, bessel_table_hi;
  FF_I1 ff_i1;
  FF_I5 ff_i5;
private:
  void freeAll ();
  void makeKernel (const double & rc);
  double integral_ff_i1_numerical (const double & k,
				   const double & rc);
  double integral_ff_i5_numerical (const double & k,
				   const double & rc);
public:
  double pxx, pyy, pzz;
  PresureCorrection (const double & rc,
		     const DensityProfile_PiecewiseConst & dp);
  ~PresureCorrection ();
  void reinit (const double & rc,
	       const DensityProfile_PiecewiseConst & dp);
  void naiveCorrection (const DensityProfile_PiecewiseConst & dp);
  void correction (const DensityProfile_PiecewiseConst & dp);
}
    ;


#include "RCutTable.h"

class PresureCorrection_RCutTable
{
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double hx, hy, hz;
  inline int  index3to1 (int  ix, int  iy, int  iz) const;
  inline void index1to3 (int& input,
			 int& ix, int& iy, int& iz) const;
  fftw_complex *rhor, *rhok;
  fftw_complex **kxxk, **kyyk, **kzzk;
  fftw_complex **convxx, **convyy, **convzz;
  fftw_plan p_forward_rho;
  fftw_plan *p_backward_kxx;
  fftw_plan *p_backward_kyy;
  fftw_plan *p_backward_kzz;
  bool malloced;
private:
  std::vector<double > rcList;
  int nrc;
private:
  double * bessel_table_1;
  double * bessel_table_5;
  int bessel_table_length;
  double integral_upper;
  double bessel_table_h, bessel_table_hi;
  FF_I1 ff_i1;
  FF_I5 ff_i5;
private:
  void freeAll ();
  void makeKernel ();
  double integral_ff_i1_numerical (const double & k,
				   const double & rc);
  double integral_ff_i5_numerical (const double & k,
				   const double & rc);
public:
  double pxx, pyy, pzz;
  PresureCorrection_RCutTable (const RCutTable & arc,
			       const DensityProfile_PiecewiseConst & dp);
  ~PresureCorrection_RCutTable ();
  void reinit (const RCutTable & arc,
	       const DensityProfile_PiecewiseConst & dp);
  void correction (const RCutTable & arc,
		   const DensityProfile_PiecewiseConst & dp);
}
    ;


int PresureCorrection::
index3to1 (int  ix, int  iy, int  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void PresureCorrection::
index1to3 (int& input,
	   int& ix, int& iy, int& iz) const
{
  int tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}

int PresureCorrection_RCutTable::
index3to1 (int  ix, int  iy, int  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void PresureCorrection_RCutTable::
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
