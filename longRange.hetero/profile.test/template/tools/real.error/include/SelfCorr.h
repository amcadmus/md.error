#ifndef __SelfCorr_wanghan_h__
#define __SelfCorr_wanghan_h__

#include <gsl/gsl_math.h>
#include <cmath>
#include <fftw3.h>
#include <vector>

struct IntVectorType
{
  int x, y, z;
};

struct VectorType
{
  double x, y, z;
};

struct MatrixType 
{
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
};

inline double 
kernel_rm1_rec_f (const double m2,
		  const double beta)
{
  double tmp = M_PI / beta;
  return exp (- tmp*tmp*m2) / m2;
}

class SelfCorr
{
  IntVectorType K;
  VectorType Ki;
  // IntVectorType Kmax;
  int kcut;
  int numk;
  double beta;
  int order;
  int nele;
  MatrixType vecA;
  MatrixType vecAStar;
  double volume;
  bool malloced;
private:
  void freeAll();
  void calV();
  void calAStar();
  void calKernel();
  // inline int  index3to1 (const IntVectorType i,
  // 			 const IntVectorType N) const;
  // inline void index1to3 (const unsigned input,
  // 			 const IntVectorType N,
  // 			 IntVectorType * i) const;
private:
  fftw_complex *selfx;
  fftw_complex *selfy;
  fftw_complex *selfz;
public:
  SelfCorr ();
  ~SelfCorr ();
  void reinit (const double & beta_,
	       const int & order_,
	       const int & nx,
	       const int & ny,
	       const int & nz,
	       const double & bx,
	       const double & by,
	       const double & bz);
  void correction (const std::vector<double > & posi,
		   std::vector<double > & force);
  void correction (const std::vector<std::vector<double > > & posi,
		   std::vector<std::vector<double > > & force);
}
    ;

// int SelfCorr::
// index3to1 (const IntVectorType i,
// 	   const IntVectorType N) const
// {
//   return i.z + N.z * (i.y + N.y * i.x);
// }

// void SelfCorr::
// index1to3 (const unsigned input,
// 	   const IntVectorType N,
// 	   IntVectorType * i) const
// {
//   unsigned tmp = input;
//   i->z = tmp % (N.z);
//   tmp = (tmp - i->z) / N.z;
//   i->y = tmp % (N.y);
//   i->x =  (tmp - i->y) / N.y;
// }

#endif
