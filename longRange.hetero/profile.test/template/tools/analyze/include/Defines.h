#ifndef __Defines_h_wanghan__
#define __Defines_h_wanghan__

#include <cmath>

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



#endif
