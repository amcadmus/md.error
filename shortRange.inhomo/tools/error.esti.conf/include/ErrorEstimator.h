#ifndef __ErrorEstimator_h_wanghan__
#define __ErrorEstimator_h_wanghan__

#include "DensityProfile.h"
#include "ForceKernel.h"
#include <vector>

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
  ErrorEstimatorConvN2 (const ForceKernel & fk);
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


#endif
