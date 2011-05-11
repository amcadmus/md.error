#ifndef __Test_h_wanghan__
#define __Test_h_wanghan__

#include "Integral1D.h"
#include "Integral3D.h"
#include <gsl/gsl_sf_bessel.h>

class tmpf3 
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class tmpf
{
public:
  double operator () (const double & r) const;
};


#endif

