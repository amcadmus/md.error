#ifndef __DensityProfile_h_wanghan__
#define __DensityProfile_h_wanghan__

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <cmath>
#include <string>

#define CPLUSPLUS
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

class DensityProfile_PiecewiseConst
{
  std::vector<double > boxsize;
  unsigned nx, ny, nz, nele;
  double   hx, hy, hz;
  std::vector<double > profile;
  int natoms;
public:
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  void reinit_xtc (const std::string & filename,
		   const double & refh);
  void deposit (const std::string & filename,
		const float & start,
		const float & end);
  void reinit_conf (const std::string & filename,
		    const double & refh);
public:
  const double & getProfile (const unsigned & ix,
			     const unsigned & iy,
			     const unsigned & iz) const
      {return profile[index3to1(ix, iy, iz)];}
  const double & getProfile (const unsigned & i) const
      {return profile[i];}
  inline const double & getValue (const double & xx,
				  const double & yy,
				  const double & zz) const;
  const unsigned & getNx () const {return nx;}
  const unsigned & getNy () const {return ny;}
  const unsigned & getNz () const {return nz;}
  const std::vector<double > & getBox () const {return boxsize;}
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
  void save (const char * file) const;
  void load (const char * file);
};

    
unsigned DensityProfile_PiecewiseConst::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void DensityProfile_PiecewiseConst::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}

const double & DensityProfile_PiecewiseConst::
getValue (const double & xx,
	  const double & yy,
	  const double & zz) const
{
  unsigned ix = unsigned (xx / hx);
  unsigned iy = unsigned (yy / hy);
  unsigned iz = unsigned (zz / hz);
  return getProfile (ix, iy, iz);
}




#endif
