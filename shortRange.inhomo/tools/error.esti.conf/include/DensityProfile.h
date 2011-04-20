#ifndef __DensityProfile_h_wanghan__
#define __DensityProfile_h_wanghan__

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>


class DensityProfile_PiecewiseConst
{
  std::vector<double > boxsize;
  unsigned nx, ny, nz;
  double   hx, hy, hz;
  std::vector<double > profile;
public:
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  DensityProfile_PiecewiseConst (const std::string & filename,
				 const double & refh);
  DensityProfile_PiecewiseConst (const std::vector<std::string> & filename,
				 const double & refh);
  const double & getProfile (const unsigned & ix,
			     const unsigned & iy,
			     const unsigned & iz) const
      {return profile[index3to1(ix, iy, iz)];}
  const unsigned & getNx () const {return nx;}
  const unsigned & getNy () const {return ny;}
  const unsigned & getNz () const {return nz;}
  const std::vector<double > & getBox () const {return boxsize;}
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
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


#endif
