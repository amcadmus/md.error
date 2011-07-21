#ifndef __DensityProfile_h_wanghan__
#define __DensityProfile_h_wanghan__

#include <vector>
#include <cmath>


class ErrorProfile_PiecewiseConst
{
  std::vector<double > boxsize;
  unsigned nx, ny, nz;
  double   hx, hy, hz;
  std::vector<double > profile;
  std::vector<int    > count;
  unsigned nfile_record;
public:
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  ErrorProfile_PiecewiseConst (const std::vector<double > &box,
			       const double & refh);
  void deposit (const std::vector<std::vector<double > > & coord,
		const std::vector<std::vector<double > > & force);
  void calculate ();
  void clear ();
public:
  const double & getProfile (const unsigned & ix,
			     const unsigned & iy,
			     const unsigned & iz) const
      {return profile[index3to1(ix, iy, iz)];}
  const unsigned & getNx () const {return nx;}
  const unsigned & getNy () const {return ny;}
  const unsigned & getNz () const {return nz;}
  const std::vector<double > & getBox () const {return boxsize;}
public:
  void print_x  (const char * filename) const;
  void print_xy (const char * filename) const;
  void print_x_avg (const char * filename) const;
};

    
unsigned ErrorProfile_PiecewiseConst::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void ErrorProfile_PiecewiseConst::
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
