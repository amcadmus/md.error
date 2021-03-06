#ifndef __DensityProfile_h_wanghan__
#define __DensityProfile_h_wanghan__

#include <vector>
#include <cmath>
#include <string>

#define CPLUSPLUS
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

#include "MDSystem.h"

class DensityProfile_PiecewiseConst
{
  std::vector<double > boxsize;
  unsigned nx, ny, nz, nele;
  double   hx, hy, hz;
  double * profile;
  int ndata;
  mutable FILE * fp_write;
public:
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
public:
  void reinit (const ScalorType & bx,
	       const ScalorType & by,
	       const ScalorType & bz,
	       const double & refh) ;
  void clearData ();
  void deposite (const CoordType * coord,
		 const IndexType numAtom);
  void calculate ();
public:
  DensityProfile_PiecewiseConst();
  ~DensityProfile_PiecewiseConst();
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
  void init_write (const char * file) const;
  void end_write () const;
  void write (const ScalorType & time) const;
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
