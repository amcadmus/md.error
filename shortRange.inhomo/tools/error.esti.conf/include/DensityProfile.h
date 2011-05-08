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

#define CPLUSPLUS
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

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
  void reinit_xtc (const std::string & filename,
		   const double & refh);
  const double & getProfile (const unsigned & ix,
			     const unsigned & iy,
			     const unsigned & iz) const
      {return profile[index3to1(ix, iy, iz)];}
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


//
// consider correlation
//
class DensityProfile_Corr_PiecewiseConst
{
  std::vector<double > boxsize;
  int nx, ny, nz;
  double hx, hy, hz;
  double zero;
  int corrBond, corrDim;
  int numCorr;
  std::vector<double > mean;
  std::vector<std::vector<double > > corr;
public:
  inline int  index3to1 (int  ix, int  iy, int  iz) const;
  inline void index1to3 (int& input,
			 int& ix, int& iy, int& iz) const;
  inline int  corrIndex3to1 (int  ix, int  iy, int  iz) const;
  inline void corrIndex1to3 (int& input,
			     int& ix, int& iy, int& iz) const;
public:
  DensityProfile_Corr_PiecewiseConst (const std::string & filename,
				      const int & corrBond,
				      const double & refh);
  void reinit_xtc (const std::string & filename,
		   const int & corrBond,
		   const double & refh);
  const double & getMean (const int & ix,
			  const int & iy,
			  const int & iz) const
      {return mean[index3to1(ix, iy, iz)];}
  const double & getMean (const int & i) const
      {return mean[i];}
  const double & getCorr (const int & ix,
			  const int & iy,
			  const int & iz,
			  const int & dx,
			  const int & dy,
			  const int & dz) const
      {return corr[index3to1(ix, iy, iz)][corrIndex3to1(dx, dy, dz)];}
  const double & getCorr (const int & i,
			  const int & j) const
      {return corr[i][j];}
  inline const double & getMean (const double & xx,
				 const double & yy,
				 const double & zz) const;
  inline const double & getCorr (const double & xx,
				 const double & yy,
				 const double & zz,
				 const double & dx,
				 const double & dy,
				 const double & dz) const;
  const int & getNx () const {return nx;}
  const int & getNy () const {return ny;}
  const int & getNz () const {return nz;}
  const std::vector<double > & getBox () const {return boxsize;}
  const int & getCorrBond () const {return corrBond;}
  const int & getCorrDim  () const {return corrDim; }
public:
  void print_x  (const std::string & filename) const;
  void print_xy (const std::string & filename) const;
};


int DensityProfile_Corr_PiecewiseConst::
corrIndex3to1 (int  ix, int  iy, int  iz) const
{
  ix += corrBond;
  iy += corrBond;
  iz += corrBond;
  return iz + corrDim * (iy + corrDim * ix);
}

void DensityProfile_Corr_PiecewiseConst::
corrIndex1to3 (int& input,
	       int& ix, int& iy, int& iz) const
{
  int tmp = input;
  iz = tmp % (corrDim);
  tmp = (tmp - iz) / corrDim;
  iy = tmp % (corrDim);
  ix =  (tmp - iy) / corrDim;
  ix -= corrBond;
  iy -= corrBond;
  iz -= corrBond;
}
    
int DensityProfile_Corr_PiecewiseConst::
index3to1 (int  ix, int  iy, int  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void DensityProfile_Corr_PiecewiseConst::
index1to3 (int& input,
	   int& ix, int& iy, int& iz) const
{
  int tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}

const double & DensityProfile_Corr_PiecewiseConst::
getMean (const double & xx,
	 const double & yy,
	 const double & zz) const
{
  int ix = int (xx / hx);
  int iy = int (yy / hy);
  int iz = int (zz / hz);
  return getMean (ix, iy, iz);
}

const double & DensityProfile_Corr_PiecewiseConst::
getCorr (const double & xx,
	 const double & yy,
	 const double & zz,
	 const double & dx,
	 const double & dy,
	 const double & dz) const
{
  int ix = int (xx / hx);
  int iy = int (yy / hy);
  int iz = int (zz / hz);
  int idx = int (dx / hx);
  int idy = int (dy / hy);
  int idz = int (dz / hz);
  if (-idx > corrBond || idx > corrBond) return zero;
  if (-idy > corrBond || idy > corrBond) return zero;
  if (-idz > corrBond || idz > corrBond) return zero;
  return getCorr (ix, iy, iz, idx, idy, idz);
}


#endif
