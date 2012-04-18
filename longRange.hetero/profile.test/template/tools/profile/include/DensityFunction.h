#ifndef __DensityFunction_h_wanghan__
#define __DensityFunction_h_wanghan__

#include <iostream>
#include <vector>
#include <cmath>

struct FunctionCoefficients
{
  double a0, a1, a2, a3;
};

class DensityFunction 
{
  unsigned mysize;
  double xlow, xup;
  std::vector<double > xx;
  std::vector<double > yy;
  std::vector<FunctionCoefficients> coeff;
public:
  void reinit (const std::vector<double > & x,
	       const std::vector<double > & y);
  double yMax () const;
  double integral () const;
  double value (const double & x) const;
  inline double Lx () const {return xup - xlow;}
  inline double operator () (const double & x) const;
}
    ;

class SystemDensityFunction
{
  unsigned np, nn, natom;
  double Lx, Ly, Lz;
  DensityFunction positive;
  DensityFunction negative;
  // DensityFunction total;
public:
  void reinit (const std::vector<double > & xp,
	       const std::vector<double > & xn,
	       const std::vector<double > & p,
	       const std::vector<double > & n,
	       const double & Ly,
	       const double & Lz);
  void genConf_gro (const char * filename);
  void genConf_rand_water_gro (const char * filename,
			       const std::vector<std::vector<double > > & posi);
  void genXtc (const char * filename,
	       const int & nframe);
  void genXtc_rand_water (const char * filename,
			  const int & nframe,
			  const std::vector<std::vector<double > > & posi);
  void genConf_native (const char * filename);
  unsigned get_natom () const {return natom;}
  unsigned get_natom_posi () const {return np;}
  unsigned get_natom_nega () const {return nn;}
};




double DensityFunction::
operator () (const double & x) const 
{
  return value (x);
}

#endif
