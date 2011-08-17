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
  DensityFunction total;
public:
  void reinit (const std::vector<double > & x,
	       const std::vector<double > & p,
	       const std::vector<double > & n,
	       const double & Ly,
	       const double & Lz);
  void genConf_gro (const char * filename);
  void genXtc (const char * filename,
	       const int & nframe);
  void genConf_native (const char * filename);
};




double DensityFunction::
operator () (const double & x) const 
{
  return value (x);
}

#endif
