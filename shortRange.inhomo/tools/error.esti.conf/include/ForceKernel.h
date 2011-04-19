#ifndef __ForceKernel_h_wanghan__
#define __ForceKernel_h_wanghan__

#include<cmath>

class ForceKernel 
{
public:
  virtual double f2 (const double & dx,
		     const double & dy,
		     const double & dz) const = 0;
  virtual void f (const double & dx,
		  const double & dy,
		  const double & dz,
		  double & fx,
		  double & fy,
		  double & fz) const = 0;
};

class Disperson6 : public ForceKernel
{
  double epsilon;
  double sigma;
  double sigma6;
  double rc2;
public:
  Disperson6 (const double & epsilon,
	      const double & sigma,
	      const double & rc);
  inline double f2 (const double & dx,
		    const double & dy,
		    const double & dz) const;
  inline void f (const double & dx,
		 const double & dy,
		 const double & dz,
		 double & fx,
		 double & fy,
		 double & fz) const;
}
    ;



double Disperson6::
f2 (const double & dx,
    const double & dy,
    const double & dz) const
{
  double diff = dx*dx + dy*dy + dz*dx;
  if (diff <= rc2) return 0.;

  double tmp = 6. * epsilon * sigma6 / (diff * diff * diff * sqrt(diff));
  return tmp * tmp;
}

void Disperson6::
f (const double & dx,
   const double & dy,
   const double & dz,
   double & fx,
   double & fy,
   double & fz) const
{
  double diff = dx*dx + dy*dy + dz*dx;
  if (diff <= rc2) {
    fx = fy = fz = 0.;
    return ;
  }
  
  double tmp = 6. * epsilon * sigma6 / (diff * diff * diff * (diff));
  fx = dx * tmp;
  fy = dy * tmp;
  fz = dz * tmp;
}

  
  



#endif


