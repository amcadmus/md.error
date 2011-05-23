#ifndef __Test_h_wanghan__
#define __Test_h_wanghan__

#include "Integral1D.h"
#include "Integral3D.h"
#include <gsl/gsl_sf_bessel.h>

void calculateGlobal ();

void inteClever (const double & epsilon,
		 const double & sigma,
		 const double & rc,
		 double & b0k,  double & b1k,  double & b2k,
		 double & a00k, double & a01k, double & a02k,
		 double & a10k, double & a11k, double & a22k);

class tmpf3 
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};


class InteNaive0
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive1
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive2
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};




class InteNaive00
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive01
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive02
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive11
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive12
{
public:
  double operator () (const double & r,
		      const double & theta,
		      const double & phi) const;
};

class InteNaive22
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

