#pragma once

inline double h00 (const double & xx)
{
  return (1 + 2 * xx) * (1 - xx) * (1 - xx);
}

inline double h10 (const double & xx) 
{
  return xx * (1 - xx) * (1 - xx);
}

inline double h01 (const double & xx) 
{
  return xx * xx * (3. - 2. * xx);
}

inline double h11 (const double & xx) 
{
  return xx * xx * (xx - 1);
}

class SymmHermiteBaseV
{
public:
  double value (const double & xx_,
		const int & knot,
		const double & hh) const 
      {
	double xx = fabs (xx_);
	if (xx < (knot - 1) * hh || xx > (knot + 1) * hh) return 0;
        if (xx < knot * hh) {
	  double x1 = xx / hh - double(knot - 1);
	  return h01 (x1);
	}
        else  {
	  double x1 = xx / hh - double(knot);
	  return h00 (x1);
	}
      }
};

class SymmHermiteBaseD
{
public:
  double value (const double & xx_,
		const int & knot,
		const double & hh) const 
      {
	double xx = fabs(xx_);
	if (xx < (knot - 1) * hh || xx > (knot + 1) * hh) return 0;
        if (xx < knot * hh) {
	  double x1 = xx / hh - double(knot - 1);
	  return h11 (x1) * hh;
	}
        else  {
	  double x1 = xx / hh - double(knot);
	  return h10 (x1) * hh;
	}
      }
};
