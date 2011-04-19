#include "GenTubeAux.h"
#include <cmath>

AcceptRatio::
AcceptRatio (const double & Lx_,
	     const double & Lh_,
	     const double & Lt_,
	     const double & rhoh_,
	     const double & rhol_)
    : Lx (Lx_), Lh(Lh_), Lt(Lt_)
      // rhoh (rhoh_), rhol(rhol_)
{
  double Ll = (Lx - Lh - Lt - Lt) * 0.5;
  x1 = Ll;
  x2 = x1 + Lt;
  x3 = x2 + Lh;
  x4 = x3 + Lt;
  v0 = rhol_ / rhoh_;
}

double AcceptRatio::
operator () (const double & x) const 
{
  if (x < x1){
    return v0;
  }
  else if (x < x2){
    return v0 + (1. - v0) * sin ((x - (x1 + 0.5 * Lt)) / Lt * M_PI);
  }
  else if (x < x3 ){
    return 1.;
  }
  else if (x < x4){
    return v0 - (1. - v0) * sin ((x - (x3 + 0.5 * Lt)) / Lt * M_PI);
  }
  else {
    return v0;
  }
}


