#ifndef __GenTubeAux_h_wanghan__
#define __GenTubeAux_h_wanghan__

struct AcceptRatio 
{
  double Lx, Lh, Lt, Ll;
  double x1, x2, x3, x4;
  // double rhol, rhoh;
  double v0;
public:
  AcceptRatio (const double & Lx,
	       const double & Lh,
	       const double & Lt,
	       const double & rhoh,
	       const double & rhol);
  double operator () (const double & x) const ;
};


#endif
