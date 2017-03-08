#pragma once

#include <cmath>

template <typename VALUETYPE>
class Basis 
{
public:
  virtual VALUETYPE value      (const VALUETYPE & xx) const = 0;
  virtual VALUETYPE derivative (const VALUETYPE & xx) const = 0;	//
};

  template <typename VALUETYPE, int ORDER>
  class CardinalBSpline : public Basis<VALUETYPE>
  {
public:
    enum {MAX_ORDER = 14};
    enum {SUPPORT_SIZE = ORDER};
    typedef VALUETYPE		ParamType;
public:
    CardinalBSpline	 (const ParamType & cc = 1.5);
    virtual VALUETYPE value      (const VALUETYPE & xx) const;
    virtual VALUETYPE derivative (const VALUETYPE & xx) const;	// 
    VALUETYPE hatValue	 (const VALUETYPE & mm) const;	// Fourier trans on [0,1)
    int order () const {return ORDER;}
    int shift () const {return ORDER/2;}
private:
    VALUETYPE aa[ORDER+2][ORDER];	// with buff avoiding round-off problem
    VALUETYPE dd[ORDER+2][ORDER-1];
private:
    void init_value_coeffs ();
    void init_deriv_coeffs ();
  };

template <typename VALUETYPE, int ORDER>
VALUETYPE
CardinalBSpline<VALUETYPE, ORDER>::
value (const VALUETYPE & xx_) const
{
  VALUETYPE xx (xx_+shift());
  const VALUETYPE * tmpa (aa[int(xx+1)]);
  VALUETYPE buff (tmpa[ORDER-1]);
  for (int ii = ORDER-2; ii >= 0; --ii){
    buff = buff * xx + tmpa[ii];
  }
//  if (ORDER == 4) cout << xx_ << " " << buff << endl;
  return buff;
}

template <typename VALUETYPE, int ORDER>
VALUETYPE
CardinalBSpline<VALUETYPE, ORDER>::
derivative (const VALUETYPE & xx_) const
{
  VALUETYPE xx (xx_+shift());
  const VALUETYPE * tmpd (dd[int(xx+1)]);
  VALUETYPE buff (tmpd[ORDER-2]);
  for (int ii = ORDER-3; ii >= 0; --ii){
    buff = buff * xx + tmpd[ii];
  }
  return buff;
}

template <typename VALUETYPE, int ORDER>
VALUETYPE
CardinalBSpline<VALUETYPE, ORDER>::
hatValue (const VALUETYPE & mm) const
{
  if (mm == 0) return 1;
  double tmp (M_PI * mm);
  tmp = sin(tmp) / tmp;
  return pow (tmp, ORDER);	// lazy implementation...
}

#include "CardinalBSpline_Impl.h"

