#ifndef __Integral3D_h__wanghan__
#define __Integral3D_h__wanghan__

#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include "ToolBox.h"

#define INIT_N 1

struct ToolBox::Integral3DInfo 
{
public:
  typedef int Method;
  static Method Simpson;
  static Method Gauss27;
  static Method Gauss4;
}
    ;

template <typename TrinaryFunction, typename ValueType >
class IntegralBox3DGauss27;
template <typename TrinaryFunction, typename ValueType >
class IntegralBox3DGauss4;


template <typename TrinaryFunction, typename ValueType >
class ToolBox::Integral3D
{
  ValueType hx;
  ValueType hy;
  ValueType hz;
  IntegralBox3DGauss27 <TrinaryFunction,ValueType >  box_Gauss27;
public:
  ValueType cal_int(const ToolBox::Integral3DInfo::Method & method,
		    const TrinaryFunction & f,
                    const ValueType & xlo, const ValueType & xup,
                    const ValueType & ylo, const ValueType & yup,
		    const ValueType & zlo, const ValueType & zup,
		    const ValueType & tol = 1e-6,
                    const int & init_nx = INIT_N,
		    const int & init_ny = INIT_N,
		    const int & init_nz = INIT_N);
}
  ;


template <typename TrinaryFunction, typename ValueType >
class IntegralBox3DGauss27
{
public:
  IntegralBox3DGauss27 () {};
  void reinit (const TrinaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & ylo, const ValueType & yup,
	       const ValueType & zlo, const ValueType & zup,
	       const ValueType & tol = 1e-6);
  ValueType & coarseInt () {return coarse_int;}
  const ValueType & coarseInt () const {return coarse_int;}
  
  void cal_int (ValueType & value);
private:
  ValueType hx;
  ValueType hy;
  ValueType hz;
  ValueType xlo, xup;
  ValueType ylo, yup;
  ValueType zlo, zup;
  ValueType tol;
  const TrinaryFunction * f;
  
  ValueType coarse_int;
};


//////////////////////////////////////////////////
// implementation
//////////////////////////////////////////////////





ToolBox::Integral3DInfo::Method ToolBox::Integral3DInfo::Simpson	= 0;
ToolBox::Integral3DInfo::Method ToolBox::Integral3DInfo::Gauss27	= 1;

static double global_hx;
static double global_hy;
static double global_hz;

#define AA (.88729833462074168851)	// 0.5 * (1 + sqrt(0.6))
#define BB (.11270166537925831148)	// 0.5 * (1 - sqrt(0.6))

template <typename TrinaryFunction, typename ValueType >
void IntegralBox3DGauss27<TrinaryFunction, ValueType>::
reinit (const TrinaryFunction & f_, 
	const ValueType & xlo_, const ValueType & xup_,
	const ValueType & ylo_, const ValueType & yup_,
	const ValueType & zlo_, const ValueType & zup_,
	const ValueType & tol_)
{
  hx = (xup_ - xlo_);
  xlo = xlo_;
  xup = xup_;
  hy = (yup_ - ylo_);
  ylo = ylo_;
  yup = yup_;
  hz = (zup_ - zlo_);
  zlo = zlo_;
  zup = zup_;
  tol = tol_;
  f = &f_;
  
  ValueType x0 = AA * xlo + BB * xup;
  ValueType x1 = 0.5 * (xlo + xup);
  ValueType x2 = BB * xlo + AA * xup;
  ValueType y0 = AA * ylo + BB * yup;
  ValueType y1 = 0.5 * (ylo + yup);
  ValueType y2 = BB * ylo + AA * yup;
  ValueType z0 = AA * zlo + BB * zup;
  ValueType z1 = 0.5 * (zlo + zup);
  ValueType z2 = BB * zlo + AA * zup;
  coarse_int = 1./(18.*18.*18.) * hx * hy * hz * (
      125 * (f_(x0, y0, z0) + f_(x0, y0, z2) + 
	     f_(x0, y2, z0) + f_(x0, y2, z2) + 
	     f_(x2, y2, z0) + f_(x2, y2, z2) + 
	     f_(x2, y0, z0) + f_(x2, y0, z2) ) +
      200 * (f_(x1, y0, z0) + f_(x1, y0, z2) + 
	     f_(x1, y2, z0) + f_(x1, y2, z2) + 
	     f_(x0, y1, z0) + f_(x0, y1, z2) + 
	     f_(x2, y1, z0) + f_(x2, y1, z2) + 
	     f_(x0, y0, z1) + f_(x0, y2, z1) + 
	     f_(x2, y0, z1) + f_(x2, y2, z1) ) +
      320 * (f_(x0, y1, z1) + f_(x2, y1, z1) +
	     f_(x1, y0, z1) + f_(x1, y2, z1) +
	     f_(x1, y1, z0) + f_(x1, y1, z2) ) +
      512 *  f_(x1, y1, z1) );
  
}

template <typename TrinaryFunction, typename ValueType >
void IntegralBox3DGauss27<TrinaryFunction, ValueType>
::cal_int (ValueType & value)
{
  IntegralBox3DGauss27 box0;
  IntegralBox3DGauss27 box1;
  IntegralBox3DGauss27 box2;
  IntegralBox3DGauss27 box3;
  IntegralBox3DGauss27 box4;
  IntegralBox3DGauss27 box5;
  IntegralBox3DGauss27 box6;
  IntegralBox3DGauss27 box7;
  ValueType xmid = 0.5 * (xup + xlo);
  ValueType ymid = 0.5 * (yup + ylo);
  ValueType zmid = 0.5 * (zup + zlo);

  ValueType newtol = tol * 0.125;
  box0.reinit (*f, xlo, xmid, ylo, ymid, zlo, zmid, newtol);
  box1.reinit (*f, xlo, xmid, ylo, ymid, zmid, zup, newtol);
  box2.reinit (*f, xlo, xmid, ymid, yup, zlo, zmid, newtol);
  box3.reinit (*f, xlo, xmid, ymid, yup, zmid, zup, newtol);
  box4.reinit (*f, xmid, xup, ylo, ymid, zlo, zmid, newtol);
  box5.reinit (*f, xmid, xup, ylo, ymid, zmid, zup, newtol);
  box6.reinit (*f, xmid, xup, ymid, yup, zlo, zmid, newtol);
  box7.reinit (*f, xmid, xup, ymid, yup, zmid, zup, newtol);
  
  value = 
      box0.coarseInt () + box4.coarseInt() +
      box1.coarseInt () + box5.coarseInt() +
      box2.coarseInt () + box6.coarseInt() +
      box3.coarseInt () + box7.coarseInt();
  
  if (fabs(coarse_int - value) < tol * 127) {
    return;
  }
  else if (xup-xlo < 1e-13 * global_hx ||
	   yup-ylo < 1e-13 * global_hy || 
	   zup-zlo < 1e-13 * global_hz){
    // // std::cout << "lala"<< std::endl;
    return;
  }
  else {
    ValueType tmp;
    value = 0;
    box0.cal_int (tmp);
    value += tmp;
    box1.cal_int (tmp);
    value += tmp;
    box2.cal_int (tmp);
    value += tmp;
    box3.cal_int (tmp);
    value += tmp;
    box4.cal_int (tmp);
    value += tmp;
    box5.cal_int (tmp);
    value += tmp;
    box6.cal_int (tmp);
    value += tmp;
    box7.cal_int (tmp);
    value += tmp;
    return;
  }
}



template <typename TrinaryFunction, typename ValueType >
ValueType ToolBox::Integral3D<TrinaryFunction, ValueType >::
cal_int (const Integral3DInfo::Method & method,
	 const TrinaryFunction & f,
	 const ValueType & xlo, const ValueType & xup,
	 const ValueType & ylo, const ValueType & yup,
	 const ValueType & zlo, const ValueType & zup,
	 const ValueType & tol,
	 const int & init_nx,
	 const int & init_ny,
	 const int & init_nz)
{
  
  hx = (xup - xlo) / ValueType (init_nx);
  global_hx = hx;
  ValueType nxi = ValueType(1)/ValueType(init_nx);
  hy = (yup - ylo) / ValueType (init_ny);
  global_hy = hy;
  ValueType nyi = ValueType(1)/ValueType(init_ny);
  hz = (zup - zlo) / ValueType (init_nz);
  global_hz = hz;
  ValueType nzi = ValueType(1)/ValueType(init_nz);
  ValueType value = 0;

  ValueType mytol = tol * nxi * nyi * nzi;
  if (method == Integral3DInfo::Gauss27) {
    for (int i = 0; i < init_nx; ++i){
      for (int j = 0; j < init_ny; ++j){
	for (int k = 0; k < init_nz; ++k){
	  ValueType x = xlo + i * hx;
	  ValueType y = ylo + j * hy;
	  ValueType z = zlo + k * hz;
	  box_Gauss27.reinit (f, 
			      x, x+hx, 
			      y, y+hy,  
			      z, z+hz, 
			      mytol);
	  ValueType tmp;
	  box_Gauss27.cal_int (tmp);
	  value += tmp;
	}
      }
    }
  }
  
  return value;
}







#endif


