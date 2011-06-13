#ifndef __Integral3D_h__wanghan__
#define __Integral3D_h__wanghan__

#include <iostream>
#include <vector>
#include <list>
#include <math.h>

#define INIT_N 1

struct Integral3DInfo 
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
class Integral3D
{
  ValueType hx;
  ValueType hy;
  ValueType hz;
  IntegralBox3DGauss27 <TrinaryFunction,ValueType >  box_Gauss27;
public:
  ValueType cal_int(const Integral3DInfo::Method & method,
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

#endif


