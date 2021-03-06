#pragma once

#include <vector>
#include <dlib/optimization.h>
#include "SimulationRegion.h"

typedef dlib::matrix<double,0,1> column_vector;
using namespace std;

class AbsLossFunc 
{
public:
  virtual ~AbsLossFunc () {};
  virtual double value (const column_vector & xx) = 0;
  virtual const column_vector deriv (const column_vector & xx) = 0;
  virtual bool evaluate (double & v,
			  column_vector & d,
			  const column_vector & xx) = 0;
  virtual int inputDim  () const = 0;
}
    ;

template <typename ErrorEsti> 
class LossFunc : public AbsLossFunc
{
public:
  LossFunc (const int & CC,
	    const int & nbins, 
	    const double & beta,
	    const vector<int> & KK,
	    const double & q2, 
	    const int & natoms, 
	    const SimulationRegion<double > & region,
	    const int & l_cut,
	    const int b_style,
	    const int numb_threads);
  virtual ~LossFunc () {};
  virtual double value (const column_vector & xx);
  virtual const column_vector deriv (const column_vector & xx);
  virtual bool evaluate (double & v,
			 column_vector & d,
			 const column_vector & xx);
  virtual int inputDim  () const {return 2 * (nbins - 1);};
public:
  int CC;
  int nbins;
  double beta;
  vector<int> KK;
  double q2;
  int natoms;
  const SimulationRegion<double > & region;
  int l_cut;
  int over_sampling;
  vector<typename ErrorEsti::HatType> hhc;
  ErrorEsti err;
};


template <typename ErrorEsti> 
LossFunc<ErrorEsti>::
LossFunc (const int & CC_,
	  const int & nbins_, 
	  const double & beta_,
	  const vector<int> & KK_,
	  const double & q2_, 
	  const int & natoms_, 
	  const SimulationRegion<double > & region_,
	  const int & l_cut_,
	  const int b_style,
	  const int numb_threads)
    : CC (CC_), 
      nbins (nbins_),
      beta (beta_),
      KK (KK_),
      q2 (q2_),
      natoms (natoms_),
      region (region_),
      l_cut (l_cut_),
      err (beta, KK, hhc, l_cut, b_style, numb_threads)
{
  over_sampling = 100 * l_cut;
  hhc.push_back (typename ErrorEsti::HatType (CC, nbins, KK[0], over_sampling));
  hhc.push_back (typename ErrorEsti::HatType (CC, nbins, KK[1], over_sampling));
  hhc.push_back (typename ErrorEsti::HatType (CC, nbins, KK[2], over_sampling));
}


template <typename ErrorEsti> 
double
LossFunc<ErrorEsti>::
value (const column_vector & xx)
{
  long size = xx.size();
  assert (size == 2 * (nbins - 1));
  
  vector<double > vv (nbins-1);
  vector<double > dd (nbins-1);
  for (long ii = 0; ii < nbins-1; ++ii){
    vv[ii] = xx(ii);
    dd[ii] = xx(ii+nbins-1);
  }
  
  hhc[0].set_value (vv, dd);
  hhc[1].set_value (vv, dd);
  hhc[2].set_value (vv, dd);

  double ret = err.estimate (q2, natoms, region);
  return ret;
}

template <typename ErrorEsti> 
const column_vector
LossFunc<ErrorEsti>::
deriv (const column_vector & xx)
{
  long size = xx.size();
  assert (size == 2 * (nbins - 1));
  
  vector<double > vv (nbins-1);
  vector<double > dd (nbins-1);
  for (long ii = 0; ii < nbins-1; ++ii){
    vv[ii] = xx(ii);
    dd[ii] = xx(ii+nbins-1);
  }

  hhc[0].set_value (vv, dd);
  hhc[1].set_value (vv, dd);
  hhc[2].set_value (vv, dd);

  vector<double > reslt;
  err.estimate (reslt, q2, natoms, region);

  assert (size == long(reslt.size()));
  
  column_vector ret (size);
  
  for (long ii = 0; ii < size; ++ii){
    ret(ii) = reslt[ii];
  }

  return ret;
}

template <typename ErrorEsti> 
bool
LossFunc<ErrorEsti>::
evaluate (double & rvv,
	  column_vector & rdd,
	  const column_vector & xx)
{
  long size = xx.size();
  assert (size == 2 * (nbins - 1));
  
  vector<double > vv (nbins-1);
  vector<double > dd (nbins-1);
  for (long ii = 0; ii < nbins-1; ++ii){
    vv[ii] = xx(ii);
    dd[ii] = xx(ii+nbins-1);
  }

  hhc[0].set_value (vv, dd);
  hhc[1].set_value (vv, dd);
  hhc[2].set_value (vv, dd);

  vector<double > reslt;
  rvv = err.estimate (reslt, q2, natoms, region);

  assert (size == long(reslt.size()));
  
  column_vector ret (size);
  
  for (long ii = 0; ii < size; ++ii){
    ret(ii) = reslt[ii];
  }
  rdd = ret;

  return true;
}




