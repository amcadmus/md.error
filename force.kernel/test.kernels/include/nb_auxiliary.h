#ifndef __nb_auxiliary_h_wanghan__
#define __nb_auxiliary_h_wanghan__

template <int MDDIM, typename ValueType>
inline void
nb_auxiliary_diff (const ValueType * dof0,
		   const ValueType * dof1,
		   ValueType * diff) 
{
  for (unsigned ii = 0; ii < MDDIM; ++ii){
    diff[ii] = dof1[ii] - dof0[ii];
  }
}

template <>
inline void
nb_auxiliary_diff<3, double> (const double * dof0,
			      const double * dof1,
			      double * diff)
{
  diff[0] = dof1[0] - dof0[0];
  diff[1] = dof1[1] - dof0[1];
  diff[2] = dof1[2] - dof0[2];
}

template <int MDDIM, typename ValueType>
inline ValueType
nb_auxiliary_dist2 (const ValueType * diff) 
{
  ValueType sum = ValueType(0);
  for (unsigned ii = 0; ii < MDDIM; ++ii){
    sum += diff[ii] * diff[ii];
  }
  return sum;
}

template <>
inline double
nb_auxiliary_dist2<3, double> (const double * diff)
{		   
  return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
}

template <int MDDIM, typename ValueType>
inline void
nb_auxiliary_scalor_multiply (const ValueType & scalor,
			      const ValueType * diff,
			      ValueType * b) 
{
  for (unsigned ii = 0; ii < MDDIM; ++ii){
    b[ii] = scalor * diff[ii];
  }
}

template <>
inline void
nb_auxiliary_scalor_multiply<3, double> (const double & scalor,
					 const double * diff,
					 double * b) 
{
  b[0] = scalor * diff[0];
  b[1] = scalor * diff[1];
  b[2] = scalor * diff[2];
}

template <int MDDIM, typename ValueType>
inline void
nb_auxiliary_scalor_multiply_two (const ValueType & scalor,
				  const ValueType * diff,
				  ValueType * b0,
				  ValueType * b1) 
{
  for (unsigned ii = 0; ii < MDDIM; ++ii){
    b1[ii] = -(b0[ii] = scalor * diff[ii]);
  }
}

template <>
inline void
nb_auxiliary_scalor_multiply_two<3, double> (const double & scalor,
					     const double * diff,
					     double * b0,
					     double * b1) 
{
  b1[0] = -(b0[0] = scalor * diff[0]);
  b1[1] = -(b0[1] = scalor * diff[1]);
  b1[2] = -(b0[2] = scalor * diff[2]);
}


#endif
