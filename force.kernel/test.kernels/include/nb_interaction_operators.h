#ifndef __nb_interaction_operators_h_wanghan__
#define __nb_interaction_operators_h_wanghan__

template<typename Acceleration>
inline void
get (const typename Acceleration::DataType in,
     typename Acceleration::ValueType * out);		// print the content of a for debug

template<typename Acceleration>
inline void
copy_index (const int   start_idx,
	    const int * nlist_data,
	    typename Acceleration::IndexType * idx);

template<typename Acceleration>
inline void
copy_index_capped (const int   start_idx,
		   const int * nlist_data,
		   typename Acceleration::IndexType * idx);

template<typename Acceleration>
inline typename Acceleration::DataType
cal_mask (const int   start_idx,
	  const int * nlist_data);

template<typename Acceleration>
inline typename Acceleration::DataType
apply_mask (const typename Acceleration::DataType mask,
	    const typename Acceleration::DataType a);

template<typename Acceleration>
inline void
index_table_trans (const typename Acceleration::IndexType in,
		   const int * table,
		   typename Acceleration::IndexType * out);

template<typename Acceleration>
inline void
index_increase (typename Acceleration::IndexType * out,
		const int inc);

template <typename Acceleration>
inline void
load_data_s2_afull (const typename Acceleration::ValueType * __restrict__ dof,
		    const typename Acceleration::IndexType input,
		    typename Acceleration::DataType * x,
		    typename Acceleration::DataType * y);

template <typename Acceleration>
inline void
load_data_s3_a1 (const typename Acceleration::ValueType * __restrict__ dof,
		 const int idx,
		 typename Acceleration::DataType * x,
		 typename Acceleration::DataType * y,
		 typename Acceleration::DataType * z);

template <typename Acceleration>
inline void
load_data_s3_afull (const typename Acceleration::ValueType * __restrict__ dof,
		    const typename Acceleration::IndexType input,
		    typename Acceleration::DataType * x,
		    typename Acceleration::DataType * y,
		    typename Acceleration::DataType * z);

template <typename Acceleration>
inline void
decrease_data_s3_afull (typename Acceleration::ValueType * __restrict__ dof,
			const typename Acceleration::IndexType input,
			typename Acceleration::DataType x,
			typename Acceleration::DataType y,
			typename Acceleration::DataType z);

template <typename Acceleration>
inline void
increase_data_s3_a1 (typename Acceleration::ValueType * __restrict__ dof,
		     typename Acceleration::ValueType * __restrict__ shift,
		     const int input,
		     typename Acceleration::DataType x,
		     typename Acceleration::DataType y,
		     typename Acceleration::DataType z);


template<typename Acceleration>
inline typename Acceleration::DataType
init (const typename Acceleration::ValueType a);			// return all elements of a

template<typename Acceleration>
inline typename Acceleration::DataType
init (void);								// return 0

template<typename Acceleration>
inline typename Acceleration::DataType
add (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b);				// a + b

template<typename Acceleration>
inline typename Acceleration::DataType
sub (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b);				// a - b

template<typename Acceleration>
inline typename Acceleration::DataType
mul (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b);				// a * b

template<typename Acceleration>
inline typename Acceleration::DataType
rsq (const typename Acceleration::DataType x,
     const typename Acceleration::DataType y,
     const typename Acceleration::DataType z);				// x^2 + y^2 + z^2

template<typename Acceleration>
inline typename Acceleration::DataType
invsqrt (const typename Acceleration::DataType a);			// 1/sqrt(a)

template<typename Acceleration>
inline typename Acceleration::DataType
sq (const typename Acceleration::DataType a);				// a*a


////////////////////////////////////////////////////////////////////////////////////////////////////
// default implementation of the operators for a single float or double
////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Acceleration>
inline void
get (const typename Acceleration::DataType in,
     typename Acceleration::ValueType * out)
{
  *out = static_cast<typename Acceleration::ValueType> (in);
}

template<typename Acceleration>
inline void
copy_index (const int   start_idx,
	    const int * nlist_data,
	    typename Acceleration::IndexType * idx)
{
  *idx = nlist_data[start_idx];
}

template<typename Acceleration>
inline void
copy_index_capped (const int   start_idx,
		   const int * nlist_data,
		   typename Acceleration::IndexType * idx)
{
  *idx = nlist_data[start_idx];		// the stride should be 1  
}

template<typename Acceleration>
inline typename Acceleration::DataType
cal_mask (const int   start_idx,
	  const int * nlist_data)
{
  return 0;
}


template<typename Acceleration>
inline typename Acceleration::DataType
apply_mask (const typename Acceleration::DataType mask,
	    const typename Acceleration::DataType a)
{
  return a;
}

template<typename Acceleration>
inline void
index_table_trans (const typename Acceleration::IndexType in,
		   const int * table,
		   typename Acceleration::IndexType * out)
{
  *out = table[in];
}

template<typename Acceleration>
inline void
index_increase (typename Acceleration::IndexType * out,
		const int inc)
{
  *out += inc;
}


template <typename Acceleration>
inline void
load_data_s2_afull (const typename Acceleration::ValueType * __restrict__ dof,
		    const typename Acceleration::IndexType idx,
		    typename Acceleration::DataType * x,
		    typename Acceleration::DataType * y)
{
  *x = dof[idx*2+0];
  *y = dof[idx*2+1];
}

template <typename Acceleration>
inline void
load_data_s3_a1 (const typename Acceleration::ValueType * __restrict__ dof,
		 const int idx,
		 typename Acceleration::DataType * x,
		 typename Acceleration::DataType * y,
		 typename Acceleration::DataType * z)
{
  *x = dof[idx*3];
  *y = dof[idx*3+1];
  *z = dof[idx*3+2];
}

template <typename Acceleration>
inline void
load_data_s3_afull (const typename Acceleration::ValueType * __restrict__ dof,
		    const typename Acceleration::IndexType idx,
		    typename Acceleration::DataType * x,
		    typename Acceleration::DataType * y,
		    typename Acceleration::DataType * z)
{
  *x = dof[idx*3];
  *y = dof[idx*3+1];
  *z = dof[idx*3+2];  
}

template <typename Acceleration>
inline void
decrease_data_s3_afull (typename Acceleration::ValueType * __restrict__ dof,
			const typename Acceleration::IndexType idx,
			typename Acceleration::DataType x,
			typename Acceleration::DataType y,
			typename Acceleration::DataType z)
{
  dof[idx*3+0] -= x;
  dof[idx*3+1] -= y;
  dof[idx*3+2] -= z;
}

template <typename Acceleration>
inline void
increase_data_s3_a1 (typename Acceleration::ValueType * __restrict__ dof,
		     typename Acceleration::ValueType * __restrict__ shift,
		     const int idx,
		     typename Acceleration::DataType x,
		     typename Acceleration::DataType y,
		     typename Acceleration::DataType z)
{
  dof[idx*3+0] += x;
  dof[idx*3+1] += y;
  dof[idx*3+2] += z;
  shift[idx*3+0] += x;
  shift[idx*3+1] += y;
  shift[idx*3+2] += z;
}


template<typename Acceleration>
inline typename Acceleration::DataType
init (const typename Acceleration::ValueType a)
{
  return static_cast <typename Acceleration::DataType> (a);
}

template<typename Acceleration>
inline typename Acceleration::DataType
init (void)
{
  return static_cast <typename Acceleration::DataType> (0.);
}

template<typename Acceleration>
inline typename Acceleration::DataType
add (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b)
{
  return a+b;
}

template<typename Acceleration>
inline typename Acceleration::DataType
sub (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b)
{
  return a-b;
}

template<typename Acceleration>
inline typename Acceleration::DataType
mul (const typename Acceleration::DataType a,
     const typename Acceleration::DataType b)
{
  return a*b;
}

template<typename Acceleration>
inline typename Acceleration::DataType
rsq (const typename Acceleration::DataType x,
     const typename Acceleration::DataType y,
     const typename Acceleration::DataType z)
{
  return x*x + y*y + z*z;
}

template<typename Acceleration>
inline typename Acceleration::DataType
invsqrt (const typename Acceleration::DataType a)
{
  return static_cast<typename Acceleration::DataType> (1./sqrt(a));
}

template<typename Acceleration>
inline typename Acceleration::DataType
sq (const typename Acceleration::DataType a)
{
  return a*a;
}





// default implementation. Acc = none, native data types
// comment out to test code

// template<typename ValueType, typename Acceleration>
// inline typename Acceleration::DataType
// init (const ValueType a)
// {
//   return typename Acceleration::DataType(a);
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// sub (const typename Acceleration::DataType a,
//      const typename Acceleration::DataType b) 
// {
//   return a-b;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// mul (const typename Acceleration::DataType a,
//      const typename Acceleration::DataType b)
// {
//   return a*b;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// rsq (const typename Acceleration::DataType x,
//      const typename Acceleration::DataType y,
//      const typename Acceleration::DataType z)
// {
//   return x*x + y*y + z*z;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// invsqrt (const typename Acceleration::DataType a)
// {
//   return 1./sqrtf(a);
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// sq (const typename Acceleration::DataType a)
// {
//   return a*a;
// }

#endif
