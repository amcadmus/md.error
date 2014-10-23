#ifndef __nb_interaction_operators_h_wanghan__
#define __nb_interaction_operators_h_wanghan__

template<typename Acceleration>
inline void
get (typename Acceleration::DataType in,
     typename Acceleration::ValueType * out);		// print the content of a for debug

template<typename Acceleration>
inline void
copy_index (const int & start_idx,
	    const int * nlist_data,
	    typename Acceleration::IndexType * idx);

template<typename Acceleration>
inline void
copy_index_capped (const int & start_idx,
		   const int * nlist_data,
		   typename Acceleration::IndexType * idx);

template<typename Acceleration>
inline typename Acceleration::DataType
cal_mask (const int & start_idx,
	  const int * nlist_data);

template<typename Acceleration>
inline typename Acceleration::DataType
apply_mask (const typename Acceleration::DataType & mask,
	    const typename Acceleration::DataType & a);

template<typename Acceleration>
inline void
index_table_trans (const typename Acceleration::IndexType * in,
		   const int * table,
		   typename Acceleration::IndexType * out);

template<typename Acceleration>
inline void
index_increase (typename Acceleration::IndexType * out,
		const int & inc);

template <typename Acceleration>
inline void
load_data_s2_afull (const typename Acceleration::ValueType * __restrict__ dof,
		    const typename Acceleration::IndexType input,
		    typename Acceleration::DataType * x,
		    typename Acceleration::DataType * y);

template <typename Acceleration>
inline void
load_data_s3_a1 (const typename Acceleration::ValueType * __restrict__ dof,
		 const int & idx,
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
		     const int &input,
		     typename Acceleration::DataType x,
		     typename Acceleration::DataType y,
		     typename Acceleration::DataType z);


template<typename Acceleration>
inline typename Acceleration::DataType
init (const typename Acceleration::ValueType & a);				// return all elements of a

template<typename Acceleration>
inline typename Acceleration::DataType
init (void);									// return 0

template<typename Acceleration>
inline typename Acceleration::DataType
add (const typename Acceleration::DataType & a,
     const typename Acceleration::DataType & b);				// a + b

template<typename Acceleration>
inline typename Acceleration::DataType
sub (const typename Acceleration::DataType & a,
     const typename Acceleration::DataType & b);				// a - b

template<typename Acceleration>
inline typename Acceleration::DataType
mul (const typename Acceleration::DataType & a,
     const typename Acceleration::DataType & b);				// a * b

template<typename Acceleration>
inline typename Acceleration::DataType
rsq (const typename Acceleration::DataType & x,
     const typename Acceleration::DataType & y,
     const typename Acceleration::DataType & z);				// x^2 + y^2 + z^2

template<typename Acceleration>
inline typename Acceleration::DataType
invsqrt (const typename Acceleration::DataType & a);				// 1/sqrt(a)

template<typename Acceleration>
inline typename Acceleration::DataType
sq (const typename Acceleration::DataType & a);					// a*a

// default implementation. Acc = none, native data types
// comment out to test code

// template<typename ValueType, typename Acceleration>
// inline typename Acceleration::DataType
// init (const ValueType & a)
// {
//   return typename Acceleration::DataType(a);
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// sub (const typename Acceleration::DataType & a,
//      const typename Acceleration::DataType & b) 
// {
//   return a-b;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// mul (const typename Acceleration::DataType & a,
//      const typename Acceleration::DataType & b)
// {
//   return a*b;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// rsq (const typename Acceleration::DataType & x,
//      const typename Acceleration::DataType & y,
//      const typename Acceleration::DataType & z)
// {
//   return x*x + y*y + z*z;
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// invsqrt (const typename Acceleration::DataType & a)
// {
//   return 1./sqrtf(a);
// }

// template<typename Acceleration>
// inline typename Acceleration::DataType
// sq (const typename Acceleration::DataType & a)
// {
//   return a*a;
// }

#endif
