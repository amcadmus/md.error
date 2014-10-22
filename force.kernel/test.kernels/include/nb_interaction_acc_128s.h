#ifndef __nb_interaction_acc_128s_h_wanghan__
#define __nb_interaction_acc_128s_h_wanghan__

#include <emmintrin.h>
#include "nb_interaction_operators.h"

struct nb_interaction_accelteration_128s_tag	{
  typedef float		ValueType;
  typedef __m128	DataType;
  typedef int*		IndexType;
  static const int	index_stride = 4;
}
    ;

// template<>
// inline typename nb_interaction_accelteration_128s_tag::DataType
// init <float, nb_interaction_accelteration_128s_tag> (const float & a)
// {
//   return _mm_set1_ps (a);
// }
// template<>
// inline void 
// copy_index <nb_interaction_accelteration_128s_tag> (const int & input,
// 						    typename nb_interaction_accelteration_128s_tag::IndexType * idx)
// {
// }

template<typename Acceleration>
inline void
copy_index (const int & start_idx,
	    const int * nlist_data,
	    typename Acceleration::IndexType * idx)
{
  idx[0] = nlist_data[start_idx+0];
  idx[1] = nlist_data[start_idx+1];
  idx[2] = nlist_data[start_idx+2];
  idx[3] = nlist_data[start_idx+3];
}

template<>
inline void
get <nb_interaction_accelteration_128s_tag> (typename nb_interaction_accelteration_128s_tag::DataType in,
					     typename nb_interaction_accelteration_128s_tag::ValueType* out)
{
  _mm_store_ps (out, in);
}

template<>
inline typename nb_interaction_accelteration_128s_tag::DataType
init <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType & a)
{
  return _mm_set1_ps (a);
}

template<>
inline typename nb_interaction_accelteration_128s_tag::DataType
init <nb_interaction_accelteration_128s_tag> ()
{
  return _mm_setzero_ps ();
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
add <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType & a,
					     const typename nb_interaction_accelteration_128s_tag::DataType & b)
{
  return _mm_add_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
sub <nb_interaction_accelteration_128s_tag> (const nb_interaction_accelteration_128s_tag::DataType & a,
					     const nb_interaction_accelteration_128s_tag::DataType & b)
{
  return _mm_sub_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
mul <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType & a,
					     const typename nb_interaction_accelteration_128s_tag::DataType & b)
{
  return _mm_mul_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType
rsq<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType & x,
					    const typename nb_interaction_accelteration_128s_tag::DataType & y,
					    const typename nb_interaction_accelteration_128s_tag::DataType & z)
{
  return _mm_add_ps( _mm_add_ps( _mm_mul_ps(x, x), _mm_mul_ps(y, y) ), _mm_mul_ps(z, z) );
}

// rough estimate of 1/sqrt(a) with an extra N-R refinement
template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
invsqrt <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType & a)
{
  const __m128 onehalf  = _mm_set_ps(0.5, 0.5, 0.5, 0.5);
  const __m128 three = _mm_set_ps(3.0, 3.0, 3.0, 3.0);  
  __m128 esti = _mm_rsqrt_ps(a);
  return _mm_mul_ps(onehalf, _mm_mul_ps(_mm_sub_ps(three, _mm_mul_ps(_mm_mul_ps(esti, esti), a)), esti));
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
sq <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType & a)
{
  return _mm_mul_ps (a, a);
}

template <>
inline void
load_data_coord_full
<3, nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					    const typename nb_interaction_accelteration_128s_tag::IndexType iidx,
					    typename nb_interaction_accelteration_128s_tag::DataType * x,
					    typename nb_interaction_accelteration_128s_tag::DataType * y,
					    typename nb_interaction_accelteration_128s_tag::DataType * z)
{
  const float * __restrict__ p0 (dof+3*iidx[0]);
  const float * __restrict__ p1 (dof+3*iidx[1]);
  const float * __restrict__ p2 (dof+3*iidx[2]);
  const float * __restrict__ p3 (dof+3*iidx[3]);
  
  __m128 t1, t2, t3, t4, t5, t6, t7, t8;
  t1 = _mm_castpd_ps(_mm_load_sd((const double *)p0));
  t2 = _mm_castpd_ps(_mm_load_sd((const double *)p1));
  t3 = _mm_castpd_ps(_mm_load_sd((const double *)p2));
  t4 = _mm_castpd_ps(_mm_load_sd((const double *)p3));
  t5 = _mm_load_ss(p0+2);
  t6 = _mm_load_ss(p1+2);
  t7 = _mm_load_ss(p2+2);
  t8 = _mm_load_ss(p3+2);
  t1 = _mm_unpacklo_ps(t1, t2);
  t3 = _mm_unpacklo_ps(t3, t4);
  *x = _mm_movelh_ps(t1, t3);
  *y = _mm_movehl_ps(t3, t1);
  t5 = _mm_unpacklo_ps(t5, t6);
  t7 = _mm_unpacklo_ps(t7, t8);
  *z = _mm_movelh_ps(t5, t7);
}

template <>
inline void
update_data_force_full
<3, nb_interaction_accelteration_128s_tag> (typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					    const typename nb_interaction_accelteration_128s_tag::IndexType iidx,
					    typename nb_interaction_accelteration_128s_tag::DataType x,
					    typename nb_interaction_accelteration_128s_tag::DataType y,
					    typename nb_interaction_accelteration_128s_tag::DataType z)
{
  float * __restrict__ p0 (dof+3*iidx[0]);
  float * __restrict__ p1 (dof+3*iidx[1]);
  float * __restrict__ p2 (dof+3*iidx[2]);
  float * __restrict__ p3 (dof+3*iidx[3]);
  
  __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  t5          = _mm_unpacklo_ps(y, z);
  t6          = _mm_unpackhi_ps(y, z);
  t7          = _mm_shuffle_ps(x, t5, _MM_SHUFFLE(1, 0, 0, 0));
  t8          = _mm_shuffle_ps(x, t5, _MM_SHUFFLE(3, 2, 0, 1));
  t9          = _mm_shuffle_ps(x, t6, _MM_SHUFFLE(1, 0, 0, 2));
  t10         = _mm_shuffle_ps(x, t6, _MM_SHUFFLE(3, 2, 0, 3));
  t1          = _mm_load_ss(p0);
  t1          = _mm_loadh_pi(t1, (__m64 *)(p0+1));
  t1          = _mm_sub_ps(t1, t7);
  _mm_store_ss(p0, t1);
  _mm_storeh_pi((__m64 *)(p0+1), t1);
  t2          = _mm_load_ss(p1);
  t2          = _mm_loadh_pi(t2, (__m64 *)(p1+1));
  t2          = _mm_sub_ps(t2, t8);
  _mm_store_ss(p1, t2);
  _mm_storeh_pi((__m64 *)(p1+1), t2);
  t3          = _mm_load_ss(p2);
  t3          = _mm_loadh_pi(t3, (__m64 *)(p2+1));
  t3          = _mm_sub_ps(t3, t9);
  _mm_store_ss(p2, t3);
  _mm_storeh_pi((__m64 *)(p2+1), t3);
  t4          = _mm_load_ss(p3);
  t4          = _mm_loadh_pi(t4, (__m64 *)(p3+1));
  t4          = _mm_sub_ps(t4, t10);
  _mm_store_ss(p3, t4);
  _mm_storeh_pi((__m64 *)(p3+1), t4);
}


template <>
inline void
load_data_coord_one
<3, nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					    const int & iidx,
					    typename nb_interaction_accelteration_128s_tag::DataType * x,
					    typename nb_interaction_accelteration_128s_tag::DataType * y,
					    typename nb_interaction_accelteration_128s_tag::DataType * z)
{
  const float * __restrict__ p0 (dof+3*iidx);
  
  __m128 t1, t2;

  t1 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p0);
  t2 = _mm_load_ss(dof+2);

  *x = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(0, 0, 0, 0));
  *y = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(1, 1, 1, 1));
  *z = _mm_shuffle_ps(t2, t2, _MM_SHUFFLE(0, 0, 0, 0));
}

// template <>
// inline void
// load_data_coord_full (const IndexType idx,
// 		      const ValueType * dof,
// 		      DataType * x,
// 		      DataType * y,
// 		      DataType * z)
// {
  
// }




#endif
