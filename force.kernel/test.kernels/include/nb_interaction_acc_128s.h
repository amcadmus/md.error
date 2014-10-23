#ifndef __nb_interaction_acc_128s_h_wanghan__
#define __nb_interaction_acc_128s_h_wanghan__

#include <emmintrin.h>
#include "nb_interaction_operators.h"

struct nb_interaction_accelteration_128s_tag	{
private:
  struct int4 {int x, y, z, w;};
public:
  typedef float		ValueType;
  typedef __m128	DataType;
  typedef int4		IndexType;
  static const int	index_stride = 4;
}
    ;

template<>
inline void
copy_index<nb_interaction_accelteration_128s_tag> (const int   start_idx,
						   const int * nlist_data,
						   typename nb_interaction_accelteration_128s_tag::IndexType * idx)
{
  (*idx).x = nlist_data[start_idx+0];
  (*idx).y = nlist_data[start_idx+1];
  (*idx).z = nlist_data[start_idx+2];
  (*idx).w = nlist_data[start_idx+3];
}

template<>
inline void
copy_index_capped <nb_interaction_accelteration_128s_tag> (const int   start_idx,
							   const int * nlist_data,
							   typename nb_interaction_accelteration_128s_tag::IndexType * idx)
{
  (*idx).x = nlist_data[start_idx+0];
  (*idx).y = (nlist_data[start_idx+1] < 0) ? nlist_data[start_idx] : nlist_data[start_idx+1] ;
  (*idx).z = (nlist_data[start_idx+2] < 0) ? nlist_data[start_idx] : nlist_data[start_idx+2] ;
  (*idx).w = (nlist_data[start_idx+3] < 0) ? nlist_data[start_idx] : nlist_data[start_idx+3] ;
}

template<>
inline nb_interaction_accelteration_128s_tag::DataType
cal_mask <nb_interaction_accelteration_128s_tag> (const int   start_idx,
						  const int * nlist_data)
{
  return (nb_interaction_accelteration_128s_tag::DataType) (
      _mm_cmplt_epi32(_mm_loadu_si128((const __m128i *)(nlist_data+start_idx)),_mm_setzero_si128()));
}

template<>
inline nb_interaction_accelteration_128s_tag::DataType
apply_mask <nb_interaction_accelteration_128s_tag> (const nb_interaction_accelteration_128s_tag::DataType mask,
						    const nb_interaction_accelteration_128s_tag::DataType a)
{
  return _mm_andnot_ps(mask, a);
}

template<>
inline void
index_table_trans<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::IndexType * in,
							  const int * table,
							  typename nb_interaction_accelteration_128s_tag::IndexType * out)
{
  out->x = table[in->x];
  out->y = table[in->y];
  out->z = table[in->z];
  out->w = table[in->w];
}

template<>
inline void
index_increase<nb_interaction_accelteration_128s_tag> (typename nb_interaction_accelteration_128s_tag::IndexType * out,
						       const int inc)
{
  out->x += inc;
  out->y += inc;
  out->z += inc;
  out->w += inc;
}

template<>
inline void
get <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType in,
					     typename nb_interaction_accelteration_128s_tag::ValueType* out)
{
  _mm_storeu_ps (out, in);
}

template<>
inline typename nb_interaction_accelteration_128s_tag::DataType
init <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType a)
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
add <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType a,
					     const typename nb_interaction_accelteration_128s_tag::DataType b)
{
  return _mm_add_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
sub <nb_interaction_accelteration_128s_tag> (const nb_interaction_accelteration_128s_tag::DataType a,
					     const nb_interaction_accelteration_128s_tag::DataType b)
{
  return _mm_sub_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
mul <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType a,
					     const typename nb_interaction_accelteration_128s_tag::DataType b)
{
  return _mm_mul_ps (a, b);
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType
rsq<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType x,
					    const typename nb_interaction_accelteration_128s_tag::DataType y,
					    const typename nb_interaction_accelteration_128s_tag::DataType z)
{
  return _mm_add_ps( _mm_add_ps( _mm_mul_ps(x, x), _mm_mul_ps(y, y) ), _mm_mul_ps(z, z) );
}

// rough estimate of 1/sqrt(a) with an extra N-R refinement
template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
invsqrt <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType a)
{
  const __m128 onehalf  = _mm_set_ps(0.5, 0.5, 0.5, 0.5);
  const __m128 three = _mm_set_ps(3.0, 3.0, 3.0, 3.0);  
  __m128 esti = _mm_rsqrt_ps(a);
  return _mm_mul_ps(onehalf, _mm_mul_ps(_mm_sub_ps(three, _mm_mul_ps(_mm_mul_ps(esti, esti), a)), esti));
}

template <>
inline typename nb_interaction_accelteration_128s_tag::DataType 
sq <nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::DataType a)
{
  return _mm_mul_ps (a, a);
}


template <>
inline void
load_data_s3_a1
<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					 const int iidx,
					 typename nb_interaction_accelteration_128s_tag::DataType * x,
					 typename nb_interaction_accelteration_128s_tag::DataType * y,
					 typename nb_interaction_accelteration_128s_tag::DataType * z)
{
  const float * __restrict__ p0 (dof+3*iidx);
  
  __m128 t1, t2;

  t1 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p0);
  t2 = _mm_load_ss(p0+2);

  *x = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(0, 0, 0, 0));
  *y = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(1, 1, 1, 1));
  *z = _mm_shuffle_ps(t2, t2, _MM_SHUFFLE(0, 0, 0, 0));

}

template <>
inline void
load_data_s3_afull
<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					    const typename nb_interaction_accelteration_128s_tag::IndexType iidx,
					    typename nb_interaction_accelteration_128s_tag::DataType * x,
					    typename nb_interaction_accelteration_128s_tag::DataType * y,
					    typename nb_interaction_accelteration_128s_tag::DataType * z)
{
  const float * __restrict__ p0 (dof+3*iidx.x);
  const float * __restrict__ p1 (dof+3*iidx.y);
  const float * __restrict__ p2 (dof+3*iidx.z);
  const float * __restrict__ p3 (dof+3*iidx.w);
  
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
load_data_s2_afull
<nb_interaction_accelteration_128s_tag> (const typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					 const typename nb_interaction_accelteration_128s_tag::IndexType iidx,
					 typename nb_interaction_accelteration_128s_tag::DataType * x,
					 typename nb_interaction_accelteration_128s_tag::DataType * y)
{
  const float * __restrict__ p0 (dof+2*iidx.x);
  const float * __restrict__ p1 (dof+2*iidx.y);
  const float * __restrict__ p2 (dof+2*iidx.z);
  const float * __restrict__ p3 (dof+2*iidx.w);
  
  __m128 t1, t2, t3, t4;

  t1 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p0);   // x x c12 c6
  t2 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p1);
  t3 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p2);
  t4 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p3);
  t1 = _mm_unpacklo_ps(t1, t2);
  t2 = _mm_unpacklo_ps(t3, t4);
  *x = _mm_movelh_ps(t1, t2);
  *y = _mm_movehl_ps(t2, t1);
}

template <>
inline void
decrease_data_s3_afull
<nb_interaction_accelteration_128s_tag> (typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof,
					 const typename nb_interaction_accelteration_128s_tag::IndexType iidx,
					 typename nb_interaction_accelteration_128s_tag::DataType x,
					 typename nb_interaction_accelteration_128s_tag::DataType y,
					 typename nb_interaction_accelteration_128s_tag::DataType z)
{
  float * __restrict__ p0 (dof+3*iidx.x);
  float * __restrict__ p1 (dof+3*iidx.y);
  float * __restrict__ p2 (dof+3*iidx.z);
  float * __restrict__ p3 (dof+3*iidx.w);
  
  __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  t5  = _mm_unpacklo_ps(y, z);
  t6  = _mm_unpackhi_ps(y, z);
  t7  = _mm_shuffle_ps(x, t5, _MM_SHUFFLE(1, 0, 0, 0));
  t8  = _mm_shuffle_ps(x, t5, _MM_SHUFFLE(3, 2, 0, 1));
  t9  = _mm_shuffle_ps(x, t6, _MM_SHUFFLE(1, 0, 0, 2));
  t10 = _mm_shuffle_ps(x, t6, _MM_SHUFFLE(3, 2, 0, 3));
  t1  = _mm_load_ss(p0);
  t1  = _mm_loadh_pi(t1, (__m64 *)(p0+1));
  t1  = _mm_sub_ps(t1, t7);
  _mm_store_ss(p0, t1);
  _mm_storeh_pi((__m64 *)(p0+1), t1);
  t2  = _mm_load_ss(p1);
  t2  = _mm_loadh_pi(t2, (__m64 *)(p1+1));
  t2  = _mm_sub_ps(t2, t8);
  _mm_store_ss(p1, t2);
  _mm_storeh_pi((__m64 *)(p1+1), t2);
  t3  = _mm_load_ss(p2);
  t3  = _mm_loadh_pi(t3, (__m64 *)(p2+1));
  t3  = _mm_sub_ps(t3, t9);
  _mm_store_ss(p2, t3);
  _mm_storeh_pi((__m64 *)(p2+1), t3);
  t4  = _mm_load_ss(p3);
  t4  = _mm_loadh_pi(t4, (__m64 *)(p3+1));
  t4  = _mm_sub_ps(t4, t10);
  _mm_store_ss(p3, t4);
  _mm_storeh_pi((__m64 *)(p3+1), t4);
}

template <>
inline void
increase_data_s3_a1
<nb_interaction_accelteration_128s_tag> (typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ dof_,
					 typename nb_interaction_accelteration_128s_tag::ValueType * __restrict__ shift_,
					 const int index,
					 typename nb_interaction_accelteration_128s_tag::DataType x,
					 typename nb_interaction_accelteration_128s_tag::DataType y,
					 typename nb_interaction_accelteration_128s_tag::DataType z)
{
  float * __restrict__ dof (dof_+3*index);
  float * __restrict__ shift (shift_+3*index);

  __m128 t1, t2, t3;

  /* transpose data */
  t1 = x;
  _MM_TRANSPOSE4_PS(x, t1, y, z);
  x = _mm_add_ps(_mm_add_ps(x, t1), _mm_add_ps(y, z));

  t2 = _mm_load_ss(dof);
  t2 = _mm_loadh_pi(t2, (__m64 *)(dof+1));
  t3 = _mm_load_ss(shift);
  t3 = _mm_loadh_pi(t3, (__m64 *)(shift+1));

  t2 = _mm_add_ps(t2, x);
  t3 = _mm_add_ps(t3, x);

  _mm_store_ss(dof, t2);
  _mm_storeh_pi((__m64 *)(dof+1), t2);
  _mm_store_ss(shift, t3);
  _mm_storeh_pi((__m64 *)(shift+1), t3);  
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
