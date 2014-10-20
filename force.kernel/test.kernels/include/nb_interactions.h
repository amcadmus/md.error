#ifndef __nb_interactions_h_wanghan__
#define __nb_interactions_h_wanghan__

#include "nb_auxiliary.h"
#include "nb_interaction_lj.h"

// accelteration options:
struct nb_interaction_geometric_none_tag	{} ;
struct nb_interaction_geometric_tip3p_tag	{} ;
struct nb_interaction_geometric_tip4p_tag	{} ;

struct nb_interaction_accelteration_none_tag	{} ;
struct nb_interaction_accelteration_sse_tag	{} ;
struct nb_interaction_accelteration_sse2_tag	{} ;
struct nb_interaction_accelteration_sse4_tag	{} ;
struct nb_interaction_accelteration_avx_128_tag	{} ;
struct nb_interaction_accelteration_avx_256_tag	{} ;

// interaction control:
struct nb_interaction_electrostatic_null_tag	{} ;
struct nb_interaction_electrostatic_rf_tag	{} ;
struct nb_interaction_electrostatic_pme_tag	{} ;
struct nb_interaction_electrostatic_zm_tag	{} ;

struct nb_interaction_vanderWaals_null_tag	{} ;
struct nb_interaction_vanderWaals_cutoff_tag	{} ;
struct nb_interaction_vanderWaals_shift_tag	{} ;
struct nb_interaction_vanderWaals_pme_tag	{} ;

struct nb_interaction_compute_f_tag		{} ;
struct nb_interaction_compute_fe_tag		{} ;
struct nb_interaction_compute_fev_tag		{} ;
struct nb_interaction_compute_fv_tag		{} ;



template <int MDDIM, typename ValueType>
inline void nb_pair_force_impl_nocut (const ValueType *		eleParam,
				      const ValueType *		vdwParam,
				      const ValueType *		dof0,
				      const ValueType *		dof1,
				      const ValueType *		charge0,
				      const ValueType *		charge1,
				      const ValueType &		dist2,
				      const ValueType *		diff,
				      ValueType *			energy,
				      ValueType *			force,
				      ValueType *			fshift,
				      nb_interaction_geometric_none_tag,
				      nb_interaction_accelteration_none_tag,
				      nb_interaction_electrostatic_null_tag,
				      nb_interaction_vanderWaals_cutoff_tag,
				      nb_interaction_compute_f_tag)
{
  ValueType ri2 = 1./dist2;
  ValueType sri2 = vdwParam[LennardJones6_12_Cutoff::sigma] * vdwParam[LennardJones6_12_Cutoff::sigma] * ri2;
  ValueType sri6 = sri2*sri2*sri2;
  ValueType scalor = 24 * vdwParam[LennardJones6_12_Cutoff::epsilon] * (ValueType(2) * (sri6*sri6) - sri6) * ri2;

  nb_auxiliary_scalor_multiply <MDDIM, ValueType> (scalor, diff, force);
}


// template <typename ValueType,
// 	  typename ElectrostaticType,
// 	  typename VanderWaalsType,
// 	  typename ComputeMode>
// inline void nb_pair_force_impl (const ValueType *		eleParam,
// 				const ValueType *		vdwParam,
// 				const ValueType *		cutoff,
// 				const int *			idx0,
// 				const int *			idx1,
// 				const ValueType *		dof0,
// 				const ValueType *		dof1,
// 				const ValueType *		charge0,
// 				const ValueType *		charge1,
// 				ValueType *			energy,
// 				ValueType *			force,
// 				ValueType *			fshift,
// 				nb_interaction_geometric_none_tag,
// 				nb_interaction_accelteration_none_tag);

template <int MDDIM,
	  typename ValueType,
	  typename Acceleration,
	  typename ElectrostaticType,
	  typename VanderWaalsType,
	  typename ComputeMode>
inline void nb_pair_force_impl (const ValueType *		eleParam,
				const ValueType *		vdwParam,
				const ValueType *		cutoff,
				const int *			idx0,
				const int *			idx1,
				const ValueType *		dof0,
				const ValueType *		dof1,
				const ValueType *		charge0,
				const ValueType *		charge1,
				ValueType *			energy,
				ValueType *			force,
				ValueType *			fshift,
				nb_interaction_geometric_none_tag)
{
  ValueType diff[MDDIM];
  nb_auxiliary_diff <MDDIM, ValueType> (dof0, dof1, diff);
  ValueType dist2 = nb_auxiliary_dist2<MDDIM, ValueType>  (diff);
  if (dist2 > cutoff[0] * cutoff[0]) return;

  nb_pair_force_impl_nocut <MDDIM, ValueType>
      (eleParam, vdwParam, dof0, dof1, charge0, charge1, dist2, diff, energy, force, fshift,
       nb_interaction_geometric_none_tag (),
       Acceleration (),
       ElectrostaticType (),
       VanderWaalsType (),
       ComputeMode() );
}
				


#endif
