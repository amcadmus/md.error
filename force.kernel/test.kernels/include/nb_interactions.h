#ifndef __nb_interactions_h_wanghan__
#define __nb_interactions_h_wanghan__

#include "nb_auxiliary.h"
#include "nb_interaction_lj.h"
#include "nb_interaction_acc_128s.h"

// accelteration options:
struct nb_interaction_geometric_none_tag	{} ;
struct nb_interaction_geometric_tip3p_tag	{} ;
struct nb_interaction_geometric_tip4p_tag	{} ;

// interaction control:
struct nb_interaction_electrostatic_null_tag	{} ;
struct nb_interaction_electrostatic_rf_tag	{} ;
struct nb_interaction_electrostatic_pme_tag	{} ;
struct nb_interaction_electrostatic_zm_tag	{} ;

struct nb_interaction_compute_f_tag		{} ;
struct nb_interaction_compute_fe_tag		{} ;
struct nb_interaction_compute_fev_tag		{} ;
struct nb_interaction_compute_fv_tag		{} ;

template <int MDDIM, typename ValueType>
inline void
nb_pair_force_lj_scale (const ValueType *		vdwParam,
			const ValueType &		dist2,
			ValueType &			fscale,
			nb_interaction_accelteration_none_tag,
			nb_interaction_vanderWaals_null_tag)
{
}

template <int MDDIM, typename ValueType>
inline void
nb_pair_force_lj_scale (const ValueType *		vdwParam,
			const ValueType &		dist2,
			ValueType &			fscale,
			nb_interaction_accelteration_none_tag,
			nb_interaction_vanderWaals_cutoff_tag)
{
  ValueType ri2 = 1./dist2;
  ValueType sri2 = vdwParam[LennardJones6_12::sigma] * vdwParam[LennardJones6_12::sigma] * ri2;
  ValueType sri6 = sri2*sri2*sri2;
  fscale += 24 * vdwParam[LennardJones6_12::epsilon] * (ValueType(2) * (sri6*sri6) - sri6) * ri2;
}

template <int MDDIM, typename ValueType>
inline void
nb_pair_force_ele_scale (const ValueType *		eleParam,
			 const ValueType &		charge0,
			 const ValueType &		charge1,
			 const ValueType &		dist2,
			 ValueType &			fscale,
			 nb_interaction_accelteration_none_tag,
			 nb_interaction_electrostatic_null_tag)
{
}

template <int MDDIM, typename ValueType>
inline void
nb_pair_force_impl_nocut (const ValueType *		eleParam,
			  const ValueType *		vdwParam,
			  const ValueType *		dof0,
			  const ValueType *		dof1,
			  const ValueType *		charge0,
			  const ValueType *		charge1,
			  const ValueType &		dist2,
			  const ValueType *		diff,
			  ValueType &			energy,
			  ValueType &			fscale,
			  ValueType *			fshift,
			  nb_interaction_geometric_none_tag,
			  nb_interaction_accelteration_none_tag,
			  nb_interaction_electrostatic_null_tag,
			  nb_interaction_vanderWaals_cutoff_tag,
			  nb_interaction_compute_f_tag)
{
  ValueType ri2 = 1./dist2;
  ValueType sri2 = vdwParam[LennardJones6_12::sigma] * vdwParam[LennardJones6_12::sigma] * ri2;
  ValueType sri6 = sri2*sri2*sri2;
  fscale = 24 * vdwParam[LennardJones6_12::epsilon] * (ValueType(2) * (sri6*sri6) - sri6) * ri2;
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
	  typename ElectrostaticType,
	  typename VanderWaalsType,
	  typename ComputeMode>
inline void
nb_pair_force_impl (const ValueType *		eleParam,
		    const ValueType *		vdwParam,
		    const ValueType *		cutoff,
		    const int *			idx0,
		    const int *			idx1,
		    const ValueType *		dof,
		    const ValueType *		charge,
		    ValueType *			force,
		    ValueType *			energy,
		    ValueType *			fshift,
		    nb_interaction_geometric_none_tag,
		    nb_interaction_accelteration_none_tag)
{
  ValueType diff[MDDIM];
  nb_auxiliary_diff <MDDIM, ValueType> (&dof[idx0[0]*MDDIM], &dof[idx1[0]*MDDIM], diff);
  ValueType dist2 = nb_auxiliary_dist2<MDDIM, ValueType>  (diff);
  if (dist2 > cutoff[0] * cutoff[0]) return;

  ValueType fscale (0);
  // nb_pair_force_impl_nocut <MDDIM, ValueType>
  //     (eleParam, vdwParam, dof0, dof1, charge0, charge1, dist2, diff, *energy, fscale, fshift,
  //      nb_interaction_geometric_none_tag (),
  //      Acceleration (),
  //      ElectrostaticType (),
  //      VanderWaalsType (),
  //      ComputeMode() );
  nb_pair_force_ele_scale <MDDIM, ValueType>
      (eleParam, charge[idx0[0]], charge[idx1[0]], dist2, fscale, nb_interaction_accelteration_none_tag(), ElectrostaticType());  
  nb_pair_force_lj_scale <MDDIM, ValueType>
      (vdwParam, dist2, fscale, nb_interaction_accelteration_none_tag(), VanderWaalsType());
  
  nb_auxiliary_scalor_multiply_two <MDDIM, ValueType> ( fscale, diff, &force[idx0[0]*MDDIM], &force[idx1[0]*MDDIM]);
}


#endif
