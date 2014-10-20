#ifndef __nb_kernels_h_wanghan__
#define __nb_kernels_h_wanghan__

#include "nb_interactions.h"

template <int MDDIM,
	  typename ValueType,
	  typename Geometric,
	  typename Acceleration,
	  typename ElectrostaticType,
	  typename VanderWaalsType,
	  typename ComputeMode> 
inline void
nb_pair_force (const ValueType *		eleParam,
	       const ValueType *		vdwParam,
	       const ValueType *		cutoff,
	       const int *		idx0,
	       const int *		idx1,
	       const ValueType *		dof0,
	       const ValueType *		dof1,
	       const ValueType *		charge0,
	       const ValueType *		charge1,
	       ValueType *			energy,
	       ValueType *			force,
	       ValueType *			fshift)
{
  nb_pair_force_impl <MDDIM, ValueType,
		      Acceleration,
		      ElectrostaticType,
		      VanderWaalsType,
		      ComputeMode>
      (eleParam, vdwParam, cutoff, idx0, idx1, dof0, dof1, charge0, charge1, energy, force, fshift, Geometric());
}

#endif

