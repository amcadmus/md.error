#ifndef __nb_interaction_vdw_h_wanghan__
#define __nb_interaction_vdw_h_wanghan__

#include "global_defines.h"

struct nb_interaction_vanderWaals_null_tag	{} ;
struct nb_interaction_vanderWaals_full_tag	{} ;
struct nb_interaction_vanderWaals_shift_tag	{} ;
struct nb_interaction_vanderWaals_pme_tag	{} ;

template <typename Acceleration,
	  typename VanderWaalsType,
	  typename ComputeMode>
class VanderWaals_Interaction
{
public:
  VanderWaals_Interaction (const VanderWaals_Control<Acceleration> & ctrl);
public:
  inline void
  force_energy (const typename Acceleration::DataType & rinv2,
		const typename Acceleration::DataType & c6,
		const typename Acceleration::DataType & c12,
		typename Acceleration::DataType * fscale,
		typename Acceleration::DataType * energy);
};

#endif

