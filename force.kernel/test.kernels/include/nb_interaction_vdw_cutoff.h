#ifndef __nb_interaction_vdw_cutoff_h_wanghan__
#define __nb_interaction_vdw_cutoff_h_wanghan__

#include "nb_interaction_acc.h"
#include "nb_interaction_vdw.h"
#include "nb_interaction_compute.h"

namespace LennardJones6_12{
  typedef enum paramIndex {
    epsilon		,
    sigma		,
    nparam
  } paramIndex_t;
}

template <typename Acceleration,
	  typename ComputeMode>
class VanderWaals_Interaction <Acceleration,
			       nb_interaction_vanderWaals_full_tag,
			       ComputeMode>
{
  typename Acceleration::DataType one_over_12;
  typename Acceleration::DataType one_over_6;
  inline void
  cal_energy (const typename Acceleration::DataType & vvdw6,
	      const typename Acceleration::DataType & vvdw12,
	      typename Acceleration::DataType * energy,
	      nb_interaction_compute_energy_no_tag);
  inline void
  cal_energy (const typename Acceleration::DataType & vvdw6,
	      const typename Acceleration::DataType & vvdw12,
	      typename Acceleration::DataType * energy,
	      nb_interaction_compute_energy_yes_tag);	  
public:
  VanderWaals_Interaction (const VanderWaals_Control<Acceleration> & ctrl) {};
public:
  inline void
  force_energy (const typename Acceleration::DataType & rinv2,
		const typename Acceleration::DataType & c6,
		const typename Acceleration::DataType & c12,
		typename Acceleration::DataType * fscale,
		typename Acceleration::DataType * energy);
}
    ;

template<typename Acceleration,
	 typename ComputeMode>
inline void
VanderWaals_Interaction<Acceleration,
			nb_interaction_vanderWaals_full_tag,
			ComputeMode>::
force_energy (const typename Acceleration::DataType & rinv2,
	      const typename Acceleration::DataType & c6,
	      const typename Acceleration::DataType & c12,
	      typename Acceleration::DataType * fscale,
	      typename Acceleration::DataType * energy)
{
  typename Acceleration::DataType rinv6, vvdw6, vvdw12;
  
  rinv6   = mul<Acceleration>(mul<Acceleration>(rinv2, rinv2), rinv2);
  vvdw6   = mul<Acceleration>(c6, rinv6);
  vvdw12  = mul<Acceleration>(c12, mul<Acceleration>(rinv6, rinv6));
  cal_energy (vvdw6, vvdw12, energy, ComputeMode());
  *fscale = mul<Acceleration>(sub<Acceleration>(vvdw12, vvdw6), rinv2);
}

template<typename Acceleration,
	 typename ComputeMode>
inline void
VanderWaals_Interaction<Acceleration,
			nb_interaction_vanderWaals_full_tag,
			ComputeMode>::
cal_energy (const typename Acceleration::DataType & vvdw6,
	    const typename Acceleration::DataType & vvdw12,
	    typename Acceleration::DataType* energy,
	    nb_interaction_compute_energy_no_tag)
{
}

template<typename Acceleration,
	 typename ComputeMode>
inline void
VanderWaals_Interaction<Acceleration,
			nb_interaction_vanderWaals_full_tag,
			ComputeMode>::
cal_energy (const typename Acceleration::DataType & vvdw6,
	    const typename Acceleration::DataType & vvdw12,
	    typename Acceleration::DataType* energy,
	    nb_interaction_compute_energy_yes_tag)
{
  *energy = sub<Acceleration> (mul<Acceleration> (vvdw12, one_over_12) , mul<Acceleration> (vvdw6, one_over_6) );
}









// template<typename Acceleration>
// inline void
// VanderWaals_Interaction<Acceleration,
// 			nb_interaction_vanderWaals_full_tag,
// 			nb_interaction_compute_energy_no_tag>::
// force_energy (const typename Acceleration::DataType & rinv2,
// 	      const typename Acceleration::DataType & c6,
// 	      const typename Acceleration::DataType & c12,
// 	      typename Acceleration::DataType * fscale,
// 	      typename Acceleration::DataType * energy)
// {
//   typename Acceleration::DataType rinv6, vvdw6, vvdw12;
  
//   rinv6   = mul<Acceleration>(mul<Acceleration>(rinv2, rinv2), rinv2);
//   vvdw6   = mul<Acceleration>(c6, rinv6);
//   vvdw12  = mul<Acceleration>(c12, mul<Acceleration>(rinv6, rinv6));
//   *fscale = mul<Acceleration>(sub<Acceleration>(vvdw12, vvdw6), rinv2);
// }

// template<typename Acceleration>
// inline void
// VanderWaals_Interaction<Acceleration,
// 			nb_interaction_vanderWaals_full_tag,
// 			nb_interaction_compute_energy_yes_tag>::
// force_energy (const typename Acceleration::DataType & rinv2,
// 	      const typename Acceleration::DataType & c6,
// 	      const typename Acceleration::DataType & c12,
// 	      typename Acceleration::DataType * fscale,
// 	      typename Acceleration::DataType * energy)
// {
//   typename Acceleration::DataType rinv6, vvdw6, vvdw12;
  
//   rinv6   = mul<Acceleration> (mul<Acceleration> (rinv2, rinv2), rinv2);
//   vvdw6   = mul<Acceleration> (c6, rinv6);
//   vvdw12  = mul<Acceleration> (c12, mul<Acceleration> (rinv6, rinv6));
//   *fscale = mul<Acceleration> (sub<Acceleration> (vvdw12, vvdw6), rinv2);
//   // extra cost for energy at the very inner loop!
//   *energy = sub<Acceleration> (mul<Acceleration> (vvdw12, one_twelfth) , mul<Acceleration> (vvdw6, one_sixth) );
// }

// template <typename Acceleration,
// 	  typename ComputeMode>
// class VanderWaals_Interaction <Acceleration,
// 			       nb_interaction_vanderWaals_full_tag,
// 			       ComputeMode>
// VanderWaals_Interaction () 
// {
//   one_over_12 = init<Acceleration> (Acceleration::ValueType(1./12.));
//   one_over_6  = init<Acceleration> (Acceleration::ValueType(1./6. ));
// }

#endif


