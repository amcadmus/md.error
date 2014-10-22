#ifndef __nb_interaction_lj_h_wanghan__
#define __nb_interaction_lj_h_wanghan__

#include "nb_interaction_acc.h"
#include "nb_interaction_vdw.h"

namespace LennardJones6_12{
  typedef enum paramIndex {
    epsilon		,
    sigma		,
    nparam
  } paramIndex_t;
}

template <typename Acceleration>
class VanderWaals_Interaction <Acceleration, nb_interaction_vanderWaals_cutoff_tag>
{
  typedef enum paramIndex {
    epsilon		,
    sigma		,
    nparam
  } paramIndex_t;
  
  static typename Acceleration::DataType one_over_12;
  static typename Acceleration::DataType one_over_6;
public:
  inline typename Acceleration::DataType
  fscale (const typename Acceleration::DataType & rinvsq,
	  const typename Acceleration::DataType & c6,
	  const typename Acceleration::DataType & c12);
}
    ;

template<typename Acceleration>
inline typename Acceleration::DataType
VanderWaals_Interaction<Acceleration, nb_interaction_vanderWaals_cutoff_tag>::
fscale (const typename Acceleration::DataType & rinvsq,
	const typename Acceleration::DataType & c6,
	const typename Acceleration::DataType & c12)
{
  typename Acceleration::DataType rinv6, vvdw6, vvdw12;
  
  rinv6  = mul<Acceleration>(mul<Acceleration>(rinvsq, rinvsq), rinvsq);
  vvdw6  = mul<Acceleration>(c6, rinv6);
  vvdw12 = mul<Acceleration>(c12, mul<Acceleration>(rinv6, rinv6));
  return  mul<Acceleration>(sub<Acceleration>(vvdw12, vvdw6), rinvsq);
}

#endif


