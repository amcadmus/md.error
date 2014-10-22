#ifndef __nb_interaction_lj_h_wanghan__
#define __nb_interaction_lj_h_wanghan__

#include "nb_interaction_acc.h"

namespace LennardJones6_12{
  typedef enum paramIndex {
    epsilon		,
    sigma		,
    nparam
  } paramIndex_t;
}

template <typename Acceleration>
class LennardJones_6_12 
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
  fscale (typename Acceleration::DataType r);
}
    ;

template<typename Acceleration>
inline typename Acceleration::DataType
LennardJones_6_12<Acceleration>::
fscale (typename Acceleration::DataType r)
{
}


#endif

