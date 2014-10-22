#ifndef __nb_interaction_vdw_h_wanghan__
#define __nb_interaction_vdw_h_wanghan__

struct nb_interaction_vanderWaals_null_tag	{} ;
struct nb_interaction_vanderWaals_cutoff_tag	{} ;
struct nb_interaction_vanderWaals_shift_tag	{} ;
struct nb_interaction_vanderWaals_pme_tag	{} ;

template <typename Acceleration, typename VanderWaalsType>
class VanderWaals_Interaction
{
public:
  inline typename Acceleration::DataType
  fscale (const typename Acceleration::DataType & rinvsq,
	  const typename Acceleration::DataType & c6,
	  const typename Acceleration::DataType & c12);  
};

#endif

