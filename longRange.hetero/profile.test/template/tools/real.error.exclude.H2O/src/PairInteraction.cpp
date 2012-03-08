#include "PairInteraction.h"
#include "VectorOperation.h"

typedef double value_type;

void LennardJones6_12::applyInteraction (StandardParticle & p0, 
					 StandardParticle & p1, 
					 const BoxGeometry & boxg)
{
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return;
  value_type r = sqrt (rr);
  value_type rsi = 1./(r - roff);
  value_type rsi2 = rsi * rsi;
  value_type term2 = sigma6 * rsi2 *  rsi2 *  rsi2 ;
  value_type scale = 24 * epsilon * (2 * term2*term2 - term2) * (rsi / r) ;
  VectorOperation::add(p0.f(), +scale, diff);
  VectorOperation::add(p1.f(), -scale, diff);
}



value_type LennardJones6_12::applyInteractionCalPotential (StandardParticle & p0, 
							   StandardParticle & p1, 
							   const BoxGeometry & boxg)
{
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return 0;
  value_type r = sqrt (rr);
  value_type rsi = 1./(r - roff);
  value_type rsi2 = rsi * rsi;
  value_type term2 = sigma6 * rsi2 *  rsi2 *  rsi2 ;
  value_type scale = 24 * epsilon * (2 * term2*term2 - term2) * (rsi / r) ;
  VectorOperation::add(p0.f(), +scale, diff);
  VectorOperation::add(p1.f(), -scale, diff);
  
  return 4 * epsilon * (term2*term2 - term2 + shift);
}

  
value_type LennardJones6_12::calPotential (const StandardParticle & p0, 
					   const StandardParticle & p1, 
					   const BoxGeometry & boxg)
{
  
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return 0;
  value_type r = sqrt (rr);
  value_type rsi = 1./(r - roff);
  value_type rsi2 = rsi * rsi;
  value_type term2 = sigma6 * rsi2 *  rsi2 *  rsi2 ;
  
  return 4 * epsilon * (term2*term2 - term2 + shift);
}
  
