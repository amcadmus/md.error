#include "ElectrostaticInteraction.h"

#include "VectorOperation.h"


ElectrostaticInteraction_dir::ElectrostaticInteraction_dir (const value_type & rcut_,
							    const value_type & beta_)
    : rcut (rcut_), beta (beta_), rcut2 (rcut*rcut), sqrtPIi (1./sqrt(M_PI))
{
}
  
value_type ElectrostaticInteraction_dir::calPotential (const StandardParticle & p0, 
						       const StandardParticle & p1,
						       const TimeConstRectangularBoxGeometry & boxg)
{
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return 0;
  value_type r = sqrt(rr);
  return p0.charge() * p1.charge() * erfc(beta * r) / r;
}

void ElectrostaticInteraction_dir::applyInteraction (StandardParticle & p0, 
						     StandardParticle & p1, 
						     const TimeConstRectangularBoxGeometry & boxg)
{
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return;
  value_type r = sqrt (rr);
  value_type ri = 1./r;
  value_type scale= (erfc(beta*r) * ri + 2*beta*sqrtPIi * exp (-beta*beta*rr)) / rr;
  scale *= p0.charge() * p1.charge();
  VectorOperation::add(p0.f(), +scale, diff);
  VectorOperation::add(p1.f(), -scale, diff);
}


value_type ElectrostaticInteraction_dir::applyInteractionCalPotential (StandardParticle & p0, 
								       StandardParticle & p1, 
								       const TimeConstRectangularBoxGeometry & boxg)
{
  std::vector<value_type > diff (NDIM, 0);
  VectorOperation::add (diff, 1., p0.r(), -1., p1.r());
  boxg.shortestImage (diff);
  value_type rr = VectorOperation::dot (diff, diff);
  if (rr > rcut2) return 0;
  value_type r = sqrt (rr);
  value_type ri = 1./r;
  value_type erfcBetaR = erfc (beta * r);
  value_type qiqj = p0.charge() * p1.charge();
  value_type scale= (erfcBetaR * ri + 2*beta*sqrtPIi*exp(-beta*beta*rr)) / rr * qiqj;
  VectorOperation::add(p0.f(), +scale, diff);
  VectorOperation::add(p1.f(), -scale, diff);
  
  return qiqj * erfcBetaR * ri;
}
