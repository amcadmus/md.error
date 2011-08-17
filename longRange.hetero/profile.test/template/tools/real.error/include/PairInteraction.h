#ifndef __wanghan_PairInteraction_h__
#define __wanghan_PairInteraction_h__

#include "Particle.h"
#include "BoxGeometry.h"

class PairInteraction 
{
public:
  virtual ~PairInteraction() {};
  virtual void applyInteraction (StandardParticle & p0, 
				 StandardParticle & p1, 
				 const BoxGeometry & boxg) = 0;
  virtual value_type calPotential (const StandardParticle & p0, 
				   const StandardParticle & p1, 
				   const BoxGeometry & boxg) = 0;
  virtual value_type applyInteractionCalPotential (StandardParticle & p0, 
						   StandardParticle & p1, 
						   const BoxGeometry & boxg) = 0;
}
    ;

typedef double value_type;
				    
class LennardJones6_12 : public PairInteraction
{
  value_type rcut;
  value_type rcut2;
  value_type epsilon;
  value_type sigma;
  value_type sigma6;
  value_type shift;
  value_type roff;
public:
  LennardJones6_12 (const value_type & rcutout = 2.5,
		    const value_type & epsilonout = 1,
		    const value_type & sigmaout = 1,
		    const value_type & shiftout = 0,
		    const value_type & roffout = 0)
      : rcut(rcutout),
	rcut2(rcut*rcut),
	epsilon(epsilonout),
	sigma(sigmaout),
	shift(shiftout),
	roff(roffout)
      {sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;};
  ~LennardJones6_12 () {};
  void applyInteraction (StandardParticle & p0, 
			 StandardParticle & p1, 
			 const BoxGeometry & boxg);
  value_type calPotential (const StandardParticle & p0, 
			   const StandardParticle & p1,
			   const BoxGeometry & boxg);
  value_type applyInteractionCalPotential (StandardParticle & p0, 
					   StandardParticle & p1, 
					   const BoxGeometry & boxg);
}
    ;


#endif
