#include "Particle.h"
#include <algorithm>

void StandardParticle::clear ()
{
  std::fill (rBuff.begin(), rBuff.end(), 0.0);
  std::fill (vBuff.begin(), vBuff.end(), 0.0);
  std::fill (fBuff.begin(), fBuff.end(), 0.0);
  std::fill (noiBuff.begin(), noiBuff.end(), 0.0);
  massBuff = 0;
  chargeBuff = 0;
  idBuff = 0;
  typeBuff = 0;
}

  
StandardParticle::StandardParticle ()
    : rBuff(NDIM, 0.0), 
      vBuff(NDIM, 0.0), 
      fBuff(NDIM, 0.0), 
      noiBuff(NDIM, 0.0), 
      massBuff(0),
      chargeBuff(0),
      idBuff(0), 
      typeBuff(0)
{
}

StandardParticle::StandardParticle(const StandardParticle & part)
    : rBuff (part.rBuff),
      vBuff (part.vBuff),
      fBuff (part.fBuff),
      noiBuff (part.noiBuff),
      massBuff (part.massBuff),
      chargeBuff (part.chargeBuff),
      idBuff (part.idBuff),
      typeBuff (part.typeBuff)
{
}

