#include "ForceKernel.h"

Disperson6::
Disperson6 (const double & epsilon_,
	      const double & sigma_,
	    const double & rc_)
 : epsilon (epsilon_),
   sigma (sigma_),
   sigma6 (sigma_ * sigma_ * sigma_ * sigma_ * sigma_ * sigma_),
   rc (rc_),
   rc2 (rc_ * rc_)
{
}

