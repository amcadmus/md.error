#include "ElectrostaticInteraction.h"
#include "Polynominal.h"
#include "VectorOperation.h"
#include <numeric>
#include <math.h>
#include <time.h>
#include "Integral1D.h"
#include "ToolBox.h"

void ElectrostaticInteraction_rec_ES::calV()
{
  V = vecA[0][0] * (vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) - 
      vecA[0][1] * (vecA[1][0]*vecA[2][2] - vecA[2][0]*vecA[1][2]) +
      vecA[0][2] * (vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]);
}
  
  
void ElectrostaticInteraction_rec_ES::calAStar ()
{
  vecAStar.resize (3);
  vecAStar[0].resize (3);
  vecAStar[1].resize (3);
  vecAStar[2].resize (3);
  vecAStar[0][0] =( vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) / V;
  vecAStar[1][1] =( vecA[0][0]*vecA[2][2] - vecA[2][0]*vecA[0][2]) / V;
  vecAStar[2][2] =( vecA[0][0]*vecA[1][1] - vecA[1][0]*vecA[0][1]) / V;
  vecAStar[1][0] =(-vecA[1][0]*vecA[2][2] + vecA[2][0]*vecA[1][2]) / V;
  vecAStar[2][0] =( vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]) / V;
  vecAStar[0][1] =(-vecA[0][1]*vecA[2][2] + vecA[2][1]*vecA[0][2]) / V;
  vecAStar[2][1] =(-vecA[0][0]*vecA[2][1] + vecA[2][0]*vecA[0][1]) / V;
  vecAStar[0][2] =( vecA[0][1]*vecA[1][2] - vecA[1][1]*vecA[0][2]) / V;
  vecAStar[1][2] =(-vecA[0][0]*vecA[1][2] + vecA[1][0]*vecA[0][2]) / V;
}


void ElectrostaticInteraction_rec_ES::init (const std::vector<std::vector<value_type > > &vecA_, 
						 const std::vector<unsigned > K_,
						 const value_type & beta_)
{
  K = K_;
  beta = beta_;
  vecA = vecA_;
  calV();
  calAStar();
}



value_type ElectrostaticInteraction_rec_ES::calPotential (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
//   std::cout << "cal poten exact with beta " << beta << std::endl;
  double poten = 0;
  for (int m0 = -int(K[0]/2); m0 <= int(K[0]/2); m0 ++){
    for (int m1 = -int(K[1]/2); m1 <= int(K[1]/2); m1 ++){
      for (int m2 = -int(K[2]/2); m2 <= int(K[2]/2); m2 ++){
	if (m0 == 0 && m1 == 0 && m2 == 0 ) continue;
	std::vector<value_type > m (3);
	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	double mm = VectorOperation::dot (m, m);
	double expp = exp (-M_PI * M_PI * mm / beta / beta) / mm;
	std::vector<value_type > sm (2);
	sm[0] = 0;
	sm[1] = 0;
	for (std::vector<StandardParticle * >::iterator ppart = beginPool;
	     ppart != endPool; ppart ++){
	  double tmp = 2 * M_PI * VectorOperation::dot (m, (*ppart)->r());
	  sm[0] += (*ppart)->charge() * cos(tmp);
	  sm[1] += (*ppart)->charge() * sin(tmp);
	}
	poten += expp * (sm[0] * sm[0] + sm[1] * sm[1]);
      }
    }
  }
  poten /= 2 * M_PI * V;
  return poten;
}



void ElectrostaticInteraction_rec_ES::applyInteraction (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{  
  for (int m0 = -int(K[0]/2); m0 <= int(K[0]/2); m0 ++){
    for (int m1 = -int(K[1]/2); m1 <= int(K[1]/2); m1 ++){
      std::cout << "m0 is " << m0 << " m1 is " << m1  <<"     \r" ;
      std::cout << std::flush;
      for (int m2 = -int(K[2]/2); m2 <= int(K[2]/2); m2 ++){
	if (m0 == 0 && m1 == 0 && m2 == 0 ) continue;
	std::vector<value_type > m (3);
	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	double mm = VectorOperation::dot (m, m);
	double expp = exp (-M_PI * M_PI * mm / beta / beta) / mm;
	double sum0 (0);
	double sum1 (0);
	for (std::vector<StandardParticle * >::iterator ppart = beginPool;
	     ppart != endPool; ppart ++){
	  double tmp = - 2 * M_PI * VectorOperation::dot (m, (*ppart)->r());
	  sum0 += (*ppart)->charge() * cos(tmp);
	  sum1 += (*ppart)->charge() * sin(tmp);
	}
	int count = 0;
	for (std::vector<StandardParticle * >::iterator ppart = beginPool;
	     ppart != endPool; ppart ++, count++){
	  double tmp = 2 * M_PI * VectorOperation::dot (m, (*ppart)->r());
	  double a ((*ppart)->charge() * cos(tmp));
	  double b ((*ppart)->charge() * sin(tmp));
	  double coeff ((sum0 * b + sum1 * a) * 2./V * expp);  
	  VectorOperation::add ((*ppart)->f(), coeff, m);
	  // if (count == 0){
	  //   std::cout << (*ppart)->f()[1] << std::endl;
	  // }
	}
// 	for (std::vector<StandardParticle * >::iterator ppart = beginPool;
// 	     ppart != endPool; ppart ++){
// 	  std::vector<double > diff = (*beginPool)->r();
// 	  VectorOperation::add (diff, -1, (*ppart)->r());
// 	  double tmp = 2 * M_PI * VectorOperation::dot (m, diff);
// 	  sum += (*ppart)->charge() * sin (tmp);
// 	}
// 	sum *= (*beginPool)->charge();
// 	force[0] += expp * sum * 2./ V * m[0];
// 	force[1] += expp * sum * 2./ V * m[1];
// 	force[2] += expp * sum * 2./ V * m[2];
      }
    }
  }
//   std::cout << "tmp force calculater: force Of part" << std::endl;
//   std::cout << (partPool[0])->r() [0] << '\t'
// 	    << (partPool[0])->r() [1] << '\t'
// 	    << (partPool[0])->r() [2] << '\n';
//   std::cout << "the force is " << std::endl;
//   std::cout << force [0] << '\t'
// 	    << force [1] << '\t'
// 	    << force [2] << '\n';
//   std::cout << std::endl;
  
}


value_type ElectrostaticInteraction_rec_ES::applyInteractionCalPotential  (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{  
  applyInteraction (beginPool, endPool);
  return calPotential (beginPool, endPool);
}

