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



  // {
  //   int count = 0;
  //   for (std::vector<StandardParticle * >::iterator ppart = beginPool;
  // 	 ppart != endPool; ppart ++, count++){
  //     std::vector<double > fcorr01(3, 0.);
  //     std::vector<double > fcorr02(3, 0.);
  //     std::vector<double > fcorr12(3, 0.);
  //     for (int m0 = -int(K[0]/2); m0 <= int(K[0]/2); m0 ++){
  // 	for (int m1 = -int(K[1]/2); m1 <= int(K[1]/2); m1 ++){
  // 	  for (int m2 = -int(K[2]/2); m2 <= int(K[2]/2); m2 ++){
  // 	    if (m0 == 0 && m1 == 0 && m2 == 0 ) continue;
  // 	    std::vector<value_type > m (3);
  // 	    m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
  // 	    m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
  // 	    m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
  // 	    double mm = VectorOperation::dot (m, m);
  // 	    double expp = exp (-M_PI * M_PI * mm / beta / beta) / mm;
  // 	    double sum0 (0);
  // 	    double sum1 (0);
  // 	    int count1 = 0;
  // 	    for (std::vector<StandardParticle * >::iterator ppart1 = beginPool;
  // 		 ppart1 != endPool && count1 < 2; ppart1 ++, ++count1){
  // 	      if (count1 == 0) continue;
  // 	      double tmp = - 2 * M_PI * VectorOperation::dot (m, (*ppart1)->r());
  // 	      sum0 += (*ppart1)->charge() * cos(tmp);
  // 	      sum1 += (*ppart1)->charge() * sin(tmp);
  // 	    }
  // 	    double tmp = 2 * M_PI * VectorOperation::dot (m, (*ppart)->r());
  // 	    double a ((*ppart)->charge() * cos(tmp));
  // 	    double b ((*ppart)->charge() * sin(tmp));
  // 	    double coeff ((sum0 * b + sum1 * a) * 2./V * expp);  
  // 	    VectorOperation::add (fcorr01, coeff, m);
  // 	  }
  // 	}
  //     }
  //     printf ("\ncount is %d,  box %f %f %f,  fcorr01: %f %f %f  fcorr02 %f %f %f\n",
  // 	      count,
  // 	      vecAStar[0][0],
  // 	      vecAStar[1][1],
  // 	      vecAStar[2][2],
  // 	      fcorr01[0], fcorr01[1], fcorr01[2],
  // 	      fcorr02[0], fcorr02[1], fcorr02[2]);
  //   }
  // }
  
  

  {    
    int count = 0;
    for (std::vector<StandardParticle * >::iterator ppart = beginPool;
  	 ppart != endPool; ++(++(++ppart)), ++count){
      std::cout << "corr " << count << " mol "  <<"     \r" ;
      std::cout << std::flush;
      std::vector<std::vector<value_type > > posi (3);
      std::vector<double > atomcharge(3);
      std::vector<StandardParticle * >::iterator patom = ppart;
      for (unsigned molc = 0; molc < 3; ++molc, ++patom){
  	posi[molc] = (*patom)->r();
  	atomcharge[molc] = (*patom)->charge();
      }
      std::vector<double > r02(3), r01(3), r12(3);
      for (unsigned dd = 0; dd < 3; ++dd){
  	r01[dd] = posi[0][dd] - posi[1][dd];
  	r02[dd] = posi[0][dd] - posi[2][dd];
  	r12[dd] = posi[1][dd] - posi[2][dd];
      }
      std::vector<double > fcorr01(3, 0.);
      std::vector<double > fcorr02(3, 0.);
      std::vector<double > fcorr12(3, 0.);
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
  	    double tmp, sum0, sum1, a, b, coeff;
	  
  	    // tmp = - 2 * M_PI * VectorOperation::dot (m, posi[1]);
  	    // sum0 = cos(tmp);
  	    // sum1 = sin(tmp);
  	    // tmp = 2 * M_PI * VectorOperation::dot (m, posi[0]);
  	    // a = cos(tmp);
  	    // b = sin(tmp);
  	    // coeff = atomcharge[0] * atomcharge[1] * ((sum0 * b + sum1 * a) * 2./V * expp);  
  	    // VectorOperation::add (fcorr01, coeff, m);

  	    // tmp = - 2 * M_PI * VectorOperation::dot (m, posi[2]);
  	    // sum0 = cos(tmp);
  	    // sum1 = sin(tmp);
  	    // tmp = 2 * M_PI * VectorOperation::dot (m, posi[0]);
  	    // a = cos(tmp);
  	    // b = sin(tmp);
  	    // coeff = atomcharge[0] * atomcharge[2] * ((sum0 * b + sum1 * a) * 2./V * expp);  
  	    // VectorOperation::add (fcorr02, coeff, m);

  	    // tmp = - 2 * M_PI * VectorOperation::dot (m, posi[2]);
  	    // sum0 = cos(tmp);
  	    // sum1 = sin(tmp);
  	    // tmp = 2 * M_PI * VectorOperation::dot (m, posi[1]);
  	    // a = cos(tmp);
  	    // b = sin(tmp);
  	    // coeff = atomcharge[1] * atomcharge[2] * ((sum0 * b + sum1 * a) * 2./V * expp);  
  	    // VectorOperation::add (fcorr12, coeff, m);


  	    tmp = 2 * M_PI * VectorOperation::dot (m, r01);
  	    b = sin(tmp);
  	    coeff = atomcharge[0] * atomcharge[1] * 2. / V * expp * b;
  	    VectorOperation::add (fcorr01, coeff, m);
  	    tmp = 2 * M_PI * VectorOperation::dot (m, r02);
  	    b = sin(tmp);
  	    coeff = atomcharge[0] * atomcharge[2] * 2. / V * expp * b;
  	    VectorOperation::add (fcorr02, coeff, m);
  	    tmp = 2 * M_PI * VectorOperation::dot (m, r12);
  	    b = sin(tmp);
  	    coeff = atomcharge[1] * atomcharge[2] * 2. / V * expp * b;
  	    VectorOperation::add (fcorr12, coeff, m);
  	  }
  	}
      }

      patom = ppart;
      (*patom)->f()[0] -= fcorr01[0];
      (*patom)->f()[1] -= fcorr01[1];
      (*patom)->f()[2] -= fcorr01[2];
      (*patom)->f()[0] -= fcorr02[0];
      (*patom)->f()[1] -= fcorr02[1];
      (*patom)->f()[2] -= fcorr02[2];
      ++patom;
      (*patom)->f()[0] -= fcorr12[0];
      (*patom)->f()[1] -= fcorr12[1];
      (*patom)->f()[2] -= fcorr12[2];
      (*patom)->f()[0] += fcorr01[0];
      (*patom)->f()[1] += fcorr01[1];
      (*patom)->f()[2] += fcorr01[2];
      ++patom;
      (*patom)->f()[0] += fcorr12[0];
      (*patom)->f()[1] += fcorr12[1];
      (*patom)->f()[2] += fcorr12[2];
      (*patom)->f()[0] += fcorr02[0];
      (*patom)->f()[1] += fcorr02[1];
      (*patom)->f()[2] += fcorr02[2];

      // printf ("\ncount is %d,  box %f %f %f,  charge0: %f, charge1 %f, charge2 %f,   fcorr01: %f %f %f  fcorr02 %f %f %f\n",
      // 	      count,
      // 	      vecAStar[0][0],
      // 	      vecAStar[1][1],
      // 	      vecAStar[2][2],
      // 	      atomcharge[0],
      // 	      atomcharge[1],
      // 	      atomcharge[2],
      // 	      fcorr01[0], fcorr01[1], fcorr01[2],
      // 	      fcorr02[0], fcorr02[1], fcorr02[2]);
    }
  }
}


value_type ElectrostaticInteraction_rec_ES::applyInteractionCalPotential  (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{  
  applyInteraction (beginPool, endPool);
  return calPotential (beginPool, endPool);
}

