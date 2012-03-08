#include "ElectrostaticInteraction.h"
#include "Polynominal.h"
#include "VectorOperation.h"
#include <numeric>
#include <math.h>
#include <time.h>
#include "Integral1D.h"
#include "ToolBox.h"

static double global_x;
static double global_u;
static unsigned global_m;

typedef double (*MFP) (double );   

double intPow (const double & x, const unsigned & alpha_);

value_type intPow (const value_type & x, const unsigned & alpha_)
{
  value_type sum = 1;
  unsigned alpha = alpha_;
  for (; alpha > 0; --alpha){
    sum *= x;
  }
  return sum;
}
    
double value0 (double t) {
  if (t == 0) {return 0;}
  return (cos(2*M_PI*global_x/t) - 1) * intPow (t, 2*global_m) / t / t
      / intPow (global_u*t + 2*M_PI, 2*global_m);
}

double value1 (double t) {
  if (t == 0) {return 0;}
  return (sin(2*M_PI*global_x/t)) * intPow (t, 2*global_m) / t / t
      / intPow (global_u*t + 2*M_PI, 2*global_m);
}

double ElectrostaticInteraction_rec_BSpline::fenmu (const value_type & u, const unsigned & _2m)
{  
  if (u == 0) return 0;
  int A = 3;
  double fenmu = 0;
  for (int i = -A+1; i < A; ++i){
    fenmu += 1./intPow (u+2*M_PI*i, _2m);
  }
  fenmu += ( 1./ (2.*M_PI * (_2m-1) * intPow(u + 2*M_PI*(A-0.5), _2m-1)) -
	     1./ (2.*M_PI * (_2m-1) * intPow(u - 2*M_PI*(A-0.5), _2m-1)) );
  return fenmu;
}
double ElectrostaticInteraction_rec_BSpline::fenzi (const value_type & u, const unsigned & _2m)
{
  if (u == 0) return 0;
  int A = 3;
  double fenmu = 0;
  for (int i = -A+1; i < A; ++i){
    if (i == 0) continue;
    fenmu += 1./intPow (u+2*M_PI*i, _2m);
  }
  fenmu += ( 1./ (2.*M_PI * (_2m-1) * intPow(u + 2*M_PI*(A-0.5), _2m-1)) -
	     1./ (2.*M_PI * (_2m-1) * intPow(u - 2*M_PI*(A-0.5), _2m-1)) );
  return fenmu;
}

double ElectrostaticInteraction_rec_BSpline::avg_e2 (const value_type & u, const unsigned & _2m)
{
  if (u == 0) return 0;
  double fenmua = fenmu (u, _2m);
  double fenzia = fenzi (u, _2m);
  return (fenzi(u, 2*_2m) + fenzia * fenzia) / fenmua / fenmua;
}

double ElectrostaticInteraction_rec_BSpline::avg_e (const value_type & u, const unsigned & _2m)
{
  if (u == 0) return 0;
  int A = 3;
  double fenmu = 0;
  for (int i = -A+1; i < A; ++i){
    fenmu += 1./intPow (u+2*M_PI*i, _2m);
  }
  fenmu += ( 1./ (2.*M_PI * (_2m-1) * intPow(u + 2*M_PI*(A-0.5), _2m-1)) -
	     1./ (2.*M_PI * (_2m-1) * intPow(u - 2*M_PI*(A-0.5), _2m-1)) );
  double fenzi = fenmu - 1./intPow (u, _2m);
  return - fenzi / fenmu;
}


std::vector<double > error (const value_type & u, const value_type & x, const unsigned & _2m)
{
  global_m = _2m/2;
  global_u = u;
  global_x = x;
  int A = 3;
  
  if ( u == 0){
    return std::vector<double > (2,0);
  }
  double fenmu = 0;
  for (int i = -A+1; i < A; ++i){
    fenmu += 1./intPow (u+2*M_PI*i, _2m);
  }
  fenmu += ( 1./ (2.*M_PI * (_2m-1) * intPow(u + 2*M_PI*(A-0.5), _2m-1)) -
	     1./ (2.*M_PI * (_2m-1) * intPow(u - 2*M_PI*(A-0.5), _2m-1)) );
  
  std::vector<double > fenzi (2, 0);
  for (int i = -A+1; i < A; ++i){
    fenzi[0] += (cos(2*M_PI*i*x) - 1) / intPow (u+2*M_PI*i, _2m);
    fenzi[1] += (sin(2*M_PI*i*x)    ) / intPow (u+2*M_PI*i, _2m);
  }
//   Integral1D <MFP, value_type > inte;
//   std::cout << "est fenzi " << fenzi[0] << '\t' << fenzi[1] << std::endl;
//   fenzi[0] += inte.cal_int (Integral1DInfo::Gauss4, value0, -1./(A-0.5), 1./(A-0.5), 1e-10, 5*(unsigned(x)+1));
//   fenzi[1] += inte.cal_int (Integral1DInfo::Gauss4, value1, -1./(A-0.5), 1./(A-0.5), 1e-10, 5*(unsigned(x)+1));
  
//   std::cout << "est fenzi " << fenzi[0] << '\t' << fenzi[1] << std::endl;

  double c = fenzi[0] / fenmu;
  double d = fenzi[1] / fenmu;
//   double a = cos(u*x);
//   double b = sin(u*x);
  
  std::vector<double > tmp(2);
  tmp[0] = c;
  tmp[1] = d;
  return tmp;
}

value_type error1 (const value_type & u, const value_type & x, const unsigned & _2m)
{
  int A = 100;

  double fenmu = 0;
  for (int k = -A+1; k < A; ++k){
    fenmu += 1./intPow (u+2*M_PI*k, _2m);
  }

  std::vector<double > fenzi (2, 0);
  for (int k = -A+1; k < A; ++k){
    fenzi[0] += (cos(2*M_PI*k*x) - 1) / intPow (u+2*M_PI*k, _2m);
    fenzi[1] += (sin(2*M_PI*k*x)    ) / intPow (u+2*M_PI*k, _2m);
  }
//   std::cout << "exa fenzi " << fenzi[0] << '\t' << fenzi[1] << std::endl;

  double c = fenzi[0] / fenmu;
  double d = fenzi[1] / fenmu;
  double a = cos(u*x);
  double b = sin(u*x);
  
  return b*c + a*d;
}
  
	
value_type ElectrostaticInteraction_rec_BSpline::errorEstimate (
    const value_type & c2, const int & N)
{
//   unsigned n = Mn->getN();
  std::vector<double > KK (3);
  KK[0] = K[0];
  KK[1] = K[1];
  KK[2] = K[2];
  
  double sum = 0;
  std::vector<double > errorForce(3, 0);
  for (int m0 = -int(KK[0]/2); m0 < KK[0]/2; ++m0){
    for (int m1 = -int(KK[1]/2); m1 < KK[1]/2; ++m1){
      for (int m2 = -int(KK[2]/2); m2 < KK[2]/2; ++m2){
	if (fabs(m0) + fabs(m1) + fabs(m2) == 0) continue;
	std::vector<value_type > m (3);
	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	double mm = VectorOperation::dot (m, m);
// 	double ees = (avg_e2(2*M_PI*m0/KK[0], n) + avg_e2 (2*M_PI*m1/KK[1], n) + avg_e2 (2*M_PI*m2/KK[2], n) +
// 		      2 * avg_e (2*M_PI*m0/KK[0], n) * avg_e (2*M_PI*m1/KK[1], n) +
// 		      2 * avg_e (2*M_PI*m1/KK[1], n) * avg_e (2*M_PI*m2/KK[2], n) + 
// 		      2 * avg_e (2*M_PI*m0/KK[0], n) * avg_e (2*M_PI*m2/KK[2], n));
// 	double eemid = ( avg_e (2*M_PI*m0/KK[0], n) + 
// 			 avg_e (2*M_PI*m1/KK[1], n) + 
// 			 avg_e (2*M_PI*m2/KK[2], n) );
// 	eemid = eemid * eemid;
	
	double ees, eemid;
	std::vector<double > E1s(3);
	std::vector<double > E2s(3);
	ecal.get_e1 (m0, m1, m2, E1s[0], E1s[1], E1s[2]);
	ecal.get_e2 (m0, m1, m2, E2s[0], E2s[1], E2s[2]);
	ees = (E2s[0] + E2s[1] + E2s[2] + 
	       2 * E1s[0] * E1s[1] + 2 * E1s[0] * E1s[2] + 2 * E1s[1] * E1s[2]);
	eemid = (E1s[0] + E1s[1] + E1s[2]) * (E1s[0] + E1s[1] + E1s[2]);

	double expp = exp (- 2*M_PI*M_PI/beta/beta*mm) / (mm) * 2/V * 2/V;
	sum += 2 * expp * (ees + eemid);

// 	if (m0 ==11 && m1 == 0 && m2 == 0){
// 	  std::cout << 2 * expp * ees << std::endl;
// 	  std::cout << 2* expp * eemid << std::endl; 
// 	  std::cout << m0 << "\t" << m1 << "\t" << m2 << std::endl;
// 	}

      }
    } 
  }
  sum *= c2*c2/double(N);
  return sqrt(sum);
}
  

value_type ElectrostaticInteraction_rec_BSpline::errorEstimate (
    std::vector<StandardParticle * >::const_iterator begin,
    std::vector<StandardParticle * >::const_iterator end) 
{
  
  std::vector<double > qs;
  for (std::vector<StandardParticle * >::const_iterator pp = begin;
       pp != end; ++pp){
    qs.push_back ((*pp)->charge());
  }
  double c2 = 0;
  for (std::vector<double >::iterator i = qs.begin();
       i != qs.end(); ++ i){
    c2 += *i * *i;
  }
  
  return errorEstimate (c2, qs.size());
}



void ElectrostaticInteraction_rec_BSpline::applyInteraction (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
  watch.start();
  calQ(beginPool, endPool);
  // FILE * fp = fopen ("tmpQ.out", "w");
  // for (unsigned i = 0; i < K[0]*K[1]*K[2]; ++i){
  //   fprintf (fp, "%.12e\n", Q[i]);
  // }
  // fclose (fp);
  watch.stop();
  loop_time += watch.user();

  watch.start();
  fftw_execute (forwardQ);
  watch.stop();
  conv_time += watch.user();

  unsigned n = Mn->getN();
  unsigned size = K[0] * K[1] * K[2];
  unsigned sizeHalf = K[0] * K[1] * (K[2]/2+1);
  value_type sizei = 1./size;
  
  watch.start();
  for (unsigned i = 0; i < sizeHalf; i ++){
    QFProdPhiF0[i][0] = QF[i][0] * phiF0[i][0] - QF[i][1] * phiF0[i][1];
    QFProdPhiF0[i][1] = QF[i][0] * phiF0[i][1] + QF[i][1] * phiF0[i][0];
    QFProdPhiF1[i][0] = QF[i][0] * phiF1[i][0] - QF[i][1] * phiF1[i][1];
    QFProdPhiF1[i][1] = QF[i][0] * phiF1[i][1] + QF[i][1] * phiF1[i][0];
    QFProdPhiF2[i][0] = QF[i][0] * phiF2[i][0] - QF[i][1] * phiF2[i][1];
    QFProdPhiF2[i][1] = QF[i][0] * phiF2[i][1] + QF[i][1] * phiF2[i][0];
  }
  fftw_execute (backwardQFProdPhiF0);
  fftw_execute (backwardQFProdPhiF1);
  fftw_execute (backwardQFProdPhiF2);
  
  for (unsigned i = 0; i < size; i ++){
    QConvPhi0[i] *= sizei;
    QConvPhi1[i] *= sizei;
    QConvPhi2[i] *= sizei;
  }
  watch.stop();
  conv_time += watch.user();

  watch.start();
  double ii0 = 1./ value_type(K[0]);
  double ii1 = 1./ value_type(K[1]);
  double ii2 = 1./ value_type(K[2]);
  bool fast = ((n < K[0]) && (n < K[1]) && (n < K[2]));
  
  for (std::vector<StandardParticle * >::iterator ppart = beginPool;
       ppart != endPool; ppart ++){
    std::vector<value_type > u(3);
    u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*ppart)->r());
    u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*ppart)->r());
    u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*ppart)->r());
    int A0 = -int(floor ((u[0]) * ii0)) ;
    int A1 = -int(floor ((u[1]) * ii1)) ;
    int A2 = -int(floor ((u[2]) * ii2)) ;
    value_type posi0 = u[0] + A0 * K[0];
    value_type posi1 = u[1] + A1 * K[1];
    value_type posi2 = u[2] + A2 * K[2];
    value_type tmp0 = 0;
    value_type tmp1 = 0;
    value_type tmp2 = 0;
    value_type tmp;
    std::vector<double > force (3, 0);
    unsigned index0, index1;
    
    if (!fast){
      int count = 0;
      for (unsigned k0 = 0; k0 < K[0]; k0 ++){
	for (unsigned k1 = 0; k1 < K[1]; k1 ++){
	  for (unsigned k2 = 0; k2 < K[2]; k2 ++, count ++){
	    value_type n0 = floor ((u[0] - k0) * ii0);
	    if (u[0] - k0 - n0 * K[0] < n){
	      Mn->value (u[0] - k0 - n0 * K[0], tmp0);
	    }
	    else{
	      continue;
	    }
	    value_type n1 = floor ((u[1] - k1) * ii1);
	    if (u[1] - k1 - n1 * K[1] < n){
	      Mn->value (u[1] - k1 - n1 * K[1], tmp1);
	    }
	    else{
	      continue;
	    }
	    value_type n2 = floor ((u[2] - k2) * ii2);
	    if (u[2] - k2 - n2 * K[2] < n){
	      Mn->value (u[2] - k2 - n2 * K[2], tmp2);
	    }
	    else{
	      continue;
	    }
	    tmp = tmp0 * tmp1 * tmp2;
	    force[0] += QConvPhi0[count] * tmp;
	    force[1] += QConvPhi1[count] * tmp;
	    force[2] += QConvPhi2[count] * tmp;
	  }
	}
      }
    }
    else{
      for (int k0 = int(ceil(posi0-n)); k0 < int(ceil(posi0)); ++ k0){
	Mn->value (posi0 - k0, tmp0);
	index0 = K[1] * (k0<0 ? k0+K[0] : k0);
	for (int k1 = int(ceil(posi1-n)); k1 < int(ceil(posi1)); ++ k1){
	  Mn->value (posi1 - k1, tmp1);
	  index1 = K[2] * ((k1<0 ? k1+K[1] : k1) + index0);
	  for (int k2 = int(ceil(posi2-n)); k2 < int(ceil(posi2)); ++ k2){
	    Mn->value (posi2 - k2, tmp2);
	    tmp = tmp0 * tmp1 * tmp2;
	    unsigned index = (k2<0 ? k2+K[2] : k2) + index1;
	    force[0] += QConvPhi0[index] * tmp;
	    force[1] += QConvPhi1[index] * tmp;
	    force[2] += QConvPhi2[index] * tmp;
	  }
	}
      }
    }

    (*ppart)->f()[0] += (*ppart)->charge() * force[0];
    (*ppart)->f()[1] += (*ppart)->charge() * force[1];
    (*ppart)->f()[2] += (*ppart)->charge() * force[2];
  }
  watch.stop();
  loop_time += watch.user();
}


value_type ElectrostaticInteraction_rec_BSpline::calPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool)
{
  watch.start();
  calQ(beginPool, endPool);
  watch.stop();
  loop_time += watch.user();
  
  watch.start();
  fftw_execute (forwardQ);
  for (unsigned i = 0; i < K[0]*K[1]*(K[2]/2+1); i ++){
    QFProdPsiF[i][0] = QF[i][0] * psiF[i][0] - QF[i][1] * psiF[i][1];
    QFProdPsiF[i][1] = QF[i][0] * psiF[i][1] + QF[i][1] * psiF[i][0];
  }
  fftw_execute (backwardQFProdPsiF);
  watch.stop();
  conv_time += watch.user();

  int size = K[0]*K[1]*K[2];
  value_type sizei = 1./size;
  value_type value = 0;
  for (int i = 0; i < size; i ++){
    value += Q[i] * QConvPsi [i] * sizei;
  }
  value *= 0.5;

  return value;
}


value_type ElectrostaticInteraction_rec_BSpline::applyInteractionCalPotential (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
  // cal potential part
  // b[0];
  watch.start();
  calQ(beginPool, endPool);
  watch.stop();
  loop_time += watch.user();

  watch.start();
  fftw_execute (forwardQ);
  for (unsigned i = 0; i < K[0]*K[1]*(K[2]/2+1); i ++){
    QFProdPsiF[i][0] = QF[i][0] * psiF[i][0] - QF[i][1] * psiF[i][1];
    QFProdPsiF[i][1] = QF[i][0] * psiF[i][1] + QF[i][1] * psiF[i][0];
  }
  fftw_execute (backwardQFProdPsiF);
  watch.stop();
  conv_time += watch.user();

  int size_cal = K[0]*K[1]*K[2];
  value_type size_cali = 1./size_cal;
  value_type value = 0;
  for (int i = 0; i < size_cal; i ++){
    value += Q[i] * QConvPsi [i] * size_cali;
  }
  value *= 0.5;

  // apply interaction part
  unsigned n = Mn->getN();
  unsigned size = K[0] * K[1] * K[2];
  unsigned sizeHalf = K[0] * K[1] * (K[2]/2+1);
  value_type sizei = 1./size;
  
  watch.start();
  for (unsigned i = 0; i < sizeHalf; i ++){
    QFProdPhiF0[i][0] = QF[i][0] * phiF0[i][0] - QF[i][1] * phiF0[i][1];
    QFProdPhiF0[i][1] = QF[i][0] * phiF0[i][1] + QF[i][1] * phiF0[i][0];
    QFProdPhiF1[i][0] = QF[i][0] * phiF1[i][0] - QF[i][1] * phiF1[i][1];
    QFProdPhiF1[i][1] = QF[i][0] * phiF1[i][1] + QF[i][1] * phiF1[i][0];
    QFProdPhiF2[i][0] = QF[i][0] * phiF2[i][0] - QF[i][1] * phiF2[i][1];
    QFProdPhiF2[i][1] = QF[i][0] * phiF2[i][1] + QF[i][1] * phiF2[i][0];
  }
  fftw_execute (backwardQFProdPhiF0);
  fftw_execute (backwardQFProdPhiF1);
  fftw_execute (backwardQFProdPhiF2);
  
  for (unsigned i = 0; i < size; i ++){
    QConvPhi0[i] *= sizei;
    QConvPhi1[i] *= sizei;
    QConvPhi2[i] *= sizei;
  }
  watch.stop();
  conv_time += watch.user();

  watch.start();
  double ii0 = 1./ value_type(K[0]);
  double ii1 = 1./ value_type(K[1]);
  double ii2 = 1./ value_type(K[2]);
  // bool fast = ((n < K[0]) && (n < K[1]) && (n < K[2]));
  
  for (std::vector<StandardParticle * >::iterator ppart = beginPool;
       ppart != endPool; ppart ++){
    std::vector<double > force (3, 0);
    std::vector<value_type > u(3);
    u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*ppart)->r());
    u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*ppart)->r());
    u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*ppart)->r());
    int A0 = -int(floor ((u[0]) * ii0)) ;
    int A1 = -int(floor ((u[1]) * ii1)) ;
    int A2 = -int(floor ((u[2]) * ii2)) ;
    value_type posi0 = u[0] + A0 * K[0];
    value_type posi1 = u[1] + A1 * K[1];
    value_type posi2 = u[2] + A2 * K[2];
    std::vector<int > start (3, 0);
    start[0] = int(posi0) - n + 1;
    start[1] = int(posi1) - n + 1;
    start[2] = int(posi2) - n + 1;
    value_type tmp0, tmp1, tmp2;
    int index0, index1, index2;
    for (int kk0 = start[0]; kk0 < start[0] + int(n); ++kk0){
      kk0 < 0 ? (index0 = kk0 + K[0]) : (index0 = kk0);
      Mn->value (posi0 - kk0, tmp0);
      for (int kk1 = start[1]; kk1 < start[1] + int(n); ++kk1){
	kk1 < 0 ? (index1 = kk1 + K[1]) : (index1 = kk1);
	Mn->value (posi1 - kk1, tmp1);
	for (int kk2 = start[2]; kk2 < start[2] + int(n); ++kk2){
	  kk2 < 0 ? (index2 = kk2 + K[2]) : (index2 = kk2);
	  Mn->value (posi2 - kk2, tmp2);
	  int index = index2 + K[2] * (index1 + K[1] * index0);
	  double tmp = tmp0 * tmp1 * tmp2;
      	  force[0] += QConvPhi0[index] * tmp;
	  force[1] += QConvPhi1[index] * tmp;
	  force[2] += QConvPhi2[index] * tmp;
	}
      }
    }
    (*ppart)->f()[0] += (*ppart)->charge() * force[0];
    (*ppart)->f()[1] += (*ppart)->charge() * force[1];
    (*ppart)->f()[2] += (*ppart)->charge() * force[2];
  }


  {
    for (std::vector<StandardParticle * >::iterator ppart = beginPool;
	 ppart != endPool; ++(++(++ppart))){
      std::vector<std::vector<value_type > > posi (3);
      std::vector<std::vector<int > > start (3);
      std::vector<double > atomcharge(3);
      std::vector<StandardParticle * >::iterator patom = ppart;
      for (unsigned molc = 0; molc < 3; ++molc, ++patom){
	atomcharge[molc] = (*patom)->charge();
	std::vector<value_type > u(3);
	u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*patom)->r());
	u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*patom)->r());
	u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*patom)->r());
	std::vector<value_type > tmpposi(3);
	tmpposi[0] = u[0] -int(floor ((u[0]) * ii0)) * K[0];
	tmpposi[1] = u[1] -int(floor ((u[1]) * ii1)) * K[1];
	tmpposi[2] = u[2] -int(floor ((u[2]) * ii2)) * K[2];
	posi[molc] = tmpposi;
	start[molc].resize (3);
	for (unsigned dd = 0; dd < 3; ++dd) start[molc][dd] = int (posi[molc][dd]) - n + 1;
      }

      // forceCorr and initialization
      std::vector<std::vector<std::vector<double > > > forceCorr (3);
      {
	for (unsigned dd = 0; dd < 3; ++dd) forceCorr[dd].resize (3);
	std::vector<double > tmp (3, 0.);
	for (unsigned ii = 0; ii < 3; ++ii){
	  for (unsigned jj = 0; jj < 3; ++jj){
	    forceCorr[ii][jj] = tmp;
	  }
	}
      }
      
      for (int mol0 = 0; mol0 < 3; ++mol0){
      for (int mol1 = mol0 + 1; mol1 < 3; ++mol1){
	std::vector<std::vector<double > > mn_mol0(3);
	std::vector<std::vector<double > > mn_mol1(3);
	for (unsigned dd = 0; dd < 3; ++dd){
	  mn_mol0[dd].resize(n, 0.);
	  mn_mol1[dd].resize(n, 0.);
	}	
	for (unsigned dd = 0; dd < 3; ++dd){
	  for (unsigned kk = start[mol0][dd]; kk < start[mol0][dd] + n; ++kk){
	    Mn->value (posi[mol0][dd] - kk, mn_mol0[dd][kk - start[mol0][dd]]);
	  }
	}
	for (unsigned dd = 0; dd < 3; ++dd){
	  for (unsigned kk = start[mol1][dd]; kk < start[mol1][dd] + n; ++kk){
	    Mn->value (posi[mol1][dd] - kk, mn_mol1[dd][kk - start[mol1][dd]]);
	  }
	}	
	  
	std::vector<int > kk(3);
	std::vector<int > ll(3);
	std::vector<int > diffkl (3);
	for (kk[0] = start[mol0][0]; kk[0] < start[mol0][0] + int(n); ++kk[0]){
	for (ll[0] = start[mol1][0]; ll[0] < start[mol1][0] + int(n); ++ll[0]){
	  diffkl[0] = kk[0] - ll[0];
	  if      (diffkl[0] <  0        ) diffkl[0] += int(K[0]);
	  else if (diffkl[0] >= int(K[0])) diffkl[0] -= int(K[0]);
	for (kk[1] = start[mol0][1]; kk[1] < start[mol0][1] + int(n); ++kk[1]){
	for (ll[1] = start[mol1][1]; ll[1] < start[mol1][1] + int(n); ++ll[1]){
	  diffkl[1] = kk[1] - ll[1];
	  if      (diffkl[1] <  0        ) diffkl[1] += int(K[1]);
	  else if (diffkl[1] >= int(K[1])) diffkl[1] -= int(K[1]);
	for (kk[2] = start[mol0][2]; kk[2] < start[mol0][2] + int(n); ++kk[2]){
	for (ll[2] = start[mol1][2]; ll[2] < start[mol1][2] + int(n); ++ll[2]){
	  diffkl[2] = kk[2] - ll[2];
	  if      (diffkl[2] <  0        ) diffkl[2] += int(K[2]);
	  else if (diffkl[2] >= int(K[2])) diffkl[2] -= int(K[2]);
	  double scalor =
	      atomcharge[mol0] * atomcharge[mol1] *
	      mn_mol0[0][kk[0] - start[mol0][0]] *
	      mn_mol0[1][kk[1] - start[mol0][1]] *
	      mn_mol0[2][kk[2] - start[mol0][2]] *
	      mn_mol1[0][ll[0] - start[mol1][0]] *
	      mn_mol1[1][ll[1] - start[mol1][1]] *
	      mn_mol1[2][ll[2] - start[mol1][2]] ;
	  forceCorr[mol0][mol1][0] += scalor *
	      phiFr0[diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0];
	  forceCorr[mol0][mol1][1] += scalor *
	      phiFr1[diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0];
	  forceCorr[mol0][mol1][2] += scalor *
	      phiFr2[diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0];
	  // std::cout << scalor << std::endl;
	  // std::cout << diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0]) << std::endl;
	  // std::cout << phiFr1[diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0] << std::endl;
	}
	}
	}
	}
	}
	}
      }
      }
      std::cout << forceCorr[0][1][0] <<  "  "
		<< forceCorr[0][1][1] <<  "  "
		<< forceCorr[0][1][2] <<  "  "
		<< std::endl;
      patom = ppart;
      (*patom)->f()[0] -= forceCorr[0][1][0];
      (*patom)->f()[1] -= forceCorr[0][1][1];
      (*patom)->f()[2] -= forceCorr[0][1][2];
      (*patom)->f()[0] -= forceCorr[0][2][0];
      (*patom)->f()[1] -= forceCorr[0][2][1];
      (*patom)->f()[2] -= forceCorr[0][2][2];
      ++patom;
      (*patom)->f()[0] -= forceCorr[1][2][0];
      (*patom)->f()[1] -= forceCorr[1][2][1];
      (*patom)->f()[2] -= forceCorr[1][2][2];
      (*patom)->f()[0] += forceCorr[0][1][0];
      (*patom)->f()[1] += forceCorr[0][1][1];
      (*patom)->f()[2] += forceCorr[0][1][2];
      ++patom;
      (*patom)->f()[0] += forceCorr[1][2][0];
      (*patom)->f()[1] += forceCorr[1][2][1];
      (*patom)->f()[2] += forceCorr[1][2][2];
      (*patom)->f()[0] += forceCorr[0][2][0];
      (*patom)->f()[1] += forceCorr[0][2][1];
      (*patom)->f()[2] += forceCorr[0][2][2];
    } // end for
    std::cout << "hit" << std::endl;
    
    //   std::vector<double > force (3, 0);
    //   std::vector<value_type > u(3);
    //   u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*ppart)->r());
    //   u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*ppart)->r());
    //   u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*ppart)->r());
    //   int A0 = -int(floor ((u[0]) * ii0)) ;
    //   int A1 = -int(floor ((u[1]) * ii1)) ;
    //   int A2 = -int(floor ((u[2]) * ii2)) ;
    //   value_type posi0 = u[0] + A0 * K[0];
    //   value_type posi1 = u[1] + A1 * K[1];
    //   value_type posi2 = u[2] + A2 * K[2];
    //   std::vector<int > start (3, 0);
    //   start[0] = int(posi0) - n + 1;
    //   start[1] = int(posi1) - n + 1;
    //   start[2] = int(posi2) - n + 1;
    //   value_type tmp0, tmp1, tmp2;
    //   int index0, index1, index2;
    //   for (int kk0 = start[0]; kk0 < start[0] + int(n); ++kk0){
    // 	kk0 < 0 ? (index0 = kk0 + K[0]) : (index0 = kk0);
    // 	Mn->value (posi0 - kk0, tmp0);
    // 	for (int kk1 = start[1]; kk1 < start[1] + int(n); ++kk1){
    // 	  kk1 < 0 ? (index1 = kk1 + K[1]) : (index1 = kk1);
    // 	  Mn->value (posi1 - kk1, tmp1);
    // 	  for (int kk2 = start[2]; kk2 < start[2] + int(n); ++kk2){
    // 	    kk2 < 0 ? (index2 = kk2 + K[2]) : (index2 = kk2);
    // 	    Mn->value (posi2 - kk2, tmp2);
    // 	    int index = index2 + K[2] * (index1 + K[1] * index0);
    // 	    double tmp = tmp0 * tmp1 * tmp2;
    // 	    force[0] += QConvPhi0[index] * tmp;
    // 	    force[1] += QConvPhi1[index] * tmp;
    // 	    force[2] += QConvPhi2[index] * tmp;
    // 	  }
    // 	}
    //   }
    //   (*ppart)->f()[0] += (*ppart)->charge() * force[0];
    //   (*ppart)->f()[1] += (*ppart)->charge() * force[1];
    //   (*ppart)->f()[2] += (*ppart)->charge() * force[2];
    // }
      
  } // end out brace


  watch.stop();
  loop_time += watch.user();

  return value;
}
  



void ElectrostaticInteraction_rec_BSpline::calQ (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
  //   std::cerr << "K0  " << K[0] << std::endl;
  
  unsigned n = Mn->getN();

  for (unsigned i = 0; i < K[0] * K[1] * K[2]; i ++){
    Q[i] = 0;
  }
  double ii0 = 1./ value_type(K[0]);
  double ii1 = 1./ value_type(K[1]);
  double ii2 = 1./ value_type(K[2]);
  bool fast = ((n < K[0]) && (n < K[1]) && (n < K[2]));
  for (std::vector<StandardParticle * >::iterator ppart = beginPool;
       ppart != endPool; ppart ++){
    std::vector<value_type > u(3);
    u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*ppart)->r());
    u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*ppart)->r());
    u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*ppart)->r());
    int A0 = -int(floor ((u[0]) * ii0)) ;
    int A1 = -int(floor ((u[1]) * ii1)) ;
    int A2 = -int(floor ((u[2]) * ii2)) ;
    value_type posi0 = u[0] + A0 * K[0];
    value_type posi1 = u[1] + A1 * K[1];
    value_type posi2 = u[2] + A2 * K[2];
    value_type tmp0 = 0;
    value_type tmp1 = 0;
    value_type tmp2 = 0;

    int count = 0;
    if (!fast){
      for (unsigned k0 = 0; k0 < K[0]; k0 ++){
	for (unsigned k1 = 0; k1 < K[1]; k1 ++){
	  for (unsigned k2 = 0; k2 < K[2]; k2 ++, count ++){
	    value_type n0 = floor ((u[0] - k0) * ii0);
	    if (u[0] - k0 - n0 * K[0] < n){
	      Mn->value (u[0] - k0 - n0 * K[0], tmp0);
	    }
	    else{
	      continue;
	    }
	    value_type n1 = floor ((u[1] - k1) * ii1);
	    if (u[1] - k1 - n1 * K[1] < n){
	      Mn->value (u[1] - k1 - n1 * K[1], tmp1);
	    }
	    else{
	      continue;
	    }
	    value_type n2 = floor ((u[2] - k2) * ii2);
	    if (u[2] - k2 - n2 * K[2] < n){
	      Mn->value (u[2] - k2 - n2 * K[2], tmp2);
	    }
	    else{
	      continue;
	    }
	    Q[count] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
	  }
	}
      }
    }
    else {
      unsigned index0, index1;
      for (int k0 = int(ceil(posi0-n)); k0 < int(ceil(posi0)); ++ k0){
	Mn->value (posi0 - k0, tmp0);
	index0 = K[1] * (k0<0 ? k0+K[0] : k0);
	for (int k1 = int(ceil(posi1-n)); k1 < int(ceil(posi1)); ++ k1){
	  Mn->value (posi1 - k1, tmp1);
	  index1 = K[2] * ((k1<0 ? k1+K[1] : k1) + index0);
	  for (int k2 = int(ceil(posi2-n)); k2 < int(ceil(posi2)); ++ k2){
	    Mn->value (posi2 - k2, tmp2);
	    Q[(k2<0 ? k2+K[2] : k2) + index1] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
	  }
	}
      }
    }
  }
}




void ElectrostaticInteraction_rec_BSpline::init (const std::vector<std::vector<value_type > > &vecA_, 
						 const std::vector<unsigned > K_,
						 const value_type & beta_,
						 const InterpolationInfo::Order & interpOrder_)
{  
  switch (interpOrder_){
  case 2 :
      Mn = new BSpline2();
      break;
  case 4 :
      Mn = new BSpline4();
      break;
//   case 5 :
//       Mn = new BSpline5();
//       break;
  case 6 :
      Mn = new BSpline6();
      break;
  case 8 :
      Mn = new BSpline8();
      break;
  case 10 :
      Mn = new BSpline10();
      break;
  case 12 :
      Mn = new BSpline12();
      break;
  case 14 :
      Mn = new BSpline14();
      break;
  default:
      std::cerr << "no such order implemented, use order 8" << std::endl;
      Mn = new BSpline8();
  }
  K = K_;
  beta = beta_;
  vecA = vecA_;
  calV();
  calAStar();
  calB ();
  loop_time = 0;
  conv_time = 0;
  ecal.reinit (K, interpOrder_);
  unBuild = true;
}

void ElectrostaticInteraction_rec_BSpline::build()
{
  int size = K[0] * K[1] * K[2];
  int sizeHalf = K[0] * K[1] * (K[2] / 2 + 1);
  Q	= (value_type *) fftw_malloc (sizeof(value_type) * size);
  psiF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  phiF0	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  phiF1	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  phiF2	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  phiFr0= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  phiFr1= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  phiFr2= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  QF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QFProdPsiF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QFProdPhiF0	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QFProdPhiF1	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QFProdPhiF2	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QConvPsi	= (value_type *) fftw_malloc (sizeof(value_type) * size);
  QConvPhi0	= (value_type *) fftw_malloc (sizeof(value_type) * size);
  QConvPhi1	= (value_type *) fftw_malloc (sizeof(value_type) * size);
  QConvPhi2	= (value_type *) fftw_malloc (sizeof(value_type) * size);

  forwardQ	= fftw_plan_dft_r2c_3d (K[0], K[1], K[2], Q  , QF  , FFTW_MEASURE);
  backwardQFProdPsiF = fftw_plan_dft_c2r_3d (K[0], K[1], K[2], QFProdPsiF, QConvPsi, FFTW_MEASURE);
  backwardQFProdPhiF0 = fftw_plan_dft_c2r_3d (K[0], K[1], K[2], QFProdPhiF0, QConvPhi0, FFTW_MEASURE);
  backwardQFProdPhiF1 = fftw_plan_dft_c2r_3d (K[0], K[1], K[2], QFProdPhiF1, QConvPhi1, FFTW_MEASURE);
  backwardQFProdPhiF2 = fftw_plan_dft_c2r_3d (K[0], K[1], K[2], QFProdPhiF2, QConvPhi2, FFTW_MEASURE);

  backward_phiF0 = fftw_plan_dft_3d (K[0], K[1], K[2], phiFr0, phiFr0, +1, FFTW_MEASURE);
  backward_phiF1 = fftw_plan_dft_3d (K[0], K[1], K[2], phiFr1, phiFr1, +1, FFTW_MEASURE);
  backward_phiF2 = fftw_plan_dft_3d (K[0], K[1], K[2], phiFr2, phiFr2, +1, FFTW_MEASURE);
  
  calPsiFPhiF();
  unBuild = false;
}



void ElectrostaticInteraction_rec_BSpline::clear()
{
  delete Mn;
  if (!unBuild){
    fftw_free (Q);
    fftw_free (QF);
    fftw_free (psiF);
    fftw_free (phiF0);
    fftw_free (phiF1);
    fftw_free (phiF2);
    fftw_free (QFProdPsiF);
    fftw_free (QFProdPhiF0);
    fftw_free (QFProdPhiF1);
    fftw_free (QFProdPhiF2);
    fftw_free (QConvPsi);
    fftw_free (QConvPhi0);
    fftw_free (QConvPhi1);
    fftw_free (QConvPhi2);
  
    fftw_destroy_plan (forwardQ);
    fftw_destroy_plan (backwardQFProdPsiF);
    fftw_destroy_plan (backwardQFProdPhiF0);
    fftw_destroy_plan (backwardQFProdPhiF1);
    fftw_destroy_plan (backwardQFProdPhiF2);
  }
}


void ElectrostaticInteraction_rec_BSpline::calB ()
{
  unsigned n = Mn->getN();
  b0.resize (K[0]);
  b1.resize (K[1]);
  b2.resize (K[2]);

  for (unsigned m = 0; m < K[0]; ++m){
    std::vector<value_type > fenzi (2);
    value_type tmp = 2 * M_PI * (n-1) * m / value_type (K[0]);
    fenzi[0] = cos(tmp);
    fenzi[1] = sin(tmp);
    std::vector<value_type > fenmu (2, 0);
    for (unsigned k = 0; k < n-1; k ++){
      value_type scale ;
      Mn->value (k+1, scale);
      tmp = 2 * M_PI * m * k / value_type (K[0]);
      fenmu[0] += scale * cos(tmp);
      fenmu[1] += scale * sin(tmp);
    }
    std::vector<value_type > btmp (2);
    value_type scale = 1./ (fenmu[0]*fenmu[0] + fenmu[1]*fenmu[1]);
    btmp[0] = scale * (fenzi[0] * fenmu[0] + fenzi[1] * fenmu[1]);
    btmp[1] = scale * (fenzi[1] * fenmu[0] - fenzi[0] * fenmu[1]);
    b0[m] = btmp[0]*btmp[0] + btmp[1]*btmp[1];
  }
  
  for (unsigned m = 0; m < K[1]; ++m){
    std::vector<value_type > fenzi (2);
    value_type tmp = 2 * M_PI * (n-1) * m / value_type (K[1]);
    fenzi[0] = cos(tmp);
    fenzi[1] = sin(tmp);
    std::vector<value_type > fenmu (2, 0);
    for (unsigned k = 0; k < n-1; k ++){
      value_type scale ;
      Mn->value (k+1, scale);
      tmp = 2 * M_PI * m * k / value_type (K[1]);
      fenmu[0] += scale * cos(tmp);
      fenmu[1] += scale * sin(tmp);
    }
    std::vector<value_type > btmp (2);
    value_type scale = 1./ (fenmu[0]*fenmu[0] + fenmu[1]*fenmu[1]);
    btmp[0] = scale * (fenzi[0] * fenmu[0] + fenzi[1] * fenmu[1]);
    btmp[1] = scale * (fenzi[1] * fenmu[0] - fenzi[0] * fenmu[1]);
    b1[m] = btmp[0]*btmp[0] + btmp[1]*btmp[1];
  }
  
  for (unsigned m = 0; m < K[2]; ++m){
    std::vector<value_type > fenzi (2);
    value_type tmp = 2 * M_PI * (n-1) * m / value_type (K[2]);
    fenzi[0] = cos(tmp);
    fenzi[1] = sin(tmp);
    std::vector<value_type > fenmu (2, 0);
    for (unsigned k = 0; k < n-1; k ++){
      value_type scale ;
      Mn->value (k+1, scale);
      tmp = 2 * M_PI * m * k / value_type (K[2]);
      fenmu[0] += scale * cos(tmp);
      fenmu[1] += scale * sin(tmp);
    }
    std::vector<value_type > btmp (2); 
    value_type scale = 1./ (fenmu[0]*fenmu[0] + fenmu[1]*fenmu[1]);
    btmp[0] = scale * (fenzi[0] * fenmu[0] + fenzi[1] * fenmu[1]);
    btmp[1] = scale * (fenzi[1] * fenmu[0] - fenzi[0] * fenmu[1]);
    b2[m] = btmp[0]*btmp[0] + btmp[1]*btmp[1];
  }
  
}


void ElectrostaticInteraction_rec_BSpline::calPsiFPhiF ()
{
  fftw_complex * psiFtmp;
  fftw_complex * phiFtmp0;
  fftw_complex * phiFtmp1;
  fftw_complex * phiFtmp2;
  
  psiFtmp= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * K[0]*K[1]*K[2]);
  phiFtmp0=(fftw_complex *) fftw_malloc (sizeof(fftw_complex) * K[0]*K[1]*K[2]);
  phiFtmp1=(fftw_complex *) fftw_malloc (sizeof(fftw_complex) * K[0]*K[1]*K[2]);
  phiFtmp2=(fftw_complex *) fftw_malloc (sizeof(fftw_complex) * K[0]*K[1]*K[2]);

  value_type oneOverPiV = 1. / (M_PI * V);
  value_type scale = M_PI * M_PI / beta / beta;
  value_type minousTwoOverV = -2./V;
  unsigned size = K[0]*K[1]*K[2];
  
  psiFtmp[0][0] = 0;
  psiFtmp[0][1] = 0;
  phiFtmp0[0][0] = 0;
  phiFtmp0[0][1] = 0;
  phiFtmp1[0][0] = 0;
  phiFtmp1[0][1] = 0;
  phiFtmp2[0][0] = 0;
  phiFtmp2[0][1] = 0;
  phiFr0[0][0] = phiFr0[0][1] = 0.;
  phiFr1[0][0] = phiFr1[0][1] = 0.;
  phiFr2[0][0] = phiFr2[0][1] = 0.;
  std::vector<value_type > m (3);
  for (unsigned i = 0; i < K[0]; i ++){
    int ip ;    
    if (i <= K[0] / 2) ip = i;
    else ip = i - K[0];
    for (unsigned j = 0; j < K[1]; j ++){
      int jp;      
      if (j <= K[1] / 2) jp = j;
      else jp = j - K[1];
      for (unsigned k = 0; k < K[2]; k ++){
	int kp ;        
	if (k <= K[2] / 2) kp = k;
	else kp = k - K[2];

	if (kp == 0 && jp == 0 && ip == 0) continue;
	m[0] = ip * vecAStar[0][0] + jp * vecAStar[1][0] + kp * vecAStar[2][0];
	m[1] = ip * vecAStar[0][1] + jp * vecAStar[1][1] + kp * vecAStar[2][1];
	m[2] = ip * vecAStar[0][2] + jp * vecAStar[1][2] + kp * vecAStar[2][2];
	value_type mm = VectorOperation::dot (m, m);
	value_type expm = exp (-scale * mm) / mm * b0[i] * b1[j] * b2[k];
	psiFtmp[k + K[2] * (j + K[1] * i)][0] = oneOverPiV * expm * size;
	psiFtmp[k + K[2] * (j + K[1] * i)][1] = 0;
	phiFtmp0[k + K[2] * (j + K[1] * i)][0] = 0;
	phiFtmp0[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[0] * size;
	phiFtmp1[k + K[2] * (j + K[1] * i)][0] = 0;
	phiFtmp1[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[1] * size;
	phiFtmp2[k + K[2] * (j + K[1] * i)][0] = 0;
	phiFtmp2[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[2] * size;
	phiFr0[k + K[2] * (j + K[1] * i)][0] = 0.;
	phiFr1[k + K[2] * (j + K[1] * i)][0] = 0.;
	phiFr2[k + K[2] * (j + K[1] * i)][0] = 0.;
	phiFr0[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[0];
	phiFr1[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[1];
	phiFr2[k + K[2] * (j + K[1] * i)][1] = minousTwoOverV * expm * m[2];
	// std::cout << phiFr1[k + K[2] * (j + K[1] * i)][1] << std::endl;
      }
    }
  }

  // for (unsigned i = 0; i < K[0]; i ++){
  //   for (unsigned j = 0; j < K[1]; j ++){
  //     for (unsigned k = 0; k < K[2]; k ++){
  // 	std::cout << phiFr1[k + K[2] * (j + K[1] * i)][1] << std::endl;
  //     }
  //   }
  // }
	
  fftw_execute (backward_phiF0);
  fftw_execute (backward_phiF1);
  fftw_execute (backward_phiF2);

  // for (unsigned i = 0; i < K[0]; i ++){
  //   for (unsigned j = 0; j < K[1]; j ++){
  //     for (unsigned k = 0; k < K[2]; k ++){
  // 	std::cout << phiFr1[k + K[2] * (j + K[1] * i)][1] << std::endl;
  //     }
  //   }
  // }

  // std::cout << phiFr1[100][0] << std::endl;
  // std::cout << phiFr1[100][1] << std::endl;

  for (unsigned i = 0; i < K[0]; ++ i){
    for (unsigned j = 0; j < K[1]; ++ j){
      for (unsigned k = 0; k < K[2]/2+1; ++ k){
	psiF[k+(K[2]/2+1)*(j+K[1]*i)][0] = psiFtmp[k+K[2]*(j+K[1]*i)][0];
	psiF[k+(K[2]/2+1)*(j+K[1]*i)][1] = psiFtmp[k+K[2]*(j+K[1]*i)][1];
	phiF0[k+(K[2]/2+1)*(j+K[1]*i)][0] = phiFtmp0[k+K[2]*(j+K[1]*i)][0];
	phiF0[k+(K[2]/2+1)*(j+K[1]*i)][1] = phiFtmp0[k+K[2]*(j+K[1]*i)][1];
	phiF1[k+(K[2]/2+1)*(j+K[1]*i)][0] = phiFtmp1[k+K[2]*(j+K[1]*i)][0];
	phiF1[k+(K[2]/2+1)*(j+K[1]*i)][1] = phiFtmp1[k+K[2]*(j+K[1]*i)][1];
	phiF2[k+(K[2]/2+1)*(j+K[1]*i)][0] = phiFtmp2[k+K[2]*(j+K[1]*i)][0];
	phiF2[k+(K[2]/2+1)*(j+K[1]*i)][1] = phiFtmp2[k+K[2]*(j+K[1]*i)][1];
      }
    }
  }

  fftw_free (psiFtmp);
  fftw_free (phiFtmp0);
  fftw_free (phiFtmp1);
  fftw_free (phiFtmp2);
  //   std::cout << "beta  " << beta << "  V  " << V << std::endl;
}




void ElectrostaticInteraction_rec_BSpline::calV()
{
  V = vecA[0][0] * (vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) - 
      vecA[0][1] * (vecA[1][0]*vecA[2][2] - vecA[2][0]*vecA[1][2]) +
      vecA[0][2] * (vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]);
}
  
void ElectrostaticInteraction_rec_BSpline::calAStar ()
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









//       if (posi0 < n){
// 	for (unsigned k0 = 0; k0 < unsigned (ceil(posi0)); k0 ++){
// 	  Mn->value (posi0 - k0, tmp0);
// 	  index0 = K[1]*k0;
// 	  if (posi1 < n){
// 	    for (unsigned k1 = 0; k1 < unsigned (ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	    for (unsigned k1 = unsigned(ceil(posi1+K[1]-n)); k1 < K[1]; k1 ++){
// 	      Mn->value (u[1] - k1 + (A1+1) * K[1], tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	  else{
// 	    for (unsigned k1 = unsigned(ceil(posi1-n)); k1 < unsigned(ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	for (unsigned k0 = unsigned(ceil(posi0+K[0]-n)); k0 < K[0]; k0 ++){
// 	  Mn->value (u[0] - k0 + (A0+1) * K[0], tmp0);
// 	  index0 = K[1]*k0;
// 	  if (posi1 < n){
// 	    for (unsigned k1 = 0; k1 < unsigned (ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	    for (unsigned k1 = unsigned(ceil(posi1+K[1]-n)); k1 < K[1]; k1 ++){
// 	      Mn->value (u[1] - k1 + (A1+1) * K[1], tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	  else{
// 	    for (unsigned k1 = unsigned(ceil(posi1-n)); k1 < unsigned(ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//       else{
// 	for (unsigned k0 = unsigned(ceil(posi0-n)); k0 < unsigned(ceil(posi0)); k0 ++){
// 	  Mn->value (posi0 - k0, tmp0);
// 	  index0 = K[1]*k0;
// 	  if (posi1 < n){
// 	    for (unsigned k1 = 0; k1 < unsigned (ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	    for (unsigned k1 = unsigned(ceil(posi1+K[1]-n)); k1 < K[1]; k1 ++){
// 	      Mn->value (u[1] - k1 + (A1+1) * K[1], tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	  else{
// 	    for (unsigned k1 = unsigned(ceil(posi1-n)); k1 < unsigned(ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);
// 	      index1 = K[2]*(k1+index0);
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// // 	  unsigned index = k2+K[2]*(k1+K[1]*k0);
// 		  unsigned index = k2 + index1;
// 		  force[0] += QConvPhi0[index] * tmp0 * tmp1 * tmp2;
// 		  force[1] += QConvPhi1[index] * tmp0 * tmp1 * tmp2;
// 		  force[2] += QConvPhi2[index] * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }



// // 	std::cout << "bububu" << std::endl;
// 	for (unsigned k0 = unsigned(ceil(posi0-n)); k0 < unsigned(ceil(posi0)); k0 ++){
// 	  Mn->value (posi0 - k0, tmp0);
// 	  if (posi1 < n){
// 	    for (unsigned k1 = 0; k1 < unsigned (ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);	      
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	    for (unsigned k1 = unsigned(ceil(posi1+K[1]-n)); k1 < K[1]; k1 ++){
// 	      Mn->value (u[1] - k1 + (A1+1) * K[1], tmp1);	      
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	  else{
// 	    for (unsigned k1 = unsigned(ceil(posi1-n)); k1 < unsigned(ceil(posi1)); k1 ++){
// 	      Mn->value (posi1 - k1, tmp1);	      
// 	      if (posi2 < n){
// 		for (unsigned k2 = 0; k2 < unsigned (ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 		for (unsigned k2 = unsigned(ceil(posi2+K[2]-n)); k2 < K[2]; k2 ++){
// 		  Mn->value (u[2] - k2 + (A2+1) * K[2], tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	      else{
// 		for (unsigned k2 = unsigned(ceil(posi2-n)); k2 < unsigned(ceil(posi2)); k2 ++){
// 		  Mn->value (posi2 - k2, tmp2);
// 		  Q[k2+K[2]*(k1+K[1]*k0)] += (*ppart)->charge() * tmp0 * tmp1 * tmp2;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
      



