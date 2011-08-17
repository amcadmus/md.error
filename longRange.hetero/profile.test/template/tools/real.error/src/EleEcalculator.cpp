#include "ElectrostaticInteraction.h"
#include <cmath>

double intPow (const double & x, const unsigned & alpha_);

//template <int N>
inline 
double intPow (int N, const double & x)
{
  double x2, x4;
  switch (N){
  case 1:
      return x;
      break;
  case 2:
      return x*x;
      break;
  case 3:
      return x*x*x;
      break;
  case 4:
      x2 = x*x;
      return x2*x2;
      break;
  case 5:
      x2 = x*x;
      return x2 * x2 * x;
      break;
  case 6:
      x2 = x*x;
      return x2 * x2 * x2;
      break;
  case 7:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x2 * x;
      break;
  case 8:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x4;
      break;
  case 9:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x4 * x;
      break;
  case 10:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x4 * x2;
      break;
  case 11:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x4 * x2 * x;
      break;
  case 12:
      x2 = x*x;
      x4 = x2 * x2;
      return x4 * x4 * x4;
      break;
  default:
      return intPow (x, unsigned(N));
  }
}

inline double sixTerms (const double & x, const int & n)
{
  double sum = 0;
  sum += 1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI);
  sum += 1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI);
  sum += 1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI);
  sum += 1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI);
  sum += 1./intPow(n, x + 10*M_PI) + 1./intPow(n, x - 10*M_PI);
  return sum;
}
inline double sixTermsp (const double & x, const int & n)
{
  double sum = 0;
  sum += (1./intPow(n, x + 2*M_PI) - 1./intPow(n, x - 2*M_PI));
  sum += (1./intPow(n, x + 4*M_PI) - 1./intPow(n, x - 4*M_PI)) * 2;
  sum += (1./intPow(n, x + 6*M_PI) - 1./intPow(n, x - 6*M_PI)) * 3;
  sum += (1./intPow(n, x + 8*M_PI) - 1./intPow(n, x - 8*M_PI)) * 4;
  sum += (1./intPow(n, x + 10*M_PI) - 1./intPow(n, x - 10*M_PI)) * 5;
  return sum ;
}
inline double sixTermspp (const double & x, const int & n)
{
//   std::cout << n << std::endl;
  double sum = 0;
  sum += (1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI));
  sum += (1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI)) * 4;
  sum += (1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI)) * 9;
  sum += (1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI)) * 16;
  sum += (1./intPow(n, x + 10*M_PI) + 1./intPow(n, x - 10*M_PI)) * 25;
  return sum ;
}

inline double sevenTerms (const double & x, const int & n)
{
  double sum = 1./intPow(n, x);
  sum += 1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI);
  sum += 1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI);
  sum += 1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI);
  sum += 1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI);
  sum += 1./intPow(n, x + 10*M_PI) + 1./intPow(n, x - 10*M_PI);
  sum += 1./intPow(n, x + 12*M_PI) + 1./intPow(n, x - 12*M_PI);
  sum += 1./intPow(n, x + 14*M_PI) + 1./intPow(n, x - 14*M_PI);
  return sum;
}

inline double denominator (const double & x, const int & n)
{
  return sevenTerms (x, n) ;
// + 
//       ( 1./ (2.*M_PI * (n-1) * intPow(x + 2*M_PI*(3-0.5), n-1)) -
// 	1./ (2.*M_PI * (n-1) * intPow(x - 2*M_PI*(3-0.5), n-1)) );
}

inline void produce4Item (const double & x, const int & n,
			  double & e1, double & e2, double & e2p, double &e2pp )
{
  if (x == 0){
    e1 = (e2 = (e2p = (e2pp = 0.)));
    return ;
  }
  double tmp6 = sixTerms (x, n);
  double de = denominator (x, n);
  e1 = - tmp6 / de;
  e2 = (tmp6 * tmp6 + sixTerms (x, n*2)) / (de * de);
  e2p = 2 * M_PI / (de*de) *  sixTermsp (x, 2*n);
  e2pp = 4 * M_PI * M_PI / (de*de) * sixTermspp (x, n*2);
}


void EleEcalculator::reinit (const std::vector<unsigned > & K_,
			     const unsigned & n_)
{
  K.resize(3, 0);
  K[0] = int(K_[0]);
  K[1] = int(K_[1]);
  K[2] = int(K_[2]);
  n = int(n_);
  sw.resize(3);
  sw[0] = K[0] / 2;
  sw[1] = K[1] / 2;
  sw[2] = K[2] / 2;
  
  e1K0.resize(K[0]+1);
  e2K0.resize(K[0]+1);
  e2pK0.resize(K[0]+1);
  e2ppK0.resize(K[0]+1);
  e1K1.resize(K[1]+1);
  e2K1.resize(K[1]+1);
  e2pK1.resize(K[1]+1);
  e2ppK1.resize(K[1]+1);
  e1K2.resize(K[2]+1);
  e2K2.resize(K[2]+1);
  e2pK2.resize(K[2]+1);
  e2ppK2.resize(K[2]+1);
  
  double u ;
  for (int m = -K[0]/2; m <= K[0]/2; ++m){
    u = 2*M_PI*m / double(K[0]);
    double a, b, c, d;
    produce4Item (u, n, a, b, c, d);
    e1K0[m+sw[0]] = a;
    e2K0[m+sw[0]] = b;
    e2pK0[m+sw[0]] = c;
    e2ppK0[m+sw[0]] = d;
  }
  
  if (K[1] == K[0]){
    e1K1 = e1K0;
    e2K1 = e2K0;
    e2pK1 = e2pK0;
    e2ppK1 = e2ppK0;
  }
  else {
    for (int m = -K[1]/2; m <= K[1]/2; ++m){
      u = 2*M_PI*m / double(K[1]);
      double a, b, c, d;
      produce4Item (u, n, a, b, c, d);
      e1K1[m+sw[1]] = a;
      e2K1[m+sw[1]] = b;
      e2pK1[m+sw[1]] = c;
      e2ppK1[m+sw[1]] = d;
    }
  }
  
  if (K[2] == K[0]){
    e1K2 = e1K0;
    e2K2 = e2K0;
    e2pK2 = e2pK0;
    e2ppK2 = e2ppK0;
  } else if (K[2] == K[1]){
    e1K2 = e1K1;
    e2K2 = e2K1;
    e2pK2 = e2pK1;
    e2ppK2 = e2ppK1;
  } else {
    for (int m = -K[2]/2; m <= K[2]/2; ++m){
      u = 2*M_PI*m / double(K[2]);
      double a, b, c, d;
      produce4Item (u, n, a, b, c, d);
      e1K2[m+sw[2]] = a;
      e2K2[m+sw[2]] = b;
      e2pK2[m+sw[2]] = c;
      e2ppK2[m+sw[2]] = d;
    }
  }
}


// void EleEcalculator::get_e1 (const int & m0, const int & m1, const int & m2,
// 			     double & e1K0_, double & e1K1_, double & e1K2_)
// {
//   e1K0_ = e1K0[m0 + sw[0]];
//   e1K1_ = e1K1[m1 + sw[1]];
//   e1K2_ = e1K2[m2 + sw[2]];
// }


// void EleEcalculator::get_e2 (const int & m0, const int & m1, const int & m2,
// 			     double & e2K0_, double & e2K1_, double & e2K2_)
// {
//   e2K0_ = e2K0[m0 + sw[0]];
//   e2K1_ = e2K1[m1 + sw[1]];
//   e2K2_ = e2K2[m2 + sw[2]];
// }


// void EleEcalculator::get_e2p (const int & m0, const int & m1, const int & m2,
// 			     double & e2pK0_, double & e2pK1_, double & e2pK2_)
// {
//   e2pK0_ = e2pK0[m0 + sw[0]];
//   e2pK1_ = e2pK1[m1 + sw[1]];
//   e2pK2_ = e2pK2[m2 + sw[2]];
// }


// void EleEcalculator::get_e2pp (const int & m0, const int & m1, const int & m2,
// 			       double & e2ppK0_, double & e2ppK1_, double & e2ppK2_)
// {
//   e2ppK0_ = e2ppK0[m0 + sw[0]];
//   e2ppK1_ = e2ppK1[m1 + sw[1]];
//   e2ppK2_ = e2ppK2[m2 + sw[2]];
// }


