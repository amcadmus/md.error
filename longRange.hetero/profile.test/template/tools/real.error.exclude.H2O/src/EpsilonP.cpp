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

static void epsilonp (double u, double x, double n, double & re, double & im)
{
  int A = 2;
  re = 0;
  im  = 0;
  if (u == 0) return ;
  double fenmu = 0;
  for (int k = -A; k <= A; ++k){
    double tmp  = 1./intPow(u + 2*M_PI*k, int(n));
    fenmu += tmp;
    re += - tmp * 2*M_PI*k * sin(2*M_PI*k*x);
//     im += tmp * 2*M_PI*k * cos(2*M_PI*k*x);
  }
  
  re /= fenmu;
  im /= fenmu;
}

void EpsilonP::reinit (const std::vector<unsigned > K_, 
		       const unsigned & n_,
		       const std::vector<unsigned > ngrid)
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
  
  ep0.resize (K[0]+1);
  for (std::vector<std::vector<double > >::iterator i = ep0.begin(); 
       i != ep0.end(); ++i){
    i->resize (ngrid[0]);
  }
  ep1.resize (K[1]+1);
  for (std::vector<std::vector<double > >::iterator i = ep1.begin(); 
       i != ep1.end(); ++i){
    i->resize (ngrid[1]);
  }
  ep2.resize (K[2]+1);
  for (std::vector<std::vector<double > >::iterator i = ep2.begin(); 
       i != ep2.end(); ++i){
    i->resize (ngrid[2]);
  }
  
  double u;
  double tmp;
  h.resize(3);
  h[0] = K[0] / double (ngrid[0]);
  h[1] = K[1] / double (ngrid[1]);
  h[2] = K[2] / double (ngrid[2]);
  hi.resize(3);
  hi[0] = 1./h[0];
  hi[1] = 1./h[1];
  hi[2] = 1./h[2];
  
  for (int m = -K[0]/2; m <= K[0]/2; ++m){
    u = 2*M_PI*m/double(K[0]);
    for (unsigned i = 0; i < ngrid[0]; ++ i){
      epsilonp (u, h[0]*i, n, ep0[m+sw[0]][i], tmp);
    }
  }
  
  for (int m = -K[1]/2; m <= K[1]/2; ++m){
    u = 2*M_PI*m/double(K[1]);
    for (unsigned i = 0; i < ngrid[1]; ++ i){
      epsilonp (u, h[1]*i, n, ep1[m+sw[1]][i], tmp);
    }
  }
  
  for (int m = -K[2]/2; m <= K[2]/2; ++m){
    u = 2*M_PI*m/double(K[2]);
    for (unsigned i = 0; i < ngrid[2]; ++ i){
      epsilonp (u, h[2]*i, n, ep2[m+sw[2]][i], tmp);
    }
  }
}
