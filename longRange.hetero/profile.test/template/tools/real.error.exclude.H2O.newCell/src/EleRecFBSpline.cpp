#include "ElectrostaticInteraction.h"
#include "Polynominal.h"
#include "VectorOperation.h"
#include <numeric>
#include <math.h>
#include <time.h>
#include "Integral1D.h"
#include "ToolBox.h"

#include <complex>
// static double global_x;
// static double global_u;
// static unsigned global_m;


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


void epsilon (double u, double x, double n, double & re, double & im)
{
  int A = 2;
  re = 0;
  im  = 0;
  if (u == 0) return ;
  double fenmu = 0;
  for (int k = -A; k <= A; ++k){
    double tmp  = 1./intPow(u + 2*M_PI*k, int(n));
    fenmu += tmp;
    re += tmp * (cos(2*M_PI*k*x) - 1);
    im += tmp * sin(2*M_PI*k*x);
  }
  re /= fenmu;
  im /= fenmu;
}

std::complex<double > epsilon (double u, double x, double n)
{
  double re, im;
  int A = 2;
  re = 0;
  im  = 0;
  if (u == 0) return std::complex<double > (0,0);
  double fenmu = 0;
  for (int k = -A; k <= A; ++k){
    double tmp  = 1./intPow(u + 2*M_PI*k, int(n));
    fenmu += tmp;
    re += tmp * (cos(2*M_PI*k*x) - 1);
    im += tmp * sin(2*M_PI*k*x);
  }
  return std::complex<double > (re, im);
}

void epsilonp (double u, double x, double n, double & re, double & im)
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

void ElectrostaticInteraction_rec_FBSpline::selfTerm(
    std::vector<double > m, std::vector<double > K, double n, std::vector<double > r,
    std::vector<double > & resultre, std::vector<double > & resultim)
{
  resultre.resize(3, 0);
  resultim.resize(3, 0);
  
//   std::vector<double > sume(2, 0);

  std::vector<double > re (3, 0), im(3, 0);
//   epsilon (2*M_PI*m[0]/K[0], K[0]*r[0], n, re[0], im[0]); 
//   epsilon (2*M_PI*m[1]/K[1], K[1]*r[1], n, re[1], im[1]); 
//   epsilon (2*M_PI*m[2]/K[2], K[2]*r[2], n, re[2], im[2]); 
//   double newre, newim;
//   newre = - 2*M_PI * (im[0] + im[1] + im[2]);
//   newim = 2 * M_PI * (re[0] + re[1] + re[2]);
//   resultre[0] += 1*m[0]*newre * vecAStar[0][0] + 1*m[0]*newre * vecAStar[1][0] + 1*m[0]*newre * vecAStar[2][0];
//   resultre[1] += 1*m[1]*newre * vecAStar[0][1] + 1*m[1]*newre * vecAStar[1][1] + 1*m[1]*newre * vecAStar[2][1];
//   resultre[2] += 1*m[2]*newre * vecAStar[0][2] + 1*m[2]*newre * vecAStar[1][2] + 1*m[2]*newre * vecAStar[2][2];
//   resultim[0] += 1*m[0]*newim * vecAStar[0][0] + 1*m[0]*newim * vecAStar[1][0] + 1*m[0]*newim * vecAStar[2][0];
//   resultim[1] += 1*m[1]*newim * vecAStar[0][1] + 1*m[1]*newim * vecAStar[1][1] + 1*m[1]*newim * vecAStar[2][1];
//   resultim[2] += 1*m[2]*newim * vecAStar[0][2] + 1*m[2]*newim * vecAStar[1][2] + 1*m[2]*newim * vecAStar[2][2];  

//   newre = 2 * M_PI * (im[0] + im[1] + im[2]);
//   newim = 2 * M_PI * (re[0] + re[1] + re[2]);  
//   resultre[0] += 1*m[0]*newre * vecAStar[0][0] + 1*m[0]*newre * vecAStar[1][0] + 1*m[0]*newre * vecAStar[2][0];
//   resultre[1] += 1*m[1]*newre * vecAStar[0][1] + 1*m[1]*newre * vecAStar[1][1] + 1*m[1]*newre * vecAStar[2][1];
//   resultre[2] += 1*m[2]*newre * vecAStar[0][2] + 1*m[2]*newre * vecAStar[1][2] + 1*m[2]*newre * vecAStar[2][2];
//   resultim[0] += 1*m[0]*newim * vecAStar[0][0] + 1*m[0]*newim * vecAStar[1][0] + 1*m[0]*newim * vecAStar[2][0];
//   resultim[1] += 1*m[1]*newim * vecAStar[0][1] + 1*m[1]*newim * vecAStar[1][1] + 1*m[1]*newim * vecAStar[2][1];
//   resultim[2] += 1*m[2]*newim * vecAStar[0][2] + 1*m[2]*newim * vecAStar[1][2] + 1*m[2]*newim * vecAStar[2][2];  

//   VectorOperation::add (resultre, 1*m[0]*newre, vecAStar[0]);
//   VectorOperation::add (resultre, 1*m[1]*newre, vecAStar[1]);
//   VectorOperation::add (resultre, 1*m[2]*newre, vecAStar[2]);
//   VectorOperation::add (resultim, 1*m[0]*newim, vecAStar[0]);
//   VectorOperation::add (resultim, 1*m[1]*newim, vecAStar[1]);
//   VectorOperation::add (resultim, 1*m[2]*newim, vecAStar[2]);
  
  epcal.get_ep (K[0]*r[0], K[1]*r[1], K[2]*r[2], 
		m[0], m[1], m[2],
		re[0], re[1], re[2]);

//   epsilonp (2*M_PI*m[0]/K[0], K[0]*r[0], n, re[0], im[0]);
//   epsilonp (2*M_PI*m[1]/K[1], K[1]*r[1], n, re[1], im[1]);
//   epsilonp (2*M_PI*m[2]/K[2], K[2]*r[2], n, re[2], im[2]);
  resultre[0] = K[0]*re[0] * vecAStar[0][0] + K[1]*re[1] * vecAStar[1][0] + K[2]*re[2] * vecAStar[2][0];
  resultre[1] = K[0]*re[0] * vecAStar[0][1] + K[1]*re[1] * vecAStar[1][1] + K[2]*re[2] * vecAStar[2][1];
  resultre[2] = K[0]*re[0] * vecAStar[0][2] + K[1]*re[1] * vecAStar[1][2] + K[2]*re[2] * vecAStar[2][2];
//   resultim[0] = K[0]*im[0] * vecAStar[0][0] + K[1]*im[1] * vecAStar[1][0] + K[2]*im[2] * vecAStar[2][0];
//   resultim[1] = K[0]*im[0] * vecAStar[0][1] + K[1]*im[1] * vecAStar[1][1] + K[2]*im[2] * vecAStar[2][1];
//   resultim[2] = K[0]*im[0] * vecAStar[0][2] + K[1]*im[1] * vecAStar[1][2] + K[2]*im[2] * vecAStar[2][2];
  

//   VectorOperation::add (resultre, K[0] * re[0], vecAStar[0]);
//   VectorOperation::add (resultim, K[0] * im[0], vecAStar[0]);
//   VectorOperation::add (resultre, K[1] * re[1], vecAStar[1]);
//   VectorOperation::add (resultim, K[1] * im[1], vecAStar[1]);
//   VectorOperation::add (resultre, K[2] * re[2], vecAStar[2]);
//   VectorOperation::add (resultim, K[2] * im[2], vecAStar[2]);
}  

  

double ElectrostaticInteraction_rec_FBSpline::gradE2 (std::vector<double > m, std::vector<double > K, double n, std::vector<double > r)
{
  double mr = VectorOperation::dot(m, r);
  std::vector<double > sume(2, 0);
  double re, im;
  epsilon (2*M_PI*m[0]/K[0], K[0]*r[0], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[1]/K[1], K[1]*r[1], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[2]/K[2], K[2]*r[2], n, re, im); sume[0] += re; sume[1] += im;
  double term0re, term0im;
  term0im = (cos(2*M_PI*mr) * sume[0] - sin(2*M_PI*mr) * sume[1]) * 2 * M_PI;
  term0re = (cos(2*M_PI*mr) * sume[1] + sin(2*M_PI*mr) * sume[0]) * 2 * M_PI * (-1);

  std::vector<double > term1re(3, 0);
  std::vector<double > term1im(3, 0);
  double re1, im1;
  epsilonp (2*M_PI*m[0]/K[0], K[0]*r[0], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[0] * re1, vecAStar[0]);
  VectorOperation::add (term1im, K[0] * im1, vecAStar[0]);
  epsilonp (2*M_PI*m[1]/K[1], K[1]*r[1], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[1] * re1, vecAStar[1]);
  VectorOperation::add (term1im, K[1] * im1, vecAStar[1]);
  epsilonp (2*M_PI*m[2]/K[2], K[2]*r[2], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[2] * re1, vecAStar[2]);
  VectorOperation::add (term1im, K[2] * im1, vecAStar[2]);
  
  std::vector<double > gradEre(3, 0);
  std::vector<double > gradEim (3, 0);
  VectorOperation::add (gradEre, m[0] * term0re, vecAStar[0]);
  VectorOperation::add (gradEre, m[1] * term0re, vecAStar[1]);
  VectorOperation::add (gradEre, m[2] * term0re, vecAStar[2]);
  VectorOperation::add (gradEim, m[0] * term0im, vecAStar[0]);
  VectorOperation::add (gradEim, m[1] * term0im, vecAStar[1]);
  VectorOperation::add (gradEim, m[2] * term0im, vecAStar[2]);
  VectorOperation::add (gradEre, 1, term1re);
  VectorOperation::add (gradEim, 1, term1im);
  
  return VectorOperation::dot(gradEre, gradEre) + VectorOperation::dot(gradEim, gradEim);
} 


void ElectrostaticInteraction_rec_FBSpline::gradEgrade (
    std::vector<double > m, std::vector<double > K, double n, std::vector<double > r,
    double & multi0re, double & multi0im, double & multi1re, double & multi1im)
{
  double mr = VectorOperation::dot(m, r);
  std::vector<double > sume(2, 0);
  double re, im;
  epsilon (2*M_PI*m[0]/K[0], K[0]*r[0], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[1]/K[1], K[1]*r[1], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[2]/K[2], K[2]*r[2], n, re, im); sume[0] += re; sume[1] += im;
  double term0re, term0im;
  term0im = (cos(2*M_PI*mr) * sume[0] - sin(2*M_PI*mr) * sume[1]) * 2 * M_PI;
  term0re = (cos(2*M_PI*mr) * sume[1] + sin(2*M_PI*mr) * sume[0]) * 2 * M_PI * (-1);

  std::vector<double > term1re(3, 0);
  std::vector<double > term1im(3, 0);
  double re1, im1;
  epsilonp (2*M_PI*m[0]/K[0], K[0]*r[0], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[0] * re1, vecAStar[0]);
  VectorOperation::add (term1im, K[0] * im1, vecAStar[0]);
  epsilonp (2*M_PI*m[1]/K[1], K[1]*r[1], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[1] * re1, vecAStar[1]);
  VectorOperation::add (term1im, K[1] * im1, vecAStar[1]);
  epsilonp (2*M_PI*m[2]/K[2], K[2]*r[2], n, re, im);
  re1 = re * cos(2*M_PI*mr) - im * sin(2*M_PI*mr);
  im1 = re * sin(2*M_PI*mr) + im * cos(2*M_PI*mr);
  VectorOperation::add (term1re, K[2] * re1, vecAStar[2]);
  VectorOperation::add (term1im, K[2] * im1, vecAStar[2]);
  
  std::vector<double > gradEre(3, 0);
  std::vector<double > gradEim (3, 0);
  VectorOperation::add (gradEre, m[0] * term0re, vecAStar[0]);
  VectorOperation::add (gradEre, m[1] * term0re, vecAStar[1]);
  VectorOperation::add (gradEre, m[2] * term0re, vecAStar[2]);
  VectorOperation::add (gradEim, m[0] * term0im, vecAStar[0]);
  VectorOperation::add (gradEim, m[1] * term0im, vecAStar[1]);
  VectorOperation::add (gradEim, m[2] * term0im, vecAStar[2]);
  VectorOperation::add (gradEre, 1, term1re);
  VectorOperation::add (gradEim, 1, term1im);
  
  std::vector<double > termGere(3, 0);
  std::vector<double > termGeim(3, 0);
  re = 2 * M_PI * sin(-2*M_PI*mr);
  im = -2* M_PI * cos(-2*M_PI*mr);
  VectorOperation::add (termGere, re * m[0], vecAStar[0]);
  VectorOperation::add (termGere, re * m[1], vecAStar[1]);
  VectorOperation::add (termGere, re * m[2], vecAStar[2]);
  VectorOperation::add (termGeim, im * m[0], vecAStar[0]);
  VectorOperation::add (termGeim, im * m[1], vecAStar[1]);
  VectorOperation::add (termGeim, im * m[2], vecAStar[2]);
  
  multi0re = 0;
  multi0im = 0;
  multi0re += gradEre[0] * termGere[0] - gradEim[0] * termGeim[0];
  multi0im += gradEre[0] * termGeim[0] + gradEim[0] * termGere[0];
  multi0re += gradEre[1] * termGere[1] - gradEim[1] * termGeim[1];
  multi0im += gradEre[1] * termGeim[1] + gradEim[1] * termGere[1];
  multi0re += gradEre[2] * termGere[2] - gradEim[2] * termGeim[2];
  multi0im += gradEre[2] * termGeim[2] + gradEim[2] * termGere[2];

  sume[0] = (sume[1] = 0);
  epsilon (2*M_PI*m[0]/K[0], K[0]*r[0], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[1]/K[1], K[1]*r[1], n, re, im); sume[0] += re; sume[1] += im;
  epsilon (2*M_PI*m[2]/K[2], K[2]*r[2], n, re, im); sume[0] += re; sume[1] += im;
  multi1re = sume[0];
  multi1im = sume[1];
//   std::cout << "re is " << multi0re * sume[0] - multi0im * sume[1] << std::endl;
//   std::cout << "im is " << multi0re * sume[1] + multi0im * sume[0] << std::endl;
//   resultre = multi0re * sume[0] - multi0im * sume[1];
//   resultim = multi0re * sume[1] + multi0im * sume[0];
} 




inline double sixTerms (const double & x, const int & n)
{
  double sum = 0;
  sum += 1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI);
  sum += 1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI);
  sum += 1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI);
//   sum += 1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI);
  return sum;
}
inline double sixTermsp (const double & x, const int & n)
{
  double sum = 0;
  sum += (1./intPow(n, x + 2*M_PI) - 1./intPow(n, x - 2*M_PI));
  sum += (1./intPow(n, x + 4*M_PI) - 1./intPow(n, x - 4*M_PI)) * 2;
  sum += (1./intPow(n, x + 6*M_PI) - 1./intPow(n, x - 6*M_PI)) * 3;
//   sum += (1./intPow(n, x + 8*M_PI) - 1./intPow(n, x - 8*M_PI)) * 4;
  return sum ;
}
inline double sixTermspp (const double & x, const int & n)
{
//   std::cout << n << std::endl;
  double sum = 0;
  sum += (1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI));
  sum += (1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI)) * 4;
  sum += (1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI)) * 9;
//   sum += (1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI)) * 16;
//   sum += (1./intPow(n, x + 10*M_PI) + 1./intPow(n, x - 10*M_PI)) * 25;
  return sum ;
}

inline double sevenTerms (const double & x, const int & n)
{
  double sum = 1./intPow(n, x);
  sum += 1./intPow(n, x + 2*M_PI) + 1./intPow(n, x - 2*M_PI);
  sum += 1./intPow(n, x + 4*M_PI) + 1./intPow(n, x - 4*M_PI);
  sum += 1./intPow(n, x + 6*M_PI) + 1./intPow(n, x - 6*M_PI);
//   sum += 1./intPow(n, x + 8*M_PI) + 1./intPow(n, x - 8*M_PI);
//   sum += 1./intPow(n, x + 10*M_PI) + 1./intPow(n, x - 10*M_PI);
  return sum;
}

inline double denominator (const double & x, const int & n)
{
  return sevenTerms (x, n) ;
// + 
//       ( 1./ (2.*M_PI * (n-1) * intPow(x + 2*M_PI*(3-0.5), n-1)) -
// 	1./ (2.*M_PI * (n-1) * intPow(x - 2*M_PI*(3-0.5), n-1)) );
}

inline double E1 (const double & x, const int & n)
{
  if (x == 0) return 0;
  return - sixTerms (x, n) / denominator (x, n);
}

inline double E2 (const double & x, const int & n)
{
  if (x == 0) return 0;
  double tmp0 = sixTerms (x, n);
  double tmp1 = denominator (x, n);
  return (tmp0 * tmp0 + sixTerms (x, n*2)) / (tmp1 * tmp1);
}

inline double Ep2 (const double & x, const int & n)
{
  if (x == 0) return 0;
  double tmp = denominator (x, n);
  return 2 * M_PI / tmp / tmp * sixTermsp (x, 2*n);
}

inline double Epp2 (const double & x, const int & n)
{
  if (x == 0) return 0;
  double tmp = denominator (x, n);
  return 4 * M_PI * M_PI / tmp / tmp * sixTermspp (x, n*2);
}


value_type ElectrostaticInteraction_rec_FBSpline::errorEstimate (
    const value_type & c2, const int & N)
{
  // unsigned n = Mn->getN();
  double sum = 0;
  std::vector<double > KK (3);
  KK[0] = K[0];
  KK[1] = K[1];
  KK[2] = K[2];  


//////////////////////////////////////////////////
// calculate self term
  double sumre = 0;
//   double sumim = 0;
  ToolBox::init_genrand (40000);
  int noc=  10;
  std::vector<double > fms ((K[0]+1)*(K[0]+1)*(K[2]+1));

  for (int count = 0; count < noc; ++count){
    std::vector<double >  tmpr(3);
    tmpr[0] = K[0]*ToolBox::genrand_real3();
    tmpr[1] = K[1]*ToolBox::genrand_real3();
    tmpr[2] = K[2]*ToolBox::genrand_real3();
    std::vector<double > Fre(3, 0);
    std::vector<double > Fim(3, 0);
    std::vector<double >::iterator i = fms.begin();
    for (int m0 = -KK[0]*0.5; m0 <= KK[0]*0.5; ++m0){
      for (int m1 = -KK[1]*0.5; m1 <= KK[1]*0.5; ++m1){
	for (int m2 = -KK[2]*0.5; m2 <= KK[2]*0.5; ++m2){
	  if (fabs(m0) + fabs(m1) + fabs(m2) == 0) continue;
	  double fm;
	  if (count == 0){
	    std::vector<double > m (3);
	    m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	    m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	    m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	    double mm = VectorOperation::dot (m, m);
	    fm = exp (-M_PI*M_PI*mm/beta/beta) / (2*M_PI*V*mm);
	    *(i++) = fm;
	  }
	  else {
	    fm = *(i++);
	  }
	  
	  std::vector<double > resultre(3);
	  std::vector<double > re(3);
	  epcal.get_ep (tmpr[0], tmpr[1], tmpr[2], m0, m1, m2, re[0], re[1], re[2]);
	  resultre[0] = KK[0]*re[0] * vecAStar[0][0] + KK[1]*re[1] * vecAStar[1][0] + KK[2]*re[2] * vecAStar[2][0];
	  resultre[1] = KK[0]*re[0] * vecAStar[0][1] + KK[1]*re[1] * vecAStar[1][1] + KK[2]*re[2] * vecAStar[2][1];
	  resultre[2] = KK[0]*re[0] * vecAStar[0][2] + KK[1]*re[1] * vecAStar[1][2] + KK[2]*re[2] * vecAStar[2][2];
	    
	  VectorOperation::add (Fre, 2 * fm * c2 / double(N), resultre);
	}
      }
    }
    sumre += VectorOperation::dot (Fre, Fre);
//     sumim += VectorOperation::dot (Fim, Fim);
////     std::cout << count << "\r ";
////     std::cout << std::flush;
  }
  sumre /= double (noc);
//   sumim /= double (noc);
////  std::cout << "\nsumre is " << sumre << std::endl;
//   std::cout << "sumim is " << sumim << std::endl;  
//////////////////////////////////////////////////


  for (int m0 = -KK[0]*0.5; m0 <= KK[0]*0.5; ++m0){
    for (int m1 = -KK[1]*0.5; m1 <= KK[1]*0.5; ++m1){
      for (int m2 = -KK[2]*0.5; m2 <= KK[2]*0.5; ++m2){
	if (fabs(m0) + fabs(m1) + fabs(m2) == 0) continue;
	std::vector<double > m (3);
	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	double mm = VectorOperation::dot (m, m);
	double fm = exp (-M_PI*M_PI*mm/beta/beta) / (2*M_PI*V*mm);
// 	std::vector<double > u(3);
// 	u[0] = 2*M_PI*m0 / KK[0];
// 	u[1] = 2*M_PI*m1 / KK[1];
// 	u[2] = 2*M_PI*m2 / KK[2];
	std::vector<double > E1s(3);
// 	E1s[0] = E1 (u[0], n);
// 	E1s[1] = E1 (u[1], n);
// 	E1s[2] = E1 (u[2], n);
	std::vector<double > E2s(3);
// 	E2s[0] = E2 (u[0], n);
// 	E2s[1] = E2 (u[1], n);
// 	E2s[2] = E2 (u[2], n);
	std::vector<double > Ep2s(3);
// 	Ep2s[0] = Ep2 (u[0], n);
// 	Ep2s[1] = Ep2 (u[1], n);
// 	Ep2s[2] = Ep2 (u[2], n);
	std::vector<double > Epp2s(3);
// 	Epp2s[0] = Epp2 (u[0], n);
// 	Epp2s[1] = Epp2 (u[1], n);
// 	Epp2s[2] = Epp2 (u[2], n);
	
	ecal.get_e1 (m0, m1, m2, E1s[0], E1s[1], E1s[2]);
	ecal.get_e2 (m0, m1, m2, E2s[0], E2s[1], E2s[2]);
	ecal.get_e2p (m0, m1, m2, Ep2s[0], Ep2s[1], Ep2s[2]);
	ecal.get_e2pp (m0, m1, m2, Epp2s[0], Epp2s[1], Epp2s[2]);
	
	double t1 = E2s[0] + E2s[1] + E2s[2] + 
	    2 * E1s[0] * E1s[1] + 2 * E1s[0] * E1s[2] + 2 * E1s[1] * E1s[2] ;
	double t2 = (E1s[0] + E1s[1] + E1s[2]) * (E1s[0] + E1s[1] + E1s[2]);
	double t30 = 4 * M_PI * M_PI * mm * t1;

	

	double t31 = 0;
	std::vector<double > tmp (3, 0);
	VectorOperation::add (tmp, Ep2s[0] * KK[0], vecAStar[0]);
	VectorOperation::add (tmp, Ep2s[1] * KK[1], vecAStar[1]);
	VectorOperation::add (tmp, Ep2s[2] * KK[2], vecAStar[2]);
	t31 += 1*4 * M_PI * VectorOperation::dot(m, tmp);
	t31 += Epp2s[0] * KK[0] * KK[0] * VectorOperation::dot(vecAStar[0], vecAStar[0]);
	t31 += Epp2s[1] * KK[1] * KK[1] * VectorOperation::dot(vecAStar[1], vecAStar[1]);
	t31 += Epp2s[2] * KK[2] * KK[2] * VectorOperation::dot(vecAStar[2], vecAStar[2]);
	double t3 = t30 + t31;


// 	std::vector<double > tmpm(3), tmpr(3);
// 	tmpm[0] = m0;
// 	tmpm[1] = m1;
// 	tmpm[2] = m2;
// 	int noc=  100;
// 	double mysum = 0;
// 	for (int count = 0; count < noc; ++count){
// 	  tmpr[0] = ToolBox::genrand_real1();
// 	  tmpr[1] = ToolBox::genrand_real1();
// 	  tmpr[2] = ToolBox::genrand_real1();
// 	  mysum += gradE2 (tmpm, KK, n, tmpr);
// 	}
// 	t3 = mysum / double(noc);
	
	sum += 4 * fm * fm * (4 * M_PI * M_PI * mm * t1 +
			      2 * 4 * M_PI * M_PI * mm * t2 +
			      t3);


// 	if (m0 == 4 && m1 == -2 && m2 == 3){
// // 	if (m1 == 0 && m2 == 0){

// 	  std::vector<double > tmpm(3), tmpr(3);
// 	  tmpm[0] = m0;
// 	  tmpm[1] = m1;
// 	  tmpm[2] = m2;
// 	  int noc=  40000;
// 	  ToolBox::init_genrand(1020);
// 	  double mysuma = 0, mysumb = 0, mysumc = 0, mysumd = 0;
// 	  double mysum3 = 0;
// 	  for (int count = 0; count < noc; ++count){
// 	    tmpr[0] = ToolBox::genrand_real1();
// 	    tmpr[1] = ToolBox::genrand_real1();
// 	    tmpr[2] = ToolBox::genrand_real1();
// 	    mysum3 += gradE2 (tmpm, KK, n, tmpr);
// 	    double a, b, c, d;
// 	    gradEgrade (tmpm, KK, n, tmpr, a, b, c, d);
// 	    mysuma += a;
// 	    mysumb += b;
// 	    mysumc += c;
// 	    mysumd += d;
// 	  }
// 	  mysuma /= double(noc);
// 	  mysumb /= double(noc);
// 	  mysumc /= double(noc);
// 	  mysumd /= double(noc);
// 	  std::cout << "MC mysum2 is " 
// 		    << 8 * fm * fm * (mysuma * mysumc - mysumb * mysumd) << std::endl;
// 	  std::cout << "MC mysum2 is "
// 		    << 8 * fm * fm * (mysuma * mysumd + mysumb * mysumc) << std::endl;
// 	  std::cout << "MC mysum3 is " << 4 * fm * fm * mysum3 / (double)(noc) << std::endl;

// // 	  std::copy (E2s.begin(), E2s.end(), std::ostream_iterator<double >(std::cout, "\t"));
// // 	  std::cout << std::endl;
// // 	  std::copy (Ep2s.begin(), Ep2s.end(), std::ostream_iterator<double >(std::cout, "\t"));
// // 	  std::cout << std::endl;
// // 	  std::copy (Epp2s.begin(), Epp2s.end(), std::ostream_iterator<double >(std::cout, "\t"));
// // 	  std::cout << std::endl;
// // 	  std::cout << mm << std::endl;
// // 	  std::cout << m0 << "\t" << m0 * m0 * VectorOperation::dot(vecAStar[0], vecAStar[0]) << std::endl;
// // 	  std::cout << KK[0] << "\t" << KK[0] * KK[0] * VectorOperation::dot(vecAStar[0], vecAStar[0]) << std::endl;
// // 	  std::cout << 4 * fm * fm * 4 * M_PI * VectorOperation::dot(m, tmp) << std::endl;
	  
// // 	  std::cout <<  mm<< std::endl;
// // 	  std::cout << t1 << std::endl;
// 	  std::cout << 4 * fm * fm * 4 * M_PI * M_PI * mm * t1 + 4 * fm * fm * t30 << "\n";
// 	  std::cout << 4 * fm * fm * 8 * M_PI * M_PI * mm * t2 << "\n";
// 	  std::cout << 4 * fm * fm * (t31+t30) << "\t"
// 		    << 4 * fm * fm * 4 * M_PI * VectorOperation::dot(m, tmp) << std::endl;
// // 	  std::cout << fm * fm << std::endl;
// // 	  std::cout << mm << std::endl;
// // 	  std::cout << m0 << "\t" << m1 << "\t" << m2 << std::endl;
// // 	  exit (0);
// 	}
      }
    }
  }
  sum *= (c2 * c2 ) / double (N);
  sum += sumre;
  
  return sqrt(sum);
}

  

value_type ElectrostaticInteraction_rec_FBSpline::errorEstimate (
    std::vector<StandardParticle * >::const_iterator beginPool,
    std::vector<StandardParticle * >::const_iterator endPool) 
{
  std::vector<double > qs;
  for (std::vector<StandardParticle * >::const_iterator pp = beginPool;
       pp != endPool; ++pp){
    qs.push_back ((*pp)->charge());
  }
  double Q = 0;
  for (std::vector<double >::iterator i = qs.begin();
       i != qs.end(); ++ i){
    Q += *i * *i;
  }

  return errorEstimate (Q, qs.size());
}


// value_type ElectrostaticInteraction_rec_FBSpline::errorEstimateDir (double rcut)
// {
  
//   std::vector<double > qs;
//   for (std::vector<StandardParticle * >::iterator pp = partPool.begin();
//        pp != partPool.end(); ++pp){
//     qs.push_back ((*pp)->charge());
//   }
//   std::sort (qs.begin(), qs.end());
//   double qmax = qs[qs.size()-1];
//   double Q = 0;
//   for (std::vector<double >::iterator i = qs.begin();
//        i != qs.end(); ++ i){
//     Q += *i * *i;
//   }  

//   return 2 * qmax * sqrt (Q / rcut / V) * exp (- beta*beta*rcut*rcut);
// }



// void ElectrostaticInteraction_rec_FBSpline::calForce_tmp (
//     std::vector<StandardParticle * >::iterator beginPool,
//     std::vector<StandardParticle * >::iterator endPool,
//     std::vector<double > & force)
// {
//   int M = 32;
//   force.clear();
//   force.resize(3, 0);
//   std::cout << "tmp force calculater" << std::endl;
  
//   for (int m0 = -M; m0 <= M; m0 ++){
//     for (int m1 = -M; m1 <= M; m1 ++){
//       for (int m2 = -M; m2 <= M; m2 ++){
// 	if (m0 == 0 && m1 == 0 && m2 == 0 ) continue;
// 	std::vector<value_type > m (3);
// 	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
// 	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
// 	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
// 	double mm = VectorOperation::dot (m, m);
// 	double expp = exp (-M_PI * M_PI * mm / beta / beta) / mm;
// 	double sum = 0;
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
//       }
//     }
//   }
// //   std::cout << "tmp force calculater: force Of part" << std::endl;
// //   std::cout << (*beginPool)->r() [0] << '\t'
// // 	    << (*beginPool)->r() [1] << '\t'
// // 	    << (*beginPool)->r() [2] << '\n';
//   std::cout << "the force is " << std::endl;
//   std::cout << force [0] << '\t'
// 	    << force [1] << '\t'
// 	    << force [2] << '\n';
//   std::cout << std::endl;
  
// }


// value_type ElectrostaticInteraction_rec_FBSpline::calPotential_tmp (
//     std::vector<StandardParticle * >::iterator beginPool,
//     std::vector<StandardParticle * >::iterator endPool)
// {
//   int M = 32; 
//   std::cout << "cal poten exact with beta " << beta << std::endl;
//   double poten = 0;
//   for (int m0 = -M; m0 <= M; m0 ++){
//     for (int m1 = -M; m1 <= M; m1 ++){
//       for (int m2 = -M; m2 <= M; m2 ++){
// 	if (m0 == 0 && m1 == 0 && m2 == 0 ) continue;
// 	std::vector<value_type > m (3);
// 	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
// 	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
// 	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
// 	double mm = VectorOperation::dot (m, m);
// 	double expp = exp (-M_PI * M_PI * mm / beta / beta) / mm;
// 	std::vector<value_type > sm (2);
// 	sm[0] = 0;
// 	sm[1] = 0;
// 	for (std::vector<StandardParticle * >::iterator ppart = beginPool;
// 	     ppart != endPool; ppart ++){
// 	  double tmp = 2 * M_PI * VectorOperation::dot (m, (*ppart)->r());
// 	  sm[0] += (*ppart)->charge() * cos(tmp);
// 	  sm[1] += (*ppart)->charge() * sin(tmp);
// 	}
// 	poten += expp * (sm[0] * sm[0] + sm[1] * sm[1]);
//       }
//     }
//   }
//   poten /= 2 * M_PI * V;
//   return poten;
// }




void ElectrostaticInteraction_rec_FBSpline::applyInteraction (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
  applyInteractionCalPotential (beginPool, endPool);
}


value_type ElectrostaticInteraction_rec_FBSpline::calPotential (
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
    QConvPsi[i] *= sizei;
    value += Q[i] * QConvPsi[i];
  }

  return value;
}

value_type ElectrostaticInteraction_rec_FBSpline::applyInteractionCalPotential (
    std::vector<StandardParticle * >::iterator beginPool,
    std::vector<StandardParticle * >::iterator endPool)
{
  // cal potential part
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
    QConvPsi[i] *= sizei;
    value += Q[i] * QConvPsi[i];
  }
  
  // apply interaction part 

  unsigned n = Mn->getN();
//   watch.start();
//   watch.stop();
//   std::cout << "time rec force conv part: " << (toc - tic) / double(CLOCKS_PER_SEC) << std::endl;

  watch.start();
  double ii0 = 1./ value_type(K[0]);
  double ii1 = 1./ value_type(K[1]);
  double ii2 = 1./ value_type(K[2]);
  bool fast = ((n < K[0]) && (n < K[1]) && (n < K[2]));
  
  int count = 0;
  for (std::vector<StandardParticle * >::iterator ppart = beginPool;
       ppart != endPool; ppart ++, count++){
    std::vector<value_type > u(3);
    u[0] = K[0] * VectorOperation::dot (vecAStar[0], (*ppart)->r());
    u[1] = K[1] * VectorOperation::dot (vecAStar[1], (*ppart)->r());
    u[2] = K[2] * VectorOperation::dot (vecAStar[2], (*ppart)->r());
    value_type tmp0 = 0;
    value_type tmp1 = 0;
    value_type tmp2 = 0;
    value_type dtmp0 = 0;
    value_type dtmp1 = 0;
    value_type dtmp2 = 0;
    std::vector<double > force (3, 0);
    
    if (!fast){
      int count = 0;
      for (unsigned k0 = 0; k0 < K[0]; k0 ++){
	for (unsigned k1 = 0; k1 < K[1]; k1 ++){
	  for (unsigned k2 = 0; k2 < K[2]; k2 ++, count ++){
	    value_type n0 = floor ((u[0] - k0) * ii0);
	    if (u[0] - k0 - n0 * K[0] < n){
	      Mn->value (u[0] - k0 - n0 * K[0], tmp0);
	      Mn->derivative (u[0] - k0 - n0 * K[0], dtmp0);
// 	      std::cout << tmp0 <<'\t' << dtmp0 << std::endl;
	    }
	    else{
	      continue;
	    }
	    value_type n1 = floor ((u[1] - k1) * ii1);
	    if (u[1] - k1 - n1 * K[1] < n){
	      Mn->value (u[1] - k1 - n1 * K[1], tmp1);
	      Mn->derivative (u[1] - k1 - n1 * K[1], dtmp1);
// 	      std::cout << tmp1 <<'\t' << dtmp1 << std::endl;
	    }
	    else{
	      continue;
	    }
	    value_type n2 = floor ((u[2] - k2) * ii2);
	    if (u[2] - k2 - n2 * K[2] < n){
	      Mn->value (u[2] - k2 - n2 * K[2], tmp2);
	      Mn->derivative (u[2] - k2 - n2 * K[2], dtmp2);
// 	      std::cout << tmp2 <<'\t' << dtmp2 << std::endl;
// 	      exit(1);
	    }
	    else{
	      continue;
	    }
	    VectorOperation::add (force, -QConvPsi[count] * dtmp0 * tmp1 * tmp2 * K[0], vecAStar[0]);
	    VectorOperation::add (force, -QConvPsi[count] * tmp0 * dtmp1 * tmp2 * K[1], vecAStar[1]);
	    VectorOperation::add (force, -QConvPsi[count] * tmp0 * tmp1 * dtmp2 * K[2], vecAStar[2]);
	  }
	}
      }
    }
    else{
      int A0 = -int(floor ((u[0]) * ii0)) ;
      int A1 = -int(floor ((u[1]) * ii1)) ;
      int A2 = -int(floor ((u[2]) * ii2)) ;
      value_type posi0 = u[0] + A0 * K[0];
      value_type posi1 = u[1] + A1 * K[1];
      value_type posi2 = u[2] + A2 * K[2];
      unsigned index0, index1;

      for (int k0 = int(ceil(posi0-n)); k0 < int(ceil(posi0)); ++ k0){
	Mn->value (posi0 - k0, tmp0);
	Mn->derivative (posi0 - k0, dtmp0);
	index0 = K[1] * (k0<0 ? k0+K[0] : k0);
	for (int k1 = int(ceil(posi1-n)); k1 < int(ceil(posi1)); ++ k1){
	  Mn->value (posi1 - k1, tmp1);
	  Mn->derivative (posi1 - k1, dtmp1);
	  index1 = K[2] * ((k1<0 ? k1+K[1] : k1) + index0);
	  for (int k2 = int(ceil(posi2-n)); k2 < int(ceil(posi2)); ++ k2){
	    Mn->value (posi2 - k2, tmp2);
	    Mn->derivative (posi2 - k2, dtmp2);
	    unsigned index = (k2<0 ? k2+K[2] : k2) + index1;
	    VectorOperation::add (force, -QConvPsi[index] * dtmp0 * tmp1 * tmp2 * K[0], vecAStar[0]);
	    VectorOperation::add (force, -QConvPsi[index] * tmp0 * dtmp1 * tmp2 * K[1], vecAStar[1]);
	    VectorOperation::add (force, -QConvPsi[index] * tmp0 * tmp1 * dtmp2 * K[2], vecAStar[2]);
	    // if (count == 4){
	    //   printf ("%d %f %f %f\n", index, force[0], force[1], force[2]);
	    // }
	  }
	}
      }
    }

    (*ppart)->f()[0] += 2 * (*ppart)->charge() * force[0];
    (*ppart)->f()[1] += 2 * (*ppart)->charge() * force[1];
    (*ppart)->f()[2] += 2 * (*ppart)->charge() * force[2];
  }



  {
    int count = 0;
    for (std::vector<StandardParticle * >::iterator ppart = beginPool;
	 ppart != endPool; ++(++(++ppart)), ++count){
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
	// std::vector<value_type > tmpposi(3);
	// tmpposi[0] = u[0] -int(floor ((u[0]) * ii0)) * K[0];
	// tmpposi[1] = u[1] -int(floor ((u[1]) * ii1)) * K[1];
	// tmpposi[2] = u[2] -int(floor ((u[2]) * ii2)) * K[2];
	if      (u[0] <  0           ) u[0] += double(K[0]);
	else if (u[0] >= double(K[0])) u[0] -= double(K[0]);
	if      (u[1] <  0           ) u[1] += double(K[1]);
	else if (u[1] >= double(K[1])) u[1] -= double(K[1]);
	if      (u[2] <  0           ) u[2] += double(K[2]);
	else if (u[2] >= double(K[2])) u[2] -= double(K[2]);
	posi[molc] = u;
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
      for (int mol1 = 0; mol1 < 3; ++mol1){
	if (mol0 == mol1) continue;
	std::vector<std::vector<double > > mn_mol0(3);
	std::vector<std::vector<double > > dmn_mol0(3);
	std::vector<std::vector<double > > mn_mol1(3);
	for (unsigned dd = 0; dd < 3; ++dd){
	  mn_mol0[dd].resize(n, 0.);
	  dmn_mol0[dd].resize(n, 0.);
	  mn_mol1[dd].resize(n, 0.);
	}	
	for (unsigned dd = 0; dd < 3; ++dd){
	  for (int kk = start[mol0][dd]; kk < start[mol0][dd] + int(n); ++kk){
	    Mn->value (posi[mol0][dd] - kk, mn_mol0[dd][kk - start[mol0][dd]]);
	    Mn->derivative (posi[mol0][dd] - kk, dmn_mol0[dd][kk - start[mol0][dd]]);
	  }
	}
	for (unsigned dd = 0; dd < 3; ++dd){
	  for (int kk = start[mol1][dd]; kk < start[mol1][dd] + int(n); ++kk){
	    Mn->value (posi[mol1][dd] - kk, mn_mol1[dd][kk - start[mol1][dd]]);
	  }
	}	
	  
	std::vector<int > kk(3);
	std::vector<int > ll(3);
	std::vector<int > diffkl (3);
	for (kk[0] = start[mol0][0]; kk[0] < start[mol0][0] + int(n); ++kk[0]){
	for (ll[0] = start[mol1][0]; ll[0] < start[mol1][0] + int(n); ++ll[0]){
	  diffkl[0] = kk[0] - ll[0];
	  if      (diffkl[0] <  0        ) {diffkl[0] += 3*int(K[0]); diffkl[0] = diffkl[0] % int(K[0]);}
	  else if (diffkl[0] >= int(K[0])) diffkl[0] -= int(K[0]);
	for (kk[1] = start[mol0][1]; kk[1] < start[mol0][1] + int(n); ++kk[1]){
	for (ll[1] = start[mol1][1]; ll[1] < start[mol1][1] + int(n); ++ll[1]){
	  diffkl[1] = kk[1] - ll[1];
	  if      (diffkl[1] <  0        ) {diffkl[1] += 3*int(K[1]); diffkl[1] = diffkl[1] % int(K[1]);}
	  else if (diffkl[1] >= int(K[1])) diffkl[1] -= int(K[1]);
	for (kk[2] = start[mol0][2]; kk[2] < start[mol0][2] + int(n); ++kk[2]){
	for (ll[2] = start[mol1][2]; ll[2] < start[mol1][2] + int(n); ++ll[2]){
	  diffkl[2] = kk[2] - ll[2];
	  if      (diffkl[2] <  0        ) {diffkl[2] += 3*int(K[2]); diffkl[2] = diffkl[2] % int(K[2]);}
	  else if (diffkl[2] >= int(K[2])) diffkl[2] -= int(K[2]);
	  // if (count == 2){
	  // printf ("%d   kk:%d %d %d   ll:%d %d %d    kl:%d %d %d    K:%d %d %d\n",
	  // 	  diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0]),
	  // 	  kk[0], kk[1], kk[2],
	  // 	  ll[0], ll[1], ll[2],
	  // 	  diffkl[0], diffkl[1], diffkl[2],
	  // 	  K[0], K[1], K[2]
	  //     );
	  // printf ("%f\n",
	  // 	  psir [diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0]);
	  // printf ("%d   %f\n",
	  // 	  diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0]),
	  // 	  psir [diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0]);
	  // }
	  double psi = psir [diffkl[2] + K[2] * (diffkl[1] + K[1] * diffkl[0])][0];
	  double Pj =
	      mn_mol1[0][ll[0] - start[mol1][0]] *
	      mn_mol1[1][ll[1] - start[mol1][1]] *
	      mn_mol1[2][ll[2] - start[mol1][2]] ;
	  double scalor0 =
	      - 2. * psi *
	      atomcharge[mol0] * atomcharge[mol1] *
	      dmn_mol0[0][kk[0] - start[mol0][0]] *
	      mn_mol0[1][kk[1] - start[mol0][1]] *
	      mn_mol0[2][kk[2] - start[mol0][2]] *
	      K[0] *
	      Pj;
	  double scalor1 =
	      - 2. * psi *
	      atomcharge[mol0] * atomcharge[mol1] *
	      mn_mol0[0][kk[0] - start[mol0][0]] *
	      dmn_mol0[1][kk[1] - start[mol0][1]] *
	      mn_mol0[2][kk[2] - start[mol0][2]] *
	      K[1] *
	      Pj;
	  double scalor2 =
	      - 2. * psi *
	      atomcharge[mol0] * atomcharge[mol1] *
	      mn_mol0[0][kk[0] - start[mol0][0]] *
	      mn_mol0[1][kk[1] - start[mol0][1]] *
	      dmn_mol0[2][kk[2] - start[mol0][2]] *
	      K[2] *
	      Pj;
	  VectorOperation::add (forceCorr[mol0][mol1], scalor0, vecAStar[0]);
	  VectorOperation::add (forceCorr[mol0][mol1], scalor1, vecAStar[1]);
	  VectorOperation::add (forceCorr[mol0][mol1], scalor2, vecAStar[2]);
	}
	}
	}
	}
	}
	}
      }
      }
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
      (*patom)->f()[0] -= forceCorr[1][0][0];
      (*patom)->f()[1] -= forceCorr[1][0][1];
      (*patom)->f()[2] -= forceCorr[1][0][2];
      ++patom;
      (*patom)->f()[0] -= forceCorr[2][0][0];
      (*patom)->f()[1] -= forceCorr[2][0][1];
      (*patom)->f()[2] -= forceCorr[2][0][2];
      (*patom)->f()[0] -= forceCorr[2][1][0];
      (*patom)->f()[1] -= forceCorr[2][1][1];
      (*patom)->f()[2] -= forceCorr[2][1][2];
    } // end for
  } // end out brace  

  watch.stop();
  loop_time += watch.user();
  
  return value;
}



void ElectrostaticInteraction_rec_FBSpline::calQ (
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
      int A0 = -int(floor ((u[0]) * ii0)) ;
      int A1 = -int(floor ((u[1]) * ii1)) ;
      int A2 = -int(floor ((u[2]) * ii2)) ;
      value_type posi0 = u[0] + A0 * K[0];
      value_type posi1 = u[1] + A1 * K[1];
      value_type posi2 = u[2] + A2 * K[2];
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


void ElectrostaticInteraction_rec_FBSpline::init (
    const std::vector<std::vector<value_type > > &vecA_, 
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
      std::cerr << "no such order: " << interpOrder_
		<< " implemented, use order 8" << std::endl;
      Mn = new BSpline8();
  }
  K = K_;
  beta = beta_;
  vecA = vecA_;
  calV();
  calAStar();
  calB ();
  reset_time ();
  ecal.reinit (K, interpOrder_);
  epcal.reinit (K, interpOrder_, std::vector<unsigned >(3, 253));
  unBuild = true;
}

void ElectrostaticInteraction_rec_FBSpline::build ()
{
  int size = K[0] * K[1] * K[2];
  int sizeHalf = K[0] * K[1] * (K[2] / 2 + 1);
  Q	= (value_type *) fftw_malloc (sizeof(value_type) * size);
  psiF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  psir	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  QF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QFProdPsiF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * sizeHalf);
  QConvPsi	= (value_type *) fftw_malloc (sizeof(value_type) * size);

  forwardQ	= fftw_plan_dft_r2c_3d (K[0], K[1], K[2], Q  , QF  , FFTW_MEASURE);
  backwardQFProdPsiF = fftw_plan_dft_c2r_3d (K[0], K[1], K[2], QFProdPsiF, QConvPsi, FFTW_MEASURE);
  backward_psi  = fftw_plan_dft_3d (K[0], K[1], K[2], psir, psir, +1, FFTW_MEASURE);

  calPsiFPhiF();
  unBuild = false;
}


void ElectrostaticInteraction_rec_FBSpline::clear()
{
  delete Mn;
  if (!unBuild){
  fftw_free (Q);
  fftw_free (QF);
  fftw_free (psiF);
  fftw_free (psir);
  fftw_free (QFProdPsiF);
  fftw_free (QConvPsi);
  
  fftw_destroy_plan (forwardQ);
  fftw_destroy_plan (backwardQFProdPsiF);
  fftw_destroy_plan (backward_psi);
  }
}


void ElectrostaticInteraction_rec_FBSpline::calB ()
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




void ElectrostaticInteraction_rec_FBSpline::calPsiFPhiF ()
{
  fftw_complex * psiFtmp;
  
  psiFtmp= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * K[0]*K[1]*K[2]);

  value_type oneOver2PiV = 1. / (2 * M_PI * V);
  value_type scale = M_PI * M_PI / beta / beta;
  unsigned size = K[0]*K[1]*K[2];
  
  psiFtmp[0][0] = 0;
  psiFtmp[0][1] = 0;
  psir[0][0] = 0;
  psir[0][1] = 0;
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
	psiFtmp[k + K[2] * (j + K[1] * i)][0] = oneOver2PiV * expm * size;
	psiFtmp[k + K[2] * (j + K[1] * i)][1] = 0;
	psir[k + K[2] * (j + K[1] * i)][0] = oneOver2PiV * expm;
	psir[k + K[2] * (j + K[1] * i)][1] = 0;
      }
    }
  }

  fftw_execute (backward_psi);
  // for (unsigned i = 0; i < K[0]; i ++){
  //   for (unsigned j = 0; j < K[1]; j ++){
  //     for (unsigned k = 0; k < K[2]; k ++){
  // 	printf ("%f %f\n", psir[k + K[2] * (j + K[1] * i)][0], psir[k + K[2] * (j + K[1] * i)][1]);
  //     }
  //   }
  // }
  
  for (unsigned i = 0; i < K[0]; ++ i){
    for (unsigned j = 0; j < K[1]; ++ j){
      for (unsigned k = 0; k < K[2]/2+1; ++ k){
	psiF[k+(K[2]/2+1)*(j+K[1]*i)][0] = psiFtmp[k+K[2]*(j+K[1]*i)][0];
	psiF[k+(K[2]/2+1)*(j+K[1]*i)][1] = psiFtmp[k+K[2]*(j+K[1]*i)][1];
      }
    }
  }

  fftw_free (psiFtmp);
}




void ElectrostaticInteraction_rec_FBSpline::calV()
{
  V = vecA[0][0] * (vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) - 
    vecA[0][1] * (vecA[1][0]*vecA[2][2] - vecA[2][0]*vecA[1][2]) +
    vecA[0][2] * (vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]);
}
  
void ElectrostaticInteraction_rec_FBSpline::calAStar ()
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





