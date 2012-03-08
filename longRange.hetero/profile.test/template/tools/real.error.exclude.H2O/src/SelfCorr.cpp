#include "SelfCorr.h"
#include <stdlib.h>

const int global_half_N_k_sum = 3;

inline double
calsum (const double & u,
	const int & order)
{
  double sum = 1./gsl_pow_int(u,order);
  double up(u), un(u);
  double incr (2. * M_PI);
  sum += 1./gsl_pow_int((up+=incr),order);
  sum += 1./gsl_pow_int((un-=incr),order);
  sum += 1./gsl_pow_int((up+=incr),order);
  sum += 1./gsl_pow_int((un-=incr),order);
  sum += 1./gsl_pow_int((up+=incr),order);
  sum += 1./gsl_pow_int((un-=incr),order);
  return sum;
}

SelfCorr::
SelfCorr ()
    : malloced (false)
{
}

SelfCorr::
~SelfCorr()
{
  freeAll();
}

void SelfCorr::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void SelfCorr::
calAStar ()
{
  double volumei = double(1.) / volume;
  vecAStar.xx = ( vecA.yy*vecA.zz - vecA.zy*vecA.yz) * volumei;
  vecAStar.yy = ( vecA.xx*vecA.zz - vecA.zx*vecA.xz) * volumei;
  vecAStar.zz = ( vecA.xx*vecA.yy - vecA.yx*vecA.xy) * volumei;
  vecAStar.yx = (-vecA.yx*vecA.zz + vecA.zx*vecA.yz) * volumei;
  vecAStar.zx = ( vecA.yx*vecA.zy - vecA.zx*vecA.yy) * volumei;
  vecAStar.xy = (-vecA.xy*vecA.zz + vecA.zy*vecA.xz) * volumei;
  vecAStar.zy = (-vecA.xx*vecA.zy + vecA.zx*vecA.xy) * volumei;
  vecAStar.xz = ( vecA.xy*vecA.yz - vecA.yy*vecA.xz) * volumei;
  vecAStar.yz = (-vecA.xx*vecA.yz + vecA.yx*vecA.xz) * volumei;
}

void SelfCorr::
freeAll ()
{
  if (malloced){
    free (selfx);
    free (selfy);
    free (selfz);
  }
}

void SelfCorr::
reinit (const double & beta_,
	const int & order_,
	const int & nx,
	const int & ny,
	const int & nz,
	const double & bx,
	const double & by,
	const double & bz)
{
  freeAll();
  
  kcut = global_half_N_k_sum;
  numk = 2 * kcut + 1;
  beta = beta_;
  order = order_;
  K.x = nx;
  K.y = ny;
  K.z = nz;
  Ki.x = 1./double(K.x);
  Ki.y = 1./double(K.y);
  Ki.z = 1./double(K.z);
  vecA.xy = vecA.yx = vecA.yz = vecA.zy = vecA.xz = vecA.zx = 0.;
  vecA.xx = bx;
  vecA.yy = by;
  vecA.zz = bz;
  calV();
  calAStar();

  nele = K.x * K.y * K.z;
  
  selfx = (fftw_complex *) malloc (sizeof(fftw_complex) * (numk));
  selfy = (fftw_complex *) malloc (sizeof(fftw_complex) * (numk));
  selfz = (fftw_complex *) malloc (sizeof(fftw_complex) * (numk));
  malloced = true;

  printf ("# reinit: start build kernel ...");
  fflush (stdout);
  calKernel();
  printf (" done\n");  
  fflush (stdout);
}


void SelfCorr::
calKernel()
{
  double scalor = 1./(2. * M_PI * volume);
  // int nele = K.x * K.y * K.z;
  IntVectorType idx;
  IntVectorType posi;
  bool cleanX = ( ((K.x >> 1) << 1) == K.x);
  bool cleanY = ( ((K.y >> 1) << 1) == K.y);
  bool cleanZ = ( ((K.z >> 1) << 1) == K.z);
  VectorType Ki;
  Ki.x = 1./double(K.x);
  Ki.y = 1./double(K.y);
  Ki.z = 1./double(K.z);
  
  for (int kk = 0; kk < numk; ++kk){
    selfx[kk][0] = selfx[kk][1] = 0.;
    selfy[kk][0] = selfy[kk][1] = 0.;
    selfz[kk][0] = selfz[kk][1] = 0.;
  }
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (idx.x = 0; idx.x < K.x; ++idx.x){
      if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
      else posi.x = idx.x;
      for (idx.y = 0; idx.y < K.y; ++idx.y){
	if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
	else posi.y = idx.y;
	for (idx.z = 0; idx.z < K.z; ++idx.z){
	  if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	  else posi.z = idx.z;
	  VectorType mm;
	  mm.x = posi.x * vecAStar.xx + posi.y * vecAStar.yx + posi.z * vecAStar.zx;
	  mm.y = posi.x * vecAStar.xy + posi.y * vecAStar.yy + posi.z * vecAStar.zy;
	  mm.z = posi.x * vecAStar.xz + posi.y * vecAStar.yz + posi.z * vecAStar.zz;
	  double m2 = (mm.x*mm.x + mm.y*mm.y + mm.z*mm.z);
	  double fm;
	  if (m2 == 0){
	    fm = 0.;
	  }
	  else {
	    fm = kernel_rm1_rec_f (m2, beta) * scalor;
	  }
	  VectorType uu;
	  uu.x = 2. * M_PI * double(posi.x) * Ki.x;
	  uu.y = 2. * M_PI * double(posi.y) * Ki.y;
	  uu.z = 2. * M_PI * double(posi.z) * Ki.z;
	  VectorType fenshu;
	  fenshu.x = 1./(gsl_pow_int(uu.x + 2.*M_PI*kk, order));
	  fenshu.x /= calsum (uu.x, order);
	  fenshu.y = 1./(gsl_pow_int(uu.y + 2.*M_PI*kk, order));
	  fenshu.y /= calsum (uu.y, order);
	  fenshu.z = 1./(gsl_pow_int(uu.z + 2.*M_PI*kk, order));
	  fenshu.z /= calsum (uu.z, order);
	  
	  // unsigned index = index3to1 (idx, K);

	  double Fmx, Fmy, Fmz;
	  Fmx = -2. * fm * fenshu.x;
	  Fmy = -2. * fm * fenshu.y;
	  Fmz = -2. * fm * fenshu.z;
	  if (idx.x == (K.x >> 1) && cleanX) Fmx = 0.;
	  if (idx.y == (K.y >> 1) && cleanY) Fmy = 0.;
	  if (idx.z == (K.z >> 1) && cleanZ) Fmz = 0.;
	  
	  selfx[myk][0] += Fmx * (-4.) * M_PI * kk * double(K.x) * vecAStar.xx;
	  selfy[myk][0] += Fmy * (-4.) * M_PI * kk * double(K.y) * vecAStar.yy;
	  selfz[myk][0] += Fmz * (-4.) * M_PI * kk * double(K.z) * vecAStar.zz;
	  // selfy[myk][0] += Fmy;
	  // selfz[myk][0] += Fmz;
	}
      }
    }
  }
}


void SelfCorr::
correction (const std::vector<double > & posi,
	    const double & charge,
	    std::vector<double > & force)
{

  for (int kk = 1; kk <= kcut; ++kk){
    int myk = kk + kcut;    
    force[0] -= selfx[myk][0] * charge * charge * sin(2. * M_PI * kk * K.x * vecAStar.xx * posi[0]);
    force[1] -= selfy[myk][0] * charge * charge * sin(2. * M_PI * kk * K.y * vecAStar.yy * posi[1]);
    force[2] -= selfz[myk][0] * charge * charge * sin(2. * M_PI * kk * K.z * vecAStar.zz * posi[2]);
  }
}


void SelfCorr::
correction (const std::vector<std::vector<double > > & posi,
	    const std::vector<double > & charge,
	    std::vector<std::vector<double > > & force)
{
  for (unsigned i = 0; i < posi.size(); ++i){
    for (int kk = 1; kk <= kcut; ++kk){
      int myk = kk + kcut;    
      force[i][0] -= selfx[myk][0] * charge[i] * charge[i] * sin(2. * M_PI * kk * K.x * vecAStar.xx * posi[i][0]);
      force[i][1] -= selfy[myk][0] * charge[i] * charge[i] * sin(2. * M_PI * kk * K.y * vecAStar.yy * posi[i][1]);
      force[i][2] -= selfz[myk][0] * charge[i] * charge[i] * sin(2. * M_PI * kk * K.z * vecAStar.zz * posi[i][2]);
    }
  }
}


