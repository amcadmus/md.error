#include <iostream>

#include "ErrorEstimate_SPME_Ik.h"

const int global_half_N_k_sum = 4;

ErrorEstimate_SPME_Ik::
ErrorEstimate_SPME_Ik ()
    : malloced (false)
{
}

ErrorEstimate_SPME_Ik::
~ErrorEstimate_SPME_Ik()
{
  freeAll();
}

void ErrorEstimate_SPME_Ik::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void ErrorEstimate_SPME_Ik::
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

void ErrorEstimate_SPME_Ik::
freeAll ()
{
  if (malloced){
    free (rho1);
    free (rho2);
    free (k1mx);
    free (k1my);
    free (k1mz);
    free (k1rx);
    free (k1ry);
    free (k1rz);
    free (k2ma);
    free (k2mb);
    free (error1x);
    free (error1y);
    free (error1z);
    free (error1);
    free (error2);
    fftw_destroy_plan (p_forward_rho1);
    fftw_destroy_plan (p_forward_rho2);
    fftw_destroy_plan (p_backward_k1mx);
    fftw_destroy_plan (p_backward_k1my);
    fftw_destroy_plan (p_backward_k1mz);
    fftw_destroy_plan (p_forward_k2a);
    fftw_destroy_plan (p_forward_k2b);
    fftw_destroy_plan (p_backward_error1x);
    fftw_destroy_plan (p_backward_error1y);
    fftw_destroy_plan (p_backward_error1z);
    fftw_destroy_plan (p_backward_error2);
  }
}


void ErrorEstimate_SPME_Ik::
reinit (const double & beta_,
	const int & order_,
	// const IntVectorType & Kmax_,
	const DensityProfile_PiecewiseConst & dp)
{
  freeAll();
  
  beta = beta_;
  order = order_;
  K.x = dp.getNx();
  K.y = dp.getNy();
  K.z = dp.getNz();
  vecA.xy = vecA.yx = vecA.yz = vecA.zy = vecA.xz = vecA.zx = 0.;
  vecA.xx = dp.getBox()[0];
  vecA.yy = dp.getBox()[1];
  vecA.zz = dp.getBox()[2];
  calV();
  calAStar();

  nele = K.x * K.y * K.z;
  rho1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rho2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1mx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1my = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1mz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1rx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1ry = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1rz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k2ma = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k2mb = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rho2 = fftw_plan_dft_3d (K.x, K.y, K.z, rho2, rho2, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_k1mx = fftw_plan_dft_3d (K.x, K.y, K.z, k1mx, k1rx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_k1my = fftw_plan_dft_3d (K.x, K.y, K.z, k1my, k1ry, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_k1mz = fftw_plan_dft_3d (K.x, K.y, K.z, k1mz, k1rz, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_k2a = fftw_plan_dft_3d (K.x, K.y, K.z, k2ma, k2ma, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_k2b = fftw_plan_dft_3d (K.x, K.y, K.z, k2mb, k2mb, FFTW_FORWARD, FFTW_MEASURE);
  p_backward_error1x = fftw_plan_dft_3d (K.x, K.y, K.z, error1x, error1x, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1y = fftw_plan_dft_3d (K.x, K.y, K.z, error1y, error1y, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1z = fftw_plan_dft_3d (K.x, K.y, K.z, error1z, error1z, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error2 = fftw_plan_dft_3d (K.x, K.y, K.z, error2, error2, FFTW_BACKWARD, FFTW_MEASURE);
  
  malloced = true;

  printf ("# reinit: start build kernel ...");
  fflush (stdout);
  calKernel();
  printf (" done\n");  
  fflush (stdout);
}

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

static void 
array_multiply (fftw_complex * a,
		const int nele,
		fftw_complex * b,
		fftw_complex * c) 
{
  for (int ii = 0; ii < nele; ++ii){
    // double tmpr, tmpi;
    a[ii][0] =
	c[ii][0] * b[ii][0] - c[ii][1] * b[ii][1];
    a[ii][1] =
	c[ii][0] * b[ii][1] + c[ii][1] * b[ii][0];
  }
}

void ErrorEstimate_SPME_Ik::
calKernel()
{
  double scalor = 1./(2. * M_PI * volume);
  int nele = K.x * K.y * K.z;
  IntVectorType idx;
  IntVectorType posi;
  bool cleanX = ( ((K.x >> 1) << 1) == K.x);
  bool cleanY = ( ((K.y >> 1) << 1) == K.y);
  bool cleanZ = ( ((K.z >> 1) << 1) == K.z);
  VectorType Ki;
  Ki.x = 1./double(K.x);
  Ki.y = 1./double(K.y);
  Ki.z = 1./double(K.z);

  for (int i = 0; i < nele; ++i){
    k2mb[i][0] = k2mb[i][1] = 0.;
  }

  VectorType * gmi = (VectorType *) malloc (sizeof(VectorType) * nele);
  
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
	unsigned index = index3to1 (idx, K);
	if (m2 == 0) {
	  gmi[index].x = gmi[index].y = gmi[index].z = 0.;
	}
	else {
	  double fm = kernel_rm1_rec_f (m2, beta) * scalor;
	  gmi[index].x = - 4. * M_PI * fm * mm.x;
	  gmi[index].y = - 4. * M_PI * fm * mm.y;
	  gmi[index].z = - 4. * M_PI * fm * mm.z;
	}
      }
    }
  }
      
  for (int myk = -global_half_N_k_sum; myk <= global_half_N_k_sum; ++myk){
    if (myk == 0) continue;
    {
      for (idx.x = 0; idx.x < K.x; ++idx.x){
	if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
	else posi.x = idx.x;
	for (idx.y = 0; idx.y < K.y; ++idx.y){
	  if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
	  else posi.y = idx.y;
	  for (idx.z = 0; idx.z < K.z; ++idx.z){
	    if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	    else posi.z = idx.z;
	    unsigned index = index3to1 (idx, K);
	    k1mx[index][0] = k1mx[index][1] = 0.;
	    k1my[index][0] = k1my[index][1] = 0.;
	    k1mz[index][0] = k1mz[index][1] = 0.;
	    if (idx.x == 0 && idx.y == 0 && idx.z == 0) continue;
	    double uu = 2. * M_PI * double(posi.x) * Ki.x;
	    double fenmu = 1./calsum (uu, order);
	    double sum = 1./gsl_pow_int(uu + 2.*M_PI * myk, order) * fenmu;

	    k1mx[index][1] = gmi[index].x * sum;
	    k1my[index][1] = gmi[index].y * sum;
	    k1mz[index][1] = gmi[index].z * sum;
	
	    if (idx.x == (K.x >> 1) && cleanX) k1mx[index][1] = 0.;
	    if (idx.y == (K.y >> 1) && cleanY) k1my[index][1] = 0.;
	    if (idx.z == (K.z >> 1) && cleanZ) k1mz[index][1] = 0.;
	  }
	}
      }
      fftw_execute (p_backward_k1mx);
      fftw_execute (p_backward_k1my);
      fftw_execute (p_backward_k1mz);
  
      for (int i = 0; i < nele; ++i){
	k2mb[i][0] += (
	    (k1rx[i][0] * k1rx[i][0] + k1rx[i][1] * k1rx[i][1]) +
	    (k1ry[i][0] * k1ry[i][0] + k1ry[i][1] * k1ry[i][1]) +
	    (k1rz[i][0] * k1rz[i][0] + k1rz[i][1] * k1rz[i][1]) ) * volume / double(nele);
	k2mb[i][1] = 0.;
      }
    }

    {
      for (idx.x = 0; idx.x < K.x; ++idx.x){
	if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
	else posi.x = idx.x;
	for (idx.y = 0; idx.y < K.y; ++idx.y){
	  if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
	  else posi.y = idx.y;
	  for (idx.z = 0; idx.z < K.z; ++idx.z){
	    if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	    else posi.z = idx.z;
	    unsigned index = index3to1 (idx, K);
	    k1mx[index][0] = k1mx[index][1] = 0.;
	    k1my[index][0] = k1my[index][1] = 0.;
	    k1mz[index][0] = k1mz[index][1] = 0.;
	    if (idx.x == 0 && idx.y == 0 && idx.z == 0) continue;
	    double uu = 2. * M_PI * double(posi.y) * Ki.y;
	    double fenmu = 1./calsum (uu, order);
	    double sum = 1./gsl_pow_int(uu + 2.*M_PI * myk, order) * fenmu;

	    k1mx[index][1] = gmi[index].x * sum;
	    k1my[index][1] = gmi[index].y * sum;
	    k1mz[index][1] = gmi[index].z * sum;
	
	    if (idx.x == (K.x >> 1) && cleanX) k1mx[index][1] = 0.;
	    if (idx.y == (K.y >> 1) && cleanY) k1my[index][1] = 0.;
	    if (idx.z == (K.z >> 1) && cleanZ) k1mz[index][1] = 0.;
	  }
	}
      }
      fftw_execute (p_backward_k1mx);
      fftw_execute (p_backward_k1my);
      fftw_execute (p_backward_k1mz);
  
      for (int i = 0; i < nele; ++i){
	k2mb[i][0] += (
	    (k1rx[i][0] * k1rx[i][0] + k1rx[i][1] * k1rx[i][1]) +
	    (k1ry[i][0] * k1ry[i][0] + k1ry[i][1] * k1ry[i][1]) +
	    (k1rz[i][0] * k1rz[i][0] + k1rz[i][1] * k1rz[i][1]) ) * volume / double(nele);
	k2mb[i][1] = 0.;
      }
    }
    
    {
      for (idx.x = 0; idx.x < K.x; ++idx.x){
	if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
	else posi.x = idx.x;
	for (idx.y = 0; idx.y < K.y; ++idx.y){
	  if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
	  else posi.y = idx.y;
	  for (idx.z = 0; idx.z < K.z; ++idx.z){
	    if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	    else posi.z = idx.z;
	    unsigned index = index3to1 (idx, K);
	    k1mx[index][0] = k1mx[index][1] = 0.;
	    k1my[index][0] = k1my[index][1] = 0.;
	    k1mz[index][0] = k1mz[index][1] = 0.;
	    if (idx.x == 0 && idx.y == 0 && idx.z == 0) continue;
	    double uu = 2. * M_PI * double(posi.z) * Ki.z;
	    double fenmu = 1./calsum (uu, order);
	    double sum = 1./gsl_pow_int(uu + 2.*M_PI * myk, order) * fenmu;

	    k1mx[index][1] = gmi[index].x * sum;
	    k1my[index][1] = gmi[index].y * sum;
	    k1mz[index][1] = gmi[index].z * sum;
	
	    if (idx.x == (K.x >> 1) && cleanX) k1mx[index][1] = 0.;
	    if (idx.y == (K.y >> 1) && cleanY) k1my[index][1] = 0.;
	    if (idx.z == (K.z >> 1) && cleanZ) k1mz[index][1] = 0.;
	  }
	}
      }
      fftw_execute (p_backward_k1mx);
      fftw_execute (p_backward_k1my);
      fftw_execute (p_backward_k1mz);
  
      for (int i = 0; i < nele; ++i){
	k2mb[i][0] += (
	    (k1rx[i][0] * k1rx[i][0] + k1rx[i][1] * k1rx[i][1]) +
	    (k1ry[i][0] * k1ry[i][0] + k1ry[i][1] * k1ry[i][1]) +
	    (k1rz[i][0] * k1rz[i][0] + k1rz[i][1] * k1rz[i][1]) ) * volume / double(nele);
	k2mb[i][1] = 0.;
      }
    }
  }
  
  for (idx.x = 0; idx.x < K.x; ++idx.x){
    if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
    else posi.x = idx.x;
    for (idx.y = 0; idx.y < K.y; ++idx.y){
      if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
      else posi.y = idx.y;
      for (idx.z = 0; idx.z < K.z; ++idx.z){
	if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	else posi.z = idx.z;
	unsigned index = index3to1 (idx, K);
	k1mx[index][0] = k1mx[index][1] = 0.;
	k1my[index][0] = k1my[index][1] = 0.;
	k1mz[index][0] = k1mz[index][1] = 0.;
	if (idx.x == 0 && idx.y == 0 && idx.z == 0) continue;
	VectorType uu;
	uu.x = 2. * M_PI * double(posi.x) * Ki.x;
	uu.y = 2. * M_PI * double(posi.y) * Ki.y;
	uu.z = 2. * M_PI * double(posi.z) * Ki.z;
	VectorType fenmu;
	fenmu.x = 1./calsum (uu.x, order);
	fenmu.y = 1./calsum (uu.y, order);
	fenmu.z = 1./calsum (uu.z, order);
	double sum = 0.;
	for (int myk = -global_half_N_k_sum; myk <= global_half_N_k_sum; ++myk){
	  if (myk == 0) continue;
	  sum += 1./gsl_pow_int(uu.x + 2.*M_PI * myk, order) * fenmu.x;
	  sum += 1./gsl_pow_int(uu.y + 2.*M_PI * myk, order) * fenmu.y;
	  sum += 1./gsl_pow_int(uu.z + 2.*M_PI * myk, order) * fenmu.z;
	}

	k1mx[index][1] = gmi[index].x * sum * 2.;
	k1my[index][1] = gmi[index].y * sum * 2.;
	k1mz[index][1] = gmi[index].z * sum * 2.;
	
	if (idx.x == (K.x >> 1) && cleanX) k1mx[index][1] = 0.;
	if (idx.y == (K.y >> 1) && cleanY) k1my[index][1] = 0.;
	if (idx.z == (K.z >> 1) && cleanZ) k1mz[index][1] = 0.;
      }
    }
  }

  fftw_execute (p_backward_k1mx);
  fftw_execute (p_backward_k1my);
  fftw_execute (p_backward_k1mz);
  
  for (int i = 0; i < nele; ++i){
    k1mx[i][1] *= volume;
    k1my[i][1] *= volume;
    k1mz[i][1] *= volume;
    k2ma[i][0] = (
	(k1rx[i][0] * k1rx[i][0] + k1rx[i][1] * k1rx[i][1]) +
	(k1ry[i][0] * k1ry[i][0] + k1ry[i][1] * k1ry[i][1]) +
	(k1rz[i][0] * k1rz[i][0] + k1rz[i][1] * k1rz[i][1]) ) * volume / double(nele);
    k2ma[i][1] = 0.;
  }

  fftw_execute (p_forward_k2a);
  fftw_execute (p_forward_k2b);

  for (int i = 0; i < nele; ++i){
    k2ma[i][0] += 2. * k2mb[i][0];
    k2ma[i][1] += 2. * k2mb[i][1];
  }

  free (gmi);
}


void ErrorEstimate_SPME_Ik::
calError (const DensityProfile_PiecewiseConst & dp)
{
  double scalor = volume/double(nele);
  for (int i = 0; i < nele; ++i){
    rho1[i][0] = dp.getProfile1(i) * scalor;
    rho1[i][1] = 0;
    rho2[i][0] = dp.getProfile2(i) * scalor;
    rho2[i][1] = 0;
  }
  fftw_execute (p_forward_rho1);
  fftw_execute (p_forward_rho2);

  array_multiply (error1x, nele, k1mx, rho1);
  array_multiply (error1y, nele, k1my, rho1);
  array_multiply (error1z, nele, k1mz, rho1);
  array_multiply (error2, nele, k2ma, rho2);

  fftw_execute (p_backward_error1x);
  fftw_execute (p_backward_error1y);
  fftw_execute (p_backward_error1z);
  fftw_execute (p_backward_error2);

  for (int i = 0; i < nele; ++i){
    error1x[i][0] /= volume;
    error1y[i][0] /= volume;
    error1z[i][0] /= volume;
    error2[i][0] /= volume;
  }
  for (int i = 0; i < nele; ++i){
    error1[i][0] =
	error1x[i][0] * error1x[i][0] + 
	error1y[i][0] * error1y[i][0] + 
	error1z[i][0] * error1z[i][0] ;
    error1[i][1] =
	error1x[i][1] * error1x[i][1] + 
	error1y[i][1] * error1y[i][1] + 
	error1z[i][1] * error1z[i][1] ;
  }
}


void ErrorEstimate_SPME_Ik::    
print_error (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  fprintf (fp, "# 1      2         3-4       5-6\n");
  fprintf (fp, "# x  RMS-E  E_inhomo^2  E_homo^2\n");
  for (int i = 0; i < K.x; ++i){
    IntVectorType idx;
    idx.x = i;
    idx.y = 0;
    idx.z = 0;
    unsigned index = index3to1 (idx, K);
    fprintf (fp, "%f %e   %e %e %e %e\n",
	     (i + 0.5) * vecA.xx / K.x,
	     sqrt(error1[index][0] + error2[index][0]),
	     error1[index][0],
	     error1[index][1],
	     error2[index][0],
	     error2[index][1]
	);
  }
  fclose (fp);
}


void ErrorEstimate_SPME_Ik::
print_meanf (const std::string & file,
	     const DensityProfile_PiecewiseConst & dp) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  fprintf (fp, "# 1       2    3-5\n");
  fprintf (fp, "# x  meanFx  meanF\n");
  for (int i = 0; i < K.x; ++i){
    IntVectorType idx;
    idx.x = i;
    idx.y = 0;
    idx.z = 0;
    unsigned index = index3to1 (idx, K);
    double scalor = 1.;
    scalor = 
    	(dp.getProfileP(index) + dp.getProfileN(index)) /
    	(dp.getProfileP(index) - dp.getProfileN(index)) ;
    scalor = 1.;
    
    fprintf (fp, "%f %e %e %e %e\n",
	     (i + 0.5) * vecA.xx / K.x,
	     error1x[index][0],
	     error1x[index][0] * scalor,
	     error1y[index][0] * scalor,
	     error1z[index][0] * scalor
	);
  }
  fclose (fp);
}



