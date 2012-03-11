#include <iostream>

#include "ErrorEstimate_SPME_St_H2O.h"

const int global_half_N_k_sum = 4;

ErrorEstimate_SPME_St_H2O::
ErrorEstimate_SPME_St_H2O ()
    : malloced (false)
{
}

ErrorEstimate_SPME_St_H2O::
~ErrorEstimate_SPME_St_H2O()
{
  freeAll();
}

void ErrorEstimate_SPME_St_H2O::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void ErrorEstimate_SPME_St_H2O::
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

void ErrorEstimate_SPME_St_H2O::
freeAll ()
{
  if (malloced){
    free (rho1);
    free (rho2);
    free (rhoO);
    free (rhoH);
    free (k1mx);
    free (k1my);
    free (k1mz);
    free (k1rx);
    free (k1ry);
    free (k1rz);
    free (k2ma);
    free (error1x);
    free (error1y);
    free (error1z);
    free (error1);
    free (error2);
    free (TO);
    free (TH);
    free (TOKx);
    free (TOKy);
    free (TOKz);
    free (THKx);
    free (THKy);
    free (THKz);
    free (coKO);
    free (coKH);
    free (coerrorO);
    free (coerrorH);
    fftw_destroy_plan (p_forward_rho1);
    fftw_destroy_plan (p_forward_rho2);
    fftw_destroy_plan (p_backward_k1mx);
    fftw_destroy_plan (p_backward_k1my);
    fftw_destroy_plan (p_backward_k1mz);
    fftw_destroy_plan (p_forward_k2a);
    fftw_destroy_plan (p_backward_error1x);
    fftw_destroy_plan (p_backward_error1y);
    fftw_destroy_plan (p_backward_error1z);
    fftw_destroy_plan (p_backward_error2);
    fftw_destroy_plan (p_backward_TOKx);
    fftw_destroy_plan (p_backward_TOKy);
    fftw_destroy_plan (p_backward_TOKz);
    fftw_destroy_plan (p_backward_THKx);
    fftw_destroy_plan (p_backward_THKy);
    fftw_destroy_plan (p_backward_THKz);
    fftw_destroy_plan (p_forward_coKO);
    fftw_destroy_plan (p_forward_coKH);
    fftw_destroy_plan (p_backward_coerrorO);
    fftw_destroy_plan (p_backward_coerrorH);
    fftw_destroy_plan (p_forward_rhoO);
    fftw_destroy_plan (p_forward_rhoH);
  }
}


void ErrorEstimate_SPME_St_H2O::
reinit (const double & beta_,
	const int & order_,
	// const IntVectorType & Kmax_,
	const DensityProfile_PiecewiseConst & dp,
	const std::vector<std::vector<double > > & moles,
	const std::vector<double > & charges)
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
  rhoO = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rhoH = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1mx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1my = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1mz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1rx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1ry = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1rz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k2ma = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  TO	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  TH	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  TOKx	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  TOKy	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  TOKz	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  THKx	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  THKy	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  THKz	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  coKO	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  coKH	= (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  coerrorO = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  coerrorH = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rho2 = fftw_plan_dft_3d (K.x, K.y, K.z, rho2, rho2, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rhoO = fftw_plan_dft_3d (K.x, K.y, K.z, rhoO, rhoO, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rhoH = fftw_plan_dft_3d (K.x, K.y, K.z, rhoH, rhoH, FFTW_FORWARD,  FFTW_MEASURE);
  
  p_backward_k1mx = fftw_plan_dft_3d (K.x, K.y, K.z, k1mx, k1rx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_k1my = fftw_plan_dft_3d (K.x, K.y, K.z, k1my, k1ry, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_k1mz = fftw_plan_dft_3d (K.x, K.y, K.z, k1mz, k1rz, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_k2a = fftw_plan_dft_3d (K.x, K.y, K.z, k2ma, k2ma, FFTW_FORWARD, FFTW_MEASURE);
  p_backward_error1x = fftw_plan_dft_3d (K.x, K.y, K.z, error1x, error1x, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1y = fftw_plan_dft_3d (K.x, K.y, K.z, error1y, error1y, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1z = fftw_plan_dft_3d (K.x, K.y, K.z, error1z, error1z, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error2 = fftw_plan_dft_3d (K.x, K.y, K.z, error2, error2, FFTW_BACKWARD, FFTW_MEASURE);

  p_backward_TOKx = fftw_plan_dft_3d (K.x, K.y, K.z, TOKx, TOKx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_TOKy = fftw_plan_dft_3d (K.x, K.y, K.z, TOKy, TOKy, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_TOKz = fftw_plan_dft_3d (K.x, K.y, K.z, TOKz, TOKz, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_THKx = fftw_plan_dft_3d (K.x, K.y, K.z, THKx, THKx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_THKy = fftw_plan_dft_3d (K.x, K.y, K.z, THKy, THKy, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_THKz = fftw_plan_dft_3d (K.x, K.y, K.z, THKz, THKz, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_coKO = fftw_plan_dft_3d (K.x, K.y, K.z, coKO, coKO, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_coKH = fftw_plan_dft_3d (K.x, K.y, K.z, coKH, coKH, FFTW_FORWARD, FFTW_MEASURE);
  p_backward_coerrorO = fftw_plan_dft_3d (K.x, K.y, K.z, coerrorO, coerrorO, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_coerrorH = fftw_plan_dft_3d (K.x, K.y, K.z, coerrorH, coerrorH, FFTW_BACKWARD, FFTW_MEASURE);  
  
  malloced = true;

  printf ("# reinit: start build TO and TH ...");
  fflush (stdout);
  calTOTH (moles, charges);
  printf (" done\n");  
  fflush (stdout);

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

void ErrorEstimate_SPME_St_H2O::
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

  // for (int i = 0; i < nele; ++i){
  //   k2ma[i][0] += 2. * k2mb[i][0];
  //   k2ma[i][1] += 2. * k2mb[i][1];
  // }

  for (int ii = 0; ii < nele; ++ii){
    TOKx[ii][0] = k1mx[ii][0] * TO[ii][0] - k1mx[ii][1] * TO[ii][1];
    TOKx[ii][1] = k1mx[ii][0] * TO[ii][1] + k1mx[ii][1] * TO[ii][0];
    TOKy[ii][0] = k1my[ii][0] * TO[ii][0] - k1my[ii][1] * TO[ii][1];
    TOKy[ii][1] = k1my[ii][0] * TO[ii][1] + k1my[ii][1] * TO[ii][0];
    TOKz[ii][0] = k1mz[ii][0] * TO[ii][0] - k1mz[ii][1] * TO[ii][1];
    TOKz[ii][1] = k1mz[ii][0] * TO[ii][1] + k1mz[ii][1] * TO[ii][0];

    THKx[ii][0] = k1mx[ii][0] * TH[ii][0] - k1mx[ii][1] * TH[ii][1];
    THKx[ii][1] = k1mx[ii][0] * TH[ii][1] + k1mx[ii][1] * TH[ii][0];
    THKy[ii][0] = k1my[ii][0] * TH[ii][0] - k1my[ii][1] * TH[ii][1];
    THKy[ii][1] = k1my[ii][0] * TH[ii][1] + k1my[ii][1] * TH[ii][0];
    THKz[ii][0] = k1mz[ii][0] * TH[ii][0] - k1mz[ii][1] * TH[ii][1];
    THKz[ii][1] = k1mz[ii][0] * TH[ii][1] + k1mz[ii][1] * TH[ii][0];
  }
  
  fftw_execute (p_backward_TOKx);
  fftw_execute (p_backward_TOKy);
  fftw_execute (p_backward_TOKz);

  fftw_execute (p_backward_THKx);
  fftw_execute (p_backward_THKy);
  fftw_execute (p_backward_THKz);
  
  for (int ii = 0; ii < nele; ++ii){
    TOKx[ii][0] /= volume;
    TOKx[ii][1] /= volume;
    TOKy[ii][0] /= volume;
    TOKy[ii][1] /= volume;
    TOKz[ii][0] /= volume;
    TOKz[ii][1] /= volume;

    THKx[ii][0] /= volume;
    THKx[ii][1] /= volume;
    THKy[ii][0] /= volume;
    THKy[ii][1] /= volume;
    THKz[ii][0] /= volume;
    THKz[ii][1] /= volume;
  }

  for (int ii = 0; ii < nele; ++ii){
    coKO[ii][0] =
	k1rx[ii][0] * TOKx[ii][0] - k1rx[ii][1] * TOKx[ii][1] +
	k1ry[ii][0] * TOKy[ii][0] - k1ry[ii][1] * TOKy[ii][1] +
	k1rz[ii][0] * TOKz[ii][0] - k1rz[ii][1] * TOKz[ii][1] ;
    coKO[ii][1] = 
	k1rx[ii][0] * TOKx[ii][1] + k1rx[ii][1] * TOKx[ii][0] +
	k1ry[ii][0] * TOKy[ii][1] + k1ry[ii][1] * TOKy[ii][0] +
	k1rz[ii][0] * TOKz[ii][1] + k1rz[ii][1] * TOKz[ii][0] ;
    coKH[ii][0] =
	k1rx[ii][0] * THKx[ii][0] - k1rx[ii][1] * THKx[ii][1] +
	k1ry[ii][0] * THKy[ii][0] - k1ry[ii][1] * THKy[ii][1] +
	k1rz[ii][0] * THKz[ii][0] - k1rz[ii][1] * THKz[ii][1] ;
    coKH[ii][1] = 
	k1rx[ii][0] * THKx[ii][1] + k1rx[ii][1] * THKx[ii][0] +
	k1ry[ii][0] * THKy[ii][1] + k1ry[ii][1] * THKy[ii][0] +
	k1rz[ii][0] * THKz[ii][1] + k1rz[ii][1] * THKz[ii][0] ;
    // printf ("%e %e    %e %e\n",
    // 	    coKO[ii][0], coKO[ii][1],
    // 	    coKH[ii][0], coKH[ii][1]);
  }

  fftw_execute (p_forward_coKO);
  fftw_execute (p_forward_coKH);
  
  for (int ii = 0; ii < nele; ++ii){
    coKO[ii][0] *= volume / double(nele);
    coKO[ii][1] *= volume / double(nele);
    coKH[ii][0] *= volume / double(nele);
    coKH[ii][1] *= volume / double(nele);
  }
  
  free (gmi);
}



void ErrorEstimate_SPME_St_H2O::
calTm (const std::vector<std::vector<double > > & dirs,
       const std::vector<double > & charges,
       fftw_complex * out)
{
  int nele = K.x * K.y * K.z;
  IntVectorType idx;
  IntVectorType posi;
  bool cleanX = ( ((K.x >> 1) << 1) == K.x);
  bool cleanY = ( ((K.y >> 1) << 1) == K.y);
  bool cleanZ = ( ((K.z >> 1) << 1) == K.z);

  for (int ii = 0; ii < nele; ++ii){
    out[ii][0] = out[ii][1] = 0.;
  }

  printf ("\n");
  for (unsigned dc = 0; dc < dirs.size(); ++dc){
    printf ("# reading %d th direction    \r", dc);
    fflush (stdout);
    for (idx.x = 0; idx.x < K.x; ++idx.x){
      if (idx.x > (K.x >> 1)) posi.x = idx.x - K.x;
      else posi.x = idx.x;
      for (idx.y = 0; idx.y < K.y; ++idx.y){
	if (idx.y > (K.y >> 1)) posi.y = idx.y - K.y;
	else posi.y = idx.y;
	for (idx.z = 0; idx.z < K.z; ++idx.z){
	  if (idx.z > (K.z >> 1)) posi.z = idx.z - K.z;
	  else posi.z = idx.z;
	  if (idx.x == 0 && idx.y == 0 && idx.z == 0) continue;
	  unsigned index = index3to1 (idx, K);
	  VectorType mm;
	  mm.x = posi.x * vecAStar.xx + posi.y * vecAStar.yx + posi.z * vecAStar.zx;
	  mm.y = posi.x * vecAStar.xy + posi.y * vecAStar.yy + posi.z * vecAStar.zy;
	  mm.z = posi.x * vecAStar.xz + posi.y * vecAStar.yz + posi.z * vecAStar.zz;
	  double value;
	  value = mm.x * dirs[dc][0] + mm.y * dirs[dc][1] + mm.z * dirs[dc][2];
	  out[index][0] += charges[dc] * cos(2. * M_PI * value);
	  out[index][1] += charges[dc] * sin(2. * M_PI * value);
	
	  if ((idx.x == (K.x >> 1) && cleanX) ||
	      (idx.y == (K.y >> 1) && cleanY) ||
	      (idx.z == (K.z >> 1) && cleanZ) ) {
	    out[index][0] = out[index][1] = 0.;
	  }
	}
      }
    }
  }
  printf ("\n");
  for (int ii = 0; ii < nele; ++ii){
    out[ii][0] /= double (dirs.size());
    out[ii][1] /= double (dirs.size());
  }
}

void ErrorEstimate_SPME_St_H2O::
calTOTH (const std::vector<std::vector<double > > & atoms,
	 const std::vector<double > & charges)
{
  int natom = atoms.size();
  if ((charges.size()) < natom){
    std::cerr << "wrong size of charge table" << std::endl;
    exit (1);
  }
  int nmol = natom / 3;
  std::vector<std::vector<double > > sh1(nmol), sh2(nmol), so(nmol), sh(nmol);
  std::vector<double > qh(nmol), qo(nmol);
  
  for (int ii = 0; ii < nmol; ++ii){
    int atomO = ii * 3;
    int atomH1 = ii * 3 + 1;
    int atomH2 = ii * 3 + 2;

    qh[ii] = charges[atomH1];
    qo[ii] = charges[atomO];

    std::vector<double > tmp(3);
    for (int dd = 0; dd < 3; ++dd){
      tmp[dd] = atoms[atomH1][dd] - atoms[atomO][dd];
    }
    sh1[ii] = tmp;
    for (int dd = 0; dd < 3; ++dd){
      tmp[dd] = atoms[atomH2][dd] - atoms[atomO][dd];
    }
    sh2[ii] = tmp;
    
    for (int dd = 0; dd < 3; ++dd){
      tmp[dd] = atoms[atomO][dd] - atoms[atomH1][dd];
    }
    so[ii] = tmp;
    for (int dd = 0; dd < 3; ++dd){
      tmp[dd] = atoms[atomH2][dd] - atoms[atomH1][dd];
    }
    sh[ii] = tmp;
  }

  for (int ii = 0; ii < nele; ++ii){
    TO[ii][0] = TO[ii][1] = 0.;
    TH[ii][0] = TH[ii][1] = 0.;
  }
  
  fftw_complex * buff = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  calTm (sh1, qh, buff);
  for (int ii = 0; ii < nele; ++ii){
    TO[ii][0] += buff[ii][0];
    TO[ii][1] += buff[ii][1];
  }
  calTm (sh2, qh, buff);
  for (int ii = 0; ii < nele; ++ii){
    TO[ii][0] += buff[ii][0];
    TO[ii][1] += buff[ii][1];
  }
  
  calTm (so, qo, buff);
  for (int ii = 0; ii < nele; ++ii){
    TH[ii][0] += buff[ii][0];
    TH[ii][1] += buff[ii][1];
  }
  calTm (sh, qh, buff);
  for (int ii = 0; ii < nele; ++ii){
    TH[ii][0] += buff[ii][0];
    TH[ii][1] += buff[ii][1];
  }

  free (buff);
}



void ErrorEstimate_SPME_St_H2O::
calError (const DensityProfile_PiecewiseConst & dp,
	  const double charge)
{
  double scalor = volume/double(nele);
  for (int i = 0; i < nele; ++i){
    rho1[i][0] = dp.getProfile1(i) * scalor;
    rho1[i][1] = 0;
    rho2[i][0] = dp.getProfile2(i) * scalor;
    rho2[i][1] = 0;
    rhoO[i][0] = dp.getProfileP(i) * scalor;
    rhoO[i][1] = 0;
    rhoH[i][0] = dp.getProfileN(i) * scalor;
    rhoH[i][1] = 0;
  }
  fftw_execute (p_forward_rho1);
  fftw_execute (p_forward_rho2);
  fftw_execute (p_forward_rhoO);
  fftw_execute (p_forward_rhoH);

  array_multiply (error1x, nele, k1mx, rho1);
  array_multiply (error1y, nele, k1my, rho1);
  array_multiply (error1z, nele, k1mz, rho1);
  array_multiply (error2, nele, k2ma, rho2);
  array_multiply (coerrorO, nele, coKO, rhoO);
  array_multiply (coerrorH, nele, coKH, rhoH);

  fftw_execute (p_backward_error1x);
  fftw_execute (p_backward_error1y);
  fftw_execute (p_backward_error1z);
  fftw_execute (p_backward_error2);
  fftw_execute (p_backward_coerrorO);
  fftw_execute (p_backward_coerrorH);

  for (int i = 0; i < nele; ++i){
    // error1x[i][0] /= volume;
    // error1y[i][0] /= volume;
    // error1z[i][0] /= volume;
    // error2[i][0] /= volume;
    error1x[i][0] *= charge / volume;
    error1y[i][0] *= charge / volume;
    error1z[i][0] *= charge / volume;
    error2[i][0] *= charge * charge / volume;
    coerrorO[i][0] *= charge * charge / volume;
    coerrorH[i][0] *= charge * charge / volume;
    coerrorO[i][1] *= charge * charge / volume;
    coerrorH[i][1] *= charge * charge / volume;
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


void ErrorEstimate_SPME_St_H2O::    
print_error (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  fprintf (fp, "# 1      2         3-4       5-6     7-8\n");
  fprintf (fp, "# x  RMS-E  E_inhomo^2  E_homo^2  E_corr\n");
  for (int i = 0; i < K.x; ++i){
    IntVectorType idx;
    idx.x = i;
    idx.y = 0;
    idx.z = 0;
    unsigned index = index3to1 (idx, K);
    fprintf (fp, "%f %e   %e %e   %e %e   %e %e\n",
	     (i + 0.5) * vecA.xx / K.x,
	     sqrt(error1[index][0] + error2[index][0] +
		  coerrorO[index][0] + coerrorH[index][0]),
	     error1[index][0],
	     error1[index][1],
	     error2[index][0],
	     error2[index][1],
	     coerrorO[index][0] + coerrorH[index][0],
	     coerrorO[index][1] + coerrorH[index][1]
	);
  }
  fclose (fp);
}


void ErrorEstimate_SPME_St_H2O::
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



