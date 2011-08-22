#include <iostream>
#include "ErrorEstimate_SPME_Ana.h"

ErrorEstimate_SPME_Ana::
ErrorEstimate_SPME_Ana ()
    : malloced (false)
{
}

ErrorEstimate_SPME_Ana::
~ErrorEstimate_SPME_Ana()
{
  freeAll();
}

fftw_complex ** ErrorEstimate_SPME_Ana::
myMalloc (int numOfk,
	  int nele)
{
  fftw_complex ** result;
  result = (fftw_complex **) malloc (sizeof(fftw_complex*) * numOfk);
  for (int i = 0; i < numOfk; ++i){
    result[i] = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  }
  return result;
}

void ErrorEstimate_SPME_Ana::
myFree (int numOfk,
	fftw_complex ** a)
{
  for (int i = 0; i < numOfk; ++i){
    free (a[i]);
  }
  free (a);
}

void ErrorEstimate_SPME_Ana::
freeAll ()
{
  if (malloced){
    free (refined_error1x);
    free (refined_error1y);
    free (refined_error1z);
    free (refined_error1);
    free (rho1);
    myFree (numAlphak, f1kx);
    myFree (numAlphak, f1ky);
    myFree (numAlphak, f1kz);
    myFree (numAlphak, error1x);
    myFree (numAlphak, error1y);
    myFree (numAlphak, error1z);
    free (selfx);
    free (selfy);
    free (selfz);
    fftw_destroy_plan (p_forward_rho1);
    for (int i = 0; i < numAlphak; ++i){
      fftw_destroy_plan (p_backward_error1x[i]);
      fftw_destroy_plan (p_backward_error1y[i]);
      fftw_destroy_plan (p_backward_error1z[i]);
    }
    free (p_backward_error1x);
    free (p_backward_error1y);
    free (p_backward_error1z);
  }
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


void ErrorEstimate_SPME_Ana::
reinit (const double & beta_,
	const int & order_,
	const DensityProfile_PiecewiseConst & dp,
	const IntVectorType refine_)
{
  freeAll();
  
  spmeik.reinit (beta_, order_, dp);

  K = spmeik.K;
  nele = K.x * K.y * K.z;
  beta = spmeik.beta;
  order = spmeik.order;
  vecA = spmeik.vecA;
  vecAStar = spmeik.vecAStar;
  volume = spmeik.volume;

  refine = refine_;
  refined_K.x = refine.x * K.x;
  refined_K.y = refine.y * K.y;
  refined_K.z = refine.z * K.z;
  kcut = 3;
  numk = 2 * kcut + 1;
  numAlphak = 1 * numk; // compact store for orthogonal box

  int refinednele = refined_K.x * refined_K.y * refined_K.z;
  refined_error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1  = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);

  rho1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  f1kx = myMalloc (numAlphak, nele);
  f1ky = myMalloc (numAlphak, nele);
  f1kz = myMalloc (numAlphak, nele);
  error1x = myMalloc (numAlphak, nele);
  error1y = myMalloc (numAlphak, nele);
  error1z = myMalloc (numAlphak, nele);
  selfx = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  selfy = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  selfz = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  
  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_error1x = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_error1y = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_error1z = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  for (int i = 0; i < numAlphak; ++i){
    p_backward_error1x[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,error1x[i],error1x[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_error1y[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,error1y[i],error1y[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_error1z[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,error1z[i],error1z[i],FFTW_BACKWARD,FFTW_MEASURE);
  }

  malloced = true;
  
  calKernel();
}

void ErrorEstimate_SPME_Ana::
calKernel ()
{
  double scalor = 1./(2. * M_PI * volume);
  //
  // clean X
  //
  // bool cleanX = ( ((K.x >> 1) << 1) == K.x);
  // bool cleanY = ( ((K.y >> 1) << 1) == K.y);
  // bool cleanZ = ( ((K.z >> 1) << 1) == K.z);
  VectorType Ki;
  Ki.x = 1./double(K.x);
  Ki.y = 1./double(K.y);
  IntVectorType idx;
  IntVectorType posi;

  for (int kk = 0; kk < numAlphak; ++kk){
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
	  VectorType fenmu;
	  fenmu.x = 1./calsum (uu.x, order);
	  fenmu.y = 1./calsum (uu.y, order);
	  fenmu.z = 1./calsum (uu.z, order);
	  unsigned index = index3to1 (idx, K);
	  
	  f1kx[myk][index][0] =
	      -2. * fm * 1./gsl_pow_int(uu.x + 2.*M_PI*kk, order) * fenmu.x;
	  f1kx[myk][index][1] = 0.;
	  f1ky[myk][index][0] =
	      -2. * fm * 1./gsl_pow_int(uu.y + 2.*M_PI*kk, order) * fenmu.y;
	  f1ky[myk][index][1] = 0.;
	  f1kz[myk][index][0] =
	      -2. * fm * 1./gsl_pow_int(uu.z + 2.*M_PI*kk, order) * fenmu.z;
	  f1kz[myk][index][1] = 0.;
	  
	  selfx[myk][0] += f1kx[myk][index][0];
	  selfy[myk][0] += f1ky[myk][index][0];
	  selfz[myk][0] += f1kz[myk][index][0];
	  
	  f1kx[myk][index][0] *= volume;
	  f1ky[myk][index][0] *= volume;
	  f1kz[myk][index][0] *= volume;
	}
      }
    }
  }
}

// static void 
// array_multiply (fftw_complex * a,
// 		const int nele,
// 		fftw_complex * b,
// 		fftw_complex * c) 
// {
//   for (int ii = 0; ii < nele; ++ii){
//     // double tmpr, tmpi;
//     a[ii][0] =
// 	c[ii][0] * b[ii][0] - c[ii][1] * b[ii][1];
//     a[ii][1] =
// 	c[ii][0] * b[ii][1] + c[ii][1] * b[ii][0];
//   }
// }


static void 
array_multiply (fftw_complex ** a,
		const int kcut,
		const int nele,
		fftw_complex ** b,
		fftw_complex * c) 
{
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (int ii = 0; ii < nele; ++ii){
      // double tmpr, tmpi;
      a[myk][ii][0] =
	  c[ii][0] * b[myk][ii][0] - c[ii][1] * b[myk][ii][1];
      a[myk][ii][1] =
	  c[ii][0] * b[myk][ii][1] + c[ii][1] * b[myk][ii][0];
    }
  }
}

static void
multiply (fftw_complex & c,
	  const fftw_complex & a,
	  const fftw_complex & b)
{
  c[0] = a[0] * b[0] - a[1] * b[1];
  c[1] = a[0] * b[1] + a[1] * b[0];
}

void ErrorEstimate_SPME_Ana::
calError (const DensityProfile_PiecewiseConst & dp)
{
  spmeik.calError(dp);

  double scalor = volume/double(nele);
  for (int i = 0; i < nele; ++i){
    rho1[i][0] = dp.getProfile1(i) * scalor;
    rho1[i][1] = 0;
    // rho2[i][0] = dp.getProfile2(i) * scalor;
    // rho2[i][1] = 0;
  }
  fftw_execute (p_forward_rho1);
  // fftw_execute (p_forward_rho2);

  array_multiply (error1x, kcut, nele, f1kx, rho1);
  array_multiply (error1y, kcut, nele, f1ky, rho1);
  array_multiply (error1z, kcut, nele, f1kz, rho1);
  // array_multiply (error2, nele, k2ma, rho2);
  
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    fftw_execute (p_backward_error1x[myk]);
    fftw_execute (p_backward_error1y[myk]);
    fftw_execute (p_backward_error1z[myk]);
  }
  // fftw_execute (p_backward_error2);

  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (int i = 0; i < nele; ++i){
      error1x[myk][i][0] /= volume;
      error1y[myk][i][0] /= volume;
      error1z[myk][i][0] /= volume;
      error1x[myk][i][1] /= volume;
      error1y[myk][i][1] /= volume;
      error1z[myk][i][1] /= volume;
      // error2[i][0] /= volume;
    }
  }
  
  IntVectorType ii;
  for (ii.x = 0; ii.x < refined_K.x; ++ ii.x){
    for (ii.y = 0; ii.y < refined_K.y; ++ ii.y){
      for (ii.z = 0; ii.z < refined_K.z; ++ ii.z){
	VectorType ir;
	ir.x = ii.x / double(refined_K.x) * vecA.xx;
	ir.y = ii.y / double(refined_K.y) * vecA.yy;
	ir.z = ii.z / double(refined_K.z) * vecA.zz;
	IntVectorType jj ;
	jj.x = ii.x / refine.x;
	jj.y = ii.y / refine.y;
	jj.z = ii.z / refine.z;
	IntVectorType jneighbor;
	VectorType dr;
	if (double(ii.x) / double(refine.x) - jj.x < 0.5){
	  jneighbor.x = jj.x - 1;
	  if (jneighbor.x < 0) jneighbor.x += K.x;
	  dr.x = -vecA.xx / double(K.x);
	}
	else {
	  jneighbor.x = jj.x + 1;
	  if (jneighbor.x >= K.x) jneighbor.x -= K.x;
	  dr.x = vecA.xx / double(K.x);
	}
	if (double(ii.y) / double(refine.y) - jj.y < 0.5){
	  jneighbor.y = jj.y - 1;
	  if (jneighbor.y < 0) jneighbor.y += K.y;
	  dr.y = -vecA.yy / double(K.y);
	}
	else {
	  jneighbor.y = jj.y + 1;
	  if (jneighbor.y >= K.y) jneighbor.y -= K.y;
	  dr.y = vecA.yy / double(K.y);
	}
	if (double(ii.z) / double(refine.z) - jj.z < 0.5){
	  jneighbor.z = jj.z - 1;
	  if (jneighbor.z < 0) jneighbor.z += K.z;
	  dr.z = -vecA.zz / double(K.z);
	}
	else {
	  jneighbor.z = jj.z + 1;
	  if (jneighbor.z >= K.z) jneighbor.z -= K.z;
	  dr.z = vecA.zz / double(K.z);
	}
	int indexii, indexjj, indexjjneighbor;
	indexii = index3to1 (ii, refined_K);
	indexjj = index3to1 (jj, K);
	indexjjneighbor = index3to1 (jneighbor, K);
	VectorType result;
	result.x =
	    (spmeik.error1x[indexjjneighbor][0] - spmeik.error1x[indexjj][0]) /
	    dr.x *
	    (ir.x - (jj.x+0.5) / double(K.x) * vecA.xx) + spmeik.error1x[indexjj][0];
	result.y =
	    (spmeik.error1y[indexjjneighbor][0] - spmeik.error1y[indexjj][0]) /
	    dr.y *
	    (ir.y - (jj.y+0.5) / double(K.y) * vecA.yy) + spmeik.error1y[indexjj][0];
	result.z =
	    (spmeik.error1z[indexjjneighbor][0] - spmeik.error1z[indexjj][0]) /
	    dr.z *
	    (ir.z - (jj.z+0.5) / double(K.z) * vecA.zz) + spmeik.error1z[indexjj][0];
	  
	// refined_error1x[indexii][0] = 0.;
	// refined_error1y[indexii][0] = 0.;
	// refined_error1z[indexii][0] = 0.;
	refined_error1x[indexii][0] = result.x;
	refined_error1y[indexii][0] = result.y;
	refined_error1z[indexii][0] = result.z;
	// refined_error1x[indexii][0] = spmeik.error1x[indexjj][0];
	// refined_error1y[indexii][0] = spmeik.error1y[indexjj][0];
	// refined_error1z[indexii][0] = spmeik.error1z[indexjj][0];
	refined_error1x[indexii][1] = 0.;
	refined_error1y[indexii][1] = 0.;
	refined_error1z[indexii][1] = 0.;
      }
    }
  }
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (ii.x = 0; ii.x < refined_K.x; ++ ii.x){
      for (ii.y = 0; ii.y < refined_K.y; ++ ii.y){
	for (ii.z = 0; ii.z < refined_K.z; ++ ii.z){
	  VectorType rr;
	  rr.x = vecA.xx / double(refined_K.x) * ii.x;
	  rr.y = vecA.yy / double(refined_K.y) * ii.y;
	  rr.z = vecA.zz / double(refined_K.z) * ii.z;
	  VectorType value ;
	  value.x = 2 * M_PI * kk * K.x * (vecAStar.xx * rr.x);
	  value.y = 2 * M_PI * kk * K.y * (vecAStar.yy * rr.y);
	  value.z = 2 * M_PI * kk * K.z * (vecAStar.zz * rr.z);
	  fftw_complex fx, fy, fz;
	  fx[0] =-2 * M_PI * kk * K.x * vecAStar.xx * sin(value.x);
	  fx[1] = 2 * M_PI * kk * K.x * vecAStar.xx * cos(value.x);
	  fy[0] =-2 * M_PI * kk * K.y * vecAStar.yy * sin(value.y);
	  fy[1] = 2 * M_PI * kk * K.y * vecAStar.yy * cos(value.y);
	  fz[0] =-2 * M_PI * kk * K.z * vecAStar.zz * sin(value.z);
	  fz[1] = 2 * M_PI * kk * K.z * vecAStar.zz * cos(value.z);
	  // fy[0] = cos(value.y);
	  // fy[1] = sin(value.y);
	  // fz[0] = cos(value.z);
	  // fz[1] = sin(value.z);

	  VectorType ir;
	  ir.x = ii.x / double(refined_K.x) * vecA.xx;
	  ir.y = ii.y / double(refined_K.y) * vecA.yy;
	  ir.z = ii.z / double(refined_K.z) * vecA.zz;
	  IntVectorType jj ;
	  jj.x = ii.x / refine.x;
	  jj.y = ii.y / refine.y;
	  jj.z = ii.z / refine.z;
	  IntVectorType jneighbor;
	  VectorType dr;
	  if (double(ii.x) / double(refine.x) - jj.x < 0.5){
	    jneighbor.x = jj.x - 1;
	    if (jneighbor.x < 0) jneighbor.x += K.x;
	    dr.x = -vecA.xx / double(K.x);
	  }
	  else {
	    jneighbor.x = jj.x + 1;
	    if (jneighbor.x >= K.x) jneighbor.x -= K.x;
	    dr.x = vecA.xx / double(K.x);
	  }
	  if (double(ii.y) / double(refine.y) - jj.y < 0.5){
	    jneighbor.y = jj.y - 1;
	    if (jneighbor.y < 0) jneighbor.y += K.y;
	    dr.y = -vecA.yy / double(K.y);
	  }
	  else {
	    jneighbor.y = jj.y + 1;
	    if (jneighbor.y >= K.y) jneighbor.y -= K.y;
	    dr.y = vecA.yy / double(K.y);
	  }
	  if (double(ii.z) / double(refine.z) - jj.z < 0.5){
	    jneighbor.z = jj.z - 1;
	    if (jneighbor.z < 0) jneighbor.z += K.z;
	    dr.z = -vecA.zz / double(K.z);
	  }
	  else {
	    jneighbor.z = jj.z + 1;
	    if (jneighbor.z >= K.z) jneighbor.z -= K.z;
	    dr.z = vecA.zz / double(K.z);
	  }
	  int indexii, indexjj, indexjjneighbor;
	  indexii = index3to1 (ii, refined_K);
	  indexjj = index3to1 (jj, K);
	  indexjjneighbor = index3to1 (jneighbor, K);
	  fftw_complex resultx, resulty, resultz;
	  resultx[0] =
	      (error1x[myk][indexjjneighbor][0] - error1x[myk][indexjj][0]) /
	      dr.x *
	      (ir.x - (jj.x+0.5) / double(K.x) * vecA.xx) + error1x[myk][indexjj][0];
	  resulty[0] =
	      (error1y[myk][indexjjneighbor][0] - error1y[myk][indexjj][0]) /
	      dr.y *
	      (ir.y - (jj.y+0.5) / double(K.y) * vecA.yy) + error1y[myk][indexjj][0];
	  resultz[0] =
	      (error1z[myk][indexjjneighbor][0] - error1z[myk][indexjj][0]) /
	      dr.z *
	      (ir.z - (jj.z+0.5) / double(K.z) * vecA.zz) + error1z[myk][indexjj][0];
	  resultx[1] =
	      (error1x[myk][indexjjneighbor][1] - error1x[myk][indexjj][1]) /
	      dr.x *
	      (ir.x - (jj.x+0.5) / double(K.x) * vecA.xx) + error1x[myk][indexjj][1];
	  resulty[1] =
	      (error1y[myk][indexjjneighbor][1] - error1y[myk][indexjj][1]) /
	      dr.y *
	      (ir.y - (jj.y+0.5) / double(K.y) * vecA.yy) + error1y[myk][indexjj][1];
	  resultz[1] =
	      (error1z[myk][indexjjneighbor][1] - error1z[myk][indexjj][1]) /
	      dr.z *
	      (ir.z - (jj.z+0.5) / double(K.z) * vecA.zz) + error1z[myk][indexjj][1];
	  // IntVectorType jj ;
	  // jj.x = ii.x / refine.x;
	  // jj.y = ii.y / refine.y;
	  // jj.z = ii.z / refine.z;
	  // int indexjj;
	  // indexjj = index3to1 (jj, K);
	  // int indexii;
	  // indexii = index3to1 (ii, refined_K);
	  fftw_complex rx, ry, rz;
	  // multiply (rx, fx, error1x[myk][indexjj]);
	  // multiply (ry, fy, error1y[myk][indexjj]);
	  // multiply (rz, fz, error1z[myk][indexjj]);
	  multiply (rx, fx, resultx);
	  multiply (ry, fy, resulty);
	  multiply (rz, fz, resultz);
	  refined_error1x[indexii][0] -= rx[0];
	  refined_error1y[indexii][0] -= ry[0];
	  refined_error1z[indexii][0] -= rz[0];
	  refined_error1x[indexii][1] -= rx[1];
	  refined_error1y[indexii][1] -= ry[1];
	  refined_error1z[indexii][1] -= rz[1];
	  multiply (rx, fx, selfx[myk]);
	  multiply (ry, fy, selfy[myk]);
	  multiply (rz, fz, selfz[myk]);
	  refined_error1x[indexii][0] -= rx[0];
	  refined_error1y[indexii][0] -= ry[0];
	  refined_error1z[indexii][0] -= rz[0];
	  refined_error1x[indexii][1] -= rx[1];
	  refined_error1y[indexii][1] -= ry[1];
	  refined_error1z[indexii][1] -= rz[1];
	}
      }
    }
  }
}


void ErrorEstimate_SPME_Ana::
print_meanf (const std::string & file,
	     const DensityProfile_PiecewiseConst & dp) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  IntVectorType idx1;
  for (idx1.x = 0; idx1.x < refined_K.x; ++idx1.x){
    double sum0 = 0.;
    double sum1 = 0.;
    double sums = 0.;
    for (idx1.y = 0; idx1.y < refined_K.y; ++idx1.y){
      for (idx1.z = 0; idx1.z < refined_K.z; ++idx1.z){
	unsigned index1 = index3to1 (idx1, refined_K);
	sum0 += refined_error1x[index1][0];
	sum1 += refined_error1x[index1][1];
	IntVectorType idx2(idx1);
	idx2.x = idx1.x/ refine.x;
	idx2.y = idx1.y/ refine.y;
	idx2.z = idx1.z/ refine.z;
	unsigned index2 = index3to1 (idx2, K);
	double scalor = 1.;
	scalor = 
	    (dp.getProfileP(index2) + dp.getProfileN(index2)) /
	    (dp.getProfileP(index2) - dp.getProfileN(index2)) ;
	sums += refined_error1x[index1][0] * scalor;
      }
    }
    sum0 /= (refined_K.y * refined_K.z);
    sum1 /= (refined_K.y * refined_K.z);
    sums /= (refined_K.y * refined_K.z);
    // fprintf (fp, "%f %e %e %e\n",
    // 	     (i + 0.5) * vecA.xx / refined_K.x,
    // 	     refined_error1x[index1][0],
    // 	     refined_error1x[index1][1],
    // 	     refined_error1x[index1][0] * scalor
    // 	);
    fprintf (fp, "%f %e %e %e\n",
    	     (idx1.x + 0.5) * vecA.xx / refined_K.x,
	     sum0, sum1, sums);
  }
  fclose (fp);
}


  // for (ii.x = 0; ii.x < spmeik.K.x; ++ii.x){
  // for (ii.y = 0; ii.y < spmeik.K.y; ++ii.y){
  // for (ii.z = 0; ii.z < spmeik.K.z; ++ii.z){
  //   IntVectorType jj (ii);
  //   jj.x *= refine.x;
  //   jj.y *= refine.y;
  //   jj.z *= refine.z;
  //   unsigned indexii, indexjj;
  //   index3to1 (ii, spmeik.K, indexii);
  //   index3to1 (jj, refined_K, indexjj);

  //   for (int ll = 0; ll < refine.x; 
  //   refined_error1x[indexjj+0]
      
