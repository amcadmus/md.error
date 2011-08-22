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
    free (rho1);
    free (refined_error1x);
    free (refined_error1y);
    free (refined_error1z);
    free (refined_error1);
    myFree (numAlphak, k1mx);
    myFree (numAlphak, k1my);
    myFree (numAlphak, k1mz);
    myFree (numAlphak, k1rx);
    myFree (numAlphak, k1ry);
    myFree (numAlphak, k1rz);
    myFree (numAlphak, error1x);
    myFree (numAlphak, error1y);
    myFree (numAlphak, error1z);
    free (selfx);
    free (selfy);
    free (selfz);
    fftw_destroy_plan (p_forward_rho1);
    for (int i = 0; i < numAlphak; ++i){
      fftw_destroy_plan (p_backward_k1mx[i]);
      fftw_destroy_plan (p_backward_k1my[i]);
      fftw_destroy_plan (p_backward_k1mz[i]);
      fftw_destroy_plan (p_backward_error1x[i]);
      fftw_destroy_plan (p_backward_error1y[i]);
      fftw_destroy_plan (p_backward_error1z[i]);
    }
    free (p_backward_k1mx);
    free (p_backward_k1my);
    free (p_backward_k1mz);
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
  beta = spmeik.beta;
  order = spmeik.order;
  K = spmeik.K;
  vecA = spmeik.vecA;
  vecAStar = spmeik.vecAStar;
  volume = spmeik.volume;
  refine = refine_;
  refined_K.x = refine.x * K.x;
  refined_K.y = refine.y * K.y;
  refined_K.z = refine.z * K.z;

  int refinednele = refined_K.x * refined_K.y * refined_K.z;
  refined_error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1  = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);

  nele = K.x * K.y * K.z;
  kcut = 3;
  numk = 2 * kcut + 1;
  numAlphak = 1 * numk; // compact store for orthogonal box
  rho1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  k1mx = myMalloc (numAlphak, nele);
  k1my = myMalloc (numAlphak, nele);
  k1mz = myMalloc (numAlphak, nele);
  k1rx = myMalloc (numAlphak, nele);
  k1ry = myMalloc (numAlphak, nele);
  k1rz = myMalloc (numAlphak, nele);
  error1x = myMalloc (numAlphak, nele);
  error1y = myMalloc (numAlphak, nele);
  error1z = myMalloc (numAlphak, nele);
  
  selfx = (VectorType *) malloc (sizeof(VectorType) * (numAlphak));
  selfy = (VectorType *) malloc (sizeof(VectorType) * (numAlphak));
  selfz = (VectorType *) malloc (sizeof(VectorType) * (numAlphak));

  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);

  p_backward_k1mx = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_k1my = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_k1mz = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_error1x = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_error1y = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_error1z = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  for (int i = 0; i < numAlphak; ++i){
    p_backward_k1mx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z, k1mx[i], k1rx[i], FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_k1my[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z, k1my[i], k1ry[i], FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_k1mz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z, k1mz[i], k1rz[i], FFTW_BACKWARD,FFTW_MEASURE);
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
  // bool cleanX = ( ((K.x >> 1) << 1) == K.x);
  // bool cleanY = ( ((K.y >> 1) << 1) == K.y);
  // bool cleanZ = ( ((K.z >> 1) << 1) == K.z);
  VectorType Ki;
  Ki.x = 1./double(K.x);
  Ki.y = 1./double(K.y);
  Ki.z = 1./double(K.z);

  for (int kk = -kcut; kk <= kcut; ++kk){
    int myk = kk + kcut;
    selfx[myk].x = selfx[myk].y = selfx[myk].z = 0.;
    selfy[myk].x = selfy[myk].y = selfy[myk].z = 0.;
    selfz[myk].x = selfz[myk].y = selfz[myk].z = 0.;
  }
  
  IntVectorType idx;
  IntVectorType posi;
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
	  if (m2 == 0) {
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
	  k1mx[myk][index][0] = 0;
	  k1mx[myk][index][1] = (
	      -2. * fm * fenmu.x / gsl_pow_int(uu.x + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.x * vecAStar.xx );
	  k1my[myk][index][0] = 0;
	  k1my[myk][index][1] = (
	      -2. * fm * fenmu.y / gsl_pow_int(uu.y + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.y * vecAStar.yy );
	  k1mz[myk][index][0] = 0;
	  k1mz[myk][index][1] = (
	      -2. * fm * fenmu.z / gsl_pow_int(uu.z + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.z * vecAStar.zz );
	  k1mx[myk][index][1] *= volume;
	  k1my[myk][index][1] *= volume;
	  k1mz[myk][index][1] *= volume;
	  
	  selfx[myk].x +=
	      2. * fm * fenmu.x / gsl_pow_int(uu.x + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.x * vecAStar.xx;
	  selfx[myk].y +=
	      2. * fm * fenmu.y / gsl_pow_int(uu.y + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.y * vecAStar.xy;
	  selfx[myk].z +=
	      2. * fm * fenmu.z / gsl_pow_int(uu.z + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.z * vecAStar.xz;
	  selfy[myk].x +=
	      2. * fm * fenmu.x / gsl_pow_int(uu.x + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.x * vecAStar.yx;
	  selfy[myk].y +=
	      2. * fm * fenmu.y / gsl_pow_int(uu.y + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.y * vecAStar.yy;
	  selfy[myk].z +=
	      2. * fm * fenmu.z / gsl_pow_int(uu.z + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.z * vecAStar.yz;
	  selfz[myk].x +=
	      2. * fm * fenmu.x / gsl_pow_int(uu.x + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.x * vecAStar.zx;
	  selfz[myk].y +=
	      2. * fm * fenmu.y / gsl_pow_int(uu.y + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.y * vecAStar.zy;
	  selfz[myk].z +=
	      2. * fm * fenmu.z / gsl_pow_int(uu.z + 2.*M_PI * kk, order) *
	      2. * M_PI * kk * K.z * vecAStar.zz;
	}
      }
    }
  }
	
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

  array_multiply (error1x, kcut, nele, k1mx, rho1);
  array_multiply (error1y, kcut, nele, k1my, rho1);
  array_multiply (error1z, kcut, nele, k1mz, rho1);
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
	IntVectorType jj ;
	jj.x = ii.x / refine.x;
	jj.y = ii.y / refine.y;
	jj.z = ii.z / refine.z;
	int indexii, indexjj;
	indexii = index3to1 (ii, refined_K);
	indexjj = index3to1 (jj, K);
	// refined_error1x[indexii][0] = 0.;
	// refined_error1y[indexii][0] = 0.;
	// refined_error1z[indexii][0] = 0.;
	refined_error1x[indexii][0] = spmeik.error1x[indexjj][0];
	refined_error1y[indexii][0] = spmeik.error1y[indexjj][0];
	refined_error1z[indexii][0] = spmeik.error1z[indexjj][0];
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
	  IntVectorType jj ;
	  jj.x = ii.x / refine.x;
	  jj.y = ii.y / refine.y;
	  jj.z = ii.z / refine.z;
	  int indexii, indexjj;
	  indexii = index3to1 (ii, refined_K);
	  indexjj = index3to1 (jj, K);
	  VectorType rr;
	  rr.x = vecA.xx / double(refined_K.x) * ii.x;
	  rr.y = vecA.yy / double(refined_K.y) * ii.y;
	  rr.z = vecA.zz / double(refined_K.z) * ii.z;
	  IntVectorType value ;
	  value.x = 2 * M_PI * kk * K.x * (vecAStar.xx * rr.x);
	  value.y = 2 * M_PI * kk * K.y * (vecAStar.yy * rr.y);
	  value.z = 2 * M_PI * kk * K.z * (vecAStar.zz * rr.z);
	  fftw_complex fx, fy, fz;
	  fx[0] = cos(value.x);
	  fx[1] = sin(value.x);
	  fy[0] = cos(value.y);
	  fy[1] = sin(value.y);
	  fz[0] = cos(value.z);
	  fz[1] = sin(value.z);
	  fftw_complex rx, ry, rz;
	  multiply (rx, fx, error1x[myk][indexjj]);
	  multiply (ry, fy, error1y[myk][indexjj]);
	  multiply (rz, fz, error1z[myk][indexjj]);
	  refined_error1x[indexii][0] -= rx[0];
	  refined_error1y[indexii][0] -= ry[0];
	  refined_error1z[indexii][0] -= rz[0];
	  refined_error1x[indexii][1] -= rx[1];
	  refined_error1y[indexii][1] -= ry[1];
	  refined_error1z[indexii][1] -= rz[1];
	  refined_error1x[indexii][0] -= selfx[myk].x * -sin(value.x);
	  refined_error1y[indexii][0] -= selfy[myk].y * -sin(value.y);
	  refined_error1z[indexii][0] -= selfz[myk].z * -sin(value.z);
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

  for (int i = 0; i < refined_K.x; ++i){
    IntVectorType idx1;
    idx1.x = i;
    idx1.y = 0;
    idx1.z = 0;
    unsigned index1 = index3to1 (idx1, refined_K);
    IntVectorType idx2(idx1);
    idx2.x = i / refine.x;
    unsigned index2 = index3to1 (idx2, K);
    double scalor = 1.;
    scalor = 
    	(dp.getProfileP(index2) + dp.getProfileN(index2)) /
    	(dp.getProfileP(index2) - dp.getProfileN(index2)) ;
    fprintf (fp, "%f %e %e %e\n",
	     (i + 0.5) * vecA.xx / refined_K.x,
	     refined_error1x[index1][0],
	     refined_error1x[index1][1],
	     refined_error1x[index1][0] * scalor
	);
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
      
