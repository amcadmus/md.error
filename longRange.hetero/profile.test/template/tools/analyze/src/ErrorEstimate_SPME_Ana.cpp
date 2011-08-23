#include <iostream>
#include "ErrorEstimate_SPME_Ana.h"

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
array_multiply_kl (fftw_complex ** a,
		   const int nk,
		   const int nele,
		   fftw_complex ** b,
		   fftw_complex * c) 
{
  for (int kk = 0; kk < nk*nk; ++kk){
    for (int ii = 0; ii < nele; ++ii){
      // double tmpr, tmpi;
      a[kk][ii][0] =
	  c[ii][0] * b[kk][ii][0] - c[ii][1] * b[kk][ii][1];
      a[kk][ii][1] =
	  c[ii][0] * b[kk][ii][1] + c[ii][1] * b[kk][ii][0];
    }
  }
}


static void
multiply (fftw_complex & c,
	  const fftw_complex & a_,
	  const fftw_complex & b_)
{
  fftw_complex a;
  fftw_complex b;
  a[0] = a_[0];
  a[1] = a_[1];
  b[0] = b_[0];
  b[1] = b_[1];
  c[0] = a[0] * b[0] - a[1] * b[1];
  c[1] = a[0] * b[1] + a[1] * b[0];
}


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

    free (refined_termIx);
    free (refined_termIy);
    free (refined_termIz);
    free (refined_termIIx);
    free (refined_termIIy);
    free (refined_termIIz);
    
    free (refined_error1);
    free (refined_error2);

    free (rho1);
    free (rho2);
    
    myFree (numAlphak, Fmx);
    myFree (numAlphak, Fmy);
    myFree (numAlphak, Fmz);
    myFree (numAlphak, FConvRho1x);
    myFree (numAlphak, FConvRho1y);
    myFree (numAlphak, FConvRho1z);
    myFree (numAlphak, Frx);
    myFree (numAlphak, Fry);
    myFree (numAlphak, Frz);
    
    myFree (numk*numk, FxFx);
    myFree (numk*numk, FyFy);
    myFree (numk*numk, FzFz);
    myFree (numk*numk, FxFxConvRho2);
    myFree (numk*numk, FyFyConvRho2);
    myFree (numk*numk, FzFzConvRho2);

    myFree (numk, Gmxx);
    myFree (numk, Gmxy);
    myFree (numk, Gmxz);
    myFree (numk, Gmyx);
    myFree (numk, Gmyy);
    myFree (numk, Gmyz);
    myFree (numk, Gmzx);
    myFree (numk, Gmzy);
    myFree (numk, Gmzz);

    myFree (numk*numk, GxxFx);
    myFree (numk*numk, GxyFy);
    myFree (numk*numk, GxzFz);
    myFree (numk*numk, GyxFx);
    myFree (numk*numk, GyyFy);
    myFree (numk*numk, GyzFz);
    myFree (numk*numk, GzxFx);
    myFree (numk*numk, GzyFy);
    myFree (numk*numk, GzzFz);
    myFree (numk*numk, GxxFxConvRho2);
    myFree (numk*numk, GxyFyConvRho2);
    myFree (numk*numk, GxzFzConvRho2);
    myFree (numk*numk, GyxFxConvRho2);
    myFree (numk*numk, GyyFyConvRho2);
    myFree (numk*numk, GyzFzConvRho2);
    myFree (numk*numk, GzxFxConvRho2);
    myFree (numk*numk, GzyFyConvRho2);
    myFree (numk*numk, GzzFzConvRho2);

    free (selfx);
    free (selfy);
    free (selfz);
    
    fftw_destroy_plan (p_forward_rho1);
    fftw_destroy_plan (p_forward_rho2);
    
    for (int i = 0; i < numAlphak; ++i){
      fftw_destroy_plan (p_backward_FConvRho1x[i]);
      fftw_destroy_plan (p_backward_FConvRho1y[i]);
      fftw_destroy_plan (p_backward_FConvRho1z[i]);

      fftw_destroy_plan (p_backward_Fmx[i]);
      fftw_destroy_plan (p_backward_Fmy[i]);
      fftw_destroy_plan (p_backward_Fmz[i]);
      fftw_destroy_plan (p_backward_Gmxx[i]);
      fftw_destroy_plan (p_backward_Gmxy[i]);
      fftw_destroy_plan (p_backward_Gmxz[i]);
      fftw_destroy_plan (p_backward_Gmyx[i]);
      fftw_destroy_plan (p_backward_Gmyy[i]);
      fftw_destroy_plan (p_backward_Gmyz[i]);
      fftw_destroy_plan (p_backward_Gmzx[i]);
      fftw_destroy_plan (p_backward_Gmzy[i]);
      fftw_destroy_plan (p_backward_Gmzz[i]);
    }
    for (int i = 0; i < numk*numk; ++i){
      fftw_destroy_plan (p_forward_FxFx[i]);
      fftw_destroy_plan (p_forward_FyFy[i]);
      fftw_destroy_plan (p_forward_FzFz[i]);
      fftw_destroy_plan (p_forward_GxxFx[i]);
      fftw_destroy_plan (p_forward_GxyFy[i]);
      fftw_destroy_plan (p_forward_GxzFz[i]);
      fftw_destroy_plan (p_forward_GyxFx[i]);
      fftw_destroy_plan (p_forward_GyyFy[i]);
      fftw_destroy_plan (p_forward_GyzFz[i]);
      fftw_destroy_plan (p_forward_GzxFx[i]);
      fftw_destroy_plan (p_forward_GzyFy[i]);
      fftw_destroy_plan (p_forward_GzzFz[i]);
    
      fftw_destroy_plan (p_backward_FxFxConvRho2[i]);
      fftw_destroy_plan (p_backward_FyFyConvRho2[i]);
      fftw_destroy_plan (p_backward_FzFzConvRho2[i]);
      fftw_destroy_plan (p_backward_GxxFxConvRho2[i]);
      fftw_destroy_plan (p_backward_GxyFyConvRho2[i]);
      fftw_destroy_plan (p_backward_GxzFzConvRho2[i]);
      fftw_destroy_plan (p_backward_GyxFxConvRho2[i]);
      fftw_destroy_plan (p_backward_GyyFyConvRho2[i]);
      fftw_destroy_plan (p_backward_GyzFzConvRho2[i]);
      fftw_destroy_plan (p_backward_GzxFxConvRho2[i]);
      fftw_destroy_plan (p_backward_GzyFyConvRho2[i]);
      fftw_destroy_plan (p_backward_GzzFzConvRho2[i]);
    }
    
    free (p_backward_FConvRho1x);
    free (p_backward_FConvRho1y);
    free (p_backward_FConvRho1z);

    free (p_backward_Fmx);
    free (p_backward_Fmy);
    free (p_backward_Fmz);
    free (p_backward_Gmxx);
    free (p_backward_Gmxy);
    free (p_backward_Gmxz);
    free (p_backward_Gmyx);
    free (p_backward_Gmyy);
    free (p_backward_Gmyz);
    free (p_backward_Gmzx);
    free (p_backward_Gmzy);
    free (p_backward_Gmzz);
    
    free (p_forward_FxFx);
    free (p_forward_FyFy);
    free (p_forward_FzFz);
    free (p_forward_GxxFx);
    free (p_forward_GxyFy);
    free (p_forward_GxzFz);
    free (p_forward_GyxFx);
    free (p_forward_GyyFy);
    free (p_forward_GyzFz);
    free (p_forward_GzxFx);
    free (p_forward_GzyFy);
    free (p_forward_GzzFz);
    
    free (p_backward_FxFxConvRho2);
    free (p_backward_FyFyConvRho2);
    free (p_backward_FzFzConvRho2);
    free (p_backward_GxxFxConvRho2);
    free (p_backward_GxyFyConvRho2);
    free (p_backward_GxzFzConvRho2);
    free (p_backward_GyxFxConvRho2);
    free (p_backward_GyyFyConvRho2);
    free (p_backward_GyzFzConvRho2);
    free (p_backward_GzxFxConvRho2);
    free (p_backward_GzyFyConvRho2);
    free (p_backward_GzzFzConvRho2);
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
  
  printf ("#@ here 0\n");
  fflush (stdout);

  spmeik.reinit (beta_, order_, dp);

  printf ("#@ here 1\n");
  fflush (stdout);
  
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
  kcut = 1;
  numk = 2 * kcut + 1;
  numAlphak = 1 * numk; // compact store for orthogonal box

  int refinednele = refined_K.x * refined_K.y * refined_K.z;
  refined_error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error1  = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_error2  = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);

  refined_termIx = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_termIy = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_termIz = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_termIIx = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_termIIy = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);
  refined_termIIz = (fftw_complex *) malloc (sizeof(fftw_complex) * refinednele);

  rho1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rho2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  Fmx = myMalloc (numAlphak, nele);
  Fmy = myMalloc (numAlphak, nele);
  Fmz = myMalloc (numAlphak, nele);
  FConvRho1x = myMalloc (numAlphak, nele);
  FConvRho1y = myMalloc (numAlphak, nele);
  FConvRho1z = myMalloc (numAlphak, nele);
  Frx = myMalloc (numAlphak, nele);
  Fry = myMalloc (numAlphak, nele);
  Frz = myMalloc (numAlphak, nele);

  FxFx = myMalloc (numk*numk, nele);
  FyFy = myMalloc (numk*numk, nele);
  FzFz = myMalloc (numk*numk, nele);
  FxFxConvRho2 = myMalloc (numk*numk, nele);
  FyFyConvRho2 = myMalloc (numk*numk, nele);
  FzFzConvRho2 = myMalloc (numk*numk, nele);

  Gmxx = myMalloc (numk, nele);
  Gmxy = myMalloc (numk, nele);
  Gmxz = myMalloc (numk, nele);
  Gmyx = myMalloc (numk, nele);
  Gmyy = myMalloc (numk, nele);
  Gmyz = myMalloc (numk, nele);
  Gmzx = myMalloc (numk, nele);
  Gmzy = myMalloc (numk, nele);
  Gmzz = myMalloc (numk, nele);

  GxxFx = myMalloc (numk*numk, nele);
  GxyFy = myMalloc (numk*numk, nele);
  GxzFz = myMalloc (numk*numk, nele);
  GyxFx = myMalloc (numk*numk, nele);
  GyyFy = myMalloc (numk*numk, nele);
  GyzFz = myMalloc (numk*numk, nele);
  GzxFx = myMalloc (numk*numk, nele);
  GzyFy = myMalloc (numk*numk, nele);
  GzzFz = myMalloc (numk*numk, nele);
  GxxFxConvRho2 = myMalloc (numk*numk, nele);
  GxyFyConvRho2 = myMalloc (numk*numk, nele);
  GxzFzConvRho2 = myMalloc (numk*numk, nele);
  GyxFxConvRho2 = myMalloc (numk*numk, nele);
  GyyFyConvRho2 = myMalloc (numk*numk, nele);
  GyzFzConvRho2 = myMalloc (numk*numk, nele);
  GzxFxConvRho2 = myMalloc (numk*numk, nele);
  GzyFyConvRho2 = myMalloc (numk*numk, nele);
  GzzFzConvRho2 = myMalloc (numk*numk, nele);
  
  selfx = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  selfy = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  selfz = (fftw_complex *) malloc (sizeof(fftw_complex) * (numAlphak));
  
  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rho2 = fftw_plan_dft_3d (K.x, K.y, K.z, rho2, rho2, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_FConvRho1x = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_FConvRho1y = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);
  p_backward_FConvRho1z = (fftw_plan *) malloc (sizeof(fftw_plan) * numAlphak);

  p_backward_Fmx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk);
  p_backward_Fmy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk);
  p_backward_Fmz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk);
  p_backward_Gmxx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmxy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmxz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmyx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmyy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmyz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmzx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmzy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk); 
  p_backward_Gmzz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk);
  
  p_forward_FxFx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_FyFy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_FzFz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GxxFx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GxyFy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GxzFz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GyxFx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GyyFy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GyzFz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GzxFx = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GzyFy = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_forward_GzzFz = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);

  p_backward_FxFxConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_FyFyConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_FzFzConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GxxFxConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GxyFyConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GxzFzConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GyxFxConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GyyFyConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GyzFzConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GzxFxConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GzyFyConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);
  p_backward_GzzFzConvRho2 = (fftw_plan *) malloc (sizeof(fftw_plan) * numk * numk);

  printf ("#@ here 2\n");
  fflush (stdout);

  for (int i = 0; i < numAlphak; ++i){
    p_backward_FConvRho1x[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FConvRho1x[i],FConvRho1x[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_FConvRho1y[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FConvRho1y[i],FConvRho1y[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_FConvRho1z[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FConvRho1z[i],FConvRho1z[i],FFTW_BACKWARD,FFTW_MEASURE);
    
    p_backward_Fmx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,Fmx[i],Frx[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Fmy[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,Fmy[i],Fry[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Fmz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,Fmz[i],Frz[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmxx[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmxx[i],Gmxx[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmxy[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmxy[i],Gmxy[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmxz[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmxz[i],Gmxz[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmyx[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmyx[i],Gmyx[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmyy[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmyy[i],Gmyy[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmyz[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmyz[i],Gmyz[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmzx[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmzx[i],Gmzx[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmzy[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmzy[i],Gmzy[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_Gmzz[i] = 
	fftw_plan_dft_3d (K.x,K.y,K.z,Gmzz[i],Gmzz[i],FFTW_BACKWARD,FFTW_MEASURE);
  }

  printf ("#@ here 3\n");
  fflush (stdout);

  for (int i = 0; i < numk*numk; ++i){
    p_forward_FxFx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FxFx[i],FxFx[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_FyFy[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FyFy[i],FyFy[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_FzFz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FzFz[i],FzFz[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GxxFx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxxFx[i],GxxFx[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GxyFy[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxyFy[i],GxyFy[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GxzFz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxzFz[i],GxzFz[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GyxFx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyxFx[i],GyxFx[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GyyFy[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyyFy[i],GyyFy[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GyzFz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyzFz[i],GyzFz[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GzxFx[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzxFx[i],GzxFx[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GzyFy[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzyFy[i],GzyFy[i],FFTW_FORWARD,FFTW_MEASURE);
    p_forward_GzzFz[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzzFz[i],GzzFz[i],FFTW_FORWARD,FFTW_MEASURE);

    p_backward_FxFxConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FxFxConvRho2[i],FxFxConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_FyFyConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FyFyConvRho2[i],FyFyConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_FzFzConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,FzFzConvRho2[i],FzFzConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GxxFxConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxxFxConvRho2[i],GxxFxConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GxyFyConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxyFyConvRho2[i],GxyFyConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GxzFzConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GxzFzConvRho2[i],GxzFzConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GyxFxConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyxFxConvRho2[i],GyxFxConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GyyFyConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyyFyConvRho2[i],GyyFyConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GyzFzConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GyzFzConvRho2[i],GyzFzConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GzxFxConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzxFxConvRho2[i],GzxFxConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GzyFyConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzyFyConvRho2[i],GzyFyConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    p_backward_GzzFzConvRho2[i] =
	fftw_plan_dft_3d (K.x,K.y,K.z,GzzFzConvRho2[i],GzzFzConvRho2[i],FFTW_BACKWARD,FFTW_MEASURE);
    
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
  Ki.z = 1./double(K.z);
  IntVectorType idx;
  IntVectorType posi;

  for (int kk = 0; kk < numAlphak; ++kk){
    selfx[kk][0] = selfx[kk][1] = 0.;
    selfy[kk][0] = selfy[kk][1] = 0.;
    selfz[kk][0] = selfz[kk][1] = 0.;
  }
  for (int index = 0; index < nele; ++index){
    Frx[kcut][index][0] = Frx[kcut][index][1] = 0.;
    Fry[kcut][index][0] = Fry[kcut][index][1] = 0.;
    Frz[kcut][index][0] = Frz[kcut][index][1] = 0.;
    Gmxx[kcut][index][0] = Gmxx[kcut][index][1] = 0.;
    Gmxy[kcut][index][0] = Gmxy[kcut][index][1] = 0.;
    Gmxz[kcut][index][0] = Gmxz[kcut][index][1] = 0.;
    Gmyx[kcut][index][0] = Gmyx[kcut][index][1] = 0.;
    Gmyy[kcut][index][0] = Gmyy[kcut][index][1] = 0.;
    Gmyz[kcut][index][0] = Gmyz[kcut][index][1] = 0.;
    Gmzx[kcut][index][0] = Gmzx[kcut][index][1] = 0.;
    Gmzy[kcut][index][0] = Gmzy[kcut][index][1] = 0.;
    Gmzz[kcut][index][0] = Gmzz[kcut][index][1] = 0.;
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
	  VectorType uu;
	  uu.x = 2. * M_PI * double(posi.x) * Ki.x;
	  uu.y = 2. * M_PI * double(posi.y) * Ki.y;
	  uu.z = 2. * M_PI * double(posi.z) * Ki.z;
	  double fm;
	  VectorType fenshu;
	  if (m2 == 0){
	    fm = 0.;
	  }
	  else {
	    fm = kernel_rm1_rec_f (m2, beta) * scalor;
	  }
	  fenshu.x = 1./(gsl_pow_int(uu.x + 2.*M_PI*kk, order) * calsum (uu.x, order));
	  fenshu.y = 1./(gsl_pow_int(uu.y + 2.*M_PI*kk, order) * calsum (uu.y, order));
	  fenshu.z = 1./(gsl_pow_int(uu.z + 2.*M_PI*kk, order) * calsum (uu.z, order));

	  // fenshu.x = fenshu.y = fenshu.z = 0;
	  // fenshu.x = 1./calsum (uu.x, order);
	  // fenshu.y = 1./calsum (uu.y, order);
	  // fenshu.z = 1./calsum (uu.z, order);
	  // fenshu.x *= 1./gsl_pow_int(uu.x + 2.*M_PI*kk, order);
	  // fenshu.y *= 1./gsl_pow_int(uu.y + 2.*M_PI*kk, order);
	  // fenshu.z *= 1./gsl_pow_int(uu.z + 2.*M_PI*kk, order);
	  
	  unsigned index = index3to1 (idx, K);
	  
	  Fmx[myk][index][0] = -2. * fm * fenshu.x;
	  Fmx[myk][index][1] = 0.;
	  Fmy[myk][index][0] = -2. * fm * fenshu.y;
	  Fmy[myk][index][1] = 0.;
	  Fmz[myk][index][0] = -2. * fm * fenshu.z;
	  Fmz[myk][index][1] = 0.;
	  printf("! %f\n", Fmz[0][0][0]);
	  
	  selfx[myk][0] += Fmx[myk][index][0];
	  selfy[myk][0] += Fmy[myk][index][0];
	  selfz[myk][0] += Fmz[myk][index][0];

	  Gmxx[myk][index][0] = 0.;
	  Gmxx[myk][index][1] =
	      -4. * M_PI * posi.x * vecAStar.xx * vecAStar.xx *
	      fm * fenshu.x;
	  Gmxy[myk][index][0] = 0.;
	  Gmxy[myk][index][1] =
	      -4. * M_PI * posi.y * vecAStar.yy * vecAStar.yy *
	      fm * fenshu.x;
	  Gmxz[myk][index][0] = 0.;
	  Gmxz[myk][index][1] =
	      -4. * M_PI * posi.z * vecAStar.zz * vecAStar.zz *
	      fm * fenshu.x;

	  Gmyx[myk][index][0] = 0.;
	  Gmyx[myk][index][1] =
	      -4. * M_PI * posi.x * vecAStar.xx * vecAStar.xx *
	      fm * fenshu.y;
	  Gmyy[myk][index][0] = 0.;
	  Gmyy[myk][index][1] =
	      -4. * M_PI * posi.y * vecAStar.yy * vecAStar.yy *
	      fm * fenshu.y;
	  Gmyz[myk][index][0] = 0.;
	  Gmyz[myk][index][1] =
	      -4. * M_PI * posi.z * vecAStar.zz * vecAStar.zz *
	      fm * fenshu.y;
	  
	  Gmzx[myk][index][0] = 0.;
	  Gmzx[myk][index][1] =
	      -4. * M_PI * posi.x * vecAStar.xx * vecAStar.xx *
	      fm * fenshu.z;
	  Gmzy[myk][index][0] = 0.;
	  Gmzy[myk][index][1] =
	      -4. * M_PI * posi.y * vecAStar.yy * vecAStar.yy *
	      fm * fenshu.z;
	  Gmzz[myk][index][0] = 0.;
	  Gmzz[myk][index][1] =
	      -4. * M_PI * posi.z * vecAStar.zz * vecAStar.zz *
	      fm * fenshu.z;
	}
      }
    }
  }

  for (int kk = -kcut; kk <= kcut; ++kk){  
    if (kk == 0) continue;
    int myk = kk + kcut;
    fftw_execute (p_backward_Fmx[myk]);   // Fmx -> Frx
    fftw_execute (p_backward_Fmy[myk]);   // 
    fftw_execute (p_backward_Fmz[myk]);   //
    fftw_execute (p_backward_Gmxx[myk]);  // Gmxx -> Grxx
    fftw_execute (p_backward_Gmxy[myk]);
    fftw_execute (p_backward_Gmxz[myk]);
    fftw_execute (p_backward_Gmyx[myk]);
    fftw_execute (p_backward_Gmyy[myk]);
    fftw_execute (p_backward_Gmyz[myk]);
    fftw_execute (p_backward_Gmzx[myk]);
    fftw_execute (p_backward_Gmzy[myk]);
    fftw_execute (p_backward_Gmzz[myk]);
  }

  printf("! %f\n", Fmz[0][0][0]);
  for (int kk = -kcut; kk <= kcut; ++kk){  
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (int index = 0; index < nele; ++index){
      Fmx[myk][index][0] *= volume;
      Fmy[myk][index][0] *= volume;
      Fmz[myk][index][0] *= volume;
    }
  }
  printf("! %f\n", Fmz[0][0][0]);

  // assemble cross terms...
  for (int kk = 0; kk < numk; ++kk){
    for (int ll = 0; ll < numk; ++ll){
      int klindex = kl2to1(kk, ll);
      for (int index = 0; index < nele; ++index){
	multiply (FxFx[klindex][index], Frx[kk][index], Frx[ll][index]);
	multiply (FyFy[klindex][index], Fry[kk][index], Fry[ll][index]);
	multiply (FzFz[klindex][index], Frz[kk][index], Frz[ll][index]);
	FxFx[klindex][index][0] *= volume / double(nele);
	FxFx[klindex][index][1] *= volume / double(nele);
	FyFy[klindex][index][0] *= volume / double(nele);
	FyFy[klindex][index][1] *= volume / double(nele);
	FzFz[klindex][index][0] *= volume / double(nele);
	FzFz[klindex][index][1] *= volume / double(nele);
	multiply (GxxFx[klindex][index], Gmxx[kk][index], Frx[ll][index]);
	multiply (GxyFy[klindex][index], Gmxy[kk][index], Fry[ll][index]);
	multiply (GxzFz[klindex][index], Gmxz[kk][index], Frz[ll][index]);
	multiply (GyxFx[klindex][index], Gmyx[kk][index], Frx[ll][index]);
	multiply (GyyFy[klindex][index], Gmyy[kk][index], Fry[ll][index]);
	multiply (GyzFz[klindex][index], Gmyz[kk][index], Frz[ll][index]);
	multiply (GzxFx[klindex][index], Gmzx[kk][index], Frx[ll][index]);
	multiply (GzyFy[klindex][index], Gmzy[kk][index], Fry[ll][index]);
	multiply (GzzFz[klindex][index], Gmzz[kk][index], Frz[ll][index]);
	GxxFx[klindex][index][0] *= volume / double(nele);
	GxyFy[klindex][index][0] *= volume / double(nele);
	GxzFz[klindex][index][0] *= volume / double(nele);
	GyxFx[klindex][index][0] *= volume / double(nele);
	GyyFy[klindex][index][0] *= volume / double(nele);
	GyzFz[klindex][index][0] *= volume / double(nele);
	GzxFx[klindex][index][0] *= volume / double(nele);
	GzyFy[klindex][index][0] *= volume / double(nele);
	GzzFz[klindex][index][0] *= volume / double(nele);
	GxxFx[klindex][index][1] *= volume / double(nele);
	GxyFy[klindex][index][1] *= volume / double(nele);
	GxzFz[klindex][index][1] *= volume / double(nele);
	GyxFx[klindex][index][1] *= volume / double(nele);
	GyyFy[klindex][index][1] *= volume / double(nele);
	GyzFz[klindex][index][1] *= volume / double(nele);
	GzxFx[klindex][index][1] *= volume / double(nele);
	GzyFy[klindex][index][1] *= volume / double(nele);
	GzzFz[klindex][index][1] *= volume / double(nele);
      }
      fftw_execute (p_forward_FxFx[klindex]);
      fftw_execute (p_forward_FyFy[klindex]);
      fftw_execute (p_forward_FzFz[klindex]);
      fftw_execute (p_forward_GxxFx[klindex]);
      fftw_execute (p_forward_GxyFy[klindex]);
      fftw_execute (p_forward_GxzFz[klindex]);
      fftw_execute (p_forward_GyxFx[klindex]);
      fftw_execute (p_forward_GyyFy[klindex]);
      fftw_execute (p_forward_GyzFz[klindex]);
      fftw_execute (p_forward_GzxFx[klindex]);
      fftw_execute (p_forward_GzyFy[klindex]);
      fftw_execute (p_forward_GzzFz[klindex]);
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


void ErrorEstimate_SPME_Ana::
interpolate (const IntVectorType ii,		// on the refined grid
	     const fftw_complex * value,	// on the coarse grid
	     fftw_complex & result)
{
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
  int indexjj, indexjjneighbor;
  indexjj = index3to1 (jj, K);
  indexjjneighbor = index3to1 (jneighbor, K);
  result[0] =
      (value[indexjjneighbor][0] - value[indexjj][0]) /
      dr.x *
      (ir.x - (jj.x+0.5) / double(K.x) * vecA.xx) + value[indexjj][0];
  result[1] =
      (value[indexjjneighbor][1] - value[indexjj][1]) /
      dr.x *
      (ir.x - (jj.x+0.5) / double(K.x) * vecA.xx) + value[indexjj][1];
}



void ErrorEstimate_SPME_Ana::
calError (const DensityProfile_PiecewiseConst & dp)
{
  spmeik.calError(dp);

  printf ("# here 0\n");
  fflush (stdout);

  double scalor = volume/double(nele);
  for (int i = 0; i < nele; ++i){
    rho1[i][0] = dp.getProfile1(i) * scalor;
    rho1[i][1] = 0;
    rho2[i][0] = dp.getProfile2(i) * scalor;
    rho2[i][1] = 0;
  }
  fftw_execute (p_forward_rho1);
  fftw_execute (p_forward_rho2);

  printf ("# here 1\n");
  fflush (stdout);
  
  array_multiply (FConvRho1x, kcut, nele, Fmx, rho1);
  array_multiply (FConvRho1y, kcut, nele, Fmy, rho1);
  array_multiply (FConvRho1z, kcut, nele, Fmz, rho1);
  printf("! %f\n", Fmz[0][0][0]);
  printf("! %f\n", Fmz[0][0][1]);
  printf("! %f\n", rho1[0][0]);
  printf("! %f\n", rho1[0][1]);
  printf("! %f\n", FConvRho1z[0][0][0]);

  printf ("# here 2\n");
  fflush (stdout);

  array_multiply_kl (FxFxConvRho2, numk, nele, FxFx, rho2);
  array_multiply_kl (FyFyConvRho2, numk, nele, FyFy, rho2);
  array_multiply_kl (FzFzConvRho2, numk, nele, FzFz, rho2);
  array_multiply_kl (GxxFxConvRho2, numk, nele, GxxFx, rho2);
  array_multiply_kl (GxyFyConvRho2, numk, nele, GxyFy, rho2);
  array_multiply_kl (GxzFzConvRho2, numk, nele, GxzFz, rho2);
  array_multiply_kl (GyxFxConvRho2, numk, nele, GyxFx, rho2);
  array_multiply_kl (GyyFyConvRho2, numk, nele, GyyFy, rho2);
  array_multiply_kl (GyzFzConvRho2, numk, nele, GyzFz, rho2);
  array_multiply_kl (GzxFxConvRho2, numk, nele, GzxFx, rho2);
  array_multiply_kl (GzyFyConvRho2, numk, nele, GzyFy, rho2);
  array_multiply_kl (GzzFzConvRho2, numk, nele, GzzFz, rho2);
  
  printf ("# here 3\n");
  fflush (stdout);

  // array_multiply (error2, nele, k2ma, rho2);
  
  printf("! %f\n", FConvRho1z[0][0][0]);
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    fftw_execute (p_backward_FConvRho1x[myk]);
    fftw_execute (p_backward_FConvRho1y[myk]);
    fftw_execute (p_backward_FConvRho1z[myk]);
  }
  printf("! %f\n", FConvRho1z[0][0][0]);
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (int ll = -kcut; ll <= kcut; ++ll){
      if (ll == 0) continue;
      int myl = ll + kcut;
      int klindex = kl2to1(myk, myl);
      fftw_execute (p_backward_FxFxConvRho2[klindex]);
      fftw_execute (p_backward_FyFyConvRho2[klindex]);
      fftw_execute (p_backward_FzFzConvRho2[klindex]);
      fftw_execute (p_backward_GxxFxConvRho2[klindex]);
      fftw_execute (p_backward_GxyFyConvRho2[klindex]);
      fftw_execute (p_backward_GxzFzConvRho2[klindex]);
      fftw_execute (p_backward_GyxFxConvRho2[klindex]);
      fftw_execute (p_backward_GyyFyConvRho2[klindex]);
      fftw_execute (p_backward_GyzFzConvRho2[klindex]);
      fftw_execute (p_backward_GzxFxConvRho2[klindex]);
      fftw_execute (p_backward_GzyFyConvRho2[klindex]);
      fftw_execute (p_backward_GzzFzConvRho2[klindex]);
      for (int index = 0; index < nele; ++index){
	FxFxConvRho2[klindex][index][0] /= volume;
	FyFyConvRho2[klindex][index][0] /= volume;
	FzFzConvRho2[klindex][index][0] /= volume;
	GxxFxConvRho2[klindex][index][0] /= volume;
	GxyFyConvRho2[klindex][index][0] /= volume;
	GxzFzConvRho2[klindex][index][0] /= volume;
	GyxFxConvRho2[klindex][index][0] /= volume;
	GyyFyConvRho2[klindex][index][0] /= volume;
	GyzFzConvRho2[klindex][index][0] /= volume;
	GzxFxConvRho2[klindex][index][0] /= volume;
	GzyFyConvRho2[klindex][index][0] /= volume;
	GzzFzConvRho2[klindex][index][0] /= volume;
      }
    }
  }
    
  printf ("# here 4\n");
  fflush (stdout);

  // fftw_execute (p_backward_error2);

  printf("! %f\n", FConvRho1z[0][0][0]);
  for (int kk = -kcut; kk <= kcut; ++kk){
    if (kk == 0) continue;
    int myk = kk + kcut;
    for (int i = 0; i < nele; ++i){
      FConvRho1x[myk][i][0] /= volume;
      FConvRho1y[myk][i][0] /= volume;
      FConvRho1z[myk][i][0] /= volume;
      FConvRho1x[myk][i][1] /= volume;
      FConvRho1y[myk][i][1] /= volume;
      FConvRho1z[myk][i][1] /= volume;
      // error2[i][0] /= volume;
    }
  }
  printf("! %f\n", FConvRho1z[0][0][0]);
  
  IntVectorType ii;
  for (ii.x = 0; ii.x < refined_K.x; ++ ii.x){
    for (ii.y = 0; ii.y < refined_K.y; ++ ii.y){
      for (ii.z = 0; ii.z < refined_K.z; ++ ii.z){
	int indexii;
	indexii = index3to1 (ii, refined_K);
	interpolate (ii, spmeik.error1x, refined_error1x[indexii]);
	interpolate (ii, spmeik.error1y, refined_error1y[indexii]);
	interpolate (ii, spmeik.error1z, refined_error1z[indexii]);
	// refined_error1x[indexii][0] = 0.;
	// refined_error1y[indexii][0] = 0.;
	// refined_error1z[indexii][0] = 0.;
	// refined_error1x[indexii][0] = result.x;
	// refined_error1y[indexii][0] = result.y;
	// refined_error1z[indexii][0] = result.z;
	// refined_error1x[indexii][0] = spmeik.error1x[indexjj][0];
	// refined_error1y[indexii][0] = spmeik.error1y[indexjj][0];
	// refined_error1z[indexii][0] = spmeik.error1z[indexjj][0];
	// refined_error1x[indexii][1] = 0.;
	// refined_error1y[indexii][1] = 0.;
	// refined_error1z[indexii][1] = 0.;
	refined_termIx[indexii][0] = refined_termIx[indexii][1] = 0.;
	refined_termIy[indexii][0] = refined_termIy[indexii][1] = 0.;
	refined_termIz[indexii][0] = refined_termIz[indexii][1] = 0.;
	refined_termIIx[indexii][0] = refined_termIIx[indexii][1] = 0.;
	refined_termIIy[indexii][0] = refined_termIIy[indexii][1] = 0.;
	refined_termIIz[indexii][0] = refined_termIIz[indexii][1] = 0.;
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

	  int indexii;
	  indexii = index3to1 (ii, refined_K);
	  fftw_complex resultx, resulty, resultz;
	  
	  // int indexii;
	  // indexii = index3to1 (ii, refined_K);
	  // fftw_complex resultx, resulty, resultz;

	  interpolate (ii, FConvRho1x[myk], resultx);
	  interpolate (ii, FConvRho1y[myk], resulty);
	  interpolate (ii, FConvRho1z[myk], resultz);
	  printf("! %f\n", FConvRho1z[myk][0][0]);
	  printf("! %f\n", FConvRho1z[myk][0][1]);
	  printf("! %f\n", FConvRho1z[myk][19682][0]);
	  printf("! %f\n", FConvRho1z[myk][19682][1]);
	  printf("! %f\n", resultz[0]);
	  printf("! %f\n", resultz[1]);
	  
	  fftw_complex rx, ry, rz;
	  multiply (rx, fx, resultx);
	  multiply (ry, fy, resulty);
	  multiply (rz, fz, resultz);
	  refined_error1x[indexii][0] -= rx[0];
	  refined_error1y[indexii][0] -= ry[0];
	  refined_error1z[indexii][0] -= rz[0];
	  refined_error1x[indexii][1] -= rx[1];
	  refined_error1y[indexii][1] -= ry[1];
	  refined_error1z[indexii][1] -= rz[1];
	  
	  printf("! %f\n", resultz[0]);
	  printf("! %f\n", resultz[1]);
	  printf("! %f\n", fz[0]);
	  printf("! %f\n", fz[1]);
	  printf("! %f\n", refined_error1z[0][0]);

	  refined_termIx[indexii][0] -= rx[0];
	  refined_termIx[indexii][1] -= rx[1];
	  refined_termIy[indexii][0] -= ry[0];
	  refined_termIy[indexii][1] -= ry[1];
	  refined_termIz[indexii][0] -= rz[0];
	  refined_termIz[indexii][1] -= rz[1];
	  
	  multiply (rx, fx, selfx[myk]);
	  multiply (ry, fy, selfy[myk]);
	  multiply (rz, fz, selfz[myk]);
	  refined_error1x[indexii][0] -= rx[0];
	  refined_error1y[indexii][0] -= ry[0];
	  refined_error1z[indexii][0] -= rz[0];
	  refined_error1x[indexii][1] -= rx[1];
	  refined_error1y[indexii][1] -= ry[1];
	  refined_error1z[indexii][1] -= rz[1];

	  refined_termIIx[indexii][0] -= rx[0];
	  refined_termIIx[indexii][1] -= rx[1];
	  refined_termIIy[indexii][0] -= ry[0];
	  refined_termIIy[indexii][1] -= ry[1];
	  refined_termIIz[indexii][0] -= rz[0];
	  refined_termIIz[indexii][1] -= rz[1];
	}
      }
    }
  }

  // calculate error 2
  for (ii.x = 0; ii.x < refined_K.x; ++ ii.x){
    for (ii.y = 0; ii.y < refined_K.y; ++ ii.y){
      for (ii.z = 0; ii.z < refined_K.z; ++ ii.z){
	int indexii;
	indexii = index3to1 (ii, refined_K);
	interpolate (ii, spmeik.error2, refined_error2[indexii]);
	printf("%f\n", refined_error2[0][0]);
	refined_error2[indexii][0] += 
	    refined_error1x[indexii][0] * refined_error1x[indexii][0] +
	    refined_error1x[indexii][1] * refined_error1x[indexii][1] + 
	    refined_error1y[indexii][0] * refined_error1y[indexii][0] + 
	    refined_error1y[indexii][1] * refined_error1y[indexii][1] + 
	    refined_error1z[indexii][0] * refined_error1z[indexii][0] + 
	    refined_error1z[indexii][1] * refined_error1z[indexii][1] ;
	printf("%f\n", refined_error1x[0][0]);
	printf("%f\n", refined_error1y[0][0]);
	printf("%f\n", refined_error1z[0][0]);
	printf("%f\n", refined_error1x[0][1]);
	printf("%f\n", refined_error1y[0][1]);
	printf("%f\n", refined_error1z[0][1]);
	printf("%f\n", refined_error2[0][0]);
	
	// refined_error2[indexii][0] += 
	//     refined_termIIx[indexii][0] * refined_termIIx[indexii][0] +
	//     refined_termIIx[indexii][1] * refined_termIIx[indexii][1] + 
	//     refined_termIIy[indexii][0] * refined_termIIy[indexii][0] + 
	//     refined_termIIy[indexii][1] * refined_termIIy[indexii][1] + 
	//     refined_termIIz[indexii][0] * refined_termIIz[indexii][0] + 
	//     refined_termIIz[indexii][1] * refined_termIIz[indexii][1] ;
	// refined_error2[indexii][0] += 
	//     2. * spmeError1x[0] * refined_termIx[indexii][0] +
	//     2. * spmeError1x[1] * refined_termIx[indexii][1] + 
	//     2. * spmeError1y[0] * refined_termIy[indexii][0] + 
	//     2. * spmeError1y[1] * refined_termIy[indexii][1] + 
	//     2. * spmeError1z[0] * refined_termIz[indexii][0] + 
	//     2. * spmeError1z[1] * refined_termIz[indexii][1] ;
	// refined_error2[indexii][0] += 
	//     2. * spmeError1x[0] * refined_termIIx[indexii][0] +
	//     2. * spmeError1x[1] * refined_termIIx[indexii][1] + 
	//     2. * spmeError1y[0] * refined_termIIy[indexii][0] + 
	//     2. * spmeError1y[1] * refined_termIIy[indexii][1] + 
	//     2. * spmeError1z[0] * refined_termIIz[indexii][0] + 
	//     2. * spmeError1z[1] * refined_termIIz[indexii][1] ;
	// refined_error2[indexii][0] += 
	//     2. * refined_termIx[indexii][0] * refined_termIIx[indexii][0] +
	//     2. * refined_termIx[indexii][1] * refined_termIIx[indexii][1] + 
	//     2. * refined_termIy[indexii][0] * refined_termIIy[indexii][0] + 
	//     2. * refined_termIy[indexii][1] * refined_termIIy[indexii][1] + 
	//     2. * refined_termIz[indexii][0] * refined_termIIz[indexii][0] + 
	//     2. * refined_termIz[indexii][1] * refined_termIIz[indexii][1] ;
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

void ErrorEstimate_SPME_Ana::
print_error (const std::string & file) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < refined_K.x; ++i){
    IntVectorType idx;
    idx.x = i;
    idx.y = 0;
    idx.z = 0;
    unsigned index = index3to1 (idx, refined_K);
    fprintf (fp, "%f %e  \n",
	     (i + 0.5) * vecA.xx / refined_K.x,
	     sqrt( refined_error2[index][0])
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
  //   refined_error1xx[indexjj+0]
      
