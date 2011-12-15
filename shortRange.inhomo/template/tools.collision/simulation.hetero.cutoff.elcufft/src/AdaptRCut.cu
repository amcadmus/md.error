#include "AdaptRCut.h"

AdaptRCut::AdaptRCut ()
    : malloced (false)
{
}

AdaptRCut::~AdaptRCut ()
{
  freeAll();
}

static void
mallocArrayComplex (fftw_complex *** a,
		    const unsigned & nrc,
		    const unsigned & nele)
{
  *a = (fftw_complex **) malloc (sizeof(fftw_complex *) * nrc);
  for (unsigned i = 0; i < nrc; ++i){
    (*a)[i] = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  }
}

static void
mallocArrayComplex (cufftComplex *** a,
		    const unsigned & nrc,
		    const unsigned & nele)
{
  *a = (cufftComplex **) malloc (sizeof(cufftComplex *) * nrc);
  for (unsigned i = 0; i < nrc; ++i){
    cudaMalloc ((void**)&((*a)[i]), sizeof(cufftComplex) * nele);
  }
}

static void
freeArrayComplex (fftw_complex *** a,
		  const unsigned & nrc)
{
  for (unsigned i = 0; i < nrc; ++i){
    free ((*a)[i]);
  }
  free (*a);
  *a = NULL;
}

static void
freeArrayComplex (cufftComplex *** a,
		  const unsigned & nrc)
{
  for (unsigned i = 0; i < nrc; ++i){
    cudaFree ((*a)[i]);
  }
  free (*a);
  *a = NULL;
}

// static void
// mallocArrayDoule (double *** a,
// 		  const unsigned & nrc,
// 		  const unsigned & nele)
// {
//   *a = (double **) malloc (sizeof(double *) * nrc);
//   for (unsigned i = 0; i < nrc; ++i){
//     (*a)[i] = (double *) malloc (sizeof(double) * nele);
//   }
// }

// static void
// freeArrayDouble (double *** a,
// 		 const unsigned & nrc)
// {
//   for (unsigned i = 0; i < nrc; ++i){
//     free ((*a)[i]);
//   }
//   free (*a);
//   *a = NULL;
// }


void AdaptRCut::
freeAll ()
{
  if (malloced) {
    // free (rhok);
    // free (rhor);
    freeArrayComplex (&s1k,  nrc);
    freeArrayComplex (&s2kx, nrc);
    freeArrayComplex (&s2ky, nrc);
    freeArrayComplex (&s2kz, nrc);
    // freeArrayComplex (&error1k,  nrc);
    // freeArrayComplex (&error2kx, nrc);
    // freeArrayComplex (&error2ky, nrc);
    // freeArrayComplex (&error2kz, nrc);
    // freeArrayComplex (&error1r,  nrc);
    // freeArrayComplex (&error2rx, nrc);
    // freeArrayComplex (&error2ry, nrc);
    // freeArrayComplex (&error2rz, nrc);
    // freeArrayComplex (&error, nrc);

    free (copyBuff);
    cudaFree (d_rhor);
    cudaFree (d_rhok);
    freeArrayComplex (&d_s1k,  nrc);
    freeArrayComplex (&d_s2kx, nrc);
    freeArrayComplex (&d_s2ky, nrc);
    freeArrayComplex (&d_s2kz, nrc);
    freeArrayComplex (&d_error1k,  nrc);
    freeArrayComplex (&d_error2kx, nrc);
    freeArrayComplex (&d_error2ky, nrc);
    freeArrayComplex (&d_error2kz, nrc);
    freeArrayComplex (&d_error1r,  nrc);
    freeArrayComplex (&d_error2rx, nrc);
    freeArrayComplex (&d_error2ry, nrc);
    freeArrayComplex (&d_error2rz, nrc);
    // freeArrayComplex (&d_error, nrc);
    cudaFree (d_error);
    for (int i = 0; i < NSTREAM; ++i){
      cufftDestroy(plan[i]);
    }
    
    free (rcut);
    free (rcutIndex);
    free (result_error);
    cudaFree (d_rcutIndex);
    cudaFree (d_result_error);
    
    // fftw_destroy_plan (p_forward_rho);
    // for (int count = 0; count < nrc; ++count){
    //   fftw_destroy_plan (p_backward_error1[count]);
    //   fftw_destroy_plan (p_backward_error2x[count]);
    //   fftw_destroy_plan (p_backward_error2y[count]);
    //   fftw_destroy_plan (p_backward_error2z[count]);
    // }
    // free (p_backward_error1);
    // free (p_backward_error2x);
    // free (p_backward_error2y);
    // free (p_backward_error2z);

    for (int i = 0; i < NSTREAM; ++i){
      cudaStreamDestroy(stream[i]);
    }
    
    malloced = false;
  }
}


void AdaptRCut::
reinit (const double & rcmin,
	const double & rcmax,
	const double & rcstep,
	DensityProfile_PiecewiseConst & dp)
{
  freeAll ();

  for (int i = 0; i < NSTREAM; ++i){
    cudaStreamCreate(&stream[i]);
  }

  boxsize = dp.getBox();
  nx = dp.getNx();
  ny = dp.getNy();
  nz = dp.getNz();
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;
  printf ("# AdaptRCut nx, ny, nz are %d %d %d\n", nx, ny, nz);
  
  rcList.clear();  
  for (double tmprc = rcmin; tmprc <= rcmax; tmprc += rcstep){
    rcList.push_back(tmprc);
  }
  nrc = rcList.size();
  
  // rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  mallocArrayComplex (&s1k,  nrc, nele);
  mallocArrayComplex (&s2kx, nrc, nele);
  mallocArrayComplex (&s2ky, nrc, nele);
  mallocArrayComplex (&s2kz, nrc, nele);
  // mallocArrayComplex (&error1k,  nrc, nele);
  // mallocArrayComplex (&error2kx, nrc, nele);
  // mallocArrayComplex (&error2ky, nrc, nele);
  // mallocArrayComplex (&error2kz, nrc, nele);
  // mallocArrayComplex (&error1r,  nrc, nele);
  // mallocArrayComplex (&error2rx, nrc, nele);
  // mallocArrayComplex (&error2ry, nrc, nele);
  // mallocArrayComplex (&error2rz, nrc, nele);
  // mallocArrayComplex (&error, nrc, nele);

  copyBuff = (cufftComplex *) malloc (sizeof(cufftComplex) * nele);
  cudaMalloc ((void**)&d_rhor, sizeof(cufftComplex) * nele);
  cudaMalloc ((void**)&d_rhok, sizeof(cufftComplex) * nele);
  mallocArrayComplex (&d_s1k,  nrc, nele);
  mallocArrayComplex (&d_s2kx, nrc, nele);
  mallocArrayComplex (&d_s2ky, nrc, nele);
  mallocArrayComplex (&d_s2kz, nrc, nele);
  mallocArrayComplex (&d_error1k,  nrc, nele);
  mallocArrayComplex (&d_error2kx, nrc, nele);
  mallocArrayComplex (&d_error2ky, nrc, nele);
  mallocArrayComplex (&d_error2kz, nrc, nele);
  mallocArrayComplex (&d_error1r,  nrc, nele);
  mallocArrayComplex (&d_error2rx, nrc, nele);
  mallocArrayComplex (&d_error2ry, nrc, nele);
  mallocArrayComplex (&d_error2rz, nrc, nele);
  // mallocArrayComplex (&d_error, nrc, nele);
  cudaMalloc ((void**)&d_error, nele*nrc*sizeof(cufftComplex));
  for (int i = 0; i < NSTREAM; ++i){
    cufftPlan3d (&plan[i], nx, ny, nz, CUFFT_C2C);
    cufftSetStream (plan[i], stream[i]);
  }
  
  rcut = (double *) malloc (sizeof(double) * nele);
  rcutIndex = (unsigned *) malloc (sizeof(unsigned) * nele);
  result_error = (float *) malloc (sizeof(float) * nele);
  cudaMalloc ((void**)&d_rcutIndex, sizeof(unsigned) * nele);
  cudaMalloc ((void**)&d_result_error, sizeof(float) * nele);

  // p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_PATIENT);
  // p_backward_error1  = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2x = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2y = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2z = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // for (int count = 0; count < nrc; ++count){
  //   p_backward_error1[count] = fftw_plan_dft_3d (nx, ny, nz, error1k[count], error1r[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2x[count] = fftw_plan_dft_3d (nx, ny, nz, error2kx[count], error2rx[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2y[count] = fftw_plan_dft_3d (nx, ny, nz, error2ky[count], error2ry[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2z[count] = fftw_plan_dft_3d (nx, ny, nz, error2kz[count], error2rz[count], FFTW_BACKWARD,  FFTW_PATIENT);
  // }
  
  malloced = true;

  printf ("# reinit: start build k mat 1 ...");
  fflush (stdout);
  makeS1k (1., 1.);
  size_t sizec = sizeof(cufftComplex) * nele;
  for (int count = 0; count < nrc; ++ count){
    for (int i = 0; i < nele; ++i){
      copyBuff[i].x = s1k[count][i][0];
      copyBuff[i].y = s1k[count][i][1];
    }
    cudaMemcpy (d_s1k[count], copyBuff, sizec, cudaMemcpyHostToDevice);
  }
  printf (" done\n");  
  fflush (stdout);
  printf ("# reinit: start build k mat 2 ...");
  fflush (stdout);
  makeS2k (1., 1.);
  for (int count = 0; count < nrc; ++ count){
    for (int i = 0; i < nele; ++i){
      copyBuff[i].x = s2kx[count][i][0];
      copyBuff[i].y = s2kx[count][i][1];
    }
    cudaMemcpy (d_s2kx[count], copyBuff, sizec, cudaMemcpyHostToDevice);
    for (int i = 0; i < nele; ++i){
      copyBuff[i].x = s2ky[count][i][0];
      copyBuff[i].y = s2ky[count][i][1];
    }
    cudaMemcpy (d_s2ky[count], copyBuff, sizec, cudaMemcpyHostToDevice);
    for (int i = 0; i < nele; ++i){
      copyBuff[i].x = s2kz[count][i][0];
      copyBuff[i].y = s2kz[count][i][1];
    }
    cudaMemcpy (d_s2kz[count], copyBuff, sizec, cudaMemcpyHostToDevice);
  }
  printf (" done\n");  
  fflush (stdout);
}

// static void 
// array_multiply (fftw_complex ** a,
// 		const int nrc,
// 		const int nele,
// 		fftw_complex ** b,
// 		fftw_complex * c) 
// {
//   for (int count = 0; count < nrc; ++count){
//     for (int ii = 0; ii < nele; ++ii){
//       // double tmpr, tmpi;
//       a[count][ii][0] =
// 	  c[ii][0] * b[count][ii][0] - c[ii][1] * b[count][ii][1];
//       a[count][ii][1] =
// 	  c[ii][0] * b[count][ii][1] + c[ii][1] * b[count][ii][0];
//     }
//   }
// }

__global__ static void 
array_multiply (cufftComplex *a,
		const int nele,
		const cufftComplex *b,
		const cufftComplex *c)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
      a[ii].x =
	  c[ii].x * b[ii].x - c[ii].y * b[ii].y;
      a[ii].y =
	  c[ii].x * b[ii].y + c[ii].y * b[ii].x;
  }
}

__global__ static void 
formErrorEr1 (const int count,
	      const int nele,
	      const float volumei,
	      const cufftComplex * error1r,
	      cufftComplex * error)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;
  unsigned start = count * nele;

  if (ii < nele){
    error[start+ii].x = error1r[ii].x * volumei;
    error[start+ii].y = error1r[ii].y * volumei;
  }
}

__global__ static void
formErrorEr2 (const int count,
	      const int nele,
	      const float volumei,
	      const cufftComplex * error2rx,
	      const cufftComplex * error2ry,
	      const cufftComplex * error2rz,
	      cufftComplex * error)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;
  unsigned start = count * nele;

  if (ii < nele){
    float tx0 = error2rx[ii].x * volumei;
    float ty0 = error2ry[ii].x * volumei;
    float tz0 = error2rz[ii].x * volumei;
    float tx1 = error2rx[ii].y * volumei;
    float ty1 = error2ry[ii].y * volumei;
    float tz1 = error2rz[ii].y * volumei;    
    error[start+ii].x += tx0*tx0 + ty0*ty0 + tz0*tz0;
    error[start+ii].y += tx1*tx1 + ty1*ty1 + tz1*tz1;
  }
}
  
__global__ static void
formErrorSqrt (const int nele,
	       cufftComplex * error)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
    error[ii].x = sqrtf(error[ii].x);
    error[ii].y = sqrtf(error[ii].y);
  }
}

#include "error.h"

void AdaptRCut::
calError (const DensityProfile_PiecewiseConst & dp)
{
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  double scale = volume / nele;
  size_t sizec = sizeof(cufftComplex) * nele;
  for (int i = 0; i < nele; ++i){
    copyBuff[i].x = dp.getProfile(i) * scale;
    copyBuff[i].y = 0.;
  }
  cudaMemcpy (d_rhor, copyBuff, sizec, cudaMemcpyHostToDevice);
  cufftExecC2C (plan[0], d_rhor, d_rhok, CUFFT_FORWARD);
  checkCUDAError ("AdaptRCut::calError: rhor -> rhok");

  unsigned blockSize = 128;
  unsigned nblock = unsigned(nele) / blockSize + 1;
  for (int count = 0; count < nrc; ++count){
    array_multiply <<<nblock, blockSize>>>
	(d_error1k [count], nele, d_s1k [count], d_rhok);
    array_multiply <<<nblock, blockSize>>>
	(d_error2kx[count], nele, d_s2kx[count], d_rhok);
    array_multiply <<<nblock, blockSize>>>
	(d_error2ky[count], nele, d_s2ky[count], d_rhok);
    array_multiply <<<nblock, blockSize>>>
	(d_error2kz[count], nele, d_s2kz[count], d_rhok);
  }
  checkCUDAError ("AdaptRCut::calError: ek = sk * rhok");

  for (int count = 0; count < nrc; ++count){
    cufftExecC2C (plan[(4*count+0)%NSTREAM], d_error1k [count], d_error1r [count], CUFFT_INVERSE);
    cufftExecC2C (plan[(4*count+1)%NSTREAM], d_error2kx[count], d_error2rx[count], CUFFT_INVERSE);
    cufftExecC2C (plan[(4*count+2)%NSTREAM], d_error2ky[count], d_error2ry[count], CUFFT_INVERSE);
    cufftExecC2C (plan[(4*count+3)%NSTREAM], d_error2kz[count], d_error2rz[count], CUFFT_INVERSE);
  }
  checkCUDAError ("AdaptRCut::calError: ek -> er");

  for (int count = 0; count < nrc; ++count){
    formErrorEr1 <<<nblock, blockSize, 0, stream[count%NSTREAM]>>>
	(count, nele, 1./volume, d_error1r[count], d_error);
  }
  checkCUDAError ("AdaptRCut::calError: cal error e1");
  for (int count = 0; count < nrc; ++count){
    formErrorEr2 <<<nblock, blockSize, 0, stream[count%NSTREAM]>>>
	(count, nele, 1./volume,
	 d_error2rx[count],
	 d_error2ry[count],
	 d_error2rz[count],
	 d_error);
  }
  checkCUDAError ("AdaptRCut::calError: cal error e2");
  dim3 tmpNBlock;
  tmpNBlock.x = nrc;
  tmpNBlock.y = nblock;
  // printf ("tmpNBlock is %d\n", tmpNBlock);
  cudaThreadSynchronize();
  formErrorSqrt <<<tmpNBlock, blockSize>>>
      (nele*nrc, d_error);
  checkCUDAError ("AdaptRCut::calError: cal error sqrt");
  
  // for (int count = 0; count < nrc; ++count){
  //   for (int i = 0; i < nele; ++i){
  //     error1r[count][i][0] /= volume;
  //     error1r[count][i][1] /= volume;
  //     error2rx[count][i][0] /= volume;
  //     error2ry[count][i][0] /= volume;
  //     error2rz[count][i][0] /= volume;
  //     error2rx[count][i][1] /= volume;
  //     error2ry[count][i][1] /= volume;
  //     error2rz[count][i][1] /= volume;
  //     error[count][i][0] =
  // 	  sqrt (
  // 	      error1r[count][i][0] +
  // 	      error2rx[count][i][0] * error2rx[count][i][0] +
  // 	      error2ry[count][i][0] * error2ry[count][i][0] +
  // 	      error2rz[count][i][0] * error2rz[count][i][0] );
  //     error[count][i][1] =
  // 	  sqrt (
  // 	      error1r[count][i][1] +
  // 	      error2rx[count][i][1] * error2rx[count][i][1] +
  // 	      error2ry[count][i][1] * error2ry[count][i][1] +
  // 	      error2rz[count][i][1] * error2rz[count][i][1] );
  //   }
  // }
}

// double AdaptRCut::
// maxError () const 
// {
//   double max = 0;
//   for (int i = 0; i < nele; ++i){
//     if (error[nrc-1][i][0] > max) max = error[nrc-1][i][0];
//   }
//   return max;
// }

double AdaptRCut::
integral_an13_numerical (const double & k,
			 const double & r1,
			 const double & r2,
			 const double & prec)
{
  f_inte13.k = k;
  return inte13.cal_int (Integral1DInfo::Gauss4,
			 f_inte13,
			 r1, r2,
			 prec);
}


double AdaptRCut::
integral_an5_numerical (const double & k,
			const double & r1,
			const double & r2,
			const double & prec)
{
  f_inte5.k = k;
  return inte5.cal_int (Integral1DInfo::Gauss4,
			f_inte5,
			r1, r2,
			prec);
}


void AdaptRCut::
makeS1k (const double & epsilon,
	 const double & sigma) 
{
  // int count = 0;
  double sigma12 = pow(sigma, 12.);
  double pref = 16. * 36. * sigma12 * epsilon;

  for (int count = nrc-1; count >= 0; --count){
    for (int ix = -nx/2; ix < nx - nx/2; ++ix){
      int locx;
      if (ix < 0) locx = ix + nx;
      else locx = ix;
      for (int iy = -ny/2; iy < ny - ny/2; ++iy){
	int locy;
	if (iy < 0) locy = iy + ny;
	else locy = iy;
	for (int iz = -nz/2; iz < nz - nz/2; ++iz){
	  int locz;
	  if (iz < 0) locz = iz + nz;
	  else locz = iz;
	  double k = sqrt(ix*ix/(boxsize[0]*boxsize[0]) +
			  iy*iy/(boxsize[1]*boxsize[1]) +
			  iz*iz/(boxsize[2]*boxsize[2]) );
	  int posi = index3to1 (locx, locy, locz);
	  if (k != 0){
	    double an;
	    if (count == nrc - 1){
	      an = integral_an13_numerical (k, rcList.back(), 30, 1e-12);
	      s1k[count][posi][0] = 2. / k * pref * an;
	      s1k[count][posi][1] = 0.;
	    }
	    else {
	      an = integral_an13_numerical (k, rcList[count], rcList[count+1], 1e-12);
	      s1k[count][posi][0] = s1k[count+1][posi][0] + 2. / k * pref * an;
	      s1k[count][posi][1] = 0.;	    
	    }
	    if (ix == -nx/2 ||
	    	iy == -ny/2 ||
	    	iz == -nz/2){
	      s1k[count][posi][0] = 0.;
	      s1k[count][posi][1] = 0.;
	    }
	  }
	  else {
	    s1k[count][posi][0] =
		4. * M_PI * pref / 11. * gsl_pow_int(rcList[count], -11);
	    s1k[count][posi][1] = 0.;
	  }
	}
      }
    }
  }
}

void AdaptRCut::
makeS2k (const double & epsilon,
	 const double & sigma) 
{
  double sigma6 = gsl_pow_6 (sigma);
  double pre = 8. * epsilon * sigma6;
  
  double kx, ky, kz;
  for (int count = nrc-1; count >= 0; --count){
    for (int ix = -nx/2; ix < nx - nx/2; ++ix){
      kx = ix / boxsize[0];
      int locx;
      if (ix < 0) locx = ix + nx;
      else locx = ix;
      for (int iy = -ny/2; iy < ny - ny/2; ++iy){
	ky = iy / boxsize[1];
	int locy;
	if (iy < 0) locy = iy + ny;
	else locy = iy;
	for (int iz = -nz/2; iz < nz - nz/2; ++iz){
	  kz = iz / boxsize[2];
	  int locz;
	  if (iz < 0) locz = iz + nz;
	  else locz = iz;
	  double k = sqrt(kx*kx + ky*ky + kz*kz);
	  int posi = index3to1 (locx, locy, locz);
	  if (k != 0){
	    double an;
	    if (count == nrc-1){
	      an = integral_an5_numerical (k, rcList[count], 30, 1e-16) * pre / k;
	      s2kx[count][posi][1] = 2 * M_PI * kx * an;
	      s2ky[count][posi][1] = 2 * M_PI * ky * an;
	      s2kz[count][posi][1] = 2 * M_PI * kz * an;
	    }
	    else {
	      an = integral_an5_numerical (k, rcList[count], rcList[count+1], 1e-16) * pre / k;
	      s2kx[count][posi][1] = s2kx[count+1][posi][1] + 2 * M_PI * kx * an;
	      s2ky[count][posi][1] = s2ky[count+1][posi][1] + 2 * M_PI * ky * an;
	      s2kz[count][posi][1] = s2kz[count+1][posi][1] + 2 * M_PI * kz * an;
	    }
	  }
	}
      }
    }
  }


  for (int count = nrc-1; count >= 0; --count){
    for (int ix = -nx/2; ix < nx - nx/2; ++ix){
      kx = ix / boxsize[0];
      int locx;
      if (ix < 0) locx = ix + nx;
      else locx = ix;
      for (int iy = -ny/2; iy < ny - ny/2; ++iy){
	ky = iy / boxsize[1];
	int locy;
	if (iy < 0) locy = iy + ny;
	else locy = iy;
	for (int iz = -nz/2; iz < nz - nz/2; ++iz){
	  kz = iz / boxsize[2];
	  int locz;
	  if (iz < 0) locz = iz + nz;
	  else locz = iz;
	  double k = sqrt(kx*kx + ky*ky + kz*kz);
	  int posi = index3to1 (locx, locy, locz);
	  double rc = rcList[count];
	  double rc6 = gsl_pow_6 (rc);
	  if (k != 0){
	    double erPiKR = 2. * M_PI * k * rc;
	    double size = - 2. * M_PI * rc * rc * (4. * epsilon * sigma6 / rc6) *
		(2. * cos(erPiKR) / erPiKR - 2. * sin(erPiKR) / (erPiKR*erPiKR));
	    size /= k;
	    s2kx[count][posi][0] = 0.;
	    s2kx[count][posi][1] += size * kx;
	    s2ky[count][posi][0] = 0.;
	    s2ky[count][posi][1] += size * ky;
	    s2kz[count][posi][0] = 0.;
	    s2kz[count][posi][1] += size * kz;
	    if (ix == -nx/2 ||
	    	iy == -ny/2 ||
	    	iz == -nz/2){
	      s2kx[count][posi][0] = 0.;
	      s2kx[count][posi][1] = 0.;
	      s2ky[count][posi][0] = 0.;
	      s2ky[count][posi][1] = 0.;
	      s2kz[count][posi][0] = 0.;
	      s2kz[count][posi][1] = 0.;
	    }
	  }
	  else {
	    s2kx[count][posi][0] = 0.;
	    s2kx[count][posi][1] = 0.;
	    s2ky[count][posi][0] = 0.;
	    s2ky[count][posi][1] = 0.;
	    s2kz[count][posi][0] = 0.;
	    s2kz[count][posi][1] = 0.;
	  }
	}
      }
    }
  }

}

// void AdaptRCut::
// calRCutOnePoint (const double & prec,
// 		 const unsigned & idx)
// {
//   unsigned posia = 0;
//   unsigned posib = nrc-1;
//   double rca = rcList[posia];
//   double rcb = rcList[posib];
//   double errora = error[posia][idx][0];
//   double errorb = error[posib][idx][0];
//   double diffa = errora - prec;
//   double diffb = errorb - prec;
//   if (diffa <= 0){
//     rcut[idx] = rca;
//     rcutIndex[idx] = posia;
//     result_error[idx] = errora;
//   }
//   else if (diffb >= 0){
//     rcut[idx] = rcb;
//     rcutIndex[idx] = posib;
//     result_error[idx] = errorb;
//   }
//   else {
//     while (posib - posia > 1) {
//       unsigned posic = (posib + posia) / 2;
//       double rcc = rcList[posic];
//       double errorc = error[posic][idx][0];
//       if (errorc > prec){
// 	posia = posic;
// 	rca = rcc;
// 	errora = errorc;
//       }
//       else {
// 	posib = posic;
// 	rcb = rcc;
// 	errorb = errorc;
//       }
//     }
//     rcut[idx] = rcb;
//     rcutIndex[idx] = posib;
//     result_error[idx] = errorb;
//   }
// }   

__global__ void static
calRCut_d (const unsigned nrc,
	   const unsigned nele,
	   const cufftComplex *error,
	   const float prec,
	   unsigned * rcutIndex,
	   float * result_error)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
    unsigned posia = 0;
    unsigned posib = nrc-1;
    float errora = error[posia * nele + ii].x;
    float errorb = error[posib * nele + ii].x;
    float diffa = errora - prec;
    float diffb = errorb - prec;
    if (diffa <= 0){
      rcutIndex[ii] = posia;
      result_error[ii] = errora;
    }
    else if (diffb >= 0){
      rcutIndex[ii] = posib;
      result_error[ii] = errorb;
    }
    else {
      while (posib - posia > 1) {
	unsigned posic = (posib + posia) / 2;
	double errorc = error[posic * nele + ii].x;
	if (errorc > prec){
	  posia = posic;
	  errora = errorc;
	}
	else {
	  posib = posic;
	  errorb = errorc;
	}
      }
      rcutIndex[ii] = posib;
      result_error[ii] = errorb;
    }
  }   
}		       


void AdaptRCut::
calRCut  (const double & prec)
{
  // for (int i = 0; i < nele; ++i){
  //   calRCutOnePoint (prec, i);
  // }
  unsigned blockSize = 128;
  unsigned nblock = unsigned(nele) / blockSize + 1;
  calRCut_d <<<nblock, blockSize>>>
      (nrc, nele, d_error, prec, d_rcutIndex, d_result_error);
  checkCUDAError ("AdaptRCut::calRCut, decide cutoff");
  cudaMemcpy (rcutIndex, d_rcutIndex, nele * sizeof(unsigned), cudaMemcpyDeviceToHost);
  cudaMemcpy (result_error, d_result_error, nele * sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError ("AdaptRCut::calRCut, copy back to host");
  for (int i = 0; i < nele; ++i){
    rcut[i] = rcList[rcutIndex[i]];
  }
}

__device__ void
index1to3_d (int input,
	     int nx, int ny, int nz,
	     int* ix, int* iy, int* iz)
{
  int tmp = input;
  *iz = tmp % (nz);
  tmp = (tmp - *iz) / nz;
  *iy = tmp % (ny);
  *ix = (tmp - *iy) / ny;
}

__device__ int
index3to1_d (int nx, int ny, int nz,
	     int ix, int iy, int iz) 
{
  return iz + nz * (iy + ny * ix);
}


// __global__ static void
// refineRCut_1 (const int nele,
// 	      unsigned * rcutIndex_bk)
// {
//   unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
//   unsigned tid = threadIdx.x;
//   unsigned ii = tid + bid * blockDim.x;

//   if (ii < nele){
//     rcutIndex_bk[ii] = 0;
//   }
// }


__global__ static void
refineRCut_d0 (const int nx,
	       const int ny,
	       const int nz,
	       const int nele,
	       const unsigned * rcutIndex,
	       unsigned * rcutIndex_bk)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;
  
  if (ii < nele){
    rcutIndex_bk[ii] = 0;
    int ix, iy, iz;
    index1to3_d (ii, nx, ny, nz, &ix, &iy, &iz);
    for (int dx = -1; dx <= 1; ++dx){
      int jx = ix + dx;
      if      (jx <  0 ) jx += nx;
      else if (jx >= nx) jx -= nx;
      for (int dy = -1; dy <= 1; ++dy){
	int jy = iy + dy;
	if      (jy <  0 ) jy += ny;
	else if (jy >= ny) jy -= ny;
	for (int dz = -1; dz <= 1; ++dz){
	  int jz = iz + dz;
	  if      (jz <  0 ) jz += nz;
	  else if (jz >= nz) jz -= nz;
	  if (rcutIndex[index3to1_d(nx, ny, nz, jx, jy, jz)] > rcutIndex_bk[ii]){
	    rcutIndex_bk[ii] = rcutIndex[index3to1_d(nx, ny, nz, jx, jy, jz)];
	  }
	}
      }
    }
  }
}

__global__ static void
refineRCut_d1 (const int nx,
	       const int ny,
	       const int nz,
	       const int nele,
	       unsigned * rcutIndex,
	       const unsigned * rcutIndex_bk)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
    rcutIndex[ii] = rcutIndex_bk[ii];
  }
}


void AdaptRCut::
refineRCut ()
{
  unsigned * rcutIndex_bk;
  cudaMalloc ((void**)&rcutIndex_bk, sizeof (unsigned) * nele);
  checkCUDAError ("AdaptRCut::refineRCut, malloc rcutIndex_bk");
  unsigned blockSize = 128;
  unsigned nblock = unsigned(nele) / blockSize + 1;
  refineRCut_d0
      <<<nblock, blockSize>>>
      (nx, ny, nz, nele, d_rcutIndex, rcutIndex_bk);
  refineRCut_d1
      <<<nblock, blockSize>>>
      (nx, ny, nz, nele, d_rcutIndex, rcutIndex_bk);
  checkCUDAError ("AdaptRCut::refineRCut, decide cutoff");
  cudaFree (rcutIndex_bk);
  checkCUDAError ("AdaptRCut::refineRCut, free rcutIndex_bk");
  cudaMemcpy (rcutIndex, d_rcutIndex, nele * sizeof(unsigned), cudaMemcpyDeviceToHost);
  checkCUDAError ("AdaptRCut::refineRCut, copy back to host");
  for (int i = 0; i < nele; ++i){
    rcut[i] = rcList[rcutIndex[i]];
  }
}

// void AdaptRCut::
// refineRCut () 
// {
//   int * rcutIndex_bk = (int *) malloc (sizeof(int) * nele);
//   // float * result_error_bk = (float*) malloc (sizeof(float) * nele);
//   for (int i = 0; i < nele; ++i){
//     rcutIndex_bk[i] = 0;
//   }

//   for (int i = 0; i < nele; ++i){
//     int ix, iy, iz;
//     index1to3 (i, ix, iy, iz);
//     for (int dx = -1; dx <= 1; ++dx){
//       int jx = ix + dx;
//       if      (jx <  0 ) jx += nx;
//       else if (jx >= nx) jx -= nx;
//       for (int dy = -1; dy <= 1; ++dy){
// 	int jy = iy + dy;
// 	if      (jy <  0 ) jy += ny;
// 	else if (jy >= ny) jy -= ny;
// 	for (int dz = -1; dz <= 1; ++dz){
// 	  int jz = iz + dz;
// 	  if      (jz <  0 ) jz += nz;
// 	  else if (jz >= nz) jz -= nz;
// 	  if (rcutIndex[index3to1(jx, jy, jz)] > rcutIndex_bk[i]){
// 	    rcutIndex_bk[i] = rcutIndex[index3to1(jx, jy, jz)];
// 	  }
// 	}
//       }
//     }
//   }

//   for (int i = 0; i < nele; ++i){
//     rcutIndex[i] = rcutIndex_bk[i];
//     rcut[i] = rcList[rcutIndex[i]];
//     // result_error[i] = error[rcutIndex[i]][i][0];
//   }

//   free (rcutIndex_bk);
// }

void AdaptRCut::    
print_x (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    // double sum = 0.;
    // for (int j = 0; j < ny; ++j){
    //   for (int k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    // fprintf (fp, "%f %e %e\n",
    // 	     (i + 0.5) * hx,
    // 	     error[4][index3to1(i,0,0)][0],
    // 	     error[4][index3to1(i,0,0)][1]
    // 	);
    fprintf (fp, "%f %e %f\n",
	     (i + 0.5) * hx,
	     result_error[index3to1(i,0,0)],
	     rcut[index3to1(i, 0, 0)]);
  }
  fclose (fp);
}


void AdaptRCut::    
print_error_avg (const DensityProfile_PiecewiseConst & dp,
		 const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    double sum = 0.;
    double sumprofile = 0.;
    for (int j = 0; j < ny; ++j){
      for (int k = 0; k < nz; ++k){
    	sum += result_error[index3to1(i,j,k)] *
	    result_error[index3to1(i,j,k)] *
	    dp.getProfile(i,j,k);
	sumprofile += dp.getProfile (i,j,k);
      }
    }
    sum /= sumprofile;
    sum = sqrt(sum);
    // fprintf (fp, "%f %e %e\n",
    // 	     (i + 0.5) * hx,
    // 	     error[4][index3to1(i,0,0)][0],
    // 	     error[4][index3to1(i,0,0)][1]
    // 	);
    fprintf (fp, "%f %e\n",
	     (i + 0.5) * hx,
	     sum);
  }
  fclose (fp);
}


  
void AdaptRCut::    
print_rc (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    // double sum = 0.;
    // for (int j = 0; j < ny; ++j){
    //   for (int k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    fprintf (fp, "%f %e %e\n",
	     (i + 0.5) * hx,
	     rcut[index3to1(i,ny/2,nz/2)], result_error[index3to1(i,ny/2,nz/2)]
	);
  }
  fclose (fp);
}


// void AdaptRCut::    
// print_rc_xy (const std::string & file) const 
// {
//   FILE * fp = fopen (file.c_str(), "w");
//   if (fp == NULL){
//     std::cerr << "cannot open file " << file << std::endl;
//     exit(1);
//   }

//   for (int i = 0; i < nx; ++i){
//     for (int k = 0; k < nz; ++k){
//       fprintf (fp, "%f %f %e %e\n",
// 	       (i + 0.5) * hx, (k + 0.5) * hz,
// 	       rcut[index3to1(i,ny/2,k)], result_error[index3to1(i,ny/2,k)]
// 	  );
//     }
//     printf (fp, "\n");
//   }
//   fclose (fp);
// }


void AdaptRCut::    
print_rc_avg (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    double sum = 0.;
    for (int j = 0; j < ny; ++j){
      for (int k = 0; k < nz; ++k){
    	sum += rcut[index3to1(i, j, k)];
      }
    }
    sum /= ny * nz;
    fprintf (fp, "%f %e\n",
	     (i + 0.5) * hx, sum
	);
  }
  fclose (fp);
}


void AdaptRCut::
load_rc (const std::string & file)
{
  FILE * fp = fopen (file.c_str(), "r");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file.c_str());
    exit (1);
  }
  fscanf (fp, "%d", &nrc);
  rcList.resize (nrc);
  for (int i = 0; i < nrc; ++i){
    fscanf (fp, "%lf", &rcList[i]);
  }
  fscanf (fp, "%lf %lf %lf", &boxsize[0], &boxsize[1], &boxsize[2]);
  fscanf (fp, "%d %d %d", &nx, &ny, &nz);
  nele = nx * ny * nz;

  freeAll();

  // rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  mallocArrayComplex (&s1k,  nrc, nele);
  mallocArrayComplex (&s2kx, nrc, nele);
  mallocArrayComplex (&s2ky, nrc, nele);
  mallocArrayComplex (&s2kz, nrc, nele);
  // mallocArrayComplex (&error1k,  nrc, nele);
  // mallocArrayComplex (&error2kx, nrc, nele);
  // mallocArrayComplex (&error2ky, nrc, nele);
  // mallocArrayComplex (&error2kz, nrc, nele);
  // mallocArrayComplex (&error1r,  nrc, nele);
  // mallocArrayComplex (&error2rx, nrc, nele);
  // mallocArrayComplex (&error2ry, nrc, nele);
  // mallocArrayComplex (&error2rz, nrc, nele);
  // mallocArrayComplex (&error, nrc, nele);
  rcut = (double *) malloc (sizeof(double) * nele);
  rcutIndex = (unsigned *) malloc (sizeof(unsigned) * nele);
  result_error = (float *) malloc (sizeof(float) * nele);

  for (int i = 0; i < nele; ++i){
    fscanf (fp, "%d ", &rcutIndex[i]);
    rcut[i] = rcList[rcutIndex[i]];
  }
  
  // p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_PATIENT);
  // p_backward_error1  = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2x = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2y = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // p_backward_error2z = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  // for (int count = 0; count < nrc; ++count){
  //   p_backward_error1[count] = fftw_plan_dft_3d (nx, ny, nz, error1k[count], error1r[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2x[count] = fftw_plan_dft_3d (nx, ny, nz, error2kx[count], error2rx[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2y[count] = fftw_plan_dft_3d (nx, ny, nz, error2ky[count], error2ry[count], FFTW_BACKWARD,  FFTW_PATIENT);
  //   p_backward_error2z[count] = fftw_plan_dft_3d (nx, ny, nz, error2kz[count], error2rz[count], FFTW_BACKWARD,  FFTW_PATIENT);
  // }
  
  malloced = true;

  printf ("# reinit: start build k mat 1 ...");
  fflush (stdout);
  makeS1k (1., 1.);
  printf (" done\n");  
  fflush (stdout);
  printf ("# reinit: start build k mat 2 ...");
  fflush (stdout);
  makeS2k (1., 1.);
  printf (" done\n");  
  fflush (stdout);
}
  

void AdaptRCut::
save_rc (const std::string & file) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file.c_str());
    exit (1);
  }

  fprintf (fp, "%d ", nrc);
  for (int i = 0; i < nrc; ++i){
    fprintf (fp, "%f ", rcList[i]);
  }
  fprintf (fp, "\n%f %f %f\n", boxsize[0], boxsize[1], boxsize[2]);
  fprintf (fp, "%d %d %d\n", nx, ny, nz);
  for (int i = 0; i < nele; ++i){
    fprintf (fp, "%d ", rcutIndex[i]);
  }
  fprintf(fp, "\n");
  
  // fwrite (&nrc, sizeof(int), 1, fp);
  // fwrite (rcList, sizeof(double), nrc, fp);
  // fwrite (boxsize, sizeof(double), 3, fp);
  // fwrite (&nx, sizeof(int), 1, fp);
  // fwrite (&ny, sizeof(int), 1, fp);
  // fwrite (&nz, sizeof(int), 1, fp);
  // fwrite (rcutIndex, sizeof(int), nele, fp);
  
  fclose (fp);
}

void AdaptRCut::
write_rc (const std::string & file) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  fprintf (fp, "%d %d %d\n", nx, ny, nz);
  for (int i = 0; i < nele; ++i){
    fprintf (fp, "%f ", rcut[i]);
  }

  fclose (fp);
}


