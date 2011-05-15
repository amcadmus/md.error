#include "ErrorEstimator.h"

ErrorEstimatorFFT_Corr::
ErrorEstimatorFFT_Corr (const ForceKernel & fk_,
			const DensityProfile_Corr_PiecewiseConst & dp)
    : fk (fk_)
{
  boxsize = dp.getBox();
  nx = dp.getNx();
  ny = dp.getNy();
  nz = dp.getNz();
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  corrBond = dp.getCorrBond();
  corrDim  = dp.getCorrDim();
  numCorr  = corrDim * corrDim * corrDim;
  corrCut = 2.;

  double prec = 1e-6;
  bond_f_i3 = 50;
  bond_ff_i1 = 50;
  bond_ff_i5 = 50;
  prec_f_i3 = prec;
  prec_ff_i1 = prec;
  prec_ff_i5 = prec;
  
  nele = nx * ny * nz;
  f2k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  b0k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  b1k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  b2k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a00k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a01k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a02k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a11k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a12k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  a22k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);

  s1r  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s1k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2r  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2rx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2ry = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2rz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2kx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2ky = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s2kz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s3r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s3k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  
  rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corrr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corrk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr0r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr1r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr2r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr0k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr1k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr2k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr00r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr01r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr02r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr11r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr12r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr22r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr00k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr01k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr02k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr11k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr12k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  corr22k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a00r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a01r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a02r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a11r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a12r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  // a22r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  
  p_forward_rho  = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_rho = fftw_plan_dft_3d (nx, ny, nz, rhok, rhor, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_corr   = fftw_plan_dft_3d (nx, ny, nz, corrr,   corrk,   FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr0  = fftw_plan_dft_3d (nx, ny, nz, corr0r,  corr0k,  FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr1  = fftw_plan_dft_3d (nx, ny, nz, corr1r,  corr1k,  FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr2  = fftw_plan_dft_3d (nx, ny, nz, corr2r,  corr2k,  FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr00 = fftw_plan_dft_3d (nx, ny, nz, corr00r, corr00k, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr01 = fftw_plan_dft_3d (nx, ny, nz, corr01r, corr01k, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr02 = fftw_plan_dft_3d (nx, ny, nz, corr02r, corr02k, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr11 = fftw_plan_dft_3d (nx, ny, nz, corr11r, corr11k, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr12 = fftw_plan_dft_3d (nx, ny, nz, corr12r, corr12k, FFTW_FORWARD, FFTW_MEASURE);
  p_forward_corr22 = fftw_plan_dft_3d (nx, ny, nz, corr22r, corr22k, FFTW_FORWARD, FFTW_MEASURE);
  
  p_backward_s1  = fftw_plan_dft_3d (nx, ny, nz, s1k,  s1r,  FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s2x = fftw_plan_dft_3d (nx, ny, nz, s2kx, s2rx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s2y = fftw_plan_dft_3d (nx, ny, nz, s2ky, s2ry, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s2z = fftw_plan_dft_3d (nx, ny, nz, s2kz, s2rz, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s3  = fftw_plan_dft_3d (nx, ny, nz, s3k,  s3r,  FFTW_BACKWARD, FFTW_MEASURE);
  
  
  for (int i = 0; i < nele; ++i){
    rhor[i][1] = 0.;
  }
}

ErrorEstimatorFFT_Corr::
~ErrorEstimatorFFT_Corr () 
{
  fftw_destroy_plan (p_backward_s1);
  fftw_destroy_plan (p_backward_s2x);
  fftw_destroy_plan (p_backward_s2y);
  fftw_destroy_plan (p_backward_s2z);
  fftw_destroy_plan (p_backward_s3);
  fftw_destroy_plan (p_forward_rho);
  fftw_destroy_plan (p_backward_rho);

  fftw_free (f2k);
  fftw_free (b0k);
  fftw_free (b1k);
  fftw_free (b2k);
  fftw_free (a00k);
  fftw_free (a01k);
  fftw_free (a02k);
  fftw_free (a11k);
  fftw_free (a12k);
  fftw_free (a22k);

  fftw_free (s1r);
  fftw_free (s1k);
  fftw_free (s2r);
  fftw_free (s2rx);
  fftw_free (s2kx);
  fftw_free (s2ry);
  fftw_free (s2ky);
  fftw_free (s2rz);
  fftw_free (s2kz);
  fftw_free (s3r);
  fftw_free (s3k);

  fftw_free (rhor);
  fftw_free (rhok);
  fftw_free (corrr);
  fftw_free (corrk);
  fftw_free (corr0r);
  fftw_free (corr1r);
  fftw_free (corr2r);
  fftw_free (corr0k);
  fftw_free (corr1k);
  fftw_free (corr2k);
  fftw_free (corr00r);
  fftw_free (corr00k);
  fftw_free (corr01r);
  fftw_free (corr01k);
  fftw_free (corr02r);
  fftw_free (corr02k);
  fftw_free (corr11r);
  fftw_free (corr11k);
  fftw_free (corr12r);
  fftw_free (corr12k);
  fftw_free (corr22r);
  fftw_free (corr22k);
}


double ErrorEstimatorFFT_Corr::
integral_an5_numerical (const double & k,
			const double & rc)
{
  f_inte5.k = k;
  return inte5.cal_int (Integral1DInfo::Gauss4,
			f_inte5,
			rc, 50,
			1e-10);
}

double ErrorEstimatorFFT_Corr::
integral_an13_numerical (const double & k,
			 const double & rc)
{
  f_inte13.k = k;
  return inte13.cal_int (Integral1DInfo::Gauss4,
			 f_inte13,
			 rc, 30,
			 1e-16);
}

double ErrorEstimatorFFT_Corr::
integral_an15_numerical (const double & k,
			 const double & rc)
{
  f_inte15.k = k;
  return inte15.cal_int (Integral1DInfo::Gauss4,
			 f_inte15,
			 rc, 30,
			 1e-17);
}


double ErrorEstimatorFFT_Corr::
integral_f_i3_numerical (const double & k,
			 const double & rc)
{
  f_i3.k = k;
  return inte_f_i3.cal_int (Integral1DInfo::Gauss4,
			    f_i3,
			    rc, bond_f_i3,
			    prec_f_i3);
}

double ErrorEstimatorFFT_Corr::
integral_ff_i1_numerical (const double & k,
			  const double & rc)
{
  ff_i1.k = k;
  return inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
			    ff_i1,
			    rc, bond_ff_i1,
			    prec_ff_i1);
}

double ErrorEstimatorFFT_Corr::
integral_ff_i5_numerical (const double & k,
			  const double & rc)
{
  ff_i5.k = k;
  return inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
			    ff_i5,
			    rc, bond_ff_i5,
			    prec_ff_i5);
}


void ErrorEstimatorFFT_Corr::
formCorrMat (const DensityProfile_Corr_PiecewiseConst & dp)
{
  corrBond = dp.getCorrBond();
  corrDim = dp.getCorrDim();
  numCorr = corrDim * corrDim * corrDim;
  
  printf ("## init\n");
  for (int ii = 0; ii < nele; ++ii){
    corrr[ii][0] = corrr[ii][1] = 0.;
    corr0r[ii][0] = corr0r[ii][1] = 0.;
    corr1r[ii][0] = corr1r[ii][1] = 0.;
    corr2r[ii][0] = corr2r[ii][1] = 0.;
    corr00r[ii][0] = corr00r[ii][1] = 0.;
    corr01r[ii][0] = corr01r[ii][1] = 0.;
    corr02r[ii][0] = corr02r[ii][1] = 0.;
    corr11r[ii][0] = corr11r[ii][1] = 0.;
    corr12r[ii][0] = corr12r[ii][1] = 0.;
    corr22r[ii][0] = corr22r[ii][1] = 0.;    
  }

  double volumeEle = hx * hy * hz;
  // for (int ii = 0; ii < nele; ++ii){
  //   for (int di = -corrBond; di <= corrBond; ++di){
  //     double dx = di * hx;
  //     for (int dj = -corrBond; dj <= corrBond; ++dj){
  // 	double dy = dj * hy;
  // 	for (int dk = -corrBond; dk <= corrBond; ++dk){
  // 	  double dz = dk * hz;
  // 	  if (di == 0 && dj == 0 && dk == 0) continue;
  // 	  double dr = sqrt (dx*dx + dy*dy + dz*dz);
  // 	  // if (dr > corrCut) continue;
  // 	  double tmp = dp.getCorr(ii, -di, -dj, -dk) * volumeEle;
  // 	  corrr[ii][0] += tmp;
  // 	  corr0r[ii][0] += dx * tmp;
  // 	  corr1r[ii][0] += dy * tmp;
  // 	  corr2r[ii][0] += dz * tmp;
  // 	  corr00r[ii][0] += dx * dx * tmp;
  // 	  corr01r[ii][0] += dx * dy * tmp;
  // 	  corr02r[ii][0] += dx * dz * tmp;
  // 	  corr11r[ii][0] += dy * dy * tmp;
  // 	  corr12r[ii][0] += dy * dz * tmp;
  // 	  corr22r[ii][0] += dz * dz * tmp;
  // 	}
  //     }
  //   }
  // }
  std::vector<double > tmpcorr (corrDim, 0.);
  int tmpcorrCount = 0;
  printf ("## integral loop\n");
  for (int ix = 0; ix < nx; ++ix) {
    // printf ("%d\n", ix);
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
	int tx, ty, tz;
	int ii = index3to1 (ix, iy, iz);
	double rRho = dp.getMean (ii);
	for (int dx = -corrBond; dx <= corrBond; ++dx){
	  double drx = dx * hx;
	  tx = ix - dx;
	  if      (tx <   0) tx += nx;
	  else if (tx >= nx) tx -= nx;		
	  for (int dy = -corrBond; dy <= corrBond; ++dy){
	    double dry = dy * hy;
	    ty = iy - dy;
	    if      (ty <   0) ty += ny;
	    else if (ty >= ny) ty -= ny;		
	    for (int dz = -corrBond; dz <= corrBond; ++dz){
	      if (dx == 0 && dy == 0 && dz == 0) continue;
	      double drz = dz * hz;
	      double dr = sqrt (drx*drx + dry*dry + drz*drz);
	      // if (dr >= corrCut) continue;
	      tz = iz - dz;
	      if      (tz <   0) tz += nz;
	      else if (tz >= nz) tz -= nz;
	      double tRho = dp.getMean (tx, ty, tz);
	      double tmp = (dp.getCorr(ii, -dx, -dy, -dz) -
			    rRho * tRho) * volumeEle;
	      if (ix == nx/2 &&
		  dy == 0 && dz == 0){
		tmpcorr[dx+corrBond] += dp.getCorr(ii, -dx, -dy, -dz) - rRho * tRho;
		tmpcorrCount ++;
	      }
	      corrr[ii][0] += tmp;
	      corr0r[ii][0] += drx * tmp;
	      corr1r[ii][0] += dry * tmp;
	      corr2r[ii][0] += drz * tmp;
	      corr00r[ii][0] += drx * drx * tmp;
	      corr01r[ii][0] += drx * dry * tmp;
	      corr02r[ii][0] += drx * drz * tmp;
	      corr11r[ii][0] += dry * dry * tmp;
	      corr12r[ii][0] += dry * drz * tmp;
	      corr22r[ii][0] += drz * drz * tmp;
	    }
	  }
	}
      }
    }
  }
  for (int i = 0;i < corrDim; ++i){
    printf ("%d %f\n", i - corrBond, tmpcorr[i] / tmpcorrCount);
  }
  
  printf ("## ffts\n");
  fftw_execute (p_forward_corr);
  fftw_execute (p_forward_corr0);
  fftw_execute (p_forward_corr1);
  fftw_execute (p_forward_corr2);
  fftw_execute (p_forward_corr00);
  fftw_execute (p_forward_corr01);
  fftw_execute (p_forward_corr02);
  fftw_execute (p_forward_corr11);
  fftw_execute (p_forward_corr12);
  fftw_execute (p_forward_corr22);
  
  printf ("## rescale\n");
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  double scale = volume / nele;
  for (int i = 0; i < nele; ++i){
    // rhok[i][0] *= scale;
    // rhok[i][1] *= scale;
    corrk[i][0] *= scale;
    corrk[i][1] *= scale;
    corr0k[i][0] *= scale;
    corr0k[i][1] *= scale;
    corr1k[i][0] *= scale;
    corr1k[i][1] *= scale;
    corr2k[i][0] *= scale;
    corr2k[i][1] *= scale;
    corr00k[i][0] *= scale;
    corr00k[i][1] *= scale;
    corr01k[i][0] *= scale;
    corr01k[i][1] *= scale;
    corr02k[i][0] *= scale;
    corr02k[i][1] *= scale;
    corr11k[i][0] *= scale;
    corr11k[i][1] *= scale;
    corr12k[i][0] *= scale;
    corr12k[i][1] *= scale;
    corr22k[i][0] *= scale;
    corr22k[i][1] *= scale; 
  }  
}

static void 
array_multiply (fftw_complex * a,
		const int n,
		const fftw_complex * b) 
{
  for (int ii = 0; ii < n; ++ii){
    // double tmpr, tmpi;
    double tmpr = a[ii][0];
    double tmpi = a[ii][1];
    a[ii][0] = tmpr * b[ii][0] - tmpi * b[ii][1];
    a[ii][1] = tmpr * b[ii][1] + tmpi * b[ii][0];
  }
}


void ErrorEstimatorFFT_Corr::
estimate (const DensityProfile_Corr_PiecewiseConst & dp,
	  const int & corrBond,
	  const double & rc)
{
  for (int i = 0; i < nele; ++i){
    rhor[i][0] = dp.getMean(i);
  }
  fftw_execute (p_forward_rho);
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  double scale = volume / nele;
  for (int i = 0; i < nele; ++i){
    rhok[i][0] *= scale;
    rhok[i][1] *= scale;
  }  

  printf ("# start build corr mat\n");
  formCorrMat (dp);
  
  printf ("# start build k mat 1\n");
  makef2k (1., rc);
  printf ("# start build k mat 2\n");
  makeS2k (1., 1., rc);
  printf ("# start build k mat 3\n");
  makeS3k (1., 1., rc);
  
  printf ("# multiply on k space\n");
  for (int i = 0; i < nele; ++i){
    s1k[i][0] = f2k[i][0];
    s1k[i][1] = f2k[i][1];    
    s3k[i][0] = f2k[i][0];
    s3k[i][1] = f2k[i][1];    
  }
  
  array_multiply (s1k,  nele, rhok);
  array_multiply (s2kx, nele, rhok);
  array_multiply (s2ky, nele, rhok);
  array_multiply (s2kz, nele, rhok);
  
  array_multiply (s3k,  nele, corrk);
  array_multiply (b0k,  nele, corr0k);
  array_multiply (b1k,  nele, corr1k);
  array_multiply (b2k,  nele, corr2k);
  array_multiply (a00k, nele, corr00k);
  array_multiply (a01k, nele, corr01k);
  array_multiply (a02k, nele, corr02k);
  array_multiply (a11k, nele, corr11k);
  array_multiply (a12k, nele, corr12k);
  array_multiply (a22k, nele, corr22k);

  for (int ii = 0; ii < nele; ++ii){
    s3k[ii][0] += (b0k[ii][0] + b1k[ii][0] + b2k[ii][0] +
		   (a00k[ii][0] + a11k[ii][0] + a22k[ii][0]) * .5 +
		   (a01k[ii][0] + a02k[ii][0] + a12k[ii][0]));
    s3k[ii][1] += (b0k[ii][1] + b1k[ii][1] + b2k[ii][1] +
		   (a00k[ii][1] + a11k[ii][1] + a22k[ii][1]) * .5 +
		   (a01k[ii][1] + a02k[ii][1] + a12k[ii][1]));
  }
  // for (int ix = -nx/2; ix < nx - nx/2; ++ix){
  //   int locx;
  //   if (ix < 0) locx = ix + nx;
  //   else locx = ix;
  //   for (int iy = -ny/2; iy < ny - ny/2; ++iy){
  //     int locy;
  //     if (iy < 0) locy = iy + ny;
  //     else locy = iy;
  //     for (int iz = -nz/2; iz < nz - nz/2; ++iz){
  // 	int locz;
  // 	if (iz < 0) locz = iz + nz;
  // 	else locz = iz;
  // 	int posi = index3to1 (locx, locy, locz);
  // 	// printf ("%d %d %d  %e %e\n",
  // 	// 	ix, iy, iz,
  // 	// 	s3k[posi][0], s3k[posi][1]);
  //     }
  //   }
  // }

  printf ("# reverse fft\n");
  fftw_execute (p_backward_s1);
  fftw_execute (p_backward_s2x);
  fftw_execute (p_backward_s2y);
  fftw_execute (p_backward_s2z);
  fftw_execute (p_backward_s3);
  
  for (int i = 0; i < nele; ++i){
    s1r[i][0] /= volume;
    s1r[i][1] /= volume;
    s2rx[i][0] /= volume;
    s2ry[i][0] /= volume;
    s2rz[i][0] /= volume;
    s2rx[i][1] /= volume;
    s2ry[i][1] /= volume;
    s2rz[i][1] /= volume;
    s2r[i][0] = s2rx[i][0] * s2rx[i][0] + s2ry[i][0] * s2ry[i][0] + s2rz[i][0] * s2rz[i][0];
    s2r[i][1] = s2rx[i][1] * s2rx[i][1] + s2ry[i][1] * s2ry[i][1] + s2rz[i][1] * s2rz[i][1];
    s3r[i][0] /= -volume;
    s3r[i][1] /= -volume;
  }
}

void ErrorEstimatorFFT_Corr::
makef2k (const double & sigma,
	 const double & rc) 
{
  // int count = 0;
  double sigma12 = pow(sigma, 12.);
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
	  // if (2*M_PI*k*rc > 2.){
	  //   an = integral_an(13, k, rc);
	  // }
	  // else {
	  //   an = integral_an1(13, k, rc);
	  // }
	  an = integral_an13_numerical (k, rc);
	  f2k[posi][0] = 16. * 36. * 2. / k * sigma12 * an;
	  f2k[posi][1] = 0.;
	}
	else {
	  f2k[posi][0] = 16. * 36. * 4. * M_PI * sigma12 / 11. / pow(rc, 11);
	  f2k[posi][1] = 0.;
	}
	// count ++;
      }
    }
  }
}

void ErrorEstimatorFFT_Corr::
makeS2k (const double & epsilon,
	 const double & sigma,
	 const double & rc) 
{
  double sigma6 = pow(sigma, 6.);
  double rc6 = pow(rc, 6.);
  double pre = 8. * epsilon * sigma6;
  // double rc1 = 9.;
  // double rc16 = pow (rc1, 6.);
  
  double kx, ky, kz;
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
	// double k = sqrt(ix*ix/(boxsize[0]*boxsize[0]) +
	// 		iy*iy/(boxsize[1]*boxsize[1]) +
	// 		iz*iz/(boxsize[2]*boxsize[2]) );
	int posi = index3to1 (locx, locy, locz);
	if (k != 0){
	  double an = integral_an5_numerical (k, rc);
	  an *= pre / k;
	  double erPiKR = 2. * M_PI * k * rc;
	  double size = - 2. * M_PI * rc * rc * (4. * epsilon * sigma6 / rc6) *
	      (2. * cos(erPiKR) / erPiKR - 2. * sin(erPiKR) / (erPiKR*erPiKR));
	  // double erPiKR1 = 2. * M_PI * k * rc1;
	  // double size1 = 2. * M_PI * rc1 * rc1 * (4. * epsilon * sigma6 / rc16) *
	  //     (2. * cos(erPiKR1) / erPiKR1 - 2. * sin(erPiKR1) / (erPiKR1*erPiKR1));
	  size /= k;
	  // size1 /= k;
	  // size += size1;
	  s2kx[posi][0] = 0.;
	  s2kx[posi][1] = 2 * M_PI * kx * an + size * kx;
	  s2ky[posi][0] = 0.;
	  s2ky[posi][1] = 2 * M_PI * ky * an + size * kx;
	  s2kz[posi][0] = 0.;
	  s2kz[posi][1] = 2 * M_PI * kz * an + size * kz;
	  if (ix == -nx/2 ||
	      iy == -ny/2 ||
	      iz == -nz/2){
	    s2kx[posi][0] = 0.;
	    s2kx[posi][1] = 0.;
	    s2ky[posi][0] = 0.;
	    s2ky[posi][1] = 0.;
	    s2kz[posi][0] = 0.;
	    s2kz[posi][1] = 0.;
	  }
	}
	else {
	  s2kx[posi][0] = 0.;
	  s2kx[posi][1] = 0.;
	  s2ky[posi][0] = 0.;
	  s2ky[posi][1] = 0.;
	  s2kz[posi][0] = 0.;
	  s2kz[posi][1] = 0.;
	}
      }
    }
  }
}
  
void ErrorEstimatorFFT_Corr::
makeS3k (const double & epsilon,
	 const double & sigma,
	 const double & rc) 
{
  double sigma6 = pow(sigma, 6.);
  double pre = epsilon * sigma6;
  pre = pre*pre;

  double kx, ky, kz;
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
	int posi = index3to1 (locx, locy, locz);

	double k = sqrt(kx*kx + ky*ky + kz*kz);
	
	if (k != 0){
	  double phi_k = acos(kz/k);
	  double theta_k;
	  if (phi_k == 0.){
	    theta_k = 0.;
	  } else {
	    double sp = sin(phi_k);
	    double tmp = kx / (k * sp);
	    if (tmp > 1.) tmp = 1.;
	    if (tmp <-1.) tmp = -1.;
	    theta_k = acos(tmp);
	    if (ky < 0) theta_k = 2. * M_PI - theta_k;
	  }
	  double i3 = integral_f_i3_numerical (k, rc);
	  b0k[posi][0] = 0;
	  b0k[posi][1] = pre * 2. * M_PI / sqrt(k) * kx / k * (-i3);
	  b1k[posi][0] = 0;
	  b1k[posi][1] = pre * 2. * M_PI / sqrt(k) * ky / k * (-i3);
	  b2k[posi][0] = 0;
	  b2k[posi][1] = pre * 2. * M_PI / sqrt(k) * kz / k * (-i3);
	  double an15 = integral_an15_numerical (k, rc);
	  double ih = 1./k * (16.) * 24. * 24. * pre * an15;
	  //double ih = 0.;
	  double i1 = integral_ff_i1_numerical (k, rc);
	  double i5 = integral_ff_i5_numerical (k, rc);
	  double tmpscale = pre * M_PI / sqrt(k);
	  double cp = cos(phi_k);
	  double sp = sin(phi_k);
	  double ct = cos(theta_k);
	  double st = sin(theta_k);
	  a00k[posi][0] = tmpscale * (
	      2./3. * i1 +
	      (-1./3. + cp*cp - sp*sp*cos(theta_k*2)) * i5) + ih;
	  a00k[posi][1] = 0;
	  a11k[posi][0] = tmpscale * (
	      2./3. * i1 +
	      (-1./3. + cp*cp + sp*sp*cos(theta_k*2)) * i5) + ih;
	  a11k[posi][1] = 0;
	  a01k[posi][0] = - tmpscale * sp*sp*sin(theta_k*2) * i5;
	  a01k[posi][1] = 0;
	  a02k[posi][0] = - tmpscale * sin(phi_k*2)*ct * i5;
	  a02k[posi][1] = 0;
	  a12k[posi][0] = - tmpscale * sin(phi_k*2)*st * i5;
	  a12k[posi][1] = 0;
	  a22k[posi][0] = tmpscale * 4./3. *
	      (0.5*i1 + (0.5 - 1.5 * cp*cp) * i5) + ih;
	  a22k[posi][1] = 0;
	  if (ix == -nx/2 ||
	      iy == -ny/2 ||
	      iz == -nz/2){
	    b0k[posi][0] = 0.;
	    b1k[posi][0] = 0.;
	    b2k[posi][0] = 0.;
	    b0k[posi][1] = 0.;
	    b1k[posi][1] = 0.;
	    b2k[posi][1] = 0.;
	    a00k[posi][0] = 0;
	    a00k[posi][1] = 0;
	    a11k[posi][0] = 0;
	    a11k[posi][1] = 0;
	    a01k[posi][0] = 0;
	    a01k[posi][1] = 0;
	    a02k[posi][0] = 0;
	    a02k[posi][1] = 0;
	    a12k[posi][0] = 0;
	    a12k[posi][1] = 0;
	    a22k[posi][0] = 0;
	    a22k[posi][1] = 0;
	  }
	}
	else {
	  b0k[posi][0] = 0;
	  b0k[posi][1] = 0;
	  b1k[posi][0] = 0;
	  b1k[posi][1] = 0;
	  b2k[posi][0] = 0;
	  b2k[posi][1] = 0;
	  
	  a00k[posi][0] = 0;
	  a00k[posi][1] = 0;
	  a11k[posi][0] = 0;
	  a11k[posi][1] = 0;
	  a01k[posi][0] = 0;
	  a01k[posi][1] = 0;
	  a02k[posi][0] = 0;
	  a02k[posi][1] = 0;
	  a12k[posi][0] = 0;
	  a12k[posi][1] = 0;
	  a22k[posi][0] = 0;
	  a22k[posi][1] = 0;
	}
      }
    }
  }
}




void ErrorEstimatorFFT_Corr::    
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
    fprintf (fp, "%f %e %e %e %e %e %e\n",
	     (i + 0.5) * hx,
	     s1r[index3to1(i,0,0)][0],
	     s1r[index3to1(i,0,0)][1],
	     s2r[index3to1(i,0,0)][0],
	     s2r[index3to1(i,0,0)][1],
	     s3r[index3to1(i,0,0)][0],
	     s3r[index3to1(i,0,0)][1]
	);
  }
  fclose (fp);
}

void ErrorEstimatorFFT_Corr::    
print_xy (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    double vx = (i + 0.5) * hx;
    for (int j = 0; j < ny; ++j){
      double vy = (j + 0.5) * hy;
      // double sum = 0.;
      // for (int k = 0; k < nz; ++k){
      // 	sum += profile[index3to1(i, j, k)];
      // }
      fprintf (fp, "%f %f %e\n", vx, vy, s1r[index3to1(i,j,0)][0]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}




// double ErrorEstimatorFFT_Corr::
// integral_an1 (const int & n_,
// 	      const double & k,
// 	      const double & rc)
// {
//   int n = (n_ >> 1);
//   double erNjiecheng = 1.;
//   double erMjiecheng = 1.;
//   double fu1Njia1 = 1.;
//   double fu1Mjia1 = -1.;
//   double erPik = 2. * M_PI * k;
//   double erPik2N = 1.;
//   double erPik2Njian2Mjian1 = 1.;
//   double x2M = 1.;
//   double cosValue = cos(erPik * rc);
//   double sinValue = sin(erPik * rc);
  
//   for (int i = 0; i < n+1; ++i){
//     fu1Njia1 *= -1.;
//   }
//   for (int i = 1; i <= 2*n; ++i){
//     erNjiecheng *= i;
//   }
//   for (int i = 1; i <= 2*n-1; ++i){
//     erPik2Njian2Mjian1 *= erPik;
//   }
//   erPik2N = erPik2Njian2Mjian1 * erPik;
  
//   double sum = 0.;
//   for (int m = 0; m < n; ++m){
//     sum += (erPik2Njian2Mjian1 * fu1Mjia1 * erMjiecheng) / x2M * cosValue;
//     erPik2Njian2Mjian1 /= erPik;
//     erMjiecheng *= (2*m+1);
//     x2M *= rc;
//     sum += (erPik2Njian2Mjian1 * fu1Mjia1 * erMjiecheng) / x2M * sinValue;
//     erPik2Njian2Mjian1 /= erPik;
//     erMjiecheng *= (2*m+2);
//     x2M *= rc;
//     fu1Mjia1 *= -1;
//   }
//   double si, ci;
//   alglib::sinecosineintegrals (2. * M_PI * k * rc, si, ci);
//   // printf ("r:%e, si:%e, ci:%e\n", 2. * M_PI * k * rc, si, ci);

//   return - fu1Njia1 * (sum + (0.5 * M_PI - si) * erPik2N) / erNjiecheng ;
// }

  

// double ErrorEstimatorFFT_Corr::
// integral_an (const int & n_,
// 	     const double & k,
// 	     const double & rc)
// {
//   int n = (n_ >> 1);
//   int nterm = 20;
  
//   double sum = 0.;

//   double erMjiecheng = 1.;
//   double erNjiecheng = 1.;
//   double x2Mjia1 = 1.;
//   double erPik = 2. * M_PI * k;
//   double erPik1jia2Mjian2N = erPik;
//   double scale = -1.;
  
//   for (int i = 1; i <= 2*n; ++i){
//     // erMjiecheng *= i;
//     x2Mjia1 *= rc;
//   }
//   erNjiecheng = erMjiecheng;
//   x2Mjia1 *= rc;
//   for (int i = 2; i <= n; ++i){
//     scale *= -1.;
//   }
//   scale *= -1.;
//   double fu1N = scale * -1.;

//   double cosValue = cos(erPik * rc);
//   double sinValue = sin(erPik * rc);
  
//   for (int i = n; i <= n + nterm; ++i){
//     sum += (scale * erMjiecheng) / (erPik1jia2Mjian2N * x2Mjia1) * cosValue;
//     erPik1jia2Mjian2N *= erPik;
//     erMjiecheng *= (2 * i + 1);
//     x2Mjia1 *= rc;
//     sum += (scale * erMjiecheng) / (erPik1jia2Mjian2N * x2Mjia1) * sinValue;
//     erPik1jia2Mjian2N *= erPik;
//     erMjiecheng *= (2 * i + 2);
//     x2Mjia1 *= rc;
//     scale *= -1.;
//   }
  
//   return - fu1N * sum / erNjiecheng ;
    
//   // if (n_ == 1){
//   //   double si, ci;
//   //   alglib::sinecosineintegrals (2. * M_PI * k * rc, si, ci);
//   //   // printf ("r:%e, si:%e, ci:%e\n", 2. * M_PI * k * rc, si, ci);
//   //   return 0.5 * M_PI - si;
//   // }
//   // else if (n_ == 2){
//   //   double twopikrc = 2 * M_PI * k * rc;
//   //   double si, ci;
//   //   alglib::sinecosineintegrals (twopikrc, si, ci);
//   //   return sin(twopikrc) / rc - 2. * M_PI * k * ci;
//   // }
//   // else {
//   //   double twopikrc = 2 * M_PI * k * rc;
//   //   double rcnm1 = pow(rc, n_-1.);
//   //   double anm2 = integral_an (n_-2, k, rc);
//   //   // printf ("n-2:%d, an-2:%e\n", n-2, anm2);
//   //   return (1. / (n_-1.) * sin(twopikrc) / rcnm1 +
//   // 	    2. * M_PI * k / ((n_-1.) * (n_-2.) * rcnm1 / rc) * cos(twopikrc) - 
//   // 	    4. * M_PI * M_PI * k * k  / ((n_-1.) * (n_-2.)) * anm2);
//   // }
// }
