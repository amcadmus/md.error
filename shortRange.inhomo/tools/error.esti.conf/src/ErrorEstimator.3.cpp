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

  nele = nx * ny * nz;
  f2k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s1r  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  s1k  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  p_forward_rho  = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_rho = fftw_plan_dft_3d (nx, ny, nz, rhok, rhor, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_s1   = fftw_plan_dft_3d (nx, ny, nz, s1r,  s1k,  FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_s1  = fftw_plan_dft_3d (nx, ny, nz, s1k,  s1r,  FFTW_BACKWARD, FFTW_MEASURE);
  for (int i = 0; i < nele; ++i){
    rhor[i][1] = 0.;
  }
}

ErrorEstimatorFFT_Corr::
~ErrorEstimatorFFT_Corr () 
{
  fftw_destroy_plan (p_forward_s1);
  fftw_destroy_plan (p_backward_s1);
  fftw_destroy_plan (p_forward_rho);
  fftw_destroy_plan (p_backward_rho);
  fftw_free (f2k);
  fftw_free (s1r);
  fftw_free (s1k);
  fftw_free (rhor);
  fftw_free (rhok);
}


double ErrorEstimatorFFT_Corr::
integral_an13_numerical (const double & k,
			 const double & rc)
{
  f_inte.k = k;
  return inte.cal_int (Integral1DInfo::Gauss4,
		       f_inte,
		       rc, 30,
		       1e-16);
}


double ErrorEstimatorFFT_Corr::
integral_an1 (const int & n_,
	      const double & k,
	      const double & rc)
{
  int n = (n_ >> 1);
  double erNjiecheng = 1.;
  double erMjiecheng = 1.;
  double fu1Njia1 = 1.;
  double fu1Mjia1 = -1.;
  double erPik = 2. * M_PI * k;
  double erPik2N = 1.;
  double erPik2Njian2Mjian1 = 1.;
  double x2M = 1.;
  double cosValue = cos(erPik * rc);
  double sinValue = sin(erPik * rc);
  
  for (int i = 0; i < n+1; ++i){
    fu1Njia1 *= -1.;
  }
  for (int i = 1; i <= 2*n; ++i){
    erNjiecheng *= i;
  }
  for (int i = 1; i <= 2*n-1; ++i){
    erPik2Njian2Mjian1 *= erPik;
  }
  erPik2N = erPik2Njian2Mjian1 * erPik;
  
  double sum = 0.;
  for (int m = 0; m < n; ++m){
    sum += (erPik2Njian2Mjian1 * fu1Mjia1 * erMjiecheng) / x2M * cosValue;
    erPik2Njian2Mjian1 /= erPik;
    erMjiecheng *= (2*m+1);
    x2M *= rc;
    sum += (erPik2Njian2Mjian1 * fu1Mjia1 * erMjiecheng) / x2M * sinValue;
    erPik2Njian2Mjian1 /= erPik;
    erMjiecheng *= (2*m+2);
    x2M *= rc;
    fu1Mjia1 *= -1;
  }
  double si, ci;
  alglib::sinecosineintegrals (2. * M_PI * k * rc, si, ci);
  // printf ("r:%e, si:%e, ci:%e\n", 2. * M_PI * k * rc, si, ci);

  return - fu1Njia1 * (sum + (0.5 * M_PI - si) * erPik2N) / erNjiecheng ;
}

  

double ErrorEstimatorFFT_Corr::
integral_an (const int & n_,
	     const double & k,
	     const double & rc)
{
  int n = (n_ >> 1);
  int nterm = 20;
  
  double sum = 0.;

  double erMjiecheng = 1.;
  double erNjiecheng = 1.;
  double x2Mjia1 = 1.;
  double erPik = 2. * M_PI * k;
  double erPik1jia2Mjian2N = erPik;
  double scale = -1.;
  
  for (int i = 1; i <= 2*n; ++i){
    // erMjiecheng *= i;
    x2Mjia1 *= rc;
  }
  erNjiecheng = erMjiecheng;
  x2Mjia1 *= rc;
  for (int i = 2; i <= n; ++i){
    scale *= -1.;
  }
  scale *= -1.;
  double fu1N = scale * -1.;

  double cosValue = cos(erPik * rc);
  double sinValue = sin(erPik * rc);
  
  for (int i = n; i <= n + nterm; ++i){
    sum += (scale * erMjiecheng) / (erPik1jia2Mjian2N * x2Mjia1) * cosValue;
    erPik1jia2Mjian2N *= erPik;
    erMjiecheng *= (2 * i + 1);
    x2Mjia1 *= rc;
    sum += (scale * erMjiecheng) / (erPik1jia2Mjian2N * x2Mjia1) * sinValue;
    erPik1jia2Mjian2N *= erPik;
    erMjiecheng *= (2 * i + 2);
    x2Mjia1 *= rc;
    scale *= -1.;
  }
  
  return - fu1N * sum / erNjiecheng ;
    
  // if (n_ == 1){
  //   double si, ci;
  //   alglib::sinecosineintegrals (2. * M_PI * k * rc, si, ci);
  //   // printf ("r:%e, si:%e, ci:%e\n", 2. * M_PI * k * rc, si, ci);
  //   return 0.5 * M_PI - si;
  // }
  // else if (n_ == 2){
  //   double twopikrc = 2 * M_PI * k * rc;
  //   double si, ci;
  //   alglib::sinecosineintegrals (twopikrc, si, ci);
  //   return sin(twopikrc) / rc - 2. * M_PI * k * ci;
  // }
  // else {
  //   double twopikrc = 2 * M_PI * k * rc;
  //   double rcnm1 = pow(rc, n_-1.);
  //   double anm2 = integral_an (n_-2, k, rc);
  //   // printf ("n-2:%d, an-2:%e\n", n-2, anm2);
  //   return (1. / (n_-1.) * sin(twopikrc) / rcnm1 +
  // 	    2. * M_PI * k / ((n_-1.) * (n_-2.) * rcnm1 / rc) * cos(twopikrc) - 
  // 	    4. * M_PI * M_PI * k * k  / ((n_-1.) * (n_-2.)) * anm2);
  // }
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
  // fftw_execute (p_backward_rho);
  // for (int i = 0; i < nele; ++i){
  //   printf ("%d %e %e %e\n", i, rhor[i][0] / nele, dp.getMean(i), rhor[i][1] / nele);
  //   if (fabs(rhor[i][0] / nele - dp.getMean(i)) > 1e-10 || fabs(rhor[i][1]) > 1e-10){
  //     printf ("%d error !!! %e %e %e\n",
  // 	      i, rhor[i][0] / nele, dp.getMean(i), fabs(rhor[i][0]/nele - dp.getMean(i)));
  //   }
  // }
  for (int i = 0; i < nele; ++i){
    rhok[i][0] *= scale;
    rhok[i][1] *= scale;
  }
  
  makef2k (1., rc);
  for (int i = 0; i < nele; ++i){
    s1k[i][0] = f2k[i][0] * rhok[i][0] - f2k[i][1] * rhok[i][1];
    s1k[i][1] = f2k[i][1] * rhok[i][0] + f2k[i][0] * rhok[i][1];
  }
  
  fftw_execute (p_backward_s1);
  for (int i = 0; i < nele; ++i){
    s1r[i][0] /= volume;
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
	  // printf ("n:%d, an:%e\n", 13, an);
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
    fprintf (fp, "%f %e %e\n",
	     (i + 0.5) * hx,
	     s1r[index3to1(i,0,0)][0],
	     s1r[index3to1(i,0,0)][1]);
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




