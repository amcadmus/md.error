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

  nele = nx * ny * nz;
  rho_real    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rho_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  p_forward_rho  = fftw_plan_dft_3d (nx, ny, nz, rho_real, rho_complex, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_rho = fftw_plan_dft_3d (nx, ny, nz, rho_real, rho_complex, FFTW_BACKWARD, FFTW_MEASURE);
  for (int i = 0; i < nele; ++i){
    rho_real[i][1] = 0.;
  }
}

ErrorEstimatorFFT_Corr::
~ErrorEstimatorFFT_Corr () 
{
  fftw_destroy_plan (p_forward_rho);
  fftw_destroy_plan (p_backward_rho);
  fftw_free (rho_real);
  fftw_free (rho_complex);
}

void ErrorEstimatorFFT_Corr::
estimate (const DensityProfile_Corr_PiecewiseConst & dp,
	  const int & corrBond)
{
  for (int i = 0; i < nele; ++i){
    rho_real[i][0] = dp.getMean(i);
  }
  fftw_execute (p_forward_rho);
  fftw_execute (p_backward_rho);
  for (int i = 0; i < nele; ++i){
    if (fabs(rho_real[i][0] - dp.getMean(i)) > 1e-10 || fabs(rho_real[i][1]) > 1e-10){
      printf ("error !!!\n");
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
    // fprintf (fp, "%f %e\n", (i + 0.5) * hx, profile[index3to1(i,0,0)]);
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
      // fprintf (fp, "%f %f %e\n", vx, vy, profile[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}




