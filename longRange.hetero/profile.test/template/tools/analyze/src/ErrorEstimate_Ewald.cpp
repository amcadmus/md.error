#include "ErrorEstimate_Ewald.h"

#include <fstream>
#include <iostream>

ErrorEstimate_Ewald::
ErrorEstimate_Ewald ()
    : malloced (false)
{
}

ErrorEstimate_Ewald::
~ErrorEstimate_Ewald()
{
  freeAll();
}

void ErrorEstimate_Ewald::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void ErrorEstimate_Ewald::
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

void ErrorEstimate_Ewald::
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
    free (k2m);
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
    fftw_destroy_plan (p_forward_k2);
    fftw_destroy_plan (p_backward_error1x);
    fftw_destroy_plan (p_backward_error1y);
    fftw_destroy_plan (p_backward_error1z);
    fftw_destroy_plan (p_backward_error2);
  }
}

    
void ErrorEstimate_Ewald::
reinit (const double & beta_,
	const IntVectorType & Kmax_,
	const DensityProfile_PiecewiseConst & dp)
{
  freeAll ();
  
  beta = beta_;
  Kmax = Kmax_;
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
  k2m  = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
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
  p_forward_k2 = fftw_plan_dft_3d (K.x, K.y, K.z, k2m, k2m, FFTW_FORWARD, FFTW_MEASURE);
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

void ErrorEstimate_Ewald::
calKernel()
{
  double scalor = 1./(2. * M_PI * volume);
  int nele = K.x * K.y * K.z;
  IntVectorType idx;
  for (idx.x = -(K.x >> 1); idx.x <= (K.x >> 1); ++idx.x){
    for (idx.y = -(K.y >> 1); idx.y <= (K.y >> 1); ++idx.y){
      for (idx.z = -(K.z >> 1); idx.z <= (K.z >> 1); ++idx.z){
	IntVectorType posiIdx (idx);
	if (posiIdx.x < 0) posiIdx.x += K.x;
	if (posiIdx.y < 0) posiIdx.y += K.y;
	if (posiIdx.z < 0) posiIdx.z += K.z;
	unsigned index = index3to1 (posiIdx, K);
	if (idx.x >= -(Kmax.x >> 1) && idx.x <= (Kmax.x >> 1) &&
	    idx.y >= -(Kmax.y >> 1) && idx.y <= (Kmax.y >> 1) &&
	    idx.z >= -(Kmax.z >> 1) && idx.z <= (Kmax.z >> 1) ){
	  k1mx[index][0] = k1mx[index][1] = 0.;
	  k1my[index][0] = k1my[index][1] = 0.;
	  k1mz[index][0] = k1mz[index][1] = 0.;
	}
	else {
	  VectorType mm;
	  mm.x = idx.x * vecAStar.xx + idx.y * vecAStar.yx + idx.z * vecAStar.zx;
	  mm.y = idx.x * vecAStar.xy + idx.y * vecAStar.yy + idx.z * vecAStar.zy;
	  mm.z = idx.x * vecAStar.xz + idx.y * vecAStar.yz + idx.z * vecAStar.zz;
	  double m2 = (mm.x*mm.x + mm.y*mm.y + mm.z*mm.z);
	  double fm = kernel_rm1_rec_f (m2, beta) * scalor;
	  k1mx[index][0] = k1my[index][0] = k1mz[index][0] = 0.;
	  k1mx[index][1] = - 4. * M_PI * fm * mm.x;
	  k1my[index][1] = - 4. * M_PI * fm * mm.y;
	  k1mz[index][1] = - 4. * M_PI * fm * mm.z;
	}
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
    k2m[i][0] = (
	(k1rx[i][0] * k1rx[i][0] - k1rx[i][1] * k1rx[i][1]) +
	(k1ry[i][0] * k1ry[i][0] - k1ry[i][1] * k1ry[i][1]) +
	(k1rz[i][0] * k1rz[i][0] - k1rz[i][1] * k1rz[i][1]) ) * volume / double(nele);
    k2m[i][1] = 0.;
  }

  fftw_execute (p_forward_k2);
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


void ErrorEstimate_Ewald::
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
  array_multiply (error2, nele, k2m, rho2);

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


// void AdaptRCut::    
// print_error_avg (const DensityProfile_PiecewiseConst & dp,
// 		 const std::string & file) const 
// {
//   FILE * fp = fopen (file.c_str(), "w");
//   if (fp == NULL){
//     std::cerr << "cannot open file " << file << std::endl;
//     exit(1);
//   }

//   for (int i = 0; i < nx; ++i){
//     double sum = 0.;
//     double sumprofile = 0.;
//     for (int j = 0; j < ny; ++j){
//       for (int k = 0; k < nz; ++k){
//     	sum +=
// 	    result_error[index3to1(i,j,k)] *
// 	    result_error[index3to1(i,j,k)] *
// 	    dp.getProfile(i,j,k);
// 	sumprofile += dp.getProfile (i,j,k);
//       }
//     }
//     sum /= sumprofile;
//     sum = sqrt(sum);
//     // fprintf (fp, "%f %e %e\n",
//     // 	     (i + 0.5) * hx,
//     // 	     error[4][index3to1(i,0,0)][0],
//     // 	     error[4][index3to1(i,0,0)][1]
//     // 	);
//     // fprintf (fp, "%f", (i + 0.5) * hx);
//     // fprintf (fp, " %e", sum);
//     fprintf (fp, "%f %e\n",
// 	     (i + 0.5) * hx,
// 	     sum);
//   }
//   fclose (fp);
// }

void ErrorEstimate_Ewald::    
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


void ErrorEstimate_Ewald::
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
    
    fprintf (fp, "%f     %e %e %e %e\n",
	     (i + 0.5) * vecA.xx / K.x,
	     error1x[index][0],
	     error1x[index][0] * scalor,
	     error1y[index][0] * scalor,
	     error1z[index][0] * scalor
	);
  }
  fclose (fp);
}
  


	
