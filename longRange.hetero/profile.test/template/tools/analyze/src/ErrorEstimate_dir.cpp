#include "ErrorEstimate_dir.h"

#include <fstream>
#include <iostream>

ErrorEstimate_dir::
ErrorEstimate_dir ()
    : malloced (false)
{
}

ErrorEstimate_dir::
~ErrorEstimate_dir()
{
  freeAll();
}

void ErrorEstimate_dir::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void ErrorEstimate_dir::
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

void ErrorEstimate_dir::
freeAll ()
{
  if (malloced){
    free (rho1);
    free (rho2);
    free (s2kx);
    free (s2ky);
    free (s2kz);
    // free (k1rx);
    // free (k1ry);
    // free (k1rz);
    free (s1k);
    free (error1x);
    free (error1y);
    free (error1z);
    free (error1);
    free (error2);
    fftw_destroy_plan (p_forward_rho1);
    fftw_destroy_plan (p_forward_rho2);
    fftw_destroy_plan (p_backward_s2kx);
    fftw_destroy_plan (p_backward_s2ky);
    fftw_destroy_plan (p_backward_s2kz);
    fftw_destroy_plan (p_forward_s1k);
    fftw_destroy_plan (p_backward_error1x);
    fftw_destroy_plan (p_backward_error1y);
    fftw_destroy_plan (p_backward_error1z);
    fftw_destroy_plan (p_backward_error2);
  }
}

    
void ErrorEstimate_dir::
reinit (const double & beta_,
	const double & rcut_,
	const DensityProfile_PiecewiseConst & dp)
{
  freeAll ();
  
  beta = beta_;
  rcut = rcut_;
  K.x = dp.getNx();
  K.y = dp.getNy();
  K.z = dp.getNz();
  boxsize = dp.getBox();
  vecA.xy = vecA.yx = vecA.yz = vecA.zy = vecA.xz = vecA.zx = 0.;
  vecA.xx = dp.getBox()[0];
  vecA.yy = dp.getBox()[1];
  vecA.zz = dp.getBox()[2];
  calV();
  calAStar();

  nele = K.x * K.y * K.z;
  rho1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rho2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  s2kx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  s2ky = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  s2kz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  // k1rx = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  // k1ry = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  // k1rz = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  s1k  = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1x = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1y = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1z = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error1 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  error2 = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  p_forward_rho1 = fftw_plan_dft_3d (K.x, K.y, K.z, rho1, rho1, FFTW_FORWARD,  FFTW_MEASURE);
  p_forward_rho2 = fftw_plan_dft_3d (K.x, K.y, K.z, rho2, rho2, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_s2kx = fftw_plan_dft_3d (K.x, K.y, K.z, s2kx, s2kx, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s2ky = fftw_plan_dft_3d (K.x, K.y, K.z, s2ky, s2ky, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_s2kz = fftw_plan_dft_3d (K.x, K.y, K.z, s2kz, s2kz, FFTW_BACKWARD, FFTW_MEASURE);
  p_forward_s1k = fftw_plan_dft_3d (K.x, K.y, K.z, s1k, s1k, FFTW_FORWARD, FFTW_MEASURE);
  p_backward_error1x = fftw_plan_dft_3d (K.x, K.y, K.z, error1x, error1x, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1y = fftw_plan_dft_3d (K.x, K.y, K.z, error1y, error1y, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error1z = fftw_plan_dft_3d (K.x, K.y, K.z, error1z, error1z, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_error2 = fftw_plan_dft_3d (K.x, K.y, K.z, error2, error2, FFTW_BACKWARD, FFTW_MEASURE);
  
  malloced = true;

  printf ("# reinit: start build k mat 1 ...");
  fflush (stdout);
  makeS1k ();
  printf (" done\n");  
  fflush (stdout);
  printf ("# reinit: start build k mat 2 ...");
  fflush (stdout);
  makeS2k ();
  printf (" done\n");  
  fflush (stdout);
}


void ErrorEstimate_dir::
makeS1k ()
{
  // int count = 0;
  // double pref = 1.;
  IntVectorType myi, loc;
  
  for (myi.x = -K.x/2; myi.x < K.x - K.x/2; ++myi.x){
    if (myi.x < 0) loc.x = myi.x + K.x;
    else loc.x = myi.x;
    for (myi.y = -K.y/2; myi.y < K.y - K.y/2; ++myi.y){
      if (myi.y < 0) loc.y = myi.y + K.y;
      else loc.y = myi.y;
      for (myi.z = -K.z/2; myi.z < K.z - K.z/2; ++myi.z){
	if (myi.z < 0) loc.z = myi.z + K.z;
	else loc.z = myi.z;
	double k = sqrt(myi.x*myi.x/(boxsize[0]*boxsize[0]) +
			myi.y*myi.y/(boxsize[1]*boxsize[1]) +
			myi.z*myi.z/(boxsize[2]*boxsize[2]) );
	int posi = index3to1 (loc, K);
	if (k != 0){
	  double an;
	  an = integral_s1k_numerical (k, rcut, 30, 1e-12);
	  s1k[posi][0] = 2. / k * an;
	  s1k[posi][1] = 0.;
	}
	else {
	  double an;
	  an = integral_s1k1_numerical (k, rcut, 30, 1e-12);
	  s1k[posi][0] = 4. * M_PI * an;
	  s1k[posi][1] = 0.;
	}
	// printf ("%d %d %d      %e\n", myi.x, myi.y, myi.z, s1k[posi][0]);
      }
    }
  }
}

double ErrorEstimate_dir::
integral_s1k_numerical (const double & k,
			 const double & r1,
			 const double & r2,
			 const double & prec)
{
  f_s1k.k = k;
  f_s1k.beta = beta;
  return inte_s1k.cal_int (Integral1DInfo::Gauss4,
			   f_s1k,
			   r1, r2,
			   prec);
}

double ErrorEstimate_dir::
integral_s1k1_numerical (const double & k,
			 const double & r1,
			 const double & r2,
			 const double & prec)
{
  f_s1k1.k = k;
  f_s1k1.beta = beta;
  return inte_s1k1.cal_int (Integral1DInfo::Gauss4,
			    f_s1k1,
			    r1, r2,
			    prec);
}


void ErrorEstimate_dir::
makeS2k ()
{
  // double sigma6 = gsl_pow_6 (sigma);
  // double pre = 8. * epsilon * sigma6;
  IntVectorType myi, loc;
  
  double kx, ky, kz;
  for (myi.x = -K.x/2; myi.x < K.x - K.x/2; ++myi.x){
    kx = myi.x / boxsize[0];
    if (myi.x < 0) loc.x = myi.x + K.x;
    else loc.x = myi.x;
    for (myi.y = -K.y/2; myi.y < K.y - K.y/2; ++myi.y){
      ky = myi.y / boxsize[1];
      if (myi.y < 0) loc.y = myi.y + K.y;
      else loc.y = myi.y;
      for (myi.z = -K.z/2; myi.z < K.z - K.z/2; ++myi.z){
	kz = myi.z / boxsize[2];
	if (myi.z < 0) loc.z = myi.z + K.z;
	else loc.z = myi.z;
	double k = sqrt(kx*kx + ky*ky + kz*kz);
	int posi = index3to1 (loc, K);
	if (k != 0){
	  double an;
	  an = integral_s2k_numerical (k, rcut, 30, 1e-16) / k;
	  // printf ("%d %d %d    %e    %e\n", myi.x, myi.y, myi.z, k, an * k);
	  s2kx[posi][1] = 2 * M_PI * kx * an;
	  s2ky[posi][1] = 2 * M_PI * ky * an;
	  s2kz[posi][1] = 2 * M_PI * kz * an;
	  double erPiKR = 2. * M_PI * k * rcut;
	  double size = 2. * M_PI * rcut * erfc(beta*rcut) *
	      (2. * cos(erPiKR) / erPiKR - 2. * sin(erPiKR) / (erPiKR*erPiKR));
	  size /= k;
	  s2kx[posi][0] = 0.;
	  s2kx[posi][1] += size * kx;
	  s2ky[posi][0] = 0.;
	  s2ky[posi][1] += size * ky;
	  s2kz[posi][0] = 0.;
	  s2kz[posi][1] += size * kz;
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

double ErrorEstimate_dir::
integral_s2k_numerical (const double & k,
			const double & r1,
			const double & r2,
			const double & prec)
{
  f_s2k.k = k;
  f_s2k.beta = beta;
  return inte_s2k.cal_int (Integral1DInfo::Gauss4,
			   f_s2k,
			   r1, r2,
			   prec);
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


void ErrorEstimate_dir::
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

  array_multiply (error1x, nele, s2kx, rho1);
  array_multiply (error1y, nele, s2ky, rho1);
  array_multiply (error1z, nele, s2kz, rho1);
  array_multiply (error2, nele, s1k, rho2);

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



void ErrorEstimate_dir::    
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


void ErrorEstimate_dir::
print_meanf (const std::string & file,
	     const DensityProfile_PiecewiseConst & dp) const
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  fprintf (fp, "# 1          2        3       4                5        6       7\n");
  fprintf (fp, "# x  meanFx(r)  meanFxR meanFxI  meanFyR meanFyI  meanFzR meanFzI\n");
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
    
    fprintf (fp, "%f %e     %e %e  %e %e  %e %e\n",
	     (i + 0.5) * vecA.xx / K.x,
	     error1x[index][0],
	     error1x[index][0] * scalor,
	     error1x[index][1] * scalor,
	     error1y[index][0] * scalor,
	     error1y[index][1] * scalor,
	     error1z[index][0] * scalor,
	     error1z[index][1] * scalor

	);
  }
  fclose (fp);
}
  


	
