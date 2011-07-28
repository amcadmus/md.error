#include "ForceCorr.h"

ForceCorr::ForceCorr ()
    : malloced (false)
{
}

ForceCorr::~ForceCorr ()
{
  freeAll();
}

static void
mallocArrayComplex (fftw_complex ** a,
		    const unsigned & nele)
{
  (*a) = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
}

static void
freeArrayComplex (fftw_complex ** a)
{
  free ((*a));
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


void ForceCorr::
freeAll ()
{
  if (malloced) {
    free (rhok);
    free (rhor);
    freeArrayComplex (&s2kx);
    freeArrayComplex (&s2ky);
    freeArrayComplex (&s2kz);
    freeArrayComplex (&error2kx);
    freeArrayComplex (&error2ky);
    freeArrayComplex (&error2kz);
    freeArrayComplex (&error2rx);
    freeArrayComplex (&error2ry);
    freeArrayComplex (&error2rz);
    fftw_destroy_plan (p_forward_rho);
    fftw_destroy_plan (p_backward_error2x);
    fftw_destroy_plan (p_backward_error2y);
    fftw_destroy_plan (p_backward_error2z);
    malloced = false;
  }
}


void ForceCorr::
reinit (const double & rc_,
	DensityProfile_PiecewiseConst & dp)
{
  freeAll ();

  rc = rc_;
  boxsize = dp.getBox();
  nx = dp.getNx();
  ny = dp.getNy();
  nz = dp.getNz();
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;
  printf ("# ForceCorr nx, ny, nz are %d %d %d\n", nx, ny, nz);

  rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  mallocArrayComplex (&s2kx, nele);
  mallocArrayComplex (&s2ky, nele);
  mallocArrayComplex (&s2kz, nele);
  mallocArrayComplex (&error2kx, nele);
  mallocArrayComplex (&error2ky, nele);
  mallocArrayComplex (&error2kz, nele);
  mallocArrayComplex (&error2rx, nele);
  mallocArrayComplex (&error2ry, nele);
  mallocArrayComplex (&error2rz, nele);

  p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_PATIENT);
  p_backward_error2x = fftw_plan_dft_3d (nx, ny, nz, error2kx, error2rx, FFTW_BACKWARD,  FFTW_PATIENT);
  p_backward_error2y = fftw_plan_dft_3d (nx, ny, nz, error2ky, error2ry, FFTW_BACKWARD,  FFTW_PATIENT);
  p_backward_error2z = fftw_plan_dft_3d (nx, ny, nz, error2kz, error2rz, FFTW_BACKWARD,  FFTW_PATIENT);
  
  malloced = true;

  printf ("# reinit: start build k mat 2 ...");
  fflush (stdout);
  makeS2k (1., 1.);
  printf (" done\n");  
  fflush (stdout);
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

void ForceCorr::
calError (const DensityProfile_PiecewiseConst & dp)
{
  for (int i = 0; i < nele; ++i){
    rhor[i][0] = dp.getProfile(i);
    rhor[i][1] = 0.;
  }
  fftw_execute (p_forward_rho);
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  double scale = volume / nele;
  for (int i = 0; i < nele; ++i){
    rhok[i][0] *= scale;
    rhok[i][1] *= scale;
  }  

  array_multiply (error2kx, nele, s2kx, rhok);
  array_multiply (error2ky, nele, s2ky, rhok);
  array_multiply (error2kz, nele, s2kz, rhok);

  fftw_execute (p_backward_error2x);
  fftw_execute (p_backward_error2y);
  fftw_execute (p_backward_error2z);

  for (int i = 0; i < nele; ++i){
    error2rx[i][0] /= volume;
    error2ry[i][0] /= volume;
    error2rz[i][0] /= volume;
    error2rx[i][1] /= volume;
    error2ry[i][1] /= volume;
    error2rz[i][1] /= volume;
  }
}



double ForceCorr::
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


double ForceCorr::
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


void ForceCorr::
makeS2k (const double & epsilon,
	 const double & sigma) 
{
  double sigma6 = gsl_pow_6 (sigma);
  double pre = 8. * epsilon * sigma6;
  
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
	int posi = index3to1 (locx, locy, locz);
	if (k != 0){
	  double an;
	  an = integral_an5_numerical (k, rc, 30, 1e-16) * pre / k;
	  s2kx[posi][1] = 2 * M_PI * kx * an;
	  s2ky[posi][1] = 2 * M_PI * ky * an;
	  s2kz[posi][1] = 2 * M_PI * kz * an;
	}
      }
    }
  }


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
	double rc6 = gsl_pow_6 (rc);
	if (k != 0){
	  double erPiKR = 2. * M_PI * k * rc;
	  double size = - 2. * M_PI * rc * rc * (4. * epsilon * sigma6 / rc6) *
	      (2. * cos(erPiKR) / erPiKR - 2. * sin(erPiKR) / (erPiKR*erPiKR));
	  size /= k;
	  s2kx[posi][0] = 0.;
	  s2kx[posi][1] += size * kx;
	  s2ky[posi][0] = 0.;
	  s2ky[posi][1] += size * ky;
	  s2kz[posi][0] = 0.;
	  s2kz[posi][1] += size * kz;
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



void ForceCorr::    
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
    fprintf (fp, "%f %e %e %e %e %e %e\n",
	     (i + 0.5) * hx,
	     error2rx[index3to1(i,0,0)][0], error2rx[index3to1(i,0,0)][1],
	     error2ry[index3to1(i,0,0)][0], error2ry[index3to1(i,0,0)][1],
	     error2rz[index3to1(i,0,0)][0], error2rz[index3to1(i,0,0)][1]
	     );
  }
  fclose (fp);
}


// void ForceCorr::    
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
//     	sum += result_error[index3to1(i,j,k)] *
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
//     fprintf (fp, "%f %e\n",
// 	     (i + 0.5) * hx,
// 	     sum);
//   }
//   fclose (fp);
// }


  


  

// void ForceCorr::
// save_rc (const std::string & file) const
// {
//   FILE * fp = fopen (file.c_str(), "w");
//   if (fp == NULL){
//     fprintf (stderr, "cannot open file %s\n", file.c_str());
//     exit (1);
//   }

//   fprintf (fp, "%d ", nrc);
//   for (int i = 0; i < nrc; ++i){
//     fprintf (fp, "%f ", rcList[i]);
//   }
//   fprintf (fp, "\n%f %f %f\n", boxsize[0], boxsize[1], boxsize[2]);
//   fprintf (fp, "%d %d %d\n", nx, ny, nz);
//   for (int i = 0; i < nele; ++i){
//     fprintf (fp, "%d ", rcutIndex[i]);
//   }
//   fprintf(fp, "\n");
  
//   // fwrite (&nrc, sizeof(int), 1, fp);
//   // fwrite (rcList, sizeof(double), nrc, fp);
//   // fwrite (boxsize, sizeof(double), 3, fp);
//   // fwrite (&nx, sizeof(int), 1, fp);
//   // fwrite (&ny, sizeof(int), 1, fp);
//   // fwrite (&nz, sizeof(int), 1, fp);
//   // fwrite (rcutIndex, sizeof(int), nele, fp);
  
//   fclose (fp);
// }

// void ForceCorr::
// write_rc (const std::string & file) const
// {
//   FILE * fp = fopen (file.c_str(), "w");
//   if (fp == NULL){
//     std::cerr << "cannot open file " << file << std::endl;
//     exit(1);
//   }

//   fprintf (fp, "%d %d %d\n", nx, ny, nz);
//   for (int i = 0; i < nele; ++i){
//     fprintf (fp, "%f ", rcut[i]);
//   }

//   fclose (fp);
// }


