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
freeArrayComplex (fftw_complex *** a,
		  const unsigned & nrc)
{
  for (unsigned i = 0; i < nrc; ++i){
    free ((*a)[i]);
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
    free (rhok);
    free (rhor);
    freeArrayComplex (&s1k,  nrc);
    freeArrayComplex (&s2kx, nrc);
    freeArrayComplex (&s2ky, nrc);
    freeArrayComplex (&s2kz, nrc);
    freeArrayComplex (&error1k,  nrc);
    freeArrayComplex (&error2kx, nrc);
    freeArrayComplex (&error2ky, nrc);
    freeArrayComplex (&error2kz, nrc);
    freeArrayComplex (&error1r,  nrc);
    freeArrayComplex (&error2rx, nrc);
    freeArrayComplex (&error2ry, nrc);
    freeArrayComplex (&error2rz, nrc);
    freeArrayComplex (&error, nrc);
    free (rcut);
    free (rcutIndex);
    free (result_error);
    fftw_destroy_plan (p_forward_rho);
    for (int count = 0; count < nrc; ++count){
      fftw_destroy_plan (p_backward_error1[count]);
      fftw_destroy_plan (p_backward_error2x[count]);
      fftw_destroy_plan (p_backward_error2y[count]);
      fftw_destroy_plan (p_backward_error2z[count]);
    }
    free (p_backward_error1);
    free (p_backward_error2x);
    free (p_backward_error2y);
    free (p_backward_error2z);
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
  
  boxsize = dp.getBox();
  nx = dp.getNx();
  ny = dp.getNy();
  nz = dp.getNz();
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;

  rcList.clear();  
  for (double tmprc = rcmin; tmprc <= rcmax; tmprc += rcstep){
    rcList.push_back(tmprc);
  }
  nrc = rcList.size();
  
  rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  mallocArrayComplex (&s1k,  nrc, nele);
  mallocArrayComplex (&s2kx, nrc, nele);
  mallocArrayComplex (&s2ky, nrc, nele);
  mallocArrayComplex (&s2kz, nrc, nele);
  mallocArrayComplex (&error1k,  nrc, nele);
  mallocArrayComplex (&error2kx, nrc, nele);
  mallocArrayComplex (&error2ky, nrc, nele);
  mallocArrayComplex (&error2kz, nrc, nele);
  mallocArrayComplex (&error1r,  nrc, nele);
  mallocArrayComplex (&error2rx, nrc, nele);
  mallocArrayComplex (&error2ry, nrc, nele);
  mallocArrayComplex (&error2rz, nrc, nele);
  mallocArrayComplex (&error, nrc, nele);
  rcut = (double *) malloc (sizeof(double) * nele);
  rcutIndex = (int *) malloc (sizeof(int) * nele);
  result_error = (double *) malloc (sizeof(double) * nele);

  p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_error1  = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2x = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2y = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2z = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  for (int count = 0; count < nrc; ++count){
    p_backward_error1[count] = fftw_plan_dft_3d (nx, ny, nz, error1k[count], error1r[count], FFTW_BACKWARD,  FFTW_MEASURE);
    p_backward_error2x[count] = fftw_plan_dft_3d (nx, ny, nz, error2kx[count], error2rx[count], FFTW_BACKWARD,  FFTW_MEASURE);
    p_backward_error2y[count] = fftw_plan_dft_3d (nx, ny, nz, error2ky[count], error2ry[count], FFTW_BACKWARD,  FFTW_MEASURE);
    p_backward_error2z[count] = fftw_plan_dft_3d (nx, ny, nz, error2kz[count], error2rz[count], FFTW_BACKWARD,  FFTW_MEASURE);
  }
  
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

static void 
array_multiply (fftw_complex ** a,
		const int nrc,
		const int nele,
		fftw_complex ** b,
		fftw_complex * c) 
{
  for (int count = 0; count < nrc; ++count){
    for (int ii = 0; ii < nele; ++ii){
      // double tmpr, tmpi;
      a[count][ii][0] =
	  c[ii][0] * b[count][ii][0] - c[ii][1] * b[count][ii][1];
      a[count][ii][1] =
	  c[ii][0] * b[count][ii][1] + c[ii][1] * b[count][ii][0];
    }
  }
}


void AdaptRCut::
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

  array_multiply (error1k,  nrc, nele, s1k,  rhok);
  array_multiply (error2kx, nrc, nele, s2kx, rhok);
  array_multiply (error2ky, nrc, nele, s2ky, rhok);
  array_multiply (error2kz, nrc, nele, s2kz, rhok);

  for (int count = 0; count < nrc; ++count){
    fftw_execute (p_backward_error1[count]);
    fftw_execute (p_backward_error2x[count]);
    fftw_execute (p_backward_error2y[count]);
    fftw_execute (p_backward_error2z[count]);
  }

  for (int count = 0; count < nrc; ++count){
    for (int i = 0; i < nele; ++i){
      error1r[count][i][0] /= volume;
      error1r[count][i][1] /= volume;
      error2rx[count][i][0] /= volume;
      error2ry[count][i][0] /= volume;
      error2rz[count][i][0] /= volume;
      error2rx[count][i][1] /= volume;
      error2ry[count][i][1] /= volume;
      error2rz[count][i][1] /= volume;
      error[count][i][0] =
	  sqrt (
	      error1r[count][i][0] +
	      error2rx[count][i][0] * error2rx[count][i][0] +
	      error2ry[count][i][0] * error2ry[count][i][0] +
	      error2rz[count][i][0] * error2rz[count][i][0] );
      error[count][i][1] =
	  sqrt (
	      error1r[count][i][1] +
	      error2rx[count][i][1] * error2rx[count][i][1] +
	      error2ry[count][i][1] * error2ry[count][i][1] +
	      error2rz[count][i][1] * error2rz[count][i][1] );
    }
  }
}



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
	    s2ky[count][posi][1] += size * kx;
	    s2kz[count][posi][0] = 0.;
	    s2kz[count][posi][1] += size * kz;
	    // if (ix == -nx/2 ||
	    // 	iy == -ny/2 ||
	    // 	iz == -nz/2){
	    //   s2kx[count][posi][0] = 0.;
	    //   s2kx[count][posi][1] = 0.;
	    //   s2ky[count][posi][0] = 0.;
	    //   s2ky[count][posi][1] = 0.;
	    //   s2kz[count][posi][0] = 0.;
	    //   s2kz[count][posi][1] = 0.;
	    // }
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

void AdaptRCut::
calRCutOnePoint (const double & prec,
		 const unsigned & idx)
{
  unsigned posia = 0;
  unsigned posib = nrc-1;
  double rca = rcList[posia];
  double rcb = rcList[posib];
  double errora = error[posia][idx][0];
  double errorb = error[posib][idx][0];
  double diffa = errora - prec;
  double diffb = errorb - prec;
  if (diffa <= 0){
    rcut[idx] = rca;
    rcutIndex[idx] = posia;
    result_error[idx] = errora;
  }
  else if (diffb >= 0){
    rcut[idx] = rcb;
    rcutIndex[idx] = posib;
    result_error[idx] = errorb;
  }
  else {
    while (posib - posia > 1) {
      unsigned posic = (posib + posia) / 2;
      double rcc = rcList[posic];
      double errorc = error[posic][idx][0];
      if (errorc > prec){
	posia = posic;
	rca = rcc;
	errora = errorc;
      }
      else {
	posib = posic;
	rcb = rcc;
	errorb = errorc;
      }
    }
    rcut[idx] = rcb;
    rcutIndex[idx] = posib;
    result_error[idx] = errorb;
  }
}   

void AdaptRCut::
calRCut  (const double & prec)
{
  for (int i = 0; i < nele; ++i){
    calRCutOnePoint (prec, i);
  }
}

void AdaptRCut::
uniformRCut (const double & rc)
{
  int minIndex = 0;
  double minDiff = fabs(rc - rcList.front());
  
  for (unsigned i = 1; i < rcList.size(); ++i){
    double diff = fabs(rc - rcList[i]);
    if (diff < minDiff){
      minDiff = diff;
      minIndex = i;
    }
  }

  for (int ii = 0; ii < nele; ++ii){
    rcut[ii] = rcList[minIndex];
    rcutIndex[ii] = minIndex;
    result_error[ii] = error[minIndex][ii][0];
  }
}


void AdaptRCut::
load_rc (const std::string & file,
	 const DensityProfile_PiecewiseConst & dp)
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
  boxsize.resize (3);
  fscanf (fp, "%lf %lf %lf", &boxsize[0], &boxsize[1], &boxsize[2]);
  fscanf (fp, "%d %d %d", &nx, &ny, &nz);
  nele = nx * ny * nz;
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;

  freeAll();

  rhor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  rhok = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nele);
  mallocArrayComplex (&s1k,  nrc, nele);
  mallocArrayComplex (&s2kx, nrc, nele);
  mallocArrayComplex (&s2ky, nrc, nele);
  mallocArrayComplex (&s2kz, nrc, nele);
  mallocArrayComplex (&error1k,  nrc, nele);
  mallocArrayComplex (&error2kx, nrc, nele);
  mallocArrayComplex (&error2ky, nrc, nele);
  mallocArrayComplex (&error2kz, nrc, nele);
  mallocArrayComplex (&error1r,  nrc, nele);
  mallocArrayComplex (&error2rx, nrc, nele);
  mallocArrayComplex (&error2ry, nrc, nele);
  mallocArrayComplex (&error2rz, nrc, nele);
  mallocArrayComplex (&error, nrc, nele);
  rcut = (double *) malloc (sizeof(double) * nele);
  rcutIndex = (int *) malloc (sizeof(int) * nele);
  result_error = (double *) malloc (sizeof(double) * nele);

  for (int i = 0; i < nele; ++i){
    fscanf (fp, "%d ", &rcutIndex[i]);
    rcut[i] = rcList[rcutIndex[i]];
  }
  
  p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_PATIENT);
  p_backward_error1  = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2x = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2y = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  p_backward_error2z = (fftw_plan*) malloc (sizeof(fftw_plan) * nrc);
  for (int count = 0; count < nrc; ++count){
    p_backward_error1[count] = fftw_plan_dft_3d (nx, ny, nz, error1k[count], error1r[count], FFTW_BACKWARD,  FFTW_PATIENT);
    p_backward_error2x[count] = fftw_plan_dft_3d (nx, ny, nz, error2kx[count], error2rx[count], FFTW_BACKWARD,  FFTW_PATIENT);
    p_backward_error2y[count] = fftw_plan_dft_3d (nx, ny, nz, error2ky[count], error2ry[count], FFTW_BACKWARD,  FFTW_PATIENT);
    p_backward_error2z[count] = fftw_plan_dft_3d (nx, ny, nz, error2kz[count], error2rz[count], FFTW_BACKWARD,  FFTW_PATIENT);
  }
  
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

  calError (dp);

  for (int ii = 0; ii < nele; ++ii){
    result_error[ii] = error[rcutIndex[ii]][ii][0];
  }
}



// void AdaptRCut::
// print_error_uni_rc (const std::string & file,
// 		    const double & rcut) const
// {
//   int minIndex = 0;
//   double minDiff = fabs(rcut - rcList.front());
  
//   for (unsigned i = 1; i < rcList.size(); ++i){
//     double diff = fabs(rcut - rcList[i]);
//     if (diff < minDiff){
//       minDiff = diff;
//       minIndex = i;
//     }
//   }

//   FILE * fp = fopen (file.c_str(), "w");
//   if (fp == NULL){
//     std::cerr << "cannot open file " << file << std::endl;
//     exit(1);
//   }
//   fprintf (fp, "# error at rcut %f\n", rcut);
//   for (int i = 0; i < nx; ++i){
//     fprintf (fp, "%f %e\n",
// 	     (i + 0.5) * hx,
// 	     error[minIndex][index3to1(i,0,0)][0]);
//   }
//   fclose (fp);  
// }


void AdaptRCut::    
print_error (const std::string & file) const 
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
    fprintf (fp, "%f", (i + 0.5) * hx);
    fprintf (fp, "%e", sum);
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
	     rcut[index3to1(i,0,0)], result_error[index3to1(i,0,0)]
	);
  }
  fclose (fp);
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


