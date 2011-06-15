#include "PresureCorrection.h"
#include "Integral1D.h"

PresureCorrection::
PresureCorrection (const DensityProfile_PiecewiseConst & dp)
    : malloced (false)
{
  reinit (dp);
}

PresureCorrection::
~PresureCorrection ()
{
  freeAll();
}

void PresureCorrection::
reinit (const DensityProfile_PiecewiseConst & dp)
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
  
  rhor = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rhok = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kxxk = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kyyk = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kzzk = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kxxr = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kyyr = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  kzzr = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);

  p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_MEASURE);
  p_backward_kxx = fftw_plan_dft_3d (nx, ny, nz, kxxk, kxxr, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_kyy = fftw_plan_dft_3d (nx, ny, nz, kyyk, kyyr, FFTW_BACKWARD, FFTW_MEASURE);
  p_backward_kzz = fftw_plan_dft_3d (nx, ny, nz, kzzk, kzzr, FFTW_BACKWARD, FFTW_MEASURE);
  
  pxx = pyy = pzz = 0.;

  integral_upper = 30.;
  bessel_table_h = 1e-3;
  bessel_table_hi = 1./bessel_table_h;
  double maxk = ((nx*nx) / (boxsize[0]*boxsize[0]) +
		 (ny*ny) / (boxsize[1]*boxsize[1]) +
		 (nz*nz) / (boxsize[2]*boxsize[2]) );
  maxk = sqrt(maxk);
  bessel_table_length = int (2 * M_PI * maxk * integral_upper / bessel_table_h);
  bessel_table_length += 2;
  bessel_table_1 = (double *) malloc (sizeof(double) * bessel_table_length);
  bessel_table_5 = (double *) malloc (sizeof(double) * bessel_table_length);
  bessel_table_1[0] = 0.;
  bessel_table_5[0] = 0.;
  for (int i = 1; i < bessel_table_length; ++i){
    bessel_table_1[i] = bessel_table_1[i-1] + bessel_table_h;
    bessel_table_5[i] = bessel_table_5[i-1] + bessel_table_h;
  }
  // bessel_table_1[0] = 1e-8;
  // bessel_table_5[0] = 1e-8;
  gsl_sf_bessel_sequence_Jnu_e (0.5,
  				GSL_PREC_DOUBLE,
  				bessel_table_length,
  				bessel_table_1);
  gsl_sf_bessel_sequence_Jnu_e (2.5,
  				GSL_PREC_DOUBLE,
  				bessel_table_length,
  				bessel_table_5);

  // double percetage = 0.;
  // double step = 0.1;
  // for (int i = 0; i < bessel_table_length; ++i){
  //   if ((double (i+1) * 100. / double(bessel_table_length)) >= percetage){
  //     percetage += step;
  //     printf ("# cal bessel table %.1f %%   \r", percetage);
  //     fflush (stdout);
  //   }
  //   double x = bessel_table_h * i;
  //   bessel_table_1[i] = gsl_sf_bessel_Jnu (0.5, x);
  //   bessel_table_5[i] = gsl_sf_bessel_Jnu (2.5, x);
  // }
  // printf ("\n");

  ff_i1.bessel_table = bessel_table_1;
  ff_i1.bessel_table_length = bessel_table_length;
  ff_i1.bessel_table_h  = bessel_table_h;
  ff_i1.bessel_table_hi = bessel_table_hi;
  ff_i5.bessel_table = bessel_table_5;
  ff_i5.bessel_table_length = bessel_table_length;
  ff_i5.bessel_table_h  = bessel_table_h;
  ff_i5.bessel_table_hi = bessel_table_hi;  

  malloced = true;
}

void PresureCorrection::
freeAll ()
{
  if (malloced){
    free (rhor);
    free (rhok);
    free (kxxk);
    free (kxxr);
    free (kyyk);
    free (kyyr);
    free (kzzk);
    free (kzzr);
    fftw_destroy_plan (p_forward_rho);
    fftw_destroy_plan (p_backward_kxx);
    fftw_destroy_plan (p_backward_kyy);
    fftw_destroy_plan (p_backward_kzz);
    free (bessel_table_1);
    free (bessel_table_5);
    malloced = false;
  }
}

void PresureCorrection::
naiveCorrection (const DensityProfile_PiecewiseConst & dp)
{
  pxx = pyy = pzz = 0.;

  double rcut1 = 5.;
  double rcut2 = 10.;
  int dnx = (rcut2 - 1e-8) / hx + 1;
  int dny = (rcut2 - 1e-8) / hy + 1;
  int dnz = (rcut2 - 1e-8) / hz + 1;
  
  double dvolume2 = hx * hy * hz;
  dvolume2 = dvolume2 * dvolume2;

  for (int i = 0; i < nele; ++i){
    if (i % 100 == 0){
      printf ("# i is %d   \r", i);
      fflush (stdout);
    }
    int ix1, iy1, iz1;
    index1to3 (i, ix1, iy1, iz1);
    double rho1 = dp.getProfile()[i];
    for (int dix = -dnx; dix <= dnx; ++ dix){
      double dx = dix * hx;
      for (int diy = -dny; diy <= dny; ++ diy){
	double dy = diy * hy;
	for (int diz = -dnz; diz <= dnz; ++ diz){
	  double dz = diz * hz;
	  double diff2 = dx*dx + dy*dy + dz*dz;
	  if (diff2 < rcut1*rcut1 || diff2 > rcut2*rcut2) continue;
	  int ix2 = ix1 + dix;
	  if      (ix2 <  0 ) ix2 += nx;
	  else if (ix2 >= nx) ix2 -= nx;
	  int iy2 = iy1 + diy;
	  if      (iy2 <  0 ) iy2 += ny;
	  else if (iy2 >= ny) iy2 -= ny;
	  int iz2 = iz1 + diz;
	  if      (iz2 <  0 ) iz2 += nz;
	  else if (iz2 >= nz) iz2 -= nz;
	  double rho2 = dp.getProfile (ix2, iy2, iz2);

	  double scalor_Kernel = -24 / gsl_pow_4 (diff2);
	  pxx += scalor_Kernel * dx * dx * rho1 * rho2 * dvolume2;
	  pyy += scalor_Kernel * dy * dy * rho1 * rho2 * dvolume2;
	  pzz += scalor_Kernel * dz * dz * rho1 * rho2 * dvolume2;
	}
      }
    }
  }

  printf ("\n");
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  pxx *= 0.5 / (volume);
  pyy *= 0.5 / (volume);
  pzz *= 0.5 / (volume);
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

void PresureCorrection::
correction (const DensityProfile_PiecewiseConst & dp)
{
  for (int i = 0; i < nele; ++i){
    rhor[i][0] = dp.getProfile()[i];
    rhor[i][1] = 0.;
  }
  fftw_execute (p_forward_rho);
  double volume = boxsize[0] * boxsize[1] * boxsize[2];
  double scale = volume / nele;
  for (int i = 0; i < nele; ++i){
    rhok[i][0] *= scale;
    rhok[i][1] *= scale;
  }  

  makeKernel ();

  array_multiply (kxxk, nele, rhok);
  array_multiply (kyyk, nele, rhok);
  array_multiply (kzzk, nele, rhok);

  fftw_execute (p_backward_kxx);
  fftw_execute (p_backward_kyy);
  fftw_execute (p_backward_kzz);
  
  double dvolume = hx * hy * hz;
  for (int i = 0; i < nele; ++i){
    pxx += 0.5 * kxxr[i][0] * rhor[i][0] / volume * dvolume;
    pyy += 0.5 * kyyr[i][0] * rhor[i][0] / volume * dvolume;
    pzz += 0.5 * kzzr[i][0] * rhor[i][0] / volume * dvolume;
  }
  pxx /= (volume);
  pyy /= (volume);
  pzz /= (volume);
}



double PresureCorrection::
integral_ff_i1_numerical (const double & k,
			  const double & rc)
{
  // double tmpvalue = 0.;
  // int ndiv = 100;
  // double h = (integral_upper - rc) / double(ndiv);
  // double x0 = rc;
  // double x1 = rc + h;
  // double v0 = ff_i1 (x0);
  // double v1;
  // for (int i = 0; i < ndiv; ++i) {
  //   v1 = ff_i1(x1);
  //   tmpvalue += 0.5 * h * (v0 + v1);
  //   v0 = v1;
  //   x0 = x1;
  //   x1 += h;
  // }
  // FF_I1 ff_i1;
  Integral1D<FF_I1, double > inte_ff_i1;
  double prec = 1e-7;
  ff_i1.k = k;
  double tmpvalue, newprec;
  tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
  				  ff_i1,
  				  rc, integral_upper,
  				  prec);
  // printf ("value: %e\n", tmpvalue);
  // newprec = fabs(tmpvalue) * 1e-4;
  // // printf ("newprec: %e\n", newprec);
  // tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
  // 				  ff_i1,
  // 				  rc, 30,
  // 				  newprec);  
  return tmpvalue;
}

double PresureCorrection::
integral_ff_i5_numerical (const double & k,
			  const double & rc)
{
  // double tmpvalue = 0.;
  // int ndiv = 100;
  // double h = (integral_upper - rc) / double(ndiv);
  // double x0 = rc;
  // double x1 = rc + h;
  // double v0 = ff_i5 (x0);
  // double v1;
  // for (int i = 0; i < ndiv; ++i) {
  //   v1 = ff_i5(x1);
  //   tmpvalue += 0.5 * h * (v0 + v1);
  //   v0 = v1;
  //   x0 = x1;
  //   x1 += h;
  // }
  // FF_I5 ff_i5;
  Integral1D<FF_I5, double > inte_ff_i5;
  double prec = 1e-7;
  ff_i5.k = k;
  double tmpvalue, newprec;
  tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
  				 ff_i5,
  				 rc, integral_upper,
  				 prec);
  // newprec = fabs(tmpvalue) * 1e-4;
  // // printf ("newprec: %e\n", newprec);
  // tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
  // 			    ff_i5,
  // 			    rc, 30,
  // 			    newprec);
  return tmpvalue;
}



void PresureCorrection::
makeKernel ()
{
  double pre = 1.;
  double rc = 7.5;

  double kx, ky, kz;
  for (int ix = -nx/2; ix < nx - nx/2; ++ix){
    kx = ix / boxsize[0];
    int locx;
    if (ix < 0) locx = ix + nx;
    else locx = ix;
    printf ("# ix is %d   \r", ix);
    fflush (stdout);
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

	  double i1 = integral_ff_i1_numerical (k, rc);
	  double i5 = integral_ff_i5_numerical (k, rc);
	  double tmpscale = pre * M_PI / sqrt(k);
	  double cp = cos(phi_k);
	  double sp = sin(phi_k);
	  // double ct = cos(theta_k);
	  // double st = sin(theta_k);
	  kxxk[posi][0] = tmpscale * (
	      2./3. * i1 +
	      (-1./3. + cp*cp - sp*sp*cos(theta_k*2)) * i5);
	  kxxk[posi][1] = 0;
	  kyyk[posi][0] = tmpscale * (
	      2./3. * i1 +
	      (-1./3. + cp*cp + sp*sp*cos(theta_k*2)) * i5);
	  kyyk[posi][1] = 0;
	  kzzk[posi][0] = tmpscale * 4./3. *
	      (0.5*i1 + (0.5 - 1.5 * cp*cp) * i5);
	  kzzk[posi][1] = 0;
	}
	else {
	  kxxk[posi][0] = 0;
	  kxxk[posi][1] = 0;
	  kyyk[posi][0] = 0;
	  kyyk[posi][1] = 0;
	  kzzk[posi][0] = 0;
	  kzzk[posi][1] = 0;
	}
      }
    }
  }

  printf ("\n");
}



// void correction (const DensityProfile_PiecewiseConst & dp)
// {
  
	  
