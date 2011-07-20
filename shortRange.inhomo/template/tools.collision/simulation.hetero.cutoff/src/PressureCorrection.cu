#include "PressureCorrection.h"
#include "Integral1D.h"

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


PressureCorrection::
PressureCorrection (const AdaptRCut & arc,
			     const DensityProfile_PiecewiseConst & dp)
    : malloced (false)
{
  reinit (arc, dp);
}

PressureCorrection::
~PressureCorrection ()
{
  freeAll();
}

void PressureCorrection::
reinit (const AdaptRCut & arc,
	const DensityProfile_PiecewiseConst & dp)
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

  // nrc = arc.getNrc();
  // rcList.resize (nrc);
  // for (int i = 0; i < nrc; ++i){
  //   rcList
  rcList = arc.getRcList ();
  nrc = rcList.size();
  
  rhor = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  rhok = (fftw_complex *) malloc (sizeof(fftw_complex) * nele);
  mallocArrayComplex (&kxxk, nrc, nele);
  mallocArrayComplex (&kyyk, nrc, nele);
  mallocArrayComplex (&kzzk, nrc, nele);
  mallocArrayComplex (&convxx, nrc, nele);
  mallocArrayComplex (&convyy, nrc, nele);
  mallocArrayComplex (&convzz, nrc, nele);
  
  p_forward_rho = fftw_plan_dft_3d (nx, ny, nz, rhor, rhok, FFTW_FORWARD,  FFTW_PATIENT);
  p_backward_kxx = (fftw_plan *) malloc (sizeof(fftw_plan) * nrc);
  p_backward_kyy = (fftw_plan *) malloc (sizeof(fftw_plan) * nrc);
  p_backward_kzz = (fftw_plan *) malloc (sizeof(fftw_plan) * nrc);
  for (int i = 0; i < nrc; ++i){
    p_backward_kxx[i] = fftw_plan_dft_3d (nx, ny, nz, convxx[i], convxx[i], FFTW_BACKWARD, FFTW_PATIENT);
    p_backward_kyy[i] = fftw_plan_dft_3d (nx, ny, nz, convyy[i], convyy[i], FFTW_BACKWARD, FFTW_PATIENT);
    p_backward_kzz[i] = fftw_plan_dft_3d (nx, ny, nz, convzz[i], convzz[i], FFTW_BACKWARD, FFTW_PATIENT);
  }
  
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
  				GSL_PREC_APPROX,
  				bessel_table_length,
  				bessel_table_1);
  gsl_sf_bessel_sequence_Jnu_e (2.5,
  				GSL_PREC_APPROX,
  				bessel_table_length,
  				bessel_table_5);

  ff_i1.bessel_table = bessel_table_1;
  ff_i1.bessel_table_length = bessel_table_length;
  ff_i1.bessel_table_h  = bessel_table_h;
  ff_i1.bessel_table_hi = bessel_table_hi;
  ff_i5.bessel_table = bessel_table_5;
  ff_i5.bessel_table_length = bessel_table_length;
  ff_i5.bessel_table_h  = bessel_table_h;
  ff_i5.bessel_table_hi = bessel_table_hi;  

  malloced = true;

  printf ("# pressure correction making kernel ...\n");
  fflush (stdout);
  makeKernel ();
  printf ("# done\n");
  fflush (stdout);
}

void PressureCorrection::
freeAll ()
{
  if (malloced){
    free (rhor);
    free (rhok);
    freeArrayComplex (&kxxk, nrc);
    freeArrayComplex (&kyyk, nrc);
    freeArrayComplex (&kzzk, nrc);
    freeArrayComplex (&convxx, nrc);
    freeArrayComplex (&convyy, nrc);
    freeArrayComplex (&convzz, nrc);
    fftw_destroy_plan (p_forward_rho);
    for (int i = 0; i < nrc; ++i){
      fftw_destroy_plan (p_backward_kxx[i]);
      fftw_destroy_plan (p_backward_kyy[i]);
      fftw_destroy_plan (p_backward_kzz[i]);
    }
    free (p_backward_kxx);
    free (p_backward_kyy);
    free (p_backward_kzz);
    free (bessel_table_1);
    free (bessel_table_5);
    malloced = false;
  }
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

void PressureCorrection::
correction (const AdaptRCut & arc,
	    const DensityProfile_PiecewiseConst & dp)
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

  array_multiply (convxx, nrc, nele, kxxk, rhok);
  array_multiply (convyy, nrc, nele, kyyk, rhok);
  array_multiply (convzz, nrc, nele, kzzk, rhok);

  for (int count = 0; count < nrc; ++count){
    fftw_execute (p_backward_kxx[count]);
    fftw_execute (p_backward_kyy[count]);
    fftw_execute (p_backward_kzz[count]);
  }
  
  double dvolume = hx * hy * hz;
  for (int i = 0; i < nele; ++i){
    int index = arc.getProfileIndex()[i];
    pxx += 0.5 * convxx[index][i][0] * rhor[i][0] / volume * dvolume;
    pyy += 0.5 * convyy[index][i][0] * rhor[i][0] / volume * dvolume;
    pzz += 0.5 * convzz[index][i][0] * rhor[i][0] / volume * dvolume;
  }
  pxx /= (volume);
  pyy /= (volume);
  pzz /= (volume);
}



double PressureCorrection::
integral_ff_i1_numerical (const double & k,
			  const double & rc,
			  const double & rc2)
{
  Integral1D<FF_I1, double > inte_ff_i1;
  double prec = 1e-7;
  ff_i1.k = k;
  double tmpvalue;
  tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
  				  ff_i1,
  				  rc, rc2,
  				  prec);
  // printf ("value: %e\n", tmpvalue);
  // double newprec = fabs(tmpvalue) * 1e-4;
  // // printf ("newprec: %e\n", newprec);
  // tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
  // 				  ff_i1,
  // 				  rc, 30,
  // 				  newprec);  
  return tmpvalue;
}

double PressureCorrection::
integral_ff_i5_numerical (const double & k,
			  const double & rc,
			  const double & rc2)
{
  Integral1D<FF_I5, double > inte_ff_i5;
  double prec = 1e-7;
  ff_i5.k = k;
  double tmpvalue;
  tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
  				 ff_i5,
  				 rc, rc2,
  				 prec);
  // double newprec = fabs(tmpvalue) * 1e-4;
  // // printf ("newprec: %e\n", newprec);
  // tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
  // 			    ff_i5,
  // 			    rc, 30,
  // 			    newprec);
  return tmpvalue;
}



void PressureCorrection::
makeKernel ()
{
  double pre = 1.;

  double kx, ky, kz;
  for (int count = nrc-1; count >= 0; --count){
    printf ("# rc is %f   \r", rcList[count]);
    for (int ix = -nx/2; ix < nx - nx/2; ++ix){
      kx = ix / boxsize[0];
      int locx;
      if (ix < 0) locx = ix + nx;
      else locx = ix;
      // printf ("# rc is %f, ix is %d   \r", rcList[count], ix);
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
	    double tmpscale = pre * M_PI / sqrt(k);
	    double cp = cos(phi_k);
	    double sp = sin(phi_k);
	    // double ct = cos(theta_k);
	    // double st = sin(theta_k);
	    if (count == nrc - 1){
	      double rc = rcList[count];
	      double i1 = integral_ff_i1_numerical (k, rc, integral_upper);
	      double i5 = integral_ff_i5_numerical (k, rc, integral_upper);
	      kxxk[count][posi][0] = tmpscale * (
		  2./3. * i1 +
		  (-1./3. + cp*cp - sp*sp*cos(theta_k*2)) * i5);
	      kyyk[count][posi][0] = tmpscale * (
		  2./3. * i1 +
		  (-1./3. + cp*cp + sp*sp*cos(theta_k*2)) * i5);
	      kzzk[count][posi][0] = tmpscale * 4./3. *
		  (0.5*i1 + (0.5 - 1.5 * cp*cp) * i5);
	    }
	    else {
	      double i1 = integral_ff_i1_numerical (k, rcList[count], rcList[count+1]);
	      double i5 = integral_ff_i5_numerical (k, rcList[count], rcList[count+1]);
	      kxxk[count][posi][0] = kxxk[count+1][posi][0] + 
		  tmpscale * (
		      2./3. * i1 +
		      (-1./3. + cp*cp - sp*sp*cos(theta_k*2)) * i5);
	      kyyk[count][posi][0] = kyyk[count+1][posi][0] +
		  tmpscale * (
		      2./3. * i1 +
		      (-1./3. + cp*cp + sp*sp*cos(theta_k*2)) * i5);
	      kzzk[count][posi][0] = kzzk[count+1][posi][0] +
		  tmpscale * 4./3. *
		  (0.5 * i1 + (0.5 - 1.5 * cp * cp) * i5);
	    }
	    kxxk[count][posi][1] = 0;
	    kyyk[count][posi][1] = 0;
	    kzzk[count][posi][1] = 0;	      
	  }
	  else {
	    kxxk[count][posi][0] = 0;
	    kxxk[count][posi][1] = 0;
	    kyyk[count][posi][0] = 0;
	    kyyk[count][posi][1] = 0;
	    kzzk[count][posi][0] = 0;
	    kzzk[count][posi][1] = 0;
	  }
	}
      }
    }
  }
  
  printf ("\n");
}


	  
