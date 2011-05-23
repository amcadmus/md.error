#include "Test.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

double global_kx = 0.1;
double global_ky = 0.2;
double global_kz = 0.3;

double global_k = 0.1;
double global_Phi = M_PI / 6.;
double global_Theta = M_PI / 6. * 5.;

double InteNaive0::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (-4032. * gsl_pow_int(r, -16)) *
	  gsl_pow_3 (r) *
	  gsl_pow_2 (sinphi) *
	  sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
	  cos(theta) )
      ;
}

double InteNaive1::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (-4032. * gsl_pow_int(r, -16)) *
	  gsl_pow_3 (r) *
	  gsl_pow_2 (sinphi) *
	  sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
	  sin(theta) )
      ;
}

double InteNaive2::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (-4032. * gsl_pow_int(r, -16)) *
	  gsl_pow_3 (r) *
	  sinphi * cos(phi) *
	  sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) 
	  )
      ;
}


double InteNaive00::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (64. * 24. * 24. * gsl_pow_int(r, -18)) *
	  gsl_pow_4 (r) *
	  gsl_pow_3 (sinphi) *
	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
	  cos(theta) * cos(theta) )
      +
      (
      	  (-8. * 24. * 24. * gsl_pow_int(r, -16)) *
      	  gsl_pow_2 (r) *
      	  sinphi *	  
      	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
      				    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) )
      ;
}

double InteNaive11::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (64. * 24. * 24. * gsl_pow_int(r, -18)) *
	  gsl_pow_4 (r) *
	  gsl_pow_3 (sinphi) *
	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
	  sin(theta) * sin(theta) )
      +
      (
      	  (-8. * 24. * 24. * gsl_pow_int(r, -16)) *
      	  gsl_pow_2 (r) *
      	  sinphi *	  
      	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
      				    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) )
      ;
}

double InteNaive22::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return
      (
	  (64. * 24. * 24. * gsl_pow_int(r, -18)) *
	  gsl_pow_4 (r) *
	  sinphi * cos(phi) * cos(phi) *
	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
	  			    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))))
      +
      (
      	  (-8. * 24. * 24. * gsl_pow_int(r, -16)) *
      	  gsl_pow_2 (r) *
      	  sinphi *	  
      	  cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
      				    sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) )
      ;
}




double InteNaive01::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return (
      (64. * 24. * 24. * gsl_pow_int(r, -18)) *
      gsl_pow_4 (r) *
      sinphi*sinphi*sinphi *
      cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
      sin(theta) * cos(theta) 
      );
}


double InteNaive02::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return (
      (64. * 24. * 24. * gsl_pow_int(r, -18)) *
      gsl_pow_4 (r) *
      sinphi*sinphi* cos(phi) *
      cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
      cos(theta)
      );
}


double InteNaive12::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);
  return (
      (64. * 24. * 24. * gsl_pow_int(r, -18)) *
      gsl_pow_4 (r) *
      sinphi*sinphi* cos(phi) *
      cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
      sin(theta)
      );
}



void
calculateGlobal ()
{
  double kx, ky, kz;
  kx = global_kx;
  ky = global_ky;
  kz = global_kz;
  global_k = sqrt(kx*kx + ky*ky + kz*kz);
	
  if (global_k != 0){
    global_Phi = acos(kz/global_k);
    if (global_Phi == 0.){
      global_Theta = 0.;
    }
    else {
      double sp = sin(global_Phi);
      double tmp = kx / (global_k * sp);
      if (tmp > 1.) tmp = 1.;
      if (tmp <-1.) tmp = -1.;
      global_Theta = acos(tmp);
      if (ky < 0) global_Theta = 2. * M_PI - global_Theta;
    }
  }
}

#include "ErrorEstimator.h"

double 
integral_an5_numerical (const double & k,
			const double & rc)
{
  f5  f_inte5;
  Integral1D<f5,  double> inte5;
  f_inte5.k = k;
  return inte5.cal_int (Integral1DInfo::Gauss4,
			f_inte5,
			rc, 50,
			1e-12);
}

double 
integral_an13_numerical (const double & k,
			 const double & rc)
{
  f13 f_inte13;
  Integral1D<f13, double> inte13;
  f_inte13.k = k;
  return inte13.cal_int (Integral1DInfo::Gauss4,
			 f_inte13,
			 rc, 30,
			 1e-16);
}

double 
integral_an15_numerical (const double & k,
			 const double & rc)
{
  f15 f_inte15;
  Integral1D<f15, double> inte15;
  f_inte15.k = k;
  return inte15.cal_int (Integral1DInfo::Gauss4,
			 f_inte15,
			 rc, 30,
			 1e-17);
}


double 
integral_f_i3_numerical (const double & k,
			 const double & rc)
{
  F_I3 f_i3;
  Integral1D<F_I3, double > inte_f_i3;
  double prec = 1e-6;
  f_i3.k = k;
  double tmpvalue, newprec;
  tmpvalue = inte_f_i3.cal_int (Integral1DInfo::Gauss4,
				f_i3,
				rc, 50,
				prec);
  newprec = fabs(tmpvalue) * 1e-3;
  printf ("newprec: %e\n", newprec);
  tmpvalue = inte_f_i3.cal_int (Integral1DInfo::Gauss4,
				f_i3,
				rc, 50,
				newprec);
  return tmpvalue;
}

double 
integral_ff_i1_numerical (const double & k,
			  const double & rc)
{
  FF_I1 ff_i1;
  Integral1D<FF_I1, double > inte_ff_i1;
  double prec = 1e-6;
  ff_i1.k = k;
  double tmpvalue, newprec;
  tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
				  ff_i1,
				  rc, 50,
				  prec);
  newprec = fabs(tmpvalue) * 1e-3;
  printf ("newprec: %e\n", newprec);
  tmpvalue =  inte_ff_i1.cal_int (Integral1DInfo::Gauss4,
				  ff_i1,
				  rc, 50,
				  newprec);  
  return tmpvalue;
}

double 
integral_ff_i5_numerical (const double & k,
			  const double & rc)
{
  FF_I5 ff_i5;
  Integral1D<FF_I5, double > inte_ff_i5;
  double prec = 1e-6;
  ff_i5.k = k;
  double tmpvalue, newprec;
  tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
			    ff_i5,
			    rc, 50,
			    prec);
  newprec = fabs(tmpvalue) * 1e-3;
  printf ("newprec: %e\n", newprec);
  tmpvalue = inte_ff_i5.cal_int (Integral1DInfo::Gauss4,
			    ff_i5,
			    rc, 50,
			    newprec);
  return tmpvalue;
}



void 
inteClever (const double & epsilon,
	    const double & sigma,
	    const double & rc,
	    double & b0k,  double & b1k,  double & b2k,
	    double & a00k, double & a01k, double & a02k,
	    double & a11k, double & a12k, double & a22k)
{
  double sigma6 = pow(sigma, 6.);
  double pre = epsilon * sigma6;
  pre = pre*pre;

  double kx, ky, kz;
  kx = global_kx;
  ky = global_ky;
  kz = global_kz;

  double i3 = integral_f_i3_numerical (global_k, rc);
  b0k = pre * 2. * M_PI / sqrt(global_k) * kx / global_k * (-i3);
  b1k = pre * 2. * M_PI / sqrt(global_k) * ky / global_k * (-i3);
  b2k = pre * 2. * M_PI / sqrt(global_k) * kz / global_k * (-i3);
  double an15 = integral_an15_numerical (global_k, rc);
  double ih = - 1./global_k * (16.) * 24. * 24. * pre * an15;
  //double ih = 0.;
  double i1 = integral_ff_i1_numerical (global_k, rc);
  double i5 = integral_ff_i5_numerical (global_k, rc);
  printf ("i1 %e, i3 %e, i5 %e, an15 %e\n", i1, i3, i5, an15);
  double tmpscale = pre * M_PI / sqrt(global_k);
  double cp = cos(global_Phi);
  double sp = sin(global_Phi);
  double ct = cos(global_Theta);
  double st = sin(global_Theta);
  a00k = tmpscale * (
      2./3. * i1 +
      (-1./3. + cp*cp - sp*sp*cos(global_Theta*2)) * i5) + ih;
  a11k = tmpscale * (
      2./3. * i1 +
      (-1./3. + cp*cp + sp*sp*cos(global_Theta*2)) * i5) + ih;
  a01k = - tmpscale * sp*sp*sin(global_Theta*2) * i5;
  a02k = - tmpscale * sin(global_Phi*2)*ct * i5;
  a12k = - tmpscale * sin(global_Phi*2)*st * i5;
  a22k = tmpscale * 4./3. *
      (0.5*i1 + (0.5 - 1.5 * cp*cp) * i5) + ih;
}


double tmpf3::
operator () (const double & r,
	     const double & theta,
	     const double & phi) const
{
  double sinphi = sin(phi);

  // // r12
  // return 1./(r*r) * sinphi*sinphi*sinphi *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     sin(theta) * cos(theta) ;
  // // r11
  // return 1./(r*r) * sinphi*sinphi*sinphi *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     cos(theta) * cos(theta) ;
  // // r22
  // return 1./(r*r) * sinphi*sinphi*sinphi *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     sin(theta) * sin(theta) ;
  // // r13
  // return 1./(r*r) * sinphi*sinphi* cos(phi) *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     cos(theta);
  // // r23
  // return 1./(r*r) * sinphi*sinphi* cos(phi) *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     sin(theta);
  // r33
  // return 1./(r*r) * sinphi* cos(phi)* cos(phi) *
  //     cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta)));
  // // r1
  // return 1./(r*r*r) * sinphi*sinphi *
  //     sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     cos(theta);
  // // r2
  // return 1./(r*r*r) * sinphi*sinphi *
  //     sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  // 				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) *
  //     sin(theta);
  // r3
  return 1./(r*r*r) * sinphi*cos(phi) *
      sin(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  				sin(phi)*sin(global_Phi) * cos(theta - global_Theta))) ;
  
}

double tmpf::
operator () (const double & r) const
{
  // // r12
  // return 
  //     1./(r*r) * (-M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi) * sin(global_Phi) *
  //     sin (2. * global_Theta) * gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r);
  // // r11
  // return
  //     1./(r*r) * ( M_PI) / sqrt(global_k*r) *
  //     (2./3. * gsl_sf_bessel_Jnu (0.5, 2 * M_PI * global_k * r) +
  //      (-1./3. + cos(global_Phi) * cos(global_Phi)) *
  //      gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r)) +
  //     1./(r*r) * (-M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi) * sin(global_Phi) *
  //     cos (2. * global_Theta) * gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r);
  // r22
  // return
  //     1./(r*r) * ( M_PI) / sqrt(global_k*r) *
  //     (2./3. * gsl_sf_bessel_Jnu (0.5, 2 * M_PI * global_k * r) +
  //      (-1./3. + cos(global_Phi) * cos(global_Phi)) *
  //      gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r)) +
  //     1./(r*r) * (M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi) * sin(global_Phi) *
  //     cos (2. * global_Theta) * gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r);
  // // r13
  // return
  //     - 1./(r*r) * (M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi * 2.) *
  //     cos(global_Theta) * gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r);
  // // r23
  // return
  //     - 1./(r*r) * (M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi * 2.) *
  //     sin(global_Theta) * gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r);
  // // r33
  // return
  //     1./(r*r) * (M_PI * 4./3.) / sqrt(global_k*r) *
  //     ( 1./2. * gsl_sf_bessel_Jnu (0.5, 2 * M_PI * global_k * r) +
  // 	(1./2. - 3./2. * cos(global_Phi) * cos(global_Phi)) *
  // 	gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r) );
  // // r1
  // return
  //     - 1./(r*r*r) * (2. * M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi) *
  //     cos(global_Theta) * gsl_sf_bessel_Jnu (1.5, 2 * M_PI * global_k * r);
  // // r2
  // return
  //     - 1./(r*r*r) * (2. * M_PI) / sqrt(global_k*r) *
  //     sin(global_Phi) *
  //     sin(global_Theta) * gsl_sf_bessel_Jnu (1.5, 2 * M_PI * global_k * r);
  // r3
  return
      - 1./(r*r*r) * (2. * M_PI) / sqrt(global_k*r) *
      cos(global_Phi) *
      gsl_sf_bessel_Jnu (1.5, 2 * M_PI * global_k * r);
}

