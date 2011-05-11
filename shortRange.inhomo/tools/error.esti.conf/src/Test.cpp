#include "Test.h"


double global_Phi = M_PI / 6.;
double global_Theta = M_PI / 6. * 5.;
double global_k = 0.1;

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
  return 1./(r*r) * sinphi* cos(phi)* cos(phi) *
      cos(-2*M_PI*global_k*r * (cos(phi)*cos(global_Phi) +
  				sin(phi)*sin(global_Phi) * cos(theta - global_Theta)));
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
  // r33
  return
      1./(r*r) * (M_PI * 4./3.) / sqrt(global_k*r) *
      ( 1./2. * gsl_sf_bessel_Jnu (0.5, 2 * M_PI * global_k * r) +
	(1./2. - 3./2. * cos(global_Phi) * cos(global_Phi)) *
	gsl_sf_bessel_Jnu (2.5, 2 * M_PI * global_k * r) );
}

