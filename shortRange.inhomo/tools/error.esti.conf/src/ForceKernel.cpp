#include "ForceKernel.h"

Disperson6::
Disperson6 (const double & epsilon_,
	    const double & sigma_,
	    const double & rc_)
    : epsilon (epsilon_),
      sigma (sigma_),
      sigma6 (sigma_ * sigma_ * sigma_ * sigma_ * sigma_ * sigma_),
      rc (rc_),
      rc2 (rc_ * rc_)
{
}

Disperson6_Taylor::
Disperson6_Taylor (const double & epsilon_,
		   const double & sigma_,
		   const double & rc_)
    : Disperson6 (epsilon_, sigma_, rc_)
{
}

double Disperson6_Taylor::
u0p (const double & r) const
{
  double r2 = r*r;
  double r4 = r2*r2;
  return - 4. * epsilon * sigma6 / (r4 * r2);
}

double Disperson6_Taylor::
u1p (const double & r) const
{
  double r2 = r*r;
  double r4 = r2*r2;
  return 24. * epsilon * sigma6 / (r4 * r2 * r);
}

double Disperson6_Taylor::
u2p (const double & r) const
{
  double r2 = r*r;
  double r4 = r2*r2;
  return -168. * epsilon * sigma6 / (r4 * r4);
}

double Disperson6_Taylor::
u3p (const double & r) const
{
  double r2 = r*r;
  double r4 = r2*r2;
  return 1344. * epsilon * sigma6 / (r4 * r4 * r);
}


void Disperson6_Taylor::
f_Taylor (const double & refx,
	  const double & refy,
	  const double & refz,
	  const double & dx,
	  const double & dy,
	  const double & dz,
	  double & fx,
	  double & fy,
	  double & fz) const
{
  double refr[3];
  refr[0] = refx;
  refr[1] = refy;
  refr[2] = refz;
  double dr[3];
  dr[0] = dx;
  dr[1] = dy;
  dr[2] = dz;

  double r1 = sqrt(refx*refx + refy*refy + refz*refz);
  double r2 = r1 * r1;
  double r3 = r2 * r1;
  double r4 = r2 * r2;
  double r5 = r3 * r2;
  double u1pr = u1p(r1);
  double u2pr = u2p(r1);
  double u3pr = u3p(r1);

  double tmpa = u2pr / r2 - u1pr / r3;
  double tmpb = u1pr / r1;
  double tmpc = u3pr / r3 - 3. * u2pr / r4 + 3. * u1pr / r5 ;
  {
    for (int ii = 0; ii < 3; ++ii){
      for (int jj = 0; jj < 3; ++jj){
	gradF[ii][jj] = tmpa * refr[ii] * refr[jj];
      }
    }
    for (int ii = 0; ii < 3; ++ii){
      gradF[ii][ii] += tmpb;
    }
  }

  {
    for (int ii = 0; ii < 3; ++ii){
      for (int jj = 0; jj < 3; ++jj){
	for (int kk = 0; kk < 3; ++kk){
	  hessianF[ii][jj][kk] = tmpc * refr[ii] * refr[jj] * refr[kk];
	}
      }
    }
    for (int ii = 0; ii < 3; ++ii){
      for (int jj = 0; jj < 3; ++jj){
	hessianF[ii][ii][jj] += tmpa * refr[jj];
	hessianF[ii][jj][ii] += tmpa * refr[jj];
	hessianF[jj][ii][ii] += tmpa * refr[jj];
      }
    }
  }

  double df[3] = {0.};
  
  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      df[ii] += gradF[ii][jj] * dr[jj];
    }
  }
  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      for (int kk = 0; kk < 3; ++kk){
	df[ii] += .5 * hessianF[ii][jj][kk] * dr[jj] * dr[kk];
      }
    }
  }

  f (refx, refy, refz, fx, fy, fz);
  fx += df[0];
  fy += df[1];
  fz += df[2];
}

double Disperson6_Taylor::
fdotf_taylor (const double & refx,
	      const double & refy,
	      const double & refz,
	      const double & dx,
	      const double & dy,
	      const double & dz) const
{
  double refr[3];
  refr[0] = refx;
  refr[1] = refy;
  refr[2] = refz;
  double dr[3];
  dr[0] = dx;
  dr[1] = dy;
  dr[2] = dz;

  double r1 = sqrt(refx*refx + refy*refy + refz*refz);
  double r2 = r1 * r1;
  double r3 = r2 * r1;
  double r4 = r2 * r2;
  double u1pr = u1p(r1);
  double u2pr = u2p(r1);
  double u3pr = u3p(r1);
  
  double tmpa = u1pr * (u3pr / r2 - u2pr / r3 + u1pr / r4);
  double tmpb = u1pr * (u2pr / r1 - u1pr / r2);
  double tmpc = u2pr * u2pr / r2 - u1pr * u1pr / r4;
  double tmpd = u1pr * u1pr / r2;

  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      matA[ii][jj] = (tmpa - tmpc) * refr[ii] * refr[jj] ;
    }
  }
  for (int ii = 0; ii < 3; ++ii){
    matA[ii][ii] += (tmpb - tmpd) ;
  }

  double sum = f2 (refx, refy, refz);
  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      sum += dr[ii] * matA[ii][jj] * dr[jj];
    }
  }

  return sum;
}


double Disperson6_Taylor::
fdotf_oneSide_taylor (const double & refx,
		      const double & refy,
		      const double & refz,
		      const double & dx,
		      const double & dy,
		      const double & dz) const
{
  double refr[3];
  refr[0] = refx;
  refr[1] = refy;
  refr[2] = refz;
  double dr[3];
  dr[0] = dx;
  dr[1] = dy;
  dr[2] = dz;
  
  double r1 = sqrt(refx*refx + refy*refy + refz*refz);
  double r2 = r1 * r1;
  double r3 = r2 * r1;
  double r4 = r2 * r2;
  double u1pr = u1p(r1);
  double u2pr = u2p(r1);
  double u3pr = u3p(r1);
  
  double tmpa = u1pr * (u3pr / r2 - u2pr / r3 + u1pr / r4);
  double tmpb = u1pr * (u2pr / r1 - u1pr / r2);
  // double tmpc = u2pr * u2pr / r2 - u1pr * u1pr / r4;
  // double tmpd = u1pr * u1pr / r2;

  for (int ii = 0; ii < 3; ++ii){
    vecB[ii] = u1pr * u2pr / r1 * refr[ii];
  }
  
  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      matA[ii][jj] = (tmpa) * refr[ii] * refr[jj] ;
    }
  }
  for (int ii = 0; ii < 3; ++ii){
    matA[ii][ii] += (tmpb) ;
  }

  double sum = f2 (refx, refy, refz);
  for (int ii = 0; ii < 3; ++ii){
    sum += vecB[ii] * dr[ii];
  }
  for (int ii = 0; ii < 3; ++ii){
    for (int jj = 0; jj < 3; ++jj){
      sum += 0.5 * dr[ii] * matA[ii][jj] * dr[jj];
    }
  }

  return sum;
}


