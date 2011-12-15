#include "ErrorEstimator.h"
#include "ForceKernel.h"

ErrorEstimatorConvN2_Corr::
ErrorEstimatorConvN2_Corr (const ForceKernel & fk_,
		      const std::vector<double > & box,
		      const double refh)
    : fk (fk_),
      boxsize (box)
{
  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  profile.resize (nx * ny * nz, 0.);
}

void ErrorEstimatorConvN2_Corr::
estimate (const DensityProfile_Corr_PiecewiseConst & dp,
	  const int & corrBond)
{
  int myCorrBond = corrBond;
  
  double dvolume = hx * hy * hz;
  printf ("# error estimator: hx %f, hy %f, hz %f\n", hx, hy, hz);  
  double convCut = double(nx/2);
  if (convCut > double(ny/2)) convCut = double(ny/2);
  if (convCut > double(nz/2)) convCut = double(nz/2);
  convCut -= 1.00001;
  printf ("# convCut: %f\n", convCut);  
  
  for (int ix = 0; ix < int(nx); ++ix){
    printf ("ix %d\n", ix);
    for (int iy = 0; iy < int(ny); ++iy){
      for (int iz = 0; iz < int(nz); ++iz){
	// if (ix != 0 || iy != 0 || iz !=0 ) continue;
	double sumx, sumy, sumz;
	sumx = sumy = sumz = 0.;
	double sum = 0.;
	double sum2 = 0.;
	double sumCorr = 0.;
	for (int jx = 0; jx < int(nx); ++jx){
	  int dx = jx - ix;
	  if      (dx < -int(nx)/2) dx += nx;
	  else if (dx >= int(nx)/2) dx -= nx;
	  if (fabs(dx) > convCut) continue;
	  for (int jy = 0; jy < int(ny); ++jy){
	    int dy = jy - iy;
	    if      (dy < -int(ny)/2) dy += ny;
	    else if (dy >= int(ny)/2) dy -= ny;
	    if (fabs(dy) > convCut) continue;
	    for (int jz = 0; jz < int(nz); ++jz){
	      int dz = jz - iz;
	      if      (dz < -int(nz)/2) dz += nz;
	      else if (dz >= int(nz)/2) dz -= nz;
	      if (fabs(dz) > convCut) continue;
	      double targetx, targety, targetz;
	      targetx = (dx) * hx;
	      targety = (dy) * hy;
	      targetz = (dz) * hz;
	      double rho = dp.getMean ((jx+0.5)*hx, (jy+0.5)*hy, (jz+0.5)*hz);
	      // double rho = 2;
	      double f2 = fk.f2 (targetx, targety, targetz);
	      double fx, fy, fz;
	      fk.f (targetx, targety, targetz, fx, fy, fz);
	      if (fabs(fx*fx + fy*fy + fz*fz - f2) > 1e-10) {
		std::cout << "force not match!!!" << std::endl;
		printf ("%.10e %.10e %e\n", fx*fx + fy*fy + fz*fz, f2, fx*fx + fy*fy + fz*fz -f2);
	      }
	      sum += f2 * rho * dvolume;
	      sum2 += f2 * rho * rho * dvolume;
	      sumx += fx * rho * dvolume;
	      sumy += fy * rho * dvolume;
	      sumz += fz * rho * dvolume;

	      // if ((ix+0.5)*hx >= 45 && (ix+0.5)*hx <= 55) {
	      // 	myCorrBond = corrBond;
	      // }
	      // else {
	      // 	myCorrBond = 0;
	      // }
	      
	      int kx, ky, kz;
	      for (int ddx = -myCorrBond; ddx <= myCorrBond; ++ddx){
	      	kx = jx + ddx;
	      	if      (kx <   0) kx += nx;
	      	else if (kx >= nx) kx -= nx;		
	      	for (int ddy = -myCorrBond; ddy <= myCorrBond; ++ddy){
	      	  ky = jy + ddy;
	      	  if      (ky <   0) ky += ny;
	      	  else if (ky >= ny) ky -= ny;		
	      	  for (int ddz = -myCorrBond; ddz <= myCorrBond; ++ddz){
	      	    if (ddx == 0 && ddy == 0 && ddz == 0) continue;
	      	    kz = jz + ddz;
	      	    if      (kz <   0) kz += nz;
	      	    else if (kz >= nz) kz -= nz;		
	      	    double gx, gy, gz;
	      	    fk.f ((dx+ddx)*hx, (dy+ddy)*hy, (dz+ddz)*hz, gx, gy, gz);
	      	    double rho2 =
	      		dp.getCorr ((jx+0.5)*hx, (jy+0.5)*hy, (jz+0.5)*hz,
	      			    ddx*hx, ddy*hy, ddz*hz);
	      	    double rho12 =
	      		dp.getMean ((kx+0.5)*hx, (ky+0.5)*hy, (kz+0.5)*hz);
	      	    sumCorr += (rho2 - rho * rho12) * (fx*gx + fy*gy + fz*gz) * dvolume * dvolume;
	      	  }
	      	}
	      }
	    }
	  }
	}
	profile[index3to1(ix, iy, iz)] = sqrt(sum +
					      sumx*sumx + sumy*sumy + sumz*sumz +
					      sumCorr);
	printf ("%d %d %d,  sum %e, sum1 %e, sumCorr %e minus: %e  profile: %e\n",
		ix, iy, iz,
		sum, sumx*sumx + sumy*sumy + sumz*sumz, sumCorr, sum2, profile[index3to1(ix, iy, iz)]);
      }
    }
  }
}


void ErrorEstimatorConvN2_Corr::    
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
    fprintf (fp, "%f %e\n", (i + 0.5) * hx, profile[index3to1(i,0,0)]);
  }

  fclose (fp);
}

void ErrorEstimatorConvN2_Corr::    
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
      fprintf (fp, "%f %f %e\n", vx, vy, profile[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}


