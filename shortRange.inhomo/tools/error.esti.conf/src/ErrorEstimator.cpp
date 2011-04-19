#include "ErrorEstimator.h"
#include "ForceKernel.h"

ErrorEstimatorConvN2::
ErrorEstimatorConvN2 (const ForceKernel & fk_)
    : fk (fk_)
{
}

void ErrorEstimatorConvN2::
estimate (const DensityProfile_PiecewiseConst & dp)
{
  boxsize = dp.getBox();
  nx = dp.getNx();
  ny = dp.getNy();
  nz = dp.getNz();
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  profile.resize (nx * ny * nz);
  double dvolume = hx * hy * hz;
  printf ("# error estimator: hx %f, hy %f, hz %f\n", hx, hy, hz);  
  int convCut = int(nx/2);
  if (convCut > int(ny/2)) convCut = int(ny/2);
  if (convCut > int(nz/2)) convCut = int(nz/2);
  printf ("# convCut: %d\n", convCut);  
  
  
  for (int ix = 0; ix < int(nx); ++ix){
    // printf ("ix %d, iy %d, iz %d\n", ix, iy, iz);
    printf ("ix %d\n", ix);
    for (int iy = 0; iy < int(ny); ++iy){
      for (int iz = 0; iz < int(nz); ++iz){
	// myx = (ix + 0.5) * hx;
	// myy = (iy + 0.5) * hy;
	// myz = (iz + 0.5) * hz;
	double sumx, sumy, sumz;
	sumx = sumy = sumz = 0.f;
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
	      // double f2 = fk.f2 (targetx, targety, targetz);
	      // double rho = dp.getProfile (jx, jy, jz);
	      // sum += f2 * rho * dvolume;
	      double rho = dp.getProfile (jx, jy, jz);
	      double fx, fy, fz;
	      fk.f (targetx, targety, targetz, fx, fy, fz);
	      sumx += fx * rho * dvolume;
	      sumy += fy * rho * dvolume;
	      sumz += fz * rho * dvolume;
	    }
	  }
	}
	profile[index3to1(ix, iy, iz)] = sqrt(sumx*sumx + sumy*sumy + sumz*sumz);
      }
    }
  }
}


void ErrorEstimatorConvN2::    
print_x (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (unsigned i = 0; i < nx; ++i){
    double sum = 0.;
    for (unsigned j = 0; j < ny; ++j){
      for (unsigned k = 0; k < nz; ++k){
	sum += profile[index3to1(i, j, k)];
      }
    }
    fprintf (fp, "%f %e\n", (i + 0.5) * hx, sum / ny / nz);
  }

  fclose (fp);
}

void ErrorEstimatorConvN2::    
print_xy (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (unsigned i = 0; i < nx; ++i){
    double vx = (i + 0.5) * hx;
    for (unsigned j = 0; j < ny; ++j){
      double vy = (j + 0.5) * hy;
      double sum = 0.;
      for (unsigned k = 0; k < nz; ++k){
	sum += profile[index3to1(i, j, k)];
      }
      fprintf (fp, "%f %f %e\n", vx, vy, sum / nz);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}


