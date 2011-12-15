#include "ErrorEstimator.h"
#include "ForceKernel.h"

ErrorEstimatorConvN2::
ErrorEstimatorConvN2 (const ForceKernel & fk_,
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

void ErrorEstimatorConvN2::
estimate (const DensityProfile_PiecewiseConst & dp)
{
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
	      double rho = dp.getValue ((jx+0.5)*hx, (jy+0.5)*hy, (jz+0.5)*hz);
	      // double rho = 2;
	      double f2 = fk.f2 (targetx, targety, targetz);
	      double fx, fy, fz;
	      fk.f (targetx, targety, targetz, fx, fy, fz);
	      if (fabs(fx*fx + fy*fy + fz*fz - f2) > 1e-10) {
		std::cout << "force not match!!!" << std::endl;
		printf ("%.10e %.10e %e\n", fx*fx + fy*fy + fz*fz, f2, fx*fx + fy*fy + fz*fz -f2);
	      }
	      sum += f2 * rho * dvolume;
	      sumx += fx * rho * dvolume;
	      sumy += fy * rho * dvolume;
	      sumz += fz * rho * dvolume;
	    }
	  }
	}
	profile[index3to1(ix, iy, iz)] = sqrt(sum + sumx*sumx + sumy*sumy + sumz*sumz);
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
    // double sum = 0.;
    // for (unsigned j = 0; j < ny; ++j){
    //   for (unsigned k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    fprintf (fp, "%f %e\n", (i + 0.5) * hx, profile[index3to1(i,0,0)]);
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
      // double sum = 0.;
      // for (unsigned k = 0; k < nz; ++k){
      // 	sum += profile[index3to1(i, j, k)];
      // }
      fprintf (fp, "%f %f %e\n", vx, vy, profile[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}


