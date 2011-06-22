#include "ErrorProfile.h"
#include "GroFileManager.h"
#include <stdlib.h>
#include <stdio.h>

ErrorProfile_PiecewiseConst::
ErrorProfile_PiecewiseConst (const std::vector<double > &box,
			     const double & refh)
{
  boxsize = box;

  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;

  profile.clear();
  profile.resize (nx * ny * nz, 0.);
  count.clear();
  count.resize (nx * ny * nz, 0);
}


void ErrorProfile_PiecewiseConst::
deposit (const std::vector<std::vector<double > > & coord,
	 const std::vector<std::vector<double > > & force)
{
  for (unsigned i = 0; i < coord.size(); ++i) {
    double tmp;
    tmp = coord[i][0];
    if      (coord[i][0] >= boxsize[0]) tmp -= boxsize[0];
    else if (coord[i][0] <  0)          tmp += boxsize[0];
    unsigned ix = unsigned (tmp / hx);
    tmp = coord[i][1];
    if      (coord[i][1] >= boxsize[1]) tmp -= boxsize[1];
    else if (coord[i][1] <  0)          tmp += boxsize[1];
    unsigned iy = unsigned (tmp / hy);
    tmp = coord[i][2];
    if      (coord[i][2] >= boxsize[2]) tmp -= boxsize[2];
    else if (coord[i][2] <  0)          tmp += boxsize[2];
    unsigned iz = unsigned (tmp / hz);
    profile[index3to1(ix, iy, iz)] +=
	( force[i][0] * force[i][0] +
	  force[i][1] * force[i][1] +
	  force[i][2] * force[i][2] );
    count[index3to1(ix, iy, iz)] += 1;
  }
}

void ErrorProfile_PiecewiseConst::
calculate ()
{
  for (unsigned i = 0; i < profile.size(); ++i){
    if (count[i] == 0) continue;
    profile[i] = sqrt (profile[i] / count[i]);
  }
}

void ErrorProfile_PiecewiseConst::
clear ()
{
  profile.clear ();
  count.clear ();
  profile.resize (nx*ny*nz, 0.);
  count.resize (nx*ny*nz, 0);
}  
  

void ErrorProfile_PiecewiseConst::    
print_x (const char * file) const 
{
  FILE * fp = fopen (file, "w");
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
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, profile[index3to1(i,0,0)]);
    // fprintf (fp, "%f %f\n", (i + 0.5) * hx, sum / ny / nz);
  }

  fclose (fp);
}

void ErrorProfile_PiecewiseConst::    
print_x_avg (const char * file) const 
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (unsigned i = 0; i < nx; ++i){
    double sum = 0.;
    int sumcount = 0;
    for (unsigned j = 0; j < ny; ++j){
      for (unsigned k = 0; k < nz; ++k){
    	sum += profile[index3to1(i, j, k)] *
	    profile[index3to1(i, j, k)] *
	    count [index3to1(i, j, k)];
	sumcount += count [index3to1(i, j, k)];
      }
    }
    if (sumcount != 0) sum /= double(sumcount);
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, sqrt(sum));
  }

  fclose (fp);
}


void ErrorProfile_PiecewiseConst::    
print_xy (const char * file) const 
{
  FILE * fp = fopen (file, "w");
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
      // fprintf (fp, "%f %f %f\n", vx, vy, sum / nz);
      fprintf (fp, "%f %f %f\n", vx, vy, profile[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}



