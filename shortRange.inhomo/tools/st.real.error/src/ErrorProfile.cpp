#include "ErrorProfile.h"
#include "GroFileManager.h"

ErrorProfile_PiecewiseConst::
ErrorProfile_PiecewiseConst (const std::string & filename,
			       const double & refh)
{
  std::vector<std::vector<double > > posi;
  std::vector<std::vector<double > > force_error;

  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    exit(1);
  }
  int natom;
  boxsize.resize(3);
  fscanf (fp, "%d%lf%lf%lf", &natom, &boxsize[0], &boxsize[1], &boxsize[2]);
  posi.resize (natom);
  force_error.resize (natom);

  for (int i = 0; i < natom; ++i){
    double xx, yy, zz;
    double fx, fy, fz;
    fscanf (fp, "%lf%lf%lf%lf%lf%lf", &xx, &yy, &zz, &fx, &fy, &fz);
    posi[i].push_back (xx);
    posi[i].push_back (yy);
    posi[i].push_back (zz);
    force_error[i].push_back (fx);
    force_error[i].push_back (fy);
    force_error[i].push_back (fz);
  }
  fclose (fp);
  
  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;

  profile.clear();
  profile.resize (nx * ny * nz, 0.);
  std::vector<int > count;
  count.resize (nx * ny * nz, 0);

  for (unsigned i = 0; i < posi.size(); ++i) {
    double tmp;
    tmp = posi[i][0];
    if      (posi[i][0] >= boxsize[0]) tmp -= boxsize[0];
    else if (posi[i][0] <  0)          tmp += boxsize[0];
    unsigned ix = unsigned (tmp / hx);
    tmp = posi[i][1];
    if      (posi[i][1] >= boxsize[1]) tmp -= boxsize[1];
    else if (posi[i][1] <  0)          tmp += boxsize[1];
    unsigned iy = unsigned (tmp / hy);
    tmp = posi[i][2];
    if      (posi[i][2] >= boxsize[2]) tmp -= boxsize[2];
    else if (posi[i][2] <  0)          tmp += boxsize[2];
    unsigned iz = unsigned (tmp / hz);
    profile[index3to1(ix, iy, iz)] +=
	( force_error[i][0] * force_error[i][0] +
	  force_error[i][1] * force_error[i][1] +
	  force_error[i][2] * force_error[i][2] );
    count[index3to1(ix, iy, iz)] += 1;
  }

  for (unsigned i = 0; i < profile.size(); ++i){
    if (count[i] == 0) continue;
    profile[i] = sqrt (profile[i] / count[i]);
  }
}

void ErrorProfile_PiecewiseConst::    
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
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, profile[index3to1(i,0,0)]);
    // fprintf (fp, "%f %f\n", (i + 0.5) * hx, sum / ny / nz);
  }

  fclose (fp);
}

void ErrorProfile_PiecewiseConst::    
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
      // fprintf (fp, "%f %f %f\n", vx, vy, sum / nz);
      fprintf (fp, "%f %f %f\n", vx, vy, profile[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}

