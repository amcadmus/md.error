#include "DensityProfile.h"
#include "GroFileManager.h"

DensityProfile_PiecewiseConst::
DensityProfile_PiecewiseConst (const std::string & filename,
			       const double & refh)
{
  std::vector<int > resdindex;
  std::vector<std::string >  resdname;
  std::vector<std::string >  atomname;
  std::vector<int > atomindex;
  std::vector<std::vector<double > > posi;
  std::vector<std::vector<double > > velo;

  GroFileManager::read (filename,
			resdindex, resdname, atomname, atomindex,
			posi, velo, boxsize);
  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;

  profile.clear();
  profile.resize (nx * ny * nz, 0.);

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
    profile[index3to1(ix, iy, iz)] += 1.;
  }

  double dvolume = hx * hy * hz;

  for (unsigned i = 0; i < profile.size(); ++i){
    profile[i] /= dvolume;
  }
}

void DensityProfile_PiecewiseConst::    
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
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, sum / ny / nz);
  }

  fclose (fp);
}

void DensityProfile_PiecewiseConst::    
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
      fprintf (fp, "%f %f %f\n", vx, vy, sum / nz);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}

