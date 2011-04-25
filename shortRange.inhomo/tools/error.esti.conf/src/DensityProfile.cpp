#include "DensityProfile.h"
#include "GroFileManager.h"
#include <string.h>

DensityProfile_PiecewiseConst::
DensityProfile_PiecewiseConst (const std::string & filename,
			       const double & refh)
{
  reinit_xtc (filename, refh);
  
  // std::vector<int > resdindex;
  // std::vector<std::string >  resdname;
  // std::vector<std::string >  atomname;
  // std::vector<int > atomindex;
  // std::vector<std::vector<double > > posi;
  // std::vector<std::vector<double > > velo;

  // GroFileManager::read (filename,
  // 			resdindex, resdname, atomname, atomindex,
  // 			posi, velo, boxsize);
  // nx = unsigned (boxsize[0] / refh);
  // ny = unsigned (boxsize[1] / refh);
  // nz = unsigned (boxsize[2] / refh);
  // hx = boxsize[0] / nx;
  // hy = boxsize[1] / ny;
  // hz = boxsize[2] / nz;

  // profile.clear();
  // profile.resize (nx * ny * nz, 0.);

  // for (unsigned i = 0; i < posi.size(); ++i) {
  //   double tmp;
  //   tmp = posi[i][0];
  //   if      (posi[i][0] >= boxsize[0]) tmp -= boxsize[0];
  //   else if (posi[i][0] <  0)          tmp += boxsize[0];
  //   unsigned ix = unsigned (tmp / hx);
  //   tmp = posi[i][1];
  //   if      (posi[i][1] >= boxsize[1]) tmp -= boxsize[1];
  //   else if (posi[i][1] <  0)          tmp += boxsize[1];
  //   unsigned iy = unsigned (tmp / hy);
  //   tmp = posi[i][2];
  //   if      (posi[i][2] >= boxsize[2]) tmp -= boxsize[2];
  //   else if (posi[i][2] <  0)          tmp += boxsize[2];
  //   unsigned iz = unsigned (tmp / hz);
  //   profile[index3to1(ix, iy, iz)] += 1.;
  // }

  // double dvolume = hx * hy * hz;

  // for (unsigned i = 0; i < profile.size(); ++i){
  //   profile[i] /= dvolume;
  // }
}


DensityProfile_PiecewiseConst::
DensityProfile_PiecewiseConst (const std::vector<std::string> & filename,
			       const double & refh)
{
  std::vector<int > resdindex;
  std::vector<std::string >  resdname;
  std::vector<std::string >  atomname;
  std::vector<int > atomindex;
  std::vector<std::vector<double > > posi;
  std::vector<std::vector<double > > velo;

  GroFileManager::read (filename[0],
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

  std::cout << "# cal profile from file " << filename[0] << std::endl;
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

  for (unsigned jj = 1; jj < filename.size(); ++jj){
    std::cout << "# cal profile from file " << filename[jj] << std::endl;
    GroFileManager::read (filename[jj],
			  resdindex, resdname, atomname, atomindex,
			  posi, velo, boxsize);
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
  }

  double dvolume = hx * hy * hz;

  for (unsigned i = 0; i < profile.size(); ++i){
    profile[i] /= dvolume * filename.size();
  }
}

void DensityProfile_PiecewiseConst::
reinit_xtc (const std::string & filename,
	    const double & refh)
{
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  XDRFILE * fp = xdrfile_open (fname, "r");
  int natoms;
  if (read_xtc_natoms (fname, &natoms) != 1){
    std::cerr << "wrong reading natoms" << std::endl;    
    exit (1);
  }
  std::cout << "# natom is " << natoms << std::endl;
  
  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  if (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) != 0){
    std::cerr << "error while loading xtc file" << std::endl;
    exit(1); 
  }
  std::cout << "loaded frame at time " << time << "ps   \r";  
  std::cout << std::flush;  

  boxsize.resize(3);
  boxsize[0] = gbox[0][0];
  boxsize[1] = gbox[1][1];
  boxsize[2] = gbox[2][2];
  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  
  profile.clear();
  profile.resize (nx * ny * nz, 0.);

  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) != 0){
    for (unsigned i = 0; i < unsigned(natoms); ++i) {
      double tmp;
      tmp = xx[i][0];
      if      (xx[i][0] >= boxsize[0]) tmp -= boxsize[0];
      else if (xx[i][0] <  0)          tmp += boxsize[0];
      unsigned ix = unsigned (tmp / hx);
      tmp = xx[i][1];
      if      (xx[i][1] >= boxsize[1]) tmp -= boxsize[1];
      else if (xx[i][1] <  0)          tmp += boxsize[1];
      unsigned iy = unsigned (tmp / hy);
      tmp = xx[i][2];
      if      (xx[i][2] >= boxsize[2]) tmp -= boxsize[2];
      else if (xx[i][2] <  0)          tmp += boxsize[2];
      unsigned iz = unsigned (tmp / hz);
      profile[index3to1(ix, iy, iz)] += 1.;
    }
  }

  free(xx);
  xdrfile_close(fp);
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
    // double sum = 0.;
    // for (unsigned j = 0; j < ny; ++j){
    //   for (unsigned k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    // fprintf (fp, "%f %f\n", (i + 0.5) * hx, sum / ny / nz);
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, profile[index3to1(i, 0, 0)]);
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

