#include "DensityProfile.h"
#include <string.h>

DensityProfile_PiecewiseConst::
DensityProfile_PiecewiseConst (const std::string & filename,
			       const double & refh)
{
  reinit_xtc (filename, refh);
}

void DensityProfile_PiecewiseConst::
reinit_xtc (const std::string & filename,
	    const double & refh)
{
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  XDRFILE * fp = xdrfile_open (fname, "r");
  if (read_xtc_natoms (fname, &natoms) != 0){
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
  // std::cout << "loaded frame at time " << time << "ps   \r";  
  // std::cout << std::flush;  

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

  // int nfile = 0;
  // while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
  //   std::cout << "loaded frame at time " << time << "ps   \r";  
  //   std::cout << std::flush;  
  //   for (unsigned i = 0; i < unsigned(natoms); ++i) {
  //     double tmp;
  //     tmp = xx[i][0];
  //     if      (xx[i][0] >= boxsize[0]) tmp -= boxsize[0];
  //     else if (xx[i][0] <  0)          tmp += boxsize[0];
  //     unsigned ix = unsigned (tmp / hx);
  //     tmp = xx[i][1];
  //     if      (xx[i][1] >= boxsize[1]) tmp -= boxsize[1];
  //     else if (xx[i][1] <  0)          tmp += boxsize[1];
  //     unsigned iy = unsigned (tmp / hy);
  //     tmp = xx[i][2];
  //     if      (xx[i][2] >= boxsize[2]) tmp -= boxsize[2];
  //     else if (xx[i][2] <  0)          tmp += boxsize[2];
  //     unsigned iz = unsigned (tmp / hz);
  //     profile[index3to1(ix, iy, iz)] += 1.;
  //   }
  //   nfile ++;
  // }
  // std::cout << std::endl;

  // double dvolume = hx * hy * hz;
  // for (unsigned i = 0; i < profile.size(); ++i){
  //   profile[i] /= dvolume * nfile;
  // }
  
  free(xx);
  xdrfile_close(fp);
}

void DensityProfile_PiecewiseConst::
deposit (const std::string & filename,
	 const float & start,
	 const float & end)
{
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  XDRFILE * fp = xdrfile_open (fname, "r");
  
  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  float time_prec = 0.001;
  
  int nfile = 0;
  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
    if (end != 0.f) {
      if (time < start - time_prec){
	continue;
      }
      else if (time > end + time_prec) {
	break;
      }	
    }
    else {
      if (time < start - time_prec) continue;
    }
    std::cout << "#! loaded frame at time " << time << "ps   \r";  
    std::cout << std::flush;  
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
    nfile ++;
  }
  std::cout << std::endl;

  double dvolume = hx * hy * hz;
  for (unsigned i = 0; i < profile.size(); ++i){
    profile[i] /= dvolume * nfile;
  }
  
  free(xx);
  xdrfile_close(fp);
}

void DensityProfile_PiecewiseConst::
clear ()
{
  for (unsigned i = 0; i < profile.size(); ++i){
    profile[i] = 0.;
  }  
}

void DensityProfile_PiecewiseConst::
deposit (const rvec * xx,
	 const int natoms)
{
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
    for (unsigned k = 0; k < nz; ++k){
      double vz = (k + 0.5) * hz;
      // double sum = 0.;
      // for (unsigned k = 0; k < nz; ++k){
      // 	sum += profile[index3to1(i, j, k)];
      // }
      // fprintf (fp, "%f %f %f\n", vx, vy, sum / nz);
      fprintf (fp, "%f %f %f\n", vx, vz, profile[index3to1(i,ny/2,k)]);
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
}

