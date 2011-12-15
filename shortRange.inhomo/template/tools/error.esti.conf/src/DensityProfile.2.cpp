#include "DensityProfile.h"
#include "GroFileManager.h"
#include <string.h>

DensityProfile_Corr_PiecewiseConst::
DensityProfile_Corr_PiecewiseConst (const std::string & filename,
				    const int & corrBond_,
				    const double & refh)
    : zero (0.)
{
  reinit_xtc (filename, corrBond_, refh);  
}

void DensityProfile_Corr_PiecewiseConst::
reinit_xtc (const std::string & filename,
	    const int & corrBond_,
	    const double & refh)
{
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  XDRFILE * fp = xdrfile_open (fname, "r");
  int natoms;
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
  double dvolume = hx * hy * hz;
  corrBond = corrBond_;
  corrDim = 2 * corrBond + 1;
  numCorr = corrDim * corrDim * corrDim;
  
  std::vector<double > profile;
  profile.clear();
  profile.resize (nx * ny * nz, 0.);
  mean.clear();
  mean.resize (nx * ny * nz, 0.);
  corr.clear();
  corr.resize (nx * ny * nz, std::vector<double >(numCorr, 0.));

  int nfile = 0;
  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
    std::cout << "loaded frame at time " << time << "ps   \r";  
    std::cout << std::flush;  
    for (int i = 0; i < nx * ny * nz; ++i){
      profile[i] = 0.;
    }
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
      profile[index3to1(ix, iy, iz)] += 1./dvolume;
    }
    nfile ++;
    for (int i = 0; i < nx * ny * nz; ++i){
      mean[i] += profile[i];
    }
    for (int i = 0; i < nx * ny * nz; ++i){
      int myi, myj, myk;
      index1to3 (i, myi, myj, myk);
      int writeIndex = 0;
      for (int di = -corrBond; di <= corrBond; ++di){
	int targeti = myi + di;
	if      (targeti <   0) targeti += nx;
	else if (targeti >= nx) targeti -= nx;
	for (int dj = -corrBond; dj <= corrBond; ++dj){
	  int targetj = myj + dj;
	  if      (targetj <   0) targetj += ny;
	  else if (targetj >= ny) targetj -= ny;
	  for (int dk = -corrBond; dk <= corrBond; ++dk){
	    int targetk = myk + dk;
	    if      (targetk <   0) targetk += nz;
	    else if (targetk >= nz) targetk -= nz;
	    // if (writeIndex != corrIndex3to1 (di, dj, dk)){
	    //   std::cerr << "wrong corrIndex3to1" << std::endl;
	    //   exit (1);
	    // }
	    // int tmpa(-1), tmpb(0), tmpc(1), tmpt;
	    // tmpt = corrIndex3to1 (tmpa, tmpb, tmpc);
	    // corrIndex1to3 (tmpt, tmpa, tmpb, tmpc);
	    corr[i][writeIndex++] +=
		profile[i] * profile[index3to1(targeti, targetj, targetk)];
	  }
	}
      }
    }	
  }
  std::cout << std::endl;

  for (unsigned i = 0; i < mean.size(); ++i){
    mean[i] /= nfile;
  }
  for (unsigned i = 0; i < corr.size(); ++i){
    for (int j = 0; j < numCorr; ++j){
      corr[i][j] /= nfile;
    }
  }
  
  free(xx);
  xdrfile_close(fp);
}


void DensityProfile_Corr_PiecewiseConst::    
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
    // fprintf (fp, "%f %f\n", (i + 0.5) * hx, sum / ny / nz);
    fprintf (fp, "%f %f\n", (i + 0.5) * hx, mean[index3to1(i, 0, 0)]);
  }

  fclose (fp);
}

void DensityProfile_Corr_PiecewiseConst::    
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
      // fprintf (fp, "%f %f %f\n", vx, vy, sum / nz);
      fprintf (fp, "%f %f %f\n", vx, vy, mean[index3to1(i,j,0)]);
    }
    fprintf (fp, "\n");
  }

  fclose (fp);
}

