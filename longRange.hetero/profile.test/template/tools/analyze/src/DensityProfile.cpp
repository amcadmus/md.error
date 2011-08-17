#include "DensityProfile.h"
#include "GroFileManager.h"
#include <fstream>
#include <iostream>


void DensityProfile_PiecewiseConst::
reinit_conf (const std::string & filename,
	     const int & nx_,
	     const int & ny_,
	     const int & nz_)
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
  nx = nx_;
  ny = ny_;
  nz = nz_;
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;
  
  profile1.clear();
  profile1.resize (nx * ny * nz, 0.);
  profile2.clear();
  profile2.resize (nx * ny * nz, 0.);
  
  {
    unsigned i = 0;
    for (; i < posi.size() / 2; ++i) {
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
      profile1[index3to1(ix, iy, iz)] += 1.;
      profile2[index3to1(ix, iy, iz)] += 1.;
    }
    for (; i < posi.size(); ++i){
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
      profile1[index3to1(ix, iy, iz)] -= 1.;
      profile2[index3to1(ix, iy, iz)] += 1.;
    }
  }

  double dvolume = hx * hy * hz;

  for (unsigned i = 0; i < profile1.size(); ++i){
    profile1[i] /= dvolume;
    profile2[i] /= dvolume;
  }
}


// void DensityProfile_PiecewiseConst::
// reinit_conf (const std::string & filename,
// 	     const double & refh)
// {
//   std::vector<int > resdindex;
//   std::vector<std::string >  resdname;
//   std::vector<std::string >  atomname;
//   std::vector<int > atomindex;
//   std::vector<std::vector<double > > posi;
//   std::vector<std::vector<double > > velo;

//   GroFileManager::read (filename,
//   			resdindex, resdname, atomname, atomindex,
//   			posi, velo, boxsize);
//   nx = unsigned (boxsize[0] / refh);
//   ny = unsigned (boxsize[1] / refh);
//   nz = unsigned (boxsize[2] / refh);
//   if ((nx - (nx/2)*2) == 0) nx++;
//   if ((ny - (ny/2)*2) == 0) ny++;
//   if ((nz - (nz/2)*2) == 0) nz++;
//   hx = boxsize[0] / nx;
//   hy = boxsize[1] / ny;
//   hz = boxsize[2] / nz;

//   profile.clear();
//   profile.resize (nx * ny * nz, 0.);

//   for (unsigned i = 0; i < posi.size(); ++i) {
//     double tmp;
//     tmp = posi[i][0];
//     if      (posi[i][0] >= boxsize[0]) tmp -= boxsize[0];
//     else if (posi[i][0] <  0)          tmp += boxsize[0];
//     unsigned ix = unsigned (tmp / hx);
//     tmp = posi[i][1];
//     if      (posi[i][1] >= boxsize[1]) tmp -= boxsize[1];
//     else if (posi[i][1] <  0)          tmp += boxsize[1];
//     unsigned iy = unsigned (tmp / hy);
//     tmp = posi[i][2];
//     if      (posi[i][2] >= boxsize[2]) tmp -= boxsize[2];
//     else if (posi[i][2] <  0)          tmp += boxsize[2];
//     unsigned iz = unsigned (tmp / hz);
//     profile[index3to1(ix, iy, iz)] += 1.;
//   }

//   double dvolume = hx * hy * hz;

//   for (unsigned i = 0; i < profile.size(); ++i){
//     profile[i] /= dvolume;
//   }
// }

// void DensityProfile_PiecewiseConst::
// reinit_xtc (const std::string & filename,
// 	    const double & refh)
// {
//   char fname [1024];
//   strncpy (fname, filename.c_str(), 1024);
//   XDRFILE * fp = xdrfile_open (fname, "r");
//   if (read_xtc_natoms (fname, &natoms) != 0){
//     std::cerr << "wrong reading natoms" << std::endl;    
//     exit (1);
//   }
//   std::cout << "# natom is " << natoms << std::endl;
  
//   int step;
//   float time, prec;
//   matrix gbox;  
//   rvec * xx;
//   xx = (rvec *) malloc (sizeof(rvec) * natoms);
//   if (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) != 0){
//     std::cerr << "error while loading xtc file" << std::endl;
//     exit(1); 
//   }
//   // std::cout << "loaded frame at time " << time << "ps   \r";  
//   // std::cout << std::flush;  

//   boxsize.resize(3);
//   boxsize[0] = gbox[0][0];
//   boxsize[1] = gbox[1][1];
//   boxsize[2] = gbox[2][2];
//   nx = unsigned (boxsize[0] / refh);
//   ny = unsigned (boxsize[1] / refh);
//   nz = unsigned (boxsize[2] / refh);
//   // nx, ny, nz should be odd
//   if ((nx - (nx/2)*2) == 0) nx++;
//   if ((ny - (ny/2)*2) == 0) ny++;
//   if ((nz - (nz/2)*2) == 0) nz++;  
//   hx = boxsize[0] / nx;
//   hy = boxsize[1] / ny;
//   hz = boxsize[2] / nz;
  
//   profile.clear();
//   profile.resize (nx * ny * nz, 0.);
  
//   free(xx);
//   xdrfile_close(fp);
// }

// void DensityProfile_PiecewiseConst::
// deposit (const std::string & filename,
// 	 const float & start,
// 	 const float & end)
// {
//   char fname [1024];
//   strncpy (fname, filename.c_str(), 1024);
//   XDRFILE * fp = xdrfile_open (fname, "r");
  
//   int step;
//   float time, prec;
//   matrix gbox;  
//   rvec * xx;
//   xx = (rvec *) malloc (sizeof(rvec) * natoms);
//   float time_prec = 0.001;
  
//   int nfile = 0;
//   while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
//     if (end != 0.f) {
//       if (time < start - time_prec){
// 	continue;
//       }
//       else if (time > end + time_prec) {
// 	break;
//       }	
//     }
//     else {
//       if (time < start - time_prec) continue;
//     }
//     std::cout << "#! loaded frame at time " << time << "ps   \r";  
//     std::cout << std::flush;  
//     for (unsigned i = 0; i < unsigned(natoms); ++i) {
//       double tmp;
//       tmp = xx[i][0];
//       if      (xx[i][0] >= boxsize[0]) tmp -= boxsize[0];
//       else if (xx[i][0] <  0)          tmp += boxsize[0];
//       unsigned ix = unsigned (tmp / hx);
//       tmp = xx[i][1];
//       if      (xx[i][1] >= boxsize[1]) tmp -= boxsize[1];
//       else if (xx[i][1] <  0)          tmp += boxsize[1];
//       unsigned iy = unsigned (tmp / hy);
//       tmp = xx[i][2];
//       if      (xx[i][2] >= boxsize[2]) tmp -= boxsize[2];
//       else if (xx[i][2] <  0)          tmp += boxsize[2];
//       unsigned iz = unsigned (tmp / hz);
//       profile[index3to1(ix, iy, iz)] += 1.;
//     }
//     nfile ++;
//   }
//   std::cout << std::endl;

//   double dvolume = hx * hy * hz;
//   for (unsigned i = 0; i < profile.size(); ++i){
//     profile[i] /= dvolume * nfile;
//   }
  
//   free(xx);
//   xdrfile_close(fp);
// }


void DensityProfile_PiecewiseConst::
reinit_xtc (const char * fname,
	    const int & nx_,
	    const int & ny_,
	    const int & nz_,
	    const float & start,
	    const float & end)
{
  XDRFILE * fp = xdrfile_open (fname, "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << fname << std::endl;
    exit(1);
  }
  char tfname [1024];
  strncpy (tfname, fname, 1024);
  if (read_xtc_natoms (tfname, &natoms) != 0){
    std::cerr << "wrong reading natoms" << std::endl;    
    exit (1);
  }
  std::cout << "# natom is " << natoms << std::endl;

  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  float time_prec = 0.001;
  if (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) != 0){
    std::cerr << "error while loading xtc file" << std::endl;
    exit(1); 
  }

  boxsize.resize(3);
  boxsize[0] = gbox[0][0];
  boxsize[1] = gbox[1][1];
  boxsize[2] = gbox[2][2];
  nx = nx_;
  ny = ny_;
  nz = nz_;
  // nx, ny, nz should be odd
  if ((nx - (nx/2)*2) == 0) nx++;
  if ((ny - (ny/2)*2) == 0) ny++;
  if ((nz - (nz/2)*2) == 0) nz++;  
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;
  profile1.clear();
  profile1.resize (nx * ny * nz, 0.);
  profile2.clear();
  profile2.resize (nx * ny * nz, 0.);
  
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
    {
      int i = 0;
      for (; i < natoms / 2; ++i) {
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
	profile1[index3to1(ix, iy, iz)] += 1.;
	profile2[index3to1(ix, iy, iz)] += 1.;
      }
      for (; i < natoms; ++i){
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
	profile1[index3to1(ix, iy, iz)] -= 1.;
	profile2[index3to1(ix, iy, iz)] += 1.;
      }
    }
    nfile ++;
  }
  std::cout << std::endl;

  double dvolume = hx * hy * hz;
  for (unsigned i = 0; i < profile1.size(); ++i){
    profile1[i] /= dvolume * nfile;
    profile2[i] /= dvolume * nfile;
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
    fprintf (fp, "%f %f %f\n",
	     (i + 0.5) * hx,
	     profile1[index3to1(i, 0, 0)],
	     profile2[index3to1(i, 0, 0)]
	);
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
      fprintf (fp, "%f %f %f %f\n",
	       vx, vy,
	       profile1[index3to1(i,j,0)],
	       profile2[index3to1(i,j,0)]
	  );
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
}

void DensityProfile_PiecewiseConst::    
print_avg_x (const std::string & file) const 
{
  FILE * fp = fopen (file.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }
  for (unsigned i = 0; i < nx; ++i){
    double sum1 = 0., sum2 = 0.;
    for (unsigned j = 0; j < ny; ++j){
      for (unsigned k = 0; k < nz; ++k){
    	sum1 += profile1[index3to1(i, j, k)];
    	sum2 += profile2[index3to1(i, j, k)];
      }
    }
    fprintf (fp, "%f %f %f\n",
	     (i + 0.5) * hx,
	     sum1 / ny / nz,
	     sum2 / ny / nz
	);
  }
  fclose (fp);
}


void DensityProfile_PiecewiseConst::
save (const char * file) const
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit (1);
  }

  fprintf (fp, "%f %f %f  ", boxsize[0], boxsize[1], boxsize[2]);
  fprintf (fp, "%d %d %d\n", nx, ny, nz);
  for (unsigned i = 0; i < nele; ++i){
    fprintf (fp, "%e %e", profile1[i], profile2[i]);
  }

  fclose(fp);
}

void DensityProfile_PiecewiseConst::
load (const char * file)
{
  FILE * fp = fopen (file, "r");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit (1);
  }

  boxsize.resize(3);
  fscanf (fp, "%lf %lf %lf  ", &boxsize[0], &boxsize[1], &boxsize[2]);
  fscanf (fp, "%d %d %d\n", &nx, &ny, &nz);
  nele = nx * ny * nz;
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  profile1.resize(nele);
  profile2.resize(nele);
  
  for (unsigned i = 0; i < nele; ++i){
    fscanf (fp, "%lf %lf", &profile1[i], &profile2[i]);
  }

  fclose(fp);
}

