#include "DensityProfile.h"
#include <fstream>
#include <iostream>

// DensityProfile_PiecewiseConst::
// DensityProfile_PiecewiseConst (const std::string & filename,
// 			       const double & refh)
// {
//   reinit_xtc (filename, refh);
// }

void DensityProfile_PiecewiseConst::
reinit (const ScalorType & bx,
	const ScalorType & by,
	const ScalorType & bz,
	const double & refh)
{
  boxsize.resize(3);
  boxsize[0] = bx;
  boxsize[1] = by;
  boxsize[2] = bz;
  nx = unsigned (boxsize[0] / refh);
  ny = unsigned (boxsize[1] / refh);
  nz = unsigned (boxsize[2] / refh);
  if ((nx - (nx/2)*2) == 0) nx++;
  if ((ny - (ny/2)*2) == 0) ny++;
  if ((nz - (nz/2)*2) == 0) nz++;
  nele = nx * ny * nz;
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  
  profile.clear();
  profile.resize (nx * ny * nz, 0.);

  clearData ();
}

void DensityProfile_PiecewiseConst::
clearData ()
{
  for (int i = 0; i < profile.size(); ++i){
    profile[i] = 0.;
  }
  ndata = 0;
}

void DensityProfile_PiecewiseConst::
deposite (const CoordType * coord,
	  const IndexType numAtom)
{
  for (IndexType i = 0; i < numAtom; ++i){
    double tmp;
    tmp = coord[i].x;
    if      (coord[i].x >= boxsize[0]) tmp -= boxsize[0];
    else if (coord[i].x <  0)          tmp += boxsize[0];
    unsigned ix = unsigned (tmp / hx);
    tmp = coord[i].y;
    if      (coord[i].y >= boxsize[1]) tmp -= boxsize[1];
    else if (coord[i].y <  0)          tmp += boxsize[1];
    unsigned iy = unsigned (tmp / hy);
    tmp = coord[i].z;
    if      (coord[i].z >= boxsize[2]) tmp -= boxsize[2];
    else if (coord[i].z <  0)          tmp += boxsize[2];
    unsigned iz = unsigned (tmp / hz);
    profile[index3to1(ix, iy, iz)] += 1.;
  }
  ndata ++;
}

void DensityProfile_PiecewiseConst::
calculate ()
{
  double dvolume = hx * hy * hz;
  for (unsigned i = 0; i < profile.size(); ++i){
    profile[i] /= dvolume * ndata;
  }
}


// void DensityProfile_PiecewiseConst::
// reinit_xtc (const std::string & filename,
// 	    const double & refh)
// {
//   char fname [1024];
//   strncpy (fname, filename.c_str(), 1024);
//   XDRFILE * fp = xdrfile_open (fname, "r");
//   int natoms;
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
//   std::cout << "loaded frame at time " << time << "ps   \r";  
//   std::cout << std::flush;  

//   boxsize.resize(3);
//   boxsize[0] = gbox[0][0];
//   boxsize[1] = gbox[1][1];
//   boxsize[2] = gbox[2][2];
//   nx = unsigned (boxsize[0] / refh);
//   ny = unsigned (boxsize[1] / refh);
//   nz = unsigned (boxsize[2] / refh);
//   hx = boxsize[0] / nx;
//   hy = boxsize[1] / ny;
//   hz = boxsize[2] / nz;
  
//   profile.clear();
//   profile.resize (nx * ny * nz, 0.);

//   int nfile = 0;
//   while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
//     std::cout << "loaded frame at time " << time << "ps   \r";  
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
  for (int i = 0; i < nele; ++i){
    fprintf (fp, "%e ", profile[i]);
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
  profile.resize(nele);
  ndata = 0;
  
  for (int i = 0; i < nele; ++i){
    fscanf (fp, "%lf", &profile[i]);
  }

  fclose(fp);
}

