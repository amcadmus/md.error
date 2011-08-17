#include "DensityFunction.h"
#include "GroFileManager.h"
#include "RandomGenerator.h"
#define CPLUSPLUS
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include <stdlib.h>

void DensityFunction::
reinit (const std::vector<double > & x,
	const std::vector<double > & y)
{
  xx = x;
  yy = y;
  mysize = y.size() - 1;
  yy[mysize] = yy[0];
  xlow = xx.front();
  xup  = xx.back();

  coeff.resize (mysize);
  for (unsigned i = 0; i < mysize; ++i){
    coeff[i].a0 = 0.5 * (xx[i] + xx[i+1]);
    coeff[i].a1 = M_PI / (xx[i+1] - xx[i]);
    coeff[i].a2 = (yy[i+1] - yy[i]) * 0.5;
    coeff[i].a3 = (yy[i+1] + yy[i]) * 0.5;
  }
}

double DensityFunction::
integral () const
{
  double sum = 0.;
  for (unsigned i = 0; i < mysize; ++i){
    sum += 0.5 * (yy[i] + yy[i+1]) * (xx[i+1] - xx[i]);
  }
  return sum;
}

double DensityFunction::
yMax() const
{
  double max = 0.;
  for (unsigned i = 0; i < mysize; ++i){
    if (yy[i] > max){
      max = yy[i];
    }
  }
  return max;
}

double DensityFunction::
value (const double & x) const 
{
  int idx;
  if (x < xx.front()){
    idx = 0;
  }
  else if (x >= xx.back()){
    idx = mysize - 1;
  }
  else {
    int ib = mysize;
    int ia = 0;
    while (ib - ia > 1){
      int ic = ((ib - ia) >> 1) + ia;
      if (xx[ic] > x){
	ib = ic;
      }
      else {
	ia = ic;
      }
    }
    idx = ia;
  }  
  return coeff[idx].a2 * sin(coeff[idx].a1 * (x - coeff[idx].a0)) + coeff[idx].a3;
}

void SystemDensityFunction::
reinit (const std::vector<double > & x,
	const std::vector<double > & p,
	const std::vector<double > & n,
	const double & Ly_,
	const double & Lz_)
{
  Ly = Ly_;
  Lz = Lz_;
  
  positive.reinit (x, p);
  negative.reinit (x, n);
  Lx = positive.Lx();

  std::vector<double > t (p);
  for (unsigned i = 0; i < p.size(); ++i){
    t[i] +=  n[i];
  }
  total.reinit (x, t);

  np = (positive.integral() * Ly * Lz);
  nn = (negative.integral() * Ly * Lz);
  natom = np + nn;
}

void SystemDensityFunction::
genConf_gro (const char * filename)
{
  std::vector<int > index (natom, 0);
  std::vector<std::string > name (natom);
  std::vector<std::vector<double > > posi (natom, std::vector<double > (3,0.));
  std::vector<std::vector<double > > velo (natom, std::vector<double > (3,0.));
  std::vector<double > box (3, 0.);
  int count = 1;

  double max = positive.yMax();
  for (unsigned i = 0; i < np; ++i, ++count){
    do{
      posi[i][0] = RandomGenerator_MT19937::genrand_real2() * Lx;
    } while (RandomGenerator_MT19937::genrand_real1() > 
	     positive (posi[i][0]) / max);
    posi[i][1] = RandomGenerator_MT19937::genrand_real2() * Ly;
    posi[i][2] = RandomGenerator_MT19937::genrand_real2() * Lz;
    index[i] = count;
    name[i] = std::string ("ljp");
  }
  max = negative.yMax();
  for (unsigned i = 0; i < nn; ++i, ++count){
    do{
      posi[i+np][0] = RandomGenerator_MT19937::genrand_real2() * Lx;
    } while (RandomGenerator_MT19937::genrand_real1() > 
	     negative (posi[i+np][0]) / max);
    posi[i+np][1] = RandomGenerator_MT19937::genrand_real2() * Ly;
    posi[i+np][2] = RandomGenerator_MT19937::genrand_real2() * Lz;
    index[i+np] = count;
    name[i+np] = std::string ("ljn");
  }

  box[0] = Lx;
  box[1] = Ly;
  box[2] = Lz;

  GroFileManager::write (filename, index, name, name, index,
			 posi, velo, box);
}

void SystemDensityFunction::
genXtc (const char * filename,
	const int & nframe)
{
  XDRFILE * xd;
  xd = xdrfile_open (filename, "w");
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    exit(1);
  }

  matrix box;
  for (unsigned i = 0; i < 3; ++i){
    for (unsigned j = 0; j < 3; ++j){
      box[i][j] = 0;
    }
  }
  box[0][0] = Lx;
  box[1][1] = Ly;
  box[2][2] = Lz;
  
  rvec * posi = (rvec *) malloc (sizeof(rvec) * natom);
  
  for (int time = 0; time < nframe; ++time){
    double max = positive.yMax();
    for (unsigned i = 0; i < np; ++i){
      do{
	posi[i][0] = RandomGenerator_MT19937::genrand_real2() * Lx;
      } while (RandomGenerator_MT19937::genrand_real1() > 
	       positive (posi[i][0]) / max);
      posi[i][1] = RandomGenerator_MT19937::genrand_real2() * Ly;
      posi[i][2] = RandomGenerator_MT19937::genrand_real2() * Lz;
    }
    max = negative.yMax();
    for (unsigned i = 0; i < nn; ++i){
      do{
	posi[i+np][0] = RandomGenerator_MT19937::genrand_real2() * Lx;
      } while (RandomGenerator_MT19937::genrand_real1() > 
	       negative (posi[i+np][0]) / max);
      posi[i+np][1] = RandomGenerator_MT19937::genrand_real2() * Ly;
      posi[i+np][2] = RandomGenerator_MT19937::genrand_real2() * Lz;
    }

    write_xtc (xd, natom, time, time, box, posi, 1000);
  }
  
  free (posi);
  xdrfile_close (xd);
}





