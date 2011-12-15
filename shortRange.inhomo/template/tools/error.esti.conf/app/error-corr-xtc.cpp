#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>

#include "ErrorEstimator.h"
#include "DensityProfile.h"

int main(int argc, char * argv[])
{
  if (argc != 6){
    std::cerr << "usage:\n"
	      << argv[0] << " xtc refh_Profile refh_Error rcut rcorr" << std::endl;
    return 1;
  }
  double refh_Profile, refh_Error, rcut;
  refh_Profile = atof(argv[2]);
  refh_Error = atof(argv[3]);
  rcut = atof(argv[4]);
  int rcorr = atoi(argv[5]);
  int ratio = refh_Profile / refh_Error;

  DensityProfile_Corr_PiecewiseConst dp (std::string(argv[1]), rcorr, refh_Profile);
  dp.print_x  (std::string("profile.x.out"));
  dp.print_xy (std::string("profile.xy.out"));

  int corrBond = 1;
  int nx = dp.getNx(), ny = dp.getNy(), nz = dp.getNz();
  int jx, jy, jz;
  int kx, ky, kz;
  double hx, hy, hz;
  hx = hy = hz = refh_Profile;
  jx = 0;
  jy = jz = 0;
  for (int ddx = -corrBond; ddx <= corrBond; ++ddx){
    kx = jx + ddx;
    if      (kx <   0) kx += nx;
    else if (kx >= nx) kx -= nx;		
    for (int ddy = -corrBond; ddy <= corrBond; ++ddy){
      ky = jy + ddy;
      if      (ky <   0) ky += ny;
      else if (ky >= ny) ky -= ny;		
      for (int ddz = -corrBond; ddz <= corrBond; ++ddz){
	// if (ddx == 0 && ddy == 0 && ddz == 0) continue;
	kz = jz + ddz;
	if      (kz <   0) kz += nz;
	else if (kz >= nz) kz -= nz;		
	double rho2 =
	    dp.getCorr ((jx+0.5)*hx, (jy+0.5)*hy, (jz+0.5)*hz,
			ddx*hx, ddy*hy, ddz*hz);
	double rho12 =
	    dp.getMean ((kx+0.5)*hx, (ky+0.5)*hy, (kz+0.5)*hz);
	double rho11 =
	    dp.getMean ((jx+0.5)*hx, (jy+0.5)*hy, (jz+0.5)*hz);
	printf ("%d %d %d    rho2: %f, rho11 %f, rho12 %f, corr %f\n",
		ddx, ddy, ddz,
		rho2, rho11, rho12,
		rho2 - rho11 * rho12);
      }
    }
  }		    

  Disperson6 lj (1., 1., rcut);
  ErrorEstimatorConvN2_Corr error (lj, dp.getBox(), refh_Error);

  error.estimate (dp, rcorr * ratio);
  error.print_x (std::string("esti.x.out"));
  error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


