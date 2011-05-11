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
  
  Disperson6 lj (1., 1., rcut);
  ErrorEstimatorFFT_Corr error (lj, dp);

  double k = 0.8660254037844386;
  double tmp1, tmp2;
  tmp1 = 0;
  tmp2 = 0;
  // double tmp1 = error.integral_an1 (13, k, rcut);
  // double tmp2 = error.integral_an (13, k, rcut);
  double tmp3 = error.integral_an13_numerical (k, rcut);
  printf ("tmp1 (small) is %e, tmp2 (large) is %e, tmp3 (numerical) is %e\n",
	  tmp1, tmp2, tmp3);
  
  error.estimate (dp, rcorr * ratio, rcut);
  error.print_x (std::string("esti.x.out"));
  error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


