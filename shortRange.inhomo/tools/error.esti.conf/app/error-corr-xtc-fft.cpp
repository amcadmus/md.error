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

  error.estimate (dp, rcorr * ratio);
  error.print_x (std::string("esti.x.out"));
  error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


